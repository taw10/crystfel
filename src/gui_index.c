/*
 * gui_index.c
 *
 * Peak search parts of GUI
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms-compat.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <datatemplate.h>
#include <peaks.h>
#include <cell.h>
#include <cell-utils.h>
#include <integration.h>
#include <predict-refine.h>

#include "crystfelimageview.h"
#include "crystfelindexingopts.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gui_peaksearch.h"
#include "gtk-util-routines.h"


void cell_explorer_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GSubprocess *sp;
	GError *error = NULL;
	const gchar *results_name;
	struct gui_indexing_result *res;
	const gchar **streams;
	int i;

	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	if ( strcmp(results_name, "crystfel-gui-internal") == 0 ) {
		STATUS("Please select results first.\n");
		return;
	}

	res = find_indexing_result_by_name(proj, results_name);
	if ( res == NULL ) {
		ERROR("Results for '%s' not found!\n", results_name);
		return;
	}

	streams = malloc((res->n_streams+2)*sizeof(gchar *));
	if ( streams == NULL ) return;

	streams[0] = get_crystfel_exe("cell_explorer");
	for ( i=0; i<res->n_streams; i++ ) {
		streams[i+1] = res->streams[i];
	}
	streams[res->n_streams+1] = NULL;

	sp = g_subprocess_newv(streams, G_SUBPROCESS_FLAGS_NONE,
	                       &error);
	free(streams);
	if ( sp == NULL ) {
		ERROR("Failed to start cell_explorer: %s\n",
		      error->message);
		g_error_free(error);
		return;
	}

	STATUS("Starting cell_explorer...\n");
}


static void get_indexing_opts(struct crystfelproject *proj,
                              CrystFELIndexingOpts *opts)
{
	/* Indexing */
	proj->indexing_params.cell_file = crystfel_indexing_opts_get_cell_file(opts);
	proj->indexing_params.indexing_methods = crystfel_indexing_opts_get_indexing_method_string(opts);
	proj->indexing_params.multi = crystfel_indexing_opts_get_multi_lattice(opts);
	proj->indexing_params.no_refine = !crystfel_indexing_opts_get_refine(opts);
	proj->indexing_params.no_retry = !crystfel_indexing_opts_get_retry(opts);
	proj->indexing_params.no_peak_check = !crystfel_indexing_opts_get_peak_check(opts);
	proj->indexing_params.no_cell_check = !crystfel_indexing_opts_get_cell_check(opts);
	proj->indexing_params.min_peaks = crystfel_indexing_opts_get_min_peaks(opts);

	/* Integration */
	proj->indexing_params.integration_method = crystfel_indexing_opts_get_integration_method_string(opts);
	proj->indexing_params.overpredict = crystfel_indexing_opts_get_overpredict(opts);
	proj->indexing_params.push_res = crystfel_indexing_opts_get_push_res(opts);
	crystfel_indexing_opts_get_integration_radii(opts,
	                                             &proj->indexing_params.ir_inn,
	                                             &proj->indexing_params.ir_mid,
	                                             &proj->indexing_params.ir_out);
	proj->indexing_params.fix_profile_radius = crystfel_indexing_opts_get_fixed_profile_radius(opts,
	                                      &proj->indexing_params.use_fix_profile_radius);
	proj->indexing_params.fix_divergence = crystfel_indexing_opts_get_fixed_divergence(opts);

	/* Stream output */
	proj->indexing_params.exclude_nonhits = crystfel_indexing_opts_get_exclude_blanks(opts);
	proj->indexing_params.exclude_peaks = crystfel_indexing_opts_get_exclude_peaks(opts);
	proj->indexing_params.exclude_refls = crystfel_indexing_opts_get_exclude_reflections(opts);
	proj->indexing_params.metadata_to_copy = crystfel_indexing_opts_get_metadata_to_copy(opts,
		               &proj->indexing_params.n_metadata);
}


static int run_indexing_all(struct crystfelproject *proj,
                            int backend_idx, const char *job_title,
                            const char *job_notes)
{
	struct crystfel_backend *be;
	void *job_priv;

	be = &proj->backends[backend_idx];
	job_priv = be->run_indexing(job_title, job_notes, proj,
	                            be->indexing_opts_priv);

	if ( job_priv != NULL ) {
		char name[256];
		snprintf(name, 255, "Indexing all frames (%s)",
		         job_title);
		add_running_task(proj, name, be, job_priv);
		return 0;
	} else {
		return 1;
	}
}


struct new_index_job_params {
	struct crystfelproject *proj;
	struct gui_job_notes_page *notes_page;
	GtkWidget *indexing_backend_combo;
	GtkWidget *indexing_backend_opts_widget;
	GtkWidget *indexing_backend_opts_box;
	GtkWidget *job_title_entry;
};


static void index_all_response_sig(GtkWidget *dialog, gint resp,
                                   struct new_index_job_params *njp)
{
	if ( resp == GTK_RESPONSE_OK ) {

		int backend_idx;
		const char *job_title;
		char *job_notes;

		get_indexing_opts(njp->proj,
		                  CRYSTFEL_INDEXING_OPTS(njp->proj->indexing_opts));

		backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(njp->indexing_backend_combo));
		if ( backend_idx < 0 ) return;

		job_title = gtk_entry_get_text(GTK_ENTRY(njp->job_title_entry));
		job_notes = get_all_text(GTK_TEXT_VIEW(njp->notes_page->textview));

		if ( job_title[0] == '\0' ) {
			ERROR("You must provide a job name.\n");
			return;
		}

		free(njp->proj->indexing_new_job_title);
		njp->proj->indexing_new_job_title = strdup(job_title);

		if ( run_indexing_all(njp->proj, backend_idx,
		                      job_title, job_notes) == 0 )
		{
			gtk_widget_destroy(dialog);
			njp->proj->indexing_opts = NULL;
		}

		free(job_notes);

	} else {
		gtk_widget_destroy(dialog);
		njp->proj->indexing_opts = NULL;
	}
}


static void indexing_backend_changed_sig(GtkWidget *combo,
                                         struct new_index_job_params *njp)
{
	int backend_idx;
	struct crystfel_backend *be;

	backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	if ( backend_idx < 0 ) return;
	njp->proj->indexing_backend_selected = backend_idx;

	be = &njp->proj->backends[backend_idx];

	if ( njp->indexing_backend_opts_widget != NULL ) {
		gtk_widget_destroy(njp->indexing_backend_opts_widget);
	}

	njp->indexing_backend_opts_widget = be->make_indexing_parameters_widget(be->indexing_opts_priv);

	gtk_box_pack_start(GTK_BOX(njp->indexing_backend_opts_box),
	                   GTK_WIDGET(njp->indexing_backend_opts_widget),
	                   FALSE, FALSE, 0);
	gtk_widget_show_all(njp->indexing_backend_opts_widget);
}


static GtkWidget *make_backend_opts(struct new_index_job_params *njp)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	int i;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Batch system:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	njp->indexing_backend_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->indexing_backend_combo),
	                   FALSE, FALSE, 0);

	for ( i=0; i<njp->proj->n_backends; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(njp->indexing_backend_combo),
		                          njp->proj->backends[i].name,
		                          njp->proj->backends[i].friendly_name);
	}

	njp->indexing_backend_opts_box = gtk_box_new(GTK_ORIENTATION_VERTICAL,
	                                              0);
	gtk_box_pack_start(GTK_BOX(box),
	                   GTK_WIDGET(njp->indexing_backend_opts_box),
	                   FALSE, FALSE, 0);
	njp->indexing_backend_opts_widget = NULL;

	/* njp->indexing_backend_opts{_box} must exist before the following */
	g_signal_connect(G_OBJECT(njp->indexing_backend_combo), "changed",
	                 G_CALLBACK(indexing_backend_changed_sig), njp);
	gtk_combo_box_set_active(GTK_COMBO_BOX(njp->indexing_backend_combo),
	                         njp->proj->indexing_backend_selected);

	return box;
}


static void set_indexing_opts(struct crystfelproject *proj,
                              CrystFELIndexingOpts *opts)
{
	/* Indexing */
	crystfel_indexing_opts_set_cell_file(opts, proj->indexing_params.cell_file);
	crystfel_indexing_opts_set_indexing_method_string(opts, proj->indexing_params.indexing_methods);
	crystfel_indexing_opts_set_multi_lattice(opts, proj->indexing_params.multi);
	crystfel_indexing_opts_set_refine(opts, !proj->indexing_params.no_refine);
	crystfel_indexing_opts_set_retry(opts, !proj->indexing_params.no_retry);
	crystfel_indexing_opts_set_peak_check(opts, !proj->indexing_params.no_peak_check);
	crystfel_indexing_opts_set_cell_check(opts, !proj->indexing_params.no_cell_check);
	crystfel_indexing_opts_set_tolerances(opts, proj->indexing_params.tols);
	crystfel_indexing_opts_set_min_peaks(opts, proj->indexing_params.min_peaks);

	/* Integration */
	crystfel_indexing_opts_set_integration_method_string(opts, proj->indexing_params.integration_method);
	crystfel_indexing_opts_set_overpredict(opts, proj->indexing_params.overpredict);
	crystfel_indexing_opts_set_push_res(opts, proj->indexing_params.push_res);
	crystfel_indexing_opts_set_integration_radii(opts,
	                                             proj->indexing_params.ir_inn,
	                                             proj->indexing_params.ir_mid,
	                                             proj->indexing_params.ir_out);
	crystfel_indexing_opts_set_fixed_profile_radius(opts,
	                                                proj->indexing_params.use_fix_profile_radius,
	                                                proj->indexing_params.fix_profile_radius);
	crystfel_indexing_opts_set_fixed_divergence(opts,
	                                            proj->indexing_params.fix_divergence);

	/* Stream output */
	crystfel_indexing_opts_set_exclude_blanks(opts,
	                                          proj->indexing_params.exclude_nonhits);
	crystfel_indexing_opts_set_exclude_peaks(opts,
	                                         proj->indexing_params.exclude_peaks);
	crystfel_indexing_opts_set_exclude_reflections(opts,
	                                               proj->indexing_params.exclude_refls);
	crystfel_indexing_opts_set_metadata_to_copy(opts,
	                                            proj->indexing_params.metadata_to_copy,
	                                            proj->indexing_params.n_metadata);
}


static void free_new_index_job_params(gpointer njp, GClosure *closure)
{
	free(njp);
}


gint index_all_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *backend_page;
	char *new_title;
	struct new_index_job_params *njp;

	if ( proj->indexing_opts != NULL ) return FALSE;

	njp = malloc(sizeof(struct new_index_job_params));
	if ( njp == NULL ) return FALSE;

	njp->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Index all frames",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect_data(G_OBJECT(dialog), "response",
	                      G_CALLBACK(index_all_response_sig),
	                      njp, free_new_index_job_params, 0);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Job/output name:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	njp->job_title_entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(njp->job_title_entry), 16);
	gtk_entry_set_placeholder_text(GTK_ENTRY(njp->job_title_entry),
	                               "indexing-trial-1");
	new_title = make_new_job_title(proj->indexing_new_job_title);
	if ( new_title != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(njp->job_title_entry), new_title);
		free(new_title);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->job_title_entry),
	                   TRUE, TRUE, 4.0);

	proj->indexing_opts = crystfel_indexing_opts_new();
	crystfel_indexing_opts_set_show_stream_opts(CRYSTFEL_INDEXING_OPTS(proj->indexing_opts),
	                                            TRUE);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(proj->indexing_opts),
	                   FALSE, FALSE, 8.0);
	set_indexing_opts(proj,
	                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));

	backend_page = make_backend_opts(njp);
	gtk_notebook_append_page(GTK_NOTEBOOK(proj->indexing_opts),
	                          backend_page,
	                          gtk_label_new("Cluster/batch system"));

	njp->notes_page = add_job_notes_page(proj->indexing_opts);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}


static IndexingFlags indexing_flags(struct index_params *params)
{
	IndexingFlags fl = 0;

	if ( !params->no_retry ) fl |= INDEXING_RETRY;
	if ( params->multi ) fl |= INDEXING_MULTI;
	if ( !params->no_refine ) fl |= INDEXING_REFINE;
	if ( !params->no_peak_check ) fl |= INDEXING_CHECK_PEAKS;
	if ( !params->no_cell_check ) fl |= INDEXING_CHECK_CELL;

	return fl;
}


static char *enter_gui_tempdir()
{
	char *tmpdir;
	struct stat s;

	tmpdir = malloc(64);
	if ( tmpdir == NULL ) {
		ERROR("Failed to allocate temporary directory name\n");
		return NULL;
	}
	snprintf(tmpdir, 63, "crystfel-gui.%i", getpid());

	if ( stat(tmpdir, &s) == -1 ) {

		int r;

		if ( errno != ENOENT ) {
			ERROR("Failed to stat temporary folder.\n");
			return NULL;
		}

		r = mkdir(tmpdir, S_IRWXU);
		if ( r ) {
			ERROR("Failed to create temporary folder: %s\n",
			      strerror(errno));
			return NULL;
		}

	}

	return tmpdir;
}


static void delete_gui_tempdir(char *tmpdir)
{
	char *path;
	int i;

	/* List of files which it's safe to delete */
	char *files[] = {"gmon.out", "mosflm.lp", "SUMMARY", "XDS.INP",
	                 "xfel_001.img", "xfel_001.spt", "xfel.drx",
	                 "xfel.felix", "xfel.gve", "xfel.ini", "xfel.log",
	                 "IDXREF.LP", "SPOT.XDS", "xfel.newmat", "XPARM.XDS"};

	/* Number of items in the above list */
	int n_files = 15;

	if ( tmpdir == NULL ) return;

	path = calloc(strlen(tmpdir)+64, 1);
	if ( path == NULL ) return;

	for ( i=0; i<n_files; i++ ) {
		snprintf(path, 127, "%s/%s", tmpdir, files[i]);
		unlink(path);
	}

	if ( rmdir(tmpdir) ) {
		ERROR("Failed to delete GUI temporary folder: %s\n",
		      strerror(errno));
	}

	free(tmpdir);
}


static void run_indexing_once(struct crystfelproject *proj)
{
	IndexingPrivate *ipriv;
	UnitCell *cell;
	IntegrationMethod int_method;
	char *methods;
	int i;
	int err;
	TakeTwoOptions *taketwoopts;
	XGandalfOptions *xgandalf_opts;
	PinkIndexerOptions *pinkIndexer_opts;
	FelixOptions *felix_opts;
	FromFileOptions *fromfile_opts;
	char *old_cwd;
	char *tmpdir;
	int r;

	if ( proj->indexing_params.cell_file != NULL ) {
		cell = load_cell_from_file(proj->indexing_params.cell_file);
	} else {
		cell = NULL;
	}

	if ( proj->cur_image->features == NULL ) {
		update_peaks(proj);
	}

	old_cwd = getcwd(NULL, 0);
	tmpdir = enter_gui_tempdir();

	if ( proj->indexing_params.indexing_methods == NULL ) {
		methods = detect_indexing_methods(cell);
		STATUS("Auto-detected indexng methods: %s\n",
		       methods);
	} else {
		methods = strdup(proj->indexing_params.indexing_methods);
	}

	/* Get default options for the indexing methods.
	 * The GUI currently does not allow them to be changed */
	default_method_options(&taketwoopts,
	                       &xgandalf_opts,
	                       &pinkIndexer_opts,
	                       &felix_opts,
	                       &fromfile_opts);

	ipriv = setup_indexing(methods, cell,
	                       proj->indexing_params.tols,
	                       indexing_flags(&proj->indexing_params),
	                       proj->cur_image->lambda,
	                       detgeom_mean_camera_length(proj->cur_image->detgeom),
	                       1,
	                       taketwoopts, xgandalf_opts,
	                       pinkIndexer_opts, felix_opts,
	                       NULL);
	free(methods);

	index_pattern(proj->cur_image, ipriv);

	for ( i=0; i<proj->cur_image->n_crystals; i++ ) {
		crystal_set_profile_radius(proj->cur_image->crystals[i], 0.02e9);
		crystal_set_mosaicity(proj->cur_image->crystals[i], 0.0);
		if ( refine_radius(proj->cur_image->crystals[i],
		                   proj->cur_image) )
		{
			ERROR("WARNING: Radius determination failed\n");
		}
	}

	r = chdir(old_cwd);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		return;
	}
	free(old_cwd);
	delete_gui_tempdir(tmpdir);

	err = 0;
	int_method = integration_method(proj->indexing_params.integration_method,
	                                &err);

	integrate_all_5(proj->cur_image, int_method, PMODEL_XSPHERE,
	                proj->indexing_params.push_res,
	                proj->indexing_params.ir_inn,
	                proj->indexing_params.ir_mid,
	                proj->indexing_params.ir_out,
	                INTDIAG_NONE, 0, 0, 0, NULL,
	                proj->indexing_params.overpredict);

	cleanup_indexing(ipriv);

	STATUS("Number of crystals: %i\n",
	       proj->cur_image->n_crystals);
	for ( i=0; i<proj->cur_image->n_crystals; i++ ) {
		cell_print(crystal_get_cell(proj->cur_image->crystals[i]));
	}

}


static void index_one_response_sig(GtkWidget *dialog, gint resp,
                                   struct crystfelproject *proj)
{

	if ( resp == GTK_RESPONSE_OK ) {
		get_indexing_opts(proj,
		                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));
		run_indexing_once(proj);

		crystfel_image_view_set_refl_box_size(CRYSTFEL_IMAGE_VIEW(proj->imageview),
		                                      proj->indexing_params.ir_inn);
		force_refls_on(proj);
		redraw_widget(proj->imageview);
	}

	gtk_widget_destroy(dialog);
	proj->indexing_opts = NULL;
}


gint index_one_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;

	if ( proj->indexing_opts != NULL ) return FALSE;

	dialog = gtk_dialog_new_with_buttons("Index one frame",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(index_one_response_sig), proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	proj->indexing_opts = crystfel_indexing_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox),
	                   GTK_WIDGET(proj->indexing_opts),
	                   FALSE, FALSE, 8.0);
	set_indexing_opts(proj,
	                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);
	crystfel_indexing_opts_set_show_stream_opts(CRYSTFEL_INDEXING_OPTS(proj->indexing_opts),
	                                            FALSE);
	return FALSE;
}


static int contains_spaces(const char *str)
{
	int i;
	size_t len;
	len = strlen(str);
	for ( i=0; i<len; i++ ) {
		if ( str[i] == ' ' ) return 1;
	}
	return 0;
}


static void add_arg(char **args, int pos, const char *label)
{
	if ( contains_spaces(label) ) {
		size_t len = strlen(label)+3;
		args[pos] = malloc(len);
		args[pos][0] = '"';
		args[pos][1] = '\0';
		strcat(args[pos], label);
		args[pos][len-2] = '"';
		args[pos][len-1] = '\0';
	} else {
		args[pos] = strdup(label);
	}
}


static void add_arg_float(char **args, int pos, const char *label,
                          float val)
{
	char *str = malloc(64);
	if ( str == NULL ) return;
	snprintf(str, 63, "--%s=%f", label, val);
	args[pos] = str;
}


static void add_arg_int(char **args, int pos, const char *label,
                        int val)
{
	char *str = malloc(64);
	if ( str == NULL ) return;
	snprintf(str, 63, "--%s=%i", label, val);
	args[pos] = str;
}


static void add_arg_string(char **args, int pos, const char *label,
                           const char *val)
{
	size_t len;
	char *str;

	len = strlen(label)+strlen(val)+4;
	str = malloc(len);
	if ( str == NULL ) return;
	snprintf(str, 63, "--%s=%s", label, val);
	args[pos] = str;
}


static char **indexamajig_command_line(const char *geom_filename,
                                       const char *n_thread_str,
                                       const char *files_list,
                                       const char *stream_filename,
                                       const char *serial_start,
                                       struct peak_params *peak_search_params,
                                       struct index_params *indexing_params)
{
	char **args;
	char tols[2048];
	char *indexamajig_path;
	int i;
	int n_args = 0;

	args = malloc(64*sizeof(char *));
	if ( args == NULL ) return NULL;

	indexamajig_path = get_crystfel_exe("indexamajig");
	if ( indexamajig_path == NULL ) return NULL;

	/* The basics */
	add_arg(args, n_args++, indexamajig_path);
	add_arg(args, n_args++, "-i");
	add_arg(args, n_args++, files_list);
	if ( geom_filename != NULL ) {
		add_arg(args, n_args++, "-g");
		add_arg(args, n_args++, geom_filename);
	}
	add_arg(args, n_args++, "-o");
	add_arg(args, n_args++, stream_filename);
	add_arg(args, n_args++, "-j");
	add_arg(args, n_args++, n_thread_str);

	/* Peak search */
	add_arg(args, n_args++, "--peaks");
	add_arg(args, n_args++, str_peaksearch(peak_search_params->method));

	if ( peak_search_params->method == PEAK_ZAEF ) {
		add_arg_float(args, n_args++, "threshold",
		              peak_search_params->threshold);
		add_arg_float(args, n_args++, "min-squared-gradient",
		              peak_search_params->min_sq_gradient);
		add_arg_float(args, n_args++, "min-snr",
		              peak_search_params->min_snr);

	} else if ( peak_search_params->method == PEAK_PEAKFINDER8 ) {
		add_arg_float(args, n_args++, "threshold",
		              peak_search_params->threshold);
		add_arg_float(args, n_args++, "min-snr",
		              peak_search_params->min_snr);
		add_arg_int(args, n_args++, "min-pix-count",
		            peak_search_params->min_pix_count);
		add_arg_int(args, n_args++, "max-pix-count",
		            peak_search_params->max_pix_count);
		add_arg_int(args, n_args++, "local-bg-radius",
		            peak_search_params->local_bg_radius);
		add_arg_int(args, n_args++, "min-res",
		            peak_search_params->min_res);
		add_arg_int(args, n_args++, "max-res",
		            peak_search_params->max_res);
	}
	snprintf(tols, 2048, "--peak-radius=%.1f,%.1f,%.1f",
	         peak_search_params->pk_inn,
	         peak_search_params->pk_mid,
	         peak_search_params->pk_out);
	add_arg(args, n_args++, tols);

	if ( indexing_params->min_peaks > 0 ) {
		add_arg_int(args, n_args++, "min-peaks",
		            indexing_params->min_peaks);
	}

	/* Indexing */
	if ( indexing_params->indexing_methods != NULL ) {
		add_arg(args, n_args++, "--indexing");
		add_arg(args, n_args++, indexing_params->indexing_methods);
	}
	if ( indexing_params->cell_file != NULL ) {
		add_arg(args, n_args++, "-p");
		add_arg(args, n_args++, indexing_params->cell_file);
	}
	/* indexing_params->tols is in frac (not %) and radians.
	 * Indexamajig command line wants percent and degrees */
	snprintf(tols, 2048, "--tolerance=%f,%f,%f,%f,%f,%f",
	         indexing_params->tols[0]*100.0,
	         indexing_params->tols[1]*100.0,
	         indexing_params->tols[2]*100.0,
	         rad2deg(indexing_params->tols[3]),
	         rad2deg(indexing_params->tols[4]),
	         rad2deg(indexing_params->tols[5]));
	add_arg(args, n_args++, tols);
	if ( indexing_params->multi ) add_arg(args, n_args++, "--multi");
	if ( indexing_params->no_refine ) add_arg(args, n_args++, "--no-refine");
	if ( indexing_params->no_retry ) add_arg(args, n_args++, "--no-retry");
	if ( indexing_params->no_peak_check ) add_arg(args, n_args++, "--no-check-peaks");
	if ( indexing_params->no_cell_check ) add_arg(args, n_args++, "--no-check-cell");

	/* Integration */
	add_arg(args, n_args++, "--integration");
	add_arg(args, n_args++, indexing_params->integration_method);
	if ( indexing_params->overpredict ) add_arg(args, n_args++, "--overpredict");
	if ( !isinf(indexing_params->push_res) ) {
		add_arg_float(args, n_args++, "push-res",
		              indexing_params->push_res);
	}
	snprintf(tols, 2048, "--int-radius=%.1f,%.1f,%.1f",
	         indexing_params->ir_inn,
	         indexing_params->ir_mid,
	         indexing_params->ir_out);
	add_arg(args, n_args++, tols);
	if ( indexing_params->use_fix_profile_radius ) {
		add_arg_float(args, n_args++, "fix-profile-radius",
		              indexing_params->fix_profile_radius);
	}
	add_arg_float(args, n_args++, "fix-divergence",
	              indexing_params->fix_divergence);

	/* Stream output */
	if ( indexing_params->exclude_nonhits ) add_arg(args, n_args++, "--no-non-hits-in-stream");
	if ( indexing_params->exclude_peaks ) add_arg(args, n_args++, "--no-peaks-in-stream");
	if ( indexing_params->exclude_refls ) add_arg(args, n_args++, "--no-refls-in-stream");
	for ( i=0; i<indexing_params->n_metadata; i++ ) {
		add_arg_string(args, n_args++, "copy-header",
		               indexing_params->metadata_to_copy[i]);
	}

	if ( serial_start != NULL ) {
		add_arg_string(args, n_args++, "serial-start", serial_start);
	}

	add_arg_string(args, n_args++, "harvest-file", "parameters.json");

	args[n_args] = NULL;
	return args;
}


int read_number_processed(const char *filename)
{
	FILE *fh = fopen(filename, "r");
	int n_proc = 0;

	/* Normal situation if SLURM job hasn't started yet */
	if ( fh == NULL ) return 0;

	/* Only look at the last part of the file */
	fseek(fh, -4096, SEEK_END);

	do {
		char line[1024];
		if ( fgets(line, 1024, fh) == NULL ) break;

		if ( strncmp(line, "Final: ", 7) == 0 ) {
			int i;
			if ( sscanf(line, "Final: %i images processed", &i) == 1 ) {
				n_proc = i;
			}
		} else if ( strstr(line, " images processed, ") != NULL ) {
			int i;
			if ( sscanf(line, "%i ", &i) == 1 ) {
				n_proc = i;
			}
		}

	} while ( 1 );

	fclose(fh);

	return n_proc;
}


int write_indexamajig_script(const char *script_filename,
                             const char *geom_filename,
                             const char *n_thread_str,
                             const char *files_list,
                             const char *stream_filename,
                             const char *serial_start,
                             int redirect_output,
                             struct peak_params *peak_search_params,
                             struct index_params *indexing_params)
{
	FILE *fh;
	int i;
	char **cmdline;

	cmdline = indexamajig_command_line(geom_filename,
	                                   n_thread_str,
	                                   files_list,
	                                   stream_filename,
	                                   serial_start,
	                                   peak_search_params,
	                                   indexing_params);
	if ( cmdline == NULL ) return 1;

	fh = fopen(script_filename, "w");
	if ( fh == NULL ) return 1;

	fprintf(fh, "#!/bin/sh\n");

	i = 0;
	while ( cmdline[i] != NULL ) {
		fprintf(fh, "%s ", cmdline[i]);
		free(cmdline[i]);
		i++;
	};
	free(cmdline);
	if ( redirect_output ) {
		fprintf(fh, ">stdout.log 2>stderr.log\n");
	}

	fclose(fh);
	return 0;
}
