/*
 * crystfel_gui.c
 *
 * CrystFEL's main graphical user interface
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#include <datatemplate.h>
#include <peaks.h>

#include "crystfelimageview.h"
#include "crystfelimageview.h"
#include "crystfel_gui.h"
#include "gui_peaksearch.h"
#include "gui_index.h"
#include "gui_backend_local.h"


static void show_help(const char *s)
{
	printf("Syntax: %s\n\n", s);
	printf(
"CrystFEL graphical user interface.\n"
"\n"
" -h, --help              Display this help message.\n"
"     --version           Print CrystFEL version number and exit.\n"

);
}


static gboolean destroy_sig(GtkWidget *da, struct crystfelproject *proj)
{
	gtk_main_quit();
	return FALSE;
}


static void add_ui_sig(GtkUIManager *ui, GtkWidget *widget,
                       GtkContainer *container)
{
	gtk_box_pack_start(GTK_BOX(container), widget, FALSE, FALSE, 0);
	if ( GTK_IS_TOOLBAR(widget) ) {
		gtk_toolbar_set_show_arrow(GTK_TOOLBAR(widget), TRUE);
	}
}


static void update_imageview(struct crystfelproject *proj)
{
	char tmp[1024];
	if ( proj->n_frames == 0 ) return;

	snprintf(tmp, 1023, "%s (frame %i of %i)",
	         proj->filenames[proj->cur_frame],
	         proj->cur_frame+1, proj->n_frames);
	gtk_label_set_text(GTK_LABEL(proj->image_info), tmp);
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              proj->filenames[proj->cur_frame],
	                              proj->events[proj->cur_frame]);
}


enum match_type_id
	{
	 MATCH_EVERYTHING,
	 MATCH_CHEETAH_LCLS_H5,
	 MATCH_CHEETAH_CXI,
	 MATCH_CBF,
	 MATCH_CBFGZ,
	};


static int match_filename(const char *fn, enum match_type_id mt)
{
	const char *ext = NULL;
	const char *ext2 = NULL;
	size_t r = strlen(fn)-1;

	while ( r > 0 ) {
		if ( fn[r] == '.' ) {
			if ( ext != NULL ) {
				ext2 = fn+r;
				break;
			} else {
				ext = fn+r;
			}
		}
		r--;
	}

	if ( mt == MATCH_EVERYTHING ) return 1;

	if ( ext == NULL ) return 0;
	if ( mt == MATCH_CHEETAH_LCLS_H5 ) {
		return ((strcmp(ext, ".h5")==0)
		        && (strncmp(fn, "LCLS", 4)==0));
	}
	if ( mt == MATCH_CHEETAH_CXI ) return strcmp(ext, ".cxi")==0;
	if ( mt == MATCH_CBF ) return strcmp(ext, ".cbf")==0;
	if ( mt == MATCH_CBFGZ  ) {
		if ( ext2 != NULL ) {
			return strcmp(ext2, ".cbf.gz")==0;
		}
	}

	return 0;
}


static void add_file_proj(struct crystfelproject *proj,
                          const char *filename)
{
	if ( proj->n_frames == proj->max_frames ) {
		int n_max = proj->max_frames + 1024;
		char **n_filenames;
		char **n_events;
		n_filenames = realloc(proj->filenames,
		                      n_max*sizeof(char *));
		n_events = realloc(proj->events,
		                   n_max*sizeof(char *));
		if ( (n_filenames == NULL) || (n_events == NULL) ) {
			ERROR("Failed to allocate new filename\n");
			return;
		}
		proj->max_frames = n_max;
		proj->filenames = n_filenames;
		proj->events = n_events;
	}

	proj->filenames[proj->n_frames] = strdup(filename);
	proj->events[proj->n_frames] = NULL;
	proj->n_frames++;
}


static void add_files(struct crystfelproject *proj, GFile *folder,
                      enum match_type_id type)
{
	GFileEnumerator *fenum;
	GFileInfo *finfo;
	GError *error = NULL;

	fenum = g_file_enumerate_children(folder, "standard::name,standard::type",
	                                  G_FILE_QUERY_INFO_NONE,
	                                  NULL, &error);

	do {

		GFile *file;

		finfo = g_file_enumerator_next_file(fenum, NULL, &error);

		if ( error != NULL ) {
			STATUS("Error!\n");
			g_object_unref(fenum);
			return;
		}

		if ( finfo == NULL ) continue;

		file = g_file_get_child(folder, g_file_info_get_name(finfo));

		if ( g_file_info_get_file_type(finfo) == G_FILE_TYPE_DIRECTORY ) {

			add_files(proj, file, type);

		} else {

			char *bn = g_file_get_basename(file);
			if ( match_filename(bn, type) ) {
				add_file_proj(proj,
				              g_file_get_path(file));
			}

		}

		g_object_unref(finfo);

	} while ( finfo != NULL );

	g_object_unref(fenum);
}


static enum match_type_id decode_matchtype(const char *type_id)
{
	if ( strcmp(type_id, "everything") == 0 ) return MATCH_EVERYTHING;
	if ( strcmp(type_id, "lcls-cheetah-hdf5") == 0 ) return MATCH_CHEETAH_LCLS_H5;
	if ( strcmp(type_id, "cheetah-cxi") == 0 ) return MATCH_CHEETAH_CXI;
	if ( strcmp(type_id, "cbf") == 0 ) return MATCH_CBF;
	if ( strcmp(type_id, "cbfgz") == 0 ) return MATCH_CBFGZ;
	ERROR("Unknown match type id '%s'\n", type_id);
	return MATCH_EVERYTHING;
}


static void finddata_response_sig(GtkWidget *dialog, gint resp,
                                  struct crystfelproject *proj)
{
	GFile *top;
	DataTemplate *dtempl;
	char *geom_filename;
	const char *type_id;
	int i;

	if ( (resp==GTK_RESPONSE_DELETE_EVENT) || (resp==GTK_RESPONSE_CANCEL) ) {
		gtk_widget_destroy(dialog);
		proj->file_chooser = NULL;
		proj->geom_chooser = NULL;
		return;
	}

	geom_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(proj->geom_chooser));
	if ( geom_filename == NULL ) return;
	dtempl = data_template_new_from_file(geom_filename);
	if ( dtempl == NULL ) return;
	proj->geom_filename = geom_filename;
	crystfel_image_view_set_datatemplate(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                     dtempl);

	top = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(proj->file_chooser));
	if ( top == NULL ) return;

	/* Totally clean up the old list */
	for ( i=0; i<proj->n_frames; i++ ) {
		free(proj->filenames[i]);
		free(proj->events[i]);
	}
	free(proj->filenames);
	free(proj->events);
	proj->n_frames = 0;
	proj->max_frames = 0;
	proj->filenames = NULL;
	proj->events = NULL;

	type_id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->type_combo));
	add_files(proj, top, decode_matchtype(type_id));

	proj->cur_frame = 0;
	update_imageview(proj);

	proj->file_chooser = NULL;
	proj->geom_chooser = NULL;
	g_object_unref(top);
	gtk_widget_destroy(dialog);
}


static gint finddata_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *chooser;
	GtkWidget *combo;

	if ( proj->file_chooser != NULL ) return FALSE;

	dialog = gtk_dialog_new_with_buttons("Find data files",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Find data", GTK_RESPONSE_ACCEPT,
	                                     NULL);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Find data in folder:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	chooser = gtk_file_chooser_button_new("Select a folder",
	                                      GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(chooser), TRUE, TRUE, 2.0);
	proj->file_chooser = chooser;

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Search pattern:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo), TRUE, TRUE, 2.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "everything",
	                "All files in folder and subfolders");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "lcls-cheetah-hdf5",
	                "Individual LCLS files from Cheetah ('LCLS*.h5')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "cheetah-cxi",
	                "Multi-event CXI files from Cheetah ('*.cxi')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "cbf",
	                "Individual CBF files ('*.cbf')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "cbfgz",
	                "Individual gzipped CBF files ('*.cbf.gz')");
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	proj->type_combo = combo;

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Geometry file:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	chooser = gtk_file_chooser_button_new("Select geometry file",
	                                      GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(chooser), TRUE, TRUE, 2.0);
	proj->geom_chooser = chooser;

	g_signal_connect(dialog, "response",
	                 G_CALLBACK(finddata_response_sig), proj);

	gtk_window_set_default_size(GTK_WINDOW(dialog), 512, 0);
	gtk_widget_show_all(dialog);
	return FALSE;
}


static gint quit_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	gtk_main_quit();
	return FALSE;
}


static gint first_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	proj->cur_frame = 0;
	update_imageview(proj);
	update_peaks(proj);
	return FALSE;
}


static gint prev_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->cur_frame == 0 ) return FALSE;
	proj->cur_frame--;
	update_imageview(proj);
	update_peaks(proj);
	return FALSE;
}


static gint next_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->cur_frame == proj->n_frames - 1 ) return FALSE;
	proj->cur_frame++;
	update_imageview(proj);
	update_peaks(proj);
	return FALSE;
}


static gint last_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	proj->cur_frame = proj->n_frames - 1;
	update_imageview(proj);
	update_peaks(proj);
	return FALSE;
}


static gint about_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *window;

	const gchar *authors[] = {
		"Thomas White <taw@physics.org>",
		NULL
	};

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(proj->window));

	gtk_about_dialog_set_program_name(GTK_ABOUT_DIALOG(window),
	        "CrystFEL graphical user interface");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window), CRYSTFEL_VERSIONSTRING);
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"© 2020 Deutsches Elektronen-Synchrotron DESY, "
		"a research centre of the Helmholtz Association.");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"https://www.desy.de/~twhite/crystfel");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static gint show_peaks_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_peaks = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	update_peaks(proj);
	return FALSE;
}


static void add_menu_bar(struct crystfelproject *proj, GtkWidget *vbox)
{
	GError *error = NULL;

	const char *ui = "<ui> <menubar name=\"mainwindow\">"
		"<menu name=\"file\" action=\"FileAction\">"
		"	<menuitem name=\"quit\" action=\"QuitAction\" />"
		"</menu>"
		"<menu name=\"view\" action=\"ViewAction\" >"
		"	<menuitem name=\"peaks\" action=\"PeaksAction\" />"
		"</menu>"
		"<menu name=\"tools\" action=\"ToolsAction\" >"
		"</menu>"
		"<menu name=\"help\" action=\"HelpAction\">"
		"	<menuitem name=\"about\" action=\"AboutAction\" />"
		"</menu>"
		"</menubar></ui>";

	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "QuitAction", GTK_STOCK_QUIT, "_Quit", NULL, NULL,
			G_CALLBACK(quit_sig) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About", NULL, NULL,
			G_CALLBACK(about_sig) },

	};

	GtkToggleActionEntry toggles[] = {
		{ "PeaksAction", NULL, "Peak detection results", NULL, NULL,
		  G_CALLBACK(show_peaks_sig), FALSE },
	};

	proj->action_group = gtk_action_group_new("cellwindow");
	gtk_action_group_add_actions(proj->action_group, entries,
	                             G_N_ELEMENTS(entries), proj);
	gtk_action_group_add_toggle_actions(proj->action_group, toggles,
	                                    G_N_ELEMENTS(toggles), proj);

	proj->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(proj->ui, proj->action_group, 0);
	g_signal_connect(proj->ui, "add_widget", G_CALLBACK(add_ui_sig), vbox);
	if ( gtk_ui_manager_add_ui_from_string(proj->ui, ui, -1, &error) == 0 )
	{
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(proj->window),
				   gtk_ui_manager_get_accel_group(proj->ui));
	gtk_ui_manager_ensure_update(proj->ui);
}


static void add_button(GtkWidget *vbox, const char *label, const char *imagen,
                       GCallback callback, struct crystfelproject *proj)
{
	GtkWidget *button;
	GtkWidget *image;

	button = gtk_button_new_with_label(label);
	g_object_set(G_OBJECT(button), "image-position", GTK_POS_TOP, NULL);
	image = gtk_image_new_from_icon_name(imagen, GTK_ICON_SIZE_DIALOG);
	g_object_set(G_OBJECT(button), "image", image, NULL);
	g_object_set(G_OBJECT(button), "always-show-image", TRUE, NULL);
	gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(button), FALSE, FALSE, 4);

	if ( callback != NULL ) {
		g_signal_connect(G_OBJECT(button), "clicked",
		                 G_CALLBACK(callback), proj);
	}
}


static void add_task_buttons(GtkWidget *vbox, struct crystfelproject *proj)
{
	add_button(vbox, "Find data", "folder-pictures",
	           G_CALLBACK(finddata_sig), proj);
	add_button(vbox, "Peak detection", "edit-find",
	           G_CALLBACK(peaksearch_sig), proj);
	add_button(vbox, "Determine unit cell", "document-page-setup",
	           G_CALLBACK(unitcell_sig), proj);
	add_button(vbox, "Index and integrate", "system-run",
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Merge", "applications-science",
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Figures of merit", "trophy-gold",
	           G_CALLBACK(NULL), proj);
}


static void add_gui_message(enum log_msg_type type, const char *msg,
                            void *vp)
{
	GtkTextBuffer *buf;
	GtkTextIter iter;
	struct crystfelproject *proj = vp;
	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(proj->report));
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, msg, -1);
}


int main(int argc, char *argv[])
{
	int c;
	struct crystfelproject proj;
	GtkWidget *vbox;
	GtkWidget *vpaned;
	GtkWidget *hpaned;
	GtkWidget *scroll;
	GtkWidget *frame;
	GtkWidget *main_vbox;
	GtkWidget *toolbar;
	GtkWidget *button;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                1 },
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "h",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 1 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			default :
			return 1;

		}

	}

	/* This isn't great, but necessary to make the command-line UI and file
	 * formats consistent with the other programs, which all use the C
	 * locale.  Better would be to have all the programs call
	 * setlocale(LC_ALL, "") and use the C locale temporarily when reading
	 * or writing a stream, reflection file, geometry file etc. */
	gtk_disable_setlocale();

	gtk_init(&argc, &argv);

	proj.file_chooser = NULL;
	proj.geom_chooser = NULL;
	proj.geom_filename = NULL;
	proj.show_peaks = 0;
	proj.n_frames = 0;
	proj.max_frames = 0;
	proj.filenames = NULL;
	proj.events = NULL;
	proj.peak_params = NULL;
	proj.unitcell_combo = NULL;
	proj.peak_search_params.method = PEAK_ZAEF;
	proj.peak_search_params.threshold = 800.0;
	proj.peak_search_params.min_sq_gradient = 100000;
	proj.peak_search_params.min_snr = 5.0;
	proj.peak_search_params.local_bg_radius = 3;
	proj.peak_search_params.min_res = 0;
	proj.peak_search_params.min_sig = 11.0;
	proj.peak_search_params.max_res = 1200;
	proj.peak_search_params.min_pix_count = 2;
	proj.peak_search_params.max_pix_count = 200;
	proj.peak_search_params.min_peak_over_neighbour = -INFINITY;
	proj.peak_search_params.pk_inn = 3.0;
	proj.peak_search_params.pk_mid = 4.0;
	proj.peak_search_params.pk_out = 5.0;
	proj.backend = backend_local;

	proj.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(proj.window), "CrystFEL");
	g_signal_connect(G_OBJECT(proj.window), "destroy", G_CALLBACK(destroy_sig),
	                 &proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_container_add(GTK_CONTAINER(proj.window), vbox);
	add_menu_bar(&proj, vbox);

	vpaned = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
	gtk_box_pack_end(GTK_BOX(vbox), vpaned, TRUE, TRUE, 0.0);

	hpaned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_paned_pack1(GTK_PANED(vpaned), hpaned, TRUE, TRUE);

	proj.imageview = crystfel_image_view_new();

	proj.cur_frame = 0;
	update_imageview(&proj);

	toolbar = gtk_hbox_new(FALSE, 0.0);

	/* First */
	button = gtk_button_new_from_icon_name("go-first", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 0.0);
	g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(first_frame_sig), &proj);

	/* Prev */
	button = gtk_button_new_from_icon_name("go-previous", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 0.0);
	g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(prev_frame_sig), &proj);

	/* Next */
	button = gtk_button_new_from_icon_name("go-next", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 0.0);
	g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(next_frame_sig), &proj);

	/* Last */
	button = gtk_button_new_from_icon_name("go-last", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 0.0);
	g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(last_frame_sig), &proj);

	proj.image_info = gtk_label_new("Ready to load images");
	gtk_label_set_selectable(GTK_LABEL(proj.image_info), TRUE);
	gtk_label_set_ellipsize(GTK_LABEL(proj.image_info),
	                        PANGO_ELLIPSIZE_START);
	gtk_box_pack_end(GTK_BOX(toolbar), proj.image_info, TRUE, TRUE, 0.0);

	main_vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(main_vbox), toolbar, FALSE, FALSE, 0.0);

	/* Main area stuff (toolbar and imageview) at right */
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_ALWAYS, GTK_POLICY_ALWAYS);
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.imageview));
	gtk_box_pack_start(GTK_BOX(main_vbox), scroll, TRUE, TRUE, 0.0);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(main_vbox));
	gtk_paned_pack2(GTK_PANED(hpaned), GTK_WIDGET(frame), TRUE, TRUE);

	/* Icon region at left */
	proj.icons = gtk_vbox_new(FALSE, 0.0);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_set_border_width(GTK_CONTAINER(proj.icons), 16);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(scroll));
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.icons));
	gtk_paned_pack1(GTK_PANED(hpaned), GTK_WIDGET(frame), FALSE, FALSE);
	add_task_buttons(proj.icons, &proj);

	/* Report (text) region at bottom */
	proj.report = gtk_text_view_new();
	gtk_text_view_set_editable(GTK_TEXT_VIEW(proj.report), FALSE);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(proj.report), FALSE);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(scroll));
	gtk_container_add(GTK_CONTAINER(scroll), GTK_WIDGET(proj.report));
	gtk_paned_pack2(GTK_PANED(vpaned), GTK_WIDGET(frame), FALSE, FALSE);

	/* Send messages to report region */
	set_log_message_func(add_gui_message, &proj);

	gtk_window_set_default_size(GTK_WINDOW(proj.window), 1024, 768);
	gtk_paned_set_position(GTK_PANED(hpaned), 172);
	gtk_paned_set_position(GTK_PANED(vpaned), 600);
	gtk_widget_show_all(proj.window);
	gtk_main();

	return 0;
}
