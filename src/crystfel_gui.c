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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <datatemplate.h>
#include <peaks.h>
#include <cell-utils.h>

#include "crystfelimageview.h"
#include "crystfelimageview.h"
#include "crystfel_gui.h"
#include "gui_peaksearch.h"
#include "gui_index.h"
#include "gui_merge.h"
#include "gui_fom.h"
#include "gui_export.h"
#include "gui_ambi.h"
#include "gui_project.h"
#include "version.h"


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

static int confirm_exit(struct crystfelproject *proj)
{
	int r;
	GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(proj->window),
	                                           0,
	                                           GTK_MESSAGE_QUESTION,
	                                           GTK_BUTTONS_NONE,
	                                           "Do you want to save the session?");
	gtk_dialog_add_buttons(GTK_DIALOG(dialog),
	                       "Save", GTK_RESPONSE_YES,
	                       "Don't save", GTK_RESPONSE_NO,
	                       "Cancel", GTK_RESPONSE_CANCEL,
	                       NULL);
	r = gtk_dialog_run(GTK_DIALOG(dialog));
	gtk_widget_destroy(dialog);
	if ( r == GTK_RESPONSE_YES ) {
		save_project(proj);
		return 1;
	}
	if ( r == GTK_RESPONSE_NO ) return 1;
	return 0;
}


/* Main window destroyed */
static gboolean delete_event_sig(GtkWidget *da, GdkEvent *event,
                                 struct crystfelproject *proj)
{
	if ( proj->unsaved ) {
		if ( !confirm_exit(proj) ) return TRUE;
	}
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


static void swap_data_arrays(struct image *a, struct image *b)
{
	float **swap;
	int **swap_bad;

	swap = a->dp;
	a->dp = b->dp;
	b->dp = swap;

	swap = a->sat;
	a->sat = b->sat;
	b->sat = swap;

	swap_bad = a->bad;
	a->bad = b->bad;
	b->bad = swap_bad;
}


void select_result(struct crystfelproject *proj,
                   const char *result_name)
{
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(proj->results_combo),
	                            result_name);
}


const char *selected_result(struct crystfelproject *proj)
{
	return gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
}


/* Return non-zero if there are any jobs running.
 * If not, it makes no sense to re-scan any result streams.
 * Possible future improvement: exclude jobs which don't produce streams
 *  (i.e. merging) */
static int have_running_jobs(struct crystfelproject *proj)
{
	int i;

	for ( i=0; i<proj->n_running_tasks; i++ ) {
		if ( proj->tasks[i].running ) return 1;
	}

	return 0;
}


/* Bring the image view up to date after changing the selected image */
void update_imageview(struct crystfelproject *proj)
{
	char tmp[1024];
	char *ev_str;
	char *ev_sep;
	struct image *image;
	const gchar *results_name;

	if ( proj->n_frames == 0 ) return;

	image = image_read(proj->dtempl,
	                   proj->filenames[proj->cur_frame],
	                   proj->events[proj->cur_frame],
	                   0, 0);

	if ( image == NULL ) {
		ERROR("Failed to load image\n");
		return;
	}

	/* Give CrystFELImageView a chance to free resources */
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              NULL);
	image_free(proj->cur_image);
	proj->cur_image = image;

	if ( proj->cur_image->ev != NULL ) {
		ev_str = proj->cur_image->ev;
		ev_sep = " ";
	} else {
		ev_str = "";
		ev_sep = "";
	}
	snprintf(tmp, 1023, "%s%s%s (frame %i of %i)",
	         proj->cur_image->filename,
	         ev_sep,
	         ev_str,
	         proj->cur_frame+1,
	         proj->n_frames);
	gtk_label_set_text(GTK_LABEL(proj->image_info), tmp);

	/* Look up results, if applicable */
	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	if ( strcmp(results_name, "crystfel-gui-internal") == 0 ) {
		update_peaks(proj);
	} else {
		struct image *res_im;

		res_im = find_indexed_image(proj,
		                            results_name,
		                            image->filename,
		                            image->ev,
		                            have_running_jobs(proj));
		if ( res_im != NULL ) {
			swap_data_arrays(image, res_im);
			image_free(proj->cur_image);
			proj->cur_image = res_im;
		}
	}

	crystfel_image_view_set_show_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->show_refls);
	crystfel_image_view_set_label_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                          proj->label_refls);
	crystfel_image_view_set_refl_box_size(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                      proj->indexing_params.ir_inn);
	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_peaks);
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              proj->cur_image);

	gtk_widget_set_sensitive(proj->next_button,
	                         !(proj->cur_frame == proj->n_frames-1));
	gtk_widget_set_sensitive(proj->last_button,
	                         !(proj->cur_frame == proj->n_frames-1));
	gtk_widget_set_sensitive(proj->prev_button,
	                         !(proj->cur_frame == 0));
	gtk_widget_set_sensitive(proj->first_button,
	                         !(proj->cur_frame == 0));
}


static void add_all_events(struct crystfelproject *proj,
                           const char *filename,
                           const DataTemplate *dtempl)
{
	char **events;
	int i;
	int n_events;

	events = image_expand_frames(dtempl, filename, &n_events);
	if ( events == NULL ) {
		ERROR("Couldn't expand event list\n");
		return;
	}

	for ( i=0; i<n_events; i++ ) {
		add_file_to_project(proj, filename, events[i]);
		free(events[i]);
	}
	free(events);
}


static void add_files(struct crystfelproject *proj, GFile *folder,
                      enum match_type_id type,
                      const DataTemplate *dtempl)
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

			add_files(proj, file, type, dtempl);

		} else {

			char *bn = g_file_get_basename(file);
			if ( match_filename(bn, type) ) {
				add_all_events(proj, g_file_get_path(file),
				               dtempl);
			}

		}

		g_object_unref(finfo);

	} while ( finfo != NULL );

	g_object_unref(fenum);
}


static void add_frames_from_stream(Stream *st,
                                   DataTemplate *dtempl,
                                   struct crystfelproject *proj)
{
	do {
		struct image *image;
		image = stream_read_chunk(st, 0);
		if ( image == NULL ) break;
		add_file_to_project(proj, image->filename, image->ev);
		image_free(image);

	} while ( 1 );
}


struct finddata_ctx
{
	struct crystfelproject *proj;
	GtkWidget *geom_file;

	/* "Select individual file" */
	GtkWidget *indiv;
	GtkWidget *indiv_chooser;

	/* Read list of files */
	GtkWidget *list;  /* "Import list" radio */
	GtkWidget *list_chooser;

	/* Search for files */
	GtkWidget *search;
	GtkWidget *search_chooser;
	GtkWidget *search_pattern;

	/* Load stream */
	GtkWidget *stream;
	GtkWidget *stream_chooser;

	GtkWidget *dump;
};

enum import_mode
{
	IMPORT_FILES,
	IMPORT_LIST,
	IMPORT_SEARCH,
	IMPORT_STREAM
};


static enum import_mode import_mode(struct finddata_ctx *ctx)
{
	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->indiv))) {
		return IMPORT_FILES;
	} else if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->list))) {
		return IMPORT_LIST;
	} else if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->search))) {
		return IMPORT_SEARCH;
	} else {
		return IMPORT_STREAM;
	}
}


static void finddata_typetoggle_sig(GtkWidget *radio,
                                    struct finddata_ctx *ctx)
{
	gtk_widget_set_sensitive(ctx->indiv_chooser, FALSE);
	gtk_widget_set_sensitive(ctx->list_chooser, FALSE);
	gtk_widget_set_sensitive(ctx->search_chooser, FALSE);
	gtk_widget_set_sensitive(ctx->search_pattern, FALSE);
	gtk_widget_set_sensitive(ctx->stream_chooser, FALSE);

	gtk_widget_set_sensitive(ctx->geom_file, TRUE);

	switch ( import_mode(ctx) ) {

		case IMPORT_FILES :
		gtk_widget_set_sensitive(ctx->indiv_chooser, TRUE);
		break;

		case IMPORT_LIST :
		gtk_widget_set_sensitive(ctx->list_chooser, TRUE);
		break;

		case IMPORT_SEARCH :
		gtk_widget_set_sensitive(ctx->search_chooser, TRUE);
		gtk_widget_set_sensitive(ctx->search_pattern, TRUE);
		break;

		case IMPORT_STREAM :
		gtk_widget_set_sensitive(ctx->geom_file, FALSE);
		gtk_widget_set_sensitive(ctx->stream_chooser, TRUE);
		break;
	}
}


static void import_via_search(struct finddata_ctx *ctx)
{
	GFile *top;
	DataTemplate *dtempl;
	char *geom_filename;
	const char *type_id;
	struct crystfelproject *proj = ctx->proj;

	geom_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->geom_file));
	if ( geom_filename == NULL ) return;

	top = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(ctx->search_chooser));
	if ( top == NULL ) return;

	dtempl = data_template_new_from_file(geom_filename);
	if ( dtempl == NULL ) return;

	type_id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(ctx->search_pattern));
	proj->data_search_pattern = decode_matchtype(type_id);

	/* Totally clean up the old list */
	clear_project_files(proj);
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              NULL);

	g_free(proj->geom_filename);
	proj->geom_filename = geom_filename;

	data_template_free(proj->dtempl);
	proj->dtempl = dtempl;

	g_free(proj->data_top_folder);
	proj->data_top_folder = g_file_get_path(top);

	add_files(proj, top, proj->data_search_pattern,
	          proj->dtempl);

	g_object_unref(top);
}


static void import_stream(struct finddata_ctx *ctx)
{
	struct crystfelproject *proj = ctx->proj;
	Stream *st;
	char *stream_filename;
	DataTemplate *dtempl;
	const char *geom_str;
	char **streams;

	stream_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->stream_chooser));
	if ( stream_filename == NULL ) return;

	st = stream_open_for_read(stream_filename);
	if ( st == NULL ) return;

	geom_str = stream_geometry_file(st);
	if ( geom_str == NULL ) {
		ERROR("No geometry file\n");
		stream_close(st);
		return;
	}

	dtempl = data_template_new_from_string(geom_str);
	if ( dtempl == NULL ) {
		stream_close(st);
		return;
	}

	clear_project_files(proj);
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              NULL);

	data_template_free(proj->dtempl);
	proj->dtempl = dtempl;

	/* Set some defaults for things we won't be using */
	g_free(proj->geom_filename);
	proj->geom_filename = NULL;
	g_free(proj->data_top_folder);
	proj->data_top_folder = NULL;
	proj->data_search_pattern = MATCH_EVERYTHING;

	add_frames_from_stream(st, proj->dtempl, proj);
	proj->stream_filename = stream_filename;
	stream_close(st);

	streams = malloc(sizeof(char *));
	if ( streams != NULL ) {
		char *result_name = safe_basename(stream_filename);
		streams[0] = strdup(stream_filename);
		add_indexing_result(proj, result_name, streams, 1);
		select_result(proj, result_name);
	}

	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   1);
}

static void finddata_response_sig(GtkWidget *dialog, gint resp,
                                  struct finddata_ctx *ctx)
{
	struct crystfelproject *proj = ctx->proj;

	if ( (resp == GTK_RESPONSE_DELETE_EVENT)
	  || (resp == GTK_RESPONSE_CANCEL) )
	{
		gtk_widget_destroy(dialog);
		free(ctx);
		return;
	}

	switch ( import_mode(ctx) ) {

		case IMPORT_FILES :
		/* FIXME */
		break;

		case IMPORT_LIST :
		/* FIXME */
		break;

		case IMPORT_SEARCH :
		import_via_search(ctx);
		break;

		case IMPORT_STREAM :
		import_stream(ctx);
		break;
	}

	proj->unsaved = 1;
	proj->cur_frame = 0;
	crystfel_image_view_reset_zoom(CRYSTFEL_IMAGE_VIEW(proj->imageview));
	update_imageview(proj);

	free(ctx);
	gtk_widget_destroy(dialog);
}


static gint finddata_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	struct finddata_ctx *ctx;

	ctx = malloc(sizeof(struct finddata_ctx));
	if ( ctx == NULL ) return FALSE;

	ctx->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Import data",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Import", GTK_RESPONSE_ACCEPT,
	                                     NULL);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Geometry file:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	ctx->geom_file = gtk_file_chooser_button_new("Select geometry file",
	                                             GTK_FILE_CHOOSER_ACTION_OPEN);
	if ( proj->geom_filename != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->geom_file),
		                              proj->geom_filename);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->geom_file), TRUE, TRUE, 2.0);

	/* Select individual files */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->indiv = gtk_radio_button_new_with_label(NULL,
	                                             "Select an individual file");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->indiv),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(ctx->indiv, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	ctx->indiv_chooser = gtk_file_chooser_button_new("Select file",
	                                                 GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->indiv_chooser),
	                   FALSE, FALSE, 4.0);

	/* Pre-prepared list of files */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->list = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(ctx->indiv),
	                                                        "Read a list of files");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->list),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(ctx->list, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	ctx->list_chooser = gtk_file_chooser_button_new("Select the list of filenames",
	                                                GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->list_chooser),
	                   FALSE, FALSE, 4.0);

	/* Search in folder */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 2.0);
	gtk_widget_set_margin_top(hbox, 6.0);
	ctx->search = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(ctx->indiv),
	                                                      "Search for files in folder");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->search),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(ctx->search, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	ctx->search_chooser = gtk_file_chooser_button_new("Select a folder",
	                                              GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
	if ( proj->data_top_folder != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->search_chooser),
		                              proj->data_top_folder);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->search_chooser),
	                   TRUE, TRUE, 2.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_widget_set_margin_bottom(hbox, 6.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 2.0);
	label = gtk_label_new("Search pattern:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_widget_set_margin_start(label, 32);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	ctx->search_pattern = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->search_pattern), TRUE, TRUE, 2.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "everything",
	                "All files in folder and subfolders");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "hdf5",
	                "All HDF5 files ('*.h5')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "lcls-cheetah-hdf5",
	                "Individual LCLS files from Cheetah ('LCLS*.h5')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "cheetah-cxi",
	                "Multi-event CXI files from Cheetah ('*.cxi')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "cbf",
	                "Individual CBF files ('*.cbf')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->search_pattern), "cbfgz",
	                "Individual gzipped CBF files ('*.cbf.gz')");
	gtk_combo_box_set_active(GTK_COMBO_BOX(ctx->search_pattern),
	                         proj->data_search_pattern);

	/* Load a stream */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->stream = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(ctx->indiv),
	                                                          "Load stream");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->stream),
	                   FALSE, FALSE, 4.0);
	ctx->stream_chooser = gtk_file_chooser_button_new("Select stream file",
	                                                  GTK_FILE_CHOOSER_ACTION_OPEN);
	if ( proj->stream_filename != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->stream_chooser),
		                              proj->stream_filename);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ctx->stream), TRUE);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->stream_chooser),
	                   TRUE, TRUE, 2.0);

	/* Replace data toggle */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->dump = gtk_check_button_new_with_label("Replace all the current data");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->dump), FALSE, FALSE, 4.0);

	g_signal_connect(dialog, "response",
	                 G_CALLBACK(finddata_response_sig), ctx);

	gtk_window_set_default_size(GTK_WINDOW(dialog), 512, 0);
	finddata_typetoggle_sig(ctx->search, ctx);
	gtk_widget_show_all(dialog);
	return FALSE;
}


/* File->Quit */
static gint quit_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->unsaved ) {
		if ( !confirm_exit(proj) ) return TRUE;
	}
	gtk_main_quit();
	return FALSE;
}


static gint save_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( save_project(proj) == 0 ) {
		STATUS("Saved project.\n");
	} else {
		ERROR("Could not save project.\n");
	}
	return FALSE;
}


static gint rescan_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	const char *results_name;

	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	if ( strcmp(results_name, "crystfel-gui-internal") != 0 ) {
		struct gui_indexing_result *res;
		res = find_indexing_result_by_name(proj, results_name);
		if ( res != NULL ) {
			update_result_index(res);
		} else {
			ERROR("Couldn't find result '%s'\n", results_name);
		}
	}
	return FALSE;
}


static gint reset_zoom_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	crystfel_image_view_reset_zoom(CRYSTFEL_IMAGE_VIEW(proj->imageview));
	return FALSE;
}


static gint first_frame_sig(GtkWidget *widget,
                            struct crystfelproject *proj)
{
	proj->cur_frame = 0;
	update_imageview(proj);
	return FALSE;
}


static gint prev_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	if ( proj->cur_frame == 0 ) return FALSE;
	proj->cur_frame--;
	update_imageview(proj);
	return FALSE;
}


static gint random_frame_sig(GtkWidget *widget,
                             struct crystfelproject *proj)
{
	proj->cur_frame = random()*proj->n_frames / RAND_MAX;
	update_imageview(proj);
	return FALSE;
}


static gint next_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	if ( proj->cur_frame == proj->n_frames - 1 ) return FALSE;
	proj->cur_frame++;
	update_imageview(proj);
	return FALSE;
}


static gint last_frame_sig(GtkWidget *widget,
                           struct crystfelproject *proj)
{
	proj->cur_frame = proj->n_frames - 1;
	update_imageview(proj);
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
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window),
	                             crystfel_version_string());
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


static gint image_info_clicked_sig(GtkWidget *widget,
                                   struct crystfelproject *proj)
{
	GtkWidget *popover;
	GtkWidget *grid;
	GtkWidget *label;
	char tmp[64];

	if ( proj->cur_image == NULL ) return FALSE;

	popover = gtk_popover_new(widget);
	gtk_popover_set_position(GTK_POPOVER(popover),
	                         GTK_POS_BOTTOM);

	grid = gtk_grid_new();
	gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 4);
	gtk_container_set_border_width(GTK_CONTAINER(grid), 6);

	label = gtk_label_new("Number of peaks:");
	gtk_grid_attach(GTK_GRID(grid), label, 0, 0, 1, 1);
	snprintf(tmp, 63, "%i", image_feature_count(proj->cur_image->features));
	label = gtk_label_new(tmp);
	gtk_grid_attach(GTK_GRID(grid), label, 1, 0, 1, 1);

	label = gtk_label_new("Number of crystals:");
	gtk_grid_attach(GTK_GRID(grid), label, 0, 1, 1, 1);
	snprintf(tmp, 63, "%i", proj->cur_image->n_crystals);
	label = gtk_label_new(tmp);
	gtk_grid_attach(GTK_GRID(grid), label, 1, 1, 1, 1);

	gtk_container_add(GTK_CONTAINER(popover), grid);
	gtk_widget_show_all(grid);

#if GTK_CHECK_VERSION(3,22,0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show_all(GTK_WIDGET(popover));
#endif

	return FALSE;
}


static gint results_combo_changed_sig(GtkComboBox *w,
                                      struct crystfelproject *proj)
{
	update_imageview(proj);
	return FALSE;
}


static gint show_peaks_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_peaks = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   proj->show_peaks);
	return FALSE;
}


static gint show_refls_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->show_refls = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_show_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                         proj->show_refls);
	return FALSE;
}


static gint label_refls_sig(GtkWidget *w, struct crystfelproject *proj)
{
	proj->label_refls = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(w));
	crystfel_image_view_set_label_reflections(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                          proj->label_refls);
	return FALSE;
}


static void add_menu_bar(struct crystfelproject *proj, GtkWidget *vbox)
{
	GError *error = NULL;

	const char *ui = "<ui> <menubar name=\"mainwindow\">"
		"<menu name=\"file\" action=\"FileAction\">"
		"	<menuitem name=\"save\" action=\"SaveAction\" />"
		"	<menuitem name=\"quit\" action=\"QuitAction\" />"
		"</menu>"
		"<menu name=\"view\" action=\"ViewAction\" >"
		"	<menuitem name=\"peaks\" action=\"PeaksAction\" />"
		"	<menuitem name=\"refls\" action=\"ReflsAction\" />"
		"	<menuitem name=\"labelrefls\" action=\"LabelReflsAction\" />"
		"       <separator />"
		"	<menuitem name=\"resetzoom\" action=\"ResetZoomAction\" />"
		"</menu>"
		"<menu name=\"tools\" action=\"ToolsAction\" >"
		"	<menuitem name=\"rescan\" action=\"RescanAction\" />"
		"</menu>"
		"<menu name=\"help\" action=\"HelpAction\">"
		"	<menuitem name=\"about\" action=\"AboutAction\" />"
		"</menu>"
		"</menubar></ui>";

	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "SaveAction", GTK_STOCK_SAVE, "_Save", NULL, NULL,
			G_CALLBACK(save_sig) },
		{ "QuitAction", GTK_STOCK_QUIT, "_Quit", NULL, NULL,
			G_CALLBACK(quit_sig) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },
		{ "ResetZoomAction", NULL, "Reset zoom", NULL, NULL,
			G_CALLBACK(reset_zoom_sig) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "RescanAction", NULL, "Rescan streams", NULL, NULL,
			G_CALLBACK(rescan_sig) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About", NULL, NULL,
			G_CALLBACK(about_sig) },

	};

	GtkToggleActionEntry toggles[] = {
		{ "PeaksAction", NULL, "Peak detection results", NULL, NULL,
		  G_CALLBACK(show_peaks_sig), FALSE },
		{ "ReflsAction", NULL, "Calculated reflection positions", NULL, NULL,
		  G_CALLBACK(show_refls_sig), FALSE },
		{ "LabelReflsAction", NULL, "Show reflection indices", NULL, NULL,
		  G_CALLBACK(label_refls_sig), FALSE },
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
	/* FIXME: All these icons are placeholders */
	add_button(vbox, "Load data", "folder-pictures",
	           G_CALLBACK(finddata_sig), proj);
	add_button(vbox, "Peak detection", "edit-find",
	           G_CALLBACK(peaksearch_sig), proj);
	add_button(vbox, "Index this frame", "system-run",
	           G_CALLBACK(index_one_sig), proj);
	add_button(vbox, "Index all frames", "view-grid",
	           G_CALLBACK(index_all_sig), proj);
	add_button(vbox, "Determine unit cell", "applications-engineering",
	           G_CALLBACK(cell_explorer_sig), proj);
	add_button(vbox, "Indexing ambiguity", "face-worried",
	           G_CALLBACK(ambi_sig), proj);
	add_button(vbox, "Merge", "applications-science",
	           G_CALLBACK(merge_sig), proj);
	add_button(vbox, "Figures of merit", "trophy-gold",
	           G_CALLBACK(fom_sig), proj);
	add_button(vbox, "Export data", "document-send",
	           G_CALLBACK(export_sig), proj);
}


static void add_gui_message(enum log_msg_type type, const char *msg,
                            void *vp)
{
	GtkTextBuffer *buf;
	GtkTextIter iter;
	GtkTextMark *mark;
	struct crystfelproject *proj = vp;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(proj->report));
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, msg, -1);

	mark = gtk_text_mark_new(NULL, FALSE);
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_add_mark(buf, mark, &iter);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(proj->report),
	                             mark, 0.0, FALSE, 0.0, 0.0);
	gtk_text_buffer_delete_mark(buf, mark);
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
	GtkWidget *results_toolbar;
	GtkWidget *button;
	GtkWidget *label;

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
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
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

	default_project(&proj);

	proj.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(proj.window), "CrystFEL");
	gtk_window_set_default_icon_name("crystfel");
	g_signal_connect(G_OBJECT(proj.window), "delete-event",
	                 G_CALLBACK(delete_event_sig),
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
	crystfel_image_view_reset_zoom(CRYSTFEL_IMAGE_VIEW(proj.imageview));

	toolbar = gtk_hbox_new(FALSE, 0.0);

	/* First */
	proj.first_button = gtk_button_new_from_icon_name("go-first", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.first_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.first_button), "clicked",
	                 G_CALLBACK(first_frame_sig), &proj);

	/* Prev */
	proj.prev_button = gtk_button_new_from_icon_name("go-previous", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.prev_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.prev_button), "clicked",
	                 G_CALLBACK(prev_frame_sig), &proj);

	/* Random */
	button = gtk_button_new_from_icon_name("media-playlist-shuffle",
	                                       GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(random_frame_sig), &proj);

	/* Next */
	proj.next_button = gtk_button_new_from_icon_name("go-next", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.next_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.next_button), "clicked",
	                 G_CALLBACK(next_frame_sig), &proj);

	/* Last */
	proj.last_button = gtk_button_new_from_icon_name("go-last", GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_start(GTK_BOX(toolbar), proj.last_button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(proj.last_button), "clicked",
	                 G_CALLBACK(last_frame_sig), &proj);

	/* Information about image */
	button = gtk_button_new_from_icon_name("document-properties",
	                                       GTK_ICON_SIZE_LARGE_TOOLBAR);
	gtk_box_pack_end(GTK_BOX(toolbar), button, FALSE, FALSE, 2.0);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(image_info_clicked_sig), &proj);

	results_toolbar = gtk_hbox_new(FALSE, 0.0);
	label = gtk_label_new("Show results from:");
	gtk_box_pack_start(GTK_BOX(results_toolbar), label,
	                   FALSE, FALSE, 4.0);
	proj.results_combo = gtk_combo_box_text_new();
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(proj.results_combo),
	                          "crystfel-gui-internal",
	                          "Calculations within GUI");
	gtk_combo_box_set_active(GTK_COMBO_BOX(proj.results_combo), 0);
	gtk_box_pack_start(GTK_BOX(results_toolbar), proj.results_combo,
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(proj.results_combo), "changed",
	                 G_CALLBACK(results_combo_changed_sig), &proj);

	/* Filename */
	proj.image_info = gtk_label_new("Ready to load images");
	gtk_label_set_selectable(GTK_LABEL(proj.image_info), TRUE);
	gtk_label_set_ellipsize(GTK_LABEL(proj.image_info),
	                        PANGO_ELLIPSIZE_START);
	gtk_box_pack_end(GTK_BOX(toolbar), proj.image_info, TRUE, TRUE, 0.0);

	main_vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(main_vbox), toolbar,
	                   FALSE, FALSE, 2.0);
	gtk_box_pack_start(GTK_BOX(main_vbox), results_toolbar,
	                   FALSE, FALSE, 2.0);

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
	proj.main_vbox = main_vbox;

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

	/* Load state from disk */
	if ( load_project(&proj) == 0 ) {
		DataTemplate *dtempl;
		GtkAction *w;
		proj.cur_frame = 0;

		if ( proj.geom_filename != NULL ) {

			dtempl = data_template_new_from_file(proj.geom_filename);
			if ( dtempl != NULL ) {
				proj.dtempl = dtempl;
			}
		} else if ( proj.stream_filename != NULL ) {

			Stream *st;
			st = stream_open_for_read(proj.stream_filename);
			if ( st != NULL ) {
				char *geom_str = stream_geometry_file(st);
				if ( geom_str != NULL ) {
					proj.dtempl = data_template_new_from_string(geom_str);
				}
			}
			stream_close(st);
		}

		w = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/peaks");
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(w),
		                             proj.show_peaks);

		w = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/refls");
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(w),
		                             proj.show_refls);

		w = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/labelrefls");
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(w),
		                             proj.label_refls);

		update_imageview(&proj);
	}

	gtk_window_set_default_size(GTK_WINDOW(proj.window), 1024, 768);
	gtk_paned_set_position(GTK_PANED(hpaned), 172);
	gtk_paned_set_position(GTK_PANED(vpaned), 600);
	gtk_widget_show_all(proj.window);
	gtk_main();

	return 0;
}


static void infobar_response_sig(GtkInfoBar *infobar, gint resp,
                                 gpointer data)
{
	struct gui_task *task = data;

	if ( resp == GTK_RESPONSE_CANCEL ) {
		task->backend->cancel_task(task->job_priv);

	} else if ( resp == GTK_RESPONSE_CLOSE ) {

		gtk_widget_destroy(GTK_WIDGET(infobar));
		/* FIXME: Remove task from list */

	} else {
		ERROR("Unrecognised infobar response!\n");
	}
}


static gboolean update_info_bar(gpointer data)
{
	struct gui_task *task = data;
	int running;
	float frac_complete;

	if ( task->backend->task_status(task->job_priv,
	                                &running, &frac_complete) ) {
		ERROR("Error retrieving task status\n");
		return G_SOURCE_CONTINUE;
	}

	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(task->progress_bar),
	                              frac_complete);

	if ( !running && task->running ) {
		/* Task is no longer running */
		task->running = 0;
		gtk_widget_destroy(task->cancel_button);
		gtk_info_bar_set_show_close_button(GTK_INFO_BAR(task->info_bar),
		                                   TRUE);
		return G_SOURCE_REMOVE;
	}

	return G_SOURCE_CONTINUE;
}


void add_running_task(struct crystfelproject *proj,
                      const char *task_desc,
                      struct crystfel_backend *backend,
                      void *job_priv)
{
	struct gui_task *task;
	GtkWidget *bar_area;

	task = &proj->tasks[proj->n_running_tasks++];
	task->job_priv = job_priv;
	task->backend = backend;
	task->running = 1;

	/* Progress info bar */
	task->info_bar = gtk_info_bar_new();
	gtk_info_bar_set_message_type(GTK_INFO_BAR(task->info_bar),
	                              GTK_MESSAGE_INFO);

	task->cancel_button = gtk_info_bar_add_button(GTK_INFO_BAR(task->info_bar),
	                                              GTK_STOCK_CANCEL,
	                                              GTK_RESPONSE_CANCEL);

	gtk_box_pack_end(GTK_BOX(proj->main_vbox), GTK_WIDGET(task->info_bar),
	                 FALSE, FALSE, 0.0);

	bar_area = gtk_info_bar_get_content_area(GTK_INFO_BAR(task->info_bar));

	/* Create progress bar */
	task->progress_bar = gtk_progress_bar_new();
	gtk_box_pack_start(GTK_BOX(bar_area),
	                   GTK_WIDGET(task->progress_bar),
	                   TRUE, TRUE, 0.0);
	gtk_progress_bar_set_text(GTK_PROGRESS_BAR(task->progress_bar),
	                          task_desc);
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(task->progress_bar),
	                               TRUE);

	g_signal_connect(G_OBJECT(task->info_bar), "response",
	                 G_CALLBACK(infobar_response_sig), task);

	gtk_widget_show_all(task->info_bar);

#if GTK_CHECK_VERSION(3,22,29)
	gtk_info_bar_set_revealed(GTK_INFO_BAR(task->info_bar), TRUE);
#endif

	g_timeout_add(500, update_info_bar, task);
}


static GFile *get_crystfel_path_gfile()
{
	GFile *self;
	GFileInfo *self_info;
	const char *self_target;
	GFile *tar;
	GFile *parent_dir;
	GError *error = NULL;

	self = g_file_new_for_path("/proc/self/exe");
	self_info = g_file_query_info(self, "standard",
	                              G_FILE_QUERY_INFO_NONE,
	                              NULL, &error);
	if ( self_info == NULL ) return NULL;

	self_target = g_file_info_get_symlink_target(self_info);
	if ( self_target == NULL ) return NULL;

	tar = g_file_new_for_path(self_target);
	if ( tar == NULL ) return NULL;

	parent_dir = g_file_get_parent(tar);
	if ( parent_dir == NULL ) return NULL;

	g_object_unref(self);
	g_object_unref(self_info);
	g_object_unref(tar);

	return parent_dir;
}


char *get_crystfel_path_str()
{
	char *path;
	GFile *crystfel_path = get_crystfel_path_gfile();
	if ( crystfel_path == NULL ) return NULL;
	path = g_file_get_path(crystfel_path);
	g_object_unref(crystfel_path);
	return path;
}


char *get_crystfel_exe(const char *program)
{
	GFile *crystfel_path;
	char *exe_path;
	GFile *exe;

	crystfel_path = get_crystfel_path_gfile();
	if ( crystfel_path == NULL ) return NULL;

	exe = g_file_get_child(crystfel_path, program);
	if ( exe == NULL ) {
		ERROR("Couldn't determine executable path. "
		      "This is OK provided the executable "
		      "path is set correctly.\n");
		exe_path = strdup(program);
	}

	exe_path = g_file_get_path(exe);
	g_object_unref(exe);
	g_object_unref(crystfel_path);

	return exe_path;
}


struct gui_job_notes_page *add_job_notes_page(GtkWidget *notebook)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *scroll;
	struct gui_job_notes_page *notes;

	notes = malloc(sizeof(struct gui_job_notes_page));
	if ( notes == NULL ) return NULL;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	label = gtk_label_new("Whatever you enter here will be placed in "
	                      "the job's folder as 'notes.txt'");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   TRUE, TRUE, 0);
	notes->textview = gtk_text_view_new();
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(notes->textview),
	                            GTK_WRAP_WORD_CHAR);
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(scroll), notes->textview);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll),
	                                    GTK_SHADOW_ETCHED_IN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(scroll),
	                   TRUE, TRUE, 2.0);

	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box,
	                         gtk_label_new("Notes"));
	return notes;
}


GFile *make_job_folder(const char *job_title, const char *job_notes)
{
	struct stat s;
	char *workdir;
	GFile *workdir_file;
	GFile *cwd_file;
	GFile *notes_file;
	char *notes_path;
	FILE *fh;

	workdir = strdup(job_title);
	if ( workdir == NULL ) return NULL;

	if ( stat(workdir, &s) != -1 ) {
		ERROR("Working directory already exists.  "
		      "Choose a different job name.\n");
		return NULL;
	}

	if ( mkdir(workdir, S_IRWXU) ) {
		ERROR("Failed to create working directory: %s\n",
		      strerror(errno));
		return NULL;
	}

	cwd_file = g_file_new_for_path(".");
	workdir_file = g_file_get_child(cwd_file, workdir);
	g_object_unref(cwd_file);

	/* Write the notes into notes.txt */
	notes_file = g_file_get_child(workdir_file, "notes.txt");
	notes_path = g_file_get_path(notes_file);
	fh = fopen(notes_path, "w");
	fputs(job_notes, fh);
	fclose(fh);
	g_free(notes_path);
	g_object_unref(notes_file);

	return workdir_file;
}
