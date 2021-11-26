/*
 * gui_import.c
 *
 * Data import to GUI
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <datatemplate.h>

#include "crystfelimageview.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"

#include "version.h"

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

	GtkWidget *replace_geom;
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
	GtkWidget *dump_results;
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
	const char *type_id;
	struct crystfelproject *proj = ctx->proj;

	top = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(ctx->search_chooser));
	if ( top == NULL ) return;

	type_id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(ctx->search_pattern));
	proj->data_search_pattern = decode_matchtype(type_id);

	g_free(proj->data_top_folder);
	proj->data_top_folder = g_file_get_path(top);

	add_files(proj, top, proj->data_search_pattern, proj->dtempl);

	g_object_unref(top);
}


/* stream_filename will be adopted */
int load_stream(struct crystfelproject *proj, char *stream_filename)
{
	Stream *st;
	DataTemplate *dtempl;
	const char *geom_str;
	char *result_name;

	st = stream_open_for_read(stream_filename);
	if ( st == NULL ) return 1;

	geom_str = stream_geometry_file(st);
	if ( geom_str == NULL ) {
		ERROR("No geometry file\n");
		stream_close(st);
		return 1;
	}

	dtempl = data_template_new_from_string(geom_str);
	if ( dtempl == NULL ) {
		stream_close(st);
		return 1;
	}

	/* If we do not yet have a DataTemplate, the one from the file
	 * becomes it.  If we already have one, it will be kept.  Note that the
	 * stream's DataTemplate will always be used for display in the GUI. */
	if ( proj->dtempl == NULL ) {
		proj->dtempl = dtempl;
	}

	/* Use the user's nominated DataTemplate over the one from the stream.
	 * If it doesn't match, better that things break earlier. */
	add_frames_from_stream(st, proj->dtempl, proj);
	stream_close(st);

	result_name = safe_basename(stream_filename);
	add_indexing_result(proj, result_name, &stream_filename, 1);
	select_result(proj, result_name);
	free(result_name);

	return 0;
}


static void import_stream(struct finddata_ctx *ctx)
{
	char *stream_filename;

	stream_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->stream_chooser));
	if ( stream_filename == NULL ) return;

	load_stream(ctx->proj, stream_filename);
}


static void import_file_list(struct finddata_ctx *ctx)
{
	struct crystfelproject *proj = ctx->proj;
	char *list_filename;
	FILE *fh;

	list_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->list_chooser));
	if ( list_filename == NULL ) return;

	fh = fopen(list_filename, "r");
	if ( fh == NULL ) return;

	do {

		char line[1024];
		char *event = "//";
		size_t n;

		if ( fgets(line, 1024, fh) == NULL ) break;
		chomp(line);

		/* Chop off event ID */
		n = strlen(line);
		while ( line[n] != ' ' && n > 2 ) n--;
		if ( n != 2 ) {
			/* Event descriptor must contain "//".
			 * If it doesn't, assume the filename just contains a
			 * space. */
			if ( strstr(&line[n], "//") != NULL ) {
				line[n] = '\0';
				event = &line[n+1];
			}
		} /* else no spaces at all */

		if ( event != NULL ) {
			/* Explicit event ID given */
			add_file_to_project(proj, line, event);
		} else {
			/* No event ID - expand (possibly 1:1) */
			add_all_events(proj, line, proj->dtempl);
		}

	} while ( 1 );

	fclose(fh);
}


static void import_file(struct finddata_ctx *ctx)
{
	struct crystfelproject *proj = ctx->proj;
	char *filename;

	filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->indiv_chooser));
	if ( filename == NULL ) return;

	add_all_events(proj, filename, proj->dtempl);
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

	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->dump)) ) {
		clear_project_files(proj);
		crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
		                              NULL);

		if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->dump_results)) ) {
			clear_indexing_results(proj);
		}
	}

	if ( import_mode(ctx) != IMPORT_STREAM ) {
		if ( (ctx->replace_geom == NULL)
		  || (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->replace_geom))) )
		{
			gchar *geom_filename;

			geom_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->geom_file));
			if ( geom_filename == NULL ) {
				error_box(proj, "Geometry file not found");
				return;
			}
			g_free(proj->geom_filename);
			proj->geom_filename = geom_filename;

			data_template_free(proj->dtempl);
			proj->dtempl = data_template_new_from_file(geom_filename);
			if ( proj->dtempl == NULL ) {
				error_box(proj, "Invalid geometry file");
				return;
			}
		}
	} /* else don't touch the geometry */

	if ( (import_mode(ctx) != IMPORT_STREAM) && (proj->dtempl == NULL) ) {
		error_box(proj, "You must specify the geometry file.");
		return;
	}

	switch ( import_mode(ctx) ) {

		case IMPORT_FILES :
		import_file(ctx);
		break;

		case IMPORT_LIST :
		import_file_list(ctx);
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


gint import_sig(GtkWidget *widget, struct crystfelproject *proj)
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

	/* Select individual files */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->indiv = gtk_radio_button_new_with_label(NULL,
	                                             "Select an individual file");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->indiv),
	                   FALSE, FALSE, 4.0);
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
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->stream_chooser),
	                   TRUE, TRUE, 2.0);

	/* Stuff at bottom */
	gtk_box_pack_start(GTK_BOX(vbox),
	                   gtk_separator_new(GTK_ORIENTATION_HORIZONTAL),
	                   FALSE, FALSE, 4.0);

	/* Geometry file */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);

	ctx->geom_file = gtk_file_chooser_button_new("Select geometry file",
	                                             GTK_FILE_CHOOSER_ACTION_OPEN);
	if ( proj->geom_filename != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->geom_file),
		                              proj->geom_filename);
	}
	if ( proj->dtempl == NULL ) {
		label = gtk_label_new("Geometry file:");
		gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
		gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
		                   FALSE, FALSE, 4.0);
		ctx->replace_geom = NULL;
	} else {
		ctx->replace_geom = gtk_check_button_new_with_label("Replace geometry file:");
		gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->replace_geom),
		                   FALSE, FALSE, 4.0);
		g_signal_connect(G_OBJECT(ctx->replace_geom), "toggled",
		                 G_CALLBACK(i_maybe_disable), ctx->geom_file);
		gtk_widget_set_sensitive(ctx->geom_file, FALSE);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->geom_file), TRUE, TRUE, 2.0);

	/* Replace data toggle */
	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	ctx->dump = gtk_check_button_new_with_label("Replace all the current data");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->dump),
	                   FALSE, FALSE, 4.0);
	ctx->dump_results = gtk_check_button_new_with_label("Forget about indexing results");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->dump_results),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(ctx->dump), "toggled",
	                 G_CALLBACK(i_maybe_disable_and_deselect),
	                 ctx->dump_results);
	gtk_widget_set_sensitive(ctx->dump_results, FALSE);

	g_signal_connect(ctx->indiv, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	g_signal_connect(ctx->list, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	g_signal_connect(ctx->search, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);
	g_signal_connect(ctx->stream, "toggled",
	                 G_CALLBACK(finddata_typetoggle_sig), ctx);

	g_signal_connect(dialog, "response",
	                 G_CALLBACK(finddata_response_sig), ctx);

	gtk_window_set_default_size(GTK_WINDOW(dialog), 512, 0);
	finddata_typetoggle_sig(ctx->search, ctx);
	gtk_widget_show_all(dialog);
	return FALSE;
}


