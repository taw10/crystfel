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
#include "gui_project.h"


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


/* Bring the image view up to date after changing the selected image */
static void update_imageview(struct crystfelproject *proj)
{
	char tmp[1024];
	char *ev_str;
	char *ev_sep;
	struct image *image;

	if ( proj->n_frames == 0 ) return;

	if ( proj->stream != NULL ) {

		if ( stream_select_chunk(proj->stream,
		                         proj->cur_frame) )
		{
			ERROR("Failed to select new chunk\n");
			return;
		}

		image = stream_read_chunk(proj->stream,
		                          proj->dtempl,
		                          STREAM_UNITCELL
		                        | STREAM_REFLECTIONS
		                        | STREAM_PEAKS
		                        | STREAM_CRYSTALS);

		if ( image == NULL ) {
			ERROR("Failed to read from stream\n");
			return;
		}

	} else {

		image = image_read(proj->dtempl,
		                   proj->filenames[proj->cur_frame],
		                   proj->events[proj->cur_frame]);

		if ( image == NULL ) {
			ERROR("Failed to load image\n");
			return;
		}

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
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              proj->cur_image);
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


struct finddata_ctx
{
	struct crystfelproject *proj;
	GtkWidget *top_folder;
	GtkWidget *geom_file;
	GtkWidget *stream;
	GtkWidget *type_combo;
	GtkWidget *fi;  /* "Find files" radio */
};


static void finddata_response_sig(GtkWidget *dialog, gint resp,
                                  struct finddata_ctx *ctx)
{
	struct crystfelproject *proj = ctx->proj;

	if ( (resp==GTK_RESPONSE_DELETE_EVENT) || (resp==GTK_RESPONSE_CANCEL) ) {
		gtk_widget_destroy(dialog);
		free(ctx);
		return;
	}

	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ctx->fi)) ) {

		GFile *top;
		DataTemplate *dtempl;
		char *geom_filename;
		const char *type_id;

		geom_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->geom_file));
		if ( geom_filename == NULL ) return;

		top = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(ctx->top_folder));
		if ( top == NULL ) return;

		dtempl = data_template_new_from_file(geom_filename);
		if ( dtempl == NULL ) return;

		type_id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(ctx->type_combo));
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

	} else {

		Stream *st;
		char *stream_filename;
		DataTemplate *dtempl;
		const char *geom_str;
		int n_chunks;

		stream_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(ctx->stream));
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

		stream_close(proj->stream);
		proj->stream = st;

		data_template_free(proj->dtempl);
		proj->dtempl = dtempl;

		/* Set some defaults for things we won't be using */
		g_free(proj->geom_filename);
		proj->geom_filename = NULL;
		g_free(proj->data_top_folder);
		proj->data_top_folder = NULL;
		proj->data_search_pattern = MATCH_EVERYTHING;

		n_chunks = stream_scan_chunks(st);

		if ( n_chunks == 0 ) {
			ERROR("No chunks found (or error reading)\n");
			stream_close(st);
			return;
		}

		proj->n_frames = n_chunks;
		proj->stream_filename = stream_filename;

		crystfel_image_view_set_show_peaks(CRYSTFEL_IMAGE_VIEW(proj->imageview),
		                                   1);

	}

	proj->unsaved = 1;
	proj->cur_frame = 0;
	crystfel_image_view_reset_zoom(CRYSTFEL_IMAGE_VIEW(proj->imageview));
	update_imageview(proj);
	update_peaks(proj);

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
	GtkWidget *ls;

	ctx = malloc(sizeof(struct finddata_ctx));
	if ( ctx == NULL ) return FALSE;

	ctx->proj = proj;

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

	ctx->fi = gtk_radio_button_new_with_label(NULL, "Load images directly");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(ctx->fi),
	                   FALSE, FALSE, 8.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Find data in folder:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	ctx->top_folder = gtk_file_chooser_button_new("Select a folder",
	                                              GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
	if ( proj->data_top_folder != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->top_folder),
		                              proj->data_top_folder);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->top_folder), TRUE, TRUE, 2.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Search pattern:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	ctx->type_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->type_combo), TRUE, TRUE, 2.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->type_combo), "everything",
	                "All files in folder and subfolders");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->type_combo), "lcls-cheetah-hdf5",
	                "Individual LCLS files from Cheetah ('LCLS*.h5')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->type_combo), "cheetah-cxi",
	                "Multi-event CXI files from Cheetah ('*.cxi')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->type_combo), "cbf",
	                "Individual CBF files ('*.cbf')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(ctx->type_combo), "cbfgz",
	                "Individual gzipped CBF files ('*.cbf.gz')");
	gtk_combo_box_set_active(GTK_COMBO_BOX(ctx->type_combo),
	                         proj->data_search_pattern);

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

	ls = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(ctx->fi),
	                                                 "Load stream");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(ls),
	                   FALSE, FALSE, 8.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Stream file:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	ctx->stream = gtk_file_chooser_button_new("Select stream file",
	                                          GTK_FILE_CHOOSER_ACTION_OPEN);
	if ( proj->stream_filename != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(ctx->stream),
		                              proj->stream_filename);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(ctx->stream), TRUE, TRUE, 2.0);

	g_signal_connect(dialog, "response",
	                 G_CALLBACK(finddata_response_sig), ctx);

	gtk_window_set_default_size(GTK_WINDOW(dialog), 512, 0);
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
		"	<menuitem name=\"save\" action=\"SaveAction\" />"
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
		{ "SaveAction", GTK_STOCK_SAVE, "_Save", NULL, NULL,
			G_CALLBACK(save_sig) },
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
	add_button(vbox, "Load data", "folder-pictures",
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
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_view_scroll_to_iter(GTK_TEXT_VIEW(proj->report),
	                             &iter, 0.0, FALSE, 0.0, 0.0);
}


static void brightness_changed_sig(GtkScaleButton *brightness,
                                   double value,
                                   struct crystfelproject *proj)
{
	crystfel_image_view_set_brightness(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                                   value);
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
	GtkWidget *brightness;

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

	default_project(&proj);

	proj.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(proj.window), "CrystFEL");
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

	/* Image view parameters */
	const gchar *icons[] = {"weather-clear", NULL};
	brightness = gtk_scale_button_new(GTK_ICON_SIZE_LARGE_TOOLBAR,
	                                  1.0, 10.0, 1.0, icons);
	gtk_box_pack_end(GTK_BOX(toolbar), brightness, FALSE, FALSE, 0.0);
	gtk_scale_button_set_value(GTK_SCALE_BUTTON(brightness), 1.0);
	g_signal_connect(G_OBJECT(brightness), "value-changed",
	                 G_CALLBACK(brightness_changed_sig), &proj);

	/* Filename */
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
		}

		w = gtk_ui_manager_get_action(proj.ui, "/mainwindow/view/peaks");
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(w),
		                             proj.show_peaks);
		update_imageview(&proj);
		update_peaks(&proj);
	}

	/* Initialise backend */
	proj.backend->init(&proj);

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
	struct crystfelproject *proj = data;

	if ( resp == GTK_RESPONSE_CANCEL ) {
		proj->backend->cancel(proj);

	} else if ( resp == GTK_RESPONSE_OK ) {
		proj->infobar_callback(proj);

	} else {
		ERROR("Unrecognised infobar response!\n");

	}
}


void remove_infobar(struct crystfelproject *proj)
{
	gtk_widget_destroy(proj->info_bar);
	proj->info_bar = NULL;
}


GtkWidget *create_infobar(struct crystfelproject *proj, const char *task,
                          const char *extra_button,
                          void (*cbfunc)(struct crystfelproject *proj))
{
	GtkWidget *info_bar;
	GtkWidget *bar_area;

	if ( proj->info_bar != NULL ) {
		STATUS("Can't create info bar - task already running\n");
		return NULL;
	}

	/* Progress info bar */
	info_bar = gtk_info_bar_new_with_buttons(GTK_STOCK_CANCEL,
	                                         GTK_RESPONSE_CANCEL,
	                                         extra_button,
	                                         GTK_RESPONSE_OK,
	                                         NULL);
	gtk_box_pack_end(GTK_BOX(proj->main_vbox), GTK_WIDGET(info_bar),
	                 FALSE, FALSE, 0.0);
	proj->info_bar = info_bar;
	proj->infobar_callback = cbfunc;

	bar_area = gtk_info_bar_get_content_area(GTK_INFO_BAR(info_bar));

	/* Create progress bar */
	proj->progressbar = gtk_progress_bar_new();
	gtk_box_pack_start(GTK_BOX(bar_area),
	                   GTK_WIDGET(proj->progressbar),
	                   TRUE, TRUE, 0.0);
	gtk_progress_bar_set_text(GTK_PROGRESS_BAR(proj->progressbar),
	                          task);
	gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(proj->progressbar),
	                               TRUE);

	g_signal_connect(G_OBJECT(info_bar), "response",
	                 G_CALLBACK(infobar_response_sig), proj);

	gtk_widget_show_all(info_bar);
	gtk_info_bar_set_revealed(GTK_INFO_BAR(info_bar), TRUE);

	return info_bar;
}
