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

#include "crystfelimageview.h"


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


struct crystfelproject {

	GtkWidget *window;
	GtkUIManager *ui;
	GtkActionGroup *action_group;

	GtkWidget *imageview;
	GtkWidget *icons;      /* Drawing area for task icons */
	GtkWidget *report;     /* Text view at the bottom for messages */

	int cur_frame;

	int n_frames;
	char **filenames;
	char **events;
};


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


static void finddata_response_sig(GtkWidget *dialog, gint resp,
                                  struct crystfelproject *proj)
{
	if ( (resp==GTK_RESPONSE_DELETE_EVENT) || (resp==GTK_RESPONSE_CANCEL) ) {
		gtk_widget_destroy(dialog);
		return;
	}
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

	dialog = gtk_dialog_new_with_buttons("Find data files",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Find data", GTK_RESPONSE_ACCEPT,
	                                     NULL);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Find data in folder:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	chooser = gtk_file_chooser_button_new("Select a folder",
	                                      GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(chooser), TRUE, TRUE, 2.0);

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
	                "LCLS, individual files from Cheetah ('LCLS*.h5')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "lcls-cheetah-cxi",
	                "Multi-event CXI files from Cheetah ('*.cxi')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "cbf",
	                "Individual CBF files ('*.cbf')");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "cbfgz",
	                "Individual gzipped CBF files ('*.cbf.gz')");
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Geometry file:");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	chooser = gtk_file_chooser_button_new("Select geometry file",
	                                      GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(chooser), TRUE, TRUE, 2.0);

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


static void update_imageview(struct crystfelproject *proj)
{
	if ( proj->n_frames == 0 ) return;
	crystfel_image_view_set_image(CRYSTFEL_IMAGE_VIEW(proj->imageview),
	                              proj->filenames[proj->cur_frame],
	                              proj->events[proj->cur_frame]);
}


static gint first_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	proj->cur_frame = 0;
	update_imageview(proj);
	return FALSE;
}


static gint prev_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->cur_frame == 0 ) return FALSE;
	proj->cur_frame--;
	update_imageview(proj);
	return FALSE;
}


static gint next_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	if ( proj->cur_frame == proj->n_frames - 1 ) return FALSE;
	proj->cur_frame++;
	update_imageview(proj);
	return FALSE;
}


static gint last_frame_sig(GtkWidget *widget, struct crystfelproject *proj)
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


static void add_menu_bar(struct crystfelproject *proj, GtkWidget *vbox)
{
	GError *error = NULL;

	const char *ui = "<ui> <menubar name=\"cellwindow\">"
		"<menu name=\"file\" action=\"FileAction\">"
		"	<menuitem name=\"quit\" action=\"QuitAction\" />"
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

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About", NULL, NULL,
			G_CALLBACK(about_sig) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	proj->action_group = gtk_action_group_new("cellwindow");
	gtk_action_group_add_actions(proj->action_group, entries, n_entries, proj);

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


static void load_image_list(struct crystfelproject *proj, const char *filename)
{
	int done = 0;
	int max_frames = 0;
	FILE *fh = fopen(filename, "r");

	proj->filenames = NULL;
	proj->events = NULL;
	proj->n_frames = 0;

	do {
		char line[1024];
		if ( fgets(line, 1024, fh) != NULL ) {

			if ( proj->n_frames == max_frames ) {
				proj->filenames = realloc(proj->filenames,
				                          (max_frames+1024)*sizeof(char *));
				proj->events = realloc(proj->events,
				                       (max_frames+1024)*sizeof(char *));
				if ( (proj->filenames == NULL)
				  || (proj->events == NULL) )
				{
					ERROR("Failed to allocate while "
					      "loading image list\n");
					proj->n_frames = 0;
					return;
				}
				max_frames += 1024;
			}

			chomp(line);
			proj->filenames[proj->n_frames] = strdup(line);
			proj->events[proj->n_frames] = NULL;
			proj->n_frames++;

		} else {
			done = 1;
		}
	} while ( !done );
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
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Determine unit cell", "document-page-setup",
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Index and integrate", "system-run",
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Merge", "applications-science",
	           G_CALLBACK(NULL), proj);
	add_button(vbox, "Figures of merit", "trophy-gold",
	           G_CALLBACK(NULL), proj);
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
	DataTemplate *dtempl;

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

	load_image_list(&proj, "files.lst");

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

	/* FIXME: Testing stuff */
	dtempl = data_template_new_from_file("5HT2b-Liu-2013.geom");
	crystfel_image_view_set_datatemplate(CRYSTFEL_IMAGE_VIEW(proj.imageview),
	                                     dtempl);

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
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(proj.report));
	gtk_paned_pack2(GTK_PANED(vpaned), GTK_WIDGET(frame), FALSE, FALSE);

	gtk_window_set_default_size(GTK_WINDOW(proj.window), 1024, 768);
	gtk_paned_set_position(GTK_PANED(hpaned), 172);
	gtk_paned_set_position(GTK_PANED(vpaned), 600);
	gtk_widget_show_all(proj.window);
	gtk_main();

	return 0;
}
