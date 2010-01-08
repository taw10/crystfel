/*
 * displaywindow.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _GNU_SOURCE
#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cairo.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

#include "displaywindow.h"
#include "render.h"
#include "hdf5-file.h"
#include "hdfsee.h"


#define INITIAL_BINNING 2


static void displaywindow_error(DisplayWindow *dw, const char *message)
{
	GtkWidget *window;

	window = gtk_message_dialog_new(GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_MESSAGE_WARNING,
					GTK_BUTTONS_CLOSE, message);

	g_signal_connect_swapped(window, "response",
				 G_CALLBACK(gtk_widget_destroy), window);
	gtk_widget_show(window);
}


static void displaywindow_update(DisplayWindow *dw)
{
	gint width;
	GdkGeometry geom;

	if ( dw->hdfile != NULL ) {
		dw->width = hdfile_get_width(dw->hdfile)/dw->binning;
		dw->height = hdfile_get_height(dw->hdfile)/dw->binning;
	} else {
		dw->width = 320;
		dw->height = 320;
	}

	width = dw->width;
	if ( dw->show_col_scale ) width += 20;

	gtk_widget_set_size_request(GTK_WIDGET(dw->drawingarea), width,
				    dw->height);
	geom.min_width = -1;
	geom.min_height = -1;
	geom.max_width = -1;
	geom.max_height = -1;
	gtk_window_set_geometry_hints(GTK_WINDOW(dw->window),
				      GTK_WIDGET(dw->drawingarea), &geom,
				      GDK_HINT_MIN_SIZE | GDK_HINT_MAX_SIZE);

	if ( dw->pixbuf != NULL ) {
		gdk_pixbuf_unref(dw->pixbuf);
	}
	if ( dw->hdfile != NULL ) {
		dw->pixbuf = render_get_image(dw->hdfile, dw->binning,
					      dw->boostint, dw->monochrome);
	} else {
		dw->pixbuf = NULL;
	}

	if ( dw->col_scale != NULL ) {
		gdk_pixbuf_unref(dw->col_scale);
	}
	dw->col_scale = render_get_colour_scale(20, dw->height, dw->monochrome);

	gdk_window_invalidate_rect(dw->drawingarea->window, NULL, FALSE);
}


/* Window closed - clean up */
static gint displaywindow_closed(GtkWidget *window, DisplayWindow *dw)
{
	if ( dw->hdfile != NULL ) {
		hdfile_close(dw->hdfile);
	}

	/* Notify 'main', so it can update the master list */
	hdfsee_window_closed(dw);

	return 0;
}


static gboolean displaywindow_expose(GtkWidget *da, GdkEventExpose *event,
				     DisplayWindow *dw)
{
	cairo_t *cr;

	cr = gdk_cairo_create(da->window);

	/* Blank white background */
	cairo_rectangle(cr, 0.0, 0.0, da->allocation.width,
			da->allocation.height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
	cairo_fill(cr);

	cairo_destroy(cr);

	if ( dw->pixbuf != NULL ) {
		gdk_draw_pixbuf(da->window,
				da->style->bg_gc[GTK_WIDGET_STATE(da)],
				dw->pixbuf,
				0, 0, 0, 0, dw->width, dw->height,
				GDK_RGB_DITHER_NONE, 0, 0);
	}

	if ( (dw->show_col_scale) && (dw->col_scale != NULL) ) {
		gdk_draw_pixbuf(da->window,
				da->style->bg_gc[GTK_WIDGET_STATE(da)],
				dw->col_scale,
				0, 0, dw->width, 0, 20, dw->height,
				GDK_RGB_DITHER_NONE, 0, 0);
	}

	return FALSE;
}


static gint displaywindow_close(GtkWidget *widget, DisplayWindow *dw)
{
	gtk_widget_destroy(dw->window);
	return 0;
}


static gint displaywindow_set_binning_response(GtkWidget *widget, gint response,
					       DisplayWindow *dw)
{
	int done = 1;

	if ( response == GTK_RESPONSE_OK ) {

		const char *sbinning;
		unsigned int binning;
		int scanval;

		sbinning = gtk_entry_get_text(
					GTK_ENTRY(dw->binning_dialog->entry));
		scanval = sscanf(sbinning, "%u", &binning);
		if ( (scanval != 1) || (binning <= 0) ) {
			displaywindow_error(dw,
				"Please enter a positive integer for the "
				"binning factor.");
			done = 0;
		} else {
			if ((binning < hdfile_get_width(dw->hdfile)/10)
			 && (binning < hdfile_get_height(dw->hdfile)/10)) {
				dw->binning = binning;
				displaywindow_update(dw);
			} else {
				displaywindow_error(dw,
					"Please enter a sensible value for "
					"the binning factor.");
				done = 0;
			}
		}
	}

	if ( done ) {
		gtk_widget_destroy(dw->binning_dialog->window);
	}

	return 0;

}


static gint displaywindow_set_binning_destroy(GtkWidget *widget,
					      DisplayWindow *dw)
{
	free(dw->binning_dialog);
	dw->binning_dialog = NULL;
	return 0;
}


static gint displaywindow_set_binning_response_ac(GtkWidget *widget,
						  DisplayWindow *dw)
{
	return displaywindow_set_binning_response(widget, GTK_RESPONSE_OK, dw);
}


/* Create a window to ask the user for a new binning factor */
static gint displaywindow_set_binning(GtkWidget *widget, DisplayWindow *dw)
{
	BinningDialog *bd;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *table;
	GtkWidget *label;
	char tmp[64];

	if ( dw->binning_dialog != NULL ) {
		return 0;
	}

	if ( dw->hdfile == NULL ) {
		return 0;
	}

	bd = malloc(sizeof(BinningDialog));
	if ( bd == NULL ) return 0;
	dw->binning_dialog = bd;

	bd->window = gtk_dialog_new_with_buttons("Set Binning",
					GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_STOCK_CANCEL, GTK_RESPONSE_CLOSE,
					GTK_STOCK_OK, GTK_RESPONSE_OK,
					NULL);

	vbox = gtk_vbox_new(FALSE, 0);
	hbox = gtk_hbox_new(TRUE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(bd->window)->vbox),
			   GTK_WIDGET(hbox), FALSE, FALSE, 7);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(vbox), FALSE, FALSE, 5);

	table = gtk_table_new(3, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(table), FALSE, FALSE, 0);

	label = gtk_label_new("Smaller numbers mean larger images on screen");
	gtk_label_set_markup(GTK_LABEL(label),
			"<span style=\"italic\" weight=\"light\">"
			"Smaller numbers mean larger images on screen</span>");
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 3, 1, 2);

	snprintf(tmp, 63, "Raw image size: %i by %i pixels",
		 (int)hdfile_get_width(dw->hdfile),
		 (int)hdfile_get_height(dw->hdfile));
	label = gtk_label_new(tmp);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 3, 2, 3);

	label = gtk_label_new("Binning Factor:");
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 2, 3, 4);

	bd->entry = gtk_entry_new();
	snprintf(tmp, 63, "%i", dw->binning);
	gtk_entry_set_text(GTK_ENTRY(bd->entry), tmp);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(bd->entry),
				  2, 3, 3, 4);

	g_signal_connect(G_OBJECT(bd->entry), "activate",
			 G_CALLBACK(displaywindow_set_binning_response_ac), dw);
	g_signal_connect(G_OBJECT(bd->window), "response",
			 G_CALLBACK(displaywindow_set_binning_response), dw);
	g_signal_connect(G_OBJECT(bd->window), "destroy",
			 G_CALLBACK(displaywindow_set_binning_destroy), dw);
	gtk_window_set_resizable(GTK_WINDOW(bd->window), FALSE);
	gtk_widget_show_all(bd->window);
	gtk_widget_grab_focus(GTK_WIDGET(bd->entry));

	return 0;
}


static gint displaywindow_set_boostint_response(GtkWidget *widget,
						gint response,
						DisplayWindow *dw)
{
	int done = 1;

	if ( response == GTK_RESPONSE_OK ) {

		const char *sboostint;
		unsigned int boostint;
		int scanval;

		sboostint = gtk_entry_get_text(
		 			GTK_ENTRY(dw->boostint_dialog->entry));
		scanval = sscanf(sboostint, "%u", &boostint);
		if ( (scanval != 1) || (boostint <= 0) ) {
			displaywindow_error(dw, "Please enter a positive "
					"integer for the intensity boost "
					"factor.");
			done = 0;
		} else {
			dw->boostint = boostint;
			displaywindow_update(dw);
		}
	}

	if ( done ) {
		gtk_widget_destroy(dw->boostint_dialog->window);
	}

	return 0;
}


static gint displaywindow_set_boostint_destroy(GtkWidget *widget,
					       DisplayWindow *dw)
{
	free(dw->boostint_dialog);
	dw->boostint_dialog = NULL;
	return 0;
}


static gint displaywindow_set_boostint_response_ac(GtkWidget *widget,
						   DisplayWindow *dw)
{
	return displaywindow_set_boostint_response(widget, GTK_RESPONSE_OK, dw);
}


/* Create a window to ask the user for a new intensity boost factor */
static gint displaywindow_set_boostint(GtkWidget *widget, DisplayWindow *dw)
{
	BoostIntDialog *bd;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *table;
	GtkWidget *label;
	char tmp[64];

	if ( dw->boostint_dialog != NULL ) {
		return 0;
	}

	if ( dw->hdfile == NULL ) {
		return 0;
	}

	bd = malloc(sizeof(BoostIntDialog));
	if ( bd == NULL ) return 0;
	dw->boostint_dialog = bd;

	bd->window = gtk_dialog_new_with_buttons("Intensity Boost",
					GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_STOCK_CANCEL, GTK_RESPONSE_CLOSE,
					GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);

	vbox = gtk_vbox_new(FALSE, 0);
	hbox = gtk_hbox_new(TRUE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(bd->window)->vbox),
			   GTK_WIDGET(hbox), FALSE, FALSE, 7);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(vbox), FALSE, FALSE, 5);

	table = gtk_table_new(3, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(table), FALSE, FALSE, 0);

	label = gtk_label_new("Boost Factor:");
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 2, 3, 4);

	bd->entry = gtk_entry_new();
	snprintf(tmp, 63, "%i", dw->boostint);
	gtk_entry_set_text(GTK_ENTRY(bd->entry), tmp);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(bd->entry),
				  2, 3, 3, 4);

	g_signal_connect(G_OBJECT(bd->entry), "activate",
			 G_CALLBACK(displaywindow_set_boostint_response_ac),
			 dw);
	g_signal_connect(G_OBJECT(bd->window), "response",
			 G_CALLBACK(displaywindow_set_boostint_response), dw);
	g_signal_connect(G_OBJECT(bd->window), "destroy",
			 G_CALLBACK(displaywindow_set_boostint_destroy), dw);
	gtk_window_set_resizable(GTK_WINDOW(bd->window), FALSE);
	gtk_widget_show_all(bd->window);
	gtk_widget_grab_focus(GTK_WIDGET(bd->entry));

	return 0;
}


static gint displaywindow_peaklist_response(GtkWidget *d, gint response,
                                            DisplayWindow *dw)
{
	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *filename;

		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		load_features_from_file(dw->image);

		g_free(filename);

	}

	gtk_widget_destroy(d);

	return 0;
}


static gint displaywindow_about(GtkWidget *widget, DisplayWindow *dw)
{
	GtkWidget *window;

	const gchar *authors[] = {
		"Thomas White <taw@physics.org>",
		"Erica Bithell <egb10@cam.ac.uk>",
		"Alex Eggeman <ase25@cam.ac.uk>",
		NULL
	};

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window),
				     GTK_WINDOW(dw->window));

	gtk_about_dialog_set_name(GTK_ABOUT_DIALOG(window), "hdfsee");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window), PACKAGE_VERSION);
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2010 Thomas White <taw@physics.org> and others");
	gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(window),
		"Quick viewer for HDF files");
	gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2009 Thomas White <taw@physics.org>\n");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"http://www.bitwiz.org.uk/");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static gint displaywindow_peak_overlay(GtkWidget *widget, DisplayWindow *dw)
{
	GtkWidget *d;

	d = gtk_file_chooser_dialog_new("Choose Peak List",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_OPEN,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
	                                NULL);

	g_signal_connect(G_OBJECT(d), "response",
	                 G_CALLBACK(displaywindow_peaklist_response), dw);

	gtk_widget_show_all(d);

	return 0;
}


static gint displaywindow_set_colscale(GtkWidget *widget, DisplayWindow *dw)
{
	dw->show_col_scale = 1 - dw->show_col_scale;
	displaywindow_update(dw);
	return 0;
}


static gint displaywindow_set_mono(GtkWidget *widget, DisplayWindow *dw)
{
	dw->monochrome = 1 - dw->monochrome;
	displaywindow_update(dw);
	return 0;
}


static gint displaywindow_numbers_response(GtkWidget *widget,
                                           gint response, DisplayWindow *dw)
{
	gtk_widget_destroy(dw->numbers_window->window);
	return 0;
}


static gint displaywindow_numbers_destroy(GtkWidget *widget, DisplayWindow *dw)
{
	free(dw->numbers_window);
	dw->numbers_window = NULL;
	return 0;
}


static gint displaywindow_show_numbers(GtkWidget *widget, DisplayWindow *dw)
{
	struct numberswindow *nw;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *table;
	unsigned int x, y;

	if ( dw->numbers_window != NULL ) {
		return 0;
	}

	if ( dw->hdfile == NULL ) {
		return 0;
	}

	nw = malloc(sizeof(struct numberswindow));
	if ( nw == NULL ) return 0;
	dw->numbers_window = nw;

	nw->window = gtk_dialog_new_with_buttons("Numbers",
					GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_STOCK_CLOSE, GTK_RESPONSE_CLOSE,
					NULL);

	vbox = gtk_vbox_new(FALSE, 0);
	hbox = gtk_hbox_new(TRUE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(nw->window)->vbox),
			   GTK_WIDGET(hbox), FALSE, FALSE, 7);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(vbox), FALSE, FALSE, 5);

	table = gtk_table_new(17, 17, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(table), FALSE, FALSE, 0);

	for ( x=0; x<17; x++ ) {
	for ( y=0; y<17; y++ ) {

		GtkWidget *label;

		label = gtk_label_new("--");
		gtk_widget_set_size_request(GTK_WIDGET(label), 40, -1);

		gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
		                          x, x+1, y, y+1);

		nw->labels[x+17*y] = label;

	}
	}

	g_signal_connect(G_OBJECT(nw->window), "response",
			 G_CALLBACK(displaywindow_numbers_response), dw);
	g_signal_connect(G_OBJECT(nw->window), "destroy",
			 G_CALLBACK(displaywindow_numbers_destroy), dw);
	gtk_window_set_resizable(GTK_WINDOW(nw->window), FALSE);

	gtk_widget_show_all(nw->window);

	return 0;
}


static void numbers_update(DisplayWindow *dw)
{
	int px, py;

	for ( px=0; px<17; px++ ) {
	for ( py=0; py<17; py++ ) {

		char s[32];
		int16_t val;
		GtkWidget *l;
		int x, y;

		x = dw->binning * dw->numbers_window->cx + (px-8);
		y = dw->binning * (dw->height-dw->numbers_window->cy) + (py-8);

		if ( (x>0) && (y>0) &&
		     !hdfile_get_unbinned_value(dw->hdfile, x, y, &val) ) {
			snprintf(s, 31, "%i", val);
		} else {
			strcpy(s, "--");
		}
		l = dw->numbers_window->labels[px+17*py];
		gtk_label_set_text(GTK_LABEL(l), s);

	}
	}
}


static void displaywindow_addui_callback(GtkUIManager *ui, GtkWidget *widget,
					 GtkContainer *container)
{
	gtk_box_pack_start(GTK_BOX(container), widget, FALSE, FALSE, 0);

	/* Enable overflow menu if this is a toolbar */
	if ( GTK_IS_TOOLBAR(widget) ) {
		gtk_toolbar_set_show_arrow(GTK_TOOLBAR(widget), TRUE);
	}
}


static void displaywindow_addmenubar(DisplayWindow *dw, GtkWidget *vbox)
{
	GError *error = NULL;
	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "ImagesAction", NULL, "Images", NULL, NULL, NULL },
		{ "CloseAction", GTK_STOCK_CLOSE, "_Close", NULL, NULL,
			G_CALLBACK(displaywindow_close) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },
		{ "BinningAction", NULL, "Set Binning...", "F3", NULL,
			G_CALLBACK(displaywindow_set_binning) },
		{ "BoostIntAction", NULL, "Boost Intensity...", "F5", NULL,
			G_CALLBACK(displaywindow_set_boostint) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "NumbersAction", NULL, "View Numbers...", "F2", NULL,
			G_CALLBACK(displaywindow_show_numbers) },
		{ "PeaksAction", NULL, "Peak Position Overlay...", NULL, NULL,
			G_CALLBACK(displaywindow_peak_overlay) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About hdfsee...",
			NULL, NULL,
			G_CALLBACK(displaywindow_about) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	GtkToggleActionEntry toggles[] = {
		{ "ColScaleAction", NULL, "Show Colour Scale", NULL, NULL,
			G_CALLBACK(displaywindow_set_colscale), FALSE },
		{ "MonoAction", NULL, "Monochrome", NULL, NULL,
			G_CALLBACK(displaywindow_set_mono), FALSE },
	};
	guint n_toggles = G_N_ELEMENTS(toggles);

	dw->action_group = gtk_action_group_new("hdfseedisplaywindow");
	gtk_action_group_add_actions(dw->action_group, entries, n_entries, dw);
	gtk_action_group_add_toggle_actions(dw->action_group, toggles,
					    n_toggles, dw);

	dw->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(dw->ui, dw->action_group, 0);
	g_signal_connect(dw->ui, "add_widget",
			 G_CALLBACK(displaywindow_addui_callback), vbox);
	if ( gtk_ui_manager_add_ui_from_file(dw->ui,
	     DATADIR"/hdfsee/displaywindow.ui", &error) == 0 ) {
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(dw->window),
				   gtk_ui_manager_get_accel_group(dw->ui));
	gtk_ui_manager_ensure_update(dw->ui);
}


struct newhdf {
	DisplayWindow *dw;
	char name[1024];
};

static gint displaywindow_newhdf(GtkMenuItem *item, struct newhdf *nh)
{
	hdfile_set_image(nh->dw->hdfile, nh->name);
	displaywindow_update(nh->dw);
	return 0;
}


static GtkWidget *displaywindow_addhdfgroup(struct hdfile *hdfile,
                                            const char *group,
                                            DisplayWindow *dw)
{
	char **names;
	int *is_group;
	int *is_image;
	GtkWidget *ms;
	GSList *rg = NULL;
	int n, i;

	names = hdfile_read_group(hdfile, &n, group, &is_group, &is_image);
	if ( n == 0 ) return NULL;

	ms = gtk_menu_new();

	for ( i=0; i<n; i++ ) {

		GtkWidget *item;
		GtkWidget *sub;

		if ( names[i] == NULL ) return NULL;

		if ( is_group[i] ) {

			item = gtk_menu_item_new_with_label(names[i]);

			sub = displaywindow_addhdfgroup(hdfile, names[i], dw);
			gtk_menu_item_set_submenu(GTK_MENU_ITEM(item), sub);

		} else if ( is_image[i] ) {

			struct newhdf *nh;

			item = gtk_radio_menu_item_new_with_label(rg, names[i]);
			rg = gtk_radio_menu_item_get_group(
		                                     GTK_RADIO_MENU_ITEM(item));

			nh = malloc(sizeof(struct newhdf));
			if ( nh != NULL ) {
				strncpy(nh->name, names[i], 1023);
				nh->dw = dw;
				g_signal_connect(G_OBJECT(item), "activate",
			                  G_CALLBACK(displaywindow_newhdf), nh);
			}

		} else {

			char *tmp;

			item = gtk_menu_item_new_with_label(names[i]);

			tmp = hdfile_get_string_value(hdfile, names[i]);
			if ( tmp != NULL ) {

				GtkWidget *ss;
				GtkWidget *mss;

				mss = gtk_menu_new();
				ss = gtk_menu_item_new_with_label(tmp);
				gtk_widget_set_sensitive(ss, FALSE);
				gtk_menu_shell_append(GTK_MENU_SHELL(mss), ss);
				gtk_menu_item_set_submenu(GTK_MENU_ITEM(item),
				                          mss);

			}


		}

		gtk_menu_shell_append(GTK_MENU_SHELL(ms), item);

		free(names[i]);


	}

	free(is_group);
	free(is_image);

	return ms;
}


static void displaywindow_update_menus(DisplayWindow *dw)
{
	GtkWidget *ms;
	GtkWidget *w;

	ms = displaywindow_addhdfgroup(dw->hdfile, "/", dw);

	if ( ms == NULL ) {

		/* Too bad.  You'd better hope that /data/data exists... */
		ERROR("Couldn't get list of images in HDF file\n");
		w = gtk_ui_manager_get_widget(dw->ui,
					      "/ui/displaywindow/file/images");
		gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

		/* Add a dummy menu so that the user knows what's going on */
		ms = gtk_menu_new();
		w = gtk_ui_manager_get_widget(dw->ui,
		                              "/ui/displaywindow/file/images");
		gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

		return;

	}

	/* Make new menu be the submenu for File->Images */
	w = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/file/images");
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

	gtk_widget_show_all(ms);
}


static void displaywindow_disable(DisplayWindow *dw)
{
	GtkWidget *w;

	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/file/images");
	gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/binning");
	gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/boostint");
	gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/tools/numbers");
	gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);
}


static gint displaywindow_release(GtkWidget *widget, GdkEventButton *event,
                                  DisplayWindow *dw)
{
	if ( (event->type == GDK_BUTTON_RELEASE) && (event->button == 1) ) {

		g_signal_handler_disconnect(GTK_OBJECT(dw->drawingarea),
		                            dw->motion_callback);
		dw->motion_callback = 0;

	}

	return 0;
}


static gint displaywindow_motion(GtkWidget *widget, GdkEventMotion *event,
                                 DisplayWindow *dw)
{
	if ( dw->numbers_window == NULL ) return 0;

	dw->numbers_window->cx = event->x;
	dw->numbers_window->cy = dw->height - 1 - event->y;

	/* Schedule redraw */
	gtk_widget_queue_draw_area(dw->drawingarea, 0, 0,
	                           dw->width, dw->height);

	/* Update numbers window */
	numbers_update(dw);

	return 0;

}

static gint displaywindow_press(GtkWidget *widget, GdkEventButton *event,
                                DisplayWindow *dw)
{
	if ( dw->motion_callback != 0 ) {
		return 0;
	}

	if ( (event->type == GDK_BUTTON_PRESS) && (event->button == 1) ) {

		dw->motion_callback = g_signal_connect(
		                               GTK_OBJECT(dw->drawingarea),
		                               "motion-notify-event",
		                               G_CALLBACK(displaywindow_motion),
		                               dw);

		if ( dw->numbers_window != NULL ) {
			dw->numbers_window->cx = event->x;
			dw->numbers_window->cy = dw->height - 1 - event->y;
			numbers_update(dw);
		}

	}

	return 0;

}


DisplayWindow *displaywindow_open(const char *filename)
{
	DisplayWindow *dw;
	char *title;
	GtkWidget *vbox;

	dw = malloc(sizeof(DisplayWindow));
	if ( dw == NULL ) return NULL;
	dw->pixbuf = NULL;
	dw->binning_dialog = NULL;
	dw->show_col_scale = 0;
	dw->col_scale = NULL;
	dw->monochrome = 0;
	dw->boostint_dialog = NULL;
	dw->boostint = 1;
	dw->motion_callback = 0;
	dw->numbers_window = NULL;

	dw->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

	if ( filename == NULL ) {
		title = strdup("No file - hdfsee");
	} else {
		title = malloc(strlen(basename(filename))+14);
		sprintf(title, "%s - hdfsee", basename(filename));
	}
	gtk_window_set_title(GTK_WINDOW(dw->window), title);
	free(title);

	g_signal_connect(G_OBJECT(dw->window), "destroy",
			 G_CALLBACK(displaywindow_closed), dw);

	vbox = gtk_vbox_new(FALSE, 0);
	gtk_container_add(GTK_CONTAINER(dw->window), vbox);
	displaywindow_addmenubar(dw, vbox);

	dw->drawingarea = gtk_drawing_area_new();
	gtk_box_pack_start(GTK_BOX(vbox), dw->drawingarea, TRUE, TRUE, 0);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "expose-event",
			 G_CALLBACK(displaywindow_expose), dw);

	/* Open the file, if any */
	if ( filename != NULL ) {

		dw->hdfile = hdfile_open(filename);
		if ( dw->hdfile == NULL ) {
			fprintf(stderr, "Couldn't open file '%s'\n", filename);
			displaywindow_disable(dw);
		} else if ( hdfile_set_first_image(dw->hdfile, "/") ) {
			fprintf(stderr, "Couldn't select path\n");
			displaywindow_disable(dw);
		}

	} else {
		dw->hdfile = NULL;
		displaywindow_disable(dw);
	}

	gtk_window_set_resizable(GTK_WINDOW(dw->window), FALSE);
	gtk_widget_show_all(dw->window);

	dw->binning = INITIAL_BINNING;
	displaywindow_update(dw);

	gtk_widget_add_events(GTK_WIDGET(dw->drawingarea),
	                      GDK_BUTTON_PRESS_MASK
	                      | GDK_BUTTON_RELEASE_MASK
	                      | GDK_BUTTON1_MOTION_MASK);
	g_object_set(G_OBJECT(dw->drawingarea), "can-focus", TRUE, NULL);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-press-event",
	                 G_CALLBACK(displaywindow_press), dw);
	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-release-event",
	                 G_CALLBACK(displaywindow_release), dw);

	displaywindow_update_menus(dw);

	return dw;
}
