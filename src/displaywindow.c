/*
 * displaywindow.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "utils.h"


static void displaywindow_error(DisplayWindow *dw, const char *message)
{
	GtkWidget *window;

	window = gtk_message_dialog_new(GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_MESSAGE_WARNING,
					GTK_BUTTONS_CLOSE, message);
	gtk_window_set_title(GTK_WINDOW(window), "Error");

	g_signal_connect_swapped(window, "response",
				 G_CALLBACK(gtk_widget_destroy), window);
	gtk_widget_show(window);
}


static void displaywindow_update(DisplayWindow *dw)
{
	gint width;
	GdkGeometry geom;

	if ( dw->image == NULL ) {
		dw->width = 320;
		dw->height = 320;
	} else {

		double min_x, min_y, max_x, max_y;

		get_pixel_extents(dw->image->det,
		                  &min_x, &min_y, &max_x, &max_y);

		if ( min_x > 0.0 ) min_x = 0.0;
		if ( max_x < 0.0 ) max_x = 0.0;
		if ( min_y > 0.0 ) min_y = 0.0;
		if ( max_y < 0.0 ) max_y = 0.0;
		dw->min_x = min_x;
		dw->max_x = max_x;
		dw->min_y = min_y;
		dw->max_y = max_y;

		dw->width = (max_x - min_x) / dw->binning;
		dw->height = (max_y - min_y) / dw->binning;

		/* Add a thin border */
		dw->width += 2.0;
		dw->height += 2.0;
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

	/* Free old pixbufs */
	if ( dw->pixbufs != NULL ) {
		int i;
		for ( i=0; i<dw->n_pixbufs; i++ ) {
			gdk_pixbuf_unref(dw->pixbufs[i]);
		}
		free(dw->pixbufs);
	}
	if ( dw->col_scale != NULL ) {
		gdk_pixbuf_unref(dw->col_scale);
	}

	if ( dw->image != NULL ) {
		dw->pixbufs = render_panels(dw->image, dw->binning,
		                            dw->scale, dw->boostint,
		                            &dw->n_pixbufs);
	} else {
		dw->pixbufs = NULL;
	}

	dw->col_scale = render_get_colour_scale(20, dw->height, dw->scale);

	/* Schedule redraw */
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


static double ring_radius(struct image *image, double d)
{
	double theta, r, r_px;

	theta = asin(image->lambda / (2.0*d));
	r = image->det->panels[0].clen * tan(2.0*theta);
	r_px = r * image->det->panels[0].res;

	return r_px;
}


static void show_ring(cairo_t *cr, DisplayWindow *dw,
                      double d, const char *label)
{
	cairo_identity_matrix(cr);
	cairo_translate(cr, -dw->min_x/dw->binning, dw->max_y/dw->binning);
	cairo_arc(cr, 0.0, 0.0, ring_radius(dw->image, d)/dw->binning,
	          0.0, 2.0*M_PI);
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);
	cairo_rotate(cr, -M_PI/4.0);
	cairo_translate(cr, 0.0,
	                ring_radius(dw->image, d)/dw->binning-5.0);
	cairo_set_font_size(cr, 20.0);
	cairo_show_text(cr, label);
	cairo_fill(cr);
}


static gboolean displaywindow_expose(GtkWidget *da, GdkEventExpose *event,
				     DisplayWindow *dw)
{
	cairo_t *cr;
	int i;
	cairo_matrix_t basic_m;
	cairo_matrix_t m;

	cr = gdk_cairo_create(da->window);

	/* Blank white background */
	cairo_rectangle(cr, 0.0, 0.0, da->allocation.width,
			da->allocation.height);
	cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	cairo_fill(cr);

	/* Set up basic coordinate system
	 *  - origin in the centre, y upwards. */
	cairo_identity_matrix(cr);
	cairo_matrix_init(&m, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0);
	cairo_translate(cr, -dw->min_x/dw->binning, dw->max_y/dw->binning);
	cairo_transform(cr, &m);
	cairo_get_matrix(cr, &basic_m);

	if ( dw->pixbufs != NULL ) {

		int i;

		for ( i=0; i<dw->image->det->n_panels; i++ ) {

			struct panel p = dw->image->det->panels[i];
			int w = gdk_pixbuf_get_width(dw->pixbufs[i]);
			int h = gdk_pixbuf_get_height(dw->pixbufs[i]);
			cairo_matrix_t m;

			/* Start with the basic coordinate system */
			cairo_set_matrix(cr, &basic_m);

			/* Move to the right location */
			cairo_translate(cr, p.cx/dw->binning, p.cy/dw->binning);

			/* Twiddle directions according to matrix */
			cairo_matrix_init(&m, p.fsx, p.fsy, p.ssx, p.ssy,
			                      0.0, 0.0);
			cairo_transform(cr, &m);

			gdk_cairo_set_source_pixbuf(cr, dw->pixbufs[i],
			                            0.0, 0.0);
			cairo_rectangle(cr, 0.0, 0.0, w, h);
			cairo_fill(cr);

		}

	}

	if ( (dw->show_col_scale) && (dw->col_scale != NULL) ) {
		cairo_identity_matrix(cr);
		cairo_translate(cr, dw->width, 0.0);
		cairo_rectangle(cr, 0.0, 0.0, 20.0, dw->height);
		gdk_cairo_set_source_pixbuf(cr, dw->col_scale, 0.0, 0.0);
		cairo_fill(cr);
	}

	if ( dw->show_rings ) {

		/* Mark the beam */
		cairo_set_matrix(cr, &basic_m);
		cairo_arc(cr, 0.0, 0.0, 5.0/dw->binning, 0.0, 2.0*M_PI);
		cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
		cairo_fill(cr);

		/* Draw resolution circles */
		show_ring(cr, dw, 5.0e-10, "5A");
		show_ring(cr, dw, 4.0e-10, "4A");
		show_ring(cr, dw, 3.0e-10, "3A");
		show_ring(cr, dw, 2.0e-10, "2A");
		show_ring(cr, dw, 1.0e-10, "1A");

	}

	if ( dw->image->features == NULL ) {
		cairo_destroy(cr);
		return FALSE;
	}

	for ( i=0; i<image_feature_count(dw->image->features); i++ ) {

		double x, y;
		struct imagefeature *f;

		f = image_get_feature(dw->image->features, i);
		if ( f == NULL ) continue;

		x = f->x / (double)dw->binning;
		y = dw->height - f->y / (double)dw->binning;

		cairo_new_path(cr);
		cairo_arc(cr, x, y, 7.0/dw->binning, 0.0, 2.0*M_PI);
		switch ( dw->scale ) {

			case SCALE_COLOUR :
			cairo_set_source_rgb(cr, 1.0, 1.0, 0.0);
			break;

			case SCALE_MONO :
			cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
			break;

			case SCALE_INVMONO:
			cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
			break;

		}

		cairo_set_line_width(cr, 0.75);
		cairo_stroke(cr);

	}

	cairo_destroy(cr);

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
			if ((binning < dw->image->width/10)
			 && (binning < dw->image->height/10)) {
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
		 dw->image->width, dw->image->height);
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
		float boostint;
		int scanval;

		sboostint = gtk_entry_get_text(
		 			GTK_ENTRY(dw->boostint_dialog->entry));
		scanval = sscanf(sboostint, "%f", &boostint);
		if ( (scanval != 1) || (boostint <= 0) ) {
			displaywindow_error(dw, "Please enter a positive "
					"number for the intensity boost "
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
	snprintf(tmp, 63, "%.2f", dw->boostint);
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


static void load_features_from_file(struct image *image, const char *filename)
{
	FILE *fh;
	char *rval;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	do {
		char line[1024];
		float x, y, df;
		int r;
		signed int h, k, l;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		/* Try long format (output of pattern_sim --near-bragg) */
		r = sscanf(line, "%i %i %i %f (at %f,%f)",
		           &h, &k, &l, &df, &x, &y);
		if ( r == 6 ) {
			char name[32];
			snprintf(name, 31, "%i %i %i", h, k, l);
			image_add_feature(image->features, x, y, image, 1.0,
			                  strdup(name));
			continue;
		}

		r = sscanf(line, "%f %f", &x, &y);
		if ( r != 2 ) continue;

		image_add_feature(image->features, x, y, image, 1.0, NULL);

	} while ( rval != NULL );
}


static gint displaywindow_peaklist_response(GtkWidget *d, gint response,
                                            DisplayWindow *dw)
{
	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *filename;

		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		load_features_from_file(dw->image, filename);
		displaywindow_update(dw);

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
		NULL
	};

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window),
				     GTK_WINDOW(dw->window));

	gtk_about_dialog_set_name(GTK_ABOUT_DIALOG(window), "hdfsee");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window), PACKAGE_VERSION);
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2011 Thomas White <taw@physics.org> and others");
	gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(window),
		"Quick viewer for HDF files");
	gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(window),
		"(c) 2006-2011 Thomas White <taw@physics.org>\n");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"http://www.bitwiz.org.uk/");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static int load_geometry_file(DisplayWindow *dw, struct image *image,
                              const char *filename)
{
	struct detector *geom;
	GtkWidget *w;

	geom = get_detector_geometry(filename);
	if ( geom == NULL ) {
		displaywindow_error(dw, "Failed to load geometry file");
		return -1;
	}

	if ( (1+geom->max_fs != dw->image->width)
	  || (1+geom->max_ss != dw->image->height) ) {

		displaywindow_error(dw, "Geometry doesn't match image.");
		return -1;

	}

	if ( dw->loaded_geom != NULL ) free_detector_geometry(dw->loaded_geom);
	dw->loaded_geom = geom;


	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/usegeom");
	gtk_widget_set_sensitive(GTK_WIDGET(w), TRUE);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), TRUE);
	dw->use_geom = 1;
	displaywindow_update(dw);

	return 0;
}


static gint displaywindow_loadgeom_response(GtkWidget *d, gint response,
                                            DisplayWindow *dw)
{
	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *filename;

		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		if ( load_geometry_file(dw, dw->image, filename) == 0 ) {
			displaywindow_update(dw);
		}

		g_free(filename);

	}

	gtk_widget_destroy(d);

	return 0;
}


static gint displaywindow_load_geom(GtkWidget *widget, DisplayWindow *dw)
{
	GtkWidget *d;

	d = gtk_file_chooser_dialog_new("Load Geometry File",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_OPEN,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
	                                NULL);

	g_signal_connect(G_OBJECT(d), "response",
	                 G_CALLBACK(displaywindow_loadgeom_response), dw);

	gtk_widget_show_all(d);

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


static gint displaywindow_set_usegeom(GtkWidget *d, DisplayWindow *dw)
{
	GtkWidget *w;

	/* Get new value */
	w =  gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/usegeom");
	dw->use_geom = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));

	if ( dw->use_geom ) {
		dw->image->det = dw->loaded_geom;
	} else {
		dw->image->det = dw->simple_geom;
	}

	displaywindow_update(dw);

	return 0;
}


static gint displaywindow_set_rings(GtkWidget *d, DisplayWindow *dw)
{
	GtkWidget *w;

	/* Get new value */
	w =  gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/rings");
	dw->show_rings = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));

	gdk_window_invalidate_rect(dw->drawingarea->window, NULL, FALSE);

	return 0;
}



struct savedialog {
	DisplayWindow *dw;
	GtkWidget *cb;
};


static gint displaywindow_save_response(GtkWidget *d, gint response,
                                        struct savedialog *cd)
{
	DisplayWindow *dw = cd->dw;
	int r;

	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *file;
		int type;

		file = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		type = gtk_combo_box_get_active(GTK_COMBO_BOX(cd->cb));

		if ( type == 0 ) {
			/* FIXME: Use Cairo or something */
			r = render_png(dw->pixbufs[0], file);
		} else if ( type == 1 ) {
			r = render_tiff_fp(dw->image, file);
		} else if ( type == 2 ) {
			r = render_tiff_int16(dw->image, file, dw->boostint);
		} else {
			r = -1;
		}

		if ( r != 0 ) {
			displaywindow_error(dw, "Unable to save the image.");
		}

		g_free(file);

	}

	gtk_widget_destroy(d);
	free(cd);

	return 0;
}


static gint displaywindow_save(GtkWidget *widget, DisplayWindow *dw)
{
	GtkWidget *d, *hbox, *l, *cb;
	struct savedialog *cd;

	d = gtk_file_chooser_dialog_new("Save Image",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_SAVE,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                NULL);

	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(d),
	                                               TRUE);

	hbox = gtk_hbox_new(FALSE, 0);
	gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(d), hbox);
	cb = gtk_combo_box_new_text();
	gtk_box_pack_end(GTK_BOX(hbox), GTK_WIDGET(cb), TRUE, TRUE, 5);
	l = gtk_label_new("Save as type:");
	gtk_box_pack_end(GTK_BOX(hbox), GTK_WIDGET(l), FALSE, FALSE, 5);

	gtk_combo_box_append_text(GTK_COMBO_BOX(cb),
	       "PNG - 8 bit RGB (colour, binned, filtered, boosted)");
	gtk_combo_box_append_text(GTK_COMBO_BOX(cb),
	       "TIFF - Floating point (mono, unbinned, filtered, not boosted)");
	gtk_combo_box_append_text(GTK_COMBO_BOX(cb),
	       "TIFF - 16 bit signed integer (mono, unbinned, filtered, boosted)");
	gtk_combo_box_set_active(GTK_COMBO_BOX(cb), 0);

	cd = malloc(sizeof(*cd));
	cd->dw = dw;
	cd->cb = cb;

	g_signal_connect(G_OBJECT(d), "response",
	                 G_CALLBACK(displaywindow_save_response), cd);

	gtk_widget_show_all(d);

	return 0;
}


static gint displaywindow_set_colscale(GtkWidget *widget, DisplayWindow *dw)
{
	dw->show_col_scale = 1 - dw->show_col_scale;
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
	GtkWidget *hbox2;
	GtkWidget *table;
	GtkWidget *label;
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

	hbox2 = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox2), FALSE, FALSE, 5);
	label = gtk_label_new("Feature:");
	gtk_box_pack_start(GTK_BOX(hbox2), GTK_WIDGET(label), FALSE, FALSE, 5);
	nw->feat = gtk_label_new("-");
	gtk_box_pack_start(GTK_BOX(hbox2), GTK_WIDGET(nw->feat), FALSE, FALSE, 5);

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
	int imin;
	double dmin;
	struct imagefeature *f;

	for ( px=0; px<17; px++ ) {
	for ( py=0; py<17; py++ ) {

		char s[32];
		float val;
		GtkWidget *l;
		int x, y;
		int valid;

		x = dw->binning * dw->numbers_window->cx + (px-8);
		y = dw->binning * dw->numbers_window->cy + (17-py-8);

		if ( (x>=dw->image->width) || (y>=dw->image->height) ) {
			valid = 0;
			val = 0;
		} else {
			val = dw->image->data[x+y*dw->image->width];
			valid = 1;
		}

		if ( (x>0) && (y>0) && valid ) {
			if ( val > 0 ) {
				if ( log(val)/log(10) < 5 ) {
					snprintf(s, 31, "%.0f", val);
				} else {
					snprintf(s, 31, "HUGE");
				}
			} else {
				if ( log(-val)/log(10) < 4 ) {
					snprintf(s, 31, "%.0f", val);
				} else {
					snprintf(s, 31, "-HUGE");
				}
			}
		} else {
			strcpy(s, "--");
		}
		l = dw->numbers_window->labels[px+17*py];
		gtk_label_set_text(GTK_LABEL(l), s);

	}
	}

	if ( dw->image->features == NULL ) return;

	f = image_feature_closest(dw->image->features,
	                          dw->binning * dw->numbers_window->cx,
	                          dw->binning * dw->numbers_window->cy,
	                          &dmin, &imin);
	if ( dmin < 20.0 ) {
		gtk_label_set_text(GTK_LABEL(dw->numbers_window->feat),
                                   f->name);
        } else {
		gtk_label_set_text(GTK_LABEL(dw->numbers_window->feat),
                                   "-");
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


static gint displaywindow_setscale(GtkWidget *widget, GtkRadioAction *action,
                                   DisplayWindow *dw)
{
	switch ( gtk_radio_action_get_current_value(action) )
	{
		case 0 : dw->scale = SCALE_COLOUR; break;
		case 1 : dw->scale = SCALE_MONO; break;
		case 2 : dw->scale = SCALE_INVMONO; break;
	}
	displaywindow_update(dw);

	return 0;
}


static void displaywindow_addmenubar(DisplayWindow *dw, GtkWidget *vbox,
                                     int colscale)
{
	GError *error = NULL;
	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "SaveAction", GTK_STOCK_SAVE, "Save Image...", NULL, NULL,
			G_CALLBACK(displaywindow_save) },
		{ "CloseAction", GTK_STOCK_CLOSE, "_Close", NULL, NULL,
			G_CALLBACK(displaywindow_close) },

		{ "ViewAction", NULL, "_View", NULL, NULL, NULL },
		{ "ImagesAction", NULL, "Images", NULL, NULL, NULL },
		{ "BinningAction", NULL, "Set Binning...", "F3", NULL,
			G_CALLBACK(displaywindow_set_binning) },
		{ "BoostIntAction", NULL, "Boost Intensity...", "F5", NULL,
			G_CALLBACK(displaywindow_set_boostint) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "NumbersAction", NULL, "View Numbers...", "F2", NULL,
			G_CALLBACK(displaywindow_show_numbers) },
		{ "PeaksAction", NULL, "Peak Position Overlay...", NULL, NULL,
			G_CALLBACK(displaywindow_peak_overlay) },
		{ "LoadGeomAction", NULL, "Load Geometry File...", NULL, NULL,
			G_CALLBACK(displaywindow_load_geom) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About hdfsee...",
			NULL, NULL,
			G_CALLBACK(displaywindow_about) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	GtkToggleActionEntry toggles[] = {
		{ "ColScaleAction", NULL, "Show Colour Scale", NULL, NULL,
			G_CALLBACK(displaywindow_set_colscale), FALSE },
		{ "GeometryAction", NULL, "Use Detector Geometry", NULL, NULL,
			G_CALLBACK(displaywindow_set_usegeom), FALSE },
		{ "RingsAction", NULL, "Show Resolution Rings", NULL, NULL,
			G_CALLBACK(displaywindow_set_rings), FALSE },
	};
	guint n_toggles = G_N_ELEMENTS(toggles);
	GtkRadioActionEntry radios[] = {
		{ "ColAction", NULL, "Colour", NULL, NULL,
			SCALE_COLOUR },
		{ "MonoAction", NULL, "Monochrome", NULL, NULL,
			SCALE_MONO },
		{ "InvMonoAction", NULL, "Inverse Monochrome", NULL, NULL,
			SCALE_INVMONO },
	};
	guint n_radios = G_N_ELEMENTS(radios);

	dw->action_group = gtk_action_group_new("hdfseedisplaywindow");
	gtk_action_group_add_actions(dw->action_group, entries, n_entries, dw);
	gtk_action_group_add_toggle_actions(dw->action_group, toggles,
					    n_toggles, dw);
	gtk_action_group_add_radio_actions(dw->action_group, radios, n_radios,
	                                   colscale,
	                                   G_CALLBACK(displaywindow_setscale),
	                                   dw);

	dw->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(dw->ui, dw->action_group, 0);
	g_signal_connect(dw->ui, "add_widget",
			 G_CALLBACK(displaywindow_addui_callback), vbox);
	if ( gtk_ui_manager_add_ui_from_file(dw->ui,
	     DATADIR"/crystfel/hdfsee.ui", &error) == 0 ) {
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(dw->window),
				   gtk_ui_manager_get_accel_group(dw->ui));
	gtk_ui_manager_ensure_update(dw->ui);
}



static int geometry_fits(struct image *image, struct detector *geom)
{
	if ( (1+geom->max_fs != image->width)
	  || (1+geom->max_ss != image->height) ) return 0;

	return 1;
}


struct newhdf {
	DisplayWindow *dw;
	char name[1024];
};

static gint displaywindow_newhdf(GtkMenuItem *item, struct newhdf *nh)
{
	hdfile_set_image(nh->dw->hdfile, nh->name);
	hdf5_read(nh->dw->hdfile, nh->dw->image, 0, 0.0);

	/* Check that the geometry still fits */
	if ( !geometry_fits(nh->dw->image, nh->dw->simple_geom) ) {
		free_detector_geometry(nh->dw->simple_geom);
		nh->dw->simple_geom = simple_geometry(nh->dw->image);
	}

	if ( (nh->dw->loaded_geom != NULL )
	  && (!geometry_fits(nh->dw->image, nh->dw->loaded_geom)) ) {

		GtkWidget *w;

		free_detector_geometry(nh->dw->loaded_geom);

		w = gtk_ui_manager_get_widget(nh->dw->ui,
				      "/ui/displaywindow/view/usegeom");
		gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), FALSE);
		nh->dw->use_geom = 0;

	}

	if ( nh->dw->use_geom ) {
		nh->dw->image->det = nh->dw->loaded_geom;
	} else {
		nh->dw->image->det = nh->dw->simple_geom;
	}

	displaywindow_update(nh->dw);
	return 0;
}


static GtkWidget *displaywindow_addhdfgroup(struct hdfile *hdfile,
                                            const char *group,
                                            DisplayWindow *dw, GSList **rgp,
                                            const char *selectme)
{
	char **names;
	int *is_group;
	int *is_image;
	GtkWidget *ms;
	int n, i;

	if ( hdfile == NULL ) return NULL;

	names = hdfile_read_group(hdfile, &n, group, &is_group, &is_image);
	if ( n == 0 ) return NULL;

	ms = gtk_menu_new();

	for ( i=0; i<n; i++ ) {

		GtkWidget *item;
		GtkWidget *sub;

		if ( names[i] == NULL ) return NULL;

		if ( is_group[i] ) {

			item = gtk_menu_item_new_with_label(names[i]);

			sub = displaywindow_addhdfgroup(hdfile, names[i],
			                                dw, rgp, selectme);
			gtk_menu_item_set_submenu(GTK_MENU_ITEM(item), sub);

		} else if ( is_image[i] ) {

			struct newhdf *nh;

			item = gtk_radio_menu_item_new_with_label(*rgp,
			                                          names[i]);

			nh = malloc(sizeof(struct newhdf));
			if ( nh != NULL ) {
				strncpy(nh->name, names[i], 1023);
				nh->dw = dw;
				g_signal_connect(G_OBJECT(item), "activate",
			                  G_CALLBACK(displaywindow_newhdf), nh);
			}

			if ( (selectme != NULL)
			  && (strcmp(names[i], selectme) == 0) ) {
				gtk_check_menu_item_set_active(
				               GTK_CHECK_MENU_ITEM(item), TRUE);
			} else {
				gtk_check_menu_item_set_active(
				              GTK_CHECK_MENU_ITEM(item), FALSE);
			}

			*rgp = gtk_radio_menu_item_get_group(
		                                     GTK_RADIO_MENU_ITEM(item));

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


static GtkWidget *displaywindow_createhdfmenus(struct hdfile *hdfile,
                                               DisplayWindow *dw,
                                               const char *selectme)
{
	GSList *rg = NULL;

	return displaywindow_addhdfgroup(hdfile, "/", dw, &rg, selectme);
}


static void displaywindow_update_menus(DisplayWindow *dw, const char *selectme)
{
	GtkWidget *ms;
	GtkWidget *w;

	ms = displaywindow_createhdfmenus(dw->hdfile, dw, selectme);

	if ( ms == NULL ) {

		/* Too bad.  You'd better hope that /data/data exists... */
		ERROR("Couldn't get list of images in HDF file\n");
		w = gtk_ui_manager_get_widget(dw->ui,
					      "/ui/displaywindow/view/images");
		gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);

		/* Add a dummy menu so that the user knows what's going on */
		ms = gtk_menu_new();
		w = gtk_ui_manager_get_widget(dw->ui,
		                              "/ui/displaywindow/view/images");
		gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

		return;

	}

	/* Make new menu be the submenu for File->Images */
	w = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/view/images");
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
				      "/ui/displaywindow/file/rings");
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

	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/tools/peaks");
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


DisplayWindow *displaywindow_open(const char *filename, const char *peaks,
                                  int boost, int binning, int cmfilter,
                                  int noisefilter, int colscale,
                                  const char *element)
{
	DisplayWindow *dw;
	char *title;
	GtkWidget *vbox;
	GtkWidget *w;

	dw = calloc(1, sizeof(DisplayWindow));
	if ( dw == NULL ) return NULL;
	dw->pixbufs = NULL;
	dw->binning_dialog = NULL;
	dw->show_col_scale = 0;
	dw->col_scale = NULL;
	dw->boostint_dialog = NULL;
	dw->boostint = 1;
	dw->motion_callback = 0;
	dw->numbers_window = NULL;
	dw->image = NULL;
	dw->use_geom = 0;
	dw->show_rings = 0;

	dw->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

	if ( filename == NULL ) {
		title = strdup("No file - hdfsee");
	} else {
		char *bn = safe_basename(filename);
		title = malloc(strlen(bn)+14);
		sprintf(title, "%s - hdfsee", bn);
		free(bn);
	}
	gtk_window_set_title(GTK_WINDOW(dw->window), title);
	free(title);

	g_signal_connect(G_OBJECT(dw->window), "destroy",
			 G_CALLBACK(displaywindow_closed), dw);

	vbox = gtk_vbox_new(FALSE, 0);
	gtk_container_add(GTK_CONTAINER(dw->window), vbox);
	displaywindow_addmenubar(dw, vbox, colscale);

	dw->drawingarea = gtk_drawing_area_new();
	gtk_box_pack_start(GTK_BOX(vbox), dw->drawingarea, TRUE, TRUE, 0);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "expose-event",
			 G_CALLBACK(displaywindow_expose), dw);

	/* Open the file, if any */
	if ( filename != NULL ) {

		dw->hdfile = hdfile_open(filename);
		if ( dw->hdfile == NULL ) {
			ERROR("Couldn't open file '%s'\n", filename);
			displaywindow_disable(dw);
		} else {
			int fail = -1;

			if ( element == NULL ) {
				fail = hdfile_set_first_image(dw->hdfile, "/");
			} else {
				fail = hdfile_set_image(dw->hdfile, element);
			}

			if ( !fail ) {
				dw->image = calloc(1, sizeof(struct image));
				hdf5_read(dw->hdfile, dw->image, 0, 0.0);
			} else {
				ERROR("Couldn't select path\n");
				displaywindow_disable(dw);
			}
		}

	} else {
		dw->hdfile = NULL;
		displaywindow_disable(dw);
	}

	gtk_window_set_resizable(GTK_WINDOW(dw->window), FALSE);
	gtk_widget_show_all(dw->window);

	/* No geometry loaded initially */
	w = gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/usegeom");
	gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);
	dw->loaded_geom = NULL;
	dw->simple_geom = simple_geometry(dw->image);
	dw->image->det = dw->simple_geom;

	dw->scale = colscale;
	dw->binning = binning;
	dw->boostint = boost;
	dw->cmfilter = cmfilter;
	dw->noisefilter = noisefilter;
	displaywindow_update(dw);

	/* Peak list provided at startup? */
	if ( (dw->hdfile != NULL) && (peaks != NULL) ) {
		load_features_from_file(dw->image, peaks);
		displaywindow_update(dw);
	}

	gtk_widget_add_events(GTK_WIDGET(dw->drawingarea),
	                      GDK_BUTTON_PRESS_MASK
	                      | GDK_BUTTON_RELEASE_MASK
	                      | GDK_BUTTON1_MOTION_MASK);
	g_object_set(G_OBJECT(dw->drawingarea), "can-focus", TRUE, NULL);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-press-event",
	                 G_CALLBACK(displaywindow_press), dw);
	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-release-event",
	                 G_CALLBACK(displaywindow_release), dw);

	if ( dw->hdfile != NULL ) displaywindow_update_menus(dw, element);

	return dw;
}
