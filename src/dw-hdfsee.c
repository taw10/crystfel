/*
 * dw-hdfsee.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2014      Takanori Nakane
 *   2012      Richard Kirian
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

#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cairo.h>
#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>

#include "dw-hdfsee.h"
#include "hdfsee-render.h"
#include "render.h"
#include "hdf5-file.h"
#include "hdfsee.h"
#include "utils.h"
#include "filters.h"
#include "events.h"


static void displaywindow_error(DisplayWindow *dw, const char *message)
{
	GtkWidget *window;

	window = gtk_message_dialog_new(GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_MESSAGE_WARNING,
					GTK_BUTTONS_CLOSE, "%s", message);
	gtk_window_set_title(GTK_WINDOW(window), "Error");

	g_signal_connect_swapped(window, "response",
				 G_CALLBACK(gtk_widget_destroy), window);
	gtk_widget_show(window);
}


/* Window closed - clean up */
static gint displaywindow_closed(GtkWidget *window, DisplayWindow *dw)
{

	if ( dw->hdfile != NULL ) {
		hdfile_close(dw->hdfile);
	}

	if ( dw->surf != NULL ) cairo_surface_destroy(dw->surf);

	if ( dw->pixbufs != NULL ) {
		int i;
		for ( i=0; i<dw->n_pixbufs; i++ ) {
			g_object_unref(dw->pixbufs[i]);
		}
		free(dw->pixbufs);
	}

	if ( dw->col_scale != NULL ) {
		g_object_unref(dw->col_scale);
	}

	if ( dw->image != NULL ) {

		free(dw->image->filename);
		free(dw->image->data);
		free(dw->image->flags);
		free(dw->image);
	}

	/* Notify 'main', so it can update the master list */
	hdfsee_window_closed(dw);

	return 0;
}


static double ring_radius(struct image *image, int p, double d)
{
	double theta, r, r_px;

	theta = asin(image->lambda / (2.0*d));
	r = image->det->panels[p].clen * tan(2.0*theta);
	r_px = r * image->det->panels[p].res;

	return r_px;
}


static void draw_panel_rectangle(cairo_t *cr, cairo_matrix_t *basic_m,
                                 DisplayWindow *dw, int i)
{
	struct panel p = dw->image->det->panels[i];
	int w = gdk_pixbuf_get_width(dw->pixbufs[i]);
	int h = gdk_pixbuf_get_height(dw->pixbufs[i]);
	cairo_matrix_t m;

	/* Start with the basic coordinate system */
	cairo_set_matrix(cr, basic_m);

	/* Move to the right location */
	cairo_translate(cr, p.cnx/dw->binning,
	                    p.cny/dw->binning);

	/* Twiddle directions according to matrix */
	cairo_matrix_init(&m, p.fsx, p.fsy, p.ssx, p.ssy,
	                      0.0, 0.0);
	cairo_transform(cr, &m);

	gdk_cairo_set_source_pixbuf(cr, dw->pixbufs[i],
	                            0.0, 0.0);
	cairo_rectangle(cr, 0.0, 0.0, w, h);
}


int render_adsc_uint16(DisplayWindow *dw, const char *filename)
{
	int x, y, fs, ss;
	double dfs, dss;
	int min_x = (int)dw->min_x;
	int max_x = (int)dw->max_x;
	int min_y = (int)dw->min_y;
	int max_y = (int)dw->max_y;
	int width = max_x - min_x;
	int height = max_y - min_y;
	unsigned short *buf;
	struct image *image = dw->image;
	FILE *fh;

	if ( image == NULL ) return 1;
	if ( image->det == NULL ) return 1;
	if ( image->det->n_panels == 0)  return 1;

	buf = calloc(width * height, sizeof(unsigned short));
	if ( buf == NULL ) return 1;

	fh = fopen(filename, "wb");
	if ( fh == NULL ) {
		free(buf);
		return 1;
	}

	fprintf(fh, "{\n"
		"HEADER_BYTES=512;\n"
		"DIM=2;\n"
		"BYTE_ORDER=little_endian;\n"
		"TYPE=unsigned_short;\n"
		"SIZE1=%d;\n"
		"SIZE2=%d;\n"
		"PIXEL_SIZE=%f;\n"
		"WAVELENGTH=%f;\n"
		"DISTANCE=%f;\n"
		"PHI=0.0;\n"
		"OSC_START=0.00;\n"
		"OSC_END=0.00;\n"
		"OSC_RANGE=0.00;\n"
		"AXIS=phi;\n"
		"BEAM_CENTER_X=%f;\n"
		"BEAM_CENTER_Y=%f;\n"
		"}\n",
		width, height, 1 / image->det->panels[0].res * 10e2,
		image->lambda * 10e9, image->det->panels[0].clen * 10e2,
		-min_x / image->det->panels[0].res * 10e2,
		-min_y / image->det->panels[0].res * 10e2);

	fseek(fh, 512, SEEK_SET);

	for ( y=min_y; y<max_y; y++ ) {
	for ( x=min_x; x<max_x; x++ ) {

		int val, invalid;
		unsigned short out;

		invalid = reverse_2d_mapping(x, y, &dfs, &dss, image->det);
		if ( invalid ) continue;

		fs = dfs;
		ss = dss;
		val = image->data[fs + image->width * ss];
		if ( val < 0 ) {
			out = 0;
		} else if ( val > 65535 ) {
			out = 65535;
		} else {
			out = val;
		}

		buf[(x - min_x) + (y - min_y) * width] = out;

	}
	}

	fwrite(buf, sizeof(unsigned short), width * height, fh);
	free(buf);
	fclose(fh);

	return 0;
}


static void draw_calib_focus_rectangle(cairo_t *cr, cairo_matrix_t *basic_m,
                                       DisplayWindow *dw, int i)
{
	struct panel p = dw->image->det->panels[i];
	int w = gdk_pixbuf_get_width(dw->pixbufs[i]);
	int h = gdk_pixbuf_get_height(dw->pixbufs[i]);
	cairo_matrix_t m;

	/* Start with the basic coordinate system */
	cairo_set_matrix(cr, basic_m);

	/* Move to the right location */
	cairo_translate(cr, p.cnx/dw->binning, p.cny/dw->binning);

	/* Twiddle directions according to matrix */
	cairo_matrix_init(&m, p.fsx, p.fsy, p.ssx, p.ssy, 0.0, 0.0);
	cairo_transform(cr, &m);

	cairo_set_line_width (cr, 3.0);
	cairo_set_source_rgb (cr, 255, 255, 255);
	cairo_rectangle(cr, 0.0, 0.0, w, h);

}


static void show_ring(cairo_t *cr, DisplayWindow *dw,
                      double d, const char *label, cairo_matrix_t *basic_m,
                      double r, double g, double b)
{
	struct detector *det;
	int i;

	det = dw->image->det;

	for ( i=0; i<det->n_panels; i++ ) {

		draw_panel_rectangle(cr, basic_m, dw, i);
		cairo_clip(cr);

		cairo_text_extents_t size;
		cairo_identity_matrix(cr);
		cairo_translate(cr, -dw->min_x/dw->binning,
		                     dw->max_y/dw->binning);
		cairo_arc(cr, 0.0, 0.0,
		          ring_radius(dw->image, i, d)/dw->binning,
			  0.0, 2.0*M_PI);
		cairo_set_source_rgb(cr, r, g, b);
		cairo_set_line_width(cr, 3.0/dw->binning);
		cairo_stroke(cr);

		cairo_reset_clip(cr);

		/* Any ideas for a better way of doing this? */
		if ( i == 0 ) {

			cairo_rotate(cr, -M_PI/4.0);
			cairo_translate(cr, 0.0,
				 ring_radius(dw->image, i, d)/dw->binning-5.0);
			cairo_set_font_size(cr, 17.0);
			cairo_text_extents(cr, label, &size);
			cairo_translate(cr, -size.width/2.0, 0.0);

			cairo_show_text(cr, label);
			cairo_fill(cr);

		}

	}
}


static void show_simple_ring(cairo_t *cr, DisplayWindow *dw,
                             double d, cairo_matrix_t *basic_m)
{
	cairo_identity_matrix(cr);
	cairo_translate(cr, -dw->min_x/dw->binning,
	                     dw->max_y/dw->binning);
	cairo_arc(cr, 0.0, 0.0, d / dw->binning, 0.0, 2.0*M_PI);
	cairo_set_source_rgb(cr, 0.3, 1.0, 0.3);
	cairo_set_line_width(cr, 1.0/dw->binning);
	cairo_stroke(cr);
}


static void maybe_draw_focus(DisplayWindow *dw, cairo_t *cr, int i,
                             cairo_matrix_t *basic_m)
{
	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		if ( &dw->image->det->panels[i] == dw->calib_mode_curr_p) {
			draw_calib_focus_rectangle(cr, basic_m, dw, i);
			cairo_stroke(cr);
		}
		break;

		case CALIBMODE_GROUPS:
		if ( dw->image->det->panels[i].rigid_group == dw->calib_mode_curr_rg )
		{
			draw_calib_focus_rectangle(cr, basic_m, dw, i);
			cairo_stroke(cr);
		}
		break;

		case CALIBMODE_ALL:
		draw_calib_focus_rectangle(cr, basic_m, dw, i);
		cairo_stroke(cr);
		break;

	}
}


static int draw_stuff(cairo_surface_t *surf, DisplayWindow *dw)
{
	cairo_t *cr;
	cairo_matrix_t basic_m;
	cairo_matrix_t m;

	cr = cairo_create(surf);

	/* Blank grey background */
	cairo_rectangle(cr, 0.0, 0.0, dw->width, dw->height);
	cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	cairo_fill(cr);

	if ( dw->image == NULL ) return 0;

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

			draw_panel_rectangle(cr, &basic_m, dw, i);
			cairo_fill(cr);

			if ( dw->calib_mode && dw->calib_mode_show_focus ) {
				maybe_draw_focus(dw, cr, i, &basic_m);
			}

		}

	}

	if ( dw->show_rings ) {

		/* Mark the beam */
		cairo_set_matrix(cr, &basic_m);
		cairo_arc(cr, 0.0, 0.0, 5.0/dw->binning, 0.0, 2.0*M_PI);
		cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
		cairo_fill(cr);

		/* Draw resolution circles */
		if ( dw->n_rings == -1 ) {
			/* n_rings == -1 means default behavior */
			show_ring(cr, dw, 10.0e-10, "10A", &basic_m,
			          1.0, 0.0, 0.0);
			show_ring(cr, dw, 9.0e-10, "9A", &basic_m,
			          1.0, 0.0, 0.0);
			show_ring(cr, dw, 8.0e-10, "8A", &basic_m,
			          1.0, 0.0, 0.0);
			show_ring(cr, dw, 7.0e-10, "7A", &basic_m,
			          1.0, 0.5, 0.0);
			show_ring(cr, dw, 6.0e-10, "6A", &basic_m,
			          1.0, 1.0, 0.0);
			show_ring(cr, dw, 5.0e-10, "5A", &basic_m,
			          0.0, 1.0, 0.0);
			show_ring(cr, dw, 4.0e-10, "4A", &basic_m,
			          0.2, 1.0, 0.2);
			show_ring(cr, dw, 3.0e-10, "3A", &basic_m,
			          0.4, 1.0, 0.4);
			show_ring(cr, dw, 2.0e-10, "2A", &basic_m,
			          0.6, 1.0, 0.6);
			show_ring(cr, dw, 1.0e-10, "1A", &basic_m,
			          0.8, 1.0, 0.8);
			show_ring(cr, dw, 1.0e-10, "0.5A", &basic_m,
			          1.0, 1.0, 1.0);
		} else {
			int i;
			for ( i=0; i<dw->n_rings; i++ ) {
				show_simple_ring(cr, dw, dw->ring_radii[i],
			                         &basic_m);
			}
		}

	}

	if ( (dw->show_col_scale) && (dw->col_scale != NULL) ) {
		cairo_identity_matrix(cr);
		cairo_translate(cr, dw->width, 0.0);
		cairo_rectangle(cr, 0.0, 0.0, 20.0, dw->height);
		gdk_cairo_set_source_pixbuf(cr, dw->col_scale, 0.0, 0.0);
		cairo_fill(cr);
	}

	/* Ensure a clean Cairo context, since the rings often cause
	 * matrix trouble */
	cairo_destroy(cr);
	cr = cairo_create(surf);

	if ( (dw->image->features != NULL) && (dw->show_peaks) ) {

		int i;

		cairo_set_matrix(cr, &basic_m);

		for ( i=0; i<image_feature_count(dw->image->features); i++ ) {

			double fs, ss;
			double x, y, xs, ys;
			struct imagefeature *f;
			struct panel *p;
			double radius = dw->ring_radius;

			f = image_get_feature(dw->image->features, i);
			if ( f == NULL ) continue;

			fs = f->fs;
			ss = f->ss;

			p = find_panel(dw->image->det, fs, ss);

			if ( p == NULL ) continue;

			xs = (fs-p->min_fs)*p->fsx + (ss-p->min_ss)*p->ssx;
			ys = (fs-p->min_fs)*p->fsy + (ss-p->min_ss)*p->ssy;
			x = xs + p->cnx;
			y = ys + p->cny;

			cairo_arc(cr, x/dw->binning, y/dw->binning,
				  radius, 0.0, 2.0*M_PI);
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

			cairo_set_line_width(cr, 0.75/dw->binning);
			cairo_stroke(cr);
;
		}

	}

	cairo_destroy(cr);

	return 0;
}


static void redraw_window(DisplayWindow *dw)
{
	int width;

	width = dw->width;
	if ( dw->show_col_scale ) width += 20;
	if ( dw->surf != NULL ) cairo_surface_destroy(dw->surf);
	dw->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
	                                      width, dw->height);
	draw_stuff(dw->surf, dw);

	gdk_window_invalidate_rect(dw->drawingarea->window, NULL, FALSE);
}


static void set_window_size(DisplayWindow *dw)
{
	gint width;

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
		dw->width += 2.0;;
		dw->height += 2.0;
	}

	width = dw->width;
	if ( dw->show_col_scale ) width += 20;

	gtk_widget_set_size_request(GTK_WIDGET(dw->drawingarea), width,
				    dw->height);

	gtk_window_resize(GTK_WINDOW(dw->window), width+30, dw->height+30);
}


static void update_colscale(DisplayWindow *dw)
{
	if ( dw->col_scale != NULL ) {
		g_object_unref(dw->col_scale);
	}
	dw->col_scale = render_get_colour_scale(20, dw->height, dw->scale);
}


static void displaywindow_update(DisplayWindow *dw)
{

	set_window_size(dw);
	update_colscale(dw);

	/* Free old pixbufs */
	if ( dw->pixbufs != NULL ) {
		int i;
		for ( i=0; i<dw->n_pixbufs; i++ ) {
			g_object_unref(dw->pixbufs[i]);
		}
		free(dw->pixbufs);
	}

	if ( dw->image != NULL ) {
		dw->pixbufs = render_panels(dw->image, dw->binning,
		                            dw->scale, dw->boostint,
		                            &dw->n_pixbufs);
	} else {
		dw->pixbufs = NULL;
	}

	redraw_window(dw);
}


static gboolean displaywindow_expose(GtkWidget *da, GdkEventExpose *event,
                                     DisplayWindow *dw)
{
	cairo_t *cr;

	cr = gdk_cairo_create(da->window);

	cairo_set_source_surface(cr, dw->surf, 0.0, 0.0);
	cairo_rectangle(cr, event->area.x, event->area.y,
	                event->area.width, event->area.height);
	cairo_fill(cr);

	cairo_destroy(cr);

	return FALSE;
}


static int write_png(const char *filename, DisplayWindow *dw)
{
	cairo_status_t r;
	cairo_surface_t *surf;

	surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
	                                  dw->width, dw->height);

	draw_stuff(surf, dw);

	r = cairo_surface_write_to_png(surf, filename);
	if ( r != CAIRO_STATUS_SUCCESS ) return 1;

	return 0;
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


static gint displaywindow_set_ringradius_response(GtkWidget *widget,
                                                  gint response,
                                                  DisplayWindow *dw)
{
	int done = 1;

	if ( response == GTK_RESPONSE_OK ) {

		const char *srad;
		float ringrad;
		int scanval;

		srad = gtk_entry_get_text(
		                       GTK_ENTRY(dw->ringradius_dialog->entry));
		scanval = sscanf(srad, "%f", &ringrad);
		if ( (scanval != 1) || (ringrad <= 0) ) {
			displaywindow_error(dw, "Please enter a positive "
					"number for the ring radius "
					"factor.");
			done = 0;
		} else {
			dw->ring_radius = ringrad;
			displaywindow_update(dw);
		}
	}

	if ( done ) {
		gtk_widget_destroy(dw->ringradius_dialog->window);
	}

	return 0;
}


static gint displaywindow_set_ringradius_destroy(GtkWidget *widget,
                                                 DisplayWindow *dw)
{
	free(dw->ringradius_dialog);
	dw->ringradius_dialog = NULL;
	return 0;
}


static gint displaywindow_set_ringradius_response_ac(GtkWidget *widget,
                                                     DisplayWindow *dw)
{
	return displaywindow_set_ringradius_response(widget, GTK_RESPONSE_OK,
	                                             dw);
}

/* Create a window to ask the user for a new ring radius */
static gint displaywindow_set_ringradius(GtkWidget *widget, DisplayWindow *dw)
{
	RingRadiusDialog *rd;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *table;
	GtkWidget *label;
	char tmp[64];

	if ( dw->ringradius_dialog != NULL ) {
		return 0;
	}

	if ( dw->hdfile == NULL ) {
		return 0;
	}

	rd = malloc(sizeof(RingRadiusDialog));
	if ( rd == NULL ) return 0;
	dw->ringradius_dialog = rd;

	rd->window = gtk_dialog_new_with_buttons("Ring Radius",
					GTK_WINDOW(dw->window),
					GTK_DIALOG_DESTROY_WITH_PARENT,
					GTK_STOCK_CANCEL, GTK_RESPONSE_CLOSE,
					GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);

	vbox = gtk_vbox_new(FALSE, 0);
	hbox = gtk_hbox_new(TRUE, 0);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(rd->window)->vbox),
			   GTK_WIDGET(hbox), FALSE, FALSE, 7);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(vbox), FALSE, FALSE, 5);

	table = gtk_table_new(3, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(table), FALSE, FALSE, 0);

	label = gtk_label_new("Ring Radius:");
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(label),
				  1, 2, 3, 4);

	rd->entry = gtk_entry_new();
	snprintf(tmp, 63, "%.2f", dw->ring_radius);
	gtk_entry_set_text(GTK_ENTRY(rd->entry), tmp);
	gtk_table_attach_defaults(GTK_TABLE(table), GTK_WIDGET(rd->entry),
				  2, 3, 3, 4);

	g_signal_connect(G_OBJECT(rd->entry), "activate",
			 G_CALLBACK(displaywindow_set_ringradius_response_ac),
			 dw);
	g_signal_connect(G_OBJECT(rd->window), "response",
			 G_CALLBACK(displaywindow_set_ringradius_response), dw);
	g_signal_connect(G_OBJECT(rd->window), "destroy",
			 G_CALLBACK(displaywindow_set_ringradius_destroy), dw);
	gtk_window_set_resizable(GTK_WINDOW(rd->window), FALSE);
	gtk_widget_show_all(rd->window);
	gtk_widget_grab_focus(GTK_WIDGET(rd->entry));

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
		float intensity, sigma, fs, ss;
		char phs[1024];
		char pn[32];
		int r;
		float cts, d;
		signed int h, k, l;
		struct panel *p = NULL;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		r = sscanf(line, "%i %i %i %f %s %f %f %f %f %s",
		                   &h, &k, &l, &intensity, phs, &sigma, &cts,
		                   &fs, &ss, pn);

		if ( r == 10 ) {

			/* Stream reflection list format 2.3 */
			char name[32];
			snprintf(name, 31, "%i %i %i", h, k, l);
			p = find_panel_by_name(image->det, pn);

			if ( p == NULL ) {

				ERROR("Unable to find panel %s "
				      "(no geometry file given?)\n", pn);

			} else {

				/* Convert coordinates to match rearranged
				 * panels in memory */
				fs = fs-p->orig_min_fs+p->min_fs;
				ss = ss-p->orig_min_ss+p->min_ss;

			}

			image_add_feature(image->features, fs, ss,
			                  image, 1.0, strdup(name));
			continue;

		} else if ( r == 9 ) {

			/* Stream reflection list format 2.2 or 2.1 */
			char name[32];

			snprintf(name, 31, "%i %i %i", h, k, l);

			p = find_orig_panel(image->det, fs, ss);

			if ( p == NULL ) {
				ERROR("Unable to find panel for %s "
				      "(no geometry file given?)\n", name);
			} else {

				/* Convert coordinates to match rearranged
				 * panels in memory */
				fs = fs - p->orig_min_fs + p->min_fs;
				ss = ss - p->orig_min_ss + p->min_ss;

			}

			image_add_feature(image->features, fs, ss,
			                  image, 1.0, strdup(name));
			continue;

		}


		/* Try long peak format from stream */
		r = sscanf(line, "%f %f %f %f %s", &fs, &ss, &d,
		           &intensity, pn);
		if ( r == 5 ) {

			p = find_panel_by_name(image->det, pn);
			if ( p == NULL ) {
				ERROR("Unable to find panel %s "
				      "(no geometry file given?)\n", pn);
			} else {

				/* Convert coordinates to match rearranged panels in
				 * memory */
				fs = fs - p->orig_min_fs + p->min_fs;
				ss = ss - p->orig_min_ss + p->min_ss;

			}

			image_add_feature(image->features, fs, ss, image, 1.0,
			                  "peak");
		} else if ( r == 4 ) {

			p = find_orig_panel(image->det, fs, ss);

			if ( p == NULL ) {
				ERROR("Unable to find panel "
				      "(no geometry file given?)\n");
			} else {

				/* Convert coordinates to match rearranged
				 * panels in memory */
				fs = fs - p->orig_min_fs + p->min_fs;
				ss = ss - p->orig_min_ss + p->min_ss;
			}
			image_add_feature(image->features, fs, ss, image, 1.0,
			                  "peak");

		}

	} while ( rval != NULL );

}


static gint displaywindow_peaklist_response(GtkWidget *d, gint response,
                                            DisplayWindow *dw)
{
	if ( response == GTK_RESPONSE_ACCEPT ) {

		char *filename;
		GtkWidget *w;

		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(d));

		load_features_from_file(dw->image, filename);
		dw->show_peaks = 1;
		w =  gtk_ui_manager_get_widget(dw->ui,
					    "/ui/displaywindow/view/showpeaks");
		gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), TRUE);

		redraw_window(dw);

		g_free(filename);

	}

	gtk_widget_destroy(d);

	return 0;;
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
		"(c) 2012-2013 Thomas White <taw@physics.org> and others");
	gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(window),
		"Quick viewer for HDF files");
	gtk_about_dialog_set_license(GTK_ABOUT_DIALOG(window),
		"(c) 2012-2013 Thomas White <taw@physics.org>\n");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"http://www.desy.de/~twhite/crystfel");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static int save_geometry_file(DisplayWindow *dw)
{
	GtkWidget *d;
	gchar *output_filename;
	int w;

	d = gtk_file_chooser_dialog_new("Save Detector Geometry",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_SAVE,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                        GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                NULL);

	gtk_dialog_run(GTK_DIALOG(d));
	output_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER (d));
	w = write_detector_geometry(dw->geom_filename, output_filename,
	                            dw->image->det);
	if ( w != 0 ) {
		ERROR("Error saving geometry!\n");
	}
	gtk_widget_destroy(d);
	g_free(output_filename);
	return w;
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


static gint displaywindow_set_calibmode(GtkWidget *d, DisplayWindow *dw)
{
	GtkWidget *w, *vbox;
	int val;

	w =  gtk_ui_manager_get_widget(dw->ui,
	                               "/ui/displaywindow/tools/calibmode");
	if ( dw->image->det == dw->simple_geom ) {
		gtk_check_menu_item_set_state(GTK_CHECK_MENU_ITEM(w), 0);
		return 0;
	}

	/* Get new value */
	val = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));

	/* When entering calibration mode */
	if ( val ) {

		guint cc;

		if ( (dw->calib_mode_curr_rg == NULL)
		  && (dw->calib_mode_curr_p == NULL) )
		{
			dw->calib_mode_curr_rg = dw->image->det->rigid_groups[0];
			dw->calib_mode_curr_p = dw->calib_mode_curr_rg->panels[0];
		}

		dw->calib_mode = CALIBMODE_PANELS;

		dw->statusbar = gtk_statusbar_new();
		gtk_widget_show(dw->statusbar);
		vbox = gtk_bin_get_child(GTK_BIN(dw->window));
		gtk_box_pack_end(GTK_BOX(vbox), dw->statusbar,
		                 FALSE, FALSE, 0);
		cc = gtk_statusbar_get_context_id(GTK_STATUSBAR(dw->statusbar),
		                                  "calibmode");
		gtk_statusbar_push(GTK_STATUSBAR(dw->statusbar), cc,
		                   "Calibration mode activated");
		displaywindow_update(dw);

	} else {

		dw->calib_mode = CALIBMODE_NONE;
		gtk_widget_destroy(dw->statusbar);
		dw->statusbar = NULL;
		displaywindow_update(dw);

	}

	return 0;
}


static gint displaywindow_set_rings(GtkWidget *d, DisplayWindow *dw)
{
	GtkWidget *w;

	/* Get new value */
	w =  gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/rings");
	dw->show_rings = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));

	redraw_window(dw);

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
			r = write_png(file, dw);
		} else if ( type == 1 ) {
			r = render_tiff_fp(dw->image, file);
		} else if ( type == 2 ) {
			r = render_tiff_int16(dw->image, file, dw->boostint);
		} else if (type == 3) {
		  r = render_adsc_uint16(dw, file);
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
	char *fn, *bfn;

	d = gtk_file_chooser_dialog_new("Save Image",
	                                GTK_WINDOW(dw->window),
	                                GTK_FILE_CHOOSER_ACTION_SAVE,
	                                GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                NULL);

	bfn = safe_basename(dw->image->filename);
	if ( bfn != NULL ) {
		fn = malloc(strlen(bfn)+10);
		if ( fn != NULL ) {
			sprintf(fn, "%s.png", bfn);
			STATUS("%s'\n", fn);
			gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(d),
			                                  fn);
			free(fn);
		}
		free(bfn);
	}

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
	       "TIFF - 16 bit signed integer "
	       "(mono, unbinned, filtered, boosted)");
	gtk_combo_box_append_text(GTK_COMBO_BOX(cb),
	       "ADSC - 16 bit unsigned integer "
	       "(unbinned, filtered, not boosted)");
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
	GtkWidget *w;

	w =  gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/colscale");
	dw->show_col_scale = gtk_check_menu_item_get_active(
	                                                GTK_CHECK_MENU_ITEM(w));
	set_window_size(dw);
	redraw_window(dw);
	return 0;
}


static gint displaywindow_set_peaks(GtkWidget *widget, DisplayWindow *dw)
{
	GtkWidget *w;

	w =  gtk_ui_manager_get_widget(dw->ui,
				      "/ui/displaywindow/view/showpeaks");
	dw->show_peaks = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(w));
	redraw_window(dw);
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
	gtk_box_pack_start(GTK_BOX(hbox2), GTK_WIDGET(nw->feat),
	                   FALSE, FALSE, 5);

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
	int ffs = 0;
	int fss = 0;
	int found = 0;

	for ( px=0; px<17; px++ ) {
	for ( py=0; py<17; py++ ) {

		char s[32];
		GtkWidget *l;
		int x, y;
		int invalid;
		double dfs, dss;
		int fs, ss;

		x = dw->binning * dw->numbers_window->cx + (px-8);
		y = dw->binning * dw->numbers_window->cy + (17-py-8);
		x += dw->min_x;
		y += dw->min_y;

		/* Map from unbinned mapped pixel coordinates to a panel */
		invalid = reverse_2d_mapping(x, y, &dfs, &dss, dw->image->det);
		fs = dfs;  ss = dss;

		if ( !invalid ) {

			float val;

			val = dw->image->data[fs+ss*dw->image->width];

			if ( (val > 0.0) && (log(val)/log(10.0) >= 5) ) {
				snprintf(s, 31, "HUGE");
			} else if ( (val < 0.0) && (log(-val)/log(10) >= 4) ) {
				snprintf(s, 31, "-HUGE");
			} else {
				size_t l, i;
				snprintf(s, 31, "%.4f", val);
				l = strlen(s);
				for ( i=l-1; i>0; i-- ) {
					if ( s[i] == '0' ) s[i] = '\0';
					if ( s[i] == '.' ) {
						s[i] = '\0';
						break;
					}
				}
			}

			if ( (px==8) && (py==8) ) {
				ffs = fs;
				fss = ss;
				found = 1;
			}

		} else {
			strcpy(s, "-");
		}
		l = dw->numbers_window->labels[px+17*py];
		gtk_label_set_text(GTK_LABEL(l), s);

	}
	}

	if ( found && (dw->image->features != NULL) ) {

		char text[64];

		f = image_feature_closest(dw->image->features, ffs, fss,
		                          &dmin, &imin, dw->image->det);
		if ( dmin < dw->ring_radius*dw->binning ) {
			strncpy(text, f->name, 32);
		} else {
			strcpy(text, "");
		}

		strcat(text, " (panel ");
		strncat(text, find_panel(dw->image->det, ffs, fss)->name, 20);
		strcat(text, ")");

		gtk_label_set_text(GTK_LABEL(dw->numbers_window->feat), text);

	} else {

		char text[64];

		text[0] = '\0';
		strcat(text, "Panel ");
		strncat(text, find_panel(dw->image->det, ffs, fss)->name, 20);

		gtk_label_set_text(GTK_LABEL(dw->numbers_window->feat), text);
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
		{ "RingRadiusAction", NULL, "Ring Radius...", "F6", NULL,
			G_CALLBACK(displaywindow_set_ringradius) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "NumbersAction", NULL, "View Numbers...", "F2", NULL,
			G_CALLBACK(displaywindow_show_numbers) },
		{ "PeaksAction", NULL, "Load Feature List...", NULL, NULL,
			G_CALLBACK(displaywindow_peak_overlay) },

		{ "EventsAction", NULL, "_Events", NULL, NULL, NULL },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About hdfsee...",
			NULL, NULL,
			G_CALLBACK(displaywindow_about) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	GtkToggleActionEntry toggles[] = {
		{ "ColScaleAction", NULL, "Colour Scale", NULL, NULL,
			G_CALLBACK(displaywindow_set_colscale), FALSE },
		{ "CalibModeAction", NULL, "Calibration Mode", NULL, NULL,
			G_CALLBACK(displaywindow_set_calibmode), FALSE },
		{ "RingsAction", NULL, "Resolution Rings", "F9", NULL,
			G_CALLBACK(displaywindow_set_rings), dw->show_rings },
		{ "ShowPeaksAction", NULL, "Features", "F8", NULL,
			G_CALLBACK(displaywindow_set_peaks), dw->show_peaks },
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


static void do_filters(DisplayWindow *dw)
{
	if ( dw->median_filter > 0 ) {
		filter_median(dw->image, dw->median_filter);
	}

	if ( dw->noisefilter ) {
		filter_noise(dw->image);
	}
}

struct newhdf {
	DisplayWindow *dw;
	GtkWidget *widget;
	char name[1024];
};

static gint displaywindow_newhdf(GtkMenuItem *item, struct newhdf *nh)
{
	gboolean a;
	int fail;

	if ( nh->dw->not_ready_yet ) return 0;

	a = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(nh->widget));
	if ( !a ) return 0;

	fail = hdf5_read(nh->dw->hdfile, nh->dw->image, nh->name, 0);
	if ( fail ) {
		ERROR("Couldn't load image");
		return 1;
	}

	nh->dw->simple_geom = simple_geometry(nh->dw->image);
	nh->dw->image->det = nh->dw->simple_geom;

	do_filters(nh->dw);
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
				nh->widget = item;
				g_signal_connect(G_OBJECT(item), "toggled",
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

			*rgp = gtk_radio_menu_item_get_group(GTK_RADIO_MENU_ITEM(item));

		} else {

			char *tmp;

			item = gtk_menu_item_new_with_label(names[i]);

			if ( hdfile_is_scalar(hdfile, names[i], 0) ) {
				tmp = hdfile_get_string_value(hdfile, names[i], NULL);
			} else {
				tmp = NULL;
			}
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
	GtkWidget *w;

	w = displaywindow_addhdfgroup(hdfile, "/", dw, &rg, selectme);

	return w;
}


struct newev {
    DisplayWindow *dw;
    GtkWidget *widget;
    int new_ev;
};


static int displaywindow_update_menus(DisplayWindow *dw, const char *selectme)
{
	GtkWidget *ms;
	GtkWidget *w;

	ms = displaywindow_createhdfmenus(dw->hdfile, dw, selectme);

	if ( ms == NULL ) return 1;

	/* Make new menu be the submenu for File->Images */
	w = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/view/images");
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ms);

	gtk_widget_show_all(ms);

	return 0;
}


static gint displaywindow_newevent(GtkMenuItem *item, struct newev *ne)
{
	gboolean a;
	int fail;
	int i;

	if ( ne->dw->not_ready_yet ) return 0;

	a = gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(ne->widget));
	if ( !a ) return 0;

	float *old_data = ne->dw->image->data;
	uint16_t *old_flags = ne->dw->image->flags;
	float **old_dp = ne->dw->image->dp;
	int **old_bad = ne->dw->image->bad;

	fail = hdf5_read2(ne->dw->hdfile, ne->dw->image,
                          ne->dw->ev_list->events[ne->new_ev], 0);
	if ( fail ) {
		ERROR("Couldn't load image");
		return 1;
	}

	ne->dw->curr_event = ne->new_ev;

	do_filters(ne->dw);
	displaywindow_update(ne->dw);
	for (i = 0; i < ne->dw->image->det->n_panels; i++) {
		free(old_dp[i]);
		free(old_bad[i]);
	}
	free(old_data);
	free(old_flags);
	free(old_dp);
	free(old_bad);
	return 0;
}


static int displaywindow_update_event_menu(DisplayWindow *dw,
                                           struct event_list *ev_list,
                                           int curr_event)
{

	int ei;
	GtkWidget *w;
	GtkWidget *ww;
	GSList *grp = NULL;

	w = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/events");
	ww = gtk_menu_new();

	for ( ei=0; ei<ev_list->num_events; ei++ ) {

		GtkWidget *www;
		struct newev *ne;
		char *ev_string;

		ev_string = get_event_string(ev_list->events[ei]);
		www = gtk_radio_menu_item_new_with_label(grp, ev_string);
		free(ev_string);

		ne = malloc(sizeof(struct newev));
		if ( ne != NULL ) {

			ne->widget = www;
			ne->dw = dw;
			ne->new_ev = ei;

			g_signal_connect(G_OBJECT(www), "toggled",
			                 G_CALLBACK(displaywindow_newevent), ne);

		}

		if ( ei == dw->curr_event ) {
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(www),
			                               TRUE);
		} else {
			gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(www),
			                               FALSE);
		}

		gtk_menu_shell_append(GTK_MENU_SHELL(ww), www);

		grp = gtk_radio_menu_item_get_group(GTK_RADIO_MENU_ITEM(www));
		gtk_widget_show_all(www);

	}

	gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), ww);

	return 0;
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


static void calibmode_press(DisplayWindow *dw, GdkEventButton *event)
{
	int x,y;
	double dfs, dss;
	int fs, ss, revmap_return_value;
	char statusbar_string[80];
	guint cc;

	cc = gtk_statusbar_get_context_id(GTK_STATUSBAR(dw->statusbar),
	                                  "calibmode");

	x = dw->binning * (event->x);
	y = dw->binning * (dw->height - 1 - event->y);
	x += dw->min_x;
	y += dw->min_x;
	snprintf(statusbar_string, 80,
	         "Last clicked position: x: %i, y: %i (Not in panel)", x, y);

	revmap_return_value = reverse_2d_mapping(x, y, &dfs, &dss,
	                                         dw->image->det);
	if ( revmap_return_value == 0 ) {
		fs = dfs;
		ss = dss;
		snprintf(statusbar_string, 80,
                        "Last clicked position: x: %i, y: %i, fs: %u, ss: %u,"
                        " (panel %s)",
                        x, y, fs, ss, find_panel(dw->image->det, fs, ss)->name);
	}
	gtk_statusbar_push(GTK_STATUSBAR(dw->statusbar), cc, statusbar_string);
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

		if ( dw->statusbar != NULL ) {
			calibmode_press(dw, event);
		}

	}

	return 0;
}


static int curr_rg_pointer_index(DisplayWindow *dw)
{
	int r;

	for ( r=0; r<dw->image->det->n_rigid_groups; ++r) {
		if ( dw->image->det->rigid_groups[r] == dw->calib_mode_curr_rg ) {
			return r;
		}
	}

	/* Never reached (we hope) */
	return 999;
}


static int curr_p_pointer_index(DisplayWindow *dw)
{
	int p;

	for ( p=0; p<dw->image->det->n_panels; ++p) {
		if ( &dw->image->det->panels[p] == dw->calib_mode_curr_p ) {
			return p;
		}
	}

	/* Never reached (we hope) */
	return 999;
}


static void select_next_group(DisplayWindow *dw, int num_rg)
{
	if ( dw->calib_mode_curr_rg == dw->image->det->rigid_groups[num_rg-1] ) {
		dw->calib_mode_curr_rg = dw->image->det->rigid_groups[0];
	} else {
		dw->calib_mode_curr_rg =
		  dw->image->det->rigid_groups[curr_rg_pointer_index(dw)+1];
	}
}


static void select_prev_group(DisplayWindow *dw, int num_rg)
{
	if ( dw->calib_mode_curr_rg == dw->image->det->rigid_groups[0] ) {
		dw->calib_mode_curr_rg = dw->image->det->rigid_groups[num_rg-1];
	} else {
		dw->calib_mode_curr_rg =
		  dw->image->det->rigid_groups[curr_rg_pointer_index(dw)-1];
	}
}


static void select_next_panel(DisplayWindow *dw, int num_p)
{
	if ( dw->calib_mode_curr_p == &dw->image->det->panels[num_p-1] ) {
		dw->calib_mode_curr_p = &dw->image->det->panels[0];
	} else {
		dw->calib_mode_curr_p =
		  &dw->image->det->panels[curr_p_pointer_index(dw)+1];
	}
}


static void select_prev_panel(DisplayWindow *dw, int num_p)
{
	if ( dw->calib_mode_curr_p == &dw->image->det->panels[0] ) {
		dw->calib_mode_curr_p = &dw->image->det->panels[num_p-1];
	} else {
		dw->calib_mode_curr_p =
		  &dw->image->det->panels[curr_p_pointer_index(dw)-1];
	}
}


static void toggle_calibmode_groupmode(DisplayWindow *dw)
{
	struct rigid_group *rg;
	struct detector *det = dw->image->det;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		if ( det->n_rigid_groups != det->n_panels ) {
			/* Only change if there are any rigid groups defined */
			dw->calib_mode = CALIBMODE_GROUPS;
			rg = dw->calib_mode_curr_p->rigid_group;
			dw->calib_mode_curr_rg = rg;
		} else {
			/* ...otherwise skip to ALL mode */
			dw->calib_mode = CALIBMODE_ALL;
		}
		break;

		case CALIBMODE_GROUPS:
		dw->calib_mode = CALIBMODE_ALL;
		break;

		case CALIBMODE_ALL:
		dw->calib_mode = CALIBMODE_PANELS;
		dw->calib_mode_curr_p = dw->calib_mode_curr_rg->panels[0];
		break;

	}
}


static void calibmode_next(DisplayWindow *dw)
{
	int n;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		n = dw->image->det->n_panels;
		select_next_panel(dw, n);
		break;

		case CALIBMODE_GROUPS:
		n = dw->image->det->n_rigid_groups;
		select_next_group(dw, n);
		break;

		case CALIBMODE_ALL:
		break;

	}
}


static void calibmode_prev(DisplayWindow *dw)
{
	int n;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		n = dw->image->det->n_panels;
		select_prev_panel(dw, n);
		break;

		case CALIBMODE_GROUPS:
		n = dw->image->det->n_rigid_groups;
		select_prev_group(dw, n);
		break;

		case CALIBMODE_ALL:
		break;

	}
}


static void calibmode_up(DisplayWindow *dw)
{
	int pi;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		dw->calib_mode_curr_p->cny += 1.0;
		break;

		case CALIBMODE_GROUPS:
		for ( pi=0; pi<dw->calib_mode_curr_rg->n_panels; ++pi ) {
			dw->calib_mode_curr_rg->panels[pi]->cny += 1.0;
		}
		break;

		case CALIBMODE_ALL:
		for ( pi=0; pi<dw->image->det->n_panels; ++pi ) {
			dw->image->det->panels[pi].cny += 1.0;
		}
		break;

	}
}


static void calibmode_down(DisplayWindow *dw)
{
	int pi;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		dw->calib_mode_curr_p->cny -= 1.0;
		break;

		case CALIBMODE_GROUPS:
		for ( pi=0; pi<dw->calib_mode_curr_rg->n_panels; ++pi ) {
			dw->calib_mode_curr_rg->panels[pi]->cny -= 1.0;
		}
		break;

		case CALIBMODE_ALL:
		for ( pi=0; pi<dw->image->det->n_panels; ++pi ) {
			dw->image->det->panels[pi].cny -= 1.0;
		}
		break;

	}
}


static void calibmode_left(DisplayWindow *dw)
{
	int pi;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		dw->calib_mode_curr_p->cnx -= 1.0;
		break;

		case CALIBMODE_GROUPS:
		for ( pi=0; pi<dw->calib_mode_curr_rg->n_panels; ++pi ) {
			dw->calib_mode_curr_rg->panels[pi]->cnx -= 1.0;
		}
		break;

		case CALIBMODE_ALL:
		for ( pi=0; pi<dw->image->det->n_panels; ++pi ) {
			dw->image->det->panels[pi].cnx -= 1.0;
		}
		break;

	}
}


static void calibmode_right(DisplayWindow *dw)
{
	int pi;

	switch ( dw->calib_mode ) {

		case CALIBMODE_NONE:
		break;

		case CALIBMODE_PANELS:
		dw->calib_mode_curr_p->cnx += 1.0;
		break;

		case CALIBMODE_GROUPS:
		for ( pi=0; pi<dw->calib_mode_curr_rg->n_panels; ++pi ) {
			dw->calib_mode_curr_rg->panels[pi]->cnx += 1.0;
		}
		break;

		case CALIBMODE_ALL:
		for ( pi=0; pi<dw->image->det->n_panels; ++pi ) {
			dw->image->det->panels[pi].cnx += 1.0;
		}
		break;

	}
}


static gint displaywindow_keypress(GtkWidget *widget, GdkEventKey *event,
                                   DisplayWindow *dw)
{
	int s;

	if ( !dw->calib_mode ) {
		return 0;
	}

	switch ( event->keyval ) {

		case GDK_Up:
		case GDK_KP_Up:
		calibmode_up(dw);
		redraw_window(dw);
		break;

		case GDK_Down:
		case GDK_KP_Down:
		calibmode_down(dw);
		redraw_window(dw);
		break;

		case GDK_Left:
		case GDK_KP_Left:
		calibmode_left(dw);
		redraw_window(dw);
		break;

		case GDK_Right:
		case GDK_KP_Right:
		calibmode_right(dw);
		redraw_window(dw);
		break;

		case GDK_plus:
		case GDK_KP_Add:
		calibmode_next(dw);
		redraw_window(dw);
		break;

		case GDK_minus:
		case GDK_KP_Subtract:
		calibmode_prev(dw);
		redraw_window(dw);
		break;

		case GDK_f:
		dw->calib_mode_show_focus = 1 - dw->calib_mode_show_focus;
		redraw_window(dw);
		break;

		case GDK_g:
		toggle_calibmode_groupmode(dw);
		redraw_window(dw);
		break;

		case GDK_s:
		s = save_geometry_file(dw);
		if ( s != 0 ) {
			if ( s != 2 ) {
				displaywindow_error(dw,
				       "Unable to save the detector geometry.");
			}
		}
		break;

	}

	return 0;
}


DisplayWindow *displaywindow_open(char *filename, char *geom_filename,
                                  const char *peaks,
                                  double boost, int binning,
                                  int noisefilter, int calibmode, int colscale,
                                  const char *element, const char *event,
                                  struct detector *det_geom,
                                  struct beam_params *beam,
                                  int show_rings, double *ring_radii,
                                  int n_rings, double ring_size,
                                  int median_filter)
{
	DisplayWindow *dw;
	char *title;
	GtkWidget *vbox;
	int check;
	GtkWidget *w;
	GtkWidget *ww;

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
	dw->simple_geom = NULL;
	dw->image = NULL;
	dw->show_rings = show_rings;
	dw->show_peaks = 0;
	dw->scale = colscale;
	dw->binning = binning;
	dw->boostint = boost;
	dw->noisefilter = noisefilter;
	dw->not_ready_yet = 1;
	dw->surf = NULL;
	dw->ring_radius = ring_size;
	dw->ring_radii = ring_radii;
	dw->n_rings = n_rings;
	dw->median_filter = median_filter;
	dw->image = calloc(1, sizeof(struct image));
	dw->calib_mode = CALIBMODE_NONE;
	dw->calib_mode_curr_rg = NULL;
	dw->calib_mode_curr_p = NULL;
	dw->calib_mode_show_focus = 1;
	dw->statusbar = NULL;
	dw->multi_event = 0;
	dw->ev_list = NULL;
	if ( geom_filename != NULL ) {
		dw->geom_filename = strdup(geom_filename);
	} else {
		dw->geom_filename = NULL;
	}
	dw->image->det = det_geom;
	dw->image->beam = beam;
	dw->image->lambda = 0.0;
	dw->image->filename = filename;

	dw->hdfile = hdfile_open(filename);
	if ( dw->hdfile == NULL ) {
		ERROR("Couldn't open file: %s\n", filename);
		free(dw->geom_filename);
		free(dw);
		return NULL;
	}

	if ( dw->image->det != NULL && ( dw->image->det->path_dim != 0 ||
	                                 dw->image->det->dim_dim != 0  )) {

		dw->multi_event = 1;

		dw->ev_list = fill_event_list(dw->hdfile, dw->image->det);

		if ( dw->ev_list == NULL ) {
			ERROR("Error while parsing file structure\n");
			free_event_list(dw->ev_list);
			free(dw->geom_filename);
			free(dw);
			return NULL;
		}
		if ( dw->ev_list->num_events == 0 ) {
			ERROR("Multi-event geometry file but no events found "
			      "in data file\n");
			free_event_list(dw->ev_list);
			free(dw->geom_filename);
			free(dw);
			return NULL;
		} else {
			dw->curr_event = 0;
		}

	}

	if ( dw->image->det != NULL ) {

		if ( dw->multi_event ) {
			struct event *ev;
			if ( event != NULL ) {
				ev = get_event_from_event_string(event);
				dw->curr_event = find_event(ev, dw->ev_list);
				if ( dw->curr_event == dw->ev_list->num_events)
				{
					ERROR("Invalid event\n");
					return NULL;
				}
			} else {
				dw->curr_event = 0;
				ev = dw->ev_list->events[dw->curr_event];
			}
			check = hdf5_read2(dw->hdfile, dw->image, ev, 0);
		} else {
			check = hdf5_read2(dw->hdfile, dw->image, NULL, 0);
		}

	} else {
		check =	hdf5_read(dw->hdfile, dw->image, element, 0);
	}
	if ( check ) {
		ERROR("Couldn't load file\n");
		hdfile_close(dw->hdfile);
		free(dw->geom_filename);
		return NULL;
	}

	dw->image->filename = strdup(filename);

	if ( dw->image->det == NULL ) {
		dw->simple_geom = simple_geometry(dw->image);
		dw->image->det = dw->simple_geom;
	}

	/* Filters need geometry */
	do_filters(dw);

	/* Peak list provided at startup? */
	if ( peaks != NULL ) {
		load_features_from_file(dw->image, peaks);
		dw->show_peaks = 1;
	}

	dw->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	char *bn = safe_basename(filename);
	title = malloc(strlen(bn)+14);
	sprintf(title, "%s - hdfsee", bn);
	free(bn);
	gtk_window_set_title(GTK_WINDOW(dw->window), title);
	free(title);

	g_signal_connect(G_OBJECT(dw->window), "destroy",
			 G_CALLBACK(displaywindow_closed), dw);

	vbox = gtk_vbox_new(FALSE, 0);
	gtk_container_add(GTK_CONTAINER(dw->window), vbox);
	displaywindow_addmenubar(dw, vbox, colscale);

	dw->drawingarea = gtk_drawing_area_new();
	dw->scrollarea = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(dw->scrollarea),
	                               GTK_POLICY_AUTOMATIC,
	                               GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(dw->scrollarea),
	                                      dw->drawingarea);
	gtk_box_pack_start(GTK_BOX(vbox), dw->scrollarea, TRUE, TRUE, 0);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "expose-event",
			 G_CALLBACK(displaywindow_expose), dw);

	gtk_window_set_resizable(GTK_WINDOW(dw->window), TRUE);
	gtk_widget_show_all(dw->window);

	w = gtk_ui_manager_get_widget(dw->ui,
	                   "/ui/displaywindow/view/images");

	if ( dw->image->det != dw->simple_geom ) {
		gtk_widget_set_sensitive(GTK_WIDGET(w), FALSE);
	}

	ww = gtk_ui_manager_get_widget(dw->ui, "/ui/displaywindow/events");

	if ( dw->image->det == dw->simple_geom || dw->multi_event == 0) {
		gtk_widget_set_sensitive(GTK_WIDGET(ww), FALSE);
	}

	displaywindow_update(dw);

	gtk_widget_add_events(GTK_WIDGET(dw->drawingarea),
	                      GDK_BUTTON_PRESS_MASK
	                      | GDK_BUTTON_RELEASE_MASK
	                      | GDK_BUTTON1_MOTION_MASK
	                      | GDK_KEY_PRESS_MASK);
	g_object_set(G_OBJECT(dw->drawingarea), "can-focus", TRUE, NULL);

	gtk_widget_grab_focus(dw->drawingarea);

	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-press-event",
	                 G_CALLBACK(displaywindow_press), dw);
	g_signal_connect(GTK_OBJECT(dw->drawingarea), "button-release-event",
	                 G_CALLBACK(displaywindow_release), dw);
	g_signal_connect(GTK_OBJECT(dw->drawingarea), "key-press-event",
	                 G_CALLBACK(displaywindow_keypress), dw);

	if ( dw->image->det == dw->simple_geom ) {
		displaywindow_update_menus(dw, element);
	} else {
		if ( dw->multi_event ) {
			displaywindow_update_event_menu(dw, dw->ev_list,
			                                dw->curr_event);
		}
	}

	dw->not_ready_yet = 0;

	return dw;
}
