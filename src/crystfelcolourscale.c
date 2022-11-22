/*
 * crystfelcolourscale.c
 *
 * CrystFEL's colour scale widget
 *
 * Copyright © 2020-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2022 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gtk/gtk.h>
#include <glib-object.h>
#include <gsl/gsl_statistics_float.h>

#include "crystfelcolourscale.h"


G_DEFINE_TYPE(CrystFELColourScale, crystfel_colour_scale,
              GTK_TYPE_DRAWING_AREA)

static void redraw(CrystFELColourScale *cs)
{
	gint w, h;
	w = gtk_widget_get_allocated_width(GTK_WIDGET(cs));
	h = gtk_widget_get_allocated_height(GTK_WIDGET(cs));
	gtk_widget_queue_draw_area(GTK_WIDGET(cs), 0, 0, w, h);
}


static gint destroy_sig(GtkWidget *window, CrystFELColourScale *cs)
{
	return FALSE;
}


static gint realise_sig(GtkWidget *window, CrystFELColourScale *cs)
{
	return FALSE;
}


static gint button_press_sig(GtkWidget *window, GdkEventButton *event,
                             CrystFELColourScale *cs)
{
	cs->drag_start_x = event->x;
	cs->drag_start_y = event->y;
	return FALSE;
}


static gint motion_sig(GtkWidget *window, GdkEventMotion *event,
                       CrystFELColourScale *cs)
{
	double ddx, ddy;
	ddx = event->x - cs->drag_start_x;
	ddy = event->y - cs->drag_start_y;
	/* FIXME: Do something */
	redraw(cs);
	return FALSE;
}


static gint configure_sig(GtkWidget *window, GdkEventConfigure *rec,
                          CrystFELColourScale *cs)
{
	cs->visible_width = rec->width;
	cs->visible_height = rec->height;
	return FALSE;
}


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELColourScale *cs)
{
	int i;
	int mx = 0;
	double max_w = cs->visible_width;
	double bin_h = cs->visible_height/COLSCALE_N_BINS;

	cairo_save(cr);

	/* Overall background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	for ( i=0; i<COLSCALE_N_BINS; i++ ) {
		if ( cs->bins[i] > mx ) mx = cs->bins[i];
	}

	/* Origin at bottom right */
	cairo_translate(cr, cs->visible_width, cs->visible_height);
	cairo_scale(cr, -1.0, -1.0);

	for ( i=0; i<COLSCALE_N_BINS; i++ ) {
		cairo_rectangle(cr, 0.0, bin_h*i, max_w*cs->bins[i]/mx, bin_h);
		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_fill(cr);
	}

	cairo_restore(cr);

	return FALSE;
}


static GtkSizeRequestMode get_request_mode(GtkWidget *widget)
{
	return GTK_SIZE_REQUEST_CONSTANT_SIZE;
}


static void get_preferred_width(GtkWidget *widget, gint *min, gint *natural)
{
	*min = 0;
	*natural = 40;
}


static void get_preferred_height(GtkWidget *widget, gint *min, gint *natural)
{
	*min = 0;
	*natural = 640;
}


static void crystfel_colour_scale_class_init(CrystFELColourScaleClass *klass)
{
	GTK_WIDGET_CLASS(klass)->get_request_mode = get_request_mode;
	GTK_WIDGET_CLASS(klass)->get_preferred_width = get_preferred_width;
	GTK_WIDGET_CLASS(klass)->get_preferred_height = get_preferred_height;
	GTK_WIDGET_CLASS(klass)->get_preferred_height_for_width = NULL;
}


static void crystfel_colour_scale_init(CrystFELColourScale *cs)
{
}


GtkWidget *crystfel_colour_scale_new()
{
	CrystFELColourScale *cs;

	cs = g_object_new(CRYSTFEL_TYPE_COLOUR_SCALE, NULL);

	g_signal_connect(G_OBJECT(cs), "destroy",
	                 G_CALLBACK(destroy_sig), cs);
	g_signal_connect(G_OBJECT(cs), "realize",
	                 G_CALLBACK(realise_sig), cs);
	g_signal_connect(G_OBJECT(cs), "button-press-event",
	                 G_CALLBACK(button_press_sig), cs);
	g_signal_connect(G_OBJECT(cs), "motion-notify-event",
	                 G_CALLBACK(motion_sig), cs);
	g_signal_connect(G_OBJECT(cs), "configure-event",
	                 G_CALLBACK(configure_sig), cs);
	g_signal_connect(G_OBJECT(cs), "draw",
	                 G_CALLBACK(draw_sig), cs);

	gtk_widget_set_can_focus(GTK_WIDGET(cs), TRUE);
	gtk_widget_add_events(GTK_WIDGET(cs),
	                      GDK_POINTER_MOTION_HINT_MASK
	                       | GDK_BUTTON1_MOTION_MASK
	                       | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK
	                       | GDK_KEY_PRESS_MASK | GDK_KEY_RELEASE_MASK);

	gtk_widget_grab_focus(GTK_WIDGET(cs));

	gtk_widget_show(GTK_WIDGET(cs));

	return GTK_WIDGET(cs);
}


static double auto_scale_top(const struct image *image)
{
	int pn;
	double total_mean = 0.0;
	double total_variance = 0.0;

	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {

		long int i, j;
		int w, h;
		float *data;
		float this_mean;

		w = image->detgeom->panels[pn].w;
		h = image->detgeom->panels[pn].h;

		data = malloc(w*h*sizeof(float));
		if ( data == NULL ) return 100.0;

		j = 0;
		for ( i=0; i<w*h; i++ ) {
			if ( !image->bad[pn][i] ) {
				data[j++] = image->dp[pn][i];
			}
		}

		this_mean = gsl_stats_float_mean(data, 1, j);

		total_mean += this_mean;
		total_variance += gsl_stats_float_variance_m(data, 1, j,
		                                             this_mean);

		free(data);
	}

	return (total_mean/image->detgeom->n_panels)
	      + 10.0*sqrt(total_variance/image->detgeom->n_panels);
}


void image_min_max(struct image *image, double *pmin, double *pmax)
{
	int pn;
	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {
		int w, h;
		long int i;
		w = image->detgeom->panels[pn].w;
		h = image->detgeom->panels[pn].h;
		for ( i=0; i<w*h; i++ ) {
			if ( !image->bad[pn][i] ) {
				double v = image->dp[pn][i];
				*pmin = fmin(v, *pmin);
				*pmax = fmax(v, *pmax);
			}
		}
	}
}


void histogram_image(struct image *image,
                     int *bins, int n_bins,
                     double min, double max)
{
	int pn;
	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {
		int w, h;
		long int i;
		w = image->detgeom->panels[pn].w;
		h = image->detgeom->panels[pn].h;
		for ( i=0; i<w*h; i++ ) {
			if ( !image->bad[pn][i] ) {
				int bin;
				double v = image->dp[pn][i];
				bin = n_bins*(v-min)/(max-min);
				if ( bin < 0 ) bin = 0;
				if ( bin >= n_bins ) bin = n_bins-1;
				bins[bin]++;
			}
		}
	}
}


void crystfel_colour_scale_scan_image(CrystFELColourScale *cs,
                                      struct image *image)
{
	double range_min, range_max;
	int i;
	int n_filled = 0;

	if ( image == NULL ) return;

	image_min_max(image, &range_min, &range_max);

	for ( i=0; i<COLSCALE_N_BINS; i++ ) {
		cs->bins[i] = 0;
	}

	histogram_image(image, cs->bins, COLSCALE_N_BINS, range_min, range_max);

	for ( i=0; i<COLSCALE_N_BINS; i++ ) {
		if ( cs->bins[i] > 0 ) n_filled++;
	}

	if ( n_filled < 3 ) {
		ERROR("WARNING: Suspicious pixel value distribution - "
		      "are there still some bad pixels to mask?\n");
	}

	redraw(cs);
}
