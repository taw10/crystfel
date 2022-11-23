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
	cs->drag_min = cs->lo;
	gtk_widget_grab_focus(GTK_WIDGET(cs));
	return TRUE;
}


static void make_histogram(CrystFELColourScale *cs)
{
	int i;
	int n_bins = COLSCALE_N_BINS;
	float min = cs->lo;
	float max = cs->hi;

	for ( i=0; i<COLSCALE_N_BINS; i++ ) {
		cs->bins[i] = 0;
	}

	for ( i=0; i<COLSCALE_SAMPLE_SIZE; i++ ) {
		int bin;
		double v = cs->sample[i];
		bin = n_bins*(v-min)/(max-min);
		if ( (bin >= 0) && (bin < n_bins) ) {
			cs->bins[bin]++;
		}
	}
}


static gint motion_sig(GtkWidget *window, GdkEventMotion *event,
                       CrystFELColourScale *cs)
{
	double span = cs->hi - cs->lo;
	//double ddx = event->x - cs->drag_start_x;
	//double ddy = event->y - cs->drag_start_y;

	cs->lo = cs->drag_min + span*(event->y - cs->drag_start_y)/cs->visible_height;
	cs->hi = cs->lo + span;

	make_histogram(cs);
	gtk_widget_queue_draw(GTK_WIDGET(cs));

	if ( event->is_hint ) {
		gdk_window_get_pointer(gtk_widget_get_window(GTK_WIDGET(cs)),
		                       NULL, NULL, NULL);
	}

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


static void handle_scroll_click(double zoom_scale, CrystFELColourScale *cs,
                                double pos)
{
	cs->lo = pos - (pos - cs->lo)*zoom_scale;
	cs->hi = pos + (cs->hi - pos)*zoom_scale;
}


static gint scroll_sig(GtkWidget *widget, GdkEventScroll *event,
                       CrystFELColourScale *cs)
{
	double xs, ys;
	double span = cs->hi - cs->lo;
	double pos = cs->lo + span*(1.0-event->y/cs->visible_height);

	switch ( event->direction ) {

		case GDK_SCROLL_UP:
		handle_scroll_click(0.9, cs, pos);
		break;

		case GDK_SCROLL_DOWN:
		handle_scroll_click(1.1, cs, pos);
		break;

		case GDK_SCROLL_SMOOTH:
		if ( gdk_event_get_scroll_deltas((GdkEvent *)event, &xs, &ys) ) {
			handle_scroll_click(1.0+ys*0.1, cs, pos);
		}
		break;

		case GDK_SCROLL_LEFT:
		case GDK_SCROLL_RIGHT:
		return FALSE;  /* Not handled here */

		default:
		STATUS("Unhandled scroll direction %i\n", event->direction);
		return FALSE;
	}

	make_histogram(cs);
	gtk_widget_grab_focus(GTK_WIDGET(cs));
	gtk_widget_queue_draw(GTK_WIDGET(cs));

	return TRUE;
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

	cs->sample = calloc(COLSCALE_SAMPLE_SIZE, sizeof(float));
	if ( cs->sample == NULL ) return NULL;

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
	g_signal_connect(G_OBJECT(cs), "scroll-event",
	                 G_CALLBACK(scroll_sig), cs);
	g_signal_connect(G_OBJECT(cs), "draw",
	                 G_CALLBACK(draw_sig), cs);

	gtk_widget_set_can_focus(GTK_WIDGET(cs), TRUE);
	gtk_widget_add_events(GTK_WIDGET(cs),
	                      GDK_POINTER_MOTION_HINT_MASK
	                       | GDK_BUTTON1_MOTION_MASK
	                       | GDK_BUTTON_PRESS_MASK
	                       | GDK_BUTTON_RELEASE_MASK
	                       | GDK_SCROLL_MASK
	                       | GDK_SMOOTH_SCROLL_MASK
	                       | GDK_KEY_PRESS_MASK
	                       | GDK_KEY_RELEASE_MASK);

	gtk_widget_grab_focus(GTK_WIDGET(cs));

	gtk_widget_show(GTK_WIDGET(cs));

	return GTK_WIDGET(cs);
}


void crystfel_colour_scale_auto_range(CrystFELColourScale *cs)
{
	double mean;
	double variance;

	mean = gsl_stats_float_mean(cs->sample, 1, cs->n_samples);
	variance = gsl_stats_float_variance_m(cs->sample, 1,
	                                      cs->n_samples, mean);

	cs->lo = 0.0;
	cs->hi = mean + 10.0*sqrt(variance);

	make_histogram(cs);
}


void crystfel_colour_scale_scan_image(CrystFELColourScale *cs,
                                      struct image *image)
{
	int i;
	int pn;
	long int n_pix;

	if ( image == NULL ) return;

	n_pix = 0;
	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {
		int w, h;
		w = image->detgeom->panels[pn].w;
		h = image->detgeom->panels[pn].h;
		for ( i=0; i<w*h; i++ ) {
			if ( !image->bad[pn][i] ) n_pix++;
		}
	}
	float p = (float)COLSCALE_SAMPLE_SIZE/n_pix;

	cs->n_samples = 0;
	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {
		int w, h;
		long int i;
		w = image->detgeom->panels[pn].w;
		h = image->detgeom->panels[pn].h;
		for ( i=0; i<w*h; i++ ) {
			if ( !image->bad[pn][i] ) {
				if ( rand() < p*RAND_MAX ) {
					cs->sample[cs->n_samples++] = image->dp[pn][i];
				}
			}
			if ( cs->n_samples == COLSCALE_SAMPLE_SIZE ) break;
		}
		if ( cs->n_samples == COLSCALE_SAMPLE_SIZE ) break;
	}

	make_histogram(cs);
}
