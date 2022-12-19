/*
 * crystfelfomgraph.c
 *
 * Figure of merit graph plot widget
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

#include <fom.h>

#include "crystfelfomgraph.h"


G_DEFINE_TYPE(CrystFELFoMGraph, crystfel_fom_graph,
              GTK_TYPE_DRAWING_AREA)

static gint destroy_sig(GtkWidget *window, CrystFELFoMGraph *fg)
{
	return FALSE;
}


static gint configure_sig(GtkWidget *window, GdkEventConfigure *rec,
                          CrystFELFoMGraph *fg)
{
	fg->visible_width = rec->width;
	fg->visible_height = rec->height;
	return FALSE;
}


static void fom_colour(enum fom_type t, double *col)
{
	switch ( t ) {
		case FOM_R1I :
		col[0] = 1.0;  col[1] = 0.0;  col[2] = 0.0; break;
		case FOM_R1F:
		col[0] = 1.0;  col[1] = 0.0;  col[2] = 0.0; break;
		case FOM_R2:
		col[0] = 1.0;  col[1] = 0.0;  col[2] = 0.0; break;
		case FOM_RSPLIT:
		col[0] = 1.0;  col[1] = 0.0;  col[2] = 0.0; break;
		case FOM_CC:
		col[0] = 0.0;  col[1] = 0.0;  col[2] = 0.3; break;
		case FOM_CCSTAR:
		col[0] = 0.0;  col[1] = 0.0;  col[2] = 0.5; break;
		case FOM_CCANO:
		col[0] = 0.0;  col[1] = 0.0;  col[2] = 1.0; break;
		case FOM_CRDANO:
		col[0] = 0.5;  col[1] = 0.0;  col[2] = 0.5; break;
		case FOM_RANO:
		col[0] = 0.5;  col[1] = 0.2;  col[2] = 0.0; break;
		case FOM_RANORSPLIT:
		col[0] = 0.8;  col[1] = 0.2;  col[2] = 0.0; break;
		case FOM_D1SIG:
		col[0] = 0.0;  col[1] = 0.3;  col[2] = 0.0; break;
		case FOM_D2SIG:
		col[0] = 1.0;  col[1] = 0.5;  col[2] = 0.0; break;
		case FOM_REDUNDANCY:
		col[0] = 0.0;  col[1] = 1.0;  col[2] = 0.0; break;
		case FOM_SNR:
		col[0] = 0.5;  col[1] = 0.5;  col[2] = 0.0; break;
		case FOM_COMPLETENESS:
		col[0] = 0.3;  col[1] = 0.3;  col[2] = 0.3; break;
		default:
		col[0] = 1.0;  col[1] = 0.0;  col[2] = 0.0; break;
	}
}


static void fom_range(enum fom_type t, double *min, double *max)
{
	switch ( t ) {
		case FOM_R1I:              *min = 0.0;   *max = 0.8;   break;
		case FOM_R1F:              *min = 0.0;   *max = 0.8;   break;
		case FOM_R2:               *min = 0.0;   *max = 0.8;   break;
		case FOM_RSPLIT:           *min = 0.0;   *max = 0.8;   break;
		case FOM_CC:               *min = 0.9;   *max = 1.0;   break;
		case FOM_CCSTAR:           *min = 0.9;   *max = 1.0;   break;
		case FOM_CCANO:            *min = 0.9;   *max = 1.0;   break;
		case FOM_CRDANO:           *min = 0.9;   *max = 1.0;   break;
		case FOM_RANO:             *min = 0.0;   *max = 0.7;   break;
		case FOM_RANORSPLIT:       *min = 0.0;   *max = 8.0;   break;
		case FOM_D1SIG:            *min = 0.0;   *max = 1.0;   break;
		case FOM_D2SIG:            *min = 0.0;   *max = 1.0;   break;
		case FOM_REDUNDANCY:       *min = 0.0;   *max = 100;   break;
		case FOM_SNR:              *min = 0.0;   *max = 10.0;  break;
		case FOM_COMPLETENESS:     *min = 0.0;   *max = 1.0;   break;
		default:                   *min = 0.0;   *max = 1.0;   break;
	}
}


static void draw_x_axis(cairo_t *cr, double *tics, int n_tics,
                        double x1, double x2, double ox, double w, double h,
                        double axsp)
{
	int i;

	cairo_new_path(cr);
	cairo_move_to(cr, ox, axsp);
	cairo_line_to(cr, w, axsp);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);

	for ( i=0; i<15; i++ ) {

		cairo_text_extents_t ext;
		char label[128];
		double x;

		if ( (1e10/tics[i] > x1) && (1e10/tics[i] < x2) ) {

			x = ox + (w-ox) * (1e10/tics[i]-x1)/(x2-x1);

			cairo_move_to(cr, x, h);
			cairo_line_to(cr, x, axsp-5.0);
			cairo_set_line_width(cr, 1.0);
			cairo_stroke(cr);

			cairo_save(cr);
			snprintf(label, 127, "%.1f A", tics[i]);
			cairo_text_extents(cr, label, &ext);
			cairo_move_to(cr, x-ext.x_advance/2, -3.0);
			cairo_scale(cr, 1.0, -1.0);
			cairo_show_text(cr, label);
			cairo_restore(cr);

		}
	}
}


static void draw_y_axis(cairo_t *cr, double h, double oy,
                        double y1, double y2, const char *name)
{
	int i;
	cairo_text_extents_t ext;

	cairo_new_path(cr);
	cairo_move_to(cr, 0, oy);
	cairo_line_to(cr, 0, h);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);

	for ( i=0; i<=5; i++ ) {

		char label[128];
		double y_pos;
		double y_val;

		y_val = y1 + i*(y2-y1)/5.0;
		y_pos = oy + (h-oy) * (y_val-y1)/(y2-y1);

		cairo_move_to(cr, 0.0, y_pos);
		cairo_line_to(cr, -5.0, y_pos);
		cairo_set_line_width(cr, 1.0);
		cairo_stroke(cr);

		cairo_save(cr);
		snprintf(label, 127, "%.2f", y_val);
		cairo_text_extents(cr, label, &ext);
		cairo_move_to(cr, -ext.width-5.0, y_pos-ext.height/2.0);
		cairo_scale(cr, 1.0, -1.0);
		cairo_show_text(cr, label);
		cairo_restore(cr);

	}

	cairo_save(cr);
	cairo_text_extents(cr, name, &ext);
	cairo_move_to(cr, -ext.height-30.0, oy+(h-oy+ext.x_advance)/2.0);
	cairo_scale(cr, 1.0, -1.0);
	cairo_rotate(cr, M_PI_2);
	cairo_show_text(cr, name);
	cairo_restore(cr);
}


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELFoMGraph *fg)
{
	int j;

	cairo_save(cr);

	/* Overall background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);

	if ( fg->n_foms == 0 ) {
		cairo_restore(cr);
		return FALSE;
	}

	double border = 10.0;
	double w = fg->visible_width;
	double h = fg->visible_height;
	double axsp = 20.0;   /* Vertical space for x-axis labels */
	double yaxsp = 50.0;  /* Horizontal space for y-axis labels */
	double x1 = fg->shell_centers[0];
	double x2 = fg->shell_centers[fg->n_shells-1];

	/* Logical coordinates */
	cairo_translate(cr, 0.0, h);
	cairo_scale(cr, 1.0, -1.0);

	/* Add empty border */
	cairo_translate(cr, border, border);
	w -= border*2.0;
	h -= border*2.0;

	/* y-axes (multiple) */
	double ox = 0.0;
	cairo_save(cr);
	cairo_translate(cr, yaxsp, 0.0);
	for ( j=0; j<fg->n_foms; j++ ) {
		double col[3];
		double y1, y2;
		fom_colour(fg->fom_types[j], col);
		fom_range(fg->fom_types[j], &y1, &y2);
		cairo_set_source_rgb(cr, col[0], col[1], col[2]);
		draw_y_axis(cr, h, axsp, y1, y2, fom_name(fg->fom_types[j]));
		ox += yaxsp;
		cairo_translate(cr, yaxsp, 0.0);
	}
	cairo_restore(cr);

	/* x-axis */
	double tics[] = {20.0, 15.0, 10.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.7, 1.5,
		         1.4, 1.3, 1.2, 1.1, 1.0};
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	draw_x_axis(cr, tics, 15, x1, x2, ox, w, h, axsp);

	cairo_translate(cr, ox, axsp);
	w -= ox;
	h -= axsp;
	cairo_rectangle(cr, 0, 0, w, h);
	cairo_set_source_rgba(cr, 0.9, 0.9, 0.9, 0.9);
	cairo_clip_preserve(cr);
	cairo_fill(cr);
	for ( j=0; j<fg->n_foms; j++ ) {
		int i;
		double col[3];
		double y1, y2;
		int split = 1;
		fom_range(fg->fom_types[j], &y1, &y2);
		for ( i=0; i<fg->n_shells; i++ ) {
			double fv = fg->fom_values[j][i];
			if ( !isnan(fv) ) {
				if ( split ) {
					cairo_move_to(cr, w*(fg->shell_centers[i]-x1)/(x2-x1),
					              h*(fv-y1)/(y2-y1));
					split = 0;
				} else {
					cairo_line_to(cr, w*(fg->shell_centers[i]-x1)/(x2-x1),
					              h*(fv-y1)/(y2-y1));
				}
			} else {
				split = 1;
			}
		}
		cairo_set_line_width(cr, 1.0);
		fom_colour(fg->fom_types[j], col);
		cairo_set_source_rgb(cr, col[0], col[1], col[2]);
		cairo_stroke(cr);
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
	*natural = 480;
}


static void get_preferred_height(GtkWidget *widget, gint *min, gint *natural)
{
	*min = 0;
	*natural = 320;
}


static void crystfel_fom_graph_class_init(CrystFELFoMGraphClass *klass)
{
	GTK_WIDGET_CLASS(klass)->get_request_mode = get_request_mode;
	GTK_WIDGET_CLASS(klass)->get_preferred_width = get_preferred_width;
	GTK_WIDGET_CLASS(klass)->get_preferred_height = get_preferred_height;
	GTK_WIDGET_CLASS(klass)->get_preferred_height_for_width = NULL;
}


static void crystfel_fom_graph_init(CrystFELFoMGraph *fg)
{
}


GtkWidget *crystfel_fom_graph_new()
{
	CrystFELFoMGraph *fg;

	fg = g_object_new(CRYSTFEL_TYPE_FOM_GRAPH, NULL);

	fg->n_shells = 0;
	fg->shell_centers = NULL;
	fg->n_foms = 0;
	fg->fom_types = NULL;
	fg->fom_values = NULL;

	g_signal_connect(G_OBJECT(fg), "destroy",
	                 G_CALLBACK(destroy_sig), fg);
	g_signal_connect(G_OBJECT(fg), "configure-event",
	                 G_CALLBACK(configure_sig), fg);
	g_signal_connect(G_OBJECT(fg), "draw",
	                 G_CALLBACK(draw_sig), fg);

	gtk_widget_set_can_focus(GTK_WIDGET(fg), FALSE);

	gtk_widget_show(GTK_WIDGET(fg));

	return GTK_WIDGET(fg);
}


void crystfel_fom_graph_set_data(CrystFELFoMGraph *fg,
                                 double *shell_centers, int n_shells,
                                 enum fom_type *fom_types, double **fom_values,
                                 int n_foms)
{
	int i;
	for ( i=0; i<fg->n_foms; i++ ) {
		free(fg->fom_values[i]);
	}
	free(fg->shell_centers);
	free(fg->fom_types);
	free(fg->fom_values);

	fg->n_shells = n_shells;
	fg->shell_centers = shell_centers;
	fg->n_foms = n_foms;
	fg->fom_types = fom_types;
	fg->fom_values= fom_values;

	gtk_widget_queue_draw(GTK_WIDGET(fg));
}
