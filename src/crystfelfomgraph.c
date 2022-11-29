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


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELFoMGraph *fg)
{
	cairo_save(cr);

	/* Overall background */
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_paint(cr);


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
