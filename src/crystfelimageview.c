/*
 * crystfelimageview.c
 *
 * CrystFEL's image viewer widget
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020 Thomas White <taw@physics.org>
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

#include "crystfelimageview.h"


static void scroll_interface_init(GtkScrollable *iface)
{
}


enum
{
	CRYSTFELIMAGEVIEW_0,
	CRYSTFELIMAGEVIEW_VADJ,
	CRYSTFELIMAGEVIEW_HADJ,
	CRYSTFELIMAGEVIEW_VPOL,
	CRYSTFELIMAGEVIEW_HPOL,
};


G_DEFINE_TYPE_WITH_CODE(CrystFELImageView, crystfel_image_view,
                        GTK_TYPE_DRAWING_AREA,
                        G_IMPLEMENT_INTERFACE(GTK_TYPE_SCROLLABLE,
                                              scroll_interface_init))


static gint destroy_sig(GtkWidget *window, CrystFELImageView *iv)
{
	return FALSE;
}


static gint realise_sig(GtkWidget *window, CrystFELImageView *iv)
{
	return FALSE;
}


static gint button_press_sig(GtkWidget *window, GdkEventButton *event,
                             CrystFELImageView *iv)
{
	return FALSE;
}


static gint motion_sig(GtkWidget *window, GdkEventMotion *event,
                       CrystFELImageView *iv)
{
	return FALSE;
}



static gint resize_sig(GtkWidget *window, GdkEventConfigure *event,
                       CrystFELImageView *iv)
{
	return FALSE;
}


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELImageView *iv)
{
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	cairo_paint(cr);
	return FALSE;
}


static void redraw(CrystFELImageView *iv)
{
}


static void horizontal_adjust(GtkAdjustment *adj, CrystFELImageView *iv)
{
	iv->x_scroll_pos = gtk_adjustment_get_value(adj);
	redraw(iv);
}


static void set_horizontal_params(CrystFELImageView *iv)
{
	if ( iv->hadj == NULL ) return;
	gtk_adjustment_configure(iv->hadj, iv->x_scroll_pos, 0, iv->w, 100,
	                         iv->visible_width, iv->visible_width);
}


static void vertical_adjust(GtkAdjustment *adj, CrystFELImageView *iv)
{
	iv->y_scroll_pos = gtk_adjustment_get_value(adj);
	redraw(iv);
}


static void set_vertical_params(CrystFELImageView *iv)
{
	if ( iv->vadj == NULL ) return;
	gtk_adjustment_configure(iv->vadj, iv->y_scroll_pos, 0, iv->w, 100,
	                         iv->visible_width, iv->visible_width);
}


static void crystfel_image_view_set_property(GObject *obj, guint id, const GValue *val,
                                             GParamSpec *spec)
{
	CrystFELImageView *iv = CRYSTFEL_IMAGE_VIEW(obj);

	switch ( id ) {

		case CRYSTFELIMAGEVIEW_VPOL :
		iv->vpol = g_value_get_enum(val);
		break;

		case CRYSTFELIMAGEVIEW_HPOL :
		iv->hpol = g_value_get_enum(val);
		break;

		case CRYSTFELIMAGEVIEW_VADJ :
		iv->vadj = g_value_get_object(val);
		set_vertical_params(iv);
		if ( iv->vadj != NULL ) {
			g_signal_connect(G_OBJECT(iv->vadj), "value-changed",
			                 G_CALLBACK(vertical_adjust), iv);
		}
		break;

		case CRYSTFELIMAGEVIEW_HADJ :
		iv->hadj = g_value_get_object(val);
		set_horizontal_params(iv);
		if ( iv->hadj != NULL ) {
			g_signal_connect(G_OBJECT(iv->hadj), "value-changed",
			                 G_CALLBACK(horizontal_adjust), iv);
		}
		break;

		default :
		printf("setting %i\n", id);
		break;

	}
}


static void crystfel_image_view_get_property(GObject *obj, guint id, GValue *val,
                                             GParamSpec *spec)
{
	CrystFELImageView *iv = CRYSTFEL_IMAGE_VIEW(obj);

	switch ( id ) {

		case CRYSTFELIMAGEVIEW_VADJ :
		g_value_set_object(val, iv->vadj);
		break;

		case CRYSTFELIMAGEVIEW_HADJ :
		g_value_set_object(val, iv->hadj);
		break;

		case CRYSTFELIMAGEVIEW_VPOL :
		g_value_set_enum(val, iv->vpol);
		break;

		case CRYSTFELIMAGEVIEW_HPOL :
		g_value_set_enum(val, iv->hpol);
		break;

		default :
		G_OBJECT_WARN_INVALID_PROPERTY_ID(obj, id, spec);
		break;

	}
}


static GtkSizeRequestMode get_request_mode(GtkWidget *widget)
{
	return GTK_SIZE_REQUEST_CONSTANT_SIZE;
}


static void get_preferred_width(GtkWidget *widget, gint *min, gint *natural)
{
	*min = 0;
	*natural = 640;
}


static void get_preferred_height(GtkWidget *widget, gint *min, gint *natural)
{
	*min = 0;
	*natural = 640;
}


static void crystfel_image_view_class_init(CrystFELImageViewClass *klass)
{
	GObjectClass *goc = G_OBJECT_CLASS(klass);
	goc->set_property = crystfel_image_view_set_property;
	goc->get_property = crystfel_image_view_get_property;
	g_object_class_override_property(goc, CRYSTFELIMAGEVIEW_VADJ, "vadjustment");
	g_object_class_override_property(goc, CRYSTFELIMAGEVIEW_HADJ, "hadjustment");
	g_object_class_override_property(goc, CRYSTFELIMAGEVIEW_VPOL, "vscroll-policy");
	g_object_class_override_property(goc, CRYSTFELIMAGEVIEW_HPOL, "hscroll-policy");

	GTK_WIDGET_CLASS(klass)->get_request_mode = get_request_mode;
	GTK_WIDGET_CLASS(klass)->get_preferred_width = get_preferred_width;
	GTK_WIDGET_CLASS(klass)->get_preferred_height = get_preferred_height;
	GTK_WIDGET_CLASS(klass)->get_preferred_height_for_width = NULL;
}


static void crystfel_image_view_init(CrystFELImageView *iv)
{
	iv->vpol = GTK_SCROLL_NATURAL;
	iv->hpol = GTK_SCROLL_NATURAL;
	iv->vadj = NULL;
	iv->hadj = NULL;
}


GtkWidget *crystfel_image_view_new()
{
	CrystFELImageView *iv;

	iv = g_object_new(CRYSTFEL_TYPE_IMAGE_VIEW, NULL);

	iv->w = 100;
	iv->h = 100;
	iv->x_scroll_pos = 0;
	iv->y_scroll_pos = 0;

	gtk_widget_set_size_request(GTK_WIDGET(iv), iv->w, iv->h);

	g_signal_connect(G_OBJECT(iv), "destroy",
	                 G_CALLBACK(destroy_sig), iv);
	g_signal_connect(G_OBJECT(iv), "realize",
	                 G_CALLBACK(realise_sig), iv);
	g_signal_connect(G_OBJECT(iv), "button-press-event",
	                 G_CALLBACK(button_press_sig), iv);
	g_signal_connect(G_OBJECT(iv), "motion-notify-event",
	                 G_CALLBACK(motion_sig), iv);
	g_signal_connect(G_OBJECT(iv), "configure-event",
	                 G_CALLBACK(resize_sig), iv);
	g_signal_connect(G_OBJECT(iv), "draw",
	                 G_CALLBACK(draw_sig), iv);

	gtk_widget_set_can_focus(GTK_WIDGET(iv), TRUE);
	gtk_widget_add_events(GTK_WIDGET(iv),
	                      GDK_POINTER_MOTION_HINT_MASK
	                       | GDK_BUTTON1_MOTION_MASK
	                       | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK
	                       | GDK_KEY_PRESS_MASK | GDK_KEY_RELEASE_MASK
	                       | GDK_SCROLL_MASK);

	gtk_widget_grab_focus(GTK_WIDGET(iv));

	gtk_widget_show(GTK_WIDGET(iv));

	return GTK_WIDGET(iv);
}
