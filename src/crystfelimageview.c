/*
 * crystfelimageview.c
 *
 * CrystFEL's image viewer widget
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2021 Thomas White <taw@physics.org>
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

#include <utils.h>
#include <detgeom.h>
#include <colscale.h>

#include "crystfelimageview.h"


static int rerender_image(CrystFELImageView *iv);


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

static void redraw(CrystFELImageView *iv)
{
	gint w, h;
	w = gtk_widget_get_allocated_width(GTK_WIDGET(iv));
	h = gtk_widget_get_allocated_height(GTK_WIDGET(iv));
	gtk_widget_queue_draw_area(GTK_WIDGET(iv), 0, 0, w, h);
}


static void cleanup_image(CrystFELImageView *iv)
{
	if ( iv->pixbufs != NULL ) {
		int i;
		for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
			if ( iv->pixbufs[i] != NULL ) {
				gdk_pixbuf_unref(iv->pixbufs[i]);
			}
		}
	}
	free(iv->pixbufs);

	iv->image = NULL;
	iv->pixbufs = NULL;
}


static gint destroy_sig(GtkWidget *window, CrystFELImageView *iv)
{
	cleanup_image(iv);
	return FALSE;
}


static gint realise_sig(GtkWidget *window, CrystFELImageView *iv)
{
	return FALSE;
}


static void configure_scroll_adjustments(CrystFELImageView *iv)
{
	if ( iv->hadj != NULL ) {
		double pos = gtk_adjustment_get_value(iv->hadj);
		double vis_size = iv->visible_width / iv->zoom;
		gtk_adjustment_configure(iv->hadj, pos,
		                         iv->min_x, iv->min_x+iv->detector_w*1.1,
		                         0.0001, 0.1, vis_size);
	}
	if ( iv->vadj != NULL ) {
		double pos = gtk_adjustment_get_value(iv->vadj);
		double vis_size = iv->visible_height / iv->zoom;
		gtk_adjustment_configure(iv->vadj, pos,
		                         -iv->max_y,
		                         -(iv->max_y-iv->detector_h*1.1),
		                         0.0001, 0.1, vis_size);
	}
}


static void handle_scroll_click(double zoom_scale,
                                CrystFELImageView *iv,
                                GdkEventScroll *event)
{
	double ratio;
	int zoom_allowed = 1;

	/* Size of a detector pixel in screen pixels */
	ratio = iv->zoom * iv->image->detgeom->panels[0].pixel_pitch;

	if ( (ratio < 0.05) && (zoom_scale < 1.0) ) zoom_allowed = 0;
	if ( (ratio > 100.0) && (zoom_scale > 1.0) ) zoom_allowed = 0;

	if ( zoom_allowed ) {

		double shift_x, shift_y;
		double scr_x, scr_y;

		scr_x = gtk_adjustment_get_value(iv->hadj);
		scr_y = gtk_adjustment_get_value(iv->vadj);

		shift_x = event->x*((1.0/(zoom_scale*iv->zoom))-(1.0/iv->zoom));
		shift_y = event->y*((1.0/(zoom_scale*iv->zoom))-(1.0/iv->zoom));
		iv->zoom *= zoom_scale;

		configure_scroll_adjustments(iv);

		gtk_adjustment_set_value(iv->hadj, scr_x - shift_x);
		gtk_adjustment_set_value(iv->vadj, scr_y - shift_y);

		redraw(iv);

	}
}


static gint scroll_sig(GtkWidget *window, GdkEventScroll *event,
                       CrystFELImageView *iv)
{
	double xs, ys;

	if ( iv->image == NULL ) return FALSE;

	switch ( event->direction ) {

		case GDK_SCROLL_SMOOTH:
		if ( gdk_event_get_scroll_deltas((GdkEvent *)event, &xs, &ys) ) {
			handle_scroll_click(1.0-ys*0.1, iv, event);
		}
		return TRUE;

		case GDK_SCROLL_UP:
		handle_scroll_click(1.1, iv, event);
		return TRUE;

		case GDK_SCROLL_DOWN:
		handle_scroll_click(0.9, iv, event);
		return TRUE;

		case GDK_SCROLL_LEFT:
		return TRUE;  /* Do not propagate further */

		case GDK_SCROLL_RIGHT:
		return TRUE;  /* Do not propagate further */

		default:
		printf("Unhandled scroll direction %i\n", event->direction);
		return FALSE;
	}
}


static gint button_press_sig(GtkWidget *window, GdkEventButton *event,
                             CrystFELImageView *iv)
{
	iv->drag_start_x = event->x;
	iv->drag_start_y = event->y;
	iv->drag_start_sp_x = gtk_adjustment_get_value(iv->hadj);
	iv->drag_start_sp_y = gtk_adjustment_get_value(iv->vadj);
	return FALSE;
}


static gint motion_sig(GtkWidget *window, GdkEventMotion *event,
                       CrystFELImageView *iv)
{
	double ddx, ddy;
	ddx = event->x - iv->drag_start_x;
	ddy = event->y - iv->drag_start_y;
	gtk_adjustment_set_value(iv->hadj, iv->drag_start_sp_x - ddx/iv->zoom);
	gtk_adjustment_set_value(iv->vadj, iv->drag_start_sp_y - ddy/iv->zoom);
	redraw(iv);
	return FALSE;
}


static gint configure_sig(GtkWidget *window, GdkEventConfigure *rec,
                          CrystFELImageView *iv)
{
	iv->visible_width = rec->width;
	iv->visible_height = rec->height;
	configure_scroll_adjustments(iv);
	return FALSE;
}


static void draw_pixel_values(cairo_t *cr,
                              int min_fs, int max_fs,
                              int min_ss, int max_ss,
                              struct detgeom_panel p, float *dp,
                              int *bad)
{
	int fs, ss;
	PangoLayout *layout;
	PangoFontDescription *fontdesc;
	double w, h;

	layout = pango_cairo_create_layout(cr);
	fontdesc = pango_font_description_from_string("Sans 10");
	pango_layout_set_font_description(layout, fontdesc);

	w = 1.0;
	h = 1.0;
	cairo_device_to_user_distance(cr, &w, &h);

	for ( fs=min_fs; fs<=max_fs; fs++ ) {
	for ( ss=min_ss; ss<=max_ss; ss++ ) {

		double x, y;
		double fsd, ssd;
		char tmp[64];
		PangoRectangle rec;
		double rw, rh;
		const char *b1;
		const char *b2;

		fsd = fs + 0.5;
		ssd = ss + 0.5;

		x = (fsd*p.fsx + ssd*p.ssx + p.cnx)*p.pixel_pitch;
		y = (fsd*p.fsy + ssd*p.ssy + p.cny)*p.pixel_pitch;

		if ( bad[fs+ss*p.w] ) {
			b1 = "(";
			b2 = ")";
		} else {
			b1 = "";
			b2 = "";
		}

		snprintf(tmp, 63, "%s%.f%s", b1, dp[fs+ss*p.w], b2);
		pango_layout_set_text(layout, tmp, -1);

		cairo_save(cr);

		cairo_translate(cr, x, y);
		cairo_scale(cr, w, h);
		pango_cairo_update_layout(cr, layout);
		pango_layout_get_extents(layout, NULL, &rec);
		rw = pango_units_to_double(rec.width);
		rh = pango_units_to_double(rec.height);

		cairo_rectangle(cr, -rw/2.0, -rh/2.0, rw, rh);
		cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
		cairo_fill(cr);

		cairo_move_to(cr, -rw/2.0, -rh/2.0);
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		pango_cairo_show_layout(cr, layout);

		cairo_restore(cr);

	}
	}

	pango_font_description_free(fontdesc);
	g_object_unref(layout);

	cairo_arc(cr, p.cnx*p.pixel_pitch, p.cny*p.pixel_pitch, 0.00002,
	          0, 2.0*M_PI);
	cairo_set_source_rgba(cr, 0.0, 1.0, 1.0, 1.0);
	cairo_fill(cr);
}


static void check_pixel_visibility(cairo_matrix_t *gtkmatrix,
                                   cairo_t *cr,
                                   struct detgeom_panel *p,
                                   double x, double y,
                                   int *min_fs, int *max_fs,
                                   int *min_ss, int *max_ss)
{
	/* Take into account the transformation which GTK already
	 * has in effect on the widget (translation from parent
	 * window's origin). */
	cairo_matrix_transform_point(gtkmatrix, &x, &y);
	cairo_device_to_user(cr, &x, &y);

	x /= p->pixel_pitch;
	y /= p->pixel_pitch;
	/* x,y is now fs,ss for the current panel */

	if ( x < 0 ) x = 0;
	if ( x > p->w ) x = p->w-1;
	if ( y < 0 ) y = 0;
	if ( y > p->h ) y = p->h-1;

	if ( (x > *max_fs) ) *max_fs = x;
	if ( (x < *min_fs) ) *min_fs = x;
	if ( (y > *max_ss) ) *max_ss = y;
	if ( (y < *min_ss) ) *min_ss = y;
}


static void draw_panel_rectangle(cairo_t *cr, CrystFELImageView *iv,
                                 int i, cairo_matrix_t *gtkmatrix)
{
	struct detgeom_panel p = iv->image->detgeom->panels[i];
	cairo_matrix_t m;
	cairo_pattern_t *patt;
	double xs, ys, pixel_size_on_screen;
	int have_pixels = 1;

	cairo_save(cr);

	/* Move to the right location */
	cairo_translate(cr, p.cnx*p.pixel_pitch, p.cny*p.pixel_pitch);

	/* Twiddle directions according to matrix */
	cairo_matrix_init(&m, p.fsx, p.fsy, p.ssx, p.ssy,
	                      0.0, 0.0);
	cairo_transform(cr, &m);

	gdk_cairo_set_source_pixbuf(cr, iv->pixbufs[i], 0.0, 0.0);
	patt = cairo_get_source(cr);

	cairo_pattern_get_matrix(patt, &m);
	cairo_matrix_scale(&m, 1.0/p.pixel_pitch, 1.0/p.pixel_pitch);
	cairo_pattern_set_matrix(patt, &m);

	cairo_pattern_set_filter(patt, CAIRO_FILTER_NEAREST);

	cairo_rectangle(cr, 0.0, 0.0, p.w*p.pixel_pitch, p.h*p.pixel_pitch);
	cairo_fill_preserve(cr);
	cairo_set_line_width(cr, 0.00001);
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	cairo_stroke(cr);

	/* Are any pixels from this panel visible? */
	int min_fs = p.w;
	int max_fs = 0;
	int min_ss = p.h;
	int max_ss = 0;
	check_pixel_visibility(gtkmatrix, cr, &p, 0.0, 0.0,
	                       &min_fs, &max_fs, &min_ss, &max_ss);
	check_pixel_visibility(gtkmatrix, cr, &p, 0.0, iv->visible_height,
	                       &min_fs, &max_fs, &min_ss, &max_ss);
	check_pixel_visibility(gtkmatrix, cr, &p, iv->visible_width, 0.0,
	                       &min_fs, &max_fs, &min_ss, &max_ss);
	check_pixel_visibility(gtkmatrix, cr, &p, iv->visible_width, iv->visible_height,
	                       &min_fs, &max_fs, &min_ss, &max_ss);

	cairo_restore(cr);

	xs = p.pixel_pitch;
	ys = p.pixel_pitch;
	cairo_user_to_device_distance(cr, &xs, &ys);
	pixel_size_on_screen = smallest(fabs(xs), fabs(ys));
	if ( (pixel_size_on_screen > 40.0) && have_pixels ) {
		draw_pixel_values(cr, min_fs, max_fs, min_ss, max_ss, p,
		                  iv->image->dp[i], iv->image->bad[i]);
	}
}


static void draw_peaks(cairo_t *cr, CrystFELImageView *iv,
                       ImageFeatureList *pks)
{
	int i, n_pks;
	double bs, lw;

	bs = 5.0;  /* Box size in pixels on screen */
	lw = 1.0;  /* Line width in pixels on screen */
	cairo_device_to_user_distance(cr, &bs, &lw);
	bs = fabs(bs);
	lw = fabs(lw);

	n_pks = image_feature_count(pks);
	for ( i=0; i<n_pks; i++ ) {

		const struct imagefeature *f;
		struct detgeom_panel *p;
		double x, y;
		double this_bs;
		double this_lw;
		int show_cen = 0;

		f = image_get_feature_const(pks, i);
		if ( f == NULL ) continue;
		p = &iv->image->detgeom->panels[f->pn];

		this_lw = biggest(0.1*p->pixel_pitch, lw);
		this_bs = biggest(iv->peak_box_size * p->pixel_pitch,
		                  bs);

		if ( this_bs > bs ) {
			show_cen = 1;
		}

		x = p->pixel_pitch*(p->cnx + p->fsx*f->fs + p->ssx*f->ss);
		y = p->pixel_pitch*(p->cny + p->fsy*f->fs + p->ssy*f->ss);
		cairo_rectangle(cr, x-this_bs, y-this_bs, 2*this_bs, 2*this_bs);
		cairo_set_line_width(cr, this_lw);
		cairo_set_source_rgba(cr, 1.0, 1.0, 0.0, 0.9);
		cairo_stroke(cr);

		if ( show_cen ) {
			cairo_move_to(cr, x-0.2*p->pixel_pitch, y);
			cairo_line_to(cr, x+0.2*p->pixel_pitch, y);
			cairo_move_to(cr, x, y-0.2*p->pixel_pitch);
			cairo_line_to(cr, x, y+0.2*p->pixel_pitch);
			cairo_set_source_rgba(cr, 0.4, 0.4, 0.0, 0.9);
			cairo_stroke(cr);
		}

	}
}


static void render_overlined_indices(cairo_t *dctx,
                                     signed int h,
                                     signed int k,
                                     signed int l)
{
	char tmp[256];
	cairo_text_extents_t size;
	double x, y;

	cairo_save(dctx);

	cairo_set_line_cap(dctx, CAIRO_LINE_CAP_ROUND);
	cairo_get_current_point(dctx, &x, &y);

	/* Draw 'h' */
	snprintf(tmp, 255, "%i", abs(h));
	cairo_text_extents(dctx, tmp, &size);
	cairo_set_line_width(dctx, size.height*0.1);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( h < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-size.height*1.1);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}
	x += size.x_advance*1.4;

	/* Draw 'k' */
	cairo_move_to(dctx, x, y);
	snprintf(tmp, 255, "%i", abs(k));
	cairo_text_extents(dctx, tmp, &size);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( k < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-size.height*1.1);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}
	x += size.x_advance*1.4;

	/* Draw 'l' */
	cairo_move_to(dctx, x, y);
	snprintf(tmp, 255, "%i", abs(l));
	cairo_text_extents(dctx, tmp, &size);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( l < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-size.height*1.1);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}

	cairo_restore(dctx);
}


static void draw_refls(cairo_t *cr,
                       CrystFELImageView *iv,
                       RefList *list,
                       int label_reflections,
                       double *colour)
{
	const Reflection *refl;
	RefListIterator *iter;
	double bs, lw;

	if ( list == NULL ) return;

	bs = 5.0;
	lw = 1.0;
	cairo_device_to_user_distance(cr, &bs, &lw);
	bs = fabs(bs);
	lw = fabs(lw);

	for ( refl = first_refl_const(list, &iter);
	      refl != NULL;
	      refl = next_refl_const(refl, iter) )
	{
		struct detgeom_panel *p;
		double fs, ss;
		int pn;
		double x, y;
		float this_bs;
		float this_lw;
		int show_cen = 0;

		get_detector_pos(refl, &fs, &ss);
		pn = get_panel_number(refl);
		p = &iv->image->detgeom->panels[pn];

		this_lw = biggest(0.1*p->pixel_pitch, lw);
		this_bs = biggest(iv->refl_box_size * p->pixel_pitch, bs);

		if ( this_bs > 2.5*bs ) {
			show_cen = 1;
		}

		x = p->pixel_pitch*(p->cnx + p->fsx*fs + p->ssx*ss);
		y = p->pixel_pitch*(p->cny + p->fsy*fs + p->ssy*ss);

		cairo_arc(cr, x, y, this_bs, 0, 2*M_PI);
		cairo_set_line_width(cr, this_lw);

		if ( get_redundancy(refl) == 0 ) {
			cairo_set_source_rgba(cr, 0.7, 0.0, 0.0, 0.9);
		} else {
			cairo_set_source_rgba(cr,
			                      colour[0],
			                      colour[1],
			                      colour[2],
			                      0.9);
		}
		cairo_stroke(cr);

		if ( show_cen ) {
			cairo_move_to(cr, x-0.2*p->pixel_pitch, y);
			cairo_line_to(cr, x+0.2*p->pixel_pitch, y);
			cairo_move_to(cr, x, y-0.2*p->pixel_pitch);
			cairo_line_to(cr, x, y+0.2*p->pixel_pitch);
			cairo_set_source_rgba(cr,
			                      colour[0]/2.0,
			                      colour[1]/2.0,
			                      colour[2]/2.0,
			                      0.9);
			cairo_stroke(cr);
		}

		if ( label_reflections ) {

			signed int h, k, l;

			get_indices(refl, &h, &k, &l);
			cairo_save(cr);
			cairo_move_to(cr, x, y);
			cairo_set_source_rgba(cr, 0.0, 0.4, 0.0, 0.9);
			cairo_set_font_size(cr, 11*p->pixel_pitch);
			cairo_scale(cr, 1.0, -1.0);
			render_overlined_indices(cr, h, k, l);
			cairo_restore(cr);

		}


	}
}


static double ring_radius(double d, double wl, double z)
{
	double theta = asin(wl / (2.0*d));
	return z * tan(2.0*theta);
}


static void show_ring(cairo_t *cr, double wl, double mean_z,
                      double d, const char *label,
                      double r, double g, double b)
{
	cairo_text_extents_t size;
	double bs, lw;
	double radius = ring_radius(d, wl, mean_z);

	if ( isnan(radius) ) return;

	bs = 17.0;
	lw = 1.0;
	cairo_device_to_user_distance(cr, &bs, &lw);
	bs = fabs(bs);
	lw = fabs(lw);

	cairo_save(cr);

	cairo_arc(cr, 0.0, 0.0, radius, 0.0, 2.0*M_PI);
	cairo_set_source_rgb(cr, r, g, b);
	cairo_set_line_width(cr, lw);
	cairo_stroke(cr);

	cairo_rotate(cr, -M_PI/4.0);
	cairo_scale(cr, 1.0, -1.0);
	cairo_set_font_size(cr, bs);
	cairo_text_extents(cr, label, &size);
	cairo_translate(cr, -size.width/2.0, radius-5.0*lw);
	cairo_show_text(cr, label);
	cairo_fill(cr);

	cairo_restore(cr);
}


static double crystal_cols[][3] =
{
	{0.0, 1.0, 0.0},   /* bright green */
	{0.0, 0.8, 0.8},   /* cyan */
	{1.0, 1.0, 0.0},   /* bright yellow */
	{1.0, 1.0, 1.0},   /* white */
};

static int n_crystal_cols = 4;


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELImageView *iv)
{
	cairo_matrix_t m;

	if ( iv->image == NULL ) return FALSE;
	if ( iv->need_rerender ) rerender_image(iv);

	cairo_save(cr);

	/* Overall background (light grey) */
	cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
	cairo_paint(cr);

	/* Get the transformation matrix before my transformations */
	cairo_get_matrix(cr, &m);

	cairo_scale(cr, iv->zoom, iv->zoom);
	cairo_scale(cr, 1.0, -1.0);

	cairo_translate(cr, -gtk_adjustment_get_value(iv->hadj),
	                     gtk_adjustment_get_value(iv->vadj));

	if ( iv->pixbufs != NULL ) {
		int i;
		for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
			cairo_save(cr);
			draw_panel_rectangle(cr, iv, i, &m);
			cairo_restore(cr);
		}
	}

	if ( iv->show_centre ) {
		cairo_arc(cr, 0.0, 0.0, 0.006, 0, 2.0*M_PI);
		cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 0.9);
		cairo_set_line_width(cr, 0.0001);
		cairo_stroke(cr);
		cairo_move_to(cr, -0.001, 0.0);
		cairo_line_to(cr, 0.001, 0.0);
		cairo_move_to(cr, 0.0, -0.001);
		cairo_line_to(cr, 0.0, 0.001);
		cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 0.9);
		cairo_stroke(cr);
	}

	if ( iv->show_peaks ) {
		draw_peaks(cr, iv, iv->image->features);
	}

	if ( iv->show_refls ) {
		int i;
		for ( i=0; i<iv->image->n_crystals; i++ ) {
			Crystal *cry = iv->image->crystals[i];
			draw_refls(cr, iv,
			           crystal_get_reflections(cry),
			           iv->label_refls,
			           crystal_cols[i % n_crystal_cols]);
		}
	}

	if ( iv->resolution_rings ) {
		double wl = iv->image->lambda;
		double mean_z = detgeom_mean_camera_length(iv->image->detgeom);
		if ( !isnan(mean_z) ) {
			show_ring(cr, wl, mean_z, 10.0e-10, "10A",   1.0, 0.0, 0.0);
			show_ring(cr, wl, mean_z, 9.0e-10,   "9A",   1.0, 0.0, 0.0);
			show_ring(cr, wl, mean_z, 8.0e-10,   "8A",   1.0, 0.0, 0.0);
			show_ring(cr, wl, mean_z, 7.0e-10,   "7A",   1.0, 0.5, 0.0);
			show_ring(cr, wl, mean_z, 6.0e-10,   "6A",   1.0, 1.0, 0.0);
			show_ring(cr, wl, mean_z, 5.0e-10,   "5A",   0.0, 1.0, 0.0);
			show_ring(cr, wl, mean_z, 4.0e-10,   "4A",   0.2, 1.0, 0.2);
			show_ring(cr, wl, mean_z, 3.0e-10,   "3A",   0.4, 1.0, 0.4);
			show_ring(cr, wl, mean_z, 2.0e-10,   "2A",   0.6, 1.0, 0.6);
			show_ring(cr, wl, mean_z, 1.0e-10,   "1A",   0.8, 1.0, 0.8);
			show_ring(cr, wl, mean_z, 0.5e-10, "0.5A",   1.0, 1.0, 1.0);
		}
	}

	cairo_restore(cr);

	return FALSE;
}


static void scroll_adjust_sig(GtkAdjustment *adj, CrystFELImageView *iv)
{
	redraw(iv);
}

static void set_adjustment_horizontal(CrystFELImageView *iv, GtkAdjustment *adj)
{
	if ( iv->hadj == adj ) return;

	if ( adj == NULL ) {
		adj = gtk_adjustment_new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}

	iv->hadj = adj;
	g_object_ref_sink(iv->hadj);

	configure_scroll_adjustments(iv);
	g_signal_connect(G_OBJECT(iv->hadj), "value-changed",
	                 G_CALLBACK(scroll_adjust_sig), iv);
}


static void set_adjustment_vertical(CrystFELImageView *iv, GtkAdjustment *adj)
{
	if ( iv->vadj == adj ) return;

	if ( adj == NULL ) {
		adj = gtk_adjustment_new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}

	iv->vadj = adj;
	g_object_ref_sink(iv->vadj);

	configure_scroll_adjustments(iv);
	g_signal_connect(G_OBJECT(iv->vadj), "value-changed",
	                 G_CALLBACK(scroll_adjust_sig), iv);
}


static void crystfel_image_view_set_property(GObject *obj, guint id,
                                             const GValue *val,
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
		set_adjustment_vertical(iv, g_value_get_object(val));
		break;

		case CRYSTFELIMAGEVIEW_HADJ :
		set_adjustment_horizontal(iv, g_value_get_object(val));
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

	/* All values initially meaningless */
	iv->detector_w = 1.0;
	iv->detector_h = 1.0;
	iv->zoom = -1.0;
	iv->image = NULL;
	iv->show_centre = 1;
	iv->show_peaks = 0;
	iv->brightness = 1.0;
	iv->pixbufs = NULL;
	iv->peak_box_size = 1.0;
	iv->refl_box_size = 1.0;
	iv->label_refls = 1;
	iv->need_rerender = 0;
	iv->need_recentre = 1;
	iv->resolution_rings = 0;
	iv->scale_lo = 0.0;
	iv->scale_hi = 100000.0;

	g_signal_connect(G_OBJECT(iv), "destroy",
	                 G_CALLBACK(destroy_sig), iv);
	g_signal_connect(G_OBJECT(iv), "realize",
	                 G_CALLBACK(realise_sig), iv);
	g_signal_connect(G_OBJECT(iv), "button-press-event",
	                 G_CALLBACK(button_press_sig), iv);
	g_signal_connect(G_OBJECT(iv), "scroll-event",
	                 G_CALLBACK(scroll_sig), iv);
	g_signal_connect(G_OBJECT(iv), "motion-notify-event",
	                 G_CALLBACK(motion_sig), iv);
	g_signal_connect(G_OBJECT(iv), "configure-event",
	                 G_CALLBACK(configure_sig), iv);
	g_signal_connect(G_OBJECT(iv), "draw",
	                 G_CALLBACK(draw_sig), iv);

	gtk_widget_set_can_focus(GTK_WIDGET(iv), TRUE);
	gtk_widget_add_events(GTK_WIDGET(iv),
	                      GDK_POINTER_MOTION_HINT_MASK
	                       | GDK_BUTTON1_MOTION_MASK
	                       | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK
	                       | GDK_KEY_PRESS_MASK | GDK_KEY_RELEASE_MASK
	                       | GDK_SCROLL_MASK | GDK_SMOOTH_SCROLL_MASK);

	gtk_widget_grab_focus(GTK_WIDGET(iv));

	gtk_widget_show(GTK_WIDGET(iv));

	return GTK_WIDGET(iv);
}


static void check_extents(struct detgeom_panel p, double *min_x, double *min_y,
                          double *max_x, double *max_y, double fs, double ss)
{
	double xs, ys;

	xs = (fs*p.fsx + ss*p.ssx + p.cnx) * p.pixel_pitch;
	ys = (fs*p.fsy + ss*p.ssy + p.cny) * p.pixel_pitch;

	if ( xs > *max_x ) *max_x = xs;
	if ( ys > *max_y ) *max_y = ys;
	if ( xs < *min_x ) *min_x = xs;
	if ( ys < *min_y ) *min_y = ys;
}


static void detgeom_pixel_extents(struct detgeom *det,
                                  double *min_x, double *min_y,
                                  double *max_x, double *max_y)
{
	int i;

	*min_x = 0.0;
	*max_x = 0.0;
	*min_y = 0.0;
	*max_y = 0.0;

	/* To determine the maximum extents of the detector, put all four
	 * corners of each panel through the transformations and watch for the
	 * biggest */

	for ( i=0; i<det->n_panels; i++ ) {

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              0.0, 0.0);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              0.0, det->panels[i].h+1);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              det->panels[i].w+1, 0.0);

		check_extents(det->panels[i], min_x, min_y, max_x, max_y,
		              det->panels[i].w+1, det->panels[i].h+1);


	}
}


static void free_pixbuf(guchar *data, gpointer p)
{
	free(data);
}


static GdkPixbuf *render_panel(float *data, int *badmap, int w, int h,
                               int scale_type, double scale_lo, double scale_hi)


{
	guchar *pixbuf_data;
	long int i;

	/* Rendered (colourful) version */
	pixbuf_data = malloc(3*w*h);
	if ( pixbuf_data == NULL ) return NULL;

	for ( i=0; i<w*h; i++ ) {

		double r, g, b;

		if ( !badmap[i] ) {

			colscale_lookup(data[i]-scale_lo, scale_hi-scale_lo,
			                scale_type, &r, &g, &b);

			pixbuf_data[3*i+0] = 255*r;
			pixbuf_data[3*i+1] = 255*g;
			pixbuf_data[3*i+2] = 255*b;

		} else {

			/* Bad pixel indicator colour */
			pixbuf_data[3*i+0] = 30;
			pixbuf_data[3*i+1] = 20;
			pixbuf_data[3*i+2] = 0;

		}

	}

	/* Create the pixbuf from the 8-bit display data */
	return gdk_pixbuf_new_from_data(pixbuf_data,
	                                GDK_COLORSPACE_RGB,
	                                FALSE, 8, w, h, w*3,
	                                free_pixbuf, NULL);

}


static void center_adjustment(GtkAdjustment *adj)
{
	double min = gtk_adjustment_get_lower(adj);
	double max = gtk_adjustment_get_upper(adj);
	double page = gtk_adjustment_get_page_size(adj);
	gtk_adjustment_set_value(adj, min+(max-min-page)/2.0);
}


static int rerender_image(CrystFELImageView *iv)
{
	int i;
	double min_x, min_y, max_x, max_y;

	if ( iv->image == NULL ) return 0;
	if ( iv->image->detgeom == NULL ) return 0;

	if ( iv->pixbufs == NULL ) {
		iv->pixbufs = calloc(iv->image->detgeom->n_panels,
		                     sizeof(GdkPixbuf *));
		if ( iv->pixbufs == NULL ) return 1;
	} else {
		for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
			gdk_pixbuf_unref(iv->pixbufs[i]);
		}
	}

	for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
		iv->pixbufs[i] = render_panel(iv->image->dp[i],
		                              iv->image->bad[i],
		                              iv->image->detgeom->panels[i].w,
		                              iv->image->detgeom->panels[i].h,
		                              SCALE_COLOUR,
		                              iv->scale_lo, iv->scale_hi);
		if ( iv->pixbufs[i] == NULL ) return 1;
	}

	detgeom_pixel_extents(iv->image->detgeom,
	                      &min_x, &min_y,
	                      &max_x, &max_y);
	iv->detector_w = max_x - min_x;
	iv->detector_h = max_y - min_y;
	iv->min_x = min_x - iv->detector_w*0.05;
	iv->max_y = max_y + iv->detector_h*0.05;
	if ( iv->zoom < 0.0 ) {
		iv->zoom = 1.0/iv->image->detgeom->panels[0].pixel_pitch;
	}
	configure_scroll_adjustments(iv);
	if ( iv->need_recentre ) {
		center_adjustment(iv->hadj);
		center_adjustment(iv->vadj);
		iv->need_recentre = 0;
	}

	iv->need_rerender = 0;
	redraw(iv);

	return 0;
}


int crystfel_image_view_set_image(CrystFELImageView *iv,
                                  const struct image *image)
{
	cleanup_image(iv);
	iv->image = image;
	iv->need_rerender = 1;
	redraw(iv);
	return 0;
}


void crystfel_image_view_reset_zoom(CrystFELImageView *iv)
{
	iv->detector_w = 1.0;
	iv->detector_h = 1.0;
	iv->zoom = -1.0;
	iv->need_recentre = 1;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_brightness(CrystFELImageView *iv,
                                        double brightness)
{
	iv->brightness = brightness;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_show_centre(CrystFELImageView *iv,
                                         int show_centre)
{
	iv->show_centre = show_centre;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_show_peaks(CrystFELImageView *iv,
                                        int show_peaks)
{
	iv->show_peaks = show_peaks;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_show_reflections(CrystFELImageView *iv,
                                              int show_refls)
{
	iv->show_refls = show_refls;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_label_reflections(CrystFELImageView *iv,
                                               int label_refls)
{
	iv->label_refls = label_refls;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_peak_box_size(CrystFELImageView *iv,
                                           float box_size)
{
	iv->peak_box_size = box_size;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_refl_box_size(CrystFELImageView *iv,
                                           float box_size)
{
	iv->refl_box_size = box_size;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_resolution_rings(CrystFELImageView *iv,
                                              int rings)
{
	iv->resolution_rings = rings;
	iv->need_rerender = 1;
	redraw(iv);
}


void crystfel_image_view_set_colour_scale(CrystFELImageView *iv,
                                          double lo, double hi)
{
	iv->scale_lo = lo;
	iv->scale_hi = hi;
	iv->need_rerender = 1;
	redraw(iv);
}
