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
#include <gsl/gsl_statistics_float.h>

#include <utils.h>
#include <detgeom.h>
#include <colscale.h>

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
		                         0.0, iv->detector_w,
		                         0.0001, 0.1, vis_size);
	}
	if ( iv->vadj != NULL ) {
		double pos = gtk_adjustment_get_value(iv->vadj);
		double vis_size = iv->visible_height / iv->zoom;
		gtk_adjustment_configure(iv->vadj, pos,
		                         0.0, iv->detector_h,
		                         0.0001, 0.1, vis_size);
	}
}


static gint scroll_sig(GtkWidget *window, GdkEventScroll *event,
                       CrystFELImageView *iv)
{
	double zoom_scale;
	int claim = FALSE;

	if ( event->direction == GDK_SCROLL_UP ) {
		zoom_scale = 1.1;
		claim = TRUE;
	}
	if ( event->direction == GDK_SCROLL_DOWN ) {
		zoom_scale = 0.9;
		claim = TRUE;
	}
	if ( event->direction == GDK_SCROLL_LEFT ) return TRUE;
	if ( event->direction == GDK_SCROLL_RIGHT ) return TRUE;

	if ( claim ) {

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

	return claim;
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


static int clamp(double val, int min, int max)
{
	if ( val < min ) return min;
	if ( val > max ) return max;
	return val;
}


static void swap(int *a, int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}


static void draw_pixel_values(cairo_t *cr,
                              double imin_fs, double imin_ss,
                              double imax_fs, double imax_ss,
                              struct detgeom_panel p, float *dp,
                              int *bad)
{
	int min_fs, max_fs, min_ss, max_ss;
	int fs, ss;
	PangoLayout *layout;
	PangoFontDescription *fontdesc;
	double w, h;

	/* FIXME: This is wrong for slanty pixels */
	min_fs = clamp(imin_fs, 0, p.w-1);
	min_ss = clamp(imin_ss, 0, p.h-1);
	max_fs = clamp(imax_fs, 0, p.w-1);
	max_ss = clamp(imax_ss, 0, p.h-1);
	if ( min_ss > max_ss ) swap(&min_ss, &max_ss);
	if ( min_fs > max_fs ) swap(&min_fs, &max_fs);

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


static void draw_panel_rectangle(cairo_t *cr, CrystFELImageView *iv,
                                 int i, cairo_matrix_t *gtkmatrix)
{
	struct detgeom_panel p = iv->image->detgeom->panels[i];
	cairo_matrix_t m;
	cairo_pattern_t *patt;
	double min_x, min_y, max_x, max_y;
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
	min_x = 0.0;
	min_y = 0.0;
	max_x = iv->visible_width;
	max_y = iv->visible_height;

	/* Take into account the transformation which GTK already
	 * has in effect on the widget (translation from parent
	 * window's origin). */
	cairo_matrix_transform_point(gtkmatrix, &min_x, &min_y);
	cairo_matrix_transform_point(gtkmatrix, &max_x, &max_y);
	cairo_device_to_user(cr, &min_x, &min_y);
	cairo_device_to_user(cr, &max_x, &max_y);

	min_x /= p.pixel_pitch;
	min_y /= p.pixel_pitch;
	max_x /= p.pixel_pitch;
	max_y /= p.pixel_pitch;
	if ( (min_x < 0.0) && (max_x < 0.0) ) have_pixels = 0;
	if ( (min_y < 0.0) && (max_y < 0.0) ) have_pixels = 0;
	if ( (min_x > p.w) && (max_x > p.w) ) have_pixels = 0;
	if ( (min_y > p.h) && (max_y > p.h) ) have_pixels = 0;

	cairo_restore(cr);

	cairo_arc(cr, 0.0, 0.0, 0.006, 0, 2.0*M_PI);
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	cairo_set_line_width(cr, 0.00001);
	cairo_stroke(cr);

	xs = p.pixel_pitch;
	ys = p.pixel_pitch;
	cairo_user_to_device_distance(cr, &xs, &ys);
	pixel_size_on_screen = smallest(fabs(xs), fabs(ys));
	if ( (pixel_size_on_screen > 40.0) && have_pixels ) {
		draw_pixel_values(cr, min_x, min_y, max_x, max_y, p,
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
		cairo_set_source_rgb(cr, 1.0, 1.0, 0.0);
		cairo_stroke(cr);

		if ( show_cen ) {
			cairo_move_to(cr, x-0.2*p->pixel_pitch, y);
			cairo_line_to(cr, x+0.2*p->pixel_pitch, y);
			cairo_move_to(cr, x, y-0.2*p->pixel_pitch);
			cairo_line_to(cr, x, y+0.2*p->pixel_pitch);
			cairo_set_source_rgb(cr, 0.4, 0.4, 0.0);
			cairo_stroke(cr);
		}

	}
}


static void draw_refls(cairo_t *cr, CrystFELImageView *iv,
                       RefList *list)
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
		this_bs = biggest(iv->refl_box_size * p->pixel_pitch,
		                  bs);

		if ( this_bs > bs ) {
			show_cen = 1;
		}

		x = p->pixel_pitch*(p->cnx + p->fsx*fs + p->ssx*ss);
		y = p->pixel_pitch*(p->cny + p->fsy*fs + p->ssy*ss);

		cairo_arc(cr, x, y, this_bs, 0, 2*M_PI);
		cairo_set_line_width(cr, this_lw);
		cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
		cairo_stroke(cr);

		if ( show_cen ) {
			cairo_move_to(cr, x-0.2*p->pixel_pitch, y);
			cairo_line_to(cr, x+0.2*p->pixel_pitch, y);
			cairo_move_to(cr, x, y-0.2*p->pixel_pitch);
			cairo_line_to(cr, x, y+0.2*p->pixel_pitch);
			cairo_set_source_rgb(cr, 0.0, 0.4, 0.0);
			cairo_stroke(cr);
		}

	}
}


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELImageView *iv)
{
	cairo_matrix_t m;

	if ( iv->image == NULL ) return FALSE;

	cairo_save(cr);

	/* Overall background (light grey) */
	cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
	cairo_paint(cr);

	/* Get the transformation matrix before my transformations */
	cairo_get_matrix(cr, &m);

	cairo_scale(cr, iv->zoom, iv->zoom);
	cairo_translate(cr, iv->offs_x, iv->offs_y);
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

	if ( iv->show_peaks ) {
		draw_peaks(cr, iv, iv->image->features);
	}

	if ( iv->show_refls ) {
		int i;
		for ( i=0; i<iv->image->n_crystals; i++ ) {
			Crystal *cry = iv->image->crystals[i];
			draw_refls(cr, iv, crystal_get_reflections(cry));
		}
	}

	cairo_restore(cr);

	return FALSE;
}


static void scroll_adjust_sig(GtkAdjustment *adj, CrystFELImageView *iv)
{
	redraw(iv);
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
		iv->vadj = g_value_get_object(val);
		configure_scroll_adjustments(iv);
		if ( iv->vadj != NULL ) {
			g_signal_connect(G_OBJECT(iv->vadj), "value-changed",
			                 G_CALLBACK(scroll_adjust_sig), iv);
		}
		break;

		case CRYSTFELIMAGEVIEW_HADJ :
		iv->hadj = g_value_get_object(val);
		configure_scroll_adjustments(iv);
		if ( iv->hadj != NULL ) {
			g_signal_connect(G_OBJECT(iv->hadj), "value-changed",
			                 G_CALLBACK(scroll_adjust_sig), iv);
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

	/* All values initially meaningless */
	iv->detector_w = 1.0;
	iv->detector_h = 1.0;
	iv->zoom = -1.0;
	iv->image = NULL;
	iv->show_peaks = 0;
	iv->brightness = 1.0;
	iv->pixbufs = NULL;
	iv->peak_box_size = 1.0;
	iv->refl_box_size = 1.0;

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
	                       | GDK_SCROLL_MASK);

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
                               int scale_type, double scale_top)


{
	guchar *pixbuf_data;
	long int i;

	/* Rendered (colourful) version */
	pixbuf_data = malloc(3*w*h);
	if ( pixbuf_data == NULL ) return NULL;

	for ( i=0; i<w*h; i++ ) {

		double r, g, b;

		if ( !badmap[i] ) {

			colscale_lookup(data[i], scale_top,
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


static int rerender_image(CrystFELImageView *iv)
{
	int i;
	double min_x, min_y, max_x, max_y;
	double border;
	double scale_top;

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

	scale_top = auto_scale_top(iv->image);

	for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
		iv->pixbufs[i] = render_panel(iv->image->dp[i],
		                              iv->image->bad[i],
		                              iv->image->detgeom->panels[i].w,
		                              iv->image->detgeom->panels[i].h,
		                              SCALE_COLOUR, scale_top);
		if ( iv->pixbufs[i] == NULL ) return 1;
	}

	detgeom_pixel_extents(iv->image->detgeom, &min_x, &min_y,
	                      &max_x, &max_y);
	iv->detector_w = max_x - min_x;
	iv->detector_h = max_y - min_y;
	border = iv->detector_w * 0.1;
	iv->detector_w += border;
	iv->detector_h += border;
	if ( iv->zoom < 0.0 ) {
		/* Set initial values */
		iv->offs_x = -min_x + border/2.0;
		iv->offs_y = max_y + border/2.0;
		iv->zoom = 1.0/iv->image->detgeom->panels[0].pixel_pitch;
	}
	configure_scroll_adjustments(iv);

	redraw(iv);

	return 0;
}


int crystfel_image_view_set_image(CrystFELImageView *iv,
                                  const struct image *image)
{
	cleanup_image(iv);
	iv->image = image;
	return rerender_image(iv);
}


void crystfel_image_view_reset_zoom(CrystFELImageView *iv)
{
	iv->detector_w = 1.0;
	iv->detector_h = 1.0;
	iv->zoom = -1.0;
}


void crystfel_image_view_set_brightness(CrystFELImageView *iv,
                                        double brightness)
{
	iv->brightness = brightness;
	rerender_image(iv);
}


void crystfel_image_view_set_show_peaks(CrystFELImageView *iv,
                                        int show_peaks)
{
	iv->show_peaks = show_peaks;
	rerender_image(iv);
}


void crystfel_image_view_set_show_reflections(CrystFELImageView *iv,
                                              int show_refls)
{
	iv->show_refls = show_refls;
	rerender_image(iv);
}


void crystfel_image_view_set_peak_box_size(CrystFELImageView *iv,
                                           float box_size)
{
	iv->peak_box_size = box_size;
}


void crystfel_image_view_set_refl_box_size(CrystFELImageView *iv,
                                           float box_size)
{
	iv->refl_box_size = box_size;
}
