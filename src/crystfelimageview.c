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

#include <utils.h>
#include <detgeom.h>
#include <render.h>
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
	int i;

	if ( iv->image != NULL ) {
		if ( (iv->image->detgeom != NULL) && (iv->pixbufs != NULL) ) {
			for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
				if ( iv->pixbufs[i] != NULL ) {
					gdk_pixbuf_unref(iv->pixbufs[i]);
				}
			}
		}
		image_free(iv->image);
	}
}


static gint destroy_sig(GtkWidget *window, CrystFELImageView *iv)
{
	cleanup_image(iv);
	free(iv->filename);
	free(iv->event);
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


static void draw_panel_rectangle(cairo_t *cr, CrystFELImageView *iv, int i)
{
	struct detgeom_panel p = iv->image->detgeom->panels[i];
	cairo_matrix_t m;
	cairo_pattern_t *patt;

	/* Move to the right location */
	cairo_translate(cr, p.cnx*p.pixel_pitch, p.cny*p.pixel_pitch);

	/* Twiddle directions according to matrix */
	cairo_matrix_init(&m, rint(p.fsx), rint(p.fsy), rint(p.ssx), rint(p.ssy),
	                      0.0, 0.0);
	cairo_transform(cr, &m);

	gdk_cairo_set_source_pixbuf(cr, iv->pixbufs[i], 0.0, 0.0);
	patt = cairo_get_source(cr);

	cairo_pattern_get_matrix(patt, &m);
	cairo_matrix_scale(&m, 1.0/p.pixel_pitch, 1.0/p.pixel_pitch);
	cairo_pattern_set_matrix(patt, &m);

	cairo_pattern_set_filter(patt, CAIRO_FILTER_NEAREST);

	cairo_rectangle(cr, 0.0, 0.0, p.w*p.pixel_pitch, p.h*p.pixel_pitch);
	cairo_fill(cr);
}


static void draw_peaks(cairo_t *cr, CrystFELImageView *iv,
                       ImageFeatureList *pks)
{
	int i, n_pks;
	double bs, lw;

	bs = 5.0;
	lw = 1.0;
	cairo_device_to_user_distance(cr, &bs, &lw);
	bs = fabs(bs);
	lw = fabs(lw);

	n_pks = image_feature_count(pks);
	for ( i=0; i<n_pks; i++ ) {
		const struct imagefeature *f;
		struct detgeom_panel *p;
		double x, y;
		f = image_get_feature_const(pks, i);
		if ( f == NULL ) continue;
		p = &iv->image->detgeom->panels[f->pn];
		x = p->pixel_pitch*(p->cnx + p->fsx*f->fs + p->ssx*f->ss);
		y = p->pixel_pitch*(p->cny + p->fsy*f->fs + p->ssy*f->ss);
		cairo_rectangle(cr, x-bs, y-bs, 2*bs, 2*bs);
		cairo_set_line_width(cr, lw);
		cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
		cairo_stroke(cr);
	}
}


static gint draw_sig(GtkWidget *window, cairo_t *cr, CrystFELImageView *iv)
{
	int i;

	cairo_save(cr);

	/* Overall background (light grey) */
	cairo_set_source_rgb(cr, 0.7, 0.7, 0.7);
	cairo_paint(cr);

	cairo_scale(cr, iv->zoom, iv->zoom);
	cairo_translate(cr, iv->offs_x, iv->offs_y);
	cairo_scale(cr, 1.0, -1.0);
	cairo_translate(cr, -gtk_adjustment_get_value(iv->hadj),
	                     gtk_adjustment_get_value(iv->vadj));

	if ( (iv->image != NULL) && (iv->pixbufs != NULL) ) {
		int i;
		for ( i=0; i<iv->image->detgeom->n_panels; i++ ) {
			cairo_save(cr);
			draw_panel_rectangle(cr, iv, i);
			cairo_restore(cr);
		}
	}

	for ( i=0; i<iv->num_peaklists; i++ ){
		if ( iv->peaklists[i] != NULL ) {
			draw_peaks(cr, iv , iv->peaklists[i]);
		}
	}

	cairo_restore(cr);
	return FALSE;
}


static void scroll_adjust_sig(GtkAdjustment *adj, CrystFELImageView *iv)
{
	redraw(iv);
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
	iv->filename = NULL;
	iv->event = NULL;
	iv->image = NULL;
	iv->num_peaklists = 0;
	iv->peaklists = NULL;

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


static int reload_image(CrystFELImageView *iv)
{
	int n_pb;
	double min_x, min_y, max_x, max_y;
	double border;

	if ( iv->dtempl == NULL ) return 0;
	if ( iv->filename == NULL ) return 0;

	cleanup_image(iv);

	iv->image = image_read(iv->dtempl, iv->filename, iv->event);
	if ( iv->image == NULL ) {
		ERROR("Failed to load image\n");
		return 1;
	}

	iv->pixbufs = render_panels(iv->image, 1, SCALE_COLOUR, 5, &n_pb);
	if ( n_pb != iv->image->detgeom->n_panels ) {
		ERROR("Wrong number of panels returned!\n");
		return 1;
	}

	detgeom_pixel_extents(iv->image->detgeom, &min_x, &min_y, &max_x, &max_y);
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


int crystfel_image_view_set_datatemplate(CrystFELImageView *iv,
                                         DataTemplate *dtempl)
{
	iv->dtempl = dtempl;
	return reload_image(iv);
}


int crystfel_image_view_set_image(CrystFELImageView *iv,
                                  const char *filename,
                                  const char *event)
{
	free(iv->filename);
	free(iv->event);
	iv->filename = safe_strdup(filename);
	iv->event = safe_strdup(event);
	return reload_image(iv);
}


struct image *crystfel_image_view_get_image_struct(CrystFELImageView *iv)
{
	return iv->image;
}


void crystfel_image_view_set_peaks(CrystFELImageView *iv,
                                   ImageFeatureList *peaks,
                                   int list_num)
{
	int i;
	if ( list_num >= iv->num_peaklists ) {
		ImageFeatureList **n_fl;
		n_fl = realloc(iv->peaklists,
		               (list_num+1)*sizeof(ImageFeatureList*));
		if ( n_fl == NULL ) return;
		for ( i=iv->num_peaklists; i<list_num+1; i++ ) {
			n_fl[i] = NULL;
		}
		iv->peaklists = n_fl;
		iv->num_peaklists = list_num+1;
	}
	if ( iv->peaklists[list_num] != NULL ) {
		image_feature_list_free(iv->peaklists[list_num]);
	}
	iv->peaklists[list_num] = image_feature_list_copy(peaks);

	redraw(iv);
}
