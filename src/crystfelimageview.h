/*
 * crystfelimageview.h
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

#ifndef CRYSTFELIMAGEVIEW_H
#define CRYSTFELIMAGEVIEW_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gdk-pixbuf/gdk-pixbuf.h>

#include <image.h>
#include <datatemplate.h>

#define CRYSTFEL_TYPE_IMAGE_VIEW (crystfel_image_view_get_type())

#define CRYSTFEL_IMAGE_VIEW(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_IMAGE_VIEW, CrystFELImageView))

#define CRYSTFEL_IS_IMAGE_VIEW(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_IMAGE_VIEW))

#define CRYSTFEL_IMAGE_VIEW_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_IMAGE_VIEW, CrystFELImageView))

#define CRYSTFEL_IS_IMAGE_VIEW_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_IMAGE_VIEW))

#define CRYSTFEL_IMAGE_VIEW_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_IMAGE_VIEW, CrystFELImageView))

struct _crystfelimageview
{
	GtkDrawingArea       parent_instance;

	/*< private >*/
	GtkIMContext        *im_context;

	/* Detector size in metres */
	double               detector_w;
	double               detector_h;

	/* Redraw/scroll stuff */
	int                  need_rerender;
	int                  need_recentre;
	GtkScrollablePolicy  hpol;
	GtkScrollablePolicy  vpol;
	GtkAdjustment       *hadj;
	GtkAdjustment       *vadj;
	double               visible_width;
	double               visible_height;
	double               zoom;  /* screen pixels per m */
	double               drag_start_x;
	double               drag_start_y;
	double               drag_start_sp_x;
	double               drag_start_sp_y;
	double               min_x;
	double               max_y;

	const struct image  *image;

	GdkPixbuf          **pixbufs;

	double               brightness;
	int                  show_centre;
	int                  show_peaks;
	int                  show_refls;
	int                  label_refls;
	float                peak_box_size;
	float                refl_box_size;
	int                  resolution_rings;
};

struct _crystfelimageviewclass
{
	GtkDrawingAreaClass parent_class;
};

typedef struct _crystfelimageview CrystFELImageView;
typedef struct _crystfelimageviewclass CrystFELImageViewClass;

extern GType crystfel_image_view_get_type(void);
extern GtkWidget *crystfel_image_view_new(void);

extern int crystfel_image_view_set_image(CrystFELImageView *iv,
                                         const struct image *image);

extern void crystfel_image_view_reset_zoom(CrystFELImageView *iv);

extern void crystfel_image_view_set_brightness(CrystFELImageView *iv,
                                               double brightness);

extern void crystfel_image_view_set_show_centre(CrystFELImageView *iv,
                                                int show_centre);

extern void crystfel_image_view_set_show_peaks(CrystFELImageView *iv,
                                               int show_peaks);

extern void crystfel_image_view_set_show_reflections(CrystFELImageView *iv,
                                                     int show_refls);

extern void crystfel_image_view_set_label_reflections(CrystFELImageView *iv,
                                                      int label_refls);

extern void crystfel_image_view_set_peak_box_size(CrystFELImageView *iv,
                                                  float box_size);

extern void crystfel_image_view_set_refl_box_size(CrystFELImageView *iv,
                                                  float box_size);

extern void crystfel_image_view_set_resolution_rings(CrystFELImageView *iv,
                                                     int rings);

#endif	/* CRYSTFELIMAGEVIEW_H */
