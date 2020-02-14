/*
 * crystfelimageview.h
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

#ifndef CRYSTFELIMAGEVIEW_H
#define CRYSTFELIMAGEVIEW_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

	int                  w;   /* Surface size in pixels */
	int                  h;

	/* Redraw/scroll stuff */
	GtkScrollablePolicy  hpol;
	GtkScrollablePolicy  vpol;
	GtkAdjustment       *hadj;
	GtkAdjustment       *vadj;
	double               x_scroll_pos;
	double               y_scroll_pos;
};

struct _crystfelimageviewclass
{
	GtkDrawingAreaClass parent_class;
};

typedef struct _crystfelimageview CrystFELImageView;
typedef struct _crystfelimageviewclass CrystFELImageViewClass;

extern GType crystfel_image_view_get_type(void);
extern GtkWidget *crystfel_image_view_new(void);


#endif	/* CRYSTFELIMAGEVIEW_H */
