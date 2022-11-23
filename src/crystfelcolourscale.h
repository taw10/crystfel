/*
 * crystfelcolourscale.h
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

#ifndef CRYSTFELCOLOURSCALE_H
#define CRYSTFELCOLOURSCALE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "image.h"

#define CRYSTFEL_TYPE_COLOUR_SCALE (crystfel_colour_scale_get_type())

#define CRYSTFEL_COLOUR_SCALE(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj), \
                                 CRYSTFEL_TYPE_COLOUR_SCALE, CrystFELColourScale))

#define CRYSTFEL_IS_COLOUR_SCALE(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj), \
                                    CRYSTFEL_TYPE_COLOUR_SCALE))

#define CRYSTFEL_COLOUR_SCALE_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((obj), \
                                         CRYSTFEL_TYPE_COLOUR_SCALE, CrystFELColourScale))

#define CRYSTFEL_IS_COLOUR_SCALE_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((obj), \
                                            CRYSTFEL_TYPE_COLOUR_SCALE))

#define CRYSTFEL_COLOUR_SCALE_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS((obj), \
                                           CRYSTFEL_TYPE_COLOUR_SCALE, CrystFELColourScale))


#define COLSCALE_N_BINS (256)
#define COLSCALE_SAMPLE_SIZE (4096)

struct _crystfelcolourscale
{
	GtkDrawingArea       parent_instance;
	double               visible_width;
	double               visible_height;
	double               drag_start_x;
	double               drag_start_y;
	double               drag_min;

	double               lo;
	double               hi;
	int                  bins[COLSCALE_N_BINS];
	float               *sample;
	int                  n_samples;
};

struct _crystfelcolourscaleclass
{
	GtkDrawingAreaClass parent_class;
};

typedef struct _crystfelcolourscale CrystFELColourScale;
typedef struct _crystfelcolourscaleclass CrystFELColourScaleClass;

extern GType crystfel_colour_scale_get_type(void);
extern GtkWidget *crystfel_colour_scale_new(void);

extern void crystfel_colour_scale_scan_image(CrystFELColourScale *cs,
                                             struct image *image);

extern void crystfel_colour_scale_get_range(CrystFELColourScale *cs,
                                            double scale_min,
                                            double scale_max);

extern void crystfel_colour_scale_auto_range(CrystFELColourScale *cs);

#endif	/* CRYSTFELCOLOURSCALE_H */
