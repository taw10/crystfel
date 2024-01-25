/*
 * peakfinder8.h
 *
 * The peakfinder8 algorithm
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
 *   2017 Valerio Mariani <valerio.mariani@desy.de>
 *   2017 Anton Barty <anton.barty@desy.de>
 *   2017 Oleksandr Yefanov <oleksandr.yefanov@desy.de>
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

#ifndef PEAKFINDER8_H
#define PEAKFINDER8_H

#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file peakfinder8.h
 * The "peakfinder8" peak search algorithm, originally implemented in Cheetah
 */

struct pf8_private_data
{
    int fast_mode;
    struct radius_maps *rmaps;
    struct radial_stats_pixels *rpixels;
};

struct pf8_private_data *prepare_peakfinder8(struct detgeom *det, int fast_mode);

void free_pf8_private_data(struct pf8_private_data *data);

extern ImageFeatureList *peakfinder8(const struct image *img, int max_n_peaks,
                                     float threshold, float min_snr,
                                     int mix_pix_count, int max_pix_count,
                                     int local_bg_radius, int min_res,
                                     int max_res, int use_saturated,
                                     int fast_mode,
                                     struct pf8_private_data *private_data);

#ifdef __cplusplus
}
#endif

#endif // PEAKFINDER8_H
