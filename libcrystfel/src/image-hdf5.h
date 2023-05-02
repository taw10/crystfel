/*
 * image-hdf5.h
 *
 * Image loading, HDF5 parts
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

/* NB This file is NOT part of the public API, and should NOT
 * be installed, but rather stays in the libcrystfel source folder. */

#ifndef IMAGE_HDF5_H
#define IMAGE_HDF5_H

#include "datatemplate_priv.h"

extern int image_hdf5_read_header_to_cache(struct image *image,
                                           const char *name);

extern int image_hdf5_read(struct image *image,
                           const DataTemplate *dtempl);

extern int image_hdf5_read_mask(struct panel_template *p,
                                const char *filename,
                                const char *event, int *bad,
                                const char *mask_location,
                                int mask_good, int mask_bad);

extern float *image_hdf5_read_satmap(struct panel_template *p,
                                     const char *filename,
                                     const char *event,
                                     const char *map_location);

extern ImageFeatureList *image_hdf5_read_peaks_cxi(const DataTemplate *dtempl,
                                                   const char *filename,
                                                   const char *event,
                                                   int half_pixel_shift);

extern ImageFeatureList *image_hdf5_read_peaks_hdf5(const DataTemplate *dtempl,
                                                    const char *filename,
                                                    const char *event,
                                                    int half_pixel_shift);

extern char **image_hdf5_expand_frames(const DataTemplate *dtempl,
                                       const char *filename,
                                       int *n_frames);

extern int is_hdf5_file(const char *filename);

#endif	/* IMAGE_HDF5_H */
