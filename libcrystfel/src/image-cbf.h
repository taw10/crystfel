/*
 * image-cbf.h
 *
 * Image loading, CBF parts
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

/* NB This file is NOT part of the public API, and should NOT
 * be installed, but rather stays in the libcrystfel source folder. */

#ifndef IMAGE_CBF_H
#define IMAGE_CBF_H

#include "datatemplate_priv.h"

extern signed int is_cbf_file(const char *filename);

extern signed int is_cbfgz_file(const char *filename);

extern int load_mask_cbf(struct panel_template *p,
                         const char *filename, const char *event,
                         int gz, int *bad, int mask_good, int mask_bad);

extern int image_cbf_read(struct image *image,
                          const DataTemplate *dtempl,
                          const char *filename,
                          const char *event,
                          int gz);

extern int image_cbf_read_mask(struct panel_template *p,
                               const char *filename, const char *event,
                               int gz, int *bad,
                               int mask_good, int mask_bad);

#endif	/* IMAGE_CBF_H */
