/*
 * datatemplate.h
 *
 * Template for loading data
 *
 * Copyright © 2019-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2020 Thomas White <taw@physics.org>
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

#ifndef DATATEMPLATE_H
#define DATATEMPLATE_H

#include "detgeom.h"

/**
 * \file datatemplate.h
 * Template for loading data.
 */

/**
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 **/
typedef struct _datatemplate DataTemplate;


#ifdef __cplusplus
extern "C" {
#endif

extern DataTemplate *data_template_new_from_file(const char *filename);
extern DataTemplate *data_template_new_from_string(const char *string_in);
extern void data_template_free(DataTemplate *dt);

extern const char *data_template_panel_number_to_name(const DataTemplate *dt,
                                                      int pn);

extern int data_template_panel_name_to_number(const DataTemplate *dt,
                                              const char *panel_name,
                                              int *pn);

extern int data_template_file_to_panel_coords(const DataTemplate *dt,
                                              float *pfs, float *pss,
                                              int pn);

extern int data_template_slabby_file_to_panel_coords(const DataTemplate *dt,
                                                     float *pfs, float *pss,
                                                     int *ppn);

extern int data_template_panel_to_file_coords(const DataTemplate *dt,
                                              int pn,
                                              float *pfs, float *pss);

extern void data_template_add_copy_header(DataTemplate *dt,
                                          const char *header);

extern int data_template_get_slab_extents(const DataTemplate *dt, int *pw, int *ph);

extern double data_template_get_wavelength_if_possible(const DataTemplate *dt);

extern double data_template_get_clen_if_possible(const DataTemplate *dt);

extern struct detgeom *data_template_get_2d_detgeom_if_possible(const DataTemplate *dt);

extern void data_template_show_hierarchy(const DataTemplate *dtempl);

#ifdef __cplusplus
}
#endif

#endif	/* DATATEMPLATE_H */
