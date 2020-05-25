/*
 * datatemplate.h
 *
 * Template for loading data
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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
extern void data_template_free(DataTemplate *dt);

extern signed int data_template_find_panel(const DataTemplate *dt,
                                           int fs, int ss);

extern void data_template_file_to_panel_coords(const DataTemplate *dt,
                                               float *pfs, float *pss);

extern void data_template_add_copy_header(DataTemplate *dt,
                                          const char *header);

#ifdef __cplusplus
}
#endif

#endif	/* DATATEMPLATE_H */
