/*
 * crystal.h
 *
 * A class representing a single crystal
 *
 * Copyright © 2013-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2013-2015 Thomas White <taw@physics.org>
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

#ifndef CRYSTAL_H
#define CRYSTAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "cell.h"
#include "reflist.h"


/**
 * Crystal:
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 **/
typedef struct _crystal Crystal;

#ifdef __cplusplus
extern "C" {
#endif

extern Crystal *crystal_new(void);
extern Crystal *crystal_copy(Crystal *cryst);
extern void crystal_free(Crystal *cryst);

extern UnitCell *crystal_get_cell(Crystal *cryst);
extern double crystal_get_profile_radius(Crystal *cryst);
extern RefList *crystal_get_reflections(Crystal *cryst);
extern double crystal_get_resolution_limit(Crystal *cryst);
extern long long int crystal_get_num_saturated_reflections(Crystal *cryst);
extern long long int crystal_get_num_implausible_reflections(Crystal *cryst);
extern int crystal_get_user_flag(Crystal *cryst);
extern double crystal_get_osf(Crystal *cryst);
extern double crystal_get_Bfac(Crystal *cryst);
extern struct image *crystal_get_image(Crystal *cryst);
extern double crystal_get_mosaicity(Crystal *cryst);
extern const char *crystal_get_notes(Crystal *cryst);

extern void crystal_set_cell(Crystal *cryst, UnitCell *cell);
extern void crystal_set_profile_radius(Crystal *cryst, double r);
extern void crystal_set_reflections(Crystal *cryst, RefList *reflist);
extern void crystal_set_resolution_limit(Crystal *cryst, double res);
extern void crystal_set_num_saturated_reflections(Crystal *cryst,
                                                  long long int n);
extern void crystal_set_num_implausible_reflections(Crystal *cryst,
                                                    long long int n);
extern void crystal_set_user_flag(Crystal *cryst, int flag);
extern void crystal_set_osf(Crystal *cryst, double osf);
extern void crystal_set_Bfac(Crystal *cryst, double B);
extern void crystal_set_image(Crystal *cryst, struct image *image);
extern void crystal_set_mosaicity(Crystal *cryst, double m);
extern void crystal_set_notes(Crystal *cryst, const char *notes);

extern void crystal_add_notes(Crystal *cryst, const char *notes_add);

#ifdef __cplusplus
}
#endif

#endif	/* CRYSTAL_H */
