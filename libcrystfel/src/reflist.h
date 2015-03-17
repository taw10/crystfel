/*
 * reflist.h
 *
 * Fast reflection/peak list
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2015 Thomas White <taw@physics.org>
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

#ifndef REFLIST_H
#define REFLIST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define SERIAL(h, k, l) ((((h)+512)<<20) + (((k)+512)<<10) + ((l)+512))
#define GET_H(serial) ((((serial) & 0x3ff00000)>>20)-512)
#define GET_K(serial) ((((serial) & 0x000ffc00)>>10)-512)
#define GET_L(serial) (((serial) & 0x000003ff)-512)

/**
 * RefList:
 *
 * A %RefList represents a list of Bragg reflections.
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 *
 **/
typedef struct _reflist RefList;

/**
 * Reflection:
 *
 * A %Reflection represents a single Bragg reflection.
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 *
 **/
typedef struct _reflection Reflection;

/**
 * RefListIterator:
 *
 * A %RefListIterator is an opaque data type used when iterating over a
 * %RefList.
 *
 **/
typedef struct _reflistiterator RefListIterator;

#ifdef __cplusplus
extern "C" {
#endif

/* Creation/deletion */
extern RefList *reflist_new(void);
extern void reflist_free(RefList *list);
extern Reflection *reflection_new(signed int h, signed int k, signed int l);
extern void reflection_free(Reflection *refl);

/* Search */
extern Reflection *find_refl(const RefList *list, signed int h, signed int k, signed int l);
extern Reflection *next_found_refl(Reflection *refl);

/* Get */
extern void get_detector_pos(const Reflection *refl, double *fs, double *ss);
extern struct panel *get_panel(const Reflection *refl);
extern double get_partiality(const Reflection *refl);
extern double get_lorentz(const Reflection *refl);
extern void get_indices(const Reflection *refl,
                        signed int *h, signed int *k, signed int *l);
extern void get_symmetric_indices(const Reflection *refl,
                                  signed int *hs, signed int *ks,
                                  signed int *ls);
extern double get_intensity(const Reflection *refl);
extern void get_partial(const Reflection *refl, double *rlow, double *rhigh,
                        double *p);
extern int get_redundancy(const Reflection *refl);
extern double get_temp1(const Reflection *refl);
extern double get_temp2(const Reflection *refl);
extern double get_esd_intensity(const Reflection *refl);
extern double get_phase(const Reflection *refl, int *have_phase);
extern double get_peak(const Reflection *refl);
extern double get_mean_bg(const Reflection *refl);

/* Set */
extern void copy_data(Reflection *to, const Reflection *from);
extern void set_detector_pos(Reflection *refl, double fs, double ss);
extern void set_panel(Reflection *refl, struct panel *p);
extern void set_partial(Reflection *refl, double rlow, double rhigh, double p);
extern void set_partiality(Reflection *refl, double p);
extern void set_lorentz(Reflection *refl, double L);
extern void set_intensity(Reflection *refl, double intensity);
extern void set_redundancy(Reflection *refl, int red);
extern void set_temp1(Reflection *refl, double temp);
extern void set_temp2(Reflection *refl, double temp);
extern void set_esd_intensity(Reflection *refl, double esd);
extern void set_phase(Reflection *refl, double phase);
extern void set_peak(Reflection *refl, double peak);
extern void set_mean_bg(Reflection *refl, double mean_bg);
extern void set_symmetric_indices(Reflection *refl,
                                  signed int hs, signed int ks, signed int ls);

/* Insertion */
extern Reflection *add_refl(RefList *list,
                            signed int h, signed int k, signed int l);
extern void add_refl_to_list(Reflection *refl, RefList *list);

/* Iteration */
extern Reflection *first_refl(RefList *list, RefListIterator **piter);
extern Reflection *next_refl(Reflection *refl, RefListIterator *iter);

/* Misc */
extern int num_reflections(RefList *list);
extern int tree_depth(RefList *list);
extern void lock_reflection(Reflection *refl);
extern void unlock_reflection(Reflection *refl);

#ifdef __cplusplus
}
#endif

#endif	/* REFLIST_H */
