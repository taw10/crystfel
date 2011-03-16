/*
 * reflist.h
 *
 * Fast reflection/peak list
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef REFLIST_H
#define REFLIST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


typedef struct _reflist RefList;
typedef struct _reflection Reflection;
typedef struct _reflistiterator RefListIterator;

#define INDICES signed int h, signed int k, signed int l

/* Creation/deletion */
extern RefList *reflist_new(void);
extern void reflist_free(RefList *list);

/* Search */
extern Reflection *find_refl(RefList *list, INDICES);
extern Reflection *next_found_refl(Reflection *refl);

/* Get */
extern double get_excitation_error(Reflection *refl);
extern void get_detector_pos(Reflection *refl, double *x, double *y);
extern void get_indices(Reflection *refl,
                        signed int *h, signed int *k, signed int *l);
extern double get_partiality(Reflection *refl);
extern double get_intensity(Reflection *refl);
extern void get_partial(Reflection *refl, double *r1, double *r2, double *p,
                        int *clamp_low, int *clamp_high);
extern int get_scalable(Reflection *refl);
extern int get_redundancy(Reflection *refl);
extern double get_sum_squared_dev(Reflection *refl);
extern double get_esd_intensity(Reflection *refl);
extern double get_phase(Reflection *refl);

/* Set */
extern void copy_data(Reflection *to, Reflection *from);
extern void set_detector_pos(Reflection *refl, double exerr,
                             double x, double y);
extern void set_partial(Reflection *refl, double r1, double r2, double p,
                        double clamp_low, double clamp_high);
extern void set_int(Reflection *refl, double intensity);
extern void set_scalable(Reflection *refl, int scalable);
extern void set_redundancy(Reflection *refl, int red);
extern void set_sum_squared_dev(Reflection *refl, double dev);
extern void set_esd_intensity(Reflection *refl, double esd);
extern void set_ph(Reflection *refl, double phase);

/* Insertion */
extern Reflection *add_refl(RefList *list, INDICES);

/* Deletion */
extern void delete_refl(Reflection *refl);

/* Iteration */
extern Reflection *first_refl(RefList *list, RefListIterator **iterator);
extern Reflection *next_refl(Reflection *refl, RefListIterator *iter);

/* Misc */
extern void optimise_reflist(RefList *list);
extern int num_reflections(RefList *list);

#endif	/* REFLIST_H */
