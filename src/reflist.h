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

/**
 * RefList:
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 */
typedef struct _reflist RefList;

/**
 * Reflection:
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 */
typedef struct _reflection Reflection;
typedef struct _reflistiterator RefListIterator;

/* Creation/deletion */
extern RefList *reflist_new(void);
extern void reflist_free(RefList *list);

/* Search */
extern Reflection *find_refl(const RefList *list, signed int h, signed int k, signed int l);
extern Reflection *next_found_refl(Reflection *refl);

/* Get */
extern double get_excitation_error(const Reflection *refl);
extern void get_detector_pos(const Reflection *refl, double *fs, double *ss);
extern double get_partiality(const Reflection *refl);
extern void get_indices(const Reflection *refl,
                        signed int *h, signed int *k, signed int *l);
extern void get_symmetric_indices(const Reflection *refl,
                                  signed int *hs, signed int *ks,
                                  signed int *ls);
extern double get_intensity(const Reflection *refl);
extern void get_partial(const Reflection *refl, double *r1, double *r2,
                        double *p, int *clamp_low, int *clamp_high);
extern int get_scalable(const Reflection *refl);
extern int get_redundancy(const Reflection *refl);
extern double get_sum_squared_dev(const Reflection *refl);
extern double get_esd_intensity(const Reflection *refl);
extern double get_phase(const Reflection *refl);

/* Set */
extern void copy_data(Reflection *to, const Reflection *from);
extern void set_detector_pos(Reflection *refl, double exerr,
                             double fs, double ss);
extern void set_partial(Reflection *refl, double r1, double r2, double p,
                        double clamp_low, double clamp_high);
extern void set_int(Reflection *refl, double intensity);
extern void set_scalable(Reflection *refl, int scalable);
extern void set_redundancy(Reflection *refl, int red);
extern void set_sum_squared_dev(Reflection *refl, double dev);
extern void set_esd_intensity(Reflection *refl, double esd);
extern void set_ph(Reflection *refl, double phase);
extern void set_symmetric_indices(Reflection *refl,
                                  signed int hs, signed int ks, signed int ls);

/* Insertion */
extern Reflection *add_refl(RefList *list,
                            signed int h, signed int k, signed int l);

/* Iteration */
extern Reflection *first_refl(RefList *list, RefListIterator **iterator);
extern Reflection *next_refl(Reflection *refl, RefListIterator *iter);

/* Misc */
extern int num_reflections(RefList *list);
extern int tree_depth(RefList *list);

#endif	/* REFLIST_H */
