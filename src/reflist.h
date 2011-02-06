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

#define INDICES signed int h, signed int k, signed int l

/* Creation/deletion */
extern RefList *reflist_new(void);
extern void reflist_free(RefList *list);

/* Search */
extern Reflection *find_refl(RefList *list, INDICES);

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

/* Set */
extern void set_detector_pos(Reflection *refl, double exerr,
                             double x, double y);
extern void set_partial(Reflection *refl, double r1, double r2, double p,
                        double clamp_low, double clamp_high);
extern void set_indices(Reflection *refl,
                        signed int h, signed int k, signed int l);
extern void set_int(Reflection *refl, double intensity);
extern void set_scalable(Reflection *refl, int scalable);

/* Insertion */
extern Reflection *add_refl(RefList *list, INDICES);
extern Reflection *add_refl_with_det_pos(RefList *list, INDICES, double exerr,
                                        double x, double y);


/* Deletion */
extern void delete_refl(Reflection *refl);

/* Iteration */
extern Reflection *first_refl(RefList *list);
extern Reflection *next_refl(Reflection *refl);

/* Voodoo */
extern void optimise_reflist(RefList *list);

#endif	/* REFLIST_H */
