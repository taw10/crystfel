/*
 * reflist.c
 *
 * Fast reflection/peak list
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>

#include "reflist.h"


struct _reflection {

	/* Listy stuff */
	unsigned int serial;          /* Unique serial number, key */
	struct _reflection *child[2]; /* Child nodes */
	struct _reflection *parent;   /* Parent node */

	signed int h;
	signed int k;
	signed int l;

	/* Partiality and related geometrical stuff */
	double r1;  /* First excitation error */
	double r2;  /* Second excitation error */
	double p;   /* Partiality */
	int clamp1; /* Clamp status for r1 */
	int clamp2; /* Clamp status for r2 */

	/* Location in image */
	double x;
	double y;

	/* The distance from the exact Bragg position to the coordinates
	 * given above. */
	double excitation_error;

	/* Non-zero if this reflection can be used for scaling */
	int scalable;

	/* Intensity */
	double intensity;
};


struct _reflist {

	struct _reflection *head;

};


/**************************** Creation / deletion *****************************/

/* Create a reflection list */
RefList *reflist_new()
{
	RefList *new;

	new = malloc(sizeof(struct _reflist));
	new->head = NULL;

	return new;
}


static void recursive_free(Reflection *refl)
{
	if ( refl->child[0] != NULL ) recursive_free(refl->child[0]);
	if ( refl->child[1] != NULL ) recursive_free(refl->child[1]);
	free(refl);
}


void reflist_free(RefList *list)
{
	recursive_free(list->head);
	free(list);
}


/********************************** Search ************************************/

Reflection *find_refl(RefList *list, INDICES)
{
}


Reflection *next_found_refl(Reflection *refl)
{
}


/********************************** Getters ***********************************/

double get_excitation_error(Reflection *refl)
{
	return refl->excitation_error;
}


void get_detector_pos(Reflection *refl, double *x, double *y)
{
	*x = refl->x;
	*y = refl->y;
}


void get_indices(Reflection *refl, signed int *h, signed int *k, signed int *l)
{
	*h = refl->h;
	*k = refl->k;
	*l = refl->l;
}


double get_partiality(Reflection *refl)
{
	return refl->p;
}


double get_intensity(Reflection *refl)
{
	return refl->intensity;
}


void get_partial(Reflection *refl, double *r1, double *r2, double *p,
                 int *clamp_low, int *clamp_high)
{
	*r1 = refl->r1;
	*r2 = refl->r2;
	*p = get_partiality(refl);
	*clamp_low = refl->clamp1;
	*clamp_high = refl->clamp2;
}


int get_scalable(Reflection *refl)
{
	return refl->scalable;
}


/********************************** Setters ***********************************/

void set_detector_pos(Reflection *refl, double exerr, double x, double y)
{
}


void set_partial(Reflection *refl, double r1, double r2, double p,
                 double clamp_low, double clamp_high)
{
}


void set_indices(Reflection *refl,
                        signed int h, signed int k, signed int l)
{
}


void set_int(Reflection *refl, double intensity)
{
}


void set_scalable(Reflection *refl, int scalable)
{
}


/********************************* Insertion **********************************/

Reflection *add_refl(RefList *list, INDICES)
{
}


Reflection *add_refl_with_det_pos(RefList *refl, INDICES, double exerr,
                                  double x, double y)
{
}


/********************************** Deletion **********************************/

void delete_refl(Reflection *refl)
{
}


/********************************* Iteration **********************************/

Reflection *first_refl(RefList *list)
{
}


Reflection *next_refl(Reflection *refl)
{
}


/*********************************** Voodoo ***********************************/

void optimise_reflist(RefList *list)
{
}
