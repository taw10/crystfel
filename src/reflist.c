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

	unsigned int serial;  /* Serial number */

	signed int h;
	signed int k;
	signed int l;

	double excitation_error;

	/* Partiality */
	double r1;  /* First excitation error */
	double r2;  /* Second excitation error */
	double p;   /* Partiality */
	int clamp1; /* Clamp status for r1 */
	int clamp2; /* Clamp status for r2 */

	/* Location in image */
	int x;
	int y;

	int scalable;

	/* Intensity */
	double intensity;


};


struct _reflist {

	struct ref *head;

};


/**************************** Creation / deletion *****************************/

/* Create a reflection list */
RefList *reflist_new()
{
	RefList *new;

	new = malloc(sizeof(struct _reflist));

	return new;
}


void reflist_free(RefList *list)
{
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
}


void get_detector_pos(Reflection *refl, double *x, double *y)
{
}


void get_indices(Reflection *refl, signed int *h, signed int *k, signed int *l)
{
}


double get_partiality(Reflection *refl)
{
}


double get_intensity(Reflection *refl)
{
}


void get_partial(Reflection *refl, double *r1, double *r2, double *p,
                 int *clamp_low, int *clamp_high)
{
}


int get_scalable(Reflection *refl)
{
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
