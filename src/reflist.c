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
#include <assert.h>

#include "reflist.h"


struct _reflection {

	/* Listy stuff */
	RefList *list;
	unsigned int serial;          /* Unique serial number, key */
	struct _reflection *child[2]; /* Child nodes */
	struct _reflection *parent;   /* Parent node */
	struct _reflection *next;     /* Another reflection with the same
	                               * indices, or NULL */

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


#define SERIAL(h, k, l) (((h)+256)*512*512 + ((k)+256)*512 + ((l)+256))

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

static Reflection *recursive_search(Reflection *refl, unsigned int search)
{
	int i;

	/* Hit the bottom of the tree? */
	if ( refl == NULL ) return NULL;

	/* Is this the correct reflection? */
	if ( refl->serial == search ) return refl;

	for ( i=0; i<2; i++ ) {

		if ( search < refl->serial ) {
			return recursive_search(refl->child[0], search);
		}

		if ( search >= refl->serial ) {
			return recursive_search(refl->child[1], search);
		}

	}

	return NULL;
}


/* Return the first reflection in 'list' with the given indices, or NULL */
Reflection *find_refl(RefList *list, INDICES)
{
	unsigned int search = SERIAL(h, k, l);
	return recursive_search(list->head, search);
}


/* Find the next reflection in 'refl's list with the same indices, or NULL */
Reflection *next_found_refl(Reflection *refl)
{
	return refl->next;  /* Well, that was easy... */
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
	refl->excitation_error = exerr;
	refl->x = x;
	refl->y = y;
}


void set_partial(Reflection *refl, double r1, double r2, double p,
                 double clamp_low, double clamp_high)
{
	refl->r1 = r1;
	refl->r2 = r2;
	refl->p = p;
	refl->clamp1 = clamp_low;
	refl->clamp2 = clamp_high;
}


void set_indices(Reflection *refl, signed int h, signed int k, signed int l)
{
	/* Tricky, because the indices determine the position in the tree */
	Reflection copy;
	Reflection *new;
	Reflection transfer;

	/* Copy all data */
	copy = *refl;

	/* Delete and re-add with new indices */
	delete_refl(refl);
	new = add_refl(copy.list, h, k, l);

	/* Transfer data back */
	transfer = *new;
	*new = copy;
	new->list = transfer.list;
	new->parent = transfer.parent;
	new->child[0] = transfer.child[0];
	new->child[1] = transfer.child[1];
	new->h = transfer.h;  new->k = transfer.k;  new->l = transfer.l;
	new->serial = transfer.serial;
}


void set_int(Reflection *refl, double intensity)
{
	refl->intensity = intensity;
}


void set_scalable(Reflection *refl, int scalable)
{
	refl->scalable = scalable;
}


/********************************* Insertion **********************************/

static int recursive_insert(Reflection *refl, Reflection *new)
{
	int i;

	if ( refl->serial == new->serial ) {

		/* Found a reflection with identical indices */
		do {
			refl = refl->next;
		} while ( refl != NULL );
		refl->next = new;
		new->parent = refl;

		return 1;
	}

	for ( i=0; i<2; i++ ) {

		if ( new->serial < refl->serial ) {
			if ( refl->child[0] != NULL ) {
				return recursive_insert(refl->child[0], new);
			} else {
				refl->child[0] = new;
				new->parent = refl;
				return 1;
			}
		}

		if ( new->serial >= refl->serial ) {
			if ( refl->child[1] != NULL ) {
				return recursive_insert(refl->child[1], new);
			} else {
				refl->child[1] = new;
				new->parent = refl;
				return 1;
			}
		}


	}

	return 0;
}


Reflection *add_refl(RefList *list, INDICES)
{
	Reflection *new;

	new = calloc(1, sizeof(struct _reflection));
	new->h = h;  new->k = k,  new->l = l;
	new->serial = SERIAL(h, k, l);
	new->next = NULL;
	new->child[0] = NULL;
	new->child[1] = NULL;

	if ( list->head == NULL ) {
		list->head = new;
		new->parent = NULL;
	} else {
		recursive_insert(list->head, new);
	}

	return new;
}


/********************************** Deletion **********************************/

void delete_refl(Reflection *refl)
{
	int i;
	Reflection **parent_pos = NULL;

	/* Remove parent's reference */
	for ( i=0; i<2; i++ ) {
		if ( refl->parent->child[i] == refl ) {
			parent_pos = &refl->parent->child[i];
			*parent_pos = NULL;
		}
	}
	assert(parent_pos != NULL);

	/* Two child nodes? */
	if ( (refl->child[0] != NULL) && (refl->child[1] != NULL ) ) {

		if ( random() > RAND_MAX/2 ) {

			*parent_pos = refl->child[0];

			/* Now sort out the right child */
			recursive_insert(refl->child[0], refl->child[1]);
		} else {

			*parent_pos = refl->child[1];

			/* Now sort out the left child */
			recursive_insert(refl->child[1], refl->child[0]);

		}

	} else if ( refl->child[0] != NULL ) {

		/* One child, left */
		*parent_pos = refl->child[0];

	} else if (refl->child[1] != NULL ) {

		/* One child, right */
		*parent_pos = refl->child[1];

	} /* else it was just a leaf node */

	free(refl);
}


/********************************* Iteration **********************************/

Reflection *first_refl(RefList *list)
{
	return list->head;
}


Reflection *next_refl(Reflection *refl)
{
	/* Does a left child exist? */
	if ( refl->child[0] != NULL ) return refl->child[0];

	/* Otherwise move up the tree to find the next right child */
	while ( refl->child[1] != NULL ) {
		refl = refl->parent;
		if ( refl == NULL ) return NULL;
	}

	return refl->child[1];
}


/*********************************** Voodoo ***********************************/

void optimise_reflist(RefList *list)
{
}
