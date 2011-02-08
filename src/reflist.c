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
	unsigned int serial;          /* Unique serial number, key */
	struct _reflection *child[2]; /* Child nodes */
	struct _reflection *parent;   /* Parent node */
	struct _reflection *next;     /* Another reflection with the same
	                               * indices, or NULL */
	struct _reflection *prev;

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

static Reflection *new_node(unsigned int serial)
{
	Reflection *new;

	new = calloc(1, sizeof(struct _reflection));
	new->serial = serial;
	new->next = NULL;
	new->prev = NULL;
	new->child[0] = NULL;
	new->child[1] = NULL;

	return new;
}


/* Create a reflection list */
RefList *reflist_new()
{
	RefList *new;

	new = malloc(sizeof(struct _reflist));

	/* Create pseudo-root with invalid indices.
	 * The "real" root will be the left child of this. */
	new->head = new_node(SERIAL(257, 257, 257));

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

/* Return the first reflection in 'list' with the given indices, or NULL */
Reflection *find_refl(RefList *list, INDICES)
{
	unsigned int search = SERIAL(h, k, l);
	Reflection *refl = list->head->child[0];

	while ( refl != NULL ) {

		if ( search < refl->serial ) {

			if ( refl->child[0] != NULL ) {
				refl = refl->child[0];
			} else {
				/* Hit the bottom of the tree */
				return NULL;
			}

		} else if ( search > refl->serial ) {

			if ( refl->child[1] != NULL ) {
				refl = refl->child[1];
			} else {
				/* Hit the bottom of the tree */
				return NULL;
			}

		} else {

			return refl;

		}

	}

	return NULL;
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


void set_int(Reflection *refl, double intensity)
{
	refl->intensity = intensity;
}


void set_scalable(Reflection *refl, int scalable)
{
	refl->scalable = scalable;
}


/********************************* Insertion **********************************/

static void insert_node(Reflection *head, Reflection *new)
{
	Reflection *refl;

	refl = head;

	while ( refl != NULL ) {

		if ( new->serial < refl->serial ) {

			if ( refl->child[0] != NULL ) {
				refl = refl->child[0];
			} else {
				refl->child[0] = new;
				new->parent = refl;
				return;
			}

		} else if ( new->serial > refl->serial ) {

			if ( refl->child[1] != NULL ) {
				refl = refl->child[1];
			} else {
				refl->child[1] = new;
				new->parent = refl;
				return;
			}

		} else {

			/* New reflection is identical to a previous one */
			do {
				if ( refl->next == NULL ) {
					refl->next = new;
					new->prev = refl;
					return;
				}
				refl = refl->next;
			} while ( 1 );

		}

	}
}


Reflection *add_refl(RefList *list, INDICES)
{
	Reflection *new;

	new = new_node(SERIAL(h, k, l));
	new->h = h;  new->k = k,  new->l = l;

	if ( list->head == NULL ) {
		list->head = new;
		new->parent = NULL;
	} else {
		insert_node(list->head, new);
	}

	return new;
}


/********************************** Deletion **********************************/

void delete_refl(Reflection *refl)
{
	int i;
	Reflection **parent_pos = NULL;

	/* Is this a duplicate, and not the first one? */
	if ( refl->prev != NULL ) {
		refl->prev->next = refl->next;
		free(refl);
		return;
	}

	/* Is this the first duplicate of many? */
	if ( refl->next != NULL ) {

		refl->next->parent = refl->parent;
		refl->next->child[0] = refl->child[0];
		refl->next->child[1] = refl->child[1];
		refl->next->prev = NULL;

		for ( i=0; i<2; i++ ) {
			if ( refl->parent->child[i] == refl ) {
				refl->parent->child[i] = refl->next;
			}
		}

		free(refl);

		return;

	}

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
			refl->child[0]->parent = refl->parent;

			/* Now sort out the right child */
			insert_node(refl->child[0], refl->child[1]);

		} else {

			*parent_pos = refl->child[1];
			refl->child[1]->parent = refl->parent;

			/* Now sort out the left child */
			insert_node(refl->child[1], refl->child[0]);

		}

	} else if ( refl->child[0] != NULL ) {

		/* One child, left */
		*parent_pos = refl->child[0];
		refl->child[0]->parent = refl->parent;

	} else if (refl->child[1] != NULL ) {

		/* One child, right */
		*parent_pos = refl->child[1];
		refl->child[1]->parent = refl->parent;

	} /* else it was just a leaf node */

	free(refl);
}


/********************************* Iteration **********************************/


struct _reflistiterator {

	int stack_size;
	int stack_ptr;
	Reflection **stack;

};

Reflection *first_refl(RefList *list, RefListIterator **piter)
{
	RefListIterator *iter;

	iter = malloc(sizeof(struct _reflistiterator));
	iter->stack_size = 32;
	iter->stack = malloc(iter->stack_size*sizeof(Reflection *));
	iter->stack_ptr = 0;
	*piter = iter;

	Reflection *refl = list->head->child[0];

	do {

		if ( refl != NULL ) {
			iter->stack[iter->stack_ptr++] = refl;
			if ( iter->stack_ptr == iter->stack_size ) {
				iter->stack_size += 32;
				iter->stack = realloc(iter->stack,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;
		}

		if ( iter->stack_ptr == 0 ) {
			free(iter->stack);
			free(iter);
			return NULL;
		}

		refl = iter->stack[iter->stack_ptr--];

		return refl;

	} while ( 1 );
}


Reflection *next_refl(Reflection *refl, RefListIterator *iter)
{
	do {

		refl = refl->child[1];;

		if ( refl != NULL ) {
			iter->stack[iter->stack_ptr++] = refl;
			if ( iter->stack_ptr == iter->stack_size ) {
				iter->stack_size += 32;
				iter->stack = realloc(iter->stack,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;
		}
		if ( iter->stack_ptr == 0 ) {
			free(iter->stack);
			free(iter);
			return NULL;
		}
		refl = iter->stack[iter->stack_ptr--];

		return refl;

	} while ( 1 );
}


/*********************************** Voodoo ***********************************/

void optimise_reflist(RefList *list)
{
}
