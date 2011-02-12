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
#include <stdio.h>

#include "reflist.h"


struct _refldata {

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


struct _reflection {

	/* Listy stuff */
	unsigned int serial;          /* Unique serial number, key */
	struct _reflection *child[2]; /* Child nodes */
	struct _reflection *parent;   /* Parent node */
	struct _reflection *next;     /* Next and previous in doubly linked */
	struct _reflection *prev;     /*  list of duplicate reflections */

	/* Payload */
	struct _refldata data;
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
	if ( list == NULL ) return;
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

			assert(search == refl->serial);
			assert(h == refl->data.h);
			assert(k == refl->data.k);
			assert(l == refl->data.l);
			return refl;

		}

	}

	return NULL;
}


/* Find the next reflection in 'refl's list with the same indices, or NULL */
Reflection *next_found_refl(Reflection *refl)
{
	if ( refl->next != NULL ) assert(refl->serial == refl->next->serial);

	return refl->next;  /* Well, that was easy... */
}


/********************************** Getters ***********************************/

double get_excitation_error(Reflection *refl)
{
	return refl->data.excitation_error;
}


void get_detector_pos(Reflection *refl, double *x, double *y)
{
	*x = refl->data.x;
	*y = refl->data.y;
}


void get_indices(Reflection *refl, signed int *h, signed int *k, signed int *l)
{
	*h = refl->data.h;
	*k = refl->data.k;
	*l = refl->data.l;
}


double get_partiality(Reflection *refl)
{
	return refl->data.p;
}


double get_intensity(Reflection *refl)
{
	return refl->data.intensity;
}


void get_partial(Reflection *refl, double *r1, double *r2, double *p,
                 int *clamp_low, int *clamp_high)
{
	*r1 = refl->data.r1;
	*r2 = refl->data.r2;
	*p = get_partiality(refl);
	*clamp_low = refl->data.clamp1;
	*clamp_high = refl->data.clamp2;
}


int get_scalable(Reflection *refl)
{
	return refl->data.scalable;
}


/********************************** Setters ***********************************/

void set_detector_pos(Reflection *refl, double exerr, double x, double y)
{
	refl->data.excitation_error = exerr;
	refl->data.x = x;
	refl->data.y = y;
}


void set_partial(Reflection *refl, double r1, double r2, double p,
                 double clamp_low, double clamp_high)
{
	refl->data.r1 = r1;
	refl->data.r2 = r2;
	refl->data.p = p;
	refl->data.clamp1 = clamp_low;
	refl->data.clamp2 = clamp_high;
}


void set_int(Reflection *refl, double intensity)
{
	refl->data.intensity = intensity;
}


void set_scalable(Reflection *refl, int scalable)
{
	refl->data.scalable = scalable;
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
			assert(refl->serial == new->serial);
			while ( refl->next != NULL ) {
				refl = refl->next;
			}
			refl->next = new;
			new->prev = refl;
			return;

		}

	}
}


Reflection *add_refl(RefList *list, INDICES)
{
	Reflection *new;

	new = new_node(SERIAL(h, k, l));
	new->data.h = h;  new->data.k = k,  new->data.l = l;

	if ( list->head == NULL ) {
		list->head = new;
		new->parent = NULL;
	} else {
		insert_node(list->head, new);
	}

	return new;
}


/********************************** Deletion **********************************/

static void lr_delete(Reflection *refl, int side)
{
	int other = 1-side;
	int i;
	Reflection *pre;

	pre = refl->child[side];
	while ( pre->child[other] != NULL ) pre = pre->child[other];

	assert(refl->next == NULL);
	assert(refl->prev == NULL); /* Should have been caught previously */

	refl->data = pre->data;
	refl->serial = pre->serial;

	/* If the predecessor node had duplicates, we need to fix things up. */
	assert(pre->prev == NULL);
	refl->next = pre->next;
	if ( pre->next != NULL ) {
		refl->next->prev = refl;
	}

	for ( i=0; i<2; i++ ) {
		if ( pre->parent->child[i] == pre ) {
			pre->parent->child[i] = pre->child[side];
		}
	}
	if ( pre->child[side] != NULL ) {
		pre->child[side]->parent = pre->parent;
	}
	free(pre);
}


void delete_refl(Reflection *refl)
{
	int i;

	/* Is this a duplicate, and not the first one? */
	if ( refl->prev != NULL ) {
		refl->prev->next = refl->next;
		if ( refl->next != NULL ) refl->next->prev = refl->prev;
		free(refl);
		return;
	}

	/* Is this the first duplicate of many? */
	if ( refl->next != NULL ) {

		assert(refl->next->prev == refl);
		assert(refl->prev == NULL);
		refl->next->parent = refl->parent;
		refl->next->prev = NULL;

		for ( i=0; i<2; i++ ) {
			refl->next->child[i] = refl->child[i];
			if ( refl->parent->child[i] == refl ) {
				refl->parent->child[i] = refl->next;
			}
			if ( refl->child[i] != NULL ) {
				refl->child[i]->parent = refl->next;
			}
		}

		free(refl);

		return;

	}

	assert(refl->next == NULL);
	assert(refl->prev == NULL);

	/* Two child nodes? */
	if ( (refl->child[0] != NULL) && (refl->child[1] != NULL ) ) {

		if ( random() > RAND_MAX/2 ) {
			lr_delete(refl, 0);
		} else {
			lr_delete(refl, 1);
		}

	} else if ( refl->child[0] != NULL ) {

		/* One child, left */
		for ( i=0; i<2; i++ ) {
			if ( refl->parent->child[i] == refl ) {
				refl->parent->child[i] = refl->child[0];
			}
		}
		refl->child[0]->parent = refl->parent;
		free(refl);

	} else if (refl->child[1] != NULL ) {

		/* One child, right */
		for ( i=0; i<2; i++ ) {
			if ( refl->parent->child[i] == refl ) {
				refl->parent->child[i] = refl->child[1];
			}
		}
		refl->child[1]->parent = refl->parent;
		free(refl);

	} else {

		/* Leaf node */
		for ( i=0; i<2; i++ ) {
			if ( refl->parent->child[i] == refl ) {
				refl->parent->child[i] = NULL;
			}
		}
		free(refl);

	}
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

		refl = iter->stack[--iter->stack_ptr];

		return refl;

	} while ( 1 );
}


Reflection *next_refl(Reflection *refl, RefListIterator *iter)
{
	int returned = 1;

	do {

		if ( returned ) refl = refl->child[1];
		returned = 0;

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

		return iter->stack[--iter->stack_ptr];

	} while ( 1 );
}


/*********************************** Voodoo ***********************************/

void optimise_reflist(RefList *list)
{
}
