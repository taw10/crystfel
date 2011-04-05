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
#include "utils.h"

/**
 * SECTION:reflist
 * @short_description: The fast reflection list
 * @title: RefList
 * @section_id:
 * @see_also:
 * @include: "reflist.h"
 * @Image:
 *
 * The fast reflection list stores reflections in a binary search tree indexed
 * by the Miller indices h, k and l.  Provided the tree has been optimised (by
 * using optimise_reflist()), any reflection can be found in a maximum length
 * of time which scales logarithmically with the number of reflections in the
 * list.
 *
 * A RefList can contain any number of reflections, and can store more than
 * one reflection with a given set of indices, for example when two distinct
 * reflections are to be stored according to their asymmetric indices.
 *
 * There are getters and setters which can be used to get and set values for an
 * individual reflection.  The reflection list does not calculate any values,
 * only stores what it was given earlier.  As such, you will need to carefully
 * examine which fields your prior processing steps have filled in.
 */


struct _refldata {

	/* Partiality and related geometrical stuff */
	double r1;  /* First excitation error */
	double r2;  /* Second excitation error */
	double p;   /* Partiality */
	int clamp1; /* Clamp status for r1 */
	int clamp2; /* Clamp status for r2 */

	/* Location in image */
	double fs;
	double ss;

	/* The distance from the exact Bragg position to the coordinates
	 * given above. */
	double excitation_error;

	/* Non-zero if this reflection can be used for scaling */
	int scalable;

	/* Intensity */
	double intensity;
	double esd_i;

	/* Phase */
	double phase;

	/* Redundancy */
	int redundancy;

	/* Total squared difference between all estimates of this reflection
	 * and the estimated mean value */
	double sum_squared_dev;
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


#define SERIAL(h, k, l) ((((h)+256)<<18) + (((k)+256)<<9) + ((l)+256))
#define GET_H(serial) ((((serial) & 0xfffc0000)>>18)-256)
#define GET_K(serial) ((((serial) & 0x0003fe00)>>9)-256)
#define GET_L(serial) (((serial) & 0x000001ff)-256)

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


/**
 * reflist_new:
 *
 * Creates a new reflection list.
 *
 * Returns: the new reflection list, or NULL on error.
 */
RefList *reflist_new()
{
	RefList *new;

	new = malloc(sizeof(struct _reflist));
	if ( new == NULL ) return NULL;

	/* Create pseudo-root with invalid indices.
	 * The "real" root will be the left child of this. */
	new->head = new_node(1<<31);

	return new;
}


static void recursive_free(Reflection *refl)
{
	if ( refl->child[0] != NULL ) recursive_free(refl->child[0]);
	if ( refl->child[1] != NULL ) recursive_free(refl->child[1]);
	free(refl);
}


/**
 * reflist_free:
 * @list: The reflection list to free.
 *
 * Destroys a reflection list.
 */
void reflist_free(RefList *list)
{
	if ( list == NULL ) return;
	recursive_free(list->head);
	free(list);
}


/********************************** Search ************************************/

/**
 * find_refl:
 * @list: The reflection list to search in
 * @h: The 'h' index to search for
 * @k: The 'k' index to search for
 * @l: The 'l' index to search for
 *
 * This function finds the first reflection in 'list' with the given indices.
 *
 * Since a %RefList can contain multiple reflections with the same indices, you
 * may need to use next_found_refl() to get the other reflections.
 *
 * Returns: The found reflection, or NULL if no reflection with the given
 * indices could be found.
 **/
Reflection *find_refl(const RefList *list, signed int h, signed int k, signed int l)
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
			assert(h == GET_H(refl->serial));
			assert(k == GET_K(refl->serial));
			assert(l == GET_L(refl->serial));
			return refl;

		}

	}

	return NULL;
}


/**
 * next_found_refl:
 * @refl: A reflection returned by find_refl() or next_found_refl()
 *
 * This function returns the next reflection in @refl's list with the same
 * indices.
 *
 * Returns: The found reflection, or NULL if there are no more reflections with
 * the same indices.
 **/
Reflection *next_found_refl(Reflection *refl)
{
	if ( refl->next != NULL ) assert(refl->serial == refl->next->serial);

	return refl->next;  /* Well, that was easy... */
}


/********************************** Getters ***********************************/

/**
 * get_excitation_error:
 * @refl: A %Reflection
 *
 * Returns: The excitation error for the reflection.
 **/
double get_excitation_error(const Reflection *refl)
{
	return refl->data.excitation_error;
}


/**
 * get_detector_pos:
 * @refl: A %Reflection
 * @fs: Location at which to store the fast scan offset of the reflection
 * @ss: Location at which to store the slow scan offset of the reflection
 *
 **/
void get_detector_pos(const Reflection *refl, double *fs, double *ss)
{
	*fs = refl->data.fs;
	*ss = refl->data.ss;
}


/**
 * get_indices:
 * @refl: A %Reflection
 * @h: Location at which to store the 'h' index of the reflection
 * @k: Location at which to store the 'k' index of the reflection
 * @l: Location at which to store the 'l' index of the reflection
 *
 **/
void get_indices(const Reflection *refl,
                 signed int *h, signed int *k, signed int *l)
{
	*h = GET_H(refl->serial);
	*k = GET_K(refl->serial);
	*l = GET_L(refl->serial);
}


/**
 * get_partiality:
 * @refl: A %Reflection
 *
 * Returns: The partiality of the reflection.
 **/
double get_partiality(const Reflection *refl)
{
	return refl->data.p;
}


/**
 * get_intensity:
 * @refl: A %Reflection
 *
 * Returns: The intensity of the reflection.
 **/
double get_intensity(const Reflection *refl)
{
	return refl->data.intensity;
}


/**
 * get_partial:
 * @refl: A %Reflection
 * @r1: Location at which to store the first excitation error
 * @r2: Location at which to store the second excitation error
 * @p: Location at which to store the partiality
 * @clamp_low: Location at which to store the first clamp status
 * @clamp_high: Location at which to store the second clamp status
 *
 * This function is used during post refinement (in conjunction with
 * set_partial()) to get access to the details of the partiality calculation.
 *
 **/
void get_partial(const Reflection *refl, double *r1, double *r2, double *p,
                 int *clamp_low, int *clamp_high)
{
	*r1 = refl->data.r1;
	*r2 = refl->data.r2;
	*p = get_partiality(refl);
	*clamp_low = refl->data.clamp1;
	*clamp_high = refl->data.clamp2;
}


/**
 * get_scalable:
 * @refl: A %Reflection
 *
 * Returns: non-zero if this reflection was marked as useful for scaling and
 * post refinement.
 *
 **/
int get_scalable(const Reflection *refl)
{
	return refl->data.scalable;
}


/**
 * get_redundancy:
 * @refl: A %Reflection
 *
 * The redundancy of the reflection is the number of measurements that have been
 * made of it.  Note that a redundancy of zero may have a special meaning, such
 * as that the reflection was impossible to integrate.  Note further that each
 * reflection in the list has its own redundancy, even if there are multiple
 * copies of the reflection in the list.  The total number of reflection
 * measurements should always be the sum of the redundancies in the entire list.
 *
 * Returns: the number of measurements of this reflection.
 *
 **/
int get_redundancy(const Reflection *refl)
{
	return refl->data.redundancy;
}


/**
 * get_sum_squared_dev:
 * @refl: A %Reflection
 *
 * The sum squared deviation is used to estimate the standard errors on the
 * intensities during 'Monte Carlo' merging.
 *
 * Returns: the sum of the squared deviations between the intensities and the
 * mean intensity from all measurements of the reflection (and probably its
 * symmetry equivalents according to some point group).
 *
 **/
double get_sum_squared_dev(const Reflection *refl)
{
	return refl->data.sum_squared_dev;
}


/**
 * get_esd_intensity:
 * @refl: A %Reflection
 *
 * Returns: the standard error in the intensity measurement (as returned by
 * get_intensity()) for this reflection.
 *
 **/
double get_esd_intensity(const Reflection *refl)
{
	return refl->data.esd_i;
}


/**
 * get_phase:
 * @refl: A %Reflection
 *
 * Returns: the phase for this reflection.
 *
 **/
double get_phase(const Reflection *refl)
{
	return refl->data.phase;
}


/********************************** Setters ***********************************/

/**
 * copy_data:
 * @to: %Reflection to copy data into
 * @from: %Reflection to copy data from
 *
 * This function is used to copy the data (which is everything listed above in
 * the list of getters and setters, apart from the indices themselves) from one
 * reflection to another.  This might be used when creating a new list from an
 * old one, perhaps using the asymmetric indices instead of the raw indices for
 * the new list.
 *
 **/
void copy_data(Reflection *to, const Reflection *from)
{
	memcpy(&to->data, &from->data, sizeof(struct _refldata));
}


void set_detector_pos(Reflection *refl, double exerr, double fs, double ss)
{
	refl->data.excitation_error = exerr;
	refl->data.fs = fs;
	refl->data.ss = ss;
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


void set_redundancy(Reflection *refl, int red)
{
	refl->data.redundancy = red;
}


void set_sum_squared_dev(Reflection *refl, double dev)
{
	refl->data.sum_squared_dev = dev;
}


void set_esd_intensity(Reflection *refl, double esd)
{
	refl->data.esd_i = esd;
}


void set_ph(Reflection *refl, double phase)
{
	refl->data.phase = phase;
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


Reflection *add_refl(RefList *list, signed int h, signed int k, signed int l)
{
	Reflection *new;

	new = new_node(SERIAL(h, k, l));

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

static int recursive_depth(Reflection *refl)
{
	int depth_left, depth_right;

	if ( refl == NULL ) return 0;

	depth_left = recursive_depth(refl->child[0]);
	depth_right = recursive_depth(refl->child[1]);

	return 1 + biggest(depth_left, depth_right);
}


static int recursive_count(Reflection *refl)
{
	int count_left, count_right;

	if ( refl == NULL ) return 0;

	count_left = recursive_count(refl->child[0]);
	count_right = recursive_count(refl->child[1]);

	return 1 + count_left + count_right;
}


static int tree_to_vine(Reflection *root)
{
	Reflection *vine_tail = root;
	Reflection *remainder = vine_tail->child[0];
	int size = 0;

	while ( remainder != NULL ) {

		if ( remainder->child[1] == NULL ) {
			vine_tail = remainder;
			remainder = remainder->child[0];
			size++;
		} else {
			Reflection *tmp = remainder->child[1];
			remainder->child[1] = tmp->child[0];
			if ( tmp->child[0] != NULL ) {
				tmp->child[0]->parent = remainder;
			}
			tmp->child[0] = remainder;
			if ( remainder != NULL ) remainder->parent = tmp;
			remainder = tmp;
			vine_tail->child[0] = tmp;
			if ( tmp != NULL ) tmp->parent = vine_tail;
		}

	}

	return size;
}


static void compress(Reflection *root, int count)
{
	Reflection *scan = root;
	int i;

	for ( i=1; i<=count; i++ ) {
		Reflection *child;
		child = scan->child[0];
		scan->child[0] = child->child[0];
		if ( child->child[0] != NULL ) {
			child->child[0]->parent = scan;
		}
		scan = scan->child[0];
		child->child[0] = scan->child[1];
		if ( scan->child[1] != NULL ) {
			scan->child[1]->parent = child;
		}
		scan->child[1] = child;
		if ( child != NULL ) {
			child->parent = scan;
		}
	}
}


static void vine_to_tree(Reflection *root, int size)
{
	int leaf_count = size + 1 - pow(2.0, floor(log(size+1)/log(2.0)));

	compress(root, leaf_count);
	size -= leaf_count;
	while ( size > 1 ) {
		compress(root, size / 2);
		size = size / 2;
	}
}


/**
 * optimise_reflist:
 * @list: The reflection list to optimise
 *
 * Optimises the ordering of reflections in the list such that the list can be
 * searched in the fastest possible way.
 *
 * This is a relatively expensive operation, so in typical usage you would call
 * it only after adding or removing many reflections from a list, when the list
 * is unlikely to be significantly modified for a long period of time.
 *
 * Note that only adding or deleting reflections may reduce the efficiency of
 * the list.  Changing the contents of the reflections (e.g. updating intensity
 * values) does not.
 **/
void optimise_reflist(RefList *list)
{
	int n_items;
	int size;
	const int verbose = 0;

	n_items = recursive_count(list->head->child[0]);
	if ( verbose ) {
		STATUS("Tree depth = %i\n",
		       recursive_depth(list->head->child[0]));
		STATUS("Number of items = %i\n", n_items);
		STATUS("Optimum depth = %5.2f\n", floor(log(n_items)/log(2.0)));
	}

	/* Now use the DSW algorithm to rebalance the tree */
	size = tree_to_vine(list->head);
	vine_to_tree(list->head, size);

	if ( verbose ) {
		STATUS("Tree depth after rebalancing = %i\n",
		       recursive_depth(list->head->child[0]));
	}
}


int num_reflections(RefList *list)
{
	return recursive_count(list->head->child[0]);
}
