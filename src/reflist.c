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
 * The fast reflection list stores reflections in an RB-tree indexed
 * by the Miller indices h, k and l.  Any reflection can be found in a
 * length of time which scales logarithmically with the number of reflections in
 * the list.
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

	/* Symmetric indices (i.e. the "real" indices) */
	signed int hs;
	signed int ks;
	signed int ls;

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


enum _nodecol {
	RED,
	BLACK
};


struct _reflection {

	/* Listy stuff */
	unsigned int serial;          /* Unique serial number, key */
	struct _reflection *child[2]; /* Child nodes */
	struct _reflection *next;     /* Next and previous in doubly linked */
	struct _reflection *prev;     /*  list of duplicate reflections */
	enum _nodecol col;            /* Colour (red or black) */

	/* Payload */
	struct _refldata data;
};


struct _reflist {

	struct _reflection *head;
	struct _reflection *tail;

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
	new->col = RED;

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
	new->head = NULL;//new_node(1<<31);

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
Reflection *find_refl(const RefList *list,
                      signed int h, signed int k, signed int l)
{
	unsigned int search = SERIAL(h, k, l);
	Reflection *refl;

	if ( list->head == NULL ) return NULL;

	refl = list->head;

	while ( refl != NULL ) {

		if ( refl->serial == search ) {

			assert(search == refl->serial);
			assert(h == GET_H(refl->serial));
			assert(k == GET_K(refl->serial));
			assert(l == GET_L(refl->serial));
			return refl;

		} else {

			int dir = search > refl->serial;
			if ( refl->child[dir] != NULL ) {
				refl = refl->child[dir];
			} else {
				/* Hit the bottom of the tree */
				return NULL;
			}

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
 * get_symmetric_indices:
 * @refl: A %Reflection
 * @hs: Location at which to store the 'h' index of the reflection
 * @ks: Location at which to store the 'k' index of the reflection
 * @ls: Location at which to store the 'l' index of the reflection
 *
 * This function gives the symmetric indices, that is, the "real" indices before
 * squashing down to the asymmetric reciprocal unit.  This may be useful if the
 * list is indexed according to the asymmetric indices, but you still need
 * access to the symmetric version.  This happens during post-refinement.
 *
 **/
void get_symmetric_indices(const Reflection *refl,
                                  signed int *hs, signed int *ks,
                                  signed int *ls)
{
	*hs = refl->data.hs;
	*ks = refl->data.ks;
	*ls = refl->data.ls;
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


/**
 * set_detector_pos:
 * @refl: A %Reflection
 * @exerr: The excitation error for this reflection
 * @fs: The fast scan offset of the reflection
 * @ss: The slow scan offset of the reflection
 *
 **/
void set_detector_pos(Reflection *refl, double exerr, double fs, double ss)
{
	refl->data.excitation_error = exerr;
	refl->data.fs = fs;
	refl->data.ss = ss;
}


/**
 * set_partial:
 * @refl: A %Reflection
 * @r1: The first excitation error
 * @r2: The second excitation error
 * @p: The partiality
 * @clamp_low: The first clamp status
 * @clamp_high: The second clamp status
 *
 * This function is used during post refinement (in conjunction with
 * get_partial()) to get access to the details of the partiality calculation.
 *
 **/
void set_partial(Reflection *refl, double r1, double r2, double p,
                 double clamp_low, double clamp_high)
{
	refl->data.r1 = r1;
	refl->data.r2 = r2;
	refl->data.p = p;
	refl->data.clamp1 = clamp_low;
	refl->data.clamp2 = clamp_high;
}


/**
 * set_int:
 * @refl: A %Reflection
 * @intensity: The intensity for the reflection.
 *
 * Set the intensity for the reflection.  Note that retrieval is done with
 * get_intensity().
 **/
void set_int(Reflection *refl, double intensity)
{
	refl->data.intensity = intensity;
}


/**
 * set_scalable:
 * @refl: A %Reflection
 * @scalable: Non-zero if this reflection was marked as useful for scaling and
 * post refinement.
 *
 **/
void set_scalable(Reflection *refl, int scalable)
{
	refl->data.scalable = scalable;
}


/**
 * set_redundancy:
 * @refl: A %Reflection
 * @red: New redundancy for the reflection
 *
 * The redundancy of the reflection is the number of measurements that have been
 * made of it.  Note that a redundancy of zero may have a special meaning, such
 * as that the reflection was impossible to integrate.  Note further that each
 * reflection in the list has its own redundancy, even if there are multiple
 * copies of the reflection in the list.  The total number of reflection
 * measurements should always be the sum of the redundancies in the entire list.
 *
 **/
void set_redundancy(Reflection *refl, int red)
{
	refl->data.redundancy = red;
}


/**
 * set_sum_squared_dev:
 * @refl: A %Reflection
 * @dev: New sum squared deviation for the reflection
 *
 * The sum squared deviation is used to estimate the standard errors on the
 * intensities during 'Monte Carlo' merging.  It is defined as the sum of the
 * squared deviations between the intensities and the mean intensity from all
 * measurements of the reflection (and probably its symmetry equivalents
 * according to some point group).
 *
 **/
void set_sum_squared_dev(Reflection *refl, double dev)
{
	refl->data.sum_squared_dev = dev;
}


/**
 * set_esd_intensity:
 * @refl: A %Reflection
 * @esd: New standard error for this reflection's intensity measurement
 *
 **/
void set_esd_intensity(Reflection *refl, double esd)
{
	refl->data.esd_i = esd;
}


/**
 * set_ph:
 * @refl: A %Reflection
 * @phase: New phase for the reflection
 *
 **/
void set_ph(Reflection *refl, double phase)
{
	refl->data.phase = phase;
}


/**
 * set_symmetric_indices:
 * @refl: A %Reflection
 * @hs: The 'h' index of the reflection
 * @ks: The 'k' index of the reflection
 * @ls: The 'l' index of the reflection
 *
 * This function gives the symmetric indices, that is, the "real" indices before
 * squashing down to the asymmetric reciprocal unit.  This may be useful if the
 * list is indexed according to the asymmetric indices, but you still need
 * access to the symmetric version.  This happens during post-refinement.
 *
 **/
void set_symmetric_indices(Reflection *refl,
                           signed int hs, signed int ks, signed int ls)
{
	refl->data.hs = hs;
	refl->data.ks = ks;
	refl->data.ls = ls;
}


/********************************* Insertion **********************************/

static Reflection *rotate_once(Reflection *refl, int dir)
{
	Reflection *s = refl->child[!dir];

	refl->child[!dir] = s->child[dir];
	s->child[dir] = refl;

	refl->col = RED;
	s->col = BLACK;

	return s;
}


static Reflection *rotate_twice(Reflection *refl, int dir)
{
	refl->child[!dir] = rotate_once(refl->child[!dir], !dir);
	return rotate_once(refl, dir);
}


static int is_red(Reflection *refl)
{
	return (refl != NULL) && (refl->col == RED);
}


static Reflection *insert_node(Reflection *refl, Reflection *new)
{
	if ( refl == NULL ) {

		refl = new;

	} else {

		int dir;
		Reflection *rcd;

		assert(new->serial != refl->serial);
		dir = new->serial > refl->serial;
		refl->child[dir] = insert_node(refl->child[dir], new);

		rcd = refl->child[dir];
		if ( is_red(rcd) ) {

			if ( is_red(refl->child[!dir]) ) {

				refl->col = RED;
				refl->child[0]->col = BLACK;
				refl->child[1]->col = BLACK;

			} else {

				if ( is_red(rcd->child[dir] ) ) {
					refl = rotate_once(refl, !dir);
				} else if ( is_red(rcd->child[!dir] ) ) {
					refl = rotate_twice(refl, !dir);
				}

			}
		}

	}

	return refl;
}


/**
 * add_refl
 * @list: A %RefList
 * @h: The 'h' index of the reflection
 * @k: The 'k' index of the reflection
 * @l: The 'l' index of the reflection
 *
 * Adds a new reflection to @list.  Note that the implementation allows there to
 * be multiple reflections with the same indices in the list, so this function
 * should succeed even if the given indices already feature in the list.
 *
 * Returns: The newly created reflection, or NULL on failure.
 *
 **/
Reflection *add_refl(RefList *list, signed int h, signed int k, signed int l)
{
	Reflection *new;
	Reflection *f;

	assert(abs(h)<256);
	assert(abs(k)<256);
	assert(abs(l)<256);

	new = new_node(SERIAL(h, k, l));
	if ( new == NULL ) return NULL;

	f = find_refl(list, h, k, l);
	if ( f == NULL ) {

		list->head = insert_node(list->head, new);
		list->head->col = BLACK;

	} else {

		/* New reflection is identical to a previous one */
		while ( f->next != NULL ) {
			f = f->next;
		}
		f->next = new;
		new->prev = f;
	}

	return new;
}


/********************************* Iteration **********************************/

struct _reflistiterator {

	int stack_size;
	int stack_ptr;
	Reflection **stack;

};


/**
 * first_refl:
 * @list: A %RefList to iterate over
 * @piter: Address at which to store a %RefListIterator
 *
 * This function sets up the state required for iteration over the entire list,
 * and then returns the first reflection in the list.  An iterator object will
 * be created and its address stored at the location given in piter.
 *
 * Returns: the first reflection in the list.
 *
 **/
Reflection *first_refl(RefList *list, RefListIterator **piter)
{
	RefListIterator *iter;

	iter = malloc(sizeof(struct _reflistiterator));
	iter->stack_size = 32;
	iter->stack = malloc(iter->stack_size*sizeof(Reflection *));
	iter->stack_ptr = 0;
	*piter = iter;

	Reflection *refl = list->head;

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


/**
 * next_refl:
 * @refl: A reflection
 * @iter: A %RefListIterator
 *
 * This function looks up the next reflection in the list that was given earlier
 * to first_refl().
 *
 * Returns: the next reflection in the list, or NULL if no more.
 *
 **/
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


/**
 * num_reflections:
 * @list: A %RefList
 *
 * Returns: the number of reflections in @list.
 *
 **/
int num_reflections(RefList *list)
{
	return recursive_count(list->head);
}


/**
 * tree_depth:
 * @list: A %RefList
 *
 * If the depth of the tree is more than about 20, access to the list will be
 * slow.  This should never happen.
 *
 * Returns: the depth of the RB-tree used internally to represent @list.
 *
 **/
int tree_depth(RefList *list)
{
	return recursive_depth(list->head);
}
