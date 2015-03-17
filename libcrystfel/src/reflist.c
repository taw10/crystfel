/*
 * reflist.c
 *
 * Fast reflection/peak list
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2014 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <pthread.h>

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
	double rlow;  /* Low excitation error */
	double rhigh;  /* High excitation error */
	double p;   /* Partiality */
	double L;   /* Lorentz factor */

	/* Location in image */
	double fs;
	double ss;
	struct panel *panel;

	/* Non-zero if this reflection can be used for scaling */
	int scalable;

	/* Non-zero if this reflection should be used as a "guide star" for
	 * post refinement */
	int refinable;

	/* Intensity */
	double intensity;
	double esd_i;

	/* Phase */
	double phase;
	int have_phase;

	/* Redundancy */
	int redundancy;

	/* Peak height and mean background */
	double peak;
	double mean_bg;

	/* User-specified temporary values */
	double temp1;
	double temp2;
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
	pthread_mutex_t lock;         /* Protects the contents of "data" */
	struct _refldata data;
};


struct _reflist {

	struct _reflection *head;
	struct _reflection *tail;

};


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
	pthread_mutex_init(&new->lock, NULL);

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

	new->head = NULL;

	return new;
}


/**
 * reflection_new:
 * @h: The h index of the new reflection
 * @k: The k index of the new reflection
 * @l: The l index of the new reflection
 *
 * Creates a new individual reflection.  You'll probably want to use
 * add_refl_to_list() at some later point.
 */
Reflection *reflection_new(signed int h, signed int k, signed int l)
{
	return new_node(SERIAL(h, k, l));
}


/**
 * reflection_free:
 * @refl: The reflection to free.
 *
 * Destroys an individual reflection.
 */
void reflection_free(Reflection *refl)
{
	pthread_mutex_destroy(&refl->lock);
	free(refl);
}


static void recursive_free(Reflection *refl)
{
	if ( refl->child[0] != NULL ) recursive_free(refl->child[0]);
	if ( refl->child[1] != NULL ) recursive_free(refl->child[1]);

	while ( refl != NULL ) {
		Reflection *next = refl->next;
		reflection_free(refl);
		refl = next;
	}
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
	if ( list->head != NULL ) {
		recursive_free(list->head);
	} /* else empty list */
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

	/* Indices greater than or equal to 512 are filtered out when
	 * reflections are added, so don't even bother looking.
	 * (also, looking for such reflections causes trouble because the search
	 * serial number would be invalid) */
	if ( abs(h) >= 512 ) return NULL;
	if ( abs(k) >= 512 ) return NULL;
	if ( abs(l) >= 512 ) return NULL;

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
 * get_panel:
 * @refl: A %Reflection
 *
 * Returns: the panel which the reflection appears on
 *
 **/
struct panel *get_panel(const Reflection *refl)
{
	return refl->data.panel;
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
 * Returns: The partiality of the reflection.  See get_lorentz().
 **/
double get_partiality(const Reflection *refl)
{
	return refl->data.p;
}


/**
 * get_lorentz:
 * @refl: A %Reflection
 *
 * Returns: The Lorentz factor for the reflection.  To "scale up" a partial
 * reflection, divide by this multiplied by the partiality.
 **/
double get_lorentz(const Reflection *refl)
{
	return refl->data.L;
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
 * @rlow: Location at which to store the "low" excitation error
 * @rhigh: Location at which to store the "high" excitation error
 * @p: Location at which to store the partiality
 *
 * This function is used during post refinement (in conjunction with
 * set_partial()) to get access to the details of the partiality calculation.
 *
 **/
void get_partial(const Reflection *refl, double *rlow, double *rhigh,
                 double *p)
{
	*rlow = refl->data.rlow;
	*rhigh = refl->data.rhigh;
	*p = get_partiality(refl);
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
 * @have_phase: Place to store a non-zero value if the phase is set, or NULL.
 *
 * Returns: the phase for this reflection.
 *
 **/
double get_phase(const Reflection *refl, int *have_phase)
{
	if ( have_phase != NULL ) *have_phase = refl->data.have_phase;
	return refl->data.phase;
}


/**
 * get_peak:
 * @refl: A %Reflection
 *
 * Returns: the peak height (value of the highest pixel, before background
 * subtraction) for this reflection.
 *
 **/
double get_peak(const Reflection *refl)
{
	return refl->data.peak;
}


/**
 * get_mean_bg:
 * @refl: A %Reflection
 *
 * Returns: the mean background level for this reflection.
 *
 **/
double get_mean_bg(const Reflection *refl)
{
	return refl->data.mean_bg;
}


/**
 * get_temp1:
 * @refl: A %Reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 * Returns: the first temporary value for this reflection.
 *
 **/
double get_temp1(const Reflection *refl)
{
	return refl->data.temp1;
}


/**
 * get_temp2:
 * @refl: A %Reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 * Returns: the second temporary value for this reflection.
 *
 **/
double get_temp2(const Reflection *refl)
{
	return refl->data.temp2;
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
 * @fs: The fast scan offset of the reflection
 * @ss: The slow scan offset of the reflection
 *
 **/
void set_detector_pos(Reflection *refl, double fs, double ss)
{
	refl->data.fs = fs;
	refl->data.ss = ss;
}


/**
 * set_panel:
 * @refl: A %Reflection
 * @panel: Pointer to the panel structure on which the reflection appears
 *
 * Note that the pointer will be stored, not the contents of the structure.
 *
 **/
void set_panel(Reflection *refl, struct panel *p)
{
	refl->data.panel = p;
}


/**
 * set_partial:
 * @refl: A %Reflection
 * @rlow: The "low" excitation error
 * @rhigh: The "high" excitation error
 * @p: The partiality
 *
 * This function is used during post refinement (in conjunction with
 * get_partial()) to get access to the details of the partiality calculation.
 *
 **/
void set_partial(Reflection *refl, double rlow, double rhigh, double p)
{
	refl->data.rlow = rlow;
	refl->data.rhigh = rhigh;
	refl->data.p = p;
}


/**
 * set_intensity:
 * @refl: A %Reflection
 * @p: The partiality for the reflection.
 *
 * Set the partiality for the reflection.  See set_lorentz().
 **/
void set_partiality(Reflection *refl, double p)
{
	refl->data.p = p;
}

/**
 * set_lorentz:
 * @refl: A %Reflection
 * @L: The Lorentz factor for the reflection.
 *
 * Set the Lorentz factor for the reflection.  To "scale up" a partial
 * reflection, divide by this multiplied by the partiality.
 **/
void set_lorentz(Reflection *refl, double L)
{
	refl->data.L = L;
}


/**
 * set_intensity:
 * @refl: A %Reflection
 * @intensity: The intensity for the reflection.
 *
 * Set the intensity for the reflection.
 **/
void set_intensity(Reflection *refl, double intensity)
{
	refl->data.intensity = intensity;
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
 * set_phase:
 * @refl: A %Reflection
 * @phase: New phase for the reflection
 *
 **/
void set_phase(Reflection *refl, double phase)
{
	refl->data.phase = phase;
	refl->data.have_phase = 1;
}


/**
 * set_peak:
 * @refl: A %Reflection
 * @peak: New peak height for the reflection
 *
 **/
void set_peak(Reflection *refl, double peak)
{
	refl->data.peak = peak;
}


/**
 * set_mean_bg:
 * @refl: A %Reflection
 * @mean_bg: New peak height for the reflection
 *
 **/
void set_mean_bg(Reflection *refl, double mean_bg)
{
	refl->data.mean_bg = mean_bg;
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


/**
 * set_temp1
 * @refl: A %Reflection
 * @temp: New temporary value for the reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 **/
void set_temp1(Reflection *refl, double temp)
{
	refl->data.temp1 = temp;
}


/**
 * set_temp2
 * @refl: A %Reflection
 * @temp: New temporary value for the reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 **/
void set_temp2(Reflection *refl, double temp)
{
	refl->data.temp2 = temp;
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


static void add_to_list(RefList *list, Reflection *new,
                        signed int h, signed int k, signed int l)
{
	Reflection *f;

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

	assert(abs(h)<512);
	assert(abs(k)<512);
	assert(abs(l)<512);

	new = new_node(SERIAL(h, k, l));
	if ( new == NULL ) return NULL;

	add_to_list(list, new, h, k, l);

	return new;
}


/**
 * add_refl_to_list
 * @refl: A %Reflection
 * @list: A %RefList
 *
 * Adds a @refl to @list.
 *
 **/
void add_refl_to_list(Reflection *refl, RefList *list)
{
	signed int h, k, l;

	get_indices(refl, &h, &k, &l);

	add_to_list(list, refl, h, k, l);
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
	Reflection *refl;
	RefListIterator *iter;

	iter = malloc(sizeof(struct _reflistiterator));
	iter->stack_size = 32;
	iter->stack = malloc(iter->stack_size*sizeof(Reflection *));
	iter->stack_ptr = 0;
	*piter = iter;

	if ( list == NULL ) return NULL;

	refl = list->head;

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
	/* Are there more reflections with the same indices? */
	if ( refl->next != NULL ) {
		return refl->next;
	} else {

		/* No, so rewind back to the head of the list */
		while ( refl->prev != NULL ) {
			refl = refl->prev;
		}

	}

	refl = refl->child[1];
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
	Reflection *probe;
	int n_this = 1;

	if ( refl == NULL ) return 0;

	probe = refl;
	while ( probe->next != NULL ) {
		probe = probe->next;
		n_this++;
	}

	count_left = recursive_count(refl->child[0]);
	count_right = recursive_count(refl->child[1]);

	return n_this + count_left + count_right;
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


/**
 * lock_reflection:
 * @refl: A %Reflection
 *
 * Acquires a lock on the reflection.
 */
void lock_reflection(Reflection *refl)
{
	pthread_mutex_lock(&refl->lock);
}


/**
 * unlock_reflection:
 * @refl: A %Reflection
 *
 * Releases a lock on the reflection.
 */
void unlock_reflection(Reflection *refl)
{
	pthread_mutex_unlock(&refl->lock);
}
