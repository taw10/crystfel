/*
 * reflist.c
 *
 * Fast reflection/peak list
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2021 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <pthread.h>

#include "reflist.h"
#include "utils.h"

/** \file reflist.h */

struct _refldata {

	/* Symmetric indices (i.e. the "real" indices) */
	signed int hs;
	signed int ks;
	signed int ls;

	/* Partiality and related geometrical stuff */
	double khalf; /* Wavenumber of middle of reflection */
	double kpred; /* Wavenumber for prediction */
	double exerr; /* Excitation error */
	double p;     /* Partiality */
	double L;     /* Lorentz factor */

	/* Location in image */
	double fs;
	double ss;
	int panel_number;

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

	/* Contributions */
	struct reflection_contributions *contribs;

	/* User-specified temporary values */
	double temp1;
	double temp2;
	int flag;
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
	int in_list;                  /* If 0, reflection is not in a list */

	/* Payload */
	pthread_mutex_t lock;         /* Protects the contents of "data" */
	struct _refldata data;
};


struct _reflist {

	struct _reflection *head;
	char *notes;

};


/**************************** Creation / deletion *****************************/

static Reflection *new_node(unsigned int serial)
{
	Reflection *new;

	new = cfcalloc(1, sizeof(struct _reflection));
	if ( new == NULL ) return NULL;
	new->in_list = 0;
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
 * Creates a new reflection list.
 *
 * \returns the new reflection list, or NULL on error.
 */
RefList *reflist_new()
{
	RefList *new;

	new = cfmalloc(sizeof(struct _reflist));
	if ( new == NULL ) return NULL;

	new->head = NULL;
	new->notes = NULL;

	return new;
}


/**
 * \param h The h index of the new reflection
 * \param k The k index of the new reflection
 * \param l The l index of the new reflection
 *
 * Creates a new individual reflection.  You'll probably want to use
 * add_refl_to_list() at some later point.
 */
Reflection *reflection_new(signed int h, signed int k, signed int l)
{
	assert(abs(h)<512);
	assert(abs(k)<512);
	assert(abs(l)<512);
	return new_node(SERIAL(h, k, l));
}


/**
 * \param refl: The reflection to free.
 *
 * Destroys an individual reflection.
 */
void reflection_free(Reflection *refl)
{
	pthread_mutex_destroy(&refl->lock);
	cffree(refl);
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
 * \param list: The reflection list to free.
 *
 * Destroys a reflection list.
 */
void reflist_free(RefList *list)
{
	if ( list == NULL ) return;
	if ( list->head != NULL ) {
		recursive_free(list->head);
	} /* else empty list */
	if ( list->notes != NULL ) cffree(list->notes);
	cffree(list);
}


/********************************** Search ************************************/

/**
 * \param list: The reflection list to search in
 * \param h: The 'h' index to search for
 * \param k: The 'k' index to search for
 * \param l: The 'l' index to search for
 *
 * This function finds the first reflection in 'list' with the given indices.
 *
 * Since a %RefList can contain multiple reflections with the same indices, you
 * may need to use next_found_refl() to get the other reflections.
 *
 * \returns The found reflection, or NULL if no reflection with the given
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
 * \param refl: A reflection returned by find_refl() or next_found_refl()
 *
 * This function returns the next reflection in \p refl's list with the same
 * indices.
 *
 * \returns The found reflection, or NULL if there are no more reflections with
 * the same indices.
 **/
Reflection *next_found_refl(Reflection *refl)
{
	if ( refl->next != NULL ) assert(refl->serial == refl->next->serial);

	return refl->next;  /* Well, that was easy... */
}


/********************************** Getters ***********************************/


/**
 * \param refl: Reflection
 * \param fs: Location at which to store the fast scan offset of the reflection
 * \param ss: Location at which to store the slow scan offset of the reflection
 *
 **/
void get_detector_pos(const Reflection *refl, double *fs, double *ss)
{
	*fs = refl->data.fs;
	*ss = refl->data.ss;
}


/**
 * \param refl: Reflection
 *
 * \returns panel number (index in detgeom/DataTemplate structure)
 *           which the reflection appears on
 *
 **/
int get_panel_number(const Reflection *refl)
{
	return refl->data.panel_number;
}


/**
 * \param refl: Reflection
 * \param h: Location at which to store the 'h' index of the reflection
 * \param k: Location at which to store the 'k' index of the reflection
 * \param l: Location at which to store the 'l' index of the reflection
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
 * \param refl: Reflection
 * \param hs: Location at which to store the 'h' index of the reflection
 * \param ks: Location at which to store the 'k' index of the reflection
 * \param ls: Location at which to store the 'l' index of the reflection
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
 * \param refl: Reflection
 *
 * \returns The partiality of the reflection.  See get_lorentz().
 **/
double get_partiality(const Reflection *refl)
{
	return refl->data.p;
}


/**
 * \param refl: Reflection
 *
 * \returns The Lorentz factor for the reflection.  To "scale up" a partial
 * reflection, divide by this multiplied by the partiality.
 **/
double get_lorentz(const Reflection *refl)
{
	return refl->data.L;
}


/**
 * \param refl: Reflection
 *
 * \returns The intensity of the reflection.
 **/
double get_intensity(const Reflection *refl)
{
	return refl->data.intensity;
}


/**
 * get_khalf
 * \param refl: Reflection
 *
 * \returns the wavenumber at the centre of the reflection
 *
 **/
double get_khalf(const Reflection *refl)
{
	return refl->data.khalf;
}




/**
 * \param refl: Reflection
 *
 * \returns the wavenumber which should be used for prediction of this reflection
 *
 **/
double get_kpred(const Reflection *refl)
{
	return refl->data.kpred;
}


/**
 * \param refl: Reflection
 *
 * \returns the excitation error (in m^-1) for this reflection
 *
 **/
double get_exerr(const Reflection *refl)
{
	return refl->data.exerr;
}


/**
 * \param refl: Reflection
 *
 * The redundancy of the reflection is the number of measurements that have been
 * made of it.  Note that a redundancy of zero may have a special meaning, such
 * as that the reflection was impossible to integrate.  Note further that each
 * reflection in the list has its own redundancy, even if there are multiple
 * copies of the reflection in the list.  The total number of reflection
 * measurements should always be the sum of the redundancies in the entire list.
 *
 * \returns the number of measurements of this reflection.
 *
 **/
int get_redundancy(const Reflection *refl)
{
	return refl->data.redundancy;
}


/**
 * \param refl: Reflection
 *
 * \returns the standard error in the intensity measurement (as returned by
 * get_intensity()) for this reflection.
 *
 **/
double get_esd_intensity(const Reflection *refl)
{
	return refl->data.esd_i;
}


/**
 * \param refl: Reflection
 * \param have_phase: Place to store a non-zero value if the phase is set, or NULL.
 *
 * \returns the phase for this reflection.
 *
 **/
double get_phase(const Reflection *refl, int *have_phase)
{
	if ( have_phase != NULL ) *have_phase = refl->data.have_phase;
	return refl->data.phase;
}


/**
 * \param refl: Reflection
 *
 * \returns the peak height (value of the highest pixel, before background
 * subtraction) for this reflection.
 *
 **/
double get_peak(const Reflection *refl)
{
	return refl->data.peak;
}


/**
 * \param refl: Reflection
 *
 * \returns the mean background level for this reflection.
 *
 **/
double get_mean_bg(const Reflection *refl)
{
	return refl->data.mean_bg;
}


/**
 * \param refl: Reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 * \returns the first temporary value for this reflection.
 *
 **/
double get_temp1(const Reflection *refl)
{
	return refl->data.temp1;
}


/**
 * \param refl: Reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 * \returns the second temporary value for this reflection.
 *
 **/
double get_temp2(const Reflection *refl)
{
	return refl->data.temp2;
}


/**
 * \param refl: Reflection
 *
 * The integer flag value can be used according to the needs of the calling
 * program.
 *
 * \returns the flag for this reflection.
 *
 **/
int get_flag(const Reflection *refl)
{
	return refl->data.flag;
}


/**
 * \param refl: Reflection
 *
 * \returns the reflection's contribution list
 *
 **/
struct reflection_contributions *get_contributions(const Reflection *refl)
{
	return refl->data.contribs;
}

/********************************** Setters ***********************************/

/**
 * \param to: %Reflection to copy data into
 * \param from: %Reflection to copy data from
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
 * \param refl: Reflection
 * \param fs: The fast scan offset of the reflection
 * \param ss: The slow scan offset of the reflection
 *
 **/
void set_detector_pos(Reflection *refl, double fs, double ss)
{
	refl->data.fs = fs;
	refl->data.ss = ss;
}


/**
 * \param refl: Reflection
 * \param pn: Panel number (index in detgeom/DataTemplate structure) of
 *      the panel on which the reflection appears.
 *
 **/
void set_panel_number(Reflection *refl, int pn)
{
	refl->data.panel_number = pn;
}


/**
 * \param refl: Reflection
 * \param khalf: The wavenumber at which the reflection should be predicted
 *
 * Sets the wavenumber at the centre of the reflection.
 **/
void set_khalf(Reflection *refl, double khalf)
{
	refl->data.khalf = khalf;
}


/**
 * \param refl: Reflection
 * \param kpred: The wavenumber at which the reflection should be predicted
 *
 * Sets the wavenumber at which the reflection should be predicted.
 * Used by predict_to_res() and update_predictions()
 **/
void set_kpred(Reflection *refl, double kpred)
{
	refl->data.kpred = kpred;
}


/**
 * \param refl: Reflection
 * \param exerr: The excitation error for the reflection
 *
 **/
void set_exerr(Reflection *refl, double exerr)
{
	refl->data.exerr = exerr;
}


/**
 * \param refl: Reflection
 * \param p: The partiality for the reflection.
 *
 * Set the partiality for the reflection.  See set_lorentz().
 **/
void set_partiality(Reflection *refl, double p)
{
	refl->data.p = p;
}

/**
 * \param refl: Reflection
 * \param L: The Lorentz factor for the reflection.
 *
 * Set the Lorentz factor for the reflection.  To "scale up" a partial
 * reflection, divide by this multiplied by the partiality.
 **/
void set_lorentz(Reflection *refl, double L)
{
	refl->data.L = L;
}


/**
 * \param refl: Reflection
 * \param intensity: The intensity for the reflection.
 *
 * Set the intensity for the reflection.
 **/
void set_intensity(Reflection *refl, double intensity)
{
	refl->data.intensity = intensity;
}


/**
 * \param refl: Reflection
 * \param red: New redundancy for the reflection
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
 * \param refl: Reflection
 * \param esd: New standard error for this reflection's intensity measurement
 *
 **/
void set_esd_intensity(Reflection *refl, double esd)
{
	refl->data.esd_i = esd;
}


/**
 * \param refl: Reflection
 * \param phase: New phase for the reflection
 *
 **/
void set_phase(Reflection *refl, double phase)
{
	refl->data.phase = phase;
	refl->data.have_phase = 1;
}


/**
 * \param refl: Reflection
 * \param peak: New peak height for the reflection
 *
 **/
void set_peak(Reflection *refl, double peak)
{
	refl->data.peak = peak;
}


/**
 * \param refl: Reflection
 * \param mean_bg: New peak height for the reflection
 *
 **/
void set_mean_bg(Reflection *refl, double mean_bg)
{
	refl->data.mean_bg = mean_bg;
}


/**
 * \param refl: Reflection
 * \param hs: The 'h' index of the reflection
 * \param ks: The 'k' index of the reflection
 * \param ls: The 'l' index of the reflection
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
 * \param refl: A \ref Reflection
 * \param temp: New temporary value for the reflection
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
 * \param refl: A \ref Reflection
 * \param temp: New temporary value for the reflection
 *
 * The temporary values can be used according to the needs of the calling
 * program.
 *
 **/
void set_temp2(Reflection *refl, double temp)
{
	refl->data.temp2 = temp;
}


/**
 * \param refl: A \ref Reflection
 * \param flag: New flag value
 *
 * \param flag is an integer value which can be used according to the needs of the
 * calling program.
 *
 **/
void set_flag(Reflection *refl, int flag)
{
	refl->data.flag = flag;
}


/**
 * \param refl: Reflection
 * \param contribs: Pointer to the contribution list
 *
 * Note that the pointer will be stored, not the contents of the structure.
 *
 **/
void set_contributions(Reflection *refl,
                       struct reflection_contributions *contribs)
{
	refl->data.contribs = contribs;
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


static void add_refl_to_list_real(RefList *list,
                                  Reflection *new,
                                  signed int h,
                                  signed int k,
                                  signed int l)
{
	Reflection *f;

	assert(!new->in_list);

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

	new->in_list = 1;
}


/**
 * \param list: A %RefList
 * \param h: The 'h' index of the reflection
 * \param k: The 'k' index of the reflection
 * \param l: The 'l' index of the reflection
 *
 * Adds a new reflection to \p list.  Note that the implementation allows there to
 * be multiple reflections with the same indices in the list, so this function
 * should succeed even if the given indices already feature in the list.
 *
 * \returns The newly created reflection, or NULL on failure.
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

	add_refl_to_list_real(list, new, h, k, l);

	return new;
}


/**
 * \param refl: Reflection
 * \param list: A %RefList
 *
 * Adds \p refl to \p list.
 *
 **/
void add_refl_to_list(Reflection *refl, RefList *list)
{
	signed int h, k, l;

	get_indices(refl, &h, &k, &l);

	add_refl_to_list_real(list, refl, h, k, l);
}


/********************************* Iteration **********************************/

struct _reflistiterator {

	int stack_size;
	int stack_ptr;
	Reflection **stack;
	const Reflection **stack_const;
	int is_const;
};


/**
 * \param list: A %RefList to iterate over
 * \param piter: Address at which to store a %RefListIterator
 *
 * This function sets up the state required for iteration over the entire list,
 * and then returns the first reflection in the list.  An iterator object will
 * be created and its address stored at the location given in piter.
 *
 * \returns the first reflection in the list.
 *
 **/
Reflection *first_refl(RefList *list, RefListIterator **piter)
{
	Reflection *refl;
	RefListIterator *iter;

	iter = cfmalloc(sizeof(struct _reflistiterator));
	iter->stack_size = 32;
	iter->stack = cfmalloc(iter->stack_size*sizeof(Reflection *));
	iter->stack_ptr = 0;
	iter->is_const = 0;
	*piter = iter;

	if ( list == NULL ) return NULL;

	refl = list->head;

	do {

		if ( refl != NULL ) {
			iter->stack[iter->stack_ptr++] = refl;
			if ( iter->stack_ptr == iter->stack_size ) {
				iter->stack_size += 32;
				iter->stack = cfrealloc(iter->stack,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;
		}

		if ( iter->stack_ptr == 0 ) {
			cffree(iter->stack);
			cffree(iter);
			return NULL;
		}

		refl = iter->stack[--iter->stack_ptr];

		return refl;

	} while ( 1 );
}


/**
 * \param list: A %RefList to iterate over
 * \param piter: Address at which to store a %RefListIterator
 *
 * As first_refl(), except returns a const %Reflection.
 * Use this when you don't need to modify any of the reflections.
 *
 * \returns the first reflection in the list.
 *
 **/
const Reflection *first_refl_const(const RefList *list, RefListIterator **piter)
{
	const Reflection *refl;
	RefListIterator *iter;

	iter = cfmalloc(sizeof(struct _reflistiterator));
	iter->stack_size = 32;
	iter->stack_const = cfmalloc(iter->stack_size*sizeof(Reflection *));
	iter->stack_ptr = 0;
	iter->is_const = 1;
	*piter = iter;

	if ( list == NULL ) return NULL;

	refl = list->head;

	do {

		if ( refl != NULL ) {
			iter->stack_const[iter->stack_ptr++] = refl;
			if ( iter->stack_ptr == iter->stack_size ) {
				iter->stack_size += 32;
				iter->stack_const = cfrealloc(iter->stack_const,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;
		}

		if ( iter->stack_ptr == 0 ) {
			cffree(iter->stack_const);
			cffree(iter);
			return NULL;
		}

		refl = iter->stack_const[--iter->stack_ptr];

		return refl;

	} while ( 1 );
}


/**
 * \param refl: A reflection
 * \param iter: A %RefListIterator
 *
 * This function looks up the next reflection in the list that was given earlier
 * to first_refl().
 *
 * \returns the next reflection in the list, or NULL if no more.
 *
 **/
Reflection *next_refl(Reflection *refl, RefListIterator *iter)
{
	assert(!iter->is_const);

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
				iter->stack = cfrealloc(iter->stack,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;

		}
		if ( iter->stack_ptr == 0 ) {
			cffree(iter->stack);
			cffree(iter);
			return NULL;
		}

		return iter->stack[--iter->stack_ptr];

	} while ( 1 );
}


/**
 * \param refl: A reflection
 * \param iter: A %RefListIterator
 *
 * As next_refl(), except returns a const %Reflection.
 * Use this when you don't need to modify any of the reflections.
 *
 * \returns the next reflection in the list, or NULL if no more.
 *
 **/
const Reflection *next_refl_const(const Reflection *refl, RefListIterator *iter)
{
	assert(iter->is_const);

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

			iter->stack_const[iter->stack_ptr++] = refl;
			if ( iter->stack_ptr == iter->stack_size ) {
				iter->stack_size += 32;
				iter->stack_const = cfrealloc(iter->stack_const,
				         iter->stack_size*sizeof(Reflection *));
			}
			refl = refl->child[0];
			continue;

		}
		if ( iter->stack_ptr == 0 ) {
			free_reflistiterator(iter);
			return NULL;
		}

		return iter->stack_const[--iter->stack_ptr];

	} while ( 1 );
}


void free_reflistiterator(RefListIterator *iter)
{
	if ( iter != NULL ) {
		if ( iter->is_const ) {
			cffree(iter->stack_const);
		} else {
			cffree(iter->stack);
		}
		cffree(iter);
	}
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
 * \param list: A %RefList
 *
 * \returns the number of reflections in \p list.
 *
 **/
int num_reflections(RefList *list)
{
	return recursive_count(list->head);
}


/**
 * \param list: A %RefList
 *
 * If the depth of the tree is more than about 20, access to the list will be
 * slow.  This should never happen.
 *
 * \returns the depth of the RB-tree used internally to represent \p list.
 *
 **/
int tree_depth(RefList *list)
{
	return recursive_depth(list->head);
}


/**
 * \param refl: Reflection
 *
 * Acquires a lock on the reflection.
 */
void lock_reflection(Reflection *refl)
{
	pthread_mutex_lock(&refl->lock);
}


/**
 * \param refl: Reflection
 *
 * Releases a lock on the reflection.
 */
void unlock_reflection(Reflection *refl)
{
	pthread_mutex_unlock(&refl->lock);
}


static void reflist_set_notes(RefList *reflist, const char *notes)
{
	cffree(reflist->notes);  /* free(NULL) is OK */
	reflist->notes = cfstrdup(notes);
}


/**
 * \param reflist: Reflection list
 *
 * \returns the notes field for \p reflist, or NULL if there are no notes.
 * See reflist_add_notes() for more details.
 */
const char *reflist_get_notes(RefList *reflist)
{
	return reflist->notes;
}


/**
 * \param reflist: Reflection list
 * \param notes_add: Notes to add
 *
 * Appends the string \p notes_add to the notes field for \p reflist.  The notes
 * will be stored in the reflection list file by, e.g.,
 * write_reflist(), and are meant to be for humans to read.
 * Possible uses include making a record of the command line arguments used to
 * create the reflection list.
 */
void reflist_add_notes(RefList *reflist, const char *notes_add)
{
	size_t len;
	char *nnotes;

	if ( reflist->notes == NULL ) {
		reflist_set_notes(reflist, notes_add);
		return;
	}

	len = strlen(notes_add) + strlen(reflist->notes) + 2;
	nnotes = cfmalloc(len);
	if ( nnotes == NULL ) {
		ERROR("Failed to add notes to crystal.\n");
		return;
	}

	strcpy(nnotes, reflist->notes);
	strcat(nnotes, "\n");
	strcat(nnotes, notes_add);
	cffree(reflist->notes);
	reflist->notes = nnotes;
}
