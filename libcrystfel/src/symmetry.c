/*
 * symmetry.c
 *
 * Symmetry
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "symmetry.h"
#include "utils.h"
#include "integer_matrix.h"
#include "symop-parse.h"
#define YYSTYPE SYMOPSTYPE
#include "symop-lex.h"


/** \file symmetry.h */

struct _symoplist
{
	IntegerMatrix **ops;
	int n_ops;
	int max_ops;
	char *name;
	int num_equivs;
};


struct _symopmask
{
	const SymOpList *list;
	int *mask;
};



static void alloc_ops(SymOpList *ops)
{
	ops->ops = cfrealloc(ops->ops, ops->max_ops*sizeof(IntegerMatrix *));
}


/**
 * new_symopmask:
 * \param list A \ref SymOpList
 *
 * \returns A new \ref SymOpMask, which you can use when filtering out special
 * reflections.
 **/
SymOpMask *new_symopmask(const SymOpList *list)
{
	SymOpMask *m;
	int i;

	m = cfmalloc(sizeof(struct _symopmask));
	if ( m == NULL ) return NULL;

	m->list = list;
	m->mask = cfmalloc(sizeof(int)*list->n_ops);
	if ( m->mask == NULL ) {
		cffree(m);
		return NULL;
	}

	for ( i=0; i<list->n_ops; i++ ) {
		m->mask[i] = 1;
	}

	return m;
}


/* Creates a new SymOpList */
static SymOpList *new_symoplist()
{
	SymOpList *new;
	new = cfmalloc(sizeof(SymOpList));
	if ( new == NULL ) return NULL;
	new->max_ops = 16;
	new->n_ops = 0;
	new->ops = NULL;
	new->name = NULL;
	new->num_equivs = 1;
	alloc_ops(new);
	return new;
}


/**
 * \param ops A \ref SymOpList to free
 *
 * Frees a \ref SymOpList and all associated resources.
 **/
void free_symoplist(SymOpList *ops)
{
	int i;

	if ( ops == NULL ) return;
	for ( i=0; i<ops->n_ops; i++ ) {
		intmat_free(ops->ops[i]);
	}
	if ( ops->ops != NULL ) cffree(ops->ops);
	if ( ops->name != NULL ) cffree(ops->name);
	cffree(ops);
}

/**
 * \param m A \ref SymOpMask to free
 *
 * Frees a \ref SymOpMask and all associated resources.
 **/
void free_symopmask(SymOpMask *m)
{
	if ( m == NULL ) return;
	cffree(m->mask);
	cffree(m);
}


/* This returns the number of operations in "ops".  This might be different
 * to num_equivs() if the point group is being constructed. */
static int num_ops(const SymOpList *ops)
{
	return ops->n_ops;
}


/**
 * \param ops A \ref SymOpList
 * \param m An \ref IntegerMatrix
 *
 * Adds \p m to \p ops.
 **/
void add_symop(SymOpList *ops, IntegerMatrix *m)
{
	if ( ops->n_ops == ops->max_ops ) {
		ops->max_ops += 16;
		alloc_ops(ops);
	}

	ops->ops[ops->n_ops++] = m;
}


/* Add a operation to a SymOpList, starting from v(..) */
static void add_symop_v(SymOpList *ops,
                        signed int *h, signed int *k, signed int *l)
{
	IntegerMatrix *m;
	int i;

	m = intmat_new(3, 3);
	assert(m != NULL);

	for ( i=0; i<3; i++ ) intmat_set(m, i, 0, h[i]);
	for ( i=0; i<3; i++ ) intmat_set(m, i, 1, k[i]);
	for ( i=0; i<3; i++ ) intmat_set(m, i, 2, l[i]);

	cffree(h);
	cffree(k);
	cffree(l);

	add_symop(ops, m);
}


/**
 * \param ops A \ref SymOpList
 * \param m A \ref SymOpMask
 * \param idx Index of the operation to get
 *
 * This function returns a pointer to an integer matrix specifying a symmetry
 * operation contained in the symmetry operator list, and identified by the
 * specified index.
 *
 * The returned IntegerMatrix is owned by the SymOpList and must not be
 * freed separately.  It is valid as long as 'ops' exists.
 **/
IntegerMatrix *get_symop(const SymOpList *ops, const SymOpMask *m, int idx)
{
	const int n = num_ops(ops);

	if ( m != NULL ) {

		int i, c;

		c = 0;
		for ( i=0; i<n; i++ ) {

			if ( (c == idx) && m->mask[i] ) {
				return ops->ops[i];
			}

			if ( m->mask[i] ) {
				c++;
			}

		}

		ERROR("Index %i out of range for point group '%s'\n",
			      idx, symmetry_name(ops));

		return NULL;

	}

	if ( idx >= n ) {

		ERROR("Index %i out of range for point group '%s'\n", idx,
		      symmetry_name(ops));

		return NULL;

	}

	return ops->ops[idx];
}

static signed int *v(signed int h, signed int k, signed int i, signed int l)
{
	signed int *vec = cfmalloc(3*sizeof(signed int));
	if ( vec == NULL ) return NULL;
	/* Convert back to 3-index form now */
	vec[0] = h-i;  vec[1] = k-i;  vec[2] = l;
	return vec;
}


/**
 * \param ops A \ref SymOpList
 * \param m A \ref SymOpMask, which has been shown to \ref special_position
 *
 * \returns The number of equivalent reflections for a general reflection
 * in point group "ops", which were not flagged by your call to
 * \ref special_position.
 **/
int num_equivs(const SymOpList *ops, const SymOpMask *m)
{
	int n = num_ops(ops);
	int i;
	int c;

	if ( m == NULL ) return n;

	c = 0;
	for ( i=0; i<n; i++ ) {
		if ( m->mask[i] ) c++;
	}

	return c;
}


static void add_identity(SymOpList *s)
{
	int i, ni;
	int found;

	found = 0;
	ni = num_ops(s);
	for ( i=0; i<ni; i++ ) {
		if ( intmat_is_identity(s->ops[i]) ) {
			found = 1;
			break;
		}
	}
	if ( !found ) {
		add_symop_v(s, v(1,0,0,0), v(0,1,0,0), v(0,0,0,1));  /* I */
	}
}


/* Fill in the other operations for a point group starting from its
 * generators */
static void expand_ops(SymOpList *s)
{
	int added;

	add_identity(s);

	do {

		int i, ni;

		added = 0;

		ni = num_ops(s);
		for ( i=0; i<ni; i++ ) {

			int j;
			IntegerMatrix *opi = s->ops[i];

			/* Apply op 'i' to all the current ops in the list */
			for ( j=0; j<ni; j++ ) {

				IntegerMatrix *opj = s->ops[j];
				IntegerMatrix *m;
				int k, nk;
				int found;

				m = intmat_times_intmat(opi, opj);
				assert(m != NULL);

				nk = num_ops(s);
				found = 0;
				for ( k=0; k<nk; k++ ) {
					if ( intmat_equals(m, s->ops[k]) ) {
						found = 1;
						intmat_free(m);
						break;
					}
				}

				if ( !found ) {
					add_symop(s, m);
					added++;
				}

			}

		}

	} while ( added );
}


/* Transform all the operations in a SymOpList by a given matrix.
 * The matrix must have a determinant of +/- 1 (otherwise its inverse would
 * not also be an integer matrix). */
static void transform_ops(SymOpList *s, IntegerMatrix *P)
{
	int n, i;
	IntegerMatrix *Pi;
	signed int det;

	det = intmat_det(P);
	if ( det == -1 ) {
		ERROR("WARNING: mirrored SymOpList.\n");
	} else if ( det != 1 ) {
		ERROR("Invalid transformation for SymOpList.\n");
		return;
	}

	Pi = intmat_inverse(P);
	if ( Pi == NULL ) {
		ERROR("Failed to invert matrix.\n");
		return;
	}

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {

		IntegerMatrix *r, *f;

		r = intmat_times_intmat(P, s->ops[i]);
		if ( r == NULL ) {
			ERROR("Matrix multiplication failed.\n");
			return;
		}

		f = intmat_times_intmat(r, Pi);
		if ( f == NULL ) {
			ERROR("Matrix multiplication failed.\n");
			return;
		}
		intmat_free(r);

		intmat_free(s->ops[i]);
		s->ops[i] = intmat_copy(f);
		intmat_free(f);

	}

	intmat_free(Pi);
}


/********************************* Triclinic **********************************/

static SymOpList *make_1bar()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1));  /* -I */
	new->name = cfstrdup("-1");
	expand_ops(new);
	return new;
}


static SymOpList *make_1()
{
	SymOpList *new = new_symoplist();
	new->name = cfstrdup("1");
	expand_ops(new);
	return new;
}


/********************************* Monoclinic *********************************/

static SymOpList *make_2m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1));  /* m -| l */
	new->name = cfstrdup("2/m");
	expand_ops(new);
	return new;
}


static SymOpList *make_2()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	new->name = cfstrdup("2");
	expand_ops(new);
	return new;
}


static SymOpList *make_m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1));  /* m -| l */
	new->name = cfstrdup("m");
	expand_ops(new);
	return new;
}


/******************************** Orthorhombic ********************************/

static SymOpList *make_mmm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	add_symop_v(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m -| k */
	new->name = cfstrdup("mmm");
	expand_ops(new);
	return new;
}


static SymOpList *make_222()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	new->name = cfstrdup("222");
	expand_ops(new);
	return new;
}


static SymOpList *make_mm2()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m -| k */
	new->name = cfstrdup("mm2");
	expand_ops(new);
	return new;
}


/********************************* Tetragonal *********************************/

static SymOpList *make_4m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	add_symop_v(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* m -| l */
	new->name = cfstrdup("4/m");
	expand_ops(new);
	return new;
}


static SymOpList *make_4()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	new->name = cfstrdup("4");
	expand_ops(new);
	return new;
}


static SymOpList *make_4mm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,1)); /* m -| l */
	new->name = cfstrdup("4mm");
	expand_ops(new);
	return new;
}


static SymOpList *make_422()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1));  /* 4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	new->name = cfstrdup("422");
	expand_ops(new);
	return new;
}


static SymOpList *make_4bar()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* -4 // l */
	new->name = cfstrdup("-4");
	expand_ops(new);
	return new;
}


static SymOpList *make_4bar2m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* -4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	new->name = cfstrdup("-42m");
	expand_ops(new);
	return new;
}


static SymOpList *make_4barm2()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* -4 // l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1)); /* 2 // h+k */
	new->name = cfstrdup("-4m2");
	expand_ops(new);
	return new;
}


static SymOpList *make_4mmm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,1)); /* m -| k */
	add_symop_v(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* m -| l */
	new->name = cfstrdup("4/mmm");
	expand_ops(new);
	return new;
}


/************************** Trigonal (Rhombohedral) ***************************/

static SymOpList *make_3_R()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0)); /* 3 // h+k+l */
	new->name = cfstrdup("3_R");
	expand_ops(new);
	return new;
}


static SymOpList *make_3bar_R()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0)); /* -3 // h+k+l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	new->name = cfstrdup("-3_R");
	expand_ops(new);
	return new;
}


static SymOpList *make_32_R()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0)); /* 3 // h+k+l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* 2 -| 3 */
	new->name = cfstrdup("32_R");
	expand_ops(new);
	return new;
}


static SymOpList *make_3m_R()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0)); /* 3 // h+k+l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1)); /* m */
	new->name = cfstrdup("3m_R");
	expand_ops(new);
	return new;
}


static SymOpList *make_3barm_R()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0)); /* -3 // h+k+l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1));    /* m */
	new->name = cfstrdup("-3m_R");
	expand_ops(new);
	return new;
}


/*************************** Trigonal (Hexagonal) *****************************/

static SymOpList *make_3_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1)); /* 3 // l */
	new->name = cfstrdup("3_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_3bar_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1));    /* 3 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	new->name = cfstrdup("-3_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_321_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1));  /* 3 // l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1)); /* 2 // h */
	new->name = cfstrdup("321_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_312_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1));    /* 3 // l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* 2 // h+k */
	new->name = cfstrdup("312_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_3m1_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1)); /* 3 // l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1)); /* m -| i */
	new->name = cfstrdup("3m1_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_31m_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1)); /* 3 // l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1)); /* m -| (k+i) */
	new->name = cfstrdup("31m_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_3barm1_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1));    /* 3 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1));   /* 2 // h */
	new->name = cfstrdup("-3m1_H");
	expand_ops(new);
	return new;
}


static SymOpList *make_3bar1m_H()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1));    /* 3 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* 2 // h+k */
	new->name = cfstrdup("-31m_H");
	expand_ops(new);
	return new;
}


/********************************** Hexgonal **********************************/

static SymOpList *make_6()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1)); /* 6 // l */
	new->name = cfstrdup("6");
	expand_ops(new);
	return new;
}


static SymOpList *make_6bar()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1)); /* -6 // l */
	new->name = cfstrdup("-6");
	expand_ops(new);
	return new;
}


static SymOpList *make_6m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1)); /* 6 // l */
	add_symop_v(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1));  /* m -| l */
	new->name = cfstrdup("6/m");
	expand_ops(new);
	return new;
}


static SymOpList *make_622()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1)); /* 6 // l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1));   /* 2 // h */
	new->name = cfstrdup("622");
	expand_ops(new);
	return new;
}


static SymOpList *make_6mm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1)); /* 6 // l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1)); /* m -| i */
	new->name = cfstrdup("6mm");
	expand_ops(new);
	return new;
}


static SymOpList *make_6barm2()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1)); /* -6 // l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1)); /* m -| i */
	new->name = cfstrdup("-6m2");
	expand_ops(new);
	return new;
}


static SymOpList *make_6bar2m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1)); /* -6 // l */
	add_symop_v(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1));  /* m -| (k+i) */
	new->name = cfstrdup("-62m");
	expand_ops(new);
	return new;
}


static SymOpList *make_6mmm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1)); /* -6 // l */
	add_symop_v(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1)); /* m -| i */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	new->name = cfstrdup("6/mmm");
	expand_ops(new);
	return new;
}


/************************************ Cubic ***********************************/

static SymOpList *make_23()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	add_symop_v(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0)); /* 3 // h+k+l */
	new->name = cfstrdup("23");
	expand_ops(new);
	return new;
}


static SymOpList *make_m3bar()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	add_symop_v(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0)); /* 3 // h+k+l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	new->name = cfstrdup("m-3");
	expand_ops(new);
	return new;
}


static SymOpList *make_432()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1));/* 2 // k */
	add_symop_v(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0));  /* 3 // h+k+l */
	new->name = cfstrdup("432");
	expand_ops(new);
	return new;
}


static SymOpList *make_4bar3m()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1)); /* -4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 // k */
	add_symop_v(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0));   /* 3 // h+k+l */
	new->name = cfstrdup("-43m");
	expand_ops(new);
	return new;
}


static SymOpList *make_m3barm()
{
	SymOpList *new = new_symoplist();
	add_symop_v(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 // l */
	add_symop_v(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1));/* 2 // k */
	add_symop_v(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0));  /* 3 // h+k+l */
	add_symop_v(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1)); /* -I */
	new->name = cfstrdup("m-3m");
	expand_ops(new);
	return new;
}


static SymOpList *getpg_uac(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "-1") == 0 ) return make_1bar();
	if ( strcmp(sym, "1") == 0 ) return make_1();

	/* Monoclinic */
	if ( strcmp(sym, "2/m") == 0 ) return make_2m();
	if ( strcmp(sym, "2") == 0 ) return make_2();
	if ( strcmp(sym, "m") == 0 ) return make_m();

	/* Orthorhombic */
	if ( strcmp(sym, "mmm") == 0 ) return make_mmm();
	if ( strcmp(sym, "222") == 0 ) return make_222();
	if ( strcmp(sym, "mm2") == 0 ) return make_mm2();

	/* Tetragonal */
	if ( strcmp(sym, "4/m") == 0 ) return make_4m();
	if ( strcmp(sym, "4") == 0 ) return make_4();
	if ( strcmp(sym, "-4") == 0 ) return make_4bar();
	if ( strcmp(sym, "4/mmm") == 0 ) return make_4mmm();
	if ( strcmp(sym, "422") == 0 ) return make_422();
	if ( strcmp(sym, "-42m") == 0 ) return make_4bar2m();
	if ( strcmp(sym, "-4m2") == 0 ) return make_4barm2();
	if ( strcmp(sym, "4mm") == 0 ) return make_4mm();

	/* Trigonal (rhombohedral) */
	if ( strcmp(sym, "3_R") == 0 ) return make_3_R();
	if ( strcmp(sym, "-3_R") == 0 ) return make_3bar_R();
	if ( strcmp(sym, "32_R") == 0 ) return make_32_R();
	if ( strcmp(sym, "3m_R") == 0 ) return make_3m_R();
	if ( strcmp(sym, "-3m_R") == 0 ) return make_3barm_R();

	/* Trigonal (hexagonal) */
	if ( strcmp(sym, "3_H") == 0 ) return make_3_H();
	if ( strcmp(sym, "-3_H") == 0 ) return make_3bar_H();
	if ( strcmp(sym, "321_H") == 0 ) return make_321_H();
	if ( strcmp(sym, "312_H") == 0 ) return make_312_H();
	if ( strcmp(sym, "3m1_H") == 0 ) return make_3m1_H();
	if ( strcmp(sym, "31m_H") == 0 ) return make_31m_H();
	if ( strcmp(sym, "-3m1_H") == 0 ) return make_3barm1_H();
	if ( strcmp(sym, "-31m_H") == 0 ) return make_3bar1m_H();

	/* Hexagonal */
	if ( strcmp(sym, "6/m") == 0 ) return make_6m();
	if ( strcmp(sym, "6") == 0 ) return make_6();
	if ( strcmp(sym, "-6") == 0 ) return make_6bar();
	if ( strcmp(sym, "6/mmm") == 0 ) return make_6mmm();
	if ( strcmp(sym, "622") == 0 ) return make_622();
	if ( strcmp(sym, "-62m") == 0 ) return make_6bar2m();
	if ( strcmp(sym, "-6m2") == 0 ) return make_6barm2();
	if ( strcmp(sym, "6mm") == 0 ) return make_6mm();

	/* Cubic */
	if ( strcmp(sym, "23") == 0 ) return make_23();
	if ( strcmp(sym, "m-3") == 0 ) return make_m3bar();
	if ( strcmp(sym, "432") == 0 ) return make_432();
	if ( strcmp(sym, "-43m") == 0 ) return make_4bar3m();
	if ( strcmp(sym, "m-3m") == 0 ) return make_m3barm();

	ERROR("Unknown point group '%s'\n", sym);
	return NULL;
}


static int char_count(const char *a, char b)
{
	size_t i;
	int n;

	i = 0;  n = 0;
	do {
		if ( a[i] == b ) n++;
		if ( a[i] == '\0' ) return n;
		i++;
	} while ( 1 );
}


static SymOpList *getpg_arbitrary_ua(const char *sym, size_t s)
{
	char ua;
	char *pg_type;
	SymOpList *pg;
	IntegerMatrix *t;
	char *new_name;

	if ( strncmp(sym+s, "ua", 2) == 0 ) {
		ua = sym[s+2];
	} else {
		ERROR("Unrecognised point group '%s'\n", sym);
		return NULL;
	}

	pg_type = cfstrndup(sym, s-1);
	if ( pg_type == NULL ) {
		ERROR("Couldn't allocate string.\n");
		return NULL;
	}

	pg = getpg_uac(pg_type);
	if ( pg == NULL ) {
		ERROR("Unrecognised point group type '%s'\n",
		      pg_type);
		return NULL;
	}
	cffree(pg_type);

	t = intmat_new(3, 3);
	if ( t == NULL ) return NULL;

	switch ( ua ) {

		case 'a' :
		intmat_set(t, 0, 2, 1);
		intmat_set(t, 1, 0, 1);
		intmat_set(t, 2, 1, 1);
		break;

		case 'b' :
		intmat_set(t, 0, 1, 1);
		intmat_set(t, 1, 2, 1);
		intmat_set(t, 2, 0, 1);

		break;

		case 'c' :
		intmat_set(t, 0, 0, 1);
		intmat_set(t, 1, 1, 1);
		intmat_set(t, 2, 2, 1);
		break;

		default :
		ERROR("Bad unique axis '%c'\n", ua);
		free_symoplist(pg);
		return NULL;

	}

	transform_ops(pg, t);
	intmat_free(t);

	new_name = cfmalloc(64);
	if ( new_name == NULL ) {
		ERROR("Couldn't allocate space for PG name\n");
		return NULL;
	}

	snprintf(new_name, 64, "%s_ua%c", pg->name, ua);
	cffree(pg->name);
	pg->name = new_name;

	return pg;
}


/**
 * \param sym A string representation of a point group
 *
 * This function parses \p sym and returns the corresponding \ref SymOpList.
 * In the string representation of the point group, use a preceding minus sign
 * for any character which would have a "bar".  Trigonal groups must be suffixed
 * with either "_H" or "_R" for a hexagonal or rhombohedral lattice
 * respectively.
 *
 * Examples: -1 1 2/m 2 m mmm 222 mm2 4/m 4 -4 4/mmm 422 -42m -4m2 4mm
 * 3_R -3_R 32_R 3m_R -3m_R 3_H -3_H 321_H 312_H 3m1_H 31m_H -3m1_H -31m_H
 * 6/m 6 -6 6/mmm 622 -62m -6m2 6mm 23 m-3 432 -43m m-3m.
 **/
SymOpList *get_pointgroup(const char *sym)
{
	int n_underscore;

	n_underscore = char_count(sym, '_');

	/* No spaces nor underscores -> old system */
	if ( n_underscore == 0 ) return getpg_uac(sym);

	/* No spaces and 1 underscore -> old system + lattice or UA */
	if ( n_underscore == 1 ) {

		const char *s;

		s = strchr(sym, '_');
		assert(s != NULL);
		s++;

		/* Old system with H/R lattice? */
		if ( (s[0] == 'H') || (s[0] == 'R') ) {
			return getpg_uac(sym);
		}

		/* Old system with unique axis */
		return getpg_arbitrary_ua(sym, s-sym);

	}

	ERROR("Unrecognised point group '%s'\n", sym);
	return NULL;
}


void pointgroup_warning(const char *sym)
{
	if ( (strcmp(sym, "m") == 0)
	  || (strcmp(sym, "2/m") == 0)
	  || (strcmp(sym, "2") == 0) )
	{
		ERROR("WARNING: You have specified a monoclinic point group "
		      "without a unique axis.  The default unique axis is 'c'. "
		      "If you want unique axis b, append '_uab' to your point "
		      "group.\n");
	}

}


static void do_op(const IntegerMatrix *op,
                  signed int h, signed int k, signed int l,
                  signed int *he, signed int *ke, signed int *le)
{
	signed int vec[3];
	signed int *ans;

	vec[0] = h;  vec[1] = k;  vec[2] = l;

	ans = transform_indices(op, vec);
	assert(ans != NULL);

	*he = ans[0];  *ke = ans[1];  *le = ans[2];
	cffree(ans);
}


/**
 * \param ops A \ref SymOpList
 * \param m A \ref SymOpMask, which has been shown to \ref special_position
 * \param idx Index of the operation to use
 * \param h index of reflection
 * \param k index of reflection
 * \param l index of reflection
 * \param he location to store h index of equivalent reflection
 * \param ke location to store k index of equivalent reflection
 * \param le location to store l index of equivalent reflection
 *
 * This function applies the \p idx-th symmetry operation from \p ops to the
 * reflection \p h, \p k, \p l, and stores the result at \p he, \p ke and \p le.
 *
 * Call this function multiple times with idx=0 .. num_equivs(ops, m) to get all
 * of the equivalent reflections in turn.
 *
 * If you don't mind that the same equivalent might appear twice, simply let
 * \p m = NULL.  Otherwise, call \ref new_symopmask and then
 * \ref special_position to set up a \ref SymOpMask appropriately.
 **/
void get_equiv(const SymOpList *ops, const SymOpMask *m, int idx,
               signed int h, signed int k, signed int l,
               signed int *he, signed int *ke, signed int *le)
{
	IntegerMatrix *op;
	op = get_symop(ops, m, idx);
	if ( op == NULL ) {
		fprintf(stderr, "Cannot proceed.\n");
		abort();
	}
	do_op(op, h, k, l, he, ke, le);
}


/**
 * \param ops A \ref SymOpList, usually corresponding to a point group
 * \param m A \ref SymOpMask created with \ref new_symopmask
 * \param h index of a reflection
 * \param k index of a reflection
 * \param l index of a reflection
 *
 * This function sets up \p m to contain information about which operations in
 * \p ops map the reflection \p h, \p k, \p l onto itself.
 **/
void special_position(const SymOpList *ops, SymOpMask *m,
                      signed int h, signed int k, signed int l)
{
	int i, n;
	signed int *htest;
	signed int *ktest;
	signed int *ltest;

	assert(m->list == ops);

	n = num_equivs(ops, NULL);
	htest = cfmalloc(n*sizeof(signed int));
	ktest = cfmalloc(n*sizeof(signed int));
	ltest = cfmalloc(n*sizeof(signed int));

	for ( i=0; i<n; i++ ) {

		signed int he, ke, le;
		int j;

		get_equiv(ops, NULL, i, h, k, l, &he, &ke, &le);

		m->mask[i] = 1;
		for ( j=0; j<i; j++ ) {
			if ( (he==htest[j]) && (ke==ktest[j])
			  && (le==ltest[j]) )
			{
				m->mask[i] = 0;
				break;  /* Only need to find one */
			}
		}

		htest[i] = he;
		ktest[i] = ke;
		ltest[i] = le;

	}

	cffree(htest);
	cffree(ktest);
	cffree(ltest);
}


static int any_negative(signed int h, signed int k, signed int l)
{
	if ( h < 0 ) return 1;
	if ( k < 0 ) return 1;
	if ( l < 0 ) return 1;
	return 0;
}


/**
 * \param h h index
 * \param k k index
 * \param l l index
 * \param ops A \ref SymOpList
 *
 * A reflection is centric if it is related by symmetry to its Friedel partner.
 *
 * \returns True if \p h \p k \p l is centric in \p ops.
 *
 **/
int is_centric(signed int h, signed int k, signed int l, const SymOpList *ops)
{
	signed int ha, ka, la;
	signed int hb, kb, lb;

	get_asymm(ops, h, k, l, &ha, &ka, &la);
	get_asymm(ops, -h, -k, -l, &hb, &kb, &lb);

	if ( ha != hb ) return 0;
	if ( ka != kb ) return 0;
	if ( la != lb ) return 0;

	return 1;
}


/**
 * \param ops A \ref SymOpList, usually corresponding to a point group
 * \param h index of a reflection
 * \param k index of a reflection
 * \param l index of a reflection
 * \param hp location for asymmetric index of reflection
 * \param kp location for asymmetric index of reflection
 * \param lp location for asymmetric index of reflection
 *
 * This function determines the asymmetric version of the reflection \p h, \p k, \p l
 * in symmetry group \p ops, and puts the result in \p hp, \p kp, \p lp.
 *
 * This is a relatively expensive operation because of its generality.
 * Therefore, if you know you'll need to make repeated use of the asymmetric
 * indices, consider creating a new \ref RefList indexed according to the asymmetric
 * indices themselves with \ref asymmetric_indices.  If you do that, you'll still
 * be able to get the original versions of the indices with
 * \ref get_symmetric_indices.
 *
 **/
void get_asymm(const SymOpList *ops,
               signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp)
{
	int nequiv;
	int p;
	signed int best_h, best_k, best_l;
	int have_negs;

	nequiv = num_equivs(ops, NULL);

	best_h = h;  best_k = k;  best_l = l;
	have_negs = any_negative(best_h, best_k, best_l);
	for ( p=0; p<nequiv; p++ ) {

		int will_have_negs;

		get_equiv(ops, NULL, p, h, k, l, hp, kp, lp);

		will_have_negs = any_negative(*hp, *kp, *lp);

		/* Don't lose "no negs" status */
		if ( !have_negs && will_have_negs ) continue;

		if ( have_negs && !will_have_negs ) {
			best_h = *hp;  best_k = *kp;  best_l = *lp;
			have_negs = 0;
			continue;
		}

		if ( *hp > best_h ) {
			best_h = *hp;  best_k = *kp;  best_l = *lp;
			have_negs = any_negative(best_h, best_k, best_l);
			continue;
		}
		if ( *hp < best_h ) continue;

		if ( *kp > best_k ) {
			best_h = *hp;  best_k = *kp;  best_l = *lp;
			have_negs = any_negative(best_h, best_k, best_l);
			continue;
		}
		if ( *kp < best_k ) continue;

		if ( *lp > best_l ) {
			best_h = *hp;  best_k = *kp;  best_l = *lp;
			have_negs = any_negative(best_h, best_k, best_l);
			continue;
		}

	}

	*hp = best_h;  *kp = best_k;  *lp = best_l;
}


/**
 * \param s A \ref SymOpList
 *
 * \returns Non-zero if \p s contains an inversion operation
 */
int is_centrosymmetric(const SymOpList *s)
{
	int i, n;

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {
		if ( intmat_is_inversion(s->ops[i]) ) return 1;
	}

	return 0;
}


/* Return true if a*b = ans */
static int check_mult(const IntegerMatrix *ans,
                      const IntegerMatrix *a, const IntegerMatrix *b)
{
	int val;
	IntegerMatrix *m;

	m = intmat_times_intmat(a, b);
	assert(m != NULL);

	val = intmat_equals(ans, m);
	intmat_free(m);

	return val;
}


/**
 * \param source A \ref SymOpList
 * \param target Another \ref SymOpList, which might be a subgroup of \p source.
 *
 * \returns Non-zero if every operation in \p target is also in \p source.
 **/
int is_subgroup(const SymOpList *source, const SymOpList *target)
{
	int n_src, n_tgt;
	int i;

	n_src = num_ops(source);
	n_tgt = num_ops(target);

	for ( i=0; i<n_tgt; i++ ) {

		int j;
		int found = 0;

		for ( j=0; j<n_src; j++ ) {

			if ( intmat_equals(target->ops[i], source->ops[j]) ) {
				found = 1;
				break;
			}

		}

		if ( !found ) return 0;

	}

	return 1;
}


/* Returns n, where m^n = I */
static int order(const IntegerMatrix *m)
{
	IntegerMatrix *a;
	int i;

	a = intmat_new(3, 3);
	assert(a != NULL);
	intmat_set(a, 0, 0, 1);
	intmat_set(a, 1, 1, 1);
	intmat_set(a, 2, 2, 1);

	i = 0;
	do {

		IntegerMatrix *anew;

		anew = intmat_times_intmat(m, a);
		assert(anew != NULL);
		intmat_free(a);
		a = anew;

		i++;

	} while ( !intmat_is_identity(a) );

	return i;
}


static SymOpList *flack_reorder(const SymOpList *source)
{
	SymOpList *src_reordered;
	SymOpMask *used;
	int i, n_src;

	src_reordered = new_symoplist();
	used = new_symopmask(source);

	n_src = num_ops(source);

	/* Find identity */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( intmat_is_identity(source->ops[i]) ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	/* Find binary options (order=2) of first kind (determinant positive) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( (order(source->ops[i]) == 2)
		  && (intmat_det(source->ops[i]) > 0) ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	/* Find other operations of first kind (determinant positive) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( intmat_det(source->ops[i]) > 0 ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	/* Find inversion */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( intmat_is_inversion(source->ops[i]) ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	/* Find binary options of second kind (determinant negative) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( (order(source->ops[i]) == 2)
		  && (intmat_det(source->ops[i]) < 0) ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	/* Find other operations of second kind (determinant negative) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( intmat_det(source->ops[i]) < 0 ) {
			add_symop(src_reordered, intmat_copy(source->ops[i]));
			used->mask[i] = 0;
		}
	}

	int n_left_over = 0;
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		n_left_over++;
	}
	if ( n_left_over != 0 ) {
		ERROR("%i operations left over after rearranging for"
		      " left coset decomposition.\n", n_left_over);
	}

	if ( num_ops(src_reordered) != num_ops(source) ) {
		ERROR("%i ops went to %i after rearranging.\n",
		      num_ops(src_reordered), num_ops(source));
	}

	free_symopmask(used);

	return src_reordered;
}


/**
 * \param source The "source" symmetry, a \ref SymOpList
 * \param target The "target" symmetry, a \ref SymOpList

 * Calculates twinning laws.  Returns a \ref SymOpList containing the twinning
 * operators, which are the symmetry operations which can be added to \p target
 * to generate \p source.  Only rotations are allowable - no mirrors nor
 * inversions.
 * To count the number of possibilities, use \ref num_equivs on the result.
 *
 * The algorithm used is "Algorithm A" from Flack (1987), Acta Cryst A43 p564.
 *
 * \returns a \ref SymOpList containing the twinning operators, or NULL if the
 * source symmetry cannot be generated from that target symmetry without using
 * mirror or inversion operations.
 */
SymOpList *get_ambiguities(const SymOpList *source, const SymOpList *target)
{
	int n_src, n_tgt;
	int i;
	SymOpList *twins;
	SymOpList *src_reordered;
	SymOpList *tgt_reordered;
	SymOpMask *used;
	char *name;
	int have_identity = 0;

	n_src = num_ops(source);
	n_tgt = num_ops(target);

	if ( !is_subgroup(source, target) ) {
		ERROR("'%s' is not a subgroup of '%s'\n",
		      symmetry_name(target), symmetry_name(source));
		return NULL;
	}

	if ( n_src % n_tgt != 0 ) {
		ERROR("Subgroup index would be fractional.\n");
		return NULL;
	}

	/* Reorder operations to prefer rotations in the output */
	src_reordered = flack_reorder(source);
	if ( src_reordered == NULL ) return NULL;

	/* Reorder the subgroup as well, but strictly speaking we only need
	 * the identity at the beginning */
	tgt_reordered = flack_reorder(target);
	if ( tgt_reordered == NULL ) return NULL;

	used = new_symopmask(src_reordered);
	for ( i=0; i<n_src; i++ ) {

		int j;
		if ( used->mask[i] == 0 ) continue;

		for ( j=1; j<n_tgt; j++ ) {

			int k;
			for ( k=i+1; k<n_src; k++ ) {
				if ( check_mult(src_reordered->ops[k],
				                src_reordered->ops[i],
				                tgt_reordered->ops[j]) )
				{
					used->mask[k] = 0;
				}
			}

		}

	}

	twins = new_symoplist();
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( intmat_det(src_reordered->ops[i]) < 0 ) {
			/* A mirror or inversion turned up in the list.
			 * That means that no pure rotational ambiguity can
			 * account for this subgroup relationship. */
			free_symoplist(twins);
			free_symopmask(used);
			free_symoplist(src_reordered);
			return NULL;
		}
		if ( !intmat_is_identity(src_reordered->ops[i]) ) {
			add_symop(twins, intmat_copy(src_reordered->ops[i]));
		} else {
			have_identity = 1;
		}

	}

	if ( !have_identity ) {
		ERROR("WARNING: Identity not found during left coset decomp\n");
	}

	free_symopmask(used);
	free_symoplist(src_reordered);
	free_symoplist(tgt_reordered);

	name = cfmalloc(64);
	snprintf(name, 63, "%s -> %s", symmetry_name(source),
	                               symmetry_name(target));
	twins->name = name;

	return twins;
}


/* Parse a single symmetry operation, e.g. 'h,-2k,(h+l)/3' */
RationalMatrix *parse_symmetry_operation(const char *s)
{
	YY_BUFFER_STATE b;
	RationalMatrix *m;
	int r;
	void *scanner;

	m = rtnl_mtx_new(3, 3);
	symoplex_init(&scanner);
	b = symop_scan_string(s, scanner);
	r = symopparse(scanner, m, NULL);
	symop_delete_buffer(b, scanner);
	symoplex_destroy(scanner);

	if ( r ) {
		ERROR("Failed to parse '%s'\n", s);
		rtnl_mtx_free(m);
		return NULL;
	}

	return m;
}


/**
 * \param s Textual representation of cell transformation
 *
 * Parses \p s, for example 'a,(b+a)/2,c', and returns the corresponding
 * \ref RationalMatrix.
 *
 * \returns a \ref RationalMatrix describing the transformation, or NULL on error.
 *
 */
RationalMatrix *parse_cell_transformation(const char *s)
{
	return parse_symmetry_operation(s);
}


/**
 * \param s Textual representation of a list of symmetry operations
 *
 * Parses \p s, for example 'h,k,l;k,h,-l', and returns the corresponding
 * \ref SymOpList
 *
 * \returns a \ref SymOpList, or NULL on error.
 *
 */
SymOpList *parse_symmetry_operations(const char *s)
{
	YY_BUFFER_STATE b;
	RationalMatrix *m;
	SymOpList *list;
	int r;
	void *scanner;

	m = rtnl_mtx_new(3, 3); /* Scratch space for parser */
	list = new_symoplist(); /* The result we want */

	symoplex_init(&scanner);
	b = symop_scan_string(s, scanner);
	r = symopparse(scanner, m, list);
	symop_delete_buffer(b, scanner);
	symoplex_destroy(scanner);
	rtnl_mtx_free(m);

	if ( r ) {
		ERROR("Failed to parse '%s'\n", s);
		free_symoplist(list);
		return NULL;
	}

	return list;
}


static void add_chars(char *t, const char *s, size_t max_len)
{
	size_t len;

	len = strlen(t) + strlen(s);
	if ( len > max_len ) {
		ERROR("get_matrix_name: String too long!\n");
		return;
	}

	strcat(t, s);
}


char *get_matrix_name(const IntegerMatrix *m, int col)
{
	char *text;
	const size_t max_len = 31;
	int i;
	int printed = 0;

	text = cfmalloc(max_len+1);
	text[0] = '\0';

	for ( i=0; i<3; i++ ) {

		signed int val;

		val = intmat_get(m, i, col);

		if ( val == 0 ) continue;

		if ( val < 0 ) {
			add_chars(text, "-", max_len);
		} else {
			if ( printed ) add_chars(text, "+", max_len);
		}

		if ( abs(val) > 1 ) {
			char num[11];
			snprintf(num, 11, "%i", abs(val));
			add_chars(text, num, max_len);
		}

		switch ( i )
		{
			case 0  : add_chars(text, "h", max_len); break;
			case 1  : add_chars(text, "k", max_len); break;
			case 2  : add_chars(text, "l", max_len); break;
			default : add_chars(text, "X", max_len); break;
		}

		printed = 1;

	}

	return text;
}


char *name_equiv(const IntegerMatrix *op)
{
	char *h, *k, *l;
	char *name;

	h = get_matrix_name(op, 0);
	k = get_matrix_name(op, 1);
	l = get_matrix_name(op, 2);
	name = cfmalloc(32);

	if ( strlen(h)+strlen(k)+strlen(l) == 3 ) {
		snprintf(name, 31, "%s%s%s", h, k, l);
	} else {
		snprintf(name, 31, "%s,%s,%s", h, k, l);
	}

	return name;
}


/**
 * \param s A \ref SymOpList
 *
 * Writes the name and a list of operations to stderr.
 */
void describe_symmetry(const SymOpList *s)
{
	int i, n;
	size_t max_len = 0;

	n = num_equivs(s, NULL);

	STATUS("%15s : ", symmetry_name(s));

	for ( i=0; i<n; i++ ) {
		size_t len;
		char *name = name_equiv(s->ops[i]);
		len = strlen(name);
		if ( len > max_len ) max_len = len;
		cffree(name);
	}
	if ( max_len < 8 ) max_len = 8;

	for ( i=0; i<n; i++ ) {

		char *name;
		size_t m, j;

		name = name_equiv(s->ops[i]);
		m = max_len - strlen(name) + 3;

		STATUS("%s", name);
		for ( j=0; j<m; j++ ) {
			STATUS(" ");
		}
		cffree(name);
		if ( (i!=0) && (i%8==0) ) STATUS("\n%15s   ", "");

	}
	STATUS("\n");
}


/**
 * \param ops A \ref SymOpList
 *
 * \returns A text description of \p ops.
 */
const char *symmetry_name(const SymOpList *ops)
{
	return ops->name;
}


/**
 * \param ops A \ref SymOpList
 * \param name New name for the \ref SymOpList
 *
 * Sets the text description of \p ops to \p name.  See \ref symmetry_name.
 * \p name will be copied, so you can safely free it after calling this function,
 * if that's otherwise appropriate.
 */
void set_symmetry_name(SymOpList *ops, const char *name)
{
	ops->name = cfstrdup(name);
}
