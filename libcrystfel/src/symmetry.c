/*
 * symmetry.c
 *
 * Symmetry
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "symmetry.h"
#include "utils.h"
#include "integer_matrix.h"


/**
 * SECTION:symmetry
 * @short_description: Point symmetry handling
 * @title: Symmetry
 * @section_id:
 * @see_also:
 * @include: "symmetry.h"
 * @Image:
 *
 * Routines to handle point symmetry.
 */


struct sym_op
{
	signed int *h;
	signed int *k;
	signed int *l;  /* Contributions to h, k and l from h, k and l */
	int order;
};


struct _symoplist
{
	struct sym_op *ops;
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
	ops->ops = realloc(ops->ops, ops->max_ops*sizeof(struct sym_op));
}


/**
 * new_symopmask:
 * @list: A %SymOpList
 *
 * Returns: a new %SymOpMask, which you can use when filtering out special
 * reflections.
 **/
SymOpMask *new_symopmask(const SymOpList *list)
{
	SymOpMask *m;
	int i;

	m = malloc(sizeof(struct _symopmask));
	if ( m == NULL ) return NULL;

	m->list = list;
	m->mask = malloc(sizeof(int)*list->n_ops);
	if ( m->mask == NULL ) {
		free(m);
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
	new = malloc(sizeof(SymOpList));
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
 * free_symoplist:
 * @ops: A %SymOpList to free
 *
 * Frees a %SymOpList and all associated resources.
 **/
void free_symoplist(SymOpList *ops)
{
	int i;

	if ( ops == NULL ) return;
	for ( i=0; i<ops->n_ops; i++ ) {
		free(ops->ops[i].h);
		free(ops->ops[i].k);
		free(ops->ops[i].l);
	}
	if ( ops->ops != NULL ) free(ops->ops);
	if ( ops->name != NULL ) free(ops->name);
	free(ops);
}

/**
 * free_symopmask:
 * @m: A %SymOpMask to free
 *
 * Frees a %SymOpMask and all associated resources.
 **/
void free_symopmask(SymOpMask *m)
{
	if ( m == NULL ) return;
	free(m->mask);
	free(m);
}


/* This returns the number of operations in "ops".  This might be different
 * to num_equivs() if the point group is being constructed. */
static int num_ops(const SymOpList *ops)
{
	return ops->n_ops;
}


/* Add a operation to a SymOpList */
static void add_symop(SymOpList *ops,
                      signed int *h, signed int *k, signed int *l,
                      int order)
{
	int n;

	if ( ops->n_ops == ops->max_ops ) {
		/* Pretty sure this never happens, but still... */
		ops->max_ops += 16;
		alloc_ops(ops);
	}

	n = ops->n_ops;
	ops->ops[n].h = h;
	ops->ops[n].k = k;
	ops->ops[n].l = l;
	ops->ops[n].order = order;
	ops->n_ops++;
}


/* Add a operation to a SymOpList */
static void add_copied_op(SymOpList *ops, struct sym_op *copyme)
{
	int n;
	signed int *h, *k, *l;

	if ( ops->n_ops == ops->max_ops ) {
		ops->max_ops += 16;
		alloc_ops(ops);
	}

	n = ops->n_ops;

	h = malloc(3*sizeof(signed int));
	k = malloc(3*sizeof(signed int));
	l = malloc(3*sizeof(signed int));

	memcpy(h, copyme->h, 3*sizeof(signed int));
	memcpy(k, copyme->k, 3*sizeof(signed int));
	memcpy(l, copyme->l, 3*sizeof(signed int));

	ops->ops[n].h = h;
	ops->ops[n].k = k;
	ops->ops[n].l = l;
	ops->ops[n].order = copyme->order;

	ops->n_ops++;
}


/**
 * num_equivs:
 * @ops: A %SymOpList
 * @m: A %SymOpMask, which has been shown to special_position()
 *
 * Returns: the number of equivalent reflections for a general reflection
 * in point group "ops", which were not flagged by your call to
 * special_position().
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


static signed int *v(signed int h, signed int k, signed int i, signed int l)
{
	signed int *vec = malloc(3*sizeof(signed int));
	if ( vec == NULL ) return NULL;
	/* Convert back to 3-index form now */
	vec[0] = h-i;  vec[1] = k-i;  vec[2] = l;
	return vec;
}


static void combine_ops(signed int *h1, signed int *k1, signed int *l1,
                        signed int *h2, signed int *k2, signed int *l2,
                        signed int *hnew, signed int *knew, signed int *lnew)
{
	/* Yay matrices */
	hnew[0] = h1[0]*h2[0] + h1[1]*k2[0] + h1[2]*l2[0];
	hnew[1] = h1[0]*h2[1] + h1[1]*k2[1] + h1[2]*l2[1];
	hnew[2] = h1[0]*h2[2] + h1[1]*k2[2] + h1[2]*l2[2];

	knew[0] = k1[0]*h2[0] + k1[1]*k2[0] + k1[2]*l2[0];
	knew[1] = k1[0]*h2[1] + k1[1]*k2[1] + k1[2]*l2[1];
	knew[2] = k1[0]*h2[2] + k1[1]*k2[2] + k1[2]*l2[2];

	lnew[0] = l1[0]*h2[0] + l1[1]*k2[0] + l1[2]*l2[0];
	lnew[1] = l1[0]*h2[1] + l1[1]*k2[1] + l1[2]*l2[1];
	lnew[2] = l1[0]*h2[2] + l1[1]*k2[2] + l1[2]*l2[2];
}


static void combine_and_add_symop(struct sym_op *opi, int oi,
                                  struct sym_op *opj,
                                  SymOpList *s)
{
	int i;
	signed int *h, *k, *l;

	h = malloc(3*sizeof(signed int));
	k = malloc(3*sizeof(signed int));
	l = malloc(3*sizeof(signed int));
	assert(h != NULL);
	assert(k != NULL);
	assert(l != NULL);

	memcpy(h, opj->h, 3*sizeof(signed int));
	memcpy(k, opj->k, 3*sizeof(signed int));
	memcpy(l, opj->l, 3*sizeof(signed int));

	for ( i=0; i<oi; i++ ) {

		signed int hfs[3], kfs[3], lfs[3];

		combine_ops(h, k, l, opi->h, opi->k, opi->l, hfs, kfs, lfs);

		memcpy(h, hfs, 3*sizeof(signed int));
		memcpy(k, kfs, 3*sizeof(signed int));
		memcpy(l, lfs, 3*sizeof(signed int));

	}

//	STATUS("Creating %3i %3i %3i\n", h[0], h[1], h[2]);
//	STATUS("         %3i %3i %3i\n", k[0], k[1], k[2]);
//	STATUS("         %3i %3i %3i\n", l[0], l[1], l[2]);

	add_symop(s, h, k, l, 1);
}


/* Fill in the other operations for a point group starting from its
 * generators */
static SymOpList *expand_ops(SymOpList *s)
{
	int n, i;
	SymOpList *e;

	e = new_symoplist();
	if ( e == NULL ) return NULL;
	e->name = strdup(symmetry_name(s));

	add_symop(e, v(1,0,0,0), v(0,1,0,0), v(0,0,0,1), 1);  /* I */

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {

		int j, nj;
		struct sym_op *opi = &s->ops[i];

		/* Apply op 'i' to all the current ops in the list */
		nj = num_ops(e);
		for ( j=0; j<nj; j++ ) {

			int oi;

			for ( oi=0; oi<opi->order-1; oi++ ) {
				combine_and_add_symop(opi, oi+1, &e->ops[j], e);
			}

		}

	}

	free_symoplist(s);

	return e;
}


/* Transform all the operations in a SymOpList by a given matrix.
 * The matrix must have a determinant of +/- 1 (otherwise its inverse would
 * not also be an integer matrix). */
static void transform_ops(SymOpList *s, signed int *na,
                                        signed int *nb,
                                        signed int *nc)
{
	int n, i;
	IntegerMatrix *t, *inv;
	signed int det;

	t = intmat_new(3, 3);
	if ( t == NULL ) {
		ERROR("Failed to allocate matrix.\n");
		return;
	}

	intmat_set(t, 0, 0, na[0]);
	intmat_set(t, 1, 0, na[1]);
	intmat_set(t, 2, 0, na[2]);
	intmat_set(t, 0, 1, nb[0]);
	intmat_set(t, 1, 1, nb[1]);
	intmat_set(t, 2, 1, nb[2]);
	intmat_set(t, 0, 2, nc[0]);
	intmat_set(t, 1, 2, nc[1]);
	intmat_set(t, 2, 2, nc[2]);

	det = intmat_det(t);
	if ( det == -1 ) {
		ERROR("Warning: mirrored SymOpList.\n");
	} else if ( det != 1 ) {
		ERROR("Invalid transformation for SymOpList.\n");
		return;
	}

	inv = intmat_inverse(t);
	if ( inv == NULL ) {
		ERROR("Failed to invert matrix.\n");
		return;
	}

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {

		IntegerMatrix *m, *r, *f;

		m = intmat_new(3, 3);
		if ( m == NULL ) {
			ERROR("Failed to allocate matrix.\n");
			return;
		}

		intmat_set(m, 0, 0, s->ops[i].h[0]);
		intmat_set(m, 1, 0, s->ops[i].h[1]);
		intmat_set(m, 2, 0, s->ops[i].h[2]);
		intmat_set(m, 0, 1, s->ops[i].k[0]);
		intmat_set(m, 1, 1, s->ops[i].k[1]);
		intmat_set(m, 2, 1, s->ops[i].k[2]);
		intmat_set(m, 0, 2, s->ops[i].l[0]);
		intmat_set(m, 1, 2, s->ops[i].l[1]);
		intmat_set(m, 2, 2, s->ops[i].l[2]);

		r = intmat_intmat_mult(m, t);
		if ( r == NULL ) {
			ERROR("Matrix multiplication failed.\n");
			return;
		}
		intmat_free(m);

		f = intmat_intmat_mult(inv, r);
		if ( f == NULL ) {
			ERROR("Matrix multiplication failed.\n");
			return;
		}
		intmat_free(r);

		s->ops[i].h[0] = intmat_get(f, 0, 0);
		s->ops[i].h[1] = intmat_get(f, 1, 0);
		s->ops[i].h[2] = intmat_get(f, 2, 0);
		s->ops[i].k[0] = intmat_get(f, 0, 1);
		s->ops[i].k[1] = intmat_get(f, 1, 1);
		s->ops[i].k[2] = intmat_get(f, 2, 1);
		s->ops[i].l[0] = intmat_get(f, 0, 2);
		s->ops[i].l[1] = intmat_get(f, 1, 2);
		s->ops[i].l[2] = intmat_get(f, 2, 2);
		intmat_free(f);

	}

	intmat_free(t);
	intmat_free(inv);
}


/********************************* Triclinic **********************************/

static SymOpList *make_1bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2);  /* -I */
	new->name = strdup("-1");
	return expand_ops(new);
}


static SymOpList *make_1()
{
	SymOpList *new = new_symoplist();
	new->name = strdup("1");
	return expand_ops(new);
}


/********************************* Monoclinic *********************************/

static SymOpList *make_2m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 // l */
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);  /* m -| l */
	new->name = strdup("2/m");
	return expand_ops(new);
}


static SymOpList *make_2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 // l */
	new->name = strdup("2");
	return expand_ops(new);
}


static SymOpList *make_m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);  /* m -| l */
	new->name = strdup("m");
	return expand_ops(new);
}


/******************************** Orthorhombic ********************************/

static SymOpList *make_mmm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 // k */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m -| k */
	new->name = strdup("mmm");
	return expand_ops(new);
}


static SymOpList *make_222()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 // k */
	new->name = strdup("222");
	return expand_ops(new);
}


static SymOpList *make_mm2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 // l */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m -| k */
	new->name = strdup("mm2");
	return expand_ops(new);
}


/********************************* Tetragonal *********************************/

static SymOpList *make_4m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* m -| l */
	new->name = strdup("4/m");
	return expand_ops(new);
}


static SymOpList *make_4()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	new->name = strdup("4");
	return expand_ops(new);
}


static SymOpList *make_4mm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,1), 2); /* m -| l */
	new->name = strdup("4mm");
	return expand_ops(new);
}


static SymOpList *make_422()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4);  /* 4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 // k */
	new->name = strdup("422");
	return expand_ops(new);
}


static SymOpList *make_4bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1), 4); /* -4 // l */
	new->name = strdup("-4");
	return expand_ops(new);
}


static SymOpList *make_4bar2m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1), 4); /* -4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 // k */
	new->name = strdup("-42m");
	return expand_ops(new);
}


static SymOpList *make_4barm2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1), 4); /* -4 // l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1), 2); /* 2 // h+k */
	new->name = strdup("-4m2");
	return expand_ops(new);
}


static SymOpList *make_4mmm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,1), 2); /* m -| k */
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* m -| l */
	new->name = strdup("4/mmm");
	return expand_ops(new);
}


/************************** Trigonal (Rhombohedral) ***************************/

static SymOpList *make_3_R()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0), 3); /* 3 // h+k+l */
	new->name = strdup("3_R");
	return expand_ops(new);
}


static SymOpList *make_3bar_R()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0), 3); /* -3 // h+k+l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	new->name = strdup("-3_R");
	return expand_ops(new);
}


static SymOpList *make_32_R()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0), 3); /* 3 // h+k+l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1), 2); /* 2 -| 3 */
	new->name = strdup("32_R");
	return expand_ops(new);
}


static SymOpList *make_3m_R()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0), 3); /* 3 // h+k+l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1), 2); /* m */
	new->name = strdup("3m_R");
	return expand_ops(new);
}


static SymOpList *make_3barm_R()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,0,1), v(1,0,0,0), v(0,1,0,0), 3); /* -3 // h+k+l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1), 2);    /* m */
	new->name = strdup("-3m_R");
	return expand_ops(new);
}


/*************************** Trigonal (Hexagonal) *****************************/

static SymOpList *make_3_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3); /* 3 // l */
	new->name = strdup("3_H");
	return expand_ops(new);
}


static SymOpList *make_3bar_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3);    /* 3 // l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	new->name = strdup("-3_H");
	return expand_ops(new);
}


static SymOpList *make_321_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3);  /* 3 // l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1), 2); /* 2 // h */
	new->name = strdup("321_H");
	return expand_ops(new);
}


static SymOpList *make_312_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3);    /* 3 // l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1), 2); /* 2 // h+k */
	new->name = strdup("312_H");
	return expand_ops(new);
}


static SymOpList *make_3m1_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3); /* 3 // l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1), 2); /* m -| i */
	new->name = strdup("3m1_H");
	return expand_ops(new);
}


static SymOpList *make_31m_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3); /* 3 // l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1), 2); /* m -| (k+i) */
	new->name = strdup("31m_H");
	return expand_ops(new);
}


static SymOpList *make_3barm1_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3);    /* 3 // l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1), 2);   /* 2 // h */
	new->name = strdup("-3m1_H");
	return expand_ops(new);
}


static SymOpList *make_3bar1m_H()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,1), 3);    /* 3 // l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,-1), 2); /* 2 // h+k */
	new->name = strdup("-31m_H");
	return expand_ops(new);
}


/********************************** Hexgonal **********************************/

static SymOpList *make_6()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1), 6); /* 6 // l */
	new->name = strdup("6");
	return expand_ops(new);
}


static SymOpList *make_6bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1), 6); /* -6 // l */
	new->name = strdup("-6");
	return expand_ops(new);
}


static SymOpList *make_6m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1), 6); /* 6 // l */
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);  /* m -| l */
	new->name = strdup("6/m");
	return expand_ops(new);
}


static SymOpList *make_622()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1), 6); /* 6 // l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,-1), 2);   /* 2 // h */
	new->name = strdup("622");
	return expand_ops(new);
}


static SymOpList *make_6mm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,-1,0), v(-1,0,0,0), v(0,0,0,1), 6); /* 6 // l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1), 2); /* m -| i */
	new->name = strdup("6mm");
	return expand_ops(new);
}


static SymOpList *make_6barm2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1), 6); /* -6 // l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1), 2); /* m -| i */
	new->name = strdup("-6m2");
	return expand_ops(new);
}


static SymOpList *make_6bar2m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1), 6); /* -6 // l */
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1), 2);  /* m -| (k+i) */
	new->name = strdup("-62m");
	return expand_ops(new);
}


static SymOpList *make_6mmm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,0,1,0), v(1,0,0,0), v(0,0,0,-1), 6); /* -6 // l */
	add_symop(new, v(0,-1,0,0), v(-1,0,0,0), v(0,0,0,1), 2); /* m -| i */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	new->name = strdup("6/mmm");
	return expand_ops(new);
}


/************************************ Cubic ***********************************/

static SymOpList *make_23()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2// l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2// k */
	add_symop(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0), 3); /* 3// h+k+l */
	new->name = strdup("23");
	return expand_ops(new);
}


static SymOpList *make_m3bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2// l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2// k */
	add_symop(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0), 3); /* 3// h+k+l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	new->name = strdup("m-3");
	return expand_ops(new);
}


static SymOpList *make_432()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);/* 2 // k */
	add_symop(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0), 3);  /* 3 // h+k+l */
	new->name = strdup("432");
	return expand_ops(new);
}


static SymOpList *make_4bar3m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,1,0,0), v(-1,0,0,0), v(0,0,0,-1), 4); /* -4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);/* 2 // k */
	add_symop(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0), 3);  /* 3 // h+k+l */
	new->name = strdup("-43m");
	return expand_ops(new);
}


static SymOpList *make_m3barm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 // l */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2);/* 2 // k */
	add_symop(new, v(0,1,0,0), v(0,0,0,1), v(1,0,0,0), 3);  /* 3 // h+k+l */
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1), 2); /* -I */
	new->name = strdup("m-3m");
	return expand_ops(new);
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

	if ( strncmp(sym+s, "ua", 2) == 0 ) {
		ua = sym[s+2];
	} else {
		ERROR("Unrecognised point group '%s'\n", sym);
		return NULL;
	}

	pg_type = strndup(sym, s-1);
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
	free(pg_type);

	switch ( ua ) {

		case 'a' :
		transform_ops(pg, v(0,0,0,1),
		                  v(0,1,0,0),
		                  v(-1,0,0,0));
		break;

		case 'b' :
		transform_ops(pg, v(1,0,0,0),
		                  v(0,0,0,1),
		                  v(0,-1,0,0));
		break;

		case 'c' :
		/* No transformation needed */
		break;

		default :
		ERROR("Bad unique axis '%c'\n", ua);
		free_symoplist(pg);
		return NULL;

	}

	return pg;
}


/**
 * get_pointgroup:
 * @sym: A string representation of a point group
 *
 * This function parses @sym and returns the corresponding %SymOpList.
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


static void do_op(const struct sym_op *op,
                  signed int h, signed int k, signed int l,
                  signed int *he, signed int *ke, signed int *le)
{
	*he = h*op->h[0] + k*op->h[1] + l*op->h[2];
	*ke = h*op->k[0] + k*op->k[1] + l*op->k[2];
	*le = h*op->l[0] + k*op->l[1] + l*op->l[2];
}


/**
 * get_equiv:
 * @ops: A %SymOpList
 * @m: A %SymOpMask, which has been shown to special_position()
 * @idx: Index of the operation to use
 * @h: index of reflection
 * @k: index of reflection
 * @l: index of reflection
 * @he: location to store h index of equivalent reflection
 * @ke: location to store k index of equivalent reflection
 * @le: location to store l index of equivalent reflection
 *
 * This function applies the @idx-th symmetry operation from @ops to the
 * reflection @h, @k, @l, and stores the result at @he, @ke and @le.
 *
 * If you don't mind that the same equivalent might appear twice, simply call
 * this function the number of times returned by num_ops(), using the actual
 * point group.  If repeating the same equivalent twice (for example, if the
 * given reflection is a special high-symmetry one), call special_position()
 * first to get a "specialised" SymOpList and use that instead.
 **/
void get_equiv(const SymOpList *ops, const SymOpMask *m, int idx,
               signed int h, signed int k, signed int l,
               signed int *he, signed int *ke, signed int *le)
{
	const int n = num_ops(ops);

	if ( m != NULL ) {

		int i, c;

		c = 0;
		for ( i=0; i<n; i++ ) {

			if ( (c == idx) && m->mask[i] ) {
				do_op(&ops->ops[i], h, k, l, he, ke, le);
				return;
			}

			if ( m->mask[i] ) {
				c++;
			}

		}

		ERROR("Index %i out of range for point group '%s' with"
		      " reflection %i %i %i\n",
		      idx, symmetry_name(ops), h, k, l);

		*he = 0;  *ke = 0;  *le = 0;

		return;

	}



	if ( idx >= n ) {

		ERROR("Index %i out of range for point group '%s'\n", idx,
		      symmetry_name(ops));

		*he = 0;  *ke = 0;  *le = 0;
		return;

	}

	do_op(&ops->ops[idx], h, k, l, he, ke, le);
}


/**
 * special_position:
 * @ops: A %SymOpList, usually corresponding to a point group
 * @m: A %SymOpMask created with new_symopmask()
 * @h: index of a reflection
 * @k: index of a reflection
 * @l: index of a reflection
 *
 * This function determines which operations in @ops map the reflection @h, @k,
 * @l onto itself, and uses @m to flag the operations in @ops which cause this.
 *
 **/
void special_position(const SymOpList *ops, SymOpMask *m,
                      signed int h, signed int k, signed int l)
{
	int i, n;
	signed int *htest;
	signed int *ktest;
	signed int *ltest;

	assert(m->list = ops);

	n = num_equivs(ops, NULL);
	htest = malloc(n*sizeof(signed int));
	ktest = malloc(n*sizeof(signed int));
	ltest = malloc(n*sizeof(signed int));

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

	free(htest);
	free(ktest);
	free(ltest);
}


static int any_negative(signed int h, signed int k, signed int l)
{
	if ( h < 0 ) return 1;
	if ( k < 0 ) return 1;
	if ( l < 0 ) return 1;
	return 0;
}


/**
 * get_asymm:
 * @ops: A %SymOpList, usually corresponding to a point group
 * @h: index of a reflection
 * @k: index of a reflection
 * @l: index of a reflection
 * @hp: location for asymmetric index of reflection
 * @kp: location for asymmetric index of reflection
 * @lp: location for asymmetric index of reflection
 *
 * This function determines the asymmetric version of the reflection @h, @k, @l
 * in symmetry group @ops, and puts the result in @hp, @kp, @lp.
 *
 * This is a relatively expensive operation because of its generality.
 * Therefore, if you know you'll need to make repeated use of the asymmetric
 * indices, consider creating a new %RefList indexed according to the asymmetric
 * indices themselves with asymmetric_indices().  If you do that, you'll still
 * be able to get the original versions of the indices with
 * get_symmetric_indices().
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


static int is_inversion(const struct sym_op *op)
{
	if ( (op->h[0]!=-1) || (op->h[1]!=0) || (op->h[2]!=0) ) return 0;
	if ( (op->k[0]!=0) || (op->k[1]!=-1) || (op->k[2]!=0) ) return 0;
	if ( (op->l[0]!=0) || (op->l[1]!=0) || (op->l[2]!=-1) ) return 0;
	return 1;
}


static int is_identity(const struct sym_op *op)
{
	if ( (op->h[0]!=1) || (op->h[1]!=0) || (op->h[2]!=0) ) return 0;
	if ( (op->k[0]!=0) || (op->k[1]!=1) || (op->k[2]!=0) ) return 0;
	if ( (op->l[0]!=0) || (op->l[1]!=0) || (op->l[2]!=1) ) return 0;
	return 1;
}


static signed int determinant(const struct sym_op *op)
{
	signed int det = 0;

	det += op->h[0] * (op->k[1]*op->l[2] - op->k[2]*op->l[1]);
	det -= op->h[1] * (op->k[0]*op->l[2] - op->k[2]*op->l[0]);
	det += op->h[2] * (op->k[0]*op->l[1] - op->k[1]*op->l[0]);

	return det;
}


/**
 * is_centrosymmetric:
 * @s: A %SymOpList
 *
 * Returns: non-zero if @s contains an inversion operation
 */
int is_centrosymmetric(const SymOpList *s)
{
	int i, n;

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {
		if ( is_inversion(&s->ops[i]) ) return 1;
	}

	return 0;
}


static int ops_equal(const struct sym_op *op,
                     signed int *h, signed int *k, signed int *l)
{
	if ( (op->h[0]!=h[0]) || (op->h[1]!=h[1]) || (op->h[2]!=h[2]) )
		return 0;
	if ( (op->k[0]!=k[0]) || (op->k[1]!=k[1]) || (op->k[2]!=k[2]) )
		return 0;
	if ( (op->l[0]!=l[0]) || (op->l[1]!=l[1]) || (op->l[2]!=l[2]) )
		return 0;
	return 1;
}


static int struct_ops_equal(const struct sym_op *op1, const struct sym_op *op2)
{
	return ops_equal(op1, op2->h, op2->k, op2->l);
}


/* Return true if a*b = ans */
static int check_mult(const struct sym_op *ans,
                      const struct sym_op *a, const struct sym_op *b)
{
	signed int *ans_h, *ans_k, *ans_l;
	int val;

	ans_h = malloc(3*sizeof(signed int));
	ans_k = malloc(3*sizeof(signed int));
	ans_l = malloc(3*sizeof(signed int));

	combine_ops(a->h, a->k, a->l, b->h, b->k, b->l, ans_h, ans_k, ans_l);
	val = ops_equal(ans, ans_h, ans_k, ans_l);

	free(ans_h);
	free(ans_k);
	free(ans_l);

	return val;
}


/**
 * is_subgroup:
 * @source: A %SymOpList
 * @target: Another %SymOpList, which might be a subgroup of @source.
 *
 * Returns: non-zero if every operation in @target is also in @source.
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

			if ( struct_ops_equal(&target->ops[i],
			                      &source->ops[j] ) )
			{
				found = 1;
				break;
			}

		}

		if ( !found ) return 0;

	}

	return 1;
}


/**
 * get_ambiguities:
 * @source: The "source" symmetry, a %SymOpList
 * @target: The "target" symmetry, a %SymOpList

 * Calculates twinning laws.  Returns a %SymOpList containing the twinning
 * operators, which are the symmetry operations which can be added to @target
 * to generate @source.  Only rotations are allowable - no mirrors nor
 * inversions.
 * To count the number of possibilities, use num_ops() on the result.
 *
 * Returns: A %SymOpList containing the twinning operators, or NULL if the
 * source symmetry cannot be generated from that target symmetry without using
 * mirror or inversion operations.
 */
SymOpList *get_ambiguities(const SymOpList *source, const SymOpList *target)
{
	int n_src, n_tgt;
	int i;
	SymOpList *twins;
	SymOpList *src_reordered;
	SymOpMask *used;
	char *name;

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

	src_reordered = new_symoplist();
	used = new_symopmask(source);

	/* Find identity */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( is_identity(&source->ops[i]) ) {
			add_copied_op(src_reordered, &source->ops[i]);
			used->mask[i] = 0;
		}
	}

	/* Find binary options (order=2) of first kind (determinant positive) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( (source->ops[i].order == 2)
		  && (determinant(&source->ops[i]) > 0) ) {
			add_copied_op(src_reordered, &source->ops[i]);
			used->mask[i] = 0;
		}
	}

	/* Find other operations of first kind (determinant positive) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( determinant(&source->ops[i]) > 0 ) {
			add_copied_op(src_reordered, &source->ops[i]);
			used->mask[i] = 0;
		}
	}

	/* Find inversion */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( is_inversion(&source->ops[i]) ) {
			add_copied_op(src_reordered, &source->ops[i]);
			used->mask[i] = 0;
		}
	}

	/* Find binary options of second kind (determinant negative) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( (source->ops[i].order == 2)
		  && (determinant(&source->ops[i]) < 0) ) {
			add_copied_op(src_reordered, &source->ops[i]);
			used->mask[i] = 0;
		}
	}

	/* Find other operations of second kind (determinant negative) */
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( determinant(&source->ops[i]) < 0 ) {
			add_copied_op(src_reordered, &source->ops[i]);
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
	used = new_symopmask(src_reordered);

	/* This is the first method from Flack (1987) */
	for ( i=0; i<n_src; i++ ) {

		int j;
		if ( used->mask[i] == 0 ) continue;

		for ( j=1; j<n_tgt; j++ ) {

			int k;
			for ( k=i+1; k<n_src; k++ ) {
				if ( check_mult(&src_reordered->ops[k],
				                &src_reordered->ops[i],
				                &target->ops[j]) )
				{
					used->mask[k] = 0;
				}
			}

		}

	}

	twins = new_symoplist();
	for ( i=0; i<n_src; i++ ) {
		if ( used->mask[i] == 0 ) continue;
		if ( determinant(&src_reordered->ops[i]) < 0 ) {
			/* A mirror or inversion turned up in the list.
			 * That means that no pure rotational ambiguity can
			 * account for this subgroup relationship. */
			free_symoplist(twins);
			free_symopmask(used);
			free_symoplist(src_reordered);
			return NULL;
		}
		add_copied_op(twins, &src_reordered->ops[i]);
	}

	free_symopmask(used);
	free_symoplist(src_reordered);

	name = malloc(64);
	snprintf(name, 63, "%s -> %s", symmetry_name(source),
	                               symmetry_name(target));
	twins->name = name;

	return twins;
}


static void add_chars(char *t, const char *s, int max_len)
{
	char *tmp;

	tmp = strdup(t);

	snprintf(t, max_len, "%s%s", tmp, s);
	free(tmp);
}


static char *get_matrix_name(signed int *v)
{
	char *text;
	const int max_len = 9;
	int i;
	int printed = 0;

	text = malloc(max_len+1);
	text[0] = '\0';

	for ( i=0; i<3; i++ ) {

		if ( v[i] == 0 ) continue;

		if ( (i==0) && (v[0]==v[1]) ) {
			if ( v[i]>0 ) add_chars(text, "-", max_len);
			add_chars(text, "i", max_len);
			v[1] -= v[0];
			continue;
		}

		if ( v[i]<0 ) {
			add_chars(text, "-", max_len);
		} else {
			if ( printed ) add_chars(text, "+", max_len);
		}

		if ( abs(v[i])>1 ) {
			char num[3];
			snprintf(num, 2, "%i", abs(v[i]));
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


static char *name_equiv(const struct sym_op *op)
{
	char *h, *k, *l;
	char *name;

	h = get_matrix_name(op->h);
	k = get_matrix_name(op->k);
	l = get_matrix_name(op->l);
	name = malloc(32);

	if ( strlen(h)+strlen(k)+strlen(l) == 3 ) {
		snprintf(name, 31, "%s%s%s", h, k, l);
	} else {
		snprintf(name, 31, "%s,%s,%s", h, k, l);
	}
	free(h);
	free(k);
	free(l);

	return name;
}


/**
 * describe_symmetry:
 * @s: A %SymOpList
 *
 * Writes the name and a list of operations to stderr.
 */
void describe_symmetry(const SymOpList *s)
{
	int i, n;

	n = num_equivs(s, NULL);

	STATUS("%15s :", symmetry_name(s));

	for ( i=0; i<n; i++ ) {
		size_t len;
		char *name = name_equiv(&s->ops[i]);
		len = strlen(name);
		if ( len > max_len ) max_len = len;
		free(name);
	}

	for ( i=0; i<n; i++ ) {

		char *name;
		size_t n, j;

		name = name_equiv(&s->ops[i]);
		n = max_len - strlen(j) + 1;

		for ( j=0; j<n; j++ ) {
			STATUS(" ");
		}
		STATUS("%s", name);
		free(name);
		if ( (i!=0) && (i%8==0) ) STATUS("\n%15s  ", "");
	}
	STATUS("\n");
}


/**
 * symmetry_name:
 * @ops: A %SymOpList
 *
 * Returns: a text description of @ops.
 */
const char *symmetry_name(const SymOpList *ops)
{
	return ops->name;
}


/**
 * add_op_intmat:
 * @s: A %SymOpList
 * @m: An %IntegerMatrix
 *
 * Adds @m to the @s.  Operations are NOT cross-multiplied, i.e. the result
 * might not be a point group.
 */

void add_op_intmat(SymOpList *s, const IntegerMatrix *m)
{
	add_symop(s,
	    v(intmat_get(m, 0, 0), intmat_get(m, 0, 1), 0, intmat_get(m, 0, 2)),
	    v(intmat_get(m, 1, 0), intmat_get(m, 1, 1), 0, intmat_get(m, 1, 2)),
	    v(intmat_get(m, 2, 0), intmat_get(m, 2, 1), 0, intmat_get(m, 2, 2)),
	          1);
}
