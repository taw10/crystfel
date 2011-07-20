/*
 * symmetry.c
 *
 * Symmetry
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
	signed int *l;  /* Contributions to h, k and l from h, k, i and l */
	int order;
};


/**
 * SECTION:symoplist
 * @short_description: A list of point symmetry operations
 * @title: SymOpList
 * @section_id:
 * @see_also:
 * @include: "symmetry.h"
 * @Image:
 *
 * The SymOpList is an opaque data structure containing a list of point symmetry
 * operations.  It could represent an point group or a list of indexing
 * ambiguities (twin laws), or similar.
 */

struct _symoplist
{
	struct sym_op *ops;
	int n_ops;
	int max_ops;
	char *name;
	int *divisors;
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
	ops->divisors = realloc(ops->divisors, ops->max_ops*sizeof(int));
}


/**
 * new_symopmask:
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
	new->divisors = NULL;
	new->name = NULL;
	new->num_equivs = 1;
	alloc_ops(new);
	return new;
}


/**
 * free_symoplist:
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


/**
 * num_equivs:
 *
 * Returns: the number of equivalent reflections for a general reflection
 * in point group "ops".
 **/
int num_equivs(const SymOpList *ops, const SymOpMask *m)
{
	return num_ops(ops);
}


static signed int *v(signed int h, signed int k, signed int i, signed int l)
{
	signed int *vec = malloc(3*sizeof(signed int));
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
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m */
	new->name = strdup("2/m");
	return expand_ops(new);
}


static SymOpList *make_2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 */
	new->name = strdup("2");
	return expand_ops(new);
}


static SymOpList *make_m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m */
	new->name = strdup("m");
	return expand_ops(new);
}


/******************************** Orthorhombic ********************************/

static SymOpList *make_mmm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m */
	new->name = strdup("mmm");
	return expand_ops(new);
}


static SymOpList *make_222()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1), 2); /* 2 */
	new->name = strdup("222");
	return expand_ops(new);
}


static SymOpList *make_mm2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1), 2);  /* m */
	new->name = strdup("mm2");
	return expand_ops(new);
}


/********************************* Tetragonal *********************************/

static SymOpList *make_4m()
{
	return NULL;
}


static SymOpList *make_4()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,-1,0,0), v(1,0,0,0), v(0,0,0,1), 4); /* 4 */
	new->name = strdup("4");
	return expand_ops(new);
}


static SymOpList *make_4bar()
{
	return NULL;
}


static SymOpList *make_4mmm()
{
	return NULL;
}


static SymOpList *make_422()
{
	return NULL;
}


static SymOpList *make_4bar2m()
{
	return NULL;
}


static SymOpList *make_4mm()
{
	return NULL;
}


/******************************** Rhombohedral ********************************/


/********************************** Hexgonal **********************************/

static SymOpList *make_6m()
{
	return NULL;
}


static SymOpList *make_6()
{
	return NULL;
}


static SymOpList *make_6bar()
{
	return NULL;
}


static SymOpList *make_6mmm()
{
	return NULL;
}


static SymOpList *make_622()
{
	return NULL;
}


static SymOpList *make_6bar2m()
{
	return NULL;
}


static SymOpList *make_6mm()
{
	return NULL;
}


/************************************ Cubic ***********************************/

SymOpList *get_pointgroup(const char *sym)
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
	if ( strcmp(sym, "4bar") == 0 ) return make_4bar();
	if ( strcmp(sym, "4/mmm") == 0 ) return make_4mmm();
	if ( strcmp(sym, "422") == 0 ) return make_422();
	if ( strcmp(sym, "4bar2m") == 0 ) return make_4bar2m();
	if ( strcmp(sym, "4mm") == 0 ) return make_4mm();

	/* Hexagonal */
	if ( strcmp(sym, "6/m") == 0 ) return make_6m();
	if ( strcmp(sym, "6") == 0 ) return make_6();
	if ( strcmp(sym, "6bar") == 0 ) return make_6bar();
	if ( strcmp(sym, "6/mmm") == 0 ) return make_6mmm();
	if ( strcmp(sym, "622") == 0 ) return make_622();
	if ( strcmp(sym, "6bar2m") == 0 ) return make_6bar2m();
	if ( strcmp(sym, "6mm") == 0 ) return make_6mm();

	ERROR("Unknown point group '%s'\n", sym);
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
	int i, n;

	n = num_ops(ops);
	for ( i=idx; i<n; i++ ) {
		if ( (m == NULL) || m->mask[i] ) {
			do_op(&ops->ops[i], h, k, l, he, ke, le);
			return;
		}
	}

	ERROR("Index %i out of range for point group '%s'\n", idx,
	      symmetry_name(ops));

	*he = 0;  *ke = 0;  *le = 0;
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
 * @l onto itself, and returns a new %SymOpList containing only the operations
 * from @ops which do not do so.
 *
 * Returns: the "specialised" %SymOpList.
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


void get_asymm(const SymOpList *ops,
               signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp)
{
	int nequiv;
	int p;
	signed int best_h, best_k, best_l;

	nequiv = num_equivs(ops, NULL);

	best_h = h;  best_k = k;  best_l = l;
	for ( p=0; p<nequiv; p++ ) {

		get_equiv(ops, NULL, p, h, k, l, hp, kp, lp);

		if ( h > best_h ) {
			best_h = h;  best_k = k;  best_l = l;
			continue;
		}

		if ( k > best_k ) {
			best_h = h;  best_k = k;  best_l = l;
			continue;
		}

		if ( l > best_l ) {
			best_h = h;  best_k = k;  best_l = l;
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


int is_centrosymmetric(const SymOpList *s)
{
	int i, n;

	n = num_ops(s);
	for ( i=0; i<n; i++ ) {
		if ( is_inversion(&s->ops[i]) ) return 1;
	}

	return 0;
}


/**
 * get_twins:
 *
 * Calculate twinning laws.
 *
 * To count the number of possibilities, use num_ops() on the result.
 */
SymOpList *get_twins(const SymOpList *source, const SymOpList *target)
{
	int n_src, n_tgt;
	int i;
	SymOpList *twins = new_symoplist();

	n_src = num_ops(source);
	n_tgt = num_ops(target);

	for ( i=0; i<n_src; i++ ) {


	}

	return twins;
}


const char *symmetry_name(const SymOpList *ops)
{
	return ops->name;
}
