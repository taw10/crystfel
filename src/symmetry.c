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


enum lattice_type
{
	L_TRICLINIC,
	L_MONOCLINIC,
	L_ORTHORHOMBIC,
	L_TETRAGONAL,
	L_RHOMBOHEDRAL,
	L_TRIGONAL,
	L_HEXAGONAL,
	L_CUBIC,
};


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
};



static void alloc_ops(SymOpList *ops)
{
	ops->ops = realloc(ops->ops, ops->max_ops*sizeof(struct sym_op));
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
	free(ops);
}


static int is_identity(signed int *h, signed int *k, signed int *l)
{
	if ( (h[0]!=1) || (h[1]!=0) || (h[2]!=0) ) return 0;
	if ( (k[0]!=0) || (k[1]!=1) || (k[2]!=0) ) return 0;
	if ( (l[0]!=0) || (l[1]!=0) || (l[2]!=1) ) return 0;
	return 1;
}


/* Calculate the order of the operation "M", which is the lowest
 * integer n such that M^n = I. */
static int order_of_op(signed int *hin, signed int *kin, signed int *lin)
{
	int n;
	signed int h[3];
	signed int k[3];
	signed int l[3];

	memcpy(h, hin, 3*sizeof(signed int));
	memcpy(k, kin, 3*sizeof(signed int));
	memcpy(l, lin, 3*sizeof(signed int));

	for ( n=1; n<6; n++ ) {

		signed int hnew[3];
		signed int knew[3];
		signed int lnew[3];

		/* Yay matrices */
		hnew[0] = h[0]*h[0] + h[1]*k[0] + h[2]*l[0];
		hnew[1] = h[0]*h[1] + h[1]*k[1] + h[2]*l[1];
		hnew[2] = h[0]*h[2] + h[1]*k[2] + h[2]*l[2];

		knew[0] = k[0]*h[0] + k[1]*k[0] + k[2]*l[0];
		knew[1] = k[0]*h[1] + k[1]*k[1] + k[2]*l[1];
		knew[2] = k[0]*h[2] + k[1]*k[2] + k[2]*l[2];

		lnew[0] = l[0]*h[0] + l[1]*k[0] + l[2]*l[0];
		lnew[1] = l[0]*h[1] + l[1]*k[1] + l[2]*l[1];
		lnew[2] = l[0]*h[2] + l[1]*k[2] + l[2]*l[2];

		if ( is_identity(hnew, knew, lnew) ) break;

		memcpy(h, hnew, 3*sizeof(signed int));
		memcpy(k, knew, 3*sizeof(signed int));
		memcpy(l, lnew, 3*sizeof(signed int));

	}

	return n;
}


/* Add a operation to a SymOpList */
static void add_symop(SymOpList *ops,
                      signed int *h, signed int *k, signed int *l)
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
	ops->ops[n].order = order_of_op(h, k, l);
	ops->n_ops++;
}


static void add_copied_symop(SymOpList *ops, struct sym_op *copyme)
{
	if ( ops->n_ops == ops->max_ops ) {
		/* Pretty sure this never happens, but still... */
		ops->max_ops += 16;
		alloc_ops(ops);
	}

	memcpy(&ops->ops[ops->n_ops], copyme, sizeof(*copyme));
	ops->n_ops++;
}


/* This returns the number of operations in "ops".  To get the number of
 * symmetric equivalents this generates, use num_equivs() instead. */
static int num_ops(const SymOpList *ops)
{
	return ops->n_ops;
}


/**
 * num_equivs:
 *
 * Returns: the number of equivalent reflections for a general reflection
 * in point group "ops".
 **/
int num_equivs(const SymOpList *ops)
{
	int i, n, tot;

	n = num_ops(ops);
	tot = 1;

	for ( i=0; i<n; i++ ) {
		tot *= ops->ops[i].order;
	}

	return tot;
}


static signed int *v(signed int h, signed int k, signed int i, signed int l)
{
	signed int *vec = malloc(3*sizeof(signed int));
	/* Convert back to 3-index form now */
	vec[0] = h-i;  vec[1] = k-i;  vec[2] = l;
	return vec;
}


/********************************* Triclinic **********************************/

static SymOpList *make_1bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1));  /* -I */
	return new;
}


static SymOpList *make_1()
{
	SymOpList *new = new_symoplist();
	return new;
}


/********************************* Monoclinic *********************************/

static SymOpList *make_2m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m */
	return NULL;
}


static SymOpList *make_2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 */
	return NULL;
}


static SymOpList *make_m()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m */
	return NULL;
}


/******************************** Orthorhombic ********************************/

static SymOpList *make_mmm()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m */
	return NULL;
}


static SymOpList *make_222()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 */
	add_symop(new, v(-1,0,0,0), v(0,1,0,0), v(0,0,0,-1)); /* 2 */
	return NULL;
}


static SymOpList *make_mm2()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,1)); /* 2 */
	add_symop(new, v(1,0,0,0), v(0,-1,0,0), v(0,0,0,1));  /* m */
	return NULL;
}


/********************************* Tetragonal *********************************/

static SymOpList *make_4m()
{
	return NULL;
}


static SymOpList *make_4()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(0,1,0,0), v(1,0,0,0), v(0,0,0,1)); /* 4 */
	return NULL;
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


static void do_op(struct sym_op *op,
                  signed int h, signed int k, signed int l,
                  signed int *he, signed int *ke, signed int *le)
{
	*he = h*op->h[0] + k*op->h[1] + l*op->h[2];
	*ke = h*op->k[0] + k*op->h[1] + l*op->k[2];
	*le = h*op->l[0] + k*op->h[1] + l*op->l[2];
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
void get_equiv(SymOpList *ops, int idx,
               signed int h, signed int k, signed int l,
               signed int *he, signed int *ke, signed int *le)
{
	int sig[32];
	int divisors[32];
	int i, n, r;

	n = num_ops(ops);
	divisors[0] = 1;
	for ( i=1; i<n; i++ ) {
		divisors[i] = divisors[i-1]*ops->ops[i].order;
	}
	r = idx;
	for ( i=n-1; i>=0; i-- ) {
		sig[i] = r / divisors[i];
		r = r % divisors[i];
		assert(sig[i] < ops->ops[i].order);
	}

	for ( i=0; i<n; i++ ) {

		int s;

		/* Do this operation "sig[i]" times */
		for ( s=0; s<sig[i]; s++ ) {
			do_op(&ops->ops[i], h, k, l, &h, &k, &l);
		}

	}
}


/**
 * special_position:
 * @ops: A %SymOpList, usually corresponding to a point group
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
SymOpList *special_position(SymOpList *ops,
                            signed int h, signed int k, signed int l)
{
	int i, n;
	SymOpList *specialised;

	n = num_ops(ops);
	specialised = new_symoplist();

	for ( i=0; i<n; i++ ) {

		signed int ht, kt, lt;

		do_op(&ops->ops[i], h, k, l, &ht, &kt, &lt);
		if ( (h==ht) || (k==kt) || (l==lt) ) continue;
		add_copied_symop(specialised, &ops->ops[i]);

	}

	return specialised;
}


void get_asymm(SymOpList *ops, int idx,
               signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp)
{
	int nequiv = num_equivs(ops);
	int p;
	signed int best_h, best_k, best_l;

	best_h = h;  best_k = k;  best_l = l;
	for ( p=0; p<nequiv; p++ ) {

		get_equiv(ops, p, h, k, l, hp, kp, lp);

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


/**
 * get_twins:
 *
 * Calculate twinning laws.
 *
 * To count the number of possibilities, use num_ops() on the result.
 */
SymOpList *get_twins(SymOpList *source, SymOpList *target)
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
