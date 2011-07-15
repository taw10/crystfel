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
	enum lattice_type latt;
	char op[6];
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


/* Add a operation to a SymOpList */
static void add_symop(SymOpList *ops,
                      signed int *h, signed int *k, signed int *l)
{
	if ( ops->n_ops == ops->max_ops ) {
		/* Pretty sure this never happens, but still... */
		ops->max_ops += 16;
		alloc_ops(ops);
	}

	ops->ops[ops->n_ops].h = h;
	ops->ops[ops->n_ops].k = k;
	ops->ops[ops->n_ops].l = l;
	ops->n_ops++;
}


int num_ops(const SymOpList *ops)
{
	return ops->n_ops;
}


static signed int *v(signed int h, signed int k, signed int i, signed int l)
{
	signed int *vec = malloc(4*sizeof(signed int));
	vec[0] = h;  vec[1] = k;  vec[2] = i;  vec[3] = l;
	return vec;
}


/********************************* Triclinic **********************************/

static SymOpList *make_1bar()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,1));
	add_symop(new, v(-1,0,0,0), v(0,-1,0,0), v(0,0,0,-1));
	return new;
}


static SymOpList *make_1()
{
	SymOpList *new = new_symoplist();
	add_symop(new, v(1,0,0,0), v(0,1,0,0), v(0,0,0,1));
	return new;
}


/********************************* Monoclinic *********************************/

static SymOpList *make_2m()
{
	return NULL;
}


static SymOpList *make_2()
{
	return NULL;
}


static SymOpList *make_m()
{
	return NULL;
}


/******************************** Orthorhombic ********************************/

static SymOpList *make_mmm()
{
	return NULL;
}


static SymOpList *make_222()
{
	return NULL;
}


static SymOpList *make_mm2()
{
	return NULL;
}


/********************************* Tetragonal *********************************/

static SymOpList *make_4m()
{
	return NULL;
}


static SymOpList *make_4()
{
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
}


/**
 * num_ops:
 * @ops: A %SymOpList
 *
 * Returns: the number of operations in @ops.
 **/
int num_ops(SymOpList *ops)
{
	return ops->n_ops;
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
	signed int i = -h-k;
	struct sym_op op = ops->ops[idx];

	*he = h*op.h[0] + k*op.h[1] + i*op.h[2] + l*op.h[3];
	*ke = h*op.k[0] + k*op.h[1] + i*op.k[2] + l*op.k[3];
	*le = h*op.l[0] + k*op.h[1] + i*op.l[2] + l*op.l[3];
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
	int n_general;
	int i;
	SymOpList *equivs;
	int n_equivs = 0;

	equivs = new_symoplist();

	for ( i=0; i<num_ops(ops); i++ ) {

		signed int ht, kt, lt;

		/* Get equivalent according to the point group */
		get_equiv(ops, i, h, k, l, &ht, &kt, &lt);

		if (

	}

	return equivs;
}


void get_asymm(SymOpList *ops, int idx,
               signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp)
{
	int nequiv = num_equivs(h, k, l, sym);
	int p;
	signed int best_h, best_k, best_l;

	best_h = h;  best_k = k;  best_l = l;
	for ( p=0; p<nequiv; p++ ) {

		get_equiv(h, k, l, hp, kp, lp, sym, p);

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
	signed int h, k, l;
	SymOpList *twins = new_symoplist();

	n_src = num_ops(source);
	n_tgt = num_ops(target);

	for ( i=0; i<n_src; i++ ) {


	}

	return twins;
}
