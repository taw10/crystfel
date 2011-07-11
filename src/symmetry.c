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


#ifdef DEBUG
#define SYM_DEBUG STATUS
#else /* DEBUG */
#define SYM_DEBUG(...)
#endif /* DEBUG */

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


struct sym_op {
	signed int h;
	signed int k;
	signed int l;
	int op;
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
 * Wibble
 */

struct _symoplist {
	struct sym_op *items;
	int n_ops;
	int max_ops;
};



static void alloc_items(SymOpList *items)
{
	items->items = realloc(items->items,
	                       items->max_ops*sizeof(struct sym_op));
}


/**
 * new_items:
 *
 * Creates a new %SymOpList.
 *
 * Returns: The new list, or NULL.
 **/
SymOpList *new_items()
{
	SymOpList *new;
	new = malloc(sizeof(SymOpList));
	if ( new == NULL ) return NULL;
	new->max_ops = 1024;
	new->n_ops = 0;
	new->items = NULL;
	alloc_items(new);
	return new;
}


void delete_items(SymOpList *items)
{
	if ( items == NULL ) return;
	if ( items->items != NULL ) free(items->items);
	free(items);
}


void add_item_with_op(SymOpList *items, signed int h, signed int k,
                      signed int l, int op)
{
	if ( items->n_ops == items->max_ops ) {
		items->max_ops += 1024;
		alloc_items(items);
	}

	items->items[items->n_ops].h = h;
	items->items[items->n_ops].k = k;
	items->items[items->n_ops].l = l;
	items->items[items->n_ops].op = op;
	items->n_ops++;
}


void add_item(SymOpList *items, signed int h, signed int k, signed int l)
{
	add_item_with_op(items, h, k, l, 0);
}


int find_item(SymOpList *items,
                     signed int h, signed int k, signed int l)
{
	int i;

	for ( i=0; i<items->n_ops; i++ ) {
		if ( items->items[i].h != h ) continue;
		if ( items->items[i].k != k ) continue;
		if ( items->items[i].l != l ) continue;
		return 1;
	}
	return 0;
}


static int find_op(SymOpList *items, int op)
{
	int i;

	for ( i=0; i<items->n_ops; i++ ) {
		if ( items->items[i].op == op ) return 1;
	}
	return 0;
}


struct sym_op *get_item(SymOpList *items, int i)
{
	if ( i >= items->n_ops ) return NULL;
	return &items->items[i];
}


int num_items(const SymOpList *items)
{
	return items->n_ops;
}


void union_op_items(SymOpList *items, SymOpList *newi)
{
	int n, i;

	n = num_items(newi);
	for ( i=0; i<n; i++ ) {

		struct sym_op *r = get_item(newi, i);
		if ( find_op(items, r->op) ) continue;

		add_item_with_op(items, r->h, r->k, r->l, r->op);

	}
}


void union_ops(SymOpList *items, SymOpList *newi)
{
	int n, i;

	n = num_items(newi);
	for ( i=0; i<n; i++ ) {

		struct sym_op *r = get_item(newi, i);
		if ( find_item(items, r->h, r->k, r->l) ) continue;

		add_item_with_op(items, r->h, r->k, r->l, r->op);

	}
}


SymOpList *intersection_ops(SymOpList *i1, SymOpList *i2)
{
	int n, i;
	SymOpList *res = new_items();

	n = num_items(i1);
	for ( i=0; i<n; i++ ) {

		struct sym_op *r = get_item(i1, i);
		if ( find_item(i2, r->h, r->k, r->l) ) {
			add_item_with_op(res, r->h, r->k, r->l, r->op);
		}

	}

	return res;
}


/* Check if a reflection is in the asymmetric unit cell */
static int check_cond(signed int h, signed int k, signed int l, const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 )
		return ( 1 );
	if ( strcmp(sym, "-1") == 0 )
		return ( (l>0)
		      || ( (l==0) && (k>0) )
		      || ( (l==0) && (k==0) && (h>=0) ) );

	/* Orthorhombic */
	if ( strcmp(sym, "mmm") == 0 )
		return ( (h>=0) && (k>=0) && (l>=0) );
	if ( strcmp(sym, "222") == 0 )
		return ( (h>=0) && (k>0) )
		    || ( (h>0) && (k>=0) )
		    || ( (h>=0) && (k==0) && (l>=0) )
		    || ( (h==0) && (k>=0) && (l>=0) );

	/* Tetragonal */
	if ( strcmp(sym, "4/mmm") == 0 )
		return ( (((h>0) && (k>=0)) || ((h==0) && (k==0))) && (l>=0)
		         && (h>=k) );  /* Like 6/mmm */
	if ( strcmp(sym, "422") == 0 )
		return ( (((h>0) && (k>=0)) || ((h==0) && (k==0)))
		         && (h>=k) );
	if ( strcmp(sym, "4/m") == 0 )
		return ( (((h>0) && (k>=0)) || ((h==0) && (k==0))) && (l>=0) );
	if ( strcmp(sym, "4") == 0 )
		return ( ((h>0) && (k>=0)) || ((h==0) && (k==0)) );

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 )
		return ( ((h>0) && (k>=0)) || ((h==0) && (k==0)) );
	if ( strcmp(sym, "6/m") == 0 )
		return ( (((h>0) && (k>=0)) || ((h==0) && (k==0))) && (l>=0) );
	if ( strcmp(sym, "6/mmm") == 0 )
		return ( (((h>0) && (k>=0)) || ((h==0) && (k==0))) && (l>=0)
		         && (h>=k) );

	/* TODO: Add more groups here */

	return 1;
}


/* Macros for checking the above conditions and returning if satisfied */
#define CHECK_COND(h, k, l, sym)                   \
	if ( check_cond((h), (k), (l), (sym)) ) {  \
		*hp = (h);  *kp = (k);  *lp = (l); \
		return;                            \
	}


int num_general_equivs(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 1;
	if ( strcmp(sym, "-1") == 0 ) return 2;

	/* Orthorhombic */
	if ( strcmp(sym, "222") == 0 ) return 4;
	if ( strcmp(sym, "mmm") == 0 ) return 8;

	/* Tetragonal */
	if ( strcmp(sym, "4") == 0 ) return 4;
	if ( strcmp(sym, "4/m") == 0 ) return 8;
	if ( strcmp(sym, "422") == 0 ) return 8;
	if ( strcmp(sym, "4/mmm") == 0 ) return 16;

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return 6;
	if ( strcmp(sym, "6/m") == 0 ) return 12;
	if ( strcmp(sym, "6/mmm") == 0 ) return 24;

	/* TODO: Add more groups here */

	return 1;
}


void get_general_equiv(signed int h, signed int k, signed int l,
                       signed int *he, signed int *ke, signed int *le,
                       const char *sym, int idx)
{
	signed int i = -h-k;

	/* The returned indices when idx=0 *must* be the same as the input.
	 * After that, the order does not matter. */

	if ( strcmp(sym, "1") == 0 ) {
		*he = h;   *ke = k;   *le = l;  return;
	}

	if ( strcmp(sym, "-1") == 0 ) {
		switch ( idx ) {
		case 0 : *he = h;   *ke = k;   *le = l;   return;
		case 1 : *he = -h;  *ke = -k;  *le = -l;  return;
		}
	}

	if ( strcmp(sym, "222") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = -h;  *ke = -k;  *le = l;  return;
		case  2 : *he = -h;  *ke = k;  *le = -l;  return;
		case  3 : *he = h;  *ke = -k;  *le = -l;  return;
		}
	}

	if ( strcmp(sym, "mmm") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;  *ke = k;  *le = l;  return;
		case  1 : *he = -h; *ke = -k; *le = l;  return;
		case  2 : *he = -h; *ke = k;  *le = -l; return;
		case  3 : *he = h;  *ke = -k; *le = -l; return;
		case  4 : *he = -h; *ke = -k; *le = -l; return;
		case  5 : *he = h;  *ke = k;  *le = -l; return;
		case  6 : *he = h;  *ke = -k; *le = l;  return;
		case  7 : *he = -h; *ke = k;  *le = l;  return;
		}
	}

	if ( strcmp(sym, "4") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = -h;  *ke = -k;  *le = l;  return;
		case  2 : *he = -k;  *ke = h;   *le = l;  return;
		case  3 : *he = k;   *ke = -h;  *le = l;  return;
		}
	}

	if ( strcmp(sym, "4/m") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = -h;  *ke = -k;  *le = l;  return;
		case  2 : *he = -k;  *ke = h;   *le = l;  return;
		case  3 : *he = k;   *ke = -h;  *le = l;  return;
		case  4 : *he = -h;  *ke = -k;  *le = -l; return;
		case  5 : *he = h;   *ke = k;   *le = -l; return;
		case  6 : *he = k;   *ke = -h;  *le = -l; return;
		case  7 : *he = -k;  *ke = h;   *le = -l; return;
		}
	}

	if ( strcmp(sym, "422") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = -h;  *ke = -k;  *le = l;  return;
		case  2 : *he = -k;  *ke = h;   *le = l;  return;
		case  3 : *he = k;   *ke = -h;  *le = l;  return;
		case  4 : *he = -h;  *ke = k;   *le = -l; return;
		case  5 : *he = h;   *ke = -k;  *le = -l; return;
		case  6 : *he = k;   *ke = h;   *le = -l; return;
		case  7 : *he = -k;  *ke = -h;  *le = -l; return;
		}
	}

	if ( strcmp(sym, "4/mmm") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = -h;  *ke = -k;  *le = l;  return;
		case  2 : *he = -k;  *ke = h;   *le = l;  return;
		case  3 : *he = k;   *ke = -h;  *le = l;  return;
		case  4 : *he = -h;  *ke = k;   *le = -l; return;
		case  5 : *he = h;   *ke = -k;  *le = -l; return;
		case  6 : *he = k;   *ke = h;   *le = -l; return;
		case  7 : *he = -k;  *ke = -h;  *le = -l; return;
		case  8 : *he = -h;  *ke = -k;  *le = -l; return;
		case  9 : *he = h;   *ke = k;   *le = -l; return;
		case 10 : *he = k;   *ke = -h;  *le = -l; return;
		case 11 : *he = -k;  *ke = h;   *le = -l; return;
		case 12 : *he = h;   *ke = -k;  *le = l;  return;
		case 13 : *he = -h;  *ke = k;   *le = l;  return;
		case 14 : *he = -k;  *ke = -h;  *le = l;  return;
		case 15 : *he = k;   *ke = h;   *le = l;  return;
		}
	}

	if ( strcmp(sym, "6") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;  return;
		case  1 : *he = i;   *ke = h;   *le = l;  return;
		case  2 : *he = k;   *ke = i;   *le = l;  return;
		case  3 : *he = -h;  *ke = -k;  *le = l;  return;
		case  4 : *he = -i;  *ke = -h;  *le = l;  return;
		case  5 : *he = -k;  *ke = -i;  *le = l;  return;
		}
	}

	if ( strcmp(sym, "6/m") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;   return;
		case  1 : *he = i;   *ke = h;   *le = l;   return;
		case  2 : *he = k;   *ke = i;   *le = l;   return;
		case  3 : *he = -h;  *ke = -k;  *le = l;   return;
		case  4 : *he = -i;  *ke = -h;  *le = l;   return;
		case  5 : *he = -k;  *ke = -i;  *le = l;   return;
		case  6 : *he = h;   *ke = k;   *le = -l;  return;
		case  7 : *he = i;   *ke = h;   *le = -l;  return;
		case  8 : *he = k;   *ke = i;   *le = -l;  return;
		case  9 : *he = -h;  *ke = -k;  *le = -l;  return;
		case 10 : *he = -i;  *ke = -h;  *le = -l;  return;
		case 11 : *he = -k;  *ke = -i;  *le = -l;  return;
		}
	}

	if ( strcmp(sym, "6/mmm") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;   *ke = k;   *le = l;   return;
		case  1 : *he = i;   *ke = h;   *le = l;   return;
		case  2 : *he = k;   *ke = i;   *le = l;   return;
		case  3 : *he = -h;  *ke = -k;  *le = l;   return;
		case  4 : *he = -i;  *ke = -h;  *le = l;   return;
		case  5 : *he = -k;  *ke = -i;  *le = l;   return;
		case  6 : *he = k;   *ke = h;   *le = -l;  return;
		case  7 : *he = h;   *ke = i;   *le = -l;  return;
		case  8 : *he = i;   *ke = k;   *le = -l;  return;
		case  9 : *he = -k;  *ke = -h;  *le = -l;  return;
		case 10 : *he = -h;  *ke = -i;  *le = -l;  return;
		case 11 : *he = -i;  *ke = -k;  *le = -l;  return;
		case 12 : *he = -h;  *ke = -k;  *le = -l;  return;
		case 13 : *he = -i;  *ke = -h;  *le = -l;  return;
		case 14 : *he = -k;  *ke = -i;  *le = -l;  return;
		case 15 : *he = h;   *ke = k;   *le = -l;  return;
		case 16 : *he = i;   *ke = h;   *le = -l;  return;
		case 17 : *he = k;   *ke = i;   *le = -l;  return;
		case 18 : *he = -k;  *ke = -h;  *le = l;   return;
		case 19 : *he = -h;  *ke = -i;  *le = l;   return;
		case 20 : *he = -i;  *ke = -k;  *le = l;   return;
		case 21 : *he = k;   *ke = h;   *le = l;   return;
		case 22 : *he = h;   *ke = i;   *le = l;   return;
		case 23 : *he = i;   *ke = k;   *le = l;   return;
		}
	}

	/* TODO: Add more groups here */

	ERROR("Unrecognised symmetry '%s'\n", sym);
	abort();
}


/* Given a reflection and a point group, this returns (by reference) the indices
 * of the "idx"th equivalent reflection, taking special positions into account.
 * It returns "idx" if successful.  Otherwise, it returns the number of
 * equivalents for the particular reflection (taking special positions into
 * account).  Therefore, set idx=-1 to get the number of equivalents. */
static int special_position(signed int hs, signed int ks, signed int ls,
                            signed int *hp, signed int *kp, signed int *lp,
                            const char *sym, signed int idx)
{
	int n_general;
	int i;
	SymOpList *equivs;
	int n_equivs = 0;

	if ( idx == 0 ) {
		/* Index zero is always the original reflection */
		*hp = hs;  *kp = ks;  *lp = ls;
		return 0;
	}

	equivs = new_items();
	n_general = num_general_equivs(sym);

	for ( i=0; i<n_general; i++ ) {

		signed int h, k, l;

		/* Get equivalent according to the holohedral group */
		get_general_equiv(hs, ks, ls, &h, &k, &l, sym, i);

		/* Already got this one? */
		if ( find_item(equivs, h, k, l) ) continue;

		if ( n_equivs == idx ) {
			*hp = h;
			*kp = k;
			*lp = l;
			delete_items(equivs);
			return n_equivs;
		}
		add_item(equivs, h, k, l);
		n_equivs++;

	}

	delete_items(equivs);
	return n_equivs;
}


void get_equiv(signed int h, signed int k, signed int l,
               signed int *he, signed int *ke, signed int *le,
               const char *sym, int idx)
{
	special_position(h, k, l, he, ke, le, sym, idx);
}


int num_equivs(signed int h, signed int k, signed int l, const char *sym)
{
	return special_position(h, k, l, NULL, NULL, NULL, sym, -1);
}


void get_asymm(signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp,
               const char *sym)
{
	int nequiv = num_equivs(h, k, l, sym);
	int p;

	SYM_DEBUG("------ %i %i %i\n", h, k, l);
	for ( p=0; p<nequiv; p++ ) {
		signed int he, ke, le;
		get_equiv(h, k, l, &he, &ke, &le, sym, p);
		SYM_DEBUG("%i : %i %i %i\n", p, he, ke, le);
		CHECK_COND(he, ke, le, sym);
	}

	/* Should never reach here */
	ERROR("No match found in %s for %i %i %i\n", sym, h, k, l);
	abort();
}


/* This is kind of like a "numerical" left coset decomposition.
 * Given a reflection index and a point group, it returns the list of twinning
 * possibilities.
 *
 * To count the number of possibilities, use num_items() on the result.
 */
static SymOpList *coset_decomp(signed int hs, signed int ks, signed int ls,
                               const char *holo, const char *mero)
{
	int n_mero, n_holo;
	int i;
	signed int h, k, l;
	SymOpList *twins = new_items();

	/* Start by putting the given reflection into the asymmetric cell
	 * for its (probably merohedral) point group. */
	get_asymm(hs, ks, ls, &h, &k, &l, mero);

	/* How many equivalents in the holohedral point group are not
	 * equivalents according to the (possibly) merohedral group? */
	n_holo = num_general_equivs(holo);
	n_mero = num_general_equivs(mero);

	for ( i=0; i<n_holo; i++ ) {

		signed int h_holo, k_holo, l_holo;
		signed int hs_holo, ks_holo, ls_holo;

		/* Get equivalent according to the holohedral group */
		get_general_equiv(h, k, l, &hs_holo, &ks_holo, &ls_holo,
		                  holo, i);

		/* Put it into the asymmetric cell for the merohedral group */
		get_asymm(hs_holo, ks_holo, ls_holo,
		          &h_holo, &k_holo, &l_holo, mero);

		/* Already got this one?
		 * Note: The list "twins" starts empty, so the first iteration
		 * (i=0) will add the original reflection to the list along with
		 * the identity operation. */
		if ( find_item(twins, h_holo, k_holo, l_holo) ) continue;

		add_item_with_op(twins, h_holo, k_holo, l_holo, i);

	}

	return twins;
}


/* Work out the twinning possibilities for this pattern.
 * To use the result, call get_general_equiv() on each reflection using
 * the holohedral point group (use get_holohedral() for this), and for "idx"
 * give each "op" field from the list returned by this function. */
SymOpList *get_twins(const char *holo, const char *mero)
{
	int i;
	SymOpList *ops = new_items();
	int expected, actual;
	SymOpList *items;

	/* Run the coset decomposition for every reflection in the "pattern",
	 * and see which gives the highest number of possibilities.  This
	 * approach takes into account that a pattern consisting entirely of
	 * special reflections might have fewer twin possibilities. */
	for ( i=0; i<num_items(items); i++ ) {

		signed int h, k, l;
		struct sym_op *item;
		SymOpList *new_ops;

		item = get_item(items, i);

		h = item->h;
		k = item->k;
		l = item->l;

		new_ops = coset_decomp(h, k, l, holo, mero);
		union_op_items(ops, new_ops);
		delete_items(new_ops);

	}

	/* Idiot check */
	actual = num_items(ops);
	expected = num_general_equivs(holo) / num_general_equivs(mero);
	if ( actual != expected ) {
		ERROR("Found %i twin possibilities, but expected %i.\n",
		       actual, expected);
		ERROR("I couldn't find the number of twin laws that I expected."
		      " This is an internal error, and shouldn't happen. "
		      "Sorry.\n");
		abort();
	}

	return ops;
}


int find_unique_equiv(SymOpList *items, signed int h, signed int k,
                      signed int l, const char *mero, signed int *hu,
                      signed int *ku, signed int *lu)
{
	int i;
	int found = 0;

	for ( i=0; i<num_equivs(h, k, l, mero); i++ ) {

		signed int he, ke, le;
		int f;
		get_equiv(h, k, l, &he, &ke, &le, mero, i);
		f = find_item(items, he, ke, le);

		/* There must only be one equivalent.  If there are more, it
		 * indicates that the user lied about the input symmetry.
		 * This situation should have been checked for earlier by
		 * calling check_symmetry() with 'items' and 'mero'. */

		if ( f && !found ) {
			*hu = he;  *ku = ke;  *lu = le;
			return 1;
		}

	}

	return 0;
}


/* Returns true if the point group is 23, m-3, 432, -43m or m-3m
 * (i.e. T, Th, O, Td or Oh) */
int is_polyhedral(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 0;
	if ( strcmp(sym, "-1") == 0 ) return 0;

	/* Orthorhombic */
	if ( strcmp(sym, "222") == 0 ) return 0;
	if ( strcmp(sym, "mmm") == 0 ) return 0;

	/* Tetragonal */
	if ( strcmp(sym, "4") == 0 ) return 0;
	if ( strcmp(sym, "4/m") == 0 ) return 0;
	if ( strcmp(sym, "422") == 0 ) return 0;
	if ( strcmp(sym, "4/mmm") == 0 ) return 0;

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return 0;
	if ( strcmp(sym, "6/m") == 0 ) return 0;
	if ( strcmp(sym, "6/mmm") == 0 ) return 0;

	/* TODO: Add more groups here */

	ERROR("Don't know if '%s' is polyhedral or not.\n", sym);
	abort();
}


/* Returns the order of the characteristic axis of proper or improper rotation
 * for the point group.  This should be the rotation about the c axis. */
int rotational_order(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 1;
	if ( strcmp(sym, "-1") == 0 ) return 2;

	/* Orthorhombic */
	if ( strcmp(sym, "222") == 0 ) return 2;
	if ( strcmp(sym, "mmm") == 0 ) return 2;

	/* Tetragonal */
	if ( strcmp(sym, "4") == 0 ) return 4;
	if ( strcmp(sym, "4/m") == 0 ) return 4;
	if ( strcmp(sym, "422") == 0 ) return 4;
	if ( strcmp(sym, "4/mmm") == 0 ) return 4;

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return 6;
	if ( strcmp(sym, "6/m") == 0 ) return 6;
	if ( strcmp(sym, "6/mmm") == 0 ) return 6;

	/* TODO: Add more groups here */

	ERROR("Couldn't find rotational order for '%s'.\n", sym);
	abort();
}


int has_perpendicular_mirror(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 0;
	if ( strcmp(sym, "-1") == 0 ) return 0;

	/* Orthorhombic */
	if ( strcmp(sym, "222") == 0 ) return 0;
	if ( strcmp(sym, "mmm") == 0 ) return 1;

	/* Tetragonal */
	if ( strcmp(sym, "4") == 0 ) return 0;
	if ( strcmp(sym, "4/m") == 0 ) return 1;
	if ( strcmp(sym, "422") == 0 ) return 0;
	if ( strcmp(sym, "4/mmm") == 0 ) return 1;

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return 0;
	if ( strcmp(sym, "6/m") == 0 ) return 1;
	if ( strcmp(sym, "6/mmm") == 0 ) return 1;

	/* TODO: Add more groups here */

	ERROR("Couldn't find mirror definition for '%s'.\n", sym);
	abort();
}


int has_bisecting_mirror_or_diad(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 0;
	if ( strcmp(sym, "-1") == 0 ) return 0;

	/* Orthorhombic */
	if ( strcmp(sym, "222") == 0 ) return 1;
	if ( strcmp(sym, "mmm") == 0 ) return 1;

	/* Tetragonal */
	if ( strcmp(sym, "4") == 0 ) return 0;
	if ( strcmp(sym, "4/m") == 0 ) return 0;
	if ( strcmp(sym, "422") == 0 ) return 0;
	if ( strcmp(sym, "4/mmm") == 0 ) return 1;

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return 0;
	if ( strcmp(sym, "6/m") == 0 ) return 1;
	if ( strcmp(sym, "6/mmm") == 0 ) return 1;

	/* TODO: Add more groups here */

	ERROR("Couldn't find mirror definition for '%s'.\n", sym);
	abort();
}


int check_symmetry(SymOpList *items, const char *sym)
{
	int i;
	unsigned char *flags;

	flags = new_list_flag();
	for ( i=0; i<num_items(items); i++ ) {
		struct sym_op *it = get_item(items, i);
		set_flag(flags, it->h, it->k, it->l, 1);
	}

	for ( i=0; i<num_items(items); i++ ) {

		int j;
		struct sym_op *it = get_item(items, i);
		int found = 0;

		for ( j=0; j<num_equivs(it->h, it->k, it->l, sym); j++ ) {

			signed int he, ke, le;
			get_equiv(it->h, it->k, it->l, &he, &ke, &le, sym, j);

			if ( abs(he) > INDMAX ) continue;
			if ( abs(le) > INDMAX ) continue;
			if ( abs(ke) > INDMAX ) continue;

			found += lookup_flag(flags, he, ke, le);

		}

		if ( found > 1 ) {
			free(flags);
			return 1;  /* Symmetry is wrong! */
		}

	}

	free(flags);

	return 0;
}
