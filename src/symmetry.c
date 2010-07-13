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
#include <math.h>

#include "utils.h"


#ifdef DEBUG
#define SYM_DEBUG STATUS
#else /* DEBUG */
#define SYM_DEBUG(...)
#endif /* DEBUG */


/* Check if a reflection is in the asymmetric unit cell */
static int check_cond(signed int h, signed int k, signed int l, const char *sym)
{
	if ( strcmp(sym, "1") == 0 )
		return ( 1 );
	if ( strcmp(sym, "-1") == 0 )
		return ( 1 );
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


static int num_general_equivs(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return 1;
	if ( strcmp(sym, "-1") == 0 ) return 2;

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

	if ( strcmp(sym, "1") == 0 ) {
		*he = h;   *ke = k;   *le = l;  return;
	}

	if ( strcmp(sym, "-1") == 0 ) {
		switch ( idx ) {
		case 0 : *he = h;   *ke = k;   *le = l;   return;
		case 1 : *he = -h;  *ke = -k;  *le = -l;  return;
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


static int special_position(signed int hs, signed int ks, signed int ls,
                            signed int *hp, signed int *kp, signed int *lp,
                            const char *sym, signed int idx)
{
	int n_general;
	int i;
	ReflItemList *equivs;
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


const char *get_holohedral(const char *sym)
{
	/* Triclinic */
	if ( strcmp(sym, "1") == 0 ) return "-1";
	if ( strcmp(sym, "1") == 0 ) return "-1";

	/* Hexagonal */
	if ( strcmp(sym, "6") == 0 ) return "6/m";
	if ( strcmp(sym, "6/m") == 0 ) return "6/mmm";
	if ( strcmp(sym, "6/mmm") == 0 ) return "6/mmm";

	/* TODO: Add more groups here */

	ERROR("Couldn't find holohedral point group for '%s'\n", sym);
	abort();
}


/* This is kind of like a "numerical" left coset decomposition.
 * Given a reflection index and a point group, it returns the list of twinning
 * possibilities.
 *
 * To count the number of possibilities, use num_items() on the result.
 */
static ReflItemList *coset_decomp(signed int hs, signed int ks, signed int ls,
                                  const char *mero)
{
	const char *holo = get_holohedral(mero);
	int n_mero, n_holo;
	int i;
	signed int h, k, l;
	ReflItemList *twins = new_items();

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
ReflItemList *get_twins(ReflItemList *items, const char *sym)
{
	int i;
	int n_twins = 1;
	ReflItemList *max_ops = NULL;

	/* Run the coset decomposition for every reflection in the "pattern",
	 * and see which gives the highest number of possibilities.  This
	 * approach takes into account that a pattern consisting entirely of
	 * special reflections might have fewer twin possibilities. */
	for ( i=0; i<num_items(items); i++ ) {

		signed int h, k, l;
		struct refl_item *item;
		ReflItemList *ops;

		item = get_item(items, i);

		h = item->h;
		k = item->k;
		l = item->l;

		ops = coset_decomp(h, k, l, sym);
		if ( num_items(ops) > n_twins ) {
			n_twins = num_items(ops);
			delete_items(max_ops);
			max_ops = ops;
		}

	}

	return max_ops;
}
