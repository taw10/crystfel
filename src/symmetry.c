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


int num_equivs(signed int h, signed int k, signed int l, const char *sym)
{
	/* 000 is always unique */
	if ( (h==0) && (k==0) && (l==0) ) return 1;

	if ( strcmp(sym, "1") == 0 ) return 1;

	if ( strcmp(sym, "-1") == 0 ) return 2;

	if ( strcmp(sym, "6") == 0 ) {
		if ( (h==0) && (k==0) ) return 1;  /* a */
		return 6;  /* b */
	}

	if ( strcmp(sym, "6/m") == 0 ) {
		if ( (h==0) && (k==0) ) return 2;  /* a */
		if ( l == 0 ) return 6;  /* b */
		return 12;  /* c */
	}

	if ( strcmp(sym, "6/mmm") == 0 ) {
		if ( (h==0) && (k==0) ) return 2;  /* a */
		if ( (l==0) && (h==k) ) return 6;  /* b */
		if ( (l==0) && (k==0) ) return 6;  /* c */
		if ( h == k ) return 12;  /* d */
		if ( k == 0 ) return 12;  /* e */
		if ( l == 0 ) return 12;  /* f */
		return 24;  /* g */
	}

	/* TODO: Add more groups here */

	return 1;
}


void get_equiv(signed int h, signed int k, signed int l,
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
		/* a */
		if ( (h==0) && (k==0) ) {
			switch ( idx ) {
			case  0 : *he = 0;   *ke = 0;   *le = l;   return;
			}
		}
		/* b */
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
		/* a */
		if ( (h==0) && (k==0) ) {
			switch ( idx ) {
			case  0 : *he = 0;   *ke = 0;   *le = l;   return;
			case  1 : *he = 0;   *ke = 0;   *le = -l;  return;
			}
		}
		/* b-c */
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
		/* a */
		if ( (h==0) && (k==0) ) {
			switch ( idx ) {
			case  0 : *he = 0;   *ke = 0;   *le = l;   return;
			case  1 : *he = 0;   *ke = 0;   *le = -l;  return;
			}
		}
		/* b-g */
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


static const char *get_holohedral(const char *sym)
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
 * Given a reflection index and a point group, it returns the "idx-th"
 * twinning possibility for the reflection.  It returns "idx" if
 * successful.  To just count the number of possibilities, set idx=-1.
 *
 * The sequence of operators producing each possibility is guaranteed to
 * be the same for any choice of indices given the same point group. */
static int coset_decomp(signed int hs, signed int ks, signed int ls,
                        signed int *hp, signed int *kp, signed int *lp,
                        const char *mero, signed int idx)
{
	const char *holo = get_holohedral(mero);
	int n_mero, n_holo;
	int i;
	signed int n_twins = 1;
	signed int h, k, l;
	ReflItemList *twins;

	twins = new_items();

	if ( idx == 0 ) {
		/* Twin index zero is always the original orientation */
		*hp = hs;  *kp = ks;  *lp = ls;
		return 0;
	}

	get_asymm(hs, ks, ls, &h, &k, &l, mero);

	/* How many equivalents in the holohedral point group are not
	 * equivalents according to the (possibly) merohedral group? */
	n_holo = num_equivs(h, k, l, holo);
	n_mero = num_equivs(h, k, l, mero);

	for ( i=0; i<n_holo; i++ ) {

		signed int h_holo, k_holo, l_holo;
		signed int hs_holo, ks_holo, ls_holo;

		/* Get equivalent according to the holohedral group */
		get_equiv(h, k, l, &hs_holo, &ks_holo, &ls_holo, holo, i);

		/* Put it into the asymmetric cell for the merohedral group */
		get_asymm(hs_holo, ks_holo, ls_holo,
		          &h_holo, &k_holo, &l_holo, mero);

		/* Is this the same reflection as we started with?
		 * If not, this reflection is 'equivalent by twinning' */
		if ( (h_holo != h) || (k_holo != k) || (l_holo != l) ) {

			if ( find_item(twins, h_holo, k_holo, l_holo) )
				continue;

			if ( n_twins == idx ) {
				*hp = h_holo;
				*kp = k_holo;
				*lp = l_holo;
				delete_items(twins);
				return n_twins;
			}
			add_item(twins, h_holo, k_holo, l_holo);
			n_twins++;
		}

	}

	delete_items(twins);
	return n_twins;
}


/* Get the number of twinned "equivalents" for this reflection */
int num_twins(signed int h, signed int k, signed int l, const char *sym)
{
	return coset_decomp(h, k, l, NULL, NULL, NULL, sym, -1);
}


void get_twins(signed int h, signed int k, signed int l,
              signed int *hp, signed int *kp, signed int *lp,
              const char *sym, int idx)
{
	if ( coset_decomp(h, k, l, hp, kp, lp, sym, idx) != idx ) {
		ERROR("Failed coset decomposition for %i %i %i in %s\n",
		      h, k, l, sym);
		abort();
	}
}
