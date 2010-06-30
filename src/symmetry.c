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


/* Conditions for a reflection to be in the asymmetric unit cell */
#define COND_1(h, k, l) (1)
#define COND_6(h, k, l) ( (h>=0) && (k>=0) )
#define COND_6M(h, k, l) ( (h>=0) && (k>=0) && (l>=0) )
#define COND_6MMM(h, k, l) ( (h>=0) && (k>=0) && (l>=0) \
                          && (h>=k) )
/* TODO: Add more groups here */


/* Macros for checking the above conditions and returning if satisfied */
#define CHECK_COND(h, k, l, cond)                  \
	if ( COND_##cond((h), (k), (l)) ) {        \
		*hp = (h);  *kp = (k);  *lp = (l); \
		return;                            \
	}


/* Abort macro if no match found */
#define SYM_ABORT                                             \
	ERROR("No match in %s for %i %i %i\n", sym, h, k, l); \
	abort();


/* FIXME: Should take into account special indices
 * e.g. l==0 has fewer equivalent reflections */
int num_equivs(signed int h, signed int k, signed int l, const char *sym)
{
	if ( strcmp(sym, "6/mmm") == 0 ) return 24;
	if ( strcmp(sym, "6") == 0 ) return 6;
	if ( strcmp(sym, "6/m") == 0 ) return 12;
	/* TODO: Add more groups here */

	return 1;
}


void get_equiv(signed int h, signed int k, signed int l,
               signed int *he, signed int *ke, signed int *le,
               const char *sym, int idx)
{
	if ( strcmp(sym, "1") == 0 ) {
		*he = h;   *ke = k;   *le = l;  return;
	}

	if ( strcmp(sym, "6") == 0 ) {
		signed int i = -h-k;
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
		signed int i = -h-k;
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
		signed int i = -h-k;
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

	/* Fallback for unrecognised groups */
	*he = h;  *ke = k;  *le = l;
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
		if ( strcmp(sym, "1") == 0 ) CHECK_COND(he, ke, le, 1);
		if ( strcmp(sym, "6") == 0 ) CHECK_COND(he, ke, le, 6);
		if ( strcmp(sym, "6/m") == 0 ) CHECK_COND(he, ke, le, 6M);
		if ( strcmp(sym, "6/mmm") == 0 ) CHECK_COND(he, ke, le, 6MMM);
	}

	SYM_ABORT;  /* Should never reach here */

	ERROR("Unknown point group '%s'\n", sym);
	abort();
}
