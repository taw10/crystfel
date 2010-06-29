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


/* Conditions for a reflection to be in the asymmetric unit cell */
#define COND_1(h, k, l) (1)
#define COND_6MMM(h, k, i, l) ( (h>=0) && (k>=0) && (l>=0) \
                             && ((h>k)||((h==0)&&(k==0))) )
/* TODO: Add more groups here */

/* Macros for checking the above conditions and returning if satisfied */
#define CHECK_COND_FOURIDX(h, k, i, l, cond)       \
	if ( COND_##cond((h), (k), (i), (l)) ) {   \
		*hp = (h);  *kp = (k);  *lp = (l); \
		return;                            \
	}

#define CHECK_COND_THREEIDX(h, k, l, cond)          \
	if ( COND_##cond((h), (k), (l)) ) {         \
		*hp = (h);  *kp = (k);  *lp = (l);  \
		return;                             \
	}


/* Abort macro if no match found */
#define SYM_ABORT                                             \
	ERROR("No match in %s for %i %i %i\n", sym, h, k, l); \
	abort();


/* FIXME: Should take into account special indices
 * e.g. l==0 has fewer equivalent reflections */
static int num_equivs4(signed int h, signed int k, signed int i, signed int l,
                       const char *sym)
{
	if ( strcmp(sym, "6/mmm") == 0 ) return 24;
	/* TODO: Add more groups here */

	return 1;
}


static void get_equiv4(signed int h, signed int k, signed int i, signed int l,
                       signed int *he, signed int *ke, signed int *ie,
                       signed int *le, const char *sym, int idx)
{
	if ( strcmp(sym, "6/mmm") == 0 ) {
		switch ( idx ) {
		case  0 : *he = h;  *ke = k;  *ie = i;  *le = l;  return;
		case  1 : *he = h;  *ke = i;  *ie = k;  *le = l;  return;
		case  2 : *he = k;  *ke = h;  *ie = i;  *le = l;  return;
		case  3 : *he = k;  *ke = i;  *ie = h;  *le = l;  return;
		case  4 : *he = i;  *ke = h;  *ie = k;  *le = l;  return;
		case  5 : *he = i;  *ke = k;  *ie = h;  *le = l;  return;
		case  6 : *he = h;  *ke = k;  *ie = i;  *le = -l;  return;
		case  7 : *he = h;  *ke = i;  *ie = k;  *le = -l;  return;
		case  8 : *he = k;  *ke = h;  *ie = i;  *le = -l;  return;
		case  9 : *he = k;  *ke = i;  *ie = h;  *le = -l;  return;
		case 10 : *he = i;  *ke = h;  *ie = k;  *le = -l;  return;
		case 11 : *he = i;  *ke = k;  *ie = h;  *le = -l;  return;
		case 12 : *he = -h;  *ke = -k;  *ie = -i;  *le = l;  return;
		case 13 : *he = -h;  *ke = -i;  *ie = -k;  *le = l;  return;
		case 14 : *he = -k;  *ke = -h;  *ie = -i;  *le = l;  return;
		case 15 : *he = -k;  *ke = -i;  *ie = -h;  *le = l;  return;
		case 16 : *he = -i;  *ke = -h;  *ie = -k;  *le = l;  return;
		case 17 : *he = -i;  *ke = -k;  *ie = -h;  *le = l;  return;
		case 18 : *he = -h;  *ke = -k;  *ie = -i;  *le = -l;  return;
		case 19 : *he = -h;  *ke = -i;  *ie = -k;  *le = -l;  return;
		case 20 : *he = -k;  *ke = -h;  *ie = -i;  *le = -l;  return;
		case 21 : *he = -k;  *ke = -i;  *ie = -h;  *le = -l;  return;
		case 22 : *he = -i;  *ke = -h;  *ie = -k;  *le = -l;  return;
		case 23 : *he = -i;  *ke = -k;  *ie = -h;  *le = -l;  return;
		}
	}

	*he = h;  *ke = k;  *ie = i;  *le = l;
}


void get_asymm(signed int h, signed int k, signed int l,
               signed int *hp, signed int *kp, signed int *lp,
               const char *sym)
{
	if ( strcmp(sym, "1") == 0 ) {
		CHECK_COND_THREEIDX(h, k, l, 1);
		SYM_ABORT;
	}

	if ( strcmp(sym, "6/mmm") == 0 ) {

		const signed int i = h+k;
		int nequiv = num_equivs4(h, k, i, l, sym);
		int p;

		for ( p=0; p<nequiv; p++ ) {
			signed int he, ke, ie, le;
			get_equiv4(h, k, i, l, &he, &ke, &ie, &le, sym, p);
			CHECK_COND_FOURIDX(he, ke, ie, le, 6MMM);
		}

		SYM_ABORT;  /* Should never reach here */

	}

	/* TODO: Add more groups here */

	ERROR("Unknown point group '%s'\n", sym);
	abort();
}
