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
#define COND_6MMM(h, k, i, l) ( (h>=0) && (k>=0) && (l>=0) && ((h>k)||((h==0)&&(k==0))) )


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

		CHECK_COND_FOURIDX(h, k, i, l, 6MMM);
		CHECK_COND_FOURIDX(h, i, k, l, 6MMM);
		CHECK_COND_FOURIDX(k, h, i, l, 6MMM);
		CHECK_COND_FOURIDX(k, i, h, l, 6MMM);
		CHECK_COND_FOURIDX(i, h, k, l, 6MMM);
		CHECK_COND_FOURIDX(i, k, h, l, 6MMM);

		CHECK_COND_FOURIDX(h, k, i, -l, 6MMM);
		CHECK_COND_FOURIDX(h, i, k, -l, 6MMM);
		CHECK_COND_FOURIDX(k, h, i, -l, 6MMM);
		CHECK_COND_FOURIDX(k, i, h, -l, 6MMM);
		CHECK_COND_FOURIDX(i, h, k, -l, 6MMM);
		CHECK_COND_FOURIDX(i, k, h, -l, 6MMM);

		CHECK_COND_FOURIDX(-h, -k, -i, l, 6MMM);
		CHECK_COND_FOURIDX(-h, -i, -k, l, 6MMM);
		CHECK_COND_FOURIDX(-k, -h, -i, l, 6MMM);
		CHECK_COND_FOURIDX(-k, -i, -h, l, 6MMM);
		CHECK_COND_FOURIDX(-i, -h, -k, l, 6MMM);
		CHECK_COND_FOURIDX(-i, -k, -h, l, 6MMM);

		CHECK_COND_FOURIDX(-h, -k, -i, -l, 6MMM);
		CHECK_COND_FOURIDX(-h, -i, -k, -l, 6MMM);
		CHECK_COND_FOURIDX(-k, -h, -i, -l, 6MMM);
		CHECK_COND_FOURIDX(-k, -i, -h, -l, 6MMM);
		CHECK_COND_FOURIDX(-i, -h, -k, -l, 6MMM);
		CHECK_COND_FOURIDX(-i, -k, -h, -l, 6MMM);

		SYM_ABORT;  /* Should never reach here */

	}

	ERROR("Unknown point group '%s'\n", sym);
	abort();
}
