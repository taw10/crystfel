/*
 * symmetry_check.c
 *
 * Check symmetry
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include "../src/symmetry.h"
#include "../src/utils.h"


static const char *maybenot(int v)
{
	if ( v ) {
		return "";
	} else {
		return " not";
	}
}


static void check_pg_props(const char *pg, int answer, int centro, int *fail)
{
	SymOpList *sym;
	int n, c;

	//STATUS("**************************************** Testing '%s'\n", pg);

	sym = get_pointgroup(pg);
	n = num_equivs(sym, NULL);

	if ( n != answer ) {
		ERROR("Number of equivalents in '%s' is %i (not %i).\n",
		      pg, n, answer);
		*fail = 1;
	}

	c = is_centrosymmetric(sym);
	if ( c != centro ) {
		ERROR("'%s' should%s be centrosymmetric, but is%s.\n",
		      pg, maybenot(centro), maybenot(c));
		*fail = 1;
	}

	free_symoplist(sym);
}


int main(int argc, char *argv[])
{
	int fail = 0;

	check_pg_props(  "1",  1, 0, &fail);
	check_pg_props( "-1",  2, 1, &fail);

	check_pg_props(  "2",  2, 0, &fail);
	check_pg_props(  "m",  2, 0, &fail);
	check_pg_props("2/m",  4, 1, &fail);

	check_pg_props("222",  4, 0, &fail);
	check_pg_props("mm2",  4, 0, &fail);
	check_pg_props("mmm",  8, 1, &fail);

	check_pg_props(  "4",  4, 0, &fail);

	return fail;
}
