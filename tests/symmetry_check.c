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


static void check_nequiv(const char *pg, int answer, int *fail)
{
	SymOpList *sym;
	int n;

	//STATUS("**************************************** Testing '%s'\n", pg);

	sym = get_pointgroup(pg);
	n = num_equivs(sym, NULL);

	if ( n != answer ) {
		ERROR("Number of equivalents in '%s' is %i (not %i)\n",
		      pg, n, answer);
		*fail = 1;
	}

	free_symoplist(sym);
}


int main(int argc, char *argv[])
{
	int fail = 0;

	check_nequiv(  "1",  1, &fail);
	check_nequiv( "-1",  2, &fail);
	check_nequiv(  "2",  2, &fail);
	check_nequiv(  "m",  2, &fail);
	check_nequiv("2/m",  4, &fail);

	return fail;
}
