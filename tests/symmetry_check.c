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

	sym = get_pointgroup(pg);

	if ( sym == NULL ) {
		*fail = 1;
		return;
	}

	n = num_equivs(sym, NULL);

	describe_symmetry(sym);

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

	check_pg_props(    "1",   1,  0, &fail);
	check_pg_props(   "-1",   2,  1, &fail);
	STATUS("\n");

	check_pg_props(    "2",   2,  0, &fail);
	check_pg_props(    "m",   2,  0, &fail);
	check_pg_props(  "2/m",   4,  1, &fail);
	STATUS("\n");

	check_pg_props(  "222",   4,  0, &fail);
	check_pg_props(  "mm2",   4,  0, &fail);
	check_pg_props(  "mmm",   8,  1, &fail);
	STATUS("\n");

	check_pg_props(     "4",  4,  0, &fail);
	check_pg_props(    "-4",  4,  0, &fail);
	check_pg_props(   "4/m",  8,  1, &fail);
	check_pg_props(   "422",  8,  0, &fail);
	check_pg_props(   "4mm",  8,  0, &fail);
	check_pg_props(  "-42m",  8,  0, &fail);
	check_pg_props(  "-4m2",  8,  0, &fail);
	check_pg_props( "4/mmm", 16,  1, &fail);
	STATUS("\n");

	check_pg_props(  "3_R",   3,  0, &fail);
	check_pg_props( "-3_R",   6,  1, &fail);
	check_pg_props( "32_R",   6,  0, &fail);
	check_pg_props( "3m_R",   6,  0, &fail);
	check_pg_props("-3m_R",  12,  1, &fail);
	STATUS("\n");

	check_pg_props(   "3_H",  3,  0, &fail);
	check_pg_props(  "-3_H",  6,  1, &fail);
	check_pg_props( "321_H",  6,  0, &fail);
	check_pg_props( "312_H",  6,  0, &fail);
	check_pg_props( "3m1_H",  6,  0, &fail);
	check_pg_props( "31m_H",  6,  0, &fail);
	check_pg_props("-3m1_H", 12,  1, &fail);
	check_pg_props("-31m_H", 12,  1, &fail);
	STATUS("\n");

	check_pg_props(     "6",  6,  0, &fail);
	check_pg_props(    "-6",  6,  0, &fail);
	check_pg_props(   "6/m", 12,  1, &fail);
	check_pg_props(   "622", 12,  0, &fail);
	check_pg_props(   "6mm", 12,  0, &fail);
	check_pg_props(  "-6m2", 12,  0, &fail);
	check_pg_props(  "-62m", 12,  0, &fail);
	check_pg_props( "6/mmm", 24,  1, &fail);
	STATUS("\n");

	check_pg_props(   "23",  12,  0, &fail);
	check_pg_props(  "m-3",  24,  1, &fail);
	check_pg_props(  "432",  24,  0, &fail);
	check_pg_props( "-432",  24,  0, &fail);
	check_pg_props( "m-3m",  48,  1, &fail);
	STATUS("\n");

	return fail;
}
