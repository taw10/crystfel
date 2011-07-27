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
#include <stdarg.h>

#include "../src/symmetry.h"
#include "../src/utils.h"


static void find_all_ambiguities(const char *first, ...)
{
	va_list vp;
	int i;
	const char *arg;
	SymOpList *test[32];
	int n = 0;

	test[n++] = get_pointgroup(first);

	va_start(vp, first);


	do {

		arg = va_arg(vp, const char *);
		if ( arg != NULL ) {
			test[n++] = get_pointgroup(arg);
		}

	} while ( arg != NULL );

	for ( i=0; i<n; i++ ) {

		SymOpList *holo;
		int j;

		holo = test[i];
		STATUS("%7s :", symmetry_name(holo));
		for ( j=0; j<n; j++ ) {

			SymOpList *twins;

			if ( i == j ) continue;

			if ( !is_subgroup(holo, test[j]) ) continue;

			twins = get_ambiguities(holo, test[j]);
			if ( twins == NULL ) continue;

			if ( num_equivs(twins, NULL) == 1 ) {
				free_symoplist(twins);
				continue;
			}

			STATUS(" %s(%i)", symmetry_name(test[j]),
			                   num_equivs(twins, NULL));

		}
		STATUS("\n");

	}

	STATUS("\n");
}


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


static void check_subgroup(const char *ssource, const char *starget,
                           int should_be_subgroup, int should_be_ambiguity,
                           int n_exp, int *fail)
{
	SymOpList *source;
	SymOpList *target;
	int sub;

	source = get_pointgroup(ssource);
	target = get_pointgroup(starget);
	if ( (source == NULL) || (target == NULL) ) {
		*fail = 1;
		return;
	}

	sub = is_subgroup(source, target);
	if ( sub != should_be_subgroup ) {
		ERROR("'%s' should%s be a subgroup of '%s', but is%s.\n",
		      starget, maybenot(should_be_subgroup),
		      ssource, maybenot(sub));
		*fail = 1;
		return;
	}

	if ( should_be_subgroup ) {

		SymOpList *twins;
		int nf;
		int amb = 1;

		twins = get_ambiguities(source, target);
		if ( twins == NULL ) amb = 0;
		if ( amb != should_be_ambiguity ) {
			ERROR("'%s' should%s be a rotational subgroup of '%s'"
			      " but is%s.\n",
			      starget, maybenot(should_be_ambiguity),
			      ssource, maybenot(amb));
			*fail = 1;
			return;
		}

		if ( amb ) {

			describe_symmetry(twins);
			nf = num_equivs(twins, NULL);
			if ( nf != n_exp ) {
				ERROR("Expected %i operations, found %i\n",
				      n_exp, nf);
				*fail = 1;
				return;
			}

		} else {

			STATUS("%15s : subgroup of %s, but no ambiguity without"
			       " inversion or mirroring\n", starget, ssource);


		}

	} else {

		STATUS("%15s : not a subgroup of %s\n", starget, ssource);

	}

	free_symoplist(target);
	free_symoplist(source);
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
	check_pg_props( "-43m",  24,  0, &fail);
	check_pg_props( "m-3m",  48,  1, &fail);
	STATUS("\n");

	/* Check some easy subgroups */
	check_subgroup("2/m",    "m",    1, 1,  2, &fail);
	check_subgroup("mmm",    "mm2",  1, 1,  2, &fail);
	check_subgroup("-4m2",   "-4",   1, 1,  2, &fail);
	check_subgroup("-42m",   "-4",   1, 1,  2, &fail);
	check_subgroup("-3m1_H", "-3_H", 1, 1,  2, &fail);
	check_subgroup("-31m_H", "-3_H", 1, 1,  2, &fail);
	check_subgroup("m-3m",   "-43m", 1, 1,  2, &fail);
	check_subgroup("m-3m",   "m-3",  1, 1,  2, &fail);
	check_subgroup("432",    "23",   1, 1,  2, &fail);
	check_subgroup("6/m",    "-3_H", 1, 1,  2, &fail);
	check_subgroup("4/m",    "-4",   1, 1,  2, &fail);

	/* Tetartohedral */
	check_subgroup("6/mmm",  "-3_H", 1, 1,  4, &fail);

	/* Check some things that are valid subgroups, but no valid ambiguities
	 * exist because inversions and mirrors are not allowed */
	check_subgroup("-1",     "1",    1, 0, -1, &fail);
	check_subgroup("4/mmm",  "4",    1, 0, -1, &fail);
	check_subgroup("m-3m",   "23",   1, 0, -1, &fail);

	/* Check some invalid combinations */
	check_subgroup("432",    "-43m", 0, 0, -1, &fail);
	check_subgroup("432",    "m-3",  0, 0, -1, &fail);

	/* Derive all merohedral ambiguities */
	STATUS("\nMerohedral ambiguities:\n\n");
	find_all_ambiguities("1", "-1", NULL);
	find_all_ambiguities("2", "m", "2/m", NULL);
	find_all_ambiguities("mm2", "mmm", "222", NULL);
	find_all_ambiguities("4", "-4", "-42m", "-4m2", "4mm",
	                     "4/m", "422", "4/mmm", NULL);
	find_all_ambiguities("3_R", "32_R", "-3_R", "3m_R", "-3m_R", NULL);
	find_all_ambiguities("6", "3_H", "312_H", "321_H", "622", "-3_H",
	                     "3m1_H", "-6", "31m_H", "-3m1_H", "-6m2",
	                     "-62m", "-31m_H", "6/mmm", "6/m", "6mm", NULL);
	find_all_ambiguities("23", "432", "-43m", "m-3", "m-3m", NULL);


	STATUS("\nPseudo-merohedral ambiguities:\n\n");
	find_all_ambiguities("1", "-1", "2", "m", "2/m", "mm2", "mmm", "222",
	                     "4", "-4", "-42m", "-4m2", "4mm", "4/m", "422",
	                     "4/mmm", "23", "432", "-43m", "m-3", "m-3m", NULL);

	find_all_ambiguities("1", "-1", "3_R", "32_R", "-3_R", "3m_R", "-3m_R",
	                     NULL);

	find_all_ambiguities("1", "-1", "2", "m", "2/m", "6", "3_H", "312_H",
	                     "321_H", "622", "-3_H", "3m1_H", "-6", "31m_H",
	                     "-3m1_H", "-6m2", "-62m", "-31m_H", "6/mmm", "6/m",
	                     "6mm", NULL);

	/* Check some pseudo-meroheral subgroups */
	check_subgroup("3_R",  "1", 1, 1,  3, &fail);
	check_subgroup("-3_R",  "-1", 1, 1,  3, &fail);
	check_subgroup("6",  "2", 1, 1,  3, &fail);

	return fail;
}
