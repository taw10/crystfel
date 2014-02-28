/*
 * ambi_check.c
 *
 * Check indexing ambiguities
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2014 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <symmetry.h>
#include <utils.h>


static int is_nonmirror_subgroup(SymOpList *holo, SymOpList *mero)
{
	SymOpList *twins;
	int index;

	if ( !is_subgroup(holo, mero) ) return 0;

	twins = get_ambiguities(holo, mero);
	if ( twins == NULL ) return 0;

	index = num_equivs(twins, NULL);
	free_symoplist(twins);

	return index;
}


static int is_maximal_nonmirror_subgroup(SymOpList *holo, SymOpList *mero,
                                         SymOpList **list, int n)
{
	int i;
	int index;

	index = is_nonmirror_subgroup(holo, mero);
	if ( index == 0 ) return 0;

	/* Try to find a group ... */
	for ( i=0; i<n; i++ ) {

		SymOpList *try_mero = list[i];

		/* ... apart from "mero" ... */
		if ( try_mero == mero ) continue;

		/* ... and apart from "holo" ... */
		if ( try_mero == holo ) continue;

		/* ... which is also a subgroup of "holo" ... */
		if ( !is_nonmirror_subgroup(holo, try_mero) ) continue;

		/* ... of which "mero" is also a subgroup. */
		if ( is_nonmirror_subgroup(try_mero, mero) ) return 0;

	}

	return index;
}


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

			int index;

			if ( i == j ) continue;

			index = is_maximal_nonmirror_subgroup(holo, test[j],
			                                      test, n);
			if ( index == 0 ) continue;

			STATUS(" %s(%i)", symmetry_name(test[j]), index);

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
	STATUS("Triclinic to monoclinic:\n");
	find_all_ambiguities("1", "-1",
	                     "2", "m", "2/m",
	                     NULL);

	STATUS("Triclinic to rhombohedral:\n");
	find_all_ambiguities("1", "-1",
	                     "3_R", "32_R", "-3_R", "3m_R", "-3m_R",
	                     NULL);

	STATUS("Triclinic to orthorhombic:\n");
	find_all_ambiguities("1", "-1",
	                     "mm2", "mmm", "222",
	                     NULL);

	STATUS("Orthorhombic to tetragonal:\n");
	find_all_ambiguities("mm2", "mmm", "222",
	                     "4", "-4", "-42m", "-4m2", "4mm", "4/m", "422",
	                          "4/mmm",
	                     NULL);

	STATUS("Monoclinic to tetragonal:\n");

	STATUS("All:\n");
	find_all_ambiguities("1", "-1", "2", "m", "2/m", "mm2", "mmm", "222",
	                     "4", "-4", "-42m", "-4m2", "4mm", "4/m", "422",
	                     "4/mmm", "23", "432", "-43m", "m-3", "m-3m", NULL);

	find_all_ambiguities("1", "-1", "3_R", "32_R", "-3_R", "3m_R", "-3m_R",
	                     NULL);

	find_all_ambiguities("1", "-1", "2", "m", "2/m", "6", "3_H", "312_H",
	                     "321_H", "622", "-3_H", "3m1_H", "-6", "31m_H",
	                     "-3m1_H", "-6m2", "-62m", "-31m_H", "6/mmm", "6/m",
	                     "6mm", NULL);

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
	check_subgroup("622",    "321_H",  1, 1,  2, &fail);

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

	/* Check some pseudo-meroheral subgroups */
	check_subgroup("3_R",  "1", 1, 1,  3, &fail);
	check_subgroup("-3_R",  "-1", 1, 1,  3, &fail);
	check_subgroup("6",  "2", 1, 1,  3, &fail);
	check_subgroup("422",  "222", 1, 1,  2, &fail);

	return fail;
}
