/*
 * symmetry_check.c
 *
 * Check symmetry
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
		STATUS("%15s : NULL!\n", pg);
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
	check_pg_props( "-43m",  24,  0, &fail);
	check_pg_props( "m-3m",  48,  1, &fail);
	STATUS("\n");

	/* Check some weird settings */
	STATUS("\nWeird settings:\n");
	check_pg_props( "2_uaa", 2, 0, &fail);
	check_pg_props( "2_uab", 2, 0, &fail);
	check_pg_props( "2_uac", 2, 0, &fail);
	check_pg_props( "4_uaa", 4, 0, &fail);
	check_pg_props( "4_uab", 4, 0, &fail);
	check_pg_props( "4_uac", 4, 0, &fail);
	check_pg_props( "4/m_uaa", 8, 1, &fail);
	check_pg_props( "4/m_uab", 8, 1, &fail);
	check_pg_props( "4/m_uac", 8, 1, &fail);
	check_pg_props( "23_uaa", 12, 0, &fail);
	check_pg_props( "23_uab", 12, 0, &fail);
	check_pg_props( "23_uac", 12, 0, &fail);
	check_pg_props( "6_uaa", 6, 0, &fail);
	check_pg_props( "6_uab", 6, 0, &fail);
	check_pg_props( "6_uac", 6, 0, &fail);

	return fail;
}
