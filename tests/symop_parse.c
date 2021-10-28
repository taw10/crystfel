/*
 * symop_parse.c
 *
 * Check that symmetry operation parsing works
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <symmetry.h>
#include <rational.h>

int main(int argc, char *argv[])
{
	int r = 0;
	RationalMatrix *mtx;
	SymOpList *sym;

	mtx = parse_symmetry_operation("h,k,l");
	if ( !rtnl_mtx_is_identity(mtx) ) {
		printf("h,k,l not an identity:\n");
		rtnl_mtx_print(mtx);
		r = 1;
	}

	mtx = parse_symmetry_operation("k,h,-l");

	mtx = parse_symmetry_operation("h,k,l");
	if ( !rtnl_mtx_is_identity(mtx) ) {
		printf("h,k,l not an identity on second attempt:\n");
		rtnl_mtx_print(mtx);
		r = 1;
	}

	sym = parse_symmetry_operations("h,k,l;k,h,-l;-h,-k,l");
	if ( sym == NULL ) r = 1;

	mtx = parse_symmetry_operation("h,k,fail");
	if ( mtx != NULL ) r = 1;

	sym = parse_symmetry_operations("k,h,-l;h,k,fail;h,k,l");
	if ( sym != NULL ) r = 1;

	return r;
}
