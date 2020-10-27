/*
 * evparse4.c
 *
 * Check that event string parsing works
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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
#include <stdarg.h>

extern int *read_dim_parts(const char *ev_orig, int *pn_dvals);

int main(int argc, char *argv[])
{
	int *dvals;
	int n_dvals = 99;
	int r = 0;

	dvals = read_dim_parts("cc/data123/bb//234/59", &n_dvals);

	if ( n_dvals != 2 ) {
		printf("Wrong number of dimension parts\n");
		r++;
	}

	if ( dvals == NULL ) {
		printf("read_dim_parts failed\n");
		return 1;
	}

	if ( dvals[0] != 234 ) {
		printf("First dimension part is wrong\n");
		r++;
	}
	if ( dvals[1] != 59 ) {
		printf("Second dimension part is wrong "
		       "(%i, should be %i)\n", dvals[1], 59);
		r++;
	}

	return r;
}
