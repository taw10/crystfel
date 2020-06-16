/*
 * evparse7.c
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdarg.h>

#include "../libcrystfel/src/image-hdf5.c"

int main(int argc, char *argv[])
{
	char *subs;

	subs = substitute_path("cc/data123/bb//234/59",
	                       "/data/%/test/%/%");

	if ( subs == NULL ) {
		printf("substitute_path() failed\n");
		return 1;
	}

	if ( strcmp(subs, "/data/cc/test/data123/bb") != 0 ) {
		printf("Wrong substitution: got '%s'\n", subs);
		return 1;
	}

	return 0;
}
