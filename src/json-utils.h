/*
 * json-utils.h
 *
 * JSON writing utilities
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

static void write_float(FILE *fh, int comma, const char *name, float val)
{
	const char *cs = comma ? "," : "";

	if ( isnan(val) ) {
		fprintf(fh, "    \"%s\": null%s\n", name, cs);
	} else if ( val > 0.0001 ) {
		fprintf(fh, "    \"%s\": %f%s\n", name, val, cs);
	} else {
		fprintf(fh, "    \"%s\": %e%s\n", name, val, cs);
	}
}


static void write_bool(FILE *fh, int comma, const char *name, int val)
{
	fprintf(fh, "    \"%s\": %s%s\n",
	        name, val?"true":"false", comma?",":"");
}


static void write_int(FILE *fh, int comma, const char *name, int val)
{
	fprintf(fh, "    \"%s\": %i%s\n", name, val, comma?",":"");
}


static void write_str(FILE *fh, int comma, const char *name, const char *val)
{
	fprintf(fh, "    \"%s\": \"%s\"%s\n", name, val, comma?",":"");
}

