/*
 * list_events.c
 *
 * Generate event lists
 *
 * Copyright © 2015 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2015 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "version.h"
#include "utils.h"
#include "detector.h"
#include "hdf5-file.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] -i files.lst -o events.lst "
	       "-g geometry.geom\n\n", s);
	printf(
"Generate event lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"\n"
"  -i, --input=<file>         Input filename (list of multi-event filenames).\n"
"  -g, --geometry=<file>      Get data layout from geometry file.\n"
"  -o, --output=<file>        Output filename (list of events).\n"
);
}


int main(int argc, char *argv[])
{
	int c;
	char *input = NULL;
	char *output = NULL;
	char *geom = NULL;
	char *rval;
	FILE *ifh;
	FILE *ofh;
	struct detector *det;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                2 },
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{"output",             1, NULL,               'o'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:o:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 2 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'o' :
			output = strdup(optarg);
			break;

			case 'i' :
			input = strdup(optarg);
			break;

			case 'g' :
			geom = strdup(optarg);
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( (input == NULL) || (output == NULL) || (geom == NULL) ) {
		ERROR("You must specify at least the input, output and geometry"
		      " filenames.\n");
		return 1;
	}

	ifh = fopen(input, "r");
	if ( ifh == NULL ) {
		ERROR("Couldn't open '%s'\n", input);
		return 1;
	}

	ofh = fopen(output, "w");
	if ( ofh == NULL ) {
		ERROR("Couldn't open '%s'\n", output);
		return 1;
	}

	det = get_detector_geometry(geom, NULL);
	if ( det == NULL ) {
		ERROR("Failed to read '%s'\n", geom);
		return 1;
	}

	if ( (det->path_dim == 0) && (det->dim_dim == 0) ) {
		ERROR("This does not look like a multi-event geometry file.\n");
		ERROR("Are you sure you need to use list_events instead of "
		      "just 'find' or 'ls'?\n");
		return 1;
	}

	do {

		char filename[1024];
		int i;

		rval = fgets(filename, 1024, ifh);
		if ( rval != NULL ) {

			struct event_list *evlist;
			struct hdfile *hdfile;

			chomp(filename);

			hdfile = hdfile_open(filename);
			if ( hdfile == NULL ) {
				ERROR("Failed to open '%s'\n", filename);
				ERROR("Aborting creation of event list.\n");
				return 1;
			}

			evlist = fill_event_list(hdfile, det);

			for ( i=0; i<evlist->num_events; i++ ) {
				char *str = get_event_string(evlist->events[i]);
				fprintf(ofh, "%s %s\n", filename, str);
				free(str);
			}

			STATUS("%i events found in %s\n", evlist->num_events,
			       filename);

			free_event_list(evlist);
			hdfile_close(hdfile);
		}

	} while ( rval != NULL );

	return 0;
}
