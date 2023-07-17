/*
 * align_detector.c
 *
 * Align detector using Millepede
 *
 * Copyright © 2023 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Thomas White <taw@physics.org>
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
#include <string.h>
#include <getopt.h>

#include <datatemplate.h>
#include <utils.h>
#include <predict-refine.h>
#include <crystfel-mille.h>

#include "version.h"


static void show_syntax(const char *s)
{
	printf("Syntax: %s [options] -g <input.geom> -o <output.geom> <mille-0.dat> [...]\n", s);
}


static void show_help(const char *s)
{
	show_syntax(s);
	printf("\nRefine detector geometry using Millepede.\n"
	       "\n"
	       "  -g, --geometry=file        Input geometry file\n"
	       "  -o, --output=file          Output geometry file\n"
	       "  -l, --level=n              Alignment hierarchy level\n"
	       "\n"
	       "  -h, --help                 Display this help message\n"
	       "      --version              Print version number and exit\n");
}


int main(int argc, char *argv[])
{
	int c;
	char *in_geom = NULL;
	char *out_geom = NULL;
	int level = 0;
	char *rval;
	int i;
	FILE *fh;
	DataTemplate *dtempl;
	struct dg_group_info *groups;
	int n_groups;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"verbose",            0, NULL,               'v'},

		{"version",            0, NULL,               'V'},
		{"input",              1, NULL,               'g'},
		{"output",             1, NULL,               'o'},
		{"level",              1, NULL,               'l'},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hVo:g:i:l:",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'V' :
			printf("CrystFEL: %s\n", crystfel_version_string());
			printf("%s\n", crystfel_licence_string());
			return 0;

			case 'g' :
			case 'i' :
			in_geom = strdup(optarg);
			break;

			case 'o' :
			out_geom = strdup(optarg);
			break;

			case 'l' :
			errno = 0;
			level = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --level.\n");
				return 1;
			}

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( (in_geom == NULL) || (out_geom == NULL) || (argc == optind) ) {
		show_syntax(argv[0]);
		return 1;
	}

	fh = fopen("millepede.txt", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open Millepede steering file\n");
		return 1;
	}

	for ( i=optind; i<argc; i++ ) {
		fprintf(fh, "%s\n", argv[i]);
	}

	dtempl = data_template_new_from_file(in_geom);
	groups = data_template_group_info(dtempl, &n_groups);

	fprintf(fh, "\nParameter\n");
	for ( i=0; i<n_groups; i++ ) {
		int f = (groups[i].hierarchy_level > level) ? -1 : 0;
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TX), f);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TY), f);
		if ( groups[i].hierarchy_level > 0 ) {
			fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TZ), f);
		} else {
			/* No overall camera length change */
			fprintf(fh, "%i 0 -1\n", mille_label(groups[i].serial, GPARAM_DET_TZ));
		}
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RX), f);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RY), f);
		if ( groups[i].hierarchy_level > 0 ) {
			fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RZ), f);
		} else {
			/* No overall rotation around beam */
			fprintf(fh, "%i 0 -1\n", mille_label(groups[i].serial, GPARAM_DET_RZ));
		}
	}
	fprintf(fh, "\n");

	fprintf(fh, "method inversion 5 0.1\n");
	fprintf(fh, "end\n");
	fclose(fh);

	system("pede millepede.txt");

	return 0;
}
