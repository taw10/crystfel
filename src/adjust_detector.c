/*
 * adjust_detector.c
 *
 * Move detector panels
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

#include "version.h"


static void show_syntax(const char *s)
{
	printf("Syntax: %s [options] -g <input.geom> -o <output.geom> [...]\n", s);
}


static void show_help(const char *s)
{
	show_syntax(s);
	printf("\nMove detector panels.\n"
	       "\n"
	       "  -g, --geometry=file        Input geometry file\n"
	       "  -o, --output=file          Output geometry file\n"
	       "  -p, --panel=p              Panel (or group) to move\n"
	       "      --mm                   Interpret shifts as mm, not px\n"
	       "\n"
	       "  -a                         Rotation around x-axis (deg)\n"
	       "  -b                         Rotation around y-axis (deg)\n"
	       "  -c                         Rotation around z-axis (deg)\n"
	       "  -x                         Shift in x direction\n"
	       "  -y                         Shift in y direction\n"
	       "  -z                         Shift in z direction\n"
	       "\n"
	       "  -h, --help                 Display this help message\n"
	       "      --version              Print version number and exit\n");
}


static double parse_double(const char *str, char complain)
{
	double v;
	char *rval;

	errno = 0;
	v = strtod(optarg, &rval);
	if ( *rval != '\0' ) {
		ERROR("Invalid value for -%c.\n", complain);
		exit(1);
	}

	return v;
}


int main(int argc, char *argv[])
{
	int c;
	DataTemplate *dtempl;
	char *in_geom = NULL;
	char *out_geom = NULL;
	double x_shift = 0.0;
	double y_shift = 0.0;
	double z_shift = 0.0;
	double x_rot = 0.0;
	double y_rot = 0.0;
	double z_rot = 0.0;
	int mm = 0;
	char *group = strdup("all");
	int r;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"verbose",            0, NULL,               'v'},

		{"version",            0, NULL,               'V'},
		{"input",              1, NULL,               'g'},
		{"output",             1, NULL,               'o'},
		{"panel",              1, NULL,               'p'},
		{"mm",                 0, NULL,                3},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hVo:g:i:l:p:x:y:z:a:b:c:",
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

			case 'p' :
			free(group);
			group = strdup(optarg);
			break;

			case 3 :
			mm = 1;
			break;

			case 'x' :
			x_shift = parse_double(optarg, 'x');
			break;

			case 'y' :
			x_shift = parse_double(optarg, 'y');
			break;

			case 'z' :
			z_shift = parse_double(optarg, 'z');
			break;

			case 'a' :
			x_rot = parse_double(optarg, 'a');
			break;

			case 'b' :
			y_rot = parse_double(optarg, 'b');
			break;

			case 'c' :
			z_rot = parse_double(optarg, 'c');
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

	if ( (in_geom == NULL) || (out_geom == NULL) ) {
		show_syntax(argv[0]);
		return 1;
	}

	dtempl = data_template_new_from_file(in_geom);
	if ( dtempl == NULL ) return 1;

	data_template_show_hierarchy(dtempl);

	if ( mm ) {
		r = data_template_translate_group_m(dtempl, group,
		                                    x_shift * 1e-3,
		                                    y_shift * 1e-3,
		                                    z_shift * 1e-3);
	} else {
		r = data_template_translate_group_px(dtempl, group,
		                                     x_shift, y_shift, z_shift);
	}

	if ( r ) {
		ERROR("Failed to translate group.\n");
		return 1;
	}

	if ( data_template_rotate_group(dtempl, group, deg2rad(x_rot), 'x') ) {
		ERROR("Failed to rotate group around x.\n");
		return 1;
	}

	if ( data_template_rotate_group(dtempl, group, deg2rad(y_rot), 'y') ) {
		ERROR("Failed to rotate group around y.\n");
		return 1;
	}

	if ( data_template_rotate_group(dtempl, group, deg2rad(z_rot), 'z') ) {
		ERROR("Failed to rotate group around z.\n");
		return 1;
	}

	if ( data_template_write_to_file(dtempl, out_geom) ) {
		ERROR("Failed to save geometry file.\n");
		return 1;
	}

	data_template_free(dtempl);

	return 0;
}
