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

#include <unistd.h>
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


static const char *str_param(enum gparam param)
{
	switch ( param ) {
		case GPARAM_DET_TX : return "x-translation";
		case GPARAM_DET_TY : return "y-translation";
		case GPARAM_DET_TZ : return "z-translation";
		case GPARAM_DET_RX : return "x-rotation";
		case GPARAM_DET_RY : return "y-rotation";
		case GPARAM_DET_RZ : return "z-rotation";
		default : return "(unknown)";
	}
}


static const char *group_serial_to_name(int serial,
                                        struct dg_group_info *groups,
                                        int n_groups)
{
	int i;

	for ( i=0; i<n_groups; i++ ) {
		if ( groups[i].serial == serial ) return groups[i].name;
	}

	return NULL;
}


static struct dg_group_info *find_group(struct dg_group_info *groups,
                                        int n_groups, const char *name)
{
	int i;

	for ( i=0; i<n_groups; i++ ) {
		if ( strcmp(name, groups[i].name) == 0 ) return &groups[i];
	}

	return NULL;
}

static int ipow(int base, int ex)
{
	int i;
	int v = 1;
	for ( i=0; i<ex; i++ ) {
		v *= base;
	}
	return v;
}


static int is_child(struct dg_group_info *parent, struct dg_group_info *child)
{
	int parent_serial;

	if ( 1+parent->hierarchy_level != child->hierarchy_level ) return 0;

	parent_serial = child->serial % ipow(100, child->hierarchy_level);
	if ( parent->serial != parent_serial ) return 0;

	return 1;
}


static void write_zero_sum(FILE *fh, struct dg_group_info *g,
                           struct dg_group_info *groups, int n_groups,
                           enum gparam p)
{
	int i;
	int found = 0;

	for ( i=0; i<n_groups; i++ ) {
		if ( is_child(g, &groups[i]) ) {
			if ( !found ) {
				fprintf(fh, "Constraint 0\n");
				found = 1;
			}
			fprintf(fh, "%i 1\n", mille_label(groups[i].serial, p));
		}
	}

	if ( found ) {
		fprintf(fh, "\n");
	}
}


static int make_zero_sum(FILE *fh, struct dg_group_info *groups, int n_groups,
                         const char *group_name, int level)
{
	int i;
	struct dg_group_info *g = find_group(groups, n_groups, group_name);

	if ( g == NULL ) {
		ERROR("Couldn't find group '%s'\n", group_name);
		return 1;
	}

	/* Millepede doesn't like excessive constraints */
	if ( g->hierarchy_level >= level ) return 0;

	fprintf(fh, "! Hierarchy constraints for group %s\n", group_name);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_TX);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_TY);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_TZ);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_RX);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_RY);
	write_zero_sum(fh, g, groups, n_groups, GPARAM_DET_RZ);
	fprintf(fh, "\n");

	for ( i=0; i<n_groups; i++ ) {
		if ( is_child(g, &groups[i]) ) {
			if ( make_zero_sum(fh, groups, n_groups, groups[i].name, level) ) return 1;
		}
	}

	return 0;
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
	int r;
	char line[256];

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

	/* Top level */
	fprintf(fh, "%i 0 0\n", mille_label(0, GPARAM_DET_TX));
	fprintf(fh, "%i 0 0\n", mille_label(0, GPARAM_DET_TY));
	fprintf(fh, "%i 0 -1\n", mille_label(0, GPARAM_DET_TZ));
	fprintf(fh, "%i 0 -1\n", mille_label(0, GPARAM_DET_RX));
	fprintf(fh, "%i 0 -1\n", mille_label(0, GPARAM_DET_RY));
	fprintf(fh, "%i 0 -1\n", mille_label(0, GPARAM_DET_RZ));

	for ( i=0; i<n_groups; i++ ) {
		int f = (groups[i].hierarchy_level > level) ? -1 : 0;
		if ( groups[i].hierarchy_level == 0 ) continue;
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TX), f);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TY), f);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_TZ), -1);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RX), -1);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RY), -1);
		fprintf(fh, "%i 0 %i\n", mille_label(groups[i].serial, GPARAM_DET_RZ), f);
	}
	fprintf(fh, "\n");

	/* All corrections must sum to zero at each level of hierarchy */
	if ( make_zero_sum(fh, groups, n_groups, "all", level) ) return 1;

	fprintf(fh, "method inversion 5 0.1\n");
	fprintf(fh, "skipemptycons\n");
	fprintf(fh, "end\n");
	fclose(fh);

	unlink("millepede.res");

	r = system("pede millepede.txt");
	if ( r == -1 ) {
		ERROR("Failed to run Millepde: %s\n", strerror(errno));
		return 1;
	}
	if ( !WIFEXITED(r) ) {
		ERROR("Millepede exited abnormally.\n");
		return 1;
	}
	if ( WEXITSTATUS(r) != 0 ) {
		ERROR("Millepede returned an error status (%i)\n", WEXITSTATUS(r));
		return 1;
	}

	STATUS("Millepede succeeded.\n");

	fh = fopen("millepede.res", "r");
	if ( fh == NULL ) {
		ERROR("Failed to open millepede.res\n");
		return 1;
	}

	if ( fgets(line, 256, fh) != line ) {
		ERROR("Failed to read first line of millepede.res\n");
		return 1;
	}
	if ( strncmp(line, " Parameter ", 11) != 0 ) {
		ERROR("First line of millepede.res is not as expected.\n");
		return 1;
	}

	int last_group_serial = -1;
	do {

		char **bits;
		int i, n;

		rval = fgets(line, 256, fh);
		if ( rval != line ) continue;

		chomp(line);
		notrail(line);
		n = assplode(line, " ", &bits, ASSPLODE_NONE);
		if ( (n != 3) && (n != 5) ) {
			ERROR("Didn't understand this line from Millepede: (%i) %s", n, line);
			return 1;
		}

		if ( n == 5 ) {

			int code;
			double shift;
			int p;
			const char *group_name;
			int group_serial;

			if ( convert_int(bits[0], &code) ) {
				ERROR("Didn't understand '%s'\n", bits[0]);
				return 1;
			}
			if ( convert_float(bits[1], &shift) ) {
				ERROR("Didn't understand '%s'\n", bits[1]);
				return 1;
			}

			p = mille_unlabel(code % 100);
			group_serial = code - (code % 100);
			group_name = group_serial_to_name(group_serial,
			                                  groups,
			                                  n_groups);

			if ( last_group_serial != group_serial ) {
				STATUS("Group %s:\n", group_name);
				last_group_serial = group_serial;
			}

			switch ( p ) {
				case GPARAM_DET_TX:
				case GPARAM_DET_TY:
				case GPARAM_DET_TZ:
				STATUS("   %14s %+f mm\n", str_param(p), 1e3*shift);
				break;

				case GPARAM_DET_RX:
				case GPARAM_DET_RY:
				case GPARAM_DET_RZ:
				STATUS("   %14s %+f deg\n", str_param(p), rad2deg(shift));
				break;
			}

			if ( group_name == NULL ) {
				ERROR("Invalid group serial number %i\n", code);
				return 1;
			}

			switch ( p ) {

				case GPARAM_DET_TX:
				data_template_translate_group_m(dtempl, group_name,
				                                -shift, 0, 0);
				break;

				case GPARAM_DET_TY:
				data_template_translate_group_m(dtempl, group_name,
				                                0, -shift, 0);
				break;

				case GPARAM_DET_TZ:
				data_template_translate_group_m(dtempl, group_name,
				                                0, 0, -shift);
				break;

				case GPARAM_DET_RX:
				data_template_rotate_group(dtempl, group_name,
				                           -shift, 'x');
				break;

				case GPARAM_DET_RY:
				data_template_rotate_group(dtempl, group_name,
				                           -shift, 'y');
				break;

				case GPARAM_DET_RZ:
				data_template_rotate_group(dtempl, group_name,
				                           -shift, 'z');
				break;

				default:
				ERROR("Invalid parameter %i\n", p);
				return 1;
			}

		}

		for ( i=0; i<n; i++ ) free(bits[i]);
		free(bits);

	} while ( rval == line );

	fclose(fh);

	data_template_write_to_file(dtempl, out_geom);

	return 0;
}
