/*
 * show_residuals.c
 *
 * Extract and display spot position residuals from Millepede files
 *
 * Copyright © 2025 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2025 Thomas White <taw@physics.org>
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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#include <datatemplate.h>
#include <utils.h>
#include <predict-refine.h>
#include <crystfel-mille.h>

#include "version.h"


static void show_syntax(const char *s)
{
	printf("Syntax: %s [options] -g <input.geom> <mille-0.dat> [...]\n", s);
}


static void show_help(const char *s)
{
	show_syntax(s);
	printf("\nDisplay spot position residuals.\n"
	       "\n"
	       "  -g, --geometry=file        Input geometry file\n"
	       "\n"
	       "  -h, --help                 Display this help message\n"
	       "      --version              Print version number and exit\n");
}


static int group_serial_to_number(int serial,
                                  struct dg_group_info *groups,
                                  int n_groups)
{
	int i;

	for ( i=0; i<n_groups; i++ ) {
		if ( groups[i].serial == serial ) return i;
	}

	ERROR("Failed to find group %i\n", serial);
	exit(1);
}


struct stable_running_mean
{
	/* I intended to use a stable average algorithm here,
	 * but there was no round-off error at all with tests up to
	 * 1.5 million records and 60 million measurements. */
	double sum;
	int n_meas;
};


static void init_ave(struct stable_running_mean *ave)
{
	ave->sum = 0.0;
	ave->n_meas = 0;
}


static void add_to_ave(struct stable_running_mean *ave, double v)
{
	ave->sum += v;
	ave->n_meas++;
}


static double calc_mean(struct stable_running_mean v)
{
	return v.sum/v.n_meas;
}


struct gm_ave
{
	struct stable_running_mean x;
	struct stable_running_mean y;
	struct stable_running_mean exerr;
};


static int readint(FILE *fh)
{
	int r;
	if ( fread(&r, 4, 1, fh) != 1 ) {
		ERROR("Failed to read.\n");
		exit(1);
	}
	return r;
}


static void find_markers(int *arri, int start, int len, int *mid, int *next)
{
	int i;
	*mid = 0;
	for ( i=start+1; i<len; i++ ) {
		if ( arri[i] == 0 ) {
			if ( *mid == 0 ) {
				*mid = i;
			} else {
				*next = i;
				return;
			}
		}
	}
	*next = len;
}


static int in_list(int ser, int *serials, int n_serials)
{
	int i;
	for ( i=0; i<n_serials; i++ ) {
		if ( serials[i] == ser ) return 1;
	}
	return 0;
}


static int read_measurement(int *arri, int cur_meas_pos, int arrlen,
                            int *serials, int *n_serials)
{
	int next_meas_pos, mid_sep;
	int lj;

	find_markers(arri, cur_meas_pos, arrlen,
	             &mid_sep, &next_meas_pos);

	for ( lj=mid_sep+1; lj<next_meas_pos; lj++ ) {
		int group_serial = arri[lj] - (arri[lj] % 100);
		if ( !in_list(group_serial, serials, *n_serials) ) {
			serials[*n_serials] = group_serial;
			(*n_serials)++;
		}
	}

	return next_meas_pos;
}


static struct detgeom_panel *which_panel(int *serials,
                                         int n_serials,
                                         struct detgeom *detgeom,
                                         struct dg_group_info *groups,
                                         int n_groups)
{
	int i;
	for ( i=0; i<n_serials; i++ ) {
		int gid = group_serial_to_number(serials[i], groups, n_groups);
		if ( groups[gid].leaf ) {
			return detgeom_find_panel(detgeom, groups[gid].name);
		}
	}
	return NULL;
}


static int read_file(const char *filename,
                     struct gm_ave *aves,
                     struct dg_group_info *groups,
                     int n_groups,
                     struct detgeom *detgeom)
{
	FILE *fh;
	int n_records = 0;

	fh = fopen(filename, "rb");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return 0;
	}

	do {

		int nw, arrlen;
		float *arrf;
		int *arri;

		if ( fread(&nw, sizeof(int), 1, fh) != 1 ) {
			if ( feof(fh) ) {
				break;
			} else {
				ERROR("Failed to read record header\n");
				return 0;
			}
		}

		arrlen = (nw-2)/2;

		arrf = malloc(arrlen*sizeof(float));
		arri = malloc(arrlen*sizeof(int));
		if ( (arri == NULL) || (arrf == NULL) ) {
			ERROR("Failed to allocate memory\n");
			return 0;
		}

		readint(fh);  /* And ignore it */
		if ( fread(arrf, 4, arrlen, fh) != arrlen ) {
			ERROR("Failed to read float array\n");
			return 0;
		}

		readint(fh);  /* And ignore it */
		if ( fread(arri, 4, arrlen, fh) != arrlen ) {
			ERROR("Failed to read integer array\n");
			return 0;
		}

		n_records++;

		int cur_meas_pos = 0;
		while ( cur_meas_pos < arrlen ) {

			double fs, ss, ex;
			int j;
			int n_serials = 0;
			int serials[32];
			struct detgeom_panel *p;
			double dx, dy;

			fs = arrf[cur_meas_pos];
			cur_meas_pos = read_measurement(arri, cur_meas_pos, arrlen,
			                                serials, &n_serials);
			ss = arrf[cur_meas_pos];
			cur_meas_pos = read_measurement(arri, cur_meas_pos, arrlen,
			                                serials, &n_serials);
			ex = arrf[cur_meas_pos];
			cur_meas_pos = read_measurement(arri, cur_meas_pos, arrlen,
			                                serials, &n_serials);

			/* Mille residuals are in pixels, we want metres */
			p = which_panel(serials, n_serials, detgeom, groups, n_groups);
			dx = (fs*p->fsx + ss*p->ssx)*p->pixel_pitch;
			dy = (fs*p->fsy + ss*p->ssy)*p->pixel_pitch;

			for ( j=0; j<n_serials; j++ ) {
				int gid;
				gid = group_serial_to_number(serials[j], groups, n_groups);
				add_to_ave(&aves[gid].x, dx);
				add_to_ave(&aves[gid].y, dy);
				add_to_ave(&aves[gid].exerr, fabs(ex/EXC_WEIGHT));
			}

		}

		free(arrf);
		free(arri);

	} while ( !feof(fh) );

	STATUS("Read %i records from %s\n", n_records, filename);
	return n_records;
}


static int read_children(const char *dir,
                         struct gm_ave *aves,
                         struct dg_group_info *groups,
                         int n_groups,
                         struct detgeom *dg)
{
	DIR *d = opendir(dir);
	struct dirent *dent = readdir(d);
	int n_records = 0;

	while ( dent != NULL ) {
		if ( dent->d_name[0] == 'm' ) {
			char *tmp = malloc(strlen(dir)+strlen(dent->d_name)+2);
			sprintf(tmp, "%s/%s\n", dir, dent->d_name);
			n_records += read_file(tmp, aves, groups, n_groups, dg);
			free(tmp);
		}
		dent = readdir(d);
	}

	closedir(d);
	return n_records;
}


int main(int argc, char *argv[])
{
	int c;
	char *in_geom = NULL;
	int i;
	DataTemplate *dtempl;
	struct dg_group_info *groups;
	int n_groups;
	int r;
	struct gm_ave *aves;
	int lvl;
	int total_records = 0;
	struct detgeom *dg;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"verbose",            0, NULL,               'v'},
		{"version",            0, NULL,               'V'},
		{"input",              1, NULL,               'g'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hVg:i:",
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

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( (in_geom == NULL) || (argc == optind) ) {
		show_syntax(argv[0]);
		return 1;
	}

	dtempl = data_template_new_from_file(in_geom);
	groups = data_template_group_info(dtempl, &n_groups);
	dg = data_template_get_2d_detgeom_if_possible(dtempl);
	if ( dg == NULL ) {
		ERROR("Can't transform residuals because geometry is not static\n");
		return 1;
	}

	aves = malloc(n_groups*sizeof(struct gm_ave));
	if ( aves == NULL ) {
		ERROR("Failed to allocate memory for residual calculation\n");
		return 1;
	}

	for ( i=0; i<n_groups; i++ ) {
		init_ave(&aves[i].x);
		init_ave(&aves[i].y);
		init_ave(&aves[i].exerr);
	}

	for ( i=optind; i<argc; i++ ) {

		struct stat statbuf;

		r = stat(argv[i], &statbuf);
		if ( r != 0 ) {
			ERROR("File/directory '%s' not found\n", argv[i]);
			return 1;
		}

		if ( is_dir(argv[i]) ) {
			total_records += read_children(argv[i], aves, groups, n_groups, dg);
		} else {
			total_records += read_file(argv[i], aves, groups, n_groups, dg);
		}
	}

	STATUS("Read %i records in total\n", total_records);

	for ( lvl=10; lvl>=0; lvl-- ) {
		int done_header = 0;
		for ( i=0; i<n_groups; i++ ) {
			if ( groups[i].hierarchy_level == lvl ) {
				if ( !done_header ) {
					STATUS("\n\n");
					STATUS("Hierarchy level %i:\n\n", lvl);
					done_header = 1;
				}
				STATUS("Group %s:\n", groups[i].name);
				STATUS("   Mean spot deviation in x-direction "
				       "(%i measurements) = %+f µm\n",
				       aves[i].x.n_meas, calc_mean(aves[i].x)*1e6);
				STATUS("   Mean spot deviation in y-direction "
				       "(%i measurements) = %+f µm\n",
				       aves[i].x.n_meas, calc_mean(aves[i].y)*1e6);
				STATUS("   Mean absolute reflection excitation error "
				       "(%i measurements) = %+f nm^-1\n",
				       aves[i].exerr.n_meas, 1e-9*calc_mean(aves[i].exerr));
			}
		}
	}

	STATUS("\n\nNote that the position residuals are NOT to be taken as the "
	       "shifts needed to refine the detector geometry.\n");
	STATUS("Run align_detector to calculate the correct shifts and update "
	       "the geometry file.\n\n");

	free(aves);

	return 0;
}
