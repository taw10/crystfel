/*
 * make_pixelmap.c
 *
 * Create a pixel map file from a CrystFEL geometry description
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2018 Thomas White <taw@physics.org>
 *   2016-2016 Omri Mor <omor1@asu.edu>
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

#define _ISOC99_SOURCE
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>

#include "detector.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <input.geom>\n\n", s);
	printf(
"Create a pixel map file.\n"
"\n"
" -h, --help                Display this help message.\n"
"\n"
" -o, --output=<filename>   Output filename.  Default: <input>.h5.\n"
"     --badmap              Generate bad pixel map instead of geometry\n"
"     --good-pixel=<n>      Value for good pixels in bad map.  Default 1.\n"
"     --bad-pixel=<n>       Value for bad pixels in bad map.  Default 0.\n"
);
}


static void create_array(hid_t gh, const char *name, void *vals,
                         hid_t type, int width, int height)
{
	hid_t dh, sh;
	herr_t r;
	hsize_t size[2];
	hsize_t max_size[2];

	/* Note the "swap" here, according to section 3.2.5,
	 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
	size[0] = height;
	size[1] = width;
	max_size[0] = height;
	max_size[1] = width;
	sh = H5Screate_simple(2, size, max_size);

	dh = H5Dcreate2(gh, name, type, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		return;
	}

	r = H5Dwrite(dh, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		return;
	}

	H5Sclose(sh);
	H5Dclose(dh);
}


static void create_scalar(hid_t gh, const char *name, float val)
{
	hid_t dh, sh;
	herr_t r;

	sh = H5Screate(H5S_SCALAR);

	STATUS("%s: %e\n", name, val);

	dh = H5Dcreate2(gh, name, H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		return;
	}

	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		return;
	}

	H5Sclose(sh);
	H5Dclose(dh);
}


static void write_pixelmap_hdf5(const char *filename,
                                float *x, float *y, float *z,
                                int width, int height, float res, float coffset)
{
	hid_t fh;

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return;
	}

	create_array(fh, "x", x, H5T_NATIVE_FLOAT, width, height);
	create_array(fh, "y", y, H5T_NATIVE_FLOAT, width, height);
	create_array(fh, "z", z, H5T_NATIVE_FLOAT, width, height);

	create_scalar(fh, "res", res);
	create_scalar(fh, "coffset", coffset);

	H5Fclose(fh);
}


static void write_badmap_hdf5(const char *filename, uint16_t *b,
                              int width, int height)
{
	hid_t fh, gh;

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return;
	}

	gh = H5Gcreate(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	create_array(gh, "/data/data", b, H5T_STD_U16LE, width, height);

	H5Fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	char *input_file = NULL;
	char *output_file = NULL;
	struct detector *det = NULL;
	int fs, ss, w, h;
	float *x, *y, *z;
	uint16_t *b;
	int i;
	float res, coffset;
	int badmap = 0;
	int good_pixel_val = 1;
	int bad_pixel_val = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"badmap",             0, &badmap,             1},
		{"good-pixel",         1, NULL,              301},
		{"bad-pixel",          1, NULL,              302},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'o' :
			output_file = strdup(optarg);
			break;

			case 301:
			if (sscanf(optarg, "%d", &good_pixel_val) != 1)
			{
				ERROR("Invalid value for --good-pixel\n");
				return 1;
			}
			break;

			case 302:
			if (sscanf(optarg, "%d", &bad_pixel_val) != 1)
			{
				ERROR("Invalid value for --bad-pixel\n");
				return 1;
			}
			break;

			case 0 :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( argc != (optind+1) ) {
		ERROR("Please provide exactly one geometry file name.\n");
		return 1;
	}

	input_file = strdup(argv[optind++]);

	if ( output_file == NULL ) {
		size_t len = strlen(input_file);
		output_file = malloc(len+4);
		strcpy(output_file, input_file);
		strip_extension(output_file);
		strcat(output_file, ".h5");
	}

	/* Load geometry */
	det = get_detector_geometry(input_file, NULL);
	if ( det == NULL ) {
		ERROR("Failed to read geometry from '%s'\n", input_file);
		return 1;
	}
	free(input_file);

	/* Determine max orig fs and ss */
	w = 0;  h = 0;
	for ( i=0; i<det->n_panels; i++ ) {
		if ( det->panels[i].orig_max_fs > w ) {
			w = det->panels[i].orig_max_fs;
		}
		if ( det->panels[i].orig_max_ss > h ) {
			h = det->panels[i].orig_max_ss;
		}
	}
	w += 1;  h += 1;  /* Inclusive -> Exclusive */
	STATUS("Data slab size: %i x %i\n", w, h);

	x = malloc(w*h*sizeof(float));
	y = malloc(w*h*sizeof(float));
	z = malloc(w*h*sizeof(float));
	b = malloc(w*h*sizeof(uint16_t));
	if ( (x==NULL) || (y==NULL) || (z==NULL) || (b==NULL) ) {
		ERROR("Failed to allocate memory.\n");
		return 1;
	}

	for ( ss=0; ss<h; ss++ ) {
		for ( fs=0; fs<w; fs++ ) {

			double rx, ry;
			struct panel *p;
			double xs, ys;
			double cfs, css;
			int nfs, nss;

			p = find_orig_panel(det, fs, ss);

			/* Add half a pixel to fs and ss to get the fs,ss
			 * coordinates of the CENTRE of the pixel */
			cfs = fs + 0.5;
			css = ss + 0.5;
			xs = (cfs - p->orig_min_fs)*p->fsx
			      + (css - p->orig_min_ss)*p->ssx;
			ys = (cfs - p->orig_min_fs)*p->fsy
			     + (css - p->orig_min_ss)*p->ssy;

			rx = (xs + p->cnx) / p->res;
			ry = (ys + p->cny) / p->res;

			x[fs + w*ss] = rx;
			y[fs + w*ss] = ry;
			z[fs + w*ss] = 0.0;  /* FIXME */

			nfs = fs - p->orig_min_fs;
			nss = ss - p->orig_min_ss;
			if ( in_bad_region(det, p, nfs, nss) ) {
				b[fs + w*ss] = bad_pixel_val;
			} else {
				b[fs + w*ss] = good_pixel_val;
			}

		}

		progress_bar(ss, h, "Converting");
	}
	STATUS("\n");

	res = det->defaults.res;
	/* use the res of the first panel if not in global definitions
	 * assumes one of these exist */
	if ( res == -1.0 ) {
		res = det->panels[0].res;
	}

	coffset = det->defaults.coffset;
	/* use the coffset of the first panel if not in global definitions
	 * assumes one of these exist*/
	if ( coffset == 0.0 ) {
		coffset = det->panels[0].coffset;
	}

	if ( badmap ) {
		write_badmap_hdf5(output_file, b, w, h);
	} else {
		write_pixelmap_hdf5(output_file, x, y, z, w, h, res, coffset);
	}

	return 0;
}
