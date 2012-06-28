/*
 * pattern_sim.c
 *
 * Simulate diffraction patterns from small crystals
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
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

#include "image.h"
#include "diffraction.h"
#include "diffraction-gpu.h"
#include "cell.h"
#include "utils.h"
#include "hdf5-file.h"
#include "detector.h"
#include "peaks.h"
#include "beam-parameters.h"
#include "symmetry.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "pattern_sim.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Simulate diffraction patterns from small crystals probed with femtosecond\n"
"pulses of X-rays from a free electron laser.\n"
"\n"
" -h, --help                Display this help message.\n"
"\n"
" -p, --pdb=<file>          PDB file from which to get the unit cell.\n"
"                            (The actual Bragg intensities come from the\n"
"                            intensities file)\n"
"     --gpu                 Use the GPU to speed up the calculation.\n"
"     --gpu-dev=<n>         Use GPU device <n>.  Omit this option to see the\n"
"                            available devices.\n"
" -g, --geometry=<file>     Get detector geometry from file.\n"
" -b, --beam=<file>         Get beam parameters from file.\n"
"\n"
" -n, --number=<N>          Generate N images.  Default 1.\n"
"     --no-images           Do not output any HDF5 files.\n"
" -o, --output=<filename>   Output HDF5 filename.  Default: sim.h5.\n"
"                            If you choose to simulate more than one pattern,\n"
"                            the filename given will be postfixed with a\n"
"                            hyphen, the image number and '.h5'.  In this\n"
"                            case, the default value is 'sim', such that the\n"
"                            files produced are sim-1.h5, sim-2.h5 and so on.\n"
" -r, --random-orientation  Use a randomly generated orientation\n"
"                            (a new orientation will be used for each image).\n"
"     --powder=<file>       Write a summed pattern of all images simulated by\n"
"                            this invocation as the given filename.\n"
" -i, --intensities=<file>  Specify file containing reflection intensities\n"
"                            (and phases) to use.\n"
" -y, --symmetry=<sym>      The symmetry of the intensities file.\n"
" -t, --gradients=<method>  Use <method> for the calculation of shape\n"
"                            transform intensities.  Choose from:\n"
"                             mosaic      : Take the intensity of the nearest\n"
"                                           Bragg position.  This is the\n"
"                                           fastest method and the only one\n"
"                                           supported on the GPU.\n"
"                             interpolate : Interpolate trilinearly between\n"
"                                           six adjacent Bragg intensities.\n"
"                                           This method has intermediate\n"
"                                           accuracy.\n"
"                             phased      : As 'interpolate', but take phase\n"
"                                           values into account.  This is the\n"
"                                           most accurate method, but the\n"
"                                           slowest.\n"
"     --really-random       Use really random numbers for the orientation and\n"
"                            crystal size.  By default, the same sequence\n"
"                            will be used for each new run.\n"
"     --min-size=<s>        Generate random crystal sizes using <s> as the\n"
"                            minimum crystal size in nm. --max-size is also\n"
"                            required.\n"
"     --max-size=<s>        Use <s> as the maximum crystal size in nm.\n"
"                            --min-size is also required.\n"
"     --no-noise            Do not calculate Poisson noise.\n"
);
}


static double *intensities_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	double *out = new_arr_intensity();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double intensity = get_intensity(refl);

		get_indices(refl, &h, &k, &l);

		set_arr_intensity(out, h, k, l, intensity);

	}

	return out;
}


static double *phases_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	double *out = new_arr_phase();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double phase = get_phase(refl, NULL);

		get_indices(refl, &h, &k, &l);

		set_arr_phase(out, h, k, l, phase);

	}

	return out;

}


static unsigned char *flags_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	unsigned char *out = new_arr_flag();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		set_arr_flag(out, h, k, l, 1);

	}

	return out;

}


static struct quaternion read_quaternion()
{
	do {

		int r;
		float w, x, y, z;
		char line[1024];
		char *rval;

		printf("Please input quaternion: w x y z\n");
		rval = fgets(line, 1023, stdin);
		if ( rval == NULL ) return invalid_quaternion();
		chomp(line);

		r = sscanf(line, "%f %f %f %f", &w, &x, &y, &z);
		if ( r == 4 ) {

			struct quaternion quat;

			quat.w = w;
			quat.x = x;
			quat.y = y;
			quat.z = z;

			return quat;

		} else {
			ERROR("Invalid rotation '%s'\n", line);
		}

	} while ( 1 );
}


static int random_ncells(double len, double min_nm, double max_nm)
{
	int min_cells, max_cells, wid;

	min_cells = min_nm*1e-9 / len;
	max_cells = max_nm*1e-9 / len;
	wid = max_cells - min_cells;

	return wid*random()/RAND_MAX + min_cells;
}


int main(int argc, char *argv[])
{
	int c;
	struct image image;
	struct gpu_context *gctx = NULL;
	double *powder;
	char *intfile = NULL;
	double *intensities;
	char *rval;
	double *phases;
	unsigned char *flags;
	int config_randomquat = 0;
	int config_noimages = 0;
	int config_nonoise = 0;
	int config_nosfac = 0;
	int config_gpu = 0;
	int config_random = 0;
	char *powder_fn = NULL;
	char *filename = NULL;
	char *grad_str = NULL;
	char *outfile = NULL;
	char *geometry = NULL;
	char *beamfile = NULL;
	GradientMethod grad;
	int ndone = 0;    /* Number of simulations done (images or not) */
	int number = 1;   /* Number used for filename of image */
	int n_images = 1; /* Generate one image by default */
	int done = 0;
	UnitCell *input_cell;
	struct quaternion orientation;
	int gpu_dev = -1;
	int random_size = 0;
	double min_size = 0.0;
	double max_size = 0.0;
	char *sym_str = NULL;
	SymOpList *sym;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"gpu",                0, &config_gpu,         1},
		{"random-orientation", 0, NULL,               'r'},
		{"number",             1, NULL,               'n'},
		{"no-images",          0, &config_noimages,    1},
		{"no-noise",           0, &config_nonoise,     1},
		{"intensities",        1, NULL,               'i'},
		{"symmetry",           1, NULL,               'y'},
		{"powder",             1, NULL,               'w'},
		{"gradients",          1, NULL,               't'},
		{"pdb",                1, NULL,               'p'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"really-random",      0, &config_random,      1},
		{"gpu-dev",            1, NULL,                2},
		{"min-size",           1, NULL,                3},
		{"max-size",           1, NULL,                4},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hrn:i:t:p:o:g:b:y:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'r' :
			config_randomquat = 1;
			break;

			case 'n' :
			n_images = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid number of images.\n");
				return 1;
			}
			break;

			case 'i' :
			intfile = strdup(optarg);
			break;

			case 't' :
			grad_str = strdup(optarg);
			break;

			case 'p' :
			filename = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'w' :
			powder_fn = strdup(optarg);
			break;

			case 'g' :
			geometry = strdup(optarg);
			break;

			case 'b' :
			beamfile = strdup(optarg);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 2 :
			gpu_dev = atoi(optarg);
			break;

			case 3 :
			min_size = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid minimum size.\n");
				return 1;
			}
			random_size++;
			break;

			case 4 :
			max_size = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid maximum size.\n");
				return 1;
			}
			random_size++;
			break;

			case 0 :
			break;

			default :
			return 1;

		}

	}

	if ( config_random ) {
		FILE *fh;
		unsigned int seed;
		fh = fopen("/dev/urandom", "r");
		fread(&seed, sizeof(seed), 1, fh);
		fclose(fh);
		srand(seed);
	}

	if ( random_size == 1 ) {
		ERROR("You must specify both --min-size and --max-size.\n");
		return 1;
	}

	if ( filename == NULL ) {
		filename = strdup("molecule.pdb");
	}

	if ( outfile == NULL ) {
		if ( n_images == 1 ) {
			outfile = strdup("sim.h5");
		} else {
			outfile = strdup("sim");
		}
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	/* sym_str is used below */

	if ( grad_str == NULL ) {
		STATUS("You didn't specify a gradient calculation method, so"
		       " I'm using the 'mosaic' method, which is fastest.\n");
		grad = GRADIENT_MOSAIC;
	} else if ( strcmp(grad_str, "mosaic") == 0 ) {
		grad = GRADIENT_MOSAIC;
	} else if ( strcmp(grad_str, "interpolate") == 0) {
		grad = GRADIENT_INTERPOLATE;
	} else if ( strcmp(grad_str, "phased") == 0) {
		grad = GRADIENT_PHASED;
	} else {
		ERROR("Unrecognised gradient method '%s'\n", grad_str);
		return 1;
	}
	free(grad_str);

	if ( config_gpu && (grad != GRADIENT_MOSAIC) ) {
		ERROR("Only the mosaic method can be used for gradients when"
		      "calculating on the GPU.\n");
		return 1;
	}

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	if ( beamfile == NULL ) {
		ERROR("You need to specify a beam parameter file"
		      " with --beam\n");
		return 1;
	}

	if ( intfile == NULL ) {

		/* Gentle reminder */
		STATUS("You didn't specify the file containing the ");
		STATUS("reflection intensities (with --intensities).\n");
		STATUS("I'll simulate a flat intensity distribution.\n");
		intensities = NULL;
		phases = NULL;
		flags = NULL;

	} else {

		RefList *reflections;

		reflections = read_reflections(intfile);
		if ( reflections == NULL ) {
			ERROR("Problem reading input file %s\n", intfile);
			return 1;
		}
		free(intfile);

		if ( grad == GRADIENT_PHASED ) {
			phases = phases_from_list(reflections);
		} else {
			phases = NULL;
		}
		intensities = intensities_from_list(reflections);
		phases = phases_from_list(reflections);
		flags = flags_from_list(reflections);

		/* Check that the intensities have the correct symmetry */
		if ( check_list_symmetry(reflections, sym) ) {
			ERROR("The input reflection list does not appear to"
			      " have symmetry %s\n", symmetry_name(sym));
			return 1;
		}

		reflist_free(reflections);

	}

	image.det = get_detector_geometry(geometry);
	if ( image.det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	image.beam = get_beam_parameters(beamfile);
	if ( image.beam == NULL ) {
		ERROR("Failed to read beam parameters from '%s'\n", beamfile);
		return 1;
	}
	free(beamfile);

	/* Define image parameters */
	image.width = image.det->max_fs + 1;
	image.height = image.det->max_ss + 1;
	image.lambda = ph_en_to_lambda(eV_to_J(image.beam->photon_energy));
	image.bw = image.beam->bandwidth;
	image.div = image.beam->divergence;

	/* Load unit cell */
	input_cell = load_cell_from_pdb(filename);
	if ( input_cell == NULL ) {
		exit(1);
	}

	/* Initialise stuff */
	image.filename = NULL;
	image.features = NULL;
	image.flags = NULL;

	powder = calloc(image.width*image.height, sizeof(*powder));

	/* Splurge a few useful numbers */
	STATUS("Wavelength is %f nm\n", image.lambda/1.0e-9);

	do {

		int na, nb, nc;
		double a, b, c, d;
		UnitCell *cell;

		if ( random_size ) {

			double alen, blen, clen, dis;

			cell_get_parameters(input_cell, &alen, &blen, &clen,
			                    &dis, &dis, &dis);

			na = random_ncells(alen, min_size, max_size);
			nb = random_ncells(blen, min_size, max_size);
			nc = random_ncells(clen, min_size, max_size);

		} else {

			na = 8;
			nb = 8;
			nc = 8;

		}

		/* Read quaternion from stdin */
		if ( config_randomquat ) {
			orientation = random_quaternion();
		} else {
			orientation = read_quaternion();
		}

		STATUS("Orientation is %5.3f %5.3f %5.3f %5.3f\n",
		       orientation.w, orientation.x,
		       orientation.y, orientation.z);

		if ( !quaternion_valid(orientation) ) {
			ERROR("Orientation modulus is not zero!\n");
			return 1;
		}

		cell = cell_rotate(input_cell, orientation);

		/* Ensure no residual information */
		image.data = NULL;
		image.twotheta = NULL;

		cell_get_parameters(cell, &a, &b, &c, &d, &d, &d);
		STATUS("Particle size = %i x %i x %i"
		       " ( = %5.2f x %5.2f x %5.2f nm)\n",
	               na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

		if ( config_gpu ) {
			if ( gctx == NULL ) {
				gctx = setup_gpu(config_nosfac,
				                 intensities, flags, sym_str,
				                 gpu_dev);
			}
			get_diffraction_gpu(gctx, &image, na, nb, nc, cell);
		} else {
			get_diffraction(&image, na, nb, nc, intensities, phases,
			                flags, cell, grad, sym);
		}
		if ( image.data == NULL ) {
			ERROR("Diffraction calculation failed.\n");
			goto skip;
		}

		record_image(&image, !config_nonoise);

		if ( powder_fn != NULL ) {

			int x, y, w;

			w = image.width;

			for ( x=0; x<image.width; x++ ) {
			for ( y=0; y<image.height; y++ ) {
				powder[x+w*y] += (double)image.data[x+w*y];
			}
			}

			if ( !(ndone % 10) ) {
				hdf5_write(powder_fn, powder,
				           image.width, image.height,
				           H5T_NATIVE_DOUBLE);
			}
		}

		if ( !config_noimages ) {

			char filename[1024];

			if ( n_images != 1 ) {
				snprintf(filename, 1023, "%s-%i.h5",
				         outfile, number);
			} else {
				strncpy(filename, outfile, 1023);
			}

			number++;

			/* Write the output file */
			hdf5_write_image(filename, &image);

		}

		/* Clean up */
		free(image.data);
		free(image.twotheta);
		cell_free(cell);

skip:
		ndone++;

		if ( n_images && (ndone >= n_images) ) done = 1;

	} while ( !done );

	if ( powder_fn != NULL ) {
		hdf5_write(powder_fn, powder,
		           image.width, image.height,
		           H5T_NATIVE_DOUBLE);
	}

	if ( gctx != NULL ) {
		cleanup_gpu(gctx);
	}

	free(image.det->panels);
	free(image.det);
	free(image.beam);
	free(powder);
	cell_free(input_cell);
	free(intensities);
	free(outfile);
	free(filename);
	free(sym_str);
	free_symoplist(sym);

	return 0;
}
