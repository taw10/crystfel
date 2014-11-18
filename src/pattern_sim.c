/*
 * pattern_sim.c
 *
 * Simulate diffraction patterns from small crystals
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
 *   2013-2014 Chun Hong Yoon <chun.hong.yoon@desy.de>
 *   2014      Valerio Mariani
 *   2013      Alexandra Tolstikova
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
#include "image.h"
#include "diffraction.h"
#include "diffraction-gpu.h"
#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "hdf5-file.h"
#include "detector.h"
#include "peaks.h"
#include "symmetry.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "pattern_sim.h"
#include "stream.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Simulate diffraction patterns from small crystals probed with femtosecond\n"
"pulses of X-rays from a free electron laser.\n"
"\n"
" -h, --help                Display this help message.\n"
"     --version             Print CrystFEL version number and exit.\n"
"\n"
" -p, --pdb=<file>          File from which to get the unit cell.\n"
"                            (The actual Bragg intensities come from the\n"
"                            intensities file)\n"
"     --gpu                 Use the GPU to speed up the calculation.\n"
"     --gpu-dev=<n>         Use GPU device <n>.  Omit this option to see the\n"
"                            available devices.\n"
" -g, --geometry=<file>     Get detector geometry from file.\n"
" -n, --number=<N>          Generate N images.  Default 1.\n"
"     --no-images           Do not output any HDF5 files.\n"
" -o, --output=<filename>   Output HDF5 filename.  Default: sim.h5.\n"
" -r, --random-orientation  Use randomly generated orientations.\n"
"     --powder=<file>       Write a summed pattern of all images simulated by\n"
"                            this invocation as the given filename.\n"
" -i, --intensities=<file>  Specify file containing reflection intensities\n"
"                            (and phases) to use.\n"
" -y, --symmetry=<sym>      The symmetry of the intensities file.\n"
" -t, --gradients=<method>  Select method for calculation of shape transforms\n"
"     --really-random       Seed the random number generator with /dev/urandom.\n"
"     --min-size=<s>        Minimum crystal size in nm.\n"
"     --max-size=<s>        Naximum crystal size in nm.\n"
"     --no-noise            Do not calculate Poisson noise.\n"
" -s, --sample-spectrum=<N> Use N samples from spectrum. Default 3.\n"
" -x, --spectrum=<type>     Type of spectrum to simulate.\n"
"     --background=<N>      Add N photons of Poisson background (default 0).\n"
"     --template=<file>     Take orientations from stream <file>.\n"
"     --no-fringes          Exclude the side maxima of Bragg peaks.\n"
"     --beam-bandwidth      Beam bandwidth as a fraction. Default 1%%.\n"
"     --photon-energy       Photon energy in eV.  Default 9000.\n"
"     --nphotons            Number of photons per X-ray pulse.  Default 1e12.\n"
);
}


static double *intensities_from_list(RefList *list, SymOpList *sym)
{
	Reflection *refl;
	RefListIterator *iter;
	double *out = new_arr_intensity();
	SymOpMask *m = new_symopmask(sym);
	int neq = num_equivs(sym, NULL);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		int eps;
		double intensity = get_intensity(refl);

		get_indices(refl, &h, &k, &l);
		special_position(sym, m, h, k, l);
		eps = neq / num_equivs(sym, m);

		set_arr_intensity(out, h, k, l, intensity / eps);

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


static int random_ncells(double len, double min_m, double max_m)
{
	int min_cells, max_cells, wid;

	min_cells = min_m / len;
	max_cells = max_m / len;
	wid = max_cells - min_cells;

	return wid*random()/RAND_MAX + min_cells;
}


int main(int argc, char *argv[])
{
	int c;
	struct image image;
	struct gpu_context *gctx = NULL;
	struct image *powder;
	float *powder_data;
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
	char *spectrum_str = NULL;
	GradientMethod grad;
	SpectrumType spectrum_type;
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
	int nsamples = 3;
	gsl_rng *rng;
	double background = 0.0;
	char *template_file = NULL;
	Stream *st = NULL;
	int no_fringes = 0;
	double nphotons = 1e12;
	double beam_radius = 1e-6;  /* metres */
	double bandwidth = 0.01;
	double photon_energy = 9000.0;
	struct beam_params beam;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               'v'},
		{"gpu",                0, &config_gpu,         1},
		{"beam",               1, NULL,               'b'},
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
		{"sample-spectrum",    1, NULL,               's'},
		{"type-spectrum",      1, NULL,               'x'},
		{"spectrum",           1, NULL,               'x'},
		{"really-random",      0, &config_random,      1},
		{"no-fringes",         0, &no_fringes,         1},

		{"gpu-dev",            1, NULL,                2},
		{"min-size",           1, NULL,                3},
		{"max-size",           1, NULL,                4},
		{"background",         1, NULL,                5},
		{"template",           1, NULL,                6},
		{"beam-bandwidth",     1, NULL,                7},
		{"photon-energy",      1, NULL,                9},
		{"nphotons",           1, NULL,               10},
		{"beam-radius",        1, NULL,               11},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hrn:i:t:p:o:g:y:s:x:vb:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'v' :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'b' :
			ERROR("WARNING: This version of CrystFEL no longer "
			      "uses beam files.  Please remove the beam file "
			      "from your pattern_sim command line.\n");
			return 1;

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

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 's' :
			nsamples = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid number of spectrum samples.\n");
				return 1;
			}
			break;

			case 'x' :
			spectrum_str = strdup(optarg);
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
			min_size /= 1e9;
			random_size++;
			break;

			case 4 :
			max_size = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid maximum size.\n");
				return 1;
			}
			max_size /= 1e9;
			random_size++;
			break;

			case 5 :
			background = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid background level.\n");
				return 1;
			}
			break;

			case 6 :
			template_file = strdup(optarg);
			break;

			case 7 :
			bandwidth = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid beam bandwidth.\n");
				return 1;
			}
			if ( bandwidth < 0.0 ) {
				ERROR("Beam bandwidth must be positive.\n");
				return 1;
			}
			break;

			case 9 :
			photon_energy = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid photon energy.\n");
				return 1;
			}
			if ( photon_energy < 0.0 ) {
				ERROR("Photon energy must be positive.\n");
				return 1;
			}
			break;

			case 10 :
			nphotons = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid number of photons.\n");
				return 1;
			}
			if ( nphotons < 0.0 ) {
				ERROR("Number of photons must be positive.\n");
				return 1;
			}
			break;

			case 11 :
			beam_radius = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid beam radius.\n");
				return 1;
			}
			if ( beam_radius < 0.0 ) {
				ERROR("Beam radius must be positive.\n");
				return 1;
			}
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

	if ( template_file != NULL ) {
		if ( config_randomquat ) {
			ERROR("You cannot use -r and --template together.\n");
			return 1;
		}
		st = open_stream_for_read(template_file);
		if ( st == NULL ) {
			ERROR("Failed to open stream.\n");
			return 1;
		}
		free(template_file);
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
	image.beam = &beam;
	image.det = get_detector_geometry(geometry, image.beam);
	if ( image.det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);
	if ( (beam.photon_energy > 0.0) && (beam.photon_energy_from == NULL) ) {
		ERROR("WARNING: An explicit photon energy was found in the "
		      "geometry file.  It will be ignored!\n");
		ERROR("The value given on the command line "
		      "(with --photon-energy) will be used instead.\n");
	}

	if ( spectrum_str == NULL ) {
		STATUS("You didn't specify a spectrum type, so"
		       " I'm using a 'tophat' spectrum.\n");
		spectrum_type = SPECTRUM_TOPHAT;
	} else if ( strcasecmp(spectrum_str, "tophat") == 0) {
		spectrum_type = SPECTRUM_TOPHAT;
	} else if ( strcasecmp(spectrum_str, "sase") == 0) {
		spectrum_type = SPECTRUM_SASE;
	} else if ( strcasecmp(spectrum_str, "twocolour") == 0 ||
                    strcasecmp(spectrum_str, "twocolor") == 0 ||
                    strcasecmp(spectrum_str, "twocolours") == 0 ||
                    strcasecmp(spectrum_str, "twocolors") == 0) {
		spectrum_type = SPECTRUM_TWOCOLOUR;
	} else {
		ERROR("Unrecognised spectrum type '%s'\n", spectrum_str);
		return 1;
	}
	free(spectrum_str);

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

		if ( grad == GRADIENT_PHASED ) {
			phases = phases_from_list(reflections);
		} else {
			phases = NULL;
		}
		intensities = intensities_from_list(reflections, sym);
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

	/* Define image parameters */
	image.width = image.det->max_fs + 1;
	image.height = image.det->max_ss + 1;

	double wl = ph_en_to_lambda(eV_to_J(photon_energy));
	image.lambda = wl;
	image.bw = bandwidth;
	image.nsamples = nsamples;

	/* Load unit cell */
	input_cell = load_cell_from_file(filename);
	if ( input_cell == NULL ) {
		exit(1);
	}

	/* Initialise stuff */
	image.filename = NULL;
	image.features = NULL;
	image.flags = NULL;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	if ( config_random ) {
		FILE *fh;
		unsigned long int seed;
		fh = fopen("/dev/urandom", "r");
		fread(&seed, sizeof(seed), 1, fh);
		fclose(fh);
		gsl_rng_set(rng, seed);
	}

	powder = calloc(1, sizeof(struct image));
	powder->width = image.width;
	powder->height = image.height;
	powder->det = image.det;
	powder_data = calloc(image.width*image.height, sizeof(float));
	powder->data = powder_data;

	/* Splurge a few useful numbers */
	STATUS("Simulation parameters:\n");
	STATUS("                  Photon energy: %.2f eV (wavelength %.5f A)\n",
	       photon_energy, image.lambda*1e10);
	STATUS("                Beam divergence: not simulated\n");
	STATUS("                     Background: %.2f photons\n", background);

	switch ( spectrum_type ) {

		case SPECTRUM_TOPHAT:
		STATUS("                 X-ray spectrum: top hat, "
		       "width %.5f %%\n", image.bw*100.0);
		break;

		case SPECTRUM_SASE:
		STATUS("                 X-ray spectrum: SASE, "
		       "bandwidth %.5f %%\n", image.bw*100.0);
		break;

		case SPECTRUM_TWOCOLOUR:
		STATUS("                 X-ray spectrum: two colour, "
		       "separation %.5f %%\n", image.bw*100.0);
		break;
	}
	if ( random_size ) {
		STATUS("                   Crystal size: random, between "
		       "%.2f and %.2f nm along each of a, b and c\n",
		       min_size*1e9, max_size*1e9);
	} else {
		STATUS("                   Crystal size: 8 unit cells along "
		       "each of a, b and c\n");
	}
	if ( intfile == NULL ) {
		STATUS("               Full intensities: all equal");
	} else {
		STATUS("               Full intensities: from %s\n", intfile);
	}
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

		if ( st == NULL ) {

			if ( config_randomquat ) {
				orientation = random_quaternion(rng);
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

		} else {

			struct image image;
			int i;
			Crystal *cr;

			image.det = NULL;

			/* Get data from next chunk */
			if ( read_chunk(st, &image) ) break;

			free(image.filename);
			image_feature_list_free(image.features);

			if ( image.n_crystals == 0 ) continue;

			cr = image.crystals[0];
			cell = crystal_get_cell(cr);

			if ( image.n_crystals > 1 ) {
				ERROR("Using the first crystal only.\n");
			}

			for ( i=1; i<image.n_crystals; i++ ) {

				Crystal *cr = image.crystals[i];
				cell = crystal_get_cell(cr);

				reflist_free(crystal_get_reflections(cr));
				cell_free(crystal_get_cell(cr));
				crystal_free(cr);

			}

			free(image.crystals);

		}

		switch ( spectrum_type ) {

			case SPECTRUM_TOPHAT :
			image.spectrum = generate_tophat(&image);
			break;

			case SPECTRUM_SASE :
			image.spectrum = generate_SASE(&image, rng);
			break;

                        case SPECTRUM_TWOCOLOUR :
                        image.spectrum = generate_twocolour(&image);
                        break;

		}

		/* Ensure no residual information */
		image.data = NULL;
		image.twotheta = NULL;

		cell_get_parameters(cell, &a, &b, &c, &d, &d, &d);
		STATUS("Particle size = %i x %i x %i"
		       " ( = %5.2f x %5.2f x %5.2f nm)\n",
	               na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

		if ( config_gpu ) {

			int err;

			if ( gctx == NULL ) {
				gctx = setup_gpu(config_nosfac,
				                 intensities, flags, sym_str,
				                 gpu_dev);
			}
			err = get_diffraction_gpu(gctx, &image, na, nb, nc,
			                          cell, no_fringes);
			if ( err ) image.data = NULL;

		} else {
			get_diffraction(&image, na, nb, nc, intensities, phases,
			                flags, cell, grad, sym, no_fringes);
		}
		if ( image.data == NULL ) {
			ERROR("Diffraction calculation failed.\n");
			done = 1;
			goto skip;
		}

		record_image(&image, !config_nonoise, background, rng,
		             beam_radius, nphotons);

		if ( powder_fn != NULL ) {

			int x, y, w;

			w = image.width;

			for ( x=0; x<image.width; x++ ) {
			for ( y=0; y<image.height; y++ ) {
				powder->data[x+w*y] += (double)image.data[x+w*y];
			}
			}

			if ( !(ndone % 10) ) {
				hdf5_write_image(powder_fn, powder, NULL);
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
			hdf5_write_image(filename, &image, NULL);

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
		hdf5_write_image(powder_fn, powder, NULL);
	}

	if ( gctx != NULL ) {
		cleanup_gpu(gctx);
	}

	free(image.det->panels);
	free(intfile);
	free(image.det);
	free(powder->data);
	free(powder);
	cell_free(input_cell);
	free(intensities);
	free(outfile);
	free(filename);
	free(sym_str);
	free_symoplist(sym);

	gsl_rng_free(rng);

	return 0;
}
