/*
 * pattern_sim.c
 *
 * Simulate diffraction patterns from small crystals
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include "sfac.h"
#include "reflections.h"
#include "beam-parameters.h"
#include "symmetry.h"


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
"     --simulation-details  Show technical details of the simulation.\n"
"     --gpu                 Use the GPU to speed up the calculation.\n"
"     --gpu-dev=<n>         Use GPU device <n>.  Omit this option to see the\n"
"                            available devices.\n"
" -g, --geometry=<file>     Get detector geometry from file.\n"
" -b, --beam=<file>         Get beam parameters from file.\n"
"\n"
"     --near-bragg          Output h,k,l,I near Bragg conditions.\n"
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
"\n"
"By default, the simulation aims to be as accurate as possible.  For greater\n"
"speed, or for testing, you can choose to disable certain things using the\n"
"following options.\n"
"\n"
"     --no-noise            Do not calculate Poisson noise.\n"
);
}


static void show_details()
{
	printf(
"This program simulates diffraction patterns from small crystals illuminated\n"
"with femtosecond X-ray pulses from a free electron laser.\n"
"\n"
"The lattice transform from the specified number of unit cells is calculated\n"
"using the closed-form solution for a truncated lattice faceted on the\n"
"(001), (010) and (100) planes:\n"
"\n"
"I_latt(q) =  sin^2(pi*na*q.a)/sin^2(pi*q.a)\n"
"           * sin^2(pi*nb*q.b)/sin^2(pi*q.b)\n"
"           * sin^2(pi*nc*q.c)/sin^2(pi*q.c)\n"
"\n"
"na = number of unit cells in 'a' direction (likewise nb, nc)\n"
" q = reciprocal vector (1/d convention, not 2pi/d)\n"
"\n"
"This is multiplied by a model of the underlying molecular transform, I_mol(q).\n"
"This can be approximated to varying levels of accuracy by the methods given by\n"
"the '--gradients' option.\n"
"\n"
"Expected intensities at the CCD are then calculated using:\n"
"\n"
"I(q) = I0 * r^2 * I_latt(q) * I_mol(q) * S\n"
"\n"
"I0 = number of photons per unit area in the incident beam\n"
" r = Thomson radius\n"
" S = solid angle of corresponding pixel\n"
"\n"
"Polarisation is not currently included in pattern_sim, although it is included\n"
"in the analysis of Bragg peaks by 'indexamajig'.\n"
"\n"
"Poisson counts are generated from the expected intensities using Knuth's\n"
"algorithm.  When the intensity is sufficiently high that Knuth's algorithm\n"
"would result in machine precision problems, a normal distribution with\n"
"standard deviation sqrt(I) is used instead.\n"
);
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
	int config_simdetails = 0;
	int config_nearbragg = 0;
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
	char *sym = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"simulation-details", 0, &config_simdetails,  1},
		{"gpu",                0, &config_gpu,         1},
		{"near-bragg",         0, &config_nearbragg,   1},
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
			sym = strdup(optarg);
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

	if ( sym == NULL ) sym = strdup("1");

	if ( config_simdetails ) {
		show_details();
		return 0;
	}

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

		int i;
		ReflItemList *items;

		if ( grad == GRADIENT_PHASED ) {
			phases = new_list_phase();
		} else {
			phases = NULL;
		}
		intensities = new_list_intensity();
		phases = new_list_phase();
		flags = new_list_flag();
		items = read_reflections(intfile, intensities, phases,
		                         NULL, NULL);
		free(intfile);

		for ( i=0; i<num_items(items); i++ ) {
			struct refl_item *it = get_item(items, i);
			set_flag(flags, it->h, it->k, it->l, 1);
		}

		/* Check that the intensities have the correct symmetry */
		if ( check_symmetry(items, sym) ) {
			ERROR("The input reflection list does not appear to"
			      " have symmetry %s\n", sym);
			return 1;
		}

		delete_items(items);
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

	/* Load unit cell */
	input_cell = load_cell_from_pdb(filename);
	if ( input_cell == NULL ) {
		exit(1);
	}

	/* Initialise stuff */
	image.filename = NULL;
	image.features = NULL;
	image.flags = NULL;
	image.i0 = 1.0;
	image.i0_available = 1;

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
		STATUS("Particle size = %i x %i x %i (=%5.2f x %5.2f x %5.2f nm)\n",
	               na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

		if ( config_gpu ) {
			if ( gctx == NULL ) {
				gctx = setup_gpu(config_nosfac,
				                 intensities, flags, sym,
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

		if ( config_nearbragg ) {

			/* Do EITHER: */

			//image.div = beam->divergence;
			//image.bw = beam->bandwidth;
			//image.profile_radius = 0.0001e9;
			//image.reflections = find_intersections(&image,
			//                               image.indexed_cell, 0);

			image.reflections = find_projected_peaks(&image,
			                            image.indexed_cell, 0, 0.1);

			integrate_reflections(&image, 0, 0);

			/* OR */

			//image.reflections = integrate_pixels(&image, 0, 0.1,
			//                                     config_polar);

		}

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
			hdf5_write(filename, image.data,
			           image.width, image.height, H5T_NATIVE_FLOAT);

		}

		/* Clean up */
		free(image.data);
		free(image.twotheta);
		cell_free(cell);

skip:
		ndone++;

		if ( n_images && (ndone >= n_images) ) done = 1;

	} while ( !done );

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
	free(sym);

	return 0;
}
