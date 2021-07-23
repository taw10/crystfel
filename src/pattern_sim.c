/*
 * pattern_sim.c
 *
 * Simulate diffraction patterns from small crystals
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
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
#include <hdf5.h>

#include <image.h>
#include <cell.h>
#include <cell-utils.h>
#include <utils.h>
#include <peaks.h>
#include <symmetry.h>
#include <reflist.h>
#include <reflist-utils.h>
#include <stream.h>
#include <datatemplate.h>

#include "diffraction.h"
#include "diffraction-gpu.h"
#include "pattern_sim.h"
#include "version.h"


enum spectrum_type {
	SPECTRUM_TOPHAT,    /**< A top hat distribution */
	SPECTRUM_GAUSSIAN,  /**< A Gaussian distribution */
	SPECTRUM_SASE,      /**< A simulated SASE spectrum */
	SPECTRUM_TWOCOLOUR, /**< A spectrum containing two peaks */
	SPECTRUM_FROMFILE   /**< An arbitrary spectrum read from a file */
};

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
"     --max-size=<s>        Maximum crystal size in nm.\n"
"     --no-noise            Do not calculate Poisson noise.\n"
" -s, --sample-spectrum=<N> Use N samples from spectrum. Default 3.\n"
" -x, --spectrum=<type>     Type of spectrum to simulate.\n"
"     --background=<N>      Add N photons of Poisson background (default 0).\n"
"     --template=<file>     Take orientations from stream <file>.\n"
"     --no-fringes          Exclude the side maxima of Bragg peaks.\n"
"     --flat                Make Bragg peaks flat.\n"
"     --beam-bandwidth      Beam bandwidth (standard deviation of wavelength as\n"
"                            a fraction of wavelenth).  Default 0.001 (1%%)\n"
"     --sase-spike-width    SASE spike width (standard deviation in m^-1)\n"
"                            Default 2e6 m^-1\n"
"     --twocol-separation   Separation between peaks in two-colour mode in m^-1\n"
"                            Default 8e6 m^-1\n"
"     --photon-energy       Photon energy in eV.  Default 9000.\n"
"     --nphotons            Number of photons per X-ray pulse.  Default 1e12.\n"
"     --beam-radius         Radius of beam in metres (default 1e-6).\n"
);
}


static void record_panel(struct detgeom_panel *p, float *dp,
                         int do_poisson,
                         gsl_rng *rng, double ph_per_e, double background,
			 double lambda,
                         int *n_neg1, int *n_inf1, int *n_nan1,
                         int *n_neg2, int *n_inf2, int *n_nan2,
			 double *max_tt)
{
	int fs, ss;

	for ( ss=0; ss<p->h; ss++ ) {
	for ( fs=0; fs<p->w; fs++ ) {

		double counts;
		double intensity, sa;
		double pix_area, Lsq;
		double xs, ys, rx, ry;
		double dsq, proj_area;
		float dval;
		double twotheta;

		intensity = (double)dp[fs + p->w*ss];
		if ( isinf(intensity) ) (*n_inf1)++;
		if ( intensity < 0.0 ) (*n_neg1)++;
		if ( isnan(intensity) ) (*n_nan1)++;

		/* Area of one pixel */
		pix_area = pow(p->pixel_pitch, 2.0);
		Lsq = pow(p->cnz*p->pixel_pitch, 2.0);

		/* Calculate distance from crystal to pixel */
		xs = fs*p->fsx + ss*p->ssx;
		ys = ss*p->fsy + ss*p->ssy;
		rx = (xs + p->cnx) * p->pixel_pitch;
		ry = (ys + p->cny) * p->pixel_pitch;
		dsq = pow(rx, 2.0) + pow(ry, 2.0);
		twotheta = atan2(sqrt(dsq), p->cnz*p->pixel_pitch);

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(twotheta);

		/* Projected area of pixel divided by distance squared */
		sa = proj_area / (dsq + Lsq);

		if ( do_poisson ) {
			counts = poisson_noise(rng, intensity * ph_per_e * sa);
		} else {
			counts = intensity * ph_per_e * sa;
		}

		/* Number of photons in pixel */
		dval = counts + poisson_noise(rng, background);

		/* Convert to ADU */
		dval *= p->adu_per_photon;

		/* Saturation */
		if ( dval > p->max_adu ) dval = p->max_adu;

		dp[fs + p->w*ss] = dval;

		/* Sanity checks */
		if ( isinf(dp[fs + p->w*ss]) ) n_inf2++;
		if ( isnan(dp[fs + p->w*ss]) ) n_nan2++;
		if ( dp[fs + p->w*ss] < 0.0 ) n_neg2++;

		if ( twotheta > *max_tt ) *max_tt = twotheta;

	}
	}
}


void record_image(struct image *image, int do_poisson, double background,
                  gsl_rng *rng, double beam_radius, double nphotons)
{
	double total_energy, energy_density;
	double ph_per_e;
	double area;
	double max_tt = 0.0;
	int n_inf1 = 0;
	int n_neg1 = 0;
	int n_nan1 = 0;
	int n_inf2 = 0;
	int n_neg2 = 0;
	int n_nan2 = 0;
	int pn;

	/* How many photons are scattered per electron? */
	area = M_PI*pow(beam_radius, 2.0);
	total_energy = nphotons * ph_lambda_to_en(image->lambda);
	energy_density = total_energy / area;
	ph_per_e = (nphotons /area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       nphotons, energy_density/1e7, total_energy*1e6);

	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {

		record_panel(&image->detgeom->panels[pn], image->dp[pn],
		             do_poisson, rng, ph_per_e, background,
			     image->lambda,
		             &n_neg1, &n_inf1, &n_nan1,
		             &n_neg2, &n_inf2, &n_nan2, &max_tt);
	}

	STATUS("Max 2theta = %.2f deg, min d = %.2f nm\n",
	        rad2deg(max_tt), (image->lambda/(2.0*sin(max_tt/2.0)))/1e-9);

	STATUS("Halve the d values to get the voxel size for a synthesis.\n");

	if ( n_neg1 + n_inf1 + n_nan1 + n_neg2 + n_inf2 + n_nan2 ) {
		ERROR("WARNING: The raw calculation produced %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg1, n_inf1, n_nan1);
		ERROR("WARNING: After processing, there were %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg2, n_inf2, n_nan2);
	}
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


static void add_metadata(const char *filename,
                         struct quaternion q,
                         UnitCell *cell)
{
	hid_t fh, gh, sh, dh;
	herr_t r;
	hsize_t size[1];
	float data[4];
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	fh = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Failed to open file to add metadata.\n");
		return;
	}

	gh = H5Gcreate2(fh, "pattern_sim", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Failed to create metadata group.\n");
		H5Fclose(fh);
		return;
	}

	size[0] = 4;
	sh = H5Screate_simple(1, size, NULL);

	dh = H5Dcreate2(gh, "quaternion", H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return;
	}

	data[0] = q.w;
	data[1] = q.x;
	data[2] = q.y;
	data[3] = q.z;
	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Failed to write quaternion to file\n");
	}
	H5Dclose(dh);

	size[0] = 3;
	sh = H5Screate_simple(1, size, NULL);

	dh = H5Dcreate2(gh, "astar", H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return;
	}
	data[0] = asx;
	data[1] = asy;
	data[2] = asz;
	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Failed to write astar to file.\n");
	}
	H5Dclose(dh);

	dh = H5Dcreate2(gh, "bstar", H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return;
	}
	data[0] = bsx;
	data[1] = bsy;
	data[2] = bsz;
	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Failed to write bstar to file.\n");
	}
	H5Dclose(dh);

	dh = H5Dcreate2(gh, "cstar", H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return;
	}
	data[0] = csx;
	data[1] = csy;
	data[2] = csz;
	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Failed to write cstar to file.\n");
	}
	H5Dclose(dh);

	H5Gclose(gh);
	H5Fclose(fh);
}


int main(int argc, char *argv[])
{
	int argn;
	struct image *image;
	DataTemplate *dtempl;
	struct gpu_context *gctx = NULL;
	struct image *powder;
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
	char *cell_filename = NULL;
	char *grad_str = NULL;
	char *outfile = NULL;
	char *geometry = NULL;
	char *spectrum_str = NULL;
	char *spectrum_fn = NULL;
	GradientMethod grad;
	enum spectrum_type spectrum_type;
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
	int flat = 0;
	double nphotons = 1e12;
	double beam_radius = 1e-6;  /* metres */
	double bandwidth = 0.01;
	double photon_energy = 9000.0;
	double sase_spike_width = 2e6;  /* m^-1 */
	double twocol_sep = 8e6;  /* m^-1 */

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
		{"flat",               0, &flat,               1},

		{"gpu-dev",            1, NULL,                2},
		{"min-size",           1, NULL,                3},
		{"max-size",           1, NULL,                4},
		{"background",         1, NULL,                5},
		{"template",           1, NULL,                6},
		{"beam-bandwidth",     1, NULL,                7},
		{"photon-energy",      1, NULL,                9},
		{"nphotons",           1, NULL,               10},
		{"beam-radius",        1, NULL,               11},
		{"spectrum-file",      1, NULL,               12},
		{"sase-spike-width",   1, NULL,               13},
		{"twocol-separation",  1, NULL,               14},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((argn = getopt_long(argc, argv, "hrn:i:t:p:o:g:y:s:x:vb:",
	                        longopts, NULL)) != -1) {

		switch (argn) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'v' :
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
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
			cell_filename = strdup(optarg);
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

			case 12 :
			spectrum_fn = strdup(optarg);
			break;

			case 13 :
			sase_spike_width = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid SASE spike width.\n");
				return 1;
			}
			if ( beam_radius < 0.0 ) {
				ERROR("SASE spike width must be positive.\n");
				return 1;
			}
			break;

			case 14 :
			twocol_sep = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid two colour separation.\n");
				return 1;
			}
			if ( beam_radius < 0.0 ) {
				ERROR("Two-colour separation must be positive.\n");
				return 1;
			}
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", argn);
			break;

		}

	}

	if ( random_size == 1 ) {
		ERROR("You must specify both --min-size and --max-size.\n");
		return 1;
	}

	if ( cell_filename == NULL ) {
		cell_filename = strdup("molecule.pdb");
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
		st = stream_open_for_read(template_file);
		if ( st == NULL ) {
			ERROR("Failed to open stream.\n");
			return 1;
		}
		free(template_file);
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
	dtempl = data_template_new_from_file(geometry);
	if ( dtempl == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n",
		      geometry);
		return 1;
	}
	free(geometry);

	if ( spectrum_str == NULL ) {
		STATUS("You didn't specify a spectrum type, so"
		       " I'm using a 'tophat' spectrum.\n");
		spectrum_type = SPECTRUM_TOPHAT;
	} else if ( strcasecmp(spectrum_str, "tophat") == 0) {
		spectrum_type = SPECTRUM_TOPHAT;
	} else if ( strcasecmp(spectrum_str, "gaussian") == 0) {
		spectrum_type = SPECTRUM_GAUSSIAN;
	} else if ( strcasecmp(spectrum_str, "sase") == 0) {
		spectrum_type = SPECTRUM_SASE;
	} else if ( strcasecmp(spectrum_str, "twocolour") == 0 ||
                    strcasecmp(spectrum_str, "twocolor") == 0 ||
                    strcasecmp(spectrum_str, "twocolours") == 0 ||
                    strcasecmp(spectrum_str, "twocolors") == 0) {
		spectrum_type = SPECTRUM_TWOCOLOUR;
	} else if ( strcasecmp(spectrum_str, "fromfile") == 0) {
		spectrum_type = SPECTRUM_FROMFILE;
	} else {
		ERROR("Unrecognised spectrum type '%s'\n", spectrum_str);
		return 1;
	}
	free(spectrum_str);

	/* Load unit cell */
	input_cell = load_cell_from_file(cell_filename);
	if ( input_cell == NULL ) {
		exit(1);
	}

	if ( intfile == NULL ) {

		/* Gentle reminder */
		STATUS("You didn't specify the file containing the ");
		STATUS("reflection intensities (with --intensities).\n");
		STATUS("I'll simulate a flat intensity distribution.\n");
		intensities = NULL;
		phases = NULL;
		flags = NULL;

		if ( sym_str == NULL ) sym_str = strdup("1");
		pointgroup_warning(sym_str);
		sym = get_pointgroup(sym_str);

	} else {

		RefList *reflections;
		char *sym_str_fromfile = NULL;

		reflections = read_reflections_2(intfile, &sym_str_fromfile);
		if ( reflections == NULL ) {
			ERROR("Problem reading input file %s\n", intfile);
			return 1;
		}

		/* If we don't have a point group yet, and if the file provides
		 * one, use the one from the file */
		if ( (sym_str == NULL) && (sym_str_fromfile != NULL) ) {
			sym_str = sym_str_fromfile;
			STATUS("Using symmetry from reflection file: %s\n",
			       sym_str);
		}

		/* If we still don't have a point group, use "1" */
		if ( sym_str == NULL ) sym_str = strdup("1");

		pointgroup_warning(sym_str);
		sym = get_pointgroup(sym_str);

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
			if ( cell_get_lattice_type(input_cell) == L_MONOCLINIC )
			{
				ERROR("You may need to specify the unique axis "
				      "in your point group.  The default is "
				      "unique axis c.\n");
				ERROR("See 'man crystfel' for more details.\n");
			}
			return 1;
		}

		reflist_free(reflections);

	}

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	if ( config_random ) {
		FILE *fh;
		unsigned long int seed;
		fh = fopen("/dev/urandom", "r");
		if ( fread(&seed, sizeof(seed), 1, fh) == 1 ) {
			gsl_rng_set(rng, seed);
		} else {
			ERROR("Failed to seed random number generator\n");
		}
		fclose(fh);
	}

	image = image_create_for_simulation(dtempl);
	powder = image_create_for_simulation(dtempl);

	/* Splurge a few useful numbers */
	STATUS("Simulation parameters:\n");
	STATUS("                     Wavelength: %.5f A (photon energy %.2f eV)\n",
	       image->lambda*1e10,
	       ph_lambda_to_eV(image->lambda));
	STATUS("              Number of photons: %.0f (%.2f mJ)\n",
	       nphotons,
	       eV_to_J(photon_energy)*nphotons*1e3);
	STATUS("                Beam divergence: not simulated\n");
	STATUS("                    Beam radius: %.2f microns\n",
	      beam_radius*1e6);
	STATUS("                     Background: %.2f photons\n", background);

	switch ( spectrum_type ) {

		case SPECTRUM_TOPHAT:
		STATUS("                 X-ray spectrum: top hat, "
		       "width %.5f %%\n", image->bw*100.0);
		break;

		case SPECTRUM_GAUSSIAN:
		STATUS("                 X-ray spectrum: Gaussian, "
		       "bandwidth %.5f %%\n", image->bw*100.0);
		break;

		case SPECTRUM_SASE:
		STATUS("                 X-ray spectrum: SASE, "
		       "bandwidth %.5f %%\n", image->bw*100.0);
		break;

		case SPECTRUM_TWOCOLOUR:
		STATUS("                 X-ray spectrum: two colour, "
		       "separation %.5f %%\n", image->bw*100.0);
		break;

		case SPECTRUM_FROMFILE:
		STATUS("                 X-ray spectrum: from %s\n",
			spectrum_fn);
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
		STATUS("               Full intensities: all equal\n");
	} else {
		STATUS("               Full intensities: from %s (symmetry %s)\n",
		       intfile, sym_str);
	}
	do {

		int na, nb, nc;
		double a, b, c, dis;
		UnitCell *cell;
		int err = 0;
		int pi;

		/* Reset image data to zero for new pattern */
		for ( pi=0; pi<image->detgeom->n_panels; pi++ ) {
			long j;
			long np = image->detgeom->panels[pi].w
			            * image->detgeom->panels[pi].h;
			for ( j=0; j<np; j++ ) {
				image->dp[pi][j] = 0.0;
			}
		}

		if ( random_size ) {

			double alen, blen, clen, discard;

			cell_get_parameters(input_cell, &alen, &blen, &clen,
			                    &discard, &discard, &discard);

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

			struct image *templ_image;
			Crystal *cr;

			/* Get data from next chunk */
			templ_image = stream_read_chunk(st, 0);
			if ( templ_image == NULL ) break;
			if ( templ_image->n_crystals == 0 ) continue;

			cr = templ_image->crystals[0];
			cell = cell_new_from_cell(crystal_get_cell(cr));

			if ( templ_image->n_crystals > 1 ) {
				ERROR("Using the first crystal only.\n");
			}

			image_free(templ_image);

		}

		switch ( spectrum_type ) {

			case SPECTRUM_TOPHAT :
			image->spectrum = spectrum_generate_tophat(image->lambda,
			                                           image->bw);
			break;

			case SPECTRUM_GAUSSIAN :
			image->spectrum = spectrum_generate_gaussian(image->lambda,
			                                             image->bw);
			break;


			case SPECTRUM_SASE :
			image->spectrum = spectrum_generate_sase(image->lambda,
			                                         image->bw,
			                                         sase_spike_width,
			                                         rng);
			break;

			case SPECTRUM_TWOCOLOUR :
			image->spectrum = spectrum_generate_twocolour(image->lambda,
			                                              image->bw,
			                                              twocol_sep);
                        break;

			case SPECTRUM_FROMFILE :
			image->spectrum = spectrum_load(spectrum_fn);
                        break;

		}

		cell_get_parameters(cell, &a, &b, &c, &dis, &dis, &dis);
		STATUS("Particle size = %i x %i x %i"
		       " ( = %5.2f x %5.2f x %5.2f nm)\n",
	               na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

		if ( config_gpu ) {

			if ( gctx == NULL ) {
				gctx = setup_gpu(config_nosfac,
				                 intensities, flags, sym_str,
				                 gpu_dev);
			}
			err = get_diffraction_gpu(gctx, image, na, nb, nc,
			                          cell, no_fringes, flat,
			                          nsamples);

		} else {
			get_diffraction(image, na, nb, nc, intensities, phases,
			                flags, cell, grad, sym, no_fringes, flat,
			                nsamples);
		}
		if ( err ) {
			ERROR("Diffraction calculation failed.\n");
			done = 1;
			goto skip;
		}

		record_image(image, !config_nonoise, background, rng,
		             beam_radius, nphotons);

		if ( powder_fn != NULL ) {

			int pn;

			for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {

				size_t w, i;

				w = image->detgeom->panels[pn].w
				  * image->detgeom->panels[pn].h;

				for ( i=0; i<w; i++ ) {
					powder->dp[pn][i] += (double)image->dp[pn][i];
				}

			}

			if ( !(ndone % 10) ) {
				image_write(powder, dtempl, powder_fn);
			}
		}

		if ( !config_noimages ) {

			char h5filename[1024];

			if ( n_images != 1 ) {
				snprintf(h5filename, 1023, "%s-%i.h5",
				         outfile, number);
			} else {
				strncpy(h5filename, outfile, 1023);
			}

			number++;

			/* Write the output file */
			image_write(image, dtempl, h5filename);

			/* Add some pattern_sim-specific metadata */
			add_metadata(h5filename, orientation, cell);

		}

		cell_free(cell);

skip:
		ndone++;

		if ( n_images && (ndone >= n_images) ) done = 1;

	} while ( !done );

	if ( powder_fn != NULL ) {
		image_write(powder, dtempl, powder_fn);
	}

	if ( gctx != NULL ) {
		cleanup_gpu(gctx);
	}

	image_free(image);
	image_free(powder);
	free(intfile);
	cell_free(input_cell);
	free(intensities);
	free(outfile);
	free(cell_filename);
	free(sym_str);
	free_symoplist(sym);

	gsl_rng_free(rng);

	return 0;
}
