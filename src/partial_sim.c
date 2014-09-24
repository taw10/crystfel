/*
 * partial_sim.c
 *
 * Generate partials for testing scaling
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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
#include <assert.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>

#include "version.h"
#include "image.h"
#include "utils.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "geometry.h"
#include "stream.h"
#include "thread-pool.h"
#include "cell-utils.h"

/* Number of bins for partiality graph */
#define NBINS 50


static void mess_up_cell(Crystal *cr, double cnoise, gsl_rng *rng)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	UnitCell *cell = crystal_get_cell(cr);

	//STATUS("Real:\n");
	//cell_print(cell);

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	ax = flat_noise(rng, ax, cnoise*fabs(ax)/100.0);
	ay = flat_noise(rng, ay, cnoise*fabs(ay)/100.0);
	az = flat_noise(rng, az, cnoise*fabs(az)/100.0);
	bx = flat_noise(rng, bx, cnoise*fabs(bx)/100.0);
	by = flat_noise(rng, by, cnoise*fabs(by)/100.0);
	bz = flat_noise(rng, bz, cnoise*fabs(bz)/100.0);
	cx = flat_noise(rng, cx, cnoise*fabs(cx)/100.0);
	cy = flat_noise(rng, cy, cnoise*fabs(cy)/100.0);
	cz = flat_noise(rng, cz, cnoise*fabs(cz)/100.0);
	cell_set_reciprocal(cell, ax, ay, az, bx, by, bz, cx, cy, cz);

	//STATUS("Changed:\n");
	//cell_print(cell);
}


/* For each reflection in "partial", fill in what the intensity would be
 * according to "full" */
static void calculate_partials(Crystal *cr,
                               RefList *full, const SymOpList *sym,
                               int random_intensities,
                               pthread_rwlock_t *full_lock,
                               unsigned long int *n_ref, double *p_hist,
                               double *p_max, double max_q, double full_stddev,
                               double noise_stddev, gsl_rng *rng)
{
	Reflection *refl;
	RefListIterator *iter;
	double res;

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		Reflection *rfull;
		double L, p, Ip, If;
		int bin;

		get_indices(refl, &h, &k, &l);
		get_asymm(sym, h, k, l, &h, &k, &l);
		p = get_partiality(refl);
		L = get_lorentz(refl);

		pthread_rwlock_rdlock(full_lock);
		rfull = find_refl(full, h, k, l);
		pthread_rwlock_unlock(full_lock);

		if ( rfull == NULL ) {
			if ( random_intensities ) {

				pthread_rwlock_wrlock(full_lock);

				/* In the gap between the unlock and the wrlock,
				 * the reflection might have been created by
				 * another thread.  So, we must check again */
				rfull = find_refl(full, h, k, l);
				if ( rfull == NULL ) {
					rfull = add_refl(full, h, k, l);
					If = fabs(gaussian_noise(rng, 0.0,
					                         full_stddev));
					set_intensity(rfull, If);
					set_redundancy(rfull, 1);
				} else {
					If = get_intensity(rfull);
				}
				pthread_rwlock_unlock(full_lock);

			} else {
				set_redundancy(refl, 0);
				If = 0.0;
			}
		} else {
			If = get_intensity(rfull);
			if ( random_intensities ) {
				lock_reflection(rfull);
				int red = get_redundancy(rfull);
				set_redundancy(rfull, red+1);
				unlock_reflection(rfull);
			}
		}

		Ip = crystal_get_osf(cr) * L * p * If;

		res = resolution(crystal_get_cell(cr), h, k, l);
		bin = NBINS*2.0*res/max_q;
		if ( (bin < NBINS) && (bin>=0) ) {
			p_hist[bin] += p;
			n_ref[bin]++;
			if ( p > p_max[bin] ) p_max[bin] = p;
		} else {
			STATUS("Reflection out of histogram range: %e %i %f\n",
			       res, bin,  p);
		}

		Ip = gaussian_noise(rng, Ip, noise_stddev);

		set_intensity(refl, Ip);
		set_esd_intensity(refl, noise_stddev);
	}
}


static void draw_and_write_image(struct image *image, RefList *reflections,
                                 gsl_rng *rng, double background)
{
	Reflection *refl;
	RefListIterator *iter;
	int i;

	image->data = calloc(image->width*image->height, sizeof(float));

	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double Ip;
		double dfs, dss;
		int fs, ss;

		Ip = get_intensity(refl);

		get_detector_pos(refl, &dfs, &dss);
		fs = nearbyint(dfs);
		ss = nearbyint(dss);
		assert(fs >= 0);
		assert(ss >= 0);
		assert(fs < image->width);
		assert(ss < image->height);
		image->data[fs + image->width*ss] += Ip;

	}

	for ( i=0; i<image->width*image->height; i++ ) {
		image->data[i] += poisson_noise(rng, background);
	}

	hdf5_write_image(image->filename, image, NULL);
	free(image->data);
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Generate a stream containing partials from a reflection list.\n"
"\n"
" -h, --help               Display this help message.\n"
"     --version            Print CrystFEL version number and exit.\n"
"\n"
"You need to provide the following basic options:\n"
" -i, --input=<file>       Read reflections from <file>.\n"
"                           Default: generate random ones instead (see -r).\n"
" -o, --output=<file>      Write partials in stream format to <file>.\n"
"     --images=<prefix>    Write images to <prefix>NNN.h5.\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -p, --pdb=<file>         PDB file from which to get the unit cell.\n"
"\n"
" -y, --symmetry=<sym>     Symmetry of the input reflection list.\n"
" -n <n>                   Simulate <n> patterns.  Default: 2.\n"
" -r, --save-random=<file> Save randomly generated intensities to file.\n"
"     --pgraph=<file>      Save a histogram of partiality values to file.\n"
" -c, --cnoise=<val>       Amount of reciprocal space cell noise, in percent.\n"
"     --osf-stddev=<val>   Standard deviation of the scaling factors.\n"
"     --full-stddev=<val>  Standard deviation of the randomly\n"
"                           generated full intensities, if not using -i.\n"
"     --noise-stddev=<val> Set the standard deviation of the noise.\n"
"     --background=<val>   Background level in photons.  Default 3000.\n"
"     --beam-divergence    Beam divergence in radians. Default 1 mrad.\n"
"     --beam-bandwidth     Beam bandwidth as a fraction. Default 1%%.\n"
"     --profile-radius     Reciprocal space reflection profile radius in m^-1.\n"
"                           Default 0.001e9 m^-1\n"
"     --photon-energy      Photon energy in eV.  Default 9000.\n"
"\n"
);
}


struct queue_args
{
	RefList *full;
	pthread_rwlock_t full_lock;

	int n_done;
	int n_started;
	int n_to_do;

	SymOpList *sym;
	int random_intensities;
	UnitCell *cell;
	double cnoise;
	double osf_stddev;
	double full_stddev;
	double noise_stddev;
	double background;
	double profile_radius;

	struct image *template_image;
	double max_q;

	char *image_prefix;

	/* The overall histogram */
	double p_hist[NBINS];
	unsigned long int n_ref[NBINS];
	double p_max[NBINS];

	Stream *stream;
	gsl_rng **rngs;
};


struct worker_args
{
	struct queue_args *qargs;
	Crystal *crystal;
	struct image image;

	/* Histogram for this image */
	double p_hist[NBINS];
	unsigned long int n_ref[NBINS];
	double p_max[NBINS];

	int n;
};


static void *create_job(void *vqargs)
{
	struct worker_args *wargs;
	struct queue_args *qargs = vqargs;

	/* All done already? */
	if ( qargs->n_started == qargs->n_to_do ) return NULL;

	wargs = malloc(sizeof(struct worker_args));

	wargs->qargs = qargs;
	wargs->image = *qargs->template_image;

	qargs->n_started++;
	wargs->n = qargs->n_started;

	return wargs;
}


static void run_job(void *vwargs, int cookie)
{
	struct quaternion orientation;
	struct worker_args *wargs = vwargs;
	struct queue_args *qargs = wargs->qargs;
	int i;
	Crystal *cr;
	RefList *reflections;
	double osf;

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to create crystal.\n");
		return;
	}
	wargs->crystal = cr;
	crystal_set_image(cr, &wargs->image);

	do {
		osf = gaussian_noise(qargs->rngs[cookie], 1.0,
		                     qargs->osf_stddev);
	} while ( osf <= 0.0 );
	crystal_set_osf(cr, osf);
	crystal_set_mosaicity(cr, 0.0);
	crystal_set_profile_radius(cr, qargs->profile_radius);

	/* Set up a random orientation */
	orientation = random_quaternion(qargs->rngs[cookie]);
	crystal_set_cell(cr, cell_rotate(qargs->cell, orientation));

	wargs->image.filename = malloc(256);
	if ( wargs->image.filename == NULL ) {
		ERROR("Failed to allocate filename.\n");
		return;
	}
	if ( qargs->image_prefix != NULL ) {
		snprintf(wargs->image.filename, 255, "%s%i.h5",
		         qargs->image_prefix, wargs->n);
	} else {
		snprintf(wargs->image.filename, 255, "dummy.h5");
	}

	reflections = find_intersections(&wargs->image, cr, PMODEL_SPHERE);
	crystal_set_reflections(cr, reflections);

	for ( i=0; i<NBINS; i++ ) {
		wargs->n_ref[i] = 0;
		wargs->p_hist[i] = 0.0;
		wargs->p_max[i] = 0.0;
	}

	calculate_partials(cr, qargs->full,
	                   qargs->sym, qargs->random_intensities,
	                   &qargs->full_lock,
	                   wargs->n_ref, wargs->p_hist, wargs->p_max,
	                   qargs->max_q, qargs->full_stddev,
	                   qargs->noise_stddev, qargs->rngs[cookie]);

	if ( qargs->image_prefix != NULL ) {
		draw_and_write_image(&wargs->image, reflections,
		                     qargs->rngs[cookie], qargs->background);
	}

	/* Give a slightly incorrect cell in the stream */
	mess_up_cell(cr, qargs->cnoise, qargs->rngs[cookie]);

	image_add_crystal(&wargs->image, cr);
}


static void finalise_job(void *vqargs, void *vwargs)
{
	struct worker_args *wargs = vwargs;
	struct queue_args *qargs = vqargs;
	int i;

	write_chunk(qargs->stream, &wargs->image, NULL, 0, 1, NULL);

	for ( i=0; i<NBINS; i++ ) {
		qargs->n_ref[i] += wargs->n_ref[i];
		qargs->p_hist[i] += wargs->p_hist[i];
		if ( wargs->p_max[i] > qargs->p_max[i] ) {
			qargs->p_max[i] = wargs->p_max[i];
		}
	}

	qargs->n_done++;
	progress_bar(qargs->n_done, qargs->n_to_do, "Simulating");

	free_all_crystals(&wargs->image);
	free(wargs->image.filename);
	free(wargs);
}


int main(int argc, char *argv[])
{
	int c;
	char *input_file = NULL;
	char *output_file = NULL;
	char *geomfile = NULL;
	char *cellfile = NULL;
	struct detector *det = NULL;
	struct beam_params beam;
	RefList *full = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	UnitCell *cell = NULL;
	Stream *stream;
	int n = 2;
	int random_intensities = 0;
	char *save_file = NULL;
	struct queue_args qargs;
	struct image image;
	int n_threads = 1;
	char *rval;
	int i;
	FILE *fh;
	char *phist_file = NULL;
	gsl_rng *rng_for_seeds;
	int config_random = 0;
	char *image_prefix = NULL;

	/* Default simulation parameters */
	double divergence = 0.001;
	double bandwidth = 0.01;
	double profile_radius = 0.001e9;
	double photon_energy = 9000.0;
	double osf_stddev = 2.0;
	double full_stddev = 1000.0;
	double noise_stddev = 20.0;
	double background = 3000.0;
	double cnoise = 0.0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               'v'},
		{"output",             1, NULL,               'o'},
		{"input",              1, NULL,               'i'},
		{"pdb",                1, NULL,               'p'},
		{"geometry",           1, NULL,               'g'},
		{"symmetry",           1, NULL,               'y'},
		{"save-random",        1, NULL,               'r'},
		{"cnoise",             1, NULL,               'c'},

		{"pgraph",             1, NULL,                2},
		{"osf-stddev",         1, NULL,                3},
		{"full-stddev",        1, NULL,                4},
		{"noise-stddev",       1, NULL,                5},
		{"images",             1, NULL,                6},
		{"background",         1, NULL,                7},
		{"beam-divergence",    1, NULL,                8},
		{"beam-bandwidth",     1, NULL,                9},
		{"profile-radius",     1, NULL,               10},
		{"photon-energy",      1, NULL,               11},

		{"really-random",      0, &config_random,      1},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:p:g:y:n:r:j:c:v",
	                        longopts, NULL)) != -1)
	{
		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'v' :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'i' :
			input_file = strdup(optarg);
			break;

			case 'o' :
			output_file = strdup(optarg);
			break;

			case 'p' :
			cellfile = strdup(optarg);
			break;

			case 'g' :
			geomfile = strdup(optarg);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'n' :
			n = atoi(optarg);
			break;

			case 'r' :
			save_file = strdup(optarg);
			break;

			case 'j' :
			n_threads = atoi(optarg);
			break;

			case 'c' :
			cnoise = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid cell noise value.\n");
				return 1;
			}
			break;

			case 2 :
			phist_file = strdup(optarg);
			break;

			case 3 :
			osf_stddev = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid OSF standard deviation.\n");
				return 1;
			}
			if ( osf_stddev < 0.0 ) {
				ERROR("Invalid OSF standard deviation.");
				ERROR(" (must be positive).\n");
				return 1;
			}
			break;

			case 4 :
			full_stddev = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid full standard deviation.\n");
				return 1;
			}
			if ( full_stddev < 0.0 ) {
				ERROR("Invalid full standard deviation.");
				ERROR(" (must be positive).\n");
				return 1;
			}
			break;

			case 5 :
			noise_stddev = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid noise standard deviation.\n");
				return 1;
			}
			if ( noise_stddev < 0.0 ) {
				ERROR("Invalid noise standard deviation.");
				ERROR(" (must be positive).\n");
				return 1;
			}
			break;

			case 6 :
			image_prefix = strdup(optarg);
			break;

			case 7 :
			background = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid background level.\n");
				return 1;
			}
			if ( background < 0.0 ) {
				ERROR("Background level must be positive.\n");
				return 1;
			}
			break;

			case 8 :
			divergence = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid beam divergence.\n");
				return 1;
			}
			if ( divergence < 0.0 ) {
				ERROR("Beam divergence must be positive.\n");
				return 1;
			}
			break;

			case 9 :
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

			case 10 :
			profile_radius = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid profile radius.\n");
				return 1;
			}
			if ( divergence < 0.0 ) {
				ERROR("Profile radius must be positive.\n");
				return 1;
			}
			break;

			case 11 :
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

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}
	}

	if ( n_threads < 1 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( (n_threads > 1) && (image_prefix != NULL) ) {
		ERROR("Option \"--images\" is incompatible with \"-j\".\n");
		return 1;
	}

	/* Load cell */
	if ( cellfile == NULL ) {
		ERROR("You need to give a PDB file with the unit cell.\n");
		return 1;
	}
	cell = load_cell_from_file(cellfile);
	if ( cell == NULL ) {
		ERROR("Failed to get cell from '%s'\n", cellfile);
		return 1;
	}
	free(cellfile);

	if ( !cell_is_sensible(cell) ) {
		ERROR("Invalid unit cell parameters:\n");
		cell_print(cell);
		return 1;
	}

	/* Load geometry */
	if ( geomfile == NULL ) {
		ERROR("You need to give a geometry file.\n");
		return 1;
	}
	det = get_detector_geometry(geomfile, &beam);
	if ( det == NULL ) {
		ERROR("Failed to read geometry from '%s'\n", geomfile);
		return 1;
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	free(sym_str);

	if ( save_file == NULL ) save_file = strdup("partial_sim.hkl");

	/* Load (full) reflections */
	if ( input_file != NULL ) {

		RefList *as;

		full = read_reflections(input_file);
		if ( full == NULL ) {
			ERROR("Failed to read reflections from '%s'\n",
			      input_file);
			return 1;
		}
		free(input_file);
		if ( check_list_symmetry(full, sym) ) {
			ERROR("The input reflection list does not appear to"
			      " have symmetry %s\n", symmetry_name(sym));
			return 1;
		}

		as = asymmetric_indices(full, sym);
		reflist_free(full);
		full = as;

	} else {
		random_intensities = 1;
	}

	if ( n < 1 ) {
		ERROR("Number of patterns must be at least 1.\n");
		return 1;
	}

	if ( output_file == NULL ) {
		ERROR("You must give a filename for the output.\n");
		return 1;
	}
	stream = open_stream_for_write_2(output_file, geomfile, argc, argv);
	if ( stream == NULL ) {
		ERROR("Couldn't open output file '%s'\n", output_file);
		return 1;
	}
	free(output_file);

	image.det = det;
	image.width = det->max_fs + 1;
	image.height = det->max_ss + 1;

	image.lambda = ph_en_to_lambda(eV_to_J(photon_energy));
	image.div = divergence;
	image.bw = bandwidth;
	image.filename = "dummy.h5";
	image.copyme = NULL;
	image.crystals = NULL;
	image.n_crystals = 0;
	image.indexed_by = INDEXING_SIMULATION;
	image.num_peaks = 0;
	image.num_saturated_peaks = 0;
	image.spectrum_size = 0;
	image.event = NULL;

	if ( random_intensities ) {
		full = reflist_new();
	}

	qargs.full = full;
	pthread_rwlock_init(&qargs.full_lock, NULL);
	qargs.n_to_do = n;
	qargs.n_done = 0;
	qargs.n_started = 0;
	qargs.sym = sym;
	qargs.random_intensities = random_intensities;
	qargs.cell = cell;
	qargs.template_image = &image;
	qargs.stream = stream;
	qargs.cnoise = cnoise;
	qargs.osf_stddev = osf_stddev;
	qargs.full_stddev = full_stddev;
	qargs.noise_stddev = noise_stddev;
	qargs.background = background;
	qargs.max_q = largest_q(&image);
	qargs.image_prefix = image_prefix;
	qargs.profile_radius = profile_radius;

	qargs.rngs = malloc(n_threads * sizeof(gsl_rng *));
	if ( qargs.rngs == NULL ) {
		ERROR("Failed to allocate RNGs\n");
		return 1;
	}

	if ( config_random ) {

		FILE *fh;

		fh = fopen("/dev/urandom", "r");
		if ( fh == NULL ) {
			ERROR("Failed to open /dev/urandom.  Try again without"
			      " --really-random.\n");
			return 1;
		}

		for ( i=0; i<n_threads; i++ ) {

			unsigned long int seed;

			fread(&seed, sizeof(seed), 1, fh);
			qargs.rngs[i] = gsl_rng_alloc(gsl_rng_mt19937);
			gsl_rng_set(qargs.rngs[i], seed);

		}

		fclose(fh);

	} else {
		rng_for_seeds = gsl_rng_alloc(gsl_rng_mt19937);
		for ( i=0; i<n_threads; i++ ) {
			qargs.rngs[i] = gsl_rng_alloc(gsl_rng_mt19937);
			gsl_rng_set(qargs.rngs[i], gsl_rng_get(rng_for_seeds));
		}
		gsl_rng_free(rng_for_seeds);
	}

	for ( i=0; i<NBINS; i++ ) {
		qargs.n_ref[i] = 0;
		qargs.p_hist[i] = 0.0;
		qargs.p_max[i] = 0.0;
	}

	run_threads(n_threads, run_job, create_job, finalise_job,
	            &qargs, n, 0, 0, 0);

	if ( random_intensities ) {
		STATUS("Writing full intensities to %s\n", save_file);
		write_reflist(save_file, full);
	}

	if ( phist_file != NULL ) {

		double overall_max = 0.0;
		double overall_mean = 0.0;
		long long int overall_total = 0;

		fh = fopen(phist_file, "w");

		if ( fh != NULL ) {

			for ( i=0; i<NBINS; i++ ) {

				double rcen;

				if ( qargs.p_max[i] > overall_max ) {
					overall_max = qargs.p_max[i];
				}

				overall_mean += qargs.p_hist[i];
				overall_total += qargs.n_ref[i];

				rcen = i/(double)NBINS*qargs.max_q
					  + qargs.max_q/(2.0*NBINS);
				fprintf(fh, "%.2f %7li %.3f %.3f\n", rcen/1.0e9,
					qargs.n_ref[i],
					qargs.p_hist[i]/qargs.n_ref[i],
					qargs.p_max[i]);

			}

			fclose(fh);

			overall_mean /= overall_total;

			STATUS("Overall max partiality = %.2f\n", overall_max);
			STATUS("Overall mean partiality = %.2f\n", overall_mean);
			STATUS("Total number of reflections = %lli\n",
			       overall_total);

		} else {
			ERROR("Failed to open file '%s' for writing.\n",
			      phist_file);
		}

	}

	for ( i=0; i<n_threads; i++ ) {
		gsl_rng_free(qargs.rngs[i]);
	}
	free(qargs.rngs);
	pthread_rwlock_destroy(&qargs.full_lock);
	close_stream(stream);
	cell_free(cell);
	free_detector_geometry(det);
	free_symoplist(sym);
	reflist_free(full);
	free(save_file);
	free(geomfile);

	return 0;
}
