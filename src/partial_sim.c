/*
 * partial_sim.c
 *
 * Generate partials for testing scaling
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2012 Thomas White <taw@physics.org>
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

#include <utils.h>
#include <reflist-utils.h>
#include <symmetry.h>
#include <beam-parameters.h>
#include <detector.h>
#include <geometry.h>
#include <stream.h>
#include <thread-pool.h>

/* Number of bins for partiality graph */
#define NBINS 50


static void mess_up_cell(UnitCell *cell, double cnoise)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	//STATUS("Real:\n");
	//cell_print(cell);

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	ax = flat_noise(ax, cnoise*fabs(ax)/100.0);
	ay = flat_noise(ay, cnoise*fabs(ay)/100.0);
	az = flat_noise(az, cnoise*fabs(az)/100.0);
	bx = flat_noise(bx, cnoise*fabs(bx)/100.0);
	by = flat_noise(by, cnoise*fabs(by)/100.0);
	bz = flat_noise(bz, cnoise*fabs(bz)/100.0);
	cx = flat_noise(cx, cnoise*fabs(cx)/100.0);
	cy = flat_noise(cy, cnoise*fabs(cy)/100.0);
	cz = flat_noise(cz, cnoise*fabs(cz)/100.0);
	cell_set_reciprocal(cell, ax, ay, az, bx, by, bz, cx, cy, cz);

	//STATUS("Changed:\n");
	//cell_print(cell);
}


/* For each reflection in "partial", fill in what the intensity would be
 * according to "full" */
static void calculate_partials(RefList *partial, double osf,
                               RefList *full, const SymOpList *sym,
                               int random_intensities,
                               pthread_mutex_t *full_lock,
                               unsigned long int *n_ref, double *p_hist,
                               double *p_max, double max_q, UnitCell *cell)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(partial, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		Reflection *rfull;
		double p, Ip, If;
		int bin;

		get_indices(refl, &h, &k, &l);
		get_asymm(sym, h, k, l, &h, &k, &l);
		p = get_partiality(refl);

		pthread_mutex_lock(full_lock);
		rfull = find_refl(full, h, k, l);
		pthread_mutex_unlock(full_lock);

		if ( rfull == NULL ) {
			if ( random_intensities ) {

				/* The full reflection is immutable (in this
				 * program) once created, but creating it must
				 * be an atomic operation.  So do the whole
				 * thing under lock. */
				pthread_mutex_lock(full_lock);
				rfull = add_refl(full, h, k, l);
				If = fabs(gaussian_noise(0.0, 1000.0));
				set_intensity(rfull, If);
				set_redundancy(rfull, 1);
				pthread_mutex_unlock(full_lock);

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

		Ip = osf * p * If;

		bin = NBINS*2.0*resolution(cell, h, k, l) / max_q;
		if ( (bin < NBINS) && (bin>=0) ) {
			p_hist[bin] += p;
			n_ref[bin]++;
			if ( p > p_max[bin] ) p_max[bin] = p;
		} else {
			STATUS("Reflection out of histogram range: %e %i %f\n",
			       resolution(cell, h, k, l), bin,  p);
		}

		Ip = gaussian_noise(Ip, 100.0);

		set_intensity(refl, Ip);
		set_esd_intensity(refl, 100.0);
	}
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Generate a stream containing partials from a reflection list.\n"
"\n"
" -h, --help              Display this help message.\n"
"\n"
"You need to provide the following basic options:\n"
" -i, --input=<file>       Read reflections from <file>.\n"
"                           Default: generate random ones instead (see -r).\n"
" -o, --output=<file>      Write partials in stream format to <file>.\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -b, --beam=<file>        Get beam parameters from file\n"
" -p, --pdb=<file>         PDB file from which to get the unit cell.\n"
"\n"
" -y, --symmetry=<sym>     Symmetry of the input reflection list.\n"
" -n <n>                   Simulate <n> patterns.  Default: 2.\n"
" -r, --save-random=<file> Save randomly generated intensities to file.\n"
" -c, --cnoise=<val>       Add random noise, with a flat distribution, to the\n"
"                          reciprocal lattice vector components given in the\n"
"                          stream, with maximum error +/- <val> percent.\n"
"\n"
);
}


struct queue_args
{
	RefList *full;
	pthread_mutex_t full_lock;

	int n_done;
	int n_started;
	int n_to_do;

	SymOpList *sym;
	int random_intensities;
	UnitCell *cell;
	double cnoise;

	struct image *template_image;
	double max_q;

	/* The overall histogram */
	double p_hist[NBINS];
	unsigned long int n_ref[NBINS];
	double p_max[NBINS];

	FILE *stream;
};


struct worker_args
{
	struct queue_args *qargs;
	struct image image;

	/* Histogram for this image */
	double p_hist[NBINS];
	unsigned long int n_ref[NBINS];
	double p_max[NBINS];
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

	return wargs;
}


static void run_job(void *vwargs, int cookie)
{
	double osf;
	struct quaternion orientation;
	struct worker_args *wargs = vwargs;
	struct queue_args *qargs = wargs->qargs;
	int i;

	osf = gaussian_noise(1.0, 0.3);

	/* Set up a random orientation */
	orientation = random_quaternion();
	wargs->image.indexed_cell = cell_rotate(qargs->cell, orientation);

	snprintf(wargs->image.filename, 255, "dummy.h5");
	wargs->image.reflections = find_intersections(&wargs->image,
	                                       wargs->image.indexed_cell);

	for ( i=0; i<NBINS; i++ ) {
		wargs->n_ref[i] = 0;
		wargs->p_hist[i] = 0.0;
		wargs->p_max[i] = 0.0;
	}

	calculate_partials(wargs->image.reflections, osf, qargs->full,
	                   qargs->sym, qargs->random_intensities,
	                   &qargs->full_lock,
	                   wargs->n_ref, wargs->p_hist, wargs->p_max,
	                   qargs->max_q, wargs->image.indexed_cell);

	/* Give a slightly incorrect cell in the stream */
	mess_up_cell(wargs->image.indexed_cell, qargs->cnoise);
}


static void finalise_job(void *vqargs, void *vwargs)
{
	struct worker_args *wargs = vwargs;
	struct queue_args *qargs = vqargs;
	int i;

	write_chunk(qargs->stream, &wargs->image, NULL, STREAM_INTEGRATED);

	for ( i=0; i<NBINS; i++ ) {
		qargs->n_ref[i] += wargs->n_ref[i];
		qargs->p_hist[i] += wargs->p_hist[i];
		if ( wargs->p_max[i] > qargs->p_max[i] ) {
			qargs->p_max[i] = wargs->p_max[i];
		}
	}

	qargs->n_done++;
	progress_bar(qargs->n_done, qargs->n_to_do, "Simulating");

	reflist_free(wargs->image.reflections);
	cell_free(wargs->image.indexed_cell);
	free(wargs);
}


int main(int argc, char *argv[])
{
	int c;
	char *input_file = NULL;
	char *output_file = NULL;
	char *beamfile = NULL;
	char *geomfile = NULL;
	char *cellfile = NULL;
	struct detector *det = NULL;
	struct beam_params *beam = NULL;
	RefList *full = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	UnitCell *cell = NULL;
	FILE *ofh;
	int n = 2;
	int random_intensities = 0;
	char *save_file = NULL;
	struct queue_args qargs;
	struct image image;
	int n_threads = 1;
	double cnoise = 0.0;
	char *rval;
	int i;
	FILE *fh;
	char *phist_file = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"input",              1, NULL,               'i'},
		{"beam",               1, NULL,               'b'},
		{"pdb",                1, NULL,               'p'},
		{"geometry",           1, NULL,               'g'},
		{"symmetry",           1, NULL,               'y'},
		{"save-random",        1, NULL,               'r'},
		{"pgraph",             1, NULL,                2},
		{"cnoise",             1, NULL,               'c'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:b:p:g:y:n:r:j:c:",
	                        longopts, NULL)) != -1)
	{
		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'o' :
			output_file = strdup(optarg);
			break;

			case 'i' :
			input_file = strdup(optarg);
			break;

			case 'b' :
			beamfile = strdup(optarg);
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

			case 0 :
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

	/* Load beam */
	if ( beamfile == NULL ) {
		ERROR("You need to provide a beam parameters file.\n");
		return 1;
	}
	beam = get_beam_parameters(beamfile);
	if ( beam == NULL ) {
		ERROR("Failed to load beam parameters from '%s'\n", beamfile);
		return 1;
	}
	free(beamfile);

	/* Load cell */
	if ( cellfile == NULL ) {
		ERROR("You need to give a PDB file with the unit cell.\n");
		return 1;
	}
	cell = load_cell_from_pdb(cellfile);
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
	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	free(sym_str);

	if ( save_file == NULL ) save_file = strdup("partial_sim.hkl");

	/* Load (full) reflections */
	if ( input_file != NULL ) {

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
	ofh = fopen(output_file, "w");
	if ( ofh == NULL ) {
		ERROR("Couldn't open output file '%s'\n", output_file);
		return 1;
	}
	free(output_file);
	write_stream_header(ofh, argc, argv);

	image.det = det;
	image.width = det->max_fs;
	image.height = det->max_ss;

	image.lambda = ph_en_to_lambda(eV_to_J(beam->photon_energy));
	image.div = beam->divergence;
	image.bw = beam->bandwidth;
	image.profile_radius = beam->profile_radius;
	image.filename = malloc(256);
	image.copyme = NULL;

	if ( random_intensities ) {
		full = reflist_new();
	}

	qargs.full = full;
	pthread_mutex_init(&qargs.full_lock, NULL);
	qargs.n_to_do = n;
	qargs.n_done = 0;
	qargs.n_started = 0;
	qargs.sym = sym;
	qargs.random_intensities = random_intensities;
	qargs.cell = cell;
	qargs.template_image = &image;
	qargs.stream = ofh;
	qargs.cnoise = cnoise;
	qargs.max_q = largest_q(&image);

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

		fh = fopen(phist_file, "w");

		if ( fh != NULL ) {

			for ( i=0; i<NBINS; i++ ) {

				double rcen;

				rcen = i/(double)NBINS*qargs.max_q
					  + qargs.max_q/(2.0*NBINS);
				fprintf(fh, "%.2f %7li %.3f %.3f\n", rcen/1.0e9,
					qargs.n_ref[i],
					qargs.p_hist[i]/qargs.n_ref[i],
					qargs.p_max[i]);

			}

			fclose(fh);

		} else {
			ERROR("Failed to open file '%s' for writing.\n",
			      phist_file);
		}

	}

	fclose(ofh);
	cell_free(cell);
	free_detector_geometry(det);
	free(beam);
	free_symoplist(sym);
	reflist_free(full);
	free(image.filename);

	return 0;
}
