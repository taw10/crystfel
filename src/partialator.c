/*
 * partialator.c
 *
 * Scaling and post refinement for coherent nanocrystallography
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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
#include <gsl/gsl_errno.h>

#include <utils.h>
#include <hdf5-file.h>
#include <symmetry.h>
#include <stream.h>
#include <geometry.h>
#include <peaks.h>
#include <thread-pool.h>
#include <beam-parameters.h>
#include <reflist.h>
#include <reflist-utils.h>

#include "post-refinement.h"
#include "hrs-scaling.h"
#include "scaling-report.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Scaling and post refinement for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"                              (must be a file, not e.g. stdin)\n"
"  -o, --output=<filename>    Output filename.  Default: partialator.hkl.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -b, --beam=<file>          Get beam parameters from file, which provides\n"
"                              initial values for parameters, and nominal\n"
"                              wavelengths if no per-shot value is found in \n"
"                              an HDF5 file.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"  -n, --iterations=<n>       Run <n> cycles of scaling and post-refinement.\n"
"      --no-scale             Fix all the scaling factors at unity.\n"
"  -r, --reference=<file>     Refine images against reflections in <file>,\n"
"                              instead of taking the mean of the intensity\n"
"                              estimates.\n"
"  -m, --model=<model>        Specify partiality model.\n"
"      --min-measurements=<n> Require at least <n> measurements before a\n"
"                             reflection appears in the output.  Default: 2\n"
"\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


struct refine_args
{
	RefList *full;
	Crystal *crystal;
	PartialityModel pmodel;
};


struct queue_args
{
	int n_started;
	int n_done;
	Crystal **crystals;
	int n_crystals;
	struct refine_args task_defaults;
};


static void refine_image(void *task, int id)
{
	struct refine_args *pargs = task;
	Crystal *cr = pargs->crystal;

	pr_refine(cr, pargs->full, pargs->pmodel);
}


static void *get_image(void *vqargs)
{
	struct refine_args *task;
	struct queue_args *qargs = vqargs;

	task = malloc(sizeof(struct refine_args));
	memcpy(task, &qargs->task_defaults, sizeof(struct refine_args));

	task->crystal = qargs->crystals[qargs->n_started];

	qargs->n_started++;

	return task;
}


static void done_image(void *vqargs, void *task)
{
	struct queue_args *qargs = vqargs;

	qargs->n_done++;

	progress_bar(qargs->n_done, qargs->n_crystals, "Refining");
	free(task);
}


static void refine_all(Crystal **crystals, int n_crystals,
                       struct detector *det,
                       RefList *full, int nthreads, PartialityModel pmodel)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	/* If the partiality model is "p=1", this refinement is really, really
	 * easy... */
	if ( pmodel == PMODEL_UNITY ) return;

	task_defaults.full = full;
	task_defaults.crystal = NULL;
	task_defaults.pmodel = pmodel;

	qargs.task_defaults = task_defaults;
	qargs.n_started = 0;
	qargs.n_done = 0;
	qargs.n_crystals = n_crystals;
	qargs.crystals = crystals;

	/* Don't have threads which are doing nothing */
	if ( n_crystals < nthreads ) nthreads = n_crystals;

	run_threads(nthreads, refine_image, get_image, done_image,
	            &qargs, n_crystals, 0, 0, 0);
}


/* Decide which reflections can be scaled */
static int select_scalable_reflections(RefList *list, RefList *reference)
{
	Reflection *refl;
	RefListIterator *iter;
	int n_acc = 0;
	int n_red = 0;
	int n_par = 0;
	int n_ref = 0;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		int sc = 1;

		/* This means the reflection was not found on the last check */
		if ( get_redundancy(refl) == 0 ) {
			sc = 0;
			n_red++;
		}

		/* Don't try to scale up reflections which are hardly there */
		if ( get_partiality(refl) < 0.1 ) {
			sc = 0;
			n_par++;
		}

		/* If we are scaling against a reference set, we additionally
		 * require that this reflection is in the reference list. */
		if ( reference != NULL ) {
			signed int h, k, l;
			get_indices(refl, &h, &k, &l);
			if ( find_refl(reference, h, k, l) == NULL ) {
				sc = 0;
				n_ref++;
			}
		}

		set_scalable(refl, sc);

		if ( sc ) n_acc++;
	}

	//STATUS("List %p: %i accepted, %i red zero, %i small part, %i no ref\n",
	//       list, n_acc, n_red, n_par, n_ref);

	return n_acc;
}


static void select_reflections_for_refinement(Crystal **crystals, int n,
                                              RefList *full, int have_reference)
{
	int i;

	for ( i=0; i<n; i++ ) {

		RefList *reflist;
		Reflection *refl;
		RefListIterator *iter;
		int n_acc = 0;
		int n_noscale = 0;
		int n_fewmatch = 0;
		int n_ref = 0;

		reflist = crystal_get_reflections(crystals[i]);
		for ( refl = first_refl(reflist, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			int sc;

			n_ref++;

			/* We require that the reflection itself is scalable
			 * (i.e. sensible partiality and intensity) and that
			 * the "full" estimate of this reflection is made from
			 * at least two parts. */
			get_indices(refl, &h, &k, &l);
			sc = get_scalable(refl);
			if ( !sc ) {

				n_noscale++;
				set_refinable(refl, 0);

			} else {

				Reflection *f = find_refl(full, h, k, l);

				if ( f != NULL ) {

					int r = get_redundancy(f);
					if ( (r >= 2) || have_reference ) {
						set_refinable(refl, 1);
						n_acc++;
					} else {
						n_fewmatch++;
					}

				} else {
					ERROR("%3i %3i %3i is scalable, but is"
					      " not in the reference list.\n",
					      h, k, l);
					abort();
				}

			}
		}

		//STATUS("Image %4i: %i guide reflections accepted "
		//       "(%i not scalable, %i few matches, %i total)\n",
		//       i, n_acc, n_noscale, n_fewmatch, n_ref);

	}
}


static void display_progress(int n_images, int n_crystals)
{
	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	pthread_mutex_lock(&stderr_lock);
	fprintf(stderr, "\r%i images loaded, %i crystals.",
	        n_images, n_crystals);
	pthread_mutex_unlock(&stderr_lock);

	fflush(stdout);
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *geomfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	int nthreads = 1;
	struct detector *det;
	int i;
	struct image *images;
	int n_iter = 10;
	struct beam_params *beam = NULL;
	RefList *full;
	int n_images = 0;
	int n_crystals = 0;
	int nobs;
	char *reference_file = NULL;
	RefList *reference = NULL;
	int have_reference = 0;
	char cmdline[1024];
	SRContext *sr;
	int noscale = 0;
	Stream *st;
	Crystal **crystals;
	char *pmodel_str = NULL;
	PartialityModel pmodel = PMODEL_SPHERE;
	int min_measurements = 2;
	char *rval;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},
		{"no-scale",           0, &noscale,            1},
		{"reference",          1, NULL,               'r'},
		{"model",              1, NULL,               'm'},
		{"min-measurements",   1, NULL,                2},

		{0, 0, NULL, 0}
	};

	cmdline[0] = '\0';
	for ( i=1; i<argc; i++ ) {
		strncat(cmdline, argv[i], 1023-strlen(cmdline));
		strncat(cmdline, " ", 1023-strlen(cmdline));
	}

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:g:b:y:n:r:j:m:",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			infile = strdup(optarg);
			break;

			case 'g' :
			geomfile = strdup(optarg);
			break;

			case 'j' :
			nthreads = atoi(optarg);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'n' :
			n_iter = atoi(optarg);
			break;

			case 'm' :
			pmodel_str = strdup(optarg);
			break;

			case 'b' :
			beam = get_beam_parameters(optarg);
			if ( beam == NULL ) {
				ERROR("Failed to load beam parameters"
				      " from '%s'\n", optarg);
				return 1;
			}
			break;

			case 'r' :
			reference_file = strdup(optarg);
			break;

			case 2 :
			errno = 0;
			min_measurements = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
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

	if ( nthreads < 1 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( infile == NULL ) {
		infile = strdup("-");
	}
	st = open_stream_for_read(infile);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", infile);
		return 1;
	}
	/* Don't free "infile", because it's needed for the scaling report */

	/* Sanitise output filename */
	if ( outfile == NULL ) {
		outfile = strdup("partialator.hkl");
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	free(sym_str);

	/* Get detector geometry */
	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	if ( beam == NULL ) {
		ERROR("You must provide a beam parameters file.\n");
		return 1;
	}

	if ( pmodel_str != NULL ) {
		if ( strcmp(pmodel_str, "sphere") == 0 ) {
			pmodel = PMODEL_SPHERE;
		} else if ( strcmp(pmodel_str, "unity") == 0 ) {
			pmodel = PMODEL_UNITY;
		} else {
			ERROR("Unknown partiality model '%s'.\n", pmodel_str);
			return 1;
		}
	}

	if ( reference_file != NULL ) {

		RefList *list;

		list = read_reflections(reference_file);
		if ( list == NULL ) {
			ERROR("Failed to read '%s'\n", reference_file);
			return 1;
		}
		free(reference_file);
		reference = asymmetric_indices(list, sym);
		reflist_free(list);
		have_reference = 1;

	}

	gsl_set_error_handler_off();

	/* Fill in what we know about the images so far */
	n_images = 0;
	n_crystals = 0;
	images = NULL;
	crystals = NULL;

	do {

		RefList *as;
		int i;
		struct image *images_new;
		struct image *cur;

		images_new = realloc(images, (n_images+1)*sizeof(struct image));
		if ( images_new == NULL ) {
			ERROR("Failed to allocate memory for image list.\n");
			return 1;
		}
		images = images_new;
		cur = &images[n_images];

		cur->det = det;
		if ( read_chunk(st, cur) != 0 ) {
			break;
		}

		/* Won't be needing this, if it exists */
		image_feature_list_free(cur->features);
		cur->features = NULL;
		cur->div = beam->divergence;
		cur->bw = beam->bandwidth;
		cur->width = det->max_fs;
		cur->height = det->max_ss;
		cur->data = NULL;
		cur->flags = NULL;
		cur->beam = NULL;

		n_images++;

		for ( i=0; i<cur->n_crystals; i++ ) {

			Crystal *cr;
			Crystal **crystals_new;
			RefList *cr_refl;

			crystals_new = realloc(crystals,
			                      (n_crystals+1)*sizeof(Crystal *));
			if ( crystals_new == NULL ) {
				ERROR("Failed to allocate memory for crystal "
				      "list.\n");
				return 1;
			}
			crystals = crystals_new;
			crystals[n_crystals] = cur->crystals[i];
			cr = crystals[n_crystals];

			/* Image pointer will change due to later reallocs */
			crystal_set_image(cr, NULL);

			/* Fill in initial estimates of stuff */
			crystal_set_osf(cr, 1.0);
			crystal_set_profile_radius(cr, beam->profile_radius);
			crystal_set_mosaicity(cr, 0.0);
			crystal_set_user_flag(cr, 0);

			/* This is the raw list of reflections */
			cr_refl = crystal_get_reflections(cr);
			as = asymmetric_indices(cr_refl, sym);
			crystal_set_reflections(cr, as);
			reflist_free(cr_refl);

			n_crystals++;

		}

		display_progress(n_images, n_crystals);

	} while ( 1 );
	fprintf(stderr, "\n");

	close_stream(st);

	/* Fill in image pointers */
	nobs = 0;
	for ( i=0; i<n_images; i++ ) {
		int j;
		for ( j=0; j<images[i].n_crystals; j++ ) {

			Crystal *cryst;
			RefList *as;
			int n_gained = 0;
			int n_lost = 0;

			cryst = images[i].crystals[j];
			crystal_set_image(cryst, &images[i]);

			/* Now it's safe to do the following */
			update_partialities_2(cryst, pmodel,
			                      &n_gained, &n_lost);
			assert(n_gained == 0);  /* That'd just be silly */
			as = crystal_get_reflections(cryst);
			nobs += select_scalable_reflections(as, reference);

		}
	}
	STATUS("%i scalable observations.\n", nobs);

	/* Make initial estimates */
	STATUS("Performing initial scaling.\n");
	if ( noscale ) STATUS("Scale factors fixed at 1.\n");
	full = scale_intensities(crystals, n_crystals, reference,
	                         nthreads, noscale, pmodel, min_measurements);

	sr = sr_titlepage(crystals, n_crystals, "scaling-report.pdf",
	                  infile, cmdline);
	sr_iteration(sr, 0, crystals, n_crystals, full);

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		int n_dud = 0;
		int j;
		RefList *comp;

		STATUS("Post refinement cycle %i of %i\n", i+1, n_iter);

		/* Refine the geometry of all patterns to get the best fit */
		comp = (reference == NULL) ? full : reference;
		select_reflections_for_refinement(crystals, n_crystals,
		                                  comp, have_reference);
		refine_all(crystals, n_crystals, det, comp, nthreads, pmodel);

		nobs = 0;
		for ( j=0; j<n_crystals; j++ ) {
			Crystal *cr = crystals[j];
			RefList *rf = crystal_get_reflections(cr);
			if ( crystal_get_user_flag(cr) ) n_dud++;
			nobs += select_scalable_reflections(rf, reference);
		}

		STATUS("%i crystals could not be refined this cycle.\n", n_dud);

		/* Re-estimate all the full intensities */
		reflist_free(full);
		full = scale_intensities(crystals, n_crystals,
		                         reference, nthreads, noscale, pmodel,
		                         min_measurements);

		sr_iteration(sr, i+1, crystals, n_crystals, full);

	}

	sr_finish(sr);

	/* Output results */
	write_reflist(outfile, full);

	/* Clean up */
	for ( i=0; i<n_crystals; i++ ) {
		reflist_free(crystal_get_reflections(crystals[i]));
		crystal_free(crystals[i]);
	}
	reflist_free(full);
	free(sym);
	free(outfile);
	free_detector_geometry(det);
	free(beam);
	free(crystals);
	if ( reference != NULL ) {
		reflist_free(reference);
	}
	for ( i=0; i<n_images; i++ ) {
		free(images[i].filename);
	}
	free(images);
	free(infile);

	return 0;
}
