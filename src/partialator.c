/*
 * partialator.c
 *
 * Scaling and post refinement for coherent nanocrystallography
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
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

#include <image.h>
#include <utils.h>
#include <symmetry.h>
#include <stream.h>
#include <geometry.h>
#include <peaks.h>
#include <thread-pool.h>
#include <reflist.h>
#include <reflist-utils.h>

#include "version.h"
#include "post-refinement.h"
#include "hrs-scaling.h"
#include "scaling-report.h"
#include "rejection.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Scaling and post refinement for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"  -o, --output=<filename>    Output filename.  Default: partialator.hkl.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"  -n, --iterations=<n>       Run <n> cycles of scaling and post-refinement.\n"
"      --no-scale             Fix all the scaling factors at unity.\n"
"  -m, --model=<model>        Specify partiality model.\n"
"      --min-measurements=<n> Minimum number of measurements to require.\n"
"      --no-polarisation      Disable polarisation correction.\n"
"      --max-adu=<n>          Saturation value of detector.\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


struct refine_args
{
	RefList *full;
	Crystal *crystal;
	PartialityModel pmodel;
	struct prdata prdata;
};


struct queue_args
{
	int n_started;
	int n_done;
	Crystal **crystals;
	int n_crystals;
	struct srdata *srdata;
	struct refine_args task_defaults;
};


static void refine_image(void *task, int id)
{
	struct refine_args *pargs = task;
	Crystal *cr = pargs->crystal;

	pargs->prdata = pr_refine(cr, pargs->full, pargs->pmodel);
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
	struct refine_args *pargs = task;

	qargs->n_done++;
	if ( pargs->prdata.refined ) {
		qargs->srdata->n_refined += pargs->prdata.refined;
		qargs->srdata->n_filtered += pargs->prdata.n_filtered;
	}

	progress_bar(qargs->n_done, qargs->n_crystals, "Refining");
	free(task);
}


static void refine_all(Crystal **crystals, int n_crystals,
                       RefList *full, int nthreads, PartialityModel pmodel,
                       struct srdata *srdata)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	/* If the partiality model is "p=1", this refinement is really, really
	 * easy... */
	if ( pmodel == PMODEL_UNITY ) return;

	task_defaults.full = full;
	task_defaults.crystal = NULL;
	task_defaults.pmodel = pmodel;
	task_defaults.prdata.refined = 0;
	task_defaults.prdata.n_filtered = 0;

	qargs.task_defaults = task_defaults;
	qargs.n_started = 0;
	qargs.n_done = 0;
	qargs.n_crystals = n_crystals;
	qargs.crystals = crystals;
	qargs.srdata = srdata;

	/* Don't have threads which are doing nothing */
	if ( n_crystals < nthreads ) nthreads = n_crystals;

	run_threads(nthreads, refine_image, get_image, done_image,
	            &qargs, n_crystals, 0, 0, 0);

	STATUS("%5.2f eigenvalues filtered on final iteration per successfully "
	       "refined crystal\n",
	       (double)srdata->n_filtered/srdata->n_refined);
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


static const char *str_flags(Crystal *cr)
{
	switch ( crystal_get_user_flag(cr) ) {

		case 0 :
		return "OK";

		case 1 :
		return "bad scaling";

		case 2 :
		return "not enough reflections";

		case 3 :
		return "PR solve failed";

		case 4 :
		return "PR lost too many reflections";

		case 5 :
		return "Early rejection";

		default :
		return "Unknown flag";
	}
}


static RefList *apply_max_adu(RefList *list, double max_adu)
{
	RefList *nlist;
	Reflection *refl;
	RefListIterator *iter;

	nlist = reflist_new();
	if ( nlist == NULL ) return NULL;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		if ( get_peak(refl) < max_adu ) {
			signed int h, k, l;
			get_indices(refl, &h, &k, &l);
			Reflection *nrefl = add_refl(nlist, h, k, l);
			copy_data(nrefl, refl);
		}
	}
	reflist_free(list);
	return nlist;
}


static void skip_to_end(FILE *fh)
{
	int c;
	do {
		c = fgetc(fh);
	} while ( (c != '\n') && (c != EOF) );
}


static int set_initial_params(Crystal *cr, FILE *fh)
{
	if ( fh != NULL ) {

		int err;
		int n;
		float osf, B;

		err = fscanf(fh, "%i %f %f", &n, &osf, &B);
		if ( err != 3 ) {
			ERROR("Failed to read parameters.\n");
			return 1;
		}

		crystal_set_osf(cr, osf);
		crystal_set_Bfac(cr, B*1e-20);

		skip_to_end(fh);

	} else {

		crystal_set_osf(cr, 1.0);
		crystal_set_Bfac(cr, 0.0);

	}

	return 0;
}


static void show_duds(Crystal **crystals, int n_crystals)
{
	int j;
	int n_dud = 0;
	int n_noscale = 0;
	int n_noref = 0;
	int n_solve = 0;
	int n_lost = 0;
	int n_early = 0;

	for ( j=0; j<n_crystals; j++ ) {
		int flag;
		flag = crystal_get_user_flag(crystals[j]);
		if ( flag != 0 ) n_dud++;
		switch ( flag ) {

			case 0:
			break;

			case 1:
			n_noscale++;
			break;

			case 2:
			n_noref++;
			break;

			case 3:
			n_solve++;
			break;

			case 4:
			n_lost++;
			break;

			case 5:
			n_early++;
			break;

			default:
			STATUS("Unknown flag %i\n", flag);
			break;
		}
	}

	if ( n_dud ) {
		STATUS("%i bad crystals:\n", n_dud);
		STATUS(" %i scaling failed.\n", n_noscale);
		STATUS(" %i not enough reflections.\n", n_noref);
		STATUS(" %i solve failed.\n", n_solve);
		STATUS(" %i lost too many reflections.\n", n_lost);
		STATUS(" %i early rejection.\n", n_early);
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	int nthreads = 1;
	int i;
	struct image *images;
	int n_iter = 10;
	RefList *full;
	int n_images = 0;
	int n_crystals = 0;
	char cmdline[1024];
	SRContext *sr;
	int noscale = 0;
	Stream *st;
	Crystal **crystals;
	char *pmodel_str = NULL;
	PartialityModel pmodel = PMODEL_SCSPHERE;
	int min_measurements = 2;
	char *rval;
	struct srdata srdata;
	int polarisation = 1;
	double max_adu = +INFINITY;
	char *sparams_fn = NULL;
	FILE *sparams_fh;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               'v'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},
		{"reference",          1, NULL,               'r'},
		{"model",              1, NULL,               'm'},

		{"min-measurements",   1, NULL,                2},
		{"max-adu",            1, NULL,                3},
		{"start-params",       1, NULL,                4},

		{"no-scale",           0, &noscale,            1},
		{"no-polarisation",    0, &polarisation,       0},
		{"no-polarization",    0, &polarisation,       0},
		{"polarisation",       0, &polarisation,       1}, /* compat */
		{"polarization",       0, &polarisation,       1}, /* compat */

		{0, 0, NULL, 0}
	};

	cmdline[0] = '\0';
	for ( i=1; i<argc; i++ ) {
		strncat(cmdline, argv[i], 1023-strlen(cmdline));
		strncat(cmdline, " ", 1023-strlen(cmdline));
	}

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:g:b:y:n:j:m:",
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
			infile = strdup(optarg);
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

			case 2 :
			errno = 0;
			min_measurements = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
				return 1;
			}
			break;

			case 3 :
			errno = 0;
			max_adu = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-adu.\n");
				return 1;
			}
			break;

			case 4 :
			sparams_fn = strdup(optarg);
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

	if ( pmodel_str != NULL ) {
		if ( strcmp(pmodel_str, "unity") == 0 ) {
			pmodel = PMODEL_UNITY;
		} else if ( strcmp(pmodel_str, "scgaussian") == 0 ) {
			pmodel = PMODEL_SCGAUSSIAN;
		} else if ( strcmp(pmodel_str, "scsphere") == 0 ) {
			pmodel = PMODEL_SCSPHERE;
		} else {
			ERROR("Unknown partiality model '%s'.\n", pmodel_str);
			return 1;
		}
	}

	gsl_set_error_handler_off();

	/* Fill in what we know about the images so far */
	n_images = 0;
	n_crystals = 0;
	images = NULL;
	crystals = NULL;
	if ( sparams_fn != NULL ) {
		char line[1024];
		sparams_fh = fopen(sparams_fn, "r");
		if ( sparams_fh == NULL ) {
			ERROR("Failed to open '%s'\n", sparams_fn);
			return 1;
		}
		fgets(line, 1024, sparams_fh);
		STATUS("Reading initial scaling factors (G,B) from '%s'\n",
		       sparams_fn);
		free(sparams_fn);
	} else {
		sparams_fh = NULL;
	}

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

		cur->div = NAN;
		cur->bw = NAN;
		cur->det = NULL;
		if ( read_chunk_2(st, cur, STREAM_READ_REFLECTIONS
		                           | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		if ( isnan(cur->div) || isnan(cur->bw) ) {
			ERROR("Chunk doesn't contain beam parameters.\n");
			return 1;
		}

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

			/* This is the raw list of reflections */
			cr_refl = crystal_get_reflections(cr);

			cr_refl = apply_max_adu(cr_refl, max_adu);

			if ( polarisation ) {
				polarisation_correction(cr_refl,
						        crystal_get_cell(cr),
						        cur);
			}

			as = asymmetric_indices(cr_refl, sym);
			crystal_set_reflections(cr, as);
			crystal_set_user_flag(cr, 0);
			reflist_free(cr_refl);

			if ( set_initial_params(cr, sparams_fh) ) {
				ERROR("Failed to set initial parameters\n");
				return 1;
			}

			n_crystals++;

		}

		if ( n_images % 100 == 0 ) {
			display_progress(n_images, n_crystals);
		}

	} while ( 1 );
	display_progress(n_images, n_crystals);
	fprintf(stderr, "\n");
	if ( sparams_fh != NULL ) fclose(sparams_fh);

	close_stream(st);

	/* Fill in image pointers */
	for ( i=0; i<n_images; i++ ) {
		int j;
		for ( j=0; j<images[i].n_crystals; j++ ) {

			Crystal *cryst;

			cryst = images[i].crystals[j];
			crystal_set_image(cryst, &images[i]);

			/* Now it's safe to do the following */
			update_partialities(cryst, pmodel);

		}
	}

	/* Make a first pass at cutting out crap */
//	STATUS("Checking patterns.\n");
//	early_rejection(crystals, n_crystals);

	/* Make initial estimates */
	STATUS("Performing initial scaling.\n");
	if ( noscale ) {
		STATUS("Skipping scaling step (--no-scale).\n");
		full = lsq_intensities(crystals, n_crystals, nthreads, pmodel,
		                       min_measurements);
	} else {
		full = scale_intensities(crystals, n_crystals, nthreads, pmodel,
		                         min_measurements);
	}

	check_rejection(crystals, n_crystals);

	srdata.crystals = crystals;
	srdata.n = n_crystals;
	srdata.full = full;
	srdata.n_filtered = 0;
	srdata.n_refined = 0;

	sr = sr_titlepage(crystals, n_crystals, "scaling-report.pdf",
	                  infile, cmdline);
	sr_iteration(sr, 0, &srdata);

	show_duds(crystals, n_crystals);

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		STATUS("Post refinement cycle %i of %i\n", i+1, n_iter);

		srdata.n_filtered = 0;

		/* Refine the geometry of all patterns to get the best fit */
		refine_all(crystals, n_crystals, full, nthreads, pmodel,
		           &srdata);

		show_duds(crystals, n_crystals);
		check_rejection(crystals, n_crystals);

		/* Re-estimate all the full intensities */
		reflist_free(full);
		if ( noscale ) {
			STATUS("Skipping scaling step (--no-scale).\n");
			full = lsq_intensities(crystals, n_crystals, nthreads,
			                       pmodel, min_measurements);
		} else {
			full = scale_intensities(crystals, n_crystals, nthreads,
			                         pmodel, min_measurements);
		}

		check_rejection(crystals, n_crystals);

		srdata.full = full;

		sr_iteration(sr, i+1, &srdata);

	}

	sr_finish(sr);

	/* Output results */
	write_reflist(outfile, full);

	/* Dump parameters */
	FILE *fh;
	fh = fopen("partialator.params", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open partialator.params!\n");
	} else {
		fprintf(fh, "  cr        OSF       relB         div flag\n");
		for ( i=0; i<n_crystals; i++ ) {
			fprintf(fh, "%4i %10.5f %10.2f %8.5e %s\n", i,
			        crystal_get_osf(crystals[i]),
				crystal_get_Bfac(crystals[i])*1e20,
			        crystal_get_image(crystals[i])->div,
			        str_flags(crystals[i]));
		}
		fclose(fh);
	}

	/* Clean up */
	for ( i=0; i<n_crystals; i++ ) {
		reflist_free(crystal_get_reflections(crystals[i]));
		crystal_free(crystals[i]);
	}
	reflist_free(full);
	free(sym);
	free(outfile);
	free(crystals);
	for ( i=0; i<n_images; i++ ) {
		free(images[i].filename);
	}
	free(images);
	free(infile);

	return 0;
}
