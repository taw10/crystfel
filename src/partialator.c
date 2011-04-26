/*
 * partialator.c
 *
 * Scaling and post refinement for coherent nanocrystallography
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include <assert.h>
#include <pthread.h>

#include "utils.h"
#include "hdf5-file.h"
#include "symmetry.h"
#include "stream.h"
#include "geometry.h"
#include "peaks.h"
#include "thread-pool.h"
#include "beam-parameters.h"
#include "post-refinement.h"
#include "hrs-scaling.h"
#include "reflist.h"
#include "reflist-utils.h"


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
"  -o, --output=<filename>    Output filename.  Default: facetron.hkl.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -b, --beam=<file>          Get beam parameters from file, which provides\n"
"                              initial values for parameters, and nominal\n"
"                              wavelengths if no per-shot value is found in \n"
"                              an HDF5 file.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"  -n, --iterations=<n>       Run <n> cycles of scaling and post-refinement.\n"
"\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


struct refine_args
{
	const char *sym;
	ReflItemList *obs;
	RefList *full;
	struct image *image;
	FILE *graph;
	FILE *pgraph;
};


struct queue_args
{
	int n;
	int n_done;
	int n_total_patterns;
	struct image *images;
	struct refine_args task_defaults;
};


static void refine_image(void *task, int id)
{
	struct refine_args *pargs = task;
	struct image *image = pargs->image;
	image->id = id;

	pr_refine(image, pargs->full, pargs->sym);
}


static void *get_image(void *vqargs)
{
	struct refine_args *task;
	struct queue_args *qargs = vqargs;

	task = malloc(sizeof(struct refine_args));
	memcpy(task, &qargs->task_defaults, sizeof(struct refine_args));

	task->image = &qargs->images[qargs->n];

	qargs->n++;

	return task;
}


static void done_image(void *vqargs, void *task)
{
	struct queue_args *qargs = vqargs;

	qargs->n_done++;

	progress_bar(qargs->n_done, qargs->n_total_patterns, "Refining");
	free(task);
}


static void refine_all(struct image *images, int n_total_patterns,
                       struct detector *det, const char *sym,
                       ReflItemList *obs, RefList *full, int nthreads,
                       FILE *graph, FILE *pgraph)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	task_defaults.sym = sym;
	task_defaults.obs = obs;
	task_defaults.full = full;
	task_defaults.image = NULL;
	task_defaults.graph = graph;
	task_defaults.pgraph = pgraph;

	qargs.task_defaults = task_defaults;
	qargs.n = 0;
	qargs.n_done = 0;
	qargs.n_total_patterns = n_total_patterns;
	qargs.images = images;

	run_threads(nthreads, refine_image, get_image, done_image,
	            &qargs, n_total_patterns, 0, 0, 0);
}


/* Decide which reflections can be scaled */
static void select_scalable_reflections(struct image *images, int n)
{
	int m;
	int n_scalable = 0;

	for ( m=0; m<n; m++ ) {

		Reflection *refl;
		RefListIterator *iter;

		for ( refl = first_refl(images[m].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			int scalable = 1;
			double v;

			if ( get_partiality(refl) < 0.1 ) scalable = 0;
			v = fabs(get_intensity(refl));
			if ( v < 0.1 ) scalable = 0;

			set_scalable(refl, scalable);
			if ( scalable ) n_scalable++;

		}

	}
	STATUS("%i reflections selected as scalable.\n", n_scalable);
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *geomfile = NULL;
	char *sym = NULL;
	FILE *fh;
	int nthreads = 1;
	struct detector *det;
	ReflItemList *obs;
	int i;
	int n_total_patterns;
	struct image *images;
	int n_iter = 10;
	struct beam_params *beam = NULL;
	RefList *full;
	int n_found = 0;
	int n_expected = 0;
	int n_notfound = 0;
	char *cref;
	int n_usable_patterns = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:x:j:y:o:b:",
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
			sym = strdup(optarg);
			break;

		case 'o' :
			outfile = strdup(optarg);
			break;

		case 'n' :
			n_iter = atoi(optarg);
			break;

		case 'b' :
			beam = get_beam_parameters(optarg);
			if ( beam == NULL ) {
				ERROR("Failed to load beam parameters"
				      " from '%s'\n", optarg);
				return 1;
			}
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	/* Sanitise input filename and open */
	if ( infile == NULL ) {
		infile = strdup("-");
	}
	if ( strcmp(infile, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(infile, "r");
	}
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", infile);
		return 1;
	}
	free(infile);

	/* Sanitise output filename */
	if ( outfile == NULL ) {
		outfile = strdup("partialator.hkl");
	}

	if ( sym == NULL ) sym = strdup("1");

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

	n_total_patterns = count_patterns(fh);
	if ( n_total_patterns == 0 ) {
		ERROR("No patterns to process.\n");
		return 1;
	}
	STATUS("There are %i patterns to process\n", n_total_patterns);

	images = malloc(n_total_patterns * sizeof(struct image));
	if ( images == NULL ) {
		ERROR("Couldn't allocate memory for images.\n");
		return 1;
	}

	/* Fill in what we know about the images so far */
	rewind(fh);
	obs = new_items();
	for ( i=0; i<n_total_patterns; i++ ) {

		images[n_usable_patterns].det = NULL;

		if ( read_chunk(fh, &images[n_usable_patterns]) != 0 ) {
			/* Should not happen, because we counted the patterns
			 * earlier. */
			ERROR("Failed to read chunk from the input stream.\n");
			return 1;
		}

		/* Won't be needing this, if it exists */
		image_feature_list_free(images[n_usable_patterns].features);
		images[n_usable_patterns].features = NULL;

		/* "n_usable_patterns" will not be incremented in this case */
		if ( images[n_usable_patterns].indexed_cell == NULL ) continue;

		/* Fill in initial estimates of stuff */
		images[n_usable_patterns].div = beam->divergence;
		images[n_usable_patterns].bw = beam->bandwidth;
		images[n_usable_patterns].det = det;
		images[n_usable_patterns].width = det->max_fs;
		images[n_usable_patterns].height = det->max_ss;
		images[n_usable_patterns].osf = 1.0;
		images[n_usable_patterns].profile_radius = 0.005e9;

		/* Muppet proofing */
		images[n_usable_patterns].data = NULL;
		images[n_usable_patterns].flags = NULL;
		images[n_usable_patterns].beam = NULL;

		/* This is the raw list of reflections */
		images[n_usable_patterns].raw_reflections =
		                          images[n_usable_patterns].reflections;
		images[n_usable_patterns].reflections = NULL;

		update_partialities_and_asymm(&images[n_usable_patterns], sym,
		                              obs, &n_expected, &n_found,
		                              &n_notfound);

		progress_bar(i, n_total_patterns-1, "Loading pattern data");
		n_usable_patterns++;

	}
	fclose(fh);
	STATUS("Found %5.2f%% of the expected peaks (missed %i of %i).\n",
	       100.0 * (double)n_found / n_expected, n_notfound, n_expected);
	STATUS("Mean measurements per unique reflection: %5.2f\n",
	       (double)n_found / num_items(obs));

	cref = find_common_reflections(images, n_usable_patterns);

	/* Make initial estimates */
	STATUS("Performing initial scaling.\n");
	select_scalable_reflections(images, n_total_patterns);
	full = scale_intensities(images, n_usable_patterns, sym, obs, cref);

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		FILE *fhg;
		FILE *fhp;
		char filename[1024];

		STATUS("Post refinement iteration %i of %i\n", i+1, n_iter);

		snprintf(filename, 1023, "p-iteration-%i.dat", i+1);
		fhg = fopen(filename, "w");
		if ( fhg == NULL ) {
			ERROR("Failed to open '%s'\n", filename);
			/* Nothing will be written later */
		}

		snprintf(filename, 1023, "g-iteration-%i.dat", i+1);
		fhp = fopen(filename, "w");
		if ( fhp == NULL ) {
			ERROR("Failed to open '%s'\n", filename);
			/* Nothing will be written later */
		}

		/* Refine the geometry of all patterns to get the best fit */
		refine_all(images, n_total_patterns, det, sym, obs, full,
		           nthreads, fhg, fhp);

		/* Re-estimate all the full intensities */
		reflist_free(full);
		select_scalable_reflections(images, n_usable_patterns);
		full = scale_intensities(images, n_usable_patterns,
		                         sym, obs, cref);

		fclose(fhg);
		fclose(fhp);

	}

	STATUS("Final scale factors:\n");
	for ( i=0; i<n_usable_patterns; i++ ) {
		STATUS("%4i : %5.2f\n", i, images[i].osf);
	}

	/* Output results */
	write_reflist(outfile, full, images[0].indexed_cell);

	/* Clean up */
	for ( i=0; i<n_usable_patterns; i++ ) {
		reflist_free(images[i].reflections);
		reflist_free(images[i].raw_reflections);
	}
	reflist_free(full);
	delete_items(obs);
	free(sym);
	free(outfile);
	free_detector_geometry(det);
	free(beam);
	free(cref);
	for ( i=0; i<n_usable_patterns; i++ ) {
		cell_free(images[i].indexed_cell);
		free(images[i].filename);
	}
	free(images);

	return 0;
}
