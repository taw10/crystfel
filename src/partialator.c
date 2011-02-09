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
#include "reflections.h"
#include "stream.h"
#include "geometry.h"
#include "peaks.h"
#include "thread-pool.h"
#include "beam-parameters.h"
#include "post-refinement.h"
#include "hrs-scaling.h"
#include "reflist.h"


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
	double *i_full;
	struct image *image;
	FILE *graph;
	FILE *pgraph;
};


static void refine_image(int mytask, void *tasks)
{
	struct refine_args *all_args = tasks;
	struct refine_args *pargs = &all_args[mytask];
	struct image *image = pargs->image;

	pr_refine(image, pargs->i_full, pargs->sym);
}


static void refine_all(struct image *images, int n_total_patterns,
                       struct detector *det, const char *sym,
                       ReflItemList *obs, double *i_full, int nthreads,
                       FILE *graph, FILE *pgraph)
{
	struct refine_args *tasks;
	int i;

	tasks = malloc(n_total_patterns * sizeof(struct refine_args));
	for ( i=0; i<n_total_patterns; i++ ) {

		tasks[i].sym = sym;
		tasks[i].obs = obs;
		tasks[i].i_full = i_full;
		tasks[i].image = &images[i];
		tasks[i].graph = graph;
		tasks[i].pgraph = pgraph;

	}

	run_thread_range(n_total_patterns, nthreads, "Refining",
	                 refine_image, tasks);

	free(tasks);
}


/* Decide which reflections can be scaled */
static void select_scalable_reflections(struct image *images, int n)
{
	int m;

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

		}

	}
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
	unsigned int *cts;
	ReflItemList *obs;
	int i;
	int n_total_patterns;
	struct image *images;
	int n_iter = 10;
	struct beam_params *beam = NULL;
	double *I_full;

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
		outfile = strdup("facetron.hkl");
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

		UnitCell *cell;
		char *filename;
		char *rval;
		char line[1024];
		RefList *peaks;
		RefList *transfer;
		Reflection *refl;
		RefListIterator *iter;
		double ph_en;

		if ( find_chunk(fh, &cell, &filename, &ph_en) == 1 ) {
			ERROR("Couldn't get all of the filenames, cells and"
			      "wavelengths from the input stream.\n");
			return 1;
		}

		images[i].indexed_cell = cell;

		images[i].filename = filename;
		images[i].div = beam->divergence;
		images[i].bw = beam->bandwidth;
		images[i].det = det;
		images[i].beam = beam;
		images[i].osf = 1.0;
		images[i].profile_radius = 0.005e9;
		images[i].reflections = reflist_new();
		images[i].lambda = ph_en_to_lambda(eV_to_J(ph_en));

		/* Muppet proofing */
		images[i].data = NULL;
		images[i].flags = NULL;

		/* Read integrated intensities from pattern */
		peaks = reflist_new();
		do {

			Reflection *refl;
			signed int h, k, l;
			float intensity;
			int r;

			rval = fgets(line, 1023, fh);
			chomp(line);

			if ( (strlen(line) == 0) || (rval == NULL) ) break;

			r = sscanf(line, "%i %i %i %f", &h, &k, &l, &intensity);
			if ( r != 4 ) continue;

			/* Not interested in the central beam */
			if ( (h==0) && (k==0) && (l==0) ) continue;

			refl = add_refl(peaks, h, k, l);
			set_int(refl, intensity);

		} while ( (strlen(line) != 0) && (rval != NULL) );

		/* Calculate initial partialities and fill in intensities from
		 * the stream */
		transfer = find_intersections(&images[i], cell, 0);
		images[i].reflections = reflist_new();

		for ( refl = first_refl(transfer, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			Reflection *peak;
			Reflection *new;
			signed int h, k, l, ha, ka, la;
			double r1, r2, p, x, y;
			int clamp1, clamp2;

			get_indices(refl, &h, &k, &l);
			STATUS("%3i %3i %3i\n", h, k, l);

			peak = find_refl(peaks, h, k, l);
			if ( peak == NULL ) {
				if ( (h==0) && (k==0) && (l==0) ) continue;
				STATUS("%3i %3i %3i not found\n", h, k, l);
				continue;
			}

			get_asymm(h, k, l, &ha, &ka, &la, sym);
			add_item(obs, ha, ka, la);
			new = add_refl(images[i].reflections, ha, ka, la);
			get_partial(refl, &r1, &r2, &p, &clamp1, &clamp2);
			get_detector_pos(refl, &x, &y);
			set_partial(new, r1, r2, p, clamp1, clamp2);
			set_detector_pos(new, 0.0, x, y);

		}
		reflist_free(peaks);
		reflist_free(transfer);

		progress_bar(i, n_total_patterns-1, "Loading pattern data");

	}
	fclose(fh);

	cts = new_list_count();

	/* Make initial estimates */
	STATUS("Performing initial scaling.\n");
	select_scalable_reflections(images, n_total_patterns);
	I_full = scale_intensities(images, n_total_patterns, sym, obs);

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
		refine_all(images, n_total_patterns, det, sym, obs, I_full,
		           nthreads, fhg, fhp);

		/* Re-estimate all the full intensities */
		free(I_full);
		select_scalable_reflections(images, n_total_patterns);
		I_full = scale_intensities(images, n_total_patterns, sym, obs);

		fclose(fhg);
		fclose(fhp);

	}

	STATUS("Final scale factors:\n");
	for ( i=0; i<n_total_patterns; i++ ) {
		STATUS("%4i : %5.2f\n", i, images[i].osf);
	}

	/* Output results */
	write_reflections(outfile, obs, I_full, NULL, NULL, cts, NULL);

	/* Clean up */
	for ( i=0; i<n_total_patterns; i++ ) {
		reflist_free(images[i].reflections);
	}
	free(I_full);
	delete_items(obs);
	free(sym);
	free(outfile);
	free(det->panels);
	free(det);
	free(beam);
	free(cts);
	for ( i=0; i<n_total_patterns; i++ ) {
		cell_free(images[i].indexed_cell);
		free(images[i].filename);
	}
	free(images);

	return 0;
}
