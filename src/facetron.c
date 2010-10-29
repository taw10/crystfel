/*
 * facetron.c
 *
 * Profile fitting for coherent nanocrystallography
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


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Post-refinement and profile fitting for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"                              (must be a file, not e.g. stdin)\n"
"  -o, --output=<filename>    Output filename.  Default: facetron.hkl.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -b, --beam=<file>          Get beam parameters from file (provides initial\n"
"                              values for parameters, and nominal wavelengths\n"
"                              if no per-shot value is found in an HDF5 file.\n"
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"  -n, --iterations=<n>       Run <n> cycles of post-refinement.\n"
"\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


struct refine_args
{
	const char *sym;
	ReflItemList *obs;
	double *i_full;
	struct image *image;
};


static void refine_image(int mytask, void *tasks)
{
	struct refine_args *all_args = tasks;
	struct refine_args *pargs = &all_args[mytask];
	/* Do, er, something. */
}


struct integrate_args
{
	const char *sym;
	ReflItemList *obs;
	double *i_full;
	unsigned int *cts;
	pthread_mutex_t *list_lock;
	struct image *image;
};


static void integrate_image(int mytask, void *tasks)
{
	struct integrate_args *all_args = tasks;
	struct integrate_args *pargs = &all_args[mytask];
	struct cpeak *spots;
	int j, n;
	struct hdfile *hdfile;
	struct image *image = pargs->image;
	double nominal_photon_energy = pargs->image->beam->photon_energy;

	hdfile = hdfile_open(image->filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open '%s'\n", image->filename);
		return;
	} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	if ( hdf5_read(hdfile, pargs->image, 0, nominal_photon_energy) ) {
		ERROR("Couldn't read '%s'\n", image->filename);
		hdfile_close(hdfile);
		return;
	}

	/* Figure out which spots should appear in this pattern */
	spots = find_intersections(image, image->indexed_cell, &n, 0);

	/* For each reflection, estimate the partiality */
	for ( j=0; j<n; j++ ) {

		signed int h, k, l;
		signed int ha, ka, la;
		float i_partial;
		float xc, yc;
		float i_full_est;

		h = spots[j].h;
		k = spots[j].k;
		l = spots[j].l;

		/* Don't attempt to use spots with very small
		 * partialities, since it won't be accurate. */
		if ( spots[j].p < 0.1 ) continue;

		/* Actual measurement of this reflection from this
		 * pattern? */
		/* FIXME: Coordinates aren't whole numbers */
		if ( integrate_peak(image, spots[j].x, spots[j].y,
		                    &xc, &yc, &i_partial, NULL, NULL, 1, 1) ) {
			continue;
		}

		i_full_est = i_partial / spots[j].p;

		get_asymm(h, k, l, &ha, &ka, &la, pargs->sym);

		pthread_mutex_lock(pargs->list_lock);
		integrate_intensity(pargs->i_full, ha, ka, la, i_full_est);
		integrate_count(pargs->cts, ha, ka, la, 1);
		if ( !find_item(pargs->obs, ha, ka, la) ) {
			add_item(pargs->obs, ha, ka, la);
		}
		pthread_mutex_unlock(pargs->list_lock);

	}

	free(image->data);
	if ( image->flags != NULL ) free(image->flags);
	hdfile_close(hdfile);
	free(spots);

	/* Muppet proofing */
	image->data = NULL;
	image->flags = NULL;
}


static void refine_all(struct image *images, int n_total_patterns,
                       struct detector *det, const char *sym,
                       ReflItemList *obs, double *i_full, int nthreads)
{
	struct refine_args *tasks;
	int i;

	tasks = malloc(n_total_patterns * sizeof(struct refine_args));
	for ( i=0; i<n_total_patterns; i++ ) {

		tasks[i].sym = sym;
		tasks[i].obs = obs;
		tasks[i].i_full = i_full;
		tasks[i].image = &images[i];

	}

	run_thread_range(n_total_patterns, nthreads, "Refining",
	                 refine_image, tasks);

	free(tasks);
}


static void estimate_full(struct image *images, int n_total_patterns,
                          struct detector *det, const char *sym,
                          ReflItemList *obs, double *i_full, unsigned int *cts,
                          int nthreads)
{
	int i;
	struct integrate_args *tasks;
	pthread_mutex_t list_lock = PTHREAD_MUTEX_INITIALIZER;

	clear_items(obs);

	tasks = malloc(n_total_patterns * sizeof(struct integrate_args));
	for ( i=0; i<n_total_patterns; i++ ) {

		tasks[i].sym = sym;
		tasks[i].obs = obs;
		tasks[i].i_full = i_full;
		tasks[i].cts = cts;
		tasks[i].list_lock = &list_lock;
		tasks[i].image = &images[i];

	}

	run_thread_range(n_total_patterns, nthreads, "Integrating",
	                 integrate_image, tasks);

	free(tasks);

	/* Divide the totals to get the means */
	for ( i=0; i<num_items(obs); i++ ) {

		struct refl_item *it;
		double total;

		it = get_item(obs, i);
		total = lookup_intensity(i_full, it->h, it->k, it->l);
		total /= lookup_count(cts, it->h, it->k, it->l);
		set_intensity(i_full, it->h, it->k, it->l, total);

	}

	free(cts);
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *geomfile = NULL;
	char *prefix = NULL;
	char *sym = NULL;
	FILE *fh;
	int nthreads = 1;
	int config_basename = 0;
	int config_checkprefix = 1;
	struct detector *det;
	double *i_full;
	unsigned int *cts;
	ReflItemList *obs;
	int i;
	int n_total_patterns;
	struct image *images;
	int n_iter = 10;
	struct beam_params *beam = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"prefix",             1, NULL,               'x'},
		{"basename",           0, &config_basename,    1},
		{"no-check-prefix",    0, &config_checkprefix, 0},
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

		case 'x' :
			prefix = strdup(optarg);
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

	/* Sanitise prefix */
	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
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

	/* Prepare for iteration */
	i_full = new_list_intensity();
	obs = new_items();

	n_total_patterns = count_patterns(fh);
	STATUS("There are %i patterns to process\n", n_total_patterns);

	images = malloc(n_total_patterns * sizeof(struct image));
	if ( images == NULL ) {
		ERROR("Couldn't allocate memory for images.\n");
		return 1;
	}

	/* Fill in what we know about the images so far */
	rewind(fh);
	for ( i=0; i<n_total_patterns; i++ ) {

		UnitCell *cell;
		char *filename;
		char *fnamereal;

		if ( find_chunk(fh, &cell, &filename) == 1 ) {
			ERROR("Couldn't get all of the filenames and cells"
			      " from the input stream.\n");
			return 1;
		}

		images[i].indexed_cell = cell;

		/* Mangle the filename now */
		if ( config_basename ) {
			char *tmp;
			tmp = strdup(basename(filename));
			free(filename);
			filename = tmp;
		}
		fnamereal = malloc(1024);
		snprintf(fnamereal, 1023, "%s%s", prefix, filename);

		images[i].filename = fnamereal;
		images[i].div = beam->divergence;
		images[i].bw = beam->bandwidth;
		images[i].orientation.w = 1.0;
		images[i].orientation.x = 0.0;
		images[i].orientation.y = 0.0;
		images[i].orientation.z = 0.0;
		images[i].det = det;
		images[i].beam = beam;

		/* Muppet proofing */
		images[i].data = NULL;
		images[i].flags = NULL;

		free(filename);

		progress_bar(i, n_total_patterns-1, "Loading pattern data");

	}
	fclose(fh);
	free(prefix);

	cts = new_list_count();

	/* Make initial estimates */
	estimate_full(images, n_total_patterns, det, sym, obs, i_full, cts,
	              nthreads);

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		STATUS("Post refinement iteration %i of %i\n", i+1, n_iter);

		/* Refine the geometry of all patterns to get the best fit */
		refine_all(images, n_total_patterns, det, sym, obs, i_full,
		           nthreads);

		/* Re-estimate all the full intensities */
		estimate_full(images, n_total_patterns, det, sym, obs, i_full,
		              cts, nthreads);

	}

	/* Output results */
	write_reflections(outfile, obs, i_full, NULL, cts, NULL, 1.0);
	STATUS("Sigma(I) values in output file are not (yet) meaningful.\n");

	/* Clean up */
	free(i_full);
	delete_items(obs);
	free(sym);
	free(outfile);
	free(det->panels);
	free(det);
	for ( i=0; i<n_total_patterns; i++ ) {
		cell_free(images[i].indexed_cell);
		free(images[i].filename);
	}
	free(images);

	return 0;
}
