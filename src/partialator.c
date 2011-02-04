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


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (100)


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
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
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
	double nominal_photon_energy = pargs->image->beam->photon_energy;
	struct hdfile *hdfile;
	struct cpeak *spots;
	int n, i;
	double dev, last_dev;

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

	double a, b, c, al, be, ga;
	cell_get_parameters(image->indexed_cell, &a, &b, &c, &al, &be, &ga);
	STATUS("Initial cell: %5.2f %5.2f %5.2f nm %5.2f %5.2f %5.2f deg\n",
	       a/1.0e-9, b/1.0e-9, c/1.0e-9,
	       rad2deg(al), rad2deg(be), rad2deg(ga));

	spots = find_intersections(image, image->indexed_cell, &n, 0);
	dev = +INFINITY;
	i = 0;
	do {
		last_dev = dev;
		dev = pr_iterate(image, pargs->i_full, pargs->sym, &spots, &n);
		STATUS("Iteration %2i: mean dev = %5.2f\n", i, dev);
		i++;
	} while ( (fabs(last_dev - dev) > 1.0) && (i < MAX_CYCLES) );
	mean_partial_dev(image, spots, n, pargs->sym,
	                 pargs->i_full, pargs->graph);
	if ( pargs->pgraph ) {
		fprintf(pargs->pgraph, "%5i %5.2f\n", mytask, dev);
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


static void uniquify(struct cpeak *spot, const char *sym)
{
	signed int ha, ka, la;

	get_asymm(spot->h, spot->k, spot->l, &ha, &ka, &la, sym);
	spot->h = ha;
	spot->k = ka;
	spot->l = la;
}


static void integrate_image(struct image *image, ReflItemList *obs,
                            const char *sym)
{
	struct cpeak *spots;
	int j, n;
	struct hdfile *hdfile;
	double nominal_photon_energy = image->beam->photon_energy;

	hdfile = hdfile_open(image->filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open '%s'\n", image->filename);
		return;
	} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	if ( hdf5_read(hdfile, image, 0, nominal_photon_energy) ) {
		ERROR("Couldn't read '%s'\n", image->filename);
		hdfile_close(hdfile);
		return;
	}

	/* Figure out which spots should appear in this pattern */
	spots = find_intersections(image, image->indexed_cell, &n, 0);

	/* For each reflection, estimate the partiality */
	for ( j=0; j<n; j++ ) {

		signed int h, k, l;
		float i_partial;
		float xc, yc;

		uniquify(&spots[j], sym);

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
		                    &xc, &yc, &i_partial, NULL, NULL, 1, 0) ) {
			spots[j].valid = 0;
			continue;
		}

		spots[j].intensity = i_partial;
		spots[j].valid = 1;

		if ( !find_item(obs, h, k, l) ) add_item(obs, h, k, l);

	}
	image->cpeaks = spots;
	image->n_cpeaks = n;

	free(image->data);
	if ( image->flags != NULL ) free(image->flags);
	hdfile_close(hdfile);

	/* Muppet proofing */
	image->data = NULL;
	image->flags = NULL;
}


/* Decide which reflections can be scaled */
static void select_scalable_reflections(struct image *images, int n)
{
	int m;

	for ( m=0; m<n; m++ ) {

		int j;

		for ( j=0; j<images[m].n_cpeaks; j++ ) {

			int scalable = 1;

			if ( images[m].cpeaks[j].p < 0.1 ) scalable = 0;
			if ( !images[m].cpeaks[j].valid ) {
				scalable = 0;
			} else {
				double v = fabs(images[m].cpeaks[j].intensity);
				if ( v < 0.1 ) scalable = 0;
			}

			images[m].cpeaks[j].scalable = scalable;

		}

	}
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
			tmp = safe_basename(filename);
			free(filename);
			filename = tmp;
		}
		fnamereal = malloc(1024);
		snprintf(fnamereal, 1023, "%s%s", prefix, filename);
		free(filename);

		images[i].filename = fnamereal;
		images[i].div = beam->divergence;
		images[i].bw = beam->bandwidth;
		images[i].det = det;
		images[i].beam = beam;
		images[i].osf = 1.0;
		images[i].profile_radius = 0.005e9;

		/* Muppet proofing */
		images[i].data = NULL;
		images[i].flags = NULL;

		/* Get reflections from this image.
		 * FIXME: Use the ones from the stream */
		integrate_image(&images[i], obs, sym);

		progress_bar(i, n_total_patterns-1, "Loading pattern data");

	}
	fclose(fh);
	free(prefix);

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
		free(images[i].cpeaks);
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
