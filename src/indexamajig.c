/*
 * indexamajig.c
 *
 * Find hits, index patterns, output hkl+intensity etc.
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
#include <hdf5.h>

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "intensities.h"
#include "ewald.h"
#include "peaks.h"
#include "diffraction.h"
#include "detector.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Process and index FEL diffraction images.\n"
"\n"
"  -h, --help              Display this help message.\n"
"\n"
"  -i, --input=<filename>  Specify file containing list of images to process.\n"
"                           '-' means stdin, which is the default.\n"
"      --no-index          Do everything else (including fine peak search and\n"
"                           writing 'xfel.drx' if DirAx is being used), but\n"
"                           don't actually index.\n"
"      --dirax             Use DirAx for indexing.\n"
"      --dump-peaks        Write the results of the peak search to stdout.\n"
"      --near-bragg        Output a list of reflection intensities to stdout.\n"
"      --simulate          Simulate the diffraction pattern using the indexed\n"
"                           unit cell.\n"
"\n");
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	FILE *fh;
	char *rval;
	int n_images;
	int n_hits;
	int config_noindex = 0;
	int config_dumpfound = 0;
	int config_dirax = 0;
	int config_nearbragg = 0;
	int config_simulate = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"no-index",           0, &config_noindex,     1},
		{"dump-peaks",         0, &config_dumpfound,   1},
		{"near-bragg",         0, &config_nearbragg,   1},
		{"dirax",              0, &config_dirax,       1},
		{"simulate",           0, &config_simulate,    1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:w", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'i' : {
			filename = strdup(optarg);
			break;
		}

		case 0 : {
			break;
		}

		default : {
			return 1;
		}
		}

	}

	if ( filename == NULL ) {
		filename = strdup("-");
	}
	if ( strcmp(filename, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(filename, "r");
	}
	free(filename);
	if ( fh == NULL ) {
		ERROR("Failed to open input file\n");
		return 1;
	}

	n_images = 0;
	n_hits = 0;
	do {

		char line[1024];
		struct hdfile *hdfile;
		struct image image;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		image.features = NULL;
		image.molecule = NULL;
		image.data = NULL;

		STATUS("Processing '%s'\n", line);

		n_images++;

		hdfile = hdfile_open(line);
		if ( hdfile == NULL ) {
			continue;
		} else if ( hdfile_set_first_image(hdfile, "/") ) {
			ERROR("Couldn't select path\n");
			continue;
		}

		hdf5_read(hdfile, &image);

		/* Perform 'fine' peak search */
		search_peaks(&image);

		if ( image_feature_count(image.features) > 5 ) {

			n_hits++;

			if ( config_dumpfound ) dump_peaks(&image);

			/* Not indexing?  Then there's nothing left to do. */
			if ( config_noindex ) goto done;

			/* Calculate orientation matrix (by magic) */
			index_pattern(&image, config_noindex,
			              config_dirax);

			if ( image.molecule == NULL ) goto done;

			if ( config_nearbragg || config_simulate ) {

				/* Simulate a diffraction pattern */
				image.sfacs = NULL;
				image.data = NULL;
				image.qvecs = NULL;
				image.twotheta = NULL;
				image.hdr = NULL;

				/* View head-on (unit cell is tilted) */
				image.orientation.w = 1.0;
				image.orientation.x = 0.0;
				image.orientation.y = 0.0;
				image.orientation.z = 0.0;
				get_ewald(&image);

			}

			if ( config_nearbragg ) {
				/* Read h,k,l,I */
				output_intensities(&image);
			}

			if ( config_simulate ) {

				get_diffraction(&image, 8, 8, 8);
				if ( image.molecule == NULL ) {
					ERROR("Couldn't open molecule.pdb\n");
					return 1;
				}
				record_image(&image, 0, 0, 0);

				hdf5_write("simulated.h5", image.data,
				           image.width, image.height);

			}

		}

done:
		free(image.data);
		image_feature_list_free(image.features);
		hdfile_close(hdfile);
		H5close();

	} while ( rval != NULL );

	fclose(fh);

	STATUS("There were %i images.\n", n_images);
	STATUS("%i hits were found.\n", n_hits);

	return 0;
}
