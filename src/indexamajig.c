/*
 * indexamajig.c
 *
 * Find hits, index patterns, output hkl+intensity etc.
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
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

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "intensities.h"
#include "ewald.h"
#include "peaks.h"


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

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"no-index",           0, &config_noindex,     1},
		{"dump-found-peaks",   0, &config_dumpfound,   1},
		{"dirax",              0, &config_dirax,       1},
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
		int fom;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		image.features = NULL;
		image.molecule = NULL;

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

		fom = image_fom(&image);
		if ( fom > 0 ) {

			/* Calculate orientation matrix (by magic) */
			index_pattern(&image, config_noindex, config_dumpfound,
			              config_dirax);

			if ( image.molecule == NULL ) continue;

			/* View head-on (unit cell is tilted) */
			image.orientation.x = 0.0;
			image.orientation.y = 0.0;
			image.orientation.z = 0.0;
			image.orientation.w = 1.0;
			get_ewald(&image);

			/* Read h,k,l,I */
			output_intensities(&image);

			n_hits++;

		}

	} while ( rval != NULL );

	fclose(fh);

	STATUS("There were %i images.\n", n_images);
	STATUS("%i hits were found.\n", n_hits);

	return 0;
}
