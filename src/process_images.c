/*
 * process_images.c
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
#include "dirax.h"


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
"\n");
}


static int sum_of_peaks(struct image *image)
{
	int x, y;
	int integr = 0;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int val;

		val = image->data[x+image->height*y];

		if ( val > 1000 ) integr+=val;

	}
	}

	return integr;
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	FILE *fh;
	char *rval;
	int n_images;
	int n_hits;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:", longopts, NULL)) != -1) {

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
		int integr;

		rval = fgets(line, 1023, fh);
		chomp(line);

		STATUS("Processing '%s'\n", line);

		hdfile = hdfile_open(line);
		if ( hdfile == NULL ) {
			ERROR("Couldn't open file '%s'\n", filename);
		} else if ( hdfile_set_first_image(hdfile, "/") ) {
			ERROR("Couldn't select path\n");
		}

		hdf5_read(hdfile, &image);

		integr = sum_of_peaks(&image);
		printf("%6i %i\n", n_images, integr);
		if ( integr > 1e7 ) {

			STATUS("Hit: %s\n", line);

			index_pattern(&image);

			n_hits++;

		}

		n_images++;

	} while ( rval != NULL );

	fclose(fh);

	STATUS("There were %i images.\n", n_images);
	STATUS("%i hits were found.\n", n_hits);

	return 0;
}
