/*
 * main.c
 *
 * "Main"
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
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

#include "cell.h"
#include "image.h"
#include "utils.h"
#include "hdf5-file.h"
#include "templates.h"


static void main_show_help(const char *s)
{
	printf("Syntax: %s <file1.h5> <file2.h5> [...]\n\n", s);
	printf("Index diffraction patterns\n\n");
	printf("  -h              Display this help message\n");
}


int main(int argc, char *argv[])
{
	int c;
	char **in_files;
	size_t nin;
	size_t i;
	UnitCell *cell;
	TemplateList *templates;
	struct image template_parameters;

	while ((c = getopt(argc, argv, "h")) != -1) {

		switch ( c ) {

			case 'h' : {
				main_show_help(argv[0]);
				return 0;
			}

		}

	}

	if ( optind < argc ) {
		nin = argc-optind;
		in_files = malloc(nin*sizeof(char *));
		for ( i=0; i<nin; i++ ) {
			in_files[i] = strdup(argv[optind+i]);
		}
	} else {
		fprintf(stderr, "No input files!\n");
		return 1;
	}

	/* Define unit cell */
	cell = cell_new_from_parameters(28.10e-9,
	                                28.10e-9,
	                                16.52e-9,
	                                deg2rad(90.0),
	                                deg2rad(90.0),
	                                deg2rad(120.0));

	/* Generate templates */
	template_parameters.width = 512;
	template_parameters.height = 512;
	template_parameters.fmode = FORMULATION_CLEN;
	template_parameters.x_centre = 255.5;
	template_parameters.y_centre = 255.5;
	template_parameters.camera_len = 0.2;  /* 20 cm */
	template_parameters.resolution = 5120; /* 512 pixels in 10 cm */
	template_parameters.lambda = 0.2e-9;   /* LCLS wavelength */
	templates = generate_templates(cell, template_parameters);

	printf("Input files (%i):\n", nin);
	for ( i=0; i<nin; i++ ) {

		struct image image;

		printf("%6i: %s  ", i+1, in_files[i]);

		image.width = 512;
		image.height = 512;
		image.fmode = FORMULATION_CLEN;
		image.x_centre = 255.5;
		image.y_centre = 255.5;
		image.camera_len = 0.2;  /* 20 cm */
		image.resolution = 5120; /* 512 pixels in 10 cm */
		image.lambda = 0.2e-9;   /* LCLS wavelength */

		if ( hdf5_read(&image, in_files[i]) ) {
			fprintf(stderr, "Couldn't read file '%s'\n",
			        in_files[i]);
			continue;
		}

		try_templates(&image, templates);

		printf("%6.2f %6.2f\n", rad2deg(image.omega),
		                        rad2deg(image.tilt));

	}

	return 0;
}
