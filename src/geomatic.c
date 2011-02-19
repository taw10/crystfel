/*
 * geomatic.c
 *
 * GUI geometry calibration
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <getopt.h>

#include "dw-geomatic.h"
#include "utils.h"
#include "render.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] image.h5\n\n", s);
	printf(
"GUI geometry calibration.\n"
"\n"
"  -h, --help                       Display this help message.\n"
"\n");
}


int main(int argc, char *argv[])
{
	int c;
	int nfiles;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{0, 0, NULL, 0}
	};

	gtk_init(&argc, &argv);

	/* Short options */
	while ((c = getopt_long(argc, argv, "hp:b:i:c:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	nfiles = argc-optind;

	if ( nfiles < 1 ) {
		ERROR("You need to give me a file to open!\n");
		return -1;
	}

	if ( geomatic_open(argv[optind]) == NULL ) {
		ERROR("Couldn't open display window\n");
		return 1;
	}
	gtk_main();

	return 0;
}
