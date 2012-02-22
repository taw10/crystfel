/*
 * hdfsee.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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

#include <gtk/gtk.h>
#include <getopt.h>

#include "dw-hdfsee.h"
#include "utils.h"
#include "render.h"


/* Global program state */
DisplayWindow *main_window_list[64];
size_t main_n_windows = 0;


static void show_help(const char *s)
{
	printf("Syntax: %s [options] image.h5\n\n", s);
	printf(
"Quick HDF5 image viewer.\n"
"\n"
"  -h, --help                       Display this help message.\n"
"\n"
"  -p, --peak-overlay=<filename>    Draw circles in positions listed in file.\n"
"      --ring-size=<n>              Set the size for those circles.\n"
"  -i, --int-boost=<n>              Multiply intensity by <n>.\n"
"  -b, --binning=<n>                Set display binning to <n>.\n"
"      --filter-cm                  Perform common-mode noise subtraction.\n"
"      --filter-noise               Apply an aggressive noise filter which\n"
"                                    sets all pixels in each 3x3 region to\n"
"                                    zero if any of them have negative\n"
"                                    values.\n"
"      --show-rings                 Overlay rings that indicate resolution.\n"
"      --simple-rings=XX,YY,...     Overlay rings at specified radii XX, YY, ...\n"
"                                    in pixel units.\n"
"  -c, --colscale=<scale>           Use the given colour scale.  Choose from:\n"
"                                    mono    : Greyscale, black is zero.\n"
"                                    invmono : Greyscale, white is zero.\n"
"                                    colour  : Colour scale:\n"
"                                               black-blue-pink-red-orange-\n"
"                                               -yellow-white.\n"
"  -e, --image=<element>            Start up displaying this image from the\n"
"                                    HDF5 file.  Example: /data/data0.\n"
"  -g, --geometry=<filename>        Use geometry from file for display.\n"
"\n");
}


/* Called to notify that an image display window has been closed */
void hdfsee_window_closed(DisplayWindow *dw)
{
	size_t i;

	for ( i=0; i<main_n_windows; i++ ) {

		if ( main_window_list[i] == dw ) {

			size_t j;

			for ( j=i+1; j<main_n_windows; j++ ) {
				main_window_list[j] = main_window_list[j+1];
			}

		}

	}

	main_n_windows--;

	if ( main_n_windows == 0 ) gtk_exit(0);

}


int main(int argc, char *argv[])
{
	int c;
	size_t i;
	int nfiles;
	char *peaks = NULL;
	int boost = 1;
	int binning = 2;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	int config_showrings = 0;
	int colscale = SCALE_COLOUR;
	char *cscale = NULL;
	char *element = NULL;
	char *geometry = NULL;
	double ring_size = 5.0;
	char *reslist = NULL;
	double ring_radii[128];
	int n_rings = -1;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"peak-overlay",       1, NULL,               'p'},
		{"int-boost",          1, NULL,               'i'},
		{"binning",            1, NULL,               'b'},
		{"filter-cm",          0, &config_cmfilter,    1},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"colscale",           1, NULL,               'c'},
		{"image",              1, NULL,               'e'},
		{"geometry",           1, NULL,               'g'},
		{"show-rings",         0, &config_showrings,   1},
		{"ring-size",          1, NULL,                2},
		{"simple-rings",       1, NULL,               'r'},
		{0, 0, NULL, 0}
	};

	gtk_init(&argc, &argv);

	/* Short options */
	while ((c = getopt_long(argc, argv, "hp:b:i:c:e:g:2:r:",
	                        longopts, NULL)) != -1) {

		char *test;

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'p' :
			peaks = strdup(optarg);
			break;

		case 'i' :
			boost = atoi(optarg);
			if ( boost < 1 ) {
				ERROR("Intensity boost must be a positive"
				      " integer.\n");
				return 1;
			}
			break;

		case 'b' :
			binning = atoi(optarg);
			if ( binning < 1 ) {
				ERROR("Binning must be a positive integer.\n");
				return 1;
			}
			break;

		case 'c' :
			cscale = strdup(optarg);
			break;

		case 'e' :
			element = strdup(optarg);
			break;

		case 'g' :
			geometry = strdup(optarg);
			break;

		case 2 :
			ring_size = strtod(optarg, &test);
			if ( test == optarg ) {
				ERROR("Ring size must be numerical.\n");
				return 1;
			}
			break;

		case 'r' :
			config_showrings = 1;
			reslist = strdup(optarg);
			int nchar = strlen(reslist);
			char thisvalue[128];
			int i;
			int j = 0;
			n_rings = 0;
			for ( i=0; i<=nchar; i++ ) {
				if ( ( reslist[i] != ',' )
				  && ( reslist[i] != '\0' ) )
				{
					thisvalue[j] = reslist[i];
					j++;
				} else {
					thisvalue[j] = '\0';
					ring_radii[n_rings] = atof(thisvalue);
					n_rings++;
					j = 0;
				}
			}
			break;

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

	if ( cscale == NULL ) cscale = strdup("colour");
	if ( strcmp(cscale, "mono") == 0 ) {
		colscale = SCALE_MONO;
	} else if ( strcmp(cscale, "invmono") == 0 ) {
		colscale = SCALE_INVMONO;
	} else if ( strcmp(cscale, "colour") == 0 ) {
		colscale = SCALE_COLOUR;
	} else if ( strcmp(cscale, "color") == 0 ) {
		colscale = SCALE_COLOUR;
	} else {
		ERROR("Unrecognised colour scale '%s'\n", cscale);
		return 1;
	}
	free(cscale);

	for ( i=0; i<nfiles; i++ ) {
		main_window_list[i] = displaywindow_open(argv[optind+i], peaks,
		                                         boost, binning,
		                                         config_cmfilter,
		                                         config_noisefilter,
		                                         colscale, element,
		                                         geometry,
		                                         config_showrings,
		                                         ring_radii,
		                                         n_rings,
		                                         ring_size);
		if ( main_window_list[i] == NULL ) {
			ERROR("Couldn't open display window\n");
		} else {
			main_n_windows++;
		}
	}

	if ( main_n_windows == 0 ) return 0;
	gtk_main();

	return 0;
}
