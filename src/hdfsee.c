/*
 * hdfsee.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2012      Richard Kirian
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

#include "version.h"
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
"      --version                    Print CrystFEL version number and exit.\n"
"\n"
"  -p, --peak-overlay=<filename>    Draw circles in positions listed in file.\n"
"      --ring-size=<n>              Set the size for those circles.\n"
"  -i, --int-boost=<n>              Multiply intensity by <n>.\n"
"  -b, --binning=<n>                Set display binning to <n>.\n"
"      --filter-noise               Apply an aggressive noise filter to the\n"
"                                    image data.\n"
"      --median-filter=<n>          Apply a median filter to the image data.\n"
"      --calibration-mode           Start in calibration mode\n"
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
"                                    HDF5 file. When this option is used,\n"
"                                    information about the data layout\n"
"                                    from the geometry file is ignored (See\n"
"                                    manual page).\n"
"                                    Example: /data/data0.\n"
"      --event=<event code>         Event to show from multi-event file.\n"
"  -g, --geometry=<filename>        Use geometry from file for display.\n"
"                                   (When this option is used, the value of\n"
"                                    of the -e parameter is ignored)"
"  -m, --beam=<filename>            Get beam parameters from <filename>.\n"
"  -o, --rigid-groups=<coll>        Use rigid group collection <coll>.\n"
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
	char *geom_filename = NULL;
	double boost = 1.0;
	int binning = 2;
	int config_noisefilter = 0;
	int config_showrings = 0;
	int config_calibmode =0;
	int colscale = SCALE_COLOUR;
	char *cscale = NULL;
	char *element = NULL;
	char *event = NULL;
	char *rgcoll_name = NULL;
	double ring_size = 5.0;
	char *reslist = NULL;
	double ring_radii[128];
	int n_rings = -1;
	int median_filter = 0;
	struct detector *det_geom = NULL;
	struct beam_params cbeam;
	struct beam_params *beam = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                4 },
		{"peak-overlay",       1, NULL,               'p'},
		{"int-boost",          1, NULL,               'i'},
		{"binning",            1, NULL,               'b'},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"colscale",           1, NULL,               'c'},
		{"image",              1, NULL,               'e'},
		{"geometry",           1, NULL,               'g'},
		{"show-rings",         0, &config_showrings,   1},
		{"ring-size",          1, NULL,                2},
		{"simple-rings",       1, NULL,               'r'},
		{"median-filter",      1, NULL,                3},
		{"calibration-mode",   0, &config_calibmode,   1},
		{"event",              1, NULL,                5},
		{"rigid-groups",       1, NULL,               'o'},
		{0, 0, NULL, 0}
	};

	/* Default beam parameters */
	cbeam.photon_energy = 0.0;
	cbeam.photon_energy_from = NULL;

	/* This isn't great, but necessary to make the command-line UI and file
	 * formats consistent with the other programs, which all use the C
	 * locale.  Better would be to have all the programs call
	 * setlocale(LC_ALL, "") and use the C locale temporarily when reading
	 * or writing a stream, reflection file, geometry file etc. */
	gtk_disable_setlocale();

	gtk_init(&argc, &argv);

	/* Short options */
	while ((c = getopt_long(argc, argv, "hp:b:i:c:e:g:2:r:m:o:",
	                        longopts, NULL)) != -1) {

		char *test;

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 4 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'p' :
			peaks = strdup(optarg);
			break;

			case 'i' :
			boost = atof(optarg);
			if ( boost <= 0 ) {
				ERROR("Intensity boost must be a positive"
				      " number.\n");
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
			geom_filename = strdup(optarg);
			det_geom = get_detector_geometry(geom_filename, &cbeam);
			if ( det_geom == NULL ) {
				ERROR("Failed to read detector geometry "
				      "from '%s'\n", optarg);
				return 1;
			}
			beam = &cbeam;
			break;

			case 'o' :
			rgcoll_name = strdup(optarg);
			break;

			case 2 :
			ring_size = strtod(optarg, &test);
			if ( test == optarg ) {
				ERROR("Ring size must be numerical.\n");
				return 1;
			}
			break;

			case 3 :
			median_filter = atoi(optarg);
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

			case 5 :
			event = strdup(optarg);
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;
		}

	}

	nfiles = argc-optind;

	if ( nfiles < 1 ) {
		ERROR("You need to give me a file to open!\n");
		return -1;
	}

	if ( (element != NULL) && (event != NULL) ) {
		ERROR("The options --event and --element are "
		      "mutually exclusive\n");
		return 1;
	}

	if ( event != NULL && geom_filename == NULL) {
		ERROR("The '--event' option requires geometry file\n");
		return 1;
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
		main_window_list[i] = displaywindow_open(argv[optind+i],
		                                         geom_filename,
		                                         peaks, boost, binning,
		                                         config_noisefilter,
		                                         config_calibmode,
		                                         colscale, element,
		                                         event, det_geom, beam,
		                                         rgcoll_name,
		                                         config_showrings,
		                                         ring_radii,
		                                         n_rings,
		                                         ring_size,
		                                         median_filter);
		if ( main_window_list[i] != NULL ) main_n_windows++;
	}

	if ( main_n_windows == 0 ) return 0;
	gtk_main();

	return 0;
}
