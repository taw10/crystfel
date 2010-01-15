/*
 * hdfsee.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <glib/gthread.h>
#include <getopt.h>

#include "displaywindow.h"
#include "utils.h"


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
"  -i, --int-boost=<n>        Multiple intensity by <n>.\n"
"  -b, --binning=<n>                Set display binning to <n>.\n"
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

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"peak-overlay",       1, NULL,               'p'},
		{"int-boost",          1, NULL,               'i'},
		{"binning",            1, NULL,               'b'},
		{0, 0, NULL, 0}
	};

	g_thread_init(NULL);
	gtk_init(&argc, &argv);

	/* Short options */
	while ((c = getopt_long(argc, argv, "hp:b:i:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'p' : {
			peaks = strdup(optarg);
			break;
		}

		case 'i' : {
			boost = atoi(optarg);
			if ( boost < 1 ) {
				ERROR("Intensity boost must be a positive"
				      " integer.\n");
			}
			break;
		}

		case 'b' : {
			binning = atoi(optarg);
			if ( boost < 1 ) {
				ERROR("Binning must be a positive integer.\n");
			}
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

	nfiles = argc-optind;

	if ( nfiles < 1 ) {
		ERROR("You need to give me a file to open!\n");
		return -1;
	}

	for ( i=0; i<nfiles; i++ ) {
		main_window_list[i] = displaywindow_open(argv[optind+i], peaks,
		                                         boost, binning);
		if ( main_window_list[i] == NULL ) {
			ERROR("Couldn't open display window\n");
		} else {
			main_n_windows++;
		}
	}
	gtk_main();

	return 0;
}
