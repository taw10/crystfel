/*
 * hdfsee.c
 *
 * Quick yet non-crappy HDF viewer
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <glib/gthread.h>

#include "displaywindow.h"
#include "utils.h"


/* Global program state */
DisplayWindow *main_window_list[64];
size_t main_n_windows = 0;


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
	g_thread_init(NULL);
	gtk_init(&argc, &argv);

	if ( argc == 1 ) {

		main_n_windows = 1;
		main_window_list[0] = displaywindow_open(NULL);
		if ( main_window_list[0] == NULL ) {
			ERROR("Couldn't open display window\n");
		}

	} else {

		size_t i;

		main_n_windows = argc - 1;
		for ( i=0; i<main_n_windows; i++ ) {
			main_window_list[i] = displaywindow_open(argv[i+1]);
			if ( main_window_list[i] == NULL ) {
				ERROR("Couldn't open display window\n");
			}
		}

	}

	gtk_main();

	return 0;
}
