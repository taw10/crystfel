/*
 * displaywindow.h
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

#ifndef DISPLAYWINDOW_H
#define DISPLAYWINDOW_H

#include <gtk/gtk.h>


typedef struct {
	GtkWidget	*window;
	GtkWidget	*entry;
} BinningDialog;


typedef struct {
	GtkWidget	*window;
	GtkWidget	*entry;
} BoostIntDialog;


struct numberswindow {
	GtkWidget *window;
	GtkWidget *labels[17*17];
	GtkWidget *feat;
	unsigned int cx;
	unsigned int cy;
};


typedef struct {

	GtkWidget	*window;
	GtkWidget	*drawingarea;
	GtkUIManager	*ui;
	GtkActionGroup	*action_group;
	int             n_pixbufs;
	GdkPixbuf	**pixbufs;
	gulong		motion_callback;

	struct detector *loaded_geom;
	struct detector *simple_geom;

	struct hdfile	*hdfile;
	struct image	*image;

	/* Dialog boxes */
	BinningDialog	*binning_dialog;
	BoostIntDialog	*boostint_dialog;
	struct numberswindow *numbers_window;

	int		width;
	int		height;		/* Size of the drawing area */
	int		binning;
	double		boostint;
	int		cmfilter;	/* Use CM subtraction */
	int		noisefilter;	/* Use aggressive noise filter */
	int             use_geom;

	int		show_col_scale;
	int		scale;
	GdkPixbuf	*col_scale;

} DisplayWindow;

/* Open an image display window showing the given filename, or NULL */
extern DisplayWindow *displaywindow_open(const char *filename,
                                         const char *peaks, int boost,
                                         int binning, int cmfilter,
                                         int noisefilter, int colscale,
                                         const char *element);


#endif	/* DISPLAYWINDOW_H */
