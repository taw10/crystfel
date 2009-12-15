/*
 * displaywindow.h
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

#ifndef DISPLAYWINDOW_H
#define DISPLAYWINDOW_H


typedef struct {
	GtkWidget	*window;
	GtkWidget	*entry;
} BinningDialog;


typedef struct {
	GtkWidget	*window;
	GtkWidget	*entry;
} BoostIntDialog;


typedef struct {

	GtkWidget	*window;
	GtkWidget	*drawingarea;
	GtkUIManager	*ui;
	GtkActionGroup	*action_group;
	GdkPixbuf	*pixbuf;

	struct hdfile	*hdfile;

	/* Dialog boxes */
	BinningDialog	*binning_dialog;
	BoostIntDialog	*boostint_dialog;

	int		width;
	int		height;		/* Size of the drawing area */
	int		binning;
	int		boostint;

	int		show_col_scale;
	int		monochrome;
	GdkPixbuf	*col_scale;

} DisplayWindow;

/* Open an image display window showing the given filename, or NULL */
extern DisplayWindow *displaywindow_open(const char *filename);


#endif	/* DISPLAYWINDOW_H */
