/*
 * dw-hdfsee.h
 *
 * Quick yet non-crappy HDF viewer
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *   2012      Richard Kirian
 *   2014      Valerio Mariani
 *   2014      Takanori Nakane <nakane.t@gmail.com>
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


typedef struct {
	GtkWidget	*window;
	GtkWidget	*entry;
} RingRadiusDialog;


struct numberswindow {
	GtkWidget *window;
	GtkWidget *labels[17*17];
	GtkWidget *feat;
	unsigned int cx;
	unsigned int cy;
};


typedef enum {
	CALIBMODE_NONE,
	CALIBMODE_PANELS,
	CALIBMODE_GROUPS,
	CALIBMODE_ALL
} CalibMode;


typedef struct {

	GtkWidget	*window;
	GtkWidget	*drawingarea;
	GtkWidget	*scrollarea;
	GtkUIManager	*ui;
	GtkActionGroup	*action_group;
	int             n_pixbufs;
	GdkPixbuf	**pixbufs;
	gulong		motion_callback;
	cairo_surface_t *surf;

	int             not_ready_yet;

	struct detector *loaded_geom;
	struct detector *simple_geom;

	struct hdfile	*hdfile;
	struct image	*image;

	/* Dialog boxes */
	BinningDialog  *binning_dialog;
	BoostIntDialog *boostint_dialog;
	RingRadiusDialog *ringradius_dialog;
	struct numberswindow *numbers_window;

	int		width;
	int		height;		/* Size of the drawing area */
	double          min_x;
	double          min_y;
	double          max_x;
	double          max_y;

	int		binning;
	double		boostint;
	int		noisefilter;	/* Use aggressive noise filter */
	int             median_filter;
	int             use_geom;
	int             show_rings;
	int		show_peaks;
	double          ring_radius;
	double          *ring_radii;
	int             n_rings;

	CalibMode       calib_mode;
	struct rigid_group *calib_mode_curr_rg;
	struct panel    *calib_mode_curr_p;
	int             calib_mode_show_focus;
	GtkWidget       *statusbar;

	int		show_col_scale;
	int		scale;
	GdkPixbuf	*col_scale;

} DisplayWindow;

/* Open an image display window showing the given filename, or NULL */
extern DisplayWindow *displaywindow_open(const char *filename,
                                         const char *peaks, double boost,
                                         int binning,
                                         int noisefilter, int calibmode, int colscale,
                                         const char *element,
                                         const char *geometry, const char *beam,
                                         int show_rings,
                                         double *ring_radii, int n_rings,
                                         double ring_size, int median_filter);


#endif	/* DISPLAYWINDOW_H */
