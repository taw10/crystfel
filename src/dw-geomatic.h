/*
 * dw-geomatic.h
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

#ifndef DW_GEOMATIC_H
#define DW_GEOMATIC_H

#include <gtk/gtk.h>

#include "cell.h"
#include "reflist.h"


struct gmdialog {
	GtkWidget	*window;
	GtkWidget	*entry;
};


typedef struct {

	GtkWidget        *window;
	GtkWidget        *drawingarea;
	GtkUIManager     *ui;
	GtkActionGroup   *action_group;
	GdkPixbuf        *pixbuf;
	gulong            motion_callback;

	struct hdfile    *hdfile;
	struct image     *image;
	int width;
	int height;
	double boostint;

	/* Dialog boxes */
	struct gmdialog   *boostint_dialog;

	int                show_col_scale;
	int                scale;
	GdkPixbuf         *col_scale;

	double             motion_origx;
	double             motion_origy;

	UnitCell          *cell;
	double             pos_x;
	double             pos_y;
	double             pos_z;

} DWGeomatic;

/* Return an image display window showing the given filename, or NULL */
extern DWGeomatic *geomatic_open(const char *filename);


#endif	/* DW_GEOMATICs_H */
