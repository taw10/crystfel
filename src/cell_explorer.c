/*
 * cell_explorer.c
 *
 * Examine cell parameter histograms
 *
 * Copyright © 2014 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <gtk/gtk.h>
#include <math.h>
#include <gdk/gdkkeysyms.h>
#include <gsl/gsl_multifit_nlin.h>

#include "version.h"
#include "stream.h"
#include "image.h"
#include "utils.h"
#include "index.h"

#include "multihistogram.h"


static void show_help(const char *s)
{
	printf("Syntax: %s <my.stream>\n\n", s);
	printf(
"Examine cell parameter histograms.\n"
"\n"
" -h, --help              Display this help message.\n"
"     --version           Print CrystFEL version number and exit.\n"

);
}


#define CAT_P (0)
#define CAT_A (1)
#define CAT_B (2)
#define CAT_C (3)
#define CAT_I (4)
#define CAT_F (5)
#define CAT_H (6)
#define CAT_R (7)
#define CAT_EXCLUDE (8)


typedef struct {

	GtkWidget *da;
	MultiHistogram *h;
	double min;
	double max;
	int n;
	const char *units;
	const char *label;

	int width;
	double dmin;  /* Display min/max */
	double dmax;

	double sel1;
	double sel2;
	int show_sel;
	int erase_sel_on_release;

	double press_x;
	double press_y;
	double press_min;
	int sel;

	int have_fit;
	double fit_a;
	double fit_b;
	double fit_c;

	struct _cellwindow *parent;

} HistoBox;


typedef struct _cellwindow {

	GtkWidget *window;
	GtkUIManager *ui;
	GtkActionGroup *action_group;

	GtkWidget *indmlist;

	UnitCell **cells;
	IndexingMethod *indms;
	int n_cells;

	IndexingMethod unique_indms[256];
	int active_indms[256];
	int n_unique_indms;

	HistoBox *hist_a;
	HistoBox *hist_b;
	HistoBox *hist_c;
	HistoBox *hist_al;
	HistoBox *hist_be;
	HistoBox *hist_ga;

	int cols_on[8];

} CellWindow;


#define P_COL 0.0, 0.0, 0.0
#define A_COL 0.0, 0.8, 0.8
#define B_COL 0.0, 0.0, 0.8
#define C_COL 0.4, 0.3, 1.0
#define I_COL 0.0, 0.8, 0.0
#define F_COL 1.0, 0.3, 1.0
#define H_COL 0.8, 0.0, 0.0
#define R_COL 0.6, 0.6, 0.0

static void set_col(cairo_t *cr, CellWindow *w, int cat)
{
	if ( w->cols_on[cat] == 0 ) {
		cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
		return;
	}

	if ( w->cols_on[cat] == 1 ) {
		cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
		return;
	}

	switch ( cat ) {

		case CAT_P :  cairo_set_source_rgb(cr, P_COL);  break;
		case CAT_A :  cairo_set_source_rgb(cr, A_COL);  break;
		case CAT_B :  cairo_set_source_rgb(cr, B_COL);  break;
		case CAT_C :  cairo_set_source_rgb(cr, C_COL);  break;
		case CAT_I :  cairo_set_source_rgb(cr, I_COL);  break;
		case CAT_F :  cairo_set_source_rgb(cr, F_COL);  break;
		case CAT_H :  cairo_set_source_rgb(cr, H_COL);  break;
		case CAT_R :  cairo_set_source_rgb(cr, R_COL);  break;

	}
}


static gboolean destroy_sig(GtkWidget *da, CellWindow *w)
{
	gtk_main_quit();
	return FALSE;
}


static void redraw_all(CellWindow *cw)
{
	gtk_widget_queue_draw(cw->hist_a->da);
	gtk_widget_queue_draw(cw->hist_b->da);
	gtk_widget_queue_draw(cw->hist_c->da);
	gtk_widget_queue_draw(cw->hist_al->da);
	gtk_widget_queue_draw(cw->hist_be->da);
	gtk_widget_queue_draw(cw->hist_ga->da);
}


/* Calculate a round interval which occurs about three times in span */
static double calc_tic_interval(double span, int *nm)
{
	if ( span / 500.0 > 2 ) { *nm = 5;  return 500.0; }
	if ( span / 100.0 > 2 ) { *nm = 10;  return 100.0; }
	if ( span / 50.0 > 2 ) { *nm = 5;  return 50.0; }
	if ( span / 10.0 > 2 ) { *nm = 10;  return 10.0; }
	if ( span / 5.0 > 2 ) { *nm = 5;  return 5.0; }
	if ( span / 1.0 > 2 ) { *nm = 10;  return 1.0; }
	if ( span / 0.5 > 2 ) { *nm = 5;  return 0.5; }
	*nm = 10;  return 0.1;
}


static void draw_axis(cairo_t *cr, HistoBox *b, int width, int height)
{
	char label[128];
	double t, ti, mt;
	int nm;
	const double ws = width / (b->dmax-b->dmin);

	/* Draw the abscissa */
	cairo_move_to(cr, 0.0, height-20.5);
	cairo_line_to(cr, width, height-20.5);
	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_stroke(cr);

	/* Draw major tics (with labels) */
	t = calc_tic_interval(b->dmax-b->dmin, &nm);

	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.0);

	ti = t*trunc(b->dmin/t);
	for ( ; ti<=b->dmax; ti+=t ) {

		cairo_text_extents_t ext;

		cairo_move_to(cr, ws*(ti-b->dmin), height-20.0);
		cairo_line_to(cr, ws*(ti-b->dmin), height-12.0);
		cairo_stroke(cr);

		snprintf(label, 127, "%.1f%s", ti, b->units);
		cairo_text_extents(cr, label, &ext);
		cairo_move_to(cr, ws*(ti-b->dmin)-ext.x_advance/2, height-3.0);
		cairo_show_text(cr, label);

	}

	mt = t / nm;

	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.0);

	ti = mt*trunc(b->dmin/mt);
	for ( ; ti<=b->dmax; ti+=mt ) {
		cairo_move_to(cr, ws*(ti-b->dmin), height-20.0);
		cairo_line_to(cr, ws*(ti-b->dmin), height-18.0);
		cairo_stroke(cr);
	}
}


static gboolean draw_sig(GtkWidget *da, GdkEventExpose *event, HistoBox *b)
{
	int width, height;
	int i, max;
	double h_height;
	cairo_t *cr;
	double gstep;
	int *data_p, *data_a, *data_b, *data_c, *data_i, *data_f;
	int *data_r, *data_h, *data_excl;
	int start, stop;

	cr = gdk_cairo_create(da->window);

	width = da->allocation.width;
	height = da->allocation.height;
	b->width = width;  /* Store for later use when dragging */

	/* Overall background */
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_fill(cr);

	cairo_save(cr);

	cairo_translate(cr, 0.0, height);
	cairo_scale(cr, 1.0, -1.0);

	data_p = multihistogram_get_data(b->h, CAT_P);
	data_a = multihistogram_get_data(b->h, CAT_A);
	data_b = multihistogram_get_data(b->h, CAT_B);
	data_c = multihistogram_get_data(b->h, CAT_C);
	data_i = multihistogram_get_data(b->h, CAT_I);
	data_f = multihistogram_get_data(b->h, CAT_F);
	data_h = multihistogram_get_data(b->h, CAT_H);
	data_r = multihistogram_get_data(b->h, CAT_R);
	data_excl = multihistogram_get_data(b->h, CAT_EXCLUDE);

	max = 0;
	for ( i=0; i<b->n; i++ ) {
		int sum;
		sum = data_p[i] + data_a[i] + data_b[i] + data_c[i]
		    + data_i[i] + data_f[i] + data_h[i] + data_r[i]
		    + data_excl[i];
		if ( sum > max ) max = sum;
	}

	h_height = height - 20.0;
	cairo_translate(cr, 0.0, 20.0);

	cairo_scale(cr, width, h_height);
	gstep = (b->max-b->min)/b->n;

	/* Start with first visible bin */
	if ( b->dmin > b->min ) {
		start = (b->dmin - b->min)/gstep;
		start--;
	} else {
		start = 0;
	}
	if ( b->dmax < b->max ) {
		stop = (b->dmax - b->min)/gstep;
		stop++;
	} else {
		stop = b->n;
	}
	for ( i=start; i<stop; i++ ) {

		double hp = (double)data_p[i] / max;
		double ha = (double)data_a[i] / max;
		double hb = (double)data_b[i] / max;
		double hc = (double)data_c[i] / max;
		double hi = (double)data_i[i] / max;
		double hf = (double)data_f[i] / max;
		double hh = (double)data_h[i] / max;
		double hr = (double)data_r[i] / max;
		double he = (double)data_excl[i] / max;
		double x = b->min+i*gstep;
		double x2, w2;
		double s;

		x2 = (x - b->dmin)/(b->dmax - b->dmin);
		w2 = gstep/(b->dmax - b->dmin);

		s = 0.0;

		cairo_rectangle(cr, x2, s, w2, hp);
		set_col(cr, b->parent, CAT_P);
		cairo_fill(cr);
		s += hp;

		cairo_rectangle(cr, x2, s, w2, ha);
		set_col(cr, b->parent, CAT_A);
		cairo_fill(cr);
		s += ha;

		cairo_rectangle(cr, x2, s, w2, hb);
		set_col(cr, b->parent, CAT_B);
		cairo_fill(cr);
		s += hb;

		cairo_rectangle(cr, x2, s, w2, hc);
		set_col(cr, b->parent, CAT_C);
		cairo_fill(cr);
		s += hc;

		cairo_rectangle(cr, x2, s, w2, hi);
		set_col(cr, b->parent, CAT_I);
		cairo_fill(cr);
		s += hi;

		cairo_rectangle(cr, x2, s, w2, hf);
		set_col(cr, b->parent, CAT_F);
		cairo_fill(cr);
		s += hf;

		cairo_rectangle(cr, x2, s, w2, hh);
		set_col(cr, b->parent, CAT_H);
		cairo_fill(cr);
		s += hh;

		cairo_rectangle(cr, x2, s, w2, hr);
		set_col(cr, b->parent, CAT_R);
		cairo_fill(cr);
		s += hr;

		cairo_rectangle(cr, x2, s, w2, he);
		cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
		cairo_fill(cr);

	}

	if ( b->show_sel ) {
		cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 0.2);
		cairo_rectangle(cr, (b->sel1-b->dmin)/(b->dmax - b->dmin), 0.0,
		                    (b->sel2-b->sel1)/(b->dmax - b->dmin), 1.0);
		cairo_fill(cr);
	}

	/* Draw fitted curve */
	if ( b->have_fit ) {

		double A, B, C;
		A = b->fit_a / max;  B = b->fit_b;  C = b->fit_c;
		cairo_new_path(cr);

		/* In the current coordinate system, 0,0 is the bottom left
		 * of the graph (coordinates b->dmin, 0), and 1,1 is the top
		 * right of the graph (coordinate b->dmax, max) */
		for ( i=0; i<300; i++ ) {
			double xd = (double)i/300.0;
			double x = b->dmin + xd*(b->dmax - b->dmin);
			cairo_line_to(cr, xd, A*exp(-(x-B)*(x-B)/(C*C)));
		}
		cairo_set_source_rgba(cr, 1.0, 0.3, 0.0, 1.0);
		cairo_set_line_width(cr, 0.005);
		cairo_stroke(cr);
	}

	cairo_restore(cr);

	draw_axis(cr, b, width, height);

	cairo_text_extents_t ext;
	char label[256];

	cairo_select_font_face(cr, "Liberation Serif", CAIRO_FONT_SLANT_ITALIC,
	                       CAIRO_FONT_WEIGHT_BOLD);
	cairo_set_font_size(cr, height/10.0);

	if ( b->have_fit ) {
		snprintf(label, 255, "%s = %.2f ± %.2f%s",
		         b->label, b->fit_b, b->fit_c/sqrt(2), b->units);
	} else {
		strncpy(label, b->label, 255);
	}

	cairo_text_extents(cr, label, &ext);
	cairo_move_to(cr, 10.0, 10.0+ext.height);
	cairo_show_text(cr, label);

	cairo_destroy(cr);

	return FALSE;
}


static void centered_text(cairo_t *cr, double x, double sq, const char *t)
{
	cairo_text_extents_t ext;

	cairo_set_font_size(cr, sq-5.0);
	cairo_text_extents(cr, t, &ext);

	cairo_move_to(cr, x+(sq-ext.x_advance)/2.0, sq-(sq-ext.height)/2.0);
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);

	cairo_show_text(cr, t);
}


static gboolean keyconf_sig(GtkWidget *key, GdkEventConfigure *event,
                            CellWindow *w)
{
	gtk_widget_set_size_request(GTK_WIDGET(key), 8*event->height, -1);
	return FALSE;
}


static gint keyclick_sig(GtkWidget *widget, GdkEventButton *event,
                         CellWindow *w)
{
	int width, cat;

	/* Ignore extra events for double click */
	if ( event->type != GDK_BUTTON_PRESS ) return FALSE;

	width = widget->allocation.width;

	cat = 8*event->x / width;

	if ( cat == 0 ) {
		/* Special handling for P so that it doesn't go
		 * black->black->grey */
		w->cols_on[cat] = (w->cols_on[cat]+1) % 2;
	} else {
		w->cols_on[cat] = (w->cols_on[cat]+1) % 3;
	}

	gtk_widget_queue_draw(widget);
	redraw_all(w);

	return TRUE;
}



static gboolean keydraw_sig(GtkWidget *da, GdkEventExpose *event, CellWindow *w)
{
	int width, height;
	cairo_t *cr;
	double x;

	cr = gdk_cairo_create(da->window);

	width = da->allocation.width;
	height = da->allocation.height;

	/* Overall background */
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_fill(cr);

	x = 0.0;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_P);
	cairo_fill(cr);
	centered_text(cr, x, height, "P");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_A);
	cairo_fill(cr);
	centered_text(cr, x, height, "A");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_B);
	cairo_fill(cr);
	centered_text(cr, x, height, "B");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_C);
	cairo_fill(cr);
	centered_text(cr, x, height, "C");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_I);
	cairo_fill(cr);
	centered_text(cr, x, height, "I");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_F);
	cairo_fill(cr);
	centered_text(cr, x, height, "F");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_H);
	cairo_fill(cr);
	centered_text(cr, x, height, "H");

	x += height;
	cairo_rectangle(cr, x, 0.0, height, height);
	set_col(cr, w, CAT_R);
	cairo_fill(cr);
	centered_text(cr, x, height, "R");

	return FALSE;
}


static int check_exclude(HistoBox *h, double v)
{
	double min, max;

	if ( !h->show_sel ) return 0;

	if ( h->sel1 > h->sel2 ) {
		min = h->sel2;
		max = h->sel1;
	} else {
		min = h->sel1;
		max = h->sel2;
	}

	if ( v < min ) return 1;
	if ( v > max ) return 1;

	return 0;
}


static void scan_cells(CellWindow *w)
{
	int i;
	UnitCell **cells = w->cells;
	int n_cells = w->n_cells;

	multihistogram_delete_all_values(w->hist_a->h);
	multihistogram_delete_all_values(w->hist_b->h);
	multihistogram_delete_all_values(w->hist_c->h);
	multihistogram_delete_all_values(w->hist_al->h);
	multihistogram_delete_all_values(w->hist_be->h);
	multihistogram_delete_all_values(w->hist_ga->h);

	multihistogram_set_num_bins(w->hist_a->h, w->hist_a->n);
	multihistogram_set_num_bins(w->hist_b->h, w->hist_b->n);
	multihistogram_set_num_bins(w->hist_c->h, w->hist_c->n);
	multihistogram_set_num_bins(w->hist_al->h, w->hist_al->n);
	multihistogram_set_num_bins(w->hist_be->h, w->hist_be->n);
	multihistogram_set_num_bins(w->hist_ga->h, w->hist_ga->n);

	for ( i=0; i<n_cells; i++ ) {

		double a, b, c, al, be, ga;
		int cat, j;
		int ignore = 0;

		for ( j=0; j<w->n_unique_indms; j++ ) {
			if ( w->unique_indms[j] == w->indms[i] ) {
				if ( !w->active_indms[j] ) ignore = 1;
			}
		}

		if ( ignore ) continue;

		cell_get_parameters(cells[i], &a, &b, &c, &al, &be, &ga);
		a *= 1e10;  b *= 1e10;  c *= 1e10;
		al = rad2deg(al);  be = rad2deg(be);  ga = rad2deg(ga);

		switch ( cell_get_centering(cells[i]) ) {
			case 'P' : cat = 1<<CAT_P; break;
			case 'A' : cat = 1<<CAT_A; break;
			case 'B' : cat = 1<<CAT_B; break;
			case 'C' : cat = 1<<CAT_C; break;
			case 'I' : cat = 1<<CAT_I; break;
			case 'F' : cat = 1<<CAT_F; break;
			case 'H' : cat = 1<<CAT_H; break;
			case 'R' : cat = 1<<CAT_R; break;
			default : abort();
		}

		if ( check_exclude(w->hist_a, a) ||
		     check_exclude(w->hist_b, b) ||
		     check_exclude(w->hist_c, c) ||
		     check_exclude(w->hist_al, al) ||
		     check_exclude(w->hist_be, be) ||
		     check_exclude(w->hist_ga, ga) )
		{
			cat = 1<<CAT_EXCLUDE;
		}

		multihistogram_add_value(w->hist_a->h, a, cat);
		multihistogram_add_value(w->hist_b->h, b, cat);
		multihistogram_add_value(w->hist_c->h, c, cat);
		multihistogram_add_value(w->hist_al->h, al, cat);
		multihistogram_add_value(w->hist_be->h, be, cat);
		multihistogram_add_value(w->hist_ga->h, ga, cat);

	}
}


static void check_minmax(HistoBox *h, double val)
{
	if ( val > h->max ) h->max = val;
	if ( val < h->min ) h->min = val;
}


static void set_minmax(HistoBox *h)
{
	multihistogram_set_min(h->h, h->min);
	multihistogram_set_max(h->h, h->max);
}


static void scan_minmax(CellWindow *w)
{
	int i;

	for ( i=0; i<w->n_cells; i++ ) {

		double a, b, c, al, be, ga;
		int j;
		int found = 0;

		cell_get_parameters(w->cells[i], &a, &b, &c, &al, &be, &ga);
		a *= 1e10;  b *= 1e10;  c *= 1e10;
		al = rad2deg(al);  be = rad2deg(be);  ga = rad2deg(ga);

		check_minmax(w->hist_a, a);
		check_minmax(w->hist_b, b);
		check_minmax(w->hist_c, c);
		check_minmax(w->hist_al, al);
		check_minmax(w->hist_be, be);
		check_minmax(w->hist_ga, ga);

		for ( j=0; j<w->n_unique_indms; j++ ) {
			if ( w->indms[i] == w->unique_indms[j] ) {
				found = 1;
				break;
			}
		}
		if ( !found ) {
			if ( w->n_unique_indms > 255 ) {
				fprintf(stderr, "Too many indexing methods\n");
			} else {
				IndexingMethod m = w->indms[i];
				w->unique_indms[w->n_unique_indms] = m;
				w->active_indms[w->n_unique_indms++] = 1;
			}
		}

	}

	set_minmax(w->hist_a);
	set_minmax(w->hist_b);
	set_minmax(w->hist_c);
	set_minmax(w->hist_al);
	set_minmax(w->hist_be);
	set_minmax(w->hist_ga);
}


static void add_ui_sig(GtkUIManager *ui, GtkWidget *widget,
                       GtkContainer *container)
{
	gtk_box_pack_start(GTK_BOX(container), widget, FALSE, FALSE, 0);
	if ( GTK_IS_TOOLBAR(widget) ) {
		gtk_toolbar_set_show_arrow(GTK_TOOLBAR(widget), TRUE);
	}
}


static gint quit_sig(GtkWidget *widget, CellWindow *w)
{
	gtk_main_quit();
	return FALSE;
}


struct gaussian_data
{
	size_t n;
	double gstep;
	double min;
	int *data;
};


static int gaussian_f(const gsl_vector *p, void *vp,  gsl_vector *f)
{
	struct gaussian_data *d = vp;
	int i;

	double a = gsl_vector_get(p, 0);
	double b = gsl_vector_get(p, 1);
	double c = gsl_vector_get(p, 2);

	for ( i=0; i<d->n; i++ ) {
		double x = i*d->gstep + d->min;
		double ycalc = a*exp(-pow(x-b, 2)/(c*c));
		gsl_vector_set(f, i, ycalc - d->data[i]);
	}

	return GSL_SUCCESS;
}


static int gaussian_df(const gsl_vector *p, void *vp, gsl_matrix *J)
{
	struct gaussian_data *d = vp;
	double a = gsl_vector_get(p, 0);
	double b = gsl_vector_get(p, 1);
	double c = gsl_vector_get(p, 2);
	int i;

	for ( i=0; i<d->n; i++ ) {
		double x = i*d->gstep + d->min;
		double ycalc = a*exp(-pow(x-b, 2)/(c*c));
		gsl_matrix_set(J, i, 0, ycalc/a);
		gsl_matrix_set(J, i, 1, 2.0*ycalc*pow(c, -2)*(x-b));
		gsl_matrix_set(J, i, 2, 2.0*ycalc*pow(c, -3)*(x-b)*(x-b));
	}

	return GSL_SUCCESS;
}


static int gaussian_fdf(const gsl_vector *x, void *data, gsl_vector *f,
                        gsl_matrix *J)
{
	gaussian_f(x, data, f);
	gaussian_df(x, data, J);
	return GSL_SUCCESS;
}


static void fit_param(HistoBox *h)
{
	gsl_multifit_fdfsolver *s;
	gsl_vector *v;
	gsl_multifit_function_fdf f;
	struct gaussian_data params;
	double min, max, gstep;
	int n, n_iter, status, i, nm;
	int lowest_bin, highest_bin;
	int *data_p, *data_a, *data_b, *data_c, *data_i, *data_f;
	int *data_r, *data_h;

	if ( h->sel1 > h->sel2 ) {
		min = h->sel2;
		max = h->sel1;
	} else {
		min = h->sel1;
		max = h->sel2;
	}
	gstep = (h->max - h->min)/h->n;
	n = (max-min)/gstep;

	/* How many bins fit entirely within the selected range? */
	lowest_bin = floor((min - h->min)/gstep)+1;
	highest_bin = floor((max - h->min)/gstep)-1;
	params.n = 1 + highest_bin - lowest_bin;
	if ( highest_bin-lowest_bin < 3 ) {
		ERROR("Not enough bins.\n");
		return;
	}
	params.min = h->min + (lowest_bin+0.5)*gstep;
	params.gstep = gstep;

	data_p = multihistogram_get_data(h->h, CAT_P);
	data_a = multihistogram_get_data(h->h, CAT_A);
	data_b = multihistogram_get_data(h->h, CAT_B);
	data_c = multihistogram_get_data(h->h, CAT_C);
	data_i = multihistogram_get_data(h->h, CAT_I);
	data_f = multihistogram_get_data(h->h, CAT_F);
	data_h = multihistogram_get_data(h->h, CAT_H);
	data_r = multihistogram_get_data(h->h, CAT_R);
	params.data = malloc(params.n*sizeof(int));
	if ( params.data == NULL ) return;
	nm = 0;
	for ( i=0; i<params.n; i++ ) {
		int j = i+lowest_bin;
		if ( (j < 0) || (j >= h->n) ) {
			params.data[i] = 0;
			continue;
		}
		params.data[i] = data_p[j] + data_a[j] + data_b[j] + data_c[j]
		               + data_i[j] + data_f[j] + data_h[j] + data_r[j];
		if ( params.data[i] > nm ) nm = params.data[i];
	}

	s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, 3);

	v = gsl_vector_alloc(3);
	gsl_vector_set(v, 0, nm);
	gsl_vector_set(v, 1, min+(max-min)/2.0);
	gsl_vector_set(v, 2, (max-min)/5.0);

	f.f = gaussian_f;
	f.df = gaussian_df;
	f.fdf = gaussian_fdf;
	f.n = n;
	f.p = 3;
	f.params = &params;
	gsl_multifit_fdfsolver_set(s, &f, v);

	n_iter = 0;
	do {
		n_iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if ( status ) break;

		status = gsl_multifit_test_delta(s->dx, s->x, 0.001, 0.001);
	} while ( (status == GSL_CONTINUE) && (n_iter < 10));

	STATUS("Fitted: %.2f %.2f %.2f after %i iterations\n",
	       gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
	       gsl_vector_get(s->x, 2), n_iter);
	h->have_fit = 1;
	h->fit_a = gsl_vector_get(s->x, 0);
	h->fit_b = gsl_vector_get(s->x, 1);
	h->fit_c = gsl_vector_get(s->x, 2);

	free(params.data);
	gsl_multifit_fdfsolver_free(s);
	gsl_vector_free(v);
}


static gint fit_sig(GtkWidget *widget, CellWindow *w)
{
	if ( w->hist_a->show_sel ) fit_param(w->hist_a);
	if ( w->hist_b->show_sel ) fit_param(w->hist_b);
	if ( w->hist_c->show_sel ) fit_param(w->hist_c);
	if ( w->hist_al->show_sel ) fit_param(w->hist_al);
	if ( w->hist_be->show_sel ) fit_param(w->hist_be);
	if ( w->hist_ga->show_sel ) fit_param(w->hist_ga);

	redraw_all(w);

	return TRUE;
}

static gint about_sig(GtkWidget *widget, CellWindow *w)
{
	GtkWidget *window;

	const gchar *authors[] = {
		"Thomas White <taw@physics.org>",
		NULL
	};

	window = gtk_about_dialog_new();
	gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(w->window));

	gtk_about_dialog_set_name(GTK_ABOUT_DIALOG(window),
	        "Unit Cell Explorer");
	gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(window), "0.0.1");
	gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(window),
		"© 2014 Deutsches Elektronen-Synchrotron DESY, "
		"a research centre of the Helmholtz Association.");
	gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(window),
		"Examine unit cell distributions");
	gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(window),
		"http://www.desy.de/~twhite/crystfel");
	gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(window), authors);

	g_signal_connect(window, "response", G_CALLBACK(gtk_widget_destroy),
			 NULL);

	gtk_widget_show_all(window);

	return 0;
}


static void add_menu_bar(CellWindow *w, GtkWidget *vbox)
{
	GError *error = NULL;

	const char *ui = "<ui> <menubar name=\"cellwindow\">"
		"<menu name=\"file\" action=\"FileAction\">"
		"	<menuitem name=\"quit\" action=\"QuitAction\" />"
		"</menu>"
		"<menu name=\"tools\" action=\"ToolsAction\" >"
		"       <menuitem name=\"fit\" action=\"FitCellAction\" />"
		"</menu>"
		"<menu name=\"help\" action=\"HelpAction\">"
		"	<menuitem name=\"about\" action=\"AboutAction\" />"
		"</menu>"
		"</menubar></ui>";

	GtkActionEntry entries[] = {

		{ "FileAction", NULL, "_File", NULL, NULL, NULL },
		{ "QuitAction", GTK_STOCK_QUIT, "_Quit", NULL, NULL,
			G_CALLBACK(quit_sig) },

		{ "ToolsAction", NULL, "_Tools", NULL, NULL, NULL },
		{ "FitCellAction", NULL, "_Fit cell", "<Control>F", NULL,
			G_CALLBACK(fit_sig) },

		{ "HelpAction", NULL, "_Help", NULL, NULL, NULL },
		{ "AboutAction", GTK_STOCK_ABOUT, "_About", NULL, NULL,
			G_CALLBACK(about_sig) },

	};
	guint n_entries = G_N_ELEMENTS(entries);

	w->action_group = gtk_action_group_new("cellwindow");
	gtk_action_group_add_actions(w->action_group, entries, n_entries, w);

	w->ui = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group(w->ui, w->action_group, 0);
	g_signal_connect(w->ui, "add_widget", G_CALLBACK(add_ui_sig), vbox);
	if ( gtk_ui_manager_add_ui_from_string(w->ui, ui, -1, &error) == 0 )
	{
		fprintf(stderr, "Error loading message window menu bar: %s\n",
			error->message);
		return;
	}

	gtk_window_add_accel_group(GTK_WINDOW(w->window),
				   gtk_ui_manager_get_accel_group(w->ui));
	gtk_ui_manager_ensure_update(w->ui);
}


static void reset_axes(HistoBox *h)
{
	/* Fudge factor makes sure that the tic interval calculation falls
	 * clearly on one side or other of its tests */
	h->dmin = (nearbyint(h->min/10.0)-1.001)*10.0;
	h->dmax = (nearbyint(h->max/10.0)+1.001)*10.0;

}


static gint press_sig(GtkWidget *widget, GdkEventButton *event, HistoBox *h)
{
	h->press_x = event->x;
	h->press_y = event->y;
	h->press_min = h->dmin;
	if ( event->state & GDK_SHIFT_MASK ) {
		h->sel = 1;
		h->show_sel = 0;
	} else {
		h->sel = 0;
	}
	gtk_widget_grab_focus(GTK_WIDGET(h->da));
	return TRUE;
}


static gint release_sig(GtkWidget *widget, GdkEventButton *event, HistoBox *h)
{
	if ( h->sel ) {
		scan_cells(h->parent);
		redraw_all(h->parent);
	}
	return TRUE;
}


static gint motion_sig(GtkWidget *da, GdkEventMotion *event, HistoBox *h)
{
	double span = h->dmax - h->dmin;

	if ( !h->sel ) {

		h->dmin = h->press_min - span*(event->x - h->press_x)/h->width;
		h->dmax = h->dmin + span;

	} else {

		h->sel1 = h->dmin + span*h->press_x / h->width;
		h->sel2 = h->dmin + span*event->x / h->width;
		h->show_sel = 1;

	}

	gtk_widget_queue_draw(h->da);

	if ( event->is_hint ) gdk_window_get_pointer(da->window,
	                                             NULL, NULL, NULL);
	return TRUE;
}


static gint scroll_sig(GtkWidget *widget, GdkEventScroll *event, HistoBox *h)
{
	double span = h->dmax - h->dmin;
	double pos = h->dmin + span*event->x/h->width;;

	if ( event->direction == GDK_SCROLL_UP ) {
		h->dmin = pos - (pos-h->dmin)*0.9;
		h->dmax = pos + (h->dmax-pos)*0.9;
	} else if ( event->direction == GDK_SCROLL_DOWN ) {
		h->dmin = pos - (pos-h->dmin)*1.1;
		h->dmax = pos + (h->dmax-pos)*1.1;
	} else {
		return FALSE;
	}

	gtk_widget_grab_focus(GTK_WIDGET(h->da));
	gtk_widget_queue_draw(h->da);

	return TRUE;
}


static gint keypress_sig(GtkWidget *widget, GdkEventKey *event, HistoBox *h)
{
	if ( (event->keyval == GDK_R) || (event->keyval == GDK_r) ) {
		reset_axes(h);
		gtk_widget_queue_draw(h->da);
	}

	/* I'm too lazy to press shift */
	if ( (event->keyval == GDK_plus) || (event->keyval == GDK_equal) ) {
		if ( h->n < 100000 ) {
			h->n *= 2;
			scan_cells(h->parent);
			gtk_widget_queue_draw(h->da);
		}
	}

	if ( (event->keyval == GDK_minus) && (h->n > 1) ) {
		h->n /= 2;
		scan_cells(h->parent);
		gtk_widget_queue_draw(h->da);
	}

	return FALSE;
}


static HistoBox *histobox_new(CellWindow *w, const char *units, const char *n)
{
	HistoBox *h;

	h = calloc(1, sizeof(HistoBox));
	if ( h == NULL ) return NULL;

	h->show_sel = 0;
	h->units = units;
	h->parent = w;
	h->min = +INFINITY;
	h->max = -INFINITY;
	h->n = 100;  /* Number of bins */
	h->label = n;

	h->h = multihistogram_new();

	h->da = gtk_drawing_area_new();

	g_object_set(G_OBJECT(h->da), "can-focus", TRUE, NULL);
	gtk_widget_set_size_request(GTK_WIDGET(h->da), 400, 200);
	gtk_widget_add_events(GTK_WIDGET(h->da),
	                      GDK_BUTTON_PRESS_MASK
	                      | GDK_BUTTON_RELEASE_MASK
	                      | GDK_BUTTON1_MOTION_MASK
	                      | GDK_POINTER_MOTION_HINT_MASK
	                      | GDK_SCROLL_MASK
	                      | GDK_KEY_PRESS_MASK);

	g_signal_connect(G_OBJECT(h->da), "expose_event", G_CALLBACK(draw_sig),
	                 h);
	g_signal_connect(G_OBJECT(h->da), "button-press-event",
	                 G_CALLBACK(press_sig), h);
	g_signal_connect(G_OBJECT(h->da), "button-release-event",
	                 G_CALLBACK(release_sig), h);
	g_signal_connect(G_OBJECT(h->da), "motion-notify-event",
	                 G_CALLBACK(motion_sig), h);
	g_signal_connect(G_OBJECT(h->da), "scroll-event",
	                 G_CALLBACK(scroll_sig), h);
	g_signal_connect(G_OBJECT(h->da), "key-press-event",
	                 G_CALLBACK(keypress_sig), h);

	return h;
}


struct toggle_method
{
	int *active;
	CellWindow *w;
};


static gint indm_toggle_sig(GtkWidget *widget, struct toggle_method *tm)
{
	*tm->active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
	scan_cells(tm->w);
	redraw_all(tm->w);
	return FALSE;
}


static void indexing_method_list(CellWindow *w, GtkWidget *vbox)
{
	GtkWidget *key;
	int j;

	w->indmlist = gtk_hbox_new(FALSE, 5.0);
	gtk_box_pack_start(GTK_BOX(vbox), w->indmlist, FALSE, FALSE, 5.0);
	gtk_box_pack_start(GTK_BOX(w->indmlist),
	                   gtk_label_new("Show results from:"),
	                   FALSE, FALSE, 5.0);

	key = gtk_drawing_area_new();
	gtk_box_pack_end(GTK_BOX(w->indmlist), key, FALSE, FALSE, 5.0);
	gtk_widget_add_events(GTK_WIDGET(key), GDK_BUTTON_PRESS_MASK);
	g_signal_connect(G_OBJECT(key), "expose_event", G_CALLBACK(keydraw_sig),
	                 w);
	g_signal_connect(G_OBJECT(key), "configure-event",
	                 G_CALLBACK(keyconf_sig), w);
	g_signal_connect(G_OBJECT(key), "button-press-event",
	                 G_CALLBACK(keyclick_sig), w);

	for ( j=0; j<w->n_unique_indms; j++ ) {

		GtkWidget *button;
		char *label;
		struct toggle_method *tm = malloc(sizeof(struct toggle_method));

		if ( tm == NULL ) {
			fprintf(stderr, "Failed to allocate toggle method\n");
			continue;
		}

		label = indexer_str(w->unique_indms[j]);
		button = gtk_toggle_button_new_with_label(label);
		free(label);

		gtk_box_pack_start(GTK_BOX(w->indmlist), button,
	                           FALSE, FALSE, 5.0);

		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

		tm->w = w;
		tm->active = &w->active_indms[j];
		g_signal_connect(G_OBJECT(button), "toggled",
		                 G_CALLBACK(indm_toggle_sig), tm);

		w->active_indms[j] = 1;

	}


}


int main(int argc, char *argv[])
{
	int c;
	Stream *st;
	int max_cells = 0;
	char *stream_filename;
	GtkWidget *box, *vbox;
	char title[1024];
	char *bn;
	CellWindow w;
	int n_chunks = 0;
	int i;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                1 },
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "h",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 1 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			default :
			return 1;

		}

	}

	/* This isn't great, but necessary to make the command-line UI and file
	 * formats consistent with the other programs, which all use the C
	 * locale.  Better would be to have all the programs call
	 * setlocale(LC_ALL, "") and use the C locale temporarily when reading
	 * or writing a stream, reflection file, geometry file etc. */
	gtk_disable_setlocale();

	gtk_init(&argc, &argv);

	if ( argc != (optind+1) ) {
		fprintf(stderr, "Please provide exactly one stream filename.\n");
		return 1;
	}
	stream_filename = strdup(argv[optind]);
	st = open_stream_for_read(stream_filename);
	if ( st == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n", stream_filename);
		return 1;
	}

	w.cells = NULL;
	w.indms = NULL;
	w.n_cells = 0;
	do {

		struct image cur;
		int i;

		cur.det = NULL;

		if ( read_chunk_2(st, &cur, STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		for ( i=0; i<cur.n_crystals; i++ ) {

			Crystal *cr;

			cr = cur.crystals[i];

			if ( w.n_cells == max_cells ) {

				UnitCell **cells_new;
				IndexingMethod *indms_new;
				size_t nsz;

				nsz = (max_cells+1024)*sizeof(UnitCell *);
				cells_new = realloc(w.cells, nsz);
				if ( cells_new == NULL ) {
					fprintf(stderr, "Failed to allocate "
					        "memory for cells.\n");
					break;
				}

				nsz = (max_cells+1024)*sizeof(IndexingMethod);
				indms_new = realloc(w.indms, nsz);
				if ( indms_new == NULL ) {
					fprintf(stderr, "Failed to allocate "
					        "memory for methods.\n");
					break;
				}

				max_cells += 1024;
				w.cells = cells_new;
				w.indms = indms_new;

			}

			w.cells[w.n_cells] = crystal_get_cell(cr);
			w.indms[w.n_cells] = cur.indexed_by;
			w.n_cells++;

			crystal_free(cr);

		}

		n_chunks++;
		if ( n_chunks % 1000 == 0 ) {
			fprintf(stderr, "Loaded %i cells from %i chunks\r",
			        w.n_cells, n_chunks);
		}

	} while ( 1 );
	fprintf(stderr, "Loaded %i cells from %i chunks\n", w.n_cells, n_chunks);

	fprintf(stderr, "\n");

	close_stream(st);

	w.cols_on[0] = 1;
	for ( i=1; i<8; i++ ) w.cols_on[i] = 2;

	w.hist_a = histobox_new(&w, " A", "a");
	w.hist_b = histobox_new(&w, " A", "b");
	w.hist_c = histobox_new(&w, " A", "c");
	w.hist_al = histobox_new(&w, "°", "α");
	w.hist_be = histobox_new(&w, "°", "β");
	w.hist_ga = histobox_new(&w, "°", "γ");

	w.n_unique_indms = 0;

	scan_minmax(&w);
	scan_cells(&w);
	reset_axes(w.hist_a);
	reset_axes(w.hist_b);
	reset_axes(w.hist_c);
	reset_axes(w.hist_al);
	reset_axes(w.hist_be);
	reset_axes(w.hist_ga);

	w.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	bn = safe_basename(stream_filename);
	snprintf(title, 1023, "%s - Unit Cell Explorer", bn);
	free(bn);
	gtk_window_set_title(GTK_WINDOW(w.window), title);
	g_signal_connect(G_OBJECT(w.window), "destroy", G_CALLBACK(destroy_sig),
	                 &w);

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_container_add(GTK_CONTAINER(w.window), vbox);
	add_menu_bar(&w, vbox);

	indexing_method_list(&w, vbox);

	box = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), box, TRUE, TRUE, 5.0);

	gtk_box_pack_start(GTK_BOX(box), w.hist_a->da, TRUE, TRUE, 5.0);
	gtk_box_pack_start(GTK_BOX(box), w.hist_b->da, TRUE, TRUE, 5.0);
	gtk_box_pack_start(GTK_BOX(box), w.hist_c->da, TRUE, TRUE, 5.0);

	box = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), box, TRUE, TRUE, 5.0);

	gtk_box_pack_start(GTK_BOX(box), w.hist_al->da, TRUE, TRUE, 5.0);
	gtk_box_pack_start(GTK_BOX(box), w.hist_be->da, TRUE, TRUE, 5.0);
	gtk_box_pack_start(GTK_BOX(box), w.hist_ga->da, TRUE, TRUE, 5.0);

	gtk_widget_show_all(w.window);
	gtk_main();

	return 0;
}
