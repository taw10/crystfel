/*
 * scaling-report.c
 *
 * Write a nice PDF of scaling parameters
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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

#include <cairo.h>
#include <cairo-pdf.h>
#include <pango/pangocairo.h>
#include <math.h>

#include "image.h"
#include "scaling-report.h"


#define PAGE_WIDTH (842.0)

enum justification
{
	J_CENTER,
	J_LEFT,
	J_RIGHT,
};


struct _srcontext
{
	cairo_surface_t *surf;
	cairo_t *cr;
	double w;
	double h;

	/* Most sampled reflections */
	signed int ms_h[9];
	signed int ms_k[9];
	signed int ms_l[9];

};


static void show_text(cairo_t *cr, const char *text, double y,
                      enum justification j, char *font)
{
	PangoLayout *layout;
	PangoFontDescription *fontdesc;
	int width, height;
	PangoAlignment just;

	if ( font == NULL ) font = "Sans 10";

	layout = pango_cairo_create_layout(cr);
	pango_layout_set_ellipsize(layout, PANGO_ELLIPSIZE_NONE);
	pango_layout_set_width(layout, PANGO_SCALE*(PAGE_WIDTH-20.0));

	switch ( j )
	{
		case J_CENTER : just = PANGO_ALIGN_CENTER; break;
		case J_LEFT : just = PANGO_ALIGN_LEFT; break;
		case J_RIGHT : just = PANGO_ALIGN_RIGHT; break;
		default: just = PANGO_ALIGN_LEFT; break;
	}

	pango_layout_set_alignment(layout, just);
	pango_layout_set_wrap(layout, PANGO_WRAP_CHAR);
	pango_layout_set_spacing(layout, 4.0*PANGO_SCALE);

	pango_layout_set_text(layout, text, -1);

	fontdesc = pango_font_description_from_string(font);
	pango_layout_set_font_description(layout, fontdesc);

	pango_cairo_update_layout(cr, layout);
	pango_layout_get_size(layout, &width, &height);

	cairo_move_to(cr, 10.0, y);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	pango_cairo_show_layout(cr, layout);
}


static void show_text_simple(cairo_t *cr, const char *text, double x, double y,
                             char *font, double rot, enum justification j)
{
	PangoLayout *layout;
	PangoFontDescription *fontdesc;
	int width, height;

	cairo_save(cr);

	if ( font == NULL ) font = "Sans 10";

	layout = pango_cairo_create_layout(cr);
	pango_layout_set_ellipsize(layout, PANGO_ELLIPSIZE_NONE);
	pango_layout_set_alignment(layout, PANGO_ALIGN_LEFT);

	pango_layout_set_text(layout, text, -1);

	fontdesc = pango_font_description_from_string(font);
	pango_layout_set_font_description(layout, fontdesc);

	pango_cairo_update_layout(cr, layout);
	pango_layout_get_size(layout, &width, &height);

	cairo_new_path(cr);
	cairo_translate(cr, x, y);
	cairo_rotate(cr, rot);
	if ( j == J_CENTER ) {
		cairo_translate(cr, -(width/2.0)/PANGO_SCALE,
		                    -(height/2.0)/PANGO_SCALE);
	} else if ( j == J_RIGHT ) {
		cairo_translate(cr, -width/PANGO_SCALE,
		                    -(height/2.0)/PANGO_SCALE);
	} else {
		cairo_translate(cr, 0.0, -(height/2.0)/PANGO_SCALE);
	}
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	pango_cairo_show_layout(cr, layout);

	cairo_restore(cr);
}


static void plot_point(cairo_t *cr, double g_width, double g_height,
                       double pcalc, double pobs)
{
	int bad = 0;

	if ( pobs > 1.0 ) {
		pobs = 1.01;
		bad = 1;

	}
	if ( pcalc > 1.0 ) {
		pcalc = 1.01;
		bad = 1;
	}
	if ( pobs < 0.0 ) {
		pobs = -0.01;
		bad = 1;
	}
	if ( pcalc < 0.0 ) {
		pobs = -0.01;
		bad = 1;
	}

	if ( bad ) {
		cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	}

	cairo_arc(cr, g_width*pcalc, g_height*(1.0-pobs),
		      1.0, 0.0, 2.0*M_PI);
	cairo_fill(cr);

	if ( bad ) {
		cairo_set_source_rgb(cr, 0.0, 0.7, 0.0);
	}
}


static void partiality_graph(cairo_t *cr, Crystal **crystals, int n,
                             RefList *full)
{
	const double g_width = 200.0;
	const double g_height = 200.0;
	int i;
	const int nbins = 25;
	double t_num[nbins];
	double t_den[nbins];
	double prob;
	double pcalcmin[nbins];
	double pcalcmax[nbins];
	int num_nondud;
	gsl_rng *rng;

	show_text_simple(cr, "Observed partiality", -20.0, g_height/2.0,
	                      NULL, -M_PI_2, J_CENTER);
	show_text_simple(cr, "Calculated partiality",
	                      g_width/2.0,g_height+20.0, NULL, 0.0, J_CENTER);

	show_text_simple(cr, "0.0", -20.0, g_height, NULL, 0.0, J_CENTER);
	show_text_simple(cr, "1.0", -20.0, 0.0, NULL, 0.0, J_CENTER);
	show_text_simple(cr, "0.0", 0.0, g_height+10.0, NULL,
	                     -M_PI/3.0, J_RIGHT);
	show_text_simple(cr, "1.0", g_width, g_height+10.0, NULL,
	                     -M_PI/3.0, J_RIGHT);

	for ( i=0; i<nbins; i++ ) {
		t_num[i] = 0.0;
		t_den[i] = 0.0;
		pcalcmin[i] = (double)i/nbins;
		pcalcmax[i] = (double)(i+1)/nbins;
	}
	pcalcmax[nbins-1] += 0.001;  /* Make sure it include pcalc = 1 */

	num_nondud = 0;
	for ( i=0; i<n; i++ ) {
		if ( crystal_get_user_flag(crystals[i]) ) continue;
		num_nondud++;
	}

	/* The reflections chosen for the graph will be the same every time
	 * (given the same sequence of input reflections, scalabilities etc) */
	rng = gsl_rng_alloc(gsl_rng_mt19937);

	cairo_set_source_rgb(cr, 0.0, 0.7, 0.0);
	prob = 1.0 / num_nondud;
	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;
		Crystal *cryst;

		cryst = crystals[i];

		if ( crystal_get_user_flag(cryst) ) continue;

		for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double Ipart, Ifull, pobs, pcalc;
			signed int h, k, l;
			Reflection *f;
			int bin;

			get_indices(refl, &h, &k, &l);
			f = find_refl(full, h, k, l);
			if ( f == NULL ) continue;
			if ( get_redundancy(f) < 2 ) continue;

			Ipart = crystal_get_osf(cryst) * get_intensity(refl);
			Ifull = get_intensity(f);

			pobs = Ipart/Ifull;
			pcalc = get_lorentz(refl) * get_partiality(refl);

			//STATUS("%4i %4i %4i : %9.6f %9.6f %e %e %e\n",
			//       h, k, l, pobs, pcalc,
			//       Ipart, Ifull, images[i].osf);

			for ( bin=0; bin<nbins; bin++ ) {
				if ( (pcalc >= pcalcmin[bin])
				  && (pcalc < pcalcmax[bin]) )
				{
					double esd_pobs, esd_Ip, esd_If;
					esd_Ip = get_esd_intensity(refl);
					esd_If = get_esd_intensity(f);
					esd_If *= crystal_get_osf(cryst);
					esd_pobs  = pow(esd_Ip/Ipart, 2.0);
					esd_pobs += pow(esd_If/Ifull, 2.0);
					esd_pobs = sqrt(esd_pobs);
					t_num[bin] += pobs / esd_pobs;
					t_den[bin] += 1.0 / esd_pobs;
				}
			}

			bin = nbins * pcalc;

			if ( random_flat(rng, 1.0) < prob ) {
				plot_point(cr, g_width, g_height, pcalc, pobs);
			}
		}

	}

	gsl_rng_free(rng);

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_clip(cr);

	cairo_new_path(cr);
	cairo_move_to(cr, 0.0, g_height);
	for ( i=0; i<nbins; i++ ) {

		double pos = pcalcmin[i] + (pcalcmax[i] - pcalcmin[i])/2.0;

		if ( t_den[i] == 0.0 ) continue;
		cairo_line_to(cr, g_width*pos,
		                  g_height - g_height*(t_num[i]/t_den[i]));

	}
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_stroke(cr);

	cairo_reset_clip(cr);

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}


static void partiality_histogram(cairo_t *cr, Crystal **crystals, int n,
                                 RefList *full, int calc, int backwards)
{
	int f_max;
	int i, b;
	const int nbins = 100;
	int counts[nbins];
	const double g_width = 200.0;
	const double g_height = 120.0;
	char tmp[32];
	double text_rot, axis_pos;

	if ( backwards ) {
		text_rot = M_PI_2;
		axis_pos = 215.0;
	} else {
		text_rot = -M_PI_2;
		axis_pos = -15.0;
	}

	show_text_simple(cr, "Frequency", axis_pos, g_height/2.0,
	                      NULL, text_rot, J_CENTER);

	for ( b=0; b<nbins; b++ ) {
		counts[b] = 0;
	}

	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;
		Crystal *cryst = crystals[i];

		if ( crystal_get_user_flag(cryst) ) continue;

		for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double Ipart, Ifull, pobs, pcalc;
			signed int h, k, l;
			Reflection *f;

			get_indices(refl, &h, &k, &l);
			f = find_refl(full, h, k, l);
			if ( f == NULL ) continue;

			Ipart = get_intensity(refl);
			Ifull = get_intensity(f);

			pobs = (Ipart * crystal_get_osf(cryst))
			           / (Ifull * get_lorentz(refl));
			pcalc = get_partiality(refl);

			if ( calc ) {
				b = pcalc*nbins;
			} else {
				b = pobs*nbins;
			}
			if ( (b>=0) && (b<nbins) ) counts[b]++;
		}
	}

	f_max = 0;
	for ( b=0; b<nbins; b++ ) {
		if ( counts[b] > f_max ) f_max = counts[b];
	}
	f_max = (f_max/10)*10 + 10;

	if ( !backwards ) {
		show_text_simple(cr, "0", axis_pos, g_height,
		                 NULL, 0.0, J_RIGHT);
	} else {
		show_text_simple(cr, "0", axis_pos, 0.0,
		                 NULL, M_PI, J_RIGHT);
	}
	snprintf(tmp, 31, "%i", f_max);
	if ( !backwards ) {
		show_text_simple(cr, tmp, axis_pos, 0.0,
		                 NULL, 0.0, J_RIGHT);
	} else {
		show_text_simple(cr, tmp, axis_pos, g_height,
		                 NULL, M_PI, J_RIGHT);
	}

	for ( b=0; b<nbins; b++ ) {

		double bar_height;

		bar_height = ((double)counts[b]/f_max)*g_height;

		cairo_new_path(cr);
		if ( !backwards ) {
			cairo_rectangle(cr, (g_width/nbins)*b, g_height,
			                    g_width/nbins, -bar_height);
		} else {
			cairo_rectangle(cr, (g_width/nbins)*b, 0.0,
			                    g_width/nbins, bar_height);
		}
		cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
		cairo_set_line_width(cr, 1.0);
		cairo_stroke(cr);

	}

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}


static void scale_factor_histogram(cairo_t *cr, Crystal **crystals, int n,
                                   const char *title)
{
	int f_max;
	int i, b;
	const int nbins = 100;
	double osf_max, osf_inc;
	double osf_low[nbins];
	double osf_high[nbins];
	int counts[nbins];
	const double g_width = 320.0;
	const double g_height = 180.0;
	char tmp[32];
	int n_zero, n_half;

	show_text_simple(cr, title, g_width/2.0, -18.0,
	                      "Sans Bold 10", 0.0, J_CENTER);

	show_text_simple(cr, "Frequency", -15.0, g_height/2.0,
	                      NULL, -M_PI_2, J_CENTER);
	show_text_simple(cr, "Scale factor", g_width/2.0, g_height+12.0,
	                      NULL, 0.0, J_CENTER);

	osf_max = 0.0;
	for ( i=0; i<n; i++ ) {
		double osf = crystal_get_osf(crystals[i]);
		if ( crystal_get_user_flag(crystals[i]) ) continue;
		if ( osf > osf_max ) osf_max = osf;
	}
	osf_max = ceil(osf_max+osf_max/10000.0);
	if ( osf_max > 1000.0 ) {
		ERROR("Silly scale factor detected.  Using 100.0 instead.\n");
		osf_max = 100.0;
	}

	do {

		osf_inc = osf_max / nbins;

		for ( b=0; b<nbins; b++ ) {
			osf_low[b] = b*osf_inc;
			osf_high[b] = (b+1)*osf_inc;
			counts[b] = 0;
		}

		for ( i=0; i<n; i++ ) {

			double osf = crystal_get_osf(crystals[i]);

			if ( crystal_get_user_flag(crystals[i]) ) continue;

			for ( b=0; b<nbins; b++ ) {
				if ( (osf >= osf_low[b])
				  && (osf < osf_high[b]) )
				{
					counts[b]++;
					break;
				}
			}
		}

		n_zero = 0;
		n_half = 0;
		if ( osf_max > 10.0 ) {

			/* Count the number of bins with no counts, subtract 1
			 * from the maximum value until this isn't the case */
			for ( b=0; b<nbins; b++ ) {
				if ( counts[b] == 0 ) n_zero++;
				n_half++;
			}
			n_half /= 2;

			if ( n_zero > n_half ) osf_max -= 1.0;

		}

	} while ( n_zero > n_half );

	f_max = 0;
	for ( b=0; b<nbins; b++ ) {
		if ( counts[b] > f_max ) f_max = counts[b];
	}
	f_max = (f_max/10)*10 + 10;

	show_text_simple(cr, "0", -10.0, g_height, NULL, 0.0, J_RIGHT);
	snprintf(tmp, 31, "%i", f_max);
	show_text_simple(cr, tmp, -10.0, 0.0, NULL, 0.0, J_RIGHT);

	show_text_simple(cr, "0.00", 0.0, g_height+10.0,
	                     NULL, -M_PI/3.0, J_RIGHT);
	snprintf(tmp, 32, "%5.2f", osf_max);
	show_text_simple(cr, tmp, g_width, g_height+10.0,
	                     NULL, -M_PI/3.0, J_RIGHT);

	for ( b=0; b<nbins; b++ ) {

		double bar_height;

		bar_height = ((double)counts[b]/f_max)*g_height;

		cairo_new_path(cr);
		cairo_rectangle(cr, (g_width/nbins)*b, g_height,
		                    g_width/nbins, -bar_height);
		cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
		cairo_set_line_width(cr, 1.0);
		cairo_stroke(cr);

	}

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);
}


static void intensity_histogram(cairo_t *cr, Crystal **crystals, int n,
                                RefList *full,
                                signed int h, signed int k, signed int l)
{
	int f_max;
	int i, b;
	const int nbins = 30;
	double int_max, int_inc;
	double int_low[nbins];
	double int_high[nbins];
	int counts[nbins];
	const double g_width = 115.0;
	const double g_height = 55.0;
	char tmp[64];
	Reflection *f;
	double Ifull, dsd_Ifull, pos, mI, bit;
	int have_full;

	f = find_refl(full, h, k, l);
	if ( f != NULL ) {
		Ifull = get_intensity(f);
		dsd_Ifull = get_esd_intensity(f) * sqrt(get_redundancy(f));
		have_full = 1;
	} else {
		Ifull = 0.0;
		dsd_Ifull = 0.0;
		have_full = 0;
	}

	snprintf(tmp, 63, "%i  %i  %i", h, k, l);
	show_text_simple(cr, tmp, g_width/2.0, -10.0,
	                      "Sans Bold 10", 0.0, J_CENTER);

	int_max = 0.0;
	int nmeas = 0;
	for ( i=0; i<n; i++ ) {

		double osf;
		Crystal *cryst = crystals[i];

		if ( crystal_get_user_flag(cryst) ) continue;

		osf = crystal_get_osf(cryst);

		for ( f = find_refl(crystal_get_reflections(cryst), h, k, l);
		      f != NULL;
		      f = next_found_refl(f) )
		{
			double Iobs, pcalc, Ifull_est;

			pcalc = get_partiality(f);
			Iobs = get_intensity(f);
			Ifull_est = Iobs / (pcalc * osf);

			if ( Ifull_est > int_max ) int_max = Ifull_est;
			nmeas++;
		}

	}
	int_max *= 1.1;
	int_inc = int_max / nbins;

	for ( b=0; b<nbins; b++ ) {
		int_low[b] = b*int_inc;
		int_high[b] = (b+1)*int_inc;
		counts[b] = 0;
	}

	for ( i=0; i<n; i++ ) {

		double osf;
		Crystal *cryst = crystals[i];

		if ( crystal_get_user_flag(cryst) ) continue;

		osf = crystal_get_osf(cryst);

		for ( f = find_refl(crystal_get_reflections(cryst), h, k, l);
		      f != NULL;
		      f = next_found_refl(f) )
		{
			double Iobs, pcalc, Ifull_est;

			pcalc = get_partiality(f);
			Iobs = get_intensity(f);
			Ifull_est = Iobs / (pcalc * osf);

			for ( b=0; b<nbins; b++ ) {
				if ( (Ifull_est >= int_low[b])
				  && (Ifull_est < int_high[b]) ) {
					counts[b]++;
					break;
				}
			}
		}

	}

	f_max = 0;
	for ( b=0; b<nbins; b++ ) {
		if ( counts[b] > f_max ) f_max = counts[b];
	}
	f_max = (f_max/10)*10 + 10;

	snprintf(tmp, 32, "Max I=%.0f", int_max);
	show_text_simple(cr, tmp, 10.0, 10.0, "Sans 9", 0.0, J_LEFT);

	for ( b=0; b<nbins; b++ ) {

		double bar_height;

		bar_height = ((double)counts[b]/f_max)*g_height;

		cairo_new_path(cr);
		cairo_rectangle(cr, (g_width/nbins)*b, g_height,
		                    g_width/nbins, -bar_height);
		cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
		cairo_set_line_width(cr, 1.0);
		cairo_stroke(cr);

	}

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);

	mI = int_high[nbins-1];
	pos = Ifull/mI;
	bit = g_height / 15.0;
	cairo_arc(cr, g_width*pos, g_height+bit, bit, 0.0, 2.0*M_PI);
	if ( have_full ) {
		cairo_set_source_rgb(cr, 0.0, 0.67, 0.45);
	} else {
		cairo_set_source_rgb(cr, 0.86, 0.0, 0.0);
	}
	cairo_fill(cr);

	if ( have_full ) {

		double eW = g_width*dsd_Ifull/mI;

		cairo_new_path(cr);
		cairo_rectangle(cr, 0.0, g_height+bit*2.0,
		                    g_width, g_height+bit*2.0);
		//cairo_clip(cr);

		cairo_new_path(cr);
		cairo_move_to(cr, g_width*pos - eW, g_height+bit);
		cairo_line_to(cr, g_width*pos + eW, g_height+bit);
		cairo_set_source_rgb(cr, 0.0, 0.67, 0.45);
		cairo_set_line_width(cr, 2.0);
		cairo_stroke(cr);

		cairo_reset_clip(cr);
	}
}


static void watermark(struct _srcontext *sr)
{
	show_text(sr->cr, "Written by partialator from CrystFEL"
	                 " version "PACKAGE_VERSION, sr->h-15.0, J_RIGHT,
	                 "Sans 7");
}


static void new_page(struct _srcontext *sr)
{
	cairo_show_page(sr->cr);
	watermark(sr);
}


static void find_most_sampled_reflections(RefList *list, int n, signed int *h,
                                          signed int *k, signed int *l)
{
	Reflection *refl;
	RefListIterator *iter;
	int *samples;
	int j;

	for ( j=0; j<n; j++ ) {
		h[j] = 0;
		k[j] = 0;
		l[j] = 0;
	}

	samples = calloc(n, sizeof(int));

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		int red;
		int i;

		red = get_redundancy(refl);

		for ( i=0; i<n; i++ ) {

			if ( red > samples[i] ) {

				int j;

				/* Shift everything down */
				for ( j=n-2; j>i; j-- ) {
					h[j+1] = h[j];
					k[j+1] = k[j];
					l[j+1] = l[j];
					samples[j+1] = samples[j];
				}

				/* Add this in its place */
				get_indices(refl, &h[i], &k[i], &l[i]);
				samples[i] = red;

				/* Don't compare against the others */
				break;

			}

		}

	}

	free(samples);
}



SRContext *sr_titlepage(Crystal **crystals, int n,
                        const char *filename, const char *stream_filename,
                        const char *cmdline)
{
	char tmp[1024];
	struct _srcontext *sr;

	sr = malloc(sizeof(*sr));
	if ( sr == NULL ) return NULL;

	sr->w = PAGE_WIDTH;
	sr->h = 595.0;

	sr->surf = cairo_pdf_surface_create(filename, sr->w, sr->h);

	if ( cairo_surface_status(sr->surf) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(sr->surf);
		free(sr);
		return NULL;
	}

	sr->cr = cairo_create(sr->surf);
	watermark(sr);

	snprintf(tmp, 1023, "%s", stream_filename);
	show_text(sr->cr, tmp, 10.0, J_CENTER, "Sans Bold 16");
	snprintf(tmp, 1023, "partialator %s", cmdline);
	show_text(sr->cr, tmp, 45.0, J_LEFT, "Mono 7");

	return sr;
}


void sr_iteration(SRContext *sr, int iteration, struct srdata *d)
{
	int i;
	char page_title[1024];
	double dash[] = {2.0, 2.0};

	if ( sr == NULL ) return;

	snprintf(page_title, 1023, "After %i iteration%s",
	         iteration, iteration==1?"":"s");

	new_page(sr);
	show_text(sr->cr, page_title, 10.0, J_CENTER, "Sans Bold 16");

	cairo_save(sr->cr);
	cairo_translate(sr->cr, 480.0, 350.0);
	scale_factor_histogram(sr->cr, d->crystals, d->n,
	                       "Distribution of overall scale factors");
	cairo_restore(sr->cr);

	/* Draw partiality plots (three graphs together) */
	cairo_save(sr->cr);

	cairo_translate(sr->cr, 70.0, 330.0);
	partiality_graph(sr->cr, d->crystals, d->n, d->full);

	cairo_save(sr->cr);
	cairo_move_to(sr->cr, 0.0, 0.0);
	cairo_line_to(sr->cr, 0.0, -30.0);
	cairo_move_to(sr->cr, 200.0, 0.0);
	cairo_line_to(sr->cr, 200.0, -30.0);
	cairo_set_dash(sr->cr, dash, 2, 0.0);
	cairo_stroke(sr->cr);
	cairo_set_dash(sr->cr, NULL, 0, 0.0);
	cairo_translate(sr->cr, 0.0, -150.0);
	partiality_histogram(sr->cr, d->crystals, d->n, d->full, 1, 0);
	cairo_restore(sr->cr);

	cairo_save(sr->cr);
	cairo_move_to(sr->cr, 200.0, 0.0);
	cairo_line_to(sr->cr, 230.0, 00.0);
	cairo_move_to(sr->cr, 200.0, 200.0);
	cairo_line_to(sr->cr, 230.0, 200.0);
	cairo_set_dash(sr->cr, dash, 2, 0.0);
	cairo_stroke(sr->cr);
	cairo_set_dash(sr->cr, NULL, 0, 0.0);
	cairo_translate(sr->cr, 230.0, 200.0);
	cairo_rotate(sr->cr, -M_PI_2);
	partiality_histogram(sr->cr, d->crystals, d->n, d->full, 0, 1);
	cairo_restore(sr->cr);

	cairo_restore(sr->cr);

	if ( iteration == 0 ) {
		find_most_sampled_reflections(d->full, 9,
		                              sr->ms_h, sr->ms_k, sr->ms_l);
	}

	for ( i=0; i<9; i++ ) {

		int x, y;

		x = i % 3;
		y = i / 3;

		cairo_save(sr->cr);
		cairo_translate(sr->cr, 400.0+140.0*x, 60.0+80.0*y);
		intensity_histogram(sr->cr, d->crystals, d->n, d->full,
		                    sr->ms_h[i], sr->ms_k[i], sr->ms_l[i]);
		cairo_restore(sr->cr);

	}

	STATUS("%i filtered total, %5.2f filtered per crystal\n",
	       d->n_filtered, (double)d->n_filtered / d->n);
}


void sr_finish(SRContext *sr)
{
	cairo_surface_finish(sr->surf);
	cairo_destroy(sr->cr);
}
