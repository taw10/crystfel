/*
 * scaling-report.c
 *
 * Write a nice PDF of scaling parameters
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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


static void partiality_graph(cairo_t *cr, const struct image *images, int n,
                             RefList *full)
{
	const double g_width = 200.0;
	const double g_height = 200.0;
	int i;
	const int nbins = 25;
	double totals[nbins];
	int counts[nbins];
	double prob;

	show_text_simple(cr, "Observed partiality", -20.0, g_height/2.0,
	                      NULL, -M_PI_2, J_CENTER);
	show_text_simple(cr, "Calculated partiality", g_width/2.0, g_height+20.0,
	                      NULL, 0.0, J_CENTER);

	show_text_simple(cr, "0.0", -20.0, g_height, NULL, 0.0, J_CENTER);
	show_text_simple(cr, "1.0", -20.0, 0.0, NULL, 0.0, J_CENTER);
	show_text_simple(cr, "0.0", 0.0, g_height+10.0, NULL,
	                     -M_PI/3.0, J_RIGHT);
	show_text_simple(cr, "1.0", g_width, g_height+10.0, NULL,
	                     -M_PI/3.0, J_RIGHT);

	for ( i=0; i<nbins; i++ ) {
		totals[i] = 0.0;
		counts[i] = 0;
	}

	cairo_set_source_rgb(cr, 0.0, 0.7, 0.0);
	prob = 1.0 / n;
	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;

		if ( images[i].pr_dud ) continue;

		for ( refl = first_refl(images[i].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double Ipart, Ifull, pobs, pcalc;
			signed int h, k, l;
			Reflection *f;
			int bin;

			if ( !get_scalable(refl) ) continue;

			get_indices(refl, &h, &k, &l);
			f = find_refl(full, h, k, l);
			if ( f == NULL ) continue;

			Ipart = get_intensity(refl);
			Ifull = get_intensity(f);

			pobs = Ipart/(images[i].osf*Ifull);
			pcalc = get_partiality(refl);

			bin = nbins * pcalc;
			totals[bin] += pobs;
			counts[bin]++;

			if ( random_flat(1.0) < prob ) {
				plot_point(cr, g_width, g_height, pcalc, pobs);
			}

		}

	}

	cairo_new_path(cr);
	cairo_rectangle(cr, 0.0, 0.0, g_width, g_height);
	cairo_clip(cr);

	cairo_new_path(cr);
	cairo_move_to(cr, 0.0, g_height);
	for ( i=0; i<nbins; i++ ) {
		if ( counts[i] == 0 ) continue;
		cairo_line_to(cr, g_width*((double)i+0.5)/nbins,
		                  g_height - g_height*(totals[i]/counts[i]));
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


static void partiality_histogram(cairo_t *cr, const struct image *images,
                                 int n, RefList *full, int calc)
{
	int f_max;
	int i, b;
	const int nbins = 100;
	int counts[nbins];
	const double g_width = 320.0;
	const double g_height = 180.0;
	char tmp[32];


	show_text_simple(cr, "Frequency", -15.0, g_height/2.0,
	                      NULL, -M_PI_2, J_CENTER);
	if ( calc ) {
		show_text_simple(cr, "Distribution of calculated partialities",
		                     g_width/2.0, -18.0, "Sans Bold 10", 0.0,
		                     J_CENTER);
		show_text_simple(cr, "Calculated partiality", g_width/2.0,
		                     g_height+12.0, NULL, 0.0, J_CENTER);
	} else {
		show_text_simple(cr, "Distribution of observed partialities",
		                     g_width/2.0, -18.0, "Sans Bold 10", 0.0,
		                     J_CENTER);
		show_text_simple(cr, "Observed partiality", g_width/2.0,
		                     g_height+12.0, NULL, 0.0, J_CENTER);
	}

	for ( b=0; b<nbins; b++ ) {
		counts[b] = 0;
	}

	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;

		if ( images[i].pr_dud ) continue;

		for ( refl = first_refl(images[i].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double Ipart, Ifull, pobs, pcalc;
			signed int h, k, l;
			Reflection *f;

			if ( !get_scalable(refl) ) continue;

			get_indices(refl, &h, &k, &l);
			f = find_refl(full, h, k, l);
			if ( f == NULL ) continue;

			Ipart = get_intensity(refl);
			Ifull = get_intensity(f);

			pobs = Ipart/(images[i].osf*Ifull);
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

	show_text_simple(cr, "0", -10.0, g_height, NULL, 0.0, J_RIGHT);
	snprintf(tmp, 31, "%i", f_max);
	show_text_simple(cr, tmp, -10.0, 0.0, NULL, 0.0, J_RIGHT);

	show_text_simple(cr, "0.0", 0.0, g_height+10.0,
	                     NULL, -M_PI/3.0, J_RIGHT);
	show_text_simple(cr, "1.0", g_width, g_height+10.0,
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
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}


static void scale_factor_histogram(cairo_t *cr, const struct image *images,
                                   int n, const char *title)
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
	show_text_simple(cr, "Inverse scale factor", g_width/2.0, g_height+12.0,
	                      NULL, 0.0, J_CENTER);

	osf_max = 0.0;
	for ( i=0; i<n; i++ ) {
		double osf = images[i].osf;
		if ( osf > osf_max ) osf_max = osf;
	}
	osf_max = ceil(osf_max);
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

			double osf = images[i].osf;

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
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}


static void intensity_histogram(cairo_t *cr, const struct image *images,
                                int n, signed int h, signed int k, signed int l)
{
	int f_max;
	int i, b;
	const int nbins = 100;
	double int_max, int_inc;
	double int_low[nbins];
	double int_high[nbins];
	int counts[nbins];
	const double g_width = 200.0;
	const double g_height = 100.0;
	char tmp[64];

	snprintf(tmp, 63, "%i  %i  %i", h, k, l);
	show_text_simple(cr, tmp, g_width/2.0, -18.0,
	                      "Sans Bold 10", 0.0, J_CENTER);

	show_text_simple(cr, "Frequency", -15.0, g_height/2.0,
	                      NULL, -M_PI_2, J_CENTER);
	show_text_simple(cr, "Full scaled intensity",
	                      g_width/2.0, g_height+12.0,
	                      NULL, 0.0, J_CENTER);

	int_max = 0.0;
	int nmeas = 0;
	for ( i=0; i<n; i++ ) {

		Reflection *f;
		double osf;

		if ( images[i].pr_dud ) continue;

		osf = images[i].osf;

		for ( f = find_refl(images[i].reflections, h, k, l);
		      f != NULL;
		      f = next_found_refl(f) )
		{
			double Iobs, pcalc, Ifull_est;

			if ( !get_scalable(f) ) continue;

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

		Reflection *f;
		double osf;

		if ( images[i].pr_dud ) continue;

		osf = images[i].osf;

		for ( f = find_refl(images[i].reflections, h, k, l);
		      f != NULL;
		      f = next_found_refl(f) )
		{
			double Iobs, pcalc, Ifull_est;

			if ( !get_scalable(f) ) continue;

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

	show_text_simple(cr, "0", -10.0, g_height, NULL, 0.0, J_RIGHT);
	snprintf(tmp, 31, "%i", f_max);
	show_text_simple(cr, tmp, -10.0, 0.0, NULL, 0.0, J_RIGHT);

	show_text_simple(cr, "0.00", 0.0, g_height+10.0,
	                     NULL, -M_PI/3.0, J_RIGHT);
	snprintf(tmp, 32, "%5.2f", int_max);
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
	cairo_set_line_width(cr, 1.5);
	cairo_stroke(cr);
}


static void watermark(struct _srcontext *sr)
{
	show_text(sr->cr, "Written by partialator from CrystFEL"
	                 " version "PACKAGE_VERSION, sr->h-15.0, J_RIGHT,
	                 "Sans 7");
}


SRContext *sr_header(const char *filename, const char *stream_filename,
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

	snprintf(tmp, 1023, "Scaling report: %s", stream_filename);
	show_text(sr->cr, tmp, 10.0, J_CENTER, "Sans Bold 16");
	snprintf(tmp, 1023, "partialator %s", cmdline);
	show_text(sr->cr, tmp, 45.0, J_LEFT, "Mono 7");
	watermark(sr);

	return sr;
}


void sr_before(SRContext *sr, struct image *images, int n, RefList *full)
{
	if ( sr == NULL ) return;

	cairo_save(sr->cr);
	cairo_translate(sr->cr, 75.0, 100.0);
	scale_factor_histogram(sr->cr, images, n, "Before refinement");
	cairo_translate(sr->cr, 60.0, 235.0);
	partiality_graph(sr->cr, images, n, full);
	cairo_restore(sr->cr);
}


static void find_most_sampled_reflections(RefList *list, signed int *h,
                                          signed int *k, signed int *l,
                                          int n)
{
	Reflection *refl;
	RefListIterator *iter;
	int *samples;

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


void sr_after(SRContext *sr, struct image *images, int n, RefList *full)
{
	int i;
	signed int h[9], k[9], l[9];

	if ( sr == NULL ) return;

	cairo_save(sr->cr);
	cairo_translate(sr->cr, 475.0, 100.0);
	scale_factor_histogram(sr->cr, images, n, "After refinement");
	cairo_translate(sr->cr, 60.0, 235.0);
	partiality_graph(sr->cr, images, n, full);
	cairo_restore(sr->cr);

	cairo_surface_show_page(sr->surf);
	watermark(sr);

	cairo_save(sr->cr);
	cairo_translate(sr->cr, 75.0, 50.0);
	partiality_histogram(sr->cr, images, n, full, 1);
	cairo_translate(sr->cr, 400.0, 0.0);
	partiality_histogram(sr->cr, images, n, full, 0);
	cairo_restore(sr->cr);

	cairo_surface_show_page(sr->surf);
	watermark(sr);

	find_most_sampled_reflections(full, h, k, l, 9);

	for ( i=0; i<9; i++ ) {

		int x, y;

		x = i % 3;
		y = i / 3;

		cairo_save(sr->cr);
		cairo_translate(sr->cr, 50.0+280.0*x, 50.0+180.0*y);
		intensity_histogram(sr->cr, images, n, h[i], k[i], l[i]);
		cairo_restore(sr->cr);

	}

	cairo_surface_finish(sr->surf);
	cairo_destroy(sr->cr);

	free(sr);
}
