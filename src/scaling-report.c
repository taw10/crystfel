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

#include "image.h"


static void write_title(cairo_t *cr, const char *filename, double w, double h)
{
	char text[1024];
	PangoLayout *layout;
	PangoFontDescription *fontdesc;
	int width, height;

	snprintf(text, 1023, "Scaling report: %s", filename);

	layout = pango_cairo_create_layout(cr);
	pango_layout_set_text(layout, text, -1);
	fontdesc = pango_font_description_from_string("Sans 14 Bold");
	pango_layout_set_font_description(layout, fontdesc);

	pango_cairo_update_layout(cr, layout);
	pango_layout_get_size(layout, &width, &height);

	cairo_move_to(cr, 0.5-width/PANGO_SCALE, 10.0);

	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	pango_cairo_show_layout(cr, layout);
}


void scaling_report(const char *filename, const struct image *images, int n,
                    const char *stream_filename)
{
	cairo_surface_t *surface;
	cairo_t *cr;
	const double w = 842.0;
	const double h = 595.0;

	surface = cairo_pdf_surface_create(filename, w, h);

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return;
	}

	cr = cairo_create(surface);

	write_title(cr, stream_filename, w, h);

	cairo_surface_finish(surface);
	cairo_destroy(cr);
}
