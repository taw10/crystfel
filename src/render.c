/*
 * render.c
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef HAVE_GTK
#include <gdk-pixbuf/gdk-pixbuf.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#ifdef HAVE_LIBPNG
#include <png.h>
#endif

#ifdef HAVE_TIFF
#include <tiffio.h>
#endif


#include "hdf5-file.h"
#include "render.h"
#include "peaks.h"
#include "filters.h"
#include "utils.h"


static void render_rgb(float val, float max, float *rp, float *gp, float *bp)
{
	int s;
	float p;
	float r, g, b;

	s = val / (max/6);
	p = fmod(val, max/6.0);
	p /= (max/6.0);

	r = 0.0;  g = 0.0;  b = 0.0;

	if ( (val < 0.0) ) {
		s = 0;
		p = 0;
	}
	if ( (val > max) ) {
		s = 6;
	}
	switch ( s ) {
		case 0 : {	/* Black to blue */
			r = 0;  g = 0;  b = p;
			break;
		}
		case 1 : {	/* Blue to pink */
			r = p;  g = 0;  b = 1.0;
			break;
		}
		case 2 : {	/* Pink to red */
			r = 1.0;  g = 0;  b = (1.0-p)*1.0;
			break;
		}
		case 3 : {	/* Red to Orange */
			r = 1.0;  g = 0.5*p;  b = 0;
			break;
		}
		case 4 : {	/* Orange to Yellow */
			r = 1.0;  g = 0.5 + 0.5*p;  b = 0;
			break;
		}
		case 5 : {	/* Yellow to White */
			r = 1.0;  g = 1.0;  b = 1.0*p;
			break;
		}
		case 6 : {	/* Pixel has hit the maximum value */
			r = 1.0;  g = 1.0;  b = 1.0;
			break;
		}
	}

	*rp = r;
	*gp = g;
	*bp = b;
}


static void render_mono(float val, float max, float *rp, float *gp, float *bp)
{
	float p;
	p = val / max;
	if ( val < 0.0 ) p = 0.0;
	if ( val > max ) p = 1.0;
	*rp = p;
	*gp = p;
	*bp = p;
}


static void render_invmono(float val, float max,
                           float *rp, float *gp, float *bp)
{
	float p;
	p = val / max;
	p = 1.0 - p;
	if ( val < 0.0 ) p = 1.0;
	if ( val > max ) p = 0.0;
	*rp = p;
	*gp = p;
	*bp = p;
}


void render_scale(float val, float max, int scale,
                  float *rp, float *gp, float *bp)
{
	switch ( scale ) {
	case SCALE_COLOUR :
		render_rgb(val, max, rp, gp, bp);
		break;
	case SCALE_MONO :
		render_mono(val, max, rp, gp, bp);
		break;
	case SCALE_INVMONO :
		render_invmono(val, max, rp, gp, bp);
		break;
	}
}


#ifdef HAVE_GTK

static float *get_binned_panel(struct image *image, int binning,
                               int min_fs, int max_fs, int min_ss, int max_ss)
{
	float *data;
	int x, y;
	int w, h;
	int fw;
	float *in;

	fw = image->width;
	in = image->data;

	/* Some pixels might get discarded */
	w = (max_fs - min_fs + 1) / binning;
	h = (max_ss - min_ss + 1) / binning;

	data = malloc(w*h*sizeof(float));

	for ( x=0; x<w; x++ ) {
	for ( y=0; y<h; y++ ) {

		double total;
		size_t xb, yb;

		total = 0;
		for ( xb=0; xb<binning; xb++ ) {
		for ( yb=0; yb<binning; yb++ ) {

			double v;
			v = in[binning*x+xb+min_fs + (binning*y+yb+min_ss)*fw];
			total += v;

		}
		}

		data[x+w*y] = total / ((double)binning * (double)binning);

	}
	}

	return data;
}


/* NB This function is shared between render_get_image() and
 * render_get_colour_scale() */
static void render_free_data(guchar *data, gpointer p)
{
	free(data);
}


static GdkPixbuf *render_panel(struct image *image,
                               int binning, int scale, double boost,
                               int min_fs, int max_fs, int min_ss, int max_ss)
{
	int w, h;
	guchar *data;
	float *hdr;
	int x, y;
	double max;
	int pw, ph;
	int i;

	/* Calculate panel width and height
	 * (add one because min and max are inclusive) */
	pw = max_fs - min_fs + 1;
	ph = max_ss - min_ss + 1;
	w = pw / binning;
	h = ph / binning;

	/* High dynamic range version */
	max = 0.0;
	for ( i=0; i<image->width*image->height; i++ ) {
		if ( image->data[i] > max ) max = image->data[i];
	}
	hdr = get_binned_panel(image, binning, min_fs, max_fs, min_ss, max_ss);
	if ( hdr == NULL ) return NULL;

	/* Rendered (colourful) version */
	data = malloc(3*w*h);
	if ( data == NULL ) {
		free(hdr);
		return NULL;
	}

	max /= boost;
	if ( max <= 6 ) { max = 10; }
	/* These x,y coordinates are measured relative to the bottom-left
	 * corner */
	for ( y=0; y<h; y++ ) {
	for ( x=0; x<w; x++ ) {

		float val;
		float r, g, b;

		val = hdr[x+w*y];
		render_scale(val, max, scale, &r, &g, &b);

		/* Stuff inside square brackets makes this pixel go to
		 * the expected location in the pixbuf (which measures
		 * from the top-left corner */
		data[3*( x+w*y )+0] = 255*r;
		data[3*( x+w*y )+1] = 255*g;
		data[3*( x+w*y )+2] = 255*b;

	}
	}

	/* Finished with this */
	free(hdr);

	/* Create the pixbuf from the 8-bit display data */
	return gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB, FALSE, 8,
					w, h, w*3, render_free_data, NULL);

}


/* Render an image into multiple pixbufs according to geometry */
GdkPixbuf **render_panels(struct image *image,
                          int binning, int scale, double boost,
                          int *n_pixbufs)
{
	int i;
	int np = image->det->n_panels;
	GdkPixbuf **pixbufs;

	pixbufs = calloc(np, sizeof(GdkPixbuf*));
	if ( pixbufs == NULL ) {
		*n_pixbufs = 0;
		return NULL;
	}

	for ( i=0; i<np; i++ ) {

		pixbufs[i] = render_panel(image, binning, scale, boost,
		                          image->det->panels[i].min_fs,
		                          image->det->panels[i].max_fs,
		                          image->det->panels[i].min_ss,
		                          image->det->panels[i].max_ss);

	}

	*n_pixbufs = np;

	return pixbufs;
}


GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int scale)
{
	guchar *data;
	size_t x, y;
	int max;

	data = malloc(3*w*h);
	if ( data == NULL ) return NULL;

	max = h;

	for ( y=0; y<h; y++ ) {

		float r, g, b;
		int val;

		val = y;

		render_scale(val, max, scale, &r, &g, &b);

		data[3*( 0+w*(h-1-y) )+0] = 0;
		data[3*( 0+w*(h-1-y) )+1] = 0;
		data[3*( 0+w*(h-1-y) )+2] = 0;
		for ( x=1; x<w; x++ ) {
			data[3*( x+w*(h-1-y) )+0] = 255*r;
			data[3*( x+w*(h-1-y) )+1] = 255*g;
			data[3*( x+w*(h-1-y) )+2] = 255*b;
		}

	}

	return gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB, FALSE, 8,
					w, h, w*3, render_free_data, NULL);
}


int render_tiff_fp(struct image *image, const char *filename)
{
#ifdef HAVE_TIFF
	TIFF *th;
	float *line;
	int y;

	th = TIFFOpen(filename, "w");
	if ( th == NULL ) return 1;

	TIFFSetField(th, TIFFTAG_IMAGEWIDTH, image->width);
	TIFFSetField(th, TIFFTAG_IMAGELENGTH, image->height);
	TIFFSetField(th, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(th, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	TIFFSetField(th, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(th, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(th, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(th, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(th, TIFFTAG_ROWSPERSTRIP,
	             TIFFDefaultStripSize(th, image->width*4));

	line = _TIFFmalloc(TIFFScanlineSize(th));
	for ( y=0; y<image->height; y++ ) {
		memcpy(line, &image->data[(image->height-1-y)*image->width],
		       image->width*4);
		TIFFWriteScanline(th, line, y, 0);
	}
	_TIFFfree(line);

	TIFFClose(th);

#else
	STATUS("No TIFF support.\n");
#endif
	return 0;
}


int render_tiff_int16(struct image *image, const char *filename, double boost)
{
#ifdef HAVE_TIFF
	TIFF *th;
	int16_t *line;
	int x, y;
	float max;

	th = TIFFOpen(filename, "w");
	if ( th == NULL ) return 1;

	TIFFSetField(th, TIFFTAG_IMAGEWIDTH, image->width);
	TIFFSetField(th, TIFFTAG_IMAGELENGTH, image->height);
	TIFFSetField(th, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(th, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT); /* (signed) */
	TIFFSetField(th, TIFFTAG_BITSPERSAMPLE, 16);
	TIFFSetField(th, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(th, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(th, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(th, TIFFTAG_ROWSPERSTRIP,
	             TIFFDefaultStripSize(th, image->width*4));

	line = _TIFFmalloc(TIFFScanlineSize(th));
	max = 0.0;
	for ( y=0; y<image->height; y++ ) {
	for ( x=0;x<image->width; x++ ) {
		float val;
		val = image->data[x+image->height*y];
		if ( val > max ) max = val;
	}
	}
	max /= 32767.0;

	for ( y=0; y<image->height; y++ ) {
		for ( x=0;x<image->width; x++ ) {

			float val;

			val = image->data[x+(image->height-1-y)*image->width];
			val *= ((float)boost/max);

			/* Clamp to 16-bit range,
			 * and work round inability of most readers to deal
			 * with signed integers. */
			val += 1000.0;
			if ( val > 32767.0 ) val = 32767.0;
			if ( val < 0.0 ) val = 0.0;

			line[x] = val;
		}

		TIFFWriteScanline(th, line, y, 0);
	}
	_TIFFfree(line);

	TIFFClose(th);
#else
	STATUS("No TIFF support.\n");
#endif
	return 0;
}

#endif /* HAVE_GTK */
