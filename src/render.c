/*
 * render.c
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <gdk-pixbuf/gdk-pixbuf.h>
#include <math.h>
#include <stdint.h>

#include "hdf5-file.h"
#include "render.h"
#include "peaks.h"
#include "filters.h"


static void *render_bin(int16_t *in, int inw, int inh,
                        int binning, int16_t *maxp)
{
	int16_t *data;
	int x, y;
	int w, h;
	int16_t max;

	w = inw / binning;
	h = inh / binning;      /* Some pixels might get discarded */

	data = malloc(w*h*sizeof(int16_t));
	max = 0;

	for ( x=0; x<w; x++ ) {
	for ( y=0; y<h; y++ ) {

		/* Big enough to hold large values */
		unsigned long long int total;
		size_t xb, yb;

		total = 0;
		for ( xb=0; xb<binning; xb++ ) {
		for ( yb=0; yb<binning; yb++ ) {

			total += in[binning*x+xb + (binning*y+yb)*inw];

		}
		}

		data[x+w*y] = total / (binning * binning);
		if ( data[x+w*y] > max ) max = data[x+w*y];

	}
	}

	*maxp = max;
	return data;
}


int16_t *render_get_image_binned(DisplayWindow *dw, int binning, int16_t *max)
{
	struct image *image;
	int16_t *data;

	if ( (dw->image == NULL) || (dw->image_dirty) ) {

		image = malloc(sizeof(struct image));
		if ( image == NULL ) return NULL;
		image->features = NULL;
		image->data = NULL;

		hdf5_read(dw->hdfile, image);
		dw->image_dirty = 0;
		clean_image(image);

		/* Deal with the old image, if existing */
		if ( dw->image != NULL ) {
			image->features = dw->image->features;
			if ( dw->image->data != NULL ) free(dw->image->data);
		}

		dw->image = image;

	}

	data = render_bin(dw->image->data, hdfile_get_width(dw->hdfile),
	                  hdfile_get_height(dw->hdfile), binning, max);

	return data;
}


#define RENDER_RGB							       \
									       \
	int s;								       \
	float p;							       \
									       \
	s = val / (max/6);						       \
	p = fmod(val, max/6);						       \
	p /= (max/6);							       \
									       \
	r = 0;	g = 0;	b = 0;						       \
									       \
	if ( (val < 0.0) || (val > max) ) {                                    \
		s = 0;                                                         \
		p = 1.0;                                                       \
	}                                                                      \
	switch ( s ) {							       \
		case 0 : {	/* Black to blue */			       \
			r = 0;		g = 0;			b = p*255;     \
			break;						       \
		}							       \
		case 1 : {	/* Blue to green */			       \
			r = 0;		g = 255*p;		b = (1-p)*255; \
			break;						       \
		}							       \
		case 2 : {	/* Green to red */			       \
			r =p*255;	g = (1-p)*255;		b = 0;	       \
			break;						       \
		}							       \
		case 3 : {	/* Red to Orange */			       \
			r = 255;	g = 127*p;		b = 0;	       \
			break;						       \
		}							       \
		case 4 : {	/* Orange to Yellow */			       \
			r = 255;	g = 127 + 127*p;	b = 0;	       \
			break;						       \
		}							       \
		case 5 : {	/* Yellow to White */			       \
			r = 255;	g = 255;		b = 255*p;     \
			break;						       \
		}							       \
		case 6 : {	/* Pixel has hit the maximum value */	       \
			r = 255;	g = 255;		b = 255;       \
			break;						       \
		}							       \
		default : {	/* Above saturation */			       \
			r = 255;	g = 255;		b = 255;       \
			break;						       \
		}							       \
	}

#define RENDER_MONO							       \
	float p;							       \
	p = (float)val / (float)max;					       \
	if ( val < 0.0 ) p = 0.0;                                              \
	if ( val > max ) p = 0.0;                                              \
	r = 255.0*p;	g = 255.0*p;	b = 255.0*p;


/* NB This function is shared between render_get_image() and
 * render_get_colour_scale() */
static void render_free_data(guchar *data, gpointer p)
{
	free(data);
}


static void show_marked_features(struct image *image, guchar *data,
                                 int w, int h, int binning)
{
	int i;

	if ( image->features == NULL ) return;

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		int x, y;
		double th;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		x = f->x;  y = f->y;

		x /= binning;
		y /= binning;

		for ( th=0; th<2*M_PI; th+=M_PI/40.0 ) {

			int nx, ny;

			nx = x + 10.0*cos(th);
			ny = y + 10.0*sin(th);

			if ( nx < 0 ) continue;
			if ( ny < 0 ) continue;
			if ( nx >= w ) continue;
			if ( ny >= h ) continue;

			data[3*( nx+w*(h-1-ny) )+0] = 255;

		}
	}
}


/* Return a pixbuf containing a rendered version of the image after binning.
 * This pixbuf might be scaled later - hopefully mostly in a downward
 * direction. */
GdkPixbuf *render_get_image(DisplayWindow *dw)
{
	int mw, mh, w, h;
	guchar *data;
	int16_t *hdr;
	size_t x, y;
	int16_t max;

	mw = hdfile_get_width(dw->hdfile);
	mh = hdfile_get_height(dw->hdfile);
	w = mw / dw->binning;
	h = mh / dw->binning;

	/* High dynamic range version */
	hdr = render_get_image_binned(dw, dw->binning, &max);
	if ( hdr == NULL ) return NULL;

	/* Rendered (colourful) version */
	data = malloc(3*w*h);
	if ( data == NULL ) {
		free(hdr);
		return NULL;
	}

	max /= dw->boostint;
	if ( max <= 6 ) { max = 10; }
	/* These x,y coordinates are measured relative to the bottom-left
	 * corner */
	for ( y=0; y<h; y++ ) {
	for ( x=0; x<w; x++ ) {

		int val;
		guchar r, g, b;

		val = hdr[x+w*y];
		if ( !dw->monochrome ) {
			RENDER_RGB
		} else {
			RENDER_MONO
		}

		/* Stuff inside square brackets makes this pixel go to
		 * the expected location in the pixbuf (which measures
		 * from the top-left corner */
		data[3*( x+w*(h-1-y) )+0] = r;
		data[3*( x+w*(h-1-y) )+1] = g;
		data[3*( x+w*(h-1-y) )+2] = b;

	}
	}

	show_marked_features(dw->image, data, w, h, dw->binning);

	/* Finished with this */
	free(hdr);

	/* Create the pixbuf from the 8-bit display data */
	return gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB, FALSE, 8,
					w, h, w*3, render_free_data, NULL);
}

GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int monochrome)
{
	guchar *data;
	size_t x, y;
	int max;

	data = malloc(3*w*h);
	if ( data == NULL ) return NULL;

	max = h;

	for ( y=0; y<h; y++ ) {

		guchar r, g, b;
		int val;

		val = y;
		if ( !monochrome ) {
			RENDER_RGB
		} else {
			RENDER_MONO
		}

		data[3*( 0+w*(h-1-y) )+0] = 0;
		data[3*( 0+w*(h-1-y) )+1] = 0;
		data[3*( 0+w*(h-1-y) )+2] = 0;
		for ( x=1; x<w; x++ ) {
			data[3*( x+w*(h-1-y) )+0] = r;
			data[3*( x+w*(h-1-y) )+1] = g;
			data[3*( x+w*(h-1-y) )+2] = b;
		}

	}

	return gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB, FALSE, 8,
					w, h, w*3, render_free_data, NULL);
}
