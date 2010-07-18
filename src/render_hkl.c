/*
 * render_hkl.c
 *
 * Draw pretty renderings of reflection lists
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#ifdef HAVE_CAIRO
#include <cairo.h>
#include <cairo-pdf.h>
#endif

#include "utils.h"
#include "reflections.h"
#include "povray.h"
#include "symmetry.h"
#include "render.h"

enum {
	WGHT_I,
	WGHT_SQRTI,
	WGHT_COUNTS,
	WGHT_RAWCOUNTS,
};

static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Render intensity lists in various ways.\n"
"\n"
"      --povray            Render a 3D animation using POV-ray.\n"
#ifdef HAVE_CAIRO
"      --zone-axis         Render a 2D zone axis pattern.\n"
#endif
"\n"
"      --boost=<val>       Squash colour scale by <val>.\n"
"  -p, --pdb=<file>        PDB file from which to get the unit cell.\n"
"  -y, --symmetry=<sym>    Expand reflections according to point group <sym>.\n"
"\n"
"  -c, --colscale=<scale>  Use the given colour scale.  Choose from:\n"
"                           mono    : Greyscale, black is zero.\n"
"                           invmono : Greyscale, white is zero.\n"
"                           colour  : Colours scale:\n"
"                                     black-blue-pink-red-orange-yellow-white\n"
"\n"
"  -w  --weighting=<wght>  Colour/shade the reciprocal lattice points\n"
"                           according to:\n"
"                            I      : the intensity of the reflection.\n"
"                            sqrtI  : the square root of the intensity.\n"
"                            count  : the number of hits for the reflection.\n"
"                                     (after correcting for 'epsilon')\n"
"                            rawcts : the raw number of hits for the\n"
"                                     reflection (no 'epsilon' correction).\n"
"\n"
"      --colour-key        Draw (only) the key for the current colour scale.\n"
"  -j <n>                  Run <n> instances of POV-ray in parallel.\n"
"  -h, --help              Display this help message.\n"
);
}


#ifdef HAVE_CAIRO
static void render_za(UnitCell *cell, double *ref, unsigned int *c,
                      double boost, const char *sym, int wght, int colscale)
{
	cairo_surface_t *surface;
	cairo_t *dctx;
	double max_u, max_v, max_res, max_val;
	double scale_u, scale_v, scale;
	double sep_u, sep_v, max_r;
	double u, v;
	signed int max_h, max_k;
	double as, bs, theta;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	signed int h, k;
	float wh, ht;

	wh = 1024;
	ht = 1024;

	surface = cairo_pdf_surface_create("za.pdf", wh, ht);

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return;
	}

	dctx = cairo_create(surface);

	/* Black background */
	cairo_rectangle(dctx, 0.0, 0.0, wh, ht);
	cairo_set_source_rgb(dctx, 0.0, 0.0, 0.0);
	cairo_fill(dctx);

	max_u = 0.0;  max_v = 0.0;  max_val = 0.0;
	max_res = 0.0;  max_h = 0;  max_k = 0;

	/* Work out reciprocal lattice spacings and angles for this cut */
	if ( cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz) ) {
		ERROR("Couldn't get reciprocal parameters\n");
		return;
	}
	theta = angle_between(asx, asy, asz, bsx, bsy, bsz);
	as = modulus(asx, asy, asz) / 1e9;
	bs = modulus(bsx, bsy, bsz) / 1e9;

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {

		double u, v, val, res;
		int ct;
		int nequiv, p;

		ct = lookup_count(c, h, k, 0);
		if ( ct < 1 ) continue;

		switch ( wght ) {
		case WGHT_I :
			val = lookup_intensity(ref, h, k, 0);
			break;
		case WGHT_SQRTI :
			val = lookup_intensity(ref, h, k, 0);
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = (float)ct;
			val /= (float)num_equivs(h, k, 0, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = (float)ct;
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		nequiv = num_equivs(h, k, 0, sym);
		for ( p=0; p<nequiv; p++ ) {

			signed int he, ke, le;
			get_equiv(h, k, 0, &he, &ke, &le, sym, p);

			u =  (double)he*as*sin(theta);
			v =  (double)he*as*cos(theta) + ke*bs;
			if ( fabs(u) > fabs(max_u) ) max_u = fabs(u);
			if ( fabs(v) > fabs(max_v) ) max_v = fabs(v);
			if ( fabs(val) > fabs(max_val) ) {
				max_val = fabs(val);
			}
			if ( fabs(he) > max_h ) max_h = fabs(he);
			if ( fabs(ke) > max_k ) max_k = fabs(ke);
			res = resolution(cell, he, ke, 0);
			if ( res > max_res ) max_res = res;

		}

	}
	}

	max_res /= 1e9;
	printf("Maximum resolution is 1/d = %5.3f nm^-1, d = %5.3f nm\n",
	       max_res*2.0, 1.0/(max_res*2.0));

	if ( max_val <= 0.0 ) {
		max_r = 4.0;
		goto out;
	}

	/* Choose whichever scaling factor gives the smallest value */
	scale_u = ((double)(wh/2.0)-50.0) / max_u;
	scale_v = ((double)(ht/2.0)-50.0) / max_v;
	scale = (scale_u < scale_v) ? scale_u : scale_v;

	sep_u = as * scale * cos(theta);
	sep_v = bs * scale;
	max_r = (sep_u < sep_v) ? sep_u : sep_v;
	max_r -= 1.0;  /* Add a tiny separation between circles */

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {

		double u, v, val;
		int ct;
		int nequiv, p;

		ct = lookup_count(c, h, k, 0);
		if ( ct < 1 ) continue;  /* Must have at least one count */

		switch ( wght ) {
		case WGHT_I :
			val = lookup_intensity(ref, h, k, 0) / (float)ct;
			break;
		case WGHT_SQRTI :
			val = lookup_intensity(ref, h, k, 0) / (float)ct;
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = (float)ct;
			val /= (float)num_equivs(h, k, 0, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = (float)ct;
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		nequiv = num_equivs(h, k, 0, sym);
		for ( p=0; p<nequiv; p++ ) {

			signed int he, ke, le;
			float r, g, b;
			get_equiv(h, k, 0, &he, &ke, &le, sym, p);

			u = (double)he*as*sin(theta);
			v = (double)he*as*cos(theta) + ke*bs;

			cairo_arc(dctx, ((double)wh/2)+u*scale,
					((double)ht/2)+v*scale, max_r,
					0, 2*M_PI);

			render_scale(val, max_val/boost, colscale, &r, &g, &b);
			cairo_set_source_rgb(dctx, r, g, b);
			cairo_fill(dctx);

		}

	}
	}

out:
	/* Centre marker */
	cairo_arc(dctx, (double)wh/2,
			(double)ht/2, max_r, 0, 2*M_PI);
	cairo_set_source_rgb(dctx, 1.0, 0.0, 0.0);
	cairo_fill(dctx);

	/* Draw indexing lines */
	cairo_set_line_width(dctx, 4.0);
	cairo_move_to(dctx, (double)wh/2, (double)ht/2);
	u = (double)max_h*as*sin(theta);
	v = (double)max_h*as*cos(theta) + 0*bs;
	cairo_line_to(dctx, ((double)wh/2)+u*scale,
	                    ((double)ht/2)+v*scale);
	cairo_set_source_rgb(dctx, 0.0, 1.0, 0.0);
	cairo_stroke(dctx);

	cairo_move_to(dctx,((double)wh/2)+u*scale-40.0,
	                   ((double)ht/2)+v*scale+40.0);
	                   cairo_set_font_size(dctx, 40.0);
	cairo_show_text(dctx, "h");
	cairo_fill(dctx);

	cairo_move_to(dctx, (double)wh/2, (double)ht/2);
	u = 0.0;
	v = (double)max_k*bs;
	cairo_line_to(dctx, ((double)wh/2)+u*scale,
	                    ((double)ht/2)+v*scale);
	cairo_set_source_rgb(dctx, 0.0, 1.0, 0.0);
	cairo_stroke(dctx);

	cairo_move_to(dctx,((double)wh/2)+u*scale-40.0,
	                   ((double)ht/2)+v*scale-40.0);
	                   cairo_set_font_size(dctx, 40.0);
	cairo_show_text(dctx, "k");
	cairo_fill(dctx);

	cairo_surface_finish(surface);
	cairo_destroy(dctx);
}

static int render_key(int colscale)
{
	cairo_surface_t *surface;
	cairo_t *dctx;
	float wh, ht;
	float y;

	wh = 128;
	ht = 1024;

	surface = cairo_pdf_surface_create("key.pdf", wh, ht);

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return 1;
	}

	dctx = cairo_create(surface);

	for ( y=0; y<ht; y++ ) {

		float r, g, b;

		cairo_rectangle(dctx, 0.0, y, wh, y+1.0);

		render_scale(ht-y, ht, colscale, &r, &g, &b);
		cairo_set_source_rgb(dctx, r, g, b);

		cairo_fill(dctx);

	}

	cairo_surface_finish(surface);
	cairo_destroy(dctx);

	return 0;
}
#endif


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *infile;
	double *ref;
	int config_povray = 0;
	int config_zoneaxis = 0;
	int config_sqrt = 0;
	int config_colkey = 0;
	unsigned int nproc = 1;
	char *pdb = NULL;
	int r = 0;
	double boost = 1.0;
	char *sym = NULL;
	char *weighting = NULL;
	int wght;
	int colscale;
	char *cscale = NULL;
	unsigned int *cts;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"povray",             0, &config_povray,      1},
		{"zone-axis",          0, &config_zoneaxis,    1},
		{"pdb",                1, NULL,               'p'},
		{"boost",              1, NULL,               'b'},
		{"symmetry",           1, NULL,               'y'},
		{"weighting",          1, NULL,               'w'},
		{"colscale",           1, NULL,               'c'},
		{"counts",             0, &config_sqrt,        1},
		{"colour-key",         0, &config_colkey,      1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:p:w:c:y:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'j' :
			nproc = atoi(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
			break;

		case 'b' :
			boost = atof(optarg);
			break;

		case 'y' :
			sym = strdup(optarg);
			break;

		case 'w' :
			weighting = strdup(optarg);
			break;

		case 'c' :
			cscale = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	if ( sym == NULL ) {
		sym = strdup("1");
	}

	if ( weighting == NULL ) {
		weighting = strdup("I");
	}

	if ( strcmp(weighting, "I") == 0 ) {
		wght = WGHT_I;
	} else if ( strcmp(weighting, "sqrtI") == 0 ) {
		wght = WGHT_SQRTI;
	} else if ( strcmp(weighting, "count") == 0 ) {
		wght = WGHT_COUNTS;
	} else if ( strcmp(weighting, "counts") == 0 ) {
		wght = WGHT_COUNTS;
	} else if ( strcmp(weighting, "rawcts") == 0 ) {
		wght = WGHT_RAWCOUNTS;
	} else if ( strcmp(weighting, "rawcount") == 0 ) {
		wght = WGHT_RAWCOUNTS;
	} else if ( strcmp(weighting, "rawcounts") == 0 ) {
		wght = WGHT_RAWCOUNTS;
	} else {
		ERROR("Unrecognised weighting '%s'\n", weighting);
		return 1;
	}
	free(weighting);

	if ( cscale == NULL ) {
		cscale = strdup("mono");
	}

	if ( strcmp(cscale, "mono") == 0 ) {
		colscale = SCALE_MONO;
	} else if ( strcmp(cscale, "invmono") == 0 ) {
		colscale = SCALE_INVMONO;
	} else if ( strcmp(cscale, "colour") == 0 ) {
		colscale = SCALE_COLOUR;
	} else if ( strcmp(cscale, "color") == 0 ) {
		colscale = SCALE_COLOUR;
	} else {
		ERROR("Unrecognised colour scale '%s'\n", cscale);
		return 1;
	}
	free(cscale);

	if ( config_colkey ) {
		return render_key(colscale);
	}

	infile = argv[optind];

	cell = load_cell_from_pdb(pdb);
	if ( cell == NULL ) {
		ERROR("Couldn't load unit cell from %s\n", pdb);
		return 1;
	}
	ref = new_list_intensity();
	cts = new_list_count();
	ReflItemList *items = read_reflections(infile, ref, NULL, cts);
	delete_items(items);
	if ( ref == NULL ) {
		ERROR("Couldn't open file '%s'\n", infile);
		return 1;
	}

	if ( config_povray ) {
		r = povray_render_animation(cell, ref, cts, nproc);
	} else if ( config_zoneaxis ) {
#ifdef HAVE_CAIRO
		render_za(cell, ref, cts, boost, sym, wght, colscale);
#else
		ERROR("This version of CrystFEL was compiled without Cairo");
		ERROR(" support, which is required to plot a zone axis");
		ERROR(" pattern.  Sorry!\n");
#endif
	} else {
		ERROR("Try again with either --povray or --zone-axis.\n");
	}

	free(pdb);
	free(sym);

	return r;
}
