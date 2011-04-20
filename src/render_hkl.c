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
#include "povray.h"
#include "symmetry.h"
#include "render.h"
#include "render_hkl.h"
#include "reflist.h"
#include "reflist-utils.h"


#define KEY_FILENAME "key.pdf"


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
"  -d, --down=<h>,<k>,<l>  Indices for the axis in the downward direction.\n"
"                           Default: 1,0,0.\n"
"  -r, --right=<h>,<k>,<l> Indices for the axis in the 'right' (roughly)\n"
"                           direction.  Default: 0,1,0.\n"
"  -o, --output=<filename> Output filename (not for POV-ray).  Default: za.pdf\n"
"      --boost=<val>       Squash colour scale by <val>.\n"
"  -p, --pdb=<file>        PDB file from which to get the unit cell.\n"
"  -y, --symmetry=<sym>    Expand reflections according to point group <sym>.\n"
"\n"
"  -c, --colscale=<scale>  Use the given colour scale.  Choose from:\n"
"                           mono    : Greyscale, black is zero.\n"
"                           invmono : Greyscale, white is zero.\n"
"                           colour  : Colour scale:\n"
"                                     black-blue-pink-red-orange-yellow-white\n"
"\n"
"  -w  --weighting=<wght>  Colour/shade the reciprocal lattice points\n"
"                           according to:\n"
"                            I      : the intensity of the reflection.\n"
"                            sqrtI  : the square root of the intensity.\n"
"                            count  : the number of measurements for the reflection.\n"
"                                     (after correcting for 'epsilon')\n"
"                            rawcts : the raw number of measurements for the\n"
"                                     reflection (no 'epsilon' correction).\n"
"\n"
"      --colour-key        Draw (only) the key for the current colour scale.\n"
"  -j <n>                  Run <n> instances of POV-ray in parallel.\n"
"  -h, --help              Display this help message.\n"
);
}


#ifdef HAVE_CAIRO


static void draw_circles(signed int xh, signed int xk, signed int xl,
                         signed int yh, signed int yk, signed int yl,
                         signed int zh, signed int zk, signed int zl,
                         RefList *list, const char *sym,
                         cairo_t *dctx, int wght, double boost, int colscale,
                         UnitCell *cell, double radius, double theta,
                         double as, double bs, double cx, double cy,
                         double scale,
                         signed int *max_ux, signed int *max_uy,
                         double *max_val, double *max_u, double *max_v,
                         double *max_res)
{
	signed int xi, yi;

	if ( dctx == NULL ) {
		*max_u = 0.0;  *max_v = 0.0;  *max_val = 0.0;
		*max_res = 0.0;  *max_ux = 0;  *max_uy = 0;
	}

	/* Iterate over possible reflections in this zone */
	for ( xi=-INDMAX; xi<INDMAX; xi++ ) {
	for ( yi=-INDMAX; yi<INDMAX; yi++ ) {

		double u, v, val, res;
		signed int h, k, l;
		signed int he, ke, le;
		Reflection *refl;
		int r;

		h = xi*xh + yi*yh;
		k = xi*xk + yi*yk;
		l = xi*xl + yi*yl;

		/* Got this reflection? */
		r = find_equiv_in_list(list, h, k, l, sym, &he, &ke, &le);
		if ( !r ) continue;
		refl = find_refl(list, he, ke, le);

		switch ( wght ) {
		case WGHT_I :
			val = get_intensity(refl);
			break;
		case WGHT_SQRTI :
			val = get_intensity(refl);
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = get_redundancy(refl);
			val /= (float)num_equivs(h, k, l, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = get_redundancy(refl);
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		/* Absolute location in image based on 2D basis */
		u = (double)xi*as*sin(theta);
		v = (double)xi*as*cos(theta) + (double)yi*bs;

		if ( dctx != NULL ) {

			double r, g, b;

			cairo_arc(dctx, ((double)cx)+u*scale,
			                ((double)cy)+v*scale,
			                radius, 0, 2*M_PI);

			render_scale(val, *max_val/boost, colscale,
			             &r, &g, &b);
			cairo_set_source_rgb(dctx, r, g, b);
			cairo_fill(dctx);

		} else {

			/* Find max vectors in plane for scaling */
			if ( fabs(u) > fabs(*max_u) ) *max_u = fabs(u);
			if ( fabs(v) > fabs(*max_v) ) *max_v = fabs(v);

			/* Find max value for colour scale */
			if ( !isnan(val) && !isinf(val)
			  && (fabs(val) > fabs(*max_val)) )
			{
				*max_val = fabs(val);
			}

			/* Find max indices */
			if ( (yi==0) && (fabs(xi) > *max_ux) )
				*max_ux = fabs(xi);
			if ( (xi==0) && (fabs(yi) > *max_uy) )
				*max_uy = fabs(yi);

			/* Find max resolution */
			res = resolution(cell, h, k, l);
			if ( res > *max_res ) *max_res = res;

		}

	}
	}
}


static void render_overlined_indices(cairo_t *dctx,
                                     signed int h, signed int k, signed int l)
{
	char tmp[256];
	cairo_text_extents_t size;
	double x, y;
	const double sh = 39.0;

	cairo_get_current_point(dctx, &x, &y);
	cairo_set_line_width(dctx, 4.0);

	/* Draw 'h' */
	snprintf(tmp, 255, "%i", abs(h));
	cairo_text_extents(dctx, tmp, &size);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( h < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-sh);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}
	x += size.x_advance;

	/* Draw 'k' */
	cairo_move_to(dctx, x, y);
	snprintf(tmp, 255, "%i", abs(k));
	cairo_text_extents(dctx, tmp, &size);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( k < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-sh);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}
	x += size.x_advance;

	/* Draw 'l' */
	cairo_move_to(dctx, x, y);
	snprintf(tmp, 255, "%i", abs(l));
	cairo_text_extents(dctx, tmp, &size);
	cairo_show_text(dctx, tmp);
	cairo_fill(dctx);
	if ( l < 0 ) {
		cairo_move_to(dctx, x+size.x_bearing, y-sh);
		cairo_rel_line_to(dctx, size.width, 0.0);
		cairo_stroke(dctx);
	}
}


static void render_za(UnitCell *cell, RefList *list,
                      double boost, const char *sym, int wght, int colscale,
                      signed int xh, signed int xk, signed int xl,
                      signed int yh, signed int yk, signed int yl,
                      const char *outfile, double scale_top)
{
	cairo_surface_t *surface;
	cairo_t *dctx;
	double max_u, max_v, max_res, max_val;
	double scale_u, scale_v, scale;
	double sep_u, sep_v, max_r;
	double u, v;
	signed int max_ux, max_uy;
	double as, bs, theta;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	float wh, ht;
	signed int zh, zk, zl;
	double xx, xy, xz;
	double yx, yy, yz;
	char tmp[256];
	cairo_text_extents_t size;
	double cx, cy;
	const double border = 200.0;
	int png;

	/* Vector product to determine the zone axis. */
	zh = xk*yl - xl*yk;
	zk = - xh*yl + xl*yh;
	zl = xh*yk - xk*yh;
	STATUS("Zone axis is %i %i %i\n", zh, zk, zl);

	/* Size of output and centre definition */
	wh = 1024;
	ht = 1024;

	/* Work out reciprocal lattice spacings and angles for this cut */
	if ( cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz) ) {
		ERROR("Couldn't get reciprocal parameters\n");
		return;
	}
	xx = xh*asx + xk*bsx + xl*csx;
	xy = xh*asy + xk*bsy + xl*csy;
	xz = xh*asz + xk*bsz + xl*csz;
	yx = yh*asx + yk*bsx + yl*csx;
	yy = yh*asy + yk*bsy + yl*csy;
	yz = yh*asz + yk*bsz + yl*csz;
	theta = angle_between(xx, xy, xz, yx, yy, yz);
	as = modulus(xx, xy, xz) / 1e9;
	bs = modulus(yx, yy, yz) / 1e9;

	scale = 1.0;
	draw_circles(xh, xk, xl, yh, yk, yl, zh, zk, zl,
	             list, sym, NULL, wght, boost, colscale, cell,
	             0.0, theta, as, bs, 0.0, 0.0, scale,
	             &max_ux, &max_uy, &max_val, &max_u, &max_v, &max_res);

	max_res /= 1e9;
	printf("Maximum resolution is 1/d = %5.3f nm^-1, d = %5.3f nm\n",
	       max_res*2.0, 1.0/(max_res*2.0));

	if ( max_val <= 0.0 ) {
		STATUS("Couldn't find max value.\n");
		return;
	}

	/* Use manual scale top if specified */
	if ( scale_top > 0.0 ) {
		max_val = scale_top;
	}

	/* Choose whichever scaling factor gives the smallest value */
	scale_u = ((double)wh-border) / (2.0*max_u);
	scale_v = ((double)ht-border) / (2.0*max_v);
	scale = (scale_u < scale_v) ? scale_u : scale_v;

	sep_u = scale*as;
	sep_v = scale*bs;
	/* We are interested in the smaller of the two separations */
	max_r = (sep_u < sep_v) ? sep_u : sep_v;
	max_r /= 2.0;  /* Max radius is half the separation */
	max_r -= 1.0;  /* Add a tiny separation between circles */
	if ( max_r < 1.0 ) {
		ERROR("Circle radius is probably too small (%f).\n", max_r);
	}

	if ( outfile == NULL ) outfile = "za.pdf";

	if ( strcmp(outfile+strlen(outfile)-4, ".png") == 0 ) {
		png = 1;
		surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
		                                     wh, ht);
	} else {
		png = 0;
		surface = cairo_pdf_surface_create(outfile, wh, ht);
	}

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		ERROR("Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return;
	}

	dctx = cairo_create(surface);
	if ( cairo_status(dctx) != CAIRO_STATUS_SUCCESS ) {
		ERROR("Couldn't create Cairo context\n");
		cairo_surface_destroy(surface);
		return;
	}

	/* Black background */
	cairo_rectangle(dctx, 0.0, 0.0, wh, ht);
	cairo_set_source_rgb(dctx, 0.0, 0.0, 0.0);
	cairo_fill(dctx);

	/* Test size of text that goes to the right(ish) */
	cairo_set_font_size(dctx, 40.0);
	snprintf(tmp, 255, "%i%i%i", abs(xh), abs(xk), abs(xl));
	cairo_text_extents(dctx, tmp, &size);

	cx = 532.0 - size.width;
	cy = 512.0 - 20.0;

	draw_circles(xh, xk, xl, yh, yk, yl, zh, zk, zl,
	             list, sym, dctx, wght, boost, colscale, cell,
	             max_r, theta, as, bs, cx, cy, scale,
	             NULL, NULL, &max_val, NULL, NULL, NULL);

	/* Centre marker */
	cairo_arc(dctx, (double)cx,
			(double)cy, max_r, 0, 2*M_PI);
	cairo_set_source_rgb(dctx, 1.0, 0.0, 0.0);
	cairo_fill(dctx);

	/* Draw indexing lines */
	cairo_set_line_cap(dctx, CAIRO_LINE_CAP_ROUND);
	cairo_set_line_width(dctx, 4.0);
	cairo_move_to(dctx, (double)cx, (double)cy);
	u = (2.0+max_ux)*as*sin(theta);
	v = (2.0+max_ux)*as*cos(theta);
	cairo_line_to(dctx, cx+u*scale, cy+v*scale);
	cairo_set_source_rgb(dctx, 0.0, 1.0, 0.0);
	cairo_stroke(dctx);

	cairo_set_font_size(dctx, 40.0);
	snprintf(tmp, 255, "%i%i%i", abs(xh), abs(xk), abs(xl));
	cairo_text_extents(dctx, tmp, &size);

	cairo_move_to(dctx, cx+u*scale + 20.0, cy+v*scale + size.height/2.0);
	render_overlined_indices(dctx, xh, xk, xl);
	cairo_fill(dctx);

	snprintf(tmp, 255, "%i%i%i", abs(yh), abs(yk), abs(yl));
	cairo_text_extents(dctx, tmp, &size);

	cairo_move_to(dctx, (double)cx, (double)cy);
	u = 0.0;
	v = (2.0+max_uy)*bs;
	cairo_line_to(dctx, cx+u*scale, cy+v*scale);
	cairo_set_source_rgb(dctx, 0.0, 1.0, 0.0);
	cairo_stroke(dctx);

	cairo_move_to(dctx, cx+u*scale - size.width/2.0,
	                    cy+v*scale + size.height + 20.0);
	render_overlined_indices(dctx, yh, yk, yl);
	cairo_fill(dctx);

	if ( png ) {
		int r = cairo_surface_write_to_png(surface, outfile);
		if ( r != CAIRO_STATUS_SUCCESS ) {
			ERROR("Failed to write PNG to '%s'\n", outfile);
		}
	}

	cairo_surface_finish(surface);
	cairo_destroy(dctx);
}


static int render_key(int colscale, double scale_top)
{
	cairo_surface_t *surface;
	cairo_t *dctx;
	double top, wh, ht, y;
	double slice;

	wh = 128.0;
	ht = 1024.0;
	slice = 1.0;

	if ( scale_top > 0.0 ) {
		top = scale_top;
	} else {
		top = 1.0;
	}

	surface = cairo_pdf_surface_create(KEY_FILENAME, wh, ht);

	if ( cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS ) {
		fprintf(stderr, "Couldn't create Cairo surface\n");
		cairo_surface_destroy(surface);
		return 1;
	}

	dctx = cairo_create(surface);

	for ( y=0.0; y<ht; y+=slice ) {

		double r, g, b;
		double val;
		double v = y;

		cairo_rectangle(dctx, 0.0, ht-y, wh/2.0, slice);

		if ( colscale == SCALE_RATIO ) {
			if ( v < ht/2.0 ) {
				val = v/(ht/2.0);
			} else {
				val = (((v-ht/2.0)/(ht/2.0))*(top-1.0))+1.0;
			}
		} else {
			val = v/ht;
		}

		render_scale(val, top, colscale, &r, &g, &b);
		cairo_set_source_rgb(dctx, r, g, b);

		cairo_stroke_preserve(dctx);
		cairo_fill(dctx);

	}

	if ( colscale == SCALE_RATIO ) {

		cairo_text_extents_t size;
		char tmp[32];

		cairo_rectangle(dctx, 0.0, ht/2.0-2.0, wh/2.0, 4.0);
		cairo_set_source_rgb(dctx, 0.0, 0.0, 0.0);
		cairo_stroke_preserve(dctx);
		cairo_fill(dctx);

		cairo_set_font_size(dctx, 20.0);
		cairo_text_extents(dctx, "1.0", &size);
		cairo_move_to(dctx, wh/2.0+5.0, ht/2.0+size.height/2.0);
		cairo_show_text(dctx, "1.0");

		cairo_set_font_size(dctx, 20.0);
		cairo_text_extents(dctx, "0.0", &size);
		cairo_move_to(dctx, wh/2.0+5.0, ht-5.0);
		cairo_show_text(dctx, "0.0");

		cairo_set_font_size(dctx, 20.0);
		snprintf(tmp, 31, "%.1f", top);
		cairo_text_extents(dctx, tmp, &size);
		cairo_move_to(dctx, wh/2.0+5.0, size.height+5.0);
		cairo_show_text(dctx, tmp);

	}


	cairo_surface_finish(surface);
	cairo_destroy(dctx);

	STATUS("Colour key written to "KEY_FILENAME"\n");

	return 0;
}


#else  /* HAVE_CAIRO */


static int render_key(int colscale, double scale_top)
{
	ERROR("This version of CrystFEL was compiled without Cairo");
	ERROR(" support, which is required to draw the colour");
	ERROR(" scale.  Sorry!\n");
	return 1;
}


static void render_za(UnitCell *cell, RefList *list,
                      double boost, const char *sym, int wght, int colscale,
                      signed int xh, signed int xk, signed int xl,
                      signed int yh, signed int yk, signed int yl,
                      const char *outfile, double scale_top)
{
	ERROR("This version of CrystFEL was compiled without Cairo");
	ERROR(" support, which is required to plot a zone axis");
	ERROR(" pattern.  Sorry!\n");
}


#endif /* HAVE_CAIRO */


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	RefList *list;
	char *infile;
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
	signed int dh=1, dk=0, dl=0;
	signed int rh=0, rk=1, rl=0;
	char *down = NULL;
	char *right = NULL;
	char *outfile = NULL;
	double scale_top = -1.0;
	char *endptr;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"povray",             0, &config_povray,      1},
		{"zone-axis",          0, &config_zoneaxis,    1},
		{"output",             1, NULL,               'o'},
		{"pdb",                1, NULL,               'p'},
		{"boost",              1, NULL,               'b'},
		{"symmetry",           1, NULL,               'y'},
		{"weighting",          1, NULL,               'w'},
		{"colscale",           1, NULL,               'c'},
		{"down",               1, NULL,               'd'},
		{"right",              1, NULL,               'r'},
		{"counts",             0, &config_sqrt,        1},
		{"colour-key",         0, &config_colkey,      1},
		{"scale-top",          1, NULL,                2},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:p:w:c:y:d:r:o:",
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

		case 'd' :
			down = strdup(optarg);
			break;

		case 'r' :
			right = strdup(optarg);
			break;

		case 'o' :
			outfile = strdup(optarg);
			break;

		case 2 :
			errno = 0;
			scale_top = strtod(optarg, &endptr);
			if ( !( (optarg[0] != '\0') && (endptr[0] == '\0') )
			   || (errno != 0) )
			{
				ERROR("Invalid scale top('%s')\n", optarg);
				return 1;
			}
			if ( scale_top < 0.0 ) {
				ERROR("Scale top must be positive.\n");
				return 1;
			}
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( (pdb == NULL) && !config_colkey ) {
		ERROR("You must specify the PDB containing the unit cell.\n");
		return 1;
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
	} else if ( strcmp(cscale, "ratio") == 0 ) {
		colscale = SCALE_RATIO;
	} else {
		ERROR("Unrecognised colour scale '%s'\n", cscale);
		return 1;
	}
	free(cscale);

	if ( config_colkey ) {
		return render_key(colscale, scale_top);
	}

	if ( config_zoneaxis ) {
		if ( (( down == NULL ) && ( right != NULL ))
		  || (( down != NULL ) && ( right == NULL )) ) {
			ERROR("Either specify both 'down' and 'right',"
			      " or neither.\n");
			return 1;
		}
		if ( down != NULL ) {
			int r;
			r = sscanf(down, "%i,%i,%i", &dh, &dk, &dl);
			if ( r != 3 ) {
				ERROR("Invalid format for 'down'\n");
				return 1;
			}
		}
		if ( right != NULL ) {
			int r;
			r = sscanf(right, "%i,%i,%i", &rh, &rk, &rl);
			if ( r != 3 ) {
				ERROR("Invalid format for 'right'\n");
				return 1;
			}
		}
	}

	infile = argv[optind];

	cell = load_cell_from_pdb(pdb);
	if ( cell == NULL ) {
		ERROR("Couldn't load unit cell from %s\n", pdb);
		return 1;
	}
	list = read_reflections(infile);
	if ( list == NULL ) {
		ERROR("Couldn't read file '%s'\n", infile);
		return 1;
	}
	if ( check_list_symmetry(list, sym) ) {
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", sym);
		return 1;
	}

	if ( config_povray ) {
		r = povray_render_animation(cell, list,
		                            nproc, sym, wght, boost, scale_top);
	} else if ( config_zoneaxis ) {
		render_za(cell, list, boost, sym, wght, colscale,
		          rh, rk, rl, dh, dk, dl, outfile, scale_top);
	} else {
		ERROR("Try again with either --povray or --zone-axis.\n");
	}

	free(pdb);
	free(sym);
	reflist_free(list);
	if ( outfile != NULL ) free(outfile);

	return r;
}
