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
/* GSL is only used if Cairo is present */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
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
"  -d, --down=<h>,<k>,<l>  Indices for the axis in the downward direction.\n"
"  -r, --right=<h>,<k>,<l> Indices for the axis in the 'right' (roughly)\n"
"                           direction.\n"
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
static int get_basis_change_coefficients(double *in, double *out)
{
	int s;
	gsl_matrix *m;
	gsl_matrix *inv;
	gsl_permutation *perm;

	m = gsl_matrix_alloc(3, 3);
	if ( m == NULL ) {
		ERROR("Couldn't allocate memory for matrix\n");
		return 1;
	}
	gsl_matrix_set(m, 0, 0, in[0]);
	gsl_matrix_set(m, 0, 1, in[1]);
	gsl_matrix_set(m, 0, 2, in[2]);
	gsl_matrix_set(m, 1, 0, in[3]);
	gsl_matrix_set(m, 1, 1, in[4]);
	gsl_matrix_set(m, 1, 2, in[5]);
	gsl_matrix_set(m, 2, 0, in[6]);
	gsl_matrix_set(m, 2, 1, in[7]);
	gsl_matrix_set(m, 2, 2, in[8]);

	/* Invert */
	perm = gsl_permutation_alloc(m->size1);
	if ( perm == NULL ) {
		ERROR("Couldn't allocate permutation\n");
		gsl_matrix_free(m);
		return 1;
	}
	inv = gsl_matrix_alloc(3, 3);
	if ( inv == NULL ) {
		ERROR("Couldn't allocate inverse\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	if ( gsl_linalg_LU_decomp(m, perm, &s) ) {
		ERROR("Couldn't decompose matrix\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	if ( gsl_linalg_LU_invert(m, perm, inv)  ) {
		ERROR("Couldn't invert matrix\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	gsl_permutation_free(perm);
	gsl_matrix_free(m);

	/* Transpose */
	gsl_matrix_transpose(inv);

	out[0] = gsl_matrix_get(inv, 0, 0);
	out[1] = gsl_matrix_get(inv, 0, 1);
	out[2] = gsl_matrix_get(inv, 0, 2);
	out[3] = gsl_matrix_get(inv, 1, 0);
	out[4] = gsl_matrix_get(inv, 1, 1);
	out[5] = gsl_matrix_get(inv, 1, 2);
	out[6] = gsl_matrix_get(inv, 2, 0);
	out[7] = gsl_matrix_get(inv, 2, 1);
	out[8] = gsl_matrix_get(inv, 2, 2);

	gsl_matrix_free(inv);

	return 0;

}


static void draw_circles(signed int xh, signed int xk, signed int xl,
                         signed int yh, signed int yk, signed int yl,
                         signed int zh, signed int zk, signed int zl,
                         double *ref, unsigned int *counts, ReflItemList *items,
                         const char *sym,
                         cairo_t *dctx, int wght, double boost, int colscale,
                         UnitCell *cell, double radius, double theta,
                         double as, double bs, double cx, double cy,
                         double scale,
                         signed int *max_ux, signed int *max_uy,
                         double *max_val, double *max_u, double *max_v,
                         double *max_res)
{
	signed int xi, yi;
	double in[9];
	double bc[9];

	if ( dctx == NULL ) {
		*max_u = 0.0;  *max_v = 0.0;  *max_val = 0.0;
		*max_res = 0.0;  *max_ux = 0;  *max_uy = 0;
	}

	in[0] = xh;  in[1] = xk;  in[2] = xl;
	in[3] = yh;  in[4] = yk;  in[5] = yl;
	in[6] = zh;  in[7] = zk;  in[8] = zl;
	if ( get_basis_change_coefficients(in, bc) ) {
		ERROR("Couldn't change basis.\n");
		return;
	}

	/* Loop across the two basis directions */
	for ( xi=-INDMAX; xi<INDMAX; xi++ ) {
	for ( yi=-INDMAX; yi<INDMAX; yi++ ) {

		double u, v, val, res;
		signed int h, k, l;
		signed int he, ke, le;

		h = xi*xh + yi*yh;
		k = xi*xk + yi*yk;
		l = xi*xl + yi*yl;

		/* Got this reflection? */
		if ( find_unique_equiv(items, h, k, l, sym,
		                       &he, &ke, &le) == 0 ) continue;

		switch ( wght ) {
		case WGHT_I :
			val = lookup_intensity(ref, he, ke, le);
			break;
		case WGHT_SQRTI :
			val = lookup_intensity(ref, he, ke, le);
			val = (val>0.0) ? sqrt(val) : 0.0;
			break;
		case WGHT_COUNTS :
			val = lookup_count(counts, he, ke, le);
			val /= (float)num_equivs(he, ke, le, sym);
			break;
		case WGHT_RAWCOUNTS :
			val = lookup_count(counts, he, ke, le);
			break;
		default :
			ERROR("Invalid weighting.\n");
			abort();
		}

		/* Absolute location in image based on 2D basis */
		u = (double)xi*as*sin(theta);
		v = (double)xi*as*cos(theta) + (double)yi*bs;

		if ( dctx != NULL ) {

			float r, g, b;

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
			if ( fabs(val) > fabs(*max_val) ) {
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


static void render_za(UnitCell *cell, ReflItemList *items,
                      double *ref, unsigned int *counts,
                      double boost, const char *sym, int wght, int colscale,
                      signed int xh, signed int xk, signed int xl,
                      signed int yh, signed int yk, signed int yl)
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
	             ref, counts, items, sym, NULL, wght, boost, colscale, cell,
	             0.0, theta, as, bs, cx, cy, scale,
	             &max_ux, &max_uy, &max_val, &max_u, &max_v, &max_res);

	max_res /= 1e9;
	printf("Maximum resolution is 1/d = %5.3f nm^-1, d = %5.3f nm\n",
	       max_res*2.0, 1.0/(max_res*2.0));

	if ( max_val <= 0.0 ) {
		max_r = 4.0;
		STATUS("Couldn't find max value.\n");
		goto out;
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

	/* Test size of text that goes to the right(ish) */
	cairo_set_font_size(dctx, 40.0);
	snprintf(tmp, 255, "%i%i%i", abs(xh), abs(xk), abs(xl));
	cairo_text_extents(dctx, tmp, &size);

	cx = 532.0 - size.width;
	cy = 512.0 - 20.0;

	draw_circles(xh, xk, xl, yh, yk, yl, zh, zk, zl,
	             ref, counts, items, sym, dctx, wght, boost, colscale, cell,
	             max_r, theta, as, bs, cx, cy, scale,
	             NULL, NULL, &max_val, NULL, NULL, NULL);

out:
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
	signed int dh=1, dk=0, dl=0;
	signed int rh=0, rk=1, rl=0;
	char *down = NULL;
	char *right = NULL;

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
		{"down",               1, NULL,               'd'},
		{"right",              1, NULL,               'r'},
		{"counts",             0, &config_sqrt,        1},
		{"colour-key",         0, &config_colkey,      1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:p:w:c:y:d:r:",
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

	if ( config_zoneaxis ) {
		if ( (( down == NULL ) && ( right != NULL ))
		  || (( down != NULL ) && ( right == NULL )) ) {
			ERROR("Either specify both 'down' and 'right', or neither.\n");
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
	ref = new_list_intensity();
	cts = new_list_count();
	ReflItemList *items = read_reflections(infile, ref, NULL, cts);
	if ( ref == NULL ) {
		ERROR("Couldn't open file '%s'\n", infile);
		return 1;
	}

	if ( config_povray ) {
		r = povray_render_animation(cell, ref, cts, nproc);
	} else if ( config_zoneaxis ) {
#ifdef HAVE_CAIRO
		render_za(cell, items, ref, cts, boost, sym, wght, colscale,
		          rh, rk, rl, dh, dk, dl);
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
	delete_items(items);

	return r;
}
