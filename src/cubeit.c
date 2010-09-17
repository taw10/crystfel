/*
 * cubeit.c
 *
 * "Full integration" of diffraction data
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
#include <errno.h>
#include <assert.h>
#include <png.h>

#include "image.h"
#include "cell.h"
#include "hdf5-file.h"
#include "diffraction.h"
#include "render.h"


#define MAX_HITS (1024)


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"'Full integration' of diffraction data.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input stream.\n"
"                              Can be '-' for stdin.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
);
}


static UnitCell *read_orientation_matrix(FILE *fh)
{
	float u, v, w;
	struct rvec as, bs, cs;
	UnitCell *cell;
	char line[1024];

	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read a-star\n");
		return NULL;
	}
	as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read b-star\n");
		return NULL;
	}
	bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read c-star\n");
		return NULL;
	}
	cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
	cell = cell_new_from_axes(as, bs, cs);

	return cell;
}


static int find_chunk(FILE *fh, UnitCell **cell, char **filename)
{
	char line[1024];
	char *rval = NULL;

	do {

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		chomp(line);

		if ( strncmp(line, "Reflections from indexing", 25) != 0 ) {
			continue;
		}

		*filename = strdup(line+29);

		/* Skip two lines (while checking for errors) */
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		*cell = read_orientation_matrix(fh);
		if ( *cell == NULL ) {
			STATUS("Got filename but no cell for %s\n", *filename);
			continue;
		}

		return 0;

	} while ( rval != NULL );

	return 1;
}


static void interpolate_linear(double *vals, double v,
                               int xs, int ys, int zs,
                               int xv, int yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (zv+0.5)*(float)zs;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < zs);

	val1 = v*(1.0-f);
	val2 = v*f;

	vals[xs*ys*k + xs*yv + xv] = val1;
	vals[xs*ys*(k+1) + xs*yv + xv] = val2;
}


static void interpolate_bilinear(double *vals, double v,
                                 int xs, int ys, int zs,
                                 int xv, double yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (yv+0.5)*(float)ys;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < ys);

	val1 = v*(1.0-f);
	val2 = v*f;

	interpolate_linear(vals, val1, xs, ys, zs, xv, k, zv);
	interpolate_linear(vals, val2, xs, ys, zs, xv, k+1, zv);
}


static void interpolate_onto_grid(double *vals, double v,
                                  int xs, int ys, int zs,
                                  double xv, double yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (xv+0.5)*(float)xs;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < xs);

	val1 = v*(1.0-f);
	val2 = v*f;

	interpolate_bilinear(vals, val1, xs, ys, zs, k, yv, zv);
	interpolate_bilinear(vals, val2, xs, ys, zs, k+1, yv, zv);
}


static void write_slice(const char *filename, double *vals, int z,
                        int xs, int ys, int zs)
{
#ifdef HAVE_LIBPNG
	FILE *fh;
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep *row_pointers;
	int x, y;
	float max = 0.0;
	int w, h;

	w = xs;
	h = ys;

	if ( max <= 6 ) { max = 10; }

	for ( y=0; y<h; y++ ) {
	for ( x=0; x<w; x++ ) {

		float val;

		val = vals[xs*ys*z + xs*y + x];

		if ( val > max ) max = val;

	}
	}

	fh = fopen(filename, "wb");
	if ( !fh ) {
		ERROR("Couldn't open output file.\n");
		return;
	}
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
	                                  NULL, NULL, NULL);
	if ( !png_ptr ) {
		ERROR("Couldn't create PNG write structure.\n");
		fclose(fh);
		return;
	}
	info_ptr = png_create_info_struct(png_ptr);
	if ( !info_ptr ) {
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		ERROR("Couldn't create PNG info structure.\n");
		fclose(fh);
		return;
	}
	if ( setjmp(png_jmpbuf(png_ptr)) ) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fh);
		ERROR("PNG write failed.\n");
		return;
	}
	png_init_io(png_ptr, fh);

	png_set_IHDR(png_ptr, info_ptr, w, h, 8,
	             PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

	row_pointers = malloc(h*sizeof(png_bytep *));

	/* Write the image data */
	max /= 1.0;
	if ( max <= 6 ) { max = 10; }

	for ( y=0; y<h; y++ ) {

		row_pointers[y] = malloc(w*3);

		for ( x=0; x<w; x++ ) {

			float r, g, b;
			float val;

			val = vals[xs*ys*z + xs*y + x];

			render_scale(val, max, SCALE_COLOUR, &r, &g, &b);
			row_pointers[y][3*x] = (png_byte)255*r;
			row_pointers[y][3*x+1] = (png_byte)255*g;
			row_pointers[y][3*x+2] = (png_byte)255*b;

		}
	}

	for ( y=0; y<h/2+1; y++ ) {
		png_bytep scratch;
		scratch = row_pointers[y];
		row_pointers[y] = row_pointers[h-y-1];
		row_pointers[h-y-1] = scratch;
	}

	png_set_rows(png_ptr, info_ptr, row_pointers);
	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

	png_destroy_write_struct(&png_ptr, &info_ptr);
	for ( y=0; y<h; y++ ) {
		free(row_pointers[y]);
	}
	free(row_pointers);
	fclose(fh);

#else
	STATUS("No PNG support.\n");
#endif
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	FILE *fh;
	UnitCell *cell;
	char *filename;
	unsigned int angles[180];
	int i;
	char *prefix = NULL;
	char *geomfile = NULL;
	int config_basename = 0;
	int config_checkprefix = 1;
	struct detector *det;
	int config_angles = 0;
	const signed int ht = 3;
	const signed int kt = 4;
	const signed int lt = 5;
	const int gs = 16;
	double *vals;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{"prefix",             1, NULL,               'x'},
		{"basename",           0, &config_basename,    1},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:x:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			infile = strdup(optarg);
			break;

		case 'g' :
			geomfile = strdup(optarg);
			break;

		case 'x' :
			prefix = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( infile == NULL ) infile = strdup("-");
	if ( strcmp(infile, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(infile, "r");
	}
	if ( fh == NULL ) {
		ERROR("Couldn't open input stream '%s'\n", infile);
		return ENOENT;
	}

	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	/* Initialise histogram */
	for ( i=0; i<180; i++ ) angles[i] = 0;

	/* Initialise shape transform array */
	vals = calloc(gs*gs*gs, sizeof(double));

	/* Loop over all successfully indexed patterns */
	while ( find_chunk(fh, &cell, &filename) == 0 ) {

		struct image image;
		struct hdfile *hdfile;
		double ang;
		double ax, ay, az;
		double bx, by, bz;
		double cx, cy, cz;
		unsigned int bin;
		struct rvec q;
		int x, y;

		STATUS("Processing '%s'\n", filename);

		hdfile = hdfile_open(filename);
		if ( hdfile == NULL ) {
			return ENOENT;
		} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
			ERROR("Couldn't select path\n");
			return ENOENT;
		}

		hdf5_read(hdfile, &image, 0);

		cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by,
		                         &bz, &cx, &cy, &cz);

		for ( x=0; x<image.width; x++ ) {
		for ( y=0; y<image.height; y++ ) {

			double hd, kd, ld;
			signed int h, k, l;
			double dh, dk, dl;

			q = get_q(&image, x, y, 1, NULL, 1.0/image.lambda);

			hd = q.u * ax + q.v * ay + q.w * az;
			kd = q.u * bx + q.v * by + q.w * bz;
			ld = q.u * cx + q.v * cy + q.w * cz;

			h = (signed int)rint(hd);
			k = (signed int)rint(kd);
			l = (signed int)rint(ld);

			if ( !((h==ht) && (k==kt) && (l==lt)) ) continue;

			dh = hd - h;
			dk = kd - k;
			dl = ld - l;

			double v = image.data[x+image.width*y];

			interpolate_onto_grid(vals, v, gs, gs, gs, dh, dk, dl);

		}
		}

		if ( config_angles ) {

			double asx, asy, asz;
			double bsx, bsy, bsz;
			double csx, csy, csz;

			cell_get_reciprocal(cell, &asx, &asy, &asz,
			                          &bsx, &bsy, &bsz,
			                          &csx, &csy, &csz);
			ang = angle_between(csx, csy, csz, 0.0, 0.0, 1.0);
			ang = rad2deg(ang);  /* 0->180 deg */
			bin = rint(ang);
			angles[bin]++;

		}

	}

	for ( i=0; i<gs; i++ ) {
		char line[64];
		snprintf(line, 63, "slice-%i.png", i);
		write_slice(line, vals, i, gs, gs, gs);
	}

	if ( config_angles ) {
		for ( i=0; i<180; i++ ) {
			STATUS("%i %i\n", i, angles[i]);
		}
	}

	return 0;
}
