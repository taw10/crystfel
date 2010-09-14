/*
 * facetron.c
 *
 * Profile fitting for coherent nanocrystallography
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

#include "image.h"
#include "cell.h"
#include "hdf5-file.h"
#include "detector.h"
#include "geometry.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Profile fitting for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"                              (must be a file, not e.g. stdin)\n"
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


static char *twiddle_filename(char *filename, int config_basename,
                              const char *prefix)
{
	char *f3;
	size_t len;

	if ( config_basename ) {
		char *f2;
		/* Happy with either the GNU or non-GNU versions of basename()
		 * here.  Because _GNU_SOURCE is defined in configure.ac, we
		 * should have the GNU version. */
		f2 = strdup(basename(filename));
		free(filename);
		filename = f2;
	}

	len = strlen(prefix)+strlen(filename)+1;
	f3 = malloc(len);
	snprintf(f3, len, "%s%s", prefix, filename);

	free(filename);
	return f3;
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	FILE *fh;
	UnitCell *cell;
	char *filename;
	struct detector *det;
	char *geometry = NULL;
	char *prefix = NULL;
	int config_basename = 0;
	int config_checkprefix = 1;

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
			geometry = strdup(optarg);
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
	fh = fopen(infile, "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open input stream '%s'\n", infile);
		return ENOENT;
	}
	free(infile);

	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	det = get_detector_geometry(geometry);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	/* Loop over all successfully indexed patterns */
	while ( find_chunk(fh, &cell, &filename) == 0 ) {

		struct image image;
		struct hdfile *hdfile;
		struct reflhit *hits;
		int np;

		filename = twiddle_filename(filename, config_basename, prefix);

		STATUS("Integrating intensities from '%s'\n", filename);

		image.det = det;

		hdfile = hdfile_open(filename);
		if ( hdfile == NULL ) {
			return ENOENT;
		} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
			ERROR("Couldn't select path\n");
			return ENOENT;
		}

		hdf5_read(hdfile, &image, 0);

		hits = find_intersections(&image, cell, 5.0e-3, 0.0001, &np, 1);

		hdfile_close(hdfile);
		cell_free(cell);
		free(filename);
		free(image.data);
		free(image.flags);

	}

	free(det->panels);
	free(det);
	fclose(fh);
	free(prefix);

	return 0;
}
