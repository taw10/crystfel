/*
 * alter_stream.c
 *
 * Do transformations on a stream
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include <assert.h>

#include "../src/utils.h"
#include "../src/stream.h"
#include "../src/cell.h"
#include "../src/image.h"
#include "../src/beam-parameters.h"
#include "../src/geometry.h"
#include "../src/peaks.h"


static void mess_up_cell(UnitCell *cell)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	/* Cell noise in percent */
	const double cnoise = 0.5;

	//STATUS("Real:\n");
	//cell_print(cell);

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	ax = gaussian_noise(ax, cnoise*fabs(ax)/100.0);
	ay = gaussian_noise(ay, cnoise*fabs(ay)/100.0);
	az = gaussian_noise(az, cnoise*fabs(az)/100.0);
	bx = gaussian_noise(bx, cnoise*fabs(bx)/100.0);
	by = gaussian_noise(by, cnoise*fabs(by)/100.0);
	bz = gaussian_noise(bz, cnoise*fabs(bz)/100.0);
	cx = gaussian_noise(cx, cnoise*fabs(cx)/100.0);
	cy = gaussian_noise(cy, cnoise*fabs(cy)/100.0);
	cz = gaussian_noise(cz, cnoise*fabs(cz)/100.0);
	cell_set_reciprocal(cell, ax, ay, az, bx, by, bz, cx, cy, cz);

	//STATUS("Changed:\n");
	//cell_print(cell);
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Alter a stream.\n"
"\n"
" -h, --help              Display this help message.\n"
"\n"
"You need to provide the following basic options:\n"
" -i, --input=<file>      Read reflections from <file>.\n"
" -o, --output=<file>     Write partials in stream format to <file>.\n"
);
}


int main(int argc, char *argv[])
{
	int c;
	char *input_file = NULL;
	char *output_file = NULL;
	struct image image;
	struct detector *det = NULL;
	struct beam_params *beam = NULL;
	char *beamfile = NULL;
	char *geomfile = NULL;
	FILE *ifh;
	FILE *ofh;
	int v;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"input",              1, NULL,               'i'},
		{"beam",               1, NULL,               'b'},
		{"geometry",           1, NULL,               'g'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:b:g:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'o' :
			output_file = strdup(optarg);
			break;

		case 'i' :
			input_file = strdup(optarg);
			break;

		case 'b' :
			beamfile = strdup(optarg);
			break;

		case 'g' :
			geomfile = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( input_file == NULL ) {
		ERROR("You must pgive a filename for the output.\n");
		return 1;
	}
	ifh = fopen(input_file, "r");
	if ( ifh == NULL ) {
		ERROR("Couldn't open input file '%s'\n", input_file);
		return 1;
	}
	free(input_file);

	if ( output_file == NULL ) {
		ERROR("You must pgive a filename for the output.\n");
		return 1;
	}
	ofh = fopen(output_file, "w");
	if ( ofh == NULL ) {
		ERROR("Couldn't open output file '%s'\n", output_file);
		return 1;
	}
	free(output_file);

	/* Load beam */
	if ( beamfile == NULL ) {
		ERROR("You need to provide a beam parameters file.\n");
		return 1;
	}
	beam = get_beam_parameters(beamfile);
	if ( beam == NULL ) {
		ERROR("Failed to load beam parameters from '%s'\n", beamfile);
		return 1;
	}
	free(beamfile);

	/* Load geometry */
	if ( geomfile == NULL ) {
		ERROR("You need to give a geometry file.\n");
		return 1;
	}
	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	write_stream_header(ofh, argc, argv);

	image.det = det;
	image.width = det->max_fs;
	image.height = det->max_ss;

	image.lambda = 0.0;
	image.div = beam->divergence;
	image.bw = beam->bandwidth;
	image.profile_radius = 0.0001e9;
	image.i0_available = 0;

	do {

		image.indexed_cell = NULL;
		v = read_chunk(ifh, &image);
		if ( (v == 0) && (image.indexed_cell != NULL) ) {

			mess_up_cell(image.indexed_cell);
			image.reflections = find_intersections(&image,
			                               image.indexed_cell);

			write_chunk(ofh, &image, STREAM_INTEGRATED);

		}

	} while ( v == 0 );


	fclose(ofh);
	fclose(ifh);

	return 0;
}
