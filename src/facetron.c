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

#include "image.h"
#include "cell.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Profile fitting for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the image to work on.\n"
"  -m, --matrix=<filename>    Specify the file which contains the initial\n"
"                              orientation matrix.  Can be '-' for stdin,\n"
"                              which is the default.  Units are nm^-1.\n"
);
}


static UnitCell *read_orientation_matrix(const char *filename)
{
	FILE *mfh;
	float u, v, w;
	struct rvec as, bs, cs;
	UnitCell *cell;

	if ( (filename == NULL) || (strcmp(filename, "-") == 0) ) {
		mfh = stdin;
	} else {
		mfh = fopen(filename, "r");
	}
	if ( mfh == NULL ) {
		ERROR("Failed to open matrix file '%s'\n", filename);
		return NULL;
	}

	if ( fscanf(mfh, "%f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read a-star\n");
		return NULL;
	}
	as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
	if ( fscanf(mfh, "%f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read b-star\n");
		return NULL;
	}
	bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
	if ( fscanf(mfh, "%f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read c-star\n");
		return NULL;
	}
	cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
	cell = cell_new_from_axes(as, bs, cs);
	fclose(mfh);

	return cell;
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *matrix = NULL;
	UnitCell *cell;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"matrix",             1, NULL,               'm'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:m:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			infile = strdup(optarg);
			break;

		case 'm' :
			matrix = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	cell = read_orientation_matrix(matrix);
	free(matrix);
	if ( cell == NULL ) {
		ERROR("Couldn't read initial orientation matrix.\n");
		return 1;
	}

	cell_print(cell);

	return 0;
}
