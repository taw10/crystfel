/*
 * main.c
 *
 * "Main"
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
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


static void main_show_help(const char *s)
{
	printf("Syntax: %s <file1.h5> <file2.h5> [...]\n\n", s);
	printf("Index diffraction patterns\n\n");
	printf("  -h              Display this help message\n");
}


int main(int argc, char *argv[])
{
	int c;
	char **in_files;
	size_t nin;
	size_t i;

	while ((c = getopt(argc, argv, "h")) != -1) {

		switch ( c ) {

			case 'h' : {
				main_show_help(argv[0]);
				return 0;
			}

		}

	}

	if ( optind < argc ) {
		nin = argc-optind;
		in_files = malloc(nin*sizeof(char *));
		for ( i=0; i<nin; i++ ) {
			in_files[i] = strdup(argv[optind+i]);
		}
	} else {
		fprintf(stderr, "No input files!\n");
		return 1;
	}

	printf("Input files (%i):\n", nin);
	for ( i=0; i<nin; i++ ) {
		printf("%6i: %s\n", i+1, in_files[i]);
	}

	return 0;
}
