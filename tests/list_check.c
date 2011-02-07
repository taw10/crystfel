/*
 * list_check.c
 *
 * Unit test for the reflection list module
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include "../src/reflist.h"


struct refltemp {
	signed int h;
	signed int k;
	signed int l;
	int del;
	int dup;
};

#define RANDOM_INDEX (128*random()/RAND_MAX - 256*random()/RAND_MAX)


static int test_lists(int num_items)
{
	struct refltemp *check;
	RefList *list;
	int i;

	check = malloc(num_items * sizeof(struct refltemp));
	list = reflist_new();

	for ( i=0; i<num_items; i++ ) {

		signed int h, k, l;
		int j;
		int duplicate = 0;

		h = RANDOM_INDEX;
		k = RANDOM_INDEX;
		l = RANDOM_INDEX;

		for ( j=0; j<i; j++ ) {
			if ( (check[j].h == h)
			  && (check[j].k == k)
			  && (check[j].l == l) ) {
				duplicate++;
			}
		}

		add_refl(list, h, k, l);
		check[i].h = h;
		check[i].k = k;
		check[i].l = l;
		check[i].del = 0;
		check[i].dup = duplicate;

	}

	for ( i=0; i<num_items; i++ ) {

		signed int h, k, l;
		Reflection *refl;

		h = check[i].h;
		k = check[i].k;
		l = check[i].l;

		refl = find_refl(list, h, k, l);
		if ( refl == NULL ) {
			fprintf(stderr, "Couldn't find %3i %3i %3i\n", h, k, l);
			return 1;
		}

		if ( i<num_items/2 ) {
			delete_refl(refl);
			check[i].del = 1;
		}

	}

	for ( i=0; i<num_items/2; i++ ) {

		signed int h, k, l;
		Reflection *refl;

		h = check[i].h;
		k = check[i].k;
		l = check[i].l;

		refl = find_refl(list, h, k, l);
		if ( refl == NULL ) {
			fprintf(stderr, "Couldn't find %3i %3i %3i\n", h, k, l);
			return 1;
		}

		if ( i<num_items/2 ) {
			delete_refl(refl);
			check[i].del = 1;
		}

	}

	free(check);

	return 0;
}

int main(int argc, char *argv[])
{
	int i;

	for ( i=0; i<100; i++ ) {
		if ( test_lists(random()) ) return 1;
	}

	return 0;
}
