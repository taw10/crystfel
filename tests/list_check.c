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
	int found;
};

#define RANDOM_INDEX (128*random()/RAND_MAX - 256*random()/RAND_MAX)


static int test_lists(int num_items)
{
	struct refltemp *check;
	RefList *list;
	int i;
	signed int h, k, l;
	Reflection *refl;
	RefListIterator *iter;

	check = malloc(num_items * sizeof(struct refltemp));
	list = reflist_new();

	printf("Testing with %i\n", num_items);

	h = RANDOM_INDEX;
	k = RANDOM_INDEX;
	l = RANDOM_INDEX;

	for ( i=0; i<num_items; i++ ) {

		int j;
		int duplicate = 0;

		//if ( random() > RAND_MAX/2 ) {
			h = RANDOM_INDEX;
			k = RANDOM_INDEX;
			l = RANDOM_INDEX;
		//} /* else use the same as last time */

		/* Count the number of times this reflection appeared before */
		for ( j=0; j<i; j++ ) {
			if ( (check[j].h == h)
			  && (check[j].k == k)
			  && (check[j].l == l) ) {
				duplicate++;
			}
		}

		/* Update all copies with this number */
		for ( j=0; j<i; j++ ) {
			if ( (check[j].h == h)
			  && (check[j].k == k)
			  && (check[j].l == l) ) {
				check[j].dup = duplicate;
			}
		}

		add_refl(list, h, k, l);
		check[i].h = h;
		check[i].k = k;
		check[i].l = l;
		check[i].del = 0;
		check[i].dup = duplicate;
		check[i].found = 0;

	}

	/* Iterate over the list and check we find everything */
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;

		get_indices(refl, &h, &k, &l);
		printf("%3i %3i %3i\n", h, k, l);

		for ( i=0; i<num_items; i++ ) {
			if ( (check[i].h == h)
			  && (check[i].k == k)
			  && (check[i].l == l)
			  && (check[i].found == 0) ) {
				check[i].found = 1;
			}
		}

	}
	for ( i=0; i<num_items; i++ ) {
		if ( check[i].found == 0 ) {

			Reflection *test;

			fprintf(stderr, "Iteration didn't find %3i %3i %3i %i\n",
			        check[i].h, check[i].k, check[i].l, i);
			test = find_refl(list, check[i].h, check[i].k,
			                 check[i].l);
			if ( test == NULL ) {
				fprintf(stderr, "Not in list\n");
			} else {
				fprintf(stderr, "But found in list.\n");
			}
			return 1;

		}
	}

	/* Check that all the reflections can be found,
	 * and delete the first few. */
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

		/* Delete some reflections */
		if ( i<num_items/2 ) {

			int j;

			delete_refl(refl);
			check[i].del = 1;

			/* Update all counts */
			for ( j=0; j<i; j++ ) {
				if ( (check[j].h == h)
				  && (check[j].k == k)
				  && (check[j].l == l) ) {
					check[j].dup--;
				}
			}

		}

	}

	/* Check that the deleted reflections can no longer be found, unless
	 * duplicated.  If duplicated, remove all the remaining copies. */
	for ( i=0; i<num_items/2; i++ ) {

		signed int h, k, l;
		Reflection *refl;

		h = check[i].h;
		k = check[i].k;
		l = check[i].l;

		refl = find_refl(list, h, k, l);
		if ( refl != NULL ) {

			/* Whoops, found it.  Was it a duplicate? */
			if ( check[i].dup == 0 ) {
				fprintf(stderr, "Found %3i %i %3i after"
				        " deletion.\n", h, k, l);
				return 1;
			} else {

				int j;
				Reflection *c;

				for ( j=0; j<check[i].dup; j++ ) {
					Reflection *r2;
					r2 = find_refl(list, h, k, l);
					if ( r2 == NULL ) {
						fprintf(stderr, "Found too few"
						        " duplicates for"
						        " %3i %3i %3i\n",
						        h, k, l);
						return 1;
					}
					delete_refl(r2);
				}

				c = find_refl(list, h, k, l);
				if ( c != NULL ) {
					fprintf(stderr, "Found too many "
					        "duplicates for %3i %3i %3i\n",
					        h, k, l);
					return 1;
				}

			}
		}

	}

	reflist_free(list);
	free(check);

	return 0;
}

int main(int argc, char *argv[])
{
	int i;

	printf("Running list test...");
	fflush(stdout);

	for ( i=0; i<1; i++ ) {
		if ( test_lists(4096*random()/RAND_MAX) ) return 1;
	}

	printf("\r");

	return 0;
}
