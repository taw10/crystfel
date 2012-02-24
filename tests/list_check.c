/*
 * list_check.c
 *
 * Unit test for the reflection list module
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include <reflist.h>


struct refltemp {
	signed int h;
	signed int k;
	signed int l;
	int num;
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

	h = RANDOM_INDEX;
	k = RANDOM_INDEX;
	l = RANDOM_INDEX;

	for ( i=0; i<num_items; i++ ) {

		int j;
		int num;

		if ( random() > RAND_MAX/2 ) {
			h = RANDOM_INDEX;
			k = RANDOM_INDEX;
			l = RANDOM_INDEX;
		} /* else use the same as last time */

		/* Count the number of times this reflection appeared before */
		num = 1;
		for ( j=0; j<i; j++ ) {
			if ( (check[j].h == h)
			  && (check[j].k == k)
			  && (check[j].l == l) ) {
				num++;
			}
		}

		/* Update all copies with this number */
		for ( j=0; j<i; j++ ) {
			if ( (check[j].h == h)
			  && (check[j].k == k)
			  && (check[j].l == l) ) {
				check[j].num = num;
			}
		}

		add_refl(list, h, k, l);
		check[i].h = h;
		check[i].k = k;
		check[i].l = l;
		check[i].num = num;
		check[i].found = 0;

	}

	printf("num_reflections is %i, tree depth is %i\n",
	       num_reflections(list), tree_depth(list));

	/* Iterate over the list and check we find everything */
	int count = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		for ( i=0; i<num_items; i++ ) {
			if ( (check[i].h == h)
			  && (check[i].k == k)
			  && (check[i].l == l)
			  && (check[i].found == 0) ) {
				check[i].found = 1;
			}
		}

		count++;

	}
	if ( count != num_reflections(list) ) {
		fprintf(stderr, "num_reflections gave %i, iteration gave %i\n",
		        num_reflections(list), count);
		return 1;
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

	/* Check that all the reflections can be found */
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

	}

	reflist_free(list);
	free(check);

	return 0;
}

int main(int argc, char *argv[])
{
	int i;

	printf("Running list test...\n");

	for ( i=0; i<100; i++ ) {
		if ( test_lists(4096*random()/RAND_MAX) ) return 1;
	}

	return 0;
}
