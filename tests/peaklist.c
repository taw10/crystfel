/*
 * peaklist.c
 *
 * Unit test for the peak ("feature") list
 *
 * Copyright © 2026 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2026 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>

#include <image.h>

static int approx(double a, int b)
{
    return (fabs(a-b) < 0.1);
}


int main(int argc, char *argv[])
{
    int i;
    ImageFeatureList *list;
    ImageFeatureList *list2;
    const int n = 100;
    int nf, nn;

    list = image_feature_list_new();

    for ( i=0; i<n-1; i++ ) {
        image_add_feature(list, 123.4, 567.8, 0, i, NULL);
    }
    image_add_feature(list, 100.0, 200.0, 0, n-1, NULL);

    nf = image_feature_count(list);
    if ( nf != n ) {
        fprintf(stderr, "Wrong number of peaks (%i, should be %i)\n", nf, n);
        return 1;
    }

    for ( i=0; i<n; i++ ) {
        const struct imagefeature *f = image_get_feature_const(list, i);
        if ( !approx(f->intensity, i) ) {
            fprintf(stderr, "Peak %.1f found in wrong position (%i)\n",
                    f->intensity, i);
            return 1;
        }
    }

    list2 = image_feature_list_copy(list);

    nn = image_feature_count(list2);
    if ( nn != n ) {
        fprintf(stderr, "Wrong number of peaks after copying: %i, should be %i\n", nn, n);
        return 1;
    }

    for ( i=0; i<n; i++ ) {
        const struct imagefeature *f = image_get_feature_const(list2, i);
        if ( !approx(f->intensity, i) ) {
            fprintf(stderr, "Peak %.1f found in wrong position in copy (%i)\n",
                    f->intensity, i);
            return 1;
        }
    }

    image_feature_list_free(list2);

    list2 = sort_peaks(list);
    nn = image_feature_count(list2);
    if ( nn != n ) {
        fprintf(stderr, "Wrong number of peaks after sorting: %i, should be %i\n", nn, n);
        return 1;
    }

    for ( i=0; i<n; i++ ) {
        const struct imagefeature *f = image_get_feature_const(list2, i);
        if ( !approx(f->intensity, n-i-1) ) {
            fprintf(stderr, "Peak %.1f found in wrong position in sorted (%i)\n",
                    f->intensity, i);
            return 1;
        }
    }

    image_feature_list_free(list2);

    image_remove_feature(list, 98);
    image_remove_feature(list, 50);
    image_remove_feature(list, 0);

    nn = image_feature_count(list);
    if ( nn != n-3 ) {
        fprintf(stderr, "Wrong number of peaks after deleting: %i, should be %i\n", nn, n);
        return 1;
    }

    for ( i=0; i<nn; i++ ) {
        const struct imagefeature *f = image_get_feature_const(list, i);
        if ( approx(f->intensity, 0) ) {
            fprintf(stderr, "Found deleted peak 0 at position %i\n", i);
            return 1;
        }
        if ( approx(f->intensity, 50) ) {
            fprintf(stderr, "Found deleted peak 50 at position %i\n", i);
            return 1;
        }
        if ( approx(f->intensity, 98) ) {
            fprintf(stderr, "Found deleted peak 98 at position %i\n", i);
            return 1;
        }
    }

    int idx;
    double d;
    struct imagefeature *c = image_feature_closest(list, 101.0, 205.0, 0,
                                                   &d, &idx);
    if ( !approx(c->intensity, 99) ) {
        fprintf(stderr, "Found wrong reflection (%.1f)\n", c->intensity);
        return 1;
    }
    if ( idx != 96 ) {
        fprintf(stderr, "Found reflection in wrong place (%i)\n", idx);
        return 1;
    }

    image_feature_list_free(list);

    return 0;
}
