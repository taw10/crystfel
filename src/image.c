/*
 * image.c
 *
 * Handle images and image features
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "image.h"
#include "utils.h"


struct _imagelist
{
	int		n_images;
	struct image	*images;
};


struct _imagefeaturelist
{
	struct imagefeature	*features;
	int			n_features;
};


int image_add(ImageList *list, struct image *image)
{
	if ( list->images ) {
		list->images = realloc(list->images,
		                       (list->n_images+1)*sizeof(struct image));
	} else {
		assert(list->n_images == 0);
		list->images = malloc(sizeof(struct image));
	}

	/* Copy the metadata */
	list->images[list->n_images] = *image;

	/* Fill out some extra fields */
	list->images[list->n_images].features = NULL;

	list->n_images++;

	return list->n_images - 1;
}


ImageList *image_list_new()
{
	ImageList *list;

	list = malloc(sizeof(ImageList));

	list->n_images = 0;
	list->images = NULL;

	return list;
}


void image_add_feature(ImageFeatureList *flist, double x, double y,
                       struct image *parent, double intensity, const char *name)
{
	if ( flist->features ) {
		flist->features = realloc(flist->features,
		                    (flist->n_features+1)
		                    *sizeof(struct imagefeature));
	} else {
		assert(flist->n_features == 0);
		flist->features = malloc(sizeof(struct imagefeature));
	}

	flist->features[flist->n_features].x = x;
	flist->features[flist->n_features].y = y;
	flist->features[flist->n_features].intensity = intensity;
	flist->features[flist->n_features].parent = parent;
	flist->features[flist->n_features].name = name;
	flist->features[flist->n_features].valid = 1;

	flist->n_features++;

}


ImageFeatureList *image_feature_list_new()
{
	ImageFeatureList *flist;

	flist = malloc(sizeof(ImageFeatureList));

	flist->n_features = 0;
	flist->features = NULL;

	return flist;
}


void image_feature_list_free(ImageFeatureList *flist)
{
	if ( !flist ) return;

	if ( flist->features ) free(flist->features);
	free(flist);
}


struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                           double x, double y,
                                           double *d, int *idx)
{
	int i;
	double dmin = +HUGE_VAL;
	int closest = 0;

	for ( i=0; i<flist->n_features; i++ ) {

		double ds;

		ds = distance(flist->features[i].x, flist->features[i].y, x, y);

		if ( ds < dmin ) {
			dmin = ds;
			closest = i;
		}

	}

	if ( dmin < +HUGE_VAL ) {
		*d = dmin;
		*idx = closest;
		return &flist->features[closest];
	}

	*d = +INFINITY;
	return NULL;
}


int image_feature_count(ImageFeatureList *flist)
{
	if ( flist == NULL ) return 0;
	return flist->n_features;
}


struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx > flist->n_features ) return NULL;

	if ( flist->features[idx].valid == 0 ) return NULL;

	return &flist->features[idx];
}


void image_remove_feature(ImageFeatureList *flist, int idx)
{
	flist->features[idx].valid = 0;
}
