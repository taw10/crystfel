/*
 * image.h
 *
 * Handle images and image features
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMAGE_H
#define IMAGE_H

#include <stdint.h>
#include <complex.h>
#include <sys/types.h>

#include "utils.h"
#include "cell.h"
#include "detector.h"
#include "reflist.h"


#define MAX_CELL_CANDIDATES (32)


/* Structure describing a feature in an image */
struct imagefeature {

	struct image                    *parent;
	double                          x;
	double                          y;
	double                          intensity;

	/* Reciprocal space coordinates (m^-1 of course) of this feature */
	double                          rx;
	double                          ry;
	double                          rz;

	/* Internal use only */
	int                             valid;

	const char                      *name;
};

/* An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;


/* Structure describing an image */
struct image {

	float                   *data;
	uint16_t                *flags;
	double                  *twotheta;
	UnitCell                *indexed_cell;
	UnitCell                *candidate_cells[MAX_CELL_CANDIDATES];
	int                     ncells;
	struct detector         *det;
	struct beam_params      *beam;  /* The nominal beam parameters */
	char                    *filename;

	int                     id;   /* ID number of the thread
	                               * handling this image */

	/* Information about the crystal */
	double                  m;  /* Mosaicity in radians */

	/* Per-shot radiation values */
	double                  lambda;        /* Wavelength in m */
	double                  div;           /* Divergence in radians */
	double                  bw;            /* Bandwidth as a fraction */
	double                  i0;            /* Incident intensity */
	int                     i0_available;  /* 0 if f0 wasn't available
	                                        * from the input. */
	double                  osf;           /* Overall scaling factor */
	double                  profile_radius; /* Radius of reflection */

	int                     width;
	int                     height;

	/* Integrated (or about-to-be-integrated) reflections */
	RefList                 *reflections;

	/* Detected peaks */
	ImageFeatureList        *features;

};

/* An opaque type representing a list of images */
typedef struct _imagelist ImageList;


/* Image lists */
extern ImageList *image_list_new(void);

extern int image_add(ImageList *list, struct image *image);


/* Feature lists */
extern ImageFeatureList *image_feature_list_new(void);

extern void image_feature_list_free(ImageFeatureList *flist);

extern void image_add_feature(ImageFeatureList *flist, double x, double y,
                              struct image *parent, double intensity,
                              const char *name);

extern void image_remove_feature(ImageFeatureList *flist, int idx);

extern struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                                  double x, double y, double *d,
                                                  int *idx);

extern int image_feature_count(ImageFeatureList *flist);
extern struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx);

#endif	/* IMAGE_H */
