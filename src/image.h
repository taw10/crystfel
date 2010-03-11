/*
 * image.h
 *
 * Handle images and image features
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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

#include "utils.h"
#include "cell.h"
#include "detector.h"


/* Structure describing a feature in an image */
struct imagefeature {

	struct image			*parent;
	double				x;
	double				y;
	double				intensity;

	/* Partner for this feature (in another feature list) or NULL */
	struct imagefeature_struct	*partner;

	/* Distance between this feature and its partner, if any. */
	double				partner_d;

	/* Reciprocal space coordinates (m^-1 of course) of this feature */
	double				rx;
	double				ry;
	double				rz;

	/* Internal use only */
	int				valid;

	const char			*name;
};

/* An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;


/* A 3D vector in reciprocal space */
struct rvec
{
	double   u;
	double   v;
	double   w;
};


/* Structure describing an image */
struct image {

	float			*data;     /* Integer counts after bloom */
	double			*twotheta;
	struct molecule		*molecule;
	UnitCell		*indexed_cell;
	struct detector		det;

	struct quaternion	orientation;

	/* Wavelength must always be given */
	double			lambda;		/* Wavelength in m */

	int			width;
	int			height;

	ImageFeatureList	*features;	/* "Experimental" features */
	ImageFeatureList	*rflist;	/* "Predicted" features */

	/* DirAx auto-indexing low-level stuff */
	int			dirax_pty;
	pid_t			dirax_pid;
	char			*dirax_rbuffer;
	int			dirax_rbufpos;
	int			dirax_rbuflen;

	/* DirAx auto-indexing high-level stuff */
	int			dirax_step;
	int			dirax_read_cell;

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
