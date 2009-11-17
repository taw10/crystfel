/*
 * image.h
 *
 * Handle images and image features
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMAGE_H
#define IMAGE_H

#include <stdint.h>
#include <complex.h>


/* How is the scaling of the image described? */
typedef enum {
	FORMULATION_CLEN,
	FORMULATION_PIXELSIZE
} FormulationMode;


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

};

/* An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;


/* A 3D vector in reciprocal space */
struct threevec
{
	double   u;
	double   v;
	double   w;
};


/* Structure describing an image */
struct image {

	uint16_t		*data;
	double complex		*sfacs;
	struct threevec		*qvecs;
	double			*twotheta;
	struct molecule		*molecule;

	/* Radians.  Defines where the pattern lies in reciprocal space */
	double			tilt;

	/* Radians.  Defines where the pattern lies in reciprocal space */
	double			omega;

	/* Image parameters can be given as camera length or pixel size.
	 * If FORMULATION_CLEN, then camera_len and resolution must be given.
	 * If FORMULATION_PIXELSIZE, then pixel_size only is needed.*/
	FormulationMode		fmode;
	double			pixel_size;
	double			camera_len;
	double			resolution;	/* pixels per metre */

	/* Wavelength must always be given */
	double			lambda;		/* Wavelength in m */
	double			xray_energy;	/* X-ray energy
	                                         * in J (per photon) */

	int			width;
	int			height;
	double			x_centre;
	double			y_centre;

	ImageFeatureList	*features;	/* "Experimental" features */
	ImageFeatureList	*rflist;	/* "Predicted" features */

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
                              struct image *parent, double intensity);

extern struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                                  double x, double y, double *d,
                                                  int *idx);

extern int image_feature_count(ImageFeatureList *flist);
extern struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx);

#endif	/* IMAGE_H */
