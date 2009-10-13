/*
 * image.h
 *
 * Handle images and image features
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMAGE_H
#define IMAGE_H

#include <stdint.h>


typedef enum {
	FORMULATION_CLEN,
	FORMULATION_PIXELSIZE
} FormulationMode;


typedef struct imagefeature_struct {

	struct image			*parent;
	double				x;
	double				y;
	double				intensity;

	/* Partner for this feature (in another feature list) or NULL */
	struct imagefeature_struct	*partner;

	/* Distance between this feature and its partner, if any. */
	double				partner_d;

	/* The reflection this was projected from, if any */
	struct reflection_struct	*reflection;

} ImageFeature;


typedef struct {

	ImageFeature		*features;
	int			n_features;

} ImageFeatureList;


struct image {

	uint16_t		*data;

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
	double			lambda;

	int			width;
	int			height;
	double			x_centre;
	double			y_centre;

	ImageFeatureList	*features;	/* "Experimental" features */
	ImageFeatureList	*rflist;	/* "Predicted" features */

};


typedef struct imagelist_struct {

	int		n_images;
	struct image	*images;

} ImageList;


extern ImageList *image_list_new(void);

extern int image_add(ImageList *list, struct image *image);

extern ImageFeatureList *image_feature_list_new(void);

extern void image_feature_list_free(ImageFeatureList *flist);

extern void image_add_feature(ImageFeatureList *flist, double x, double y,
                              struct image *parent, double intensity);

extern void image_add_feature_reflection(ImageFeatureList *flist,
                                         double x, double y,
                                         struct image *parent,
                                         double intensity);

extern ImageFeature *image_feature_closest(ImageFeatureList *flist,
                                           double x, double y, double *d,
                                           int *idx);

extern ImageFeature *image_feature_second_closest(ImageFeatureList *flist,
                                                  double x, double y, double *d,
                                                  int *idx);


#endif	/* IMAGE_H */
