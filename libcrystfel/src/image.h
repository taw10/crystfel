/*
 * image.h
 *
 * Handle images and image features
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *
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

#ifndef IMAGE_H
#define IMAGE_H

#include <stdint.h>
#include <complex.h>
#include <sys/types.h>

struct imagefeature;
struct image;

#include "utils.h"
#include "cell.h"
#include "reflist.h"
#include "crystal.h"
#include "index.h"
#include "spectrum.h"
#include "datatemplate.h"

/**
 * \file image.h
 *
 * Information about images
 */

/** Represents a peak in an image. */
struct imagefeature {

	/** \name Coordinates on panel (fast scan, slow scan)
	 *  Note carefully that these are the distances, measured in pixels,
	 *  from the corner of the panel.  They are NOT pixel indices.
	 *  If the peak is in the middle of the first pixel, its coordinates would be
	 *  0.5,0.5. */
	/**@{*/
	double                          fs;
	double                          ss;
	/**@}*/

	int                             pn;         /**< Panel number */
	double                          intensity;  /**< Intensity */

	/** \name Reciprocal space coordinates (m^-1) of this feature */
	/** @{ */
	double                          rx;
	double                          ry;
	double                          rz;
	/** @} */

	const char                      *name;  /**< Text name, e.g. "5,3,-1" */
};


/** An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;


struct image
{
	/** The image data, by panel */
	float                   **dp;

	/** The bad pixel mask, by panel */
	int                     **bad;

	/** The per-pixel saturation values, by panel */
	float                   **sat;

	/** Non-zero if the frame was determined to be a "hit" */
	int                     hit;

	/**Array of crystals in the image */
	Crystal                 **crystals;

	/** The number of crystals in the image (size of \p crystals) */
	int                     n_crystals;

	/** Indexing method which indexed this pattern */
	IndexingMethod          indexed_by;

	/** Number of times the indexer was tried before succeeding */
	int                     n_indexing_tries;

	/** The detector structure */
	struct detgeom          *detgeom;

	/** \name The filename and event ID for the image
	 * @{ */
	char                    *filename;
	char                    *ev;
	/** @} */

	/** A list of metadata read from the stream */
	char                    *copied_headers;

	/** Mean of the camera length values for all panels */
	double                  avg_clen;

	/** ID number of the worker processing handling this image */
	int                     id;

	/** Monotonically increasing serial number for this image */
	int                     serial;

	/** Spectrum information */
	Spectrum               *spectrum;

	/** Wavelength of the incident radiation, in metres */
	double                  lambda;

	/** Convergence angle of the incident ration, in radians (full angle) */
	double                  div;

	/** Full-width half-maximum bandwidth as a fraction, applied to wavelength */
	double                  bw;

	/** Resolution estimate based on peaks */
	double                  peak_resolution;

	/** List of peaks found in the image */
	ImageFeatureList        *features;

};

#ifdef __cplusplus
extern "C" {
#endif

/* Feature lists */
extern ImageFeatureList *image_feature_list_new(void);

extern void image_feature_list_free(ImageFeatureList *flist);

extern void image_add_feature(ImageFeatureList *flist, double x, double y,
                              int pn,
                              struct image *parent, double intensity,
                              const char *name);

extern void image_remove_feature(ImageFeatureList *flist, int idx);

extern struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                                  double fs, double ss,
                                                  int pn,
                                                  double *d, int *idx);

extern int image_feature_count(ImageFeatureList *flist);
extern struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx);
extern const struct imagefeature *image_get_feature_const(const ImageFeatureList *flist,
                                                          int idx);
extern ImageFeatureList *sort_peaks(ImageFeatureList *flist);
extern ImageFeatureList *image_feature_list_copy(const ImageFeatureList *flist);

extern void image_add_crystal(struct image *image, Crystal *cryst);
extern int remove_flagged_crystals(struct image *image);
extern void free_all_crystals(struct image *image);

extern void mark_resolution_range_as_bad(struct image *image,
                                         double min, double max);

extern struct image *image_new(void);
extern struct image *image_read(DataTemplate *dtempl, const char *filename,
                                const char *event);
extern void image_free(struct image *image);

extern ImageFeatureList *image_read_peaks(const DataTemplate *dtempl,
                                          const char *filename,
                                          const char *event,
                                          int half_pixel_shift);

extern char **image_expand_frames(const DataTemplate *dtempl,
                                  const char *filename, int *nframes);

/* The following functions are not part of the public API -
 * use within libcrystfel only */
extern void create_detgeom(struct image *image,
                           const DataTemplate *dtempl);

extern int create_blank_arrays(struct image *image);

#ifdef __cplusplus
}
#endif

#endif	/* IMAGE_H */
