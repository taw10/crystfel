/*
 * image.h
 *
 * Handle images and image features
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2021 Thomas White <taw@physics.org>
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

	const char                      *name;  /**< Text name, e.g. "5,3,-1" */
};


/** An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;

typedef struct _image_data_arrays ImageDataArrays;


#define HEADER_CACHE_SIZE (128)

typedef enum
{
	HEADER_FLOAT,
	HEADER_INT,
	HEADER_STR
} HeaderCacheType;

struct header_cache_entry {
	char *header_name;
	HeaderCacheType type;
	union {
		long long int val_int;
		double val_float;
		char *val_str;
	};
};


typedef enum
{
	DATA_SOURCE_TYPE_UNKNOWN,
	DATA_SOURCE_TYPE_NONE,
	DATA_SOURCE_TYPE_HDF5,
	DATA_SOURCE_TYPE_CBF,
	DATA_SOURCE_TYPE_CBFGZ,
	DATA_SOURCE_TYPE_MSGPACK,
	DATA_SOURCE_TYPE_SEEDEE,
	DATA_SOURCE_TYPE_TIFF
} DataSourceType;


struct crystal_refls
{
	Crystal *cr;
	RefList *refls;
	int image_owns_crystal;
	int image_owns_refls;
};

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

	/** Array of crystals (with reflection lists) in the image */
	struct crystal_refls   *crystals;

	/** The number of crystals in the image (size of \p crystals) */
	int                     n_crystals;

	/** Indexing method which indexed this pattern */
	IndexingMethod          indexed_by;

	/** Number of times the indexer was tried before succeeding */
	int                     n_indexing_tries;

	/** The detector structure */
	struct detgeom          *detgeom;

	DataSourceType           data_source_type;

	/** \name The filename and event ID for the image
	 * @{ */
	char                    *filename;
	char                    *ev;
	/** @} */

	/** The data block, e.g. received over ZMQ, for the image.
	 * filenename/ev OR this should be filled in, but not both */
	void                    *data_block;
	size_t                   data_block_size;
	char                    *meta_data;

	/** A list of metadata read from the stream */
	struct header_cache_entry *header_cache[HEADER_CACHE_SIZE];
	int                        n_cached_headers;

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

	/** Re-usable data array structure, or NULL if not used */
	ImageDataArrays         *ida;

	/** If set, then 'features' should be freed with the image.
	 * Otherwise, it is managed externally (e.g. by Julia) */
	int                      owns_peaklist;

};

#ifdef __cplusplus
extern "C" {
#endif

/* File classifiers */
extern int is_hdf5_file(const char *filename, int *err);
extern int is_cbf_file(const char *filename, int *err);
extern int is_cbfgz_file(const char *filename, int *err);

/* Feature lists */
extern ImageFeatureList *image_feature_list_new(void);

extern void image_feature_list_free(ImageFeatureList *flist);

extern void image_add_feature(ImageFeatureList *flist, double x, double y,
                              int pn, double intensity, const char *name);

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
extern void image_add_crystal_refls(struct image *image,
                                    Crystal *cryst, RefList *reflist);
extern int remove_flagged_crystals(struct image *image);
extern void free_all_crystals(struct image *image);

extern void mark_resolution_range_as_bad(struct image *image,
                                         double min, double max);

extern struct image *image_new(void);
extern struct image *image_read(const DataTemplate *dtempl,
                                const char *filename,
                                const char *event,
                                int no_image_data,
                                int no_mask_data,
                                ImageDataArrays *ida);

extern struct image *image_create_for_simulation(const DataTemplate *dtempl);
extern struct image *image_read_data_block(const DataTemplate *dtempl,
                                           void *data_block,
                                           size_t data_block_size,
                                           char *meta_data,
                                           DataSourceType type,
                                           int serial,
                                           int no_image_data,
                                           int no_mask_data,
                                           ImageDataArrays *ida);
extern void image_free(struct image *image);

extern int image_read_header_float(struct image *image, const char *from,
                                   double *val);

extern int image_read_header_int(struct image *image, const char *from,
                                 long long int *val);

/* NB image_read_header_str should exist, but there is currently no use for it.
 * Get in touch if you disagree! */

extern void image_cache_header_float(struct image *image,
                                     const char *header_name,
                                     double header_val);

extern void image_cache_header_int(struct image *image,
                                   const char *header_name,
                                   long long int header_val);

extern void image_cache_header_str(struct image *image,
                                   const char *header_name,
                                   const char *header_val);

extern ImageFeatureList *image_read_peaks(const DataTemplate *dtempl,
                                          const char *filename,
                                          const char *event,
                                          int half_pixel_shift);

extern char **image_expand_frames(const DataTemplate *dtempl,
                                  const char *filename, int *nframes);

extern ImageDataArrays *image_data_arrays_new(void);

extern void image_data_arrays_free(ImageDataArrays *ida);

extern int image_create_dp_bad(struct image *image,
                               const DataTemplate *dtempl);

extern int image_set_zero_data(struct image *image,
                               const DataTemplate *dtempl);

#ifdef __cplusplus
}
#endif

#endif	/* IMAGE_H */
