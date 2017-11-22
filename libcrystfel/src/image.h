/*
 * image.h
 *
 * Handle images and image features
 *
 * Copyright © 2012-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2018 Thomas White <taw@physics.org>
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

struct detector;

#include <stdint.h>
#include <complex.h>
#include <sys/types.h>

struct imagefeature;
struct sample;
struct image;
struct imagefile;
struct imagefile_field_list;

#include "utils.h"
#include "cell.h"
#include "detector.h"
#include "reflist.h"
#include "crystal.h"
#include "index.h"
#include "events.h"

/**
 * SpectrumType:
 * @SPECTRUM_TOPHAT: A top hat distribution of wavelengths
 * @SPECTRUM_SASE: A simulated SASE spectrum
 * @SPECTRUM_TWOCOLOUR: A spectrum containing two peaks
 * @SPECTRUM_FROMFILE: An arbitrary spectrum read from input file
 *
 * A %SpectrumType represents a type of X-ray energy spectrum to use for
 * generating simulated data.
 **/
typedef enum {
	SPECTRUM_TOPHAT,
	SPECTRUM_SASE,
	SPECTRUM_TWOCOLOUR,
	SPECTRUM_FROMFILE
} SpectrumType;


/**
 * imagefeature:
 *  @parent: Image this feature belongs to
 *  @fs: Fast scan coordinate
 *  @ss: Slow scan coordinate
 *  @p: Pointer to panel
 *  @intensity: Intensity of peak
 *  @rx: Reciprocal x coordinate in m^-1
 *  @ry: Reciprocal y coordinate in m^-1
 *  @rz: Reciprocal z coordinate in m^-1
 *  @name: Text name for feature
 *
 *  Represents a peak in an image.
 *
 *  Note carefully that the @fs and @ss coordinates are the distances, measured
 *  in pixels, from the corner of the panel.  They are NOT pixel indices.
 *  If the peak is in the middle of the first pixel, its coordinates would be
 *  0.5,0.5.
 */
struct imagefeature {

	struct image                    *parent;
	double                          fs;
	double                          ss;
	struct panel                    *p;
	double                          intensity;

	/* Reciprocal space coordinates (m^-1 of course) of this feature */
	double                          rx;
	double                          ry;
	double                          rz;

	const char                      *name;

	/*< private >*/
	int                             valid;
};


/* An enum representing the image file formats we can handle */
enum imagefile_type
{
	IMAGEFILE_HDF5,
	IMAGEFILE_CBF
};


/* An opaque type representing a list of image features */
typedef struct _imagefeaturelist ImageFeatureList;


struct spectrum
{
	int n;
	double *ks; /* 1/m */
	double *weights;
};


/* Structure describing a wavelength sample from a spectrum */
struct sample
{
	double k; /* 1/m */
	double weight;
};


/**
 * beam_params:
 *  @photon_energy: eV per photon
 *  @photon_energy_from: HDF5 dataset name
 *  @photon_energy_scale: Scale factor for photon energy, if it comes from HDF5
 */
struct beam_params
{
	double photon_energy;
	char  *photon_energy_from;
	double photon_energy_scale;
};


/**
 * image:
 *   @crystals: Array of crystals in the image
 *   @n_crystals: The number of crystals in the image
 *   @indexed_by: Indexing method which indexed this pattern
 *   @n_indexing_tries: Number of times the indexer was tried before indexing
 *   @det: Detector structure
 *   @beam: Beam parameters structure
 *   @filename: Filename for the image file
 *   @copyme: Fields to copy from the image file to the stream
 *   @id: ID number of the thread handling this image
 *   @serial: Serial number for this image
 *   @lambda: Wavelength
 *   @div: Divergence
 *   @bw: Bandwidth
 *   @num_peaks: The number of peaks
 *   @num_saturated_peaks: The number of saturated peaks
 *   @features: The peaks found in the image
 *   @dp: The image data, by panel
 *   @bad: The bad pixel mask, array by panel
 *   @sat: The per-pixel saturation mask, array by panel
 *   @event: Event ID for the image
 *   @stuff_from_stream: Items read back from the stream
 *   @avg_clen: Mean of camera length values for all panels
 *   @spectrum: Spectrum information
 *   @nsamples: Number of spectrum samples
 *   @spectrum_size: Size of spectrum array
 *
 * The field <structfield>data</structfield> contains the raw image data, if it
 * is currently available.  The data might be available throughout the
 * processing of an experimental pattern, but it might not be available when
 * simulating, scaling or merging patterns.
 *
 * <structfield>crystals</structfield> is an array of %Crystal directly
 * returned by the low-level indexing system. <structfield>n_crystals</structfield>
 * is the number of crystals which were found in the image.
 *
 * <structfield>copyme</structfield> represents a list of fields in the image
 * file (e.g. HDF5 fields or CBF headers) to copy to the output stream.
 **/
struct image;

struct image {

	float                   **dp;    /* Data in panel */
	int                     **bad;   /* Bad pixels by panel */
	float                   **sat;   /* Per-pixel saturation values */

	Crystal                 **crystals;
	int                     n_crystals;
	IndexingMethod          indexed_by;
	int                     n_indexing_tries;

	struct detector         *det;
	struct beam_params      *beam;  /* The nominal beam parameters */
	char                    *filename;
	struct event            *event;
	const struct imagefile_field_list *copyme;
	struct stuff_from_stream *stuff_from_stream;

	double                  avg_clen;  /* Average camera length extracted
	                                    * from stuff_from_stream */

	int                     id;   /* ID number of the thread
	                               * handling this image */
	int                     serial;  /* Monotonically ascending serial
	                                  * number for this image */

	struct spectrum *spectrum; /* Beam spectrum for pink beam data */

	// These only used in pattern_sim, to be changed to struct spectrum from above later...
	struct sample           *spectrum0;
	int                     nsamples; /* Number of wavelengths */
	int                     spectrum_size;  /* Size of "spectrum" */

	/* Per-shot radiation values */
	double                  lambda;        /* Wavelength in m */
	double                  div;           /* Divergence in radians */
	double                  bw;            /* FWHM bandwidth as a fraction */

	/* Detected peaks */
	long long               num_peaks;
	long long               num_saturated_peaks;
	double                  peak_resolution;  /* Estimate of resolution
	                                           * based on peaks only */
	ImageFeatureList        *features;

};

#ifdef __cplusplus
extern "C" {
#endif

/* Feature lists */
extern ImageFeatureList *image_feature_list_new(void);

extern void image_feature_list_free(ImageFeatureList *flist);

extern void image_add_feature(ImageFeatureList *flist, double x, double y,
                              struct panel *p,
                              struct image *parent, double intensity,
                              const char *name);

extern void image_remove_feature(ImageFeatureList *flist, int idx);

extern struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                                  double fs, double ss,
                                                  struct panel *p,
                                                  double *d, int *idx);

extern Reflection *image_reflection_closest(RefList *rlist,
                                            double fs, double ss,
                                            struct panel *p,
                                            struct detector *det,
                                            double *d);

extern int image_feature_count(ImageFeatureList *flist);
extern struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx);
extern ImageFeatureList *sort_peaks(ImageFeatureList *flist);

extern void image_add_crystal(struct image *image, Crystal *cryst);
extern int remove_flagged_crystals(struct image *image);
extern void free_all_crystals(struct image *image);

/* Image files */
extern struct imagefile *imagefile_open(const char *filename);
extern int imagefile_read(struct imagefile *f, struct image *image,
                          struct event *event);
extern int imagefile_read_simple(struct imagefile *f, struct image *image);
extern struct hdfile *imagefile_get_hdfile(struct imagefile *f);
extern enum imagefile_type imagefile_get_type(struct imagefile *f);
extern void imagefile_copy_fields(struct imagefile *f,
                                  const struct imagefile_field_list *copyme,
                                  FILE *fh, struct event *ev);
extern void imagefile_close(struct imagefile *f);
extern signed int is_cbf_file(const char *filename);

/* Field lists */
extern struct imagefile_field_list *new_imagefile_field_list(void);
extern void free_imagefile_field_list(struct imagefile_field_list *f);

extern void add_imagefile_field(struct imagefile_field_list *copyme,
                                const char *name);

#ifdef __cplusplus
}
#endif

#endif	/* IMAGE_H */
