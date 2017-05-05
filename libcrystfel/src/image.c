/*
 * image.c
 *
 * Handle images and image features
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011-2017 Thomas White <taw@physics.org>
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
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <hdf5.h>
#include <cbflib/cbf.h>

#include "image.h"
#include "utils.h"
#include "events.h"
#include "hdf5-file.h"
#include "detector.h"

/**
 * SECTION:image
 * @short_description: Data structure representing an image
 * @title: Image
 * @section_id:
 * @see_also:
 * @include: "image.h"
 * @Image:
 *
 * The <structname>image</structname> structure represents an image, usually one
 * frame from a large series of diffraction patterns, which might be from the
 * same or different crystals.
 */


struct imagefile
{
	enum imagefile_type type;
	char *filename;
	struct hdfile *hdfile;
};


struct _imagefeaturelist
{
	struct imagefeature	*features;
	int			n_features;
};


void image_add_feature(ImageFeatureList *flist, double fs, double ss,
                       struct panel *p,
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

	flist->features[flist->n_features].fs = fs;
	flist->features[flist->n_features].ss = ss;
	flist->features[flist->n_features].p = p;
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


static int comp(const void *a, const void *b)
{
	const struct imagefeature *ap = a;
	const struct imagefeature *bp = b;

	return ap->intensity < bp->intensity;
}


/* Strongest first.  Returned list is guaranteed not to have any holes
 * (feature->valid = 0) */
ImageFeatureList *sort_peaks(ImageFeatureList *flist)
{
	ImageFeatureList *n = image_feature_list_new();
	int nf, i;

	if ( n == NULL ) return NULL;

	n->features = malloc(flist->n_features*sizeof(struct imagefeature));
	if ( n->features == NULL ) {
		free(n);
		return NULL;
	}

	nf = 0;
	for ( i=0; i<flist->n_features; i++ ) {
		struct imagefeature *f;
		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;
		n->features[nf++] = flist->features[i];
	}
	n->n_features = nf;

	qsort(n->features, nf, sizeof(struct imagefeature), comp);

	return n;
}


void image_feature_list_free(ImageFeatureList *flist)
{
	if ( !flist ) return;

	if ( flist->features ) free(flist->features);
	free(flist);
}


struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                           double fs, double ss,
                                           struct panel *p, double *d, int *idx)
{
	int i;
	double dmin = +HUGE_VAL;
	int closest = 0;

	for ( i=0; i<flist->n_features; i++ ) {

		double ds;

		if ( p != flist->features[i].p ) continue;

		ds = distance(flist->features[i].fs, flist->features[i].ss,
		              fs, ss);

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


Reflection *image_reflection_closest(RefList *rlist,
                                     double fs, double ss, struct panel *p,
                                     struct detector *det,
                                     double *d)
{

	double dmin = HUGE_VAL;
	Reflection *closest = NULL;
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(rlist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double ds;
		struct panel *p2;
		double rfs, rss;

		get_detector_pos(refl, &rfs, &rss);
		p2 = get_panel(refl);

		if ( p != p2 ) continue;

		ds = distance(rfs, rss, fs, ss);

		if ( ds < dmin ) {
			dmin = ds;
			closest = refl;
		}

	}

	if ( dmin < +HUGE_VAL ) {
		*d = dmin;
		return closest;
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
	if ( idx >= flist->n_features ) return NULL;

	if ( flist->features[idx].valid == 0 ) return NULL;

	return &flist->features[idx];
}


void image_remove_feature(ImageFeatureList *flist, int idx)
{
	flist->features[idx].valid = 0;
}


void image_add_crystal(struct image *image, Crystal *cryst)
{
	Crystal **crs;
	int n;

	n = image->n_crystals;
	crs = realloc(image->crystals, (n+1)*sizeof(Crystal *));
	if ( crs == NULL ) {
		ERROR("Failed to allocate memory for crystals.\n");
		return;
	}

	crs[n] = cryst;
	image->crystals = crs;
	image->n_crystals = n+1;
}


void remove_flagged_crystals(struct image *image)
{
	int i;

	for ( i=0; i<image->n_crystals; i++ ) {
		if ( crystal_get_user_flag(image->crystals[i]) ) {
			int j;
			Crystal *deleteme = image->crystals[i];
			cell_free(crystal_get_cell(deleteme));
			crystal_free(deleteme);
			for ( j=i; j<image->n_crystals-1; j++ ) {
				image->crystals[j] = image->crystals[j+1];
			}
			image->n_crystals--;
		}
	}

}


/* Free all crystals, including their RefLists and UnitCells */
void free_all_crystals(struct image *image)
{
	int i;
	if ( image->crystals == NULL ) return;
	for ( i=0; i<image->n_crystals; i++ ) {
		Crystal *cr = image->crystals[i];
		reflist_free(crystal_get_reflections(cr));
		cell_free(crystal_get_cell(cr));
		crystal_free(image->crystals[i]);
	}
	free(image->crystals);
	image->n_crystals = 0;
}


/**************************** Image field lists *******************************/

struct imagefile_field_list
{
	char **fields;
	int n_fields;
	int max_fields;
};


struct imagefile_field_list *new_imagefile_field_list()
{
	struct imagefile_field_list *n;

	n = calloc(1, sizeof(struct imagefile_field_list));
	if ( n == NULL ) return NULL;

	n->max_fields = 32;
	n->fields = malloc(n->max_fields*sizeof(char *));
	if ( n->fields == NULL ) {
		free(n);
		return NULL;
	}

	return n;
}


void free_imagefile_field_list(struct imagefile_field_list *n)
{
	int i;
	for ( i=0; i<n->n_fields; i++ ) {
		free(n->fields[i]);
	}
	free(n->fields);
	free(n);
}


void add_imagefile_field(struct imagefile_field_list *copyme, const char *name)
{
	int i;

	/* Already on the list?   Don't re-add if so. */
	for ( i=0; i<copyme->n_fields; i++ ) {
		if ( strcmp(copyme->fields[i], name) == 0 ) return;
	}

	/* Need more space? */
	if ( copyme->n_fields == copyme->max_fields ) {

		char **nfields;
		int nmax = copyme->max_fields + 32;

		nfields = realloc(copyme->fields, nmax*sizeof(char *));
		if ( nfields == NULL ) {
			ERROR("Failed to allocate space for new HDF5 field.\n");
			return;
		}

		copyme->max_fields = nmax;
		copyme->fields = nfields;

	}

	copyme->fields[copyme->n_fields] = strdup(name);
	if ( copyme->fields[copyme->n_fields] == NULL ) {
		ERROR("Failed to add field for copying '%s'\n", name);
		return;
	}

	copyme->n_fields++;
}


/******************************* CBF files ************************************/


static char *cbf_strerr(int e)
{
	char *err;

	err = malloc(1024);
	if ( err == NULL ) return NULL;

	err[0] = '\0';

	/* NB Sum of lengths of all strings must be less than 1024 */
	if ( e & CBF_FORMAT ) strcat(err, "Invalid format");
	if ( e & CBF_ALLOC ) strcat(err, "Memory allocation failed");
	if ( e & CBF_ARGUMENT ) strcat(err, "Invalid argument");
	if ( e & CBF_ASCII ) strcat(err, "Value is ASCII");
	if ( e & CBF_BINARY ) strcat(err, "Value is binary");
	if ( e & CBF_BITCOUNT ) strcat(err, "Wrong number of bits");
	if ( e & CBF_ENDOFDATA ) strcat(err, "End of data");
	if ( e & CBF_FILECLOSE ) strcat(err, "File close error");
	if ( e & CBF_FILEOPEN ) strcat(err, "File open error");
	if ( e & CBF_FILEREAD ) strcat(err, "File read error");
	if ( e & CBF_FILETELL ) strcat(err, "File tell error");
	if ( e & CBF_FILEWRITE ) strcat(err, "File write error");
	if ( e & CBF_IDENTICAL ) strcat(err, "Name already exists");
	if ( e & CBF_NOTFOUND ) strcat(err, "Not found");
	if ( e & CBF_OVERFLOW ) strcat(err, "Overflow");
	if ( e & CBF_UNDEFINED ) strcat(err, "Number undefined");
	if ( e & CBF_NOTIMPLEMENTED ) strcat(err, "Not implemented");

	return err;
}


static int unpack_panels(struct image *image, signed int *data, int data_width)
{
	int pi;

	/* FIXME: Load these masks from an HDF5 file, if filenames are
	 * given in the geometry file */
	uint16_t *flags = NULL;
	float *sat = NULL;

	image->dp = malloc(image->det->n_panels * sizeof(float *));
	image->bad = malloc(image->det->n_panels * sizeof(int *));
	image->sat = malloc(image->det->n_panels * sizeof(float *));
	if ( (image->dp == NULL) || (image->bad == NULL)
	  || (image->sat == NULL) )
	{
		ERROR("Failed to allocate panels.\n");
		return 1;
	}

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		struct panel *p;
		int fs, ss;

		p = &image->det->panels[pi];
		image->dp[pi] = malloc(p->w*p->h*sizeof(float));
		image->bad[pi] = calloc(p->w*p->h, sizeof(int));
		image->sat[pi] = malloc(p->w*p->h*sizeof(float));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL)
		  || (image->sat[pi] == NULL) )
		{
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		if ( p->mask != NULL ) {
			ERROR("WARNING: Bad pixel masks do not currently work "
			      "with CBF files\n");
			ERROR(" (bad pixel regions specified in the geometry "
			      "file will be used, however)\n");
		}

		if ( p->satmap != NULL ) {
			ERROR("WARNING: Saturation maps do not currently work "
			      "with CBF files\n");
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			int cfs, css;
			int bad = 0;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p->w*ss] = data[idx];

			if ( sat != NULL ) {
				image->sat[pi][fs+p->w*ss] = sat[idx];
			} else {
				image->sat[pi][fs+p->w*ss] = INFINITY;
			}

			if ( p->no_index ) bad = 1;

			if ( in_bad_region(image->det, p, cfs, css) ) {
				bad = 1;
			}

			if ( flags != NULL ) {

				int f;

				f = flags[idx];

				/* Bad if it's missing any of the "good" bits */
				if ( (f & image->det->mask_good)
				       != image->det->mask_good ) bad = 1;

				/* Bad if it has any of the "bad" bits. */
				if ( f & image->det->mask_bad ) bad = 1;

			}
			image->bad[pi][fs+p->w*ss] = bad;

		}
		}

	}

	return 0;
}


static void cbf_fill_in_beam_parameters(struct beam_params *beam,
                                        struct imagefile *f,
                                        struct image *image)
{
	double eV;

	if ( beam->photon_energy_from == NULL ) {

		/* Explicit value given */
		eV = beam->photon_energy;

	} else {

		ERROR("Can't get photon energy from CBF yet.\n");
		eV = 0.0;

	}

	image->lambda = ph_en_to_lambda(eV_to_J(eV))*beam->photon_energy_scale;
}


static void cbf_fill_in_clen(struct detector *det, struct imagefile *f)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];

		if ( p->clen_from != NULL ) {

			ERROR("Can't get clen from CBF yet.\n");

		}

		adjust_centering_for_rail(p);

	}
}


static int read_cbf(struct imagefile *f, struct image *image)
{
	cbf_handle cbfh;
	FILE *fh;
	int r;
	unsigned int compression;
	int binary_id, minelement, maxelement, elsigned, elunsigned;
	size_t elsize, elements, elread, dimfast, dimmid, dimslow, padding;
	const char *byteorder;
	signed int *data;

	if ( cbf_make_handle(&cbfh) ) {
		ERROR("Failed to allocate CBF handle\n");
		return 1;
	}

	fh = fopen(f->filename, "rb");
	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", f->filename);
		return 1;
	}
	/* CBFlib calls fclose(fh) when it's ready */

	if ( cbf_read_widefile(cbfh, fh, 0) ) {
		ERROR("Failed to read CBF file\n");
		cbf_free_handle(cbfh);
		return 1;
	}

	/* Select row 0 in data column inside array_data */
	cbf_find_category(cbfh, "array_data");
	cbf_find_column(cbfh, "data");
	cbf_select_row(cbfh, 0);

	/* Get parameters for array read */
	r = cbf_get_integerarrayparameters_wdims(cbfh, &compression, &binary_id,
	                                         &elsize, &elsigned, &elunsigned,
	                                         &elements,
	                                         &minelement, &maxelement,
	                                         &byteorder,
	                                         &dimfast, &dimmid, &dimslow,
	                                         &padding);
	if ( r ) {
		char *err = cbf_strerr(r);
		ERROR("Failed to read CBF array parameters: %s\n", err);
		free(err);
		cbf_free_handle(cbfh);
		return 1;
	}

	if ( dimslow != 0 ) {
		ERROR("CBF data array is 3D - don't know what to do with it\n");
		cbf_free_handle(cbfh);
		return 1;
	}

	if ( dimfast*dimmid*elsize > 10e9 ) {
		ERROR("CBF data is far too big (%i x %i x %i bytes).\n",
		      (int)dimfast, (int)dimmid, (int)elsize);
		cbf_free_handle(cbfh);
		return 1;
	}

	if ( elsize != 4 ) {
		STATUS("Don't know what to do with element size %i\n",
		       (int)elsize);
		cbf_free_handle(cbfh);
		return 1;
	}

	if ( strcmp(byteorder, "little_endian") != 0 ) {
		STATUS("Don't know what to do with non-little-endian datan\n");
		cbf_free_handle(cbfh);
		return 1;
	}

	data = malloc(elsize*dimfast*dimmid);
	if ( data == NULL ) {
		ERROR("Failed to allocate memory for CBF data\n");
		cbf_free_handle(cbfh);
		return 1;
	}

	r = cbf_get_integerarray(cbfh, &binary_id, data, elsize, 1,
	                         elsize*dimfast*dimmid, &elread);
	if ( r ) {
		char *err = cbf_strerr(r);
		ERROR("Failed to read CBF array: %s\n", err);
		free(err);
		cbf_free_handle(cbfh);
		return 1;
	}

	unpack_panels(image, data, dimfast);
	free(data);

	if ( image->beam != NULL ) {
		cbf_fill_in_beam_parameters(image->beam, f, image);
		if ( image->lambda > 1000 ) {
			ERROR("WARNING: Missing or nonsensical wavelength "
			      "(%e m) for %s.\n",
			      image->lambda, image->filename);
		}
	}
	cbf_fill_in_clen(image->det, f);
	fill_in_adu(image);

	cbf_free_handle(cbfh);
	return 0;
}


/****************************** Image files ***********************************/


static signed int is_cbf_file(const char *filename)
{
	FILE *fh;
	char line[1024];

	fh = fopen(filename, "r");
	if ( fh == NULL ) return -1;

	if ( fgets(line, 1024, fh) == NULL ) return -1;
	fclose(fh);

	if ( strstr(line, "CBF") == NULL ) {
		return 0;
	}

	return 1;
}


struct imagefile *imagefile_open(const char *filename)
{
	struct imagefile *f;

	f = malloc(sizeof(struct imagefile));
	if ( f == NULL ) return NULL;

	if ( H5Fis_hdf5(filename) > 0 ) {

		/* This is an HDF5, pass through to HDF5 layer */
		f->type = IMAGEFILE_HDF5;
		f->hdfile = hdfile_open(filename);

		if ( f->hdfile == NULL ) {
			free(f);
			return NULL;
		}

	} else if ( is_cbf_file(filename) > 0 ) {

		f->type = IMAGEFILE_CBF;

	}

	f->filename = strdup(filename);
	return f;
}


int imagefile_read(struct imagefile *f, struct image *image,
                   struct event *event)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		return hdf5_read2(f->hdfile, image, event, 0);
	} else if ( f->type == IMAGEFILE_CBF ) {
		return read_cbf(f, image);
	} else {
		ERROR("Unknown file type %i\n", f->type);
		return 1;
	}
}


/* Read a simple file, no multi-event, no prior geometry etc, and
 * generate a geometry for it */
int imagefile_read_simple(struct imagefile *f, struct image *image)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		return hdf5_read(f->hdfile, image, NULL, 0);
	} else {
		STATUS("Mock CBF simple read\n");
		return 0;
	}
}


enum imagefile_type imagefile_get_type(struct imagefile *f)
{
	assert(f != NULL);
	return f->type;
}


struct hdfile *imagefile_get_hdfile(struct imagefile *f)
{
	if ( f->type != IMAGEFILE_HDF5 ) {
		ERROR("Not an HDF5 file!\n");
		return NULL;
	}

	return f->hdfile;
}


void imagefile_copy_fields(struct imagefile *f,
                           const struct imagefile_field_list *copyme,
                           FILE *fh, struct event *ev)
{
	int i;

	if ( copyme == NULL ) return;

	for ( i=0; i<copyme->n_fields; i++ ) {

		char *val;
		char *field;

		field = copyme->fields[i];

		if ( f->type == IMAGEFILE_HDF5 ) {
			val = hdfile_get_string_value(f->hdfile, field, ev);
			if ( field[0] == '/' ) {
				fprintf(fh, "hdf5%s = %s\n", field, val);
			} else {
				fprintf(fh, "hdf5/%s = %s\n", field, val);
			}
			free(val);

		} else {
			STATUS("Mock CBF variable\n");
			fprintf(fh, "cbf/%s = %s\n", field, "(FIXME)");
		}


	}
}


void imagefile_close(struct imagefile *f)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		hdfile_close(f->hdfile);
	}
	free(f->filename);
	free(f);
}
