/*
 * hdf5-file.c
 *
 * Read/write HDF5 data files
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <hdf5.h>
#include <assert.h>

#include "image.h"
#include "hdf5-file.h"
#include "utils.h"


struct hdfile {

	const char      *path;  /* Current data path */

	size_t          nx;  /* Image width */
	size_t          ny;  /* Image height */

	hid_t           fh;  /* HDF file handle */
	hid_t           dh;  /* Dataset handle */

	int             data_open;  /* True if dh is initialised */
};


struct hdfile *hdfile_open(const char *filename)
{
	struct hdfile *f;

	f = malloc(sizeof(struct hdfile));
	if ( f == NULL ) return NULL;

	/* Please stop spamming my terminal */
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

	f->fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( f->fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		free(f);
		return NULL;
	}

	f->data_open = 0;

	return f;
}


int hdfile_set_image(struct hdfile *f, const char *path)
{
	hsize_t size[2];
	hsize_t max_size[2];
	hid_t sh;

	f->dh = H5Dopen2(f->fh, path, H5P_DEFAULT);
	if ( f->dh < 0 ) {
		ERROR("Couldn't open dataset\n");
		return -1;
	}
	f->data_open = 1;

	sh = H5Dget_space(f->dh);
	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Dataset is not two-dimensional\n");
		return -1;
	}
	H5Sget_simple_extent_dims(sh, size, max_size);
	H5Sclose(sh);

	f->nx = size[0];
	f->ny = size[1];

	return 0;
}


int get_peaks(struct image *image, struct hdfile *f, const char *p)
{
	hid_t dh, sh;
	hsize_t size[2];
	hsize_t max_size[2];
	int i;
	float *buf;
	herr_t r;
	int tw;

	dh = H5Dopen2(f->fh, p, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Peak list (%s) not found.\n", p);
		return 1;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		H5Dclose(dh);
		ERROR("Couldn't get dataspace for peak list.\n");
		return 1;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Peak list has the wrong dimensionality (%i).\n",
		      H5Sget_simple_extent_ndims(sh));
		H5Sclose(sh);
		H5Dclose(dh);
		return 1;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	tw = size[1];
	if ( (tw != 3) && (tw != 4) ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Peak list has the wrong dimensions.\n");
		return 1;
	}

	buf = malloc(sizeof(float)*size[0]*size[1]);
	if ( buf == NULL ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Couldn't reserve memory for the peak list.\n");
		return 1;
	}
	r = H5Dread(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read peak list.\n");
		free(buf);
		return 1;
	}

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( i=0; i<size[0]; i++ ) {

		float fs, ss, val;
		struct panel *p;

		fs = buf[tw*i+0];
		ss = buf[tw*i+1];
		val = buf[tw*i+2];

		p = find_panel(image->det, fs, ss);
		if ( p == NULL ) continue;
		if ( p->no_index ) continue;

		image_add_feature(image->features, fs, ss, image, val, NULL);

	}

	free(buf);
	H5Sclose(sh);
	H5Dclose(dh);

	return 0;
}


static void cleanup(hid_t fh)
{
	int n_ids, i;
	hid_t ids[256];

	n_ids = H5Fget_obj_ids(fh, H5F_OBJ_ALL, 256, ids);
	for ( i=0; i<n_ids; i++ ) {

		hid_t id;
		H5I_type_t type;

		id = ids[i];
		type = H5Iget_type(id);

		if ( type == H5I_GROUP ) H5Gclose(id);
		if ( type == H5I_DATASET ) H5Dclose(id);
		if ( type == H5I_DATATYPE ) H5Tclose(id);
		if ( type == H5I_DATASPACE ) H5Sclose(id);
		if ( type == H5I_ATTR ) H5Aclose(id);

	}
}


void hdfile_close(struct hdfile *f)
{
	if ( f->data_open ) {
		H5Dclose(f->dh);
	}
	cleanup(f->fh);

	H5Fclose(f->fh);
	free(f);
}


/* Deprecated */
int hdf5_write(const char *filename, const void *data,
               int width, int height, int type)
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	hid_t ph;  /* Property list */
	herr_t r;
	hsize_t size[2];

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return 1;
	}

	gh = H5Gcreate2(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
		return 1;
	}

	/* Note the "swap" here, according to section 3.2.5,
	 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
	size[0] = height;
	size[1] = width;
	sh = H5Screate_simple(2, size, NULL);

	/* Set compression */
	ph = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(ph, 2, size);
	H5Pset_deflate(ph, 3);

	dh = H5Dcreate2(gh, "data", type, sh,
	                H5P_DEFAULT, ph, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return 1;
	}

	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, NULL);

	r = H5Dwrite(dh, type, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}
	H5Dclose(dh);
	H5Gclose(gh);
	H5Pclose(ph);
	H5Fclose(fh);

	return 0;
}


int hdf5_write_image(const char *filename, struct image *image)
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	hid_t ph;  /* Property list */
	herr_t r;
	hsize_t size[2];
	double lambda, eV;
	double *arr;
	int i;

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return 1;
	}

	gh = H5Gcreate2(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
		return 1;
	}

	/* Note the "swap" here, according to section 3.2.5,
	 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
	size[0] = image->height;
	size[1] = image->width;
	sh = H5Screate_simple(2, size, NULL);

	/* Set compression */
	ph = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(ph, 2, size);
	H5Pset_deflate(ph, 3);

	dh = H5Dcreate2(gh, "data", H5T_NATIVE_FLOAT, sh,
	                H5P_DEFAULT, ph, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return 1;
	}

	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, NULL);

	r = H5Dwrite(dh, H5T_NATIVE_FLOAT, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, image->data);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}
	H5Dclose(dh);

	H5Gclose(gh);

	gh = H5Gcreate2(fh, "LCLS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		printf("Couldn't create group\n");
		H5Fclose(fh);
		return 1;
	}

	size[0] = 1;
	sh = H5Screate_simple(1, size, NULL);

	dh = H5Dcreate2(gh, "photon_energy_eV", H5T_NATIVE_DOUBLE, sh,
	                H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		H5Fclose(fh);
		return 1;
	}
	eV = ph_lambda_to_eV(image->lambda);
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, &eV);

	H5Dclose(dh);

	dh = H5Dcreate2(fh, "/LCLS/photon_wavelength_A", H5T_NATIVE_DOUBLE, sh,
	                H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		H5Fclose(fh);
		return 1;
	}
	lambda = image->lambda * 1e10;
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, &lambda);
	if ( r < 0 ) {
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}

	H5Dclose(dh);

	if ( image->spectrum_size > 0 ) {

		arr = malloc(image->spectrum_size*sizeof(double));
		if ( arr == NULL ) {
			H5Fclose(fh);
			return 1;
		}
		for ( i=0; i<image->spectrum_size; i++ ) {
			arr[i] = 1.0e10/image->spectrum[i].k;
		}

		size[0] = image->spectrum_size;
		sh = H5Screate_simple(1, size, NULL);

		dh = H5Dcreate2(gh, "spectrum_wavelengths_A", H5T_NATIVE_DOUBLE,
		                sh, H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
		if ( dh < 0 ) {
			H5Fclose(fh);
			return 1;
		}
		r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
			     H5S_ALL, H5P_DEFAULT, arr);
		H5Dclose(dh);

		for ( i=0; i<image->spectrum_size; i++ ) {
			arr[i] = image->spectrum[i].weight;
		}
		dh = H5Dcreate2(gh, "spectrum_weights", H5T_NATIVE_DOUBLE, sh,
			        H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
		if ( dh < 0 ) {
			H5Fclose(fh);
			return 1;
		}
		r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
			     H5S_ALL, H5P_DEFAULT, arr);

		H5Dclose(dh);
		free(arr);

		size[0] = 1;
		sh = H5Screate_simple(1, size, NULL);

		dh = H5Dcreate2(gh, "number_of_samples", H5T_NATIVE_INT, sh,
			        H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
		if ( dh < 0 ) {
			H5Fclose(fh);
			return 1;
		}

		r = H5Dwrite(dh, H5T_NATIVE_INT, H5S_ALL,
			     H5S_ALL, H5P_DEFAULT, &image->nsamples);

		H5Dclose(dh);
	}

	H5Gclose(gh);

	H5Pclose(ph);

	H5Fclose(fh);

	return 0;
}


static void debodge_saturation(struct hdfile *f, struct image *image)
{
	hid_t dh, sh;
	hsize_t size[2];
	hsize_t max_size[2];
	int i;
	float *buf;
	herr_t r;

	dh = H5Dopen2(f->fh, "/processing/hitfinder/peakinfo_saturated",
	              H5P_DEFAULT);

	if ( dh < 0 ) {
		/* This isn't an error */
		return;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		H5Dclose(dh);
		ERROR("Couldn't get dataspace for saturation table.\n");
		return;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		H5Sclose(sh);
		H5Dclose(dh);
		return;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	if ( size[1] != 3 ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Saturation table has the wrong dimensions.\n");
		return;
	}

	buf = malloc(sizeof(float)*size[0]*size[1]);
	if ( buf == NULL ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Couldn't reserve memory for saturation table.\n");
		return;
	}
	r = H5Dread(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read saturation table.\n");
		free(buf);
		return;
	}

	for ( i=0; i<size[0]; i++ ) {

		unsigned int x, y;
		float val;

		x = buf[3*i+0];
		y = buf[3*i+1];
		val = buf[3*i+2];

		image->data[x+image->width*y] = val / 5.0;
		image->data[x+1+image->width*y] = val / 5.0;
		image->data[x-1+image->width*y] = val / 5.0;
		image->data[x+image->width*(y+1)] = val / 5.0;
		image->data[x+image->width*(y-1)] = val / 5.0;

	}

	free(buf);
	H5Sclose(sh);
	H5Dclose(dh);
}


static int unpack_panels(struct image *image, struct detector *det)
{
	int pi;

	image->dp = malloc(det->n_panels * sizeof(float *));
	image->bad = malloc(det->n_panels * sizeof(int *));
	if ( (image->dp == NULL) || (image->bad == NULL) ) {
		ERROR("Failed to allocate panels.\n");
		return 1;
	}

	for ( pi=0; pi<det->n_panels; pi++ ) {

		struct panel *p;
		int fs, ss;

		p = &det->panels[pi];
		image->dp[pi] = malloc(p->w*p->h*sizeof(float));
		image->bad[pi] = calloc(p->w*p->h, sizeof(int));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL) ) {
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			int cfs, css;
			int bad = 0;

			cfs = fs+p->min_fs;
			css = ss+p->min_ss;
			idx = cfs + css*image->width;

			image->dp[pi][fs+p->w*ss] = image->data[idx];

			if ( p->no_index ) bad = 1;

			if ( in_bad_region(det, cfs, css) ) {
				bad = 1;
			}

			if ( image->flags != NULL ) {

				int flags;

				flags = image->flags[idx];

				/* Bad if it's missing any of the "good" bits */
				if ( !((flags & image->det->mask_good)
			                   == image->det->mask_good) ) bad = 1;

				/* Bad if it has any of the "bad" bits. */
				if ( flags & image->det->mask_bad ) bad = 1;

			}

			image->bad[pi][fs+p->w*ss] = bad;

		}
		}

	}

	return 0;
}


int hdf5_read(struct hdfile *f, struct image *image, int satcorr)
{
	herr_t r;
	float *buf;
	uint16_t *flags;
	hid_t mask_dh;

	/* Note the "swap" here, according to section 3.2.5,
	 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
	image->width = f->ny;
	image->height = f->nx;

	buf = malloc(sizeof(float)*f->nx*f->ny);

	r = H5Dread(f->dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read data\n");
		free(buf);
		return 1;
	}
	image->data = buf;

	if ( (image->det != NULL) && (image->det->mask != NULL) ) {

		mask_dh = H5Dopen2(f->fh, image->det->mask, H5P_DEFAULT);
		if ( mask_dh <= 0 ) {
			ERROR("Couldn't open flags\n");
			image->flags = NULL;
		} else {
			flags = malloc(sizeof(uint16_t)*f->nx*f->ny);
			r = H5Dread(mask_dh, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL,
				    H5P_DEFAULT, flags);
			if ( r < 0 ) {
				ERROR("Couldn't read flags\n");
				free(flags);
				image->flags = NULL;
			} else {
				image->flags = flags;
			}
			H5Dclose(mask_dh);
		}

	}

	if ( satcorr ) debodge_saturation(f, image);

	if ( image->det != NULL ) {

		if ( (image->width != image->det->max_fs + 1 )
		  || (image->height != image->det->max_ss + 1))
		{
			ERROR("Image size doesn't match geometry size"
				" - rejecting image.\n");
			ERROR("Image size: %i,%i.  Geometry size: %i,%i\n",
			      image->width, image->height,
			      image->det->max_fs + 1, image->det->max_ss + 1);
			return 1;
		}

		fill_in_values(image->det, f);

		unpack_panels(image, image->det);

	}

	if ( image->beam != NULL ) {

		fill_in_beam_parameters(image->beam, f);
		image->lambda = ph_en_to_lambda(eV_to_J(image->beam->photon_energy));

		if ( (image->beam->photon_energy < 0.0)
		  || (image->lambda > 1000) ) {
			/* Error message covers a silly value in the beam file
			 * or in the HDF5 file. */
			ERROR("Nonsensical wavelength (%e m or %e eV) value "
			      "for %s.\n",
			      image->lambda, image->beam->photon_energy,
			      image->filename);
			return 1;
		}

	}

	return 0;
}


static int looks_like_image(hid_t h)
{
	hid_t sh;
	hsize_t size[2];
	hsize_t max_size[2];

	sh = H5Dget_space(h);
	if ( sh < 0 ) return 0;

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		return 0;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	if ( ( size[0] > 64 ) && ( size[1] > 64 ) ) return 1;

	return 0;
}


int hdfile_is_scalar(struct hdfile *f, const char *name, int verbose)
{
	hid_t dh;
	hid_t sh;
	hsize_t size[3];
	hid_t type;
	int ndims;
	int i;

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("No such field '%s'\n", name);
		return 0;
	}

	type = H5Dget_type(dh);

	/* Get the dimensionality.  We have to cope with scalars expressed as
	 * arrays with all dimensions 1, as well as zero-d arrays. */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > 3 ) {
		if ( verbose ) {
			ERROR("Too many dimensions (%i).\n", ndims);
		}
		H5Tclose(type);
		H5Dclose(dh);
		return 0;
	}

	/* Check that the size in all dimensions is 1 */
	H5Sget_simple_extent_dims(sh, size, NULL);
	for ( i=0; i<ndims; i++ ) {
		if ( size[i] != 1 ) {
			if ( verbose ) {
				ERROR("%s not a scalar value (ndims=%i,"
				      "size[%i]=%i)\n",
				      name, ndims, i, (int)size[i]);
			}
			H5Tclose(type);
			H5Dclose(dh);
			return 0;
		}
	}

	H5Tclose(type);
	H5Dclose(dh);

	return 1;
}


static int get_f_value(struct hdfile *f, const char *name, double *val)
{
	hid_t dh;
	hid_t type;
	hid_t class;
	herr_t r;
	double buf;

	if ( !hdfile_is_scalar(f, name, 1) ) return 1;

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("No such field '%s'\n", name);
		return 1;
	}

	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( class != H5T_FLOAT ) {
		ERROR("Not a floating point value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	r = H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, &buf);
	if ( r < 0 )  {
		ERROR("Couldn't read value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	*val = buf;
	return 0;
}


static int get_i_value(struct hdfile *f, const char *name, int *val)
{
	hid_t dh;
	hid_t type;
	hid_t class;
	herr_t r;
	int buf;

	if ( !hdfile_is_scalar(f, name, 1) ) return 1;

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("No such field '%s'\n", name);
		return 1;
	}

	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( class != H5T_INTEGER ) {
		ERROR("Not an integer value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	r = H5Dread(dh, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, &buf);
	if ( r < 0 )  {
		ERROR("Couldn't read value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	*val = buf;
	return 0;
}


double get_value(struct hdfile *f, const char *name)
{
	double val;
	get_f_value(f, name, &val);
	return val;
}


struct copy_hdf5_field
{
	char **fields;
	int n_fields;
	int max_fields;
};


struct copy_hdf5_field *new_copy_hdf5_field_list()
{
	struct copy_hdf5_field *n;

	n = calloc(1, sizeof(struct copy_hdf5_field));
	if ( n == NULL ) return NULL;

	n->max_fields = 32;
	n->fields = malloc(n->max_fields*sizeof(char *));
	if ( n->fields == NULL ) {
		free(n);
		return NULL;
	}

	return n;
}


void free_copy_hdf5_field_list(struct copy_hdf5_field *n)
{
	int i;
	for ( i=0; i<n->n_fields; i++ ) {
		free(n->fields[i]);
	}
	free(n->fields);
	free(n);
}


void add_copy_hdf5_field(struct copy_hdf5_field *copyme,
                         const char *name)
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


void copy_hdf5_fields(struct hdfile *f, const struct copy_hdf5_field *copyme,
                      FILE *fh)
{
	int i;

	if ( copyme == NULL ) return;

	for ( i=0; i<copyme->n_fields; i++ ) {

		char *val;
		char *field;

		field = copyme->fields[i];
		val = hdfile_get_string_value(f, field);

		if ( field[0] == '/' ) {
			fprintf(fh, "hdf5%s = %s\n", field, val);
		} else {
			fprintf(fh, "hdf5/%s = %s\n", field, val);
		}

		free(val);

	}
}


char *hdfile_get_string_value(struct hdfile *f, const char *name)
{
	hid_t dh;
	hsize_t size;
	hid_t type;
	hid_t class;
	int buf_i;
	double buf_f;
	char *tmp;

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
	if ( dh < 0 ) return NULL;

	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( class == H5T_STRING ) {

		herr_t r;
		char *tmp;
		hid_t sh;

		size = H5Tget_size(type);
		tmp = malloc(size+1);

		sh = H5Screate(H5S_SCALAR);

		r = H5Dread(dh, type, sh, sh, H5P_DEFAULT, tmp);
		if ( r < 0 ) goto fail;

		/* Two possibilities:
		 *   String is already zero-terminated
		 *   String is not terminated.
		 * Make sure things are done properly... */
		tmp[size] = '\0';
		chomp(tmp);

		return tmp;

	}

	switch ( class ) {

		case H5T_FLOAT :
		if ( get_f_value(f, name, &buf_f) ) goto fail;
		tmp = malloc(256);
		snprintf(tmp, 255, "%f", buf_f);
		return tmp;

		case H5T_INTEGER :
		if ( get_i_value(f, name, &buf_i) ) goto fail;
		tmp = malloc(256);
		snprintf(tmp, 255, "%d", buf_i);
		return tmp;

		default :
		goto fail;

	}

fail:
	H5Tclose(type);
	H5Dclose(dh);
	return NULL;
}


char **hdfile_read_group(struct hdfile *f, int *n, const char *parent,
                        int **p_is_group, int **p_is_image)
{
	hid_t gh;
	hsize_t num;
	char **res;
	int i;
	int *is_group;
	int *is_image;
	H5G_info_t ginfo;

	gh = H5Gopen2(f->fh, parent, H5P_DEFAULT);
	if ( gh < 0 ) {
		*n = 0;
		return NULL;
	}

	if ( H5Gget_info(gh, &ginfo) < 0 ) {
		/* Whoopsie */
		*n = 0;
		return NULL;
	}
	num = ginfo.nlinks;
	*n = num;
	if ( num == 0 ) return NULL;

	res = malloc(num*sizeof(char *));
	is_image = malloc(num*sizeof(int));
	is_group = malloc(num*sizeof(int));
	*p_is_image = is_image;
	*p_is_group = is_group;

	for ( i=0; i<num; i++ ) {

		char buf[256];
		hid_t dh;
		H5I_type_t type;

		H5Lget_name_by_idx(gh, ".", H5_INDEX_NAME, H5_ITER_NATIVE,
		                   i, buf, 255, H5P_DEFAULT);
		res[i] = malloc(256);
		if ( strlen(parent) > 1 ) {
			snprintf(res[i], 255, "%s/%s", parent, buf);
		} else {
			snprintf(res[i], 255, "%s%s", parent, buf);
		} /* ick */

		is_image[i] = 0;
		is_group[i] = 0;
		dh = H5Oopen(gh, buf, H5P_DEFAULT);
		if ( dh < 0 ) continue;
		type = H5Iget_type(dh);

		if ( type == H5I_GROUP ) {
			is_group[i] = 1;
		} else if ( type == H5I_DATASET ) {
			is_image[i] = looks_like_image(dh);
		}
		H5Oclose(dh);

	}

	return res;
}


int hdfile_set_first_image(struct hdfile *f, const char *group)
{
	char **names;
	int *is_group;
	int *is_image;
	int n, i, j;

	names = hdfile_read_group(f, &n, group, &is_group, &is_image);
	if ( n == 0 ) return 1;

	for ( i=0; i<n; i++ ) {

		if ( is_image[i] ) {
			hdfile_set_image(f, names[i]);
			for ( j=0; j<n; j++ ) free(names[j]);
			free(is_image);
			free(is_group);
			free(names);
			return 0;
		} else if ( is_group[i] ) {
			if ( !hdfile_set_first_image(f, names[i]) ) {
				for ( j=0; j<n; j++ ) free(names[j]);
				free(is_image);
				free(is_group);
				free(names);
				return 0;
			}
		}

	}

	for ( j=0; j<n; j++ ) free(names[j]);
	free(is_image);
	free(is_group);
	free(names);

	return 1;
}
