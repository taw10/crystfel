/*
 * hdf5.c
 *
 * Read/write HDF5 data files
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <hdf5.h>

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


int hdfile_get_width(struct hdfile *f)
{
	return f->ny;
}


int hdfile_get_height(struct hdfile *f)
{
	return f->nx;
}


int get_peaks(struct image *image, struct hdfile *f)
{
	hid_t dh, sh;
	hsize_t size[2];
	hsize_t max_size[2];
	int i;
	float *buf;
	herr_t r;

	dh = H5Dopen2(f->fh, "/processing/hitfinder/peakinfo", H5P_DEFAULT);

	if ( dh < 0 ) dh = H5Dopen2(f->fh, "/data/peakinfo", H5P_DEFAULT);

	if ( dh < 0 ) {
		ERROR("No peak list found!\n");
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

	if ( size[1] != 3 ) {
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

		unsigned int x, y;
		float val;

		x = buf[3*i+0];
		y = buf[3*i+1];
		val = buf[3*i+2];

		image_add_feature(image->features, x, y, image, val, NULL);

	}
	STATUS("Got %i peaks from peak list.\n", (int)size[0]);

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


int hdf5_write(const char *filename, const void *data,
               int width, int height, int type)
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	herr_t r;
	hsize_t size[2];
	hsize_t max_size[2];

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

	size[0] = height;
	size[1] = width;
	max_size[0] = height;
	max_size[1] = width;
	sh = H5Screate_simple(2, size, max_size);

	dh = H5Dcreate2(gh, "data", type, sh,
	                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return 1;
	}

	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, max_size);

	r = H5Dwrite(dh, type, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}

	H5Gclose(gh);
	H5Dclose(dh);
	H5Fclose(fh);

	return 0;
}


static double get_wavelength(struct hdfile *f)
{
	herr_t r;
	hid_t dh;
	double lambda;
	int nm = 1;

	dh = H5Dopen2(f->fh, "/LCLS/photon_wavelength_nm", H5P_DEFAULT);
	if ( dh < 0 ) {
		dh = H5Dopen2(f->fh, "/LCLS/photon_wavelength_A", H5P_DEFAULT);
		if ( dh < 0 ) return -1.0;
		nm = 0;

	}

	r = H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, &lambda);
	H5Dclose(dh);

	if ( r < 0 ) return -1.0;

	/* Convert nm -> m */
	if ( nm ) return lambda / 1.0e9;
	return lambda / 1.0e10;
}


static double get_f0(struct hdfile *f)
{
	herr_t r;
	hid_t dh;
	double f0;

	dh = H5Dopen2(f->fh, "/LCLS/f_11_ENRC", H5P_DEFAULT);
	if ( dh < 0 ) return -1.0;

	r = H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, &f0);
	H5Dclose(dh);
	if ( r < 0 ) return -1.0;

	return f0;
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
		STATUS("No saturation table found.\n");
		return;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		H5Dclose(dh);
		ERROR("Couldn't get dataspace for saturation table.\n");
		return;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Saturation table has the wrong dimensionality (%i).\n",
		      H5Sget_simple_extent_ndims(sh));
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
	STATUS("Corrected %i saturation values\n", (int)size[0]);

	free(buf);
	H5Sclose(sh);
	H5Dclose(dh);
}


int hdf5_read(struct hdfile *f, struct image *image, int satcorr,
              double nominal_photon_energy)
{
	herr_t r;
	float *buf;
	uint16_t *flags;
	hid_t mask_dh;

	image->width = hdfile_get_width(f);
	image->height = hdfile_get_height(f);
	STATUS("%i, %i\n", image->width, image->height);

	buf = malloc(sizeof(float)*f->nx*f->ny);

	r = H5Dread(f->dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read data\n");
		free(buf);
		return 1;
	}
	image->data = buf;

	mask_dh = H5Dopen2(f->fh, "/processing/hitfinder/masks", H5P_DEFAULT);
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

	/* Read wavelength from file */
	image->lambda = get_wavelength(f);
	if ( image->lambda < 0.0 ) {
		ERROR("Couldn't read wavelength - using nominal value.\n");
		image->lambda = ph_en_to_lambda(eV_to_J(nominal_photon_energy));
	}

	image->f0 = get_f0(f);
	if ( image->f0 < 0.0 ) {
		ERROR("Couldn't read incident intensity - using 1.0.\n");
		image->f0 = 1.0;
		image->f0_available = 0;
	} else {
		image->f0_available = 1;
	}

	if ( satcorr ) debodge_saturation(f, image);

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


char *hdfile_get_string_value(struct hdfile *f, const char *name)
{
	hid_t dh;
	hid_t sh;
	hsize_t size;
	hsize_t max_size;
	hid_t type;
	hid_t class;

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

	sh = H5Dget_space(dh);
	if ( H5Sget_simple_extent_ndims(sh) != 1 ) goto fail;

	H5Sget_simple_extent_dims(sh, &size, &max_size);
	if ( size != 1 ) {
		H5Dclose(dh);
		goto fail;
	}

	switch ( class ) {
	case H5T_FLOAT : {

		herr_t r;
		double buf;
		char *tmp;

		r = H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, &buf);
		if ( r < 0 ) goto fail;

		tmp = malloc(256);
		snprintf(tmp, 255, "%f", buf);

		return tmp;

	}
	case H5T_INTEGER : {

		herr_t r;
		int buf;
		char *tmp;

		r = H5Dread(dh, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, &buf);
		if ( r < 0 ) goto fail;

		tmp = malloc(256);
		snprintf(tmp, 255, "%d", buf);

		return tmp;
	}
	default : {
		goto fail;
	}
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

	gh = H5Gopen2(f->fh, parent, H5P_DEFAULT);
	if ( gh < 0 ) {
		*n = 0;
		return NULL;
	}

	if ( H5Gget_num_objs(gh, &num) < 0 ) {
		/* Whoopsie */
		*n = 0;
		return NULL;
	}
	*n = num;
	if ( num == 0 ) return NULL;  /* Bail out now */

	res = malloc(num*sizeof(char *));
	is_image = malloc(num*sizeof(int));
	is_group = malloc(num*sizeof(int));
	*p_is_image = is_image;
	*p_is_group = is_group;

	for ( i=0; i<num; i++ ) {

		char buf[256];
		int type;

		H5Gget_objname_by_idx(gh, i, buf, 255);
		res[i] = malloc(256);
		if ( strlen(parent) > 1 ) {
			snprintf(res[i], 255, "%s/%s", parent, buf);
		} else {
			snprintf(res[i], 255, "%s%s", parent, buf);
		} /* ick */

		type = H5Gget_objtype_by_idx(gh, i);
		is_image[i] = 0;
		is_group[i] = 0;
		if ( type == H5G_GROUP ) {
			is_group[i] = 1;
		} else if ( type == H5G_DATASET ) {
			hid_t dh;
			dh = H5Dopen2(gh, res[i], H5P_DEFAULT);
			is_image[i] = looks_like_image(dh);
			H5Dclose(dh);
		}

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
