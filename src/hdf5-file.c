/*
 * hdf5.c
 *
 * Read/write HDF5 data files
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
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


struct hdfile {

	const char      *path;  /* Current data path */

	struct image    *image;

	size_t          nx;  /* Image width */
	size_t          ny;  /* Image height */

	hid_t           fh;  /* HDF file handle */
	hid_t           dh;  /* Dataset handle */
	hid_t           sh;  /* Dataspace handle */
};


struct hdfile *hdfile_open(const char *filename)
{
	struct hdfile *f;

	f = malloc(sizeof(struct hdfile));
	if ( f == NULL ) return NULL;

	f->fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( f->fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		free(f);
		return NULL;
	}

	return f;
}


int hdfile_set_image(struct hdfile *f, const char *path)
{
	hsize_t size[2];
	hsize_t max_size[2];

	f->dh = H5Dopen(f->fh, path, H5P_DEFAULT);
	if ( f->dh < 0 ) {
		ERROR("Couldn't open dataset\n");
		return -1;
	}

	f->sh = H5Dget_space(f->dh);
	if ( H5Sget_simple_extent_ndims(f->sh) != 2 ) {
		ERROR("Dataset is not two-dimensional\n");
		return -1;
	}

	H5Sget_simple_extent_dims(f->sh, size, max_size);

	f->nx = size[0];
	f->ny = size[1];

	return 0;
}


int hdfile_get_width(struct hdfile *f)
{
	return f->nx;
}


int hdfile_get_height(struct hdfile *f)
{
	return f->ny;
}


void hdfile_close(struct hdfile *f)
{
	H5Fclose(f->fh);
	free(f->image);
	free(f);
}


static void *hdfile_bin(uint16_t *in, int inw, int inh,
                        int binning, uint16_t *maxp)
{
	uint16_t *data;
	int x, y;
	int w, h;
	uint16_t max;

	w = inw / binning;
	h = inh / binning;      /* Some pixels might get discarded */

	data = malloc(w*h*sizeof(uint16_t));
	max = 0;

	for ( x=0; x<w; x++ ) {
	for ( y=0; y<h; y++ ) {

		/* Big enough to hold large values */
		unsigned long long int total;
		size_t xb, yb;

		total = 0;
		for ( xb=0; xb<binning; xb++ ) {
		for ( yb=0; yb<binning; yb++ ) {

			total += in[inh*(binning*x+xb)+binning*y+yb];

		}
		}

		data[y+h*(w-1-x)] = total / (binning * binning);
		if ( data[y+h*(w-1-x)] > max ) max = data[y+h*(w-1-x)];

	}
	}

	*maxp = max;
	return data;
}


uint16_t *hdfile_get_image_binned(struct hdfile *f, int binning, uint16_t *max)
{
	struct image *image;
	uint16_t *data;

	image = malloc(sizeof(struct image));
	if ( image == NULL ) return NULL;

	hdf5_read(f, image);
	f->image = image;

	data = hdfile_bin(image->data, f->nx, f->ny, binning, max);

	return data;
}


int hdf5_write(const char *filename, const uint16_t *data,
               int width, int height)
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

	gh = H5Gcreate(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
		return 1;
	}

	size[0] = width;
	size[1] = height;
	max_size[0] = width;
	max_size[1] = height;
	sh = H5Screate_simple(2, size, max_size);

	dh = H5Dcreate(gh, "data", H5T_NATIVE_UINT16, sh,
	               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return 1;
	}

	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, max_size);

	r = H5Dwrite(dh, H5T_NATIVE_UINT16, H5S_ALL,
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


int hdf5_read(struct hdfile *f, struct image *image)
{
	herr_t r;
	uint16_t *buf;

	buf = malloc(sizeof(float)*f->nx*f->ny);

	r = H5Dread(f->dh, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read data\n");
		H5Dclose(f->dh);
		return 1;
	}

	image->data = buf;
	image->height = f->nx;
	image->width = f->ny;
	image->x_centre = image->width/2;
	image->y_centre = image->height/2;

	return 0;
}
