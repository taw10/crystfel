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


int hdf5_write(const char *filename, const uint16_t *data,
               int width, int height)
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	herr_t r;
	hsize_t size[2];
	hsize_t max_size[2];

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		fprintf(stderr, "Couldn't create file: %s\n", filename);
		return 1;
	}

	gh = H5Gcreate(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		fprintf(stderr, "Couldn't create group\n");
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
		fprintf(stderr, "Couldn't create dataset\n");
		H5Fclose(fh);
		return 1;
	}

	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, max_size);
	printf("Data dimensions %i %i (max %i %i)\n",
	                           (int)size[1], (int)size[0],
	                           (int)max_size[1], (int)max_size[0]);

	r = H5Dwrite(dh, H5T_NATIVE_UINT16, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		fprintf(stderr, "Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}

	H5Fclose(fh);

	return 0;
}


int hdf5_read(struct image *image, const char *filename)
{
	hid_t fh, sh, dh;	/* File, dataspace and data handles */
	herr_t r;
	hsize_t size[2];
	hsize_t max_size[2];
	uint16_t *buf;

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		/* TODO: Try other formats here. */
		fprintf(stderr, "Couldn't open file: %s\n", filename);
		return 1;
	}

	dh = H5Dopen(fh, "/data/data", H5P_DEFAULT);
	if ( dh < 0 ) {
		fprintf(stderr, "Couldn't open dataset\n");
		H5Fclose(fh);
		return 1;
	}

	sh = H5Dget_space(dh);
	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		fprintf(stderr, "Dataset is not two-dimensional\n");
		H5Fclose(fh);
		return 1;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	buf = malloc(sizeof(float)*size[0]*size[1]);

	r = H5Dread(dh, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if ( r < 0 ) {
		fprintf(stderr, "Couldn't read data\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return 1;
	}

	image->data = buf;
	image->height = size[0];
	image->width = size[1];
	image->x_centre = image->width/2;
	image->y_centre = image->height/2;

	H5Fclose(fh);

	return 0;
}
