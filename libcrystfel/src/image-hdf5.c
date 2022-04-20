/*
 * image-hdf5.c
 *
 * Image loading, HDF5 parts
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
 *   2014-2018 Valerio Mariani
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include "image.h"
#include "utils.h"
#include "detgeom.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"


/* Get the path parts of the event ID
 * e.g. ev_orig = abc/def/ghi//5/2/7
 *     -> [abc, def, ghi], with *pn_plvals=3.
 *
 * Not part of public API.  Not "static" for testing. */
char **read_path_parts(const char *ev_orig, int *pn_plvals)
{
	char **plvals;
	char *ev;
	int n_plvals = 0;
	char *start;

	plvals = malloc(MAX_PATH_PARTS*sizeof(char *));
	if ( plvals == NULL ) return NULL;

	if ( ev_orig == NULL ) {
		/* No ev -> no path parts */
		*pn_plvals = 0;
		return plvals;
	}

	ev = strdup(ev_orig);
	if ( ev == NULL ) {
		free(plvals);
		return NULL;
	}

	start = ev;
	do {

		char *sep;

		sep = strchr(start, '/');

		if ( sep == NULL ) {
			/* This would be very strange, because it
			 * must at least have // */
			ERROR("Couldn't read path parts ('%s')\n",
			      start);
			free(ev);
			free(plvals);
			return NULL;
		}

		/* Remaining string starts with '/' is end condition */
		if ( sep == start ) break;

		if ( n_plvals == MAX_PATH_PARTS ) {
			ERROR("Too many path parts: %s\n", ev_orig);
			free(ev);
			free(plvals);
			return NULL;
		}

		sep[0] = '\0';
		plvals[n_plvals++] = strdup(start);

		start = sep+1;

	} while ( 1 );

	free(ev);
	*pn_plvals = n_plvals;
	return plvals;
}


/* Get the dimension parts of the event ID
 * e.g. ev_orig = abc/def/ghi//5/2/7
 *     -> [5, 2, 7], with *pn_dvals=3
 *
 * Not part of public API.  Not "static" for testing. */
int *read_dim_parts(const char *ev_orig, int *pn_dvals)

{
	char *ev;
	int n_dvals = 0;
	int *dvals;
	char *start;
	int done;

	if ( ev_orig == NULL ) ev_orig = "//";

	/* Valid event ID? (Just the part after //, please) */
	ev = strstr(ev_orig, "//");
	if ( ev == NULL ) return NULL;

	dvals = malloc(MAX_DIMS*sizeof(int));
	if ( dvals == NULL ) return NULL;

	if ( ev[2] == '\0' ) {
		/* No dimension parts - early bailout */
		*pn_dvals = 0;
		return dvals;  /* NB Not NULL */
	}

	ev = strdup(ev+2);
	if ( ev == NULL ) {
		free(dvals);
		return NULL;
	}

	start = ev;
	done = 0;
	do {

		char *sep = strchr(start, '/');

		if ( sep != NULL ) {
			sep[0] = '\0';
		} else {
			done = 1;
		}

		if ( n_dvals == MAX_PATH_PARTS ) {
			ERROR("Too many path parts: %s\n", ev_orig);
			free(ev);
			free(dvals);
			return NULL;
		}

		if ( start[0] == '\0' ) {
			ERROR("Missing dimension: %s\n", ev_orig);
			free(ev);
			free(dvals);
			return NULL;
		}

		dvals[n_dvals++] = atoi(start);

		start = sep+1;

	} while ( !done );

	free(ev);
	*pn_dvals = n_dvals;
	return dvals;
}


static int imh_num_path_placeholders(const char *pattern)
{
	size_t l, i;
	int n_pl_exp = 0;

	l = strlen(pattern);
	for ( i=0; i<l; i++ ) {
		if ( pattern[i] == '%' ) n_pl_exp++;
	}
	return n_pl_exp;
}


/* ev = abc/def/ghi//5/2/7
 * pattern = /data/%/somewhere/%/%/data
 * output = /data/abc/somewhere/def/ghi/data
 *
 * Not part of public API.  Not "static" for testing.
 */
char *substitute_path(const char *ev, const char *pattern, int skip_ok)
{
	char **plvals;
	int n_plvals;
	int n_pl_exp;
	size_t total_len;
	int i;
	char *subs;
	const char *start;
	const char *pl_pos;

	if ( pattern == NULL ) {
		ERROR("Pattern cannot be NULL\n");
		return NULL;
	}

	/* tag_1334/wibble1//3/6 */
	plvals = read_path_parts(ev, &n_plvals);
	if ( plvals == NULL ) return NULL;

	/* /data/%/image_data/%/rawpixels */
	n_pl_exp = imh_num_path_placeholders(pattern);

	if ( n_plvals < n_pl_exp ) {
		ERROR("Event ID does not have enough path placeholder values: "
		      "event ID '%s' (%i) into pattern '%s' (%i)\n",
		      ev, n_plvals, pattern, n_pl_exp);
		return NULL;
	}

	if ( (n_plvals > n_pl_exp) && !skip_ok ) {
		ERROR("Event ID has too many path placeholder values: "
		      "event ID '%s' (%i) into pattern '%s' (%i)\n",
		      ev, n_plvals, pattern, n_pl_exp);
		return NULL;
	}

	if ( n_pl_exp == 0 ) {
		/* No placeholders in path */
		for ( i=0; i<n_plvals; i++ ) {
			free(plvals[i]);
		}
		free(plvals);
		return strdup(pattern);
	}

	total_len = strlen(pattern) - n_pl_exp;
	for ( i=0; i<n_plvals; i++ ) {
		total_len += strlen(plvals[i]);
	}
	subs = malloc(total_len+1);
	if ( subs == NULL ) {
		free(plvals);
		return NULL;
	}

	pl_pos = strchr(pattern, '%');
	if ( pl_pos == NULL ) {
		ERROR("Expected a placeholder char (%): '%s'\n",
		      pattern);
		return NULL;
	}
	strncpy(subs, pattern, pl_pos-pattern);
	subs[pl_pos-pattern] = '\0';

	start = pl_pos+1;
	for ( i=0; i<n_pl_exp; i++ ) {

		/* Add the placeholder's value */
		strcat(subs, plvals[i]);
		free(plvals[i]);

		/* Add the chars up to the next placeholder... */
		pl_pos = strchr(start, '%');
		if ( pl_pos == NULL ) {
			/* ... or the end */
			pl_pos = start+strlen(start);
		}
		strncat(subs, start, pl_pos-start);
		start = pl_pos+1;
	}

	free(plvals);

	return subs;
}


#ifdef HAVE_HDF5

static void make_placeholder_skip(signed int *dt_dims,
                                  signed int *panel_dims)
{
	int i;
	int n_dt = 0;
	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( panel_dims[i] != DIM_PLACEHOLDER ) {
			dt_dims[n_dt++] = panel_dims[i];
		}
	}
}


static int imh_num_placeholders(const struct panel_template *p)
{
	int i;
	int n_pl = 0;
	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( p->dims[i] == DIM_PLACEHOLDER ) n_pl++;
	}
	return n_pl;
}


static int total_dimensions(const struct panel_template *p)
{
	int i;
	int n_dim = 0;
	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( p->dims[i] != DIM_UNDEFINED ) n_dim++;
	}
	return n_dim;
}


static void close_hdf5(hid_t fh)
{
        int n_ids, i;
        hid_t ids[2048];

        n_ids = H5Fget_obj_ids(fh, H5F_OBJ_ALL, 2048, ids);

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

        H5Fclose(fh);
}




static int load_hdf5_hyperslab(struct panel_template *p,
                               hid_t fh,
                               const char *event,
                               void **pdata,
                               hid_t el_type, size_t el_size,
                               int skip_placeholders_ok,
                               const char *path_spec)
{
	int total_dt_dims;
	int plh_dt_dims;
	int dt_dims[MAX_DIMS];
	int n_dt_dims;
	herr_t r;
	hsize_t *f_offset, *f_count;
	hid_t dh;
	herr_t check;
	hid_t dataspace, memspace;
	hsize_t dims[2];
	char *panel_full_path;
	void *data;
	int ndims;
	int dim;
	int *dim_vals;
	int n_dim_vals;
	int pl_pos;

	panel_full_path = substitute_path(event, path_spec,
	                                  skip_placeholders_ok);
	if ( panel_full_path == NULL ) {
		ERROR("Invalid path substitution: '%s' '%s'\n",
		      event, path_spec);
		return 1;
	}

	dh = H5Dopen2(fh, panel_full_path, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Cannot open data for panel %s (%s)\n",
		      p->name, panel_full_path);
		free(panel_full_path);
		return 1;
	}

	free(panel_full_path);

	/* Set up dataspace for file
	 * (determine where to read the data from) */
	dataspace = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(dataspace);
	if ( ndims < 0 ) {
		ERROR("Failed to get number of dimensions for panel %s\n",
		      p->name);
		return 1;
	}

	/* Does the array have the expected number of dimensions? */
	total_dt_dims = total_dimensions(p);
	plh_dt_dims = imh_num_placeholders(p);
	if ( ndims != total_dt_dims ) {
		/* If the dimensions match after excluding
		 * placeholders, it's OK - probably a static mask
		 * in a multi-event file. */
		if ( skip_placeholders_ok
		  && (ndims == total_dt_dims - plh_dt_dims) )
		{
			make_placeholder_skip(dt_dims, p->dims);
			n_dt_dims = total_dt_dims - plh_dt_dims;
		} else {
			ERROR("Unexpected number of dimensions for "
			      "panel %s (%i, but expected %i or %i)\n",
			      p->name, ndims, total_dt_dims,
			      total_dt_dims - plh_dt_dims);
			return 1;
		}
	} else {
		int i;
		for ( i=0; i<MAX_DIMS; i++ ) {
			dt_dims[i] = p->dims[i];
		}
		n_dt_dims = total_dt_dims;
	}

	f_offset = malloc(ndims*sizeof(hsize_t));
	f_count = malloc(ndims*sizeof(hsize_t));
	if ( (f_offset == NULL) || (f_count == NULL ) ) {
		ERROR("Failed to allocate offset or count.\n");
		return 1;
	}

	/* Get those placeholder values from the event ID */
	dim_vals = read_dim_parts(event, &n_dim_vals);

	pl_pos = 0;
	for ( dim=0; dim<n_dt_dims; dim++ ) {

		switch ( dt_dims[dim] ) {

			case DIM_FS:
			f_offset[dim] = p->orig_min_fs;
			f_count[dim] = p->orig_max_fs - p->orig_min_fs+1;
			break;

			case DIM_SS:
			f_offset[dim] = p->orig_min_ss;
			f_count[dim] = p->orig_max_ss - p->orig_min_ss+1;
			break;

			case DIM_PLACEHOLDER:
			f_offset[dim] = dim_vals[pl_pos++];
			f_count[dim] = 1;
			break;

			case DIM_UNDEFINED:
			ERROR("Undefined dimension found!\n");
			break;

			default:
			/* Fixed value */
			f_offset[dim] = dt_dims[dim];
			f_count[dim] = 1;
			break;

		}
	}

	free(dim_vals);

	check = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
	                            f_offset, NULL, f_count, NULL);
	if ( check < 0 ) {
		ERROR("Error selecting file dataspace for panel %s\n",
		      p->name);
		free(f_offset);
		free(f_count);
		return 1;
	}

	dims[0] = p->orig_max_ss - p->orig_min_ss + 1;
	dims[1] = p->orig_max_fs - p->orig_min_fs + 1;
	memspace = H5Screate_simple(2, dims, NULL);

	data = malloc(dims[0]*dims[1]*el_size);
	if ( data == NULL ) {
		ERROR("Failed to allocate panel %s\n", p->name);
		free(f_offset);
		free(f_count);
		return 1;
	}

	r = H5Dread(dh, el_type, memspace, dataspace, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't read data for panel %s\n",
		      p->name);
		free(f_offset);
		free(f_count);
		free(data);
		return 1;
	}

	free(f_offset);
	free(f_count);

	*pdata = data;
	return 0;
}


static hid_t open_hdf5_file(const char *filename)
{
	hid_t fh;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return -1;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open HDF5 file: %s\n", filename);
		return -1;
	}

	return fh;
}


static hid_t open_hdf5(struct image *image)
{
	if ( image->data_block == NULL ) {

		return open_hdf5_file(image->filename);

	} else {

		hid_t fh;

		fh = H5LTopen_file_image(image->data_block,
		                         image->data_block_size,
		                         H5LT_FILE_IMAGE_DONT_COPY
		                         | H5LT_FILE_IMAGE_DONT_RELEASE);

		if ( fh < 0 ) {
			ERROR("Couldn't open HDF5 image (%p %lli)\n",
			      image->data_block, image->data_block_size);
			return -1;
		}

		return fh;
	}
}


int image_hdf5_read(struct image *image,
                    const DataTemplate *dtempl)
{
	int i;
	hid_t fh;

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate data array.\n");
		return 1;
	}

	if ( image->ev == NULL ) {
		image->ev = "//";
	}

	/* Set all pointers to NULL for easier clean-up */
	for ( i=0; i<dtempl->n_panels; i++ ) image->dp[i] = NULL;

	fh = open_hdf5(image);
	if ( fh < 0 ) {
		ERROR("Failed to open file\n");
		return 1;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( load_hdf5_hyperslab(&dtempl->panels[i], fh,
		                         image->ev, (void *)&image->dp[i],
		                         H5T_NATIVE_FLOAT,
		                         sizeof(float), 0,
		                         dtempl->panels[i].data) )
		{
			ERROR("Failed to load panel data\n");
			close_hdf5(fh);
			return 1;
		}
	}

	close_hdf5(fh);
	return 0;
}


float *image_hdf5_read_satmap(struct panel_template *p,
                              const char *filename,
                              const char *event,
                              const char *map_location)
{
	float *map = NULL;
	hid_t fh;

	fh = open_hdf5_file(filename);
	if ( fh < 0 ) return NULL;

	if ( load_hdf5_hyperslab(p, fh, event,
	                         (void *)&map, H5T_NATIVE_FLOAT,
	                         sizeof(float), 1, map_location) )
	{
		ERROR("Failed to load saturation map data\n");
		free(map);
		return NULL;
	}

	close_hdf5(fh);

	return map;
}



int image_hdf5_read_mask(struct panel_template *p,
                         const char *filename, const char *event,
                         int *bad, const char *mask_location,
                         int mask_good, int mask_bad)
{
	int p_w, p_h;
	int *mask = NULL;
	long unsigned int j;
	hid_t fh;

	p_w = p->orig_max_fs - p->orig_min_fs + 1;
	p_h = p->orig_max_ss - p->orig_min_ss + 1;

	fh = open_hdf5_file(filename);
	if ( fh < 0 ) return 1;

	if ( load_hdf5_hyperslab(p, fh, event,
	                         (void *)&mask, H5T_NATIVE_INT,
	                         sizeof(int), 1, mask_location) )
	{
		ERROR("Failed to load mask data\n");
		free(mask);
		return 1;
	}

	close_hdf5(fh);

	for ( j=0; j<p_w*p_h; j++ ) {

		/* Bad if it's missing any of the "good" bits */
		if ( (mask[j] & mask_good) != mask_good ) bad[j] = 1;

		/* Bad if it has any of the "bad" bits. */
		if ( mask[j] & mask_bad ) bad[j] = 1;

	}

	free(mask);
	return 0;
}


static char *read_single_fixed_string(hid_t dh)
{
	hid_t sh, type;
	herr_t r;
	size_t size;
	char *tmp;

	sh = H5Dget_space(dh);
	if ( H5Sget_simple_extent_ndims(sh) ) {
		ERROR("Non-scalar string\n");
		return NULL;
	}

	sh = H5Screate(H5S_SCALAR);
	type = H5Dget_type(dh);
	size = H5Tget_size(type);
	tmp = malloc(size+1);
	if ( tmp == NULL ) {
		H5Tclose(type);
		return NULL;
	}
	r = H5Dread(dh, type, sh, H5S_ALL, H5P_DEFAULT, tmp);
	H5Sclose(sh);
	H5Tclose(type);
	if ( r < 0 ) {
		free(tmp);
		ERROR("Couldn't read scalar string\n");
		return NULL;
	} else {
		tmp[size] = '\0';
		chomp(tmp);
		return tmp;
	}
}


static char *read_single_vlen_string(hid_t dh)
{
	herr_t r;
	char *tmp;
	hid_t type;

	type = H5Dget_type(dh);
	r = H5Dread(dh, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
	H5Tclose(type);
	if ( r < 0 ) {
		ERROR("Couldn't read vlen string\n");
		return NULL;
	}

	chomp(tmp);
	return tmp;
}


int image_hdf5_read_header_to_cache(struct image *image, const char *name)
{
	hid_t dh;
	hid_t type;
	hid_t class;
	hid_t sh;
	hid_t ms;
	hsize_t *f_offset = NULL;
	hsize_t *f_count = NULL;
	hsize_t m_offset[1];
	hsize_t m_count[1];
	hsize_t msdims[1];
	hsize_t size[64];
	herr_t r;
	herr_t check;
	int ndims;
	int i;
	char *subst_name = NULL;
	hid_t fh;
	int *dim_vals;
	int n_dim_vals;
	int dim_val_pos;

	fh = open_hdf5(image);
	if ( fh < 0 ) {
		ERROR("Couldn't open file (header): %s\n", image->filename);
		return 1;
	}

	subst_name = substitute_path(image->ev, name, 1);
	if ( subst_name == NULL ) {
		ERROR("Invalid event ID '%s'\n", image->ev);
		close_hdf5(fh);
		return 1;
	}

	dh = H5Dopen2(fh, subst_name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("No such numeric field '%s'\n", subst_name);
		close_hdf5(fh);
		return 1;
	}

	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	switch ( class ) {

		/* Acceptable types */
		case H5T_FLOAT:
		case H5T_INTEGER:
		case H5T_STRING:
		break;

		default:
		ERROR("HDF5 header is not a recognised type (%s).\n",
		      subst_name);
		close_hdf5(fh);
		return 1;
	}

	/* Get the dimensionality.  We have to cope with scalars expressed as
	 * arrays with all dimensions 1, as well as zero-d arrays. */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > 64 ) {
		ERROR("Too many dimensions for numeric value\n");
		close_hdf5(fh);
		return 1;
	}
	H5Sget_simple_extent_dims(sh, size, NULL);

	/* We want to read the value as a scalar */
	m_offset[0] = 0;
	m_count[0] = 1;
	msdims[0] = 1;
	ms = H5Screate_simple(1, msdims, NULL);

	if ( ndims == 0 ) {

		/* Easy case, because value is a scalar */
		if ( class == H5T_FLOAT ) {

			double val;
			r = H5Dread(dh, H5T_NATIVE_DOUBLE, ms, sh, H5P_DEFAULT,
			            &val);
			if ( r < 0 )  {
				ERROR("Couldn't read scalar value from %s.\n",
				      subst_name);
				free(subst_name);
				close_hdf5(fh);
				return 1;
			}
			image_cache_header_float(image, name, val);
			return 0;

		} else if ( class == H5T_INTEGER ) {

			int val;
			r = H5Dread(dh, H5T_NATIVE_INT, ms, sh, H5P_DEFAULT,
			            &val);
			if ( r < 0 )  {
				ERROR("Couldn't read scalar value from %s.\n",
				      subst_name);
				free(subst_name);
				close_hdf5(fh);
				return 1;
			}
			image_cache_header_int(image, name, val);
			return 0;

		} else if ( class == H5T_STRING ) {

			htri_t v;
			hid_t stype;
			char *val;
			int rv;

			stype = H5Dget_type(dh);
			v = H5Tis_variable_str(stype);
			H5Tclose(stype);

			if ( v == 0 ) {
				val = read_single_fixed_string(dh);
			} else if ( v > 0 ) {
				val = read_single_vlen_string(dh);
			} else {
				ERROR("Unrecognised string type: %s\n",
				      subst_name);
				val = NULL;
			}

			if ( val != NULL ) {
				image_cache_header_str(image, name, val);
				free(val);
				rv = 0;
			} else {
				ERROR("Failed to read string '%s'\n",
				      subst_name);
				rv = 1;
			}

			free(subst_name);
			close_hdf5(fh);
			return rv;

		} else {
			/* Should never be reached */
			ERROR("Invalid HDF5 class %i\n", class);
			return 1;
		}
	}

	dim_vals = read_dim_parts(image->ev, &n_dim_vals);
	if ( dim_vals == NULL ) {
		ERROR("Couldn't parse event '%s'\n");
		close_hdf5(fh);
		return 1;
	}

	f_offset = malloc(ndims*sizeof(hsize_t));
	f_count = malloc(ndims*sizeof(hsize_t));
	if ( (f_offset == NULL) || (f_count == NULL) ) {
		ERROR("Couldn't allocate dimension arrays\n");
		close_hdf5(fh);
		return 1;
	}

	/* Every dimension of the dataset must either be size 1 or
	 * large enough to contain the next value from the event ID */
	dim_val_pos = 0;
	for ( i=0; i<ndims; i++ ) {

		if ( size[i] != 1 ) {

			if ( size[i] <= dim_vals[dim_val_pos] ) {
				ERROR("Array of scalar values is too "
				      "small (%s, dim %i, ev value %i,"
				      " size %i)\n",
				      subst_name, i,
				      dim_vals[dim_val_pos], size[i]);
				close_hdf5(fh);
				return 1;
			}

			f_offset[i] = dim_vals[dim_val_pos];
			f_count[i] = 1;
			dim_val_pos++;

		} else {

			f_offset[i] = 0;
			f_count[i] = 1;

		}

	}

	check = H5Sselect_hyperslab(sh, H5S_SELECT_SET,
	                            f_offset, NULL, f_count, NULL);
	if ( check < 0 ) {
		ERROR("Error selecting dataspace for header value\n");
		free(f_offset);
		free(f_count);
		free(subst_name);
		close_hdf5(fh);
		return 1;
	}

	free(f_offset);
	free(f_count);

	ms = H5Screate_simple(1,msdims,NULL);
	check = H5Sselect_hyperslab(ms, H5S_SELECT_SET,
	                            m_offset, NULL, m_count, NULL);
	if ( check < 0 ) {
		ERROR("Error selecting memory dataspace for header value\n");
		close_hdf5(fh);
		free(subst_name);
		return 1;
	}

	if ( class == H5T_FLOAT ) {

		double val;
		r = H5Dread(dh, H5T_NATIVE_DOUBLE, ms, sh, H5P_DEFAULT, &val);
		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			close_hdf5(fh);
			free(subst_name);
			return 1;
		}

		image_cache_header_float(image, name, val);
		close_hdf5(fh);
		free(subst_name);
		return 0;

	} else if ( class == H5T_INTEGER ) {

		int val;
		r = H5Dread(dh, H5T_NATIVE_INT, ms, sh, H5P_DEFAULT, &val);
		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			close_hdf5(fh);
			free(subst_name);
			return 1;
		}

		image_cache_header_int(image, name, val);
		close_hdf5(fh);
		free(subst_name);
		return 0;

	} else if ( class == H5T_STRING ) {

		hid_t stype;

		stype = H5Dget_type(dh);
		if ( H5Tis_variable_str(stype) ) {

			/* Vlen string from array */

			herr_t rv;
			char *val;

			rv = H5Dread(dh, stype, ms, sh, H5P_DEFAULT, &val);
			if ( rv < 0 ) {
				ERROR("Can't read HDF5 vlen string from array - %s\n",
				      subst_name);
				free(subst_name);
				close_hdf5(fh);
				return 1;
			} else {

				chomp(val);
				image_cache_header_str(image, name, val);
				free(val);
				close_hdf5(fh);
				free(subst_name);
				return 0;
			}

		} else {

			/* Fixed-length string from array */

			herr_t rv;
			char *val;
			size_t ssize;

			ssize = H5Tget_size(stype);
			val = malloc(ssize+1);
			if ( val == NULL ) {
				close_hdf5(fh);
				free(subst_name);
				return 1;
			}
			rv = H5Dread(dh, stype, ms, sh, H5P_DEFAULT, val);
			H5Tclose(stype);
			if ( rv < 0 ) {
				ERROR("Couldn't read HDF5 fixed string from array - %s\n",
				      subst_name);
				close_hdf5(fh);
				free(subst_name);
				return 1;
			} else {

				val[ssize] = '\0';
				chomp(val);
				image_cache_header_str(image, name, val);
				free(val);
				close_hdf5(fh);
				free(subst_name);
				return 0;

			}

		}  /* I feel sick. */

	} else {
		/* Should never be reached */
		ERROR("Invalid HDF5 class %i\n", class);
		close_hdf5(fh);
		free(subst_name);
		return 1;
	}
}


static int read_peak_count(hid_t fh, char *path, int line,
                           int *num_peaks)
{

	hid_t dh, sh, mh;
	hsize_t size[1];
	hsize_t max_size[1];
	hsize_t offset[1], count[1];
	hsize_t m_offset[1], m_count[1], dimmh[1];
	int tw, r;

	dh = H5Dopen2(fh, path, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Data block %s not found.\n", path);
		return 1;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		H5Dclose(dh);
		ERROR("Couldn't get dataspace for data.\n");
		return 1;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 1 ) {
		ERROR("Data block %s has the wrong dimensionality (%i).\n",
		      path, H5Sget_simple_extent_ndims(sh));
		H5Sclose(sh);
		H5Dclose(dh);
		return 1;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	tw = size[0];

	if ( line > tw-1 ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Data block %s does not contain data for required event.\n",
		      path);
		return 1;
	}

	offset[0] = line;
	count[0] = 1;

	r = H5Sselect_hyperslab(sh, H5S_SELECT_SET,
	                        offset, NULL, count, NULL);
	if ( r < 0 ) {
		ERROR("Error selecting file dataspace "
		      "for data block %s\n", path);
		H5Dclose(dh);
		H5Sclose(sh);
		return 1;
	}

	m_offset[0] = 0;
	m_count[0] = 1;
	dimmh[0] = 1;
	mh = H5Screate_simple(1, dimmh, NULL);
	r = H5Sselect_hyperslab(mh, H5S_SELECT_SET,
	                        m_offset, NULL, m_count, NULL);
	if ( r < 0 ) {
		ERROR("Error selecting memory dataspace "
		      "for data block %s\n", path);
		H5Dclose(dh);
		H5Sclose(sh);
		H5Sclose(mh);
		return 1;
	}

	r = H5Dread(dh, H5T_NATIVE_INT, mh,
	            sh, H5P_DEFAULT, num_peaks);
	if ( r < 0 ) {
		ERROR("Couldn't read data for block %s, line %i\n", path, line);
		H5Dclose(dh);
		H5Sclose(sh);
		H5Sclose(mh);
		return 1;
	}

	H5Dclose(dh);
	H5Sclose(sh);
	H5Sclose(mh);
	return 0;
}


static float *read_peak_line(hid_t fh, char *path, int line,
                             int num_peaks)
{

	hid_t dh, sh, mh;
	hsize_t size[2];
	hsize_t max_size[2];
	hsize_t offset[2], count[2];
	hsize_t m_offset[2], m_count[2], dimmh[2];
	float *buf;
	int tw, r;

	dh = H5Dopen2(fh, path, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Data block (%s) not found.\n", path);
		return NULL;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		H5Dclose(dh);
		ERROR("Couldn't get dataspace for data.\n");
		return NULL;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Data block %s has the wrong dimensionality (%i).\n",
		      path, H5Sget_simple_extent_ndims(sh));
		H5Sclose(sh);
		H5Dclose(dh);
		return NULL;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);

	tw = size[0];
	if ( line > tw-1 ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Data block %s does not contain data for required event.\n",
		      path);
		return NULL;
	}

	/* NB The array might be bigger - Cheetah allocates in blocks of 2048 */
	if ( size[1] < num_peaks ) {
		ERROR("Data block %s is too small for the specified number of "
		      "peaks (has %i, expected %i)\n", path, size[1], num_peaks);
		H5Sclose(sh);
		H5Dclose(dh);
		return NULL;
	}

	offset[0] = line;
	offset[1] = 0;
	count[0] = 1;
	count[1] = size[1];

	r = H5Sselect_hyperslab(sh, H5S_SELECT_SET, offset, NULL, count, NULL);
	if ( r < 0 ) {
	    ERROR("Error selecting file dataspace "
	          "for data block %s\n", path);
	    H5Dclose(dh);
	    H5Sclose(sh);
	    return NULL;
	}

	m_offset[0] = 0;
	m_offset[1] = 0;
	m_count[0] = 1;
	m_count[1] = size[1];
	dimmh[0] = 1;
	dimmh[1] = size[1];

	mh = H5Screate_simple(2, dimmh, NULL);
	r = H5Sselect_hyperslab(mh, H5S_SELECT_SET,
	                        m_offset, NULL, m_count, NULL);
	if ( r < 0 ) {
		ERROR("Error selecting memory dataspace "
		      "for data block %s\n", path);
		H5Dclose(dh);
		H5Sclose(sh);
		H5Sclose(mh);
		return NULL;
	}

	buf = malloc(size[1]*sizeof(float));
	if ( buf == NULL ) return NULL;
	r = H5Dread(dh, H5T_NATIVE_FLOAT, mh, sh, H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read data for block %s, line %i\n", path, line);
		H5Dclose(dh);
		H5Sclose(sh);
		H5Sclose(mh);
		return NULL;
	}

	H5Dclose(dh);
	H5Sclose(sh);
	H5Sclose(mh);
	return buf;
}


ImageFeatureList *image_hdf5_read_peaks_cxi(const DataTemplate *dtempl,
                                            const char *filename,
                                            const char *event,
                                            int half_pixel_shift)
{
	ImageFeatureList *features;
	hid_t fh;
	char path_n[1024];
	char path_x[1024];
	char path_y[1024];
	char path_i[1024];
	int r;
	int pk;
	char *subst_name;
	int line;
	int num_peaks;
	float *buf_x;
	float *buf_y;
	float *buf_i;
	int *dim_vals;
	int n_dim_vals;

	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return NULL;
	}

	subst_name = substitute_path(event, dtempl->peak_list, 0);
	if ( subst_name == NULL ) {
		ERROR("Invalid peak path %s\n", subst_name);
		return NULL;
	}

	dim_vals = read_dim_parts(event, &n_dim_vals);
	if ( dim_vals == NULL ) {
		ERROR("Couldn't parse event '%s'\n");
		return NULL;
	}

	if ( n_dim_vals < 1 ) {
		ERROR("Not enough dimensions in event ID to use CXI "
		      "peak lists (%i)\n", n_dim_vals);
		return NULL;
	}

	line = dim_vals[0];
	free(dim_vals);

	snprintf(path_n, 1024, "%s/nPeaks", subst_name);
	snprintf(path_x, 1024, "%s/peakXPosRaw", subst_name);
	snprintf(path_y, 1024, "%s/peakYPosRaw", subst_name);
	snprintf(path_i, 1024, "%s/peakTotalIntensity", subst_name);

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file (peaks/cxi): %s\n", filename);
		return NULL;
	}

	r = read_peak_count(fh, path_n, line, &num_peaks);
	if ( r != 0 ) {
		close_hdf5(fh);
		return NULL;
	}

	buf_x = read_peak_line(fh, path_x, line, num_peaks);
	if ( buf_x == NULL ) {
		close_hdf5(fh);
		return NULL;
	}

	buf_y = read_peak_line(fh, path_y, line, num_peaks);
	if ( buf_y == NULL ) {
		close_hdf5(fh);
		return NULL;
	}

	buf_i = read_peak_line(fh, path_i, line, num_peaks);
	if ( buf_i == NULL ) {
		close_hdf5(fh);
		return NULL;
	}

	features = image_feature_list_new();

	for ( pk=0; pk<num_peaks; pk++ ) {

		float fs, ss, val;
		int pn;

		fs = buf_x[pk] + peak_offset;
		ss = buf_y[pk] + peak_offset;
		val = buf_i[pk];

		if ( data_template_slabby_file_to_panel_coords(dtempl,
		                                               &fs, &ss, &pn) )
		{
			ERROR("Failed to convert %i,%i to "
			      "panel-relative coordinates\n", fs, ss);
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
		}

	}

	close_hdf5(fh);

	return features;
}


ImageFeatureList *image_hdf5_read_peaks_hdf5(const DataTemplate *dtempl,
                                             const char *filename,
                                             const char *event,
                                             int half_pixel_shift)
{
	hid_t fh, dh, sh;
	hsize_t size[2];
	hsize_t max_size[2];
	int i;
	float *buf;
	herr_t r;
	int tw;
	char *subst_name;
	ImageFeatureList *features;
	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( dtempl->peak_list == NULL ) {
		ERROR("Peak location is not given in geometry file.\n");
		return NULL;
	}

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return NULL;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file (peaks/hdf5): %s\n", filename);
		return NULL;
	}

	subst_name = substitute_path(event, dtempl->peak_list, 0);
	if ( subst_name == NULL ) {
		ERROR("Invalid peak path: '%s' '%s'\n",
		      event, dtempl->peak_list);
		close_hdf5(fh);
		return NULL;
	}

	dh = H5Dopen2(fh, subst_name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Peak list (%s) not found.\n", subst_name);
		free(subst_name);
		close_hdf5(fh);
		return NULL;
	}
	free(subst_name);

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		ERROR("Couldn't get dataspace for peak list.\n");
		close_hdf5(fh);
		return NULL;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Peak list has the wrong dimensionality (%i).\n",
		      H5Sget_simple_extent_ndims(sh));
		close_hdf5(fh);
		return NULL;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);
	H5Sclose(sh);

	tw = size[1];
	if ( (tw != 3) && (tw != 4) ) {
		ERROR("Peak list has the wrong dimensions.\n");
		close_hdf5(fh);
		return NULL;
	}

	buf = malloc(sizeof(float)*size[0]*size[1]);
	if ( buf == NULL ) {
		ERROR("Couldn't reserve memory for the peak list.\n");
		close_hdf5(fh);
		return NULL;
	}
	r = H5Dread(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read peak list.\n");
		close_hdf5(fh);
		return NULL;
	}

	features = image_feature_list_new();
	if ( features == NULL ) {
		ERROR("Failed to allocate peak list\n");
		close_hdf5(fh);
		return NULL;
	}

	for ( i=0; i<size[0]; i++ ) {

		float fs, ss, val;
		int pn;

		fs = buf[tw*i+0] + peak_offset;
		ss = buf[tw*i+1] + peak_offset;
		val = buf[tw*i+2];

		if ( data_template_slabby_file_to_panel_coords(dtempl,
		                                               &fs, &ss, &pn) )
		{
			ERROR("Failed to convert %i,%i to "
			      "panel-relative coordinates\n", fs, ss);
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
		}

	}

	free(buf);
	close_hdf5(fh);

	return features;
}


/* This could be extended, later, to include patterns other than just
 * a literal string (no placeholders) and just %.  However, pattern
 * matching is in general not that easy. */
static char *matches_pattern(const char *name, const char *pattern,
                             const char *ev_str_old)
{
	if ( strcmp(pattern, "%") == 0 ) {
		char *nstr = malloc(strlen(ev_str_old)+strlen(name)+2);
		if ( nstr == NULL ) {
			ERROR("Couldn't allocate memory\n");
			return NULL;
		}
		strcpy(nstr, ev_str_old);
		strcat(nstr, "/");
		strcat(nstr, name);
		return nstr;
	} else {
		if ( strcmp(name, pattern) == 0 ) {
			return strdup(ev_str_old);
		} else {
			return NULL;
		}
	}
}


/* Private structure, just to avoid passing char *** around */
struct ev_list
{
	char **events;
	int n_events;
	int max_events;
};


static int add_ev_to_list(struct ev_list *list, char *ev_str)
{
	if ( list->n_events == list->max_events ) {
		char **new_events = realloc(list->events,
		                            (list->max_events+128)*sizeof(char *));
		if ( new_events == NULL ) return 1;
		list->max_events += 128;
		list->events = new_events;
	}

	list->events[list->n_events++] = strdup(ev_str);

	return 0;
}


static char *demunge_event(const char *orig)
{
	size_t len = strlen(orig);
	char *slash;

	if ( len == 0 ) return strdup("//");

	slash = malloc(len+3);
	if ( slash == NULL ) return NULL;
	strcpy(slash, orig+1);
	strcat(slash, "//");
	return slash;
}


static int rec_expand_paths(hid_t gh, struct ev_list *list,
                            const char *ev_str,
                            char **pattern_bits, int n_pattern_bits)
{
	int i;
	H5G_info_t group_info;

	if ( H5Gget_info(gh, &group_info) < 0 ) {
		ERROR("Couldn't get group info\n");
		return 1;
	}

	for ( i=0; i<group_info.nlinks; i++ ) {

		ssize_t size;
		char *name;
		H5O_info_t obj_info;
		char *ev_str_new;

		size = H5Lget_name_by_idx(gh, ".", H5_INDEX_NAME,
		                          H5_ITER_INC, i, NULL, 0,
		                          H5P_DEFAULT);
		if ( (size < 0) || (size > 20000) ) {
			ERROR("Couldn't get link name\n");
			return 1;
		}

		name = malloc(size+1);
		if ( name == NULL ) {
			ERROR("Couldn't allocate memory\n");
			return 1;
		}

		if ( H5Lget_name_by_idx(gh, ".", H5_INDEX_NAME,
		                        H5_ITER_INC, i, name, size+1,
		                        H5P_DEFAULT) < 0 )
		{
			ERROR("Couldn't get name\n");
			return 1;
		}

		ev_str_new = matches_pattern(name, pattern_bits[0],
		                             ev_str);
		if ( ev_str_new == NULL ) {
			free(name);
			continue;
		}

		if ( H5Oget_info_by_idx(gh, ".", H5_INDEX_NAME,
		                        H5_ITER_INC, i, &obj_info, 0) < 0 )
		{
			ERROR("Couldn't get info\n");
			free(name);
			free(ev_str_new);
			return 1;
		}

		if ( obj_info.type == H5O_TYPE_GROUP ) {

			hid_t child_gh;

			if ( n_pattern_bits == 1 ) {
				ERROR("Pattern doesn't match file"
				      " (too short)\n");
				free(name);
				free(ev_str_new);
				return 1;
			}

			child_gh = H5Gopen1(gh, name);
			if ( child_gh < 0 ) {
				ERROR("Couldn't open '%s'\n", name);
				free(name);
				free(ev_str_new);
				return 1;
			}

			if ( rec_expand_paths(child_gh, list,
			                      ev_str_new,
			                      &pattern_bits[1],
			                      n_pattern_bits - 1) )
			{
				free(name);
				free(ev_str_new);
				return 1;
			}

			free(ev_str_new);
			H5Gclose(child_gh);

		} else if ( obj_info.type == H5O_TYPE_DATASET ) {

			char *addme;

			if ( n_pattern_bits != 1 ) {
				ERROR("Pattern doesn't match file"
				      " (too long by %i)\n",
				      n_pattern_bits);
				free(name);
				free(ev_str_new);
				return 1;
			}

			addme = demunge_event(ev_str_new);
			if ( addme != NULL ) {
				add_ev_to_list(list, addme);
				free(addme);
			}
			free(ev_str_new);

		}

		free(name);

	}

	return 0;
}


/* Not "static" so that ev_enumX can test it.
 * Not part of public API! */
char **expand_paths(hid_t fh, char *pattern, int *n_evs)
{
	int n_sep;
	size_t len;
	char **pattern_bits;
	struct ev_list list;
	int i;
	char *start;

	if ( pattern == NULL ) return NULL;
	if ( pattern[0] != '/' ) return NULL;

	/* Chop up the pattern into path bits */
	len = strlen(pattern);
	n_sep = 0;
	for ( i=0; i<len; i++ ) {
		if ( pattern[i] == '/' ) n_sep++;
	}

	pattern_bits = malloc(n_sep*sizeof(char *));
	if ( pattern_bits == NULL ) return NULL;

	start = pattern+1;
	for ( i=0; i<n_sep; i++ ) {
		char *sep = strchr(start, '/');
		if ( sep == NULL ) {
			sep = start+strlen(start);
		}
		pattern_bits[i] = strndup(start, sep-start);
		if ( pattern_bits[i] == NULL ) return NULL;
		start = sep+1;
	}

	list.n_events = 0;
	list.max_events = 0;
	list.events = NULL;

	rec_expand_paths(fh, &list, "", pattern_bits, n_sep);

	for ( i=0; i<n_sep; i++ ) {
		free(pattern_bits[i]);
	}
	free(pattern_bits);

	*n_evs = list.n_events;
	return list.events;
}


static int rec_expand_dims(struct ev_list *list,
                           int *placeholder_sizes,
                           int n_placeholder_dims,
                           char *path_ev)
{
	int i;
	char *dim_ev;
	size_t len;

	len = strlen(path_ev);
	dim_ev = malloc(len+16);
	if ( dim_ev == NULL ) return 1;

	if ( n_placeholder_dims == 1 ) {
		for ( i=0; i<placeholder_sizes[0]; i++ ) {
			snprintf(dim_ev, 16, "%s/%i", path_ev, i);
			if ( add_ev_to_list(list, dim_ev) ) return 1;
		}
	} else {

		for ( i=0; i<placeholder_sizes[0]; i++ ) {
			snprintf(dim_ev, 16, "%s/%i", path_ev, i);
			if ( rec_expand_dims(list,
			                     &placeholder_sizes[1],
			                     n_placeholder_dims - 1,
			                     dim_ev) ) return 1;
		}

	}

	free(dim_ev);
	return 0;
}


static char **expand_dims(int *placeholder_sizes,
                          int n_placeholder_dims,
                          char *path_ev,
                          int *n_evs)
{
	struct ev_list list;

	list.n_events = 0;
	list.max_events = 0;
	list.events = NULL;

	if ( rec_expand_dims(&list, placeholder_sizes,
	                     n_placeholder_dims, path_ev) )
	{
		*n_evs = 0;
		return NULL;
	}

	*n_evs = list.n_events;
	return list.events;
}


static int n_dims_expected(struct panel_template *p)
{
	int i;
	int n_dims = 0;
	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( p->dims[i] != DIM_UNDEFINED ) n_dims++;
	}
	return n_dims;
}


char **image_hdf5_expand_frames(const DataTemplate *dtempl,
                                const char *filename,
                                int *pn_frames)
{
	char **path_evs;
	int n_path_evs;
	hid_t fh;
	int i;
	int dims_expected;
	struct ev_list full_evs;

	if ( dtempl->n_panels == 0 ) return NULL;

	full_evs.events = NULL;
	full_evs.n_events = 0;
	full_evs.max_events = 0;

	/* If the DataTemplate already says that one frame will be
	 * found per file, short-circuit this whole affair */
	if ( (imh_num_placeholders(&dtempl->panels[0]) == 0)
	  && (imh_num_path_placeholders(dtempl->panels[0].data) == 0) )
	{
		add_ev_to_list(&full_evs, "//");
		*pn_frames = full_evs.n_events;
		return full_evs.events;
	}

	if ( !file_exists(filename) ) {
		ERROR("File not found: %s (hdf5_expand_frames)\n", filename);
		return NULL;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file (hdf5_expand_frames): %s\n", filename);
		return NULL;
	}

	/* First, expand placeholders in the HDF5 paths.
	 *
	 * Since we require the number of placeholders to be the same
	 * for all panels, and the placeholders will be substituted
	 * with the same values for each panel (since they come from
	 * the same event ID), this only needs to be done for the
	 * first panel. */
	path_evs = expand_paths(fh, dtempl->panels[0].data,
	                        &n_path_evs);
	if ( path_evs == NULL ) {
		ERROR("Failed to enumerate paths.\n");
		close_hdf5(fh);
		return NULL;
	}

	dims_expected = n_dims_expected(&dtempl->panels[0]);

	/* For each expanded path, enumerate the placeholder
	 * dimensions.  Once again, since the number of placeholders
	 * must be the same for each panel, and the substituted values
	 * will be the same, this only needs to be done for one panel.
	 */
	for ( i=0; i<n_path_evs; i++ ) {

		hid_t dh, sh;
		char *path;
		hsize_t *size;
		int dims;
		int *placeholder_sizes;
		int n_placeholder_dims;
		int j;
		struct panel_template *p = &dtempl->panels[0];

		path = substitute_path(path_evs[i], p->data, 0);
		if ( path == NULL ) {
			ERROR("Path substitution failed during "
			      "expansion of '%s' with partial event "
			      "ID '%s'\n",
			      p->data, path_evs[i]);
			return NULL;
		}

		dh = H5Dopen2(fh, path, H5P_DEFAULT);
		if ( dh < 0 ) {
			ERROR("Error opening '%s'\n", path);
			ERROR("Failed to enumerate events.  "
			      "Check your geometry file.\n");
			close_hdf5(fh);
			return NULL;
		}

		sh = H5Dget_space(dh);
		dims = H5Sget_simple_extent_ndims(sh);
		if ( dims != dims_expected ) {
			ERROR("Unexpected number of dimensions"
			      "(%s has %i, expected %i)\n",
			      path, dims, dims_expected);
			close_hdf5(fh);
			return NULL;
		}

		size = malloc(dims*sizeof(hsize_t));
		placeholder_sizes = malloc(dims*sizeof(int));
		if ( (size == NULL) || (placeholder_sizes == NULL) ) {
			ERROR("Failed to allocate dimensions\n");
			close_hdf5(fh);
			return NULL;
		}

		if ( H5Sget_simple_extent_dims(sh, size, NULL) < 0 ) {
			ERROR("Failed to get size\n");
			close_hdf5(fh);
			return NULL;
		}

		n_placeholder_dims = 0;
		for ( j=0; j<dims; j++ ) {
			if ( p->dims[j] == DIM_PLACEHOLDER ) {
				placeholder_sizes[n_placeholder_dims++] = size[j];
			}
		}
		free(size);

		/* Path event ID ends with //, but expand_dims will
		 * add a slash.  So, remove one slash */
		if ( n_placeholder_dims > 0 ) {

			char **evs_this_path;
			int n_evs_this_path;

			path_evs[i][strlen(path_evs[i])-1] = '\0';
			evs_this_path = expand_dims(placeholder_sizes,
			                            n_placeholder_dims,
			                            path_evs[i],
			                            &n_evs_this_path);

			for ( j=0; j<n_evs_this_path; j++ ) {
				add_ev_to_list(&full_evs, evs_this_path[j]);
				free(evs_this_path[j]);
			}

			free(evs_this_path);

		} else {

			/* Easy case with no dims to expand */
			add_ev_to_list(&full_evs, path_evs[i]);

		}

		free(placeholder_sizes);
		free(path);
		free(path_evs[i]);

	}

	close_hdf5(fh);
	free(path_evs);
	*pn_frames = full_evs.n_events;
	return full_evs.events;
}


int is_hdf5_file(const char *filename)
{
	const char *ext = NULL;

	ext = filename_extension(filename, NULL);
	if ( ext == NULL ) return 0;

	return ( (strcmp(ext, ".h5") == 0)
	      || (strcmp(ext, ".cxi") == 0)
	      || (strcmp(ext, ".nx5") == 0) );
}


/***************************** Writing *****************************/

struct hdf5_write_location {

	const char      *location;
	int              n_panels;
	int             *panel_idxs;

	int              max_ss;
	int              max_fs;

};


static void add_panel_to_location(struct hdf5_write_location *loc,
                                  struct panel_template *p, int pi)
{
	int *new_panel_idxs;

	new_panel_idxs = realloc(loc->panel_idxs,
	                         (loc->n_panels+1)*sizeof(int));
	if ( new_panel_idxs == NULL ) {
		ERROR("Error while managing write location list.\n");
		return;
	}
	loc->panel_idxs = new_panel_idxs;
	loc->panel_idxs[loc->n_panels] = pi;
	loc->n_panels += 1;
	if ( p->orig_max_fs > loc->max_fs ) {
		loc->max_fs = p->orig_max_fs;
	}
	if ( p->orig_max_ss > loc->max_ss ) {
		loc->max_ss = p->orig_max_ss;
	}
}


static struct hdf5_write_location *add_panel_location(struct panel_template *p,
                                                      const char *p_location, int pi,
                                                      struct hdf5_write_location *locations,
                                                      int *pnum_locations)
{
	int li;
	int num_locations = *pnum_locations;
	int done = 0;

	/* Does this HDF5 path already exist in the location list?
	 * If so, add the new panel to it (with a unique index, we hope) */
	for ( li=0; li<num_locations; li++ ) {
		if ( strcmp(p_location, locations[li].location) == 0 ) {
			add_panel_to_location(&locations[li], p, pi);
			done = 1;
		}
	}

	/* If not, add a new location to ths list */
	if ( !done ) {

		struct hdf5_write_location *new_locations;
		size_t nsz;

		nsz = (num_locations+1)*sizeof(struct hdf5_write_location);
		new_locations = realloc(locations, nsz);
		if ( new_locations == NULL ) {
			ERROR("Failed to grow location list.\n");
			return NULL;
		}
		locations = new_locations;

		locations[num_locations].max_ss = p->orig_max_ss;
		locations[num_locations].max_fs = p->orig_max_fs;
		locations[num_locations].location = p_location;
		locations[num_locations].panel_idxs = malloc(sizeof(int));
		if ( locations[num_locations].panel_idxs == NULL ) {
			ERROR("Failed to allocate single idx (!)\n");
			return NULL;
		}
		locations[num_locations].panel_idxs[0] = pi;
		locations[num_locations].n_panels = 1;

		num_locations += 1;

	}

	*pnum_locations = num_locations;
	return locations;
}


static struct hdf5_write_location *make_location_list(const DataTemplate *dtempl,
                                                      int *pnum_locations)
{
	int pi;
	struct hdf5_write_location *locations = NULL;
	int num_locations = 0;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		const char *p_location;

		p = &dtempl->panels[pi];

		assert(p->data != NULL);
		p_location = p->data;

		locations = add_panel_location(p, p_location, pi,
		                               locations,
		                               &num_locations);
		if ( locations == NULL ) {
			ERROR("Failed to add location for panel %s\n",
			      p->name);
			return NULL;
		}

	}

	*pnum_locations = num_locations;
	return locations;
}


static void write_location(hid_t fh, const DataTemplate *dtempl,
                           float **dp,
                           struct hdf5_write_location *loc)
{
	hid_t sh, dh, ph;
	hid_t dh_dataspace;
	hsize_t size[2];
	int pi;

	/* Note the "swap" here, according to section 3.2.5,
	 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
	size[0] = loc->max_ss+1;
	size[1] = loc->max_fs+1;
	sh = H5Screate_simple(2, size, NULL);

	ph = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(ph, 1);

	dh = H5Dcreate2(fh, loc->location, H5T_NATIVE_FLOAT, sh,
	                ph, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
		return;
	}

	H5Sget_simple_extent_dims(sh, size, NULL);

	for ( pi=0; pi<loc->n_panels; pi++ ) {

		hsize_t f_offset[2], f_count[2], dims[2];
		hid_t memspace;
		struct panel_template *p;
		int r;

		p = &dtempl->panels[loc->panel_idxs[pi]];

		f_offset[0] = p->orig_min_ss;
		f_offset[1] = p->orig_min_fs;
		f_count[0] = p->orig_max_ss - p->orig_min_ss +1;
		f_count[1] = p->orig_max_fs - p->orig_min_fs +1;

		dh_dataspace = H5Dget_space(dh);
		r = H5Sselect_hyperslab(dh_dataspace, H5S_SELECT_SET,
		                        f_offset, NULL, f_count, NULL);
		if ( r < 0 ) {
			ERROR("Error selecting file dataspace "
			      "for panel %s\n", p->name);
			H5Pclose(ph);
			H5Dclose(dh);
			H5Sclose(dh_dataspace);
			H5Sclose(sh);
			H5Fclose(fh);
			return;
		}

		dims[0] = PANEL_HEIGHT(p);
		dims[1] = PANEL_WIDTH(p);
		memspace = H5Screate_simple(2, dims, NULL);

		r = H5Dwrite(dh, H5T_NATIVE_FLOAT, memspace, dh_dataspace,
		             H5P_DEFAULT, dp[loc->panel_idxs[pi]]);
		if ( r < 0 ) {
			ERROR("Couldn't write data\n");
			H5Pclose(ph);
			H5Dclose(dh);
			H5Sclose(dh_dataspace);
			H5Sclose(memspace);
			H5Sclose(sh);
			H5Fclose(fh);
			return;
		}

		H5Sclose(dh_dataspace);
		H5Sclose(memspace);
	}
	H5Pclose(ph);
	H5Sclose(sh);
	H5Dclose(dh);
}


static void write_wavelength(hid_t fh, double wl,
                             const DataTemplate *dtempl)
{
	hid_t ph, sh, dh;
	hsize_t size1d[1];
	int r;

	ph = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(ph, 1);

	size1d[0] = 1;
	sh = H5Screate_simple(1, size1d, NULL);

	dh = H5Dcreate2(fh, "/wavelength", H5T_NATIVE_DOUBLE, sh,
	                ph, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset for photon energy.\n");
		return;
	}
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
	             H5P_DEFAULT, &wl);
	if ( r < 0 ) {
		ERROR("Couldn't write photon energy.\n");
		/* carry on */
	}

	H5Pclose(ph);
	H5Dclose(dh);
}


static void write_spectrum(hid_t fh, Spectrum *s,
                           const DataTemplate *dtempl)
{
	herr_t r;
	double *arr;
	int i;
	hid_t sh, dh, ph;
	double kmin, kmax, step;
	const hsize_t n = 1024;

	ph = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(ph, 1);

	arr = malloc(n*sizeof(double));
	if ( arr == NULL ) {
		ERROR("Failed to allocate memory for spectrum.\n");
		return;
	}

	/* Save the wavelength values */
	spectrum_get_range(s, &kmin, &kmax);
	step = (kmax-kmin)/n;
	for ( i=0; i<n; i++ ) {
		arr[i] = 1.0e10/(kmin+i*step);
	}

	sh = H5Screate_simple(1, &n, NULL);

	dh = H5Dcreate2(fh, "/spectrum/wavelengths_A", H5T_NATIVE_DOUBLE,
	                sh, ph, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Failed to create dataset for spectrum wavelengths.\n");
		return;
	}
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
		     H5S_ALL, H5P_DEFAULT, arr);
	if ( r < 0 ) {
		ERROR("Failed to write spectrum wavelengths.\n");
		return;
	}
	H5Dclose(dh);

	/* Save the probability density values */
	for ( i=0; i<n; i++ ) {
		arr[i] = spectrum_get_density_at_k(s, kmin+i*step);
	}

	dh = H5Dcreate2(fh, "/spectrum/pdf", H5T_NATIVE_DOUBLE, sh,
		        H5P_DEFAULT, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Failed to create dataset for spectrum p.d.f.\n");
		return;
	}
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL,
		     H5S_ALL, H5P_DEFAULT, arr);
	if ( r < 0 ) {
		ERROR("Failed to write spectrum p.d.f.\n");
		return;
	}

	H5Dclose(dh);
	H5Pclose(ph);
	free(arr);
}


int image_hdf5_write(const struct image *image,
                     const DataTemplate *dtempl,
                     const char *filename)
{
	hid_t fh;
	int li;
	struct hdf5_write_location *locations;
	int num_locations;

	if ( dtempl == NULL ) return 1;

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return 1;
	}

	locations = make_location_list(dtempl, &num_locations);
	if ( locations == NULL ) {
		ERROR("Failed to create location list\n");
		H5Fclose(fh);
		return 1;
	}

	for ( li=0; li<num_locations; li++ ) {
		write_location(fh, dtempl, image->dp, &locations[li]);
	}

	write_wavelength(fh, image->lambda, dtempl);

	if ( image->spectrum != NULL ) {
		write_spectrum(fh, image->spectrum, dtempl);
	}

	H5Fclose(fh);
	for ( li=0; li<num_locations; li ++ ) {
		free(locations[li].panel_idxs);
	}
	free(locations);
	return 0;
}

#endif  /* HAVE_HDF5 */
