/*
 * image-hdf5.c
 *
 * Image loading, HDF5 parts
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#include <config.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <hdf5.h>
#include <unistd.h>

#include "image.h"
#include "utils.h"
#include "events.h"
#include "detgeom.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"


static int load_hdf5_hyperslab(struct panel_template *p,
                               const char *filename,
                               const char *event,
                               void **pdata,
                               hid_t el_type, size_t el_size,
                               int skip_placeholders_ok,
                               const char *path_spec)
{
	struct event *ev;
	hid_t fh;
	herr_t r;
	hsize_t *f_offset, *f_count;
	hid_t dh;
	int hsi;
	herr_t check;
	hid_t dataspace, memspace;
	hsize_t dims[2];
	char *panel_full_path;
	void *data;
	int ndims, dpos;
	int skip_placeholders = 0;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return 1;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		return 1;
	}

	ev = get_event_from_event_string(event);
	if ( (ev == NULL) && (event != NULL) ) {
		ERROR("Invalid event identifier '%s'\n", event);
		H5Fclose(fh);
		return 1;
	}

	panel_full_path = retrieve_full_path(ev, path_spec);

	dh = H5Dopen2(fh, panel_full_path, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Cannot open data for panel %s (%s)\n",
		      p->name, panel_full_path);
		free(panel_full_path);
		free_event(ev);
		H5Fclose(fh);
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
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}

	if ( ndims != p->dim_structure->num_dims ) {
		/* Dimensionality doesn't match */
		int n_nonplaceholder = 0;
		for ( hsi=0; hsi<p->dim_structure->num_dims; hsi++ ) {
			if ( p->dim_structure->dims[hsi] != HYSL_PLACEHOLDER ) {
				n_nonplaceholder++;
			}
		}

		/* If the dimensions match after excluding
		 * placeholders, it's OK - probably a static mask
		 * in a multi-event file. */
		if ( ndims == n_nonplaceholder ) {
			skip_placeholders = 1;
		}
	}

	f_offset = malloc(ndims*sizeof(hsize_t));
	f_count = malloc(ndims*sizeof(hsize_t));
	if ( (f_offset == NULL) || (f_count == NULL ) ) {
		ERROR("Failed to allocate offset or count.\n");
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}

	dpos = 0;
	for ( hsi=0; hsi<p->dim_structure->num_dims; hsi++ ) {

		switch ( p->dim_structure->dims[hsi] ) {

			case HYSL_FS:
			f_offset[dpos] = p->orig_min_fs;
			f_count[dpos] = p->orig_max_fs - p->orig_min_fs+1;
			dpos++;
			break;

			case HYSL_SS:
			f_offset[dpos] = p->orig_min_ss;
			f_count[dpos] = p->orig_max_ss - p->orig_min_ss+1;
			dpos++;
			break;

			case HYSL_PLACEHOLDER:
			if ( !skip_placeholders ) {
				f_offset[dpos] = ev->dim_entries[0];
				f_count[dpos] = 1;
				dpos++;
			}
			break;

			default:
			f_offset[dpos] = p->dim_structure->dims[hsi];
			f_count[dpos] = 1;
			dpos++;
			break;

		}

	}

	check = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
	                            f_offset, NULL, f_count, NULL);
	if ( check < 0 ) {
		ERROR("Error selecting file dataspace for panel %s\n",
		      p->name);
		free(f_offset);
		free(f_count);
		free_event(ev);
		H5Fclose(fh);
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
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}

	r = H5Dread(dh, el_type, memspace, dataspace,
	            H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't read data for panel %s\n",
		      p->name);
		free(f_offset);
		free(f_count);
		free(data);
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}

	H5Dclose(dh);
	H5Sclose(dataspace);
	free(f_offset);
	free(f_count);
	free_event(ev);
	H5Fclose(fh);

	*pdata = data;
	return 0;
}


struct image *image_hdf5_read(DataTemplate *dtempl,
                              const char *filename, const char *event)
{
	struct image *image;
	int i;

	image = image_new();
	if ( image == NULL ) {
		ERROR("Couldn't allocate image structure.\n");
		return NULL;
	}

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate data array.\n");
		image_free(image);
		return NULL;
	}

	/* Set all pointers to NULL for easier clean-up */
	for ( i=0; i<dtempl->n_panels; i++ ) image->dp[i] = NULL;

	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( load_hdf5_hyperslab(&dtempl->panels[i], filename,
		                         event, (void *)&image->dp[i],
		                         H5T_NATIVE_FLOAT,
		                         sizeof(float), 0,
		                         dtempl->panels[i].data) )
		{
			ERROR("Failed to load panel data\n");
			image_free(image);
			return NULL;
		}
	}

	image->filename = strdup(filename);
	image->ev = safe_strdup(event);

	return image;
}


int image_hdf5_read_mask(struct panel_template *p,
                         const char *filename, const char *event,
                         int *bad, int mask_good, int mask_bad)
{
	int p_w, p_h;
	int *mask = NULL;
	long unsigned int j;

	p_w = p->orig_max_fs - p->orig_min_fs + 1;
	p_h = p->orig_max_ss - p->orig_min_ss + 1;

	if ( load_hdf5_hyperslab(p, filename, event,
	                         (void *)&mask, H5T_NATIVE_INT,
	                         sizeof(int), 1, p->mask) )
	{
		ERROR("Failed to load mask data\n");
		free(mask);
		return 1;
	}

	for ( j=0; j<p_w*p_h; j++ ) {

		/* Bad if it's missing any of the "good" bits */
		if ( (mask[j] & mask_good) != mask_good ) bad[j] = 1;

		/* Bad if it has any of the "bad" bits. */
		if ( mask[j] & mask_bad ) bad[j] = 1;

	}

	free(mask);
	return 0;
}


double image_hdf5_get_value(const char *name, const char *filename,
                            const char *event)
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
	int dim_flag;
	int ndims;
	int i;
	char *subst_name = NULL;
	struct event *ev;
	hid_t fh;
	double val;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n", filename);
		return NAN;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		return NAN;
	}
	ev = get_event_from_event_string(event);
	if ( (ev == NULL) && (event != NULL) ) {
		ERROR("Invalid event identifier '%s'\n", event);
		H5Fclose(fh);
		return NAN;
	}

	subst_name = retrieve_full_path(ev, name);

	dh = H5Dopen2(fh, subst_name, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("No such numeric field '%s'\n", subst_name);
		H5Fclose(fh);
		return NAN;
	}

	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( (class != H5T_FLOAT) && (class != H5T_INTEGER) ) {
		ERROR("Not a floating point or integer value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return NAN;
	}

	/* Get the dimensionality.  We have to cope with scalars expressed as
	 * arrays with all dimensions 1, as well as zero-d arrays. */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > 64 ) {
		ERROR("Too many dimensions for numeric value\n");
		H5Tclose(type);
		H5Dclose(dh);
		return NAN;
	}
	H5Sget_simple_extent_dims(sh, size, NULL);

	m_offset[0] = 0;
	m_count[0] = 1;
	msdims[0] = 1;
	ms = H5Screate_simple(1,msdims,NULL);

	/* Check that the size in all dimensions is 1
	 * or that one of the dimensions has the same
	 * size as the hyperplane events */

	dim_flag = 0;

	for ( i=0; i<ndims; i++ ) {
		if ( size[i] == 1 ) continue;
		if ( ( i==0 ) && (ev != NULL) && (size[i] > ev->dim_entries[0]) ) {
			dim_flag = 1;
		} else {
			H5Tclose(type);
			H5Dclose(dh);
			return NAN;
		}
	}

	if ( dim_flag == 0 ) {

		if ( H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		             H5P_DEFAULT, &val) < 0 )
		{
			ERROR("Couldn't read value.\n");
			H5Tclose(type);
			H5Dclose(dh);
			return NAN;
		}

	} else {

		f_offset = malloc(ndims*sizeof(hsize_t));
		f_count = malloc(ndims*sizeof(hsize_t));

		for ( i=0; i<ndims; i++ ) {

			if ( i == 0 ) {
				f_offset[i] = ev->dim_entries[0];
				f_count[i] = 1;
			} else {
				f_offset[i] = 0;
				f_count[i] = 0;
			}

		}

		check = H5Sselect_hyperslab(sh, H5S_SELECT_SET,
		                            f_offset, NULL, f_count, NULL);
		if ( check <0 ) {
			ERROR("Error selecting dataspace for float value");
			free(f_offset);
			free(f_count);
			return NAN;
		}

		ms = H5Screate_simple(1,msdims,NULL);
		check = H5Sselect_hyperslab(ms, H5S_SELECT_SET,
		                            m_offset, NULL, m_count, NULL);
		if ( check < 0 ) {
			ERROR("Error selecting memory dataspace for float value");
			free(f_offset);
			free(f_count);
			return NAN;
		}

		r = H5Dread(dh, H5T_NATIVE_DOUBLE, ms, sh, H5P_DEFAULT, &val);
		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			H5Tclose(type);
			H5Dclose(dh);
			return NAN;
		}

	}

	free_event(ev);
	free(subst_name);
	H5Fclose(fh);

	return val;
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


static float *read_hdf5_data(hid_t fh, char *path, int line)
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
	if ( line> tw-1 ) {
		H5Sclose(sh);
		H5Dclose(dh);
		ERROR("Data block %s does not contain data for required event.\n",
		      path);
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
	struct event *ev;
	char *subst_name;

	int line = 0;
	int num_peaks;

	float *buf_x;
	float *buf_y;
	float *buf_i;

	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return NULL;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		return NULL;
	}

	ev = get_event_from_event_string(event);
	if ( (ev == NULL) && (event != NULL) ) {
		ERROR("Invalid event identifier '%s'\n", event);
		H5Fclose(fh);
		return NULL;
	}

	if ( ev->dim_entries == NULL ) {
		ERROR("CXI format peak list format selected,"
		      "but file has no event structure");
		return NULL;
	}
	line = ev->dim_entries[0];

	subst_name = retrieve_full_path(ev, dtempl->peak_list);
	free_event(ev);
	if ( subst_name == NULL ) {
		ERROR("Invalid peak path %s\n", subst_name);
		H5Fclose(fh);
		return NULL;
	}

	snprintf(path_n, 1024, "%s/nPeaks", subst_name);
	snprintf(path_x, 1024, "%s/peakXPosRaw", subst_name);
	snprintf(path_y, 1024, "%s/peakYPosRaw", subst_name);
	snprintf(path_i, 1024, "%s/peakTotalIntensity", subst_name);

	r = read_peak_count(fh, path_n, line, &num_peaks);
	if ( r != 0 ) return NULL;

	buf_x = read_hdf5_data(fh, path_x, line);
	if ( r != 0 ) return NULL;

	buf_y = read_hdf5_data(fh, path_y, line);
	if ( r != 0 ) return NULL;

	buf_i = read_hdf5_data(fh, path_i, line);
	if ( r != 0 ) return NULL;

	features = image_feature_list_new();

	for ( pk=0; pk<num_peaks; pk++ ) {

		float fs, ss, val;
		int pn;

		fs = buf_x[pk] + peak_offset;
		ss = buf_y[pk] + peak_offset;
		val = buf_i[pk];

		if ( data_template_file_to_panel_coords(dtempl,
		                                        &fs, &ss,
		                                        &pn) )
		{
			ERROR("Failed to convert %i,%i to "
			      "panel-relative coordinates\n", fs, ss);
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
		}

	}

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
	struct event *ev;
	char *subst_name;
	ImageFeatureList *features;
	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
		return NULL;
	}

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		return NULL;
	}

	ev = get_event_from_event_string(event);
	if ( (ev == NULL) && (event != NULL) ) {
		ERROR("Invalid event identifier '%s'\n", event);
		H5Fclose(fh);
		return NULL;
	}

	subst_name = retrieve_full_path(ev, dtempl->peak_list);
	free_event(ev);
	if ( subst_name == NULL ) {
		ERROR("Invalid peak path %s\n", subst_name);
		free_event(ev);
		H5Fclose(fh);
		return NULL;
	}

	dh = H5Dopen2(fh, subst_name, H5P_DEFAULT);
	free(subst_name);
	if ( dh < 0 ) {
		ERROR("Peak list (%s) not found.\n", subst_name);
		H5Fclose(fh);
		return NULL;
	}

	sh = H5Dget_space(dh);
	if ( sh < 0 ) {
		ERROR("Couldn't get dataspace for peak list.\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}

	if ( H5Sget_simple_extent_ndims(sh) != 2 ) {
		ERROR("Peak list has the wrong dimensionality (%i).\n",
		      H5Sget_simple_extent_ndims(sh));
		H5Sclose(sh);
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}

	H5Sget_simple_extent_dims(sh, size, max_size);
	H5Sclose(sh);

	tw = size[1];
	if ( (tw != 3) && (tw != 4) ) {
		ERROR("Peak list has the wrong dimensions.\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}

	buf = malloc(sizeof(float)*size[0]*size[1]);
	if ( buf == NULL ) {
		ERROR("Couldn't reserve memory for the peak list.\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}
	r = H5Dread(dh, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	            H5P_DEFAULT, buf);
	if ( r < 0 ) {
		ERROR("Couldn't read peak list.\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}

	features = image_feature_list_new();
	if ( features == NULL ) {
		ERROR("Failed to allocate peak list\n");
		H5Dclose(dh);
		H5Fclose(fh);
		return NULL;
	}

	for ( i=0; i<size[0]; i++ ) {

		float fs, ss, val;
		int pn;

		fs = buf[tw*i+0] + peak_offset;
		ss = buf[tw*i+1] + peak_offset;
		val = buf[tw*i+2];

		if ( data_template_file_to_panel_coords(dtempl,
		                                        &fs, &ss,
		                                        &pn) )
		{
			ERROR("Failed to convert %i,%i to "
			      "panel-relative coordinates\n", fs, ss);
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
		}

	}

	free(buf);
	H5Dclose(dh);
	H5Fclose(fh);

	return features;
}


struct parse_params {
	hid_t fh;
	int path_dim;
	const char *path;
	struct event *curr_event;
	struct event_list *ev_list;
	int top_level;
};


static herr_t parse_file_event_structure(hid_t loc_id, char *name,
                                         const H5L_info_t *info,
                                         struct parse_params *pp)

{
	char *substituted_path;
	char *ph_loc;
	char *truncated_path;
	herr_t herrt_iterate, herrt_info;
	struct H5O_info_t object_info;

	if ( !pp->top_level ) {

		int fail_push;

		fail_push = push_path_entry_to_event(pp->curr_event, name);
		if ( fail_push ) {
			return -1;
		}

		substituted_path = event_path_placeholder_subst(name, pp->path);

	} else {
		substituted_path = strdup(pp->path);
	}

	if ( pp->top_level == 1 ) {
		pp->top_level = 0;
	}

	truncated_path = strdup(substituted_path);
	ph_loc = strstr(substituted_path,"%");
	if ( ph_loc != NULL) {
		truncated_path[ph_loc-substituted_path] = '\0';
	}

	herrt_iterate = 0;
	herrt_info = 0;

	herrt_info = H5Oget_info_by_name(pp->fh, truncated_path,
                                 &object_info, H5P_DEFAULT);
	if ( herrt_info < 0 ) {
		free(truncated_path);
		free(substituted_path);
		return -1;
	}

	if ( pp->curr_event->path_length == pp->path_dim
	 &&  object_info.type == H5O_TYPE_DATASET )
	{

		int fail_append;

		fail_append = append_event_to_event_list(pp->ev_list,
		                                         pp->curr_event);
		if ( fail_append ) {
			free(truncated_path);
			free(substituted_path);
			return -1;
		}

		pop_path_entry_from_event(pp->curr_event);
		return 0;

	} else {

		pp->path = substituted_path;

		if ( object_info.type == H5O_TYPE_GROUP ) {

			herrt_iterate = H5Literate_by_name(pp->fh,
			      truncated_path, H5_INDEX_NAME,
			      H5_ITER_NATIVE, NULL,
			      (H5L_iterate_t)parse_file_event_structure,
			      (void *)pp, H5P_DEFAULT);
		}
	}

	pop_path_entry_from_event(pp->curr_event);

	free(truncated_path);
	free(substituted_path);

	return herrt_iterate;
}


static int fill_paths(hid_t fh, const DataTemplate *dtempl, int pi,
                      struct event_list *master_el)
{
	struct parse_params pparams;
	struct event *empty_event;
	struct event_list *panel_ev_list;
	int ei;
	int check;

	empty_event = initialize_event();
	panel_ev_list = initialize_event_list();
	if ( (empty_event == NULL) || (panel_ev_list == NULL) )
	{
		ERROR("Failed to allocate memory for event list.\n");
		return 1;
	}

	pparams.path = dtempl->panels[pi].data;
	pparams.fh = fh;
	pparams.path_dim = dtempl->path_dim;
	pparams.curr_event = empty_event;
	pparams.top_level = 1;
	pparams.ev_list = panel_ev_list;

	check = parse_file_event_structure(fh, NULL, NULL, &pparams);
	if ( check < 0 ) {
		free_event(empty_event);
		free_event_list(panel_ev_list);
		return 1;
	}

	for ( ei=0; ei<panel_ev_list->num_events; ei++ ) {

		int fail_add;

		fail_add = add_non_existing_event_to_event_list(master_el,
		                                     panel_ev_list->events[ei]);
		if ( fail_add ) {
			free_event(empty_event);
			free_event_list(panel_ev_list);
			return 1;
		}

	}

	free_event(empty_event);
	free_event_list(panel_ev_list);

	return 0;
}


static int check_dims(hid_t fh, struct panel_template *p,
                      struct event *ev,
                      struct event_list *events, int *global_path_dim)
{
	char *full_panel_path;
	hid_t dh;
	hid_t sh;
	int dims;
	hsize_t *size;
	hsize_t *max_size;
	int hsdi;
	int panel_path_dim = 0;
	struct dim_structure *panel_dim_structure;

	/* Get the full path for this panel in this event */
	full_panel_path = retrieve_full_path(ev, p->data);

	dh = H5Dopen2(fh, full_panel_path, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Error opening '%s'\n", full_panel_path);
		ERROR("Failed to enumerate events.  "
		      "Check your geometry file.\n");
		return 1;
	}

	sh = H5Dget_space(dh);
	dims = H5Sget_simple_extent_ndims(sh);
	size = malloc(dims*sizeof(hsize_t));
	max_size = malloc(dims*sizeof(hsize_t));
	if ( (size==NULL) || (max_size==NULL) ) {
		ERROR("Failed to allocate memory for dimensions\n");
		return 1;
	}

	dims = H5Sget_simple_extent_dims(sh, size, max_size);

	panel_dim_structure = p->dim_structure;
	for ( hsdi=0; hsdi<panel_dim_structure->num_dims; hsdi++ ) {
		if ( panel_dim_structure->dims[hsdi] == HYSL_PLACEHOLDER ) {
			panel_path_dim = size[hsdi];
			break;
		}
	}

	if ( *global_path_dim == -1 ) {

		*global_path_dim = panel_path_dim;

	} else if ( panel_path_dim != *global_path_dim ) {

		ERROR("All panels must have the same number of frames\n");
		ERROR("Panel %s has %i frames in one dimension, but the first "
		      "panel has %i.\n",
		      p->name, panel_path_dim, *global_path_dim);
		free(size);
		free(max_size);
		return 1;
	}

	H5Sclose(sh);
	H5Dclose(dh);

	return 0;
}


struct event_list *image_hdf5_expand_frames(const DataTemplate *dtempl,
                                            const char *filename)
{
	struct event_list *master_el;
	hid_t fh;

	fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file '%s'\n", filename);
		return NULL;
	}

	master_el = initialize_event_list();
	if ( master_el == NULL ) {
		ERROR("Failed to allocate event list.\n");
		H5Fclose(fh);
		return NULL;
	}

	/* First expand any placeholders in the HDF5 paths */
	if ( dtempl->path_dim != 0 ) {
		int pi;
		for ( pi=0; pi<dtempl->n_panels; pi++ ) {
			if ( fill_paths(fh, dtempl, pi, master_el) ) {
				ERROR("Failed to enumerate paths.\n");
				H5Fclose(fh);
				return NULL;
			}
		}
	}

	/* Now enumerate the placeholder dimensions */
	if ( dtempl->dim_dim > 0 ) {

		struct event_list *master_el_with_dims;
		int evi;

		/* If there were no HDF5 path placeholders, add a dummy event */
		if ( master_el->num_events == 0 ) {
			struct event *empty_ev;
			empty_ev = initialize_event();
			append_event_to_event_list(master_el, empty_ev);
			free(empty_ev);
		}

		master_el_with_dims = initialize_event_list();

		/* For each event so far, expand the dimensions */
		for ( evi=0; evi<master_el->num_events; evi++ ) {

			int pi;
			int global_path_dim = -1;
			int mlwd;

			/* Check the dimensionality of each panel */
			for ( pi=0; pi<dtempl->n_panels; pi++ ) {
				if ( check_dims(fh, &dtempl->panels[pi],
				                master_el->events[evi],
				                master_el_with_dims,
				                &global_path_dim) )
				{
					ERROR("Failed to enumerate dims.\n");
					H5Fclose(fh);
					return NULL;
				}
			}

			/* Add this dimensionality to all events */
			for ( mlwd=0; mlwd<global_path_dim; mlwd++ ) {

				struct event *mlwd_ev;

				mlwd_ev = copy_event(master_el->events[evi]);
				push_dim_entry_to_event(mlwd_ev, mlwd);
				append_event_to_event_list(master_el_with_dims,
				                           mlwd_ev);
				free_event(mlwd_ev);
			}

		}

		free_event_list(master_el);
		H5Fclose(fh);
		return master_el_with_dims;

	} else {

		H5Fclose(fh);
		return master_el;

	}
}


int is_hdf5_file(const char *filename)
{
	return H5Fis_hdf5(filename);
}
