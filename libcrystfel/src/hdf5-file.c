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
 *   2014      Valerio Mariani
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
#include <unistd.h>

#include "events.h"
#include "image.h"
#include "hdf5-file.h"
#include "utils.h"


struct hdf5_write_location {

	char            *location;
	int              n_panels;
	int             *panel_idxs;

	int              max_ss;
	int              max_fs;

};


int split_group_and_object(char *path, char **group, char **object)
{
	char 		*sep;
	char 		*store;

	sep = path;
	store = sep;
	sep = strpbrk(sep + 1, "/");
	if ( sep != NULL ) {
		while ( 1 ) {
			store = sep;
			sep = strpbrk(sep + 1, "/");
			if ( sep == NULL ) {
				break;
			}
		}
	}
	if ( store == path ) {
		*group = NULL;
		*object = strdup(path);
	} else {
		*group = strndup(path, store - path);
		*object = strdup(store+1);
	}
	return 0;
};


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

	if ( access( filename, R_OK ) == -1 ) {
		ERROR("File does not exists or it cannot be read: %s\n", filename);
		free(f);
		return NULL;
	}

	f->fh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( f->fh < 0 ) {
		ERROR("Couldn't open file: %s\n", filename);
		free(f);
		return NULL;
	}

	f->data_open = 0;
	return f;
}


int hdfile_set_image(struct hdfile *f, const char *path,
                     struct panel *p)
{
	hsize_t *size;
	hsize_t *max_size;
	hid_t sh;
	int sh_dim;
	int di;

	f->dh = H5Dopen2(f->fh, path, H5P_DEFAULT);
	if ( f->dh < 0 ) {
		ERROR("Couldn't open dataset\n");
		return -1;
	}
	f->data_open = 1;
	sh = H5Dget_space(f->dh);
	sh_dim = H5Sget_simple_extent_ndims(sh);

	if ( p == NULL ) {

		if ( sh_dim != 2 ) {
			ERROR("Dataset is not two-dimensional\n");
			return -1;
		}

	} else {

		if ( sh_dim != p->dim_structure->num_dims ) {
			ERROR("Dataset dimensionality does not match geometry file\n");
			return -1;
		}

	}

	size = malloc(sh_dim*sizeof(hsize_t));
	max_size = malloc(sh_dim*sizeof(hsize_t));

	H5Sget_simple_extent_dims(sh, size, max_size);
	H5Sclose(sh);

	if ( p == NULL ) {

		f->nx = size[0];
		f->ny = size[1];

	} else {

		for (di=0; di<p->dim_structure->num_dims; di++ ) {

			if (p->dim_structure->dims[di] == HYSL_SS ) {
				f->ny = size[di];
			}
			if (p->dim_structure->dims[di] == HYSL_FS ) {
				f->nx = size[di];
			}

		}

	}

	free(size);
	free(max_size);

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


int hdf5_write_image(const char *filename, struct image *image, char *element)
{
	double lambda, eV;
	double *arr;
	hsize_t size1d[1];
	herr_t r;
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	int i, pi, li;
	char *default_location;
	struct hdf5_write_location *locations;
	struct hdf5_write_location *new_location;
	int num_locations;

	if ( image->det == NULL ) {
		ERROR("Geometry not available\n");
		return 1;
	}

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return 1;
	}

	if ( element != NULL) {
		default_location = strdup(element);
	} else {
		default_location = strdup("/data/data");
	}

	locations = malloc(sizeof(struct hdf5_write_location));
	if ( locations  == NULL ) {
		ERROR("Couldn't create write location list for file: %s\n",
		       filename);
		return 1;
	}
	locations[0].max_ss = 0;
	locations[0].max_fs = 0;
	if ( image->det->panels[0].data != NULL ) {
		locations[0].location = image->det->panels[0].data;
	} else {
		locations[0].location = default_location;
	}
	locations[0].panel_idxs = NULL;
	locations[0].n_panels = 0;
	num_locations = 1;

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		struct panel p;
		int li;
		int panel_processed;
		char *p_location;

		p = image->det->panels[pi];

		if ( p.data == NULL ) {
			p_location = default_location;
		} else {
			p_location = p.data;
		}

		panel_processed = 0;

		for ( li=0; li<num_locations; li++ ) {

			if ( strcmp(p_location, locations[li].location) == 0 ) {

				int *new_panel_idxs;

				new_panel_idxs = realloc(locations[li].panel_idxs,
				                        (locations[li].n_panels+1)*sizeof(int));
				if ( new_panel_idxs == NULL ) {
					ERROR("Error while managing write location list for"
					      "file: %s\n", filename);
					return 1;
				}
				locations[li].panel_idxs = new_panel_idxs;
				locations[li].panel_idxs[locations[li].n_panels] = pi;
				locations[li].n_panels += 1;
				if ( p.orig_max_fs > locations[li].max_fs ) {
					locations[li].max_fs = p.orig_max_fs;
				}
				if ( p.orig_max_ss > locations[li].max_ss ) {
					locations[li].max_ss = p.orig_max_ss;
				}
				panel_processed = 1;
			}
		}

		if ( panel_processed == 0) {

			struct hdf5_write_location *new_locations;
			new_locations = realloc(locations,
			                       (num_locations+1)*
			                       sizeof(struct hdf5_write_location));
			if ( new_locations == NULL ) {
				ERROR("Error while managing write location list for "
				       "file: %s\n",
				       filename);
				return 1;
			}
			locations = new_locations;
			new_location = &locations[num_locations];
			new_location = malloc(sizeof(struct hdf5_write_location));
			if ( new_location == NULL ) {
				ERROR("Error while managing write location list for "
				      "file: %s\n", filename);
			      return 1;
			}
			locations[num_locations].max_ss = p.orig_max_ss;
			locations[num_locations].max_fs = p.orig_max_fs;
			locations[num_locations].location = p_location;
			locations[num_locations].panel_idxs = malloc(sizeof(int));
			if ( locations[num_locations].panel_idxs == NULL ) {
				ERROR("Error while managing write location list for "
				      "file: %s\n", filename);
			      return 1;
			}
			locations[num_locations].panel_idxs[0] = pi;
			locations[num_locations].n_panels = 1;

			num_locations += 1;
		}

	}

	for ( li=0; li<num_locations; li++ ) {

		hid_t ph;
		hid_t dh_dataspace;
		hsize_t size[2];

		char *path, *group =  NULL, *object = NULL;

		path = locations[li].location;
		if ( split_group_and_object(path, &group, &object) ) {
			ERROR("Error while determining write locations "
			      "for file: %s\n", filename);
			return 1;
		}

		if ( group != NULL ) {

			htri_t exists;
			hid_t gph = H5Pcreate(H5P_LINK_CREATE);
			H5Pset_create_intermediate_group(gph, 1);

			exists = H5Lexists(fh, group, H5P_DEFAULT);

			if ( !exists ) {

				gh = H5Gcreate2(fh, group, gph, H5P_DEFAULT,
				                H5P_DEFAULT);
				H5Pclose(gph);
				if ( gh < 0 ) {
					ERROR("Couldn't create group\n");
					H5Fclose(fh);
					return 1;
				}

			} else {
				gh = H5Gopen2(fh, group, H5P_DEFAULT);
			}

		} else {
			gh = -1;
		}

		/* Note the "swap" here, according to section 3.2.5,
		 * "C versus Fortran Dataspaces", of the HDF5 user's guide. */
		size[0] = locations[li].max_ss+1;
		size[1] = locations[li].max_fs+1;
		sh = H5Screate_simple(2, size, NULL);

		/* Set compression */
		ph = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(ph, 2, size);
		H5Pset_deflate(ph, 3);

		if ( group != NULL ) {
			dh = H5Dcreate2(gh, object, H5T_NATIVE_FLOAT, sh,
			                H5P_DEFAULT, ph, H5P_DEFAULT);
		} else {
			dh = H5Dcreate2(fh, object, H5T_NATIVE_FLOAT, sh,
			                H5P_DEFAULT, ph, H5P_DEFAULT);
		}

		if ( dh < 0 ) {
			ERROR("Couldn't create dataset\n");
			H5Fclose(fh);
			return 1;
		}

		/* Muppet check */
		H5Sget_simple_extent_dims(sh, size, NULL);

		for ( pi=0; pi<locations[li].n_panels; pi ++ ) {

			hsize_t f_offset[2], f_count[2];
			hsize_t m_offset[2], m_count[2];
			hsize_t dimsm[2];
			hid_t memspace;
			struct panel p;
			int check;

			p = image->det->panels[locations[li].panel_idxs[pi]];

			f_offset[0] = p.orig_min_ss;
			f_offset[1] = p.orig_min_fs;
			f_count[0] = p.orig_max_ss - p.orig_min_ss +1;
			f_count[1] = p.orig_max_fs - p.orig_min_fs +1;

			dh_dataspace = H5Dget_space(dh);
			check = H5Sselect_hyperslab(dh_dataspace,
			                            H5S_SELECT_SET,
			                            f_offset, NULL,
			                            f_count, NULL);
			if ( check <0 ) {
				ERROR("Error selecting file dataspace "
				      "for panel %s\n", p.name);
				free(group);
				free(object);
				H5Pclose(ph);
				H5Dclose(dh);
				H5Sclose(dh_dataspace);
				H5Sclose(sh);
				if ( gh != -1 ) H5Gclose(gh);
				H5Fclose(fh);
				return 1;
			}

			m_offset[0] = p.min_ss;
			m_offset[1] = p.min_fs;
			m_count[0] = p.max_ss - p.min_ss +1;
			m_count[1] = p.max_fs - p.min_fs +1;
			dimsm[0] = image->height;
			dimsm[1] = image->width;
			memspace = H5Screate_simple(2, dimsm, NULL);
			check = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			                            m_offset, NULL,
			                            m_count, NULL);

			r = H5Dwrite(dh, H5T_NATIVE_FLOAT, memspace,
			             dh_dataspace, H5P_DEFAULT, image->data);
			if ( r < 0 ) {
				ERROR("Couldn't write data\n");
				free(group);
				free(object);
				H5Pclose(ph);
				H5Dclose(dh);
				H5Sclose(dh_dataspace);
				H5Sclose(sh);
				H5Sclose(memspace);
				if ( gh != -1 ) H5Gclose(gh);
				H5Fclose(fh);
				return 1;
			}

			H5Sclose(dh_dataspace);
			H5Sclose(memspace);
		}

		free(group);
		free(object);
		H5Pclose(ph);
		H5Sclose(sh);
		H5Dclose(dh);
		if ( gh != -1 ) H5Gclose(gh);

	}

	gh = H5Gcreate2(fh, "LCLS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
		return 1;
	}

	size1d[0] = 1;
	sh = H5Screate_simple(1, size1d, NULL);

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

		size1d[0] = image->spectrum_size;
		sh = H5Screate_simple(1, size1d, NULL);

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

		size1d[0] = 1;
		sh = H5Screate_simple(1, size1d, NULL);

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
	H5Fclose(fh);
	free(default_location);
	for ( li=0; li<num_locations; li ++ ) {
		free(locations[li].panel_idxs);
	}
	free(locations);

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


int hdf5_read(struct hdfile *f, struct image *image,
			  const char *element, int satcorr)
{
	herr_t r;
	float *buf;
	int fail;

	if ( element == NULL ) {
		fail = hdfile_set_first_image(f, "/");
	} else {
		fail = hdfile_set_image(f, element, NULL);
	}

	if ( fail ) {
		ERROR("Couldn't select path\n");
		return 1;
	}

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

	if ( satcorr ) debodge_saturation(f, image);

	if ( image->beam != NULL ) {

		fill_in_beam_parameters(image->beam, f, NULL);
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


int hdf5_read2(struct hdfile *f, struct image *image,
               struct event *ev, int satcorr)
{
	herr_t r;
	float *buf;
	uint16_t *flags;
	int sum_p_h;
	int p_w;
	int m_min_fs, curr_ss, m_max_fs;
	int mask_is_present;
	int no_mask_loaded;
	int pi;
	hid_t mask_dh = 0;

	if ( image->det == NULL ) {
		ERROR("Geometry not available\n");
		return 1;
	}

	p_w = image->det->panels[0].w;
	sum_p_h = 0;

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		if ( image->det->panels[pi].w != p_w ) {
			ERROR("Panels must have the same width.");
			return 1;
		}

		if ( image->det->panels[pi].mask != NULL ) mask_is_present = 1;

		sum_p_h += image->det->panels[pi].h;

	}

	buf = malloc(sizeof(float)*p_w*sum_p_h);
	if ( mask_is_present ) {
		flags = calloc(p_w*sum_p_h,sizeof(uint16_t));
	}
	image->width = p_w;
	image->height = sum_p_h;

	m_min_fs = 0;
	m_max_fs = p_w-1;
	curr_ss = 0;
	no_mask_loaded = 1;

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		int data_width, data_height;
		hsize_t *f_offset, *f_count;
		int hsi;
		struct dim_structure *hsd;
		hsize_t m_offset[2], m_count[2];
		hsize_t dimsm[2];
		herr_t check;
		hid_t dataspace, memspace, mask_dataspace;
		int fail;

		struct panel *p;
		p=&image->det->panels[pi];

		if ( ev != NULL ) {

			int exists;
			char *panel_full_path;

			panel_full_path = retrieve_full_path(ev, p->data);

			exists = check_path_existence(f->fh, panel_full_path);
			if ( !exists ) {
				ERROR("Cannot find data for panel %s\n",
				      p->name);
				return 1;
			}

			fail = hdfile_set_image(f, panel_full_path, p);

			free(panel_full_path);

		} else {

			if ( p->data == NULL ) {

				fail = hdfile_set_first_image(f, "/");

			} else {

				int exists;
				exists = check_path_existence(f->fh, p->data);
				if ( !exists ) {
					ERROR("Cannot find data for panel %s\n",
					      p->name);
					return 1;
				}
				fail = hdfile_set_image(f, p->data, p);

			}

		}
		if ( fail ) {
			ERROR("Couldn't select path for panel %s\n",
			      p->name);
			return 1;
		}

		data_width = f->ny;
		data_height = f->nx;

		if ( (data_width < p->w )
		  || (data_height < p->h) )
		{
			ERROR("Data size doesn't match panel geometry size"
			      " - rejecting image.\n");
			ERROR("Panel name: %s.  Data size: %i,%i. "
			      "Geometry size: %i,%i\n",
			      p->name, data_width, data_height, p->w, p->h);
			hdfile_close(f);
			return 1;
		}

		hsd = image->det->panels[pi].dim_structure;

		f_offset = malloc(hsd->num_dims*sizeof(hsize_t));
		f_count = malloc(hsd->num_dims*sizeof(hsize_t));

		for ( hsi=0; hsi<hsd->num_dims; hsi++ ) {

			if ( hsd->dims[hsi] == HYSL_FS ) {
				f_offset[hsi] = p->orig_min_fs;
				f_count[hsi] = p->orig_max_fs - p->orig_min_fs+1;
			} else if ( hsd->dims[hsi] == HYSL_SS ) {
				f_offset[hsi] = p->orig_min_ss;
				f_count[hsi] = p->orig_max_ss - p->orig_min_ss+1;
			} else if (hsd->dims[hsi] == HYSL_PLACEHOLDER ) {
				f_offset[hsi] = ev->dim_entries[0];
				f_count[hsi] = 1;
			} else {
				f_offset[hsi] = hsd->dims[hsi];
				f_count[hsi] = 1;
			}

		}

		dataspace = H5Dget_space(f->dh);
		check = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		                            f_offset, NULL, f_count, NULL);
		if ( check <0 ) {
			ERROR("Error selecting file dataspace for panel %s\n",
			      p->name);
			free(buf);
			return 1;
		}

		m_offset[0] = curr_ss;
		m_offset[1] = 0;
		m_count[0] = p->orig_max_ss - p->orig_min_ss +1;
		m_count[1] = m_max_fs - m_min_fs +1;
		dimsm[0] = sum_p_h;
		dimsm[1] = p_w;
		memspace = H5Screate_simple(2,dimsm,NULL);
		check = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
		                            m_offset, NULL, m_count, NULL);
		if ( check < 0 ) {
			ERROR("Error selecting memory dataspace for panel %s\n",
			      p->name);
			free(buf);
			return 1;
		}

		r = H5Dread(f->dh, H5T_NATIVE_FLOAT, memspace, dataspace,
		            H5P_DEFAULT, buf);
		if ( r < 0 ) {
			ERROR("Couldn't read data for panel %s\n",
			      p->name);
			free(buf);
			return 1;
		}
		H5Dclose(f->dh);
		f->data_open = 0;
		H5Sclose(dataspace);
		H5Sclose(memspace);

		if ( p->mask != NULL ) {

			if ( ev != NULL ) {

				int exists;
				char *mask_full_path;

				mask_full_path = retrieve_full_path(ev, p->mask);

				exists = check_path_existence(f->fh,
				                              mask_full_path);
				if ( !exists ) {
					ERROR("Cannot find flags for panel %s\n",
					      p->name);
					return 1;
				}

				mask_dh = H5Dopen2(f->fh, mask_full_path,
				                   H5P_DEFAULT);

				if ( mask_dh <= 0 ) {
					ERROR("Couldn't open flags for panel %s\n",
					      p->name);
					image->flags = NULL;
				} else {

					mask_dataspace = H5Dget_space(mask_dh);
					check = H5Sselect_hyperslab(mask_dataspace, H5S_SELECT_SET,
					                            f_offset, NULL, f_count, NULL);
					if ( check < 0 ) {
						ERROR("Error selecting mask dataspace for panel %s\n",
						      p->name);
					}
					r = H5Dread(mask_dh, H5T_NATIVE_UINT16, memspace,
					        mask_dataspace, H5P_DEFAULT, flags);
					if ( r < 0 ) {
						ERROR("Couldn't read flags for panel %s\n",
						      p->name);
					} else {
						no_mask_loaded = 0;
					}

					H5Sclose(mask_dataspace);
					H5Dclose(mask_dh);

				}

			} else {

				int exists;
				exists = check_path_existence(f->fh, p->mask);
				if ( !exists ) {
					ERROR("Cannot find flags for panel %s\n",
					      p->name);
					return 1;
				}

				mask_dh = H5Dopen2(f->fh, p->mask, H5P_DEFAULT);
				if ( mask_dh <= 0 ) {
					ERROR("Couldn't open flags for panel %s\n",
					      p->name);
					image->flags = NULL;
				} else {

					mask_dataspace = H5Dget_space(mask_dh);
					check = H5Sselect_hyperslab(mask_dataspace, H5S_SELECT_SET,
					                            f_offset, NULL, f_count, NULL);
					if ( check < 0 ) {
						ERROR("Error selecting mask dataspace for panel %s\n",
						      p->name);
					}
					r = H5Dread(mask_dh, H5T_NATIVE_UINT16, memspace,
					            mask_dataspace, H5P_DEFAULT, flags);
					if ( r < 0 ) {
						ERROR("Couldn't read flags for panel %s\n",
						      p->name);
					} else {
						no_mask_loaded = 0;
					}

					H5Sclose(mask_dataspace);
					H5Dclose(mask_dh);

				}

			}

		}

		p->min_fs = m_min_fs;
		p->max_fs = m_max_fs;
		p->min_ss = curr_ss;
		p->max_ss = curr_ss + p->h-1;
		curr_ss += p->h;

		free(f_offset);
		free(f_count);

	}

	image->data = buf;

	if ( no_mask_loaded ) {
		free(flags);
	} else {
		image->flags = flags;
	}

	if ( satcorr ) debodge_saturation(f, image);

	fill_in_values(image->det, f, ev);

	unpack_panels(image, image->det);

	if ( image->beam != NULL ) {

		fill_in_beam_parameters(image->beam, f, ev);

		image->lambda = ph_en_to_lambda(eV_to_J(image->beam->photon_energy));

		if ( (image->beam->photon_energy < 0.0)
		  || (image->lambda > 1000) ) {
			/* Error message covers a silly value in the beam file
			 * or in the HDF5 file. */
			ERROR("Nonsensical wavelength (%e m or %e eV) value "
			      "for file: %s, event: %s.\n",
			      image->lambda, image->beam->photon_energy,
			      image->filename, get_event_string(image->event));
			hdfile_close(f);
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
	int check;

	check = check_path_existence(f->fh, name);
	if ( check == 0 ) {
		ERROR("No such field '%s'\n", name);
		return 0;
	}

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
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
	int check;

	if ( !hdfile_is_scalar(f, name, 1) ) return 1;

	check = check_path_existence(f->fh, name);
	if ( check == 0 ) {
		ERROR("No such field '%s'\n", name);
		return 1;
	}

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);

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


static int get_ev_based_f_value(struct hdfile *f, const char *name,
                                struct event *ev, double *val)
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
	hsize_t size[3];
	herr_t r;
	herr_t check;
	double buf;
	int check_pe;
	int dim_flag;
	int ndims;
	int i;
	char *subst_name = NULL;

	if ( ev->path_length != 0 ) {
		subst_name = partial_event_substitution(ev, name);
	} else {
		subst_name = strdup(name);
	}

	check_pe = check_path_existence(f->fh, subst_name);
	if ( check_pe == 0 ) {
		ERROR("No such field '%s'\n", subst_name);
		return 1;
	}

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( class != H5T_FLOAT ) {
		ERROR("Not a floating point value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	/* Get the dimensionality.  We have to cope with scalars expressed as
	 * arrays with all dimensions 1, as well as zero-d arrays. */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > 3 ) {
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}
	H5Sget_simple_extent_dims(sh, size, NULL);

	m_offset[0] = 0;
	m_count[0] = 1;
	msdims[0] = 1;
	ms = H5Screate_simple(1,msdims,NULL);

	/* Check that the size in all dimensions is 1 */
	/* or that one of the dimensions has the same */
	/* size as the hyperplane events              */

	dim_flag = 0;

	for ( i=0; i<ndims; i++ ) {
		if ( size[i] != 1 ) {
			if ( i == 0 && size[i] > ev->dim_entries[0] ) {
				dim_flag = 1;
			} else {
				H5Tclose(type);
				H5Dclose(dh);
				return 1;
			}
		}
	}

	if ( dim_flag == 0  ) {

		r = H5Dread(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		            H5P_DEFAULT, &buf);

		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			H5Tclose(type);
			H5Dclose(dh);
			return 1;
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
			return 1;
		}

		ms = H5Screate_simple(1,msdims,NULL);
		check = H5Sselect_hyperslab(ms, H5S_SELECT_SET,
		                            m_offset, NULL, m_count, NULL);
		if ( check < 0 ) {
			ERROR("Error selecting memory dataspace for float value");
			free(f_offset);
			free(f_count);
			return 1;
		}

		r = H5Dread(dh, H5T_NATIVE_DOUBLE, ms, sh,
		            H5P_DEFAULT, &buf);
		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			H5Tclose(type);
			H5Dclose(dh);
			return 1;
		}

	}

	free(subst_name);
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
	int check;

	if ( !hdfile_is_scalar(f, name, 1) ) return 1;

	check = check_path_existence(f->fh, name);
	if ( check == 0 ) {
		ERROR("No such field '%s'\n", name);
		return 1;
	}

	dh = H5Dopen2(f->fh, name, H5P_DEFAULT);
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

double get_ev_based_value(struct hdfile *f, const char *name,
                          struct event *ev)
{
	double val = -1;
	get_ev_based_f_value(f, name, ev, &val);
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
                      FILE *fh, struct event *ev)
{
	int i;

	if ( copyme == NULL ) return;

	for ( i=0; i<copyme->n_fields; i++ ) {

		char *val;
		char *field;

		field = copyme->fields[i];
		val = hdfile_get_string_value(f, field, ev);

		if ( field[0] == '/' ) {
			fprintf(fh, "hdf5%s = %s\n", field, val);
		} else {
			fprintf(fh, "hdf5/%s = %s\n", field, val);
		}

		free(val);

	}
}


char *hdfile_get_string_value(struct hdfile *f, const char *name,
                              struct event *ev)
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
		if ( ev != NULL ) {
			if ( get_ev_based_f_value(f, name, ev, &buf_f) ) goto fail;
		} else {
			if ( get_f_value(f, name, &buf_f) ) goto fail;
		}
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

	H5Gclose(gh);

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
			hdfile_set_image(f, names[i], NULL);
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


struct parse_params {
	struct hdfile *hdfile;
	int path_dim;
	const char *path;
	struct event *curr_event;
	struct event_list *ev_list;
	int top_level;
};


int check_path_existence(hid_t fh, const char *path)
{

	char buffer[256];
	char buffer_full_path[2048];
	herr_t herrt;
	struct H5O_info_t ob_info;
	char *path_copy = strdup(path);
	char *start = path_copy;
	char *sep = NULL;

	strncpy(buffer, "\0",1);
	strncpy(buffer_full_path, "\0", 1);

	if ( strcmp(path_copy, "/" ) == 0 ) {
		return 1;
	}

	do {

		int check;

		sep = strstr(start, "/");

		if ( sep != NULL ) {

			if ( sep == start ) {
				start = sep+1;
				strcat(buffer_full_path, "/");
				continue;
			}

			strncpy(buffer, start, sep-start);
			buffer[sep-start]='\0';
			strcat(buffer_full_path, buffer);

			check = H5Lexists(fh, buffer_full_path, H5P_DEFAULT);
			if ( check == 0 ) {
				return 0;
			} else {
				herrt = H5Oget_info_by_name(fh, buffer_full_path,
				                            &ob_info,
				                            H5P_DEFAULT);
				if ( herrt < 0 ) {
					return -1;
				}
				if ( ob_info.type != H5O_TYPE_GROUP ) {
					return 0;
				}

				start = sep+1;
				strcat(buffer_full_path, "/");

			}

		} else {

			strncpy(buffer, start, strlen(start)+1);
			strcat(buffer_full_path, buffer);

			check = H5Lexists(fh, buffer_full_path, H5P_DEFAULT);
			if ( check == 0 ) {
				return 0;
			}

		}
	} while (sep);

	free(path_copy);
	return 1;

}


static herr_t parse_file_event_structure(hid_t loc_id, char *name,
                                         const H5L_info_t *info,
                                         void *operator_data)

{
	struct parse_params *pp;
	char *substituted_path;
	char *ph_loc;
	char *truncated_path;
	htri_t check;
	herr_t herrt_iterate, herrt_info;
	struct H5O_info_t object_info;
	pp = (struct parse_params *)operator_data;

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
		strncpy(&truncated_path[ph_loc-substituted_path],"\0",1);
	}

	herrt_iterate = 0;
	herrt_info = 0;

	check = check_path_existence(pp->hdfile->fh, truncated_path);
	if ( check == 0 ) {
			pop_path_entry_from_event(pp->curr_event);
			return 0;
	} else {

		herrt_info = H5Oget_info_by_name(pp->hdfile->fh, truncated_path,
                                         &object_info, H5P_DEFAULT);
		if ( herrt_info < 0 ) {
			free(truncated_path);
			free(substituted_path);
			return -1;
		}

		if ( pp->curr_event->path_length == pp->path_dim &&
			object_info.type == H5O_TYPE_DATASET ) {

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

				herrt_iterate = H5Literate_by_name(pp->hdfile->fh,
				      truncated_path, H5_INDEX_NAME,
				      H5_ITER_NATIVE, NULL,
				      (H5L_iterate_t)parse_file_event_structure,
				      (void *)pp, H5P_DEFAULT);
			}
		}
	}

	pop_path_entry_from_event(pp->curr_event);

	free(truncated_path);
	free(substituted_path);

	return herrt_iterate;
}


struct event_list *fill_event_list(struct hdfile *hdfile, struct detector *det)
{
	int pi;
	int evi;
	herr_t check;
	struct event_list *master_el;
	struct event_list *master_el_with_dims;

	master_el = initialize_event_list();

	if ( det->path_dim != 0 ) {

		for ( pi=0; pi<det->n_panels; pi++ ) {

			struct parse_params pparams;
			struct event *empty_event;
			struct event_list *panel_ev_list;
			int ei;

			empty_event = initialize_event();
			panel_ev_list = initialize_event_list();

			pparams.path = det->panels[pi].data;
			pparams.hdfile = hdfile;
			pparams.path_dim = det->path_dim;
			pparams.curr_event = empty_event;
			pparams.top_level = 1;
			pparams.ev_list = panel_ev_list;

			check = parse_file_event_structure(hdfile->fh, NULL,
				                           NULL,
			                                    (void *)&pparams);

			if ( check < 0 ) {
				free_event(empty_event);
				free_event_list(panel_ev_list);
				return NULL;
			}

			for ( ei=0; ei<panel_ev_list->num_events; ei++ ) {

				int fail_add;

				fail_add = add_non_existing_event_to_event_list(master_el,
				                     panel_ev_list->events[ei]);
				if ( fail_add ) {

					free_event(empty_event);
					free_event_list(panel_ev_list);
					return NULL;
				}
			}

			free_event(empty_event);
			free_event_list(panel_ev_list);
		}

	}

	if ( det->dim_dim > 0 ) {

		if ( master_el->num_events == 0 ) {

			struct event *empty_ev;
			empty_ev = initialize_event();
			append_event_to_event_list(master_el, empty_ev);
			free(empty_ev);

		}

		master_el_with_dims = initialize_event_list();

		for (evi=0; evi<master_el->num_events; evi++ ) {

			int global_path_dim = -1;
			int pai;
			int mlwd;

			for ( pai=0; pai<det->n_panels; pai++ ) {

				char *full_panel_path;
				hid_t dh;
				hid_t sh;
				int dims;
				hsize_t *size;
				hsize_t *max_size;
				int hsdi;
				int panel_path_dim = 0;

				full_panel_path = retrieve_full_path(master_el->events[evi],
				                                     det->panels[pai].data);

				dh = H5Dopen2(hdfile->fh, full_panel_path, H5P_DEFAULT);
				sh = H5Dget_space(dh);
				dims = H5Sget_simple_extent_ndims(sh);

				size = malloc(dims*sizeof(hsize_t));
				max_size = malloc(dims*sizeof(hsize_t));

				dims = H5Sget_simple_extent_dims(sh, size, max_size);

				for ( hsdi=0; hsdi<det->panels[pai].dim_structure->num_dims;
				      hsdi++ ) {
					if (det->panels[pai].dim_structure->dims[hsdi] ==
				    HYSL_PLACEHOLDER ) {
						panel_path_dim = size[hsdi];
						break;
					}
				}


				if ( global_path_dim == -1 ) {

					global_path_dim = panel_path_dim;

				} else if ( panel_path_dim != global_path_dim ) {

					ERROR("Data blocks paths for panels must have the same "
					      "number of placeholders");
					free(size);
					free(max_size);
					return NULL;
				}

			}

			for ( mlwd=0; mlwd<global_path_dim; mlwd++ ) {

				struct event *mlwd_ev;

				mlwd_ev = copy_event(master_el->events[evi]);
				push_dim_entry_to_event(mlwd_ev, mlwd);
				append_event_to_event_list(master_el_with_dims,
										mlwd_ev);
				free(mlwd_ev);
			}

		}

		free_event_list(master_el);
		return master_el_with_dims;

	}

	return master_el;
}
