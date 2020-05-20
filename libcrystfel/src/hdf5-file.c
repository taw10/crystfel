/*
 * hdf5-file.c
 *
 * Read/write HDF5 data files
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2016 Thomas White <taw@physics.org>
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

/** \file hdf5-file.h */

struct hdf5_write_location {

	const char      *location;
	int              n_panels;
	int             *panel_idxs;

	int              max_ss;
	int              max_fs;

};


struct hdfile {

	const char      *path;  /* Current data path */

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
		ERROR("File does not exist or cannot be read: %s\n",
		      filename);
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


int hdfile_set_image(struct hdfile *f, const char *path)
{
	f->dh = H5Dopen2(f->fh, path, H5P_DEFAULT);
	if ( f->dh < 0 ) {
		ERROR("Couldn't open dataset\n");
		return -1;
	}
	f->data_open = 1;
	return 0;
}


static void add_panel_to_location(struct hdf5_write_location *loc,
                                  struct panel *p, int pi)
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


static void add_panel_location(struct panel *p, const char *p_location, int pi,
                               struct hdf5_write_location **plocations,
                               int *pnum_locations)
{
	int li;
	int num_locations = *pnum_locations;
	struct hdf5_write_location *locations = *plocations;
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
			return;
		}
		locations = new_locations;

		locations[num_locations].max_ss = p->orig_max_ss;
		locations[num_locations].max_fs = p->orig_max_fs;
		locations[num_locations].location = p_location;
		locations[num_locations].panel_idxs = malloc(sizeof(int));
		if ( locations[num_locations].panel_idxs == NULL ) {
			ERROR("Failed to allocate single idx (!)\n");
			return;
		}
		locations[num_locations].panel_idxs[0] = pi;
		locations[num_locations].n_panels = 1;

		num_locations += 1;

	}

	*plocations = locations;
	*pnum_locations = num_locations;
}


static struct hdf5_write_location *make_location_list(struct detector *det,
                                                      const char *def_location,
                                                      int *pnum_locations)
{
	int pi;
	struct hdf5_write_location *locations = NULL;
	int num_locations = 0;

	for ( pi=0; pi<det->n_panels; pi++ ) {

		struct panel *p;
		const char *p_location;

		p = &det->panels[pi];

		if ( p->data == NULL ) {
			p_location = def_location;
		} else {
			p_location = p->data;
		}

		add_panel_location(p, p_location, pi,
		                   &locations, &num_locations);

	}

	*pnum_locations = num_locations;
	return locations;
}


static void write_location(hid_t fh, struct detector *det, float **dp,
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
		struct panel p;
		int r;

		p = det->panels[loc->panel_idxs[pi]];

		f_offset[0] = p.orig_min_ss;
		f_offset[1] = p.orig_min_fs;
		f_count[0] = p.orig_max_ss - p.orig_min_ss +1;
		f_count[1] = p.orig_max_fs - p.orig_min_fs +1;

		dh_dataspace = H5Dget_space(dh);
		r = H5Sselect_hyperslab(dh_dataspace, H5S_SELECT_SET,
		                        f_offset, NULL, f_count, NULL);
		if ( r < 0 ) {
			ERROR("Error selecting file dataspace "
			      "for panel %s\n", p.name);
			H5Pclose(ph);
			H5Dclose(dh);
			H5Sclose(dh_dataspace);
			H5Sclose(sh);
			H5Fclose(fh);
			return;
		}

		dims[0] = p.h;
		dims[1] = p.w;
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


static void write_photon_energy(hid_t fh, double eV, const char *ph_en_loc)
{
	hid_t ph, sh, dh;
	hsize_t size1d[1];
	int r;

	ph = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(ph, 1);

	size1d[0] = 1;
	sh = H5Screate_simple(1, size1d, NULL);

	dh = H5Dcreate2(fh, ph_en_loc, H5T_NATIVE_DOUBLE, sh,
	                ph, H5S_ALL, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset for photon energy.\n");
		return;
	}
	r = H5Dwrite(dh, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &eV);
	if ( r < 0 ) {
		ERROR("Couldn't write photon energy.\n");
		/* carry on */
	}

	H5Pclose(ph);
	H5Dclose(dh);
}


static void write_spectrum(hid_t fh, Spectrum *s)
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


int hdf5_write_image(const char *filename, const struct image *image,
                     char *element)
{
	hid_t fh;
	int li;
	char *default_location;
	struct hdf5_write_location *locations;
	int num_locations;
	const char *ph_en_loc;

	if ( image->det == NULL ) {
		ERROR("Geometry not available\n");
		return 1;
	}

	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
		return 1;
	}

	if ( element != NULL ) {
		default_location = strdup(element);
	} else {
		default_location = strdup("/data/data");
	}

	locations = make_location_list(image->det, default_location,
	                               &num_locations);

	for ( li=0; li<num_locations; li++ ) {
		write_location(fh, image->det, image->dp, &locations[li]);
	}

	if ( image->beam == NULL
	 || (image->beam != NULL && image->beam->photon_energy_from == NULL) ) {
		ph_en_loc = "photon_energy_eV";
	} else {
		ph_en_loc = image->beam->photon_energy_from;
	}

	write_photon_energy(fh, ph_lambda_to_eV(image->lambda), ph_en_loc);

	if ( image->spectrum != NULL ) {
		write_spectrum(fh, image->spectrum);
	}

	H5Fclose(fh);
	free(default_location);
	for ( li=0; li<num_locations; li ++ ) {
		free(locations[li].panel_idxs);
	}
	free(locations);
	return 0;
}


int hdfile_get_value(struct hdfile *f, const char *name, struct event *ev,
                     void *val, hid_t memtype)
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
	int check_pe;
	int dim_flag;
	int ndims;
	int i;
	char *subst_name = NULL;

	if ( (ev != NULL) && (ev->path_length != 0) ) {
		subst_name = retrieve_full_path(ev, name);
	} else {
		subst_name = strdup(name);
	}

	check_pe = check_path_existence(f->fh, subst_name);
	if ( check_pe == 0 ) {
		ERROR("No such event-based float field '%s'\n", subst_name);
		return 1;
	}

	dh = H5Dopen2(f->fh, subst_name, H5P_DEFAULT);
	type = H5Dget_type(dh);
	class = H5Tget_class(type);

	if ( (class != H5T_FLOAT) && (class != H5T_INTEGER) ) {
		ERROR("Not a floating point or integer value.\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
	}

	/* Get the dimensionality.  We have to cope with scalars expressed as
	 * arrays with all dimensions 1, as well as zero-d arrays. */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > 64 ) {
		ERROR("Too many dimensions for hdfile_get_value\n");
		H5Tclose(type);
		H5Dclose(dh);
		return 1;
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
		if ( ( i==0 ) && (ev != NULL) && (ev->dim_length > 0)
		  && (size[i] > ev->dim_entries[0]) )
		{
			dim_flag = 1;
		} else {
			H5Tclose(type);
			H5Dclose(dh);
			return 1;
		}
	}

	if ( dim_flag == 0 ) {

		if ( H5Dread(dh, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, val) < 0 ) {
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

		r = H5Dread(dh, memtype, ms, sh, H5P_DEFAULT, val);
		if ( r < 0 )  {
			ERROR("Couldn't read value.\n");
			H5Tclose(type);
			H5Dclose(dh);
			return 1;
		}

	}

	free(subst_name);

	return 0;
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


static int make_dataspaces(hid_t dh, struct event *ev, hid_t *memspace,
                           hid_t *filespace)
{
	hsize_t *f_offset, *f_count;
	hsize_t *m_offset, *m_count;
	hid_t sh, mh;
	int ndims;
	int i;

	if ( ev == NULL ) {
		ERROR("Can't make dataspaces: no event ID\n");
		return 1;
	}

	/* Check that there are at least as many dim entries as dimensions */
	sh = H5Dget_space(dh);
	ndims = H5Sget_simple_extent_ndims(sh);
	if ( ndims > ev->dim_length ) {
		return 1;
	}

	/* Now set up arrays of offsets and counts in files and memory */
	f_offset = malloc(sizeof(hsize_t)*ndims);
	f_count = malloc(sizeof(hsize_t)*ndims);
	m_offset = malloc(sizeof(hsize_t)*ndims);
	m_count = malloc(sizeof(hsize_t)*ndims);
	if ( (f_offset == NULL) || (f_count == NULL)
	  || (m_offset == NULL) || (m_count == NULL) ) return 1;

	for ( i=0; i<ev->dim_length; i++ ) {
		f_offset[i] = ev->dim_entries[i];
		f_count[i] = 1;
		m_offset[i] = 0;
		m_count[i] = 1;
	}

	if ( H5Sselect_hyperslab(sh, H5S_SELECT_SET,
	                         f_offset, NULL, f_count, NULL) ) return 1;

	free(f_offset);
	free(f_count);

	mh = H5Screate_simple(ndims, m_count, NULL);
	if ( H5Sselect_hyperslab(mh, H5S_SELECT_SET,
	                         m_offset, NULL, m_count, NULL) ) return 1;
	free(m_offset);
	free(m_count);

	*memspace = mh;
	*filespace = sh;

	return 0;
}


static char *read_vlen_string(hid_t dh, struct event *ev)
{
	hid_t memspace, filespace;
	herr_t r;
	char *tmp;
	hid_t type;

	if ( make_dataspaces(dh, ev, &memspace, &filespace) ) {
		return strdup("[couldn't make dataspaces - variable len]");
	}

	type = H5Dget_type(dh);
	r = H5Dread(dh, type, memspace, filespace, H5P_DEFAULT, &tmp);
	H5Tclose(type);
	if ( r < 0 ) {
		return strdup("[couldn't read vlen string]");
	}

	H5Sclose(memspace);
	H5Sclose(filespace);

	/* Variable strings are 0-terminated */
	chomp(tmp);
	return tmp;
}


static char *read_fixed_string(hid_t dh, struct event *ev)
{
	hid_t memspace, filespace;
	herr_t r;
	hid_t sh, type;
	size_t size;
	char *tmp;

	type = H5Dget_type(dh);
	size = H5Tget_size(type);
	tmp = malloc(size+1);
	if ( tmp == NULL ) {
		H5Tclose(type);
		return strdup("[couldn't allocate string]");
	}

	if ( ev == NULL ) {

		/* Try a simple fixed-length string */
		sh = H5Dget_space(dh);
		if ( H5Sget_simple_extent_ndims(sh) ) {
			H5Tclose(type);
			return strdup("[non-scalar string]");
		}

		sh = H5Screate(H5S_SCALAR);
		r = H5Dread(dh, type, sh, H5S_ALL, H5P_DEFAULT, tmp);
		H5Sclose(sh);
		if ( r < 0 ) {
			free(tmp);
			H5Tclose(type);
			return strdup("[couldn't read scalar string]");
		} else {
			H5Tclose(type);
			tmp[size] = '\0';
			chomp(tmp);
			return tmp;
		}
	}

	if ( make_dataspaces(dh, ev, &memspace, &filespace) ) {
		H5Tclose(type);
		return strdup("[couldn't make dataspaces - fixed len]");
	}

	r = H5Dread(dh, type, memspace, filespace, H5P_DEFAULT, tmp);
	if ( r < 0 ) {
		H5Tclose(type);
		return strdup("[couldn't read string]");
	}

	H5Tclose(type);
	H5Sclose(memspace);
	H5Sclose(filespace);

	tmp[size] = '\0';
	chomp(tmp);
	return tmp;
}


static char *read_general_string(hid_t dh, struct event *ev)
{
	htri_t v;
	hid_t type;

	type = H5Dget_type(dh);
	v = H5Tis_variable_str(type);
	H5Tclose(type);

	if ( v < 0 ) {
		return strdup("[unrecognised string type]");

	} else if ( v > 0 ) {
		/* Variable length string */
		return read_vlen_string(dh, ev);

	} else {
		/* Fixed-length string */
		return read_fixed_string(dh, ev);

	}
}


char *hdfile_get_string_value(struct hdfile *f, const char *name,
                              struct event *ev)
{
	hid_t dh;
	hid_t type;
	hid_t class;
	int buf_i;
	double buf_f;
	char *tmp = NULL;
	char *subst_name = NULL;

	if ( (ev != NULL) && (ev->path_length != 0) ) {
		subst_name = retrieve_full_path(ev, name);
	} else {
		subst_name = strdup(name);
	}

	dh = H5Dopen2(f->fh, subst_name, H5P_DEFAULT);
	if ( dh < 0 ) {
		free(subst_name);
		return strdup("[couldn't read string]");
	}

	type = H5Dget_type(dh);
	class = H5Tget_class(type);
	H5Tclose(type);

	if ( class == H5T_STRING ) {

		free(subst_name);
		tmp = read_general_string(dh, ev);
		H5Dclose(dh);
		return tmp;

	} else {

		int r;

		H5Dclose(dh);

		switch ( class ) {

			case H5T_FLOAT :
			r = hdfile_get_value(f, subst_name, ev, &buf_f,
			                     H5T_NATIVE_DOUBLE);
			free(subst_name);
			if ( r == 0 ) {
				tmp = malloc(256);
				if ( tmp == NULL ) {
					ERROR("Failed to allocate float\n");
					return NULL;
				}
				snprintf(tmp, 255, "%f", buf_f);
				return tmp;
			} else {
				return NULL;
			}
			break;

			case H5T_INTEGER :
			r = hdfile_get_value(f, subst_name, ev, &buf_i,
			                     H5T_NATIVE_INT);
			free(subst_name);
			if ( r == 0 ) {
				tmp = malloc(256);
				if ( tmp == NULL ) {
					ERROR("Failed to allocate int buf!\n");
					return NULL;
				}
				snprintf(tmp, 255, "%d", buf_i);
				return tmp;

			} else {
				return NULL;
			}
			break;

			default :
			ERROR("Unrecognised type: %s\n", subst_name);
			free(subst_name);
			return NULL;
		}

	}
}


int check_path_existence(hid_t fh, const char *path)
{

	char buffer[256];
	char buffer_full_path[2048];
	herr_t herrt;
	struct H5O_info_t ob_info;
	char *path_copy = strdup(path);
	char *start = path_copy;
	char *sep = NULL;

	buffer[0] = '\0';
	buffer_full_path[0] = '\0';

	if ( strcmp(path_copy, "/" ) == 0 ) {
		return 1;
	}

	do {

		int check;

		sep = strstr(start, "/");
		if ( sep != NULL && strlen(sep) == 1 ) {
			ERROR("Error: Data path ends with a / symbol\n");
			free(path_copy);
			return 1;
		}

		if ( sep != NULL ) {

			if ( sep == start ) {
				start = sep+1;
				strcat(buffer_full_path, "/");
				continue;
			}

			strncpy(buffer, start, sep-start);
			buffer[sep-start] = '\0';
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

			strcpy(buffer, start);
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
