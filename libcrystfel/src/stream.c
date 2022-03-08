/*
 * stream.c
 *
 * Stream tools
 *
 * Copyright © 2013-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
 *   2014-2016 Valerio Mariani
 *   2011      Richard Kirian
 *   2011      Andrew Aquila
 *   2014      Takanori Nakane <nakane.t@gmail.com>
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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "image.h"
#include "stream.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "datatemplate.h"
#include "detgeom.h"
#include "libcrystfel-version.h"


/** \file stream.h */

#define LATEST_MAJOR_VERSION (2)
#define LATEST_MINOR_VERSION (3)

struct _stream
{
	FILE *fh;

	int major_version;
	int minor_version;
	char *audit_info;
	char *geometry_file;

	/* The DataTemplate provided to us for writing things.
	 * We don't own this. */
	const DataTemplate *dtempl_write;

	/* The DataTemplate we got from the stream itself, used for
	 * reading things.  We own this. */
	DataTemplate *dtempl_read;

	long long int ln;

	int old_indexers;  /* True if the stream reader encountered a deprecated
	                    * indexing method */

	long *chunk_offsets;
	int n_chunks;
};


int stream_has_old_indexers(Stream *st)
{
	return st->old_indexers;
}


static ImageFeatureList *read_peaks(Stream *st,
                                    struct image *image)
{
	char *rval = NULL;
	int first = 1;
	ImageFeatureList *features;

	features = image_feature_list_new();

	do {

		char line[1024];
		float x, y, d, intensity;
		int r, exp_n;
		char panel_name[1024];

		rval = fgets(line, 1023, st->fh);
		st->ln++;
		if ( rval == NULL ) {
			image_feature_list_free(features);
			return NULL;
		}
		chomp(line);

		if ( strcmp(line, STREAM_PEAK_LIST_END_MARKER) == 0 ) {
			return features;
		}

		if ( first ) {
			first = 0;
			continue;
		}

		if ( AT_LEAST_VERSION(st, 2, 3) ) {
			r = sscanf(line, "%f %f %f %f %64s",
			           &x, &y, &d, &intensity, panel_name);
			exp_n = 5;
		} else {
			r = sscanf(line, "%f %f %f %f",
			           &x, &y, &d, &intensity);
			exp_n = 4;
		}

		if ( r != exp_n ) {
			ERROR("Failed to parse peak list line.\n");
			ERROR("The failed line was: '%s'\n", line);
			image_feature_list_free(features);
			return NULL;
		}

		if ( (panel_name[0] != '\0') && (st->dtempl_read != NULL) ) {

			int pn;

			if ( data_template_panel_name_to_number(st->dtempl_read,
			                                        panel_name,
			                                        &pn) )
			{
				ERROR("No such panel '%s'\n",
				      panel_name);
			} else {

				data_template_file_to_panel_coords(st->dtempl_read,
				                                   &x, &y, pn);

				image_add_feature(features, x, y,
				                  pn, image, intensity,
				                  NULL);

			}

		} else {

			/* Either it's an old format stream (in which
			 * case the data is probably "slabby", so no
			 * coordinate conversion is needed), or
			 * the caller isn't interested in panel
			 * locations */
			image_add_feature(features, x, y, 0,
			                  image, intensity, NULL);

		}

	} while ( rval != NULL );

	return features;
}


static int write_peaks(const struct image *image,
                       const DataTemplate *dtempl, FILE *ofh)
{
	int i;

	fprintf(ofh, STREAM_PEAK_LIST_START_MARKER"\n");
	fprintf(ofh, "  fs/px   ss/px (1/d)/nm^-1   Intensity  Panel\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double r[3];
		double q;
		float write_fs, write_ss;
		struct detgeom_panel *p;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		p = &image->detgeom->panels[f->pn];
		detgeom_transform_coords(p, f->fs, f->ss,
		                         image->lambda, 0.0, 0.0, r);
		q = modulus(r[0], r[1], r[2]);

		write_fs = f->fs;
		write_ss = f->ss;
		data_template_panel_to_file_coords(dtempl, f->pn,
		                                   &write_fs, &write_ss);

		fprintf(ofh, "%7.2f %7.2f %10.2f  %10.2f   %s\n",
		        write_fs, write_ss, q/1.0e9, f->intensity,
		        data_template_panel_number_to_name(dtempl, f->pn));

	}

	fprintf(ofh, STREAM_PEAK_LIST_END_MARKER"\n");
	return 0;
}


static RefList *read_stream_reflections_2_3(Stream *st)
{
	char *rval = NULL;
	int first = 1;
	RefList *out;

	out = reflist_new();
	if ( out == NULL ) {
		ERROR("Failed to allocate reflection list\n");
		return NULL;
	}

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, sigma, fs, ss, pk, bg;
		char pname[64];
		int r;

		rval = fgets(line, 1023, st->fh);
		st->ln++;
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, STREAM_REFLECTION_END_MARKER) == 0 ) return out;

		r = sscanf(line, "%i %i %i %f %f %f %f %f %f %63s",
		           &h, &k, &l, &intensity, &sigma, &pk, &bg,
			   &fs, &ss, pname);

		if ( (r != 10) && (!first) ) {
			reflist_free(out);
			return NULL;
		}

		first = 0;

		if ( r == 10 ) {

			Reflection *refl;
			refl = add_refl(out, h, k, l);
			if ( refl == NULL ) {
				ERROR("Failed to add reflection\n");
				return NULL;
			}
			set_intensity(refl, intensity);
			if ( st->dtempl_read != NULL ) {
				int pn;

				if ( data_template_panel_name_to_number(st->dtempl_read,
			                                                pname, &pn) )
				{
					ERROR("No such panel '%s'\n", pname);
				} else {
					if ( data_template_file_to_panel_coords(st->dtempl_read, &fs, &ss, pn) ) {
						ERROR("Failed to convert\n");
					} else {
						set_detector_pos(refl, fs, ss);
						set_panel_number(refl, pn);
					}
				}
			}
			set_esd_intensity(refl, sigma);
			set_peak(refl, pk);
			set_mean_bg(refl, bg);
			set_redundancy(refl, 1);
			set_symmetric_indices(refl, h, k, l);
		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding STREAM_PEAK_LIST_END_MARKER */
	return NULL;
}


static int write_stream_reflections(FILE *fh, RefList *list,
                                    const DataTemplate *dtempl)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "   h    k    l          I   sigma(I)       "
	            "peak background  fs/px  ss/px panel\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		signed int h, k, l;
		double intensity, esd_i, pk, bg;
		double dfs, dss;
		float fs, ss;
		int pn;

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &dfs, &dss);
		fs = dfs;  ss = dss;
		pn = get_panel_number(refl);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		pk = get_peak(refl);
		bg = get_mean_bg(refl);

		/* Reflections with redundancy = 0 are not written */
		if ( get_redundancy(refl) == 0 ) continue;

		data_template_panel_to_file_coords(dtempl, pn,
		                                   &fs, &ss);

		fprintf(fh, "%4i %4i %4i %10.2f %10.2f %10.2f %10.2f "
		        "%6.1f %6.1f %s\n",
		        h, k, l, intensity, esd_i, pk, bg,
		        fs, ss, data_template_panel_number_to_name(dtempl, pn));

	}
	return 0;
}


static int num_integrated_reflections(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	int n = 0;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		if ( get_redundancy(refl) > 0 ) n++;
	}

	return n;
}


static int write_crystal(Stream *st, Crystal *cr,
                         int include_reflections)
{
	UnitCell *cell;
	RefList *reflist;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;
	double rad;
	double det_shift_x, det_shift_y;
	int ret = 0;

	fprintf(st->fh, STREAM_CRYSTAL_START_MARKER"\n");

	cell = crystal_get_cell(cr);
	assert(cell != NULL);

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	fprintf(st->fh, "Cell parameters %7.5f %7.5f %7.5f nm,"
			" %7.5f %7.5f %7.5f deg\n",
			a*1.0e9, b*1.0e9, c*1.0e9,
			rad2deg(al), rad2deg(be), rad2deg(ga));

	cell_get_reciprocal(cell, &asx, &asy, &asz,
				  &bsx, &bsy, &bsz,
				  &csx, &csy, &csz);
	fprintf(st->fh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	        asx/1e9, asy/1e9, asz/1e9);
	fprintf(st->fh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	        bsx/1e9, bsy/1e9, bsz/1e9);
	fprintf(st->fh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
		csx/1e9, csy/1e9, csz/1e9);

	fprintf(st->fh, "lattice_type = %s\n",
		str_lattice(cell_get_lattice_type(cell)));
	fprintf(st->fh, "centering = %c\n", cell_get_centering(cell));
	fprintf(st->fh, "unique_axis = %c\n", cell_get_unique_axis(cell));

	rad = crystal_get_profile_radius(cr);
	fprintf(st->fh, "profile_radius = %.5f nm^-1\n", rad/1e9);

	if ( crystal_get_notes(cr) != NULL ) {
		fprintf(st->fh, "%s\n", crystal_get_notes(cr));
	}

	crystal_get_det_shift(cr, &det_shift_x, &det_shift_y);

	fprintf(st->fh, "predict_refine/det_shift x = %.3f y = %.3f mm\n",
	        det_shift_x*1e3, det_shift_y*1e3);

	reflist = crystal_get_reflections(cr);
	if ( reflist != NULL ) {

		fprintf(st->fh, "diffraction_resolution_limit"
				" = %.2f nm^-1 or %.2f A\n",
				crystal_get_resolution_limit(cr)/1e9,
				1e10 / crystal_get_resolution_limit(cr));

		fprintf(st->fh, "num_reflections = %i\n",
		                num_integrated_reflections(reflist));
		fprintf(st->fh, "num_saturated_reflections = %lli\n",
		                crystal_get_num_saturated_reflections(cr));
		fprintf(st->fh, "num_implausible_reflections = %lli\n",
		                crystal_get_num_implausible_reflections(cr));

	}

	if ( include_reflections ) {

		if ( reflist != NULL ) {

			fprintf(st->fh, STREAM_REFLECTION_START_MARKER"\n");
			ret = write_stream_reflections(st->fh, reflist,
			                               st->dtempl_write);
			fprintf(st->fh, STREAM_REFLECTION_END_MARKER"\n");

		} else {

			fprintf(st->fh, "No integrated reflections.\n");

		}
	}

	fprintf(st->fh, STREAM_CRYSTAL_END_MARKER"\n");

	return ret;
}


/**
 * \param st A \ref Stream
 * \param i An \ref image structure
 * \param srf A \ref StreamFlags enum saying what to write
 *
 * Writes a new chunk to \p st.
 *
 * \returns non-zero on error.
 */
int stream_write_chunk(Stream *st, const struct image *i,
                       StreamFlags srf)
{
	int j;
	char *indexer;
	int ret = 0;

	fprintf(st->fh, STREAM_CHUNK_START_MARKER"\n");

	fprintf(st->fh, "Image filename: %s\n", i->filename);
	fprintf(st->fh, "Event: %s\n", i->ev);
	fprintf(st->fh, "Image serial number: %i\n", i->serial);

	fprintf(st->fh, "hit = %i\n", i->hit);
	indexer = indexer_str(i->indexed_by);
	fprintf(st->fh, "indexed_by = %s\n", indexer);
	free(indexer);
	if ( i->indexed_by != INDEXING_NONE ) {
		fprintf(st->fh, "n_indexing_tries = %i\n", i->n_indexing_tries);
	}

	fprintf(st->fh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(i->lambda)));

	fprintf(st->fh, "beam_divergence = %.2e rad\n", i->div);
	fprintf(st->fh, "beam_bandwidth = %.2e (fraction)\n", i->bw);

	for ( j=0; j<i->n_cached_headers; j++ ) {
		struct header_cache_entry *ce = i->header_cache[j];
		switch ( ce->type ) {

			case HEADER_FLOAT:
			fprintf(st->fh, "header/float/%s = %f\n",
			        ce->header_name, ce->val_float);
			break;

			case HEADER_INT:
			fprintf(st->fh, "header/int/%s = %i\n",
			        ce->header_name, ce->val_int);
			break;

			case HEADER_STR:
			fprintf(st->fh, "header/str/%s = %s\n",
			        ce->header_name, ce->val_str);
			break;

			default:
			ERROR("Unrecognised header cache type %i\n", ce->type);
			break;

		}
	}

	if ( i->detgeom != NULL ) {

		int pn;
		double tclen = 0.0;

		for ( pn=0; pn<i->detgeom->n_panels; pn++ ) {
			tclen += i->detgeom->panels[pn].cnz
				* i->detgeom->panels[pn].pixel_pitch;
		}
		fprintf(st->fh, "average_camera_length = %f m\n",
		        tclen / i->detgeom->n_panels);

	}

	fprintf(st->fh, "num_peaks = %i\n", image_feature_count(i->features));
	fprintf(st->fh, "peak_resolution = %f nm^-1 or %f A\n",
	        i->peak_resolution/1e9, 1e10/i->peak_resolution);
	if ( srf & STREAM_PEAKS ) {
		ret = write_peaks(i, st->dtempl_write, st->fh);
	}

	for ( j=0; j<i->n_crystals; j++ ) {
		if ( crystal_get_user_flag(i->crystals[j]) ) {
			continue;
		}
		ret = write_crystal(st, i->crystals[j],
		                    srf & STREAM_REFLECTIONS);
	}

	fprintf(st->fh, STREAM_CHUNK_END_MARKER"\n");

	fflush(st->fh);

	return ret;
}


static int find_start_of_chunk(Stream *st)
{
	char *rval = NULL;
	char line[1024];

	do {

		rval = fgets(line, 1023, st->fh);
		st->ln++;

		/* Trouble? */
		if ( rval == NULL ) return 1;

		chomp(line);

	} while ( strcmp(line, STREAM_CHUNK_START_MARKER) != 0 );

	return 0;
}


static void read_crystal(Stream *st, struct image *image,
                         StreamFlags srf)
{
	char line[1024];
	char *rval = NULL;
	struct rvec as, bs, cs;
	int have_as = 0;
	int have_bs = 0;
	int have_cs = 0;
	int have_latt = 0;
	int have_cen = 0;
	int have_ua = 0;
	char centering = 'P';
	char unique_axis = '*';
	LatticeType lattice_type = L_TRICLINIC;
	Crystal *cr;
	int n;
	Crystal **crystals_new;
	double shift_x, shift_y;

	as.u = 0.0;  as.v = 0.0;  as.w = 0.0;
	bs.u = 0.0;  bs.v = 0.0;  bs.w = 0.0;
	cs.u = 0.0;  cs.v = 0.0;  cs.w = 0.0;

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal!\n");
		return;
	}

	do {

		float u, v, w, lim, rad;
		char c;

		rval = fgets(line, 1023, st->fh);
		st->ln++;

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);
		if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) == 3 )
		{
			as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
			have_as = 1;
		}

		if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) == 3 )
		{
			bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
			have_bs = 1;
		}

		if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) == 3 )
		{
			cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
			have_cs = 1;
		}

		if ( sscanf(line, "centering = %c", &c) == 1 )
		{
			if ( !have_cen ) {
				centering = c;
				have_cen = 1;
			} else {
				ERROR("Duplicate centering (line %lli) - "
				      "stream may be corrupted!\n", st->ln);
			}
		}

		if ( sscanf(line, "unique_axis = %c", &c) == 1 )
		{
			if ( !have_ua ) {
				unique_axis = c;
				have_ua = 1;
			} else {
				ERROR("Duplicate unique axis (line %lli) - "
				      "stream may be corrupted!\n", st->ln);
			}
		}

		if ( strncmp(line, "lattice_type = ", 15) == 0 )
		{
			if ( !have_latt ) {
				lattice_type = lattice_from_str(line+15);
				have_latt = 1;
			} else {
				ERROR("Duplicate lattice type (line %lli) - "
				      "stream may be corrupted!\n", st->ln);
			}
		}

		if ( strncmp(line, "num_saturated_reflections = ", 28) == 0 ) {
			int nsat = atoi(line+28);
			crystal_set_num_saturated_reflections(cr, nsat);
		}

		if ( sscanf(line, "diffraction_resolution_limit = %f nm^-1",
		            &lim) == 1 ) {
			crystal_set_resolution_limit(cr, lim*1e9);
		}

		if ( sscanf(line, "profile_radius = %e nm^-1", &rad) == 1 ) {
			crystal_set_profile_radius(cr, rad*1e9);
		}

		if ( sscanf(line, "predict_refine/det_shift x = %lf "
		                  "y = %lf mm\n", &shift_x, &shift_y ) == 2 ) {
			crystal_set_det_shift(cr, shift_x*1e-3, shift_y*1e-3);
		}


		if ( (strcmp(line, STREAM_REFLECTION_START_MARKER) == 0)
		  && (srf & STREAM_REFLECTIONS) )
		{

			RefList *reflist;
			reflist = read_stream_reflections_2_3(st);
			if ( reflist == NULL ) {
				ERROR("Failed while reading reflections\n");
				ERROR("Filename = %s\n", image->filename);
				ERROR("Event = %s\n", image->ev);
				break;
			}

			crystal_set_reflections(cr, reflist);

		}

		if ( strcmp(line, STREAM_CRYSTAL_END_MARKER) == 0 ) break;

	} while ( 1 );

	if ( have_as && have_bs && have_cs ) {

		UnitCell *cell;

		cell = crystal_get_cell(cr);

		if ( cell != NULL ) {
			ERROR("Duplicate cell found in stream!\n");
			ERROR("I'll use the most recent one.\n");
			cell_free(cell);
		}

		cell = cell_new_from_reciprocal_axes(as, bs, cs);
		if ( cell == NULL ) {
			ERROR("Failed to allocate cell\n");
			return;
		}

		if ( have_cen && have_ua && have_latt ) {
			cell_set_centering(cell, centering);
			cell_set_unique_axis(cell, unique_axis);
			cell_set_lattice_type(cell, lattice_type);
		} /* else keep default triclinic P */

		crystal_set_cell(cr, cell);

	}

	/* Unused at the moment */
	crystal_set_mosaicity(cr, 0.0);

	/* Add crystal to the list for this image */
	n = image->n_crystals+1;
	crystals_new = realloc(image->crystals, n*sizeof(Crystal *));

	if ( crystals_new == NULL ) {
		ERROR("Failed to expand crystal list!\n");
	} else {
		image->crystals = crystals_new;
		image->crystals[image->n_crystals++] = cr;
	}

}


static void parse_header(const char *line_in, struct image *image,
                         HeaderCacheType type)
{
	char *line;
	char *pos;

	line = strdup(line_in);
	chomp(line);

	pos = strchr(line, ' ');
	if ( pos == NULL ) {
		ERROR("Invalid header line '%s' (no space)\n", line);
		return;
	}
	pos[0] = '\0';

	if ( strlen(line_in) < strlen(line) + 3 ) {
		ERROR("Invalid header line '%s' (too short)\n", line);
		return;
	}

	if ( (pos[1] != '=') || (pos[2] != ' ') ) {
		ERROR("Invalid header line '%s' (wrong separator)\n", line);
		return;
	}

	switch ( type ) {

		double vf;
		int vi;

		case HEADER_FLOAT:
		if ( convert_float(pos+3, &vf) != 0 ) {
			ERROR("Invalid header line '%s' (invalid value)\n", line);
			return;
		}
		image_cache_header_float(image, line, vf);
		break;

		case HEADER_INT:
		if ( convert_int(pos+3, &vi) != 0 ) {
			ERROR("Invalid header line '%s' (invalid value)\n", line);
			return;
		}
		image_cache_header_int(image, line, vi);
		break;

		case HEADER_STR:
		image_cache_header_str(image, line, pos+3);
		break;

		default:
		ERROR("Unrecognised header cache type %i (from stream)\n", type);
		break;
	}
}


/**
 * Read the next chunk from a stream and return an image structure
 */
struct image *stream_read_chunk(Stream *st, StreamFlags srf)
{
	char line[1024];
	char *rval = NULL;
	int have_filename = 0;
	int have_ev = 0;
	struct image *image;

	if ( find_start_of_chunk(st) ) return NULL;

	image = image_new();
	if ( image == NULL ) return NULL;

	image->data_source_type = DATA_SOURCE_TYPE_NONE;

	do {
		int ser;
		float div, bw;

		rval = fgets(line, 1023, st->fh);
		st->ln++;

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);

		if ( strncmp(line, "Image filename: ", 16) == 0 ) {
			image->filename = strdup(line+16);
			have_filename = 1;
		}

		if ( strncmp(line, "Event: ", 7) == 0 ) {
			image->ev = strdup(line+7);
		}

		if ( strncmp(line, "hdf5/", 5) == 0 ) {
			parse_header(line+4, image, HEADER_FLOAT);
		}

		if ( strncmp(line, "header/int/", 11) == 0 ) {
			parse_header(line+11, image, HEADER_INT);
		}

		if ( strncmp(line, "header/float/", 13) == 0 ) {
			parse_header(line+13, image, HEADER_FLOAT);
		}

		if ( strncmp(line, "header/str/", 13) == 0 ) {
			parse_header(line+11, image, HEADER_STR);
		}

		if ( strncmp(line, "indexed_by = ", 13) == 0 ) {
			int err = 0;
			image->indexed_by = get_indm_from_string_2(line+13, &err);
			if ( image->indexed_by == INDEXING_ERROR ) {
				ERROR("Failed to read indexer list\n");
			}
			if ( err ) {
				st->old_indexers = 1;
			}
		}

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			image->lambda = ph_en_to_lambda(eV_to_J(atof(line+19)));
			have_ev = 1;
		}

		if ( sscanf(line, "beam_divergence = %e rad", &div) == 1 ) {
			image->div = div;
		}

		if ( sscanf(line, "beam_bandwidth = %f", &bw) == 1 ) {
			image->bw = bw;
		}

		if ( sscanf(line, "Image serial number: %i", &ser) == 1 ) {
			image->serial = ser;
		}


		if ( (srf & STREAM_PEAKS)
		    && strcmp(line, STREAM_PEAK_LIST_START_MARKER) == 0 ) {

			ImageFeatureList *peaks;
			peaks = read_peaks(st, image);

			if ( peaks == NULL ) {
				ERROR("Failed while reading peaks\n");
				image_free(image);
				return NULL;
			}

			image->features = peaks;


		}

		if ( strcmp(line, STREAM_CRYSTAL_START_MARKER) == 0 ) {
			read_crystal(st, image, srf);
		}

		/* A chunk must have at least a filename and a wavelength,
		 * otherwise it's incomplete */
		if ( strcmp(line, STREAM_CHUNK_END_MARKER) == 0 ) {
			if ( have_filename && have_ev ) {
				/* Success */
				if ( srf & STREAM_DATA_DETGEOM ) {
					if ( create_detgeom(image, st->dtempl_read) ) {
						image_free(image);
						return NULL;
					}
					image_set_zero_data(image, st->dtempl_read);
					image_set_zero_mask(image, st->dtempl_read);
				}
				image->spectrum = spectrum_generate_gaussian(image->lambda,
				                                             image->bw);
				return image;
			}
			ERROR("Incomplete chunk found in input file.\n");
			image_free(image);
			return NULL;
		}

	} while ( 1 );

	if ( !feof(st->fh) ) {
		ERROR("Error reading stream.\n");
	}

	image_free(image);
	return NULL;  /* Either error or EOF, don't care because we will complain
	               * on the terminal if it was an error. */
}


char *stream_audit_info(Stream *st)
{
	if ( st->audit_info == NULL ) return NULL;
	return strdup(st->audit_info);
}


char *stream_geometry_file(Stream *st)
{
	return st->geometry_file;
}


static void read_geometry_file(Stream *st)
{
	int done = 0;
	size_t len = 0;
	const size_t max_geom_len = 64*1024;
	char *geom;

	geom = malloc(max_geom_len);
	if ( geom == NULL ) {
		ERROR("Failed to allocate memory for geometry file\n");
		return;
	}
	geom[0] = '\0';

	do {

		char line[1024];
		char *rval;

		rval = fgets(line, 1023, st->fh);
		st->ln++;
		if ( rval == NULL ) {
			ERROR("Failed to read stream geometry file.\n");
			stream_close(st);
			free(geom);
			return;
		}

		if ( strcmp(line, STREAM_GEOM_END_MARKER"\n") == 0 ) {
			done = 1;
			continue;
		}

		len += strlen(line);
		if ( len > max_geom_len-1 ) {
			ERROR("Stream's geometry file is too long (%li > %i).\n",
			      (long)len, (int)max_geom_len);
			free(geom);
			return;
		} else {
			strcat(geom, line);
		}

	} while  ( !done );

	st->geometry_file = geom;
	st->dtempl_read = data_template_new_from_string(geom);
}


static void read_headers(Stream *st)
{
	int done = 0;
	size_t len = 0;

	st->audit_info = malloc(4096);
	if ( st->audit_info == NULL ) {
		ERROR("Failed to allocate memory for audit information\n");
		return;
	}
	st->audit_info[0] = '\0';

	/* Read lines from stream until one of them starts with "-----",
	 * then rewind to the start of that line */
	do {

		char line[1024];
		char *rval;

		rval = fgets(line, 1023, st->fh);
		st->ln++;
		if ( rval == NULL ) {
			ERROR("Failed to read stream audit info.\n");
			stream_close(st);
			return;
		}

		if ( strcmp(line, STREAM_GEOM_START_MARKER"\n") == 0 ) {
			read_geometry_file(st);
			done = 1;
		} else {
			len += strlen(line);
			if ( len > 4090 ) {
				ERROR("Too much audit information.\n");
				return;
			} else {
				strcat(st->audit_info, line);
			}
		}

	} while  ( !done );
}


Stream *stream_open_for_read(const char *filename)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;
	st->old_indexers = 0;
	st->audit_info = NULL;
	st->geometry_file = NULL;
	st->n_chunks = 0;
	st->chunk_offsets = NULL;
	st->dtempl_read = NULL;
	st->dtempl_write = NULL;

	if ( strcmp(filename, "-") == 0 ) {
		st->fh = stdin;
	} else {
		st->fh = fopen(filename, "r");
	}

	if ( st->fh == NULL ) {
		free(st);
		return NULL;
	}

	char line[1024];
	char *rval;

	rval = fgets(line, 1023, st->fh);
	if ( rval == NULL ) {
		ERROR("Failed to read stream version.\n");
		stream_close(st);
		return NULL;
	}

	if ( strncmp(line, "CrystFEL stream format 2.3", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 3;
	} else {
		ERROR("Invalid stream, or stream format is not understood.\n");
		stream_close(st);
		return NULL;
	}

	st->ln = 1;

	read_headers(st);

	return st;
}


/**
 * \param fd File descriptor (e.g. from open()) to use for stream data.
 *
 * Creates a new \ref Stream from \p fd, so that stream data can be written to \p fd
 * using \ref stream_write_chunk.
 *
 * In contrast to \ref stream_open_for_write, this function does not write any of
 * the usual headers.  This function is mostly for use when multiple substreams
 * need to be multiplexed into a single master stream.  The master would be
 * opened using \ref stream_open_for_write, and the substreams using this function.
 *
 * \returns A \ref Stream, or NULL on failure.
 */
Stream *stream_open_fd_for_write(int fd, const DataTemplate *dtempl)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;
	st->old_indexers = 0;
	st->audit_info = NULL;
	st->geometry_file = NULL;
	st->n_chunks = 0;
	st->chunk_offsets = NULL;
	st->dtempl_read = NULL;
	st->dtempl_write = NULL;

	st->fh = fdopen(fd, "w");
	if ( st->fh == NULL ) {
		free(st);
		return NULL;
	}

	st->dtempl_write = dtempl;
	st->major_version = LATEST_MAJOR_VERSION;
	st->minor_version = LATEST_MINOR_VERSION;

	return st;
}


void stream_write_target_cell(Stream *st, UnitCell *cell)
{
	if ( cell == NULL ) return;
	fprintf(st->fh, STREAM_CELL_START_MARKER"\n");
	write_cell(cell, st->fh);
	fprintf(st->fh, "; Please note: this is the target unit cell.\n");
	fprintf(st->fh, "; The actual unit cells produced by indexing "
	                "depend on many other factors.\n");
	fprintf(st->fh, STREAM_CELL_END_MARKER"\n");
	fflush(st->fh);
}


/**
 * \param filename Filename of new stream
 * \param dtempl A DataTemplate
 *
 * Creates a new stream with name \p filename.  If \p filename already
 * exists, it will be overwritten.
 *
 * The CrystFEL version number will be written, but you should call
 * stream_write_geometry_file, stream_write_target_cell,
 * stream_write_commandline_args and stream_write_Indexing_methods to write
 * extended audit information.
 *
 * \returns A \ref Stream, or NULL on failure.
 */
Stream *stream_open_for_write(const char *filename,
                              const DataTemplate *dtempl)

{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;
	st->old_indexers = 0;
	st->audit_info = NULL;
	st->geometry_file = NULL;
	st->n_chunks = 0;
	st->chunk_offsets = NULL;
	st->dtempl_write = dtempl;
	st->dtempl_read = NULL;

	st->fh = fopen(filename, "w");
	if ( st->fh == NULL ) {
		ERROR("Failed to open stream.\n");
		free(st);
		return NULL;
	}

	st->major_version = LATEST_MAJOR_VERSION;
	st->minor_version = LATEST_MINOR_VERSION;

	fprintf(st->fh, "CrystFEL stream format %i.%i\n",
	        st->major_version, st->minor_version);
	fprintf(st->fh, "Generated by CrystFEL %s\n",
	        libcrystfel_version_string());
	fflush(st->fh);

	return st;
}


int stream_get_fd(Stream *st)
{
	return fileno(st->fh);
}


/**
 * \param st A \ref Stream
 *
 * Closes the stream
 */
void stream_close(Stream *st)
{
	if ( st == NULL ) return;
	free(st->audit_info);
	free(st->geometry_file);
	data_template_free(st->dtempl_read);
	fclose(st->fh);
	free(st);
}


/**
 * \param st A \ref Stream
 * \param argc number of arguments
 * \param argv command-line arguments
 *
 * Writes the command line to \p st.  \p argc and \p argv should be
 * exactly as were given to main().  This should usually be called
 * immediately after \ref stream_open_for_write.
 */
void stream_write_commandline_args(Stream *st, int argc, char *argv[])
{
	int i;

	if ( argc == 0 ) return;

	for ( i=0; i<argc; i++ ) {
		if ( i > 0 ) fprintf(st->fh, " ");
		fprintf(st->fh, "%s", argv[i]);
	}
	fprintf(st->fh, "\n");
	fflush(st->fh);
}


void stream_write_indexing_methods(Stream *st, const char *indm_str)
{
	fprintf(st->fh, "Indexing methods selected: %s\n", indm_str);
	fflush(st->fh);
}


/**
 * \param st A \ref Stream
 * \param geom_filename geomtry file name
 *
 * Writes the content of the geometry file to \p st. This should usually be
 * called immediately after \ref write_command.
 */
void stream_write_geometry_file(Stream *st, const char *geom_filename)
{
	char line[2014];
	FILE *geom_fh;
	char *rval;
	int eol;

	if ( geom_filename == NULL ) return;

	geom_fh = fopen(geom_filename, "r");
	if ( geom_fh == NULL ) {
		ERROR("Failed to read detector geometry from "
		      "'%s'\n", geom_filename);
		return;
	}
	fprintf(st->fh, STREAM_GEOM_START_MARKER"\n");

	do {
		rval = fgets(line, 1023, geom_fh);
		if ( rval != NULL ) fputs(line, st->fh);
		eol = ( line[strlen(line)-1] == '\n' );
	} while ( rval != NULL );


	if ( !eol ) {
		fprintf(st->fh, "\n");
	}

	fclose(geom_fh);

	fprintf(st->fh, STREAM_GEOM_END_MARKER"\n");
	fflush(st->fh);
}


/**
 * \param st A \ref Stream
 *
 * Attempts to set the file pointer for \p st to the start of the stream, so that
 * later calls to \ref stream_read_chunk will repeat the sequence of chunks from the
 * start.
 *
 * Programs must not assume that this operation always succeeds!
 *
 * \returns Non-zero if the stream could not be rewound.
 */
int stream_rewind(Stream *st)
{
	st->ln = 0;
	return fseek(st->fh, 0, SEEK_SET);
}


struct _streamindex
{
	char **keys;
	long int *ptrs;
	int n_keys;
	int max_keys;
};


void stream_index_free(StreamIndex *index)
{
	if ( index == NULL ) return;
	free(index->keys);
	free(index->ptrs);
	free(index);
}


static char *make_key(const char *filename,
                      const char *ev)
{
	char *key;

	if ( ev == NULL ) ev = "//";

	key = malloc(strlen(filename)+strlen(ev)+2);
	if ( key == NULL ) return NULL;

	strcpy(key, filename);
	strcat(key, " ");
	strcat(key, ev);

	return key;
}


int stream_select_chunk(Stream *st,
                        StreamIndex *index,
                        const char *filename,
                        const char *ev)
{
	int i;
	char *key;

	if ( index == NULL ) return 1;

	key = make_key(filename, ev);
	for ( i=0; i<index->n_keys; i++ ) {
		if ( strcmp(index->keys[i], key) == 0 ) {
			if ( st != NULL ) {
				fseek(st->fh, index->ptrs[i], SEEK_SET);
			}
			return 0;
		}
	}
	return 1;
}


static void add_index_record(StreamIndex *index,
                             long int ptr,
                             const char *filename,
                             const char *ev)
{
	char *key;

	if ( index->n_keys == index->max_keys ) {

		int new_max_keys = index->max_keys + 256;
		char **new_keys;
		long int *new_ptrs;

		new_keys = realloc(index->keys,
		                   new_max_keys*sizeof(char *));
		if ( new_keys == NULL ) return;

		new_ptrs = realloc(index->ptrs,
		                   new_max_keys*sizeof(long int));
		if ( new_ptrs == NULL ) {
			free(new_keys);
			return;
		}

		index->keys = new_keys;
		index->ptrs = new_ptrs;
		index->max_keys = new_max_keys;

	}

	key = make_key(filename, ev);
	if ( key == NULL ) return;

	index->keys[index->n_keys] = key;
	index->ptrs[index->n_keys] = ptr;
	index->n_keys++;
}


StreamIndex *stream_make_index(const char *filename)
{
	FILE *fh;
	StreamIndex *index;
	long int last_start_pos = 0;
	char *last_filename = NULL;
	char *last_ev = NULL;
	int done = 0;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	index = malloc(sizeof(StreamIndex));
	if ( index == NULL ) {
		fclose(fh);
		return NULL;
	}

	index->keys = NULL;
	index->ptrs = NULL;
	index->n_keys = 0;
	index->max_keys = 0;

	STATUS("Scanning %s\n", filename);

	do {

		char *rval;
		char line[1024];
		long int pos;

		pos = ftell(fh);
		rval = fgets(line, 1024, fh);
		if ( rval == NULL ) break;
		chomp(line);

		if ( strcmp(line, STREAM_CHUNK_START_MARKER) == 0 ) {
			last_start_pos = pos;
			last_filename = NULL;
			last_ev = NULL;
		}

		if ( strncmp(line, "Image filename: ", 16) == 0 ) {
			last_filename = strdup(line+16);
		}

		if ( strncmp(line, "Event: ", 7) == 0 ) {
			last_ev = strdup(line+7);
		}

		if ( strcmp(line, STREAM_CHUNK_END_MARKER) == 0 ) {
			if ( (last_start_pos != 0)
			     && (last_filename != NULL) )
			{
				add_index_record(index,
				                 last_start_pos,
				                 last_filename,
				                 last_ev);
			}
			free(last_filename);
			free(last_ev);
			last_start_pos = 0;
			last_filename = NULL;
			last_ev = NULL;
		}

	} while ( !done );

	fclose(fh);
	return index;
}
