/*
 * stream.c
 *
 * Stream tools
 *
 * Copyright © 2013-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "version.h"
#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "image.h"
#include "stream.h"
#include "reflist.h"
#include "reflist-utils.h"

#define LATEST_MAJOR_VERSION (2)
#define LATEST_MINOR_VERSION (3)

#define AT_LEAST_VERSION(st, a, b) ((st->major_version>=(a)) \
                                    && (st->minor_version>=(b)))

struct _stream
{
	FILE *fh;

	int major_version;
	int minor_version;
};

static int read_peaks(FILE *fh, struct image *image)
{
	char *rval = NULL;
	int first = 1;

	image->features = image_feature_list_new();

	do {

		char line[1024];
		float x, y, d, intensity;
		int r;
		struct panel *p = NULL;
		float add_x, add_y;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, PEAK_LIST_END_MARKER) == 0 ) return 0;

		r = sscanf(line, "%f %f %f %f", &x, &y, &d, &intensity);
		if ( (r != 4) && (!first) ) {
			ERROR("Failed to parse peak list line.\n");
			ERROR("The failed line was: '%s'\n", line);
			return 1;
		}

		first = 0;
		if ( r == 4 ) {

			if ( image->det != NULL ) {

				p = find_orig_panel(image->det, x, y);
				if ( p == NULL ) {
					ERROR("Panel not found\n");
					return 1;
				}

				add_x = x-p->orig_min_fs+p->min_fs;
				add_y = y-p->orig_min_ss+p->min_ss;

				image_add_feature(image->features, add_x, add_y,
				                  image, intensity, NULL);

			} else {

				image_add_feature(image->features, x, y,
				image, intensity, NULL);
			}
		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return 1;
}


static int read_peaks_2_3(FILE *fh, struct image *image)
{
	char *rval = NULL;
	int first = 1;

	image->features = image_feature_list_new();

	do {

		char line[1024];
		char pn[32];
		float x, y, d, intensity;
		int r;
		struct panel *p = NULL;
		float add_x, add_y;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, PEAK_LIST_END_MARKER) == 0 ) return 0;

		r = sscanf(line, "%f %f %f %f %s", &x, &y, &d, &intensity, pn);
		if ( (r != 5) && (!first) ) {
			ERROR("Failed to parse peak list line.\n");
			ERROR("The failed line was: '%s'\n", line);
			return 1;
		}

		first = 0;

		if ( r == 5 ) {

			p = find_panel_by_name(image->det, pn);
			if ( p == NULL ) {
				ERROR("Panel not found: %s\n", pn);
				return 1;
			}

			add_x = x-p->orig_min_fs+p->min_fs;
			add_y = y-p->orig_min_ss+p->min_ss;

			image_add_feature(image->features, add_x, add_y,
			                  image, intensity, NULL);

		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return 1;
}


static int write_peaks(struct image *image, FILE *ofh)
{
	int i;

	fprintf(ofh, PEAK_LIST_START_MARKER"\n");
	fprintf(ofh, "  fs/px   ss/px  (1/d)/nm^-1   Intensity\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;
		double q;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		if ( image->det != NULL ) {

			struct panel *p;
			double write_fs, write_ss;

			p = find_orig_panel(image->det, f->fs, f->ss);
			if ( p == NULL ) {
				ERROR("Panel not found\n");
				return 1;
			}

			/* Convert coordinates to match arrangement of panels in
			 * HDF5 file */
			write_fs = f->fs - p->min_fs + p->orig_min_fs;
			write_ss = f->ss - p->min_ss + p->orig_min_ss;

			fprintf(ofh, "%7.2f %7.2f   %10.2f  %10.2f\n",
			        write_fs, write_ss, q/1.0e9, f->intensity);

		} else {

			fprintf(ofh, "%7.2f %7.2f   %10.2f  %10.2f\n",
			        f->fs, f->ss, q/1.0e9, f->intensity);

		}

	}

	fprintf(ofh, PEAK_LIST_END_MARKER"\n");
	return 0;
}


static int write_peaks_2_3(struct image *image, FILE *ofh)
{
	int i;

	fprintf(ofh, PEAK_LIST_START_MARKER"\n");
	fprintf(ofh, "  fs/px   ss/px (1/d)/nm^-1   Intensity  Panel\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;
		double q;
		struct panel *p;
		double write_fs, write_ss;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		p = find_panel(image->det, f->fs, f->ss);
		if ( p == NULL ) {
			ERROR("Panel not found\n");
			return 1;
		}

		/* Convert coordinates to match arrangement of panels in HDF5
		 * file */
		write_fs = f->fs - p->min_fs + p->orig_min_fs;
		write_ss = f->ss - p->min_ss + p->orig_min_ss;

		fprintf(ofh, "%7.2f %7.2f %10.2f  %10.2f   %s\n",
		        write_fs, write_ss, q/1.0e9, f->intensity, p->name);

	}

	fprintf(ofh, PEAK_LIST_END_MARKER"\n");
	return 0;
}


static RefList *read_stream_reflections_2_3(FILE *fh, struct detector *det)
{
	char *rval = NULL;
	int first = 1;
	RefList *out;

	out = reflist_new();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, sigma, fs, ss, pk, bg;
		char pn[32];
		int r;
		Reflection *refl;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) return out;

		r = sscanf(line, "%i %i %i %f %f %f %f %f %f %s",
		           &h, &k, &l, &intensity, &sigma, &pk, &bg,
			   &fs, &ss, pn);

		if ( (r != 10) && (!first) ) {
			reflist_free(out);
			return NULL;
		}

		first = 0;

		if ( r == 10 ) {

			struct panel *p;

			refl = add_refl(out, h, k, l);
			set_intensity(refl, intensity);
			if ( det != NULL ) {
				double write_fs, write_ss;
				p = find_panel_by_name(det,pn);
				write_fs = fs - p->orig_min_fs + p->min_fs;
				write_ss = ss - p->orig_min_ss + p->min_ss;
				set_detector_pos(refl, write_fs, write_ss);
			}
			set_esd_intensity(refl, sigma);
			set_peak(refl, pk);
			set_mean_bg(refl, bg);
			set_redundancy(refl, 1);
		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return NULL;
}


static RefList *read_stream_reflections_2_1(FILE *fh, struct detector *det)
{
	char *rval = NULL;
	int first = 1;
	RefList *out;

	out = reflist_new();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, sigma, fs, ss;
		char phs[1024];
		int cts;
		int r;
		Reflection *refl;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) return out;

		r = sscanf(line, "%i %i %i %f %s %f %i %f %f",
		                  &h, &k, &l, &intensity, phs, &sigma, &cts,
		                   &fs, &ss);
		if ( (r != 9) && (!first) ) {
			reflist_free(out);
			return NULL;
		}

		first = 0;
		if ( r == 9 ) {

			double ph;
			char *v;

			refl = add_refl(out, h, k, l);
			set_intensity(refl, intensity);

			if ( det != NULL ) {

				double write_fs, write_ss;
				struct panel *p;

				p = find_orig_panel(det, fs, ss);
				write_fs = fs - p->orig_min_fs + p->min_fs;
				write_ss = ss - p->orig_min_ss + p->min_ss;
				set_detector_pos(refl, write_fs, write_ss);

			} else {

				set_detector_pos(refl, fs, ss);

			}
			set_esd_intensity(refl, sigma);
			set_redundancy(refl, cts);

			ph = strtod(phs, &v);
			if ( v != phs ) set_phase(refl, deg2rad(ph));

		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return NULL;
}


static RefList *read_stream_reflections_2_2(FILE *fh, struct detector *det)
{
	char *rval = NULL;
	int first = 1;
	RefList *out;

	out = reflist_new();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, sigma, fs, ss, pk, bg;
		int r;
		Reflection *refl;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) return out;

		r = sscanf(line, "%i %i %i %f %f %f %f %f %f",
		           &h, &k, &l, &intensity, &sigma, &pk, &bg, &fs, &ss);
		if ( (r != 9) && (!first) ) {
			reflist_free(out);
			return NULL;
		}

		first = 0;
		if ( r == 9 ) {

			refl = add_refl(out, h, k, l);
			set_intensity(refl, intensity);

			if ( det != NULL ) {

				double write_fs, write_ss;
				struct panel *p;

				p = find_orig_panel(det, fs, ss);
				write_fs = fs - p->orig_min_fs + p->min_fs;
				write_ss = ss - p->orig_min_ss + p->min_ss;
				set_detector_pos(refl, write_fs, write_ss);

			} else {

				set_detector_pos(refl, fs, ss);

			}

			set_esd_intensity(refl, sigma);
			set_redundancy(refl, 1);
			set_peak(refl, pk);
			set_mean_bg(refl, bg);

		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding REFLECTION_END_MARKER */
	return NULL;
}


static int write_stream_reflections_2_1(FILE *fh, RefList *list,
                                        struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "  h   k   l          I    phase   sigma(I) "
			 " counts  fs/px  ss/px\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		signed int h, k, l;
		double intensity, esd_i, ph;
		int red;
		double fs, ss;
		char phs[16];
		int have_phase;

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &fs, &ss);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		red = get_redundancy(refl);
		ph = get_phase(refl, &have_phase);

		/* Reflections with redundancy = 0 are not written */
		if ( red == 0 ) continue;

		if ( have_phase ) {
			snprintf(phs, 16, "%8.2f", rad2deg(ph));
		} else {
			strncpy(phs, "       -", 15);
		}

		if ( image->det != NULL ) {

			struct panel *p;
			double write_fs, write_ss;

			p = find_orig_panel(image->det, fs, ss);
			if ( p == NULL ) {
				ERROR("Panel not found\n");
				return 1;
			}

			/* Convert coordinates to match arrangement of panels
			 * in HDF5 file */
			write_fs = fs - p->min_fs + p->orig_min_fs;
			write_ss = ss - p->min_ss + p->orig_min_ss;

			fprintf(fh, "%3i %3i %3i %10.2f %s %10.2f %7i "
			            "%6.1f %6.1f\n",
			             h, k, l, intensity, phs, esd_i, red,
				     write_fs, write_ss);

		} else {

			fprintf(fh, "%3i %3i %3i %10.2f %s %10.2f %7i "
			            "%6.1f %6.1f\n",
			            h, k, l, intensity, phs, esd_i, red,
				    fs, ss);

		}
	}
	return 0;
}


static int write_stream_reflections_2_2(FILE *fh, RefList *list,
                                        struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "   h    k    l          I   sigma(I)       "
	            "peak background  fs/px  ss/px\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		signed int h, k, l;
		double intensity, esd_i, bg, pk;
		double fs, ss;

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &fs, &ss);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		pk = get_peak(refl);
		bg = get_mean_bg(refl);

		/* Reflections with redundancy = 0 are not written */
		if ( get_redundancy(refl) == 0 ) continue;

		if ( image->det != NULL ) {

			struct panel *p;
			double write_fs, write_ss;

			p = find_orig_panel(image->det, fs, ss);
			if ( p == NULL ) {
				ERROR("Panel not found\n");
				return 1;
			}

			/* Convert coordinates to match arrangement of panels in HDF5
			 * file */
			write_fs = fs - p->min_fs + p->orig_min_fs;
			write_ss = ss - p->min_ss + p->orig_min_ss;

			fprintf(fh, "%4i %4i %4i %10.2f %10.2f %10.2f %10.2f"
			            " %6.1f %6.1f\n",
			        h, k, l, intensity, esd_i, pk, bg, write_fs,
			        write_ss);

		} else {

			fprintf(fh, "%4i %4i %4i %10.2f %10.2f %10.2f %10.2f"
			            " %6.1f %6.1f\n",
			        h, k, l, intensity, esd_i, pk, bg, fs, ss);
		}
	}
	return 0;
}


static int write_stream_reflections_2_3(FILE *fh, RefList *list,
                                        struct image *image)
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
		double fs, ss;
		double write_fs, write_ss;
		struct panel *p = NULL;

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &fs, &ss);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		pk = get_peak(refl);
		bg = get_mean_bg(refl);

		/* Reflections with redundancy = 0 are not written */
		if ( get_redundancy(refl) == 0 ) continue;

		p = find_panel(image->det,fs,ss);
		if ( p == NULL ) {
			ERROR("Panel not found\n");
			return 1;
		}

		write_fs = fs-p->min_fs+p->orig_min_fs;
		write_ss = ss-p->min_ss+p->orig_min_ss;

		fprintf(fh,
                          "%4i %4i %4i %10.2f %10.2f %10.2f %10.2f "
                          "%6.1f %6.1f %s\n",
                           h, k, l, intensity, esd_i, pk, bg,
                           write_fs, write_ss, p->name);

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


static int write_crystal(Stream *st, Crystal *cr, int include_reflections)
{
	UnitCell *cell;
	RefList *reflist;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;
	double rad;
	int ret = 0;

	fprintf(st->fh, CRYSTAL_START_MARKER"\n");

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

			struct image *image;

			image = crystal_get_image(cr);

			fprintf(st->fh, REFLECTION_START_MARKER"\n");
			if ( AT_LEAST_VERSION(st, 2, 3) ) {
				ret = write_stream_reflections_2_3(st->fh,
				                                   reflist,
				                                   image);
			} else if ( AT_LEAST_VERSION(st, 2, 2) ) {
				ret = write_stream_reflections_2_2(st->fh,
			                                           reflist,
			                                           image);
			} else {
				/* This function writes like a normal reflection
				 * list was written in stream 2.1 */
				ret = write_stream_reflections_2_1(st->fh,
				                                   reflist,
				                                   image);
			}
			fprintf(st->fh, REFLECTION_END_MARKER"\n");

		} else {

			fprintf(st->fh, "No integrated reflections.\n");

		}
	}

	fprintf(st->fh, CRYSTAL_END_MARKER"\n");

	return ret;
}


int write_chunk(Stream *st, struct image *i, struct hdfile *hdfile,
                int include_peaks, int include_reflections, struct event* ev)
{
	int j;
	char *indexer;
	int ret = 0;

	fprintf(st->fh, CHUNK_START_MARKER"\n");

	fprintf(st->fh, "Image filename: %s\n", i->filename);
	if ( i->event != NULL ) {
		fprintf(st->fh, "Event: %s\n", get_event_string(i->event));
	}

	fprintf(st->fh, "Image serial number: %i\n", i->serial);

	indexer = indexer_str(i->indexed_by);
	fprintf(st->fh, "indexed_by = %s\n", indexer);
	free(indexer);

	fprintf(st->fh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(i->lambda)));

	fprintf(st->fh, "beam_divergence = %.2e rad\n", i->div);
	fprintf(st->fh, "beam_bandwidth = %.2e (fraction)\n", i->bw);

	copy_hdf5_fields(hdfile, i->copyme, st->fh, ev);

	if ( i->det != NULL ) {

		int j;
		double tclen = 0.0;

		for ( j=0; j<i->det->n_panels; j++ ) {
			tclen += i->det->panels[j].clen;
		}
		fprintf(st->fh, "average_camera_length = %f m\n",
		        tclen / i->det->n_panels);

		for ( j=0; j<i->det->n_rigid_groups; j++ ) {

			struct rigid_group *rg = i->det->rigid_groups[j];

			if ( !rg->have_deltas ) continue;

			fprintf(st->fh, "rg_delta_%s_fsx = %f\n",
			        rg->name, rg->d_fsx);
			fprintf(st->fh, "rg_delta_%s_ssx = %f\n",
			        rg->name, rg->d_ssx);
			fprintf(st->fh, "rg_delta_%s_cnx = %f\n",
			        rg->name, rg->d_cnx);

			fprintf(st->fh, "rg_delta_%s_fsy = %f\n",
			        rg->name, rg->d_fsy);
			fprintf(st->fh, "rg_delta_%s_ssy = %f\n",
			        rg->name, rg->d_ssy);
			fprintf(st->fh, "rg_delta_%s_cny = %f\n",
			        rg->name, rg->d_cny);

		}

	}

	fprintf(st->fh, "num_peaks = %lli\n", i->num_peaks);
	fprintf(st->fh, "num_saturated_peaks = %lli\n", i->num_saturated_peaks);
	if ( include_peaks ) {
		if ( AT_LEAST_VERSION(st, 2, 3) ) {
			ret = write_peaks_2_3(i, st->fh);
		} else {
			ret = write_peaks(i, st->fh);
		}
	}

	for ( j=0; j<i->n_crystals; j++ ) {
		if ( crystal_get_user_flag(i->crystals[j]) == 0 ) {
			ret = write_crystal(st, i->crystals[j],
			                    include_reflections);
		}
	}

	fprintf(st->fh, CHUNK_END_MARKER"\n");

	fflush(st->fh);

	return ret;
}


static int find_start_of_chunk(FILE *fh)
{
	char *rval = NULL;
	char line[1024];

	do {

		rval = fgets(line, 1023, fh);

		/* Trouble? */
		if ( rval == NULL ) return 1;

		chomp(line);

	} while ( strcmp(line, CHUNK_START_MARKER) != 0 );

	return 0;
}


static void read_crystal(Stream *st, struct image *image, StreamReadFlags srf)
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

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);
		if ( (srf & STREAM_READ_UNITCELL)
		  && (sscanf(line, "astar = %f %f %f", &u, &v, &w) == 3) )
		{
			as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
			have_as = 1;
		}

		if ( (srf & STREAM_READ_UNITCELL)
		  && (sscanf(line, "bstar = %f %f %f", &u, &v, &w) == 3) )
		{
			bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
			have_bs = 1;
		}

		if ( (srf & STREAM_READ_UNITCELL)
		  && (sscanf(line, "cstar = %f %f %f", &u, &v, &w) == 3) )
		{
			cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
			have_cs = 1;
		}

		if ( (srf & STREAM_READ_UNITCELL)
		  && (sscanf(line, "centering = %c", &c) == 1) )
		{
			if ( !have_cen ) {
				centering = c;
				have_cen = 1;
			} else {
				ERROR("Duplicate centering ignored.\n");
			}
		}

		if ( (srf & STREAM_READ_UNITCELL)
		  && (sscanf(line, "unique_axis = %c", &c) == 1) )
		{
			if ( !have_ua ) {
				unique_axis = c;
				have_ua = 1;
			} else {
				ERROR("Duplicate unique axis ignored.\n");
			}
		}

		if ( (srf & STREAM_READ_UNITCELL)
		  && (strncmp(line, "lattice_type = ", 15) == 0) )
		{
			if ( !have_latt ) {
				lattice_type = lattice_from_str(line+15);
				have_latt = 1;
			} else {
				ERROR("Duplicate lattice type ignored.\n");
			}
		}

		if ( strncmp(line, "num_saturated_reflections = ", 28) == 0 ) {
			int n = atoi(line+28);
			crystal_set_num_saturated_reflections(cr, n);
		}

		if ( sscanf(line, "diffraction_resolution_limit = %f nm^-1",
		            &lim) == 1 ) {
			crystal_set_resolution_limit(cr, lim*1e9);
		}

		if ( sscanf(line, "profile_radius = %e nm^-1", &rad) == 1 ) {
			crystal_set_profile_radius(cr, rad*1e9);
		}

		if ( (strcmp(line, REFLECTION_START_MARKER) == 0)
		  && (srf & STREAM_READ_REFLECTIONS) )
		{

			RefList *reflist;

			/* The reflection list format in the stream diverges
			 * after 2.2 */
			if ( AT_LEAST_VERSION(st, 2, 3) ) {
				reflist = read_stream_reflections_2_3(st->fh,
                                                      image->det);
			} else if ( AT_LEAST_VERSION(st, 2, 2) ) {
				reflist = read_stream_reflections_2_2(st->fh,
				          image->det);
			} else {
				reflist = read_stream_reflections_2_1(st->fh,
				          image->det);
			}
			if ( reflist == NULL ) {
				ERROR("Failed while reading reflections\n");
				break;
			}

			crystal_set_reflections(cr, reflist);

		}

		if ( strcmp(line, CRYSTAL_END_MARKER) == 0 ) break;

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

		if ( have_cen && have_ua && have_latt ) {
			cell_set_centering(cell, centering);
			cell_set_unique_axis(cell, unique_axis);
			cell_set_lattice_type(cell, lattice_type);
		} /* else keep default triclinic P */

		crystal_set_cell(cr, cell);

		have_as = 0;  have_bs = 0;  have_cs = 0;
		have_latt = 0;  have_ua = 0;  have_cen = 0;

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


static int read_and_store_hdf5_field(struct image *image, const char *line)
{

	char **new_fields;

	if ( image->stuff_from_stream == NULL ) {
		image->stuff_from_stream =
		       malloc(sizeof(struct stuff_from_stream));
		if ( image->stuff_from_stream == NULL) {
			ERROR("Failed reading hdf5 entries from "
			      "stream\n");
			return 1;
		}
		image->stuff_from_stream->fields = NULL;
		image->stuff_from_stream->n_fields = 0;
	}

	new_fields = realloc(image->stuff_from_stream->fields,
			     (1+image->stuff_from_stream->n_fields)*
			     sizeof(char *));
	if ( new_fields == NULL ) {
		ERROR("Failed reading hdf5 entries from stream\n");
		return 1;
	}
	image->stuff_from_stream->fields = new_fields;
	image->stuff_from_stream->fields[image->stuff_from_stream->n_fields]
							     = strdup(line);
	image->stuff_from_stream->n_fields++;

	return 0;
}


/* Read the next chunk from a stream and fill in 'image' */
int read_chunk_2(Stream *st, struct image *image,  StreamReadFlags srf)
{
	char line[1024];
	char *rval = NULL;
	int have_filename = 0;
	int have_ev = 0;

	if ( find_start_of_chunk(st->fh) ) return 1;

	image->lambda = -1.0;
	image->features = NULL;
	image->crystals = NULL;
	image->n_crystals = 0;
	image->event = NULL;
	image->stuff_from_stream = NULL;

	if ( (srf & STREAM_READ_REFLECTIONS) || (srf & STREAM_READ_UNITCELL) ) {
		srf |= STREAM_READ_CRYSTALS;
	}

	do {
		long long num_peaks;
		int ser;
		float div, bw;

		rval = fgets(line, 1023, st->fh);

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);

		if ( strncmp(line, "Image filename: ", 16) == 0 ) {
			image->filename = strdup(line+16);
			have_filename = 1;
		}

		if ( strncmp(line, "Event: ", 7) == 0 ) {
			image->event = get_event_from_event_string(line+7);
		}

		if ( strncmp(line, "indexed_by = ", 13) == 0 ) {
			IndexingMethod *list;
			list = build_indexer_list(line+13);
			image->indexed_by = list[0];
			free(list);
			have_filename = 1;
		}

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			image->lambda = ph_en_to_lambda(eV_to_J(atof(line+19)));
			have_ev = 1;
		}

		if ( sscanf(line, "beam_divergence = %e rad", &div) == 1 ) {
			image->div = div;
		}

		if ( sscanf(line, "beam_bandwidth = %f %%", &bw) == 1 ) {
			image->bw = bw/100.0;
		}

		if ( sscanf(line, "num_peaks = %lld %%", &num_peaks) == 1 ) {
			image->num_peaks = num_peaks;
		}

		if ( sscanf(line, "Image serial number: %i", &ser) == 1 ) {
			image->serial = ser;
		}

		if ( strncmp(line, "camera_length_", 14) == 0 ) {
			if ( image->det != NULL ) {

				int k;
				char name[1024];
				struct panel *p;

				for ( k=0; k<strlen(line)-14; k++ ) {
					char ch = line[k+14];
					name[k] = ch;
					if ( (ch == ' ') || (ch == '=') ) {
						name[k] = '\0';
						break;
					}
				}

				p = find_panel_by_name(image->det, name);
				if ( p == NULL ) {
					ERROR("No panel '%s'\n", name);
				} else {
					p->clen = atof(line+14+k+3);
				}

			}
		}

		if ( strncmp(line, "hdf5", 3) == 0 ) {

			int fail;

			fail = read_and_store_hdf5_field(image, line);
			if ( fail ) {
				ERROR("Failed to read hd5 fields from stream.\n");
				return 1;
			}
		}

		if ( (srf & STREAM_READ_PEAKS)
		    && strcmp(line, PEAK_LIST_START_MARKER) == 0 ) {

			int fail;

			if ( AT_LEAST_VERSION(st, 2, 3) ) {
				fail = read_peaks_2_3(st->fh, image);
			} else {
				fail = read_peaks(st->fh, image);
			}
			if ( fail ) {
				ERROR("Failed while reading peaks\n");
				return 1;
			}
		}

		if ( (srf & STREAM_READ_CRYSTALS)
		  && (strcmp(line, CRYSTAL_START_MARKER) == 0) ) {
			read_crystal(st, image, srf);
		}

		/* A chunk must have at least a filename and a wavelength,
		 * otherwise it's incomplete */
		if ( strcmp(line, CHUNK_END_MARKER) == 0 ) {
			if ( have_filename && have_ev ) return 0;
			ERROR("Incomplete chunk found in input file.\n");
			return 1;
		}

	} while ( 1 );

	if ( !feof(st->fh) ) {
		ERROR("Error reading stream.\n");
	}

	return 1;  /* Either error or EOF, don't care because we will complain
	            * on the terminal if it was an error. */
}


int read_chunk(Stream *st, struct image *image)
{
	return read_chunk_2(st, image, STREAM_READ_UNITCELL
	                               | STREAM_READ_REFLECTIONS
	                               | STREAM_READ_PEAKS);
}


void write_stream_header(FILE *ofh, int argc, char *argv[])
{
	int i;

	fprintf(ofh, "Command line:");
	for ( i=0; i<argc; i++ ) {
		fprintf(ofh, " %s", argv[i]);
	}
	fprintf(ofh, "\n");
	fflush(ofh);
}


Stream *open_stream_for_read(const char *filename)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;

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
		close_stream(st);
		return NULL;
	}

	if ( strncmp(line, "CrystFEL stream format 2.0", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 0;
	} else if ( strncmp(line, "CrystFEL stream format 2.1", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 1;
	} else if ( strncmp(line, "CrystFEL stream format 2.2", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 2;
	} else if ( strncmp(line, "CrystFEL stream format 2.3", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 3;
	} else {
		ERROR("Invalid stream, or stream format is too new.\n");
		close_stream(st);
		return NULL;
	}

	return st;
}


/**
 * open_stream_fd_for_write
 * @fd: File descriptor (e.g. from open()) to use for stream data.
 *
 * Creates a new %Stream from @fd, so that stream data can be written to @fd
 * using write_chunk().
 *
 * In contrast to open_stream_for_write(), this function does not write any of
 * the usual headers.  This function is mostly for use when multiple substreams
 * need to be multiplexed into a single master stream.  The master would be
 * opened using open_stream_for_write(), and the substreams using this function.
 *
 * Returns: a %Stream, or NULL on failure.
 */
Stream *open_stream_fd_for_write(int fd)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;

	st->fh = fdopen(fd, "w");
	if ( st->fh == NULL ) {
		free(st);
		return NULL;
	}

	st->major_version = LATEST_MAJOR_VERSION;
	st->minor_version = LATEST_MINOR_VERSION;

	return st;
}


/**
 * open_stream_for_write_2
 * @filename: Filename of new stream
 * @geom_filename: The geometry filename to copy
 * @argc: The number of arguments to the program
 * @argv: The arguments to the program
 *
 * Creates a new stream with name @filename, and adds the stream format
 * and version header, plus a verbatim copy of the geometry file
 *
 * You may want to follow this with a call to write_command() to record the
 * command line.
 *
 * Returns: a %Stream, or NULL on failure.
 */
Stream *open_stream_for_write_2(const char *filename,
                                const char *geom_filename, int argc,
                                char *argv[])

{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;

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
	fprintf(st->fh, "Generated by CrystFEL "CRYSTFEL_VERSIONSTRING"\n");
	fflush(st->fh);

	if ( (argc > 0) && (argv != NULL) ) {
		write_command(st, argc, argv);
	}
	if ( geom_filename != NULL ) {
		write_geometry_file(st, geom_filename);
	}

	return st;
}


/**
 * open_stream_for_write
 * @filename: Filename of new stream
 *
 * Creates a new stream with name @filename, and adds the stream format
 * and version headers.
 *
 * You may want to follow this with a call to write_command() to record the
 * command line.
 *
 * Returns: a %Stream, or NULL on failure.
 */
Stream *open_stream_for_write(const char *filename)
{
	return open_stream_for_write_2(filename, NULL, 0, NULL);
}


/**
 * get_stream_fd
 * @st: A %Stream
 *
 * This function gets the integer file descriptor for @st, a bit like fileno().
 *
 * This is useful in conjunction with open_stream_fd_for_write(), to get the
 * underlying file descriptor to which the multiplexed stream data should be
 * written.  In this case, the only other operations you should ever do (or have
 * done) on @st are open_stream_for_write() and close_stream().
 *
 * Returns: an integer file descriptor
 */
int get_stream_fd(Stream *st)
{
	return fileno(st->fh);
}


void close_stream(Stream *st)
{
	fclose(st->fh);
	free(st);
}


int is_stream(const char *filename)
{
	FILE *fh;
	char line[1024];
	char *rval;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return 0;

	rval = fgets(line, 1023, fh);
	fclose(fh);
	if ( rval == NULL ) return 0;

	if ( strncmp(line, "CrystFEL stream format 2.0", 26) == 0 ) return 1;
	if ( strncmp(line, "CrystFEL stream format 2.1", 26) == 0 ) return 1;
	if ( strncmp(line, "CrystFEL stream format 2.2", 26) == 0 ) return 1;

	return 0;
}


/**
 * write_command
 * @st: A %Stream
 * @argc: number of arguments
 * @argv: command-line arguments
 *
 * Writes the command line to @st.  @argc and @argv should be exactly as were
 * given to main().  This should usually be called immediately after
 * open_stream_for_write().
 */
void write_command(Stream *st, int argc, char *argv[])
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


/**
 * write_geometry_file
 * @st: A %Stream
 * @geom_filename: geomtry file name
 *
 * Writes the content of the geometry file to @st. This should usually be
 * called immediately after write_command().
 */
void write_geometry_file(Stream *st, const char *geom_filename) {

	char line[2014];
	FILE *geom_fh;
	char *rval;

	if ( geom_filename == NULL ) return;

	geom_fh = fopen(geom_filename, "r");
	if ( geom_fh == NULL ) {
		ERROR("Failed to read detector geometry from "
		      "'%s'\n", geom_filename);
		return;
	}
	fprintf(st->fh, GEOM_START_MARKER"\n");

	do {
		rval = fgets(line, 1023, geom_fh);
		if ( rval != NULL ) fputs(line, st->fh);
	} while ( rval != NULL );

	fclose(geom_fh);

	fprintf(st->fh, GEOM_END_MARKER"\n");
	fflush(st->fh);
}


/**
 * rewind_stream:
 * @st: A %Stream
 *
 * Attempts to set the file pointer for @st to the start of the stream, so that
 * later calls to read_chunk() will repeat the sequence of chunks from the
 * start.
 *
 * Programs must not assume that this operation always succeeds!
 *
 * Returns: non-zero if the stream could not be rewound.
 */
int rewind_stream(Stream *st)
{
	return fseek(st->fh, 0, SEEK_SET);
}

