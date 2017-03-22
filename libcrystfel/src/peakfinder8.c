/*
 * peakfinder8.c
 *
 * The peakfinder8 algorithm
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017  Valerio Mariani <valerio.mariani@desy.de>
 *   2017  Anton Barty <anton.barty@desy.de>
 *   2017  Oleksandr Yefanov <oleksandr.yefanov@desy.de>
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


#include "peakfinder8.h"


struct radius_maps
{
	float **r_maps;
	int n_rmaps;
};


struct peakfinder_mask
{
	char **masks;
	int n_masks;
};


struct radial_stats
{
	float *roffset;
	float *rthreshold;
	float *rsigma;
	int *rcount;
	int n_rad_bins;
};


struct peakfinder_panel_data
{
	float ** panel_data;
	int *panel_h;
	int *panel_w;
	int num_panels;
};


struct peakfinder_intern_data
{
	char *pix_in_peak_map;
	int *infs;
	int *inss;
	int *peak_pixels;
};


struct peakfinder_peak_data
{
	int num_found_peaks;
	int *npix;
	float *com_fs;
	float *com_ss;
	int *com_index;
	float *tot_i;
	float *max_i;
	float *sigma;
	float *snr;
};


static double compute_r_sigma(float *rsigma, int *rcount, float* roffset,
                              int i)
{
	return sqrt(rsigma[i] / rcount[i] -
	            ((roffset[i] / rcount[i]) *
	             (roffset[i] / rcount[i])));
}


static void update_radial_stats (float *roffset, float *rsigma, int *rcount,
                                 float value, int curr_radius)
{
	roffset[curr_radius] += value;
	rsigma[curr_radius] += (value * value);
	rcount[curr_radius] += 1;
}


static float get_radius(struct panel p, float fs, float ss)
{
	float x, y;

	/* Calculate 3D position of given position, in m */
	x = (p.cnx  + fs * p.fsx + ss * p.ssx);
	y = (p.cny  + fs * p.fsy + ss * p.ssy);

	return sqrt(x * x + y * y);
}


static struct radius_maps* compute_radius_maps(struct detector *det)
{
	int i, u, iss, ifs;
	struct panel p;
	struct radius_maps *rm = NULL;

	rm = (struct radius_maps *)malloc(sizeof(struct radius_maps));
	if ( rm == NULL ) {
		return NULL;
	}

	rm->r_maps = (float **)malloc(det->n_panels*sizeof(float*));
	if ( rm->r_maps == NULL ) {
		free(rm);
		return NULL;
	}

	rm->n_rmaps = det->n_panels;

	for( i=0 ; i<det->n_panels ; i++ ) {

		p = det->panels[i];
		rm->r_maps[i] = (float *)malloc(p.h*p.w*sizeof(float));

		if ( rm->r_maps[i] == NULL ) {
			for ( u = 0; u<i; u++ ) {
				free(rm->r_maps[u]);
			}
			free(rm);
			return NULL;
		}

		for ( iss=0 ; iss<p.h ; iss++ ) {
			for ( ifs=0; ifs<p.w; ifs++ ) {

				int rmi;

				rmi = ifs + p.w * iss;

				rm->r_maps[i][rmi] = get_radius(p, ifs, iss);
			}
		}

	}

	return rm;
}


static void free_radius_maps(struct radius_maps *r_maps)
{
	int i;

	for ( i=0 ; i<r_maps->n_rmaps ; i++ ) {
		free(r_maps->r_maps[i]);
	}
	free(r_maps->r_maps);
	free(r_maps);
}


static struct peakfinder_mask *create_peakfinder_mask(struct image *img,
                                                      struct radius_maps *rmps,
                                                      int min_res,
                                                      int max_res)
{
	int i;
	struct peakfinder_mask *msk;

	msk = (struct peakfinder_mask *)malloc(sizeof(struct peakfinder_mask));
	msk->masks =(char **) malloc(img->det->n_panels*sizeof(char*));
	msk->n_masks = img->det->n_panels;
	for ( i=0; i<img->det->n_panels; i++) {

		struct panel p;
		int iss, ifs;

		p = img->det->panels[i];

		msk->masks[i] = (char *)calloc(p.w*p.h,sizeof(char));

		for ( iss=0 ; iss<p.h ; iss++ ) {
			for ( ifs=0 ; ifs<p.w ; ifs++ ) {

				int idx;

				idx = ifs + iss*p.w;

				if ( rmps->r_maps[i][idx] < max_res
				  && rmps->r_maps[i][idx] > min_res ) {

					if  (! ( ( img->bad != NULL )
					      && ( img->bad[i] != NULL )
					      && ( img->bad[i][idx] != 0 ) ) ) {
						msk->masks[i][idx] = 1;
					}

				}

			}
		}
	}
	return msk;
}


static void free_peakfinder_mask(struct peakfinder_mask * pfmask)
{
	int i;

	for ( i=0 ; i<pfmask->n_masks ; i++ ) {
		free(pfmask->masks[i]);
	}
	free(pfmask->masks);
	free(pfmask);
}


static int compute_num_radial_bins(int num_panels, int *w, int *h,
                                   float **r_maps )
{

	float max_r;
	int max_r_int;
	int m;

	max_r = -1e9;

	for ( m=0 ; m<num_panels ; m++ ) {

		int ifs, iss;
		int pidx;

		for ( iss=0 ; iss<h[m] ; iss++ ) {
			for ( ifs=0 ; ifs<w[m] ; ifs++ ) {
				pidx = iss * w[m] + ifs;
				if ( r_maps[m][pidx] > max_r ) {
					max_r = r_maps[m][pidx];
				}
			}
		}
	}

	max_r_int = (int)ceil(max_r) + 1;

	return max_r_int;
}


static struct peakfinder_panel_data* allocate_panel_data(int num_panels)
{

	struct peakfinder_panel_data *pfdata;


	pfdata = (struct peakfinder_panel_data *)malloc(
	                        sizeof(struct peakfinder_panel_data));
	if ( pfdata == NULL ) {
		return NULL;
	}

	pfdata->panel_h = (int *)malloc(num_panels*sizeof(int));
	if ( pfdata->panel_h == NULL ) {
		free(pfdata);
		return NULL;
	}

	pfdata->panel_w = (int *)malloc(num_panels*sizeof(int));
	if ( pfdata->panel_w == NULL ) {
		free(pfdata->panel_h);
		free(pfdata);
		return NULL;
	}

	pfdata->panel_data = (float **)malloc(num_panels*sizeof(float*));
	if ( pfdata->panel_data == NULL ) {
		free(pfdata->panel_w);
		free(pfdata->panel_h);
		free(pfdata);
		return NULL;
	}

	pfdata->num_panels = num_panels;

	return pfdata;
}


static void free_panel_data(struct peakfinder_panel_data *pfdata)
{
	free(pfdata->panel_data);
	free(pfdata->panel_w);
	free(pfdata->panel_h);
	free(pfdata);
}


static struct radial_stats* allocate_radial_stats(int num_rad_bins)
{
	struct radial_stats* rstats;

	rstats = (struct radial_stats *)malloc(sizeof(struct radial_stats));
	if ( rstats == NULL ) {
		return NULL;
	}

	rstats->roffset = (float *)malloc(num_rad_bins*sizeof(float));
	if ( rstats->roffset == NULL ) {
		free(rstats);
		return NULL;
	}

	rstats->rthreshold = (float *)malloc(num_rad_bins*sizeof(float));
	if ( rstats->rthreshold == NULL ) {
		free(rstats);
		free(rstats->roffset);
		return NULL;
	}

	rstats->rsigma = (float *)malloc(num_rad_bins*sizeof(float));
	if ( rstats->rsigma == NULL ) {
		free(rstats);
		free(rstats->roffset);
		free(rstats->rthreshold);
		return NULL;
	}

	rstats->rcount = (int *)malloc(num_rad_bins*sizeof(int));
	if ( rstats->rcount == NULL ) {
		free(rstats);
		free(rstats->roffset);
		free(rstats->rthreshold);
		free(rstats->rsigma);
		return NULL;
	}

	rstats->n_rad_bins = num_rad_bins;

	return rstats;
}


static void free_radial_stats(struct radial_stats *rstats)
{
	free(rstats->roffset);
	free(rstats->rthreshold);
	free(rstats->rsigma);
	free(rstats->rcount);
	free(rstats);
}


static void fill_radial_bins(float *data,
                             int w,
                             int h,
                             float *r_map,
                             char *mask,
                             float *rthreshold,
                             float *roffset,
                             float *rsigma,
                             int *rcount)
{
	int iss, ifs;
	int pidx;

	int curr_r;
	float value;

	for ( iss = 0; iss<h ; iss++ ) {

		for ( ifs = 0; ifs<w ; ifs++ ) {

			pidx = iss * w + ifs;

			if ( mask[pidx] != 0 ) {

				curr_r = (int)rint(r_map[pidx]);

				value = data[pidx];

				if ( value < rthreshold[curr_r ] ) {

					update_radial_stats(roffset,
					                    rsigma,
					                    rcount,
					                    value,
					                    curr_r);
				}
			}
		}
	}
}


static void compute_radial_stats(float *rthreshold,
                                 float *roffset,
                                 float *rsigma,
                                 int *rcount,
                                 int num_rad_bins,
                                 float min_snr,
                                 float acd_threshold)
{

	int ri;
	float this_offset, this_sigma;

	for ( ri=0 ; ri<num_rad_bins ; ri++ ) {

		if ( rcount[ri] == 0 ) {
			roffset[ri] = 0;
			rsigma[ri] = 0;
			rthreshold[ri] = 1e9;
		} else {
			this_offset  = roffset[ri]/rcount[ri];

			this_sigma = compute_r_sigma(rsigma,
			                             rcount,
			                             roffset,
			                             ri);

			roffset[ri] = this_offset;
			rsigma[ri] = this_sigma;

			rthreshold[ri] = roffset[ri] +
			                min_snr*rsigma[ri];

			if ( rthreshold[ri] < acd_threshold ) {
				rthreshold[ri] = acd_threshold;
			}
		}
	}

}


struct peakfinder_peak_data *allocate_peak_data(int max_num_peaks)
{
	struct peakfinder_peak_data *pkdata;

	pkdata = (struct peakfinder_peak_data*)malloc(
	                        sizeof(struct peakfinder_peak_data));
	if ( pkdata == NULL ) {
		return NULL;
	}

	pkdata->npix = (int *)malloc(max_num_peaks*sizeof(int));
	if ( pkdata->npix == NULL ) {
		free(pkdata);
		free(pkdata->npix);
		return NULL;
	}

	pkdata->com_fs = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->com_fs == NULL ) {
		free(pkdata->npix);
		free(pkdata);
		return NULL;
	}

	pkdata->com_ss = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->com_ss == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata);
		return NULL;
	}

	pkdata->com_index = (int *)malloc(max_num_peaks*sizeof(int));
	if ( pkdata->com_ss == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata->com_ss);
		free(pkdata);
		return NULL;
	}


	pkdata->tot_i = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->tot_i == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata->com_ss);
		free(pkdata->com_index);
		free(pkdata);
		return NULL;
	}

	pkdata->max_i = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->max_i == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata->com_ss);
		free(pkdata->com_index);
		free(pkdata->tot_i);
		free(pkdata);
		return NULL;
	}

	pkdata->sigma = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->sigma == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata->com_ss);
		free(pkdata->com_index);
		free(pkdata->tot_i);
		free(pkdata->max_i);
		free(pkdata);
		return NULL;
	}

	pkdata->snr = (float *)malloc(max_num_peaks*sizeof(float));
	if ( pkdata->snr == NULL ) {
		free(pkdata->npix);
		free(pkdata->com_fs);
		free(pkdata->com_ss);
		free(pkdata->com_index);
		free(pkdata->tot_i);
		free(pkdata->max_i);
		free(pkdata->sigma);
		free(pkdata);
		return NULL;
	}

	return pkdata;
}


static void free_peak_data(struct peakfinder_peak_data *pkdata) {
	free(pkdata->npix);
	free(pkdata->com_fs);
	free(pkdata->com_ss);
	free(pkdata->com_index);
	free(pkdata->tot_i);
	free(pkdata->max_i);
	free(pkdata->sigma);
	free(pkdata->snr);
	free(pkdata);
}


static struct peakfinder_intern_data *allocate_peakfinder_intern_data(
                                     int data_size, int max_pix_count)
{

	struct peakfinder_intern_data *intern_data;

	intern_data = (struct peakfinder_intern_data *)malloc(
	              sizeof(struct peakfinder_intern_data));
	if ( intern_data == NULL ) {
		return NULL;
	}

	intern_data->pix_in_peak_map =(char *)calloc(data_size, sizeof(char));
	if ( intern_data->pix_in_peak_map == NULL ) {
		free(intern_data);
		return NULL;
	}

	intern_data->infs =(int *)calloc(data_size, sizeof(int));
	if ( intern_data->infs == NULL ) {
		free(intern_data->pix_in_peak_map);
		free(intern_data);
		return NULL;
	}

	intern_data->inss =(int *)calloc(data_size, sizeof(int));
	if ( intern_data->inss == NULL ) {
		free(intern_data->pix_in_peak_map);
		free(intern_data->infs);
		free(intern_data);
		return NULL;
	}

	intern_data->peak_pixels =(int *)calloc(max_pix_count, sizeof(int));
	if ( intern_data->peak_pixels == NULL ) {
		free(intern_data->pix_in_peak_map);
		free(intern_data->infs);
		free(intern_data->inss);
		free(intern_data);
		return NULL;
	}

	return intern_data;
}


static void free_peakfinder_intern_data(struct peakfinder_intern_data *pfid)
{
	free(pfid->peak_pixels);
	free(pfid->pix_in_peak_map);
	free(pfid->infs);
	free(pfid->inss);
	free(pfid);
}



static void peak_search(int p,
                        struct peakfinder_intern_data *pfinter,
                        float *copy, char *mask, float *r_map,
                        float *rthreshold, float *roffset,
                        int *num_pix_in_peak, int asic_size_fs,
                        int asic_size_ss, int aifs, int aiss,
                        int num_pix_fs, float *sum_com_fs,
                        float *sum_com_ss, float *sum_i, int max_pix_count)
{

	int k, pi;
	int curr_radius;
	float curr_threshold;
	int curr_fs;
	int curr_ss;
	float curr_i;

	int search_fs[9] = { 0, -1, 0, 1, -1, 1, -1, 0, 1 };
	int search_ss[9] = { 0, -1, -1, -1, 0, 0, 1, 1, 1 };
	int search_n = 9;

	for ( k=0; k<search_n; k++ ) {

		if ( (pfinter->infs[p] + search_fs[k]) < 0 )
			continue;
		if ( (pfinter->infs[p] + search_fs[k]) >= asic_size_fs )
			continue;
		if ( (pfinter->inss[p] + search_ss[k]) < 0 )
			continue;
		if ( (pfinter->inss[p] + search_ss[k]) >= asic_size_ss )
			continue;

		curr_fs = pfinter->infs[p] + search_fs[k] + aifs * asic_size_fs;
		curr_ss = pfinter->inss[p] + search_ss[k] + aiss * asic_size_ss;
		pi = curr_fs + curr_ss * num_pix_fs;

		curr_radius = (int)rint(r_map[pi]);
		curr_threshold = rthreshold[curr_radius];

		if ( copy[pi] > curr_threshold
		  && pfinter->pix_in_peak_map[pi] == 0
		  && mask[pi] != 0 ) {


			curr_i = copy[pi] - roffset[curr_radius];

			*sum_i += curr_i;
			*sum_com_fs += curr_i * ((float)curr_fs);
			*sum_com_ss += curr_i * ((float)curr_ss);

			pfinter->inss[*num_pix_in_peak] = pfinter->inss[p] +
			                                  search_ss[k];
			pfinter->infs[*num_pix_in_peak] = pfinter->infs[p] +
			                                  search_fs[k];

			pfinter->pix_in_peak_map[pi] = 1;
			if ( *num_pix_in_peak < max_pix_count ) {
				pfinter->peak_pixels[*num_pix_in_peak] = pi;
			}
			*num_pix_in_peak = *num_pix_in_peak+1;
		}
	}
}


static void search_in_ring(int ring_width, int com_fs_int, int com_ss_int,
                           float *copy, float *r_map,
                           float *rthreshold, float *roffset,
                           char *pix_in_peak_map, char *mask, int asic_size_fs,
                           int asic_size_ss, int aifs, int aiss,
                           int num_pix_fs,float *local_sigma,
                           float *local_offset, float *background_max_i,
                           int com_idx, int local_bg_radius)
{
	int ssj, fsi;
	float pix_radius;
	int curr_fs, curr_ss;
	int pi;
	int curr_radius;
	float curr_threshold;
	float curr_i;

	int np_sigma;
	int np_counted;
	int local_radius;

	float sum_i;
	float sum_i_squared;

	ring_width = 2 * local_bg_radius;

	sum_i = 0;
	sum_i_squared = 0;
	np_sigma = 0;
	np_counted = 0;
	local_radius = 0;

	for ( ssj = -ring_width ; ssj<ring_width ; ssj++ ) {
		for ( fsi = -ring_width ; fsi<ring_width ; fsi++ ) {

			if ( (com_fs_int + fsi) < 0 )
				continue;
			if ( (com_fs_int + fsi) >= asic_size_fs )
				continue;
			if ( (com_ss_int + ssj) < 0 )
				continue;
			if ( (com_ss_int + ssj) >= asic_size_ss )
				continue;

			pix_radius = sqrt(fsi * fsi + ssj * ssj);
			if ( pix_radius>ring_width )
				continue;

			curr_fs = com_fs_int +fsi + aifs * asic_size_fs;
			curr_ss = com_ss_int +ssj + aiss * asic_size_ss;
			pi = curr_fs + curr_ss * num_pix_fs;

			curr_radius = rint(r_map[pi]);
			curr_threshold = rthreshold[curr_radius];
			curr_i = copy[pi];

			if ( copy[pi] < curr_threshold
			                && pix_in_peak_map[pi] == 0
			                && mask[pi] != 0 ) {

				np_sigma++;
				sum_i += curr_i;
				sum_i_squared += (curr_i * curr_i);

				if ( curr_i > *background_max_i ) {
					*background_max_i = curr_i;
				}
			}
			np_counted += 1;
		}
	}

	if ( np_sigma != 0 ) {
		*local_offset = sum_i / np_sigma;
		*local_sigma = sqrt(sum_i_squared / np_sigma -
		                    ((sum_i / np_sigma) * (sum_i / np_sigma)));
	} else {

		local_radius = rint(r_map[(int)rint(com_idx)]);
		*local_offset = roffset[local_radius];
		*local_sigma = 0.01;
	}
}


static void process_panel(int asic_size_fs, int asic_size_ss, int num_pix_fs,
                          int aiss, int aifs, float *rthreshold,
                          float *roffset, int *peak_count,
                          float *copy, struct peakfinder_intern_data *pfinter,
                          float *r_map, char *mask, int *npix, float *com_fs,
                          float *com_ss, int *com_index, float *tot_i,
                          float *max_i, float *sigma, float *snr,
                          int min_pix_count, int max_pix_count,
                          int local_bg_radius, float min_snr, int max_n_peaks)
{
	int pxss, pxfs;
	int num_pix_in_peak;

	for ( pxss=1 ; pxss<asic_size_ss-1 ; pxss++ ) {
		for ( pxfs=1 ; pxfs<asic_size_fs-1 ; pxfs++ ) {

			float curr_thresh;
			int pxidx;
			int curr_rad;

			pxidx = (pxss + aiss * asic_size_ss) *
			                num_pix_fs + pxfs +
			                aifs * asic_size_fs;

			curr_rad = (int)rint(r_map[pxidx]);
			curr_thresh = rthreshold[curr_rad];

			if ( copy[pxidx] > curr_thresh
			  && pfinter->pix_in_peak_map[pxidx] == 0 ) {

				float sum_com_fs, sum_com_ss;
				float sum_i;
				float peak_com_fs, peak_com_ss;
				float peak_com_fs_int, peak_com_ss_int;
				float peak_tot_i, pk_tot_i_raw;
				float peak_max_i, pk_max_i_raw;
				float peak_snr;
				float local_sigma, local_offset;
				float background_max_i;
				float f_background_thresh;
				int lt_num_pix_in_pk;
				int ring_width;
				int peak_idx;
				int com_idx;
				int p;

				pfinter->infs[0] = pxfs;
				pfinter->inss[0] = pxss;
				pfinter->peak_pixels[0] = pxidx;
				num_pix_in_peak = 1;

				sum_i = 0;
				sum_com_fs = 0;
				sum_com_ss = 0;

				do {
					lt_num_pix_in_pk = num_pix_in_peak;

					for ( p=0; p<num_pix_in_peak; p++ ) {
						peak_search(p,
						            pfinter, copy, mask,
						            r_map,
						            rthreshold,
						            roffset,
						            &num_pix_in_peak,
						            asic_size_fs,
						            asic_size_ss,
						            aifs, aiss,
						            num_pix_fs,
						            &sum_com_fs,
						            &sum_com_ss,
						            &sum_i,
						            max_pix_count);
					}

				} while ( lt_num_pix_in_pk != num_pix_in_peak );

				if ( num_pix_in_peak < min_pix_count
				  || num_pix_in_peak > max_pix_count) {
					continue;
				}

				peak_com_fs = sum_com_fs / fabs(sum_i);
				peak_com_ss = sum_com_ss / fabs(sum_i);

				com_idx = rint(peak_com_fs) +
				                rint(peak_com_ss) * num_pix_fs;

				peak_com_fs_int = (int)rint(peak_com_fs) -
				                  aifs * asic_size_fs;
				peak_com_ss_int = (int)rint(peak_com_ss) -
				                  aiss * asic_size_ss;

				local_sigma = 0.0;
				local_offset = 0.0;
				background_max_i = 0.0;

				ring_width = 2 * local_bg_radius;

				search_in_ring(ring_width, peak_com_fs_int,
				               peak_com_ss_int,
				               copy, r_map, rthreshold,
				               roffset,
				               pfinter->pix_in_peak_map,
				               mask, asic_size_fs,
				               asic_size_ss,
				               aifs, aiss,
				               num_pix_fs,
				               &local_sigma,
				               &local_offset,
				               &background_max_i,
				               com_idx, local_bg_radius);

				peak_tot_i = 0;
				pk_tot_i_raw = 0;
				peak_max_i = 0;
				pk_max_i_raw = 0;
				sum_com_fs = 0;
				sum_com_ss = 0;

				for ( peak_idx = 1 ;
				      peak_idx < num_pix_in_peak
				   && peak_idx <= max_pix_count ;
				      peak_idx++ ) {

					int curr_idx;
					float curr_i;
					float curr_i_raw;
					int curr_fs, curr_ss;

					// BUG HERE, I THINK. PEAK_PIXELS
					// HAS SIZE MAX_PIX_COUNT, BUT
					// IN THE FOLLOWING LINE PEAK_IDX
					// CAN HAVE VALUE EXACTLY MAX_PEAK_COUNT
					// (SEE THE FOR LOOP ABOVE)
					curr_idx =
					         pfinter->peak_pixels[peak_idx];

					curr_i_raw = copy[curr_idx];
					curr_i = curr_i_raw - local_offset;

					peak_tot_i += curr_i;
					pk_tot_i_raw += curr_i_raw;

					curr_fs = curr_idx % asic_size_fs;
					curr_ss = curr_idx / asic_size_fs;
					sum_com_fs += curr_i * ((float)curr_fs);
					sum_com_ss += curr_i * ((float)curr_ss);

					if ( curr_i_raw > pk_max_i_raw )
						pk_max_i_raw = curr_i_raw;
					if ( curr_i > peak_max_i )
						peak_max_i = curr_i;
				}


				peak_com_fs = sum_com_fs / fabs(peak_tot_i);
				peak_com_ss = sum_com_ss / fabs(peak_tot_i);

				peak_snr = peak_tot_i / local_sigma;

				if (peak_snr < min_snr) {
					continue;
				}

				f_background_thresh = 1;
				f_background_thresh *=
				                background_max_i - local_offset;
				if (peak_max_i < f_background_thresh) {
					continue;
				}

				if ( num_pix_in_peak >= min_pix_count
				  && num_pix_in_peak <= max_pix_count ) {

					int peak_com_idx;

					if ( peak_tot_i == 0 ) {
						continue;
					}

					peak_com_idx = rint(peak_com_fs) +
					               rint(peak_com_ss) *
					               asic_size_fs;

					if ( *peak_count < max_n_peaks ) {

						int pidx;
						pidx = *peak_count;

						npix[pidx] = num_pix_in_peak;
						com_fs[pidx] = peak_com_fs;
						com_ss[pidx] = peak_com_ss;
						com_index[pidx] = peak_com_idx;
						tot_i[pidx] = peak_tot_i;
						max_i[pidx] = peak_max_i;
						sigma[pidx] = local_sigma;
						snr[pidx] = peak_snr;
						*peak_count = *peak_count + 1;
					} else {
						*peak_count = *peak_count + 1;
					}
				}
			}
		}
	}
}


static int peakfinder8_base(float *roffset, float *rthreshold,
                            float *data, char *mask, float *r_map,
                            int asic_size_fs, int num_asics_fs,
                            int asic_size_ss, int num_asics_ss,
                            int max_n_peaks, int *num_found_peaks,
                            int *npix, float *com_fs,
                            float *com_ss, int *com_index, float *tot_i,
                            float *max_i, float *sigma, float *snr,
                            int min_pix_count, int max_pix_count,
                            int local_bg_radius, float min_snr)
{
	int num_pix_fs, num_pix_ss, num_pix_tot;
	int i, aifs, aiss;
	int peak_count;
	float *copy;
	struct peakfinder_intern_data *pfinter;

	num_pix_fs = asic_size_fs * num_asics_fs;
	num_pix_ss = asic_size_ss * num_asics_ss;
	num_pix_tot = num_pix_fs * num_pix_ss;

	copy = (float *)calloc(num_pix_tot, sizeof(float));
	if ( copy == NULL ) {
		return 1;
	}

	memcpy(copy, data, num_pix_tot*sizeof(float));

	for (i = 0; i < num_pix_tot; i++) {
		copy[i] *= mask[i];
	}

	pfinter = allocate_peakfinder_intern_data(num_pix_tot, max_pix_count);
	if ( pfinter == NULL ) {
		free(copy);
		return 1;
	}

	peak_count = 0;

	for ( aiss=0 ; aiss<num_asics_ss ; aiss++ ) {
		for ( aifs=0 ; aifs<num_asics_fs ; aifs++ ) {

			process_panel(asic_size_fs, asic_size_ss, num_pix_fs,
			              aiss, aifs, rthreshold, roffset,
			              &peak_count, copy, pfinter, r_map, mask,
			              npix, com_fs, com_ss, com_index, tot_i,
			              max_i, sigma, snr, min_pix_count,
			              max_pix_count, local_bg_radius, min_snr,
			              max_n_peaks);
		}
	}
	*num_found_peaks = peak_count;

	free_peakfinder_intern_data(pfinter);
	free(copy);

	return 0;
}


int peakfinder8(struct image *img, int max_n_peaks,
                float threshold, float min_snr,
                int min_pix_count, int max_pix_count,
                int local_bg_radius, int min_res,
                int max_res, int use_saturated)
{
	struct radius_maps *rmaps;
	struct peakfinder_mask *pfmask;
	struct peakfinder_panel_data *pfdata;
	struct radial_stats *rstats;
	struct peakfinder_peak_data *pkdata;
	int num_rad_bins;
	int pi;
	int i, it_counter;
	int num_found_peaks;
	int remaining_max_num_peaks;
	int iterations;

	iterations = 5;

	if ( img-> det == NULL) {
		return 1;
	}

	rmaps = compute_radius_maps(img->det);
	if ( rmaps == NULL ) {
		return 1;
	}

	pfmask = create_peakfinder_mask(img, rmaps, min_res, max_res);
	if ( pfmask == NULL ) {
		free_radius_maps(rmaps);
		return 1;
	}

	pfdata = allocate_panel_data(img->det->n_panels);
	if ( pfdata == NULL) {
		free_radius_maps(rmaps);
		free_peakfinder_mask(pfmask);
		return 1;
	}

	for ( pi=0 ; pi<img->det->n_panels ; pi++ ) {
		pfdata->panel_h[pi] = img->det->panels[pi].h;
		pfdata->panel_w[pi] = img->det->panels[pi].w;
		pfdata->panel_data[pi] = img->dp[pi];
		pfdata->num_panels = img->det->n_panels;
	}

	num_rad_bins = compute_num_radial_bins(img->det->n_panels,
	                                       pfdata->panel_w,
	                                       pfdata->panel_h,
	                                       rmaps->r_maps);

	rstats = allocate_radial_stats(num_rad_bins);
	if ( rstats == NULL ) {
		free_radius_maps(rmaps);
		free_peakfinder_mask(pfmask);
		free_panel_data(pfdata);
		return 1;
	}

	for ( i=0 ; i<rstats->n_rad_bins ; i++) {
		rstats->rthreshold[i] = 1e9;
	}

	for ( it_counter=0 ; it_counter<iterations ; it_counter++ ) {

		for ( i=0; i<num_rad_bins; i++ ) {
			rstats->roffset[i] = 0;
			rstats->rsigma[i] = 0;
			rstats->rcount[i] = 0;
		}

		for ( pi=0 ; pi<pfdata->num_panels ; pi++ ) {

			fill_radial_bins(pfdata->panel_data[pi],
			                 pfdata->panel_w[pi],
			                 pfdata->panel_h[pi],
			                 rmaps->r_maps[pi],
			                 pfmask->masks[pi],
			                 rstats->rthreshold,
			                 rstats->roffset,
			                 rstats->rsigma,
			                 rstats->rcount);

		}

		compute_radial_stats(rstats->rthreshold,
		                     rstats->roffset,
		                     rstats->rsigma,
		                     rstats->rcount,
		                     num_rad_bins,
		                     min_snr,
		                     threshold);

	}

	pkdata = allocate_peak_data(max_n_peaks);
	if ( pkdata == NULL ) {
		free_radius_maps(rmaps);
		free_peakfinder_mask(pfmask);
		free_panel_data(pfdata);
		free_radial_stats(rstats);
		return 1;
	}

	remaining_max_num_peaks = max_n_peaks;

	for ( pi=0 ; pi<img->det->n_panels ; pi++) {

		int peaks_to_add;
		int pki;
		int ret;

		num_found_peaks = 0;

		if ( img->det->panels[pi].no_index ) {
			continue;
		}

		ret = peakfinder8_base(rstats->roffset,
		                       rstats->rthreshold,
		                       pfdata->panel_data[pi],
		                       pfmask->masks[pi],
		                       rmaps->r_maps[pi],
		                       pfdata->panel_w[pi], 1,
		                       pfdata->panel_h[pi], 1,
		                       max_n_peaks,
		                       &num_found_peaks,
		                       pkdata->npix,
		                       pkdata->com_fs,
		                       pkdata->com_ss,
		                       pkdata->com_index,
		                       pkdata->tot_i,
		                       pkdata->max_i,
		                       pkdata->sigma,
		                       pkdata->snr,
		                       min_pix_count,
		                       max_pix_count,
		                       local_bg_radius,
		                       min_snr);

		if ( ret != 0 ) {
			free_radius_maps(rmaps);
			free_peakfinder_mask(pfmask);
			free_panel_data(pfdata);
			free_radial_stats(rstats);
			return 1;
		}

		peaks_to_add = num_found_peaks;

		if ( num_found_peaks > remaining_max_num_peaks ) {
			peaks_to_add = remaining_max_num_peaks;
		}

		remaining_max_num_peaks -= peaks_to_add;

		for ( pki=0 ; pki<peaks_to_add ; pki++ ) {

			struct panel *p;

			p = &img->det->panels[pi];

			img->num_peaks += 1;
			if ( pkdata->max_i[pki] > p->max_adu ) {
				img->num_saturated_peaks++;
				if ( !use_saturated ) {
					continue;
				}
			}

			image_add_feature(img->features,
			                  pkdata->com_fs[pki],
			                  pkdata->com_ss[pki],
			                  p,
			                  img,
			                  pkdata->tot_i[pki],
			                  NULL);
		}
	}

	free_radius_maps(rmaps);
	free_peakfinder_mask(pfmask);
	free_panel_data(pfdata);
	free_radial_stats(rstats);
	free_peak_data(pkdata);
	return 0;
}
