/*
 * peakfinder8.c
 *
 * The peakfinder8 algorithm
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017-2021 Thomas White <taw@physics.org>
 *   2017      Valerio Mariani <valerio.mariani@desy.de>
 *   2017      Anton Barty <anton.barty@desy.de>
 *   2017      Oleksandr Yefanov <oleksandr.yefanov@desy.de>
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

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <profile.h>

#include "peakfinder8.h"
#include "detgeom.h"
#include "image.h"


/** \file peakfinder8.h */

// CrystFEL-only block 1

struct radius_maps
{
	float **r_maps;
	int *n_pixels;
	int n_rmaps;
};


struct radial_stats_pixels
{
	int n_panels;
	int *n_pixels; // n_pixels[panel]
	int **pidx;	   // pixel_index[panel][0..n_pixels]
	int **radius;  // pixel_radius[panel][0..n_pixels]
};


struct peakfinder_mask
{
	char **masks;
	int n_masks;
};


struct peakfinder_panel_data
{
	float **panel_data;
	int *panel_h;
	int *panel_w;
	int num_panels;
};
// End of CrystFEL-only block 1


struct radial_stats
{
	float *roffset;
	float *rthreshold;
	float *lthreshold;
	float *rsigma;
	int *rcount;
	int n_rad_bins;
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


static struct radial_stats_pixels *compute_rstats_pixels(struct radius_maps *rmaps)
{
	int p;
	int i;

	struct radial_stats_pixels *rsp = NULL;
	rsp = (struct radial_stats_pixels *)cfmalloc(sizeof(struct radial_stats_pixels));
	if ( rsp == NULL ) {
		return NULL;
	}
	rsp->n_pixels = (int *)cfmalloc(rmaps->n_rmaps * sizeof(int));
	if ( rsp->n_pixels == NULL ) {
		cffree(rsp);
		return NULL;
	}
	rsp->pidx = (int **)cfmalloc(rmaps->n_rmaps * sizeof(int *));
	if ( rsp->pidx == NULL ) {
		cffree(rsp->n_pixels);
		cffree(rsp);
		return NULL;
	}
	rsp->radius = (int **)cfmalloc(rmaps->n_rmaps * sizeof(int *));
	if ( rsp->radius == NULL ) {
		cffree(rsp->n_pixels);
		cffree(rsp->pidx);
		cffree(rsp);
		return NULL;
	}
	srand(0);

	int n_pixels_per_bin = 100; // Can make this a parameter

	// Assuming 5000 is the maximum possible radius
	int n_bins = 5000;
	int *n_pixels = (int *)cfmalloc(n_bins * sizeof(int)); // selected pixels per bin
	int *n_tot_pixels = (int *)cfmalloc(n_bins * sizeof(int));; // total pixels per bin
	int **panel = (int **)cfmalloc(n_bins * sizeof(int *)); // panel ID of selected pixels
	int **idx = (int **)cfmalloc(n_bins * sizeof(int *)); // index of selected pixels

	for ( i = 0; i < n_bins; i++ ) {
		n_pixels[i] = 0;
		n_tot_pixels[i] = 0;
		panel[i] = (int *)cfmalloc(n_pixels_per_bin * sizeof(int));
		idx[i] = (int *)cfmalloc(n_pixels_per_bin * sizeof(int));
	}
	int radius;

	for ( p = 0; p < rmaps->n_rmaps; p++ ) {
		rsp->n_pixels[p] = 0;
		for ( i = 0; i < rmaps->n_pixels[p]; i++ ) {
			// Reservoir sampling:
			radius = (int)rint(rmaps->r_maps[p][i]);
			n_tot_pixels[radius] += 1;

			if ( n_pixels[radius] < n_pixels_per_bin ) {
				panel[radius][n_pixels[radius]] = p;
				idx[radius][n_pixels[radius]] = i;

				n_pixels[radius] += 1;
				rsp->n_pixels[p] += 1;
			} else {
				int rand_i = rand() % n_tot_pixels[radius];
				if ( rand_i < n_pixels_per_bin ) {
					rsp->n_pixels[panel[radius][rand_i]] -= 1;
					rsp->n_pixels[p] += 1;

					panel[radius][rand_i] = p;
					idx[radius][rand_i] = i;
				}
			}
		}
	}

	int *sidx = (int *)cfmalloc(rmaps->n_rmaps * sizeof(int));
	if ( sidx == NULL ) {
		cffree(rsp->n_pixels);
		cffree(rsp->pidx);
		cffree(rsp->radius);
		cffree(rsp);
		return NULL;
	}
	for ( p = 0; p < rmaps->n_rmaps; p++ ) {
		rsp->pidx[p] = (int *)cfmalloc(rsp->n_pixels[p] * sizeof(int));
		if ( rsp->pidx[p] == NULL ) {
			for ( i = 0; i < p; i++ ) {
				cffree(rsp->pidx[i]);
				cffree(rsp->radius[i]);
			}
			cffree(rsp->pidx);
			cffree(rsp->radius);
			cffree(rsp->n_pixels);
			cffree(rsp);
			cffree(sidx);
			return NULL;
		}
		rsp->radius[p] = (int *)cfmalloc(rsp->n_pixels[p] * sizeof(int));
		if ( rsp->radius[p] == NULL ) {
			for ( i = 0; i < p; i++ ) {
				cffree(rsp->pidx[i]);
				cffree(rsp->radius[i]);
			}
			cffree(rsp->pidx[p]);
			cffree(rsp->pidx);
			cffree(rsp->radius);
			cffree(rsp->n_pixels);
			cffree(rsp);
			cffree(sidx);
			return NULL;
		}
		sidx[p] = 0;
	}

	for ( radius = 0; radius < n_bins; radius++ ) {
		for ( i = 0; i < n_pixels[radius]; i++ ) {
			p = panel[radius][i];
			rsp->pidx[p][sidx[p]] = idx[radius][i];
			rsp->radius[p][sidx[p]] = radius;
			sidx[p] += 1;
		}
	}
	cffree(sidx);
	for ( i = 0; i < n_bins; i++ ) {
		cffree(panel[i]);
		cffree(idx[i]);
	}
	cffree(panel);
	cffree(idx);
	cffree(n_pixels);
	cffree(n_tot_pixels);

	rsp->n_panels = rmaps->n_rmaps;
	return rsp;
}

static void free_rstats_pixels(struct radial_stats_pixels *rsp)
{
	int i;
	for ( i = 0; i < rsp->n_panels; i++ ) {
		cffree(rsp->pidx[i]);
		cffree(rsp->radius[i]);
	}
	cffree(rsp->pidx);
	cffree(rsp->radius);
	cffree(rsp->n_pixels);
	cffree(rsp);
}


static struct radius_maps *compute_radius_maps(struct detgeom *det)
{
	int i, u, iss, ifs;
	struct detgeom_panel p;
	struct radius_maps *rm = NULL;

	rm = (struct radius_maps *)cfmalloc(sizeof(struct radius_maps));
	if ( rm == NULL ) {
		return NULL;
	}

	rm->r_maps = (float **)cfmalloc(det->n_panels*sizeof(float*));
	if ( rm->r_maps == NULL ) {
		cffree(rm);
		return NULL;
	}

	rm->n_pixels = (int *)cfmalloc(det->n_panels*sizeof(int*));
	if ( rm->r_maps == NULL ) {
		cffree(rm);
		return NULL;
	}

	rm->n_rmaps = det->n_panels;

	for( i=0 ; i<det->n_panels ; i++ ) {

		p = det->panels[i];
		rm->r_maps[i] = (float *)cfmalloc(p.h*p.w*sizeof(float));

		if ( rm->r_maps[i] == NULL ) {
			for ( u = 0; u<i; u++ ) {
				cffree(rm->r_maps[u]);
			}
			cffree(rm);
			return NULL;
		}
		rm->n_pixels[i] = p.h * p.w;
		for ( iss=0 ; iss<p.h ; iss++ ) {
			for ( ifs=0; ifs<p.w; ifs++ ) {

				int rmi;
				int x,y;

				rmi = ifs + p.w * iss;

				x = (p.cnx  + ifs * p.fsx + iss * p.ssx);
				y = (p.cny  + ifs * p.fsy + iss * p.ssy);

				rm->r_maps[i][rmi] = sqrt(x * x + y * y);
			}
		}
	}
	return rm;
}


static void free_radius_maps(struct radius_maps *r_maps)
{
	int i;

	for ( i=0 ; i<r_maps->n_rmaps ; i++ ) {
		cffree(r_maps->r_maps[i]);
	}
	cffree(r_maps->r_maps);
	cffree(r_maps->n_pixels);
	cffree(r_maps);
}


// CrystFEL-only block 2
struct pf8_private_data *prepare_peakfinder8(struct detgeom *det, int fast_mode)
{
	struct pf8_private_data *data = NULL;
	if ( det == NULL ) {
		return NULL;
	}

	data = (struct pf8_private_data *)cfmalloc(sizeof(struct pf8_private_data));
	if ( data == NULL ) {
		return NULL;
	}
	data->rmaps = compute_radius_maps(det);
	if ( data->rmaps == NULL ) {
		cffree(data);
		return NULL;
	}
	if ( fast_mode ) {
		data->rpixels = compute_rstats_pixels(data->rmaps);
		if ( data->rpixels == NULL ) {
			free_radius_maps(data->rmaps);
			free(data);
			return NULL;
		}
	} else {
		data->rpixels = NULL;
	}
	data->fast_mode = fast_mode;
	return data;
}


void free_pf8_private_data(struct pf8_private_data *data)
{
	free_radius_maps(data->rmaps);
	if ( data->fast_mode ) {
		free_rstats_pixels(data->rpixels);
	}
	cffree(data);
}


static struct peakfinder_mask *create_peakfinder_mask(const struct image *img,
                                                      struct radius_maps *rmps,
                                                      int min_res,
                                                      int max_res)
{
	int i;
	struct peakfinder_mask *msk;

	msk = (struct peakfinder_mask *)cfmalloc(sizeof(struct peakfinder_mask));
	msk->masks =(char **) cfmalloc(img->detgeom->n_panels*sizeof(char*));
	msk->n_masks = img->detgeom->n_panels;
	for ( i=0; i<img->detgeom->n_panels; i++) {

		struct detgeom_panel p;
		int iss, ifs;

		p = img->detgeom->panels[i];

		msk->masks[i] = (char *)cfcalloc(p.w*p.h,sizeof(char));

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
		cffree(pfmask->masks[i]);
	}
	cffree(pfmask->masks);
	cffree(pfmask);
}


static struct peakfinder_panel_data *allocate_panel_data(int num_panels)
{

	struct peakfinder_panel_data *pfdata;

	pfdata = (struct peakfinder_panel_data *)cfmalloc(sizeof(struct peakfinder_panel_data));
	if ( pfdata == NULL ) {
		return NULL;
	}

	pfdata->panel_h = (int *)cfmalloc(num_panels*sizeof(int));
	if ( pfdata->panel_h == NULL ) {
		cffree(pfdata);
		return NULL;
	}

	pfdata->panel_w = (int *)cfmalloc(num_panels*sizeof(int));
	if ( pfdata->panel_w == NULL ) {
		cffree(pfdata->panel_h);
		cffree(pfdata);
		return NULL;
	}

	pfdata->panel_data = (float **)cfmalloc(num_panels*sizeof(float*));
	if ( pfdata->panel_data == NULL ) {
		cffree(pfdata->panel_w);
		cffree(pfdata->panel_h);
		cffree(pfdata);
		return NULL;
	}

	pfdata->num_panels = num_panels;

	return pfdata;
}


static void free_panel_data(struct peakfinder_panel_data *pfdata)
{
	cffree(pfdata->panel_data);
	cffree(pfdata->panel_w);
	cffree(pfdata->panel_h);
	cffree(pfdata);
}


static void compute_num_radial_bins(int w, int h, float *r_map, float *max_r)
{
	int ifs, iss;
	int pidx;

	for ( iss=0 ; iss<h ; iss++ ) {
		for ( ifs=0 ; ifs<w ; ifs++ ) {
			pidx = iss * w + ifs;
			if ( r_map[pidx] > *max_r ) {
				*max_r = r_map[pidx];
			}
		}
	}
}
// End of CrystFEL-only block 2


static struct radial_stats* allocate_radial_stats(int num_rad_bins)
{
	struct radial_stats* rstats;

	rstats = (struct radial_stats *)cfmalloc(sizeof(struct radial_stats));
	if ( rstats == NULL ) {
		return NULL;
	}

	rstats->roffset = (float *)cfmalloc(num_rad_bins*sizeof(float));
	if ( rstats->roffset == NULL ) {
		cffree(rstats);
		return NULL;
	}

	rstats->rthreshold = (float *)cfmalloc(num_rad_bins*sizeof(float));
	if ( rstats->rthreshold == NULL ) {
		cffree(rstats->roffset);
		cffree(rstats);
		return NULL;
	}

	rstats->lthreshold = (float *)cfmalloc(num_rad_bins*sizeof(float));
	if ( rstats->lthreshold == NULL ) {
		cffree(rstats->rthreshold);
		cffree(rstats->roffset);
		cffree(rstats);
		return NULL;
	}

	rstats->rsigma = (float *)cfmalloc(num_rad_bins*sizeof(float));
	if ( rstats->rsigma == NULL ) {
		cffree(rstats->roffset);
		cffree(rstats->rthreshold);
		cffree(rstats->lthreshold);
		cffree(rstats);
		return NULL;
	}

	rstats->rcount = (int *)cfmalloc(num_rad_bins*sizeof(int));
	if ( rstats->rcount == NULL ) {
		cffree(rstats->roffset);
		cffree(rstats->rthreshold);
		cffree(rstats->lthreshold);
		cffree(rstats->rsigma);
		cffree(rstats);
		return NULL;
	}

	rstats->n_rad_bins = num_rad_bins;

	return rstats;
}


static void free_radial_stats(struct radial_stats *rstats)
{
	cffree(rstats->roffset);
	cffree(rstats->rthreshold);
	cffree(rstats->lthreshold);
	cffree(rstats->rsigma);
	cffree(rstats->rcount);
	cffree(rstats);
}


static void fill_radial_bins(float *data,
                             int w,
                             int h,
                             float *r_map,
                             char *mask,
                             float *rthreshold,
                             float *lthreshold,
                             float *roffset,
                             float *rsigma,
                             int *rcount)
{
	int iss, ifs;
	int pidx;

	int curr_r;
	float value;

	for ( iss=0; iss<h; iss++ ) {
		for ( ifs=0; ifs<w; ifs++ ) {
			pidx = iss * w + ifs;
			if ( mask[pidx] != 0 ) {
				curr_r = (int)rint(r_map[pidx]);
				value = data[pidx];
				if ( value < rthreshold[curr_r]
				  && value > lthreshold[curr_r] )
				{
					roffset[curr_r] += value;
					rsigma[curr_r] += (value * value);
					rcount[curr_r] += 1;
				}
			}
		}
	}
}

static void fill_radial_bins_fast(float *data, int w, int h, int n_pixels,
				  int *pidx, int *radius, char *mask,
				  float *rthreshold, float *lthreshold,
				  float *roffset, float *rsigma, int *rcount)
{
	int i;

	int curr_r;
	float value;

	for (i = 0; i < n_pixels; i++)
	{
		if (mask[pidx[i]] != 0)
		{
			curr_r = radius[i];
			value = data[pidx[i]];
			if (value < rthreshold[curr_r] && value > lthreshold[curr_r])
			{
				roffset[curr_r] += value;
				rsigma[curr_r] += (value * value);
				rcount[curr_r] += 1;
			}
		}
	}
}

static void compute_radial_stats(float *rthreshold,
                                 float *lthreshold,
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
			rthreshold[ri] = FLT_MAX;
			lthreshold[ri] = FLT_MIN;
		} else {
			this_offset = roffset[ri] / rcount[ri];
			this_sigma = rsigma[ri] / rcount[ri] - (this_offset * this_offset);
			if ( this_sigma >= 0 ) {
				this_sigma = sqrt(this_sigma);
			}

			roffset[ri] = this_offset;
			rsigma[ri] = this_sigma;
			rthreshold[ri] = roffset[ri] + min_snr*rsigma[ri];
			lthreshold[ri] = roffset[ri] - min_snr*rsigma[ri];

			if ( rthreshold[ri] < acd_threshold ) {
				rthreshold[ri] = acd_threshold;
			}
		}
	}

}


struct peakfinder_peak_data *allocate_peak_data(int max_num_peaks)
{
	struct peakfinder_peak_data *pkdata;

	pkdata = (struct peakfinder_peak_data*)cfmalloc(sizeof(struct peakfinder_peak_data));
	if ( pkdata == NULL ) {
		return NULL;
	}

	pkdata->npix = (int *)cfmalloc(max_num_peaks*sizeof(int));
	if ( pkdata->npix == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata);
		return NULL;
	}

	pkdata->com_fs = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->com_fs == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata);
		return NULL;
	}

	pkdata->com_ss = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->com_ss == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata);
		return NULL;
	}

	pkdata->com_index = (int *)cfmalloc(max_num_peaks*sizeof(int));
	if ( pkdata->com_ss == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata->com_ss);
		cffree(pkdata);
		return NULL;
	}

	pkdata->tot_i = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->tot_i == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata->com_ss);
		cffree(pkdata->com_index);
		cffree(pkdata);
		return NULL;
	}

	pkdata->max_i = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->max_i == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata->com_ss);
		cffree(pkdata->com_index);
		cffree(pkdata->tot_i);
		cffree(pkdata);
		return NULL;
	}

	pkdata->sigma = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->sigma == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata->com_ss);
		cffree(pkdata->com_index);
		cffree(pkdata->tot_i);
		cffree(pkdata->max_i);
		cffree(pkdata);
		return NULL;
	}

	pkdata->snr = (float *)cfmalloc(max_num_peaks*sizeof(float));
	if ( pkdata->snr == NULL ) {
		cffree(pkdata->npix);
		cffree(pkdata->com_fs);
		cffree(pkdata->com_ss);
		cffree(pkdata->com_index);
		cffree(pkdata->tot_i);
		cffree(pkdata->max_i);
		cffree(pkdata->sigma);
		cffree(pkdata);
		return NULL;
	}

	return pkdata;
}


static void free_peak_data(struct peakfinder_peak_data *pkdata) {
	cffree(pkdata->npix);
	cffree(pkdata->com_fs);
	cffree(pkdata->com_ss);
	cffree(pkdata->com_index);
	cffree(pkdata->tot_i);
	cffree(pkdata->max_i);
	cffree(pkdata->sigma);
	cffree(pkdata->snr);
	cffree(pkdata);
}


static struct peakfinder_intern_data *allocate_peakfinder_intern_data(int data_size,
                                                                      int max_pix_count)
{

	struct peakfinder_intern_data *intern_data;

	intern_data = (struct peakfinder_intern_data *)cfmalloc(sizeof(struct peakfinder_intern_data));
	if ( intern_data == NULL ) {
		return NULL;
	}

	intern_data->pix_in_peak_map =(char *)cfcalloc(data_size, sizeof(char));
	if ( intern_data->pix_in_peak_map == NULL ) {
		cffree(intern_data);
		return NULL;
	}

	intern_data->infs =(int *)cfcalloc(data_size, sizeof(int));
	if ( intern_data->infs == NULL ) {
		cffree(intern_data->pix_in_peak_map);
		cffree(intern_data);
		return NULL;
	}

	intern_data->inss =(int *)cfcalloc(data_size, sizeof(int));
	if ( intern_data->inss == NULL ) {
		cffree(intern_data->pix_in_peak_map);
		cffree(intern_data->infs);
		cffree(intern_data);
		return NULL;
	}

	intern_data->peak_pixels =(int *)cfcalloc(max_pix_count, sizeof(int));
	if ( intern_data->peak_pixels == NULL ) {
		cffree(intern_data->pix_in_peak_map);
		cffree(intern_data->infs);
		cffree(intern_data->inss);
		cffree(intern_data);
		return NULL;
	}

	return intern_data;
}


static void free_peakfinder_intern_data(struct peakfinder_intern_data *pfid)
{
	cffree(pfid->peak_pixels);
	cffree(pfid->pix_in_peak_map);
	cffree(pfid->infs);
	cffree(pfid->inss);
	cffree(pfid);
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

	// Loop through search pattern
	for ( k=0; k<search_n; k++ ) {

		if ( (pfinter->infs[p] + search_fs[k]) < 0 ) continue;
		if ( (pfinter->infs[p] + search_fs[k]) >= asic_size_fs ) continue;
		if ( (pfinter->inss[p] + search_ss[k]) < 0 ) continue;
		if ( (pfinter->inss[p] + search_ss[k]) >= asic_size_ss ) continue;

		// Neighbour point in big array
		curr_fs = pfinter->infs[p] + search_fs[k] + aifs * asic_size_fs;
		curr_ss = pfinter->inss[p] + search_ss[k] + aiss * asic_size_ss;
		pi = curr_fs + curr_ss * num_pix_fs;

		curr_radius = (int)rint(r_map[pi]);
		curr_threshold = rthreshold[curr_radius];

		// Above threshold?
		if ( copy[pi] > curr_threshold
		  && pfinter->pix_in_peak_map[pi] == 0
		  && mask[pi] != 0 ) {

			curr_i = copy[pi] - roffset[curr_radius];
			*sum_i += curr_i;
			*sum_com_fs += curr_i * ((float)curr_fs);  // for center of mass x
			*sum_com_ss += curr_i * ((float)curr_ss);  // for center of mass y

			pfinter->inss[*num_pix_in_peak] = pfinter->inss[p] + search_ss[k];
			pfinter->infs[*num_pix_in_peak] = pfinter->infs[p] + search_fs[k];
			pfinter->pix_in_peak_map[pi] = 1;
			if ( *num_pix_in_peak < max_pix_count ) {
				  pfinter->peak_pixels[*num_pix_in_peak] = pi;
			}
			*num_pix_in_peak = *num_pix_in_peak + 1;
		}
	}
}


static void search_in_ring(int ring_width, int com_fs_int, int com_ss_int,
                           float *copy, float *r_map,
                           float *rthreshold, float *roffset,
                           char *pix_in_peak_map, char *mask, int asic_size_fs,
                           int asic_size_ss, int aifs, int aiss,
                           int num_pix_fs,float *local_sigma, float *local_offset,
                           float *background_max_i, int com_idx,
                           int local_bg_radius)
{
	int ssj, fsi;
	float pix_radius;
	int curr_fs, curr_ss;
	int pi;
	int curr_radius;
	float curr_threshold;
	float curr_i;

	int np_sigma;
	int local_radius;

	float sum_i;
	float sum_i_squared;

	ring_width = 2 * local_bg_radius;

	sum_i = 0;
	sum_i_squared = 0;
	np_sigma = 0;
	local_radius = 0;

	for ( ssj = -ring_width ; ssj<ring_width ; ssj++ ) {
		for ( fsi = -ring_width ; fsi<ring_width ; fsi++ ) {

			// Within-ASIC check
			if ( (com_fs_int + fsi) < 0 ) continue;
			if ( (com_fs_int + fsi) >= asic_size_fs ) continue;
			if ( (com_ss_int + ssj) < 0 ) continue;
			if ( (com_ss_int + ssj) >= asic_size_ss )
			continue;

			// Within outer ring check
			pix_radius = sqrt(fsi * fsi + ssj * ssj);
			if ( pix_radius>ring_width ) continue;

			// Position of this point in data stream
			curr_fs = com_fs_int + fsi + aifs * asic_size_fs;
			curr_ss = com_ss_int + ssj + aiss * asic_size_ss;
			pi = curr_fs + curr_ss * num_pix_fs;

			curr_radius = (int)rint(r_map[pi]);
			curr_threshold = rthreshold[curr_radius];

			// Intensity above background ??? just intensity?
			curr_i = copy[pi];

			// Keep track of value and value-squared for offset and sigma calculation
			if ( curr_i < curr_threshold && pix_in_peak_map[pi] == 0 && mask[pi] != 0 ) {

				np_sigma++;
				sum_i += curr_i;
				sum_i_squared += (curr_i * curr_i);

				if ( curr_i > *background_max_i ) {
					*background_max_i = curr_i;
				}
			}
		}
	}

	// Calculate local background and standard deviation
	if ( np_sigma != 0 ) {
		*local_offset = sum_i / np_sigma;
		*local_sigma = sum_i_squared / np_sigma - (*local_offset * *local_offset);
		if (*local_sigma >= 0) {
			*local_sigma = sqrt(*local_sigma);
		} else {
			*local_sigma = 0.01;
		}
	} else {
		local_radius = (int)rint(r_map[(int)rint(com_idx)]);
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

	// Loop over pixels within a module
	for ( pxss=1 ; pxss<asic_size_ss-1 ; pxss++ ) {
		for ( pxfs=1 ; pxfs<asic_size_fs-1 ; pxfs++ ) {

			float curr_thresh;
			int pxidx;
			int curr_rad;

			pxidx = (pxss + aiss * asic_size_ss) * num_pix_fs +
			pxfs + aifs * asic_size_fs;

			curr_rad = (int)rint(r_map[pxidx]);
			curr_thresh = rthreshold[curr_rad];

			if ( copy[pxidx] > curr_thresh
			  && pfinter->pix_in_peak_map[pxidx] == 0
			  && mask[pxidx] != 0 ) {   //??? not sure if needed

				// This might be the start of a new peak - start searching
				float sum_com_fs, sum_com_ss;
				float sum_i;
				float peak_com_fs, peak_com_ss;
				float peak_com_fs_int, peak_com_ss_int;
				float peak_tot_i, pk_tot_i_raw;
				float peak_max_i, pk_max_i_raw;
				float peak_snr;
				float local_sigma, local_offset;
				float background_max_i;
				int lt_num_pix_in_pk;
				int ring_width;
				int peak_idx;
				int com_idx;
				int p;

				pfinter->infs[0] = pxfs;
				pfinter->inss[0] = pxss;
				pfinter->peak_pixels[0] = pxidx;
				num_pix_in_peak = 0; //y 1;

				sum_i = 0;
				sum_com_fs = 0;
				sum_com_ss = 0;

				// Keep looping until the pixel count within this peak does not change
				do {
					lt_num_pix_in_pk = num_pix_in_peak;

					// Loop through points known to be within this peak
					for ( p=0; p<=num_pix_in_peak; p++ ) { //changed from 1 to 0 by O.Y.
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

				// Too many or too few pixels means ignore this 'peak'; move on now
				if ( num_pix_in_peak < min_pix_count || num_pix_in_peak > max_pix_count ) continue;

				// If for some reason sum_i is 0 - it's better to skip
				if ( fabs(sum_i) < 1e-10 ) continue;

				// Calculate center of mass for this peak from initial peak search
				peak_com_fs = sum_com_fs / fabs(sum_i);
				peak_com_ss = sum_com_ss / fabs(sum_i);

				com_idx = (int)rint(peak_com_fs) + (int)rint(peak_com_ss) * num_pix_fs;

				peak_com_fs_int = (int)rint(peak_com_fs) - aifs * asic_size_fs;
				peak_com_ss_int = (int)rint(peak_com_ss) - aiss * asic_size_ss;

				// Calculate the local signal-to-noise ratio and local background in an annulus around
				// this peak (excluding pixels which look like they might be part of another peak)
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

				// Re-integrate (and re-centroid) peak using local background estimates
				peak_tot_i = 0;
				pk_tot_i_raw = 0;
				peak_max_i = 0;
				pk_max_i_raw = 0;
				sum_com_fs = 0;
				sum_com_ss = 0;

				for ( peak_idx = 0 ;
					peak_idx < num_pix_in_peak && peak_idx < max_pix_count ;
					peak_idx++ ) {

					int curr_idx;
					float curr_i;
					float curr_i_raw;
					int curr_fs, curr_ss;

					curr_idx = pfinter->peak_pixels[peak_idx];
					curr_i_raw = copy[curr_idx];
					curr_i = curr_i_raw - local_offset;
					peak_tot_i += curr_i;
					pk_tot_i_raw += curr_i_raw;

					// Remember that curr_idx = curr_fs + curr_ss*num_pix_fs
					curr_fs = curr_idx % num_pix_fs;
					curr_ss = curr_idx / num_pix_fs;
					sum_com_fs += curr_i_raw * ((float)curr_fs);
					sum_com_ss += curr_i_raw * ((float)curr_ss);

					if ( curr_i_raw > pk_max_i_raw ) pk_max_i_raw = curr_i_raw;
					if ( curr_i > peak_max_i ) peak_max_i = curr_i;
				}


				// This CAN happen! Better to skip...
				if ( fabs(pk_tot_i_raw) < 1e-10 ) continue;

				peak_com_fs = sum_com_fs / fabs(pk_tot_i_raw);
				peak_com_ss = sum_com_ss / fabs(pk_tot_i_raw);

				// Calculate signal-to-noise and apply SNR criteria
				if ( fabs(local_sigma) > 1e-10 ) {
					peak_snr = peak_tot_i / local_sigma;
				} else {
					peak_snr = 0;
				}

				if (peak_snr < min_snr) continue;

				// Is the maximum intensity in the peak enough above intensity in background region to
				// be a peak and not noise? The more pixels there are in the peak, the more relaxed we
				// are about this criterion
				//f_background_thresh = background_max_i - local_offset; //!!! Ofiget'!  If I uncomment
				// if (peak_max_i < f_background_thresh) {               // these lines the result is
				// different!
				if (peak_max_i < background_max_i - local_offset) continue;

				if ( peak_com_fs < aifs*asic_size_fs
				  || peak_com_fs > (aifs+1)*asic_size_fs-1
				  || peak_com_ss < aiss*asic_size_ss
				  || peak_com_ss > (aiss+1)*asic_size_ss-1)
				{
					continue;
				}

				// This is a peak? If so, add info to peak list
				if ( num_pix_in_peak >= min_pix_count
				  && num_pix_in_peak <= max_pix_count ) {

					// Bragg peaks in the mask
					for ( peak_idx = 0 ;
					      peak_idx < num_pix_in_peak &&
					      peak_idx < max_pix_count ;
					      peak_idx++ ) {
						pfinter->pix_in_peak_map[pfinter->peak_pixels[peak_idx]] = 2;
					}

					int peak_com_idx;
					peak_com_idx = (int)rint(peak_com_fs) + (int)rint(peak_com_ss) *
						                num_pix_fs;
					// Remember peak information
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
					}
					*peak_count += 1;
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
                            int local_bg_radius, float min_snr,
                            char* outliersMask)
{

	int num_pix_fs, num_pix_ss, num_pix_tot;
	int aifs, aiss;
	int peak_count;
	struct peakfinder_intern_data *pfinter;

	num_pix_fs = asic_size_fs * num_asics_fs;
	num_pix_ss = asic_size_ss * num_asics_ss;
	num_pix_tot = num_pix_fs * num_pix_ss;

	pfinter = allocate_peakfinder_intern_data(num_pix_tot, max_pix_count);
	if ( pfinter == NULL ) {
		return 1;
	}

	peak_count = 0;

	// Loop over modules (nxn array)
	for ( aiss=0 ; aiss<num_asics_ss ; aiss++ ) {
		for ( aifs=0 ; aifs<num_asics_fs ; aifs++ ) {                 // ??? to change to proper panels need
			process_panel(asic_size_fs, asic_size_ss, num_pix_fs, // change copy, mask, r_map
			              aiss, aifs, rthreshold, roffset,
			              &peak_count, data, pfinter, r_map, mask,
			              npix, com_fs, com_ss, com_index, tot_i,
			              max_i, sigma, snr, min_pix_count,
			              max_pix_count, local_bg_radius, min_snr,
			              max_n_peaks);
		}
	}
	*num_found_peaks = peak_count;

	if (outliersMask != NULL) {
		memcpy(outliersMask, pfinter->pix_in_peak_map, num_pix_tot*sizeof(char));
	}

	free_peakfinder_intern_data(pfinter);

	return 0;
}


/**
 * \param img An \ref image structure
 * \param max_n_peaks The maximum number of peaks to be searched for
 * \param threshold The image threshold value, in detector units
 * \param min_snr The minimum signal to noise ratio for a peak
 * \param min_pix_count The minimum number of pixels in a peak
 * \param max_pix_count The maximum number of pixels in a peak
 * \param local_bg_radius The averaging radius for background calculation
 * \param min_res The minimum number of pixels out from the center
 * \param max_res The maximum number of pixels out from the center
 * \param use_saturated Whether saturated peaks should be considered
 *
 * Runs the peakfinder8 peak search algorithm, and returns an \ref ImageFeatureList,
 * or NULL on error.
 */
ImageFeatureList *peakfinder8(const struct image *img, int max_n_peaks,
                              float threshold, float min_snr,
                              int min_pix_count, int max_pix_count,
                              int local_bg_radius, int min_res,
                              int max_res, int use_saturated,
                              int fast_mode, struct pf8_private_data *private_data)
{
	struct pf8_private_data *geomdata;
	struct radius_maps *rmaps;
	struct radial_stats_pixels *rspixels;

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
	float max_r;
	ImageFeatureList *peaks;

	iterations = 5;

	if ( img->detgeom == NULL) return NULL;

	profile_start("pf8-rmaps");
	if ( private_data == NULL ) {
		geomdata = prepare_peakfinder8(img->detgeom, fast_mode);
	} else {
		geomdata = private_data;
	}
	rmaps = geomdata->rmaps;
	rspixels = geomdata->rpixels;
	profile_end("pf8-rmaps");
	if (geomdata == NULL) return NULL;

	profile_start("pf8-mask");
	pfmask = create_peakfinder_mask(img, rmaps, min_res, max_res);
	profile_end("pf8-mask");
	if ( pfmask == NULL ) {
		if ( private_data == NULL ) free_pf8_private_data(geomdata);
		return NULL;
	}

	pfdata = allocate_panel_data(img->detgeom->n_panels);
	if ( pfdata == NULL) {
		if ( private_data == NULL ) free_pf8_private_data(geomdata);
		free_peakfinder_mask(pfmask);
		return NULL;
	}

	for ( pi=0 ; pi<img->detgeom->n_panels ; pi++ ) {
		pfdata->panel_h[pi] = img->detgeom->panels[pi].h;
		pfdata->panel_w[pi] = img->detgeom->panels[pi].w;
		pfdata->panel_data[pi] = img->dp[pi];
		pfdata->num_panels = img->detgeom->n_panels;
	}

	max_r = -1e9;

	for ( pi=0 ; pi<pfdata->num_panels ; pi++ ) {

		compute_num_radial_bins(pfdata->panel_w[pi],
		                        pfdata->panel_h[pi],
		                        rmaps->r_maps[pi],
		                        &max_r);
	}

	num_rad_bins = (int)ceil(max_r) + 1;

	rstats = allocate_radial_stats(num_rad_bins);
	if ( rstats == NULL ) {
		if ( private_data == NULL ) free_pf8_private_data(geomdata);
		free_peakfinder_mask(pfmask);
		free_panel_data(pfdata);
		return NULL;
	}

	for ( i=0 ; i<rstats->n_rad_bins ; i++) {
		rstats->rthreshold[i] = 1e9;
		rstats->lthreshold[i] = -1e9;
	}
	profile_start("pf8-rstats");
	for ( it_counter=0 ; it_counter<iterations ; it_counter++ ) {

		for ( i=0; i<num_rad_bins; i++ ) {
			rstats->roffset[i] = 0;
			rstats->rsigma[i] = 0;
			rstats->rcount[i] = 0;
		}

		for ( pi=0 ; pi<pfdata->num_panels ; pi++ ) {
			if ( fast_mode ) {
				fill_radial_bins_fast(pfdata->panel_data[pi],
						      pfdata->panel_w[pi],
						      pfdata->panel_h[pi],
						      rspixels->n_pixels[pi],
						      rspixels->pidx[pi],
						      rspixels->radius[pi],
						      pfmask->masks[pi],
						      rstats->rthreshold,
						      rstats->lthreshold,
						      rstats->roffset,
						      rstats->rsigma,
						      rstats->rcount);
			} else {
				fill_radial_bins(pfdata->panel_data[pi],
						 pfdata->panel_w[pi],
						 pfdata->panel_h[pi],
						 rmaps->r_maps[pi],
						 pfmask->masks[pi],
						 rstats->rthreshold,
						 rstats->lthreshold,
						 rstats->roffset,
						 rstats->rsigma,
						 rstats->rcount);
			}
		}

		compute_radial_stats(rstats->rthreshold,
		                     rstats->lthreshold,
		                     rstats->roffset,
		                     rstats->rsigma,
		                     rstats->rcount,
		                     num_rad_bins,
		                     min_snr,
		                     threshold);

	}
	profile_end("pf8-rstats");

	pkdata = allocate_peak_data(max_n_peaks);
	if ( pkdata == NULL ) {
		if ( private_data == NULL ) free_pf8_private_data(geomdata);
		free_peakfinder_mask(pfmask);
		free_panel_data(pfdata);
		free_radial_stats(rstats);
		return NULL;
	}

	remaining_max_num_peaks = max_n_peaks;
	peaks = image_feature_list_new();
	profile_start("pf8-search");
	for ( pi=0 ; pi<img->detgeom->n_panels ; pi++) {

		int peaks_to_add;
		int pki;
		int ret;

		num_found_peaks = 0;

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
		                       min_snr,
		                       NULL);

		if ( ret != 0 ) {
			if ( private_data == NULL ) free_pf8_private_data(geomdata);
			free_peakfinder_mask(pfmask);
			free_panel_data(pfdata);
			free_radial_stats(rstats);
			image_feature_list_free(peaks);
			profile_end("pf8-search");
			return NULL;
		}

		peaks_to_add = num_found_peaks;

		if ( num_found_peaks > remaining_max_num_peaks ) {
			peaks_to_add = remaining_max_num_peaks;
		}

		remaining_max_num_peaks -= peaks_to_add;

		for ( pki=0 ; pki<peaks_to_add ; pki++ ) {

			struct detgeom_panel *p;

			p = &img->detgeom->panels[pi];

			if ( pkdata->max_i[pki] > p->max_adu ) {
				if ( !use_saturated ) {
					continue;
				}
			}

			image_add_feature(peaks,
			                  pkdata->com_fs[pki]+0.5,
			                  pkdata->com_ss[pki]+0.5,
			                  pi, pkdata->tot_i[pki], NULL);
		}
	}
	profile_end("pf8-search");

	if ( private_data == NULL ) free_pf8_private_data(geomdata);
	free_peakfinder_mask(pfmask);
	free_panel_data(pfdata);
	free_radial_stats(rstats);
	free_peak_data(pkdata);
	return peaks;
}
