/*
 * powder_plot.c
 *
 * Plot powder patterns
 *
 * (c) 2011 Andrew Aquila <andrew.aquila@cfel.de>
 * (c) 2006-2011 Thomas White  <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "stream.h"
#include "reflist.h"
#include "utils.h"
#include "image.h"
#include "detector.h"
#include "index.h"
#include "hdf5-file.h"
#include "beam-parameters.h"
#include "reflist-utils.h"
#include "symmetry.h"


struct bin_stats {
	unsigned int N;
	double total;
	double mean;
	double std_dev;
	double q_min;
	double q_max;
	double q_value;
	double fract;
};

struct histogram_info {
	double q_max;
	double q_min;
	double q_delta;
	unsigned int histsize;
	int spacing; //linear, q^2, & equal volume
};

enum {
	PLOT_PEAKS,
	PLOT_HKL,
	PLOT_REFL,
	PLOT_H5,
	PLOT_D
};

enum {
	FILE_STREAM,
	FILE_HKL,
	FILE_H5
};

enum {
	LINEAR,
	q2,
	VOLUME
};


static int find_q_bin_index(double q, struct histogram_info *info,
                            struct bin_stats *hist)
{
	/* bisection search alg. find q_bin index of order Log(n) time */
	int mid;
	int min = 0;
	int max = info->histsize-1;
	if (q < hist[min].q_max) {return min;}
	if (q > hist[max].q_min) {return max;}
	do {
		mid = (min + max) / 2;
		if (q < hist[mid].q_min) {
			max = mid;
		} else if (q > hist[mid].q_max){
			min = mid ;
		} else {
			return mid;
		}
	} while(max - min > 1);
	return mid;
}


/* Used for HDF5 files, peak list and stream positions */
static int add_peak_to_histogram(double fs, double ss, double intensity,
                                 struct image *image,
                                 struct histogram_info *info,
                                 struct bin_stats *hist)
{
	struct rvec r;
	double q, delta;
	int i;

	r = get_q(image, fs, ss, NULL, 1.0/ image->lambda);
	q = modulus(r.u, r.v, r.w);

	/* Ignore q value if outside of range */
	if ( (q<info->q_min) || (q>info->q_max) ) {
		return 1;
	}
	i = find_q_bin_index(q, info, hist);

	/* See Knuth TAOCP vol 2, 3rd ed, pg 232 for running variance */
	delta = intensity - hist[i].mean;
	hist[i].N++;
	hist[i].total += intensity;
	hist[i].mean = hist[i].mean + delta /hist[i].N;
	hist[i].std_dev = hist[i].std_dev + (delta *(intensity - hist[i].mean));

	return 0;
}


/* Used for d and hkl of stream files where redundancy = 1 */
static int add_d_to_histogram(double q, double intensity,
                              struct histogram_info *info,
                              struct bin_stats *hist)
{
	double delta;
	int i;

	/* Ignore q value if outside of range */
	if ( (q<info->q_min) || (q>info->q_max) ) {
		return 1;
	}
	i = find_q_bin_index(q, info, hist);

	delta = intensity - hist[i].mean;
	hist[i].N++;
	hist[i].total += intensity;
	hist[i].mean = hist[i].mean + delta /hist[i].N;
	hist[i].std_dev = hist[i].std_dev + (delta *(intensity - hist[i].mean));

	return 0;
}


static int add_hkl_to_histogram(double q, double intensity, int redundancy,
                                int q_scaling, struct histogram_info *info,
                                struct bin_stats *hist)
{
	int i = 0;

	/* Ignore q value if outside of range */
	if ( (q<info->q_min) || (q>info->q_max) ) {
		return 1;
	}

	/* The accounting is the intensity of the reflection times the
	 * number of occurance of that reflection smeared out over the
	 * surface area which is 4*pi*q^2  the 4*pi is left out since it is a
	 * common constant and the total is in arbitrary units.
	 */
	for ( i=0; i<redundancy; i++ ) {
		if ( q_scaling ) {
			add_d_to_histogram(q, intensity/(q*q), info, hist);
		} else {
			add_d_to_histogram(q, intensity, info, hist);
		}
	}

	return 0;
}


static int histogram_setup(struct histogram_info *info,
                           struct bin_stats *histdata)
{
	int i;
	double x;

	if ( info->spacing == LINEAR ) {
		x = 1.0;
	} else if ( info->spacing == q2 ) {
		x = 2.0;
	} else {
		x = 3.0;
	}

	for ( i=0; i<info->histsize; i++ ) {

		double qd, qm;

		histdata[i].N  = 0;
		histdata[i].total = 0.0;
		histdata[i].mean = 0.0;
		histdata[i].std_dev = 0.0;
		histdata[i].fract = 0.0;

		qd = info->q_delta;
		qm = info->q_min;

		histdata[i].q_min = pow( (i*qd) + pow(qm, x), 1.0/x);

		histdata[i].q_max = pow( ((i+1.0)*qd) + pow(qm, x), 1.0/x);

		histdata[i].q_value= pow( ((i+0.5)*qd) + pow(qm, x), 1.0/x);

	}
	return 0;
}


static int ring_fraction_calc(struct histogram_info *info,
                              struct bin_stats *hist, struct image *image)
{
	int fs,ss;
	int bin;

	/* Check that detector geometry is present and wavelength is valid */
	if ( (image->det == NULL) || (image->lambda < 0.0) ) return 1;

	/* Loop over all pixels */
	for ( ss=0; ss<image->height; ss++ ) {
	for ( fs=0; fs<image->width;  fs++ ) {

		struct panel *p;
		struct rvec r;
		int i;
		double q, q_fs, q_ss;

		r = get_q(image, fs, ss, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		/* If pixel is valid (not a bad pixel and not out of range) */
		if ( (q>info->q_min) && (q<info->q_max) &&
		     (in_bad_region(image->det,fs,ss) == 0) ) {

			/* Select the panel, then (sometimes) ask for the q
			 * of the (corner of the) pixel one step beyond the
			 * edge, to get the exact size of the required pixel.
			 */
			p = find_panel(image->det, fs, ss);

			r = get_q_for_panel(p, (fs+1)-(double)p->min_fs,
			                       ss-(double)p->min_ss,
			                       NULL, 1.0/image->lambda);
			q_fs =  modulus(r.u, r.v, r.w);

			r = get_q_for_panel(p, fs-(double)p->min_fs,
			                       (ss+1) -(double)p->min_ss,
			                       NULL, 1.0/image->lambda);
			q_ss =  modulus(r.u, r.v, r.w);

			i = find_q_bin_index(q, info, hist);

			hist[i].fract = hist[i].fract + fabs((q_fs-q)*(q_ss-q));

		}

	}
	}

	/* Divide measured area by ring area */
	for ( bin=0; bin<info->histsize; bin++ ) {

		double inner_area, outer_area, ring_area;

		outer_area = pow(hist[bin].q_max, 2.0);
		inner_area = pow(hist[bin].q_min, 2.0);
		ring_area = M_PI*(outer_area - inner_area);

		hist[bin].fract = hist[bin].fract / ring_area;

	}

	return 0;
}


static unsigned int process_h5(struct image *image, struct histogram_info *info,
                               struct bin_stats *histdata)
{
	int fs, ss;
	double intensity;

	for ( ss=0; ss<image->height; ss++ ) {
	for ( fs=0; fs<image->width; fs++ ) {

		intensity = image->data[fs + image->width*ss];
		if ( !in_bad_region(image->det,fs,ss) ) {
			add_peak_to_histogram(fs, ss, intensity,
		                              image, info, histdata);
		}
	}
	progress_bar(ss, image->height, "Processing");
	}
	return 0;
}


static unsigned int process_hkl(struct image *image, const SymOpList *sym,
                                UnitCell *cell,
                                struct histogram_info *info,
                                struct bin_stats *histdata,
                                int q_scaling, int use_redundancy)
{
	Reflection *refl;
	RefListIterator *iter;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	int h, k, l, redundancy;
	double q, intensity;
	unsigned int nref;
	SymOpMask *m;

	m = new_symopmask(sym);

	nref = num_reflections(image->reflections);

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		get_indices(refl, &h, &k, &l);
		intensity = get_intensity(refl);
		if ( use_redundancy ) {
			redundancy = get_redundancy(refl);
		} else {
			special_position(sym, m, h, k, l);
			redundancy = num_equivs(sym, m);
		}

		/* Multiply by 2 to get 1/d (in m^-1) */
		q = 2.0 * resolution(cell, h, k, l);

		add_hkl_to_histogram(q, intensity, redundancy, q_scaling,
		                     info, histdata);

		n_peaks += redundancy;

		i++;
		progress_bar(i, nref, "Processing");

	}

	free_symopmask(m);

	return n_peaks;
}


static unsigned int process_stream_reflection(FILE *fh, struct image *image,
                                              struct histogram_info *info,
                                              struct bin_stats *histdata,
                                              unsigned int *n_patterns)
{
	int rval;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	Reflection *refl;
	RefListIterator *iter;
	double intensity, fs_double, ss_double;
	unsigned int processing_total;

	processing_total = count_patterns(fh);
	rewind(fh);

	do {
		/* Get data from next chunk */
		rval = read_chunk(fh, image);
		if ( rval ) continue;

		/* Check if the pattern indexed, if so use those peaks */
		if ( image->reflections != NULL ) {

			(*n_patterns)++;
			for ( refl = first_refl(image->reflections, &iter);
			      refl != NULL;
		      	      refl = next_refl(refl, iter) ) {
				/* note added fs_double as fs is an int */
				intensity = get_intensity(refl);
				get_detector_pos(refl, &fs_double, &ss_double);

				if ( !add_peak_to_histogram(fs_double,
				                            ss_double,
				                            intensity,
				                            image, info,
				                            histdata) )
				{
					n_peaks++;
				}

			}

		}

		free(image->filename);
		reflist_free(image->reflections);
		image_feature_list_free(image->features);
		cell_free(image->indexed_cell);

		i++;
		progress_bar(i, processing_total, "Processing");

	} while ( rval == 0 );

	return n_peaks;
}


static unsigned int process_stream_d(FILE *fh, struct image *image,
                                     struct histogram_info *info,
                                     struct bin_stats *histdata,
                                     unsigned int *n_patterns)
{
	int h, k, l, rval;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	Reflection *refl;
	RefListIterator *iter;
	double intensity, q;
	unsigned int processing_total;

	processing_total = count_patterns(fh);
	rewind(fh);

	do {

		/* Get data from next chunk */
		rval = read_chunk(fh, image);
		if ( rval ) continue;

		if ( image->reflections != NULL ) {

			(*n_patterns)++;

			for ( refl = first_refl(image->reflections, &iter);
			      refl != NULL;
			      refl = next_refl(refl, iter) )
			{
				get_indices(refl, &h, &k, &l);
				intensity = get_intensity(refl);
				q = 2.0 * resolution(image->indexed_cell,
				                     h, k, l);
				if ( !add_d_to_histogram(q, intensity, info,
				                        histdata)) n_peaks++;
			}
		}

		free(image->filename);
		reflist_free(image->reflections);
		image_feature_list_free(image->features);
		cell_free(image->indexed_cell);

		i++;
		progress_bar(i, processing_total, "Processing");

	} while ( rval == 0 );

	return n_peaks;
}


static unsigned int process_stream_hkl(FILE *fh, struct image *image,
                                       struct histogram_info *info,
                                       struct bin_stats *histdata,
                                       UnitCell *cell, unsigned int *n_patterns)
{
	int rval;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	Reflection *refl;
	RefListIterator *iter;
	double intensity, q;
	unsigned int processing_total;

	processing_total = count_patterns(fh);
	rewind(fh);

	do {

		/* Get data from next chunk */
		rval = read_chunk(fh, image);
		if ( rval ) continue;
		if ( image->reflections != NULL ) {

			(*n_patterns)++;

			for ( refl = first_refl(image->reflections, &iter);
		      	      refl != NULL;
			      refl = next_refl(refl, iter) )
			{
				int h, k, l;

				get_indices(refl, &h, &k, &l);
				intensity = get_intensity(refl);
				q = 2.0 * resolution(cell, h, k, l);

				if ( !add_d_to_histogram(q, intensity, info,
				                         histdata) ) n_peaks++;
			}
		}

		free(image->filename);
		reflist_free(image->reflections);
		image_feature_list_free(image->features);
		cell_free(image->indexed_cell);

		i++;
		progress_bar(i, processing_total, "Processing");

	} while ( rval == 0 );

	return n_peaks;
}


static int add_features_to_histogram(struct image *image,
                                     struct histogram_info *info,
                                     struct bin_stats *histdata)
{
	int j;
	int n_peaks;

	n_peaks = 0;
	for ( j=0; j<image_feature_count(image->features); j++) {

		struct imagefeature *f;

		f = image_get_feature(image->features, j);
		if ( !f->valid ) continue;

		if ( !in_bad_region(image->det, f->fs,f->ss) ) {

			int r;

			r = add_peak_to_histogram(f->fs, f->ss, f->intensity,
			                          image, info, histdata);

			if ( !r ) n_peaks++;

		}

	}

	return n_peaks;
}


static unsigned int process_stream_peaks(FILE *fh, struct image *image,
                                         struct histogram_info *info,
                                         struct bin_stats *histdata,
                                         unsigned int *n_patterns,
                                         int only_indexed)
{
	int rval;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	unsigned int processing_total;

	processing_total = count_patterns(fh);
	rewind(fh);

	do {

		/* Get data from next chunk */
		rval = read_chunk(fh, image);
		if ( rval ) continue;

		if ( image->features != NULL ) {

			if ( (!only_indexed)
			  || ( only_indexed && (image->reflections != NULL)) )
			{
				(*n_patterns)++;
				n_peaks += add_features_to_histogram(image,
				                                     info,
				                                     histdata);
			}
		}

		free(image->filename);
		reflist_free(image->reflections);
		image_feature_list_free(image->features);
		cell_free(image->indexed_cell);

		i++;
		progress_bar(i, processing_total, "Processing");

	} while ( rval == 0 );

	return n_peaks;
}


static unsigned int process_stream_h5(FILE *fh, struct image *image,
                                      struct histogram_info *info,
                                      struct bin_stats *histdata,
                                      int config_satcorr, int only_indexed,
                                      unsigned int *n_patterns,
                                      const char *element)
{
	int fs, ss, rval;
	double intensity;
	unsigned int i = 0;
	unsigned int n_peaks = 0;
	struct hdfile *hdfile = NULL;
	unsigned int processing_total;

	processing_total = count_patterns(fh);
	rewind(fh);

	do {

		/* Get data from next chunk */
		rval = read_chunk(fh, image);
		if ( rval ) continue;
		if ( !only_indexed ||
		   ( only_indexed && (image->reflections != NULL) ) )
		{
			hdfile = hdfile_open(image->filename);
			if ( element != NULL ) {
				int r;
				r = hdfile_set_image(hdfile, element);
				if ( r ) {
					ERROR("Couldn't select path '%s'\n",
					     element);
					hdfile_close(hdfile);
					return 0;
				}
			} else {
				int r;
				r = hdfile_set_first_image(hdfile, "/");
				if ( r ) {
					ERROR("Couldn't select first path\n");
					hdfile_close(hdfile);
					return 0;
				}

			}
			hdf5_read(hdfile, image, config_satcorr);
			hdfile_close(hdfile);
			n_patterns++;
			for ( ss=0; ss<image->height; ss++ ) {
			for ( fs=0; fs<image->width;  fs++ ) {

				intensity = image->data[fs + image->width*ss];
				if ( !in_bad_region(image->det,fs,ss) ) {

					add_peak_to_histogram(fs, ss, intensity,
					                      image, info,
					                      histdata);

				}

			}
			}
		}

		free(image->data);
		free(image->filename);
		reflist_free(image->reflections);
		image_feature_list_free(image->features);
		cell_free(image->indexed_cell);

		i++;
		progress_bar(i, processing_total, "Processing");

	} while ( rval == 0 );

	return n_peaks;
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Plot a powder pattern as a 1D graph using the detector geometry.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. (*.stream, *.hkl, or *.h5)\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"  -g. --geometry=<file>   Get detector geometry from <file>.\n"
"  -b, --beam=<file>       Get beam parameters (wavelength) from <file>.\n"
"  -p, --pdb=<file>        Get unit cell from PDB file. (.hkl files only)\n"
"  -y, --symmetry=<sym>    The symmetry of crystal (.hkl files only)\n"
"  -s, --bins=n            Makes histogram with n bins (default is 100).\n"
"      --spacing=<type>    Use 'type' to select the q spacing.\n"
"                          Choose from:\n"
"                            linear      : linear (default)\n"
"                            q2          : even spacing in Wilson plots\n"
"                            volume      : constant volume\n"
"      --q-max=n           The maximum q to be considered in plot.\n"
"      --q-min=n           The minimum q to be considered in plot.\n"
"  -d, --data=<type>       Use to select the kind of stream data in histogram.\n"
"                          Choose from:\n"
"                            reflection  : uses peak positons from indexed\n"
"                                          reflection \n"
"                            hkl         : uses the hkl list from indexed\n"
"                                          reflections (requires pdb file)\n"
"                            d           : uses the 1/d list from indexed\n"
"                                          reflections (default)\n"
"                            peaks       : uses all peaks found from peak\n"
"                                          search\n"
"                            h5          : all points in h5, excluding bad\n"
"                                          regions\n"
"     --no-sat-corr        Don't correct values of saturated peaks using a\n"
"                            table included in the HDF5 file.\n"
"     --only-indexed       Use with -data=peaks or h5 if you want to use the\n"
"                            peak list of only indexed patterns\n"
"     --no-q-scaling       Use with .hkl files if you want to not scale the\n"
"                            powder by 1/q^2\n"
"     --ring-corr          Use if you want to scale the powder plot to\n"
"                            correct for the fractional area sampled of the\n"
"                            powder ring\n"
"     --use-redundancy     Use with .hkl files if you want to use the number\n"
"                            of measurements and not the number of symetrical\n"
"                            equivelent reflections as the number of time a\n"
"                            reflection occurs in the powder\n"
" -e, --image=<element>    Use this image when reading an HDF5 file.\n"
"                           Example: /data/data0.\n"
"                           Default: The first one found.\n"

"\n");
}


int main(int argc, char *argv[])
{
	FILE *fh = NULL;
	UnitCell *cell = NULL;
	struct image image;
	struct hdfile *hdfile = NULL;
	struct bin_stats *histdata = NULL;
	struct histogram_info hist_info;
	SymOpList *sym;

	/* Default settings */
	hist_info.histsize =  100;
	hist_info.q_min    = -1.0;
	hist_info.q_max    = -1.0;
	hist_info.spacing  = LINEAR;
	image.lambda = -1.0;
	image.beam = NULL;
	image.det  = NULL;
	char *element = NULL;

	unsigned int n_patterns = 0;
	unsigned int n_peaks    = 0;

	int c, rval, file_type, data_type;
	int config_satcorr = 1;
	int need_geometry  = 0;
	int need_beam      = 0;
	int need_pdb       = 0;
	int only_indexed   = 0;
	int q_scaling      = 1;
	int ring_corr      = 0;
	int use_redundancy = 0;
	unsigned int i;

	char *filename = NULL;
	char *geometry = NULL;
	char *beamf    = NULL;
	char *pdb      = NULL;
	char *output   = NULL;
	char *datatype = NULL;
	char *sym_str  = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"pdb",                1, NULL,               'p'},
		{"symmetry",           1, NULL,               'y'},
		{"bins",               1, NULL,               's'},
		{"q-max",              1, NULL,                1 },
		{"q-min",              1, NULL,                2 },
		{"spacing",            1, NULL,                3 },
		{"no-sat-corr",        0, &config_satcorr,     0 },
		{"sat-corr",           0, &config_satcorr,     1 },
		{"only-indexed",       0, &only_indexed,       1 },
		{"no-q-scaling",       0, &q_scaling,          0 },
		{"ring-corr",          0, &ring_corr,          1 },
		{"use-redundancy",     0, &use_redundancy,     1 },
		{"data",               1, NULL,               'd'},
		{"image",              1, NULL,               'e'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:g:b:p:s:d:y:",
	        longopts, NULL)) != -1)
	{

		switch (c) {

		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			filename = strdup(optarg);
			break;

		case 'o' :
			output = strdup(optarg);
			break;

		case 'g' :
			geometry = strdup(optarg);
			break;

		case 'b' :
			beamf = strdup(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
			break;

		case 'y' :
			sym_str = strdup(optarg);
			break;

		case 's' :
			hist_info.histsize = atoi(optarg);
			break;

		case 'e' :
			element = strdup(optarg);
			break;

		case 1 :
			hist_info.q_max = atof(optarg);
			break;

		case 2 :
			hist_info.q_min = atof(optarg);
			break;

		case 3 :
			if (strcmp(optarg, "linear") == 0 ) {
				hist_info.spacing = LINEAR;
			} else if (strcmp(optarg, "q2") == 0 ) {
				hist_info.spacing = q2;
			} else if (strcmp(optarg, "volume") == 0) {
				hist_info.spacing = VOLUME;
			} else {
				ERROR("Invalid spacing type: '%s'\n", optarg);
				return 1;
			}
			break;

		case 'd' :
			datatype = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	/* Process input file type */
	if ( filename == NULL ) {

		ERROR("You must specify the input filename with -i\n");
		return 1;

	}

	if ( is_stream(filename) == 1 ) {

		file_type = FILE_STREAM;
		
	} else if ( H5Fis_hdf5(filename) > 0 ) {

		file_type = FILE_H5;
		need_geometry = 1;

	} else {

		image.reflections = read_reflections(filename);

		if ( image.reflections != NULL ) {
			file_type = FILE_HKL;
			need_pdb = 1;
			need_geometry = 0;
			need_beam = 0;
			image.lambda = 0.0;
		} else {
			ERROR("Couldn't recognise %s as reflection list,"
			      " stream or image.\n", filename);
			return 1;
		}

	}

	if ( datatype == NULL ) {
		data_type = PLOT_D;
		if ((hist_info.q_min < 0.0) || (hist_info.q_max < 0.0)) {
			need_geometry = 1;
		}

	} else if ( strcmp(datatype, "reflection") == 0 ) {
		data_type = PLOT_REFL;
		need_geometry = 1;

	} else if ( strcmp(datatype, "hkl") == 0 ) {
		data_type = PLOT_HKL;
		need_pdb = 1;

	} else if ( strcmp(datatype, "d") == 0 ) {
		data_type = PLOT_D;
		if ((hist_info.q_min < 0.0) || (hist_info.q_max < 0.0)) {
			need_geometry = 1;
		}

	} else if ( strcmp(datatype, "peaks") == 0 ) {
		data_type = PLOT_PEAKS;
		need_geometry = 1;

	} else if ( strcmp(datatype, "h5") == 0 ) {
		data_type = PLOT_H5;
		need_geometry = 1;

	} else {

		ERROR("Failed to read data plot type: '%s'\n", datatype);
		return 1;
	}

	/* doubt this is needed, but double check just in case */
	if ( file_type == FILE_HKL ) {
		need_geometry = 0;
		need_beam = 0;
	}

	/* Logic checks */
	if ( need_geometry && (image.lambda < 0.0) ) {
		need_beam = 1;
	}
	if ( hist_info.histsize <= 0 ) {
		ERROR("You need to specify a histogram with more then 0 "
                      "bins\n");
		return 1;
	}

	/* Get geometry, beam and pdb files and parameters as needed */
	if ( need_geometry ) {
		if ( geometry == NULL ) {
			ERROR("You need to specify a geometry file with "
		              "--geometry\n");
			return 1;
		} else {
			image.det = get_detector_geometry(geometry);
			if ( image.det == NULL ) {
				ERROR("Failed to read detector geometry "
				      "from '%s'\n", geometry);
				return 1;
			}
			image.width  = image.det->max_fs;
			image.height = image.det->max_ss;
		}
	}
	free(geometry);

	/* Open files to get wavelength if it exists & camera length 
	   if they are not found in the geometry file */
	if (file_type == FILE_STREAM) {
		fh = fopen(filename, "r");
		if ( fh == NULL ) {
			ERROR("Failed to open input file\n");
			return 1;
		}
		/* Use wavelength from first chunk */
		rval = read_chunk(fh, &image);
		rewind(fh);
	} else if (file_type == FILE_H5) {
		hdfile = hdfile_open(filename);
		if ( element != NULL ) {
			int r;
			r = hdfile_set_image(hdfile, element);
			if ( r ) {
				ERROR("Couldn't select path '%s'\n",
				     element);
				hdfile_close(hdfile);
				return 0;
			}
		} else {
			int r;
			r = hdfile_set_first_image(hdfile, "/");
			if ( r ) {
				ERROR("Couldn't select first path\n");
				hdfile_close(hdfile);
				return 0;
			}

		}
		hdf5_read(hdfile, &image, config_satcorr);
		hdfile_close(hdfile);
	}
	free(filename);

	/* Logic checks */
	if ( need_geometry && (image.lambda < 0.0) ) {
		need_beam = 1;
	}
	if ( hist_info.histsize <= 0 ) {
		ERROR("You need to specify a histogram with more then 0 "
                      "bins\n");
		return 1;
	}

	if ( need_beam ) {

		if ( beamf == NULL ) {
			ERROR("No wavelength in file, so you need to specify "
			      "a beam parameters file with --beam\n");
			return 1;
		} else {
			image.beam = get_beam_parameters(beamf);
			if ( image.beam == NULL ) {
				ERROR("Failed to read beam from '%s'\n",
				      beamf);
				return 1;
			}
			image.lambda = ph_en_to_lambda(eV_to_J(
			                 image.beam->photon_energy));
		}

	}
	free(beamf);

	if ( need_pdb ) {

		if (pdb == NULL) {
			ERROR("You need to specify a pdb file with --pdb.\n");
			return 1;
		} else {
			cell = load_cell_from_pdb(pdb);
			if ( cell == NULL ) {
				ERROR("Couldn't read unit cell (from %s)\n",
				       pdb);
				return 1;
			}
		}
	}
	free(pdb);

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	/* Set up histogram info*/
	if ( file_type == FILE_HKL || data_type == PLOT_HKL) {
		/* get q range from Miller indices in hkl
		   file. */
		if ((hist_info.q_min < 0.0) && (hist_info.q_max < 0.0)) {
			resolution_limits(image.reflections, cell,
		                  &hist_info.q_min, &hist_info.q_max);
		} else if (hist_info.q_min < 0.0) {
			double dummy;
			resolution_limits(image.reflections, cell,
		                  &hist_info.q_min, &dummy);
		} else if (hist_info.q_max < 0.0) {
			double dummy;
			resolution_limits(image.reflections, cell,
		                  &dummy, &hist_info.q_max);
		}
	} else {
		if ( hist_info.q_min < 0.0 ) {
			hist_info.q_min = smallest_q(&image);
		}
		if ( hist_info.q_max < 0.0 ) {
			hist_info.q_max = largest_q(&image);
		}
	}

	if ( hist_info.q_min >= hist_info.q_max ) {
		ERROR("the minimum q value of: %e "
	              "is greator then your max q value of: %e\n",
                      hist_info.q_min, hist_info.q_max);
		return 1;
	}

	if ( hist_info.spacing == LINEAR) {
		hist_info.q_delta = (hist_info.q_max - hist_info.q_min)/
		                    hist_info.histsize;
	} else if ( hist_info.spacing == q2) {
		hist_info.q_delta = (pow(hist_info.q_max, 2.0) -
		                     pow(hist_info.q_min, 2.0)) /
		                     hist_info.histsize;
	} else {  //by default must be in VOLUME
		hist_info.q_delta = (pow(hist_info.q_max, 3.0) -
		                     pow(hist_info.q_min, 3.0)) /
		                     hist_info.histsize;
	}
	/* Set up histogram data */
	histdata = malloc((hist_info.histsize) * sizeof(struct bin_stats));
	histogram_setup(&hist_info, histdata);

	/* Set up ring scaling */
	if ( ring_corr ) {

		if ( ring_fraction_calc(&hist_info, histdata, &image) ) {

			ERROR("Detector geometry is required to correct for"
			      " finite ring area.\n");
			return 1;
		}
	}

	/* Process reflections based on file type and data type */
	switch (file_type) {
	case FILE_H5 :
		n_patterns++;
		n_peaks = process_h5(&image, &hist_info, histdata);
		free(image.data);
		break;

	case FILE_HKL :
		n_patterns++; //inc number of patterns used
		n_peaks = process_hkl(&image, sym, cell, &hist_info, histdata,
                                      q_scaling, use_redundancy);
		break;

	case FILE_STREAM :
		switch (data_type) {
		case PLOT_REFL :
			n_peaks = process_stream_reflection(fh, &image,
			             &hist_info, histdata, &n_patterns);
			break;
		case PLOT_D :
			n_peaks = process_stream_d(fh, &image, &hist_info,
				     histdata, &n_patterns);
			break;
		case PLOT_HKL :
			n_peaks = process_stream_hkl(fh, &image, &hist_info,
				     histdata, cell, &n_patterns);
			break;
		case PLOT_PEAKS :
			n_peaks = process_stream_peaks(fh, &image, &hist_info,
			             histdata, &n_patterns, only_indexed);
			break;
		case PLOT_H5 :
			n_peaks =  process_stream_h5(fh, &image, &hist_info,
			             histdata, config_satcorr, only_indexed,
			             &n_patterns, element);
			break;
		default :
			break;
		}
		fclose(fh);
		break;

	default :
		break;

	}

	/* Sqrt the variance to get the std_dev */
	for( i=0; i<hist_info.histsize; i++ ) {
		if (histdata[i].N > 1) {
			histdata[i].std_dev = sqrt(histdata[i].std_dev/
			                          (histdata[i].N-1));
		}
	}
	if ( ring_corr ) {
		for( i=0; i<hist_info.histsize; i++ ) {
			histdata[i].N       = histdata[i].N      /
			                      histdata[i].fract;
			histdata[i].total   = histdata[i].total  /
			                      histdata[i].fract;
			histdata[i].mean    = histdata[i].mean   /
			                      histdata[i].fract;
			histdata[i].std_dev = histdata[i].total  /
			                      histdata[i].fract;
		}
	}

	/* Print out the results */
	if ( output != NULL ) {
		fh = fopen(output, "w");
		if ( fh == NULL ) {
			ERROR("Failed to open output file\n");
			return 1;
		}
	} else {
		fh = stdout;
	}

	fprintf(fh, "I read %i patterns with %i peaks\n", n_patterns, n_peaks);
	fprintf(fh, "q\tN\ttotal\tmean\tstd dev\t std dev of mean\n");

	for( i=0; i<hist_info.histsize; i++ ) {
		fprintf(fh, "%5e\t%i\t%5e\t%5e\t%5e\t%5e\n",
		        histdata[i].q_min, histdata[i].N,
		        histdata[i].total, histdata[i].mean,
		        histdata[i].std_dev,
		        histdata[i].std_dev/sqrt(histdata[i].N));
	}

	if ( cell != NULL ) cell_free(cell);
	if ( image.det  != NULL ) free(image.det);
	if ( image.beam != NULL ) free(image.beam);
	fclose(fh);
	free(histdata);

	return 0;
}
