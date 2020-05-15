/*
 * image.c
 *
 * Handle images and image features
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011-2017 Thomas White <taw@physics.org>
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
#include <zlib.h>

#include "image.h"
#include "utils.h"
#include "events.h"
#include "hdf5-file.h"
#include "detector.h"
#include "detgeom.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"

/** \file image.h */

struct imagefile
{
	enum imagefile_type type;
	char *filename;
	struct hdfile *hdfile;
};


struct _imagefeaturelist
{
	struct imagefeature *features;
	int                  max_features;
	int                  n_features;
};


void image_add_feature(ImageFeatureList *flist, double fs, double ss,
                       int pn,
                       struct image *parent, double intensity, const char *name)
{
	if ( flist->n_features == flist->max_features ) {
		struct imagefeature *nf;
		int nmf = flist->max_features + 128;
		nf = realloc(flist->features, nmf*sizeof(struct imagefeature));
		if ( nf == NULL ) return;
		flist->features = nf;
		flist->max_features = nmf;
	}

	flist->features[flist->n_features].fs = fs;
	flist->features[flist->n_features].ss = ss;
	flist->features[flist->n_features].pn = pn;
	flist->features[flist->n_features].intensity = intensity;
	flist->features[flist->n_features].name = name;

	flist->n_features++;
}


ImageFeatureList *image_feature_list_new()
{
	ImageFeatureList *flist;

	flist = malloc(sizeof(ImageFeatureList));

	flist->n_features = 0;
	flist->max_features = 0;
	flist->features = NULL;

	return flist;
}


static int comp(const void *a, const void *b)
{
	const struct imagefeature *ap = a;
	const struct imagefeature *bp = b;

	return ap->intensity < bp->intensity;
}


ImageFeatureList *image_feature_list_copy(const ImageFeatureList *flist)
{
	ImageFeatureList *n;
	int nf, i;

	if ( flist == NULL ) return NULL;

	n = image_feature_list_new();
	if ( n == NULL ) return NULL;

	n->features = malloc(flist->n_features*sizeof(struct imagefeature));
	if ( n->features == NULL ) {
		free(n);
		return NULL;
	}

	nf = 0;
	for ( i=0; i<flist->n_features; i++ ) {
		const struct imagefeature *f;
		f = image_get_feature_const(flist, i);
		if ( f == NULL ) continue;
		n->features[nf++] = flist->features[i];
	}
	n->n_features = nf;

	return n;
}


/**
 * Strongest first.
 */
ImageFeatureList *sort_peaks(ImageFeatureList *flist)
{
	ImageFeatureList *n = image_feature_list_copy(flist);
	qsort(n->features, image_feature_count(n),
	      sizeof(struct imagefeature), comp);
	return n;
}


void image_feature_list_free(ImageFeatureList *flist)
{
	if ( flist == NULL ) return;
	free(flist->features);
	free(flist);
}


struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                           double fs, double ss,
                                           int pn, double *d, int *idx)
{
	int i;
	double dmin = +HUGE_VAL;
	int closest = 0;

	for ( i=0; i<flist->n_features; i++ ) {

		double ds;

		if ( pn != flist->features[i].pn ) continue;

		ds = distance(flist->features[i].fs, flist->features[i].ss,
		              fs, ss);

		if ( ds < dmin ) {
			dmin = ds;
			closest = i;
		}

	}

	if ( dmin < +HUGE_VAL ) {
		*d = dmin;
		*idx = closest;
		return &flist->features[closest];
	}

	*d = +INFINITY;
	return NULL;
}


Reflection *image_reflection_closest(RefList *rlist,
                                     double fs, double ss, struct panel *p,
                                     struct detector *det,
                                     double *d)
{

	double dmin = HUGE_VAL;
	Reflection *closest = NULL;
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(rlist, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double ds;
		struct panel *p2;
		double rfs, rss;

		get_detector_pos(refl, &rfs, &rss);
		p2 = get_panel(refl);

		if ( p != p2 ) continue;

		ds = distance(rfs, rss, fs, ss);

		if ( ds < dmin ) {
			dmin = ds;
			closest = refl;
		}

	}

	if ( dmin < +HUGE_VAL ) {
		*d = dmin;
		return closest;
	}

	*d = +INFINITY;
	return NULL;
}


int image_feature_count(ImageFeatureList *flist)
{
	if ( flist == NULL ) return 0;
	return flist->n_features;
}


const struct imagefeature *image_get_feature_const(const ImageFeatureList *flist,
                                                   int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx >= flist->n_features ) return NULL;

	return &flist->features[idx];
}


struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx >= flist->n_features ) return NULL;

	return &flist->features[idx];
}


void image_remove_feature(ImageFeatureList *flist, int idx)
{
	memmove(&flist->features[idx], &flist->features[idx+1],
	        (flist->n_features-idx-1)*sizeof(struct imagefeature));
	flist->n_features--;
}


void image_add_crystal(struct image *image, Crystal *cryst)
{
	Crystal **crs;
	int n;

	n = image->n_crystals;
	crs = realloc(image->crystals, (n+1)*sizeof(Crystal *));
	if ( crs == NULL ) {
		ERROR("Failed to allocate memory for crystals.\n");
		return;
	}

	crs[n] = cryst;
	image->crystals = crs;
	image->n_crystals = n+1;
}


int remove_flagged_crystals(struct image *image)
{
	int i;
	int n_bad = 0;

	for ( i=0; i<image->n_crystals; i++ ) {
		if ( crystal_get_user_flag(image->crystals[i]) ) {
			int j;
			Crystal *deleteme = image->crystals[i];
			cell_free(crystal_get_cell(deleteme));
			crystal_free(deleteme);
			for ( j=i; j<image->n_crystals-1; j++ ) {
				image->crystals[j] = image->crystals[j+1];
			}
			image->n_crystals--;
			n_bad++;
			i--;
		}
	}

	return n_bad;
}


/* Free all crystals, including their RefLists and UnitCells */
void free_all_crystals(struct image *image)
{
	int i;
	if ( image->crystals == NULL ) return;
	for ( i=0; i<image->n_crystals; i++ ) {
		Crystal *cr = image->crystals[i];
		reflist_free(crystal_get_reflections(cr));
		cell_free(crystal_get_cell(cr));
		crystal_free(image->crystals[i]);
	}
	free(image->crystals);
	image->n_crystals = 0;
}


/**************************** Image field lists *******************************/

struct imagefile_field_list
{
	char **fields;
	int n_fields;
	int max_fields;
};


struct imagefile_field_list *new_imagefile_field_list()
{
	struct imagefile_field_list *n;

	n = calloc(1, sizeof(struct imagefile_field_list));
	if ( n == NULL ) return NULL;

	n->max_fields = 32;
	n->fields = malloc(n->max_fields*sizeof(char *));
	if ( n->fields == NULL ) {
		free(n);
		return NULL;
	}

	return n;
}


void free_imagefile_field_list(struct imagefile_field_list *n)
{
	int i;
	for ( i=0; i<n->n_fields; i++ ) {
		free(n->fields[i]);
	}
	free(n->fields);
	free(n);
}


void add_imagefile_field(struct imagefile_field_list *copyme, const char *name)
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


/******************************* CBF files ************************************/

static int unpack_panels(struct image *image, float *data, int data_width,
                         int data_height)
{
	int pi;

	/* FIXME: Load these masks from an HDF5 file, if filenames are
	 * given in the geometry file */
	uint16_t *flags = NULL;
	float *sat = NULL;

	image->dp = malloc(image->det->n_panels * sizeof(float *));
	image->bad = malloc(image->det->n_panels * sizeof(int *));
	image->sat = malloc(image->det->n_panels * sizeof(float *));
	if ( (image->dp == NULL) || (image->bad == NULL)
	  || (image->sat == NULL) )
	{
		ERROR("Failed to allocate panels.\n");
		return 1;
	}

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		struct panel *p;
		int fs, ss;

		p = &image->det->panels[pi];
		image->dp[pi] = malloc(p->w*p->h*sizeof(float));
		image->bad[pi] = calloc(p->w*p->h, sizeof(int));
		image->sat[pi] = malloc(p->w*p->h*sizeof(float));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL)
		  || (image->sat[pi] == NULL) )
		{
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		if ( p->mask != NULL ) {
			ERROR("WARNING: Bad pixel masks do not currently work "
			      "with CBF files\n");
			ERROR(" (bad pixel regions specified in the geometry "
			      "file will be used, however)\n");
		}

		if ( p->satmap != NULL ) {
			ERROR("WARNING: Saturation maps do not currently work "
			      "with CBF files\n");
		}

		if ( (p->orig_min_fs + p->w > data_width)
		  || (p->orig_min_ss + p->h > data_height) )
		{
			ERROR("Panel %s is outside range of data in CBF file\n",
			      p->name);
			return 1;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int idx;
			int cfs, css;
			int bad = 0;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p->w*ss] = data[idx];

			if ( sat != NULL ) {
				image->sat[pi][fs+p->w*ss] = sat[idx];
			} else {
				image->sat[pi][fs+p->w*ss] = INFINITY;
			}

			if ( p->no_index ) bad = 1;

			if ( in_bad_region(image->det, p, cfs, css) ) {
				bad = 1;
			}

			if ( isnan(data[idx]) || isinf(data[idx]) ) bad = 1;

			if ( flags != NULL ) {

				int f;

				f = flags[idx];

				/* Bad if it's missing any of the "good" bits */
				if ( (f & image->det->mask_good)
				       != image->det->mask_good ) bad = 1;

				/* Bad if it has any of the "bad" bits. */
				if ( f & image->det->mask_bad ) bad = 1;

			}
			image->bad[pi][fs+p->w*ss] = bad;

		}
		}

	}

	return 0;
}


static void cbf_fill_in_beam_parameters(struct beam_params *beam,
                                        struct imagefile *f,
                                        struct image *image)
{
	double eV;

	if ( beam->photon_energy_from == NULL ) {

		/* Explicit value given */
		eV = beam->photon_energy;

	} else {

		ERROR("Can't get photon energy from CBF yet.\n");
		eV = 0.0;

	}

	image->lambda = ph_en_to_lambda(eV_to_J(eV))*beam->photon_energy_scale;
}


static void cbf_fill_in_clen(struct detector *det, struct imagefile *f)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];

		if ( p->clen_from != NULL ) {

			ERROR("Can't get clen from CBF yet.\n");

		}

		adjust_centering_for_rail(p);

	}
}


static void add_out(float val, float *data_out, int nmemb_out,
                    int *outpos, int *nrej)
{
	if ( *outpos < nmemb_out ) {
		data_out[(*outpos)++] = val;
	} else {
		(*nrej)++;
	}
}


/* Reverses byte offset compression and converts to single precision float.
 * Note that this compression scheme specifies the data format of the input
 * data, therefore the X-Binary-Element-Type is completely ignored. */
static void decode_cbf_byte_offset(float *data_out, int nmemb_out,
                                   const int8_t *data_in, const size_t n)
{
	int inpos = 0;
	int outpos = 0;
	int nrej = 0;
	float val = 0.0;

	while ( inpos < n ) {

		int64_t delta = data_in[inpos++];

		if ( (delta >= -127) && (delta <= 127) ) {
			val += delta;
			add_out(val, data_out, nmemb_out, &outpos, &nrej);
			continue;
		}

		delta = *(int16_t *)(data_in+inpos);
		inpos += 2;

		if ( (delta >= -32767) && (delta <= 32767) ) {
			val += delta;
			add_out(val, data_out, nmemb_out, &outpos, &nrej);
			continue;
		}

		delta = *(int32_t *)(data_in+inpos);
		inpos += 4;

		if ( (delta >= -2147483647) && (delta <= 2147483647) ) {
			val += delta;
			add_out(val, data_out, nmemb_out, &outpos, &nrej);
			continue;
		}

		delta = *(int64_t *)(data_in+inpos);
		inpos += 8;
		val += delta;
		add_out(val, data_out, nmemb_out, &outpos, &nrej);

	}

	if ( nrej > 0 ) {
		STATUS("%i elements rejected\n", nrej);
	}
}


static int binary_start(char *data)
{
	char *datac = data;
	if ( (datac[0] == (char)0x0c) && (datac[1] == (char)0x1a)
	  && (datac[2] == (char)0x04) && (datac[3] == (char)0xd5) ) return 1;
	return 0;
}


enum cbf_data_conversion
{
	CBF_NO_CONVERSION,
	CBF_BYTE_OFFSET,
	CBF_PACKED,
	CBF_CANONICAL
};

enum cbf_data_type
{
	CBF_NO_TYPE,
	CBF_ELEMENT_U8,
	CBF_ELEMENT_S8,
	CBF_ELEMENT_U16,
	CBF_ELEMENT_S16,
	CBF_ELEMENT_U32,
	CBF_ELEMENT_S32,
	CBF_ELEMENT_F32,
	CBF_ELEMENT_F64,
};


static enum cbf_data_type parse_element_type(const char *t)
{
	if ( strstr(t, "signed 8-bit integer") != NULL )
	{
		return CBF_ELEMENT_S8;
	}

	if ( strstr(t, "unsigned 8-bit integer") != NULL )
	{
		return CBF_ELEMENT_U8;
	}

	if ( strstr(t, "signed 16-bit integer") != NULL )
	{
		return CBF_ELEMENT_S16;
	}

	if ( strstr(t, "unsigned 16-bit integer") != NULL )
	{
		return CBF_ELEMENT_U16;
	}

	if ( strstr(t, "signed 32-bit integer") != NULL )
	{
		return CBF_ELEMENT_S32;
	}

	if ( strstr(t, "unsigned 32-bit integer") != NULL )
	{
		return CBF_ELEMENT_U32;
	}

	if ( strstr(t, "signed 32-bit real IEEE") != NULL )
	{
		return CBF_ELEMENT_F32;
	}

	if ( strstr(t, "signed 64-bit real IEEE") != NULL )
	{
		return CBF_ELEMENT_F64;
	}

	/* complex type is unsupported */

	return CBF_NO_TYPE;
}


static size_t element_size(enum cbf_data_type t)
{
	switch ( t ) {
		case CBF_ELEMENT_S8  : return 1;
		case CBF_ELEMENT_U8  : return 1;
		case CBF_ELEMENT_S16 : return 2;
		case CBF_ELEMENT_U16 : return 2;
		case CBF_ELEMENT_S32 : return 4;
		case CBF_ELEMENT_U32 : return 4;
		case CBF_ELEMENT_F32 : return 4;
		case CBF_ELEMENT_F64 : return 8;
		default : return 0;
	}
}



static int convert_type(float *data_out, long nmemb_exp,
                        enum cbf_data_type eltype,
                        void *data_in, size_t data_in_len)
{
	long int i;
	long int o = 0;
	size_t elsize = element_size(eltype);

	if ( elsize == 0 ) return 1;

	if ( nmemb_exp * elsize > data_in_len ) {
		ERROR("Not enough CBF data for image size/type!\n");
		return 1;
	}

	for ( i=0; i<nmemb_exp; i++ ) {
		switch ( eltype ) {

			case CBF_ELEMENT_S8:
			data_out[o++] = ((int8_t *)data_in)[i];
			break;

			case CBF_ELEMENT_U8:
			data_out[o++] = ((uint8_t *)data_in)[i];
			break;

			case CBF_ELEMENT_S16:
			data_out[o++] = ((int16_t *)data_in)[i];
			break;

			case CBF_ELEMENT_U16:
			data_out[o++] = ((uint16_t *)data_in)[i];
			break;

			case CBF_ELEMENT_S32:
			data_out[o++] = ((int32_t *)data_in)[i];
			break;

			case CBF_ELEMENT_U32:
			data_out[o++] = ((uint32_t *)data_in)[i];
			break;

			case CBF_ELEMENT_F32:
			data_out[o++] = ((float *)data_in)[i];
			break;

			case CBF_ELEMENT_F64:
			data_out[o++] = ((double *)data_in)[i];
			break;

			case CBF_NO_TYPE:
			break;
		}
	}

	return 0;
}


static float *read_cbf_data(const char *filename, int gz, int *w, int *h)
{
	FILE *fh;
	void *buf = NULL;
	char *rval;
	size_t data_compressed_len = 0;
	float *data_out = NULL;
	enum cbf_data_conversion data_conversion = CBF_NO_CONVERSION;
	enum cbf_data_type data_type = CBF_ELEMENT_U32;  /* ITG (2006) 2.3.3.3 */
	int in_binary_section = 0;

	*w = 0;
	*h = 0;

	if ( !gz ) {

		fh = fopen(filename, "rb");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", filename);
			return NULL;
		}

	} else {

		gzFile gzfh;
		size_t len, len_read;
		const size_t bufinc = 8*1024*1024;  /* Allocate buffer in 8Mb chunks */
		size_t bufsz = bufinc;

		gzfh = gzopen(filename, "rb");
		if ( gzfh == NULL ) return NULL;

		#ifdef HAVE_GZBUFFER
		/* Set larger buffer size for hopefully faster uncompression */
		gzbuffer(gzfh, 128*1024);
		#endif

		buf = malloc(bufsz);
		if ( buf == NULL ) return NULL;

		len = 0;
		do {

			len_read = gzread(gzfh, buf+len, bufinc);
			if ( len_read == -1 ) return NULL;
			len += len_read;

			if ( len_read == bufinc ) {
				bufsz += bufinc;
				buf = realloc(buf, bufsz);
				if ( buf == NULL ) return NULL;
			}

		} while ( len_read == bufinc );

		fh = fmemopen(buf, len, "rb");
		if ( fh == NULL ) return NULL;

		gzclose(gzfh);

	}

	/* This is really horrible, but there are at least three different types
	 * of header mingled together (CIF, MIME, DECTRIS), so a real parser
	 * would be very complicated and much more likely to have weird bugs. */
	do {

		char line[1024];
		long line_start;

		line_start = ftell(fh);
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;
		chomp(line);

		if ( strcmp(line, "--CIF-BINARY-FORMAT-SECTION--") == 0 ) {
			in_binary_section = 1;
		}

		if ( strcmp(line, "--CIF-BINARY-FORMAT-SECTION----") == 0 ) {
			in_binary_section = 0;
		}

		if ( in_binary_section ) {

			if ( strncmp(line, "X-Binary-Size: ", 15) == 0 ) {
				data_compressed_len = atoi(line+15);
			}

			if ( strncmp(line, "X-Binary-Element-Byte-Order: ", 29) == 0 ) {
				const char *elbo = line+29;
				if ( strcmp(elbo, "LITTLE_ENDIAN") != 0 ) {
					ERROR("Unsupported endianness: %s\n", elbo);
					free(buf);
					fclose(fh);
					return NULL;
				}
			}

			/* Try to spot compression algorithm */
			if ( strstr(line, "conversions=\"x-CBF_BYTE_OFFSET\"") != NULL ) {
				data_conversion = CBF_BYTE_OFFSET;
			} else if ( strstr(line, "conversions=\"x-CBF_CANONICAL\"") != NULL ) {
				data_conversion = CBF_CANONICAL;
			} else if ( strstr(line, "conversions=\"x-CBF_PACKED\"") != NULL ) {
				data_conversion = CBF_PACKED;
			} else if ( strstr(line, "conversions=") != NULL ) {
				ERROR("Unrecognised CBF content conversion: %s\n", line);
				free(buf);
				fclose(fh);
				return NULL;
			}

			/* Likewise, element type */
			if ( strncmp(line, "X-Binary-Element-Type: ", 23) == 0 )
			{
				const char *eltype = (line+23);
				data_type = parse_element_type(eltype);
				if ( data_type == CBF_NO_TYPE ) {
					ERROR("Unrecognised element type: %s\n",
					      eltype);
					free(buf);
					fclose(fh);
					return NULL;
				}
			}

			if ( strncmp(line, "X-Binary-Size-Fastest-Dimension: ", 33) == 0 ) {
				*w = atoi(line+33);
			}

			if ( strncmp(line, "X-Binary-Size-Second-Dimension: ", 32) == 0 ) {
				*h = atoi(line+32);
			}

		}

		if ( in_binary_section && binary_start(line) ) {

			size_t len_read;
			int nmemb_exp;
			void *data_compressed;
			int r = 0;

			if ( data_compressed_len == 0 ) {
				ERROR("Found CBF data before X-Binary-Size!\n");
				free(buf);
				fclose(fh);
				return NULL;
			}

			if ( (*w == 0) || (*h == 0) ) {
				ERROR("Found CBF data before dimensions!\n");
				free(buf);
				fclose(fh);
				return NULL;
			}

			if ( data_compressed_len > 100*1024*1024 ) {
				ERROR("Stated CBF data size too big\n");
				free(buf);
				fclose(fh);
				return NULL;
			}

			data_compressed = malloc(data_compressed_len);
			if ( data_compressed == NULL ) {
				ERROR("Failed to allocate memory for CBF data\n");
				free(buf);
				fclose(fh);
				return NULL;
			}

			fseek(fh, line_start+4, SEEK_SET);
			len_read = fread(data_compressed, 1, data_compressed_len, fh);
			if ( len_read < data_compressed_len ) {
				ERROR("Couldn't read entire CBF data\n");
				free(buf);
				free(data_compressed);
				fclose(fh);
				return NULL;
			}

			nmemb_exp = (*w) * (*h);
			data_out = malloc(nmemb_exp*sizeof(float));
			if ( data_out == NULL ) {
				ERROR("Failed to allocate memory for CBF data\n");
				free(buf);
				free(data_compressed);
				fclose(fh);
				return NULL;
			}

			switch ( data_conversion ) {

				case CBF_NO_CONVERSION:
				r = convert_type(data_out, nmemb_exp, data_type,
				                 data_compressed,
				                 data_compressed_len);
				break;

				case CBF_BYTE_OFFSET:
				decode_cbf_byte_offset(data_out, nmemb_exp,
				                       data_compressed,
				                       data_compressed_len);
				break;

				case CBF_PACKED:
				case CBF_CANONICAL:
				ERROR("Don't yet know how to decompress "
				      "CBF_PACKED or CBF_CANONICAL\n");
				free(buf);
				free(data_compressed);
				fclose(fh);
				return NULL;

			}

			free(data_compressed);

			if ( r ) {
				free(buf);
				free(data_out);
				fclose(fh);
				return NULL;
			}

			free(buf);
			fclose(fh);
			return data_out;

		}

	} while ( rval != NULL );

	ERROR("Reached end of CBF file before finding data.\n");
	free(buf);  /* might be NULL */
	return NULL;
}


static int read_cbf(struct imagefile *f, struct image *image)
{
	float *data;
	int w, h;

	data = read_cbf_data(f->filename, f->type == IMAGEFILE_CBF,
	                     &w, &h);
	if ( data == NULL ) {
		ERROR("Failed to read CBF data\n");
		return 1;
	}

	unpack_panels(image, data, w, h);
	free(data);

	if ( image->beam != NULL ) {
		cbf_fill_in_beam_parameters(image->beam, f, image);
		if ( image->lambda > 1000 ) {
			ERROR("WARNING: Missing or nonsensical wavelength "
			      "(%e m) for %s.\n",
			      image->lambda, image->filename);
		}
	}
	cbf_fill_in_clen(image->det, f);
	fill_in_adu(image);

	return 0;
}


static int read_cbf_simple(struct imagefile *f, struct image *image)
{
	float *data;
	int w, h;

	data = read_cbf_data(f->filename, f->type == IMAGEFILE_CBF,
	                     &w, &h);
	if ( data == NULL ) {
		ERROR("Failed to read CBF data\n");
		return 1;
	}

	image->det = simple_geometry(image, w, h);
	image->dp = malloc(sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate dp array\n");
		return 1;
	}

	image->dp[0] = data;

	if ( image->beam != NULL ) {
		cbf_fill_in_beam_parameters(image->beam, f, image);
		if ( image->lambda > 1000 ) {
			ERROR("WARNING: Missing or nonsensical wavelength "
			      "(%e m) for %s.\n",
			      image->lambda, image->filename);
		}
	}
	cbf_fill_in_clen(image->det, f);
	fill_in_adu(image);

	return 0;
}


/****************************** Image files ***********************************/


signed int is_cbf_file(const char *filename)
{
	FILE *fh;
	char line[1024];

	fh = fopen(filename, "r");
	if ( fh == NULL ) return -1;

	if ( fgets(line, 1024, fh) == NULL ) return -1;
	fclose(fh);

	if ( strstr(line, "CBF") == NULL ) {
		return 0;
	}

	return 1;
}


signed int is_cbfgz_file(const char *filename)
{
	gzFile gzfh;
	char line[1024];

	gzfh = gzopen(filename, "rb");
	if ( gzfh == NULL ) return -1;
	if ( gzgets(gzfh, line, 1024) == NULL ) return -1;
	gzclose(gzfh);

	if ( strstr(line, "CBF") == NULL ) {
		return 0;
	}

	return 1;
}


struct imagefile *imagefile_open(const char *filename)
{
	struct imagefile *f;

	f = malloc(sizeof(struct imagefile));
	if ( f == NULL ) return NULL;

	if ( H5Fis_hdf5(filename) > 0 ) {

		/* This is an HDF5, pass through to HDF5 layer */
		f->type = IMAGEFILE_HDF5;
		f->hdfile = hdfile_open(filename);

		if ( f->hdfile == NULL ) {
			free(f);
			return NULL;
		}

	} else if ( is_cbf_file(filename) > 0 ) {

		f->type = IMAGEFILE_CBF;

	} else if ( is_cbfgz_file(filename) ) {

		f->type = IMAGEFILE_CBFGZ;

	} else {

		ERROR("Unrecognised file type: %s\n", filename);
		return NULL;

	}

	f->filename = strdup(filename);
	return f;
}


int imagefile_read(struct imagefile *f, struct image *image,
                   struct event *event)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		return hdf5_read2(f->hdfile, image, event, 0);
	} else if ( f->type == IMAGEFILE_CBF ) {
		return read_cbf(f, image);
	} else if ( f->type == IMAGEFILE_CBFGZ ) {
		return read_cbf(f, image);
	} else {
		ERROR("Unknown file type %i\n", f->type);
		return 1;
	}
}


/* Read a simple file, no multi-event, no prior geometry etc, and
 * generate a geometry for it */
int imagefile_read_simple(struct imagefile *f, struct image *image)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		return hdf5_read(f->hdfile, image, NULL, 0);
	} else if ( f->type == IMAGEFILE_CBF ) {
		return read_cbf_simple(f, image);
	} else if ( f->type == IMAGEFILE_CBFGZ ) {
		return read_cbf_simple(f, image);
	} else {
		ERROR("Unknown file type %i\n", f->type);
		return 1;
	}
}


enum imagefile_type imagefile_get_type(struct imagefile *f)
{
	assert(f != NULL);
	return f->type;
}


struct hdfile *imagefile_get_hdfile(struct imagefile *f)
{
	if ( f == NULL ) return NULL;

	if ( f->type != IMAGEFILE_HDF5 ) {
		ERROR("Not an HDF5 file!\n");
		return NULL;
	}

	return f->hdfile;
}


void imagefile_copy_fields(struct imagefile *f,
                           const struct imagefile_field_list *copyme,
                           FILE *fh, struct event *ev)
{
	int i;

	if ( copyme == NULL ) return;

	for ( i=0; i<copyme->n_fields; i++ ) {

		char *val;
		char *field;

		field = copyme->fields[i];

		if ( f->type == IMAGEFILE_HDF5 ) {
			val = hdfile_get_string_value(f->hdfile, field, ev);
			if ( field[0] == '/' ) {
				fprintf(fh, "hdf5%s = %s\n", field, val);
			} else {
				fprintf(fh, "hdf5/%s = %s\n", field, val);
			}
			free(val);

		} else {
			STATUS("Mock CBF variable\n");
			fprintf(fh, "cbf/%s = %s\n", field, "(FIXME)");
		}


	}
}


void imagefile_close(struct imagefile *f)
{
	if ( f->type == IMAGEFILE_HDF5 ) {
		hdfile_close(f->hdfile);
	}
	free(f->filename);
	free(f);
}


/************************** New API (DataTemplate) ****************************/

static struct image *image_new()
{
	struct image *image;

	image = malloc(sizeof(struct image));
	if ( image == NULL ) return NULL;

	image->dp = NULL;
	image->bad = NULL;
	image->sat = NULL;
	image->hit = 0;
	image->crystals = NULL;
	image->n_crystals = 0;
	image->indexed_by = INDEXING_NONE;
	image->detgeom = NULL;
	image->filename = NULL;
	image->ev = NULL;
	image->copyme = NULL;
	image->stuff_from_stream = NULL;
	image->avg_clen = -1.0;
	image->id = 0;
	image->serial = 0;
	image->spectrum = NULL;
	image->lambda = -1.0;
	image->div = -1.0;
	image->bw = -1.0;
	image->peak_resolution = -1.0;
	image->features = NULL;

	/* Deprecated stuff */
	image->beam = NULL;
	image->event = NULL;
	image->det = NULL;

	return image;
}


static int unpack_panels_dtempl(struct image *image,
                                DataTemplate *dtempl,
                                float *data,
                                int data_width, int data_height)
{
	int pi;

	image->dp = malloc(dtempl->n_panels * sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate panels.\n");
		return 1;
	}

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int fs, ss;
		int p_w, p_h;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;
		image->dp[pi] = malloc(p_w*p_h*sizeof(float));
		if ( image->dp[pi] == NULL ) {
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		if ( (p->orig_min_fs + p_w > data_width)
		  || (p->orig_min_ss + p_h > data_height) )
		{
			ERROR("Panel %s is outside range of data in CBF file\n",
			      p->name);
			return 1;
		}

		for ( ss=0; ss<p_h; ss++ ) {
		for ( fs=0; fs<p_w; fs++ ) {

			int idx;
			int cfs, css;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p_w*ss] = data[idx];

		}
		}

	}

	return 0;
}


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
	int ndims;
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

	if ( !check_path_existence(fh, panel_full_path) ) {
		ERROR("Cannot find data for panel %s (%s)\n",
		      p->name, panel_full_path);
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}

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

	f_offset = malloc(p->dim_structure->num_dims*sizeof(hsize_t));
	f_count = malloc(p->dim_structure->num_dims*sizeof(hsize_t));
	if ( (f_offset == NULL) || (f_count == NULL ) ) {
		ERROR("Failed to allocate offset or count.\n");
		free_event(ev);
		H5Fclose(fh);
		return 1;
	}
	for ( hsi=0; hsi<p->dim_structure->num_dims; hsi++ ) {

		if ( p->dim_structure->dims[hsi] == HYSL_FS ) {
			f_offset[hsi] = p->orig_min_fs;
			f_count[hsi] = p->orig_max_fs - p->orig_min_fs+1;
		} else if ( p->dim_structure->dims[hsi] == HYSL_SS ) {
			f_offset[hsi] = p->orig_min_ss;
			f_count[hsi] = p->orig_max_ss - p->orig_min_ss+1;
		} else if (p->dim_structure->dims[hsi] == HYSL_PLACEHOLDER ) {
			if ( !skip_placeholders ) {
				f_offset[hsi] = ev->dim_entries[0];
				f_count[hsi] = 1;
			}
		} else {
			f_offset[hsi] = p->dim_structure->dims[hsi];
			f_count[hsi] = 1;
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


static struct image *image_read_hdf5(DataTemplate *dtempl,
                                     const char *filename,
                                     const char *event)
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


static int load_mask_cbf(struct panel_template *p,
                         const char *filename, const char *event,
                         int gz, int *bad, int mask_good, int mask_bad)
{
	ERROR("Mask loading from CBF not yet supported\n");
	return 1;
}


/* Load bad pixels for this panel from given filename/event, and merge
 * with (already allocated/initialised) mask "bad" */
static int load_mask_hdf5(struct panel_template *p,
                          const char *filename, const char *event,
                          int *bad, int mask_good, int mask_bad)
{
	int p_w, p_h;
	int *mask;
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


struct image *image_read_cbf(DataTemplate *dtempl, const char *filename,
                             const char *event, int gz)
{
	struct image *image;
	float *data;
	int w, h;

	if ( access(filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n", filename);
		return NULL;
	}

	image = image_new();
	if ( image == NULL ) {
		ERROR("Couldn't allocate image structure.\n");
		return NULL;
	}

	data = read_cbf_data(filename, gz, &w, &h);
	if ( data == NULL ) {
		ERROR("Failed to read CBF data\n");
		return NULL;
	}

	unpack_panels_dtempl(image, dtempl, data, w, h);
	free(data);

	//cbf_fill_in_beam_parameters(image->beam, f, image);
	//cbf_fill_in_clen(image->det, f);
	//fill_in_adu(image);

	return image;
}


static double get_value_hdf5(const char *name, const char *filename,
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
	int check_pe;
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

	check_pe = check_path_existence(fh, subst_name);
	if ( check_pe == 0 ) {
		ERROR("No such event-based numeric field '%s'\n", subst_name);
		return NAN;
	}

	dh = H5Dopen2(fh, subst_name, H5P_DEFAULT);
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


static double get_value(struct image *image, const char *from)
{
	double val;
	char *rval;

	if ( from == NULL ) return NAN;

	val = strtod(from, &rval);
	if ( *rval == '\0' ) return val;

	if ( H5Fis_hdf5(image->filename) > 0 ) {
		return get_value_hdf5(from, image->filename, image->ev);

	} else if ( is_cbf_file(image->filename) > 0 ) {
		/* FIXME: From headers */
		return NAN;

	} else if ( is_cbfgz_file(image->filename) ) {
		/* FIXME: From headers */
		return NAN;

	} else {
		ERROR("Unrecognised file type: %s\n", image->filename);
		return NAN;
	}
}


static void create_detgeom(struct image *image, DataTemplate *dtempl)
{
	struct detgeom *detgeom;
	int i;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return;
	}

	detgeom = malloc(sizeof(struct detgeom));
	if ( detgeom == NULL ) return;

	detgeom->panels = malloc(dtempl->n_panels*sizeof(struct detgeom_panel));
	if ( detgeom->panels == NULL ) return;

	detgeom->n_panels = dtempl->n_panels;

	for ( i=0; i<dtempl->n_panels; i++ ) {

		detgeom->panels[i].name = safe_strdup(dtempl->panels[i].name);

		detgeom->panels[i].cnx = dtempl->panels[i].cnx;
		detgeom->panels[i].cny = dtempl->panels[i].cny;
		detgeom->panels[i].cnz = get_value(image, dtempl->panels[i].cnz_from)
		                                + dtempl->panels[i].cnz_offset;

		detgeom->panels[i].pixel_pitch = dtempl->panels[i].pixel_pitch;
		detgeom->panels[i].max_adu = dtempl->panels[i].max_adu;
		detgeom->panels[i].adu_per_photon = 1.0;  /* FIXME ! */

		detgeom->panels[i].w = dtempl->panels[i].orig_max_fs
		                        - dtempl->panels[i].orig_min_fs + 1;
		detgeom->panels[i].h = dtempl->panels[i].orig_max_ss
		                        - dtempl->panels[i].orig_min_ss + 1;

		detgeom->panels[i].fsx = dtempl->panels[i].fsx;
		detgeom->panels[i].fsy = dtempl->panels[i].fsy;
		detgeom->panels[i].fsz = dtempl->panels[i].fsz;
		detgeom->panels[i].ssx = dtempl->panels[i].ssx;
		detgeom->panels[i].ssy = dtempl->panels[i].ssy;
		detgeom->panels[i].ssz = dtempl->panels[i].ssz;

	}

	image->lambda = get_value(image, dtempl->wavelength_from);
	image->detgeom = detgeom;
	/* FIXME: spectrum */
}


/* Return non-zero if pixel fs,ss on panel p is in a bad region
 * as specified in the geometry file (regions only, not including
 * masks, NaN/inf, no_index etc */
static int in_bad_region_dtempl(DataTemplate *dtempl,
                                struct panel_template *p,
                                double fs, double ss)
{
	double rx, ry;
	double xs, ys;
	int i;

	/* Convert xs and ys, which are in fast scan/slow scan coordinates,
	 * to x and y */
	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	rx = xs + p->cnx;
	ry = ys + p->cny;

	for ( i=0; i<dtempl->n_bad; i++ ) {

		struct dt_badregion *b = &dtempl->bad[i];

		if ( (b->panel != NULL)
		  && (strcmp(b->panel, p->name) != 0) ) continue;

		if ( b->is_fsss ) {

			int nfs, nss;

			/* fs/ss bad regions are specified according
			 * to the original coordinates */
			nfs = fs + p->orig_min_fs;
			nss = ss + p->orig_min_ss;

			if ( nfs < b->min_fs ) continue;
			if ( nfs > b->max_fs ) continue;
			if ( nss < b->min_ss ) continue;
			if ( nss > b->max_ss ) continue;

		} else {

			if ( rx < b->min_x ) continue;
			if ( rx > b->max_x ) continue;
			if ( ry < b->min_y ) continue;
			if ( ry > b->max_y ) continue;

		}

		return 1;
	}

	return 0;
}


struct image *image_read(DataTemplate *dtempl, const char *filename,
                         const char *event)
{
	struct image *image;
	int i;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	if ( H5Fis_hdf5(filename) > 0 ) {
		image = image_read_hdf5(dtempl, filename, event);

	} else if ( is_cbf_file(filename) > 0 ) {
		image = image_read_cbf(dtempl, filename, event, 0);

	} else if ( is_cbfgz_file(filename) ) {
		image = image_read_cbf(dtempl, filename, event, 1);

	} else {
		ERROR("Unrecognised file type: %s\n", filename);
		return NULL;
	}

	if ( image == NULL ) return NULL;

	create_detgeom(image, dtempl);

	image->bad = malloc(dtempl->n_panels * sizeof(int *));
	if ( image->bad == NULL ) {
		ERROR("Failed to allocate bad pixel mask\n");
		return NULL;
	}

	for ( i=0; i<dtempl->n_panels; i++ ) {

		const char *mask_fn;
		int p_w, p_h;
		struct panel_template *p = &dtempl->panels[i];

		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->bad[i] = calloc(p_w*p_h, sizeof(int));
		if ( image->bad[i] == NULL ) {
			ERROR("Failed to allocate bad pixel mask\n");
			return NULL;
		}

		/* Panel marked as bad? */
		if ( p->bad ) {
			/* NB this sets every element to 0x1111,
			 * but that's OK - value is still 'true'. */
			memset(image->bad[i], 1, p_w*p_h);
		}

		/* Add bad regions (skip if panel is bad anyway) */
		if ( !p->bad ) {
			int fs, ss;
			for ( fs=0; fs<p_w; fs++ ) {
			for ( ss=0; ss<p_h; ss++ ) {
				if ( in_bad_region_dtempl(dtempl, p, fs, ss)
				     || isnan(image->dp[i][fs+ss*p_w])
				     || isinf(image->dp[i][fs+ss*p_w]) )
				{
					image->bad[i][fs+ss*p_w] = 1;
				}
			}
			}
		}

		/* Load mask (skip if panel is bad anyway) */
		if ( (!p->bad) && (p->mask != NULL) ) {
			if ( p->mask_file == NULL ) {
				mask_fn = filename;
			} else {
				mask_fn = p->mask_file;
			}
			if ( H5Fis_hdf5(mask_fn) > 0 ) {
				load_mask_hdf5(p, mask_fn, event,
				               image->bad[i],
				               dtempl->mask_good,
				               dtempl->mask_bad);

			} else if ( is_cbf_file(filename) > 0 ) {
				load_mask_cbf(p, mask_fn, event,
				              0, image->bad[i],
				              dtempl->mask_good,
				              dtempl->mask_bad);

			} else if ( is_cbfgz_file(filename) ) {
				load_mask_cbf(p, mask_fn, event,
				              1, image->bad[i],
				              dtempl->mask_good,
				              dtempl->mask_bad);

			} else {
				ERROR("Unrecognised mask file type"
				      " (%s)\n", filename);
				return NULL;
			}
		}
	}

	/* FIXME: Load saturation map */

	return image;
}


void image_free(struct image *image)
{
	int i, np;

	if ( image == NULL ) return;
	image_feature_list_free(image->features);
	free_all_crystals(image);
	free(image->filename);
	free(image->ev);

	if ( image->detgeom != NULL ) {
		np = image->detgeom->n_panels;
	} else if ( image->det != NULL ) {
		np = image->det->n_panels;
	} else {
		np = 0;
	}

	for ( i=0; i<np; i++ ) {
		if ( image->dp != NULL ) free(image->dp[i]);
		if ( image->sat != NULL ) free(image->sat[i]);
		if ( image->bad != NULL ) free(image->bad[i]);
	}

	free(image->dp);
	free(image->sat);
	free(image->bad);

	free(image);
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


ImageFeatureList *get_peaks_cxi_dtempl(const DataTemplate *dtempl,
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

	subst_name = retrieve_full_path(ev, dtempl->peak_list);
	free_event(ev);
	if ( subst_name == NULL ) {
		ERROR("Invalid peak path %s\n", subst_name);
		free_event(ev);
		H5Fclose(fh);
		return NULL;
	}

	if ( check_path_existence(fh, subst_name) == 0 ) {
		ERROR("Peak path not found: %s:%s",
		      filename, subst_name);
		free(subst_name);
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

		pn = data_template_find_panel(dtempl, fs, ss);
		if ( pn < -1 ) {
			ERROR("Peak not in panel!\n");
			continue;
		}

		/* Convert coordinates to panel-relative */
		data_template_file_to_panel_coords(dtempl, &fs, &ss);

		image_add_feature(features, fs, ss, pn,
		                  NULL, val, NULL);

	}

	return NULL;
}


ImageFeatureList *get_peaks_hdf5_dtempl(const DataTemplate *dtempl,
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

	if ( check_path_existence(fh, subst_name) == 0 ) {
		ERROR("Peak path not found: %s:%s",
		      filename, subst_name);
		free(subst_name);
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

		pn = data_template_find_panel(dtempl, fs, ss);
		if ( pn < -1 ) {
			ERROR("Peak not in panel!\n");
			continue;
		}

		/* Convert coordinates to panel-relative */
		data_template_file_to_panel_coords(dtempl, &fs, &ss);

		image_add_feature(features, fs, ss, pn,
		                  NULL, val, NULL);

	}

	free(buf);
	H5Dclose(dh);
	H5Fclose(fh);

	return features;
}


ImageFeatureList *image_read_peaks(const DataTemplate *dtempl,
                                   const char *filename,
                                   const char *event,
                                   int half_pixel_shift)
{
	if ( H5Fis_hdf5(filename) > 0 ) {

		const char *ext;
		ext = filename_extension(filename, NULL);
		if ( strcmp(ext, ".cxi") == 0 ) {
			return get_peaks_cxi_dtempl(dtempl,
			                            filename,
			                            event,
			                            half_pixel_shift);

		} else {
			return get_peaks_hdf5_dtempl(dtempl,
			                             filename,
			                             event,
			                             half_pixel_shift);
		}

	} else  {
		ERROR("Peak lists can only be read from HDF5 files\n");
		return NULL;
	}
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
	htri_t check;
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

	check = check_path_existence(pp->fh, truncated_path);
	if ( check == 0 ) {
			pop_path_entry_from_event(pp->curr_event);
			return 0;
	} else {

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


struct event_list *image_expand_frames(const DataTemplate *dtempl,
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
