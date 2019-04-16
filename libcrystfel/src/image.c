/*
 * image.c
 *
 * Handle images and image features
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
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

/**
 * SECTION:image
 * @short_description: Data structure representing an image
 * @title: Image
 * @section_id:
 * @see_also:
 * @include: "image.h"
 * @Image:
 *
 * The <structname>image</structname> structure represents an image, usually one
 * frame from a large series of diffraction patterns, which might be from the
 * same or different crystals.
 */


struct imagefile
{
	enum imagefile_type type;
	char *filename;
	struct hdfile *hdfile;
};


struct _imagefeaturelist
{
	struct imagefeature	*features;
	int			n_features;
};


void image_add_feature(ImageFeatureList *flist, double fs, double ss,
                       struct panel *p,
                       struct image *parent, double intensity, const char *name)
{
	if ( flist->features ) {
		flist->features = realloc(flist->features,
		                    (flist->n_features+1)
		                    *sizeof(struct imagefeature));
	} else {
		assert(flist->n_features == 0);
		flist->features = malloc(sizeof(struct imagefeature));
	}

	flist->features[flist->n_features].fs = fs;
	flist->features[flist->n_features].ss = ss;
	flist->features[flist->n_features].p = p;
	flist->features[flist->n_features].intensity = intensity;
	flist->features[flist->n_features].parent = parent;
	flist->features[flist->n_features].name = name;
	flist->features[flist->n_features].valid = 1;

	flist->n_features++;

}


ImageFeatureList *image_feature_list_new()
{
	ImageFeatureList *flist;

	flist = malloc(sizeof(ImageFeatureList));

	flist->n_features = 0;
	flist->features = NULL;

	return flist;
}


static int comp(const void *a, const void *b)
{
	const struct imagefeature *ap = a;
	const struct imagefeature *bp = b;

	return ap->intensity < bp->intensity;
}


/* Strongest first.  Returned list is guaranteed not to have any holes
 * (feature->valid = 0) */
ImageFeatureList *sort_peaks(ImageFeatureList *flist)
{
	ImageFeatureList *n = image_feature_list_new();
	int nf, i;

	if ( n == NULL ) return NULL;

	n->features = malloc(flist->n_features*sizeof(struct imagefeature));
	if ( n->features == NULL ) {
		free(n);
		return NULL;
	}

	nf = 0;
	for ( i=0; i<flist->n_features; i++ ) {
		struct imagefeature *f;
		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;
		n->features[nf++] = flist->features[i];
	}
	n->n_features = nf;

	qsort(n->features, nf, sizeof(struct imagefeature), comp);

	return n;
}


void image_feature_list_free(ImageFeatureList *flist)
{
	if ( !flist ) return;

	if ( flist->features ) free(flist->features);
	free(flist);
}


struct imagefeature *image_feature_closest(ImageFeatureList *flist,
                                           double fs, double ss,
                                           struct panel *p, double *d, int *idx)
{
	int i;
	double dmin = +HUGE_VAL;
	int closest = 0;

	for ( i=0; i<flist->n_features; i++ ) {

		double ds;

		if ( p != flist->features[i].p ) continue;

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


struct imagefeature *image_get_feature(ImageFeatureList *flist, int idx)
{
	/* Sanity check */
	if ( flist == NULL ) return NULL;
	if ( idx >= flist->n_features ) return NULL;

	if ( flist->features[idx].valid == 0 ) return NULL;

	return &flist->features[idx];
}


void image_remove_feature(ImageFeatureList *flist, int idx)
{
	flist->features[idx].valid = 0;
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


static float *read_cbf_data(struct imagefile *f, int *w, int *h)
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

	if ( f->type == IMAGEFILE_CBF ) {

		fh = fopen(f->filename, "rb");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", f->filename);
			return NULL;
		}

	} else if ( f->type == IMAGEFILE_CBFGZ ) {

		gzFile gzfh;
		size_t len, len_read;
		const size_t bufinc = 8*1024*1024;  /* Allocate buffer in 8Mb chunks */
		size_t bufsz = bufinc;

		gzfh = gzopen(f->filename, "rb");
		if ( gzfh == NULL ) return NULL;

		/* Set larger buffer size for hopefully faster uncompression */
		gzbuffer(gzfh, 128*1024);

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

	} else {
		/* Don't know how we ended up here */
		return NULL;
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

	data = read_cbf_data(f, &w, &h);
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

	data = read_cbf_data(f, &w, &h);
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
