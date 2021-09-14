/*
 * image-cbf.c
 *
 * Image loading, CBF parts
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>

#include "image.h"
#include "utils.h"
#include "detgeom.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"

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
	char *buf = NULL;
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
		int len_read;
		size_t len;
		const size_t bufinc = 8*1024*1024;  /* Allocate buffer in 8Mb chunks */
		size_t bufsz = bufinc;

		gzfh = gzopen(filename, "rb");
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
				buf = srealloc(buf, bufsz);
				if ( buf == NULL ) return NULL;
			}

		} while ( len_read == bufinc );

		fh = fmemopen(buf, len, "rb");
		if ( fh == NULL ) {
			free(buf);
			return NULL;
		}

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


signed int is_cbf_file(const char *filename)
{
	FILE *fh;
	char line[1024];

	fh = fopen(filename, "r");
	if ( fh == NULL ) return -1;

	if ( fgets(line, 1024, fh) == NULL ) {
		fclose(fh);
		return 0;
	}

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
	if ( gzgets(gzfh, line, 1024) == NULL ) return 0;
	gzclose(gzfh);

	if ( strstr(line, "CBF") == NULL ) {
		return 0;
	}

	return 1;
}


int image_cbf_read_mask(struct panel_template *p,
                        const char *filename, const char *event,
                        int gz, int *bad, int mask_good, int mask_bad)
{
	ERROR("Mask loading from CBF not yet supported\n");
	return 1;
}


static int unpack_panels(struct image *image,
                         const DataTemplate *dtempl,
                         float *data, int data_width, int data_height)
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


int image_cbf_read(struct image *image,
                   const DataTemplate *dtempl,
                   int gz)
{
	float *data;
	int w, h;

	if ( image->data_block != NULL ) {
		ERROR("In-memory CBF not (yet!) implemented.\n");
		return 1;
	}

	if ( access(image->filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      image->filename);
		return 1;
	}

	data = read_cbf_data(image->filename, gz, &w, &h);
	if ( data == NULL ) {
		ERROR("Failed to read CBF data\n");
		return 1;
	}

	unpack_panels(image, dtempl, data, w, h);
	free(data);

	//cbf_fill_in_beam_parameters(image->beam, f, image);
	//cbf_fill_in_clen(image->det, f);
	//fill_in_adu(image);

	return 0;
}


int image_cbf_read_header_to_cache(struct image *image,
                                   const char *from)
{
	/* FIXME: Implementation (GitLab #10) */
	ERROR("Reading headers from CBF files is not currently supported.\n");
	return 1;
}
