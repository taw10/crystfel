/*
 * image-tiff.c
 *
 * Image loading, TIFF parts
 *
 * Copyright © 2025 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2025 Thomas White <taw@physics.org>
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
#include <unistd.h>

#ifdef HAVE_LIBTIFF

#include <tiffio.h>

#include "image.h"
#include "utils.h"
#include "detgeom.h"

#include "datatemplate.h"
#include "datatemplate_priv.h"


static int unpack_panels(struct image *image,
                         const DataTemplate *dtempl,
                         float *data, int data_width, int data_height)
{
	int pi;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int fs, ss;
		int p_w, p_h;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		if ( (p->orig_min_fs + p_w > data_width)
		  || (p->orig_min_ss + p_h > data_height) )
		{
			ERROR("Panel %s is outside range of data in TIFF file\n",
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
			if ( !isfinite(data[idx]) ) {
				image->bad[pi][fs+p_w*ss] = 1;
			}

		}
		}

	}

	return 0;
}


static int min(int x, int y)
{
	if ( x < y ) return x;
	return y;
}


struct sis_md
{
	uint32_t magic;
	uint32_t ignore1;
	uint16_t ignore2;
	uint16_t min;
	uint16_t hour;
	uint16_t day;
	uint16_t month;
	uint16_t year;
} __attribute__((packed));


struct sis_md2
{
	char ignore[10];
	int16_t unitExp;
	double sizeX;
	double sizeY;
	char ignore2[8];
	double mag;
	uint16_t ignore3;
	char cname[34];
	char ptype[32];
} __attribute__((packed));

#define TIFF_SISOFFSET (33560)

static const TIFFFieldInfo sis_fields[] = {
	{TIFF_SISOFFSET, 1, 1, TIFF_LONG, FIELD_CUSTOM, 0, 0, "Olympus SIS Metadata Offset"}
};


static void merge_tiff_tags(TIFF *fh)
{
	TIFFMergeFieldInfo(fh, sis_fields, 1);
}


int image_tiff_read(struct image *image, const DataTemplate *dtempl)
{
	TIFF *fh;
	float *data;
	uint32_t w, h ;
	uint16_t sampleformat, bps, spp;
	uint32_t sis_offs;
	int have_sis_metadata;

	if ( image->data_block != NULL ) {
		ERROR("In-memory TIFF not (yet!) implemented.\n");
		return 1;
	}

	if ( access(image->filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      image->filename);
		return 1;
	}

	TIFFSetTagExtender(merge_tiff_tags);

	fh = TIFFOpen(image->filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open TIFF file: %s\n", image->filename);
		return 1;
	}

	TIFFGetField(fh, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(fh, TIFFTAG_IMAGELENGTH, &h);
        if ( TIFFGetField(fh, TIFFTAG_SAMPLEFORMAT, &sampleformat) != 1 ) {
	        sampleformat = SAMPLEFORMAT_UINT;
	}
        TIFFGetField(fh, TIFFTAG_BITSPERSAMPLE, &bps);
        if ( TIFFGetField(fh, TIFFTAG_SAMPLESPERPIXEL, &spp) == 1 ) {
	        if ( spp != 1 ) {
		        ERROR("Unable to read TIFF data - more than one sample per pixel\n");
			TIFFClose(fh);
			return 1;
		}
	}

	if ( sampleformat != SAMPLEFORMAT_UINT ) {
		ERROR("Unable to read TIFF data - not unsigned int format\n");
		TIFFClose(fh);
		return 1;
	}

	if ( bps != 16 ) {
		ERROR("Unable to read TIFF data - not 16 bit format\n");
		TIFFClose(fh);
		return 1;
	}

	uint32_t *bc;
	int nstrips = TIFFNumberOfStrips(fh);
	size_t stripsize = 0;
	int s, rps;
	int y = 0;
	TIFFGetField(fh, TIFFTAG_STRIPBYTECOUNTS, &bc);
	TIFFGetField(fh, TIFFTAG_ROWSPERSTRIP, &rps);
	for ( s=0; s<nstrips; s++ ) {
		if ( bc[s] > stripsize ) stripsize = bc[s];
	}

	data = cfmalloc(w*h*sizeof(float));
	uint16_t *strip = cfmalloc(stripsize);
	if ( (data == NULL) || (strip == NULL) ) {
		ERROR("Failed to allocate memory for TIFF read\n");
		TIFFClose(fh);
		return 1;
	}

	for ( s=0; s<TIFFNumberOfStrips(fh); s++ ) {
		int i, rows_in_strip;
		TIFFReadEncodedStrip(fh, s, strip, stripsize);
		rows_in_strip = min(h-y, rps);
		for ( i=0; i<rows_in_strip; i++ ) {
			int x;
			for ( x=0; x<w; x++ ) {
				data[x+(y+i)*w] = strip[x+i*w];
			}
		}
		y += rows_in_strip;
	}

	cffree(strip);

	have_sis_metadata = TIFFGetField(fh, TIFF_SISOFFSET, &sis_offs);

	TIFFClose(fh);

	if ( have_sis_metadata ) {

		FILE *ifh = fopen(image->filename, "rb");

		struct sis_md buf;
		fseek(ifh, sis_offs, SEEK_SET);
		fread(&buf, sizeof(buf), 1, ifh);

		if ( buf.magic != 810764627 ) {
			ERROR("Invalid SIS metadata\n");
		} else {

			uint32_t offs;
			fseek(ifh, sis_offs+64, SEEK_SET);
			fread(&offs, 4, 1, ifh);

			struct sis_md2 buf2;
			fseek(ifh, offs, SEEK_SET);
			fread(&buf2, sizeof(buf2), 1, ifh);
			printf("unitexp = %i\n", buf2.unitExp);
			printf("sizeX = %e\n", buf2.sizeX);
			printf("sizeY = %e\n", buf2.sizeY);
			printf("mag = %e\n", buf2.mag);
			printf("camname = %s\n", buf2.cname);
			printf("pictype = %s\n", buf2.ptype);

		}

		fclose(ifh);

	}

	if ( unpack_panels(image, dtempl, data, w, h) ) {
		ERROR("Failed to read TIFF data\n");
		cffree(data);
		return 1;
	}
	cffree(data);

	return 0;
}


int image_tiff_read_header_to_cache(struct image *image,
                                   const char *from)
{
	/* FIXME: Implementation (GitLab #10) */
	ERROR("Reading headers from TIFF files is not currently supported.\n");
	return 1;
}


int image_tiff_read_mask(struct panel_template *p,
                         const char *filename, const char *event,
                         int *bad, int mask_good, int mask_bad)
{
	ERROR("Mask loading from TIFF not yet supported\n");
	return 1;
}

#endif /* HAVE_LIBTIFF */
