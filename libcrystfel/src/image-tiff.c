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


int image_tiff_read(struct image *image, const DataTemplate *dtempl)
{
	TIFF *fh;
	float *data;
	uint32_t w, h ;
	uint16_t sampleformat, bps;

	if ( image->data_block != NULL ) {
		ERROR("In-memory TIFF not (yet!) implemented.\n");
		return 1;
	}

	if ( access(image->filename, R_OK) == -1 ) {
		ERROR("File does not exist or cannot be read: %s\n",
		      image->filename);
		return 1;
	}

	fh = TIFFOpen(image->filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open TIFF file: %s\n", image->filename);
		return 1;
	}

	TIFFPrintDirectory(fh, stdout, TIFFPRINT_STRIPS);

	TIFFGetField(fh, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(fh, TIFFTAG_IMAGELENGTH, &h);
        TIFFGetField(fh, TIFFTAG_SAMPLEFORMAT, &sampleformat); //SAMPLEFORMAT_IEEEFP);
        TIFFGetField(fh, TIFFTAG_BITSPERSAMPLE, &bps);
        //TIFFGetField(fh, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(fh, 0));

	printf("got %i x %i,  %i bits per sample, format=%i\n", w, h, bps, sampleformat);

	TIFFClose(fh);

	if ( unpack_panels(image, dtempl, data, w, h) ) {
		ERROR("Failed to read TIFF data\n");
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
