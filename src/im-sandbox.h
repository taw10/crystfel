/*
 * im-sandbox.h
 *
 * Sandbox for indexing
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
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

enum {
	PEAK_ZAEF,
	PEAK_HDF5,
};


/* Information about the indexing process which is common to all patterns */
struct index_args
{
	UnitCell *cell;
	int config_cmfilter;
	int config_noisefilter;
	int config_verbose;
	int stream_flags;         /* What goes into the output? */
	int config_satcorr;
	int config_closer;
	int config_insane;
	int config_bgsub;
	float threshold;
	float min_gradient;
	float min_snr;
	double min_int_snr;
	struct detector *det;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int peaks;                /* Peak detection method */
	int cellr;
	float tols[4];
	struct beam_params *beam;
	char *element;
	char *hdf5_peak_path;
	double ir_inn;
	double ir_mid;
	double ir_out;
	struct copy_hdf5_field *copyme;
	int integrate_saturated;
	int use_saturated;
	int no_revalidate;
	int integrate_found;
};




/* Information about the indexing process for one pattern */
struct pattern_args
{
	/* "Input" */
	char *filename;

	/* "Output" */
	int indexable;
};

extern void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                           int config_basename, FILE *fh,
                           char *use_this_one_instead, FILE *ofh);
