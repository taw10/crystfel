/*
 * im-argparse.h
 *
 * Command line argument parsing for indexamajig
 *
 * Copyright © 2023-2024 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023-2024 Thomas White <taw@physics.org>
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

#ifndef IM_ARGPARSE_H
#define IM_ARGPARSE_H

struct indexamajig_arguments;

#include <index.h>

#include "process_image.h"

struct indexamajig_arguments
{
	struct index_args iargs;  /* These are the options that will be
	                           * given to process_image */
	char *filename;
	char *geom_filename;
	char *outfile;
	char *prefix;
	int check_prefix;
	int n_proc;
	char *cellfile;
	char *indm_str;
	int basename;
	struct im_zmq_params zmq_params;
	struct im_asapo_params asapo_params;
	int serial_start;
	char *temp_location;
	int if_refine;
	int if_checkcell;
	int if_peaks;
	int if_multi;
	int if_retry;
	int profile;  /* Whether to do wall-clock time profiling */
	int no_data_timeout;
	char **copy_headers;
	int n_copy_headers;
	char *harvest_file;
	char *milledir;
	char *millefile;
	int cpu_pin;
	int worker;
	int worker_id;
	char *worker_tmpdir;
	int fd_stream;
	int fd_mille;
	char *queue_sem;
	char *shm_name;

	struct taketwo_options **taketwo_opts_ptr;
	struct felix_options **felix_opts_ptr;
	struct xgandalf_options **xgandalf_opts_ptr;
	struct pinkindexer_options **pinkindexer_opts_ptr;
	struct fromfile_options **fromfile_opts_ptr;
	struct asdf_options **asdf_opts_ptr;
};

extern struct indexamajig_arguments *parse_indexamajig_args(int argc, char *argv[]);

#endif /* IM_ARGPARSE_H */
