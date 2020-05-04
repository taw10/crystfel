/*
 * gui_project.c
 *
 * GUI project persistence
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "crystfel_gui.h"
#include "gui_backend_local.h"
#include "gui_project.h"

static double parse_float(const char *val)
{
	float v;

	if (sscanf(val, "%f", &v) != 1) {
		ERROR("Invalid float value '%s'\n", val);
		return NAN;
	}

	return v;
}


static int parse_int(const char *val)
{
	int v;

	if (sscanf(val, "%i", &v) != 1) {
		ERROR("Invalid int value '%s'\n", val);
		return 0;
	}

	return v;
}


static struct crystfel_backend *parse_backend(const char *val)
{
	if ( strcmp(val, "local") == 0 ) {
		return backend_local;
	}

	ERROR("Invalid backend '%s'\n", val);
	return NULL;
}


enum match_type_id decode_matchtype(const char *type_id)
{
	if ( strcmp(type_id, "everything") == 0 ) return MATCH_EVERYTHING;
	if ( strcmp(type_id, "lcls-cheetah-hdf5") == 0 ) return MATCH_CHEETAH_LCLS_H5;
	if ( strcmp(type_id, "cheetah-cxi") == 0 ) return MATCH_CHEETAH_CXI;
	if ( strcmp(type_id, "cbf") == 0 ) return MATCH_CBF;
	if ( strcmp(type_id, "cbfgz") == 0 ) return MATCH_CBFGZ;
	ERROR("Unknown match type id '%s'\n", type_id);
	return MATCH_EVERYTHING;
}


int match_filename(const char *fn, enum match_type_id mt)
{
	const char *ext = NULL;
	const char *ext2 = NULL;

	if ( mt == MATCH_EVERYTHING ) return 1;

	ext = filename_extension(fn, &ext2);
	if ( ext == NULL ) return 0;

	if ( mt == MATCH_CHEETAH_LCLS_H5 ) {
		return ((strcmp(ext, ".h5")==0)
		        && (strncmp(fn, "LCLS", 4)==0));
	}
	if ( mt == MATCH_CHEETAH_CXI ) return strcmp(ext, ".cxi")==0;
	if ( mt == MATCH_CBF ) return strcmp(ext, ".cbf")==0;
	if ( mt == MATCH_CBFGZ  ) {
		if ( ext2 != NULL ) {
			return strcmp(ext2, ".cbf.gz")==0;
		}
	}

	return 0;
}


static void handle_var(const char *key, const char *val,
                       struct crystfelproject *proj)
{
	if ( strcmp(key, "peak_search_params.method") == 0 ) {
		proj->peak_search_params.method = parse_peaksearch(val);
		if ( proj->peak_search_params.method == PEAK_ERROR ) {
			ERROR("Unrecognised peak search method '%s'\n",
			      val);
		}
	}

	if ( strcmp(key, "peak_search_params.threshold") == 0 ) {
		proj->peak_search_params.threshold = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.min_sq_gradient") == 0 ) {
		proj->peak_search_params.min_sq_gradient = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.min_snr") == 0 ) {
		proj->peak_search_params.min_snr = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.local_bg_radius") == 0 ) {
		proj->peak_search_params.local_bg_radius = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.min_res") == 0 ) {
		proj->peak_search_params.min_res = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.min_sig") == 0 ) {
		proj->peak_search_params.min_sig = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.max_res") == 0 ) {
		proj->peak_search_params.max_res = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.min_pix_count") == 0 ) {
		proj->peak_search_params.min_pix_count = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.max_pix_count") == 0 ) {
		proj->peak_search_params.max_pix_count = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.min_peak_over_neighbour") == 0 ) {
		proj->peak_search_params.min_peak_over_neighbour = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.pk_inn") == 0 ) {
		proj->peak_search_params.pk_inn = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.pk_mid") == 0 ) {
		proj->peak_search_params.pk_mid = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.pk_out") == 0 ) {
		proj->peak_search_params.pk_out = parse_float(val);
	}

	if ( strcmp(key, "peak_search_params.half_pixel_shift") == 0 ) {
		proj->peak_search_params.half_pixel_shift = parse_int(val);
	}

	if ( strcmp(key, "peak_search_params.revalidate") == 0 ) {
		proj->peak_search_params.revalidate = parse_int(val);
	}

	if ( strcmp(key, "backend") == 0 ) {
		proj->backend = parse_backend(val);
	}

	if ( strcmp(key, "geom") == 0 ) {
		proj->geom_filename = strdup(val);
	}

	if ( strcmp(key, "data_folder") == 0 ) {
		proj->data_top_folder = strdup(val);
	}

	if ( strcmp(key, "search_pattern") == 0 ) {
		proj->data_search_pattern = decode_matchtype(val);
	}
}


void clear_project_files(struct crystfelproject *proj)
{
	int i;

	for ( i=0; i<proj->n_frames; i++ ) {
		free(proj->filenames[i]);
		free(proj->events[i]);
	}
	free(proj->filenames);
	free(proj->events);
	proj->n_frames = 0;
	proj->max_frames = 0;
	proj->filenames = NULL;
	proj->events = NULL;
}


void add_file_to_project(struct crystfelproject *proj,
                         const char *filename, const char *event)
{
	if ( proj->n_frames == proj->max_frames ) {
		int n_max = proj->max_frames + 1024;
		char **n_filenames;
		char **n_events;
		n_filenames = realloc(proj->filenames,
		                      n_max*sizeof(char *));
		n_events = realloc(proj->events,
		                   n_max*sizeof(char *));
		if ( (n_filenames == NULL) || (n_events == NULL) ) {
			ERROR("Failed to allocate new filename\n");
			return;
		}
		proj->max_frames = n_max;
		proj->filenames = n_filenames;
		proj->events = n_events;
	}

	proj->filenames[proj->n_frames] = strdup(filename);
	proj->events[proj->n_frames] = safe_strdup(event);
	proj->n_frames++;
}


int load_project(struct crystfelproject *proj)
{
	FILE *fh;
	char line[1024];
	char *rval;
	int image_list_mode = 0;

	fh = fopen("crystfel.project", "r");
	if ( fh == NULL ) return 1;

	do {

		char *sp;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		chomp(line);

		if ( line[0] == '\0' ) continue;

		if ( image_list_mode ) {
			char *ev = NULL;
			size_t n = strlen(line)-1;
			for ( ; n>0; n-- ) {
				if ( line[n] == ' ' ) {
					line[n] = '\0';
					ev = &line[n+1];
					break;
				}
			}
			add_file_to_project(proj, line, ev);
			continue;
		}

		if ( strcmp(line, "-----") == 0 ) {
			image_list_mode = 1;
			continue;
		}

		sp = strchr(line, ' ');
		if ( sp == NULL ) {
			ERROR("Unrecognised line in crystfel.project "
			      "file: '%s'\n", line);
			continue;
		}

		sp[0] = '\0';

		handle_var(line, sp+1, proj);

	} while ( rval != NULL );

	fclose(fh);

	return 0;
}
