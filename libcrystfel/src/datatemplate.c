/*
 * datatemplate.c
 *
 * Data template structure
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2020 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "utils.h"
#include "datatemplate.h"

#include "datatemplate_priv.h"


/**
 * \file datatemplate.h
 */

struct rg_definition {
	char *name;
	char *pns;
};


struct rgc_definition {
	char *name;
	char *rgs;
};


static struct panel_template *new_panel(DataTemplate *det,
                                        const char *name,
                                        struct panel_template *defaults)
{
	struct panel_template *new;

	det->n_panels++;
	det->panels = realloc(det->panels,
	                      det->n_panels*sizeof(struct panel_template));

	new = &det->panels[det->n_panels-1];
	memcpy(new, defaults, sizeof(struct panel_template));

	/* Set name */
	new->name = strdup(name);

	/* Copy strings */
	new->cnz_from = safe_strdup(defaults->cnz_from);
	new->data = safe_strdup(defaults->data);
	new->mask = safe_strdup(defaults->mask);
	new->mask_file = safe_strdup(defaults->mask_file);
	new->satmap = safe_strdup(defaults->satmap);
	new->satmap_file = safe_strdup(defaults->satmap_file);

	return new;
}


static struct dt_badregion *new_bad_region(DataTemplate *det, const char *name)
{
	struct dt_badregion *new;

	det->n_bad++;
	det->bad = realloc(det->bad, det->n_bad*sizeof(struct dt_badregion));

	new = &det->bad[det->n_bad-1];
	new->min_x = NAN;
	new->max_x = NAN;
	new->min_y = NAN;
	new->max_y = NAN;
	new->min_fs = 0;
	new->max_fs = 0;
	new->min_ss = 0;
	new->max_ss = 0;
	new->is_fsss = 99; /* Slightly nasty: means "unassigned" */
	new->panel_name = NULL;
	new->panel_number = 0;  /* Needs to be set after loading */
	strcpy(new->name, name);

	return new;
}


static struct panel_template *find_panel_by_name(DataTemplate *det,
                                                 const char *name)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {
		if ( strcmp(det->panels[i].name, name) == 0 ) {
			return &det->panels[i];
		}
	}

	return NULL;
}


static struct dt_badregion *find_bad_region_by_name(DataTemplate *det,
                                                    const char *name)
{
	int i;

	for ( i=0; i<det->n_bad; i++ ) {
		if ( strcmp(det->bad[i].name, name) == 0 ) {
			return &det->bad[i];
		}
	}

	return NULL;
}


static struct rigid_group *find_or_add_rg(DataTemplate *det,
                                          const char *name)
{
	int i;
	struct rigid_group **new;
	struct rigid_group *rg;

	for ( i=0; i<det->n_rigid_groups; i++ ) {

		if ( strcmp(det->rigid_groups[i]->name, name) == 0 ) {
			return det->rigid_groups[i];
		}

	}

	new = realloc(det->rigid_groups,
	              (1+det->n_rigid_groups)*sizeof(struct rigid_group *));
	if ( new == NULL ) return NULL;

	det->rigid_groups = new;

	rg = malloc(sizeof(struct rigid_group));
	if ( rg == NULL ) return NULL;

	det->rigid_groups[det->n_rigid_groups++] = rg;

	rg->name = strdup(name);
	rg->panel_numbers = NULL;
	rg->n_panels = 0;

	return rg;
}


static struct rg_collection *find_or_add_rg_coll(DataTemplate *det,
                                                 const char *name)
{
	int i;
	struct rg_collection **new;
	struct rg_collection *rgc;

	for ( i=0; i<det->n_rg_collections; i++ ) {
		if ( strcmp(det->rigid_group_collections[i]->name, name) == 0 )
		{
			return det->rigid_group_collections[i];
		}
	}

	new = realloc(det->rigid_group_collections,
	              (1+det->n_rg_collections)*sizeof(struct rg_collection *));
	if ( new == NULL ) return NULL;

	det->rigid_group_collections = new;

	rgc = malloc(sizeof(struct rg_collection));
	if ( rgc == NULL ) return NULL;

	det->rigid_group_collections[det->n_rg_collections++] = rgc;

	rgc->name = strdup(name);
	rgc->rigid_groups = NULL;
	rgc->n_rigid_groups = 0;

	return rgc;
}


static void add_to_rigid_group(struct rigid_group *rg, int panel_number)
{
	int *pn;

	pn = realloc(rg->panel_numbers, (1+rg->n_panels)*sizeof(int));
	if ( pn == NULL ) {
		ERROR("Couldn't add panel to rigid group.\n");
		return;
	}

	rg->panel_numbers = pn;
	rg->panel_numbers[rg->n_panels++] = panel_number;
}


static void add_to_rigid_group_coll(struct rg_collection *rgc,
                                    struct rigid_group *rg)
{
	struct rigid_group **r;

	r = realloc(rgc->rigid_groups, (1+rgc->n_rigid_groups)*
	            sizeof(struct rigid_group *));
	if ( r == NULL ) {
		ERROR("Couldn't add rigid group to collection.\n");
		return;
	}

	rgc->rigid_groups = r;
	rgc->rigid_groups[rgc->n_rigid_groups++] = rg;
}


/* Free all rigid groups in detector */
static void free_all_rigid_groups(DataTemplate *det)
{
	int i;

	if ( det->rigid_groups == NULL ) return;
	for ( i=0; i<det->n_rigid_groups; i++ ) {
		free(det->rigid_groups[i]->name);
		free(det->rigid_groups[i]->panel_numbers);
		free(det->rigid_groups[i]);
	}
	free(det->rigid_groups);
}


/* Free all rigid groups in detector */
static void free_all_rigid_group_collections(DataTemplate *det)
{
	int i;

	if ( det->rigid_group_collections == NULL ) return;
	for ( i=0; i<det->n_rg_collections; i++ ) {
		free(det->rigid_group_collections[i]->name);
		free(det->rigid_group_collections[i]->rigid_groups);
		free(det->rigid_group_collections[i]);
	}
	free(det->rigid_group_collections);
}


static struct rigid_group *find_rigid_group_by_name(DataTemplate *det,
                                                    char *name)
{
	int i;

	for ( i=0; i<det->n_rigid_groups; i++ ) {
		if ( strcmp(det->rigid_groups[i]->name, name) == 0 ) {
			return det->rigid_groups[i];
		}
	}

	return NULL;
}


static int atob(const char *a)
{
	if ( strcasecmp(a, "true") == 0 ) return 1;
	if ( strcasecmp(a, "false") == 0 ) return 0;
	return atoi(a);
}


static int assplode_algebraic(const char *a_orig, char ***pbits)
{
	int len, i;
	int nexp;
	char **bits;
	char *a;
	int idx, istr;

	len = strlen(a_orig);

	/* Add plus at start if no sign there already */
	if ( (a_orig[0] != '+') && (a_orig[0] != '-') ) {
		len += 1;
		a = malloc(len+1);
		snprintf(a, len+1, "+%s", a_orig);
		a[len] = '\0';

	} else {
		a = strdup(a_orig);
	}

	/* Count the expressions */
	nexp = 0;
	for ( i=0; i<len; i++ ) {
		if ( (a[i] == '+') || (a[i] == '-') ) nexp++;
	}

	bits = calloc(nexp, sizeof(char *));

	/* Break the string up */
	idx = -1;
	istr = 0;
	assert((a[0] == '+') || (a[0] == '-'));
	for ( i=0; i<len; i++ ) {

		char ch;

		ch = a[i];

		if ( ch == ' ' ) continue;

		if ( (ch == '+') || (ch == '-') ) {
			if ( idx >= 0 ) bits[idx][istr] = '\0';
			idx++;
			bits[idx] = malloc(len+1);
			istr = 0;
		}

		if ( !isdigit(ch) && (ch != '.') && (ch != '+') && (ch != '-')
		  && (ch != 'x') && (ch != 'y') && (ch != 'z') )
		{
			ERROR("Invalid character '%c' found.\n", ch);
			return 0;
		}

		assert(idx >= 0);
		bits[idx][istr++] = ch;

	}
	if ( idx >= 0 ) bits[idx][istr] = '\0';

	*pbits = bits;
	free(a);

	return nexp;
}


/* Parses the scan directions (accounting for possible rotation)
 * Assumes all white spaces have been already removed */
static int dir_conv(const char *a, double *sx, double *sy, double *sz)
{
	int n;
	char **bits;
	int i;

	*sx = 0.0;  *sy = 0.0;  *sz = 0.0;

	n = assplode_algebraic(a, &bits);

	if ( n == 0 ) {
		ERROR("Invalid direction '%s'\n", a);
		return 1;
	}

	for ( i=0; i<n; i++ ) {

		int len;
		double val;
		char axis;
		int j;

		len = strlen(bits[i]);
		assert(len != 0);
		axis = bits[i][len-1];
		if ( (axis != 'x') && (axis != 'y') && (axis != 'z') ) {
			ERROR("Invalid symbol '%c' - must be x, y or z.\n",
			      axis);
			return 1;
		}

		/* Chop off the symbol now it's dealt with */
		bits[i][len-1] = '\0';

		/* Check for anything that isn't part of a number */
		for ( j=0; j<strlen(bits[i]); j++ ) {
			if ( isdigit(bits[i][j]) ) continue;
			if ( bits[i][j] == '+' ) continue;
			if ( bits[i][j] == '-' ) continue;
			if ( bits[i][j] == '.' ) continue;
			ERROR("Invalid coefficient '%s'\n", bits[i]);
		}

		if ( strlen(bits[i]) == 0 ) {
			val = 1.0;
		} else {
			val = atof(bits[i]);
		}
		if ( strlen(bits[i]) == 1 ) {
			if ( bits[i][0] == '+' ) val = 1.0;
			if ( bits[i][0] == '-' ) val = -1.0;
		}
		switch ( axis ) {

			case 'x' :
			*sx += val;
			break;

			case 'y' :
			*sy += val;
			break;

			case 'z' :
			*sz += val;
			break;
		}

		free(bits[i]);

	}
	free(bits);

	return 0;
}


int set_dim(struct panel_template *panel, int dimension,
            const char *val)
{
	if ( dimension >= MAX_DIMS ) {
		ERROR("Too many dimensions!\n");
		return 1;
	}

	if ( strcmp(val, "fs") == 0 ) {
		panel->dims[dimension] = DIM_FS;
	} else if ( strcmp(val, "ss") == 0 ) {
		panel->dims[dimension] = DIM_SS;
	} else if ( strcmp(val, "%") == 0 ) {
		panel->dims[dimension] = DIM_PLACEHOLDER;
	} else {
		char *endptr;
		unsigned long int fix_val = strtoul(val, &endptr, 10);
		if ( endptr[0] != '\0' ) {
			ERROR("Invalid dimension value '%s'\n", val);
			return 1;
		} else {
			panel->dims[dimension] = fix_val;
		}
	}
	return 0;
}


static int add_flag_value(struct panel_template *p,
                          float val,
                          enum flag_value_type type)
{
	int i;

	for ( i=0; i<MAX_FLAG_VALUES; i++ ) {
		if ( p->flag_types[i] == FLAG_NOTHING ) {
			p->flag_types[i] = type;
			p->flag_values[i] = val;
			return 0;
		}
	}

	ERROR("Too many flag values.\n");
	return 1;
}


static int parse_field_for_panel(struct panel_template *panel, const char *key,
                                 const char *val, DataTemplate *det)
{
	int reject = 0;

	if ( strcmp(key, "min_fs") == 0 ) {
		panel->orig_min_fs = atof(val);
	} else if ( strcmp(key, "max_fs") == 0 ) {
		panel->orig_max_fs = atof(val);
	} else if ( strcmp(key, "min_ss") == 0 ) {
		panel->orig_min_ss = atof(val);
	} else if ( strcmp(key, "max_ss") == 0 ) {
		panel->orig_max_ss = atof(val);
	} else if ( strcmp(key, "corner_x") == 0 ) {
		panel->cnx = atof(val);
	} else if ( strcmp(key, "corner_y") == 0 ) {
		panel->cny = atof(val);
	} else if ( strcmp(key, "rail_direction") == 0 ) {
		if ( dir_conv(val, &panel->rail_x,
		                   &panel->rail_y,
		                   &panel->rail_z) )
		{
			ERROR("Invalid rail direction '%s'\n", val);
			reject = 1;
		}
	} else if ( strcmp(key, "clen_for_centering") == 0 ) {
		panel->clen_for_centering = atof(val);
	} else if ( strcmp(key, "adu_per_eV") == 0 ) {
		panel->adu_scale = atof(val);
		panel->adu_scale_unit = ADU_PER_EV;
	} else if ( strcmp(key, "adu_per_photon") == 0 ) {
		panel->adu_scale = atof(val);
		panel->adu_scale_unit = ADU_PER_PHOTON;
	} else if ( strcmp(key, "clen") == 0 ) {
		/* Gets expanded when image is loaded */
		panel->cnz_from = strdup(val);

	} else if ( strcmp(key, "data") == 0 ) {
		if ( strncmp(val,"/",1) != 0 ) {
			ERROR("Invalid data location '%s'\n", val);
			reject = -1;
		}
		free(panel->data);
		panel->data = strdup(val);

	} else if ( strcmp(key, "mask") == 0 ) {
		if ( strncmp(val,"/",1) != 0 ) {
			ERROR("Invalid mask location '%s'\n", val);
			reject = -1;
		}
		panel->mask = strdup(val);

	} else if ( strcmp(key, "mask_file") == 0 ) {
		panel->mask_file = strdup(val);

	} else if ( strcmp(key, "saturation_map") == 0 ) {
		panel->satmap = strdup(val);
	} else if ( strcmp(key, "saturation_map_file") == 0 ) {
		panel->satmap_file = strdup(val);

	} else if ( strcmp(key, "coffset") == 0) {
		panel->cnz_offset = atof(val);
	} else if ( strcmp(key, "res") == 0 ) {
		panel->pixel_pitch = 1.0/atof(val);
	} else if ( strcmp(key, "max_adu") == 0 ) {
		panel->max_adu = atof(val);
		ERROR("WARNING: It's usually better not to set max_adu "
		      "in the geometry file.  Use --max-adu during "
		      "merging instead.\n");

	} else if ( strcmp(key, "flag_equal") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_EQUAL) ) {
			reject = -1;
		}
	} else if ( strcmp(key, "flag_lessthan") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_LESSTHAN) ) {
			reject = -1;
		}
	} else if ( strcmp(key, "flag_morethan") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_MORETHAN) ) {
			reject = -1;
		}

	} else if ( strcmp(key, "badrow_direction") == 0 ) {
		ERROR("WARNING 'badrow_direction' is ignored in this version.\n");
	} else if ( strcmp(key, "no_index") == 0 ) {
		panel->bad = atob(val);
	} else if ( strcmp(key, "fs") == 0 ) {
		if ( dir_conv(val, &panel->fsx, &panel->fsy, &panel->fsz) != 0 )
		{
			ERROR("Invalid fast scan direction '%s'\n", val);
			reject = 1;
		}
	} else if ( strcmp(key, "ss") == 0 ) {
		if ( dir_conv(val, &panel->ssx, &panel->ssy, &panel->ssz) != 0 )
		{
			ERROR("Invalid slow scan direction '%s'\n", val);
			reject = 1;
		}
	} else if ( strncmp(key, "dim", 3) == 0) {
		int dim_entry;
		char *endptr;
		if ( key[3] != '\0' ) {
			dim_entry = strtoul(key+3, &endptr, 10);
			if ( endptr[0] != '\0' ) {
				ERROR("Invalid dimension number %s\n",
				      key+3);
			} else {
				if ( set_dim(panel, dim_entry, val) ) {
					ERROR("Failed to set dim structure entry\n");
				}
			}
		} else {
			ERROR("'dim' must be followed by a number, e.g. 'dim0'\n");
		}
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static int check_badr_fsss(struct dt_badregion *badr, int is_fsss)
{
	/* First assignment? */
	if ( badr->is_fsss == 99 ) {
		badr->is_fsss = is_fsss;
		return 0;
	}

	if ( is_fsss != badr->is_fsss ) {
		ERROR("You can't mix x/y and fs/ss in a bad region.\n");
		return 1;
	}

	return 0;
}


static int parse_field_bad(struct dt_badregion *badr, const char *key,
                           const char *val)
{
	int reject = 0;

	if ( strcmp(key, "min_x") == 0 ) {
		badr->min_x = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "max_x") == 0 ) {
		badr->max_x = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "min_y") == 0 ) {
		badr->min_y = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "max_y") == 0 ) {
		badr->max_y = atof(val);
		reject = check_badr_fsss(badr, 0);
	} else if ( strcmp(key, "min_fs") == 0 ) {
		badr->min_fs = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "max_fs") == 0 ) {
		badr->max_fs = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "min_ss") == 0 ) {
		badr->min_ss = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "max_ss") == 0 ) {
		badr->max_ss = atof(val);
		reject = check_badr_fsss(badr, 1);
	} else if ( strcmp(key, "panel") == 0 ) {
		badr->panel_name = strdup(val);
	} else {
		ERROR("Unrecognised field '%s'\n", key);
	}

	return reject;
}


static int parse_electron_voltage(const char *val,
                                  char **p_from,
                                  enum wavelength_unit *punit)
{
	char *valcpy;
	char *sp;

	valcpy = strdup(val);
	if ( valcpy == NULL ) return 1;

	/* "electron_voltage" directive must have explicit units */
	sp = strchr(valcpy, ' ');
	if ( sp == NULL ) {
		free(valcpy);
		return 1;
	}

	if ( strcmp(sp+1, "V") == 0 ) {
		*punit = WAVELENGTH_ELECTRON_V;
	} else if ( strcmp(sp+1, "kV") == 0 ) {
		*punit = WAVELENGTH_ELECTRON_KV;
	} else {
		free(valcpy);
		return 1;
	}

	sp[0] = '\0';
	*p_from = valcpy;
	return 0;
}


static int parse_wavelength(const char *val,
                            char **p_from,
                            enum wavelength_unit *punit)
{
	char *valcpy;
	char *sp;

	valcpy = strdup(val);
	if ( valcpy == NULL ) return 1;

	/* "wavelength" directive must have explicit units */
	sp = strchr(valcpy, ' ');
	if ( sp == NULL ) {
		free(valcpy);
		return 1;
	}

	if ( strcmp(sp+1, "m") == 0 ) {
		*punit = WAVELENGTH_M;
	} else if ( strcmp(sp+1, "A") == 0 ) {
		*punit = WAVELENGTH_A;
	} else {
		free(valcpy);
		return 1;
	}

	sp[0] = '\0';
	*p_from = valcpy;
	return 0;
}


static int parse_photon_energy(const char *val,
                               char **p_from,
                               enum wavelength_unit *punit)
{
	char *valcpy;
	char *sp;

	valcpy = strdup(val);
	if ( valcpy == NULL ) return 1;

	/* "photon_energy" is the only one of the wavelength
	 * directives which is allowed to not have units */
	sp = strchr(valcpy, ' ');
	if ( sp == NULL ) {
		*punit = WAVELENGTH_PHOTON_EV;
	} else if ( strcmp(sp+1, "eV") == 0 ) {
		*punit = WAVELENGTH_PHOTON_EV;
		sp[0] = '\0';
	} else if ( strcmp(sp+1, "keV") == 0 ) {
		*punit = WAVELENGTH_PHOTON_KEV;
		sp[0] = '\0';
	} else {
		/* Unit specified, but unrecognised */
		free(valcpy);
		return 1;
	}

	*p_from = valcpy;
	return 0;
}


static int parse_peak_layout(const char *val,
                             enum peak_layout *layout)
{
	if ( strcmp(val, "auto") == 0 ) {
		*layout = PEAK_LIST_AUTO;
		return 0;
	}

	if ( strcmp(val, "cxi") == 0 ) {
		*layout = PEAK_LIST_CXI;
		return 0;
	}

	if ( (strcmp(val, "list3") == 0) ) {
		*layout = PEAK_LIST_LIST3;
		return 0;
	}

	return 1;
}


static int parse_toplevel(DataTemplate *dt,
                          const char *key,
                          const char *val,
                          struct rg_definition ***rg_defl,
                          struct rgc_definition ***rgc_defl,
                          int *n_rg_defs,
                          int *n_rgc_defs,
                          struct panel_template *defaults)
{

	if ( strcmp(key, "mask_bad") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			dt->mask_bad = v;
		} else {
			return 1;
		}

	} else if ( strcmp(key, "mask_good") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			dt->mask_good = v;
		} else {
			return 1;
		}

	} else if ( strcmp(key, "detector_shift_x") == 0 ) {
		dt->shift_x_from = strdup(val);

	} else if ( strcmp(key, "detector_shift_y") == 0 ) {
		dt->shift_y_from = strdup(val);

	} else if ( strcmp(key, "photon_energy") == 0 ) {
		return parse_photon_energy(val,
		                           &dt->wavelength_from,
		                           &dt->wavelength_unit);

	} else if ( strcmp(key, "electron_voltage") == 0 ) {
		return parse_electron_voltage(val,
		                              &dt->wavelength_from,
		                              &dt->wavelength_unit);

	} else if ( strcmp(key, "wavelength") == 0 ) {
		return parse_wavelength(val,
		                        &dt->wavelength_from,
		                        &dt->wavelength_unit);

	} else if ( strcmp(key, "peak_list") == 0 ) {
		dt->peak_list = strdup(val);

	} else if ( strcmp(key, "peak_list_type") == 0 ) {
		return parse_peak_layout(val, &dt->peak_list_type);

	} else if ( strcmp(key, "bandwidth") == 0 ) {
		double v;
		char *end;
		v = strtod(val, &end);
		if ( (val[0] != '\0') && (end[0] == '\0') ) {
			dt->bandwidth = v;
		} else {
			ERROR("Invalid value for bandwidth\n");
		}

	} else if (strncmp(key, "rigid_group", 11) == 0
	        && strncmp(key, "rigid_group_collection", 22) != 0 ) {

		struct rg_definition **new;

		new = realloc(*rg_defl,
		             ((*n_rg_defs)+1)*sizeof(struct rg_definition*));
		*rg_defl = new;

		(*rg_defl)[*n_rg_defs] = malloc(sizeof(struct rg_definition));
		(*rg_defl)[*n_rg_defs]->name = strdup(key+12);
		(*rg_defl)[*n_rg_defs]->pns = strdup(val);
		*n_rg_defs = *n_rg_defs+1;

	} else if ( strncmp(key, "rigid_group_collection", 22) == 0 ) {

		struct rgc_definition **new;

		new = realloc(*rgc_defl, ((*n_rgc_defs)+1)*
		      sizeof(struct rgc_definition*));
		*rgc_defl = new;

		(*rgc_defl)[*n_rgc_defs] =
		                   malloc(sizeof(struct rgc_definition));
		(*rgc_defl)[*n_rgc_defs]->name = strdup(key+23);
		(*rgc_defl)[*n_rgc_defs]->rgs = strdup(val);
		*n_rgc_defs = *n_rgc_defs+1;

	} else if ( parse_field_for_panel(defaults, key, val, dt) ) {
		return 1;
	}

	return 0;
}


static int dt_num_path_placeholders(const char *str)
{
	size_t i, len;
	int n_pl = 0;

	if ( str == NULL ) return 0;

	len = strlen(str);
	for ( i=0; i<len; i++ ) {
		if ( str[i] == '%' ) n_pl++;
	}

	return n_pl;
}


signed int find_dim(signed int *dims, int which)
{
	int i;

	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( dims[i] == DIM_UNDEFINED ) break;
		if ( dims[i] == which ) return i;
	}

	return -1;
}


static int lookup_panel(const char *panel_name,
                        const DataTemplate *dt,
                        int *res)
{
	int i;

	/* If there is exactly one panel, you can get away without
	 * specifying the panel name */
	if ( (panel_name == NULL) && (dt->n_panels == 1) ) {
		*res = 0;
		return 0;
	}

	for ( i=0; i<dt->n_panels; i++ ) {
		if ( strcmp(dt->panels[i].name, panel_name) == 0 ) {
			*res = i;
			return 0;
		}
	}

	return 1;
}


DataTemplate *data_template_new_from_string(const char *string_in)
{
	DataTemplate *dt;
	char **bits;
	int done = 0;
	int i;
	int rgi, rgci;
	int reject = 0;
	struct rg_definition **rg_defl = NULL;
	struct rgc_definition **rgc_defl = NULL;
	int n_rg_definitions = 0;
	int n_rgc_definitions = 0;
	char *string;
	char *string_orig;
	size_t len;
	int num_data_pl;
	int num_mask_pl;
	int num_satmap_pl;
	struct panel_template defaults;

	dt = calloc(1, sizeof(DataTemplate));
	if ( dt == NULL ) return NULL;

	dt->n_panels = 0;
	dt->panels = NULL;
	dt->n_bad = 0;
	dt->bad = NULL;
	dt->mask_good = 0;
	dt->mask_bad = 0;
	dt->n_rigid_groups = 0;
	dt->rigid_groups = NULL;
	dt->n_rg_collections = 0;
	dt->rigid_group_collections = NULL;
	dt->bandwidth = 0.00000001;
	dt->peak_list = NULL;
	dt->shift_x_from = NULL;
	dt->shift_y_from = NULL;

	/* The default defaults... */
	defaults.orig_min_fs = -1;
	defaults.orig_min_ss = -1;
	defaults.orig_max_fs = -1;
	defaults.orig_max_ss = -1;
	defaults.cnx = NAN;
	defaults.cny = NAN;
	defaults.cnz_from = NULL;
	defaults.cnz_offset = 0.0;
	defaults.pixel_pitch = -1.0;
	defaults.bad = 0;
	defaults.fsx = 1.0;
	defaults.fsy = 0.0;
	defaults.fsz = 0.0;
	defaults.ssx = 0.0;
	defaults.ssy = 1.0;
	defaults.ssz = 0.0;
	defaults.rail_x = NAN;  /* The actual default rail direction */
	defaults.rail_y = NAN;  /*  is below */
	defaults.rail_z = NAN;
        defaults.clen_for_centering = NAN;
	defaults.adu_scale = NAN;
	defaults.adu_scale_unit = ADU_PER_PHOTON;
	for ( i=0; i<MAX_FLAG_VALUES; i++ ) defaults.flag_values[i] = 0;
	for ( i=0; i<MAX_FLAG_VALUES; i++ ) defaults.flag_types[i] = FLAG_NOTHING;
	defaults.max_adu = +INFINITY;
	defaults.mask = NULL;
	defaults.mask_file = NULL;
	defaults.satmap = NULL;
	defaults.satmap_file = NULL;
	defaults.data = strdup("/data/data");
	defaults.name = NULL;
	defaults.dims[0] = DIM_SS;
	defaults.dims[1] = DIM_FS;
	for ( i=2; i<MAX_DIMS; i++ ) defaults.dims[i] = DIM_UNDEFINED;

	string = strdup(string_in);
	if ( string == NULL ) return NULL;
	len = strlen(string);
	for ( i=0; i<len; i++ ) {
		if ( string_in[i] == '\r' ) string[i] = '\n';
	}

	/* Becaue 'string' will get modified */
	string_orig = string;

	do {

		char *line;
		struct dt_badregion *badregion = NULL;
		struct panel_template *panel = NULL;

		/* Copy the next line from the big string */
		const char *nl = strchr(string, '\n');
		if ( nl != NULL ) {
			size_t len = nl - string;
			line = strndup(string, nl-string);
			line[len] = '\0';
			string += len+1;
		} else {
			line = strdup(string);
			done = 1;
		}

		/* Trim leading spaces */
		i = 0;
		char *line_orig = line;
		while ( (line_orig[i] == ' ')
		     || (line_orig[i] == '\t') ) i++;
		line = strdup(line+i);
		free(line_orig);

		/* Stop at comment symbol */
		char *comm = strchr(line, ';');
		if ( comm != NULL ) comm[0] = '\0';

		/* Nothing left? Entire line was commented out,
		 * and can be silently ignored */
		if ( line[0] == '\0' ) {
			free(line);
			continue;
		}

		/* Find the equals sign */
		char *eq = strchr(line, '=');
		if ( eq == NULL ) {
			ERROR("Bad line in geometry file: '%s'\n", line);
			free(line);
			reject = 1;
			continue;
		}

		/* Split into two strings */
		eq[0] = '\0';
		char *val = eq+1;

		/* Trim leading and trailing spaces in value */
		while ( (val[0] == ' ') || (val[0] == '\t') ) val++;
	        notrail(val);

	        /* Trim trailing spaces in key
	         * (leading spaces already done above) */
	        notrail(line);

	        /* Find slash after panel name */
	        char *slash = strchr(line, '/');
	        if ( slash == NULL ) {

			/* Top-level option */
		        if ( parse_toplevel(dt, line, val,
		                            &rg_defl,
		                            &rgc_defl,
		                            &n_rg_definitions,
		                            &n_rgc_definitions,
		                            &defaults) )
			{
				ERROR("Invalid top-level line '%s'\n",
				      line);
				reject = 1;
			}
			free(line);
			continue;
		}

	        slash[0] = '\0';
	        char *key = slash+1;
	        /* No space trimming this time - must be "panel/key" */

	        /* Find either panel or bad region */
		if ( strncmp(line, "bad", 3) == 0 ) {
			badregion = find_bad_region_by_name(dt, line);
			if ( badregion == NULL ) {
				badregion = new_bad_region(dt, line);
			}
		} else {
			panel = find_panel_by_name(dt, line);
			if ( panel == NULL ) {
				panel = new_panel(dt, line, &defaults);
			}
		}

		if ( panel != NULL ) {
			if ( parse_field_for_panel(panel, key, val,
			                           dt) ) reject = 1;
		} else {
			if ( parse_field_bad(badregion, key,
			                     val) ) reject = 1;
		}

		free(line);

	} while ( !done );

	if ( dt->n_panels == 0 ) {
		ERROR("No panel descriptions in geometry file.\n");
		free(dt);
		return NULL;
	}

	num_data_pl = dt_num_path_placeholders(dt->panels[0].data);
	num_mask_pl = dt_num_path_placeholders(dt->panels[0].mask);
	num_satmap_pl = dt_num_path_placeholders(dt->panels[0].satmap);

	/* This is because the "data" path will be used to expand
	 * the path to generate the event list */
	if ( (num_mask_pl > num_data_pl)
	  || (num_satmap_pl > num_data_pl) )
	{
		ERROR("Mask and saturation map paths must have fewer "
		      "placeholders than image data path.\n");
		reject = 1;
	}

	for ( i=0; i<dt->n_panels; i++ ) {

		struct panel_template *p = &dt->panels[i];
		signed int dim_fs = find_dim(p->dims, DIM_FS);
		signed int dim_ss = find_dim(p->dims, DIM_SS);

		if ( (dim_fs<0) || (dim_ss<0) ) {
			ERROR("Panel %s does not have dimensions "
			      "assigned to both fs and ss.\n",
			      p->name);
			reject = 1;
		}

		if ( dim_ss > dim_fs ) {
			ERROR("Fast scan dimension must be lower than "
			      "slow scan (panel %s)\n", p->name);
			reject = 1;
		}

		if ( p->orig_min_fs < 0 ) {
			ERROR("Please specify the minimum FS coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( p->orig_max_fs < 0 ) {
			ERROR("Please specify the maximum FS coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( p->orig_min_ss < 0 ) {
			ERROR("Please specify the minimum SS coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( p->orig_max_ss < 0 ) {
			ERROR("Please specify the maximum SS coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( isnan(p->cnx) ) {
			ERROR("Please specify the corner X coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( isnan(p->cny) ) {
			ERROR("Please specify the corner Y coordinate for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( p->cnz_from == NULL ) {
			ERROR("Please specify the camera length for panel %s\n",
			      dt->panels[i].name);
			reject = 1;
		}
		if ( p->pixel_pitch < 0 ) {
			ERROR("Please specify the pixel size for"
			      " panel %s\n", dt->panels[i].name);
			reject = 1;
		}
		if ( p->data == NULL ) {
			ERROR("Please specify the data location for panel %s\n",
			      p->name);
			reject = 1;
		}
		if ( isnan(p->adu_scale) ) {
			ERROR("Please specify either adu_per_eV or "
			      "adu_per_photon for panel %s\n",
			      dt->panels[i].name);
			reject = 1;
		}

		if ( isnan(p->clen_for_centering) && !isnan(p->rail_x) )
		{
			ERROR("You must specify clen_for_centering if you "
			      "specify the rail direction (panel %s)\n",
			      p->name);
			reject = 1;
		}

		if ( (p->mask_file != NULL) && (p->mask == NULL) ) {
			ERROR("You have specified 'mask_file' but not 'mask'.  "
			      "'mask_file' will therefore have no effect.  "
			      "(panel %s)\n", p->name);
			reject = 1;
		}

		if ( dt_num_path_placeholders(p->data) != num_data_pl ) {
			ERROR("Data locations for all panels must "
			      "have the same number of placeholders\n");
			reject = 1;
		}

		if ( dt_num_path_placeholders(p->mask) != num_mask_pl ) {
			ERROR("Mask locations for all panels must "
			      "have the same number of placeholders\n");
			reject = 1;
		}

		if ( dt_num_path_placeholders(p->satmap) != num_satmap_pl ) {
			ERROR("Satmap locations for all panels must "
			      "have the same number of placeholders\n");
			reject = 1;
		}

		/* The default rail direction */
		if ( isnan(p->rail_x) ) {
			p->rail_x = 0.0;
			p->rail_y = 0.0;
			p->rail_z = 1.0;
		}
		if ( isnan(p->clen_for_centering) ) p->clen_for_centering = 0.0;

	}

	for ( i=0; i<dt->n_bad; i++ ) {

		if ( dt->bad[i].is_fsss == 99 ) {
			ERROR("Please specify the coordinate ranges for"
			      " bad region %s\n", dt->bad[i].name);
			reject = 1;
		}

		if ( dt->bad[i].is_fsss ) {
			if ( dt->bad[i].panel_name == NULL ) {

				ERROR("Panel not specified for bad region '%s'\n",
				      dt->bad[i].name);
				reject = 1;

			} else if ( lookup_panel(dt->bad[i].panel_name, dt,
			                         &dt->bad[i].panel_number) )
			{
				ERROR("No such panel '%s' for bad region %s\n",
				      dt->bad[i].panel_name, dt->bad[i].name);
				reject = 1;

			} else {
				struct panel_template *p;
				struct dt_badregion *bad;
				int r = 0;
				p = &dt->panels[dt->bad[i].panel_number];
				bad = &dt->bad[i];
				if ( bad->min_fs < p->orig_min_fs ) r = 1;
				if ( bad->min_ss < p->orig_min_ss ) r = 1;
				if ( bad->max_fs > p->orig_max_fs ) r = 1;
				if ( bad->max_ss > p->orig_max_ss ) r = 1;
				if ( r ) {
					ERROR("Bad region '%s' is outside the "
					      "panel bounds (%s) as presented "
					      "in data (%i %i, %i %i inclusive): "
					      "Bad region %i,%i to %i, %i "
					      "inclusive\n",
					      bad->name, p->name,
					      p->orig_min_fs, p->orig_min_ss,
					      p->orig_max_fs, p->orig_max_ss,
					      bad->min_fs, bad->min_ss,
					      bad->max_fs, bad->max_ss);
					reject = 1;
				}
				bad->min_fs -= p->orig_min_fs;
				bad->max_fs -= p->orig_min_fs;
				bad->min_ss -= p->orig_min_ss;
				bad->max_ss -= p->orig_min_ss;
			}
		}
	}

	free(defaults.cnz_from);
	free(defaults.data);
	free(defaults.mask);
	free(defaults.mask_file);

	for ( rgi=0; rgi<n_rg_definitions; rgi++) {

		int pi, n1;
		struct rigid_group *rigidgroup = NULL;

		rigidgroup = find_or_add_rg(dt, rg_defl[rgi]->name);

		n1 = assplode(rg_defl[rgi]->pns, ",", &bits, ASSPLODE_NONE);

		for ( pi=0; pi<n1; pi++ ) {

			int panel_number;
			if ( data_template_panel_name_to_number(dt,
			                                        bits[pi],
			                                        &panel_number) )
			{
				ERROR("Cannot add panel to rigid group\n");
				ERROR("Panel not found: %s\n", bits[pi]);
				return NULL;
			}
			add_to_rigid_group(rigidgroup, panel_number);
			free(bits[pi]);

		}
		free(bits);
		free(rg_defl[rgi]->name);
		free(rg_defl[rgi]->pns);
		free(rg_defl[rgi]);
	}
	free(rg_defl);

	for ( rgci=0; rgci<n_rgc_definitions; rgci++ ) {

		int rgi, n2;
		struct rg_collection *rgcollection = NULL;

		rgcollection = find_or_add_rg_coll(dt, rgc_defl[rgci]->name);

		n2 = assplode(rgc_defl[rgci]->rgs, ",", &bits, ASSPLODE_NONE);

		for ( rgi=0; rgi<n2; rgi++ ) {

			struct rigid_group *r;

			r = find_rigid_group_by_name(dt, bits[rgi]);
			if ( r == NULL ) {
				ERROR("Cannot add rigid group to collection\n");
				ERROR("Rigid group not found: %s\n", bits[rgi]);
				return NULL;
			}
			add_to_rigid_group_coll(rgcollection, r);
			free(bits[rgi]);
		}
		free(bits);
		free(rgc_defl[rgci]->name);
		free(rgc_defl[rgci]->rgs);
		free(rgc_defl[rgci]);

	}
	free(rgc_defl);

	free(string_orig);

	if ( reject ) return NULL;

	return dt;
}


DataTemplate *data_template_new_from_file(const char *filename)
{
	char *contents;
	DataTemplate *dt;

	contents = load_entire_file(filename);
	if ( contents == NULL ) {
		ERROR("Failed to load geometry file '%s'\n", filename);
		return NULL;
	}

	dt = data_template_new_from_string(contents);
	free(contents);
	return dt;
}


void data_template_free(DataTemplate *dt)
{
	int i;

	if ( dt == NULL ) return;

	free_all_rigid_groups(dt);
	free_all_rigid_group_collections(dt);

	for ( i=0; i<dt->n_panels; i++ ) {
		free(dt->panels[i].name);
		free(dt->panels[i].data);
		free(dt->panels[i].mask);
		free(dt->panels[i].mask_file);
		free(dt->panels[i].satmap);
		free(dt->panels[i].satmap_file);
		free(dt->panels[i].cnz_from);
	}

	free(dt->wavelength_from);
	free(dt->peak_list);

	free(dt->panels);
	free(dt->bad);
	free(dt);
}


static int data_template_find_panel(const DataTemplate *dt,
                                    int fs, int ss, int *ppn)
{
	int p;

	for ( p=0; p<dt->n_panels; p++ ) {
		if ( (fs >= dt->panels[p].orig_min_fs)
		  && (fs < dt->panels[p].orig_max_fs+1)
		  && (ss >= dt->panels[p].orig_min_ss)
		  && (ss < dt->panels[p].orig_max_ss+1) ) {
			*ppn = p;
			return 0;
		}
	}

	return 1;
}


int data_template_file_to_panel_coords(const DataTemplate *dt,
                                       float *pfs, float *pss,
                                       int *ppn)
{
	int pn;

	if ( data_template_find_panel(dt, *pfs, *pss, &pn) ) {
		return 1;
	}

	*ppn = pn;
	*pfs = *pfs - dt->panels[pn].orig_min_fs;
	*pss = *pss - dt->panels[pn].orig_min_ss;
	return 0;
}


int data_template_panel_to_file_coords(const DataTemplate *dt,
                                       int pn, float *pfs, float *pss)
{
	if ( pn >= dt->n_panels ) return 1;
	*pfs = *pfs + dt->panels[pn].orig_min_fs;
	*pss = *pss + dt->panels[pn].orig_min_ss;
	return 0;
}


const char *data_template_panel_number_to_name(const DataTemplate *dt,
                                               int pn)
{
	if ( pn >= dt->n_panels ) return NULL;
	return dt->panels[pn].name;
}


int data_template_panel_name_to_number(const DataTemplate *dt,
                                       const char *panel_name,
                                       int *pn)
{
	int i;

	if ( panel_name == NULL ) return 1;

	for ( i=0; i<dt->n_panels; i++ ) {
		if ( strcmp(panel_name, dt->panels[i].name) == 0 ) {
			*pn = i;
			return 0;
		}
	}

	return 1;
}


void data_template_add_copy_header(DataTemplate *dt,
                                   const char *header)
{
	/* FIXME: Add "header" to list of things to copy */
	STATUS("Adding %s\n", header);
}


static int dt_num_placeholders(const struct panel_template *p)
{
	int i;
	int n_pl = 0;
	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( p->dims[i] == DIM_PLACEHOLDER ) n_pl++;
	}
	return n_pl;
}


int data_template_get_slab_extents(const DataTemplate *dt,
                                   int *pw, int *ph)
{
	int w, h;
	char *data_from;
	int i;

	data_from = dt->panels[0].data;

	w = 0;  h = 0;
	for ( i=0; i<dt->n_panels; i++ ) {

		struct panel_template *p = &dt->panels[i];

		if ( strcmp(data_from, p->data) != 0 ) {
			/* Not slabby */
			return 1;
		}

		if ( dt_num_placeholders(p) > 0 ) {
			/* Not slabby */
			return 1;
		}

		if ( p->orig_max_fs > w ) {
			w = p->orig_max_fs;
		}
		if ( p->orig_max_ss > h ) {
			h = p->orig_max_ss;
		}

	}

	/* Inclusive -> exclusive */
	*pw = w + 1;
	*ph = h + 1;
	return 0;
}
