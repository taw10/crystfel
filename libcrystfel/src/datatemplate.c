/*
 * datatemplate.c
 *
 * Data template structure
 *
 * Copyright © 2019-2024 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2024 Thomas White <taw@physics.org>
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
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "utils.h"
#include "datatemplate.h"
#include "image.h"

#include "datatemplate_priv.h"


/**
 * \file datatemplate.h
 */

static struct panel_group_template *find_group(const DataTemplate *dt, const char *name)
{
	int i;

	for ( i=0; i<dt->n_groups; i++ ) {
		if ( strcmp(dt->groups[i]->name, name) == 0 ) {
			return dt->groups[i];
		}
	}

	return NULL;
}


static struct panel_group_template *add_group(const char *name, DataTemplate *dt)
{
	struct panel_group_template *gt;

	if ( find_group(dt, name) != NULL ) {
		ERROR("Duplicate panel group '%s'\n", name);
		return NULL;
	}

	if ( dt->n_groups >= MAX_PANEL_GROUPS ) {
		ERROR("Too many panel groups\n");
		return NULL;
	}

	gt = cfmalloc(sizeof(struct panel_group_template));
	if ( gt == NULL ) return NULL;

	gt->name = cfstrdup(name);
	gt->n_children = 0;

	if ( gt->name == NULL ) {
		cffree(gt);
		return NULL;
	}

	dt->groups[dt->n_groups++] = gt;

	return gt;
}


static int add_group_members(const char *name, DataTemplate *dt,
                             char **members, int n_members)
{
	struct panel_group_template *gt;
	int i;
	int fail = 0;

	if ( n_members == 0 ) {
		ERROR("Panel group '%s' has no members\n", name);
		fail = 1;
	}

	if ( n_members > MAX_PANEL_GROUP_CHILDREN ) {
		ERROR("Panel group '%s' has too many members\n", name);
		fail = 1;
	}


	/* A simple typo in the geometry file can segfault other
	* stuff, so check */
	for ( i=0; i<n_members; i++ ) {
		int j;
		for ( j=0; j<i; j++ ) {
			if ( strcmp(members[i], members[j]) == 0 ) {
				ERROR("Duplicate member '%s' in group '%s'\n",
				      members[i], name);
				fail = 1;
			}
		}
	}

	if ( fail ) return fail;

	gt = add_group(name, dt);
	if ( gt == NULL ) {
		ERROR("Failed to add group\n");
		return 1;
	}

	for ( i=0; i<n_members; i++ ) {
		gt->children[i] = find_group(dt, members[i]);
		if ( gt->children[i] == NULL ) {
			ERROR("Unknown panel group '%s'\n", members[i]);
			ERROR("Make sure the hierarchy groups definitions are AFTER the "
			      "panel definitions in the geometry file, and start from "
			      "the lowest hierachy level.\n");
			fail = 1;
		}
	}

	gt->n_children = n_members;

	return fail;
}


static int parse_group(const char *name, DataTemplate *dt, const char *val)
{
	int n_members;
	char **members;
	int i;
	int fail = 0;

	n_members = assplode(val, ",",  &members, ASSPLODE_NONE);

	fail = add_group_members(name, dt, members, n_members);

	for ( i=0; i<n_members; i++ ) cffree(members[i]);
	cffree(members);

	return fail;
}


static struct panel_template *new_panel(DataTemplate *det,
                                        const char *name,
                                        struct panel_template *defaults)
{
	struct panel_template *new;
	int i;

	det->n_panels++;
	det->panels = cfrealloc(det->panels,
	                        det->n_panels*sizeof(struct panel_template));

	new = &det->panels[det->n_panels-1];
	memcpy(new, defaults, sizeof(struct panel_template));

	/* Set name */
	new->name = cfstrdup(name);

	/* Copy strings */
	new->data = safe_strdup(defaults->data);
	new->satmap = safe_strdup(defaults->satmap);
	new->satmap_file = safe_strdup(defaults->satmap_file);
	for ( i=0; i<MAX_MASKS; i++ ) {
		new->masks[i].data_location = safe_strdup(defaults->masks[i].data_location);
		new->masks[i].filename = safe_strdup(defaults->masks[i].filename);
	}

	/* Create a new group just for this panel */
	add_group(name, det);

	return new;
}


static struct dt_badregion *new_bad_region(DataTemplate *det, const char *name)
{
	struct dt_badregion *new;

	det->n_bad++;
	det->bad = cfrealloc(det->bad, det->n_bad*sizeof(struct dt_badregion));

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
		a = cfmalloc(len+1);
		snprintf(a, len+1, "+%s", a_orig);
		a[len] = '\0';

	} else {
		a = cfstrdup(a_orig);
	}

	/* Count the expressions */
	nexp = 0;
	for ( i=0; i<len; i++ ) {
		if ( (a[i] == '+') || (a[i] == '-') ) nexp++;
	}

	bits = cfcalloc(nexp, sizeof(char *));

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
			bits[idx] = cfmalloc(len+1);
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
	cffree(a);

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

		cffree(bits[i]);

	}
	cffree(bits);

	return 0;
}


static int set_dim(struct panel_template *panel, int dimension,
                   const char *val, int def)
{
	if ( dimension >= MAX_DIMS ) {
		ERROR("Too many dimensions!\n");
		return 1;
	}

	if ( strcmp(val, "fs") == 0 ) {
		panel->dims[dimension] = DIM_FS;
		panel->dims_default[dimension] = def;
	} else if ( strcmp(val, "ss") == 0 ) {
		panel->dims[dimension] = DIM_SS;
		panel->dims_default[dimension] = def;
	} else if ( strcmp(val, "%") == 0 ) {
		panel->dims[dimension] = DIM_PLACEHOLDER;
		panel->dims_default[dimension] = def;
	} else {
		char *endptr;
		unsigned long int fix_val = strtoul(val, &endptr, 10);
		if ( endptr[0] != '\0' ) {
			ERROR("Invalid dimension value '%s'\n", val);
			return 1;
		} else {
			panel->dims[dimension] = fix_val;
			panel->dims_default[dimension] = def;
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


static int parse_mask(struct panel_template *panel,
                      const char *key_orig,
                      const char *val,
                      int def)
{
	int n;
	char *key;

	if ( sscanf(key_orig, "mask%d_", &n) != 1 ) {
		ERROR("Invalid mask directive '%s'\n", key_orig);
		return 1;
	}

	key = cfstrdup(key_orig);
	if ( key == NULL ) return 1;

	key[4] = '_';

	/* The mask number has been replaced with '_'.
	 * Double underscore is deliberate! */
	if ( strcmp(key, "mask__file") == 0 ) {

		panel->masks[n].filename = cfstrdup(val);

	} else if ( strcmp(key, "mask__data") == 0 ) {

		if ( strncmp(val, "/", 1) != 0 ) {
			ERROR("Invalid mask location '%s'\n", val);
			cffree(key);
			return 1;
		}
		panel->masks[n].data_location = cfstrdup(val);

	} else if ( strcmp(key, "mask__goodbits") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			panel->masks[n].good_bits = v;
		} else {
			cffree(key);
			return 1;
		}

	} else if ( strcmp(key, "mask__badbits") == 0 ) {

		char *end;
		double v = strtod(val, &end);

		if ( end != val ) {
			panel->masks[n].bad_bits = v;
		} else {
			cffree(key);
			return 1;
		}

	} else {

		ERROR("Invalid mask directive '%s'\n", key_orig);
		cffree(key);
		return 1;
	}

	panel->masks[n].mask_default = def;
	cffree(key);
	return 0;
}


static int parse_field_for_panel(struct panel_template *panel, const char *key,
                                 const char *val, DataTemplate *det,
                                 int def)
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
	} else if ( strcmp(key, "adu_per_eV") == 0 ) {
		panel->adu_scale = atof(val);
		panel->adu_scale_unit = ADU_PER_EV;
		panel->adu_scale_default = def;
	} else if ( strcmp(key, "adu_per_photon") == 0 ) {
		panel->adu_scale = atof(val);
		panel->adu_scale_unit = ADU_PER_PHOTON;
		panel->adu_scale_default = def;
	} else if ( strcmp(key, "clen") == 0 ) {
		ERROR("'clen' is a top-level property in this version of CrystFEL.\n");
		reject = 1;

	} else if ( strcmp(key, "data") == 0 ) {
		cffree(panel->data);
		panel->data = cfstrdup(val);
		panel->data_default = def;

	} else if ( strcmp(key, "mask_edge_pixels") == 0 ) {
		if ( convert_int(val, &panel->mask_edge_pixels) ) {
			ERROR("Invalid value for %s/mask_edge_pixels (%s)\n",
			      panel->name, val);
			reject = 1;
		}
		panel->mask_edge_pixels_default = def;

	} else if ( strcmp(key, "mask_bad") == 0 ) {
		parse_field_for_panel(panel, "mask0_badbits", val, det, def);
	} else if ( strcmp(key, "mask_good") == 0 ) {
		parse_field_for_panel(panel, "mask0_goodbits", val, det, def);
	} else if ( strcmp(key, "mask") == 0 ) {
		parse_field_for_panel(panel, "mask0_data", val, det, def);
	} else if ( strcmp(key, "mask_file") == 0 ) {
		parse_field_for_panel(panel, "mask0_file", val, det, def);

	} else if ( strncmp(key, "mask", 4) == 0 ) {
		reject = parse_mask(panel, key, val, def);

	} else if ( strcmp(key, "saturation_map") == 0 ) {
		panel->satmap = cfstrdup(val);
		panel->satmap_default = def;
	} else if ( strcmp(key, "saturation_map_file") == 0 ) {
		panel->satmap_file = cfstrdup(val);
		panel->satmap_file_default = def;

	} else if ( strcmp(key, "coffset") == 0) {
		panel->cnz_offset = atof(val);
	} else if ( strcmp(key, "res") == 0 ) {
		panel->pixel_pitch = 1.0/atof(val);
		panel->pixel_pitch_default = def;
	} else if ( strcmp(key, "max_adu") == 0 ) {
		panel->max_adu = atof(val);
		panel->max_adu_default = def;
		ERROR("WARNING: It's usually better not to set max_adu "
		      "in the geometry file.  Use --max-adu during "
		      "merging instead.\n");

	} else if ( strcmp(key, "flag_equal") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_EQUAL) ) {
			reject = -1;
		}
		panel->flag_values_default = def;
	} else if ( strcmp(key, "flag_lessthan") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_LESSTHAN) ) {
			reject = -1;
		}
		panel->flag_values_default = def;
	} else if ( strcmp(key, "flag_morethan") == 0 ) {
		if ( add_flag_value(panel, atof(val), FLAG_MORETHAN) ) {
			reject = -1;
		}
		panel->flag_values_default = def;

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
		char *endptr;
		if ( key[3] != '\0' ) {
			int dim_entry;
			dim_entry = strtoul(key+3, &endptr, 10);
			if ( endptr[0] != '\0' ) {
				ERROR("Invalid dimension number %s\n",
				      key+3);
			} else {
				if ( set_dim(panel, dim_entry, val, def) ) {
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
		badr->panel_name = cfstrdup(val);
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

	valcpy = cfstrdup(val);
	if ( valcpy == NULL ) return 1;

	/* "electron_voltage" directive must have explicit units */
	sp = strchr(valcpy, ' ');
	if ( sp == NULL ) {
		cffree(valcpy);
		return 1;
	}

	if ( strcmp(sp+1, "V") == 0 ) {
		*punit = WAVELENGTH_ELECTRON_V;
	} else if ( strcmp(sp+1, "kV") == 0 ) {
		*punit = WAVELENGTH_ELECTRON_KV;
	} else {
		cffree(valcpy);
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

	valcpy = cfstrdup(val);
	if ( valcpy == NULL ) return 1;

	/* "wavelength" directive must have explicit units */
	sp = strchr(valcpy, ' ');
	if ( sp == NULL ) {
		cffree(valcpy);
		return 1;
	}

	if ( strcmp(sp+1, "m") == 0 ) {
		*punit = WAVELENGTH_M;
	} else if ( strcmp(sp+1, "A") == 0 ) {
		*punit = WAVELENGTH_A;
	} else {
		cffree(valcpy);
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

	valcpy = cfstrdup(val);
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
		cffree(valcpy);
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

#define MAX_FOR_LATER (1024)

struct forlater
{
	char *keys[MAX_FOR_LATER];
	char *vals[MAX_FOR_LATER];
	int n_forlater;
};


static void store_for_later(struct forlater *fl, const char *key, const char *val)
{
	if ( fl->n_forlater >= MAX_FOR_LATER ) {
		ERROR("Too many lines stored.\n");
		return;
	}

	fl->keys[fl->n_forlater] = cfstrdup(key);
	fl->vals[fl->n_forlater] = cfstrdup(val);
	fl->n_forlater++;
}


static int parse_toplevel(DataTemplate *dt,
                          const char *key,
                          const char *val,
                          struct panel_template *defaults,
                          int *defaults_updated,
                          struct forlater *for_later)
{
	if ( strcmp(key, "detector_shift_x") == 0 ) {
		dt->shift_x_from = cfstrdup(val);

	} else if ( strcmp(key, "detector_shift_y") == 0 ) {
		dt->shift_y_from = cfstrdup(val);

	} else if ( strcmp(key, "clen") == 0 ) {
		dt->cnz_from = cfstrdup(val);

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
		dt->peak_list = cfstrdup(val);

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

	} else if ( strncmp(key, "rigid_group", 11) == 0 ) {

		/* Rigid group lines are ignored in this version */

	} else if ( strncmp(key, "group_", 6) == 0 ) {

		if ( for_later != NULL ) {
			store_for_later(for_later, key, val);
		} else {
			if ( parse_group(key+6, dt, val) ) {
				return 1;
			}
		}

	} else {

		/* If there are any panels, the value in 'defaults' gets marked
		 * as "not default".  This will cause it to be written out for
		 * each subsequent panel. */
		if ( parse_field_for_panel(defaults, key, val, dt, (dt->n_panels==0)) == 0 ) {
			*defaults_updated = 1;
		} else {
			return 1;
		}
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

	if ( panel_name == NULL ) {
		ERROR("Panel name must be specified.\n");
		return 1;
	}

	for ( i=0; i<dt->n_panels; i++ ) {
		if ( strcmp(dt->panels[i].name, panel_name) == 0 ) {
			*res = i;
			return 0;
		}
	}

	return 1;
}


static int check_mask_and_satmap_placeholders(const DataTemplate *dt)
{
	int i;

	for ( i=0; i<dt->n_panels; i++ ) {

		int num_data_pl;
		int num_satmap_pl;
		int j;

		num_data_pl = dt_num_path_placeholders(dt->panels[i].data);
		num_satmap_pl = dt_num_path_placeholders(dt->panels[i].satmap);

		if ( num_satmap_pl > num_data_pl ) return 1;

		for ( j=0; j<MAX_MASKS; j++ ) {

			int num_mask_pl;

			/* Unused slot? */
			if ( dt->panels[i].masks[j].data_location == NULL ) continue;

			num_mask_pl = dt_num_path_placeholders(dt->panels[i].masks[j].data_location);
			if ( num_mask_pl > num_data_pl ) return 1;
		}
	}

	return 0;
}


static int try_guess_panel(struct dt_badregion *bad, DataTemplate *dt)
{
	if ( dt->n_panels == 1 ) {
		bad->panel_name = dt->panels[0].name;
		ERROR("WARNING: Assuming bad_%s/panel = %s\n",
		      bad->name, dt->panels[0].name);
		return 1;
	}

	return 0;
}


static void show_group(const struct panel_group_template *gt, int level)
{
	int i;

	for ( i=0; i<level; i++ ) STATUS("  ");

	if ( gt == NULL ) {
		STATUS("!!!\n");
		return;
	}

	STATUS("%s\n", gt->name);

	for ( i=0; i<gt->n_children; i++ ) {
		show_group(gt->children[i], level+1);
	}
}


void data_template_show_hierarchy(const DataTemplate *dtempl)
{
	STATUS("Hierarchy:\n");
	show_group(find_group(dtempl, "all"), 0);
}


DataTemplate *data_template_new_from_string(const char *string_in)
{
	DataTemplate *dt;
	int done = 0;
	int i;
	int reject = 0;
	char *string;
	char *string_orig;
	size_t len;
	struct panel_template defaults;
	int have_unused_defaults = 0;
	struct forlater lines_for_later;

	lines_for_later.n_forlater = 0;

	dt = cfcalloc(1, sizeof(DataTemplate));
	if ( dt == NULL ) return NULL;

	dt->n_panels = 0;
	dt->panels = NULL;
	dt->n_bad = 0;
	dt->bad = NULL;
	dt->bandwidth = 0.00000001;
	dt->peak_list = NULL;
	dt->shift_x_from = NULL;
	dt->shift_y_from = NULL;
	dt->cnz_from = NULL;
	dt->n_headers_to_copy = 0;
	dt->n_groups = 0;

	/* The default defaults... */
	defaults.orig_min_fs = -1;
	defaults.orig_min_ss = -1;
	defaults.orig_max_fs = -1;
	defaults.orig_max_ss = -1;
	defaults.cnx = NAN;
	defaults.cny = NAN;
	defaults.cnz_offset = 0.0;
	defaults.pixel_pitch = -1.0;
	defaults.pixel_pitch_default = 1;
	defaults.bad = 0;
	defaults.mask_edge_pixels = 0;
	defaults.mask_edge_pixels_default = 1;
	defaults.fsx = NAN;
	defaults.fsy = NAN;
	defaults.fsz = NAN;
	defaults.ssx = NAN;
	defaults.ssy = NAN;
	defaults.ssz = NAN;
	defaults.adu_scale = NAN;
	defaults.adu_scale_unit = ADU_PER_PHOTON;
	defaults.adu_scale_default = 1;
	for ( i=0; i<MAX_FLAG_VALUES; i++ ) defaults.flag_values[i] = 0;
	for ( i=0; i<MAX_FLAG_VALUES; i++ ) defaults.flag_types[i] = FLAG_NOTHING;
	defaults.flag_values_default = 1;
	for ( i=0; i<MAX_MASKS; i++ ) {
		defaults.masks[i].data_location = NULL;
		defaults.masks[i].filename = NULL;
		defaults.masks[i].good_bits = 0;
		defaults.masks[i].bad_bits = 0;
		defaults.masks[i].mask_default = 1;
	}
	defaults.max_adu = +INFINITY;
	defaults.max_adu_default = 1;
	defaults.satmap = NULL;
	defaults.satmap_default = 1;
	defaults.satmap_file = NULL;
	defaults.satmap_file_default = 1;
	defaults.data = cfstrdup("/data/data");
	defaults.data_default = 1;
	defaults.name = NULL;
	defaults.dims[0] = DIM_SS;
	defaults.dims[1] = DIM_FS;
	for ( i=2; i<MAX_DIMS; i++ ) defaults.dims[i] = DIM_UNDEFINED;
	for ( i=0; i<MAX_DIMS; i++ ) defaults.dims_default[i] = 1;

	string = cfstrdup(string_in);
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
			size_t nlen = nl - string;
			line = cfstrndup(string, nlen);
			line[nlen] = '\0';
			string += nlen+1;
		} else {
			line = cfstrdup(string);
			done = 1;
		}

		/* Trim leading spaces */
		i = 0;
		char *line_orig = line;
		while ( (line_orig[i] == ' ')
		     || (line_orig[i] == '\t') ) i++;
		line = cfstrdup(line+i);
		cffree(line_orig);

		/* Stop at comment symbol */
		char *comm = strchr(line, ';');
		if ( comm != NULL ) comm[0] = '\0';

		/* Nothing left? Entire line was commented out,
		 * and can be silently ignored */
		if ( line[0] == '\0' ) {
			cffree(line);
			continue;
		}

		/* Find the equals sign */
		char *eq = strchr(line, '=');
		if ( eq == NULL ) {
			ERROR("Bad line in geometry file: '%s'\n", line);
			cffree(line);
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
		                            &defaults,
		                            &have_unused_defaults,
			                    &lines_for_later) )
			{
				ERROR("Invalid top-level line '%s'\n", line);
				reject = 1;
			}
			cffree(line);
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
				have_unused_defaults = 0;
			}
		}

		if ( panel != NULL ) {
			if ( parse_field_for_panel(panel, key, val, dt, 0) ) reject = 1;
		} else {
			if ( parse_field_bad(badregion, key, val) ) reject = 1;
		}

		cffree(line);

	} while ( !done );

	for ( i=0; i<lines_for_later.n_forlater; i++ ) {
		char *key = lines_for_later.keys[i];
		char *val = lines_for_later.vals[i];
		if ( parse_toplevel(dt, key, val,
		                    &defaults,
	                            &have_unused_defaults,
		                    NULL) )
		{
			ERROR("Invalid top-level line '%s' = '%s'\n", key, val);
			reject = 1;
		}
		cffree(key);
		cffree(val);
	}

	if ( dt->n_panels == 0 ) {
		ERROR("No panel descriptions in geometry file.\n");
		cffree(dt);
		return NULL;
	}

	if ( check_mask_and_satmap_placeholders(dt) ) {
		ERROR("Mask and saturation map paths must have fewer "
		      "placeholders than image data path.\n");
		reject = 1;
	}

	if ( have_unused_defaults ) {
		ERROR("WARNING: There are statements at the end of the geometry "
		      "file which have no effect because they only apply to "
		      "subsequently-defined panels.\n");
		reject = 1;
	}

	if ( dt->wavelength_from == NULL ) {
		ERROR("Geometry file must specify the wavelength "
		      "(value or location)\n");
		reject = 1;
	}

	if ( dt->cnz_from == NULL ) {
		ERROR("Geometry file must specify the camera length\n");
		reject = 1;
	}

	for ( i=0; i<dt->n_panels; i++ ) {

		int j;
		struct panel_template *p = &dt->panels[i];
		signed int dim_fs = find_dim(p->dims, DIM_FS);
		signed int dim_ss = find_dim(p->dims, DIM_SS);

		if ( (dim_fs<0) || (dim_ss<0) ) {
			ERROR("Panel %s does not have dimensions "
			      "assigned to both fs and ss.\n",
			      p->name);
			reject = 1;
		}

		if ( dim_ss >= dim_fs ) {
			ERROR("Fast scan dimension must be higher than "
			      "slow scan (panel %s)\n", p->name);
			reject = 1;
		}

		if ( isnan(p->fsx) ) {
			ERROR("Please specify the FS direction for panel %s\n",
			      dt->panels[i].name);
			reject = 1;
		}

		if ( isnan(p->ssx) ) {
			ERROR("Please specify the SS direction for panel %s\n",
			      dt->panels[i].name);
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

		for ( j=0; j<MAX_MASKS; j++ ) {

			struct mask_template *mt;
			mt = &p->masks[j];

			if ( (mt->filename != NULL)
			  && (mt->data_location == NULL) )
			{
				ERROR("You have specified filename but not data"
				      " location for mask %i of panel %s\n",
				      j, p->name);
				reject = 1;
			}

			if ( (mt->good_bits || mt->bad_bits)
			  && (mt->filename == NULL)
			  && (mt->data_location == NULL) )
			{
				ERROR("You have specified good/bad bits for "
				      "mask %i of panel %s, but not the mask "
				      "location.\n", j, p->name);
				reject = 1;
			}
		}

	}

	for ( i=0; i<dt->n_bad; i++ ) {

		if ( dt->bad[i].is_fsss == 99 ) {
			ERROR("Please specify the coordinate ranges for"
			      " bad region %s\n", dt->bad[i].name);
			reject = 1;
		}

		if ( dt->bad[i].is_fsss ) {
			if ( dt->bad[i].panel_name == NULL ) {

				if ( !try_guess_panel(&dt->bad[i], dt) ) {
					ERROR("Panel not specified for bad "
					      "region '%s'\n", dt->bad[i].name);
					reject = 1;
				}

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

	cffree(defaults.data);
	for ( i=0; i<MAX_MASKS; i++ ) {
		cffree(defaults.masks[i].data_location);
		cffree(defaults.masks[i].filename);
	}

	/* If no groups are defined, put everything in one group.
	 * This allows at least basic geometry refinement to work. */
	if ( dt->n_groups == dt->n_panels ) {
		char **allg = cfmalloc(dt->n_groups*sizeof(char *));
		if ( allg == NULL ) {
			ERROR("Failed to create top group\n");
		} else {
			int i;
			for ( i=0; i<dt->n_groups; i++ ) {
				allg[i] = dt->groups[i]->name;
			}
			add_group_members("all", dt, allg, dt->n_groups);
			cffree(allg);
		}
	}

	cffree(string_orig);

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
	cffree(contents);
	return dt;
}


void data_template_free(DataTemplate *dt)
{
	int i;

	if ( dt == NULL ) return;

	for ( i=0; i<dt->n_panels; i++ ) {

		int j;

		cffree(dt->panels[i].name);
		cffree(dt->panels[i].data);
		cffree(dt->panels[i].satmap);
		cffree(dt->panels[i].satmap_file);

		for ( j=0; j<MAX_MASKS; j++ ) {
			cffree(dt->panels[i].masks[j].filename);
			cffree(dt->panels[i].masks[j].data_location);
		}
	}

	for ( i=0; i<dt->n_headers_to_copy; i++ ) {
		cffree(dt->headers_to_copy[i]);
	}

	for ( i=0; i<dt->n_groups; i++ ) {
		cffree(dt->groups[i]->name);
		cffree(dt->groups[i]);
	}

	cffree(dt->wavelength_from);
	cffree(dt->peak_list);
	cffree(dt->cnz_from);

	cffree(dt->panels);
	cffree(dt->bad);
	cffree(dt);
}


int data_template_file_to_panel_coords(const DataTemplate *dt,
                                       float *pfs, float *pss,
                                       int pn)
{
	*pfs = *pfs - dt->panels[pn].orig_min_fs;
	*pss = *pss - dt->panels[pn].orig_min_ss;
	return 0;
}


/**
 * Convert image-data-space fs/ss coordinates to panel-relative fs/ss
 * coordinates and panel number, assuming that the data is all in one slab.
 *
 * WARNING: This is probably not the routine you are looking for!
 *   If you use this routine, your code will only work with 'slabby' data, and
 *   will break for (amongst others) EuXFEL data.  Use
 *   data_template_file_to_panel_coords instead, and provide the panel number.
 *
 * \returns 0 on success, 1 on failure
 *
 */
int data_template_slabby_file_to_panel_coords(const DataTemplate *dt,
                                              float *pfs, float *pss, int *ppn)
{
	int p;
	int found = 0;

	for ( p=0; p<dt->n_panels; p++ ) {
		if ( (*pfs >= dt->panels[p].orig_min_fs)
		  && (*pfs < dt->panels[p].orig_max_fs+1)
		  && (*pss >= dt->panels[p].orig_min_ss)
		  && (*pss < dt->panels[p].orig_max_ss+1) )
		{
			if ( found ) {
				ERROR("Panel is ambiguous for fs,ss %f,%f\n");
				return 1;
			}
			*ppn = p;
			found = 1;
		}
	}

	if ( !found ) {
		ERROR("Couldn't find panel for fs,ss %f,%f\n", *pfs, *pss);
		return 1;
	}

	return data_template_file_to_panel_coords(dt, pfs, pss, *ppn);
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
	if ( dt->n_headers_to_copy >= MAX_COPY_HEADERS ) {
		ERROR("Too many extra headers to copy\n");
		return;
	}

	dt->headers_to_copy[dt->n_headers_to_copy++] = cfstrdup(header);
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


double convert_to_m(double val, int units)
{
	switch ( units ) {

		case WAVELENGTH_M :
		return val;

		case WAVELENGTH_A :
		return val * 1e-10;

		case WAVELENGTH_PHOTON_EV :
		return ph_eV_to_lambda(val);

		case WAVELENGTH_PHOTON_KEV :
		return ph_eV_to_lambda(val*1e3);

		case WAVELENGTH_ELECTRON_V :
		return el_V_to_lambda(val);

		case WAVELENGTH_ELECTRON_KV :
		return el_V_to_lambda(val*1e3);

	}

	return NAN;
}


/**
 * Get the wavelength from a DataTemplate, if possible.
 *
 * WARNING: This is probably not the routine you are looking for!
 *  See the disclaimer for image_create_for_simulation(), which applies
 *  equally to this routine.
 *
 * \returns the wavelength, in metres, or NAN if impossible.
 */
double data_template_get_wavelength_if_possible(const DataTemplate *dt)
{
	float val;
	char *rval;

	if ( dt->wavelength_from == NULL ) return NAN;

	val = strtod(dt->wavelength_from, &rval);
	if ( (*rval == '\0') && (rval != dt->wavelength_from) ) {
		return convert_to_m(val, dt->wavelength_unit);
	} else {
		return NAN;
	}
}


static int separate_value_and_units(const char *from,
                                    char **pvalue,
                                    char **punits)
{
	char *sp;
	char *fromcpy;
	char *unitscpy;

	if ( from == NULL ) return 1;

	fromcpy = cfstrdup(from);
	if ( fromcpy == NULL ) return 1;

	sp = strchr(fromcpy, ' ');
	if ( sp == NULL ) {
		unitscpy = NULL;
	} else {
		unitscpy = cfstrdup(sp+1);
		sp[0] = '\0';
	}

	*pvalue = fromcpy;
	*punits = unitscpy;
	return 0;
}


/* default_scale is a value to be used if both of the following
 * conditions are met:
 *
 *  1. The value is a reference to image headers/metadata,
 *      rather than a literal number.
 *  2. No units are specified in the number.
 *
 * This is totally horrible.  Sorry.  Blame history.
 */
static int im_get_length(struct image *image, const char *from,
                         double default_scale, double *pval)
{
	char *value_str;
	char *units;

	if ( from == NULL ) return 1;

	if ( separate_value_and_units(from, &value_str, &units) ) return 1;

	if ( units == NULL ) {

		/* No units given */

		if ( convert_float(value_str, pval) == 0 ) {

			/* Literal value with no units */
			cffree(value_str);
			return 0;

		} else {

			int r;
			r = image_read_header_float(image, value_str, pval);
			cffree(value_str);

			if ( r == 0 ) {
				/* Value read from headers with no units */
				*pval *= default_scale;
				return 0;
			} else {
				/* Failed to read value from headers */
				return 1;
			}
		}

	} else {

		/* Units are specified */

		double scale;

		if ( strcmp(units, "mm") == 0 ) {
			scale = 1e-3;
		} else if ( strcmp(units, "m") == 0 ) {
			scale = 1.0;
		} else {
			ERROR("Invalid length unit '%s'\n", units);
			cffree(value_str);
			cffree(units);
			return 1;
		}

		if ( convert_float(value_str, pval) == 0 ) {

			/* Literal value, units specified */
			cffree(value_str);
			cffree(units);
			*pval *= scale;
			return 0;

		} else {

			int r;
			r = image_read_header_float(image, value_str, pval);
			cffree(value_str);

			if ( r == 0 ) {
				/* Value read from headers, units specified */
				*pval *= scale;
				return 0;
			} else {
				/* Failed to read value from headers */
				return 1;
			}
		}
	}
}


static int all_panels_same_coffset(const DataTemplate *dtempl)
{
	int i;
	double total;
	double mean;

	total = 0.0;
	for ( i=0; i<dtempl->n_panels; i++ ) {
		total += dtempl->panels[i].cnz_offset;
	}
	mean = total/dtempl->n_panels;

	for ( i=0; i<dtempl->n_panels; i++ ) {
		struct panel_template *p = &dtempl->panels[i];
		if ( fabs(dtempl->panels[i].cnz_offset - mean) > 10.0*p->pixel_pitch ) return 0;
	}

	return 1;
}


static int all_panels_perpendicular_to_beam(const DataTemplate *dtempl)
{
	int i;

	for ( i=0; i<dtempl->n_panels; i++ ) {
		double z_diff;
		struct panel_template *p = &dtempl->panels[i];
		z_diff = p->fsz*PANEL_WIDTH(p) + p->ssz*PANEL_HEIGHT(p);
		if ( z_diff > 10.0 ) return 0;
	}
	return 1;
}


static int detector_flat(const DataTemplate *dtempl)
{
	return all_panels_perpendicular_to_beam(dtempl)
	    && all_panels_same_coffset(dtempl);
}


static void add_dg_point(const struct detgeom_panel *p,
                         int fs, int ss,
                         double *tx, double *ty, double *tz)
{
	*tx += (p->cnx + fs*p->fsx + ss*p->ssx) * p->pixel_pitch;
	*ty += (p->cny + fs*p->fsy + ss*p->ssy) * p->pixel_pitch;
	*tz += (p->cnz + fs*p->fsz + ss*p->ssz) * p->pixel_pitch;
}


static struct detgeom_panel_group *walk_group(const DataTemplate *dtempl,
                                              struct panel_group_template *gt,
                                              struct detgeom *detgeom,
                                              int serial, int c_mul)
{
	struct detgeom_panel_group *gr;

	if ( gt == NULL ) return NULL;

	gr = cfmalloc(sizeof(struct detgeom_panel_group));
	if ( gr == NULL ) return NULL;

	gr->name = cfstrdup(gt->name);
	gr->n_children = gt->n_children;

	if ( gr->n_children == 0 ) {

		/* Leaf node */
		gr->children = NULL;
		gr->panel = detgeom_find_panel(detgeom, gr->name);
		if ( gr->panel == NULL ) {
			ERROR("Couldn't find panel %s for leaf group\n", gr->name);
			return NULL;
		}
		gr->panel->group = gr;

		/* Calculate and make a note of the panel center */
		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		add_dg_point(gr->panel, 0, 0, &tx, &ty, &tz);
		add_dg_point(gr->panel, gr->panel->w, 0, &tx, &ty, &tz);
		add_dg_point(gr->panel, 0, gr->panel->h, &tx, &ty, &tz);
		add_dg_point(gr->panel, gr->panel->w, gr->panel->h, &tx, &ty, &tz);

		gr->cx = tx / 4.0;
		gr->cy = ty / 4.0;
		gr->cz = tz / 4.0;

	} else {

		int i;
		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		gr->panel = NULL;
		gr->children = cfmalloc(gt->n_children*sizeof(struct detgeom_panel_group *));
		if ( gr->children == NULL ) {
			cffree(gr);
			return NULL;
		}

		for ( i=0; i<gt->n_children; i++ ) {
			gr->children[i] = walk_group(dtempl, gt->children[i], detgeom,
			                             serial + c_mul*(i+1), c_mul*100);
			if ( gr->children[i] == NULL ) return NULL;
			gr->children[i]->parent = gr;
			tx += gr->children[i]->cx;
			ty += gr->children[i]->cy;
			tz += gr->children[i]->cz;
		}

		gr->cx = tx / gt->n_children;
		gr->cy = ty / gt->n_children;
		gr->cz = tz / gt->n_children;

	}

	gr->serial = serial;
	return gr;
}


struct detgeom *create_detgeom(struct image *image,
                               const DataTemplate *dtempl,
                               int two_d_only)
{
	struct detgeom *detgeom;
	int i;
	double clen;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	detgeom = cfmalloc(sizeof(struct detgeom));
	if ( detgeom == NULL ) return NULL;

	detgeom->top_group = NULL;

	detgeom->panels = cfmalloc(dtempl->n_panels*sizeof(struct detgeom_panel));
	if ( detgeom->panels == NULL ) {
		cffree(detgeom);
		return NULL;
	}

	detgeom->n_panels = dtempl->n_panels;

	if ( two_d_only ) {
		if ( !detector_flat(dtempl)
		  || (dtempl->shift_x_from != NULL)
		  || (dtempl->shift_y_from != NULL) )
		{
			cffree(detgeom->panels);
			cffree(detgeom);
			return NULL;
		}
	}

	if ( im_get_length(image, dtempl->cnz_from, 1e-3, &clen) )
	{
		if ( two_d_only ) {
			clen = NAN;
		} else {
			ERROR("Failed to read length from '%s'\n", dtempl->cnz_from);
			return NULL;
		}
	}


	for ( i=0; i<dtempl->n_panels; i++ ) {

		struct detgeom_panel *p = &detgeom->panels[i];
		struct panel_template *tmpl = &dtempl->panels[i];
		double shift_x, shift_y;

		p->name = safe_strdup(tmpl->name);

		p->pixel_pitch = tmpl->pixel_pitch;

		/* NB cnx,cny are in pixels, cnz is in m */
		p->cnx = tmpl->cnx;
		p->cny = tmpl->cny;

		/* Apply offset (in m) and then convert cnz from m to pixels */
		p->cnz = (clen + tmpl->cnz_offset) / p->pixel_pitch;

		/* Apply overall shift (already in m) */
		if ( dtempl->shift_x_from != NULL ) {
			if ( im_get_length(image, dtempl->shift_x_from, 1.0, &shift_x) ) {
				ERROR("Failed to read length from '%s'\n",
				      dtempl->shift_x_from);
				return NULL;
			}
			if ( im_get_length(image, dtempl->shift_y_from, 1.0, &shift_y) ) {
				ERROR("Failed to read length from '%s'\n",
				      dtempl->shift_y_from);
				return NULL;
			}
		} else {
			shift_x = 0.0;
			shift_y = 0.0;
		}

		if ( !isnan(shift_x) ) {
			p->cnx += shift_x / p->pixel_pitch;
		}
		if ( !isnan(shift_y) ) {
			p->cny += shift_y / p->pixel_pitch;
		}

		p->max_adu = tmpl->max_adu;

		switch ( tmpl->adu_scale_unit ) {

			case ADU_PER_PHOTON:
			p->adu_per_photon = tmpl->adu_scale;
			break;

			case ADU_PER_EV:
			if ( image == NULL ) {
				p->adu_per_photon = NAN;
			} else {
				p->adu_per_photon = tmpl->adu_scale
					* ph_lambda_to_eV(image->lambda);
			}
			break;

			default:
			p->adu_per_photon = 1.0;
			ERROR("Invalid ADU/ph scale unit (%i)\n",
			      tmpl->adu_scale_unit);
			break;

		}

		p->w = tmpl->orig_max_fs - tmpl->orig_min_fs + 1;
		p->h = tmpl->orig_max_ss - tmpl->orig_min_ss + 1;

		p->fsx = tmpl->fsx;
		p->fsy = tmpl->fsy;
		p->fsz = tmpl->fsz;
		p->ssx = tmpl->ssx;
		p->ssy = tmpl->ssy;
		p->ssz = tmpl->ssz;

	}

	detgeom->top_group = walk_group(dtempl, find_group(dtempl, "all"), detgeom, 0, 100);
	if ( detgeom->top_group != NULL ) {
		detgeom->top_group->parent = NULL;
	}

	return detgeom;
}


/**
 * Create a detgeom structure from the DataTemplate, if possible, and ignoring
 * 3D information.
 *
 * This procedure will create a detgeom structure provided that the detector
 *  is close to lying in a single flat plane perpendicular to the beam
 *  direction.  If certain things (e.g. panel z-positions) refer to headers,
 *  it might not be possible to determine that the detector is really flat
 *  until an image is loaded.  Therefore you must gracefully handle a NULL
 *  return value from this routine.
 *
 * \returns the detgeom structure, or NULL if impossible.
 */
struct detgeom *data_template_get_2d_detgeom_if_possible(const DataTemplate *dt)
{
	return create_detgeom(NULL, dt, 1);
}


/**
 * Returns the mean clen in m, or NAN in the following circumstances:
 * 1. If the individual panel distances vary by more than 10% of the average
 * 2. If the tilt of the panel creates a distance variation of more than 10%
 *    of the corner value over the extent of the panel
 * 3. If the detector geometry is not static (per-frame clen)
 *
 * \returns the mean camera length, or NAN if impossible.
 */
double data_template_get_clen_if_possible(const DataTemplate *dt)
{
	struct detgeom *dg;
	double clen;
	dg = data_template_get_2d_detgeom_if_possible(dt);
	if ( dg == NULL ) return NAN;
	clen = detgeom_mean_camera_length(dg);
	detgeom_free(dg);
	return clen;
}


static int translate_group_contents(DataTemplate *dtempl,
                                    const struct panel_group_template *group,
                                    double x, double y, double z,
                                    int is_metres)
{
	int i;

	if ( group->n_children == 0 ) {

		struct panel_template *p = find_panel_by_name(dtempl, group->name);
		if ( p == NULL ) return 1;

		if ( is_metres ) {
			p->cnx += x/p->pixel_pitch;
			p->cny += y/p->pixel_pitch;
			p->cnz_offset += z;
		} else {
			p->cnx += x;
			p->cny += y;
			p->cnz_offset += z*p->pixel_pitch;
		}

	} else {
		for ( i=0; i<group->n_children; i++ ) {
			translate_group_contents(dtempl, group->children[i],
			                         x, y, z, is_metres);
		}
	}

	return 0;
}


/**
 * Alters dtempl by shifting the named panel group by x,y,z in the CrystFEL
 * coordinate system.  x,y,z are in pixels, and all panels in the group must
 * have the same pixel size (but, this will not be checked).
 *
 * \returns zero for success, non-zero on error
 */
int data_template_translate_group_px(DataTemplate *dtempl, const char *group_name,
                                     double x, double y, double z)
{
	const struct panel_group_template *group = find_group(dtempl, group_name);
	if ( group == NULL ) return 1;
	return translate_group_contents(dtempl, group, x, y, z, 0);
}


/**
 * Alters dtempl by shifting the named panel group by x,y,z in the CrystFEL
 * coordinate system.  x,y,z are in metres.
 *
 * \returns zero for success, non-zero on error
 */
int data_template_translate_group_m(DataTemplate *dtempl, const char *group_name,
                                    double x, double y, double z)
{
	const struct panel_group_template *group = find_group(dtempl, group_name);
	if ( group == NULL ) return 1;
	return translate_group_contents(dtempl, group, x, y, z, 1);
}


static void add_point(const struct panel_template *p,
                      int fs, int ss,
                      double *tx, double *ty, double *tz)
{
	*tx += (p->cnx + fs*p->fsx + ss*p->ssx) * p->pixel_pitch;
	*ty += (p->cny + fs*p->fsy + ss*p->ssy) * p->pixel_pitch;
	*tz += p->cnz_offset + (fs*p->fsz + ss*p->ssz) * p->pixel_pitch;
}


static int group_center(DataTemplate *dtempl,
                        const struct panel_group_template *group,
                        double *cx, double *cy, double *cz)
{
	if ( group->n_children == 0 ) {

		const struct panel_template *p = find_panel_by_name(dtempl, group->name);
		if ( p == NULL ) return 1;

		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		add_point(p, 0, 0, &tx, &ty, &tz);
		add_point(p, PANEL_WIDTH(p), 0, &tx, &ty, &tz);
		add_point(p, 0, PANEL_HEIGHT(p), &tx, &ty, &tz);
		add_point(p, PANEL_WIDTH(p), PANEL_HEIGHT(p), &tx, &ty, &tz);

		*cx = tx / 4.0;
		*cy = ty / 4.0;
		*cz = tz / 4.0;

		return 0;

	} else {

		int i;
		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		for ( i=0; i<group->n_children; i++ ) {
			double gcx, gcy, gcz;
			group_center(dtempl, group->children[i], &gcx, &gcy, &gcz);
			tx += gcx;
			ty += gcy;
			tz += gcz;
		}

		*cx = tx / group->n_children;
		*cy = ty / group->n_children;
		*cz = tz / group->n_children;

		return 0;

	}
}


static int rotate_all_panels(DataTemplate *dtempl,
                             struct panel_group_template *group,
                             char axis, double ang,
                             double cx, double cy, double cz)
{
	if ( group->n_children == 0 ) {

		double cnz_px;
		struct panel_template *p = find_panel_by_name(dtempl, group->name);
		if ( p == NULL ) return 1;

		cx /= p->pixel_pitch;
		cy /= p->pixel_pitch;
		cz /= p->pixel_pitch;
		cnz_px = p->cnz_offset / p->pixel_pitch;

		switch ( axis )
		{
			case 'x':
			rotate2d(&p->cny, &cnz_px, cy, cz, ang);
			rotate2d(&p->fsy, &p->fsz, 0, 0, ang);
			rotate2d(&p->ssy, &p->ssz, 0, 0, ang);
			p->cnz_offset = cnz_px * p->pixel_pitch;
			break;

			case 'y':
			rotate2d(&cnz_px, &p->cnx, cz, cx, ang);
			rotate2d(&p->fsz, &p->fsx, 0, 0, ang);
			rotate2d(&p->ssz, &p->ssx, 0, 0, ang);
			p->cnz_offset = cnz_px * p->pixel_pitch;
			break;

			case 'z':
			rotate2d(&p->cnx, &p->cny, cx, cy, ang);
			rotate2d(&p->fsx, &p->fsy, 0, 0, ang);
			rotate2d(&p->ssx, &p->ssy, 0, 0, ang);
			break;

			default:
			ERROR("Invalid rotation axis '%c'\n", axis);
			return 1;
		}

		return 0;

	} else {

		int i;

		for ( i=0; i<group->n_children; i++ ) {
			rotate_all_panels(dtempl, group->children[i],
			                  axis, ang, cx, cy, cz);
		}

		return 0;

	}
}

/**
 * Alters dtempl by rotating the named panel group by ang (radians) about the
 * specified axis (char 'x', 'y' or 'z'), around the center of the group.
 *
 * \returns zero for success, non-zero on error
 */
int data_template_rotate_group(DataTemplate *dtempl, const char *group_name,
                               double ang, char axis)
{
	struct panel_group_template *group;
	double cx, cy, cz;

	group = find_group(dtempl, group_name);
	if ( group == NULL ) return 1;

	if ( group_center(dtempl, group, &cx, &cy, &cz) ) return 1;

	return rotate_all_panels(dtempl, group, axis, ang, cx, cy, cz);
}


static const char *str_dim(int dim)
{
	switch ( dim ) {
		case DIM_FS: return "fs";
		case DIM_SS: return "ss";
		case DIM_PLACEHOLDER: return "%";
		default: return NULL;
	}
}


int data_template_write_to_fh(const DataTemplate *dtempl, FILE *fh)
{
	int i;

	/* Basic top-level parameters */
	switch ( dtempl->wavelength_unit ) {

		case WAVELENGTH_M:
		fprintf(fh, "wavelength = %s m\n", dtempl->wavelength_from);
		break;

		case WAVELENGTH_A:
		fprintf(fh, "wavelength = %s A\n", dtempl->wavelength_from);
		break;

		case WAVELENGTH_ELECTRON_KV:
		fprintf(fh, "electron_voltage = %s kV\n", dtempl->wavelength_from);
		break;

		case WAVELENGTH_ELECTRON_V:
		fprintf(fh, "electron_voltage = %s V\n", dtempl->wavelength_from);
		break;

		case WAVELENGTH_PHOTON_KEV:
		fprintf(fh, "photon_energy = %s keV\n", dtempl->wavelength_from);
		break;

		case WAVELENGTH_PHOTON_EV:
		fprintf(fh, "photon_energy = %s eV\n", dtempl->wavelength_from);
		break;

		default:
		ERROR("Unknown wavelength unit (%i)\n", dtempl->wavelength_unit);
		return 1;

	}

	fprintf(fh, "clen = %s\n", dtempl->cnz_from);

	if ( dtempl->peak_list != NULL ) {
		fprintf(fh, "peak_list = %s\n", dtempl->peak_list);
	}
	switch ( dtempl->peak_list_type ) {
		case PEAK_LIST_AUTO:
		break;

		case PEAK_LIST_CXI:
		fprintf(fh, "peak_list_type = cxi\n");
		break;

		case PEAK_LIST_LIST3:
		fprintf(fh, "peak_list_type = list3\n");
		break;

		default:
		ERROR("Unknown peak list type (%i)\n", dtempl->peak_list_type);
		return 1;
	}

	fprintf(fh, "bandwidth = %e\n", dtempl->bandwidth);

	if ( dtempl->shift_x_from != NULL ) {
		fprintf(fh, "detector_shift_x = %s\n", dtempl->shift_x_from);
	}
	if ( dtempl->shift_y_from != NULL ) {
		fprintf(fh, "detector_shift_y = %s\n", dtempl->shift_y_from);
	}

	/* Other top-levels */
	int mask_done[MAX_MASKS] = {0};
	int satmap_done = 0;
	int satmap_file_done = 0;
	int mask_edge_pixels_done = 0;
	int pixel_pitch_done = 0;
	int adu_scale_done = 0;
	int max_adu_done = 0;
	int flag_values_done = 0;
	int data_done = 0;
	int dims_done[MAX_DIMS] = {0};
	for ( i=0; i<dtempl->n_panels; i++ ) {

		const struct panel_template *p = &dtempl->panels[i];
		int j;

		for ( j=0; j<MAX_MASKS; j++ ) {
			if ( p->masks[j].data_location == NULL ) continue;
			if ( !p->masks[j].mask_default ) continue;
			if ( mask_done[j] ) continue;
			fprintf(fh, "mask%i_data = %s\n",
			        j, p->masks[j].data_location);
			if ( p->masks[j].filename != NULL ) {
				fprintf(fh, "mask%i_file = %s\n",
				        j, p->masks[j].filename);
			}
			fprintf(fh, "mask%i_goodbits = 0x%x\n",
			        j, p->masks[j].good_bits);
			fprintf(fh, "mask%i_badbits = 0x%x\n",
			        j, p->masks[j].bad_bits);
			mask_done[j] = 1;
		}

		if ( p->satmap_default && !satmap_done && (p->satmap != NULL) ) {
			fprintf(fh, "saturation_map = %s\n", p->satmap);
			satmap_done = 1;
		}

		if ( p->satmap_file_default && !satmap_file_done && (p->satmap_file != NULL) ) {
			fprintf(fh, "saturation_map_file = %s\n", p->satmap);
			satmap_file_done = 1;
		}

		if ( p->mask_edge_pixels_default && !mask_edge_pixels_done && (p->mask_edge_pixels != 0) ) {
			fprintf(fh, "mask_edge_pixels = %i\n", p->mask_edge_pixels);
			mask_edge_pixels_done = 1;
		}

		if ( p->pixel_pitch_default && !pixel_pitch_done ) {
			fprintf(fh, "res = %f\n", 1.0/p->pixel_pitch);
			pixel_pitch_done = 1;
		}

		if ( p->max_adu_default && !max_adu_done && !isinf(p->max_adu) ) {
			fprintf(fh, "max_adu = %f\n", p->max_adu);
			max_adu_done = 1;
		}

		if ( p->data_default && !data_done ) {
			fprintf(fh, "data = %s\n", p->data);
			data_done = 1;
		}

		if ( p->flag_values_default && !flag_values_done ) {
			for ( j=0; j<MAX_FLAG_VALUES; j++ ) {
				switch ( p->flag_types[j] ) {
					case FLAG_NOTHING :
					break;

					case FLAG_EQUAL:
					fprintf(fh, "flag_equal = %i\n",
					        p->flag_values[j]);
					break;

					case FLAG_MORETHAN:
					fprintf(fh, "flag_morethan = %i\n",
					        p->flag_values[j]);
					break;

					case FLAG_LESSTHAN:
					fprintf(fh, "flag_lessthan = %i\n",
					        p->flag_values[j]);
					break;
				}
			}
			flag_values_done = 1;
		}

		if ( p->adu_scale_default && !adu_scale_done ) {
			switch ( p->adu_scale_unit ) {

				case ADU_PER_EV:
				fprintf(fh, "adu_per_eV = %f\n", p->adu_scale);
				break;

				case ADU_PER_PHOTON:
				fprintf(fh, "adu_per_photon = %f\n", p->adu_scale);
				break;
			}
			adu_scale_done = 1;
		}

		for ( j=0; j<MAX_DIMS; j++ ) {
			if ( p->dims_default[j] && !dims_done[j] && p->dims[j] != DIM_UNDEFINED ) {
				if ( p->dims[j] < 0 ) {
					fprintf(fh, "dim%i = %s\n", j, str_dim(p->dims[j]));
				} else {
					fprintf(fh, "dim%i = %i\n", j, p->dims[j]);
				}
				dims_done[j] = 1;
			}
		}
	}

	fprintf(fh, "\n");

	/* Bad regions */
	for ( i=0; i<dtempl->n_bad; i++ ) {
		const struct dt_badregion *bad = &dtempl->bad[i];
		assert(strncmp(bad->name, "bad", 3) == 0);
		if ( bad->is_fsss ) {
			const struct panel_template *p = &dtempl->panels[bad->panel_number];
			fprintf(fh, "%s/panel = %s\n", bad->name, p->name);
			fprintf(fh, "%s/min_fs = %i\n", bad->name, bad->min_fs+p->orig_min_fs);
			fprintf(fh, "%s/max_fs = %i\n", bad->name, bad->max_fs+p->orig_min_fs);
			fprintf(fh, "%s/min_ss = %i\n", bad->name, bad->min_ss+p->orig_min_ss);
			fprintf(fh, "%s/max_ss = %i\n", bad->name, bad->max_ss+p->orig_min_ss);
		} else {
			fprintf(fh, "%s/min_x = %f\n", bad->name, bad->min_x);
			fprintf(fh, "%s/max_x = %f\n", bad->name, bad->max_x);
			fprintf(fh, "%s/min_y = %f\n", bad->name, bad->min_y);
			fprintf(fh, "%s/max_y = %f\n", bad->name, bad->max_y);
		}
		fprintf(fh, "\n");
	}

	/* Panels */
	for ( i=0; i<dtempl->n_panels; i++ ) {

		int j;
		const struct panel_template *p = &dtempl->panels[i];

		fprintf(fh, "%s/min_fs = %i\n", p->name, p->orig_min_fs);
		fprintf(fh, "%s/max_fs = %i\n", p->name, p->orig_max_fs);
		fprintf(fh, "%s/min_ss = %i\n", p->name, p->orig_min_ss);
		fprintf(fh, "%s/max_ss = %i\n", p->name, p->orig_max_ss);
		fprintf(fh, "%s/corner_x = %f\n", p->name, p->cnx);
		fprintf(fh, "%s/corner_y = %f\n", p->name, p->cny);
		fprintf(fh, "%s/fs = %fx %+fy %+fz\n", p->name,
		        p->fsx, p->fsy, p->fsz);
		fprintf(fh, "%s/ss = %fx %+fy %+fz\n", p->name,
		        p->ssx, p->ssy, p->ssz);

		fprintf(fh, "%s/coffset = %f\n", p->name, p->cnz_offset);

		for ( j=0; j<MAX_MASKS; j++ ) {
			if ( p->masks[j].data_location == NULL ) continue;
			if ( p->masks[j].mask_default ) continue;
			fprintf(fh, "%s/mask%i_data = %s\n",
			        p->name, j, p->masks[j].data_location);
			if ( p->masks[j].filename != NULL ) {
				fprintf(fh, "%smask%i_file = %s\n",
				        p->name, j, p->masks[j].filename);
			}
			fprintf(fh, "%s/mask%i_goodbits = 0x%x\n",
			        p->name, j, p->masks[j].good_bits);
			fprintf(fh, "%s/mask%i_badbits = 0x%x\n",
			        p->name, j, p->masks[j].bad_bits);
		}

		if ( !p->satmap_default && (p->satmap != NULL) ) {
			fprintf(fh, "%s/saturation_map = %s\n", p->name, p->satmap);
		}

		if ( !p->satmap_file_default && (p->satmap_file != NULL) ) {
			fprintf(fh, "%s/saturation_map_file = %s\n", p->name, p->satmap_file);
		}

		if ( !p->mask_edge_pixels_default && (p->mask_edge_pixels != 0) ) {
			fprintf(fh, "%s/mask_edge_pixels = %i\n", p->name, p->mask_edge_pixels);
		}

		if ( !p->pixel_pitch_default ) {
			fprintf(fh, "%s/res = %f\n", p->name, 1.0/p->pixel_pitch);
		}

		if ( !p->adu_scale_default ) {
			switch ( p->adu_scale_unit ) {

				case ADU_PER_EV:
				fprintf(fh, "%s/adu_per_eV = %f\n", p->name, p->adu_scale);
				break;

				case ADU_PER_PHOTON:
				fprintf(fh, "%s/adu_per_photon = %f\n", p->name, p->adu_scale);
				break;
			}
		}

		if ( !p->max_adu_default ) {
			fprintf(fh, "%s/max_adu = %f\n", p->name, p->max_adu);
		}

		if ( !p->flag_values_default ) {
			for ( j=0; j<MAX_FLAG_VALUES; j++ ) {
				switch ( p->flag_types[j] ) {
					case FLAG_NOTHING :
					break;

					case FLAG_EQUAL:
					fprintf(fh, "%s/flag_equal = %i\n",
					        p->name, p->flag_values[j]);
					break;

					case FLAG_MORETHAN:
					fprintf(fh, "%s/flag_morethan = %i\n",
					        p->name, p->flag_values[j]);
					break;

					case FLAG_LESSTHAN:
					fprintf(fh, "%s/flag_lessthan = %i\n",
					        p->name, p->flag_values[j]);
					break;
				}
			}
		}

		if ( !p->data_default ) {
			fprintf(fh, "%s/data = %s\n", p->name, p->data);
		}

		for ( j=0; j<MAX_DIMS; j++ ) {
			if ( !p->dims_default[j] && (p->dims[j] != DIM_UNDEFINED) ) {
				if ( p->dims[j] < 0 ) {
					fprintf(fh, "%s/dim%i = %s\n", p->name, j, str_dim(p->dims[j]));
				} else {
					fprintf(fh, "%s/dim%i = %i\n", p->name, j, p->dims[j]);
				}
				dims_done[j] = 1;
			}
		}

		if ( p->bad ) {
			fprintf(fh, "%s/no_index = 1\n", p->name);
		}

		fprintf(fh, "\n");
	}

	/* Groups */
	for ( i=0; i<dtempl->n_groups; i++ ) {
		int j;
		if ( dtempl->groups[i]->n_children == 0 ) continue;
		fprintf(fh, "group_%s = ", dtempl->groups[i]->name);
		for ( j=0; j<dtempl->groups[i]->n_children; j++ ) {
			if ( j > 0 ) fprintf(fh, ",");
			fprintf(fh, "%s", dtempl->groups[i]->children[j]->name);
		}
		fprintf(fh, "\n");
	}

	return 0;
}


int data_template_write_to_file(const DataTemplate *dtempl, const char *filename)
{
	FILE *fh;
	int r;
	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;
	r = data_template_write_to_fh(dtempl, fh);
	fclose(fh);
	return r;
}


static void add_group_info(struct dg_group_info *ginfo, int *ppos,
                           struct panel_group_template *group,
                           int serial, int level, int c_mul)
{
	int j;
	int i = *ppos;
	(*ppos)++;

	ginfo[i].name = group->name;
	ginfo[i].serial = serial;
	ginfo[i].hierarchy_level = level;

	for ( j=0; j<group->n_children; j++ ) {
		add_group_info(ginfo, ppos, group->children[j],
		               serial+c_mul*(j+1), level+1, c_mul*100);
	}
}


struct dg_group_info *data_template_group_info(const DataTemplate *dtempl, int *n)
{
	struct dg_group_info *ginfo;
	int i;
	struct panel_group_template *group;

	ginfo = cfmalloc(sizeof(struct dg_group_info)*dtempl->n_groups);
	if ( ginfo == NULL ) return NULL;

	group = find_group(dtempl, "all");
	i = 0;
	add_group_info(ginfo, &i, group, 0, 0, 100);

	*n = dtempl->n_groups;
	return ginfo;
}
