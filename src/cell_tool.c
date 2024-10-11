/*
 * cell_tool.c
 *
 * Unit cell tool
 *
 * Copyright © 2018-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2020 Thomas White <taw@physics.org>
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>

#include <cell.h>
#include <cell-utils.h>
#include <reflist-utils.h>
#include <reflist.h>

#include "version.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Unit cell manipulation tool.\n"
"\n"
" -h, --help                 Display this help message.\n"
" -p, --pdb=<file>           Get unit cell from <file> (PDB or CrystFEL format).\n"
" -o <file>                  Output unit cell file.\n"
"     --version              Print CrystFEL version number and exit.\n"
"\n"
"  Actions:\n"
"     --find-ambi            Find indexing ambiguities for the cell.\n"
"     --uncenter             Calculate a primitive cell.\n"
"     --rings                Calculate powder ring positions.\n"
"     --compare-cell <file>  Compare unit cell with cell from <file>.\n"
"     --cell-choices         Calculate all three cell choices for monoclinic C cell.\n"
"     --transform=<op>       Transform unit cell.\n"
"\n"
" -y <pointgroup>            Real point group of the structure.\n"
"     --tolerance=<tol>      Set the tolerances for cell comparison.\n"
"                             Default: 5,1.5 (axis percentage, angle deg).\n"
"     --highres=n            Resolution limit (Angstroms) for --rings\n"
);
}


static int comparecells(UnitCell *cell, const char *comparecell,
                        double ltl, double atl)
{
	UnitCell *cell2;
	RationalMatrix *m;
	double tolerance[6];

	cell2 = load_cell_from_file(comparecell);
	if ( cell2 == NULL ) {
		ERROR("Failed to load unit cell from '%s'\n", comparecell);
		return 1;
	}
	if ( validate_cell(cell2) > 1 ) {
		ERROR("Comparison cell is invalid.\n");
		return 1;
	}
	STATUS("------------------> The reference unit cell:\n");
	cell_print(cell2);

	tolerance[0] = ltl;
	tolerance[1] = ltl;
	tolerance[2] = ltl;
	tolerance[3] = atl;
	tolerance[4] = atl;
	tolerance[5] = atl;

	STATUS("------------------> Reindexed (strictly the same lattice):\n");
	STATUS("Tolerances applied directly to the unit cells\n");
	if ( !compare_reindexed_cell_parameters(cell, cell2, tolerance, &m) ) {
		STATUS("No relationship found between lattices.\n");
	} else {
		UnitCell *trans;
		STATUS("Relationship found.  To become similar to the reference"
		       " cell, the input cell should be transformed by:\n");
		rtnl_mtx_print(m);
		STATUS("Transformed version of input unit cell:\n");
		trans = cell_transform_rational(cell, m);
		cell_print(trans);
		cell_free(trans);
		STATUS("NB transformed cell might not really be triclinic, "
		       "it's just that I don't (yet) know how to work out what "
		       "it is.\n");

	}

	STATUS("------------------> Derivative lattice  "
	       "(strictly the same lattice):\n");
	STATUS("Tolerances applied to primitive versions of the unit cells\n");
	if ( !compare_derivative_cell_parameters(cell, cell2, tolerance, 0, &m) ) {
		STATUS("No relationship found between lattices.\n");
	} else {
		UnitCell *trans;
		STATUS("Relationship found.  To become similar to the reference"
		       " cell, the input cell should be transformed by:\n");
		rtnl_mtx_print(m);
		STATUS("Transformed version of input unit cell:\n");
		trans = cell_transform_rational(cell, m);
		cell_print(trans);
		cell_free(trans);
		STATUS("NB transformed cell might not really be triclinic, "
		       "it's just that I don't (yet) know how to work out what "
		       "it is.\n");
		rtnl_mtx_free(m);
	}

	STATUS("------------------> Coincidence site lattice "
	       "(not strictly the same lattice):\n");
	STATUS("Tolerances applied to primitive versions of the unit cells\n");
	if ( !compare_derivative_cell_parameters(cell, cell2, tolerance, 1, &m) ) {
		STATUS("No relationship found between lattices.\n");
		return 0;
	} else {
		UnitCell *trans;
		STATUS("Relationship found.  To become similar to the reference"
		       " cell, the input cell should be transformed by:\n");
		rtnl_mtx_print(m);
		STATUS("Transformed version of input unit cell:\n");
		trans = cell_transform_rational(cell, m);
		cell_print(trans);
		cell_free(trans);
		STATUS("NB transformed cell might not really be triclinic, "
		       "it's just that I don't (yet) know how to work out what "
		       "it is.\n");
		rtnl_mtx_free(m);
	}

	return 0;
}


static int all_rings(UnitCell *cell, SymOpList *sym, double mres)
{
	struct powder_ring *rings;
	int i, n;

	rings = powder_rings(cell, sym, mres, &n);

	STATUS("\nAll powder rings up to %f Ångstrøms.\n", 1e+10/mres);
	STATUS("Note that screw axis or glide plane absences are not "
	       "omitted from this list.\n");
	STATUS("\n   d (Å)   1/d (m^-1)    h    k    l    multiplicity\n");
	STATUS("------------------------------------------------------\n");
	for ( i=0; i<n; i++ ) {
		printf("%10.3f %10.3e %4i %4i %4i    m = %i\n",
		       1e10/rings[i].resolution, rings[i].resolution,
		       rings[i].h, rings[i].k, rings[i].l,
		       rings[i].multi);
	}

	free(rings);

	return 0;
}


static int find_ambi(UnitCell *cell, SymOpList *sym, double ltl, double atl)
{
	SymOpList *amb;
	SymOpList *ops;
	signed int i[9];
	const int maxorder = 3;
	double tolerance[6];

	tolerance[0] = ltl;
	tolerance[1] = ltl;
	tolerance[2] = ltl;
	tolerance[3] = atl;
	tolerance[4] = atl;
	tolerance[5] = atl;

	ops = get_pointgroup("1");
	if ( ops == NULL ) return 1;
	set_symmetry_name(ops, "Observed");

	if ( sym == NULL ) {
		sym = get_pointgroup("1");
	}

	STATUS("Looking for ambiguities up to %ix each lattice length.\n", maxorder);
	STATUS("This will take about 30 seconds.  Please wait...\n");

	for ( i[0]=-maxorder; i[0]<=+maxorder; i[0]++ ) {
	for ( i[1]=-maxorder; i[1]<=+maxorder; i[1]++ ) {
	for ( i[2]=-maxorder; i[2]<=+maxorder; i[2]++ ) {
	for ( i[3]=-maxorder; i[3]<=+maxorder; i[3]++ ) {
	for ( i[4]=-maxorder; i[4]<=+maxorder; i[4]++ ) {
	for ( i[5]=-maxorder; i[5]<=+maxorder; i[5]++ ) {
	for ( i[6]=-maxorder; i[6]<=+maxorder; i[6]++ ) {
	for ( i[7]=-maxorder; i[7]<=+maxorder; i[7]++ ) {
	for ( i[8]=-maxorder; i[8]<=+maxorder; i[8]++ ) {

		UnitCell *nc;
		IntegerMatrix *m;
		int j, k;
		int l = 0;

		m = intmat_new(3, 3);
		for ( j=0; j<3; j++ ) {
		for ( k=0; k<3; k++ ) {
			intmat_set(m, j, k, i[l++]);
		}
		}

		if ( intmat_det(m) != +1 ) {
			intmat_free(m);
			continue;
		}

		nc = cell_transform_intmat(cell, m);

		if ( compare_cell_parameters(cell, nc, tolerance) ) {
			if ( !intmat_is_identity(m) ) add_symop(ops, m);
			STATUS("-----------------------------------------------"
			       "-------------------------------------------\n");
			cell_print(nc);
			intmat_print(m);
		} else {
			intmat_free(m);
		}

		cell_free(nc);

	}
	}
	}
	}
	}
	}
	}
	}
	}

	STATUS("Observed symmetry operations:\n");
	describe_symmetry(ops);

	amb = get_ambiguities(ops, sym);
	if ( amb == NULL ) {
		STATUS("No ambiguities (or error calculating them)\n");
	} else {
		STATUS("Ambiguity operations:\n");
		describe_symmetry(amb);
		free_symoplist(amb);
	}

	free_symoplist(ops);

	return 0;
}


static int uncenter(UnitCell *cell, const char *out_file)
{
	UnitCell *cnew;
	IntegerMatrix *C;
	RationalMatrix *Ci;

	cnew = uncenter_cell(cell, &C, &Ci);

	STATUS("------------------> The primitive unit cell:\n");
	cell_print(cnew);

	STATUS("------------------> The centering transformation:\n");
	intmat_print(C);

	STATUS("------------------> The un-centering transformation:\n");
	rtnl_mtx_print(Ci);

	if ( out_file != NULL ) {
		FILE *fh = fopen(out_file, "w");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", out_file);
			return 1;
		}
		write_cell(cnew, fh);
		fclose(fh);
	}

	return 0;
}


static int transform(UnitCell *cell, const char *trans_str,
                     const char *out_file)
{
	RationalMatrix *trans;
	Rational det;
	UnitCell *nc;

	trans = parse_cell_transformation(trans_str);
	if ( trans == NULL ) {
		ERROR("Invalid cell transformation '%s'\n", trans_str);
		return 1;
	}

	nc = cell_transform_rational(cell, trans);

	STATUS("------------------> The transformation matrix:\n");
	rtnl_mtx_print(trans);
	det = rtnl_mtx_det(trans);
	STATUS("Determinant = %s\n", rtnl_format(det));
	if ( rtnl_cmp(det, rtnl_zero()) == 0 ) {
		ERROR("Singular transformation matrix - cannot transform.\n");
		return 1;
	}

	STATUS("------------------> The transformed unit cell:\n");
	cell_print(nc);
	STATUS("NB transformed cell might not really be triclinic, "
	       "it's just that I don't (yet) know how to work out what "
	       "it is.\n");

	if ( out_file != NULL ) {
		FILE *fh = fopen(out_file, "w");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", out_file);
			return 1;
		}
		write_cell(nc, fh);
		fclose(fh);
	}

	return 0;
}


static int cell_choices(UnitCell *cell)
{
	if ( cell_get_lattice_type(cell) != L_MONOCLINIC ) {
		ERROR("Cell must be monoclinic to use --cell-choices\n");
		return 1;
	}

	if ( cell_get_unique_axis(cell) == 'b' ) {
		transform(cell, "-a-c,b,a", NULL);
		transform(cell, "c,b,-a-c", NULL);
	} else {
		ERROR("Sorry, --cell-choices only supports unique axis b.\n");
		return 1;
	}

	return 0;
}


enum {
	CT_NOTHING,
	CT_FINDAMBI,
	CT_UNCENTER,
	CT_RINGS,
	CT_COMPARE,
	CT_CHOICES,
	CT_TRANSFORM,
};


int main(int argc, char *argv[])
{
	int c;
	char *cell_file = NULL;
	UnitCell *cell;
	char *toler = NULL;
	float ltl = 0.05;         /* fraction */
	float atl = deg2rad(1.5); /* radians */
	char *sym_str = NULL;
	SymOpList *sym = NULL;
	int mode = CT_NOTHING;
	char *comparecell = NULL;
	char *out_file = NULL;
	float highres;
	double rmax = 1/(2.0e-10);
	char *trans_str = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"pdb",                1, NULL,               'p'},
		{"tolerance",          1, NULL,                2},
		{"output",             1, NULL,               'o'},
		{"version",            0, NULL,                6},

		/* Modes of operation */
		{"find-ambi",          0, &mode,               CT_FINDAMBI},
		{"uncenter",           0, &mode,               CT_UNCENTER},
		{"uncentre",           0, &mode,               CT_UNCENTER},
		{"rings",              0, &mode,               CT_RINGS},
		{"compare-cell",       1, NULL,                3},
		{"cell-choices",       0, &mode,               CT_CHOICES},

		{"transform",          1, NULL,                4},
		{"highres",            1, NULL,                5},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hp:y:o:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'p' :
			cell_file = strdup(optarg);
			break;

			case 'o' :
			out_file = strdup(optarg);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 2 :
			toler = strdup(optarg);
			break;

			case 3 :
			comparecell = strdup(optarg);
			mode = CT_COMPARE;
			break;

			case 4 :
			trans_str = strdup(optarg);
			mode = CT_TRANSFORM;
			break;

			case 5 :
			if ( sscanf(optarg, "%e", &highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			rmax = 1.0 / (highres/1e10);
			break;

			case 6 :
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
			return 0;

			case 0 :
			break;

			default :
			return 1;

		}

	}

	/* If there's a parameter left over, we assume it's the unit cell */
	if ( (argc > optind) && (cell_file == NULL) ) {
		cell_file = strdup(argv[optind++]);
	}

	/* If there's STILL a parameter left over, complain*/
	if ( argc > optind ) {
		ERROR("Excess command-line arguments:\n");
		do {
			ERROR("'%s'\n", argv[optind++]);
		} while ( argc > optind );
		return 1;
	}

	if ( cell_file == NULL ) {
		ERROR("You must give a filename for the unit cell PDB file.\n");
		return 1;
	}
	STATUS("Input unit cell: %s\n", cell_file);
	cell = load_cell_from_file(cell_file);
	if ( cell == NULL ) {
		ERROR("Failed to load cell from '%s'\n", cell_file);
		return 1;
	}
	free(cell_file);

	if ( toler != NULL ) {
		int i;
		int ncomma = 0;
		size_t l = strlen(toler);
		for ( i=0; i<l; i++ ) if ( toler[i] == ',' ) ncomma++;
		if ( ncomma != 1 ) {
			ERROR("Invalid parameters for --tolerance.  "
			      "Should be: --tolerance=lengthtol,angtol "
			      "(percent,degrees)\n");
			return 1;
		}
		if ( sscanf(toler, "%f,%f", &ltl, &atl) != 2 ) {
			ERROR("Invalid parameters for --tolerance\n");
			return 1;
		}
		ltl /= 100.0;  /* percent to fraction */
		atl = deg2rad(atl);
		free(toler);
	}

	STATUS("------------------> The input unit cell:\n");
	cell_print(cell);

	if ( validate_cell(cell) > 1 ) {
		ERROR("Cell is invalid.\n");
		return 1;
	}

	if ( sym_str != NULL ) {
		sym = get_pointgroup(sym_str);
		if ( sym == NULL ) return 1;
		free(sym_str);
	}

	if ( mode == CT_NOTHING ) {
		ERROR("Please specify mode of operation (see --help)\n");
		return 1;
	}

	if ( mode == CT_FINDAMBI ) return find_ambi(cell, sym, ltl, atl);
	if ( mode == CT_UNCENTER ) return uncenter(cell, out_file);
	if ( mode == CT_RINGS ) return all_rings(cell, sym, rmax);
	if ( mode == CT_COMPARE ) return comparecells(cell, comparecell, ltl, atl);
	if ( mode == CT_TRANSFORM ) return transform(cell, trans_str, out_file);
	if ( mode == CT_CHOICES ) return cell_choices(cell);

	return 1;
}
