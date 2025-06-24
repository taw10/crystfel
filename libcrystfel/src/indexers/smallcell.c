/*
 * smallcell.c
 *
 * Re-implementation of graph theory indexing algorithm for small unit cells
 *  borrowed from cctbx.small_cell
 *
 * Copyright © 2024 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2024 Isabel Costello <isabel.costello@desy.de>
 *   2024 Thomas White <thomas.white@desy.de>
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>
#include <unistd.h>
#include <argp.h>

#include "image.h"
#include "index.h"
#include "cell-utils.h"
#include "symmetry.h"
#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "geometry.h"
#include "detgeom.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include "cell.h"
#include "taketwo.h"

/** \file smallcell.h */

#define MAX_NEIGH (50)
#define MAX_NODES (8072)
#define MAX_CLIQUES (128)


struct g_matrix
{
	double A;
	double B;
	double C;
	double D;
	double E;
	double F;
	double G;
	double H;
	double J;
};

struct smallcell_private
{
	struct powder_ring *powderrings;
	int num_rings;
	SymOpList *sym;
	struct g_matrix g9;
	UnitCell *template;
};


void *smallcell_prepare(IndexingMethod *indm, struct smallcell_options *opts,
                        UnitCell *cell)
{
	struct smallcell_private *dp;
	double asx, bsx, csx;
	double asy, bsy, csy;
	double asz, bsz, csz;

	dp = cfmalloc(sizeof(struct smallcell_private));
	if ( dp == NULL ) return NULL;

	dp->sym = get_lattice_symmetry(cell);
	dp->powderrings = powder_rings(cell, dp->sym, 1/1e-10, &dp->num_rings);
	dp->template = cell_new_from_cell(cell);

	/* Get reciprocal unit cell elements in order to create G* matrix
	 * and store in private */

	cell_get_reciprocal(cell,
	                    &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);

	dp->g9.A = (asx * asx) + (asy * asy) + (asz * asz);
	dp->g9.B = (asx * bsx) + (asy * bsy) + (asz * bsz);
	dp->g9.C = (asx * csx) + (asy * csy) + (asz * csz);
	dp->g9.D = (bsx * asx) + (bsy * asy) + (bsz * asz);
	dp->g9.E = (bsx * bsx) + (bsy * bsy) + (bsz * bsz);
	dp->g9.F = (bsx * csx) + (bsy * csy) + (bsz * csz);
	dp->g9.G = (csx * asx) + (csy * asy) + (csz * asz);
	dp->g9.H = (csx * bsx) + (csy * bsy) + (csz * bsz);
	dp->g9.J = (csx * csx) + (csy * csy) + (csz * csz);

	return dp;
}


/* PeakInfo structure for storing peak-related information for every match found */
typedef struct PeakInfo
{
	int peak_number;
	double peak_res;
	double x, y, z;
	int h, k, l;
	struct PeakInfo *neigh[MAX_NEIGH];
	int n_neigh;
	double weight_list[MAX_NEIGH];
} PeakInfo;


struct Nodelist
{
	int n_mem;
	struct PeakInfo *mem[MAX_NODES];
};


struct Cliquelist
{
	int n;
	struct Nodelist *list[MAX_CLIQUES];
};


static struct Nodelist *new_nodelist()
{
	struct Nodelist *r = cfmalloc(sizeof(struct Nodelist));
	if ( r == NULL ) return NULL;
	r->n_mem = 0;
	return r;
}


static struct Nodelist *CopyRlist(struct Nodelist *R)
{
	int i;
	struct Nodelist *Rcopy = new_nodelist();
	for ( i=0; i<R->n_mem; i++ ) {
		Rcopy->mem[i] = R->mem[i];
		Rcopy->n_mem++;
	}
	return Rcopy;
}


/* To make a list of neighbour nodes */
static struct Nodelist *neighbours(struct PeakInfo *in)
{
	int i;
	struct Nodelist *list = new_nodelist();
	for ( i=0; i<in->n_neigh; i++ ) {
		list->mem[i] = in->neigh[i];
	}
	list->n_mem = in->n_neigh;
	return list;
}


/* Function to see if a node is in a list */
static int isin(struct Nodelist *b, struct PeakInfo *test)
{
	int i;
	for ( i=0; i<b->n_mem; i++ ) {
		if ( test == b->mem[i] ) return 1;
	}
	return 0;
}


static void add(struct Nodelist *c, struct PeakInfo *test)
{
	if (isin(c, test)) return;
	assert(c->n_mem < MAX_NODES);
	c->mem[c->n_mem++] = test;
}


/* To create a new list of nodes from the union between two node lists
 *   i.e. c = a U b */
static struct Nodelist *Union(struct Nodelist *a, struct Nodelist *b)
{
	int i, j;
	struct Nodelist *c = new_nodelist();
	for ( i=0; i<a->n_mem; i++ ) {
		add(c, a->mem[i]);
	}
	for ( j=0; j<b->n_mem; j++ ) {
		add(c, b->mem[j]);
	}
	return c;
}


/* To create a new list of nodes from the intersection of two lists
 *  i.e. c = a intersection b */
static struct Nodelist *intersection(struct Nodelist *a, struct Nodelist *b)
{
	int j;
	struct Nodelist *c = new_nodelist();
	for ( j=0; j<a->n_mem; j++ ) {
		if ( isin(b, a->mem[j]) == 1 ) {
			add(c, a->mem[j]);
		}
	}
	return c;
}


/* To create a new list of nodes (c) from the exclusion of list b from list a
 *  i.e. c = a\b */
static struct Nodelist *exclusion(struct Nodelist *a, struct Nodelist *b)
{
	int j;
	struct Nodelist *c = new_nodelist();
	for ( j=0; j<a->n_mem; j++ ) {
		if ( isin(b, a->mem[j]) == 0 ) {
			add(c, a->mem[j]);
		}
	}
	return c;
}


/* Function to exclude a single node */
static struct Nodelist *exclunode(struct Nodelist *a, struct PeakInfo *v)
{
	int j;
	for ( j=0; j<a->n_mem; j++ ) {
		if (a->mem[j] == v) {
			int i;
			int hole_index = j;
			for ( i=hole_index; i<a->n_mem - 1; i++ ) {
				a->mem[i] = a->mem[i+1];
			}
		}
	}
	a->n_mem--;
	return a;
}


/* Function to append a node to a list (creating a new list) */
static struct Nodelist *append(struct Nodelist *a, struct PeakInfo *v)
{
	int j;
	struct Nodelist *c = new_nodelist();
	for ( j=0; j<a->n_mem; j++ ) {
		add(c, a->mem[j]);
	}

	add(c, v);
	return c;
}


/* Average weight function */
static double avg_weight(const struct PeakInfo *pk)
{
	int j;
	double avg_add = 0.0;
	for ( j=0; j<pk->n_neigh; j++ ) {
		avg_add += pk->weight_list[j];
	}
	return avg_add / pk->n_neigh;
}


static struct PeakInfo *find_pivot(struct Nodelist *piv_pool)
{
	int i;
	struct PeakInfo *piv = NULL;
	double min_weight = +INFINITY;
	int max_neigh = 0;

	for ( i=0; i<piv_pool->n_mem; i++ ) {

		if ( piv_pool->mem[i]->n_neigh > max_neigh ) {
			min_weight = avg_weight(piv_pool->mem[i]);
			piv = piv_pool->mem[i];

		} else if ( piv_pool->mem[i]->n_neigh == max_neigh ) {
			if ( avg_weight(piv_pool->mem[i]) < min_weight ) {
				min_weight = avg_weight(piv_pool->mem[i]);
				piv = piv_pool->mem[i];
			}
		}
	}
	return piv;
}


static int smallest_clique(struct Cliquelist *cliques)
{
	int i;
	int mini = 0;
	int minc = cliques->list[0]->n_mem;
	for ( i=1; i<cliques->n; i++ ) {
		if ( cliques->list[i]->n_mem < minc ) {
			minc = cliques->list[i]->n_mem;
			mini = i;
		}
	}
	return mini;
}


/* Bron-Kerbosch algorithm, with pivoting */
static void BK(struct Nodelist *R,
               struct Nodelist *P,
               struct Nodelist *X,
               struct Cliquelist *Max_cliques)
{
	/* If P & X are empty -> add R to the un-mapped clique array */
	if ( P->n_mem == 0 && X->n_mem == 0 ) {
		if ( Max_cliques->n >= MAX_CLIQUES ) {
			int mini = smallest_clique(Max_cliques);
			if ( R->n_mem > Max_cliques->list[mini]->n_mem ) {
				cffree(Max_cliques->list[mini]);
				Max_cliques->list[mini] = CopyRlist(R);
			}
		} else {
			Max_cliques->list[Max_cliques->n] = CopyRlist(R);
			Max_cliques->n++;
		}
		return;
	}

	/* Find pivot u from set of nodes in P U X */
	struct Nodelist *piv_pool = Union(P, X);
	struct PeakInfo *piv = find_pivot(piv_pool);
	cffree(piv_pool);
	if ( piv == NULL ) {
		ERROR("Couldn't find pivot\n");
		return;
	}

	/* Get list of neighbours of the pivot */
	struct Nodelist *piv_neighbours = neighbours(piv);

	/* Remove pivot neighbours from P list */
	struct Nodelist *P_excl = exclusion(P, piv_neighbours);

	cffree(piv_neighbours);

	int u;
	for ( u=0; u<P_excl->n_mem; u++ ) {

		/* Set up recursive call:
		 * BK(R <union> {v},
		 *    P <intersect> N(v),
		 *    X <intersect> N(v)
		 */
		struct PeakInfo *v = P_excl->mem[u];
		struct Nodelist *v_neighs = neighbours(v);
		struct Nodelist *R_new = append(R, v);
		struct Nodelist *P_new = intersection(P, v_neighs);
		struct Nodelist *X_new = intersection(X, v_neighs);
		cffree(v_neighs);

		BK(R_new, P_new, X_new, Max_cliques);

		cffree(R_new);
		cffree(P_new);
		cffree(X_new);

		/* Redefine P and X as P\v and X <union> {v} */
		exclunode(P, v);
		add(X, v);
	}
	cffree(P_excl);
	return;
}


static double calc_d2(struct PeakInfo a, struct PeakInfo b, struct g_matrix g9)
{
	/* d_2 = sqrt(Transpose(MI_b - MI_a).G*.(MI_b - MI_a)) */
	return sqrt((b.h - a.h) * (g9.A * (b.h - a.h) + g9.D * (b.k - a.k) + g9.G * (b.l - a.l))
	          + (b.k - a.k) * (g9.B * (b.h - a.h) + g9.E * (b.k - a.k) + g9.H * (b.l - a.l))
	          + (b.l - a.l) * (g9.C * (b.h - a.h) + g9.F * (b.k - a.k) + g9.J * (b.l - a.l)));
}


static struct PeakInfo *associate_to_rings(ImageFeatureList *peaks,
                                           struct powder_ring *rings,
                                           int n_rings,
                                           SymOpList *sym,
                                           struct detgeom *det,
                                           double lambda,
                                           int *pn_matches)
{
	const double tol = 6 / (lambda * det->panels[0].cnz);  /* 6 pixels at center */
	int peak_infos_size = 100;      /*  Arbitrary initial size allocation */
	PeakInfo *peak_infos = cfmalloc(peak_infos_size * sizeof(PeakInfo));
	int num_peak_infos = 0;
	int i;
	int npk = image_feature_count(peaks);
	SymOpMask *msk = new_symopmask(sym);

	/* Loop through each peak, calculate d, then 1/d value
	 * (based on estimate_peak_resolution from peak.c), then use match_rings
	 * function to create structs for this peak with all possible h,k,l values */
	for ( i=0; i<npk; i++ ) {

		struct imagefeature *f = image_get_feature(peaks, i);
		double r[3];
		detgeom_transform_coords(&det->panels[f->pn], f->fs, f->ss,
		                         lambda, 0.0, 0.0, r);
		double rns = modulus(r[0], r[1], r[2]);
		int j;

		for ( j=0; j<n_rings; j++ ) {

			if ( fabs(rings[j].resolution - rns) <= tol ) {

				signed int h = rings[j].h;
				signed int k = rings[j].k;
				signed int l = rings[j].l;

				/* Looking for symmetries and creating more
				 * PeakInfo structs with these symmetry indices */
				special_position(sym, msk, h, k, l);
				int n = num_equivs(sym, msk);
				int y;
				for ( y=0; y<n; y++ ) {

					if ( num_peak_infos >= peak_infos_size ) {
						peak_infos_size *= 2;
						peak_infos = cfrealloc(peak_infos,
						                       peak_infos_size * sizeof(PeakInfo));
					}

					signed int ha, ka, la;
					get_equiv(sym, msk, y, h, k, l, &ha, &ka, &la);

					/*  Add ha, ka, la to list of reflections  */
					peak_infos[num_peak_infos].peak_number = i;
					peak_infos[num_peak_infos].peak_res = rns;
					peak_infos[num_peak_infos].h = ha;
					peak_infos[num_peak_infos].k = ka;
					peak_infos[num_peak_infos].l = la;
					peak_infos[num_peak_infos].x = r[0];
					peak_infos[num_peak_infos].y = r[1];
					peak_infos[num_peak_infos].z = r[2];
					peak_infos[num_peak_infos].n_neigh = 0;
					num_peak_infos++;
				}
			}

		}
	}

	free_symopmask(msk);

	*pn_matches = num_peak_infos;
	return peak_infos;
}


static void link_nodes(struct PeakInfo *peak_infos, int num_peak_infos, struct g_matrix g9)
{
	const double dtol = 1e8;
	int j;

	/* Loop through peak numbers */
	for ( j=0; j<num_peak_infos; j++ ) {

		int y;

		/* Loop through the rest of the peak infos */
		for ( y=j+1; y<num_peak_infos-1; y++ ) {

			if ( peak_infos[y].peak_number == peak_infos[j].peak_number ) continue;

			/* Observed dist. d_1 = fabs(res.vec_j - res.vec_y) */
			double d_1 = modulus(peak_infos[j].x - peak_infos[y].x,
			                     peak_infos[j].y - peak_infos[y].y,
			                     peak_infos[j].z - peak_infos[y].z);

			/* Predicted d_2 */
			double d_2 = calc_d2(peak_infos[j], peak_infos[y], g9);

			double diff = fabs(d_2 - d_1);

			/* Now test difference */
			if ( diff <= dtol ) {

				/* Connect nodes */
				if ( peak_infos[j].n_neigh <= MAX_NEIGH ) {

					peak_infos[j].neigh[peak_infos[j].n_neigh] = &peak_infos[y];
					peak_infos[j].weight_list[peak_infos[j].n_neigh] = diff;
					peak_infos[j].n_neigh++;

				} else {
					ERROR("Too many neighbours.\n");
				}

				if ( peak_infos[y].n_neigh <= MAX_NEIGH ) {

					peak_infos[y].neigh[peak_infos[y].n_neigh] = &peak_infos[j];
					peak_infos[y].weight_list[peak_infos[y].n_neigh] = diff;
					peak_infos[y].n_neigh++;

				} else {
					ERROR("Too many neighbours.\n");
				}
			}
		}
	}
}


static struct Cliquelist *find_max_cliques(struct PeakInfo *peak_infos,
                                           int num_peak_infos)
{
	struct Nodelist *R;
	struct Nodelist *X;
	struct Nodelist *P;
	struct Cliquelist *Max_cliques;
	int i;

	/*  R: array of nodes forming a clique
	 *  P: array of all prospective nodes that are connected to R which
	 *     may be added to R. To begin, this is all nodes i.e all peak_infos
	 *  X: exculsion set (same form as R but nodes that are NOT candidates for
	 *     the max. clique, were originaly in P) */
	P = new_nodelist();

	/* To make P; create nodelist of all peak_infos */
	for ( i=0; i<num_peak_infos; i++ ) {
		if ( peak_infos[i].n_neigh != 0 ) {
			add(P, &peak_infos[i]);
		}
	}

	if ( P->n_mem <= 2 ) {
		cffree(P);
		return NULL;
	}

	R = new_nodelist();
	X = new_nodelist();

	Max_cliques = cfmalloc(sizeof(struct Cliquelist));
	Max_cliques->n = 0;

	BK(R, P, X, Max_cliques);

	cffree(X);
	cffree(P);
	cffree(R);

	return Max_cliques;
}


static int compare_cliques(const void *av, const void *bv)
{
	struct Nodelist * const *a = av;
	struct Nodelist * const *b = bv;
	if ( (*a)->n_mem < (*b)->n_mem ) return 1;
	if ( (*a)->n_mem > (*b)->n_mem ) return -1;
	return 0;
}


static UnitCell *fit_cell(struct Nodelist *clique)
{
	UnitCell *uc;
	gsl_matrix *h_mat;
	gsl_vector *h_vec;
	gsl_vector *cell_vec;
	double chisq;
	gsl_matrix *cov;
	gsl_multifit_linear_workspace *work;
	int i;
	int r;
	int have_a = 0, have_b = 0, have_c = 0;

	h_mat = gsl_matrix_calloc(3*clique->n_mem, 9);
	h_vec = gsl_vector_alloc(3*clique->n_mem);

	for ( i=0; i<clique->n_mem; i++) {

		struct PeakInfo *node = clique->mem[i];

		gsl_matrix_set(h_mat, 3*i, 0, node->h);
		gsl_matrix_set(h_mat, 3*i, 3, node->k);
		gsl_matrix_set(h_mat, 3*i, 6, node->l);
		gsl_vector_set(h_vec, 3*i, node->x);

		gsl_matrix_set(h_mat, 3*i+1, 1, node->h);
		gsl_matrix_set(h_mat, 3*i+1, 4, node->k);
		gsl_matrix_set(h_mat, 3*i+1, 7, node->l);
		gsl_vector_set(h_vec, 3*i+1, node->y);

		gsl_matrix_set(h_mat, 3*i+2, 2, node->h);
		gsl_matrix_set(h_mat, 3*i+2, 5, node->k);
		gsl_matrix_set(h_mat, 3*i+2, 8, node->l);
		gsl_vector_set(h_vec, 3*i+2, node->z);

		if ( node->h != 0 ) have_a = 1;
		if ( node->k != 0 ) have_b = 1;
		if ( node->l != 0 ) have_c = 1;

	}

	/* Must have all three axes represented */
	if ( !(have_a && have_b && have_c) ) {
		gsl_vector_free(h_vec);
		gsl_matrix_free(h_mat);
		return NULL;
	}

	cell_vec = gsl_vector_alloc(9);
	cov = gsl_matrix_alloc(3*clique->n_mem, 9);
	work = gsl_multifit_linear_alloc(3*clique->n_mem, 9);
	r = gsl_multifit_linear(h_mat, h_vec, cell_vec, cov, &chisq, work);

	gsl_multifit_linear_free(work);
	gsl_matrix_free(cov);
	gsl_vector_free(h_vec);
	gsl_matrix_free(h_mat);

	if ( r ) return NULL;

	/* cell_vec = [a*x a*y a*z b*x b*y b*z c*x c*y c*z]  */
	uc = cell_new();
	cell_set_reciprocal(uc,
	                    gsl_vector_get(cell_vec, 0),
	                    gsl_vector_get(cell_vec, 1),
	                    gsl_vector_get(cell_vec, 2),
	                    gsl_vector_get(cell_vec, 3),
	                    gsl_vector_get(cell_vec, 4),
	                    gsl_vector_get(cell_vec, 5),
	                    gsl_vector_get(cell_vec, 6),
	                    gsl_vector_get(cell_vec, 7),
	                    gsl_vector_get(cell_vec, 8));

	gsl_vector_free(cell_vec);

	return uc;
}


static void free_cliquelist(struct Cliquelist *cliques)
{
	int i;

	for ( i=0; i<cliques->n; i++ ) {
		cffree(cliques->list[i]);
	}

	cffree(cliques);
}


int smallcell_index(struct image *image, void *mpriv)
{
	struct PeakInfo *peak_infos;
	int num_peak_infos;
	int i;
	struct Cliquelist *cliques;
	struct smallcell_private *priv = (struct smallcell_private *)mpriv;

	peak_infos = associate_to_rings(image->features,
	                                priv->powderrings,
	                                priv->num_rings,
	                                priv->sym,
	                                image->detgeom,
	                                image->lambda,
					&num_peak_infos);

	link_nodes(peak_infos, num_peak_infos, priv->g9);

	cliques = find_max_cliques(peak_infos, num_peak_infos);
	if ( cliques == NULL ) {
		cffree(peak_infos);
		return 0;
	}

	/* Sort the list of cliques by number of members */
	qsort(cliques->list, cliques->n, sizeof(struct Nodelist *),
	      compare_cliques);

	/* Go down the list until we find an acceptable solution */
	for ( i=0; i<cliques->n; i++ ) {

		if ( cliques->list[i]->n_mem < 5 ) continue;

		UnitCell *uc = fit_cell(cliques->list[i]);
		if ( uc == NULL ) continue;
		cell_set_lattice_type(uc, cell_get_lattice_type(priv->template));
		cell_set_centering(uc, cell_get_centering(priv->template));
		cell_set_unique_axis(uc, cell_get_unique_axis(priv->template));

		if ( right_handed(uc) && !validate_cell(uc) ) {

			Crystal *cr = crystal_new();
			if ( cr == NULL ) {
				ERROR("Failed to allocate crystal.\n");
				continue;
			}
			crystal_set_cell(cr, uc);
			image_add_crystal(image, cr);
			free_cliquelist(cliques);
			cffree(peak_infos);
			return 1;
		}

		cell_free(uc);

	}

	free_cliquelist(cliques);
	cffree(peak_infos);
	return 0;
}


void smallcell_cleanup(void *mpriv)
{
	struct smallcell_private *dp = mpriv;
	cffree(dp->powderrings);
	free_symoplist(dp->sym);
	cffree(dp);
}


static void smallcell_show_help()
{
}


int smallcell_default_options(struct smallcell_options **opts_ptr)
{
	struct smallcell_options *opts;
	opts = cfmalloc(sizeof(struct smallcell_options));
	if ( opts == NULL ) return ENOMEM;
	opts->dummy = 0;
	*opts_ptr = opts;
	return 0;
}


static error_t smallcell_parse_arg(int key, char *arg, struct argp_state *state)
{
	struct smallcell_options **opts_ptr = state->input;
	int r;

	switch (key) {

		case ARGP_KEY_INIT:
		r = smallcell_default_options(opts_ptr);
		if (r) return r;
		break;

		case 1:
		smallcell_show_help();
		return EINVAL;

		default:
		return ARGP_ERR_UNKNOWN;

	}

	return 0;
}


static struct argp_option smallcell_options[] = {

	{"help-smallcell", 1, NULL, OPTION_NO_USAGE,
	 "Show options for 'smallcell' indexing", 99},

	{0}
};


struct argp smallcell_argp = { smallcell_options, smallcell_parse_arg,
	NULL, NULL, NULL, NULL, NULL
};
