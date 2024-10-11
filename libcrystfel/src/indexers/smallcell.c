/*
 * smallcell.c
 *
 * Perform indexing from solution file
 *
 * Copyright © 2020-2021 Max-Planck-Gesellschaft
 *                       zur Förderung der Wissenschaften e.V.
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>
 *   2021 Thomas White <thomas.white@desy.de>
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
#define DIFF_TOL (1e8)
#define PIXEL_RING_TOL (6)


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


void *smallcell_prepare(IndexingMethod * indm, struct smallcell_options *opts,
                        UnitCell * cell)
{
	STATUS("\n");
	STATUS("*******************************************************************\n");
	STATUS("****                    Welcome to SmallCell                  ****\n");
	STATUS("*******************************************************************\n");
	STATUS("\n");


	struct smallcell_private *dp;
	double asx, bsx, csx;
	double asy, bsy, csy;
	double asz, bsz, csz;

	dp = cfmalloc(sizeof(struct smallcell_private));

	dp->sym = get_lattice_symmetry(cell);
	dp->powderrings = powder_rings(cell, dp->sym, 1/1e-10, &dp->num_rings);

	/* Get reciprocal unit cell elements in order to create G* matrix
	 * and store in private */

	cell_get_reciprocal(cell,
	                    &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz, &csx, &csy, &csz);

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


/* PeakInfo structure for storing peak-related information for every match found*/
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
	struct Nodelist *list[MAX_NEIGH];
};

static struct Nodelist *CopyRlist(struct Nodelist *R)
{
	struct Nodelist *Rcopy = cfmalloc(sizeof(struct Nodelist));
	int i;
	Rcopy->n_mem = 0;
	for (i = 0; i < R->n_mem; i++) {
		Rcopy->mem[i] = R->mem[i];
		Rcopy->n_mem++;
	}
	return Rcopy;
};

//To make a list of neighbour nodes
static struct Nodelist *neighbours(struct PeakInfo *in)
{
	int i;
	struct Nodelist *list = cfmalloc(sizeof(struct Nodelist));
	for (i = 0; i < in->n_neigh; i++) {
		list->mem[i] = in->neigh[i];
	}
	list->n_mem = in->n_neigh;
	return list;
};

//Function to see if a node is in a list
static int isin(struct Nodelist *b, struct PeakInfo *test)
{
	int i;
	for (i = 0; i < b->n_mem; i++) {
		if (test == b->mem[i]) {
			return 1;       //isin
		}
	}
	return 0;               //is notin
};

static void add(struct Nodelist *c, struct PeakInfo *test)
{

	if (isin(c, test))
		return;
	assert(c->n_mem < MAX_NODES);
	c->mem[c->n_mem] = test;
	c->n_mem++;
};

//To create a new list of nodes from the union between two node lists (i.e. c = a U b)
static struct Nodelist *Union(struct Nodelist *a, struct Nodelist *b)
{
	struct Nodelist *c = cfmalloc(sizeof(struct Nodelist));
	c->n_mem = 0;
	int i;
	int j;
	for (i = 0; i < a->n_mem; i++) {
		add(c, a->mem[i]);
	}
	for (j = 0; j < b->n_mem; j++) {
		add(c, b->mem[j]);
	}
	return c;
};

//To create a new list of nodes from the intersection of two lists (i.e. c = a intersection b)
static struct Nodelist *intersection(struct Nodelist *a, struct Nodelist *b)
{
	struct Nodelist *c = cfmalloc(sizeof(struct Nodelist));
	c->n_mem = 0;
	int j;
	for (j = 0; j < a->n_mem; j++) {
		if (isin(b, a->mem[j]) == 1) {
			add(c, a->mem[j]);
		}
	}
	return c;
};

//To create a new list of nodes (c) from the exclusion of list b from list a (i.e. c = a\b)
static struct Nodelist *exclusion(struct Nodelist *a, struct Nodelist *b)
{
	struct Nodelist *c = cfmalloc(sizeof(struct Nodelist));
	c->n_mem = 0;
	int j;
	for (j = 0; j < a->n_mem; j++) {
		if (isin(b, a->mem[j]) == 0) {
			add(c, a->mem[j]);
		}
	}
	return c;
};

//Function to exclude a single node
static struct Nodelist *exclunode(struct Nodelist *a, struct PeakInfo *v)
{
	int j;
	for (j = 0; j < a->n_mem; j++) {
		if (a->mem[j] == v) {
			int i;
			int hole_index = j;
			for (i = hole_index; i < a->n_mem - 1; i++) {
				a->mem[i] = a->mem[i + 1];
			}
		}
	}
	a->n_mem = a->n_mem - 1;
	return a;
};

//Function to append a node to a list (creating a new list)
static struct Nodelist *append(struct Nodelist *a, struct PeakInfo *v)
{
	struct Nodelist *c = cfmalloc(sizeof(struct Nodelist));
	c->n_mem = 0;
	int j;
	for (j = 0; j < a->n_mem; j++) {
		add(c, a->mem[j]);
	}

	add(c, v);
	return c;
};

//Average weight function
static double avg_weight(int num_neigh, double *weights)
{

	double avg;
	double avg_add = 0.0;
	int j;
	for (j = 0; j < num_neigh; j++) {
		avg_add = avg_add + weights[j];
	}

	avg = avg_add / num_neigh;

	return avg;
};

//Bron-Kerbosch algorithm, with pivoting
void BK(struct Nodelist *R, struct Nodelist *P, struct Nodelist *X,
        struct Cliquelist *Max_cliques)
{

	//If P & X are empty -> add R to the un-mapped clique array
	if (P->n_mem == 0 && X->n_mem == 0) {
		Max_cliques->list[Max_cliques->n] = CopyRlist(R);
		Max_cliques->n++;


		return;
	}
	//Find pivot u from set of nodes in P U X
	struct Nodelist *piv_pool = Union(P, X);
	struct PeakInfo *piv = NULL;
	piv = piv_pool->mem[0];
	int n_max = piv_pool->mem[0]->n_neigh;
	double n_max_weight =
	        avg_weight(piv_pool->mem[0]->n_neigh,
	                   piv_pool->mem[0]->weight_list);
	int i;
	for (i = 1; i < piv_pool->n_mem; i++) {
		if (piv_pool->mem[i]->n_neigh > n_max) {
			n_max = piv_pool->mem[i]->n_neigh;
			piv = piv_pool->mem[i];
			n_max_weight =
			        avg_weight(piv_pool->mem[i]->n_neigh,
			                   piv_pool->mem[i]->weight_list);
		} else if (piv_pool->mem[i]->n_neigh == n_max) {
			if (piv_pool->n_mem == 1) {
				piv = piv_pool->mem[0];
			} else if ((n_max > 0)
			           &&
			           (avg_weight
			            (piv_pool->mem[i]->n_neigh,
			             piv_pool->mem[i]->weight_list) <
			            n_max_weight)) {
				piv = piv_pool->mem[i];
			}
		}
	}

	if (piv == NULL) {
		printf("Couldn't find pivot\n");
		return;
	}
	//Get list of neighbours of the pivot
	struct Nodelist *piv_neighbours = neighbours(piv);
	//Remove pivot neighbours from P list
	struct Nodelist *P_excl = exclusion(P, piv_neighbours);
	cffree(piv_pool);
	cffree(piv_neighbours);

	int u;
	for (u = 0; u < P_excl->n_mem; u++) {
		//For all members of the subset of P we are concered with (P_excl = P\N(pivot))
		struct PeakInfo *v = P_excl->mem[u];

		//create Nodelist for v using neighbour function
		struct Nodelist *v_neighs = neighbours(v);


		//add v to R
		struct Nodelist *R_new = append(R, v);


		//Find intersection of P_excl with N(v)
		struct Nodelist *P_new = intersection(P, v_neighs);


		//''                 '' X with N(v)    
		struct Nodelist *X_new = intersection(X, v_neighs);
		BK(R_new, P_new, X_new, Max_cliques);

		cffree(R_new);
		cffree(P_new);
		cffree(X_new);
		cffree(v_neighs);


		//Redefine P and X as P\v and XUv

		exclunode(P, v);

		add(X, v);
	}
	cffree(P_excl);
	return;
};


int smallcell_index(struct image *image, void *mpriv)
{
	struct smallcell_private *priv = (struct smallcell_private *) mpriv;
	struct powder_ring *powderrings = priv->powderrings;
	int num_rings = priv->num_rings;


	//Assigning image data
	ImageFeatureList *peaks = image->features;
	double lambda = image->lambda;
	struct detgeom *det = image->detgeom;

	//Calculating a tolerance value of PIXEL_RING_TOL  pixel lengths from a peak
	double pixel_metres = det->panels[0].pixel_pitch;
	double ten_pix_len = pixel_metres * PIXEL_RING_TOL;
	double detector_len = det->panels[0].cnz;
	double detector_len_m = pixel_metres * detector_len;
	double tol = ten_pix_len / (lambda * detector_len_m);

	int npk = image_feature_count(peaks);

	// No peaks -> no resolution! 
	if (npk <= 3) {
		printf("not enough peaks\n");
		return 0;
	}

	printf("Current Image %s, %s\n", image->filename, image->ev);

	int peak_infos_size = 100;      // Arbitrary initial size allocation
	PeakInfo *peak_infos = cfmalloc(peak_infos_size * sizeof(PeakInfo));

	int num_peak_infos = 0;
	int peaks_with_matches = 0;
	int i;

	//Loop through each peak, calculate d, then 1/d value (based on estimate_peak_resolution from peak.c), then use match_rings function to create structs for this peak with all possible h,k,l values
	for (i = 0; i < npk; i++) {
		struct imagefeature *f = image_get_feature(peaks, i);
		double r[3];
		detgeom_transform_coords(&det->panels[f->pn], f->fs, f->ss,
		                         lambda, 0.0, 0.0, r);
		double x_res = r[0];
		double y_res = (r[1]);
		double z_res = (r[2]);
		double rns = modulus(r[0], r[1], r[2]);
		int init_num_peak_infos = num_peak_infos;
		int j;

		for (j = 0; j < num_rings; j++) {
			if (fabs(powderrings[j].resolution - rns) <= tol) {

				signed int h = powderrings[j].h;
				signed int k = powderrings[j].k;
				signed int l = powderrings[j].l;

				//Looking for symmetries and creating more PeakInfo structs with these symmetry indices
				SymOpMask *m = new_symopmask(priv->sym);
				special_position(priv->sym, m, h, k, l);
				int n = num_equivs(priv->sym, m);
				int y;
				for (y = 0; y < n; y++) {

					if (num_peak_infos >= peak_infos_size) {        //Adding more structs if necessary
						peak_infos_size *= 2;
						peak_infos =
						        cfrealloc(peak_infos,
						                  peak_infos_size
						                  *
						                  sizeof
						                  (PeakInfo));
					}

					signed int ha, ka, la;
					get_equiv(priv->sym, m, y, h, k, l, &ha,
					          &ka, &la);
					// Add ha, ka, la to list of reflections 
					peak_infos[num_peak_infos].peak_number =
					        i;
					peak_infos[num_peak_infos].peak_res =
					        rns;
					peak_infos[num_peak_infos].h = ha;
					peak_infos[num_peak_infos].k = ka;
					peak_infos[num_peak_infos].l = la;
					peak_infos[num_peak_infos].x = x_res;
					peak_infos[num_peak_infos].y = y_res;
					peak_infos[num_peak_infos].z = z_res;
					int w;
					for (w = 0; w < MAX_NEIGH; w++) {
						peak_infos[num_peak_infos].
						        neigh[w] = NULL;
					}
					peak_infos[num_peak_infos].n_neigh = 0; //Initialise
					(num_peak_infos)++;
				}
				free_symopmask(m);
			}

		}

		if (num_peak_infos != init_num_peak_infos) {
			peaks_with_matches++;
		}
	}

	printf("The number of matches in this image (including symmetric indices) is %d for %d/%d peaks\n", num_peak_infos, peaks_with_matches, npk);


	// Now to connect the nodes using calculated and measured reciprocal distance


	struct g_matrix g9 = priv->g9;
	double dtol = DIFF_TOL;
	int n_connected_nodes = 0;
	int j;

	//Loop through peak numbers
	for (j = 0; j < num_peak_infos; j++) {

		int node_a_h = peak_infos[j].h;
		int node_a_k = peak_infos[j].k;
		int node_a_l = peak_infos[j].l;
		int y;

		//Loop through the rest of the peak infos
		for (y = j + 1; y < num_peak_infos - 1; y++) {

			if (peak_infos[y].peak_number ==
			    peak_infos[j].peak_number)
				continue;

			//Observed dist. d_1 = fabs(res.vec_j - res.vec_y)
			double d_1 =
			        modulus(peak_infos[j].x - peak_infos[y].x,
			                peak_infos[j].y - peak_infos[y].y,
			                peak_infos[j].z - peak_infos[y].z);

			//Predicted d_2
			int node_b_h = peak_infos[y].h;
			int node_b_k = peak_infos[y].k;
			int node_b_l = peak_infos[y].l;

			//d_2 = sqrt(Transpose(MI_b - MI_a).G*.(MI_b - MI_a))
			double d_2 =
			        sqrt((node_b_h -
			              node_a_h) * (g9.A * (node_b_h -
			                                   node_a_h) +
			                           g9.D * (node_b_k -
			                                   node_a_k) +
			                           g9.G * (node_b_l - node_a_l))
			             + (node_b_k -
			                node_a_k) * (g9.B * (node_b_h -
			                                     node_a_h) +
			                             g9.E * (node_b_k -
			                                     node_a_k) +
			                             g9.H * (node_b_l -
			                                     node_a_l))
			             + (node_b_l -
			                node_a_l) * (g9.C * (node_b_h -
			                                     node_a_h) +
			                             g9.F * (node_b_k -
			                                     node_a_k) +
			                             g9.J * (node_b_l -
			                                     node_a_l)));


			double diff = fabs(d_2 - d_1);

			//Now test difference
			if (diff <= dtol) {

				//Connect nodes

				if (peak_infos[j].n_neigh <= MAX_NEIGH) {

					peak_infos[j].neigh[peak_infos[j].n_neigh] = &peak_infos[y];    //Pointing to the info of the connected peak (adding node_b as neighbour due to its connection to node_a)
					peak_infos[j].weight_list[peak_infos[j].n_neigh] = diff;        //add weight for later
					peak_infos[j].n_neigh++;        //increasing number of connection/neighbours that node_a has

				} else {

					printf("The number of neighbours for this node has exceeded MAX_NEIGH, need to manage memory\n");
				}

				if (peak_infos[y].n_neigh <= MAX_NEIGH) {

					peak_infos[y].neigh[peak_infos[y].n_neigh] = &peak_infos[j];    //Also pointing other way(i.e. adding node_a as a neghbour to node_b as we wont loop back through the upper part of the list when j = peak_no of node_b)
					peak_infos[y].weight_list[peak_infos[y].n_neigh] = diff;        //add weight for later
					peak_infos[y].n_neigh++;        //increasing number of connection/neighbours that node_b has
				} else {
					printf("The number of neighbours for this node has exceeded MAX_NEIGH, need to manage memory\n");
				}

			}

		}

		//Counting number of nodes with 1 or more connection (may be unnecessary)...
		if (peak_infos[j].n_neigh != 0) {
			n_connected_nodes++;
		}

	}

	//Store the max.cliques found
	struct Cliquelist *Max_cliques =
	        cfmalloc(sizeof(struct Nodelist) * sizeof(struct Cliquelist));
	Max_cliques->n = 0;
	// R: array of nodes forming a clique (array of pointers to node infos)
	// P: array of all prosepctive nodes that are connected to R which may be added to R. To begin, this is all nodes i.e all peak_infos
	// X: exculsion set (same form as R but nodes that are NOT candidats for the max. clique, were originaly in P)
	struct Nodelist *R = cfmalloc(sizeof(struct Nodelist));
	struct Nodelist *X = cfmalloc(sizeof(struct Nodelist));
	struct Nodelist *P = cfmalloc(sizeof(struct Nodelist));
	R->n_mem = 0;
	X->n_mem = 0;
	P->n_mem = 0;
	//To make P; create nodelist of all peak_infos
	for (i = 0; i < num_peak_infos; i++) {
		if (peak_infos[i].n_neigh != 0) {
			add(P, &peak_infos[i]);
		}
	}

	if (P->n_mem <= 2) {
		printf("no peaks with neighbours\n");
		return 0;
	}
	//Call BK using current peak info nodes for this image

	printf("running BK\n");
	BK(R, P, X, Max_cliques);
	printf("done\n");
	printf("The number of cliques found = %d.\n", Max_cliques->n);

	//get the max. clique from list of Maximal cliques found
	int Max_clique_len = Max_cliques->list[0]->n_mem;
	struct Cliquelist *Max =
	        cfmalloc(sizeof(struct Nodelist) * sizeof(struct Cliquelist));
	Max->n = 0;
	Max->list[0] = Max_cliques->list[0];
	Max->n++;
	int m;
	for (m = 1; m < Max_cliques->n; m++) {
		if (Max_cliques->list[m]->n_mem > Max_clique_len) {
			Max_clique_len = Max_cliques->list[m]->n_mem;
			int t;
			for (t = 0; t < Max->n; t++) {
				Max->list[t] = NULL;
			}
			Max->n = 0;
			Max->list[0] = Max_cliques->list[m];
			Max->n++;
		} else if (Max_cliques->list[m]->n_mem == Max_clique_len) {
			Max->list[Max->n] = Max_cliques->list[m];
			Max->n++;
		}
	}
	//If more than one max_clique with the same number of nodes is found, take only the right-handed solution
	//This requires first getting the unit cell for each max_clique, and then using the right_handed function from cell-utils
	for (m = 0; m < Max->n; m++) {
		if (Max->list[m]->n_mem < 5 && m == (Max->n) - 1)
			return 0;
		if (Max->list[m]->n_mem < 5)
			continue;
		gsl_matrix *h_mat =
		        gsl_matrix_calloc(3 * (Max->list[m]->n_mem), 9);
		gsl_vector *h_vec = gsl_vector_alloc(3 * (Max->list[m]->n_mem));

		int count_node = 0;
		int col_count = 0;
		int have_a = 0, have_b = 0, have_c = 0;
		int i;
		for (i = 0; i < 3 * (Max->list[m]->n_mem); i++) {
			if (i > 0 && i % 3 == 0) {
				count_node++;
				col_count = 0;
			}

			if (Max->list[m]->mem[count_node]->h != 0)
				have_a = 1;
			if (Max->list[m]->mem[count_node]->k != 0)
				have_b = 1;
			if (Max->list[m]->mem[count_node]->l != 0)
				have_c = 1;

			gsl_matrix_set(h_mat, i, col_count,
			               Max->list[m]->mem[count_node]->h);
			gsl_matrix_set(h_mat, i, col_count + 3,
			               Max->list[m]->mem[count_node]->k);
			gsl_matrix_set(h_mat, i, col_count + 6,
			               Max->list[m]->mem[count_node]->l);
			col_count++;
		}

		int count_mem = 0;
		int count_comp = 0;
		for (i = 0; i < 3 * (Max->list[m]->n_mem); i++) {

			gsl_vector_set(h_vec, i,
			               Max->list[m]->mem[count_mem]->x);
			if (count_comp == 0) {
				gsl_vector_set(h_vec, i,
				               Max->list[m]->mem[count_mem]->x);
				count_comp++;
			} else if (count_comp == 1) {
				gsl_vector_set(h_vec, i,
				               Max->list[m]->mem[count_mem]->y);
				count_comp++;
			} else if (count_comp == 2) {
				gsl_vector_set(h_vec, i,
				               Max->list[m]->mem[count_mem]->z);
				count_comp = 0;
				count_mem++;
			}
		}
		//Solve matrix-vector equation for unit-cell for this clique m
		gsl_vector *cell_vecs = gsl_vector_alloc(9);
		//cell_vec = [a*x a*y a*z b*x b*y b*z c*x c*y c*z] 
		double chisq;
		gsl_matrix *cov =
		        gsl_matrix_alloc(3 * (Max->list[m]->n_mem), 9);
		gsl_multifit_linear_workspace *work =
		        gsl_multifit_linear_alloc(3 * (Max->list[m]->n_mem), 9);
		if (gsl_multifit_linear
		    (h_mat, h_vec, cell_vecs, cov, &chisq, work)) {
			ERROR("Multifit failed\n");
		}
		//Use the following function to make a unit cell file then can use the checks directly to see if it's a viable solution

		UnitCell *uc;

		uc = cell_new();
/*	UnitCell *cell = priv->template;
	cell_set_lattice_type(uc, cell_get_lattice_type(cell));*/
		cell_set_reciprocal(uc,
		                    gsl_vector_get(cell_vecs, 0),
		                    gsl_vector_get(cell_vecs, 1),
		                    gsl_vector_get(cell_vecs, 2),
		                    gsl_vector_get(cell_vecs, 3),
		                    gsl_vector_get(cell_vecs, 4),
		                    gsl_vector_get(cell_vecs, 5),
		                    gsl_vector_get(cell_vecs, 6),
		                    gsl_vector_get(cell_vecs, 7),
		                    gsl_vector_get(cell_vecs, 8));


		if (uc == NULL) {
			printf("Unit Cell not created.. returned NULL\n");
			continue;
		}

		//Free up matrix and vector memeories
		gsl_multifit_linear_free(work);
		gsl_vector_free(cell_vecs);
		gsl_vector_free(h_vec);
		gsl_matrix_free(cov);
		gsl_matrix_free(h_mat);

		if (!(have_a && have_b && have_c))
			return 0;

		printf("Unit Cell created, testing if valid solution..\n");
		if (validate_cell(uc) == 0) {

			Crystal *cr;
			printf("unit cell valid!\n");
			cell_print(uc);
			cr = crystal_new();
			if (cr == NULL) {
				ERROR("Failed to allocate crystal.\n");
				continue;
			}
			crystal_set_cell(cr, uc);
			image_add_crystal(image, cr);
			return 1;
		}

		cell_free(uc);

	}


	for (i = 0; i < Max_cliques->n; i++) {
		cffree(Max_cliques->list[i]);
	}
	for (i = 0; i < Max->n; i++) {
		cffree(Max->list[i]);
	}

	cffree(Max_cliques);
	cffree(Max);
	cffree(X);
	cffree(P);
	cffree(R);
	cffree(peak_infos);
	return 0;
}


void smallcell_cleanup(void *mpriv)
{
/*	struct smallcell_private *dp = mpriv;
	struct smallcell_entry *item, *tmp;

	HASH_ITER(hh, dp->sol_hash, item, tmp) {
		int i;
		HASH_DEL(dp->sol_hash, item);
		for ( i=0; i<item->n_crystals; i++ ) {
			Crystal *cr = item->crystals[i];
			cell_free(crystal_get_cell(cr));
			crystal_free(cr);
		}
	}

	cffree(dp);*/
}


static void smallcell_show_help()
{
	printf("Parameters for 'smallcell' indexing:\n"
	       "     --smallcell-input-file\n"
	       "                           Filename of indexing solution file\n");
}


int smallcell_default_options(struct smallcell_options **opts_ptr)
{
	struct smallcell_options *opts;
	opts = cfmalloc(sizeof(struct smallcell_options));
	if (opts == NULL)
		return ENOMEM;
	opts->filename = NULL;
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
			if (r)
				return r;
			break;

		case 1:
			smallcell_show_help();
			return EINVAL;

		case 2:
			(*opts_ptr)->filename = cfstrdup(arg);
			break;

		default:
			return ARGP_ERR_UNKNOWN;

	}

	return 0;
}


static struct argp_option smallcell_options[] = {

	{"help-smallcell", 1, NULL, OPTION_NO_USAGE,
	 "Show options for 'from file' indexing", 99},

	{"smallcell-input-file", 2, "filename", OPTION_HIDDEN, NULL},
	{0}
};


struct argp smallcell_argp = { smallcell_options, smallcell_parse_arg,
	NULL, NULL, NULL, NULL, NULL
};
