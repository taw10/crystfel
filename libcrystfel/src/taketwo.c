/*
 * taketwo.c
 *
 * Rewrite of TakeTwo algorithm (Acta D72 (8) 956-965) for CrystFEL
 *
 * Copyright © 2016-2017 Helen Ginn
 * Copyright © 2016-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2016-2017 Helen Ginn <helen@strubi.ox.ac.uk>
 *   2016-2017 Thomas White <taw@physics.org>
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

/**
 * \class TakeTwo
 * Code outline.
 * --- Get ready for calculation ---
 * Pre-calculate symmetry operations (generate_rotation_symops())
 * Pre-calculate theoretical vectors from unit cell dimensions
 * (gen_theoretical_vecs())
 * Generate observed vectors from data (gen_observed_vecs())
 * Match observed vectors to theoretical vectors (match_obs_to_cell_vecs())
 *
 * --- Business bit ---
 * ... n.b. rearranging to find all seeds in advance.
 *
 * Find starting seeds (find_seeds()):
 *   - Loop through pairs of observed vectors
 *   - If they share a spot, find matching pairs of theoretical vectors
 *   - Remove all duplicate matches due to symmetry operations
 *   - For the remainder, loop through the matches and extend the seeds
 *   (start_seed()).
 *   - If it returns a membership greater than the highest member threshold,
 *   return the matrix to CrystFEL.
 *
 * Extending a seed (start_seed()):
 *   - Generate a rotation matrix which matches the chosen start seed.
 *   - Loop through all observed vectors starting from 0.
 *   - Find another vector to add to the network, if available
 *   (find_next_index()).
 *   - If the index is not available, then give up if there were too many dead
 *   ends. Otherwise, remove the last member and pretend like that didn't
 *   happen, so it won't happen again.
 *   - Add the vector to increment the membership list.
 *   - If the membership number exceeds the maximum, tidy up the solution and
 *   return a success.
 *   - If the membership does not, then resume the loop and search for the
 *   next vector.
 *
 * Finding the next member (find_next_index()):
 *   - Go through the observed vectors, starting from the last index + 1 to
 *   explore only the "new" vectors.
 *   - If the vector does not share a spot with the current array of vectors,
 *   then skip it.
 *   - We must loop through all the current vectors in the network, and see if
 *   they match the newcomer for a given matching theoretical vector.
 *   - We only accept a match if the rotation matrix matches the seed matrix
 *   for a single matching theoretical vector.
 *   - If it does match, we can return a success.
 *
 * Tidying the solution (finish_solution()):
 *   - This chooses the most common rotation matrix of the bunch to choose to
 *   send to CrystFEL. But this should probably take the average instead,
 *   which is very possible.
 *
 * Clean up the mess (cleanup_taketwo_obs_vecs())
 */

/**
 * Helen's to-do list
 * -
 *
 *
 * - Improve the final solution
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "cell-utils.h"
#include "index.h"
#include "taketwo.h"
#include "peaks.h"
#include "symmetry.h"

/**
 * spotvec
 * @obsvec: an observed vector between two spots
 * @matches: array of matching theoretical vectors from unit cell
 * @match_num: number of matches
 * @distance: length of obsvec (do I need this?)
 * @her_rlp: pointer to first rlp position for difference vec
 * @his_rlp: pointer to second rlp position for difference vec
 *
 * Structure representing 3D vector between two potential Bragg peaks
 * in reciprocal space, and an array of potential matching theoretical
 * vectors from unit cell/centering considerations.
 **/
struct SpotVec
{
	struct rvec obsvec;
	struct TheoryVec *matches;
	int match_num;
	int assignment;
	int in_network;
	double distance;
	struct rvec *her_rlp;
	struct rvec *his_rlp;
};

/**
 * theoryvec
 */

struct TheoryVec
{
	struct rvec vec;
	int asym;
};


/**
 * seed
 */

struct Seed
{
	int obs1;
	int obs2;
	int idx1;
	int idx2;
	double score;
};

struct taketwo_private
{
	IndexingMethod indm;
	UnitCell       *cell;
	int   serial_num; /* -1 if unassigned */
	unsigned int attempts;

	gsl_matrix     **prevSols; /**< Previous solutions to be ignored */
	unsigned int   numPrevs; /**< Previous solution count */

};

/**
 * Internal structure which gets passed the various functions looking after
 * the indexing bits and bobs. */
struct TakeTwoCell
{
	UnitCell       *cell; /**< Contains unit cell dimensions */
	gsl_matrix     **rotSymOps;
	unsigned int   numOps;

	struct SpotVec *obs_vecs;
	struct Seed *seeds;
	int seed_count;
	int obs_vec_count;

	/* Options */
	int member_thresh;
	double len_tol;   /**< In reciprocal metres */
	double angle_tol; /**< In radians */
	double trace_tol; /**< Contains sqrt(4*(1-cos(angle))) */

	/** A given solution to refine */
	gsl_matrix *solution;

	double x_ang; /**< Rotations in radians to apply to x axis of solution */
	double y_ang; /**< Rotations in radians to apply to y axis of solution */
	double z_ang; /**< Rotations in radians to apply to z axis of solution */

	/**< Temporary memory always allocated for calculations */
	gsl_matrix *twiz1Tmp;
	/**< Temporary memory always allocated for calculations */
	gsl_matrix *twiz2Tmp;
	/**< Temporary memory always allocated for calculations */
	gsl_vector *vec1Tmp;
	/**< Temporary memory always allocated for calculations */
	gsl_vector *vec2Tmp;
};


/* Maximum distance between two rlp sizes to consider info for indexing */
#define MAX_RECIP_DISTANCE (0.15*1e10)

/* Tolerance for two lengths in reciprocal space to be considered the same */
#define RECIP_TOLERANCE (0.0010*1e10)

/* Threshold for network members to consider a potential solution */
#define NETWORK_MEMBER_THRESHOLD (20)

/* Minimum for network members to consider a potential solution */
#define MINIMUM_MEMBER_THRESHOLD (5)

/* Maximum dead ends for a single branch extension during indexing */
#define MAX_DEAD_ENDS (10)

/* Maximum observed vectors before TakeTwo gives up and deals with
 * what is already there. */
#define MAX_OBS_VECTORS 100000

/* Tolerance for two angles to be considered the same */
#define ANGLE_TOLERANCE (deg2rad(0.6))

/* Tolerance for rot_mats_are_similar */
#define TRACE_TOLERANCE (deg2rad(3.0))

/* Initial step size for refinement of solutions */
#define ANGLE_STEP_SIZE (deg2rad(0.5))

/* Final required step size for refinement of solutions */
#define ANGLE_CONVERGE_SIZE (deg2rad(0.01))

/* TODO: Multiple lattices */


/* ------------------------------------------------------------------------
 * apologetic function
 * ------------------------------------------------------------------------*/

void apologise()
{
	printf("Error - could not allocate memory for indexing.\n");
}

/* ------------------------------------------------------------------------
 * functions concerning aspects of rvec which are very likely to be
 * incorporated somewhere else in CrystFEL and therefore may need to be
 * deleted and references connected to a pre-existing function. (Lowest level)
 * ------------------------------------------------------------------------*/

static struct rvec new_rvec(double new_u, double new_v, double new_w)
{
	struct rvec new_rvector;
	new_rvector.u = new_u;
	new_rvector.v = new_v;
	new_rvector.w = new_w;

	return new_rvector;
}

static struct rvec rvec_add_rvec(struct rvec first, struct rvec second)
{
	struct rvec diff = new_rvec(second.u + first.u,
				    second.v + first.v,
				    second.w + first.w);

	return diff;
}


static struct rvec diff_vec(struct rvec from, struct rvec to)
{
	struct rvec diff = new_rvec(to.u - from.u,
				    to.v - from.v,
				    to.w - from.w);

	return diff;
}

static double sq_length(struct rvec vec)
{
	double sqlength = (vec.u * vec.u + vec.v * vec.v + vec.w * vec.w);

	return sqlength;
}


static double rvec_length(struct rvec vec)
{
	return sqrt(sq_length(vec));
}


static void normalise_rvec(struct rvec *vec)
{
	double length = rvec_length(*vec);
	vec->u /= length;
	vec->v /= length;
	vec->w /= length;
}


static double rvec_cosine(struct rvec v1, struct rvec v2)
{
	double dot_prod = v1.u * v2.u + v1.v * v2.v + v1.w * v2.w;
	double v1_length = rvec_length(v1);
	double v2_length = rvec_length(v2);

	double cos_theta = dot_prod / (v1_length * v2_length);

	return cos_theta;
}


static double rvec_angle(struct rvec v1, struct rvec v2)
{
	double cos_theta = rvec_cosine(v1, v2);
	double angle = acos(cos_theta);

	return angle;
}


static struct rvec rvec_cross(struct rvec a, struct rvec b)
{
	struct rvec c;

	c.u = a.v*b.w - a.w*b.v;
	c.v = -(a.u*b.w - a.w*b.u);
	c.w = a.u*b.v - a.v*b.u;

	return c;
}

/*
static void show_rvec(struct rvec r2)
{
	struct rvec r = r2;
	normalise_rvec(&r);
	STATUS("[ %.3f %.3f %.3f ]\n", r.u, r.v, r.w);
}
*/


/* ------------------------------------------------------------------------
 * functions called under the core functions, still specialised (Level 3)
 * ------------------------------------------------------------------------*/

static void rotation_around_axis(struct rvec c, double th,
				 gsl_matrix *res)
{
	double omc = 1.0 - cos(th);
	double s = sin(th);
	gsl_matrix_set(res, 0, 0, cos(th) + c.u*c.u*omc);
	gsl_matrix_set(res, 0, 1, c.u*c.v*omc - c.w*s);
	gsl_matrix_set(res, 0, 2, c.u*c.w*omc + c.v*s);
	gsl_matrix_set(res, 1, 0, c.u*c.v*omc + c.w*s);
	gsl_matrix_set(res, 1, 1, cos(th) + c.v*c.v*omc);
	gsl_matrix_set(res, 1, 2, c.v*c.w*omc - c.u*s);
	gsl_matrix_set(res, 2, 0, c.w*c.u*omc - c.v*s);
	gsl_matrix_set(res, 2, 1, c.w*c.v*omc + c.u*s);
	gsl_matrix_set(res, 2, 2, cos(th) + c.w*c.w*omc);
}

/** Rotate GSL matrix by three angles along x, y and z axes */
static void rotate_gsl_by_angles(gsl_matrix *sol, double x, double y,
                                 double z, gsl_matrix *result)
{
	gsl_matrix *x_rot = gsl_matrix_alloc(3, 3);
	gsl_matrix *y_rot = gsl_matrix_alloc(3, 3);
	gsl_matrix *z_rot = gsl_matrix_alloc(3, 3);
	gsl_matrix *xy_rot = gsl_matrix_alloc(3, 3);
	gsl_matrix *xyz_rot = gsl_matrix_alloc(3, 3);

	struct rvec x_axis = new_rvec(1, 0, 0);
	struct rvec y_axis = new_rvec(1, 0, 0);
	struct rvec z_axis = new_rvec(1, 0, 0);

	rotation_around_axis(x_axis, x, x_rot);
	rotation_around_axis(y_axis, y, y_rot);
	rotation_around_axis(z_axis, z, z_rot);

	/* Collapse the rotations in x and y onto z */
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, x_rot,
	               y_rot, 0.0, xy_rot);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, xy_rot,
	               z_rot, 0.0, xyz_rot);

	/* Apply the whole rotation offset to the solution */
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, xyz_rot,
	               sol, 0.0, result);
}


/* Rotate vector (vec1) around axis (axis) by angle theta. Find value of
 * theta for which the angle between (vec1) and (vec2) is minimised. */
static void closest_rot_mat(struct rvec vec1, struct rvec vec2,
			    struct rvec axis, gsl_matrix *twizzle)
{
	/* Let's have unit vectors */
	normalise_rvec(&vec1);
	normalise_rvec(&vec2);
	normalise_rvec(&axis);

	/* Redeclaring these to try and maintain readability and
	 * check-ability against the maths I wrote down */
	double a = vec2.u; double b = vec2.v; double c = vec2.w;
	double p = vec1.u; double q = vec1.v; double r = vec1.w;
	double x = axis.u; double y = axis.v; double z = axis.w;

	/* Components in handwritten maths online when I upload it */
	double A = a*(p*x*x - p + x*y*q + x*z*r) +
	b*(p*x*y + q*y*y - q + r*y*z) +
	c*(p*x*z + q*y*z + r*z*z - r);

	double B = a*(y*r - z*q) + b*(p*z - r*x) + c*(q*x - p*y);

	double tan_theta = - B / A;
	double theta = atan(tan_theta);

	/* Now we have two possible solutions, theta or theta+pi
	 * and we need to work out which one. This could potentially be
	 * simplified - do we really need so many cos/sins? maybe check
	 * the 2nd derivative instead? */
	double cc = cos(theta);
	double C = 1 - cc;
	double s = sin(theta);
	double occ = -cc;
	double oC = 1 - occ;
	double os = -s;

	double pPrime = (x*x*C+cc)*p + (x*y*C-z*s)*q + (x*z*C+y*s)*r;
	double qPrime = (y*x*C+z*s)*p + (y*y*C+cc)*q + (y*z*C-x*s)*r;
	double rPrime = (z*x*C-y*s)*p + (z*y*C+x*s)*q + (z*z*C+cc)*r;

	double pDbPrime = (x*x*oC+occ)*p + (x*y*oC-z*os)*q + (x*z*oC+y*os)*r;
	double qDbPrime = (y*x*oC+z*os)*p + (y*y*oC+occ)*q + (y*z*oC-x*os)*r;
	double rDbPrime = (z*x*oC-y*os)*p + (z*y*oC+x*os)*q + (z*z*oC+occ)*r;

	double cosAlpha = pPrime * a + qPrime * b + rPrime * c;
	double cosAlphaOther = pDbPrime * a + qDbPrime * b + rDbPrime * c;

	int addPi = (cosAlphaOther > cosAlpha);
	double bestAngle = theta + addPi * M_PI;

	/* Don't return an identity matrix which has been rotated by
	 * theta around "axis", but do assign it to twizzle. */
	rotation_around_axis(axis, bestAngle, twizzle);
}

/*
static double matrix_angle(gsl_matrix *m)
{
	double a = gsl_matrix_get(m, 0, 0);
	double b = gsl_matrix_get(m, 1, 1);
	double c = gsl_matrix_get(m, 2, 2);

	double cos_t = (a + b + c - 1) / 2;
	double theta = acos(cos_t);

	return theta;
}

static struct rvec matrix_axis(gsl_matrix *a)
{
	double ang = matrix_angle(a);
	double cos_t = cos(ang);
	double p = gsl_matrix_get(a, 0, 0);
	double q = gsl_matrix_get(a, 1, 1);
	double r = gsl_matrix_get(a, 2, 2);
	double x = sqrt((p - cos_t) / (1 - cos_t));
	double y = sqrt((q - cos_t) / (1 - cos_t));
	double z = sqrt((r - cos_t) / (1 - cos_t));

	struct rvec v = new_rvec(x, y, z);
	return v;
}
*/

static double matrix_trace(gsl_matrix *a)
{
	int i;
	double tr = 0.0;

	assert(a->size1 == a->size2);
	for ( i=0; i<a->size1; i++ ) {
		tr += gsl_matrix_get(a, i, i);
	}
	return tr;
}

static char *add_ua(const char *inp, char ua)
{
	char *pg = malloc(64);
	if ( pg == NULL ) return NULL;
	snprintf(pg, 63, "%s_ua%c", inp, ua);
	return pg;
}


static char *get_chiral_holohedry(UnitCell *cell)
{
	LatticeType lattice = cell_get_lattice_type(cell);
	char *pg = malloc(64);
	char *pgout = 0;

	if ( pg == NULL ) return NULL;

	switch (lattice)
	{
		case L_TRICLINIC:
		pg = "1";
		break;

		case L_MONOCLINIC:
		pg = "2";
		break;

		case L_ORTHORHOMBIC:
		pg = "222";
		break;

		case L_TETRAGONAL:
		pg = "422";
		break;

		case L_RHOMBOHEDRAL:
		pg = "3_R";
		break;

		case L_HEXAGONAL:
		if ( cell_get_centering(cell) == 'H' ) {
			pg = "3_H";
		} else {
			pg = "622";
		}
		break;

		case L_CUBIC:
		pg = "432";
		break;

		default:
		pg = "error";
		break;
	}

	switch (lattice)
	{
		case L_TRICLINIC:
		case L_ORTHORHOMBIC:
		case L_RHOMBOHEDRAL:
		case L_CUBIC:
		pgout = strdup(pg);
		break;

		case L_MONOCLINIC:
		case L_TETRAGONAL:
		case L_HEXAGONAL:
		pgout = add_ua(pg, cell_get_unique_axis(cell));
		break;

		default:
		break;
	}

	return pgout;
}


static SymOpList *sym_ops_for_cell(UnitCell *cell)
{
	SymOpList *rawList;

	char *pg = get_chiral_holohedry(cell);
	rawList = get_pointgroup(pg);
	free(pg);

	return rawList;
}

static int rot_mats_are_similar(gsl_matrix *rot1, gsl_matrix *rot2,
                                gsl_matrix *sub, gsl_matrix *mul,
                                double *score, struct TakeTwoCell *cell)
{
	double tr;

	gsl_matrix_memcpy(sub, rot1);
	gsl_matrix_sub(sub, rot2);  /* sub = rot1 - rot2 */

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, sub, sub, 0.0, mul);

	tr = matrix_trace(mul);
	if (score != NULL) *score = tr;

	return (tr < cell->trace_tol);
}

static int symm_rot_mats_are_similar(gsl_matrix *rot1, gsl_matrix *rot2,
                                     struct TakeTwoCell *cell)
{
	int i;

	gsl_matrix *sub = gsl_matrix_calloc(3, 3);
	gsl_matrix *mul = gsl_matrix_calloc(3, 3);

	for (i = 0; i < cell->numOps; i++) {
		gsl_matrix *testRot = gsl_matrix_alloc(3, 3);
		gsl_matrix *symOp = cell->rotSymOps[i];

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rot1, symOp,
		0.0, testRot);

		if (rot_mats_are_similar(testRot, rot2, sub, mul, NULL, cell)) {
			gsl_matrix_free(testRot);
			gsl_matrix_free(sub);
			gsl_matrix_free(mul);
			return 1;
		}

		gsl_matrix_free(testRot);
	}

	gsl_matrix_free(sub);
	gsl_matrix_free(mul);

	return 0;
}

static void rotation_between_vectors(struct rvec a, struct rvec b,
				     gsl_matrix *twizzle)
{
	double th = rvec_angle(a, b);
	struct rvec c = rvec_cross(a, b);
	normalise_rvec(&c);
	rotation_around_axis(c, th, twizzle);
}


static void rvec_to_gsl(gsl_vector *newVec, struct rvec v)
{
	gsl_vector_set(newVec, 0, v.u);
	gsl_vector_set(newVec, 1, v.v);
	gsl_vector_set(newVec, 2, v.w);
}


struct rvec gsl_to_rvec(gsl_vector *a)
{
	struct rvec v;
	v.u = gsl_vector_get(a, 0);
	v.v = gsl_vector_get(a, 1);
	v.w = gsl_vector_get(a, 2);
	return v;
}


/* Code me in gsl_matrix language to copy the contents of the function
 * in cppxfel (IndexingSolution::createSolution). This function is quite
 * intensive on the number crunching side so simple angle checks are used
 * to 'pre-scan' vectors beforehand. */
static gsl_matrix *generate_rot_mat(struct rvec obs1, struct rvec obs2,
				    struct rvec cell1, struct rvec cell2,
				    struct TakeTwoCell *cell)
{
	gsl_matrix *fullMat;
	rvec_to_gsl(cell->vec1Tmp, cell2);

//	gsl_vector *cell2v = rvec_to_gsl(cell2);
//	gsl_vector *cell2vr = gsl_vector_calloc(3);

	normalise_rvec(&obs1);
	normalise_rvec(&obs2);
	normalise_rvec(&cell1);
	normalise_rvec(&cell2);

	/* Rotate reciprocal space so that the first simulated vector lines up
	 * with the observed vector. */
	rotation_between_vectors(cell1, obs1, cell->twiz1Tmp);

	normalise_rvec(&obs1);

	/* Multiply cell2 by rotateSpotDiffMatrix --> cell2vr */
	gsl_blas_dgemv(CblasNoTrans, 1.0, cell->twiz1Tmp, cell->vec1Tmp,
		       0.0, cell->vec2Tmp);

	/* Now we twirl around the firstAxisUnit until the rotated simulated
	 * vector matches the second observed vector as closely as possible. */

	closest_rot_mat(gsl_to_rvec(cell->vec2Tmp), obs2, obs1, cell->twiz2Tmp);

	/* We want to apply the first matrix and then the second matrix,
	 * so we multiply these. */
	fullMat = gsl_matrix_calloc(3, 3);
	gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0,
		       cell->twiz1Tmp, cell->twiz2Tmp, 0.0, fullMat);
	gsl_matrix_transpose(fullMat);

	return fullMat;
}


static int obs_vecs_share_spot(struct SpotVec *her_obs, struct SpotVec *his_obs)
{
	if ( (her_obs->her_rlp == his_obs->her_rlp) ||
	    (her_obs->her_rlp == his_obs->his_rlp) ||
	    (her_obs->his_rlp == his_obs->her_rlp) ||
	    (her_obs->his_rlp == his_obs->his_rlp) ) {
		return 1;
	}

	return 0;
}


static int obs_shares_spot_w_array(struct SpotVec *obs_vecs, int test_idx,
				   int *members, int num)
{
	int i;

	struct SpotVec *her_obs = &obs_vecs[test_idx];

	for ( i=0; i<num; i++ ) {
		struct SpotVec *his_obs = &obs_vecs[members[i]];

		int shares = obs_vecs_share_spot(her_obs, his_obs);

		if ( shares ) return 1;
	}

	return 0;
}


static int obs_vecs_match_angles(int her, int his,
                                 struct Seed **seeds, int *match_count,
                                 struct TakeTwoCell *cell)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;
	struct SpotVec *her_obs = &obs_vecs[her];
	struct SpotVec *his_obs = &obs_vecs[his];

	*match_count = 0;

	double min_angle = deg2rad(2.5);
	double max_angle = deg2rad(187.5);

	/* calculate angle between observed vectors */
	double obs_angle = rvec_angle(her_obs->obsvec, his_obs->obsvec);

	/* calculate angle between all potential theoretical vectors */

	int i, j;
	for ( i=0; i<her_obs->match_num; i++ ) {
		for ( j=0; j<his_obs->match_num; j++ ) {
			double score = 0;

			struct rvec *her_match = &her_obs->matches[i].vec;
			struct rvec *his_match = &his_obs->matches[j].vec;

			double her_dist = rvec_length(*her_match);
			double his_dist = rvec_length(*his_match);

			double theory_angle = rvec_angle(*her_match,
			                                 *his_match);

			/* is this angle a match? */

			double angle_diff = fabs(theory_angle - obs_angle);

			/* within basic tolerance? */
			if ( angle_diff != angle_diff ||
			     angle_diff > cell->angle_tol ) {
				continue;
			}

			double add = angle_diff;
			if (add == add) {
				score += add * her_dist * his_dist;
			}

			/* If the angles are too close to 0
			or 180, one axis ill-determined */
			if (theory_angle < min_angle ||
			    theory_angle > max_angle) {
				continue;
			}

			/* check that third vector adequately completes
			*  triangle */

			struct rvec theory_diff = diff_vec(*his_match, *her_match);
			struct rvec obs_diff = diff_vec(his_obs->obsvec,
			                                her_obs->obsvec);

			theory_angle = rvec_angle(*her_match,
			                          theory_diff);
			obs_angle = rvec_angle(her_obs->obsvec, obs_diff);
			angle_diff = fabs(obs_angle - theory_angle);

			double diff_dist = rvec_length(obs_diff);

			if (angle_diff > ANGLE_TOLERANCE) {
				continue;
			}

			add = angle_diff;
			if (add == add) {
				score += add * her_dist * diff_dist;
			}

			theory_angle = rvec_angle(*his_match,
			                          theory_diff);
			obs_angle = rvec_angle(his_obs->obsvec, obs_diff);

			if (fabs(obs_angle - theory_angle) > ANGLE_TOLERANCE) {
				continue;
			}

			add = angle_diff;
			if (add == add) {
				score += add * his_dist * diff_dist;
			}

			/* we add a new seed to the array */
			size_t new_size = (*match_count + 1);
			new_size *= sizeof(struct Seed);

			/* Reallocate the array to fit in another match */
			struct Seed *tmp_seeds = realloc(*seeds, new_size);

			if ( tmp_seeds == NULL ) {
				apologise();
			}

			(*seeds) = tmp_seeds;

			(*seeds)[*match_count].obs1 = her;
			(*seeds)[*match_count].obs2 = his;
			(*seeds)[*match_count].idx1 = i;
			(*seeds)[*match_count].idx2 = j;
			(*seeds)[*match_count].score = score * 1000;

                        (*match_count)++;
                }
        }

        return (*match_count > 0);
}

/* ------------------------------------------------------------------------
 * core functions regarding the meat of the TakeTwo algorithm (Level 2)
 * ------------------------------------------------------------------------*/

static signed int finish_solution(gsl_matrix *rot, struct SpotVec *obs_vecs,
                                    int *obs_members, int *match_members,
                                    int member_num, struct TakeTwoCell *cell)
{
	gsl_matrix *sub = gsl_matrix_calloc(3, 3);
	gsl_matrix *mul = gsl_matrix_calloc(3, 3);

	gsl_matrix **rotations = malloc(sizeof(*rotations)* pow(member_num, 2)
	                                - member_num);

	int i, j, count;

	count = 0;
	for ( i=0; i<1; i++ ) {
		for ( j=0; j<member_num; j++ ) {
			if (i == j) continue;
			struct SpotVec i_vec = obs_vecs[obs_members[i]];
			struct SpotVec j_vec = obs_vecs[obs_members[j]];

			struct rvec i_obsvec = i_vec.obsvec;
			struct rvec j_obsvec = j_vec.obsvec;
			struct rvec i_cellvec = i_vec.matches[match_members[i]].vec;
			struct rvec j_cellvec = j_vec.matches[match_members[j]].vec;

			rotations[count] = generate_rot_mat(i_obsvec, j_obsvec,
							    i_cellvec, j_cellvec,
							    cell);

			count++;
		}
	}

	double min_score = FLT_MAX;
	int min_rot_index = 0;

	for (i=0; i<count; i++) {
		double current_score = 0;
		for (j=0; j<count; j++) {
			double addition;
			rot_mats_are_similar(rotations[i], rotations[j],
					     sub, mul,
					     &addition, cell);

			current_score += addition;
		}

		if (current_score < min_score) {
			min_score = current_score;
			min_rot_index = i;
		}
	}

	gsl_matrix_memcpy(rot, rotations[min_rot_index]);

	for (i=0; i<count; i++) {
		gsl_matrix_free(rotations[i]);
	}

	free(rotations);
	gsl_matrix_free(sub);
	gsl_matrix_free(mul);

	return 1;
}

gsl_matrix *rot_mat_from_indices(int her, int his,
                                 int her_match, int his_match,
                                 struct TakeTwoCell *cell)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;
	struct SpotVec *her_obs = &obs_vecs[her];
	struct SpotVec *his_obs = &obs_vecs[his];
	struct rvec i_obsvec = her_obs->obsvec;
	struct rvec j_obsvec = his_obs->obsvec;
	struct rvec i_cellvec = her_obs->matches[her_match].vec;
	struct rvec j_cellvec = his_obs->matches[his_match].vec;

	gsl_matrix *mat = generate_rot_mat(i_obsvec, j_obsvec,
	                                   i_cellvec, j_cellvec, cell);

	return mat;
}

static int weed_duplicate_matches(struct Seed **seeds,
                                  int *match_count, struct TakeTwoCell *cell)
{
	int num_occupied = 0;
	gsl_matrix **old_mats = calloc(*match_count, sizeof(gsl_matrix *));

	if (old_mats == NULL)
	{
		apologise();
		return 0;
	}

	signed int i, j;

	int duplicates = 0;

	/* Now we weed out the self-duplicates from the remaining batch */

	for (i = *match_count - 1; i >= 0; i--) {
		int her_match = (*seeds)[i].idx1;
		int his_match = (*seeds)[i].idx2;

		gsl_matrix *mat;
		mat = rot_mat_from_indices((*seeds)[i].obs1, (*seeds)[i].obs2,
		                           her_match, his_match, cell);

		int found = 0;

		for (j = 0; j < num_occupied; j++) {
			if (old_mats[j] && mat &&
			    symm_rot_mats_are_similar(old_mats[j], mat, cell))
			{
				// we have found a duplicate, so flag as bad.
				(*seeds)[i].idx1 = -1;
				(*seeds)[i].idx2 = -1;
				found = 1;

				duplicates++;

				gsl_matrix_free(mat);
				break;
			}
		}

		if (!found) {
			// we have not found a duplicate, add to list.
			old_mats[num_occupied] = mat;
			num_occupied++;
		}
	}

	for (i = 0; i < num_occupied; i++) {
		if (old_mats[i]) {
			gsl_matrix_free(old_mats[i]);
		}
	}

	free(old_mats);

	return 1;
}

static signed int find_next_index(gsl_matrix *rot, int *obs_members,
				  int *match_members, int start, int member_num,
				  int *match_found, struct TakeTwoCell *cell)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;
	int obs_vec_count = cell->obs_vec_count;
	gsl_matrix *sub = gsl_matrix_calloc(3, 3);
	gsl_matrix *mul = gsl_matrix_calloc(3, 3);

	int i, j, k;

	for ( i=start; i<obs_vec_count; i++ ) {

		/* If we've considered this vector before, ignore it */
		if (obs_vecs[i].in_network == 1)
		{
			continue;
		}

		/* first we check for a shared spot - harshest condition */
		int shared = obs_shares_spot_w_array(obs_vecs, i, obs_members,
						     member_num);

		if ( !shared ) continue;

		int all_ok = 1;
		int matched = -1;

		/* Check all existing members are happy to let in the newcomer */
		for ( j=0; j<member_num && all_ok; j++ ) {

			struct SpotVec *me = &obs_vecs[i];
			struct SpotVec *you = &obs_vecs[obs_members[j]];
			struct rvec you_cell;
			you_cell = you->matches[match_members[j]].vec;

			struct rvec me_obs = me->obsvec;
			struct rvec you_obs = you->obsvec;

			int one_is_okay = 0;

			/* Loop through all possible theoretical vector
			* matches for the newcomer.. */

			for ( k=0; k<me->match_num; k++ ) {

				gsl_matrix *test_rot;

				struct rvec me_cell = me->matches[k].vec;

				test_rot = generate_rot_mat(me_obs,
					    you_obs, me_cell, you_cell,
					    cell);

				double trace = 0;
				int ok = rot_mats_are_similar(rot, test_rot,
						       sub, mul, &trace, cell);

				gsl_matrix_free(test_rot);

				if (ok) {
					one_is_okay = 1;

					/* We are only happy if the vector
					* matches for only one kind of
					* theoretical vector. We don't want to
					* accept mixtures of theoretical vector
					* matches. */
					if (matched >= 0 && k == matched) {
						*match_found = k;
					} else if (matched < 0) {
						matched = k;
					} else {
						one_is_okay = 0;
						break;
					}
				}
			}

			if (!one_is_okay) {
				all_ok = 0;
				break;
			}
		}


		if (all_ok) {
			return i;
		}
	}

	/* give up. */

	return -1;
}

/**
 * Reward target function for refining solution to all vectors in a
 * given image. Sum of exponentials of the negative distances, which
 * means that the reward decays as the distance from the nearest
 * theoretical vector reduces. */
static double obs_to_sol_score(struct TakeTwoCell *ttCell)
{
	double total = 0;
	int count = 0;
	int i;
	gsl_matrix *solution = ttCell->solution;
	gsl_matrix *full_rot = gsl_matrix_alloc(3, 3);
	rotate_gsl_by_angles(solution, ttCell->x_ang, ttCell->y_ang,
	                     ttCell->z_ang, full_rot);

	for (i = 0; i < ttCell->obs_vec_count; i++)
	{
		struct rvec *obs = &ttCell->obs_vecs[i].obsvec;
		rvec_to_gsl(ttCell->vec1Tmp, *obs);

		/* Rotate all the observed vectors by the modified soln */
		/* ttCell->vec2Tmp = 1.0 * full_rot * ttCell->vec1Tmp */
		gsl_blas_dgemv(CblasTrans, 1.0, full_rot, ttCell->vec1Tmp,
		               0.0, ttCell->vec2Tmp);
		struct rvec rotated = gsl_to_rvec(ttCell->vec2Tmp);

		int j = ttCell->obs_vecs[i].assignment;

		if (j < 0) continue;

		struct rvec *match = &ttCell->obs_vecs[i].matches[j].vec;
		struct rvec diff = diff_vec(rotated, *match);

		double length = rvec_length(diff);

		double addition = exp(-(1 / RECIP_TOLERANCE) * length);
		total += addition;
		count++;
	}

	total /= (double)-count;

	gsl_matrix_free(full_rot);

	return total;
}

/**
 * Matches every observed vector in the image to its closest theoretical
 * neighbour after applying the rotation matrix, in preparation for
 * refining the rotation matrix to match these. */
static void match_all_obs_to_sol(struct TakeTwoCell *ttCell)
{
	int i, j;
	double total = 0;
	int count = 0;
	gsl_matrix *solution = ttCell->solution;

	for (i = 0; i < ttCell->obs_vec_count; i++)
	{
		struct rvec *obs = &ttCell->obs_vecs[i].obsvec;
		rvec_to_gsl(ttCell->vec1Tmp, *obs);

		/* ttCell->vec2Tmp = 1.0 * solution * ttCell->vec1Tmp */
		gsl_blas_dgemv(CblasTrans, 1.0, solution, ttCell->vec1Tmp,
		               0.0, ttCell->vec2Tmp);
		struct rvec rotated = gsl_to_rvec(ttCell->vec2Tmp);

		double smallest = FLT_MAX;
		int assigned = -1;

		for (j = 0; j < ttCell->obs_vecs[i].match_num; j++)
		{
			struct rvec *match = &ttCell->obs_vecs[i].matches[j].vec;
			struct rvec diff = diff_vec(rotated, *match);

			double length = rvec_length(diff);
			if (length < smallest)
			{
				smallest = length;
				assigned = j;
			}
		}

		ttCell->obs_vecs[i].assignment = assigned;

		if (smallest != FLT_MAX)
		{
			double addition = exp(-(1 / RECIP_TOLERANCE) * smallest);
			total += addition;
			count++;

		}
	}

	total /= (double)count;
}

/**
 * Refines a matrix against all of the observed vectors against their
 * closest theoretical neighbour, by perturbing the matrix along the principle
 * axes until it maximises a reward function consisting of the sum of
 * the distances of individual observed vectors to their closest
 * theoretical neighbour. Reward function means that noise and alternative
 * lattices do not dominate the equation!
 **/
static void refine_solution(struct TakeTwoCell *ttCell)
{
	match_all_obs_to_sol(ttCell);

	int i, j, k;
	const int total = 3 * 3 * 3;
	const int middle = (total - 1) / 2;

	struct rvec steps[total];
	double start = obs_to_sol_score(ttCell);
	const int max_tries = 100;

	int count = 0;
	double size = ANGLE_STEP_SIZE;

	/* First we create our combinations of steps */
	for (i = -1; i <= 1; i++) {
	for (j = -1; j <= 1; j++) {
	for (k = -1; k <= 1; k++) {
		struct rvec vec = new_rvec(i, j, k);
		steps[count] = vec;
		count++;
	}
	}
	}

	struct rvec current = new_rvec(ttCell->x_ang, ttCell->y_ang,
	                               ttCell->z_ang);

	double best = start;
	count = 0;

	while (size > ANGLE_CONVERGE_SIZE && count < max_tries)
	{
		struct rvec sized[total];

		int best_num = middle;
		for (i = 0; i < total; i++)
		{
			struct rvec sized_step = steps[i];
			sized_step.u *= size;
			sized_step.v *= size;
			sized_step.w *= size;

			struct rvec new_angles = rvec_add_rvec(current,
			                                       sized_step);

			sized[i] = new_angles;

			ttCell->x_ang = sized[i].u;
			ttCell->y_ang = sized[i].v;
			ttCell->z_ang = sized[i].w;

			double score = obs_to_sol_score(ttCell);

			if (score < best)
			{
				best = score;
				best_num = i;
			}
		}

		if (best_num == middle)
		{
			size /= 2;
		}

		current = sized[best_num];
		count++;
	}

	ttCell->x_ang = 0;
	ttCell->y_ang = 0;
	ttCell->z_ang = 0;

	gsl_matrix *tmp = gsl_matrix_alloc(3, 3);
	rotate_gsl_by_angles(ttCell->solution, current.u,
	                     current.v, current.w, tmp);
	gsl_matrix_free(ttCell->solution);
	ttCell->solution = tmp;
}


static unsigned int grow_network(gsl_matrix *rot, int obs_idx1, int obs_idx2,
                                 int match_idx1, int match_idx2,
			         struct TakeTwoCell *cell)
{

	struct SpotVec *obs_vecs = cell->obs_vecs;
	int obs_vec_count = cell->obs_vec_count;
	int *obs_members;
	int *match_members;

	/* Clear the in_network status of all vectors to start */
	int i;
	for (i = 0; i < obs_vec_count; i++)
	{
		obs_vecs[i].in_network = 0;
	}

	/* indices of members of the self-consistent network of vectors */
	obs_members = malloc((cell->member_thresh+3)*sizeof(int));
	match_members = malloc((cell->member_thresh+3)*sizeof(int));
	if ( (obs_members == NULL) || (match_members == NULL) ) {
		apologise();
		return 0;
	}

	/* initialise the ones we know already */
	obs_members[0] = obs_idx1;
	obs_members[1] = obs_idx2;
	match_members[0] = match_idx1;
	match_members[1] = match_idx2;
	int member_num = 2;

	/* counter for dead ends which must not exceed MAX_DEAD_ENDS
	 * before it is reset in an additional branch */
	int dead_ends = 0;

	/* we start from 0 */
	int start = 0;

	while ( 1 ) {

		if (start > obs_vec_count) {
			free(obs_members);
			free(match_members);
			return 0;
		}

		int match_found = -1;

		signed int next_index = find_next_index(rot, obs_members,
							match_members,
							0, member_num,
							&match_found, cell);

		if ( member_num < 2 ) {
			free(obs_members);
			free(match_members);
			return 0;
		}

		if ( next_index < 0 ) {
			/* If there have been too many dead ends, give up
			 * on indexing altogether.
			 **/
			if ( dead_ends > MAX_DEAD_ENDS ) {
				break;
			}

			/* We have not had too many dead ends. Try removing
			 the last member and continue. */
			member_num--;
			dead_ends++;

			continue;
		}

		/* Elongation of the network was successful */
		obs_vecs[next_index].in_network = 1;
		obs_members[member_num] = next_index;
		match_members[member_num] = match_found;

		member_num++;

		/* If member_num is high enough, we want to return a yes */
		if ( member_num > cell->member_thresh ) break;

	}

	finish_solution(rot, obs_vecs, obs_members,
	                match_members, member_num, cell);

	free(obs_members);
	free(match_members);

	return ( member_num );
}


static unsigned int start_seed(int i, int j, int i_match, int j_match,
                               gsl_matrix **rotation, struct TakeTwoCell *cell)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;

	gsl_matrix *rot_mat;

	rot_mat = generate_rot_mat(obs_vecs[i].obsvec,
				   obs_vecs[j].obsvec,
				   obs_vecs[i].matches[i_match].vec,
				   obs_vecs[j].matches[j_match].vec,
				   cell);

	/* Try to expand this rotation matrix to a larger network */

	int member_num = grow_network(rot_mat, i, j, i_match, j_match,
	                              cell);

	/* return this matrix and if it was immediately successful */
	*rotation = rot_mat;

	return member_num;
}

static int sort_seed_by_score(const void *av, const void *bv)
{
	struct Seed *a = (struct Seed *)av;
	struct Seed *b = (struct Seed *)bv;
	return a->score > b->score;
}

static void remove_old_solutions(struct TakeTwoCell *cell,
                                 struct taketwo_private *tp)
{
	int duplicates = 0;
	struct Seed *seeds = cell->seeds;
	unsigned int total = cell->seed_count;

	/* First we remove duplicates with previous solutions */

	int i, j;
	for (i = total - 1; i >= 0; i--) {
		int her_match = seeds[i].idx1;
		int his_match = seeds[i].idx2;

		gsl_matrix *mat;
		mat = rot_mat_from_indices(seeds[i].obs1, seeds[i].obs2,
		                           her_match, his_match, cell);

		if (mat == NULL)
		{
			continue;
		}

		for (j = 0; j < tp->numPrevs; j++)
		{
			int sim = symm_rot_mats_are_similar(tp->prevSols[j],
			                                    mat, cell);

			/* Found a duplicate with a previous solution */
			if (sim)
			{
				seeds[i].idx1 = -1;
				seeds[i].idx2 = -1;
				duplicates++;
				break;
			}
		}

		gsl_matrix_free(mat);
	}

//	STATUS("Removing %i duplicates due to prev solutions.\n", duplicates);
}

static int find_seeds(struct TakeTwoCell *cell, struct taketwo_private *tp)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;
	int obs_vec_count = cell->obs_vec_count;

	/* loop round pairs of vectors to try and find a suitable
	 * seed to start building a self-consistent network of vectors
	 */
	int i, j;

	for ( i=1; i<obs_vec_count; i++ ) {

		for ( j=0; j<i; j++ ) {

			/** Only check distances which are accumulatively less
			* than the limit */
			if (obs_vecs[j].distance + obs_vecs[i].distance >
			    MAX_RECIP_DISTANCE) {
				continue;
			}

			/** Check to see if there is a shared spot - opportunity
			 * for optimisation by generating a look-up table
			 * by spot instead of by vector.
			 */
			int shared = obs_vecs_share_spot(&obs_vecs[i],
			                                 &obs_vecs[j]);
			if ( !shared ) continue;

			/* cell vector index matches stored in i, j and total
			 * number stored in int matches.
			 */
			int seed_num = 0;
			struct Seed *seeds = NULL;

			/* Check to see if any angles match from the cell
			 * vectors */
			obs_vecs_match_angles(i, j, &seeds, &seed_num, cell);

			if (seed_num == 0)
			{
				/* Nothing to clean up here */
				continue;
			}

			/* Weed out the duplicate seeds (from symmetric
			 * reflection pairs) */
			weed_duplicate_matches(&seeds, &seed_num, cell);

			/* Add all the new seeds to the full list */

			size_t new_size = cell->seed_count + seed_num;
			new_size *= sizeof(struct Seed);

			struct Seed *tmp = realloc(cell->seeds, new_size);

			if (tmp == NULL) {
				apologise();
			}

			cell->seeds = tmp;

			for (int i = 0; i < seed_num; i++)
			{
				if (seeds[i].idx1 < 0 || seeds[i].idx2 < 0)
				{
					continue;
				}

				cell->seeds[cell->seed_count] = seeds[i];
				cell->seed_count++;
			}
		}
	}

	qsort(cell->seeds, cell->seed_count, sizeof(struct Seed),
	      sort_seed_by_score);


	return 1;
}

static int start_seeds(gsl_matrix **rotation, struct TakeTwoCell *cell)
{
	struct Seed *seeds = cell->seeds;
	int seed_num = cell->seed_count;
	int member_num = 0;

	/* We have seeds! Pass each of them through the seed-starter  */
	/* If a seed has the highest achieved membership, make note...*/
	int k;
	for ( k=0; k<seed_num; k++ ) {
		int seed_idx1 = seeds[k].idx1;
		int seed_idx2 = seeds[k].idx2;

		if (seed_idx1 < 0 || seed_idx2 < 0) {
			continue;
		}

		int seed_obs1 = seeds[k].obs1;
		int seed_obs2 = seeds[k].obs2;

		member_num = start_seed(seed_obs1, seed_obs2, seed_idx1,
		                        seed_idx2, rotation, cell);

		if (member_num >= NETWORK_MEMBER_THRESHOLD) {
			free(seeds);
			return 1;
		}
	}

	free(seeds);

	return (member_num > MINIMUM_MEMBER_THRESHOLD && rotation != NULL);
}


static void set_gsl_matrix(gsl_matrix *mat, double asx, double asy, double asz,
                           double bsx, double bsy, double bsz,
                           double csx, double csy, double csz)
{
	gsl_matrix_set(mat, 0, 0, asx);
	gsl_matrix_set(mat, 0, 1, asy);
	gsl_matrix_set(mat, 0, 2, asz);
	gsl_matrix_set(mat, 1, 0, bsx);
	gsl_matrix_set(mat, 1, 1, bsy);
	gsl_matrix_set(mat, 1, 2, bsz);
	gsl_matrix_set(mat, 2, 0, csx);
	gsl_matrix_set(mat, 2, 1, csy);
	gsl_matrix_set(mat, 2, 2, csz);
}

static int generate_rotation_sym_ops(struct TakeTwoCell *ttCell)
{
	SymOpList *rawList = sym_ops_for_cell(ttCell->cell);

	/* Now we must convert these into rotation matrices rather than hkl
	 * transformations (affects triclinic, monoclinic, rhombohedral and
	 * hexagonal space groups only) */

	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	gsl_matrix *recip = gsl_matrix_alloc(3, 3);
	gsl_matrix *cart = gsl_matrix_alloc(3, 3);
	cell_get_reciprocal(ttCell->cell, &asx, &asy, &asz, &bsx, &bsy,
						&bsz, &csx, &csy, &csz);

	set_gsl_matrix(recip, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	cell_get_cartesian(ttCell->cell, &asx, &asy, &asz, &bsx, &bsy,
						&bsz, &csx, &csy, &csz);

	set_gsl_matrix(cart, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	int i, j, k;
	int numOps = num_equivs(rawList, NULL);

	ttCell->rotSymOps = malloc(numOps * sizeof(gsl_matrix *));
	ttCell->numOps = numOps;

	if (ttCell->rotSymOps == NULL) {
		apologise();
		return 0;
	}

	for (i = 0; i < numOps; i++)
	{
		gsl_matrix *symOp = gsl_matrix_alloc(3, 3);
		IntegerMatrix *op = get_symop(rawList, NULL, i);

		for (j = 0; j < 3; j++) {
		for (k = 0; k < 3; k++) {
			gsl_matrix_set(symOp, j, k, intmat_get(op, j, k));
		}
		}

		gsl_matrix *first = gsl_matrix_calloc(3, 3);
		gsl_matrix *second = gsl_matrix_calloc(3, 3);

		/* Each equivalence is of the form:
		   cartesian * symOp * reciprocal.
		   First multiplication: symOp * reciprocal */

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		               1.0, symOp, recip,
		               0.0, first);

		/* Second multiplication: cartesian * first */

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		               1.0, cart, first,
		               0.0, second);

		ttCell->rotSymOps[i] = second;

		gsl_matrix_free(symOp);
		gsl_matrix_free(first);
	}

	gsl_matrix_free(cart);
	gsl_matrix_free(recip);

	free_symoplist(rawList);

	return 1;
}

struct sortme
{
	struct TheoryVec v;
	double dist;
};

static int sort_theory_distances(const void *av, const void *bv)
{
	struct sortme *a = (struct sortme *)av;
	struct sortme *b = (struct sortme *)bv;
	return a->dist > b->dist;
}

static int match_obs_to_cell_vecs(struct TheoryVec *cell_vecs, int cell_vec_count,
				  struct TakeTwoCell *cell)
{
	struct SpotVec *obs_vecs = cell->obs_vecs;
	int obs_vec_count = cell->obs_vec_count;
	int i, j;

	for ( i=0; i<obs_vec_count; i++ ) {

		int count = 0;
		struct sortme *for_sort = NULL;

		for ( j=0; j<cell_vec_count; j++ ) {
			/* get distance for unit cell vector */
			double cell_length = rvec_length(cell_vecs[j].vec);
			double obs_length = obs_vecs[i].distance;

			/* check if this matches the observed length */
			double dist_diff = fabs(cell_length - obs_length);

			if ( dist_diff > cell->len_tol ) continue;

			/* we have a match, add to array! */

			size_t new_size = (count+1)*sizeof(struct sortme);
			for_sort = realloc(for_sort, new_size);

			if ( for_sort == NULL ) return 0;

			for_sort[count].v = cell_vecs[j];
			for_sort[count].dist = dist_diff;
			count++;

		}

		/* Pointers to relevant things */

		struct TheoryVec **match_array;
		int *match_count;

		match_array = &(obs_vecs[i].matches);
		match_count = &(obs_vecs[i].match_num);

		/* Sort in order to get most agreeable matches first */
		qsort(for_sort, count, sizeof(struct sortme), sort_theory_distances);
		*match_array = malloc(count*sizeof(struct TheoryVec));
		*match_count = count;
		for ( j=0; j<count; j++ ) {
			(*match_array)[j] = for_sort[j].v;

		}
		free(for_sort);
	}

	return 1;
}

static int compare_spot_vecs(const void *av, const void *bv)
{
	struct SpotVec *a = (struct SpotVec *)av;
	struct SpotVec *b = (struct SpotVec *)bv;
	return a->distance > b->distance;
}

static int gen_observed_vecs(struct rvec *rlps, int rlp_count,
			     struct TakeTwoCell *cell)
{
	int i, j;
	int count = 0;

	/* maximum distance squared for comparisons */
	double max_sq_length = pow(MAX_RECIP_DISTANCE, 2);

	for ( i=0; i<rlp_count-1 && count < MAX_OBS_VECTORS; i++ ) {
		for ( j=i+1; j<rlp_count; j++ ) {

			/* calculate difference vector between rlps */
			struct rvec diff = diff_vec(rlps[i], rlps[j]);

			/* are these two far from each other? */
			double sqlength = sq_length(diff);

			if ( sqlength > max_sq_length ) continue;

			count++;

			struct SpotVec *temp_obs_vecs;
			temp_obs_vecs = realloc(cell->obs_vecs,
						count*sizeof(struct SpotVec));

			if ( temp_obs_vecs == NULL ) {
				return 0;
			} else {
				cell->obs_vecs = temp_obs_vecs;

				/* initialise all SpotVec struct members */

				struct SpotVec spot_vec;
				spot_vec.obsvec = diff;
				spot_vec.distance = sqrt(sqlength);
				spot_vec.matches = NULL;
				spot_vec.assignment = -1;
				spot_vec.match_num = 0;
				spot_vec.her_rlp = &rlps[i];
				spot_vec.his_rlp = &rlps[j];

				cell->obs_vecs[count - 1] = spot_vec;
			}
		}
	}

	/* Sort such that the shortest distances are searched first. */
	qsort(cell->obs_vecs, count, sizeof(struct SpotVec), compare_spot_vecs);

	cell->obs_vec_count = count;

	return 1;
}


static int gen_theoretical_vecs(UnitCell *cell, struct TheoryVec **cell_vecs,
				int *vec_count)
{
	double a, b, c, alpha, beta, gamma;
	int h_max, k_max, l_max;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
			    &bsx, &bsy, &bsz,
			    &csx, &csy, &csz);

	SymOpList *rawList = sym_ops_for_cell(cell);

	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

	/* find maximum Miller (h, k, l) indices for a given resolution */
	h_max = MAX_RECIP_DISTANCE * a;
	k_max = MAX_RECIP_DISTANCE * b;
	l_max = MAX_RECIP_DISTANCE * c;

	int h, k, l;
	int _h, _k, _l;
	int count = 0;

	for ( h=-h_max; h<=+h_max; h++ ) {
	for ( k=-k_max; k<=+k_max; k++ ) {
	for ( l=-l_max; l<=+l_max; l++ ) {

		struct rvec cell_vec;

		/* Exclude systematic absences from centering concerns */
		if ( forbidden_reflection(cell, h, k, l) ) continue;

		int asymmetric = 0;
		get_asymm(rawList, h, k, l, &_h, &_k, &_l);

		if (h == _h && k == _k && l == _l) {
			asymmetric = 1;
		}

		cell_vec.u = h*asx + k*bsx + l*csx;
		cell_vec.v = h*asy + k*bsy + l*csy;
		cell_vec.w = h*asz + k*bsz + l*csz;

		struct TheoryVec theory;
		theory.vec = cell_vec;
		theory.asym = asymmetric;

		/* add this to our array - which may require expanding */
		count++;

		struct TheoryVec *temp_cell_vecs;
		temp_cell_vecs = realloc(*cell_vecs,
		                         count*sizeof(struct TheoryVec));

		if ( temp_cell_vecs == NULL ) {
			return 0;
		} else {
			*cell_vecs = temp_cell_vecs;
			(*cell_vecs)[count - 1] = theory;
		}
	}
	}
	}

	*vec_count = count;

	free_symoplist(rawList);

	return 1;
}


/* ------------------------------------------------------------------------
 * cleanup functions - called from run_taketwo().
 * ------------------------------------------------------------------------*/

static void cleanup_taketwo_obs_vecs(struct SpotVec *obs_vecs,
				     int obs_vec_count)
{
	int i;
	for ( i=0; i<obs_vec_count; i++ ) {
		free(obs_vecs[i].matches);
	}

	free(obs_vecs);
}

static void cleanup_taketwo_cell(struct TakeTwoCell *ttCell)
{
	/* n.b. solutions in ttCell are taken care of in the
	* partial taketwo cleanup. */
	int i;
	for ( i=0; i<ttCell->numOps; i++ ) {
		gsl_matrix_free(ttCell->rotSymOps[i]);
	}

	cleanup_taketwo_obs_vecs(ttCell->obs_vecs,
	                         ttCell->obs_vec_count);

	free(ttCell->vec1Tmp);
	free(ttCell->vec2Tmp);
	free(ttCell->twiz1Tmp);
	free(ttCell->twiz2Tmp);
	free(ttCell->rotSymOps);
}


/* ------------------------------------------------------------------------
 * external functions - top level functions (Level 1)
 * ------------------------------------------------------------------------*/

/**
 * @cell: target unit cell
 * @rlps: spot positions on detector back-projected into recripocal
 * space depending on detector geometry etc.
 * @rlp_count: number of rlps in rlps array.
 * @rot: pointer to be given an assignment (hopefully) if indexing is
 * successful.
 **/
static UnitCell *run_taketwo(UnitCell *cell, const struct taketwo_options *opts,
                             struct rvec *rlps, int rlp_count,
                             struct taketwo_private *tp)
{
	int cell_vec_count = 0;
	struct TheoryVec *theory_vecs = NULL;
	UnitCell *result;
	int success = 0;
	gsl_matrix *solution = NULL;

	/* Initialise TakeTwoCell */
	struct TakeTwoCell ttCell;
	ttCell.cell = cell;
	ttCell.seeds = NULL;
	ttCell.seed_count = 0;
	ttCell.rotSymOps = NULL;
	ttCell.obs_vecs = NULL;
	ttCell.twiz1Tmp = gsl_matrix_calloc(3, 3);
	ttCell.twiz2Tmp = gsl_matrix_calloc(3, 3);
	ttCell.vec1Tmp = gsl_vector_calloc(3);
	ttCell.vec2Tmp = gsl_vector_calloc(3);
	ttCell.numOps = 0;
	ttCell.obs_vec_count = 0;
	ttCell.solution = NULL;
	ttCell.x_ang = 0;
	ttCell.y_ang = 0;
	ttCell.z_ang = 0;

	success = generate_rotation_sym_ops(&ttCell);

	success = gen_theoretical_vecs(cell, &theory_vecs, &cell_vec_count);
	if ( !success ) return NULL;

	success = gen_observed_vecs(rlps, rlp_count, &ttCell);
	if ( !success ) return NULL;

	if ( opts->member_thresh < 0 ) {
		ttCell.member_thresh = NETWORK_MEMBER_THRESHOLD;
	} else {
		ttCell.member_thresh = opts->member_thresh;
	}

	if ( opts->len_tol < 0.0 ) {
		ttCell.len_tol = RECIP_TOLERANCE;
	} else {
		ttCell.len_tol = opts->len_tol;  /* Already in m^-1 */
	}

	if ( opts->angle_tol < 0.0 ) {
		ttCell.angle_tol = ANGLE_TOLERANCE;
	} else {
		ttCell.angle_tol = opts->angle_tol;  /* Already in radians */
	}

	if ( opts->trace_tol < 0.0 ) {
		ttCell.trace_tol = sqrt(4.0*(1.0-cos(TRACE_TOLERANCE)));
	} else {
		ttCell.trace_tol = sqrt(4.0*(1.0-cos(opts->trace_tol)));
	}

	success = match_obs_to_cell_vecs(theory_vecs, cell_vec_count,
					 &ttCell);

	free(theory_vecs);

	if ( !success ) return NULL;

	/* Find all the seeds, then take each one and extend them, returning
	* a solution if it exceeds the NETWORK_MEMBER_THRESHOLD. */
	find_seeds(&ttCell, tp);
	remove_old_solutions(&ttCell, tp);
	start_seeds(&solution, &ttCell);

	if ( solution == NULL ) {
		return NULL;
	}

	/* If we have a solution, refine against vectors in the entire image */
	ttCell.solution = solution;
	refine_solution(&ttCell);
	solution = ttCell.solution;

	int i;
	for (i = 0; i < tp->numPrevs; i++)
	{
		gsl_matrix *sol = tp->prevSols[i];

		int sim = symm_rot_mats_are_similar(sol, solution, &ttCell);
		if (sim)
		{
//			STATUS("Warning! Returning previous solution.\n");
		}
	}

	int new_size = (tp->numPrevs + 1) * sizeof(gsl_matrix *);
	gsl_matrix **tmp = realloc(tp->prevSols, new_size);

	if (!tmp) {
		apologise();
	}

	tp->prevSols = tmp;

	tp->prevSols[tp->numPrevs] = solution;
	tp->numPrevs++;

	result = transform_cell_gsl(cell, solution);
	cleanup_taketwo_cell(&ttCell);

	return result;
}

/* Cleans up the per-image information for taketwo */

static void partial_taketwo_cleanup(struct taketwo_private *tp)
{
	if (tp->prevSols != NULL)
	{
		for (int i = 0; i < tp->numPrevs; i++)
		{
			gsl_matrix_free(tp->prevSols[i]);
		}

		free(tp->prevSols);
	}

	tp->attempts = 0;
	tp->numPrevs = 0;
	tp->prevSols = NULL;

}

/* CrystFEL interface hooks */

int taketwo_index(struct image *image, const struct taketwo_options *opts,
                  void *priv)
{
	Crystal *cr;
	UnitCell *cell;
	struct rvec *rlps;
	int n_rlps = 0;
	int i;
	struct taketwo_private *tp = (struct taketwo_private *)priv;

	/* Check serial number against previous for solution tracking */
	int this_serial = image->serial;

	if (tp->serial_num == this_serial)
	{
		tp->attempts++;
	}
	else
	{
		partial_taketwo_cleanup(tp);
		tp->serial_num = this_serial;
	}

	/*
	STATUS("Indexing %i with %i attempts, %i crystals\n", this_serial, tp->attempts,
	       image->n_crystals);
	*/

	rlps = malloc((image_feature_count(image->features)+1)*sizeof(struct rvec));
	for ( i=0; i<image_feature_count(image->features); i++ ) {
		struct imagefeature *pk = image_get_feature(image->features, i);
		if ( pk == NULL ) continue;
		rlps[n_rlps].u = pk->rx;
		rlps[n_rlps].v = pk->ry;
		rlps[n_rlps].w = pk->rz;
		n_rlps++;
	}
	rlps[n_rlps].u = 0.0;
	rlps[n_rlps].v = 0.0;
	rlps[n_rlps++].w = 0.0;

	cell = run_taketwo(tp->cell, opts, rlps, n_rlps, tp);
	free(rlps);
	if ( cell == NULL ) return 0;

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}

	crystal_set_cell(cr, cell);

	image_add_crystal(image, cr);

	return 1;
}


void *taketwo_prepare(IndexingMethod *indm, UnitCell *cell)
{
	struct taketwo_private *tp;

	/* Flags that TakeTwo knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_USE_LATTICE_TYPE
	                              | INDEXING_USE_CELL_PARAMETERS;

	if ( !( (*indm & INDEXING_USE_LATTICE_TYPE)
	       && (*indm & INDEXING_USE_CELL_PARAMETERS)) )
	{
		ERROR("TakeTwo indexing requires cell and lattice type "
		      "information.\n");
		return NULL;
	}

	if ( cell == NULL ) {
		ERROR("TakeTwo indexing requires a unit cell.\n");
		return NULL;
	}

	STATUS("*******************************************************************\n");
	STATUS("*****                    Welcome to TakeTwo                   *****\n");
	STATUS("*******************************************************************\n");
	STATUS("      If you use these indexing results, please keep a roof\n");
	STATUS("           over the author's head by citing this paper.\n\n");

	STATUS("o     o     o     o     o     o     o     o     o     o     o     o\n");
	STATUS("   o     o     o     o     o     o     o     o     o     o     o   \n");
	STATUS("o                                                                 o\n");
	STATUS("   o                      The citation is:                     o   \n");
	STATUS("o           Ginn et al., Acta Cryst. (2016). D72, 956-965         o\n");
	STATUS("   o                         Thank you!                        o   \n");
	STATUS("o                                                                 o\n");
	STATUS("   o     o     o     o     o     o     o     o     o     o     o   \n");
	STATUS("o     o     o     o     o     o     o     o     o     o     o     o\n");


	STATUS("\n");

	tp = malloc(sizeof(struct taketwo_private));
	if ( tp == NULL ) return NULL;

	tp->cell = cell;
	tp->indm = *indm;
	tp->serial_num = -1;
	tp->attempts = 0;
	tp->prevSols = NULL;
	tp->numPrevs = 0;

	return tp;
}

void taketwo_cleanup(IndexingPrivate *pp)
{
	struct taketwo_private *tp = (struct taketwo_private *)pp;

	partial_taketwo_cleanup(tp);
	free(tp);
}


const char *taketwo_probe(UnitCell *cell)
{
	if ( cell_has_parameters(cell) ) return "taketwo";
	return NULL;
}
