/*
 * taketwo.c
 *
 * Rewrite of TakeTwo algorithm (Acta D72 (8) 956-965) for CrystFEL
 *
 * Copyright © 2016 Helen Ginn
 * Copyright © 2016 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2016 Helen Ginn <helen@strubi.ox.ac.uk>
 *   2016 Thomas White <taw@physics.org>
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>

#include "cell-utils.h"
#include "index.h"
#include "taketwo.h"

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
	struct rvec *matches;
	int match_num;
	double distance;
	struct rvec *her_rlp;
	struct rvec *his_rlp;
};


struct taketwo_private
{
	IndexingMethod indm;
	float          *ltl;
	UnitCell       *cell;
};


/* Maximum distance between two rlp sizes to consider info for indexing */
#define MAX_RECIP_DISTANCE 0.15

/* Tolerance for two lengths in reciprocal space to be considered the same */
#define RECIP_TOLERANCE 0.001

/* Threshold for network members to consider a potential solution */
#define NETWORK_MEMBER_THRESHOLD 20

/* Maximum network members (obviously a solution so should stop) */
#define MAX_NETWORK_MEMBERS 100

/* Maximum dead ends for a single branch extension during indexing */
#define MAX_DEAD_ENDS 5

/* Tolerance for two angles to be considered the same */
#define ANGLE_TOLERANCE (deg2rad(1.0))

/** TODO:
 *
 * - May need to be capable of playing with the tolerances/#defined stuff.
 */


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
	c.w = a.u*b.v - a.v-b.u;

	return c;
}


/* ------------------------------------------------------------------------
 * functions called under the core functions, still specialised (Level 3)
 * ------------------------------------------------------------------------*/

static gsl_matrix *rotation_around_axis(struct rvec c, double th)
{
	double omc = 1.0 - cos(th);
	double s = sin(th);
	gsl_matrix *res = gsl_matrix_alloc(3, 3);

	gsl_matrix_set(res, 0, 0, cos(th) + c.u*c.u*omc);
	gsl_matrix_set(res, 0, 1, c.u*c.v*omc - c.w*s);
	gsl_matrix_set(res, 0, 2, c.u*c.w*omc + c.v*s);
	gsl_matrix_set(res, 1, 0, c.u*c.v*omc + c.w*s);
	gsl_matrix_set(res, 1, 1, cos(th) + c.v*c.v*omc);
	gsl_matrix_set(res, 1, 2, c.v*c.w*omc - c.u*s);
	gsl_matrix_set(res, 2, 0, c.w*c.u*omc - c.v*s);
	gsl_matrix_set(res, 2, 1, c.w*c.v*omc + c.u*s);
	gsl_matrix_set(res, 2, 2, cos(th) + c.w*c.w*omc);

	return res;
}


/* Rotate vector (vec1) around axis (axis) by angle theta. Find value of
 * theta for which the angle between (vec1) and (vec2) is minimised.
 * Behold! Finally an analytical solution for this one. Assuming
 * that @result has already been allocated. Will upload the maths to the
 * shared Google drive. */
static gsl_matrix *closest_rot_mat(struct rvec vec1, struct rvec vec2,
                                   struct rvec axis)
{
	/* Let's have unit vectors */
	normalise_rvec(&vec1);
	normalise_rvec(&vec2);
	normalise_rvec(&axis);

	/* Redeclaring these to try and maintain readability and
	 * check-ability against the maths I wrote down */
	double a = vec2.u; double b = vec2.w; double c = vec2.v;
	double p = vec1.u; double q = vec1.w; double r = vec1.v;
	double x = axis.u; double y = axis.w; double z = axis.v;

	/* Components in handwritten maths online when I upload it */
	double A = a*(p*x*x - p + x*y*q + x*z*r) +
	           b*(p*x*y + q*y*y - q + r*y*z) +
	           c*(p*x*z + q*y*z + r*z*z - r);

	double B = a*(y*r - z*q) + b*(p*z - r*x) + c*(q*x - p*y);

	double tan_theta = - B / A;

	/* Not sure why I always have to take the + M_PI solution. Work
	 * this one out. But it always works!? */
	double theta = atan(tan_theta) + M_PI;

	/* Return an identity matrix which has been rotated by
	 * theta around "axis" */
	return rotation_around_axis(axis, theta);
}


static double matrix_trace(gsl_matrix *a)
{
	int i;
	double tr = 0.0;

	assert(a->size1 == a->size2);
	for ( i=0; i<a->size1; a++ ) {
		tr += gsl_matrix_get(a, i, i);
	}
	return tr;
}


static int rot_mats_are_similar(gsl_matrix *rot1, gsl_matrix *rot2)
{
	double tr;
	gsl_matrix *sub;
	gsl_matrix *mul;

	sub = gsl_matrix_calloc(3, 3);
	gsl_matrix_memcpy(sub, rot1);
	gsl_matrix_sub(sub, rot2);  /* sub = rot1 - rot2 */

	mul = gsl_matrix_calloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, sub, sub, 0.0, mul);

	tr = matrix_trace(mul);
	gsl_matrix_free(mul);

	return tr < sqrt(4.0*(1.0-cos(ANGLE_TOLERANCE)));;
}


static gsl_matrix *rotation_between_vectors(struct rvec a, struct rvec b)
{
	double th = rvec_angle(a, b);
	struct rvec c = rvec_cross(a, b);
	return rotation_around_axis(c, th);
}


static gsl_vector *rvec_to_gsl(struct rvec v)
{
	gsl_vector *a = gsl_vector_alloc(3);
	gsl_vector_set(a, 0, v.u);
	gsl_vector_set(a, 1, v.v);
	gsl_vector_set(a, 2, v.w);
	return a;
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
                                    struct rvec cell1, struct rvec cell2)
{
	gsl_matrix *rotateSpotDiffMatrix;
	gsl_matrix *secondTwizzleMatrix;
	gsl_matrix *fullMat;
	gsl_vector *cell2v = rvec_to_gsl(cell2);
	gsl_vector *cell2vr = gsl_vector_calloc(3);

	/* Rotate reciprocal space so that the first simulated vector lines up
	 * with the observed vector. */
	rotateSpotDiffMatrix = rotation_between_vectors(cell1, obs1);

	normalise_rvec(&obs1);

	/* Multiply cell2 by rotateSpotDiffMatrix --> cell2vr */
	gsl_blas_dgemv(CblasNoTrans, 1.0,  rotateSpotDiffMatrix, cell2v,
	               0.0, cell2vr);

	/* Now we twirl around the firstAxisUnit until the rotated simulated
	 * vector matches the second observed vector as closely as possible. */
	secondTwizzleMatrix = closest_rot_mat(gsl_to_rvec(cell2vr),
	                                      obs2, obs1);

	/* We want to apply the first matrix and then the second matrix,
	 * so we multiply these. */
	fullMat = gsl_matrix_calloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
	               secondTwizzleMatrix, rotateSpotDiffMatrix, 0.0, fullMat);

	return fullMat;
}


static int obs_vecs_share_spot(struct SpotVec *her_obs, struct SpotVec *his_obs)
{
	/* FIXME: Disgusting... can I tone this down a bit? */
	if ( (her_obs->her_rlp == his_obs->her_rlp) ||
	     (her_obs->her_rlp == his_obs->his_rlp) ||
	     (her_obs->his_rlp == his_obs->her_rlp) ||
	     (her_obs->his_rlp == his_obs->his_rlp) ) {
		return 1;
	}

	return 0;
}


static int obs_shares_spot_w_array(struct SpotVec *obs_vecs, int test_idx,
                                   int *members[MAX_NETWORK_MEMBERS], int num)
{
	int i;
	struct SpotVec *her_obs = &obs_vecs[test_idx];

	for ( i=0; i<num; i++ ) {
		struct SpotVec *his_obs = &obs_vecs[i];

		int shares = obs_vecs_share_spot(her_obs, his_obs);

		if ( shares ) {
			return 1;
		}
	}

	return 0;
}


/** Note: this could potentially (and cautiously) converted to comparing
 * cosines rather than angles, to lose an "acos" but different parts of the
 * cosine graph are more sensitive than others, so may be a trade off... or not.
 */
static int obs_vecs_match_angles(struct SpotVec *her_obs,
                          struct SpotVec *his_obs, int *her_match_idx,
                          int *his_match_idx)
{
	int i, j;

	/* calculate angle between observed vectors */
	double obs_angle = rvec_angle(her_obs->obsvec, his_obs->obsvec);

	/* calculate angle between all potential theoretical vectors */

	for ( i=0; i<her_obs->match_num; i++ ) {
	for ( j=0; j<his_obs->match_num; j++ ) {

		struct rvec *her_match = &her_obs->matches[i];
		struct rvec *his_match = &his_obs->matches[j];

		double theory_angle = rvec_angle(*her_match, *his_match);

		/* is this angle a match? */

		double angle_diff = fabs(theory_angle - obs_angle);

		if ( angle_diff < ANGLE_TOLERANCE ) {
			*her_match_idx = i;
			*his_match_idx = j;

			return 1;
		}
	}
	}

	return 0;
}


static int obs_angles_match_array(struct SpotVec *obs_vecs, int test_idx,
                                  int *members[MAX_NETWORK_MEMBERS], int num)
{
	/* note: this is just a preliminary check to reduce unnecessary
	 * computation later down the line, but is not entirely accurate.
	 * For example, I have not checked that the 'matching cell vector'
	 * is identical - too much faff.
	 **/

	int i;
	struct SpotVec *her_obs = &obs_vecs[test_idx];

	for ( i=0; i<num; i++ ) {
		struct SpotVec *his_obs = &obs_vecs[i];

		/* placeholders, but results are ignored */
		int idx1, idx2;

		/* check our test vector matches existing network member */

		int matches = obs_vecs_match_angles(her_obs, his_obs,
						    &idx1, &idx2);

		if ( !matches )
		{
			return 0;
		}
	}

	return 1;
}


/* ------------------------------------------------------------------------
 * core functions regarding the meat of the TakeTwo algorithm (Level 2)
 * ------------------------------------------------------------------------*/

static int find_next_index(gsl_matrix *rot, struct SpotVec *obs_vecs,
                           int obs_vec_count, int **members,
                           int start, int member_num)
{
	int i;

	for ( i=start; i<obs_vec_count; i++ ) {

		/* first we check for a shared spot - harshest condition */
		int shared = obs_shares_spot_w_array(obs_vecs, i, members,
		                                     member_num);

		if ( !shared ) {
			continue;
		}

		/* now we check that angles between all vectors match */

		int matches = obs_angles_match_array(obs_vecs, i, members,
		                                     member_num);

		if ( !matches ) {
			continue;
		}

		/* final test: does the corresponding rotation matrix
		 * match the others? NOTE: have not tested to see if
		 * every combination of test/existing member has to be checked
		 * so for now, just sticking to the first two...
		 */

		/* need to grab the matching vector index */

		int member_idx, test_idx;

		obs_vecs_match_angles(&obs_vecs[(*members)[0]], &obs_vecs[i],
		                      &member_idx, &test_idx);

		struct rvec *test_match = &obs_vecs[i].matches[test_idx];
		struct rvec *member_match;
		member_match = &obs_vecs[(*members)[0]].matches[member_idx];

		int j;

		/* if ok is set to 0, give up on this vector before
		 * checking the next value of j
		 */
		int ok = 1;

		for ( j=0; j<2 && ok; j++ ) {
			gsl_matrix *test_rot = gsl_matrix_calloc(3, 3);

			int j_idx = (*members)[j];
			test_rot = generate_rot_mat(obs_vecs[j_idx].obsvec,
			                            obs_vecs[i].obsvec,
			                            *member_match,
			                            *test_match);

			ok = rot_mats_are_similar(rot, test_rot);
		}

		if ( !ok ) {
			continue;
		}

		/* successful branch - return to calling function...*/

		return i;
	}

	/* give up. */

	return -1;
}


static int grow_network(gsl_matrix *rot, struct SpotVec *obs_vecs,
                        int obs_vec_count, int obs_idx1, int obs_idx2)
{
	/* indices of members of the self-consistent network of vectors */
	int members[MAX_NETWORK_MEMBERS];

	/* initialise the ones we know already */
	members[0] = obs_idx1;
	members[1] = obs_idx2;
	int member_num = 2;

	/* counter for dead ends which must not exceed MAX_DEAD_ENDS
	 * before it is reset in an additional branch */
	int dead_ends = 0;

	/* we can start from after the 2nd observed vector in the seed */
	int start = obs_idx2 + 1;

	while ( 1 ) {
		/* There must be a better way of passing int *[100] pointer
		 * to an int ** parameter, but I don't know...
		 **/
		int next_index = find_next_index(rot, obs_vecs, obs_vec_count,
	                                         (int **)&members,
	                                         start, member_num);

		if ( member_num < 2 ) return 0;

		if ( next_index < 0 ) {
			/* If there have been too many dead ends, give up
			 * on indexing altogether.
			 **/
			if ( dead_ends > MAX_DEAD_ENDS ) break;

			/* We have not had too many dead ends. Try removing
			   the last member and continue. */
			start = members[member_num - 1] + 1;
			member_num--;
			dead_ends++;

			continue;
		}

		/* we have elongated membership - so reset dead_ends counter */
		dead_ends = 0;

		members[member_num] = next_index;
		start = next_index + 1;
		member_num++;

		/* If member_num is high enough, we want to return a yes */

		if ( member_num > NETWORK_MEMBER_THRESHOLD ) {
			break;
		}
	}

	/* Deal with this shit after coffee */

	return ( member_num > NETWORK_MEMBER_THRESHOLD );
}


static int find_seed_and_network(struct SpotVec *obs_vecs, int obs_vec_count,
                                 gsl_matrix **rotation)
{
	/* loop round pairs of vectors to try and find a suitable
	 * seed to start building a self-consistent network of vectors
	 */
	int i, j;

	for ( i=0; i<obs_vec_count-1; i++ ) {
	for ( j=i+1; j<obs_vec_count; j++ ) {

		/** Check to see if there is a shared spot - opportunity
		  * for optimisation by generating a look-up table
		  * by spot instead of by vector.
		  */
		int shared = obs_vecs_share_spot(&obs_vecs[i], &obs_vecs[j]);

		if ( !shared ) {
			continue;
		}

		/* cell vector "matches" index for i, j respectively */
		int i_idx = -1;
		int j_idx = -1;

		/* Check to see if any angles match from the cell vectors */

		int match = 0;
		match = obs_vecs_match_angles(&obs_vecs[i], &obs_vecs[j],
					     &i_idx, &j_idx);

		if ( !match ) {
			continue;
		}

		/* We have a seed! Generate a matrix based on this solution */

		gsl_matrix *rot_mat = gsl_matrix_calloc(3, 3);

		rot_mat = generate_rot_mat(obs_vecs[i].obsvec,
		                           obs_vecs[j].obsvec,
		                           obs_vecs[i].matches[i_idx],
		                           obs_vecs[j].matches[j_idx]);

		/* try to expand this rotation matrix to a larger network */

		int success = grow_network(rot_mat, obs_vecs, obs_vec_count,
		                           i_idx, j_idx);

		/* return this matrix or free it and try again */

		if ( success ) {
			*rotation = rot_mat;
			return 1;
		} else {
			gsl_matrix_free(rot_mat);
		}
	}
	} /* yes this } is meant to be here */

	return 0;
}


static int match_obs_to_cell_vecs(struct rvec *cell_vecs, int cell_vec_count,
                                  struct SpotVec *obs_vecs, int obs_vec_count)
{
	int i, j;

	/* Now I'm definitely bending the indentation rules! */
	for ( i=0; i<obs_vec_count; i++ ) {
		int count = 0;

	for ( j=0; j<cell_vec_count; j++ ) {

		/* get distance for unit cell vector */
		double cell_length = rvec_length(cell_vecs[j]);
		double obs_length = obs_vecs[i].distance;

		/* check if this matches the observed length */
		double dist_diff = fabs(cell_length - obs_length);

		if ( dist_diff > RECIP_TOLERANCE ) {
			continue;
		}

		/* we have a match, add to array! */

		count++;
		int new_size = count*sizeof(struct rvec);
		struct rvec *temp_matches;
		temp_matches = realloc(obs_vecs[i].matches, new_size);

		if ( temp_matches == NULL ) {
			return 0;
		} else {
			obs_vecs[i].matches = temp_matches;
			temp_matches[count - 1] = cell_vecs[j];
			obs_vecs[i].match_num = count;
		}
	}
	}

	return 1;
}


static int gen_observed_vecs(struct rvec *rlps, int rlp_count,
                             struct SpotVec **obs_vecs, int *obs_vec_count)
{
	int i, j;
	int count = 0;

	/* maximum distance squared for comparisons */
	double max_sq_length = pow(MAX_RECIP_DISTANCE, 2);

	/* Indentation... bending the rules a bit? */
	for ( i=0; i<rlp_count; i++ ) {
	for ( j=0; j<rlp_count; j++ ) {

		/* calculate difference vector between rlps */
		struct rvec diff = diff_vec(rlps[i], rlps[j]);

		/* are these two far from each other? */
		double sqlength = sq_length(diff);

		if ( sqlength > max_sq_length ) {
			continue;
		}

		count++;

		struct SpotVec *temp_obs_vecs;
		temp_obs_vecs = realloc(*obs_vecs,
		                        count*sizeof(struct SpotVec));

		if ( temp_obs_vecs == NULL ) {
			return 0;
		} else {
			*obs_vecs = temp_obs_vecs;

			/* initialise all SpotVec struct members */

			struct SpotVec spot_vec;
			spot_vec.obsvec = diff;
			spot_vec.distance = sqrt(sqlength);
			spot_vec.matches = NULL;
			spot_vec.match_num = 0;
			spot_vec.her_rlp = &rlps[i];
			spot_vec.his_rlp = &rlps[j];

			(*obs_vecs)[count - 1] = spot_vec;
		}
	}
	}

	*obs_vec_count = count;

	return 1;
}


static int gen_theoretical_vecs(UnitCell *cell, struct rvec **cell_vecs,
                                int *vec_count)
{
	double a, b, c, alpha, beta, gamma;
	int h_max, k_max, l_max;

	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

	/* find maximum Miller (h, k, l) indices for a given resolution */
	h_max = MAX_RECIP_DISTANCE / a + 1;
	k_max = MAX_RECIP_DISTANCE / b + 1;
	l_max = MAX_RECIP_DISTANCE / c + 1;

	int h, k, l;
	int count = 0;

	for ( h=-h_max; h<=+h_max; h++ ) {
	for ( k=-k_max; k<=+k_max; k++ ) {
	for ( l=-l_max; l<=+l_max; l++ ) {

		/* Exclude systematic absences from centering concerns */
		if ( forbidden_reflection(cell, h, k, l) ) {
			continue;
		}

		struct rvec cell_vec = new_rvec(h, k, l);

		/* FIXME: transform int (h, k, l) to reciprocal coordinates.
		 * Don't want to do this manually if there is already a
		 * function in libcrystfel to do this. Would like to map
		 * Miller index (5, 13, 2) onto (0.05, 0.13, 0.02) for example,
		 * if I had a 100 x 100 x 100 Å cubic cell. */

		cell_vec = cell_vec;

		/* Assumption that "cell_vec" now has transformed coords */

		/* add this to our array - which may require expanding */
		count++;

		struct rvec *temp_cell_vecs;
		temp_cell_vecs = realloc(*cell_vecs, count*sizeof(struct rvec));

		if ( temp_cell_vecs == NULL ) {
			return 0;
		} else {
			*cell_vecs = temp_cell_vecs;
			(*cell_vecs)[count - 1] = cell_vec;
		}
	}
	}
	}

	*vec_count = count;

	return 1;
}


static void generate_basis_vectors(UnitCell *cell, gsl_matrix *rot,
                                   struct rvec *a_star, struct rvec *b_star,
                                   struct rvec *c_star)
{
        /* FIXME: more matrix stuff - multiply cell matrix by rotation matrix
         * and extract the reciprocal axes from the definition of the matrix.
         */
}

/* ------------------------------------------------------------------------
 * cleanup functions - called from run_taketwo().
 * ------------------------------------------------------------------------*/

static void cleanup_taketwo_cell_vecs(struct rvec *cell_vecs)
{
	free(cell_vecs);
}


static void cleanup_taketwo_obs_vecs(struct SpotVec *obs_vecs,
                                    int obs_vec_count)
{
	int i;
	for ( i=0; i<obs_vec_count; i++ ) {
		free(obs_vecs[i].matches);
	}

	free(obs_vecs);
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
int run_taketwo(UnitCell *cell, struct rvec *rlps,
                int rlp_count, struct rvec *a_star, struct rvec *b_star,
                struct rvec *c_star)
{
	int cell_vec_count = 0;
	struct rvec *cell_vecs = NULL;

	int success = 0;
	success = gen_theoretical_vecs(cell, &cell_vecs, &cell_vec_count);

	if ( !success ) {
		apologise();
		return 0;
	}

	int obs_vec_count = 0;
	struct SpotVec *obs_vecs = NULL;

	success = gen_observed_vecs(rlps, rlp_count, &obs_vecs, &obs_vec_count);

	if ( !success ) {
		apologise();
		return 0;
	}

	success = match_obs_to_cell_vecs(cell_vecs, cell_vec_count,
	                                 obs_vecs, obs_vec_count);

	if ( !success ) {
		apologise();
		return 0;
	}

	cleanup_taketwo_cell_vecs(cell_vecs);

	gsl_matrix *solution = NULL;

	find_seed_and_network(obs_vecs, obs_vec_count, &solution);

        generate_basis_vectors(cell, solution, a_star, b_star, c_star);

	cleanup_taketwo_obs_vecs(obs_vecs, obs_vec_count);

	return (solution != NULL);
}


/* CrystFEL interface hooks */

int taketwo_index(struct image *image, IndexingPrivate *ipriv)
{
	struct taketwo_private *tp = (struct taketwo_private *)ipriv;
	return 1;
}


IndexingPrivate *taketwo_prepare(IndexingMethod *indm, UnitCell *cell,
                                 struct detector *det, float *ltl)
{
	struct taketwo_private *tp;
	int need_cell = 0;

	if ( *indm & INDEXING_CHECK_CELL_COMBINATIONS ) need_cell = 1;
	if ( *indm & INDEXING_CHECK_CELL_AXES ) need_cell = 1;

	if ( need_cell && !cell_has_parameters(cell) ) {
		ERROR("Altering your TakeTwo flags because cell parameters were"
		      " not provided.\n");
		*indm &= ~INDEXING_CHECK_CELL_COMBINATIONS;
		*indm &= ~INDEXING_CHECK_CELL_AXES;
	}

	/* Flags that TakeTwo knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_CELL_COMBINATIONS
	       | INDEXING_CHECK_CELL_AXES | INDEXING_CHECK_PEAKS
	       | INDEXING_CONTROL_FLAGS;

	tp = malloc(sizeof(struct taketwo_private));
	if ( tp == NULL ) return NULL;

	tp->ltl = ltl;
	tp->cell = cell;
	tp->indm = *indm;

	return (IndexingPrivate *)tp;
}


void taketwo_cleanup(IndexingPrivate *pp)
{
	struct taketwo_private *tp = (struct taketwo_private *)pp;
	free(tp);
}
