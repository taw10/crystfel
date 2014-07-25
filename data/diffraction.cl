/*
 * diffraction.cl
 *
 * GPU calculation kernel for truncated lattice diffraction
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
 *   2013      Alexandra Tolstikova
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


/* Maxmimum index to hold values up to (can be increased if necessary)
 * WARNING: Altering this value constitutes an ABI change, and means you must
 * update src/pattern_sim.h then recompile and reinstall everything. */
#define INDMAX 130
#define IDIM (INDMAX*2 +1)

#ifndef M_PI
#define M_PI ((float)(3.14159265))
#endif

const sampler_t sampler_a = CLK_NORMALIZED_COORDS_TRUE
                             | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;
const sampler_t sampler_b = CLK_NORMALIZED_COORDS_TRUE
                             | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;
const sampler_t sampler_c = CLK_NORMALIZED_COORDS_TRUE
                             | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;


float4 get_q(float fs, float ss, float res, float clen, float k,
             float corner_x, float corner_y,
             float fsx, float fsy, float ssx, float ssy)
{
	float rx, ry, r;
	float az, tt;
	float4 q;
	float xs, ys;
	float kx, ky, kz;

	xs = fs*fsx + ss*ssx;
	ys = fs*fsy + ss*ssy;

	rx = (xs + corner_x) / res;
	ry = (ys + corner_y) / res;

	r = sqrt(pow(rx, 2.0f) + pow(ry, 2.0f));

	tt = atan2(r, clen);

	az = atan2(ry, rx);

	kx = k*native_sin(tt)*native_cos(az);
	ky = k*native_sin(tt)*native_sin(az);
	kz = k*(native_cos(tt)-1.0);

	q = (float4)(kx, ky, kz, 0.0);

	return q;
}


float lattice_factor(float16 cell, float4 q,
                     read_only image2d_t func_a,
                     read_only image2d_t func_b,
                     read_only image2d_t func_c)
{
	float f1, f2, f3, v;
	float4 Udotq;

	Udotq.x = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;
	Udotq.y = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;
	Udotq.z = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;

	/* Look up values from precalculated sinc() table */
	f1 = read_imagef(func_a, sampler_a, (float2)(Udotq.x, 0.0)).s0;
	f2 = read_imagef(func_b, sampler_b, (float2)(Udotq.y, 0.0)).s0;
	f3 = read_imagef(func_c, sampler_c, (float2)(Udotq.z, 0.0)).s0;

	return f1 * f2 * f3;
}


float lookup_intensity(global float *intensities,
                       signed int h, signed int k, signed int l)
{
	int idx;

	/* Out of range? */
	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		return 0.0;
	}

	h = (h>=0) ? h : h+IDIM;
	k = (k>=0) ? k : k+IDIM;
	l = (l>=0) ? l : l+IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);

	return intensities[idx];
}


float lookup_flagged_intensity(global float *intensities, global float *flags,
                               signed int h, signed int k, signed int l)
{
	return lookup_intensity(intensities, h, k, l)
	        * lookup_intensity(flags, h, k, l);
}


float molecule_factor(global float *intensities, global float *flags,
                      float16 cell, float4 q)
{
	float hf, kf, lf;
	int h, k, l;
	float val = 0.0;

	#ifdef FLAT_INTENSITIES
	return 100.0;
	#else

	hf = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;  /* h */
	kf = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;  /* k */
	lf = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;  /* l */

	h = round(hf);
	k = round(kf);
	l = round(lf);

	/* Symmetry stuff goes here */
	INSERT_HERE

	return val;
	#endif /* FLAT_INTENSITIIES */
}


kernel void diffraction(global float *diff, float k, float weight,
                        int w, float corner_x, float corner_y,
                        float fsx, float fsy, float ssx, float ssy,
                        float res, float clen, float16 cell,
                        global float *intensities, global float *flags,
                        read_only image2d_t func_a,
                        read_only image2d_t func_b,
                        read_only image2d_t func_c,
                        local float *tmp)
{
	float fs, ss;
	float f_lattice, I_lattice;
	float I_molecule;
	float4 q;
	const int ls0 = get_local_size(0);
	const int ls1 = get_local_size(1);
	const int li0 = get_local_id(0);
	const int li1 = get_local_id(1);
	const int ls = ls0 * ls1;

	/* Calculate fractional coordinates in fs/ss */
	fs = convert_float(get_global_id(0)) / convert_float(ls0);
	ss = convert_float(get_global_id(1)) / convert_float(ls1);

	/* Get the scattering vector */
	q = get_q(fs, ss, res, clen, k,
	          corner_x, corner_y, fsx, fsy, ssx, ssy);

	/* Calculate the diffraction */
	f_lattice = lattice_factor(cell, q, func_a, func_b, func_c);
	I_molecule = molecule_factor(intensities, flags, cell, q);
	I_lattice = pow(f_lattice, 2.0f);

	tmp[li0 + ls0*li1] = I_molecule * I_lattice;

	barrier(CLK_LOCAL_MEM_FENCE);

	/* First thread in group sums the samples */
	if ( li0 + li1 == 0 ) {

		int i;
		float sum = 0.0;
		float val;
		int idx;

		idx = convert_int_rtz(fs) + w*convert_int_rtz(ss);

		for ( i=0; i<ls; i++ ) sum += tmp[i];

		val = weight * sum / convert_float(ls);
		diff[idx] = val;

	}

}
