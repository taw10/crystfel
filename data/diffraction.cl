/*
 * diffraction.cl
 *
 * GPU calculation kernel for truncated lattice diffraction
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <defs.h>
#ifndef M_PI
#define M_PI ((float)(3.14159265))
#endif


const sampler_t sampler_a = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;
const sampler_t sampler_b = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;
const sampler_t sampler_c = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT
                             | CLK_FILTER_LINEAR;


float4 quat_rot(float4 q, float4 z)
{
	float4 res;
	float t01, t02, t03, t11, t12, t13, t22, t23, t33;

	t01 = z.x*z.y;
	t02 = z.x*z.z;
	t03 = z.x*z.w;
	t11 = z.y*z.y;
	t12 = z.y*z.z;
	t13 = z.y*z.w;
	t22 = z.z*z.z;
	t23 = z.z*z.w;
	t33 = z.w*z.w;

	res.x = (1.0 - 2.0 * (t22 + t33)) * q.x
	            + (2.0 * (t12 + t03)) * q.y
	            + (2.0 * (t13 - t02)) * q.z;

	res.y =       (2.0 * (t12 - t03)) * q.x
	      + (1.0 - 2.0 * (t11 + t33)) * q.y
	            + (2.0 * (t01 + t23)) * q.z;

	res.z =       (2.0 * (t02 + t13)) * q.x
	            + (2.0 * (t23 - t01)) * q.y
	      + (1.0 - 2.0 * (t11 + t22)) * q.z;

	return res;
}


float4 get_q(int x, int y, float cx, float cy, float res, float clen, float k,
             float *ttp, float4 z, int sampling)
{
	float rx, ry, r;
	float az, tt;
	float4 q;

	rx = ((float)x - sampling*cx)/(res*sampling);
	ry = ((float)y - sampling*cy)/(res*sampling);

	r = sqrt(pow(rx, 2.0f) + pow(ry, 2.0f));

	tt = atan2(r, clen);
	*ttp = tt;

	az = atan2(ry, rx);

	q = (float4)(k*native_sin(tt)*native_cos(az),
	             k*native_sin(tt)*native_sin(az),
	             k-k*native_cos(tt), 0.0);

	return quat_rot(q, z);
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


float get_intensity(global float *intensities, float16 cell, float4 q)
{
	float hf, kf, lf;
	int h, k, l;
	int idx;

	hf = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;  /* h */
	kf = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;  /* k */
	lf = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;  /* l */

	h = round(hf);
	k = round(kf);
	l = round(lf);

	/* Return a silly value if indices are out of range */
	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		return 0.0;
	}

	h = (h>=0) ? h : h+IDIM;
	k = (k>=0) ? k : k+IDIM;
	l = (l>=0) ? l : l+IDIM;

	if ( (h>=IDIM) || (k>=IDIM) || (l>=IDIM) ) return 0.0;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);

	return intensities[idx];
}


kernel void diffraction(global float *diff, global float *tt, float klow,
                       int w, float cx, float cy,
                       float res, float clen, float16 cell,
                       global float *intensities, float4 z,
                       int xmin, int ymin, int sampling, local float *tmp,
                       float kstep,
                       read_only image2d_t func_a,
                       read_only image2d_t func_b,
                       read_only image2d_t func_c)
{
	float ttv;
	const int x = get_global_id(0) + (xmin*sampling);
	const int y = get_global_id(1) + (ymin*sampling);
	float f_lattice, I_lattice;
	float I_molecule;
	float4 q;
	const int lx = get_local_id(0);
	const int ly = get_local_id(1);
	const int lb = get_local_id(2);
	float k = klow + kstep * get_local_id(2);
	const int ax = x / sampling;
	const int ay = y / sampling;
	float intensity;

	/* Calculate value */
	q = get_q(x, y, cx, cy, res, clen, k, &ttv, z, sampling);
	f_lattice = lattice_factor(cell, q, func_a, func_b, func_c);
	I_molecule = get_intensity(intensities, cell, q);
	I_lattice = pow(f_lattice, 2.0f);

	/* Write the value to local memory */
	intensity = I_molecule * I_lattice;
	tmp[lx+sampling*ly+sampling*sampling*lb] = intensity;

	/* Memory fence */
	barrier(CLK_LOCAL_MEM_FENCE);

	/* Leader thread sums the values */
	if ( lx + ly + lb == 0 ) {

		int i;
		float sum = 0.0;
		float val;

		for ( i=0; i<sampling*sampling*get_local_size(2); i++ )
			sum += tmp[i];

		val = sum / (float)(sampling*sampling*get_local_size(2));
		diff[ax+w*ay] = val;

		/* Leader thread also records the 2theta value.
		 * This should really be averaged across all pixels, but
		 * I strongly suspect this would be a waste of time. */
		tt[ax+w*ay] = ttv;
	}
}
