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


float range(float a)
{
	a = a >= 0.0 ? a : 2.0*M_PI + a;
	a = a < 2.0*M_PI ? a : a - 2.0*M_PI;

	return a;
}


float lattice_factor(float16 cell, float4 q, int4 ncells)
{
	float f1, f2, f3, v;
	float4 Udotq;
	const int na = ncells.s0;
	const int nb = ncells.s1;
	const int nc = ncells.s2;

	Udotq.x = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;
	Udotq.y = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;
	Udotq.z = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;

	/* At exact Bragg condition, f1 = na */
	v = M_PI*Udotq.x;
	f1 = sin(v*(float)na) / sin(v);
	f1 = isnan(f1) ? na : f1;

	/* At exact Bragg condition, f2 = nb */
	v = M_PI*Udotq.y;
	f2 = sin(v*(float)nb) / sin(v);
	f2 = isnan(f2) ? nb : f2;

	/* At exact Bragg condition, f3 = nc */
	v = M_PI*Udotq.z;
	f3 = sin(v*(float)nc) / sin(v);
	f3 = isnan(f3) ? nc : f3;

	/* At exact Bragg condition, this will multiply the molecular
	 * part of the structure factor by the number of unit cells,
	 * as desired (more scattering from bigger crystal!) */
	return f1 * f2 * f3;
}


float2 get_sfac(global float2 *sfacs, float16 cell, float4 q)
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
		return 100000.0;
	}

	h = (h>=0) ? h : h+IDIM;
	k = (k>=0) ? k : k+IDIM;
	l = (l>=0) ? l : l+IDIM;

	if ( (h>=IDIM) || (k>=IDIM) || (l>=IDIM) ) return 100000.0;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);

	return sfacs[idx];
}


kernel void diffraction(global float *diff, global float *tt, float klow,
                       int w, float cx, float cy,
                       float res, float clen, float16 cell,
                       global float2 *sfacs, float4 z, int4 ncells,
                       int xmin, int ymin, int sampling, local float *tmp,
                       float kstep)
{
	float ttv;
	const int x = get_global_id(0) + (xmin*sampling);
	const int y = get_global_id(1) + (ymin*sampling);
	float f_lattice;
	float2 f_molecule;
	float4 q;
	const int lx = get_local_id(0);
	const int ly = get_local_id(1);
	const int lb = get_local_id(2);
	float k = klow + kstep * get_local_id(2);
	const int ax = x / sampling;
	const int ay = y / sampling;
	float intensity;
	float2 val;

	/* Calculate value */
	q = get_q(x, y, cx, cy, res, clen, k, &ttv, z, sampling);
	f_lattice = lattice_factor(cell, q, ncells);
	f_molecule = get_sfac(sfacs, cell, q);

	/* Write the value to local memory */
	val = f_molecule * f_lattice;
	intensity = pow(val.x, 2.0f) + pow(val.y, 2.0f);
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
