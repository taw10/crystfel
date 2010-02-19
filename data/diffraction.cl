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


#define INDMAX 70
#define IDIM (INDMAX*2 +1)


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
             float *ttp, float4 z)
{
	float rx, ry, r;
	float ttx, tty, tt;
	float4 q;

	rx = ((float)x - cx)/res;
	ry = ((float)y - cy)/res;

	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	ttx = atan2(rx, clen);
	tty = atan2(ry, clen);
	tt = atan2(r, clen);

	*ttp = tt;

	q = (float4)(k*sin(ttx), k*sin(tty), k-k*cos(tt), 0.0);

	return quat_rot(q, z);
}


float lattice_factor(float16 cell, float4 q)
{
	float f1, f2, f3;
	float4 Udotq;
	const int na = 8;
	const int nb = 8;
	const int nc = 8;

	Udotq.x = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;
	Udotq.y = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;
	Udotq.z = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;

	/* At exact Bragg condition, f1 = na */
	f1 = sin(M_PI*(float)na*Udotq.x) / sin(M_PI*Udotq.x);

	/* At exact Bragg condition, f2 = nb */
	f2 = sin(M_PI*(float)nb*Udotq.y) / sin(M_PI*Udotq.y);

	/* At exact Bragg condition, f3 = nc */
	f3 = sin(M_PI*(float)nc*Udotq.z) / sin(M_PI*Udotq.z);

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

	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		return 100.0;
	}

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	if ( (h>=IDIM) || (k>=IDIM) || (l>=IDIM) ) return 100.0;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);

	return sfacs[idx];
}


kernel void diffraction(global float2 *diff, global float *tt, float k,
                       int w, float cx, float cy,
                       float res, float clen, float16 cell,
                       global float2 *sfacs, float4 z)
{
	float ttv;
	const int x = get_global_id(0);
	const int y = get_global_id(1);
	float f_lattice;
	float2 f_molecule;

	float4 q = get_q(x, y, cx, cy, res, clen, k, &ttv, z);

	f_lattice = lattice_factor(cell, q);
	f_molecule = get_sfac(sfacs, cell, q);

	diff[x+w*y] = f_molecule * f_lattice;
	tt[x+w*y] = ttv;
}
