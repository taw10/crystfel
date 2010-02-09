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


#define INDMAX 30
#define IDIM (INDMAX*2 +1)


float4 get_q(int x, int y, float cx, float cy, float res, float clen, float k,
             float *ttp)
{
	float rx, ry, r;
	float ttx, tty, tt;

	rx = ((float)x - cx)/res;
	ry = ((float)y - cy)/res;

	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

	ttx = atan2(rx, clen);
	tty = atan2(ry, clen);
	tt = atan2(r, clen);

	*ttp = tt;

	return (float4)(k*sin(ttx), k*sin(tty), k-k*cos(tt), 0.0);
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
	signed int h, k, l;
	int idx;

	h = rint(cell.s0*q.x + cell.s1*q.y + cell.s2*q.z);  /* h */
	k = rint(cell.s3*q.x + cell.s4*q.y + cell.s5*q.z);  /* k */
	l = rint(cell.s6*q.x + cell.s7*q.y + cell.s8*q.z);  /* l */

	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		return 100.0;
	}

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	return sfacs[idx];
}


kernel void diffraction(global float2 *diff, global float *tt, float k,
                       int w, float cx, float cy,
                       float res, float clen, float16 cell,
                       global float2 *sfacs)
{
	float ttv;
	const int x = get_global_id(0);
	const int y = get_global_id(1);
	float f_lattice;
	float2 f_molecule;

	float4 q = get_q(x, y, cx, cy, res, clen, k, &ttv);

	f_lattice = lattice_factor(cell, q);
	f_molecule = get_sfac(sfacs, cell, q);

	diff[x+w*y] = f_molecule*f_lattice;
	tt[x+w*y] = ttv;
}
