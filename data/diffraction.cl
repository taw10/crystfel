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


/* Maxmimum index to hold values up to (can be increased if necessary)
 * WARNING: Altering this value constitutes an ABI change, and means you must
 * update src/pattern_sim.h then recompile and reinstall everything. */
#define INDMAX 200

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
             float *ttp, float corner_x, float corner_y,
             float fsx, float fsy, float ssx, float ssy,
             float xdiv, float ydiv)
{
	float rx, ry, r;
	float az, tt;
	float4 q;
	float xs, ys;
	float kx, ky, kz;
	float kxn, kyn, kzn;

	xs = fs*fsx + ss*ssx;
	ys = fs*fsy + ss*ssy;

	rx = (xs + corner_x) / res;
	ry = (ys + corner_y) / res;

	r = sqrt(pow(rx, 2.0f) + pow(ry, 2.0f));

	tt = atan2(r, clen);
	*ttp = tt;

	az = atan2(ry, rx);

	kxn = k*native_sin(tt)*native_cos(az);
	kyn = k*native_sin(tt)*native_sin(az);
	kzn = k*(native_cos(tt)-1.0);

	/* x divergence */
	kx =  kxn*native_cos(xdiv) +kzn*native_sin(xdiv);
	ky =  kyn;
	kz = -kxn*native_sin(xdiv) +kzn*native_cos(xdiv);
	kxn = kx;  kyn = ky;  kzn = kz;

	/* y divergence */
	kx =  kxn;
	ky =  kyn*cos(ydiv) +kzn*sin(ydiv);
	kz = -kyn*sin(ydiv) +kzn*cos(ydiv);

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
	signed int i;

	#ifdef FLAT_INTENSITIES
	return 1.0e5;
	#else

	hf = cell.s0*q.x + cell.s1*q.y + cell.s2*q.z;  /* h */
	kf = cell.s3*q.x + cell.s4*q.y + cell.s5*q.z;  /* k */
	lf = cell.s6*q.x + cell.s7*q.y + cell.s8*q.z;  /* l */

	h = round(hf);
	k = round(kf);
	l = round(lf);
	i = -h-k;

	#ifdef PG1
	val += lookup_flagged_intensity(intensities, flags, h, k, l);
	#endif /* PG1 */

	#ifdef PG1BAR
	val += lookup_flagged_intensity(intensities, flags,  h,  k,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k, -l);
	#endif /* PG1BAR */

	#ifdef PG6
	val += lookup_flagged_intensity(intensities, flags,  h,  k,  l);
	val += lookup_flagged_intensity(intensities, flags,  i,  h,  l);
	val += lookup_flagged_intensity(intensities, flags,  k,  i,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k,  l);
	val += lookup_flagged_intensity(intensities, flags, -i, -h,  l);
	val += lookup_flagged_intensity(intensities, flags, -k, -i,  l);
	#endif /* PG6 */

	#ifdef PG6M
	val += lookup_flagged_intensity(intensities, flags,  h,  k,  l);
	val += lookup_flagged_intensity(intensities, flags,  i,  h,  l);
	val += lookup_flagged_intensity(intensities, flags,  k,  i,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k,  l);
	val += lookup_flagged_intensity(intensities, flags, -i, -h,  l);
	val += lookup_flagged_intensity(intensities, flags, -k, -i,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k, -l);
	val += lookup_flagged_intensity(intensities, flags, -i, -h, -l);
	val += lookup_flagged_intensity(intensities, flags, -k, -i, -l);
	val += lookup_flagged_intensity(intensities, flags,  h,  k, -l);
	val += lookup_flagged_intensity(intensities, flags,  i,  h, -l);
	val += lookup_flagged_intensity(intensities, flags,  k,  i, -l);
	#endif /* PG6M */

	#ifdef PG6MMM
	val += lookup_flagged_intensity(intensities, flags,  h,  k,  l);
	val += lookup_flagged_intensity(intensities, flags,  i,  h,  l);
	val += lookup_flagged_intensity(intensities, flags,  k,  i,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k,  l);
	val += lookup_flagged_intensity(intensities, flags, -i, -h,  l);
	val += lookup_flagged_intensity(intensities, flags, -k, -i,  l);
	val += lookup_flagged_intensity(intensities, flags,  k,  h, -l);
	val += lookup_flagged_intensity(intensities, flags,  h,  i, -l);
	val += lookup_flagged_intensity(intensities, flags,  i,  k, -l);
	val += lookup_flagged_intensity(intensities, flags, -k, -h, -l);
	val += lookup_flagged_intensity(intensities, flags, -h, -i, -l);
	val += lookup_flagged_intensity(intensities, flags, -i, -k, -l);
	val += lookup_flagged_intensity(intensities, flags, -h, -k, -l);
	val += lookup_flagged_intensity(intensities, flags, -i, -h, -l);
	val += lookup_flagged_intensity(intensities, flags, -k,  i, -l);
	val += lookup_flagged_intensity(intensities, flags,  h,  k, -l);
	val += lookup_flagged_intensity(intensities, flags,  i,  h, -l);
	val += lookup_flagged_intensity(intensities, flags,  k,  i, -l);
	val += lookup_flagged_intensity(intensities, flags, -k, -h,  l);
	val += lookup_flagged_intensity(intensities, flags, -h, -i,  l);
	val += lookup_flagged_intensity(intensities, flags, -i, -k,  l);
	val += lookup_flagged_intensity(intensities, flags,  k,  h,  l);
	val += lookup_flagged_intensity(intensities, flags,  h,  i,  l);
	val += lookup_flagged_intensity(intensities, flags,  i,  k,  l);
	#endif /* PG6MMM */

	return val;
	#endif /* FLAT_INTENSITIIES */
}


kernel void diffraction(global float *diff, global float *tt, float klow,
                        int w, float corner_x, float corner_y,
                        float res, float clen, float16 cell,
                        global float *intensities,
                        int min_fs, int min_ss, int sampling, local float *tmp,
                        float kstep,
                        read_only image2d_t func_a,
                        read_only image2d_t func_b,
                        read_only image2d_t func_c,
                        global float *flags,
                        float fsx, float fsy, float ssx, float ssy,
                        float divxlow, float divxstep, int divxsamp,
                        float divylow, float divystep, int divysamp)
{
	float ttv;
	float fs, ss;
	float f_lattice, I_lattice;
	float I_molecule;
	float4 q;
	float k = klow + kstep * get_local_id(2);
	float intensity;
	const int ls0 = get_local_size(0);
	const int ls1 = get_local_size(1);
	const int ls2 = get_local_size(2) / (divxsamp*divysamp);
	const int ls3 = divxsamp;
	const int ls4 = divysamp;
	const int li0 = get_local_id(0);
	const int li1 = get_local_id(1);
	const int li234 = get_local_id(2);
	const int li2 = li234 / (ls3*ls4);
	const int li234leftover = li234 % (ls3*ls4);
	const int li3 = li234leftover / ls4;
	const int li4 = li234leftover % ls4;
	const int ls = ls0 * ls1 * ls2 * ls3 * ls4;
	float xdiv = divxlow + divxstep*ls4;
	float ydiv = divylow + divystep*ls3;

	/* Calculate fractional coordinates in fs/ss */
	fs = convert_float(get_global_id(0)) / convert_float(sampling);
	ss = convert_float(get_global_id(1)) / convert_float(sampling);

	/* Get the scattering vector */
	q = get_q(fs, ss, res, clen, k, &ttv,
	          corner_x, corner_y, fsx, fsy, ssx, ssy, xdiv, ydiv);

	/* Calculate the diffraction */
	f_lattice = lattice_factor(cell, q, func_a, func_b, func_c);
	I_molecule = molecule_factor(intensities, flags, cell, q);
	I_lattice = pow(f_lattice, 2.0f);

	/* Write the value to local memory */
	intensity = I_molecule * I_lattice;
	tmp[li0 + ls0*li1 + ls0*ls1*li2 + ls0*ls1*ls2*li3
	        + ls0*ls1*ls2*ls3*li4] = intensity;

	/* Memory fence */
	barrier(CLK_LOCAL_MEM_FENCE);

	/* Leader thread sums the values */
	if ( li0 + li1 + li2 + li3 + li4 == 0 ) {

		int i;
		float sum = 0.0;
		float val;
		int idx;

		idx = convert_int_rtz(fs) + w*convert_int_rtz(ss);

		for ( i=0; i<ls; i++ ) sum += tmp[i];

		val = sum / convert_float(ls);
		diff[idx] = val;

		/* Leader thread also records the 2theta value */
		tt[idx] = ttv;

	}
}
