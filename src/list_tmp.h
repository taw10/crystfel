static inline void LABEL(integrate)(TYPE *ref, signed int h,
                                    signed int k, signed int l,
                                    TYPE i)
{
	int idx;

	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		printf("\nReflection %i %i %i is out of range!\n", h, k, l);
		printf("You need to re-configure INDMAX, delete the reflection"
		       " cache file and re-run.\n");
		exit(1);
	}

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	ref[idx] += i;
}


static inline void LABEL(set)(TYPE *ref, signed int h,
                              signed int k, signed int l,
                              TYPE i)
{
	int idx;

	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		printf("\nReflection %i %i %i is out of range!\n", h, k, l);
		printf("You need to re-configure INDMAX, delete the reflection"
		       " cache file and re-run.\n");
		exit(1);
	}

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	ref[idx] = i;
}


static inline TYPE LABEL(lookup)(TYPE *ref, signed int h,
                                 signed int k, signed int l)
{
	int idx;

	if ( (abs(h) > INDMAX) || (abs(k) > INDMAX) || (abs(l) > INDMAX) ) {
		printf("\nReflection %i %i %i is out of range!\n", h, k, l);
		printf("You need to re-configure INDMAX, delete the reflection"
		       " cache file and re-run.\n");
		exit(1);
	}

	if ( h < 0 ) h += IDIM;
	if ( k < 0 ) k += IDIM;
	if ( l < 0 ) l += IDIM;

	idx = h + (IDIM*k) + (IDIM*IDIM*l);
	return ref[idx];
}


static inline TYPE *LABEL(new_list)(void)
{
	TYPE *r;
	size_t r_size;
	r_size = IDIM*IDIM*IDIM*sizeof(TYPE);
	r = malloc(r_size);
	memset(r, 0, r_size);
	return r;
}
