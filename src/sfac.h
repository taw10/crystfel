/*
 * sfac.h
 *
 * Scattering factors
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SFAC_H
#define SFAC_H

#include <complex.h>


struct mol_species
{
	char species[4];    /* Species name */
	int n_atoms;        /* Number of atoms of this species */

	double x[32*1024];
	double y[32*1024];
	double z[32*1024];
	double occ[32*1024];
	double B[32*1024];
};


struct molecule
{
	int n_species;
	struct mol_species *species[32];

	double xc;
	double yc;
	double zc;
};


extern double complex get_sfac(const char *n, double s, double en);
extern struct molecule *load_molecule(void);

#endif	/* SFAC_H */
