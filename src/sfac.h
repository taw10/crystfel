/*
 * sfac.h
 *
 * Scattering factors
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SFAC_H
#define SFAC_H

#include <complex.h>

#include "cell.h"


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

	/* Unit cell */
	UnitCell *cell;

	/* Reflection intensities at Bragg positions */
	double complex *reflections;

	/* Offset to molecule's centre of scattering power */
	double xc;
	double yc;
	double zc;
};


extern double complex get_sfac(const char *n, double s, double en);
extern struct molecule *load_molecule(void);
extern double complex *get_reflections(struct molecule *mol, double en);
extern void get_reflections_cached(struct molecule *mol, double en);

#endif	/* SFAC_H */
