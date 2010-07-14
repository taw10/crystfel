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
#include "utils.h"


#define MAX_ATOMS (128*1024)


struct mol_species
{
	char species[4];    /* Species name */
	int n_atoms;        /* Number of atoms of this species */

	double x[MAX_ATOMS];
	double y[MAX_ATOMS];
	double z[MAX_ATOMS];
	double occ[MAX_ATOMS];
	double B[MAX_ATOMS];
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


/* This is so that the water background calculation can use it */
extern double complex get_sfac(const char *n, double s, double en);

extern struct molecule *load_molecule(const char *filename);
extern void free_molecule(struct molecule *mol);
extern double *get_reflections(struct molecule *mol, double en, double res,
                               double *phases, ReflItemList *items);


#endif	/* SFAC_H */
