/*
 * utils.h
 *
 * Utility stuff
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifndef UTILS_H
#define UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

extern unsigned int biggest(signed int a, signed int b);
extern unsigned int smallest(signed int a, signed int b);
extern double distance(double x1, double y1, double x2, double y2);
extern double modulus(double x, double y, double z);
extern double modulus_squared(double x, double y, double z);
extern double angle_between(double x1, double y1, double z1,
                            double x2, double y2, double z2);
extern double angle_between_d(double x1, double y1, double z1,
                              double x2, double y2, double z2);
extern double lambda(double voltage);
extern double distance3d(double x1, double y1, double z1,
                         double x2, double y2, double z2);
extern size_t skipspace(const char *s);
extern void chomp(char *s);
extern int sign(double a);

#define rad2deg(a) ((a)*180/M_PI)
#define deg2rad(a) ((a)*M_PI/180)

#define is_odd(a) ((a)%2==1)

#endif	/* UTILS_H */
