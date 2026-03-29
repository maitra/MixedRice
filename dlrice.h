#include "logBessel.h"

#ifndef __DLRICE_H__
#define __DLRICE_H__

/**
 * Returns the log-density of the distribution Rice(v, a)
 * at x. For x > 0, this is the logarithm of 
 x/v^2 * exp(-(x^2+a^2)/(2*v^2)) * besselI(x*a/v^2, 0)
 *
 * For a = 0. this is the log-density of the Rayleigh distribution with
 Rayleigh(v). For x > 0, this is the logarithm of 
 x / v^2 * exp(-x^2/(2*v^2))

*/

double dlrice(double x, double v, double a);

#endif

