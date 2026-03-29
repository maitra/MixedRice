#ifndef __BESSEL_RATIO_H__
#define __BESSEL_RATIO_H__

/**
 * Returns an approximation to I(x, 1) / I(x, 0)
 * where I(x, v) is the modified Bessel function
 * of the first kind with order v.
 * See p. 378 Abramowitz.
 *
 * Assumes x >= 0.
 */
double besselI1_I0(double x);

#endif

