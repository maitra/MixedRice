#include "bessel_ratio.h"

#if USE_GSL

#include <gsl/gsl_sf_bessel.h>

double besselI1_I0(double x)
{
  return gsl_sf_bessel_I1_scaled(x) / gsl_sf_bessel_I0_scaled(x);
}

#else

/**
 * Returns an approximation to I(x, 1) / I(x, 0)
 * where I(x, v) is the modified Bessel function
 * of the first kind with order v.
 * See p. 378 Abramowitz.
 *
 * Assumes x >= 0.
 */
double besselI1_I0(double x)
{ 
  if (x > 3.75) {
    /*
     * Based off of 9.8.2 on p. 378 of Abramowitz and Stegun,
     * currently available at http://www.math.sfu.ca/~cbm/aands/page_378.htm
     * If less precision is required, we can simply throw out
     * some of the v terms, starting with v8 and working down.
     */
    const double tRecip = 3.75 / x;

    const double v8 = 0.00392377 * tRecip;
    const double v7 = (v8 - 0.01647633) * tRecip;
    const double v6 = (0.02635537 + v7) * tRecip;
    const double v5 = (-0.02057706 + v6) * tRecip;
    const double v4 = (0.00916281 + v5) * tRecip;
    const double v3 = (-0.00157565 + v4) * tRecip;
    const double v2 = (0.00225319 + v3) * tRecip;
    const double v1 = 0.39894228 + (0.01328592 + v2) * tRecip;

    /**
     * Based off of 9.8.4 on p. 378 of Abramowitz and Stegun
     */
    const double w8 = -0.00420059 * tRecip;
    const double w7 = (w8 + 0.01787654) * tRecip;
    const double w6 = (w7 - 0.02895312) *tRecip;
    const double w5 = (w6 + 0.02282967) * tRecip;
    const double w4 = (w5 - 0.01031555) * tRecip;
    const double w3 = (w4 + 0.00163801) * tRecip;
    const double w2 = (w3 - 0.00362018) * tRecip;
    const double w1 = 0.39894228 + (w2 - 0.03988024) * tRecip;

    return w1 / v1;
  }
  else if (x == 0.0) {
    return 0.0;
  }
  else {
    /*
     * Based off of 9.8.1 on p. 378 of Abramowitz and Stegun.
     */
    const double t = x/3.75;
    const double tSq = t * t;
    const double v12 = 0.0045813 * tSq;
    const double v10 = (0.0360768 + v12) * tSq;
    const double v8 = (0.2659732 + v10) * tSq;
    const double v6 = (1.2067492 + v8) * tSq;
    const double v4 = (3.0899424 + v6) * tSq;
    const double v2 = 1.0 + (3.5156229 + v4) * tSq;

    /**
     * Based off of 9.8.3 on p. 378 of Abramowitz and Stegun.
     */
    const double w12 = 0.00032411 * tSq;
    const double w10 = (w12 + 0.00301532) * tSq;
    const double w8 = (w10 + 0.02658733) * tSq;
    const double w6 = (w8 + 0.15084934) * tSq;
    const double w4 = (w6 + 0.51498869) * tSq;
    const double w2 = 0.5 + (w4 + 0.87890594) * tSq;
    
    return x * (w2 / v2);
  }

}

#endif
