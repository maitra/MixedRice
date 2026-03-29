#include "logBessel.h"
#include <math.h>

double logBesselI0(double x) {
  if (x >= 3.75) {
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
    const double v1 = (0.01328592 + v2) * tRecip;
    /*
     * Is sqrt cheaper than log?
     * -0.5 log(x) + log(c) = log(c/sqrt(x)) -- inconclusive so far
     */ 
    return x - 0.5 * log(x) + log(v1 + 0.39894228);
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
    const double v2 = (3.5156229 + v4) * tSq;
    return log(1.0 + v2);
  }
}
