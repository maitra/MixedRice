#include "logBessel.h"
#include <math.h>

double riceLogLikelihood(const double* y, int numPoints,
			     int numClusters, const double* mu, const double* logP,
			     double sigmaSq, const double* const* z)
{
  double result = 0.0;
  int i;
  for (i = 0; i < numPoints; ++i) {
    const double* zi = z[i];
    const double yi = y[i];
    const double yv = yi / sigmaSq;
    const double logYs1 =  log(yi / sigmaSq);
    const double twoSigmaSq = 2.0 * sigmaSq;
    const double ys = (yi / twoSigmaSq) * yi;
    int k;

    for (k = 0; k < numClusters; ++k) {
      const double ak = mu[k];
      const double xk = yv * ak;
      result += zi[k] * (logP[k] + logYs1 + logBesselI0(xk) - ys - (ak / twoSigmaSq) * ak);
    }
  }
  return result;
}

double observedDataLogLikelihood(const double* y, const int numPoints,
				 const double* p, const double* mu,
				 const int numClusters,
				 const double sigma)
{
  int i;
  double result = 0.0;
  const double sigmaSq = sigma * sigma;
  const double twoSigmaSqR = 1.0 / (2 * sigmaSq);

  for (i = 0; i < numPoints; ++i) {
    int j;
    const double yi = y[i];
    const double yiScaled = yi / sigmaSq;
    const double logYiScaled = log(yiScaled);
    const double ys = yi * yi * twoSigmaSqR;
    double sum = 0.0;

    /*
     * XXX! How can we make this safer numerically?
     */

    for (j = 0; j < numClusters; ++j) {
      const double muj = mu[j];
      sum += p[j] * exp(logYiScaled + logBesselI0(muj * yiScaled)
			- ys - muj * twoSigmaSqR * muj);
    }
    result += log(sum);
  }
  return result;
}

/**
 * Computes the log-density of component Rice(mu, sigma) with mixture
 * proportion p at y.
 * Assumes y > 0.
 * not used at all
 */
double logDensityComponent(double *y, double p, double mu, double sigma)
{
  return observedDataLogLikelihood(y, 1, &p, &mu, 1, sigma);
}

