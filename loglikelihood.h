#ifndef __LOGLIKELIHOOD_H__
#define  __LOGLIKELIHOOD_H__

/**
 * Estimated complete data log-likelihood.
 */
double riceLogLikelihood(const double* y, int numPoints,
			 int numClusters, const double* mu, const double* logP,
			 double sigmaSq, const double* const* z);

/**
 * Computes the log-likelihood of the specified mixture of Rician
 * components at the given data points y.
 */
double observedDataLogLikelihood(const double* y, const int numPoints,
				 const double* p, const double* mu,
				 const int numClusters,
				 const double sigma);

/**
 * Computes the log-density of component Rice(mu, sigma) with mixture
 * proportion p at y.
 * Assumes y > 0.
 */
double logDensityComponent(double y, double p, double mu, double sigma);

#endif

