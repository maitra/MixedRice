#include <assert.h>
#include <stdlib.h>
#include "sorted.h"

#ifndef USING_RLIB
#define MATHLIB_STANDALONE 1 /*It is essential to have this before the call 
                               to the Rmath's header file because this decides
                               the definitions to be set. */
#endif
#include <Rmath.h>
#define SQ(x) ((x) *(x))
#define MIN(a, b) ((a) < (b) ? a : b)
#define MAX(a, b) ((a) < (b) ? b : a)

double AndersonDarlingTest(int n, double *x, double (*teststat));

//double K(int n,double d) ;                        //Kolmogorov distribution

/* This function calculates the p-value of the Kolmogorov test statistic
   obtained by fitting the Rician mixture models at the MLE as defined 
   in the Marsaglia, et al (2003, JSS) paper. The code provided in that paper
   is also used in the calculation (function K(samplesize, test statistic)
   and is in the function Kolmogorov.c.
   
   The basic idea of this test is as follows. Since the observations are from
   a Rice-Rayleigh mixture with same noise parameter, the squared observations
   are from a mixture of scaled chi-squared distributions with different 
   non-centrality parameters and one central chi-squared distribution. These
   chi-squared distributions have 2 degrees of freedom.

   Note that if X \sim R(\sigma, \mu), then X^2 \sim \sigma^2 \chisq(2, \mu^2).
   
   So F(x) = P[X <= x] = P[X^2 <= x^2] = \sum_{i = 1}^K pi_k * F(x^2/\sigma^2; df = 2, ncp = \mu_k^2/\sigma^2) + pi_{K+1} * F(x^2/sigma^2; df = 2, ncp = 0)
*/

double prayleighricemix(double x, int nclus, double sigma, double *mu, 
			double *pi)
{
	int k;
	double sum = 0;
/* Usage:
   
   central: double pchisq(double x, double df, int lower_tail, int log_p)
   noncentral: double pnchisq(double x, double df, double ncp, int lower_tail, 
                              int log_p) 
*/

	for (k = 0; k < nclus - 1; k++) 
		sum += pi[k] * pnchisq(SQ(x / sigma), 2, SQ(mu[k] / sigma), 1,
				       0);
	sum += pi[nclus-1] * pchisq(SQ(x / sigma), 2, 1, 0);

/*

	if (sum > 1) {
		double tmp = 0;
		printf("%e\n", sum - 1);
		for (k = 0; k < nclus;  k++) 
		{
			tmp += pi[k] * pnchisq(SQ(x / sigma), 2, SQ(mu[k] / sigma), 1, 0);
			printf("%i %i sigma = %f mu = %f pi = %f x = %f cdf = %e prod = %e tmp = %e\n",
			       nclus, k, sigma, mu[k], pi[k], x,
			       pnchisq(SQ(x / sigma), 2, SQ(mu[k] / sigma), 1,
				       0),
			       pi[k] * pnchisq(SQ(x / sigma), 2, SQ(mu[k] / sigma), 1,
					       0), tmp);


		}
		tmp += pi[nclus] * pchisq(SQ(x / sigma), 2, 1, 0);

		printf("%i sigma = %f pi = %f x = %f cdf = %e prod = %e tmp = %e\n",
		       nclus, sigma, pi[nclus], x,
		       pchisq(SQ(x / sigma), 2, 1, 0), pi[nclus] * pchisq(SQ(x / sigma), 2, 1, 0), tmp);
		printf("%e\n", sum);
	}
*/

	return MIN(1, MAX(sum, 0));
}

double findmax(int n, double *x)
{
	int i;
	double max = x[0];

	for (i = 1; i < n; i ++) {
		if (max < x[i]) max = x[i];
	}
	return max;	
}

double gofsignif(double *x, int n, int nclus, double sigma, double *mu, 
		 double *pi, double *teststat)
{
	double *u, pval;
	int i;

	MAKE_VECTOR(u, n);

	for (i = 0; i < n; i++) u[i] = prayleighricemix(x[i], nclus, sigma, 
							mu, pi);
	
	sort(n, u);

/*	for (i = 0; i < n; i++) printf("x[i] = %f F(x) = %f", x[i],
				       prayleighricemix(x[i], nclus, sigma, 
								mu, pi));

	for (i = 0; i < n; i++) printf("%f ", u[i]);
	printf("\n");*/
					       


/*
	MAKE_VECTOR(y, 2 * n);

	for (i = 0; i < n; i++) y[i] = u[i] - i / (double)n;
	for (i = 0; i < n; i++) y[i + n] = (i + 1) / (double) n - u[i];

	(*teststat) =  findmax(2 * n, y);

	FREE_VECTOR(y);

	printf("\n %d %f %f \n", n, *teststat, K(n, *teststat));

	return K(n, (*teststat)); */
	
	pval = 1 - AndersonDarlingTest(n, u, teststat);

	if (isnan(pval)) {
		for (i = 0; i < n; i++) printf("x[i] = %f F(x) = %f", x[i],
		prayleighricemix(x[i], nclus, sigma, 								mu, pi));
		for (i = 0; i < n; i++) printf("%f ", u[i]);
		printf("\n");
	}


	FREE_VECTOR(u);
	
	return pval;
}
