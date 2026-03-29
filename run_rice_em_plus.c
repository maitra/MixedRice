#include<stdio.h>
#include<stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#define SQ(x) ((x) * (x)) 
#define Inf 1e+140

double **run_rice_em_plus(int nvoxels, int minclus, int maxclus, double *x,
 			  double *llhval, double *icl, double *var, 
			  double *varq)
/* similar to run_rice_em, but returns parameter estimates in a two-dimensional
   ragged array. In addition, arrays var, varq, llhval, icl return values of
   the variance estimate of sigma and varq as well as llhval and icl. */ 
{
	int i, nclus, *iclass, *noparameters;
	double sigma, *mu, *pi, **parests;

	MAKE_VECTOR(noparameters, maxclus - minclus + 1);
	for (nclus = minclus; nclus <= maxclus; nclus++) 
		noparameters[nclus - minclus] = 2 * nclus - 1;
	MAKE_2RAGGEDARRAY(parests, maxclus - minclus + 1, noparameters);

	for (nclus = minclus; nclus <= maxclus; nclus++) {

		MAKE_VECTOR(mu, nclus);
		MAKE_VECTOR(pi, nclus);
		
		rice_em_rndinit(nvoxels, nclus, x, pi, mu, &sigma, 200);
//				3 * nclus + 50);
		iclass = rice_emcluster(nvoxels, nclus, pi, x, mu, &sigma, 
					1000, 0.0001, &llhval[nclus - minclus], &icl[nclus - minclus]);
		FREE_VECTOR(iclass);

		if (nclus == 1) {
			double sum = 0.;
			for (i = 0; i < nvoxels; i++) sum += SQ((-2 + SQ(x[i]/sigma))/sigma);
			var[nclus - minclus] = 1/sum;
			sum = 0.;
			for (i = 0; i < nvoxels; i++) sum -= SQ(x[i]);
			sum *= 3/SQ(SQ(sigma));
			sum += 2 * nvoxels / SQ(sigma);
//			printf("sum = %f\n", sum);
			varq[nclus - minclus] = var * SQ(sum) / 2;
		}
		else {
			double **Gamma, **RGamma;

			MAKE_MATRIX(Gamma, nvoxels, nclus);
			MAKE_MATRIX(RGamma, nvoxels, nclus);
						
			rice_estep(nvoxels, nclus, x, Gamma, RGamma, pi, mu, 
				   sigma);
			varq[nclus - minclus] = VarQ(nvoxels, nclus, x, Gamma, 
						     RGamma, mu, sigma, pi, 
						     &(var[nclus - minclus]));
			
			FREE_MATRIX(Gamma);	
			FREE_MATRIX(RGamma);
		}

		m = 0;
		for (j = 0; j < nclus - 1; j++) 
			parests[nclus - minclus][m++] = mu[j];
		for (j = 0; j < nclus - 1; j++)
			parests[nclus - minclus][m++] = pi[j];
		parests[nclus-miclus][m] = sigma;
		
		FREE_VECTOR(pi);
		FREE_VECTOR(mu);
	}
		
	return parests;
}

double rrice(double mu, double sigma)
{
	double x = rnorm(mu, sigma), y = rnorm(0, sigma);
	return sqrt(SQ(x) + SQ(y));
}

double *rricemix(int n, int k, double *mu, double *pi, double sigma)
{
	int i, j, *nclass, l = 0;
	double *y;
	
	MAKE_VECTOR(nclass, k);

	rmultinom(n, pi, k, nclass);
	/* note that rmultinom only returns the number of observations in each
	   class */

	MAKE_VECTOR(y, n);
	for (i = 0; i < k; i++) {
		for (j = 0; j < nclass[i]; j++) y[l++] = rrice(mu[i], sigma);
	}
	FREE_VECTOR(nclass);

	return y;
}

double expected_llhd(int nvoxels, int nclus, double *pi, double *mu, 
		     double sigma, int numbest)
/* function to estimate the expected value of loglikelihood of theta by simple
   parametric bootstrap */
{
	double *x, llh, sum = 0.0;
	int i, j;
	
	for (i = 0; i < numbest; i++) {
		x = rricemix(nvoxels, nclus, mu, pi, sigma);
		llh = observedDataLogLikelihood(x, n, pi, mu, nclus, (*sigma));
		sum += llh;
	}

	FREE_VECTOR(x);
	return sum;
}

double generate_rice_data(int nvoxels, int nclus, int altnclus, double *x, 
			  double *pi, double *mu, double sigma, int numbest)
{
	double *x, llh, icl, var, varq, **Gamma, **RGamma, *pitry, *mutry;
	int i, j;
	
	x = rricemix(nvoxels, nclus, mu, pi, sigma);

	MAKE_VECTOR(pitry, nclus);
	MAKE_VECTOR(mutry, nclus);

	rice_emcluster(nvoxels, k, pitry, x, mutry, &sigma, 1000, 1e-04, &llh, 
		       &icl);

	MAKE_MATRIX(Gamma, nvoxels, nclus);
	MAKE_MATRIX(RGamma, nvoxels, nclus);
	
	rice_estep(nvoxels, nclus, x, Gamma, RGamma, pitry, mutry, sigma);
	varq = VarQ(nvoxels, nclus, x, Gamma, RGamma, mu, sigma, pi, &var);


	rice_em_rndinit(n, altclus, x, pi, mu, sigma, numbest);


	

	FREE_MATRIX(Gamma);
	FREE_MATRIX(RGamma);

	FREE_VECTOR(x);
}

