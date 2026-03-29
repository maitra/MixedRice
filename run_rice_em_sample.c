#include<stdio.h>
#include<stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#include "loglikelihood.h"
#include "order.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#define SQ(x) ((x) * (x)) 
#define Inf 1e+140
#define MIN(a, b) ((a) < (b) ? a : b)
#define MAX(a, b) ((a) < (b) ? b : a)

double gofsignif(double *x, int n, int nclus, double sigma, double *mu, 
		 double *pi, double *teststat);

double *bh_fdr(int n, double *p)
{
	int i;
	size_t *orderedp = orderDouble(p, n);
	double *q;

	MAKE_VECTOR(q, n);

	for (i = 0; i < n; i++) q[i] = p[i] * n / (orderedp[i] + 1.0); 

	free(orderedp);

	return q;

}


double rrice(double mu, double sigma)
{
	double x = rnorm(mu, sigma), y = rnorm(0, sigma);
	return sqrt(SQ(x) + SQ(y));
}

void rricemix(double *y, int n, int k, double *mu, double *pi, double sigma)
{
	int i, j, *nclass, l = 0;
	
	MAKE_VECTOR(nclass, k);

	rmultinom(n, pi, k, nclass);
	/* note that rmultinom only returns the number of observations in each
	   class */

	for (i = 0; i < k; i++) {
		for (j = 0; j < nclass[i]; j++) y[l++] = rrice(mu[i], sigma);
	}
	FREE_VECTOR(nclass);

	return;
}

double expected_llhd(int nvoxels, int nclus, double *pi, double *mu, 
		     double sigma, int numbest)
/* function to estimate the expected value of loglikelihood of theta by simple
   parametric bootstrap */
{
	double *x, llh, sum = 0.0;
	int i;
	
	MAKE_VECTOR(x, nvoxels);

	for (i = 0; i < numbest; i++) {
		rricemix(x, nvoxels, nclus, mu, pi, sigma);
		llh = observedDataLogLikelihood(x, nvoxels, pi, mu, nclus, 
						sigma); 
		sum += llh;
	}
	FREE_VECTOR(x);
	return sum;
}

int run_rice_em_sample(int nx, int ny, int nz, int gapgrid, int minclus, 
		       int maxclus, double *X, double *BIC, double *ICL,  
		       double *Sigma, double *SD_sigma, double *SD_q, 
		       double *pval, int *K) 
{
	int i, nclus, *iclass, j, l, ri, rj, rl, inum = 0, nvoxels, nvox;
	double sigma, *mu, *pi, llhval, curricl, currbic, icl, 
		teststat, currllhval = -Inf, var, varq, *x;

	nvox = ceil(nx / (double)gapgrid) * ceil(ny / (double) gapgrid) * ceil(nz / (double)gapgrid);

	MAKE_VECTOR(x, nvox);

/*	get the sample */

	ri = MIN(gapgrid, nx) * runif(0, 1);
	rj = MIN(gapgrid, ny) * runif(0, 1);
	rl = MIN(gapgrid, nz) * runif(0, 1); 

//	printf("nvoxels = %i\n", nvox);
	
	nvoxels = 0;
	for (i = 0; (gapgrid * i + ri) < nx; i++) {
		for (j = 0; (gapgrid * j + rj) < ny; j++) {
			for (l = 0; (gapgrid * l + rl) < nz; l++) {
				x[nvoxels++] = X[(gapgrid * l + rl)* nx * ny + (gapgrid * i + ri) * ny + (gapgrid * j + rj)];
			}
		}
	}	
	
//	printf("new-nvoxels = %i ", nvoxels);


	for (nclus = minclus; nclus <= maxclus; nclus++) {
		MAKE_VECTOR(mu, nclus);
		MAKE_VECTOR(pi, nclus);
		rice_em_rndinit(nvoxels, nclus, x, pi, mu, &sigma, 
				500 + 50 * nclus);
		iclass = rice_emcluster(nvoxels, nclus, pi, x, mu, &sigma, 
					1000, 0.0001, &llhval, &icl);
		llhval = observedDataLogLikelihood(x, nvoxels, pi, mu,
						   nclus, sigma);
		FREE_VECTOR(iclass);
		
		if (nclus == 1) {
			double sum = 0.;
			for (i = 0; i < nvoxels; i++) sum += SQ((-2 + SQ(x[i]/sigma))/sigma);
			var = 1/sum;
			sum = 0.;
			for (i = 0; i < nvoxels; i++) sum -= SQ(x[i]);
			sum *= 3/SQ(SQ(sigma));
			sum += 2 * nvoxels / SQ(sigma);
//			printf("sum = %f\n", sum);
			varq = var * SQ(sum) / 2;
		}
		else {
			double **Gamma, **RGamma;
			
			MAKE_MATRIX(Gamma, nvoxels, nclus);
			MAKE_MATRIX(RGamma, nvoxels, nclus);
			
			rice_estep(nvoxels, nclus, x, Gamma, RGamma, pi, mu, 
				   sigma);
			
			varq = VarQ(nvoxels, nclus, x, Gamma, RGamma, mu, 
				   sigma, pi, &var);

			FREE_MATRIX(Gamma);	
			FREE_MATRIX(RGamma);
		}

		if (llhval > currllhval) {
			pval[inum] = gofsignif(x, nvoxels, nclus, 
					       sigma, mu, pi, &teststat);
			currbic = -2 * llhval + (2 * nclus - 1) * 
				log((double)nvoxels);
			curricl = -2 * icl + (2 * nclus - 1) * 
				log((double)nvoxels);
			currllhval = llhval;
		
			//			printf("K = %i sigma = %2.3f llhval = %f BIC = %f ICL-BIC = %f SD = %f SDq = %f ICL = %f, pval = %0.3f \n", nclus, sigma, llhval, currbic, curricl, sqrt(var), sqrt(varq), icl, pval[inum]);
			
			BIC[inum] = currbic;
			ICL[inum] = curricl;
			SD_sigma[inum] = sqrt(var);
			SD_q[inum] = sqrt(varq);
			Sigma[inum] = sigma;
			K[inum] = nclus;
			inum++;
		}
		FREE_VECTOR(pi);
		FREE_VECTOR(mu);
	}
		
	FREE_VECTOR(x);
	
	return inum;
}

double remainingaveragedSD(int i, int n, double *x)
{
	int j;
	double sum = 0;
	
	for (j = (i + 1); j <= (i + 2); j++) sum += (i + 2 - j) * SQ(x[j]);

	return sqrt(sum / 2);
} 



int get_estimates(int totalKs, double *Sigma, double *BIC, double *ICL, 
		   double *SD_sigma, double *SD_q, double *ests,
		   double *pvalues, double alpha, int *K)
{
	/* totalKs is the number of cases for which we hope that we have got 
	   valid optimized loglikelihoods (means that we have nondecreasing 
	   loglikelihoods) */
	
	int i;
	double curricl, currbic, maxdip_q = -Inf, maxdip_s = -Inf, dip, dipq;
	double *qvalues;

	/* get the BIC and ICL estimates*/
	
	for (i = 0; i < 6; i++) ests[i] = Sigma[0];
	curricl = ICL[0];
	currbic = BIC[0];

	for (i = 0; i < totalKs; i++) {
		if (BIC[i] < currbic) {
			ests[0] = Sigma[i];
			currbic = BIC[i];
		}		
		if (ICL[i] < curricl) {
			ests[1] = Sigma[i];
			curricl = ICL[i];
		}
	}
	
	for (i = 0; i < (totalKs - 2); i++) {
		dip = remainingaveragedSD(i, totalKs, SD_sigma) - SD_sigma[i];///SD_sigma[i];
		dipq = remainingaveragedSD(i, totalKs, SD_q) - SD_q[i];///SD_q[i];
		if (dip > maxdip_s) {
			ests[2] = Sigma[i];
			maxdip_s = dip;
		}
		if (dipq > maxdip_q) {
			ests[3] = Sigma[i];
			maxdip_q = dipq;
		}
	}

	qvalues = bh_fdr(totalKs, pvalues);

	for (i = 0; ((i < totalKs)  && (pvalues[i] < alpha)); i++);


	if (i < totalKs) ests[4] = Sigma[i];
	else ests[4] = Sigma[totalKs - 1];

	for (i = 0; ((i < totalKs)  && (qvalues[i] < alpha)); i++);

	if (i < totalKs) ests[5] = Sigma[i];
	else ests[5] = Sigma[totalKs - 1];

	FREE_VECTOR(qvalues);

	if (i < totalKs) return K[i];
	else return K[totalKs - 1];
}


