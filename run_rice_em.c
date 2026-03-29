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

void run_rice_em(int nvoxels, int minclus, int maxclus, double *x, 
		 double *BIC_sigma, double *ICL_sigma, double *CV_sigma,
		 double *Varq_sigma)
{
	int i, nclus, *iclass;
	double sigma, *mu, *pi, llhval, curricl, currbic, currcv, bic = Inf, 
		icl, iclbic = Inf, cv = Inf, var, varq, currvarq = Inf,
		currllhval = -Inf;

	for (nclus = minclus; nclus <= maxclus; nclus++) {

		MAKE_VECTOR(mu, nclus);
		MAKE_VECTOR(pi, nclus);
		
		rice_em_rndinit(nvoxels, nclus, x, pi, mu, &sigma, 200);
//				3 * nclus + 50);
		iclass = rice_emcluster(nvoxels, nclus, pi, x, mu, &sigma, 
					1000, 0.0001, &llhval, &icl);
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

		currcv = sqrt(var);
		currbic = -2 * llhval + (2 * nclus - 1) * log((double)nvoxels);
		curricl = -2 * icl + (2 * nclus - 1) * log((double)nvoxels);

//		printf("k=%i %f %f %f %e %e %e %e \n", nclus, sigma, currbic, curricl, currcv, varq, llhval, icl);

		if (currllhval < llhval) {
			currllhval = llhval;
			if (currcv < cv) {
				cv = currcv;
				(*CV_sigma) = sigma;
			}
			if (varq < currvarq) {
				currvarq = varq;
				(*Varq_sigma) = sigma;
			}
			if (currbic < bic) {
				bic = currbic;
				(*BIC_sigma) = sigma;
			}
			if (curricl < iclbic) {
				iclbic = curricl;
				(*ICL_sigma) = sigma;
			}
		}
		FREE_VECTOR(pi);
		FREE_VECTOR(mu);

	}
//	printf("%f %f %f %f\n", *CV_sigma, *Varq_sigma, *BIC_sigma, *ICL_sigma);
		
	return;
}
