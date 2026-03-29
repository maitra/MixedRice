#include<stdio.h>
#include<math.h>
#include "loglikelihood.h"
#include "dlrice.h"
#include "array.h"
#include "matvecops.h"
#include "bessel_ratio.h"
#include "rice_emcluster.h"

#define PI 3.141593
#define Inf 1e+140
#define SQ(x) ((x) * (x))

void rice_estep(int n, int k, double *X, double **Gamma, double **RGamm,
		double *pi, double *Mu, double Sigma)
{
  /* This is the E-Step: assigns the responsibilities of each of the k groups 
     to each of the n observations. 
     The inputs are:
     n     = number of Ricean observations
     k     = number of components (all Rice; Rayleigh incorporated within, in 
             the last k'th component if present)
     X     = vector of Ricean observations
     pi    = vector of prior probabilities
     Mu    = vector of Ricean-distributed means (last zero, if Rayleigh 
             component is present)     
     Sigma = (common) noise parameter of Rice-distributed components
     Gamma = n x k matrix of posterior probability of observation X[i] 
             belonging to the kth group. 
     RGamm = n x k matrix of posterior expectation of the product of the
             indicator (that the ith observation is in kth group) and 
	     the cosine of angle given ith magnitude observation belongs to 
	     the kth group
  */

	int i, l;
	double *temp, sum, *minmax;
	
	MAKE_VECTOR(temp, k);
	for (i = 0; i < n; i++) {
		sum = 0.;
		for (l = 0; l < k; l++) 
			temp[l] = log(pi[l]) + dlrice(X[i], Sigma, Mu[l]);
		minmax = range(temp, k); /* get the minimum and maximum values*/
		for (l = 0; l < k; l++) {
			temp[l] -= minmax[1]; /*reduce to get smaller values*/
			sum += exp(temp[l]);
		}				
		FREE_VECTOR(minmax);
		
		for (l = 0; l < k; l++) {
			Gamma[i][l] = exp(temp[l])/sum;
			RGamm[i][l] = Gamma[i][l] * besselI1_I0( X[i] * Mu[l] / SQ(Sigma));
		}
	}
	FREE_VECTOR(temp);
	return;
}

void rice_mstep(double *X, int n, int k, double *pi, double *Mu, 
		double *Sigma, double **Gamma, double **RGamm)  
{
	int i, ll, iRayleigh = 0;
	double *sum, sums = 0;
	
	MAKE_VECTOR(sum, k);

	if (Mu[k - 1] == 0) iRayleigh = 1; /*this means the Rayleigh component*/
		


	for (ll = 0; ll < k; ll++) {
		sum[ll]=0.;
		Mu[ll] = 0;
		for (i = 0; i < n; i++) {
			sum[ll] += Gamma[i][ll];
			Mu[ll] += X[i] * RGamm[i][ll];
		}
		Mu[ll] /= sum[ll];
		pi[ll] = sum[ll]/n;
	}
	FREE_VECTOR(sum);

	if (iRayleigh) Mu[k - 1] = 0;
	
	for (i = 0; i < n; i ++) {
		sums += SQ(X[i]);
		for (ll = 0; ll < k; ll++) {
			sums -= 2 * Mu[ll] * X[i] * RGamm[i][ll];
			sums += Gamma[i][ll] * SQ(Mu[ll]);
		}
	}
	sums /= 2 * n;

	(*Sigma) = sqrt(sums);
	
/*	printf("sums = %f Sigma = %f\n", sums, *Sigma);*/
	
	return;
}

int *classify(int n, int k, double **Gamma)
{
	int i, j, *clas;
	double rowmax = -Inf;
	
	MAKE_VECTOR(clas, n);

	for (i = 0; i < n; i++) {
		rowmax = -Inf;
		for (j = 0; j < k; j++) {
			if (Gamma[i][j] > rowmax) {
				clas[i] = j;
				rowmax = Gamma[i][j];
			}
		}
	}
	return clas;
}



double icl(int n, int k, double *X, int *class, double *Mu, double *pi, 
	   double sigma)
{
	int i;
	double sum = 0.;
	for (i = 0; i < n; i++) {
		sum += log(pi[class[i]]) - 2*log(sigma) 
			- (SQ(X[i]) - 2 * X[i] * Mu[class[i]] + SQ(Mu[class[i]]) )/(2*SQ(sigma));
	}
	return sum;
}


int *rice_emcluster(int n, int k, double *pi, double *X, double *Mu, 
		    double *Sigma, int maxiter, double eps, double *llhdval,
		    double *ICL)
{
	int iter = 0, *class;
	double **gamm, **rgamm, oldllhd;

	MAKE_MATRIX(gamm, n, k);
	MAKE_MATRIX(rgamm, n, k);
	
	*llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, (*Sigma));

	do  {
		rice_estep(n, k, X, gamm, rgamm, pi, Mu, (*Sigma));
		rice_mstep(X, n, k, pi, Mu, &(*Sigma), gamm, rgamm);
		oldllhd = *llhdval;
		*llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, (*Sigma));
		iter++;
	} 
	while (((*llhdval - oldllhd) > eps*fabs(oldllhd)) && (iter < maxiter));

	*llhdval = observedDataLogLikelihood(X, n, pi, Mu, k, (*Sigma));
	class = classify(n, k, gamm);
	*ICL = icl(n, k, X, class, Mu, pi, (*Sigma));

	FREE_MATRIX(rgamm);
	FREE_MATRIX(gamm);
	
	return class;
}
