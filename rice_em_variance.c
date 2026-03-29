#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "inverse.h"
#include "eigens.h"
#include "rice_em_variance.h"
#define SQ(x) ((x) * (x))


double delqdeltheta(int i, int l, int K, double *X, double **Gamma, 
		    double **RGamma, double *Mu, double Sigma, double *pi)
{
	if (l < (K - 1)) 
		return (X[i] * RGamma[i][l] - Gamma[i][l] * Mu[l])/SQ(Sigma);
	/* this is the part for the derivative with respect to Mu */
	else {
		if (l < (2 * K - 2)) 
			return Gamma[i][l - (K - 1)]/pi[l - (K - 1)] - Gamma[i][K-1]/pi[K-1];
		/* partial derivative wrt pi_k*/
		else { /*partial derivative wrt Sigma */
			int j;
			double sum = 0.;
			
			for (j = 0; j < (K - 1); j++)
				sum += 2 * X[i] * RGamma[i][j] * Mu[j] - Gamma[i][j] * SQ(Mu[j]);
			return (-2 + (SQ(X[i]) - sum)/SQ(Sigma))/Sigma;
		}
	}
}

double VarSigma(int n, int k, double *X, double **Gamma, double **RGamma,
		double *Mu, double Sigma, double *pi, double **sqrtmat)  
{
	double sum, *varmat, *evec, *eval, cumsum = 0;
	int i, j, l, m = 0, p = 2 * k - 1;

	MAKE_VECTOR(varmat, p * (p + 1) / 2);

//	printf("p = %d ", p);

	for (l = 0; l < p; l++) {
		for (j = 0; j <= l; j++) {
			sum = 0;
			for (i = 0; i < n; i++) 
				sum += delqdeltheta(i, l, k, X, Gamma, RGamma, Mu, Sigma, pi) * delqdeltheta(i, j, k, X, Gamma, RGamma, Mu, Sigma, pi);
/*			if (j == l) printf("k = %d p = %d l = %d j= %d m = %d  sum = %e  \n", k, p, l, j, m, sum);*/
			varmat[m++] = sum;
//			printf("%f ", sum);
		}
	}

	MAKE_VECTOR(evec, p * p);
	MAKE_VECTOR(eval, p);
	eigend(varmat, evec, eval, p);

/*	printf("\n");
	for (i = 0; i < p; i++) printf("%f ", eval[i]);
	printf("\n");
	for (i = 0; i < p * p; i++) printf("%f ", evec[i]);
	printf("\n");*/

	for (i = 0; i < p; i++) if (eval[i] < 0) eval[i] = 0;
	
	sum = 0; 
	for (i = 0; i < p; i++) sum += eval[i];
	
	for (i = 0; (i < p) && (cumsum < 0.999*sum); i++) {
		cumsum += eval[i];
	}
	for (j = i; j < p; j++) eval[j] = 0.; 

	for (i = 0; i < p; i ++) {
		for (j = 0; j < p; j++) {
			sqrtmat[i][j] = 0;
			for (l = 0; l < p; l++) {
				if (eval[l] > 0) {
				sqrtmat[i][j] += sqrt(1/eval[l]) * evec[l * p + i] * evec[l * p + j];
				}
			}
		}
	}

	sum = 0.;
	for (l = 0; l < p; l++) {
		if (eval[l] > 0) {
			sum += (1/eval[l]) * evec[l * p + p - 1] * evec[l * p + p - 1];
		}
	}
	FREE_VECTOR(evec);
	FREE_VECTOR(eval);
	FREE_VECTOR(varmat);

	return sum;	
}

double** Hessian(int K, int n, double *X, double **Gamma, double **RGamma, 
		 double *Mu, double Sigma, double *pi)
{
	int i, l, j, p = 2 * K - 1;
	double **H;
	
	MAKE_MATRIX(H, p, p);
	
	for (i = 0; i < p; i++) {
		for (l = 0; l < p; l++) H[i][l] =0.;
	}

	for (i = 0; i < (K - 1); i++) {
		for (l = 0; l < n; l++) {
			H[i][i] -= Gamma[l][i];
			H[i][2 * K - 2] -= 2 * (RGamma[l][i] * X[l] - Gamma[l][i] * Mu[i]);
		}
		H[i][2 * K - 2] /= SQ(Sigma) * Sigma;
//		H[i][2 * K - 2] -= 2 * n / Sigma;
		H[2 * K - 2][i] = H[i][2 * K - 2];
		H[i][i] /= SQ(Sigma); 
	}               /* done with block for Mu's */
		
	for (i = (K - 1); i < (2 * K - 2); i++) {
		for (l = 0; l < n; l++) {
			H[i][i] -= Gamma[l][i - K + 1] / SQ(pi[i - K + 1])
				+ Gamma[l][K - 1] / SQ(pi[K - 1]);
			for (j = (K - 1); j < i; j++) 
				H[i][j] -= Gamma[l][K - 1] / SQ(pi[K - 1]);
		}
		for (j = (K - 1); j < i; j++) H[j][i] = H[i][j];
	}               /* done with block for pi's */
	
	for (l = 0; l < n; l++) {
		H[2 * K - 2][2 * K - 2] += SQ(X[l]);
		for (i = 0; i < (K - 1); i++) { 
			H[2 * K - 2][2 * K - 2] -= 2 * RGamma[l][i] * X[l] * Mu[i];
			H[2 * K - 2][2 * K - 2] += Gamma[l][i] * SQ(Mu[i]);

			if (Gamma[l][i] < RGamma[l][i]) printf("Yelp! %d %d %lf %lf", l, i, Gamma[l][i], RGamma[l][i]);
		} 
	}
	H[2 * K - 2][2 * K - 2] *= -3./SQ(SQ(Sigma));
	H[2 * K - 2][2 * K - 2] += 2 * n / SQ(Sigma);

	return H;
}

double VarQ(int n, int K, double *X, double **Gamma, double **RGamma,
	    double *Mu, double Sigma, double *pi, double *varsigma)  
{
	int i, j, p = 2 * K - 1;
	double **hess, **sqrtmat, **a, sum = 0;

	MAKE_MATRIX(sqrtmat, p, p);

	*varsigma = VarSigma(n, K, X, Gamma, RGamma, Mu, Sigma, pi, sqrtmat);
	
	hess = Hessian(K, n, X, Gamma, RGamma, Mu, Sigma, pi);

/*	for (i = 0; i < p; i++) {
		for (j = 0; j < p; j++) printf("%f ", hess[i][j]);
		printf("\n");
		} */

	MAKE_MATRIX(a, p, p);

	i = multiply(hess, p, p, sqrtmat, p, p, a);

/*
	for (i = 0; i < p; i++) {
		for (j = 0; j < p; j++) printf("%f ", sqrtmat[i][j]);
		printf("\n");
		}
*/
	i = multiply(sqrtmat, p, p, a, p, p, hess); /* no longer the hessian */
	
	FREE_MATRIX(sqrtmat);
	FREE_MATRIX(a);

	for (i = 0; i < p; i++) {
		for (j = 0; j < p; j++) sum +=SQ(hess[i][j]);
	}

	FREE_MATRIX(hess);
		
//	printf("sum = %f\n", sum);
				
	return 2 * sum / SQ(n);
}

