#ifndef __RICE_EM_VARIANCE_H__
#define __RICE_EM_VARIANCE_H__

double delqdeltheta(int i, int l, int K, double *X, double **Gamma, 
		    double **RGamma, double *Mu, double Sigma, double *pi);

double VarSigma(int n, int k, double *X, double **Gamma, double **RGamma,
		double *Mu, double Sigma, double *pi, double **sqrtmat);

double VarQ(int n, int K, double *X, double **Gamma, double **RGamma,
	    double *Mu, double Sigma, double *pi, double *varsigma);  

#endif
