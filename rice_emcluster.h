#ifndef __RICE_EMCLUSTER_H__
#define __RICE_EMCLUSTER_H__

int *rice_emcluster(int n, int k, double *pi, double *X, double *Mu, 
		    double *Sigma, int maxiter, double eps, double *llhdval,
		    double *ICL);

void rice_estep(int n, int k, double *X, double **Gamma, double **RGamm,
		double *pi, double *Mu, double Sigma);

#endif
