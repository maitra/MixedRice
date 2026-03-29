#include "stdlib.h"
#include "rice_em_rndinit.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#include "srswor.h"
#include "array.h"
#include "rice_emcluster.h"
#define SQ(x) ((x) * (x))
#define Inf 1e+140

double getdistWSSassignclass(int n, double *X, double newMu, int *currclass, 
			     double *currdist, int k)
{
	/* X - vector of n observations
	   newMu - new center from which distances are to be calculated
	   currclass - current class to which observation is closest 
	               (input, updated, if reassigned)
	   currdist - vector of squared distances of each observation to 
	              closest mean (input, updated if reassigned) 
	   k - group id to which newMu belongs
	 */

	int i;
	double dissq, sum = 0;

	for (i = 0; i < n; i++) {
		dissq = SQ(X[i] - newMu);
		if (dissq < currdist[i]) {
			currdist[i] = dissq;
			currclass[i] = k;
		}
		sum += currdist[i];
	}
	return sum;
}

void random_inits(int n, int k, double* X, double* pi, double* Mu, 
		  double* Sigma)
{
	double *W, *prob, rnf;
	int i, j, *class, *nk;

	MAKE_VECTOR(class, n);		    
	MAKE_VECTOR(prob, n);
	MAKE_VECTOR(W, n);   /*vector of current squared distances to closest
			       Mu*/
	Mu[k-1] = 0.;
	
	(*Sigma) = 0.;
	for (i = 0; i < n; i++) { /* get an initial estimate of sigma */
		class[i] = k - 1;
		W[i] = SQ(X[i]);
		(*Sigma) += W[i];
	}
	(*Sigma) /= 2*n; /* Note that this is actually Sigma-squared. Also, 
			    the division by 2 is because of the assumption that
			    we have observations in the background from the 
			    Rayleigh distribution. */

	for (j = 0; j < (k - 1); j++) {
		for (i = 0; i < n; i++) prob[i] = 1 - exp(-W[i]/(2*(*Sigma)));
                /* these are the probabilities of inclusion of each X[i]*/
		for (i = 1; i < n; i++)  prob[i] += prob[i-1]; 

		for (i = 0; i < n; i++) prob[i] /= prob[n-1];

		/* this is the cdf */
		rnf = unif_rand();
		for (i = 0; (i < n) && (prob[i] < rnf); i++);

		Mu[j] = X[i];
		
		(*Sigma) = getdistWSSassignclass(n, X, Mu[j], class, W, j);
		(*Sigma) /= 2*n;
	}
	 
	(*Sigma) = (sqrt(*Sigma)); /* get to the correct scale */
 
	FREE_VECTOR(W);
	FREE_VECTOR(prob);
	
	MAKE_VECTOR(nk, k);
	for (i = 0; i < k; i++) nk[i] = 0;
	for (i = 0; i < n; i++) nk[class[i]]++;
	for (i = 0; i < k; i++) pi[i] = nk[i]/(double) n;
	FREE_VECTOR(nk);

	FREE_VECTOR(class);
	return;
}

void rice_em_rndinit(int n, int k, double *X, double *pi, double *Mu, 
		     double* Sigma, int numbest)
{
	int i;

	(*Sigma) = 0;
	for (i = 0; i < n; i++) (*Sigma) += SQ(X[i]);
	(*Sigma) /= 2*n;

	(*Sigma) = (sqrt(*Sigma));
          /* a possibility is that we take just some initial value of sigma 
	     which is between 0 and the correct estimate for 
	     Rayleigh-distributed data, randomly chosen */
/*	printf("k = %i Sigma = %f runif(0,1) = %f ", k, *Sigma, rnf);*/

	if (k == 1) {
		/* this is the case where we only have Rayleigh-distributed 
		   background noise, nothing else. */
		Mu[0] = 0.;
		pi[0] = 1.;
	}
	else {
		int ii, *iclass;
		double llhval, oldlval = -Inf, Sigmadum, *pidum, icl,
			*Mudum;

//	   	(*Sigma) *= rnf;
//		Sigmadum = (*Sigma);

		MAKE_VECTOR(pidum, k);
		MAKE_VECTOR(Mudum, k);
		
		for (ii = 0; ii < numbest; ii++) {

			random_inits(n, k, X, pidum, Mudum, &Sigmadum);

//			Mudum = sample(X, n, k - 1);
//			Mudum = realloc(Mudum, k * sizeof(double));
			/*In general, not a good idea to directly enlarge Mudum,
			  but probably ok since we are only really talking 
			  about adding one extra cell */ 
//			Mudum[k-1] = 0.; /*allocate new value */

//			for (i = 0; i < k; i++) pidum[i] = 1./k;

			iclass = rice_emcluster(n, k, pidum, X, Mudum, 
						&Sigmadum, 5, 0.01, &llhval, 
						&icl);

			FREE_VECTOR(iclass);

			if (llhval > oldlval) {
				for (i = 0; i < k; i++) {
					Mu[i] = Mudum[i];
					pi[i] = pidum[i];
				}
				(*Sigma) = Sigmadum;
				oldlval = llhval;
			}
			//	free(Mudum);
			llhval = oldlval;
/*			printf("k = %i Sigma = %f  llhval = %f\n", k, *Sigma, llhval);*/
		}
		FREE_VECTOR(pidum);
		FREE_VECTOR(Mudum);
	}
	return;	
}
