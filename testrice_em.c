#include<stdio.h>
#include<stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#define SQ(x) ((x) * (x))
#define MAX(a, b) ((a) < (b) ? b : a)

int run_rice_em_sample(int nx, int ny, int nz, int gapgrid, int minclus, 
		       int maxclus, double *X, double *BIC,
		       double *ICL,  double *Sigma, double *SD_sigma, 
		       double *SD_q, double *pval, int *K);

int get_estimates(int totalKs, double *Sigma, double *BIC, double *ICL, 
		  double *SD_sigma, double *SD_q, double *ests,
		  double *pvalues, double alpha, int *K);

int main(void) 
{
	int i, ifl, nx, ny, nz, gridgap = 4, minclus = 1, totalKs, *Ks, 
		maxclus = 50;
	double *bic, *sigma, *icl, *sd_sigma, *sd_q, *X, *ests, *pvals;
	unsigned int seed1, seed2;
	char str[50];
	FILE *fout, *ffile;

	ffile = fopen("random.seed", "r");
	fscanf(ffile, "%u %u", &seed1, &seed2);
	fclose(ffile);

	set_seed(seed1, seed2);
	
	fout = fopen("phantom.results", "w");


	MAKE_VECTOR(bic, maxclus - minclus + 1);
	MAKE_VECTOR(icl, maxclus - minclus + 1);
	MAKE_VECTOR(sigma, maxclus - minclus + 1);
	MAKE_VECTOR(sd_sigma, maxclus - minclus + 1);
	MAKE_VECTOR(sd_q, maxclus - minclus + 1);
	MAKE_VECTOR(pvals, maxclus - minclus + 1);
	MAKE_VECTOR(Ks, maxclus - minclus + 1);
	MAKE_VECTOR(ests, 7);

        /* allocate space for X */

	ffile = fopen("../../dfaden/final_phantom/phantom/phantom_1.txt", "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);
	fclose(ffile);

	MAKE_VECTOR(X, nx * ny * nz);

	for (ifl = 1; ifl <= 12; ifl++) {
		sprintf(str,"../../dfaden/final_phantom/phantom/phantom_%i.txt",ifl);
		ffile = fopen(str, "r");
		fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);
		for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &X[i]);
		fclose(ffile);
		
		totalKs = run_rice_em_sample(nx, ny, nz, gridgap, minclus, 
					     maxclus, X, bic, icl, sigma,
					     sd_sigma, sd_q, pvals, Ks);
		
		ests[6] = get_estimates(totalKs, sigma, bic, icl, sd_sigma, 
					sd_q, ests, pvals, 0.05, Ks);
		for (i = 0; i < 7; i++) fprintf(fout, "%f ",ests[i]);
		fprintf(fout, "\n");
		
	}

	FREE_VECTOR(X);
	FREE_VECTOR(pvals);
	FREE_VECTOR(Ks);
	FREE_VECTOR(bic);
	FREE_VECTOR(icl);
	FREE_VECTOR(sigma);
	FREE_VECTOR(sd_sigma);
	FREE_VECTOR(sd_q);
	FREE_VECTOR(ests);

	
	fclose(fout);

	get_seed(&seed1, &seed2);
	
	ffile = fopen("random.seed", "w");
	fprintf(ffile, "%u %u\n", seed1, seed2);
	fclose(ffile);
	
	return EXIT_SUCCESS;
}
	
/*
true       BIC      ICL-BIC  CV
1.190374  1.292863 1.117615 1.072686
1.426326  1.230822 1.230822 1.202513
1.002233  0.927902 1.110726 1.273852
0.696797  0.719855 0.426593 0.817575
0.8185849 0.834701 0.814483 1.291181
0.5980409 0.629673 0.836614 0.836614
1.254058  1.409012 1.244281 1.244281
1.316242  1.533286 1.604344 1.390389
0.9265976 1.136157 1.103105 1.451422
0.6965587 0.881380 0.703792 1.554320
0.900108  1.040337 0.830191 1.393987
0.69413   0.734279 0.734279 1.199623

*/
