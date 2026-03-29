#include<stdio.h>
#include<stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#define SQ(x) ((x) * (x))
#define CUBE(x) ((x) * SQ(x))

int run_rice_em_sample(int nx, int ny, int nz, int gapgrid, int minclus, 
		       int maxclus, double *X, double *BIC,
		       double *ICL,  double *Sigma, double *SD_sigma, 
		       double *SD_q, double *pval, int *K);

int get_estimates(int totalKs, double *Sigma, double *BIC, double *ICL, 
		  double *SD_sigma, double *SD_q, double *ests,
		  double *pvalues, double alpha, int *K);

double rrice(double mu, double sigma);

int main(void) 
{
	int i, j, nx, ny, nz, minclus = 1, maxclus = 15, *Ks, totalKs,
		gridgap = 8, inu = 40;
	double *mu, *bic, *sigma, *icl, *sd_sigma, *sd_q, *X, *ests, *pvals, 
	  trsig = 50;
	char str[50], outstr[20];

	unsigned int seed1, seed2;
	FILE *ffile, *fout;

	ffile = fopen("random.seed", "r");
	fscanf(ffile, "%u %u", &seed1, &seed2);
	fclose(ffile);

	set_seed(seed1, seed2);

	sprintf(str, "../../dfaden/bwrf%i_180x216x180.txt", inu);
	
	ffile = fopen(str, "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);

	MAKE_VECTOR(X, nx * ny * nz);
	MAKE_VECTOR(mu, nx *ny *nz);	

	for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &mu[i]);
	fclose(ffile);

	sprintf(outstr, "result_%i_%i.txt", inu, (short)trsig); 

	fout = fopen(outstr, "w");
	
	MAKE_VECTOR(bic, maxclus - minclus);
	MAKE_VECTOR(icl, maxclus - minclus);
	MAKE_VECTOR(sigma, maxclus - minclus);
	MAKE_VECTOR(sd_sigma, maxclus - minclus);
	MAKE_VECTOR(sd_q, maxclus - minclus);
	MAKE_VECTOR(pvals, maxclus - minclus + 1);
	MAKE_VECTOR(Ks, maxclus - minclus + 1);
	
	MAKE_VECTOR(ests, 7);

	for(i = 0; i < 50; i++) {
		for (j = 0; j < nx * ny * nz; j++) X[j] = rrice(mu[j], trsig); 

		totalKs = run_rice_em_sample(nx, ny, nz, gridgap, minclus, 
					     maxclus, X, bic, icl, sigma, 
					     sd_sigma, sd_q, pvals, Ks);

		ests[6] = get_estimates(totalKs, sigma, bic,icl, sd_sigma, 
					sd_q, ests, pvals, 0.05, Ks);
	
//		printf("rep = %i", i);
//		for (j = 0; j < 7; j++) printf("%f ",ests[j]);
		for (j = 0; j < 7; j++) fprintf(fout, "%f ",ests[j]);
//		printf("\n");
		fprintf(fout, "\n");
	}

	FREE_VECTOR(pvals);
	FREE_VECTOR(Ks);
	FREE_VECTOR(bic);
	FREE_VECTOR(icl);
	FREE_VECTOR(sigma);
	FREE_VECTOR(sd_sigma);
	FREE_VECTOR(sd_q);
	FREE_VECTOR(ests);
	
	FREE_VECTOR(mu);
	FREE_VECTOR(X);


	fclose(fout);
	
	get_seed(&seed1, &seed2);
	
	ffile = fopen("random.seed", "w");
	fprintf(ffile, "%u %u\n", seed1, seed2);
	fclose(ffile);

	return EXIT_SUCCESS;
}


