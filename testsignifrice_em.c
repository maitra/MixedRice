#include<stdio.h>
#include<stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>
#define SQ(x) ((x) * (x))

double **run_rice_em_plus(int nvoxels, int minclus, int maxclus, double *x,
 			  double *llhval, double *icl, double *var, double *varq);

int main(void) 
{
	int i, j, l, nvoxels = 0, nx, ny, nz, minclus = 10, maxclus = 30;
	double *x, sigmab, sigmai, sigmac, sigmaq, *X, **parmests,
		*llh, *icl, *var, *varq;
	unsigned int seed1, seed2;
	FILE *ffile;

	ffile = fopen("random.seed", "r");
	fscanf(ffile, "%u %u", &seed1, &seed2);
	fclose(ffile);

	set_seed(seed1, seed2);
 	
	ffile = fopen("../../dfaden/final_phantom/phantom/phantom_3.txt", "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);

	MAKE_VECTOR(X, nx * ny * nz);

	for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &X[i]);
	fclose(ffile);

	MAKE_VECTOR(x, nx * ny * nz / 64);

	for (i = 0; i < nx/8; i++) {
		for (j = 0; j < ny/8.; j++) {
			for (l = 0; l < nz/8.; l++) {
				x[nvoxels++] = X[8 * i * nx + 8 * j];
			}
		}
		}	
/*
	MAKE_VECTOR(x, nx * ny * nz / 16);

	for (i = 0; i < nx/4; i++) {
		for (j = 0; j < ny/4.; j++) {
			for (l = 0; l < nz/4.; l++) {
				x[nvoxels++] = X[4 * i * nx + 4 * j];
			}
		}
		}	
*/
/*	printf("%d\n", nvoxels);*/

	FREE_VECTOR(X);

	MAKE_VECTOR(llh, maxclus - minclus + 1);
	MAKE_VECTOR(icl, maxclus - minclus + 1);
	MAKE_VECTOR(var, maxclus - minclus + 1);
	MAKE_VECTOR(varq, maxclus - minclus + 1);

	parmests = run_rice_em_plus(nvoxels, minclus, maxclus, x, llh, icl, 
				    var, varq);

//	printf("BIC-estimated sigma = %f ICL-estimated sigma = %f Minimum CV-estimated sigma = %f\n", sigmab, sigmai, sigmac);
//	printf("%f %f %f %f\n", sigmab, sigmai, sigmac, sigmaq);	
	
	for (i = minclus; i < (maxclus - 1); i++) {
		
	}

	FREE_VECTOR(varq);
	FREE_VECTOR(var);
	FREE_VECTOR(icl);
	FREE_VECTOR(llh);
	

	FREE_VECTOR(x);

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
