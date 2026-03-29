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

void overall_local_estimates(int nx, int ny, int nz, double *x, double ests[4],
	int plane_z);

int main(void) 
{
	int i, ifl, nx, ny, nz;
	double *X, est[4];
	char str[50];
	FILE *fout, *ffile;

        /* allocate space for X */

	ffile = fopen("../../dfaden/final_phantom/phantom/phantom_1.txt", "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);
	fclose(ffile);

	MAKE_VECTOR(X, nx * ny * nz);

	fout = fopen("phantom-local.results", "w");
	for (ifl = 1; ifl <= 12; ifl++) {
		sprintf(str,"../../dfaden/final_phantom/phantom/phantom_%i.txt",ifl);
		ffile = fopen(str, "r");
		fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);
		for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &X[i]);
		fclose(ffile);

		overall_local_estimates(nx, ny, nz, X, est, 0);

		for (i = 0; i < 4; i++) printf("%f ",est[i]);
		printf("\n");
		
		for (i = 0; i < 4; i++) fprintf(fout, "%f ",est[i]);
		fprintf(fout, "\n");
		
	}

	FREE_VECTOR(X);

//	fclose(fout);

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
