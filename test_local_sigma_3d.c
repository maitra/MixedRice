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
			     int z_plane);

int main(int argc, char *argv[]) 
{
        int i, nx, ny, nz, zplane = 0;
	double *X, est[4], y;
	unsigned int seed1, seed2;
	FILE *ffile, *fout;

	if (argc != 4) {
		fprintf(stderr, "Usage:\n\t ./test_local_sigma_3d [input_file] [output_file] [zplane]\nwhere zplane is 0 if only 2d slices are considered and 1 if 3d slices are considered\n");
		exit(1);
	}

	ffile = fopen("random.seed", "r");
	fscanf(ffile, "%u %u", &seed1, &seed2);
	fclose(ffile);

	set_seed(seed1, seed2);


        /* allocate space for X */

	/*	ffile = fopen("../../dfaden/final_breast/breast.txt", "r");
		zplane = 1;
		ffile = fopen("../../dfaden/final_MYDATA/T1-weighted.txt", "r");*/
	zplane = atoi(argv[3]);

	ffile = fopen(argv[1], "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);

	MAKE_VECTOR(X, nx * ny * nz);

	for (i = 0; i < nx * ny * nz; i++) 
	  { 
	    fscanf(ffile, "%lf ", &y);
	    if (y == 0) X[i] = runif(0, 0.01); 
	    else X[i] = y;
	  }

	fclose(ffile);

	overall_local_estimates(nx, ny, nz, X, est, zplane);
	
	fout = fopen(argv[2], "w");
	for (i = 0; i < 4; i++) fprintf(fout, "%f ",est[i]);
	fclose(fout);

	FREE_VECTOR(X);
	
	get_seed(&seed1, &seed2);
	
	ffile = fopen("random.seed", "w");
	fprintf(ffile, "%u %u\n", seed1, seed2);
	fclose(ffile);
	
	return EXIT_SUCCESS;
}
	
/*
./test_local_sigma_3d ../../dfaden/final_MYDATA/T1-weighted.txt T1-local-sds.txt 0
./test_local_sigma_3d ../../dfaden/final_MYDATA/T2-weighted.txt T2-local-sds.txt 0
./test_local_sigma_3d ../../dfaden/final_MYDATA/rho-weighted.txt rho-local-sds.txt 0
./test_local_sigma_3d ../../dfaden/final_breast/breast.txt breast-local-sds.txt 1
*/
