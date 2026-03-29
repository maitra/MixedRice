#include <stdio.h>
#include <stdlib.h>
#include "rice_emcluster.h"
#include "rice_em_rndinit.h"
#include "rice_em_variance.h"
#include "array.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>
#define SQ(x) ((x) * (x))
#define CUBE(x) ((x) * SQ(x))

void overall_local_estimates(int nx, int ny, int nz, double *x, double ests[4],
			     int z_plane);

double rrice(double mu, double sigma);

int main(int argc, char *argv[]) 
{
	int i, j, nx, ny, nz, nsample = 20, inu, zplane = 1;
	double *X, *mu, est[4], trsig;
	char str[50], outstr[20], *end;
	unsigned int seed1, seed2;
	FILE *ffile, *fout;

	if (argc != 5) {
		fprintf(stderr, "Incorrect usage!\nUsage is as follows:\n  ./testevaluate_local_sigma [inu] [true sigma] [random seed 1] [random seed 2]\n  where \"inu\" is one of 0, 10 and 20, and the seed 1 and seed 2 are long ints.\nNote that both seeds are between 0 and 2^31-1.\t Example Usage:\n\t ./testevaluate_local_sigma 0 10 2088683337 1317056914\n");
		exit(1);
	}

	inu = atoi(argv[1]);
	trsig = atof(argv[2]);
	seed1 = strtoul(argv[3], &end, 0);
	seed2 = strtoul(argv[4], &end, 0);

	printf("%d %f %d %d \n", inu, trsig, seed1, seed2);

	set_seed(seed1, seed2);

	sprintf(outstr, "start-local_%i_%i.txt", inu, (short)trsig); 

	ffile = fopen(outstr, "w");
	fprintf(ffile, "%u %u\n", seed1, seed2);
	fclose(ffile);

	sprintf(str, "../../dfaden/bwrf%i_180x216x180.txt", inu);
	
	ffile = fopen(str, "r");
	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);

	MAKE_VECTOR(X, nx * ny * nz);
	MAKE_VECTOR(mu, nx * ny * nz);	

	for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &mu[i]);
	fclose(ffile);

	sprintf(outstr, "result-local_%i_%i.txt", inu, (short)trsig); 

	fout = fopen(outstr, "w");
		
	for(i = 0; i < nsample; i++) {
		for (j = 0; j < nx * ny * nz; j++) X[j] = rrice(mu[j], trsig); 

		overall_local_estimates(nx, ny, nz, X, est, zplane);
			
//		printf("rep = %i", i);
//		for (j = 0; j < 7; j++) printf("%f ",ests[j]);
		for (j = 0; j < 4; j++) fprintf(fout, "%f ",est[j]);
//		printf("\n");
		fprintf(fout, "\n");
		fflush(fout);
	}
	
	FREE_VECTOR(mu);
	FREE_VECTOR(X);

	fclose(fout);
	
	return EXIT_SUCCESS;
}
