#include<stdio.h>
#include<stdlib.h>
#include "array.h"

#define CUBE(x) ((x) * SQ(x))

void get_voxel_coordinates(int i, int *temp, int nx, int ny, int nz)
{
	int ii = i;
	if (nz > 1) {
		temp[2] = ii / (nx * ny);
		ii %= nx * ny;
	}
	temp[1] = ii / nx;
	temp[0] = ii % nx;
	return;
}

int main(void) 
{
	int i, j, k, nx, ny, nz, inu = 40, *tmp;
	double *mu, ***Mu;
	char str[50];

	FILE *ffile;

	sprintf(str, "../../dfaden/bwrf%i_180x216x180.txt", inu);
	
	ffile = fopen(str, "r");
	ffile = fopen("../../dfaden/final_MYDATA/T1-weighted.txt", "r");

	fscanf(ffile, "%d %d %d\n", &nx, &ny, &nz);

	MAKE_VECTOR(mu, nx * ny * nz);	

	MAKE_3ARRAY(Mu, nx, ny, nz); 

	for (i = 0; i < nx * ny * nz; i++) fscanf(ffile, "%lf ", &mu[i]);
	
	fclose(ffile);

	MAKE_VECTOR(tmp, 3);
	for (i = 0; i < nx * ny *nz; i++) {
		get_voxel_coordinates(i, tmp, nx, ny, nz);
		Mu[tmp[0]][tmp[1]][tmp[2]] = mu[i];
	}
	FREE_VECTOR(tmp);
	ffile = fopen("tmp.txt", "w");

	k = 90;
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			fprintf(ffile, "%f ", Mu[i][j][k]);
			}
	}
	
	fclose(ffile);
	
	FREE_VECTOR(mu);
	FREE_3ARRAY(Mu);

	return EXIT_SUCCESS;
}


