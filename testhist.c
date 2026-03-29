/* compile with:
   gcc -o testhist testhist.c -std=c11 -Wall -pedantic -lm
*/

#include<stdio.h>
#include<stdlib.h>
#include "array.h"
#include "math.h"

double get_mode_from_binned_data(double *x, int n)
{
	double min=x[0], max=x[0], ra, *endpts;
	int i, j, k, nbins_Sturges = ceil(log2((double)n)) + 1, *bins;
	
	MAKE_VECTOR(bins, nbins_Sturges);
	MAKE_VECTOR(endpts, nbins_Sturges + 1);
	
	for (i = 1; i < n; i++) {
		if (x[i] < min) 
			min = x[i];
		else {
			if (x[i] > max)
				max = x[i];
		}
	}
	ra = (max - min)/nbins_Sturges;
	
	printf("%d %f \n", nbins_Sturges, ra);

	endpts[0] = min;
	for (i = 1; i < nbins_Sturges + 1; i++) {
		endpts[i] = min + ra*i;
		bins[i-1] = 0;
	}
	
	for (i = 0; i < n; i++) {
		for (j = 0 ; x[i] > endpts[j + 1]; j++);
		bins[j]++;
	}
	
	for (i = 0; i < nbins_Sturges; i++) 
		printf("%d ", bins[i]);

	/* get the index of the bin with the most occupancy */
	
	k = 0;
	j = 0;
	for (i = 0; i < nbins_Sturges; i++) {
		if (bins[i] > k) {
			k = bins[i];
			j = i;
		}
	}
		
	ra = (endpts[j] + endpts[j + 1])/2;
	
	FREE_VECTOR(bins);
	FREE_VECTOR(endpts);
	
	return ra;
}

int main(void)
{
	int i;
	double *x;
	FILE *ffile;
	
ffile = fopen("temp.txt", "r");

	MAKE_VECTOR(x, 1000);
	for (i = 0; i < 1000; i++) 
		fscanf(ffile, "%lf ", &x[i]);
	fclose(ffile);
	
	printf("%f \n", get_mode_from_binned_data(x, 1000));

	FREE_VECTOR(x);
	
	return EXIT_SUCCESS;
}
