#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "dlrice.h"
#include "nelder_mead_min.h"

#define SQ(x) ((x) * (x)) 
#define CUBE(x) ((x) * SQ(x))
#define Inf 1e+140
#define MIN(a, b) ((a) < (b) ? a : b)
#define MAX(a, b) ((a) < (b) ? b : a)

typedef struct realdata {
	int n;
	int p;
	double *x;
} RealData;

int findMode(int *array, int elems)
{
	int i, lastnum = array[0], curlen = 1, mode = array[0], modelen = 1;
	for (i = 1; i < elems; i++)	   {
		if (array[i] == lastnum)
			curlen++;
		else	{
			if (curlen > modelen)	{
				modelen = curlen;
				mode = lastnum;
			}
			lastnum = array[i];
			curlen = 1;
		}
	}
	return mode;
}



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

void getnhbrs(double ***X, int *vox, int win_length, double *x, int nz)
{
	int i, j, k, l = 0;
	
	for (i = -1 * win_length; i < win_length; i++) {
		for (j = -1 * win_length; j < win_length; j++) {
			if (nz > win_length) {
				for (k = -1 * win_length; k < win_length; k++) {
					x[l++] = X[vox[0] + i][vox[1] + j][vox[2] + k];
				}
			}
			else x[l++] = X[vox[0] + i][vox[1] + j][vox[2]];
		}
	}

	return;
}

void calculate3moments(double *x, double *mean, double *var, double *m3, int n)
{
	double sum1 = 0, sum2 = 0, sum3 = 0, tmp; 
	int i;
	
	for (i = 0; i < n; i++) 
		sum1 += x[i];
	(*mean) = sum1/n;

	for (i = 0; i < n; i++) {
		tmp = (x[i] - (*mean));
		sum2 += SQ(tmp);
		sum3 += CUBE(tmp);
	}
	*var = sum2/(n - 1);
	*m3 = sum3/n;

//	printf("%f %f %f\n", *mean, *var, *m3);

	return;
}

double negllhd(int p, const double *theta, const void *ext)
{
	const RealData *xx = ext;
	int i;
	double sum = 0.;

//	printf("%f %f %i\n", theta[0], theta[1], xx -> n);

/*	for (i = 0; i < xx -> n; i++) 
		printf("%f ", xx -> x[i]);
		printf("\n");*/

	for (i = 0; i < xx -> n; i++) {
		if (xx -> x[i] > 0)
			sum -= dlrice(xx -> x[i], theta[1], theta[0]);
	}
	return sum;
}



double mle(double *x, double mu, double sigma, int n)
{
	const double abstol = 1e-20, reltol = 1e-20;
	const double alpha  = 1.0, beta   = 0.5, gamma  = 2.0;
	const int trc    = 0;   

	int p = 2, fail, fncnt, maxit=10000;
	double *init, *final, fmin, sigmahat;
	RealData *Y;
	
	Y = malloc(sizeof(RealData));
	Y -> n = n;
	Y -> p = 2;
	MAKE_VECTOR(Y -> x, Y -> n);
	MAKE_VECTOR(init, p);
	MAKE_VECTOR(final, p);

	for (p = 0; p < Y -> n; p++)
		Y -> x[p] = x[p];
	
	init[0] = mu;
	init[1] = sigma;
	
//	printf("%f %f\n", init[0], init[1]);

	nelder_mead_min(Y -> p, init, final, &fmin, negllhd, &fail, abstol, 
			reltol, Y, alpha, beta, gamma, trc, &fncnt, maxit);
	
	sigmahat = final[1];

//	printf("%f %f %d %f\n", init[0], init[1], fail, fmin);

/*	for (p = 0; p < 2; p++)
		printf("%f %f ", init[p], final[p]);
	printf("\n");

*/
	FREE_VECTOR(init);
	FREE_VECTOR(final);
	FREE_VECTOR(Y -> x);
	free(Y);

//	printf("%f ", sigmahat);
	
	return sigmahat;

}


double integerpower(double x, int n) 
{
	double tmp = 1;
	int i;
	
	for (i = 0; i < n; i++) 
		tmp *= x;
	return tmp;
}

double polyexp(double x) 
{
	double cr = 0, p[10] = {1.0007570413, 2.8981188340, -72.9432278777,  
				   1162.6792136360, -9838.85598962208, 
				   47813.9607638493, -137448.5785417688, 
				   230670.4056296062, -208666.38136498138, 
				   78562.5551923769};
	int i;

	for (i = 0; i < 10; i++) 
		cr += p[i] * integerpower(x, i);
	return cr;
}

int rounded(double number)
{
	return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

double local_skewness_estimate(double var, double m3)
{
	double skewness = m3 / (var * sqrt(var));
	double cr = polyexp(skewness);

//	printf("skewness = %f, cr = %f, m3 = %f, sigma = %f ", skewness, cr, m3, sqrt(cr * var));
	
	return sqrt(cr * var);
}

double calculate_range(double *x, int n)
{
	double min=x[0], max=x[0];
	int i;

	for (i = 1; i < n; i++) {
		if (x[i] < min) 
			min = x[i];
		else {
			if (x[i] > max)
				max = x[i];
		}
	}
	return max - min;
}


void overall_local_estimates(int nx, int ny, int nz, double *x, double ests[2])
{

	int winlength, n, nlocests, *locests, *posests, *tmp, i, ll = 0, l = 0,
		*mlests, *posmlests, nmlests;
	double *y, ***X, mean, var, m3, ra = 256/calculate_range(x, nx * ny * nz);
	
	if ((nx >=256) && (ny >= 256))
		winlength = 4;
	else 
		winlength = 3;

	if (nz >= 2*winlength + 1) {
		n = CUBE(2*winlength + 1);
		nlocests = (nx - 2 * winlength) * (ny - 2 * winlength) * (nz - 2 * winlength);
	}
	else {
		n = SQ(2*winlength + 1);
		nlocests = (nx - 2 * winlength) * (ny - 2 * winlength);
	}
	MAKE_VECTOR(y, n);
	MAKE_3ARRAY(X, nx, ny, nz);
	MAKE_VECTOR(locests, nlocests);
	MAKE_VECTOR(mlests, nlocests);

	MAKE_VECTOR(tmp, 3);
	for (i = 0; i < nx * ny *nz; i++) {
		get_voxel_coordinates(i, tmp, nx, ny, nz);
		X[tmp[0]][tmp[1]][tmp[2]] = x[i];
	}

	if (nz > winlength) {
		for (tmp[2] = winlength; tmp[2] < nz - winlength; tmp[2]++) {
			for (tmp[1] = winlength; tmp[1] < ny - winlength; tmp[1]++) {
				for (tmp[0] = winlength; tmp[0] < nx - winlength; tmp[0]++) {
					getnhbrs(X, tmp, winlength, y, nz);
					calculate3moments(y, &mean, &var, &m3, n);
					locests[l++] = rounded(ra * local_skewness_estimate(var, m3));
//					printf("%f %f %d\n", mean, var, n);
					mlests[ll++] = rounded(ra * mle(y, mean, sqrt(var), n));
				}
			}
		}
	}
	else {
		tmp[2] = nz  - 1;
		for (tmp[1] = winlength; tmp[1] < ny - winlength; tmp[1]++) {
			for (tmp[0] = winlength; tmp[0] < nx - winlength; tmp[0]++) {
				getnhbrs(X, tmp, winlength, y, nz);
				calculate3moments(y, &mean, &var, &m3, n);
				locests[l++] = local_skewness_estimate(var, m3);
				mlests[ll++] = rounded(mle(y, mean, sqrt(var), n));
			}
		}
	}

	FREE_VECTOR(tmp);
	FREE_VECTOR(y);

	nlocests = l;
	nmlests = ll;
	MAKE_VECTOR(posests, nlocests);
	MAKE_VECTOR(posmlests, nmlests);

	/*remove zero local estimates */

	l = 0; 
	for (i = 0; i < nlocests; i++) {
		if (locests[i] > 0) 
			posests[l++] = locests[i];
	}
	ll = 0;
	for (i = 0; i < nmlests; i++) {
		if (mlests[i] > 0) 
			posmlests[ll++] = mlests[i];
	}
	
	FREE_3ARRAY(X);
	
	ests[0] = findMode(posests, l) / ra;
	ests[1] = findMode(posmlests, ll) / ra;

	FREE_VECTOR(locests);
	FREE_VECTOR(posests);
	FREE_VECTOR(mlests);
	FREE_VECTOR(posmlests);

	return;
}
