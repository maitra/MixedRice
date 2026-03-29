#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "dlrice.h"
#include "nelder_mead_min.h"
#define MATHLIB_STANDALONE
#include<Rmath.h>

#define SQ(x) ((x) * (x)) 
#define CUBE(x) ((x) * SQ(x))
#define Inf 1e+140
#define MIN(a, b) ((a) < (b) ? a : b)
#define MAX(a, b) ((a) < (b) ? b : a)


typedef struct brentdata {
	int n;
	double *x;
} BrentData;

int quantile(int n,double *x,double *p,double *q, int numqs);

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double *fmin, double tol);

typedef struct nelder_meaddata {
	int n;
	int p;
	double *x;
} Nelder_MeadData;

double AndersonDarlingTest(int n, double *x, double (*teststat));

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

void get_voxel_coordinates(int i, int temp[3], int nx, int ny, int nz)
{
	int ii = i;
	if (nz > 1) {
		temp[2] = ii / (nx * ny);
		ii %= nx * ny;
	}
	else 
		temp[2] = 0;
	temp[1] = ii / nx;
	temp[0] = ii % nx;
	return;
}

void getnhbrs(double ***X, int vox[3], int win_length, int winz_length,
	      double *x, int nz)
{
	int i, j, k, l = 0;

	for (i = -1 * win_length; i <= win_length; i++) {
		for (j = -1 * win_length; j <= win_length; j++) {
			if (nz > winz_length) {
				for (k = -1 * winz_length; k <= winz_length; k++)
					x[l++] = X[vox[0] + i][vox[1] + j][vox[2] + k];
			}
			else
				x[l++] = X[vox[0] + i][vox[1] + j][vox[2]];
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
	const Nelder_MeadData *xx = ext;
	int i;
	double sum = 0.;

	for (i = 0; i < xx -> n; i++) {
		if (xx -> x[i] > 0)
			sum -= dlrice(xx -> x[i], theta[1], theta[0]);
	}
	return sum;
}

double mle(double *x, double mu, double sigma, int n)
{
	const double abstol = 1e-20, reltol = 1e-20;
	const double alpha = 1.0, beta = 0.5, gamma = 2.0;
	const int trc = 0;   

	int p, fail, fncnt, maxit=10000;
	double *init, *final, fmin, sigmahat;
	Nelder_MeadData *Y;
	
	Y = malloc(sizeof(Nelder_MeadData));
	Y -> n = n;
	Y -> p = 2;
	MAKE_VECTOR(Y -> x, Y -> n);
	MAKE_VECTOR(init, Y -> p);
	MAKE_VECTOR(final, Y -> p);

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

	if (skewness > 1.5) return 0; 
	else 
//		printf("variance = %f, skewness = %f, cr = %f, m3 = %f, sigma = %f \n", var, skewness, cr, m3, sqrt(cr * var));
		return sqrt(cr * var);
}

double boxcox(double x, double lambda)
{
	if (SQ(lambda) < 1e-06) 
		return log(x);
	else 
		return (pow(x, lambda) - 1)/lambda;
}

double inverseboxcox(double x, double lambda)
{
	if (SQ(lambda) < 1e-06) 
		return exp(x);
	else 
		return pow(lambda * x + 1, 1/lambda);
}

double funcboxcox(double lambda, void *ext)
{
	BrentData *xx = ext;
	int i;
	double *y, *u, v, m, m3, teststat, pv;
  
	MAKE_VECTOR(y, xx -> n);
	MAKE_VECTOR(u, xx -> n);

	for (i = 0; i < xx -> n; i++) 
		y[i] = boxcox(xx -> x[i], lambda);
	calculate3moments(y, &m, &v, &m3, xx -> n);
	for (i = 0; i < xx -> n; i++)
			u[i] = pnorm(y[i], m, sqrt(v), 1, 0);

	pv = AndersonDarlingTest(xx -> n, u, &teststat);

	FREE_VECTOR(y);
	FREE_VECTOR(u);
	
	if (isfinite(teststat)) {
		if (pv >= 0.05) 
			return (pv - 0.05);
		else 
			return 1 + (0.05 - pv);
	}
	else 
		return Inf;
}

double findbestBoxCox(double *xx, double *var, int nn)
{
	double xmin, optlambda;
	double m, v, m3;
	int i;
	BrentData *Y;
	
	Y = malloc(sizeof(BrentData));
	Y -> n = nn;

	MAKE_VECTOR(Y -> x, Y -> n);
	for (i = 0; i < Y -> n; i++) 
		Y -> x[i] = xx[i];
	
	optlambda = Brent_fmin(-4.0, 1.0, funcboxcox, Y, &xmin, 0.01);
	/* no crying need to be exceptionally accurate so end @ 0.01 accuracy */

	FREE_VECTOR(Y -> x);
	free(Y);

//	printf("%f \n", optlambda);

	for (i = 0; i < nn; i++) 
		xx[i] = boxcox(xx[i], optlambda);
	calculate3moments(xx, &m, &v, &m3, nn);
	(*var) = v; 

	return optlambda;
}
	
double get_mode_from_binned_data(double *y, int nn)
{
	double min=Inf, max=-Inf, ra, *endpts, *x, var, lambda, *p, *q;
	int i, j, k = 0, nbins, *bins, n = 0;
	double h; /* for Scott's nbins formula */

	MAKE_VECTOR(x, nn);

	MAKE_VECTOR(p, 2);
	MAKE_VECTOR(q, 2);
	p[0] = 0.25;
	p[1] = 0.75;

	i = quantile(nn, y, p, q, 2);

//	printf("median = %f IQR = %f\n", q[1], q[2]-q[0]);

	for (i = 0; i < nn; i++) {
		if (y[i] <= (q[1] + 4.5 * (q[1]-q[0]))) 
			x[n++] = y[i];
	}

	FREE_VECTOR(q);
	FREE_VECTOR(p);

	lambda = findbestBoxCox(x, &var, n);
	h = 3.5 * sqrt(var) * pow((double)n, -1.0/3); /* for Scott's nbins formula */

	for (i = 0; i < n; i++) {
		if (x[i] < min) 
			min = x[i];
		else {
			if (x[i] > max)
				max = x[i];
		}
	}

	nbins = ceil((max - min)/h); /* Scott's formula for number of bins */

//	nbins = ceil(log2((double)n) + 1); /* Sturges formula for nbins*/
	
	MAKE_VECTOR(bins, nbins);
	MAKE_VECTOR(endpts, nbins + 1);

	ra = (max - min)/nbins;

//	printf("\n min = %f, max = %f, ra = %f nbins = %d, v = %f, h = %f\n", min, max, ra, nbins, var, h);
	
	endpts[0] = min;
	for (i = 1; i < nbins + 1; i++) {
		endpts[i] = min + ra*i;
		bins[i-1] = 0;
	}
	
	for (i = 0; i < n; i++) {
//		fprintf(stderr, "%f %f %f\n", x[i], max, endpts[nbins]);
		//for (j = 0 ; x[i] > endpts[j + 1]; j++);
		for (j = 1 ; j <= nbins ; j++) {
			if (x[i] < endpts[j])
				break;
		}
		if (j >= nbins)	// check!
			j--;
		bins[--j]++;
	}
	FREE_VECTOR(x);

	k = 0;
	j = 0;
	for (i = 0; i < nbins; i++) {
		if (bins[i] > k) {
			k = bins[i];
			j = i;
		}
	}

	ra = inverseboxcox((endpts[j] + endpts[j + 1])/2, lambda);

	FREE_VECTOR(bins);
	FREE_VECTOR(endpts);
	
	return ra;
}


void overall_local_estimates(int nx, int ny, int nz, double *x, double ests[4],
			     int plane_z)
{

	int winlength, winzlength = 0, n, nlocests, *locests, *posests, tmp[3], 
		i, ll = 0, l = 0, *mlests, *posmlests, nmlests, myl = 0, 
		myll = 0;
	double *y, ***X, mean, var, m3, *mylocests, *myposests, *mymlests, 
		*myposmlests, temp;

	if ((nz > 1) && (plane_z)) /* slice thickness < 2mm */  {
		if ((nx >=256) && (ny >= 256)) 
			winlength = 3;
		else 
			winlength = 2;
		winzlength = 1;
	}
	else /* only 2d-slices to be calculated on */ {
		if ((nx >=256) && (ny >= 256)) 
			winlength = 4;
		else 
			winlength = 3;
	}
	
	n = SQ(2*winlength + 1) * (2*winzlength + 1);

	/*	nlocests = (nx - 2*winlength) * (ny - 2*winlength) * (nz - 2*winzlength); // I think this is correct but let us allocate more to be on the safe side */
	nlocests = nx * ny * nz;

	MAKE_VECTOR(y, n);
	MAKE_3ARRAY(X, nx, ny, nz);
	MAKE_VECTOR(locests, nlocests);
	MAKE_VECTOR(mlests, nlocests);
	MAKE_VECTOR(mylocests, nlocests);
	MAKE_VECTOR(mymlests, nlocests);

	for (i = 0; i < nx * ny * nz; i++) {
		get_voxel_coordinates(i, tmp, nx, ny, nz);
		X[tmp[0]][tmp[1]][tmp[2]] = x[i];
	}

	for (tmp[2] = winzlength; tmp[2] < nz - winzlength; tmp[2]++) {
			for (tmp[1] = winlength; tmp[1] < ny - winlength; tmp[1]++) {
				for (tmp[0] = winlength; tmp[0] < nx - winlength; tmp[0]++) {
					getnhbrs(X, tmp, winlength, winzlength, y, nz);
					calculate3moments(y, &mean, &var, &m3, n);
					if ((var > 0) && (m3 >= 0)) {
							temp = local_skewness_estimate(var, m3);
							mylocests[l] = temp;
							locests[l++] = rounded(temp);
							temp = mle(y, mean, sqrt(var), n);
							mymlests[ll] = temp;
							mlests[ll++] = rounded(temp);
					}
				}
			}
	}
	FREE_VECTOR(y);

//	printf("%d %d\n", l, ll);

	nlocests = l;
	nmlests = ll;
	MAKE_VECTOR(posests, nlocests);
	MAKE_VECTOR(posmlests, nmlests);
	MAKE_VECTOR(myposests, nlocests);
	MAKE_VECTOR(myposmlests, nmlests);

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

	for (i = 0; i < nlocests; i++) {
		if (mylocests[i] > 0) 
			myposests[myl++] = mylocests[i];
	}
	for (i = 0; i < nmlests; i++) {
		if (mymlests[i] > 0) 
			myposmlests[myll++] = mymlests[i];
	}
//	printf("%d %d\n", myl, myll);

	FREE_3ARRAY(X);

	if (l == 0) 
		ests[0] = 0;
	else 
		ests[0] = findMode(posests, l);
	if (ll == 0) 
		ests[1] = 0;
	else
		ests[1] = findMode(posmlests, ll);
	if (myl == 0) 
		ests[2] = 0;
	else 
		ests[2] = get_mode_from_binned_data(myposests, myl);
	if (myll == 0) 
		ests[3] = 0;
	else 
		ests[3] = get_mode_from_binned_data(myposmlests, myll);

	FREE_VECTOR(mylocests);
	FREE_VECTOR(myposests);
	FREE_VECTOR(mymlests);
	FREE_VECTOR(myposmlests);
	FREE_VECTOR(locests);
	FREE_VECTOR(posests);
	FREE_VECTOR(mlests);
	FREE_VECTOR(posmlests);

	return;
}
