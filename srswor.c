#include <assert.h>
#include <stdlib.h>
#ifndef USING_RLIB
#define MATHLIB_STANDALONE 1 /*It is essential to have this before the call 
                               to the Rmath's header file because this decides
                               the definitions to be set. */
#endif
#include <Rmath.h>

/* Note that use of this function involves a prior call to the Rmath library to
   get the seeds in place. It is assumed that this is done in the calling 
   function. */

/* Equal probability sampling; without-replacement case */
/* Adapted from the R function called SampleNoReplace */

/**
 * Stores k out of n indices sampled at random without replacement
 * in y.
 */
int srswor(int n, int k, int *y)
{ 
	if (k > n) {
		return 1;
	}
	else {
		const double len = (double) n;
		int i;
		int* x = malloc(n * sizeof(int));
		if (!x) {
			return 2;
		}
		
		for (i = 0; i < n; ++i)	x[i] = i;
    
		for (i = 0; i < k; ++i) {
			const int j = (int)(len * runif(0.0, 1.0));
			y[i] = x[j];
			x[j] = x[--n];
		}
		free(x);
	}
	return 0;
}

#define SWAP(y, i, j) do {			\
		const double oldYi = (y)[(i)];	\
		(y)[(i)] = (y)[(j)];		\
		(y)[(j)] = oldYi;		\
	} while(0)

/**
 * Swaps sampleSize elements to the range y[0:(sampleSize-1)]
 * so that this is a random sample without replacement of
 * y (and so that y[sampleSize:len] is the remainder).
 */
void swapSample(double* y, int len, int sampleSize)
{
	assert(len >= sampleSize);
	assert(sampleSize >= 0);

	if (sampleSize != len) {
		int remainder = len - sampleSize;
		if (sampleSize >= remainder) {
			int i;
			for (i = 0; i != sampleSize; ++i, --remainder) {
				const int j = i + (int)(remainder * runif(0.0, 1.0));
				SWAP(y, i, j);
			}
		}
		else {
			int i;
			for (i = len - 1; i != sampleSize; --i) {
				const int j = (int)(i * runif(0.0, 1.0));
				SWAP(y, i, j);
			}
		}
	}
	/* There is nothing to do if sampleSize == len */
}

/**
 * Takes a random sample without replacement from y of size sampleSize
 * and places it in sample.
 *
 * y -- n length array
 * sample -- sampleSize length array
 */
void sampleInto(const double* y, int n, int sampleSize, double* sample)
{
	int i;
	int* x = malloc(sampleSize * sizeof(int));

	assert(sampleSize <= n);
	assert(sampleSize >= 0);

	srswor(n, sampleSize, x);
	for (i = 0; i < sampleSize; ++i) {
		sample[i] = y[x[i]];
	}
	free(x);
}

/**
 * Returns a dynamically-allocated array containing a sample
 * without replacement from y of size sampleSize.
 */
double* sample(const double* y, int n, int sampleSize)
{
	double* result = malloc(sizeof(double) * sampleSize);
	sampleInto(y, n, sampleSize, result);
	return result;
}

