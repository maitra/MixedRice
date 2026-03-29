#ifndef __SRSWOR_H__
#define __SRSWOR_H__

/**
 * Stores k out of n indices sampled at random without replacement
 * in y.
 */
int srswor(int n, int k, int *y);

/**
 * Swaps sampleSize elements to the range y[0:(sampleSize-1)]
 * so that this is a random sample without replacement of
 * y (and so that y[sampleSize:len] is the remainder).
 */
void swapSample(double* y, int len, int sampleSize);

/**
 * Takes a random sample without replacement from y of size sampleSize
 * and places it in sample.
 *
 * y -- n length array
 * sample -- sampleSize length array
 */
void sampleInto(const double* y, int n, int sampleSize, double* sample);

/**
 * Returns a dynamically-allocated array containing a sample
 * without replacement from y of size sampleSize.
 */
double* sample(const double* y, int n, int sampleSize);

#endif

