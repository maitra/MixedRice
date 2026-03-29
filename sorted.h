#include <stdlib.h>
#include "order.h"
#include "array.h"

#define RETURN_CMP(a,b) if (a < b) { return -1; } \
	else if (a == b) { return 0; }		  \
	else { return 1; }

int compareDouble(const void* v1, const void* v2);
/* written by David Faden. All rights reserved. */

int sort(int n, double *x);
/* This function sorts the n-dimensional array x in increasing order. It uses
   the standard library function qsort().
   The input array is replaced by the sorted array on output. 
   
   Written by Ranjan Maitra, Ames, IA 50014.
   2005/09/16. All rights reserved. */

int mdimsort(int n,int p,double **x,int sortdim);
/* This function sorts the columns of the elements of an n x p matrix, in
   increasing order in the sortdim'th dimension. The input array is replaced 
   by the sorted array on output. 
   
   Written by Ranjan Maitra, Ames, IA 50014.
   2005/09/16. All rights reserved. */
