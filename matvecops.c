#include "array.h"

#include "matvecops.h"

double *range(double *x, int n)
{
	double *ra;
	int i;

	MAKE_VECTOR(ra, 2);

	ra[0] = x[0];
	ra[1] = x[0];

	for (i = 1; i < n; i++) {
		if (x[i] < ra[0]) ra[0] = x[i];
		else {
			if (x[i] > ra[1]) ra[1] = x[i];
		}
	}
	
	return ra;
}

