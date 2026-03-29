#include "dlrice.h"
#include <math.h>
#define SQ(x) ((x) * (x))
#define Inf 1e+140

double dlrice(double x, double sigma, double mu) 
{
/* This density also has as its special case (mu = 0), the Rayleigh density */

	double sigmaSq = SQ(sigma);

	if (x > 0) {
		const double xv = x / sigmaSq;
		if (mu > 0) {
			return log(xv) - (SQ(x) + SQ(mu))/(2 * sigmaSq) + logBesselI0(xv * mu);
		}
		else {
			if (mu == 0) return log(xv) - SQ(x)/(2 * sigmaSq);
			else return -Inf;
		}
	}
	else return -Inf;
}

