/*
 * This is copied here from the source code for version
 * 2.6 of R as it is not part of the public API.
 * I have changed the name of the function to avoid
 * conflicts with the symbol in the R library.
 *
 *
 * Added modification from Numerical Recipes in C which performs minimization 
 * with first derivatives (Ranjan Maitra, Ames, IA, 2008/04/16)
 *
 *
 * Added modification such that the value at the minimum is also returned
 *
 */

#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))  /*for one-dimensional
						      Brent minimization 
						      with derivatives*/
#define MOV3(a, b, c,  d, e, f) (a) = (d); (b) = (e); (c) = (f); 
                    /*for one-dimensional Brent minimization with derivatives*/

/*
 *  R : A Computer Language for Statistical Data Analysis

 *  Copyright (C) 1997-2001  The R Development Core Team
 *  Copyright (C) 2003-2004  The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* fmin.f -- translated by f2c (version 19990503).
*/

/* R's  optimize() :   function	fmin(ax,bx,f,tol)
   =    ==========		~~~~~~~~~~~~~~~~~

        an approximation  x  to the point where  f  attains a minimum  on
    the interval  (ax,bx)  is determined.

    INPUT..

    ax    left endpoint of initial interval
    bx    right endpoint of initial interval
    f     function which evaluates  f(x, info)  for any  x
          in the interval  (ax,bx)
    fmin  minimized function value
    tol   desired length of the interval of uncertainty of the final
          result ( >= 0.)

    OUTPUT..

    fmin  abcissa approximating the point where  f  attains a minimum

        The method used is a combination of  golden  section  search  and
    successive parabolic interpolation.  convergence is never much slower
    than  that  for  a  Fibonacci search.  If  f  has a continuous second
    derivative which is positive at the minimum (which is not  at  ax  or
    bx),  then  convergence  is  superlinear, and usually of the order of
    about  1.324....
        The function  f  is never evaluated at two points closer together
    than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
    root  of  the  relative  machine  precision.   if   f   is a unimodal
    function and the computed values of   f   are  always  unimodal  when
    separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
    the abcissa of the global minimum of  f  on the interval  ax,bx  with
    an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
    then fmin may approximate a local, but perhaps non-global, minimum to
    the same accuracy.
        This function subprogram is a slightly modified  version  of  the
    Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
    Minimization without Derivatives, Prentice-Hall, Inc. (1973).
*/
#include <math.h>
#include<stdio.h>
#include <float.h> /* DBL_EPSILON */

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double *fmin, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    (*fmin) = fx;

    return x;
}

double Brent_fderfmin(double ax, double bx, double (*f)(double, void *),
		      double (*df)(double, void *), void *info, double *fmin, 
		      double tol)
{
	/* same variables as previously, only exception is that there is a 
	   function pointer for the derivative (*df)(double, void *) */

	int ok1, ok2;
	int iter, ITMAX = 1000;
	double a, b, d, d1, d2, du, dv, dw, dx, e=0.0, fu, fv, fw, fx, olde, 
		tol1, tol2, u, u1, u2, v, w, x, xm, ZEPS = DBL_EPSILON;
	const double c = (3. - sqrt(5.)) * .5;

	

	a = (ax < bx ? ax : bx);
	b = (ax > bx ? ax : bx);

	v = a + c * (b - a);
	x = v;
	w = v;
	fx = (*f)(x, info);
	fv = fx;
	fw = fx;
	dx = (*df)(x, info);
	dv = dx;
	dw = dx;

	d = 0.;/* -Wall */

	for (iter = 1; iter <= ITMAX; iter++) {
//	for (;;) {
		xm = 0.5*(a+b);
		tol1 = tol*fabs(x)+ZEPS;
		tol2 = 2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) return x;
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx)  d1=(w-x)*dx/(dx-dw);
			if (dv != dx)  d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u, info);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*f)(u, info);
			if (fu > fx) return x;
		}
		du=(*df)(u, info);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	printf("Too many iterations in routine DBRENT");
	
	(*fmin) = fx;

	return x;
}

#undef SIGN
#undef MOV3
