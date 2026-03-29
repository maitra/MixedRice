#ifndef __NMMIN_H__
#define __NMMIN_H__

/* int n               = number of parameters
   double *Bvec        = initial values (pointer to 1-D array of p elements)
   double *X           = final value (pointer to 1-D array of p elements)
   double *Fmin        = value at located minimum
   double fminfn       = objective function
   int fail            = 1 if no convergence in maximum number of iterations
   const double abstol = 1e-16;
   const double reltol = 1e-8;
   void *ex            = pointer to external data (perhaps as a struct)
   const double alpha  = 1.0 default  reflection factor
   const double beta   = 0.5 default contraction factor 
   const double gamm   = 2.0 default expansion factor 
   const int trace     = 0; default tracing on if 1.
   int fncount         = number of function evaluations/iterations
   int maxit           = maximum number of iterations  */
void nelder_mead_min(int n, double *Bvec, double *X, double *Fmin,
		     double (*fminfn)(int, const double *, const void *), 
		     int *fail, double abstol, 
		     double intol, const void *ex, double alpha, double bet, 
		     double gamm, int trace, int *fncount, int maxit);

/* Nelder-Mead optimization:
   Modified from R's source code (version 2.4.0). Please see at the end of 
   these comments for copyright information.

   Modification Date: 2006/10/21 

   Changes: - removed headers from Rmath.h (and made standalone) 
            - brought in some sort of documentation to help in calling.
            - renamed function to nelder_mead_min to avoid future 
	      potential conflicts in common calls to R.
	    - freed the matrix of vertices P (for some reason, not in the code)
	    
   This program is modified from the source code for the R statistical software
   package, simply so that we can minimize objective functions with external
   fields such as data in them. The modification is minor and listed above.
   	    
   Therefore, the following copyright probably applies. Regardless, the credits
   should go the R team.

 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2006  the R Development Core Team
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street Fifth Floor, Boston, MA 02110-1301  USA
 */

#endif

