#ifndef __EIGENS_H__
#define __EIGENS_H__

#include <stdio.h> 
#include "array.h"

void dsyevr_(char *jobz, char *range, char *uplo,
	     int *n, double *a, int *lda,
	     double *vl, double *vu, int *il, int *iu,
	     double *abstol, int *m, double *w,
	     double *z, int *ldz, int *isuppz,
	     double *work, int *lwork, int *iwork, int *liwork,
	     int *info);

void dspevd_(char *jobz,char *UPLO,int *n,double *ap,double *w,double *z,
	     int *ldz,double *work,int *lwork,int *iwork,int *liwork,int *info);


/* C front end to LAPACK's dspevd_() routine.
 *
 * Calculates the eigenvalues and the eigenvectors of the nxn symmetric matrix 
 * A and returns them in the vector E, and the vector EV, respectively.
 * The i'the eigenvector corresponds to the i'th eigenvalue.
 *
 * Note that the eigenvalues are returned in the descending order.
 * 
 * */

int eigend(double *A, double *EV, double *E, int n);
/* E = eigenvalues
   EV = eigenvectors
   A = real symmetric matrix in packed storage format, of dimension n x n.
*/


/* C front end to LINPACK's dsyevr_() routine.
 *
 * Calculates the eigenvalues and eigenvectors of the nxn symmetric matrix A.
 * The eigenvalues are returned in the vector w.
 * The (orthonormal) eigenvectors are returned in the matrix z.
 * The ith column of z holds the eigenvector associated with w[i].
 * Written by Rouben Rostamian and Ranjan Maitra
 * */
int LP_sym_eigvecs(double *a, int n, double *w, double *z);

/* Front end to LP_sym_eigvecs
 *
 * Calculates the eigenvalues and eigenvectors of the nxn symmetric matrix A.
 * The eigenvalues are returned in the vector w.
 * The (orthonormal) eigenvectors are returned in the matrix z.
 * The ith column of z holds the eigenvector associated with w[i].
 * Written by Ranjan Maitra
 * Note that one major task done here is to reverse the order of the 
 * eigenvalues (to something that makes more sense) and to put them in
 * decreasing order and the corresponding eigenvectors
 * */
int symeigens(double *a, int n, double *w, double *z);

/*
  Calculates eigevectors and eigenvalues for a packed symmetric matrix. 
  The eigenvalues returned are in descending order.
  The returned eigenvectors correspond to the eigenvalues.
 */

int eigens(double *A, double *EVec, double *EVal, int n);

#endif
