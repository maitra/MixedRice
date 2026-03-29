/* Matrix inverse. Use posymatinv for inversion of a positive-definite 
   matrix. pposymatinv for the same in packed format
  
   Author: Ranjan Maitra <maitra@iastate.edu>
   Date:   03/07/2005
   Uses:   LAPACK 
*/  

#ifndef __INVERSE_H__
#define __INVERSE_H__


#include <stdio.h>              /* for fprintf() */
#include <stdlib.h>
#include <string.h>
#include "array.h"
#include "mat_vec.h"

void dgetrf_(int *Mp, int *Np, double *A, int *LDA, int *PIVOT, int *INFOp);
void dgetri_(int *Np, double *A, int *LDA, int *PIVOT, double *WORK, 
	     int *LWORK, int *INFOp);

void  dpotrf_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);
void  dpotri_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);

void dpptrf_(char *UPLOp,int *Np,double *A,int *INFOp);
void dpptri_(char *UPLOp,int *Np,double *A,int *INFOp);

int matinv(int sizeA,double **A,double (*determinant));

int posymatinv(int size,double **A,double (*determinant));

int pposymatinv(int N,double *A, char UPLO, double *determinant);

double pposymatdet(int N,double *A, char UPLO);

#endif

