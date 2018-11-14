/* matrix_subs.c
*                          HydroBase3
*		           Ruth Curry
*		           WHOI
*		           March 2012
***************************************************************
*
*    Contains libray functions adapted from Numerical Recipes in C
*    to handle matrix inversion
* 
*  int invert(float **a, int n)
*   Uses LU decomposition and backsubstitution to find the inverse
*   of an n-by-n matrix, a.  The original matrix is destroyed, and the
*   inverse is returned in a.
* 
*  int ludcmp(float **a, int n, int *indx, float *d)
*  Replaces matrix n-by-n matrix "a" with the LU decomposition of a 
*  rowwise  permutation of itself.  a and n are input, indx is an 
*  output vector that records the row permutation resulting from partial
*  pivoting. d is output as +/-1 denoting whether the number of row 
*  interchanges was even or odd. Returns 0 if successful, or 1 if the 
*  matrix is singular.
*        
*  void lubksb(float **a, int n, int *indx, float b[])
*  Solves the set of n linear equation A*X = B,  The n-by-n matrix a 
*  is the LU decomposition of A. indx is the permutatio vector returned 
*  by ludcmp.  b is input as the right-hand side vector B, and returned 
*  as the vector solution X.
*  This routine takes into account the possibility that B will begin with many 
*  zero elements (partial pivoting using Crout's algorithm) for efficient 
*  matrix inversion.
*     		    
*/


#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "hb_nr.h"

/*************************************************/
int invert(float **a, int n)
/*
*   Uses LU decomposition and backsubstitution to find the inverse
*   of an n-by-n matrix, a.  The original matrix is destroyed, and the
*   inverse is returned in its place.

*/
{
  int i,j, error, *indx;
  float **y, d, *col;
  
  /* allocate memory */
  
  indx = (int *)calloc(n,sizeof(int));
  col = (float *) calloc(n,sizeof(float));
  y = (float **) malloc(n * sizeof(float *));
  for (i = 0; i < n; ++i) 
      y[i] = (float *) malloc(n * sizeof(float));
      
/* Decompose the matrix once, find inverse by columns */
  
  if (error = ludcmp(a, n, indx, &d))
     return(error);
  
  for (j = 0; j < n; ++j) {
     for (i = 0; i < n; ++i) 
        col[i] = 0.0;
	
     col[j] = 1.0;
     lubksb(a, n, indx, col);
     for (i = 0; i < n; ++i)
        y[i][j] = col[i];
  }
  
/* copy the inverted matrix into the old matrix */

  for (i=0; i<n; ++i) {
     for (j=0; j<n; ++j)
        a[i][j] = y[i][j];
  }
  
  free(indx);
  free(col);
  for (i=0; i<n; ++i)
     free(y[i]);
  free(y);
  
  return(0);
}
/*************************************************/
int ludcmp(float **a, int n, int *indx, float *d)
 /*  Replaces matrix n-by-n matrix "a" with the LU decomposition of a rowwise
    permutation of itself.  a and n are input, indx is an output vector that
    records the row permutation resulting from partial pivoting. d is output 
    as +/-1 denoting whether the number of row interchanges was even or odd.
    Returns 0 if successful, or 1 if the matrix is singular */
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv;

	vv = (float *) malloc(n * sizeof(float));
	*d = 1.0;
	for (i=0; i<n; i++) {
		big = 0.0;
		for (j=0; j<n; j++)
		   if ((temp=fabs(a[i][j])) > big) big=temp;
		   
		if (big == 0.0) {
		   fprintf(stderr, "\nERROR: Singular matrix in function ludcmp()\n");
		   free(vv);
		   return(1);
		}   
		vv[i]=1.0/big;
	}
	for (j=0; j<n; j++) {
		for (i=0; i<j; i++) {
			sum=a[i][j];
			for (k=0; k<i;k++) 
			    sum -= a[i][k]*a[k][j];
			    
			a[i][j] = sum;
		}
		big=0.0;
		for (i=j; i<n; i++) {
			sum=a[i][j];
			for (k=0; k<j; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum = vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k=0; k<n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]= vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j]= TINY;
		if (j != (n-1)) {
			dum= 1.0/(a[j][j]);
			for (i=j+1; i<n; i++) 
			   a[i][j] *= dum;
		}
	}
	free(vv);
	return(0);
	
} /* end ludcmp() */
/*****************************************************/
void lubksb(float **a, int n, int *indx, float b[])
/*   Solves the set of n linear equation A*X = B,  The n-by-n matrix a is the LU decomposition of A. indx is the permutation vector returned by ludcmp().  b is input as the right-hand side vector B, and returned as the vector solution X.  This routine takes into account the possibility that b will begin with many zero elements to be efficient in matrix inversion */
  
{
	int i,iflag,ip,j;
	float sum;
	
        iflag = -1;
	for (i=0; i<n; i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (iflag >= 0)
		   for (j=iflag;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) 
		   iflag=i;
		   
		b[i]=sum;
	}
	for (i=n-1; i>=0; i--) {
		sum=b[i];
		for (j=i+1; j<n; j++) 
		   sum -= a[i][j]*b[j];
		   
		b[i]=sum/a[i][i];
	}
}  /* end lubksb() */

/*****************************************************/
