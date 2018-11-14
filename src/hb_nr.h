/* hb_nr.h
..............................................................
                      HydroBase3 
		      Ruth Curry
		      March 2012
..............................................................
  Definitions for matrix operations adapted from
  Numerical Recipes in C
  
*/

#ifndef _NR_H_
#define _NR_H_

#define TINY 1.0e-20

int invert(float **a, int n);
void lubksb(float **a, int n, int *indx, float b[]);
int ludcmp(float **a, int n, int *indx, float *d);

#endif /* _NR_H_ */
