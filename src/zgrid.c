/* zgrid.c
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
*
*  Parts of these functions were adapted from GMT-SYSTEM
*	Copyright (c) 1991-2000 by P. Wessel and W. H. F. Smith
*  and from PLOT+ Scientific Graphics System software
*
................................................................................
 * Modules in this file:
 *
 *      zgrid()       Laplacian gridding/interpolating algorithm
 ................................................................................
*/ 
 
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <float.h>
#include <stddef.h>
#include "hb_grids.h"
#include "zgrid.h"
#include "hb_memory.h"
/****************************************************************************/

double * zgrid (double converge_lim, int itmax, int N, int xrad, int yrad, double tension, struct GRID_INFO *hptr, int *knxt, int *imnew, double *zpij, double *oldz, struct POINT *point, int *nmask_addr, int *nempty_addr, int *nset_addr)

/*
converge_lim:  Convergence criteria. 
       itmax:  Maximum iteration for convergence. 
           N:  number of data pts starting at point. 
        xrad:  search radius in x-direction [number of gridpts]
        yrad:  search radius in y-direction [number of gridpts]
     tension:  0 -> 1 tension parameter (0 gives pure spline, 1 gives pure laplacian)
        hptr:  pointer to gridinfo
        knxt:  work array 
       imnew:  work array
        zpij:  work array
        oldz:  input grid -- contains 0's or DBL_MAX for masked nodes.
       point:  pointer to start of input data array 
  nmask_addr:  returned number of masked gridnodes
 nempty_addr:  returned number of empty gridnodes 
   nset_addr:  returned number of gridnodes with data values 

>>> On input, zgrid() expects masked values to be DBL_MAX, all other 
>>> points to be zero.
>>> On return, empty nodes will be set to the value 1.0E+35
>>> Masked nodes will be set to the value DBL_MAX (machine dependent).
>>>  The calling routine can reset these values to something else.
>>>  The number of grid nodes set, empty, and masked are returned.
*******************************************************************************
>  The guts of this procedure were adapted from a subroutine in PLOT+ Scientific
>  Graphics System (see below).

>     1)  Masking information is considered within the interpolating procedure.
          If a masked node is encountered, information beyond that node is
	  not incorporated into the computation.  
>     2)  Also crucial is the ordering of Z-array elements: first point 
>         at [xmin,ymin], stored by row such that sq = row * ncols + col

>  The relative weighting of laplacian/biharmonic equations (cayin in the
>  original code) was adapted to reflect a tension parameter 
>  (between 0 and 1) where tension = 0 gives pure spline and 
>  tension = 1 gives pure laplacian solution.  It is a closer parallel to
>  the way the GMT module "surface" (and hb_fitsurf2d) works and I have 
>  implemented the tension parameter to avoid conflict.
>        The original cayin = (1-tension)/ tension.
>
>  Initialization procedure is as follows:
>  NRNG (original search radius) was differentiated into x and y radii. 
>  Empty nodes are initialized as weighted average of immediate neighbors within 
>  x-,y-radius. 
>  JUNE2005:  Updated to search diagonally (northwest/northeast/southwest/southeast) 
>  as well as along rows/cols (west/north/east/south). If masked node is encountered 
>  in search for neighbor, search is 
>  suspended in that direction.  The initialized values are stored in a new
>  array (newz[]) so that the old array values can be preserved.  Those old
>  values are iteratively restored to the newz array.
>  
  
>   -R.Curry
*********************************************************************/

/*********************************************************************
**  Routine zgrid.c is a C language version of the FORTRAN language
**  subroutine of the same name found in the PLOT+ Scientific Graphics 
**  System.  The original boilerplate follows.
**    @(#)zgrid.f	1.3    3/30/93
**
***********************************************************************
**
**                 PLOT+ Scientific Graphics System
**
***********************************************************************
**
**
**     Sets up square grid for contouring, given arbitrarily placed 
**     data points. Laplace interpolation is used. 
**     The method used here was lifted directly from notes left by 
**     Mr Ian Crain formerly with the comp.science div. 
**     Info on relaxation solutionn of laplace equations  supplied by 
**     Dr T Murty.   FORTRAN II   Oceanography/EMR   Dec 68   JDT 
** 
**     z = 2-d array of nodes to be gridded/interpolate points to be masked
**      should be initialized to 10**35 . The rest should be 0.0 
**     nx,ny = max subscripts of z in x and y directions. 
**     x1,y1 = coordinates of z(1,1) 
**     dx,dy = x and y increments . 
**     xp,yp,zp = arrays giving position and hgt of each data point. 
**     n = size of arrays xp,yp and zp . 
** 
**     modification feb/69   to get smoother results a portion of the 
**     beam eqn  was added to the laplace eqn giving 
**     delta2x(z)+delta2y(z)) - k((delta4x(z)+delta4y(z)) = 0 . 
**     k=0 gives pure laplace solution; k=inf gives pure spline solution.
**     cayin = k = amount of spline eqn (between 0 and inf). 
**
**     NRNG...grid points more than NRNG grid spaces from the nearest 
**            data point are set to undefined. 
** 
**     Modification Dec 23 69   data pts no longer moved to grid pts. 
** 
**     Modification may 5 79  common blocks work1 and work2 must 
**     be dimension at least n points long by the user.  Common 
**     block work3 must be dimensioned at least ny points long. 
** 
**	Modification Jun 17 1985 - handles data values of 1e35.  If at
**	least one data value near a grid point is equal to 1e35, the z
**	array is initialized to 1e35 at that grid point.
**	- by G.R. Halliwell
**********************************************************************
 *  The program was modified to use single dimension subscript for the
 *  output array Z.  The subscripting was modified to reflect the
 *  0 base C language convention but the first index, column based
 *  array FORTRAN language convention was maintained in the gridding
 *  procedure .
 *
 *   Adapted by:  R. Goldsmith
 *   Date:  June, 1999
***********************************************************************/
{
int     i;
int     im, imask;
int     iter;
int     j, jm, jmnew;
int     k, kk;
int     kksav[SSIZE];
int     nnew, npg, npt;
int     NX, NY;
int     row, col, sq, nsq;
int     found, index, rad;
int     jrow, icol;

double  DX, DY, X1, Y1;
double  a, b, c, d;
double  abz, cayin;
double  big, bigger, biggest, zbig;
double checkmin, checkmax;
double  delz, delzm;
double  derzm;
double  dz, dzmax, dzmaxf, dzrms, dzrmsp, dzrms8;
double  hrange;
double  relax, relaxn;
double  root, rootgs;
double  tpy;
double  weight, maxdist;
double  weightsum;
double  x;
double  y;
double  z00;
double  zbase, zinit;
double  zijn;
double  zim, zimm, zip, zipp;
double  zjm, zjmm, zjp, zjpp;
double  zmin, zmax;
double  zpxy;
double  zrange, zrms;
double  zsum;
double  zw, ze, zs, zn;
double  zxy, *newz;

/* use these original fortran variables to avoid introducing a 
   transcription error */
   
  NX = hptr->nx;    
  NY = hptr->ny;
  DX = hptr->x_inc;
  DY = hptr->y_inc;
  X1 = hptr->x_min;
  Y1 = hptr->y_min;
  if (hptr->node_offset) {
     X1 = hptr->x_min + 0.5 * hptr->x_inc;
     Y1 = hptr->y_min + 0.5 * hptr->y_inc;
  }
  nsq = NX * NY;
  newz = (double *) get_memory((void *)NULL, (size_t)nsq, sizeof(double));
  
/* adapt the tension parameter to the PLOT+ method of weighting 
   the laplacian and biharmonic components */

  cayin = 1.0;   
  if (tension > 0.0)   
     cayin *= (1 - tension) ;


/*  Set up some internal representations for masked and empty nodes.
    We assume no hydrographic property has a value exceeding 0.9E+35 and that
   the max double value (machine dependent)is something on the order 1.0E+308 */


big = 0.9E+35;
bigger = 1.0E+35;
biggest = DBL_MAX;

 
     /* Get zbase which will make all zp values positive by 20*(zmax-zmin) */
zmin =  big;
zmax = -big;
i = 0;
for (k = 0; k < N; k++)  {
  if (point[k].ZP < big)  {
    if (point[k].ZP > zmax)  zmax = point[k].ZP;
    if (point[k].ZP < zmin)  zmin = point[k].ZP;
    zsum += point[k].ZP * point[k].ZP;
    ++i;
  }
}
zrange = zmax - zmin;
zbase = zrange * 20.0 - zmin;
hrange = DX*(NX-1);
if (DY*(NY-1) < hrange)  hrange = DY*(NY-1);
derzm = 2.0 * zrange / hrange;

zrms = sqrt(zsum / i);   /* zsum is sum of squared z-values */
   

/* printf ("ZRANGE: %lf %lf %lf %lf\n", 
         zrange, zbase, hrange, derzm); */
     
  /* Set pointer array knxt.  Test that subscripts are in range. */
     
for (kk=0; kk < N; kk++)  {
  k =  N - 1 - kk;
  knxt[k] = 0;
  col = (int) NINT((point[k].XP - X1) / DX);
  if ((col < NX) && (col >= 0))  {
    row = (int) NINT((point[k].YP - Y1) / DY );
    if ((row < NY) && (row >= 0))  {
      sq = col + NX * row;
      if (oldz[sq] <  big)  {    /* check for mask */
        knxt[k] = N;
        if (oldz[sq] > 0)        /* check for multiple entries at this gridnode */
	  knxt[k] = (int) NINT(oldz[sq]);
        oldz[sq] = k;
      }
    }
  }
}

 
     /* Affix each data point ZP to its nearby grid point.  Take avg ZP if 
     ** more than one ZP nearby the grid point.  Add zbase and complement. 
     */
for (k = 0; k < N; k++)  {
  if (knxt[k] > 0)  {
    npt = 0;
    imask = 0;

    zsum = 0.0;
    col = (int) NINT((point[k].XP - X1) / DX);
    row = (int) NINT((point[k].YP - Y1) / DY);
    sq = col + NX * row;
    kk = k;

    while (kk < N  ||  npt == 0)  {
      if (npt >= SSIZE)  {
        fprintf (stderr, "zgrid error: internal array size SSIZE (%d) not\n", 
                 SSIZE);
        fprintf (stderr, "  large enough to handle initialization.\n", N);
        exit(1);
        }
      kksav[npt++] = kk;
      if (point[kk].ZP > big)  imask = 1;

      zsum = zsum + point[kk].ZP;
      knxt[kk] = -knxt[kk];
      kk = -knxt[kk];
      }  /* end while */

    if (!imask )  {
      oldz[sq] = -zsum / (double) npt - zbase;
    }
    else  {
      oldz[sq] = biggest;
      for (i=0; i < npt; i++)  {
        knxt[kksav[i]] = 0;
      }
    }  /* end if */
  }  /* end if */
}  /* end for */

  checkmin = fabs(-zmin - zbase);
  checkmax = fabs(-zmax - zbase);

/* Initially set each empty grid point to a large negative value.... 
   and copy the old grid values into the new grid. */ 

for (i = 0; i < NX; i++)  {
  for (j=0; j < NY; j++)  {
    sq = i+NX*j;
    if (oldz[sq] == 0)  oldz[sq] = -bigger;
    newz[sq] = oldz[sq];
  }
}
 
/* ... now fill in empty nodes of new grid with average of nearest neighbors.
   As each column is evaluated, update oldz so next column can use those values.
   After completing this loop, oldz will no longer be used. */

maxdist = (double)  xrad;
if (xrad < yrad)
    maxdist = yrad;
    
maxdist *= 1.4;   /* max radius * sqrt(2) */ 
    
for (iter = 0; iter < itmax; ++iter) {
  for (i = 0; i < NX; i++)  {
    for (j = 0; j < NY; j++)  {
  
      sq = i+NX*j;
      if (oldz[sq] + big < 0) {  /* very negative value marks an empty node */
       
	/* search west */
	
        rad = 0;
        weightsum = 0.0;
        zsum = 0.0;
	icol = i;
	jrow = j;
	found = 0;
	while ( (--icol >= 0) && (jrow < NY) && (++rad <= xrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = xrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 -  (double) rad / maxdist;      /* weight by ratio of distance to maxdist in search radius */
	      zsum += fabs(oldz[index])  * weight;
	      weightsum += weight;
	    }
	  }
	}
	
	/* search east */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while (  (++icol < NX) && (jrow < NY) && (jrow >= 0)   && (++rad <= xrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = xrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 -  (double) rad / maxdist;      /* weight by ratio of distance to maxdist in search radius */
	      zsum += fabs(oldz[index])  * weight;
	      weightsum += weight;
	    }
	  }
	}
	
	/* search north */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while (  (icol < NX) && (icol >= 0) && (++jrow < NY) && (jrow >= 0)   && (++rad <= yrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 -  (double) rad / maxdist;      /* weight by ratio of distance to maxdist in search radius */
	      zsum += fabs(oldz[index])  * weight;
	      weightsum += weight;
	    }
	  }
	}
	
	
	/* search south */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while (  (icol < NX) &&  (icol >= 0)  && (--jrow >= 0)  && (++rad <= yrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 -  (double) rad / maxdist;      /* weight by ratio of distance to maxdist in search radius */
	      zsum += fabs(oldz[index])  * weight;
	      weightsum += weight;
	    }
	  }
	}
	
	/* search southeast */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while ( (++icol < NX) &&  (icol >= 0)  && (--jrow >= 0)   && (++rad <= yrad)  && (rad <= xrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 -  (double) rad * 1.4 / maxdist;      /* weight by ratio of distance to maxdist in search radius */
	      weightsum += weight; 
	      zsum += fabs(oldz[index]) * weight;
	    }
	  }
	}
	
	
	/* search southwest */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while ( (--icol >= 0)  && (--jrow >= 0)  && (++rad <= yrad)  &&( rad <= xrad) && !found) {  
	  index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1.0 - (double) rad  * 1.4  / maxdist;   /* weight by ratio of distance to maxdist in search radius */
	      weightsum += weight; 
	      zsum += fabs(oldz[index]) * weight;
	    }
	  }
	}
	
	/* search northwest */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while (  (--icol >= 0)  && (++jrow < NY)  && (++rad <= yrad)  &&( rad <= xrad) && !found) {  	 
	 index = icol + NX * jrow;

	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1. /  ((double) rad  * 1.4  / maxdist);   /* weight by ratio of distance to maxdist in search radius */
	      weightsum += weight; 
	      zsum += fabs(oldz[index]) * weight;
	    }
	  }
	}
	
	/* search northeast */
        rad = 0;
	found = 0;
	icol = i;
	jrow = j;
	while (  (++icol < NX)  && (++jrow < NY)  && (++rad <= yrad)  &&( rad <= xrad) && !found) {  	 
	 index = icol + NX * jrow;
	  if (oldz[index] + big > 0) {  /* not another empty node */
	    if (oldz[index] > big) {  /* masked val -- stop search */
	      rad = yrad;
	    }
	    else {
	      found = 1;
	      weight = 1. /  ((double) rad  * 1.4  / maxdist);   /* weight by ratio of distance to maxdist in search radius */
	      weightsum += weight; 
	      zsum += fabs(oldz[index]) * weight;
	    }
	  }
	}
	
	
	if (weightsum > 0) 
	   newz[sq] = zsum / weightsum;
	
	
     }/* end if oldz[sq] + big */

   } /* end for j */
    
 } /* end for i */

 for (k=0; k < NY; k++)  {
   for (i=0; i < NX; ++i) {
        sq = i+NX*k;
        oldz[sq] = newz[sq];
   }
 } 
  
} /* end for iter */


/*
  At this point, nodes initialized with real data have negative values: -(zreal +zbase).
  Nodes initialized as weighted sum are now positive:  zmean + zbase.
  Masked nodes are very positive = DBL_MAX; 
  Now, change nodes that are still empty from very negative value to very positive value
  (but still smaller than the  value of masked nodes). 
  
  Only real data nodes should be negative after this.
  
  */

for (i = 0; i< NX; i++)  {
  for (j=0; j < NY; j++)  {
    sq = i+NX*j;
    abz = fabs (newz[sq]);
    if (abz >= big)  {
      newz[sq] = abz;
    } 
  }
}



/* Improve the non-data points by applying point over-relaxation 
      using the Laplace-spline equation  (Carres method) */
      
dzrmsp = zrange;
relax = 1.0;
for (iter=1; iter < itmax; ++iter)  {
  dzrms = 0.0;
  dzmax = 0.0;
  npg = 0;
  for (i=0; i<NX; i++)  {
    for (j=0; j<NY; j++)  {
      sq = i+NX*j;
      z00 = newz[sq];
      if (z00 < big)  {
        if (z00 >= 0)  {   /* this tests for non-data node */
          weight = 0.0;
          zsum = 0.0;

          im = 0;
          if (i > 0)  {
            zim = fabs (newz[(i-1)+NX*j]);
            if (zim < big)  {
              im = 1;
              weight = weight + 1.0;
              zsum = zsum + zim;
              if (i > 1)  {
                zimm = fabs(newz[(i-2)+NX*j]);
                if (zimm < big)  {
                  weight = weight + cayin;
                  zsum = zsum - cayin * (zimm - 2.0*zim);
                  }
                }
              }
            }  /* end if i > 0 */
          if (i < NX-1)  {
            zip = fabs (newz[(i+1)+NX*j]);
            if (zip < big)  {
              weight = weight + 1.0;
              zsum = zsum + zip;
              if (im > 0)  {
                weight = weight + 4. * cayin;
                zsum = zsum + 2. * cayin * (zim + zip);
                }
              if (i < NX-2)  {
                zipp = fabs (newz[(i+2)+NX*j]);
                if (zipp < big)  {
                  weight = weight + cayin;
                  zsum = zsum - cayin * (zipp - 2.0*zip);
                  }
                }
              }
            }  /* end if i < NX-1 */

          jm = 0;
          if (j > 0)  {
            zjm = fabs (newz[i+NX*(j-1)]);
            if (zjm < big)  {
              jm = 1;
              weight = weight + 1.0;
              zsum = zsum + zjm;
              if (j > 1)  {
                zjmm = fabs (newz[i+NX*(j-2)]);
                if (zjmm < big)  {
                  weight = weight + cayin;
                  zsum = zsum - cayin * (zjmm - 2.0*zjm);
                  }
                }
              }
            }  /* end if j > 1 */
          if (j < NY-1)  {
            zjp = fabs (newz[i+NX*(j+1)]);
            if (zjp < big)  {
              weight = weight + 1.0;
              zsum = zsum + zjp;
              if (jm > 0)  {
                weight = weight + 4.0 * cayin;
                zsum = zsum + 2.0 * cayin * (zjm + zjp);
                }
              if (j < NY-2)  {
                zjpp = fabs (newz[i+NX*(j+2)]);
                if (zjpp < big)  {
                  weight = weight + cayin;
                  zsum = zsum - cayin * (zjpp - 2.0*zjp);
                  }
                }
              }
            }  /* end if j < NY-1  */
	    
          if (weight != 0.0) {
            dz = zsum / weight - z00;
            npg = npg + 1;
            dzrms = dzrms + dz*dz;
            if (fabs (dz) > dzmax)  dzmax = fabs (dz);
            newz[sq] = z00 + dz * relax;
	    if (newz[sq] > checkmax)
	        newz[sq] = checkmax;
	    if (newz[sq] < checkmin)   
	       newz[sq] = checkmin;
	  }
        }  /* end if z00 >= 0  */
      }  /* end if z00 < big */
    }  /* end for j loop */
  }  /* end for i loop */

 /* Shift data points zp progressively back to their proper places as 
      the shape of surface z becomes evident.  */
if (iter % 10 == 0) {
    for (k=0; k<N; k++)  {
      if (knxt[k] < 0)  knxt[k] = -knxt[k];
      if (knxt[k] > 0)  {
        x = (point[k].XP - X1) / DX;
	col = (int) NINT(x);
        x = x - (double) col;
        y = (point[k].YP - Y1) / DY;
	row = (int) NINT(y);
        y = y  - (double)row;
        zpxy = point[k].ZP + zbase;
        z00 = fabs (newz[col+NX*row]);

        zw = bigger;
        if (col > 0)  {
          zw = fabs (newz[(col-1)+NX*row]);
          }
        ze = bigger;
        if (col < NX-1)  {
          ze = fabs (newz[(col+1)+NX*row]);
          }
        if (ze >= big)  {
          if (zw < big)  {
            ze = 2.0 * z00 - zw;
            }
           else  {
            ze = z00;
            zw = z00;
            }
          }
         else  {  /* if ze < big  */
          if (zw >= big)  {
            zw = 2.0 * z00 - ze;
            }
          }

        zs = bigger;
        if (row > 0)  {
          zs = fabs (newz[col+NX*(row-1)]);
          }
        zn = bigger;
        if (row < NY-1)  {
          zn = fabs (newz[col+NX*(row+1)]);
          }
        if (zn >= big)  {
          if (zs < big)  {
            zn = 2.0 * z00 - zs;
            }
           else  {
            zn = z00;
            zs = z00;
            }
          }
         else  {  /* zn < big */
          if (zs >= big)  {
            zs = 2.0 * z00 - zn;
            }
          }

        a = (ze - zw) * 0.5;
        b = (zn - zs) * 0.5;
        c = (ze + zw) * 0.5 - z00;
        d = (zn + zs) * 0.5 - z00;
        zxy = z00 + a*x + b*y + c*x*x + d*y*y;
        delz = z00 - zxy;
        delzm = derzm * (fabs (x)*DX + fabs (y)*DY) * 0.80;
        if (delz > delzm)  {
          delz = delzm;
          }
        if (delz+delzm < 0)  {
          delz = -delzm;
          }
        zpij[k] = zpxy + delz;
        }  /* end if knxt[k] > 0  */
      }  /* end for k loop */

 
    for (k=0; k < N; k++)  {
      if (knxt[k] > 0)  {
        npt = 0;
        zsum = 0.0;
        col = (int) NINT((point[k].XP-X1) / DX);
        row = (int) NINT((point[k].YP-Y1) / DY);
        kk = k;
        while (kk < N  ||  npt == 0)  {
          ++npt;
          zsum = zsum + zpij[kk];
          knxt[kk] = -knxt[kk];
          kk = -knxt[kk];
          }  /* end while  (kk < N) */
        if (npt > 0)  newz[col+NX*row] = -zsum / (double) npt;
        }
      }  /* end for k loop */
    }  /* end if (iter % 10 == 0)  */


     /* Test for convergence */
 
  if (npg == 0)  break;
  dzrms = sqrt(dzrms/npg);
/*
  fprintf(stderr,"  convergence (npg, dzrms): %ld %lf\n", npg, dzrms);
*/
  root = dzrms / dzrmsp;
  dzrmsp = dzrms;
  dzmaxf = dzmax / zrange;
  if (iter % 10 == 2)  {
    dzrms8 = dzrms;
    }
  if (iter % 10 == 0)  {
    root = sqrt(sqrt(sqrt(dzrms/dzrms8)));
    if (root < 0.9999)  {
      if (dzmaxf/(1.0 - root) - converge_lim <= 0)  break;
 
          /* Improve the relaxation factor. */
 
      if ((iter % 20)*(iter % 40)*(iter % 60) == 0)  {
        if (relax - 1.0 - root < 0)  {
          tpy = (root + relax - 1.0) / relax;
          rootgs = tpy * tpy / root;
          relaxn = 2.0 / (1. + sqrt(1.0 - rootgs));
          if (iter != 60)  {
            relaxn = relaxn - 0.25 * (2.0 - relaxn);
            } 
          if (relaxn > relax)  relax = relaxn;
          }
        }
      }
    }
  }  /* end for iter loop */
 

  if (iter > itmax)  iter = itmax;
/*  printf ("Convergence = %lg after %ld iterations.\n", dzmaxf/(1.0-root), iter); */

     /* Remove zbase from array z, count empty and masked nodes 
        and return. */
*nset_addr = 0;
*nmask_addr = 0;
*nempty_addr = 0;
for (i = 0; i< NX; i++)  {
  for (j = 0; j < NY; j++)  {
    sq = i+NX*j;
    if (newz[sq] < big)  {
      newz[sq] = fabs (newz[sq]) - zbase;
      ++(*nset_addr);
      }
    else if (newz[sq] < biggest-1) {
      ++(*nempty_addr);
    }
    else {
      ++(*nmask_addr);
    }
  } 
}

return (newz);
}  /* end zgrid() */

/************************************************************** */
