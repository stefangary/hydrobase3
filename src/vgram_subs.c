/*  vgram_subs.c 
................................................................................
                          *******  HydroBase3  *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             March 2011
			     
			     
................................................................................
.
.  Functions for working with semi-variograms
.  used in HydroBase 3D optimal interpolation routines.  
.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hb_variograms.h"
#include "hb_grids.h"
#include "hb_memory.h"

/***************************************************************************/

double *vgram_vals(int modelCode, double *h, int nh, double *parms)
   /* Returns a vector of  semivariogram values  associated with the input vector of lags
    given the parameters supplied in the 3-element input array, parms.  */
{
   int i;
   double c1,c0,range;
   double *Vh;
   
   if (nh <= 0)  {
      fprintf(stderr, " WARNING:  No lags supplied to vgram_vals() ");
      return(NULL);
   }
      
      
   Vh = (double *)get_memory(NULL, nh, sizeof(double));
   
   c0 = parms[0];
   c1 = parms[1];
   range = parms[2];
  
   
   switch (modelCode)  {
      case 1:   /* exponential model */
         for (i = 0; i < nh; ++i) {
	     if (h[i] == 0)
	          Vh[i] = 0.0;
	     else
	         Vh[i] = c0 + c1 * (1 - exp(-h[i]/range ) );
	 }
         break;
	 
      case 2:    /* gaussian model */
         for (i = 0; i < nh; ++i) {
	     if (h[i] == 0)
	          Vh[i] = 0.0;
	     else
	       Vh[i] = c0 + c1 * (1 - exp(-pow(h[i]/range,2) ) );
	 }
	 break;
	 
      case 3:   /* power model */
         for (i = 0; i < nh; ++i) {
	     if (h[i] == 0)
	          Vh[i] = 0.0;
	     else
	          Vh[i] = c1 * pow(h[i],range);
	 }
           
	 break;
	 
      default:
         fprintf(stderr, "WARNING: Unknown modelCode passed to vgram_vals() ");
         free(Vh);
	 return(NULL);
   
   }  /* end switch */
   
    return(Vh);
    
}   /* end vgram_vals() */


/***************************************************************************/

int vgram_vals2(int modelCode, double *h, int nh, double *parms, double *Vh)

/* Same as vgram_vals() except that the returned array, Vh, is already 
   allocated by calling program.  Returns 0 for successful completion or
   1 for an error*/

{
   int i;
   double c1,c0,range;
   
   c0 = parms[0];   /* nugget */
   c1 = parms[1];    /* sill - nugget */
   range = parms[2];  /* range */
   
   
   if (nh <= 0)  {
      fprintf(stderr, " WARNING:  No lags supplied to vgram_vals() ");
      return(1);
   }
   switch (modelCode)  {
      case 1:   /* exponential model */
         for (i = 0; i < nh; ++i) {
	   if (h[i] == 0)
	          Vh[i] = 0.0;
	   else
	      Vh[i] = c0 + c1 * (1 - exp(-h[i]/range ) );
	 }
         break;
	 
      case 2:    /* gaussian model */
         for (i = 0; i < nh; ++i) {
	   if (h[i] == 0)
	          Vh[i] = 0.0;
	   else
	       Vh[i] = c0 + c1 * (1 - exp(-pow(h[i]/range,2) ) );
	 }
	 break;
	 
      case 3:   /* power model */
         for (i = 0; i < nh; ++i) {
	   if (h[i] == 0)
	          Vh[i] = 0.0;
	   else
	       Vh[i] = c1 * pow(h[i],range);
	 }
	 break;
	 
      default:
         fprintf(stderr, "WARNING: Unknown modelCode passed to vgram_vals() ");
	 return(1);
   
   }  /* end switch */
   
    return(0);
    
}   /* end vgram_vals2() */


/***************************************************************************/
int get_parm_range(int modelCode, struct VARIOGRAM *Vptr, double *p1, double *p2)
/* Returns min/max/incr for 2 parameters (p1=sill, p2=length) appropriate to
   the estimated variogram values in Vptr for use in fitting to a model.
   Memory for p1[3], p2[3]  must by allocated by calling function.
   Returns 0 if successful, or -1 if an error occurred.

*/   

{
   int ilag, ipwr;
   double pmin, pmax, x;
   
 /* Range of sill vals determined in same way for all models */  
 
    pmin = 999999999;
    pmax = 0;
    for (ilag = 0; ilag < Vptr->nlags; ++ilag) {
       if (Vptr->Nh[ilag] > 0) {
          x = Vptr->Vh[ilag];
             if (x < pmin)
                pmin = x;
             if (x > pmax)
                pmax = x;
       }
    }
  
        p1[0] = pmin;          /* set sill values: min, max, incr */
        p1[1] = 1.5 * pmax;
        ipwr = NINT(log10(p1[1]));
        p1[2] = pow(10, ipwr-2);
	
	if (pmin == 0)   /* don't allow zero values in parm range*/
	   p1[0] += p1[2];
   
   /* p2 differs for power model */
   
   if (modelCode == 3) {  /* power model  */
       p2[0] = 1;   /* power is set between 1-2 */
       p2[1] = 2;   /* max */
       p2[2] = 0.1;  /* incr */
       return(0);
   }
   
   
    pmin = 999999999;
    pmax = 0;
    for (ilag = 1; ilag < Vptr->nlags; ++ilag) {
       if (Vptr->Nh[ilag] > 0) {
           x = Vptr->h[ilag];
             if (x < pmin)
                pmin = x;
             if (x > pmax)
                pmax = x;
       }
    }
    
    /* limit grid search to 400 km, and incr to 10 km */  
    p2[0] = 100;     /* arbitrary lower limit for decorrelation scale*/
    p2[1] = pmax < 400 ? pmax: 400;
    p2[2] = 10  ;
    
    if (pmin < p2[2])
       p2[0] = p2[2];
   
     return(0);
}   /* end get_parm_range() */
/***************************************************************************/

struct VARIOGRAM *vgram_estimate(double **e, int ngrids, float *lat_e, float *lon_e, int nrows, int ncols,double maxdist, double binsize)
/* Estimates the semi-variogram V(h) where h is the lag (distance in km).

   INPUT:  
     e:  arrays of residuals (bin3d - FG) on a gridded surface (row order)
     ngrids:  number of arrays (1st dimension of e)
     lat_e, lon_e : coordinate vectors for e
     nrows, ncols :  dimension of coordinate vectors
     maxdist : max distance between points to consider
     binsize : (km) for computing lags
   OUTPUT
      Vptr->nlags :  dimension of output variogram
          ->Vh  : array -- variogram values for each bin
	  ->h  : array -- average distance for each bin
	  ->Nh : array -- number of pts in each bin
	  
    NOTE:  e is either a value, masked, or missing 
    
    Added code to use Cressie and Hawkins version of computing RSS (residual sum of squares)
    but decided to stick with original based on Zimmerman and Zimmerman (1989) evaluation
    
*/

{
   int error, i,j, igrid, npts, row, col, ilag, nbins, sq, icntr;
   double mask_val, diff, hdg, dist;
   struct VARIOGRAM *Vptr;
   double *sse, *rms_sum, *dsum, *lat, *lon;
   int *count;
   
   npts = ncols * nrows;
   mask_val = TOOBIG * 0.1;
   nbins = (int) ceil(maxdist / binsize);

/* allocate memory */
   
    sse = (double *) calloc( nbins, sizeof(double));
/*    rms_sum = (double *) calloc(nbins, sizeof(double)); */
    dsum = (double *) calloc(nbins, sizeof(double));
    count  = (int *) calloc(nbins, sizeof(int));
    lat = (double *) calloc(npts, sizeof(double));
    lon = (double *) calloc(npts, sizeof(double));

 /* create lat/lon arrays for each point in grid */
     
    sq = 0;
    for (row = 0; row < nrows; ++row) {
       for (col = 0; col < ncols; ++ col) {
          lat[sq] = (double) lat_e[row];
	  lon[sq] = (double) lon_e[col];
	  ++sq;
       }
    }

   icntr = npts / 2;   /* index to central node */
   
/* Evaluate all pairs of points */
   for (igrid = 0; igrid < ngrids; ++igrid) {
      if (e[igrid][icntr] < mask_val) {
         for (i = 0; i < npts; ++i) {
            if (ABS(e[igrid][i]) < mask_val) {
         
	       for (j = i; j < npts; ++j) {
	          if ( j != i && ABS(e[igrid][j]) < mask_val) {
	    
	          diff = ABS(e[igrid][j] - e[igrid][i]);
	          dist = distance_c(lat[j], lon[j], lat[i],lon[i], 1, &hdg);

	             if (dist < maxdist) {
	             ilag = NINT(dist/binsize);
		        if (ilag >= 0 && ilag < nbins) {
		        sse[ilag] += diff*diff;
			/* rms_sum += sqrt(diff); */
		        dsum[ilag]+= dist;
		        ++count[ilag];
		         }  
	             }
	    
	          } /* end if */
	       }  /* end for j*/
            } /* end if */
         } /* end for i */
      } /* end if e[igrid] */ 
   } /* end for igrid */

   free(lat);
   free(lon);
   
   for (ilag = 0; ilag < nbins; ++ilag) {
      if (count[ilag] > 0) {
         sse[ilag] = sse[ilag]  / (count[ilag] + count[ilag]); 
/*	 rms_sum[ilag] = pow(rms_sum[ilag]/count[ilag],4) / (0.914 + 0.988/count[ilag]); */
	 dsum[ilag] /= count[ilag];
      }
   }
   
   Vptr = (struct VARIOGRAM *) get_memory(NULL, 1, sizeof(struct VARIOGRAM ));
   Vptr->nlags = nbins;
   Vptr->Vh = sse; 
   Vptr->h = dsum;
   Vptr->Nh = count;
   return(Vptr);
}   /* end vgram_estimate() */
/***************************************************************************/

int vgram_check(struct VARIOGRAM *vgram_ptr, int nbins, int minpts, int maxlag)
/* Determines whether the variogram is valid or not based on distribution
   of observations in each lag bin <= maxlag:  allows up to nbins to have less 
   than minpts.   Returns 1 if valid, 0 if not. */
   
{
   int i, nbad;
   
   nbad = 0;
   for (i = 0; i < vgram_ptr->nlags; ++i) {
      if (vgram_ptr->h[i] > maxlag)
           return(1);
	   
      if (vgram_ptr->Nh[i] < minpts)
           ++nbad;
	   
      if (nbad > nbins)
         return(0);
	 
    }
    return(1);

}  /* end vgram_check() */
/***************************************************************************/
