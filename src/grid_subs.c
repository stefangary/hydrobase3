/* grid_subs.c
................................................................................
                          *******  HydroBase3  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2010
................................................................................
................................................................................
 * Modules in this file:
 *
 *   ij2xy()   determines x,y (lon,lat) values corresponding to index i (col),j (row) based on GRID_INFO
 *   xy2ij()   determines index (col,row) values corresponding to x(lon),y(lat) values based on GRID_INFO

 *  ij2xy_nochk
 *  xy2ij_nochk()
 
*  sq2rc()  returns row, col associated with sq (1-dimensional array index)
*  sq2latlon()  returns lat, lon associated with sq
* distance_c() determines cartesian distance and direction between 2 points
*  pointB() determines lat/lon of pointB given starting point, bearing and distance
 *  inside ()     Determines whether gridnode sits within polygon.
 *  get_mask()    Masks/unmasks gridnodes based on one or more polygons.
 *  zero_weights_2d()  sets weights to zero based on distribution of empty, masked nodes relative to center gridnode;
 *  set_to_zero() sets array elements to zero based on azimuth and distance criteria.
 *  set_node() sets array elements to value based on azimuth and distance criteria.

 *  nn_interp2d()    Weighted nearest neighbor interpolation algorithm      
 *  weighted_mean() returns a weighted average for an m x n matrix, weights supplied by calling routine
 *  weighted_variance() returns a weighted variance for an m x n matrix,and mean value, weights supplied by calling routine
*
   void get_weights_c(double *weight, struct GRID_INFO *hptr, double alpha, double xdist, double ydist, int icntr, int jcntr)
   sets weight based on cartesian distance of grid node from center node (located at icntr,jcntr)
*  
    void get_weights_g(double *weight, struct GRID_INFO *hptr, double alpha, int xcntr, int ycntr)
  /* sets weight based on gridpoint distance between node and center node: 
 ................................................................................
*/ 
 
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <float.h>
#include <stddef.h>
#include "hb_grids.h"
#include "hb_memory.h"


/****************************************************************************/
int ij2xy(struct GRID_INFO *hptr, int i, int j, double *xptr, double *yptr)
  /* translates index i,j (col,row) into x,y (lon,lat) values based on GRID_INFO.
     Returns 0 for successful conversion, -1 if the point is out of the gridbounds*/
{
  if ((j >= hptr->ny)  || (j < 0)  || (i >= hptr->nx)  || (i < 0) )
     return(-1);
     
  if (hptr->node_offset == 0) {  /* for node grids */
      *xptr = hptr->x_min + i * hptr->x_inc;
      *yptr = hptr->y_min + j * hptr->y_inc;
      return (0);
  }
    
  /* for pixel grids */
  *xptr = hptr->x_min + (i + 0.5) * hptr->x_inc;
   *yptr = hptr->y_min + (j + 0.5) * hptr->y_inc;
  return (0);
}  /* end ij2xy */
/****************************************************************************/
void ij2xy_nochk(struct GRID_INFO *hptr, int i, int j, double *xptr, double *yptr)
  /* translates index i,j (col,row) into x,y (lon,lat) values based on GRID_INFO.  Does not check boundaries.*/
{
     
  if (hptr->node_offset == 0) {  /* for node grids */
      *xptr = hptr->x_min + i * hptr->x_inc;
      *yptr = hptr->y_min + j * hptr->y_inc;
      return;
  }
    
  /* for pixel grids */
  *xptr = hptr->x_min + (i + 0.5) * hptr->x_inc;
  *yptr = hptr->y_min + (j + 0.5) * hptr->y_inc;

}  /* end ij2xy_nochk() */
/****************************************************************************/
int xy2ij(struct GRID_INFO *hptr, double x, double y, int *iptr, int *jptr)
  /* translates x,y values into index i,j based on GRID_INFO.
     Returns 0 for successful conversion, -1 if the point is out of the gridbounds*/
{
 /* the extra smidgeon puts points that fall right on a border into the proper box */
   x +=  0.0000001;
   y +=  0.0000001;
   
  if (hptr->node_offset == 0) {  /* for node grids */
      if (  (x >= hptr->x_max + (0.5 * hptr->x_inc)) || (x  < hptr->x_min - 0.5 * hptr->x_inc))
         return (-1);
      if (  (y  >= hptr->y_max + (0.5 * hptr->y_inc)) || (y  < hptr->y_min - 0.5 * hptr->y_inc))
         return (-1);
	 
      *iptr = NINT((x  - hptr->x_min) / hptr->x_inc); 
      *jptr = NINT((y  - hptr->y_min) / hptr->y_inc); 
      return (0);
  }
 
    
  /* for pixel grids */
  
  if ( (x >= hptr->x_max) || (x < hptr-> x_min))
         return (-1);
   if ( (y >= hptr->y_max) || (y < hptr-> y_min))
         return (-1);
 
      *iptr = NINT((x  - hptr->x_min) / hptr->x_inc - 0.5); 
      *jptr = NINT((y  - hptr->y_min) / hptr->y_inc - 0.5); 
  return (0);
}  /* end xy2ij */
/****************************************************************************/
void xy2ij_nochk(struct GRID_INFO *hptr, double x, double y, int *iptr, int *jptr)
  /* translates x,y values into index i,y based on GRID_INFO.
    Does not check grid bounds   */
{
 /* the extra smidgeon puts points that fall right on a border into the proper box */
   x +=  0.0000001;
   y +=  0.0000001;
   
  if (hptr->node_offset == 0) {  /* for node grids */
	 
      *iptr = NINT((x  - hptr->x_min) / hptr->x_inc); 
      *jptr = NINT((y  - hptr->y_min) / hptr->y_inc); 
      return ;
  }
 
  /* for pixel grids */
 
   *iptr = NINT((x  - hptr->x_min) / hptr->x_inc - 0.5); 
   *jptr = NINT((y  - hptr->y_min) / hptr->y_inc - 0.5); 
  return;
}  /* end xy2ij_nochk() */
/****************************************************************************/
int sq2rc(int sq, struct GRID_INFO *hptr, int *rowptr, int *colptr)
  /* translates index sq into row,col values based on GRID_INFO.
     Returns 0 for successful conversion, -1 if the point is out of the gridbounds*/
{
   *colptr = sq % hptr->nx;
   *rowptr = (sq - *colptr) / hptr->nx;
   
   if ((*rowptr < 0) || (*rowptr >= hptr->ny))
       return(-1);
       
   return(0);
     
}  /* end sq2rc */

/****************************************************************************/
int sq2latlon(int sq, struct GRID_INFO *hptr, double *latptr, double *lonptr)
  /* translates index sq into lat, lon values based on GRID_INFO.
     Returns 0 for successful conversion, -1 if the point is out of the gridbounds*/
{
   int row, col;
   col = sq % hptr->nx;
   row = (sq - col) / hptr->nx;
   
   if ((row < 0) || (row >= hptr->ny))
       return(-1);
       
   if ((col < 0) || (col >= hptr->nx))
       return(-1);
       
     if (hptr->node_offset) {   /* pixel grids */
       *lonptr = hptr->x_min + (col + 0.5) * hptr->x_inc;
       *latptr = hptr->y_min + (row + 0.5) * hptr->y_inc;
       return(0);
    }  
    
    /* node grids */
    *lonptr = hptr->x_min + col  * hptr->x_inc;
    *latptr = hptr->y_min + row  * hptr->y_inc;
       
    return(0);   
     
}  /* end sq2latlon */
/****************************************************************************/
double distance_c(double lat0, double lon0, double lat1, double lon1, int kilo, double *dir_addr)
/* Returns the cartesian distance between 2 points using haversine (great circle) 
  formula and direction in degrees from north. Distance is returned 
  in kilometers (kilo!=0) or nautical miles (kilo == 0).
*/

{
   double dlat, dlon, phi, a, c, R;
   
   R = EarthRadius;
   if (kilo)
      R = EarthRadiusKm;
      
  /* convert decimal degrees to radians */
  
   lat0 *= RADperDEG;
   lon0 *= RADperDEG;
   lat1 *= RADperDEG;
   lon1 *= RADperDEG;
   
   dlat = lat1 - lat0;
   dlon = lon1 - lon0;
   
   
   if (dlat == 0) {  /* course is zonal */
      
      *dir_addr = 90;
      if (dlon < 0)
         *dir_addr = 270;
      
      return(abs(dlon * cos(0.5*(lat0+lat1))* R));
   }
   
   if (dlon == 0 ) {  /* meridional */
   
      *dir_addr = 0.0;
       if (dlat < 0)
         *dir_addr = 180.0;
	 
      return(abs(dlat*R));
   }
   
   phi = atan(dlon/dlat);
   
   if (phi > 0) {
      *dir_addr = phi * DEGperRAD;  /* 0 - 90 */
      
      if (dlon < 0)
         *dir_addr += 180.0;     /* 180 - 270 */
   }
   else {
   
       *dir_addr = 180 + phi * DEGperRAD;   /* 90 - 180 */
       
       if (dlon < 0)
         *dir_addr = 360 + phi * DEGperRAD;   /* 270 - 360 */
   }   
   
   a =  pow(sin(dlat * 0.5),2) + cos(lat0) * cos(lat1) * pow(sin(dlon *0.5),2);
   c = 2.0 * atan2(sqrt(a), sqrt(1-a));
   
   return(R*c);

} /* end distance_c() */  
/****************************************************************************/
double distance_g(struct GRID_INFO *ginfo, int sq1, int sq2, double *phi_addr)
  /* Returns the gridpoint distance between 2 gridnodes. 
     The direction (radians) is returned at phi_addr*/
{

}
/****************************************************************************/
 int pointB(double latA, double lonA, double bearing, double dist, int kilo, double *latB, double *lonB)
 /* Computes position of pointB given starting point, bearing and distance.
    If kilo is set (non-zero), distance is kilometers, otherwise nautical miles.
    Returns 1 if the pole was crossed, or 0 otherwise. */
{
   int cross;
   double R, rr, dx, dy, phi, distrad;
   double rlat, rlon;
   
   cross = 0;
   if (latA == 90.)
      latA = 89.9999;
   if (latA == -90.)
      latA = -89.9999;
      
   rlat = latA * RADperDEG;
   rlon = lonA * RADperDEG;
   phi = bearing * RADperDEG;
   R = EarthRadius;
   if (kilo) 
       R = EarthRadiusKm;
   
       
    dx = abs(dist * sin(phi));
    dy = abs(dist * cos(phi));
    if (bearing > 180)
       dx = -dx;
       
    if (bearing > 90 && bearing < 270 )
       dy = -dy;
     
     rr = RADperDEG * R;
       
    *latB = latA + dy /  rr;
    *lonB = lonA + dx / (rr * cos(rlat)); 
      
      if (*latB > 90) {    /* crosses north pole */
         dy = *latB - 90;
	 *latB = 90 - dy;
	 *lonB += 180;
	 cross = 1;
      }
       if (*latB < -90) {    /* crosses south pole */
         dy = *latB + 90;
	 *latB = -90 - dy;
	 *lonB += 180;
	 cross = 1;
      }
     
      if (*lonB > 180)
         *lonB -= 360;
	 
      if (*lonB < -180)
         *lonB += 360;
	 
     return(cross);
} /* end pointB() */

/***************************************************/
int pointB_hav(double latA, double lonA, double bearing, double dist, int kilo, double *latB, double *lonB)
 /* Computes position of pointB given starting point, bearing and distance using Haversine formula (great circle distance).
    If kilo is set (non-zero), distance is kilometers, otherwise nautical miles.
    Returns  0 for normal */
{
   double R, rr, phi, dist_rad;
   double rlat, rlon;
   
   rlat = latA * RADperDEG;
   rlon = lonA * RADperDEG;
   phi = bearing * RADperDEG;
   R = EarthRadius;
   if (kilo) 
       R = EarthRadiusKm;
      
   dist_rad = dist / R;
       
   *latB = asin(sin(rlat) * cos(dist_rad) + cos(rlat) * sin(dist_rad) * cos(phi)); 
   
   *lonB =  DEGperRAD * (rlon + atan2(sin(phi)*sin(dist_rad) * cos(rlat), cos(dist_rad) - sin(rlat) * sin(*latB))); 
   
   if (*lonB > 180)
      *lonB -= 360;
     
   if (*lonB < -180)
         *lonB += 360;
	 
   *latB *= DEGperRAD;      
   return(0);
} /* end pointB_hav() */

int inside (double x, double y,  double *xb,  double *yb, int nb) 
/* ***********************************************************
/* 
/*    Given a point x,y and the series xb(k),yb(k) (k=1...nb) defining 
/*    vertices of a closed polygon.  Ind is set to 1 if the point is in 
/*    the polygon and 0 if outside.  Each time a new set of bound points 
/*    is introduced ind should be set to 999 on input. 
/*    It is best to do a series of y for a single fixed x. 
/*    Method ... a count is made of the no. of times the boundary cuts 
/*    the meridian thru (x,y) south of (x,y).   An odd count indicates 
/*    the point is inside , even indicates outside. 
/*    See a long way from euclid by Constance Reid  p 174 . 
/*    Oceanography emr   oct/69 
/*
/*    As should be readily apparent, this C language version is a direct
/*    conversion from the original PlotPlus FORTRAN language version.
/*    Any other changes were avoided for the sake of checkout and
/*    compatibility.  That does not imply an approval for efficiency
/*    but rather expediency.
/*       Version 1.0          R. Goldsmith          May, 1999
************************************************************* */
{
int    ind;
int    k, ke, kw;
int    kp1;
int    nc;

double slope;
double xprev;
double yc[NCROSS];

if (nb <= 0)  {
  ind = 1;
  return ind;
}
  
xprev = x;
nc = 0;
for (k=0; k < nb; k++)  {
  if (k == nb - 1)  
        kp1 = 0;
  else  
        kp1 = k + 1;
  kw = k; 
  if (xb[k] == xb[kp1])  continue;
  if (xb[k] > xb[kp1])   kw = kp1;
  ke = k + kp1 - kw;
  if (x > xb[ke])  continue;
  if (x < xb[ke]  &&  x <= xb[kw])  continue;
  ++nc;
  if (nc == NCROSS)  {
    fprintf (stderr, "function inside(): internal crossings array overflow.\n");
    ind = 0;
    exit (-1);
  }
  slope = (yb[ke] - yb[kw])/(xb[ke] - xb[kw]);
  yc[nc-1] = yb[kw] + (x - xb[kw])*slope; 
}

ind = 0; 
if (nc <= 0)  return (ind);
for (k=0; k<nc; k++)  {
  if (yc[k] < y)  ind = 1 - ind;
}
return (ind);
}   /* end inside() */



/***********************************************************/
int get_mask(FILE *maskfile, struct GRID_INFO *hptr, char *mask, BOOLEAN latfirst)

   /* Reads multi-segment {lon/lat} (or general [x,y]) values and 
      sets appropriate grid nodes in array mask to 1 or 0.
      
      Multiple segments are separated by a line starting with the '>' character
      >O means mask outside the polygon specified by following points
      >I means mask inside the polygon 
      >C (OR '>' alone OR no '>' lines) means mask is specified as 
         individual cells.
      
      This method is being used to replicate the original zgrid procedure.  
      It is tortuous to implement but let's get the original working
      before I slip in some more efficient code.  Hopefully
      the nodes and polygons are few or the machine is fast.  In the 
      future one might look at the gmt/grdmask code and the
      use of the non-zero_winding/out_edge_in routines to do the
      polygon check.  
      
      The array order is generated to coordinate with the zgrid function 
      (it stores elements row-by-row  with first point at xmin,ymin:
           sq = row * ncols  + col).  grid dimensions are passed in the
	   struct GRID_INFO hptr.  

 */
{
int nrows, ncols, ix, iy;
int n_alloc, n, nread;
int n_mask, row, col, sq;
int cell_zval;

BOOLEAN cell_mode;
BOOLEAN mask_outside;
BOOLEAN poly_inside;
BOOLEAN multi_segments;
BOOLEAN pt_okay;

double  DX, DY, X1, Y1;
double *xmask, *ymask;
double xnode, ynode;
double xy[2];
double xnode_offset, ynode_offset;

char line[BUFSIZ];
char *moremask;
char EOL_flag = '>';

  xnode_offset = hptr->node_offset ? .5 * hptr->x_inc : 0.0;
  ynode_offset = hptr->node_offset ? .5 * hptr->y_inc : 0.0;
  
  ncols  = hptr->nx;
  nrows  = hptr->ny;
  X1 = hptr->x_min + xnode_offset;
  Y1 = hptr->y_min + ynode_offset;
  DX = hptr->x_inc;
  DY = hptr->y_inc;


  ix = 0;
  if (latfirst) ix = 1;
  iy = 1 - ix;

  cell_mode = TRUE;
  mask_outside = FALSE;
  n_mask = 0;
  
/* If first character of first line is '>', then it is a multi-segment
   polygon mask */
 
  moremask = fgets (line, BUFSIZ, maskfile);
  multi_segments = line[0] == '>';
  
  if (multi_segments) {
  
    /* default masking is by cell_mode (point by point) */

     if (line[1] == 'I'  ||  line[1] == 'i')  {
       cell_mode = FALSE;
       mask_outside = FALSE;
     }
     if (line[1] == 'O'  ||  line[1] == 'o')  {
        cell_mode = FALSE;
        mask_outside = TRUE;
     }
     moremask = fgets (line, BUFSIZ, maskfile);
  }
             
  while (moremask) {
    n_alloc = MEM_CHUNK;
    xmask = (double *) get_memory ((void *)NULL, n_alloc, sizeof(double));
    ymask = (double *) get_memory ((void *)NULL, n_alloc, sizeof(double));

    n = 0;
    while ((moremask && !multi_segments) || (moremask && multi_segments && line[0] != EOL_flag)) {

        if (line[0] == EOL_flag) {
          moremask = fgets (line, BUFSIZ, maskfile);
          continue;
        }
        
          /* Cell mode assumes points are at nodes. */
	  
        if (cell_mode)  {
	
	  cell_zval = 1;
          nread = sscanf (line, "%lf %lf %d", &xy[ix], &xy[iy], &cell_zval);
	  
          if (nread < 2)  {
            fprintf (stderr, "get_mask(): < 2 fields in masking file cell spec.\n");
            exit(1);
          }
	  
	  /* check for out-of-bounds points */
	  
          col  = (int)(NINT((xy[0] - X1) / DX ));
          row  = (int)(NINT((xy[1] - Y1) / DY ));
	  
	  pt_okay = TRUE;
	  if (col < 0 || col >= hptr->nx) pt_okay = FALSE;
	  if (row < 0 || row >= hptr->ny) pt_okay = FALSE;   
	  
	  if (pt_okay) {
	     sq = col + ncols * row;   
	       
             if (cell_zval)  {   /* mask it */
		 if (!mask[sq])  ++n_mask;
	         mask[sq] = '1';
             }  
             else  {   /* Or unmask it -- set to 0 */
               if (mask[sq] != 0)  --n_mask;
               mask[sq] = 0;  
             }  
           } /* end if pt_okay */
       }  /* end if cell node. */
	  
       else  {     /* Not cell_mode. */
          sscanf (line, "%lf %lf", &xy[ix], &xy[iy]);
          xmask[n] = xy[0];
          ymask[n] = xy[1];
          n++;
          if (n == n_alloc) {
            n_alloc += MEM_CHUNK;
            xmask = (double *) get_memory ((char *)xmask, n_alloc, sizeof (double));
            ymask = (double *) get_memory ((char *)ymask, n_alloc, sizeof (double));
          }
        }  /* end if not cell_mode. */
	
        moremask = fgets(line, BUFSIZ, maskfile);
    }  /* end while same segments. */

    if (!cell_mode)  {
        xmask = (double *) get_memory ((char *)xmask, n, sizeof (double));
        ymask = (double *) get_memory ((char *)ymask, n, sizeof (double));
        n_alloc = n;

       /* visit each gridsquare */
       sq = -1;
       for (row = 0; row < nrows; row++)  {
         ynode = Y1 + (double) row * DY;
	
         for (col = 0; col < ncols; col++)  {
	   ++sq;
	   
           if (mask[sq])  continue;  /* already masked? */
	  
           xnode = X1 + (double) col * DX;
           poly_inside = inside (xnode, ynode, xmask, ymask, n);
	   
           if (mask_outside)  {
              mask[sq] = (char)(1 - poly_inside);
              if (poly_inside == 0)  {
                n_mask++;
              }
           }
           else  {
             mask[sq] =  (char) (poly_inside);
             if (poly_inside == 1)  {
                n_mask++;
             }
           }

          }  /* end  for col  */
        }  /* end  for row  */
     }  /* end if not cell_mode. */
      
     free ((char *) xmask);
     free ((char *) ymask);
      
      
     if (line[0] == EOL_flag) {
        cell_mode = TRUE;
        mask_outside = FALSE;
        if (line[1] == 'I'  ||  line[1] == 'i') 
	  cell_mode = FALSE;
	
	if (line[1] == 'O'  ||  line[1] == 'o') {
	  cell_mode = FALSE;
	  mask_outside = TRUE;
	}
     }
      
     if (moremask) moremask = fgets (line, BUFSIZ, maskfile);
 }  /* end while moremask segments. */
      
 return(n_mask);

}  /* end get_mask() */

/****************************************************************************/
void zero_weights_2d(double *x, double *wghts,double empty_val, double mask_val, struct GRID_INFO *hptr, int xmid, int ymid)
/*  Zeros weights corresponding to empty and masked nodes in array x. 
    Additionally, when a masked node is encountered, all weights further distant along the azimuth
    from the center of the 2D array are set to zero.  

    x:         input array representing 2D matrix, parameters defined in hptr.
    wghts:   input array of weights (preset in calling routine) with same dimensions as x
    empty_val:
    mask_val:  values corresponding to flags in x
    hptr:      info about the grid parameters
    xmid, ymid:  index to column and row at center of 2D array
*/
{
   int row, col, mask_rest, sq, nsq;
   int row2, col2, sq2, rad;
   int dx, dy;
   int *mask;
   double   toosmall, toobig;
   double phi_min, phi_max, dist_min;  
   double *phi, *dist;
   
   toobig = mask_val / 1.1;    /* generate values for testing */
   toosmall = empty_val / 1.1;
   nsq = hptr->nx * hptr->ny;
   
   /* create and populate arrays to store distance, azimuth, mask flag info */
   
   phi = (double *)calloc(nsq, sizeof(double));
   dist = (double *)calloc(nsq, sizeof(double));
   mask = (int *) calloc(nsq, sizeof(int));
   
   for (row = 0; row < hptr->ny; ++row) {
      for (col = 0; col < hptr->nx; ++col) {
         sq = row * hptr->nx + col;
	 dy = row - ymid;
	 dx = col - xmid;
	 dist[sq] = sqrt(dx*dx + dy*dy);
	 phi[sq] = atan2(dy, dx) * DEGperRAD;  /* -180 < phi < 180 */
	 mask[sq] = 1;   
         if (x[sq] < toosmall)   /*empty node */
            mask[sq] = -1;
         else if (x[sq] > toobig)  /* masked node */
            mask[sq] = 0;
      }
   } 
   
 /* at this point, every element of mask is set to either -1, 0 or 1 */
 
   
/* search west for masked nodes */
   row = ymid;
   col = xmid;
   dx = 0;
   mask_rest = 0;
   while (++dx <= xmid) {
      --col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
      
   } /* end while */
  
/*search east */
   row = ymid;
   col = xmid;
   dx = 0;
   mask_rest = 0;
   while (++dx <= xmid) {
      ++col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
      
   } /* end while */

/*search north */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      ++row;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
   } /* end while */
   
/*search south */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      --row;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
   } /* end while */

/*search southwest */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      --row;
      --col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
      
   } /* end while */

/*search northwest */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      ++row;
      --col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
  } /* end while */

 /*search northeast */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      ++row;
      ++col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
   } /* end while */
  
 /*search southeast */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      --row;
      ++col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 0;
      else if (mask[sq] == 0)  
	 mask_rest = 1;
   } /* end while */

  /* Check for consecutive masked nodes along a row or a column (a topographic barrier),
     and, where found, mask all nodes at greater distance between the 2 azimuths connecting the 
    center node and those masked nodes */
     
 /* first do quadrants east and west of center node*/ 
      
   for (rad=1; rad <= xmid; ++rad) {
     col = xmid + rad;
     col2 = xmid - rad;
     row = ymid + rad;
     row2 = ymid - rad;
     if (row >= hptr->ny) {
        row = hptr->ny - 1;
	row2 = 0;
     }
         /* work top to bottom */
     while ( row > row2) {   
           /* eastern quadrant */
        sq = row * hptr->nx + col;
        sq2 = (row-1) * hptr->nx + col;
        if ( (mask[sq]==0) && (mask[sq2]==0) ) {
          phi_max = phi[sq];
	  phi_min = phi[sq2];
	  dist_min = dist[sq] > dist[sq2] ? dist[sq2]: dist[sq];
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist_min);
        }
	 /* western quadrant */
        sq = row * hptr->nx + col2;
        sq2 = (row-1) * hptr->nx + col2;
        if ( (mask[sq]==0) && (mask[sq2]==0) ) {
          phi_max = phi[sq2];
	  phi_min = phi[sq];
	  
	  if (phi_max > 0 && phi_min < 0)  /* special case for phi = +/-180 */
	      phi_max = -180.0;
	      
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist_min);
        }
	
	--row;
     } /* end while */
     
   } /* end for rad */
   
   /* quadrants north and south of center */
   
   for (rad=1; rad <= ymid; ++rad) {
     col = xmid - rad;
     col2 = xmid + rad;
     row = ymid + rad;
     row2 = ymid - rad;
     if (col < 0) {
	col = 0;
        col2 = hptr->nx - 1;
     }
        /* work left to right */
     while ( col < col2) {
        /* north quadrant */
        sq = row * hptr->nx + col;
        sq2 = sq+1;
        if ( (mask[sq]==0) && (mask[sq2]==0) ) {
          phi_max = phi[sq];
	  phi_min = phi[sq2];
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist_min);
        }
        /* south quadrant */
        sq = row2 * hptr->nx + col;
        sq2 = sq + 1;
        if ( (mask[sq]==0) && (mask[sq2]==0) ) {
          phi_max = phi[sq2];
	  phi_min = phi[sq];
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist_min);	  
        }
	++col;
     } /* end while */
     
   } /* end for rad */
   
   
   /*  find 0s and -1s in mask array and zero the corresponding elements in wghts array */
   for (sq = 0; sq < nsq; ++sq) {
     if (mask[sq] <= 0)
        wghts[sq] = 0.0;
   }
   
   free(dist);
   free(phi);
   free (mask);
   return;

}  /* end zero_weights_2d() */
/****************************************************************************/
void set_to_zero(int *m, double *a, double *dist, int npts, double a_min, double a_max, double d_min)
   /* determines indices where elements of array a between a_min and a_max AND elements of dist > d_min
      and sets corresponding elements of m to zero */
{
   int sq;
   
   for (sq = 0; sq < npts; ++sq) {
      if ((a[sq] >= a_min) && (a[sq] <= a_max) && (dist[sq] >= d_min))
         m[sq] = 0;
   }
   return;

} /* end set_to_zero() */
/****************************************************************************/
void set_node(int *m, double *a, double *dist, int npts, double a_min, double a_max, double d_min, int val)
   /* determines indices where elements of array a between a_min and a_max AND elements of dist > d_min
      and sets corresponding elements of m to val  */
{
   int sq;
   
   for (sq = 0; sq < npts; ++sq) {
      if ((a[sq] >= a_min) && (a[sq] <= a_max) && (dist[sq] >= d_min))
         m[sq] = val;
   }
   return;

} /* end set_node() */
/****************************************************************************/
void nn_interp2d(double *xin, double *wghts, double empty_val, double mask_val, int xmid, int ymid, struct GRID_INFO *hptr, double *xoutptr, double *sdptr, double *avg_wgt_ptr, int *cptr)
   /*  Interpolates value in middle of array xin (at index xmid, ymid) using nearest neighbor algorithm.
       2D array has dimensions 
                ncols = xmid * 2 + 1 
		nrows = ymid * 2 + 1 
       Array element 0 corresponds to row0, col0;
       elements are stored in row order (column varying fastest). 

      Corresponding weights stored in similar grid (wghts) are applied for interpolation.
      Searches in all directions from central (missing) gridnode and uses 2 closest points 
      in each of 4 sectors to compute a weighted average.. 
      Missing values are denoted by empty_val (large negative), masked values (below seafloor) 
      are large positive values (mask_val).
      If mask_val is encountered, the search is suspended in that direction to prevent
      mixing watermasses across topographic barriers.  
      Returned values:
          xoutptr:   weighted mean
	  cptr:   number of gridnodes contributing to mean
	  sdptr:  std deviation of xvalues contributing to mean
	  avg_wgt_ptr:  average weight  (weightsum/number of gridnodes) 
   */
{
   int row, col, dy, dx, sq, npts, nobs;
   int rad, row2, col2, nsq;
   int *mask;
   double  wsum, xwsum, x2sum, xsum; 
   double phi_min, phi_max;
   double *phi, *dist;
   
   
   /* set weights to zero for missing nodes, masked nodes and topographic barriers */
   zero_weights_2d(xin, wghts, empty_val, mask_val, hptr, xmid, ymid);

      /* create and populate arrays to store distance, azimuth, mask flag info */
   nsq = hptr->nx * hptr->ny;   
   phi = (double *)calloc(nsq, sizeof(double));
   dist = (double *)calloc(nsq, sizeof(double));
   mask = (int *) calloc(nsq, sizeof(int));
   
   for (row = 0; row < hptr->ny; ++row) {
      for (col = 0; col < hptr->nx; ++col) {
         sq = row * hptr->nx + col;
	 dy = row - ymid;
	 dx = col - xmid;
	 dist[sq] = sqrt(dx*dx + dy*dy);
	 phi[sq] = atan2(dy, dx) * DEGperRAD;  /* -180 < phi < 180 */
	 mask[sq] = 1;   
         if (wghts[sq] == 0.0)  
            mask[sq] = 0;
      }
   } 

   xwsum = 0.0;
   x2sum = 0.0;
   xsum = 0.0;
   wsum = 0.0;
   nobs = 0;
  
 /* search each quadrant around center node*/ 
 
   /* eastern quadrant, work north to south */
   npts = 0;  
   rad = 1;  
   while ((rad <= xmid) && npts < 2) {
     col = xmid + rad;
     row = ymid + rad;
     row2 = ymid - rad;
     if (row >= hptr->ny) {
        row = hptr->ny - 1;
	row2 = 0;
     }
     
     while (( row > row2) && (npts < 2)) {   
        sq = row * hptr->nx + col;
        if ( mask[sq] > 0 ) {
          phi_max = phi[sq] > 0 ? phi[sq] + 1.0 : phi[sq]-1.0;
	  phi_min = phi[sq] > 0 ? phi[sq] - 1.0 : phi[sq]+1.0;
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist[sq]);
	   xwsum += xin[sq] * wghts[sq];
	   x2sum += xin[sq] * xin[sq] ;
	   xsum += xin[sq];
	   wsum += wghts[sq];
	   ++npts;
        }
	--row;
     } /* end while */
     
     ++rad;
   } /* end while rad */
   
	 /* western quadrant, work south to north */
   nobs += npts;
   npts = 0;
   rad = 1;  
  while ((rad <= xmid) && npts < 2) {
     col = xmid - rad;
     row = ymid - rad;
     row2 = ymid + rad;
     if (row2 >= hptr->ny) {
        row2 = hptr->ny - 1;
	row = 0;
     }
     while ((row < row2) && (npts < 2)) {   
        sq = row * hptr->nx + col;
        if ( mask[sq] > 0 ) {
          phi_max = phi[sq] > 0 ? phi[sq] + 1.0 : phi[sq]-1.0;
	  phi_min = phi[sq] > 0 ? phi[sq] - 1.0 : phi[sq]+1.0;
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist[sq]);
	   xwsum += xin[sq] * wghts[sq];
	   x2sum += xin[sq] * xin[sq] ;
	   xsum += xin[sq];
	   wsum += wghts[sq];
	   ++npts;
        }
	
	++row;
     } /* end while */
     ++rad;
   } /* end for rad */
   

   /* northern quadrant, search west to east*/
   nobs += npts;
   npts = 0;
   rad = 1;
   while ( (rad <= ymid) && npts < 2) {
     col = xmid - rad;
     col2 = xmid + rad;
     row = ymid + rad;
     if (col < 0) {
	col = 0;
        col2 = hptr->nx - 1;
     }
     while ( (col < col2) && (npts < 2)) {
        sq = row * hptr->nx + col;
        if ( mask[sq] > 0 ) {
          phi_max = phi[sq] > 0 ? phi[sq] + 1.0 : phi[sq]-1.0;
	  phi_min = phi[sq] > 0 ? phi[sq] - 1.0 : phi[sq]+1.0;
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist[sq]);
	   xwsum += xin[sq] * wghts[sq];
	   x2sum += xin[sq] * xin[sq] ;
	   xsum += xin[sq];
	   wsum += wghts[sq];
	   ++npts;
        }
	
	++col;
     } /* end while */
     ++rad;
   } /* end while rad */
   

   /* southern quadrant, search east to west*/
   nobs += npts;
   npts = 0;
   rad = 1;
   while ( (rad <= ymid) && npts < 2) {
     col2 = xmid - rad;
     col = xmid + rad;
     row = ymid - rad;
     if (col2 < 0) {
	col2 = 0;
        col = hptr->nx - 1;
     }
     while ( (col > col2) && (npts < 2)) {
        sq = row * hptr->nx + col;
        if ( mask[sq] > 0 ) {
          phi_max = phi[sq] > 0 ? phi[sq] + 1.0 : phi[sq]-1.0;
	  phi_min = phi[sq] > 0 ? phi[sq] - 1.0 : phi[sq]+1.0;
	  set_to_zero(mask, phi, dist, nsq, phi_min, phi_max, dist[sq]);
	   xwsum += xin[sq] * wghts[sq];
	   x2sum += xin[sq] * xin[sq] ;
	   xsum += xin[sq];
	   wsum += wghts[sq];
	   ++npts;
        }
	--col;
     } /* end while */
     ++rad;
   } /* end while rad */
 
    free(phi);
    free(dist);
    free(mask);
   
   nobs += npts;
   if (nobs == 0) {
      *cptr = 0;
      *xoutptr = empty_val;
      *avg_wgt_ptr = empty_val;
      *sdptr = empty_val;
      return;   
   }
   
   *cptr = nobs;
   *xoutptr = xwsum / wsum;
   *avg_wgt_ptr = wsum / nobs;
   *sdptr = 0.0;
   if (nobs > 1) {
      *sdptr = ABS( (x2sum -  *xoutptr * *xoutptr *nobs ) / (nobs-1) );
      *sdptr = sqrt(*sdptr);
    }
    
   return;
   
}  /* end nn_interp2d() */
/****************************************************************************/
double weighted_mean(double *xin, double *weights, double empty_val, double mask_val, int npts)
/* returns the weighted average of all values in xin.  Nodes set to empty-val or mask_val are
   assigned zero-weight by this routine. */
{
   int i, nobs;
   double wsum, xwsum, toosmall, toobig;

   xwsum = 0.0;
   wsum = 0.0;
   nobs = 0;
   toobig = mask_val / 1.1;    /* for testing */
   toosmall = empty_val / 1.1;
   
   
   for (i = 0; i < npts; ++i ) {
     if (weights[i] > 0) {
         if (xin[i] > toosmall)  {  /* either a real value or a mask */
            if (xin[i] < toobig) { 
              xwsum += xin[i] * weights[i];
	      wsum += weights[i];
	      ++nobs;
	    }
         }
      }  
   }
   
   if (nobs == 0) 
   return (empty_val);   
   
   return (xwsum / wsum);    

} /*end weighted_mean() */

/****************************************************************************/
double weighted_variance(double *xin, double *weights, double xbar, double empty_val, double mask_val, int npts)
/* returns the weighted variance of xin with mean value, xbar. 
   Weights are supplied by the calling routine.  Nodes set to 
   empty-val or mask_val are assigned zero-weight by this routine. */
{
   int i, nobs;
   double wsum, delta, dprop, dpropsq, toosmall, toobig;

   dprop = 0.0;
   dpropsq = 0.0;
   wsum = 0.0;
   nobs = 0;
   toobig = mask_val / 1.1;    /* for testing */
   toosmall = empty_val / 1.1;
   
   
   for (i = 0; i < npts; ++i ) {
     if (weights[i] > 0) {
         if (xin[i] > toosmall)  {  /* either a real value or a mask */
            if (xin[i] < toobig) { 
	      delta = xin[i] - xbar;
              dprop += delta * weights[i];
	      dpropsq += delta * delta *  weights[i];
	      wsum += weights[i];
	      ++nobs;
	    }
         }
      }  
   }
   
   if (nobs < 3) 
   return (empty_val);   
   
   dprop /= wsum;
   return (dpropsq / wsum - dprop * dprop) * nobs / (nobs-1);

} /*end weighted_variance() */
/**********************************************************/
void get_weights_c(double *weight, struct GRID_INFO *hptr, double alpha, double xdist, double ydist, int icntr, int jcntr)
  /* sets weight based on cartesian distance of grid node from center node (located at icntr,jcntr)
  
        weight[n] = e^[- alpha * dist^2 ].  
	where  alpha = (pi/halfwidth)^2 /4.5 unless a non-zero value is passed as an argument.
	
	weight = matrix of weights,stored as vector (row order)
	hptr:  struct containing info about matrix
	alpha: weighting factor for computing Gaussian weights 
	xdist: x radius (km) of search ellipse
	ydist: y-radius (km) of search ellipse
	icntr, jcntr : index of central node
*/
{

int ix, iy, sq, error;
double dist, halfwidth, dx, dy, dx2, dy2;
double xlen2, ylen2, lon_c, lat_c, dlat, dlon;
double *xvec, *yvec;

   xvec = (double *) calloc(hptr->nx, sizeof(double));
   yvec = (double *) calloc(hptr->ny, sizeof(double));
   
   for (ix = 0; ix < hptr->nx; ++ix) {
     error = ij2xy(hptr, ix, 0, &dlon, &dlat); 
     if (error) {
       fprintf(stderr,"FATAL ERROR computing weights in work grid\n");
       exit(1);
     }
     xvec[ix] = dlon;
   }
   for (iy = 0; iy < hptr->ny; ++iy) {
     error = ij2xy(hptr, 0, iy, &dlon, &dlat); 
     if (error) {
       fprintf(stderr,"FATAL ERROR computing weights in get_weights_c()\n");
       exit(1);
     }
     yvec[iy] = dlat;
   }
   
   
   xlen2 = xdist * xdist;
   ylen2 = ydist * ydist;
   lat_c = yvec[jcntr];
   lon_c = xvec[icntr];
   
   
   if (alpha == 0.0) {   
      halfwidth = xdist > ydist ? xdist: ydist;
      alpha = PI / halfwidth;
      alpha  *= alpha;
      alpha /= 4.5;
   }
   
   for (iy = 0; iy < hptr->ny; ++iy) {
      dlat = abs(yvec[iy] - lat_c);
      dy = dlat * RADperDEG * EarthRadius * KMperNM;
      dy2 = dy * dy;   /* in km^2 */
      
      for (ix = 0; ix < hptr->nx; ++ix) {
         sq = iy * hptr->nx + ix;
         dlon = abs(xvec[ix] - lon_c);
	 dx = dlon * RADperDEG * cos(RADperDEG *lat_c) * EarthRadius * KMperNM;
	 dx2 = dx * dx;
	 if ((dx2 / xlen2 + dy2 / ylen2) > 1.0) { /* outside search ellipse? */
	    weight[sq] = 0.0;
	 }
	 else {
	    dist = sqrt(dx2 +dy2);
	    weight[sq] = exp( -alpha * dist * dist);
         }
       }
   }
   free(xvec);
   free(yvec);
   return;
} /* end get_weights_c() */
/**********************************************************/
void get_weights_g(double *weight, struct GRID_INFO *hptr, double alpha, int xcntr, int ycntr)
  /* sets weight based on gridpoint distance between node and center node: 
  
        weight[n] = e^[- alpha * dist^2 ].  
	where  alpha = (pi/halfwidth)^2 /4.5 unless a non-zero value is passed as an argument.
	
	weight: matrix of weights,stored as vector (row order)
	hptr:  struct containing info about matrix
	alpha: alternate factor for computing Gaussian weights 
 xcntr, ycntr: index (i,j) to central node and also x- and y-radii of the search ellipse.
*/
{
int ix, iy, sq, dx, dy, x2, y2;
double  dist, halfwidth;

   x2 = xcntr * xcntr;
   y2 = ycntr * ycntr;
   if (alpha == 0.0) {   
      halfwidth = xcntr > ycntr ? xcntr: ycntr;
      alpha = PI / halfwidth;
      alpha  *= alpha;
      alpha /= 4.5;
   }
   
   for (ix = 0; ix < hptr->nx; ++ix) {
      for (iy = 0; iy < hptr->ny; ++iy) {
          sq = iy * hptr->nx + ix;
	  dy = iy-ycntr;
	  dx = ix-xcntr;
	  
	  if ((dx * dx / x2 + dy *dy/y2) > 1.0)  /* outside search ellipse? */
	       weight[sq] = 0.0;
	  else {
	     dist = sqrt(dy*dy + dx*dx);
	     weight[sq] = exp( -alpha * dist * dist);
	  }
      }
   }

   return;
}/* end get_weights_g() */
/**********************************************************/

