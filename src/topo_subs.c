/*  topo_subs.c


 ................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             March 2010
................................................................................

 functions for reading HydroBase3 topography files (COARDS-compliant netcdf format)
*/

#include <stdio.h> 
#include <stdlib.h>
#include "netcdf.h"
#include "hb_grids.h"

/* 
short *hb_get_topo(char *fname, struct GRID_INFO *hptr, double **latvec_addr, double **lonvec_addr, int depths_neg, int lon0to360, short *missing)

short find_nearest_topo_val(double lat, double lon, short *topo, struct GRID_INFO *gptr);
Returns bottom depth for lat/lon specified 

*/


short *hb_get_topo(char *fname, struct GRID_INFO *hptr, double **latvec_addr, double **lonvec_addr, int depths_neg, int lon0to360, short *missing)

/*  reads topography data from specified file (netcdf format) and returns a pointer to an 
     array of seafloor depths.  The geographic range of output values 
     is specified upon input in hptr (fields: x_min,x_max,y_min, y_max). Other grid params 
     are set to match the input grid. Depths below sea level will be positive if depths_neg 
     is 0 (FALSE), or negative if set to a non-zero value. Longitudes will be in the range
      -180 to +180 if lon0to360 is 0.  Arrays of latitude and longitude are created 
     (memory is allocated in this function) and returned.  The calling function needs to 
     pass the ADDRESS of a double * variable for each:
      example:   
	 double *latvector, *lonvector;
	 struct GRID_INFO *ginfo;
	 seafloor = hb_get_topo("/bathpath/etopo1.grd", ginfo, &latvector, &lonvector, 0, 0)
	 
     Missing value is returned as the last argument. The function returns NULL if an error
     occurred, and writes an error description to stderr device */

{
   int ncid, error;
   int lat_dim, lon_dim;
   int start_row, end_row, start_col, end_col;
   int strow_out, endrow_out, stcol_out, endcol_out;
   int indx, i, j, irow, icol;
   int split_region;
   char varname[MAX_NC_NAME];
   int  varid, nsqout, nsqin, z_is_neg;
   size_t len;
   double dval2[2];
   double lon1min, lon1max, lon2min, lon2max;
   double *latptr, *lonptr, *lonout, *latout;
   int lval;
   short *zin, *seafloor, no_val;
   struct GRID_INFO z_info;
   nc_type vtype;


/*  Open file */
   
   error = nc_open(fname, NC_NOWRITE, &ncid);

   if (error == NC_NOERR)  
      fprintf(stderr,"\nOpened %s ...\n", fname);
   else {
      fprintf(stderr,"\n Unable to open %s.\n", fname);
      fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
      return(NULL);
   }

/*  Read grid dimensions, node offset */

   if ((error = nc_inq_dimid(ncid, "lat", &lat_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo latitude dimension.\n");
       return (NULL);
   }
   error = nc_inq_dimlen(ncid, lat_dim, &len );
   z_info.ny = (int) len;

   if ((error = nc_inq_dimid(ncid, "lon", &lon_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo longitude dimension.\n");
       return (NULL);
   }
   error = nc_inq_dimlen(ncid, lon_dim, &len );
   z_info.nx = (int) len;
   nsqin = z_info.nx * z_info.ny;
  
   if ((error = nc_get_att_int(ncid, NC_GLOBAL, "node_offset", &lval)) != NC_NOERR ) {
      fprintf(stderr,"\nError reading cdf file attribute: node_offset.\n");
       return (NULL);
   }
   z_info.node_offset = lval;
   hptr->node_offset = lval;
   
   /* get lat and lon vectors */
   
   
   latptr = (double *) calloc(z_info.ny, sizeof(double));
   lonptr = (double *) calloc(z_info.nx, sizeof(double));
   
   if ((error = nc_inq_varid(ncid, "lon", &varid)) != NC_NOERR ) {
      fprintf(stderr,"\nNo variable 'lon' in file.\n");
      return(NULL);
   }
   
   if ( (error = nc_inq_vartype(ncid, varid, &vtype)) == NC_NOERR ) {
      if (vtype != NC_DOUBLE) {
        fprintf(stderr,"\nExpecting longitude variable to be of type NC_DOUBLE.\n");
        fprintf(stderr,"Not a valid topography file\n");
        return(NULL);
      }
   }
   
    if ((error = nc_get_var_double(ncid, varid, lonptr)) != NC_NOERR) {
      fprintf(stderr,"\nError reading longitude values.\n");
      return(NULL);
   }
   z_info.x_min= lonptr[0];
   z_info.x_max =lonptr[z_info.nx-1];
  
   if ((error = nc_inq_varid(ncid, "lat", &varid)) != NC_NOERR ) {
      fprintf(stderr,"\nNo variable 'lat' in file.\n");
      return(NULL);
   }
   if ( (error = nc_inq_vartype(ncid, varid, &vtype)) == NC_NOERR ) {
      if (vtype != NC_DOUBLE) {
        fprintf(stderr,"\nExpecting latitude variable to be of type NC_DOUBLE.\n");
        fprintf(stderr,"Not a valid topography file\n");
        return(NULL);
      }
   }
    if ((error = nc_get_var_double(ncid, varid, latptr)) != NC_NOERR) {
      fprintf(stderr,"\nError reading latitude values.\n");
      return(NULL);
   }

   z_info.y_min = latptr[0];
   z_info.y_max = latptr[z_info.ny-1];

      
   if ( z_info.node_offset == 0)  { 
      z_info.x_inc = (z_info.x_max - z_info.x_min) / (double)(z_info.nx - 1);
      z_info.y_inc = (z_info.y_max - z_info.y_min) / (double) (z_info.ny - 1);
   }
   else {
      z_info.x_inc = (z_info.x_max - z_info.x_min) / (double) z_info.nx ;
      z_info.y_inc = (z_info.y_max - z_info.y_min) / (double) z_info.ny;
   }

    z_info.x_inc =  ((double) (NINT(z_info.x_inc * 1000000000))) / 1000000000.0;
    z_info.y_inc =  ((double) (NINT(z_info.y_inc * 1000000000))) / 1000000000.0;
    
    z_info.lon0to360 = 1;
    if (z_info.x_min < 0)
       z_info.lon0to360 = 0;
       
    z_info.xgreenwich = 0;
    if ((z_info.x_min < 0) && (z_info.x_max >= 0))
         z_info.xgreenwich = 1;
      

/* Set fields of output grid   */
   
   hptr->x_inc = z_info.x_inc;
   hptr->y_inc = z_info.y_inc;

   hptr->lon0to360 = 1;
   if (hptr->x_min < 0) {
      hptr->lon0to360 = 0; 
   }
   lon0to360 = hptr->lon0to360;

   hptr->xgreenwich = 0;
   if (hptr->x_min < 0 && hptr->x_max >= 0)
      hptr->xgreenwich = 1;
      
 /*  Determine if output longitude ranges are continuous or split regions in the input grid 
     After this, lon1min, lon1max, lon2min, lon2max will specify longitude regions needed from input topography 
     Fields of hptr-> x_min, x_max, y_min, y_max will reflect output topography grid values */
   
 
   if (lon0to360 == z_info.lon0to360) {    /* input/output grids have same longitude convention */

       split_region = 0;      
       lon1min = hptr->x_min;
       lon1max = hptr->x_max;
       
       error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
       error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
       error += xy2ij(&z_info, lon1max, hptr->y_max, &icol, &irow);
       error += ij2xy(&z_info, icol, irow, &lon1max, &hptr->y_max);
       
       hptr->x_min = lon1min - hptr->node_offset * hptr->x_inc;
       hptr->x_max = lon1max + hptr->node_offset * hptr->x_inc;
       
       if (!lon0to360) {            /* both are mixed sign, but requested area crosses dateline */
          if (hptr->x_max > 180) {
	      split_region = 1;
	      lon1min = hptr->x_min;
	      lon1max = 180;
	      lon2min = -180;
	      lon2max = hptr->x_max - 360;
	      
              error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
              error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
	      
              error += xy2ij(&z_info, lon2max, hptr->y_max, &icol, &irow);
              error += ij2xy(&z_info, icol, irow, &lon2max, &hptr->y_max);
	      
	      hptr->x_min = lon1min - hptr->node_offset * hptr->x_inc * 0.5;
	      hptr->x_max = (lon2max + 360) + hptr->node_offset * hptr->x_inc * 0.5;
	  }
	  
	  else {
	     if (hptr->x_min < -180) {
	       split_region = 1;
	       lon1min = hptr->x_min + 360;
	       lon1max = 180;
	       lon2min = -180;
	       lon2max = hptr->x_max;
	  
               error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
               error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
               error += xy2ij(&z_info, lon2max, hptr->y_max, &icol, &irow);
               error += ij2xy(&z_info, icol, irow, &lon2max, &hptr->y_max);
	      
	       hptr->x_min = (lon1min - 360) - hptr->node_offset * hptr->x_inc * 0.5;
	       hptr->x_max = lon2max  + hptr->node_offset * hptr->x_inc * 0.5;
	    }
	 }
       }
   }  /* end if same conventions */
   
   
   if (lon0to360 != z_info.lon0to360) {  /* input/output grids have different conventions */
   
      if (lon0to360)  {  /* input grid is -180 -> +180, output is all positive lons */ 
	  
	  /* no dateline crossing */
	  
	  if (hptr->x_max <= 180) {  
	     split_region = 0;      
             lon1min = hptr->x_min;
             lon1max = hptr->x_max;
	  
             error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
             error += xy2ij(&z_info, lon1max, hptr->y_max, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon1max, &hptr->y_max);
       
             hptr->x_min = lon1min - hptr->node_offset * hptr->x_inc * 0.5;
             hptr->x_max = lon1max + hptr->node_offset * hptr->x_inc * 0.5;
	  }
	  
	  /* case where we cross dateline */
	  
	  if (hptr->x_max > 180) {
	     split_region = 1;        
             lon1min = hptr->x_min;
	     lon1max = 180;
	     lon2min = -180;
	     lon2max = hptr->x_max - 360;

             error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
             error += xy2ij(&z_info, lon2max, hptr->y_max, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon2max, &hptr->y_max);
	      
	     hptr->x_min = lon1min - hptr->node_offset * hptr->x_inc * 0.5;
	     hptr->x_max = (lon2max+360) + hptr->node_offset * hptr->x_inc * 0.5;
	  }
       }
       
       /* input grid is 0->360, output has negative longitudes */
       
       if (!lon0to360) {   
       
          /* case where output region specified as all negative lons */
	 if (hptr->x_max < 0) {
	    split_region = 0;   
	    lon1min = hptr->x_min + 360;
 	    lon1max = hptr->x_max + 360;
	  
             error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
             error += xy2ij(&z_info, lon1max, hptr->y_max, &icol, &irow);
             error += ij2xy(&z_info, icol, irow, &lon1max, &hptr->y_max);
             
	     hptr->x_min = (lon1min-360) - hptr->node_offset * hptr->x_inc * 0.5;
             hptr->x_max = (lon1max-360) + hptr->node_offset * hptr->x_inc * 0.5;
	 }  
	 
	 if (hptr->xgreenwich) {                
	    split_region = 1;
	    lon1min = hptr->x_min + 360;
	    lon1max = 360;
	    lon2min = 0;
	    lon2max = hptr->x_max;
	    
            error = xy2ij(&z_info, lon1min, hptr->y_min, &icol, &irow);
            error += ij2xy(&z_info, icol, irow, &lon1min, &hptr->y_min);
            error += xy2ij(&z_info, lon2max, hptr->y_max, &icol, &irow);
            error += ij2xy(&z_info, icol, irow, &lon2max, &hptr->y_max);
	      
	    hptr->x_min = (lon1min-360) - hptr->node_offset * hptr->x_inc * 0.5;
	    hptr->x_max =  lon2max + hptr->node_offset * hptr->x_inc * 0.5;
	 }
      }
   } /* end if different conventions */
   
   /* adjust output latitude ranges to actual topography grid */

   if (hptr->node_offset) {
        hptr->y_min -=  hptr->y_inc * 0.5;
        hptr->y_max +=  hptr->y_inc * 0.5;
   }
  
    hptr->nx = 1 + NINT((hptr->x_max - hptr->x_min) / hptr->x_inc);
    hptr->ny = 1 + NINT((hptr->y_max - hptr->y_min) / hptr->y_inc);
    
   if (hptr->node_offset)  {
    --hptr->nx;
    --hptr->ny;
   }
   
/******************************/
/* read in topography values */
   
   zin = (short *) calloc(nsqin, sizeof(short));
   if (zin == NULL) {
      fprintf(stderr,"\nUnable to allocate memory for topography data.\n");
      return(NULL);
   }
   
   if ((error = nc_inq_varid(ncid, "z", &varid)) != NC_NOERR ) {
      fprintf(stderr,"\nNo variable 'z' in file.\n");
      return(NULL);
   }
   
   if ( (error = nc_inq_vartype(ncid, varid, &vtype)) == NC_NOERR ) {
      if (vtype != NC_SHORT) {
        fprintf(stderr,"\nExpecting z-values to be of type NC_SHORT.\n");
        fprintf(stderr,"Not a valid topography file\n");
        return(NULL);
      }
   }
   if ( (error = nc_get_att_short(ncid, varid, "_FillValue", &no_val)) != NC_NOERR)
        no_val = NC_FILL_SHORT;  

   *missing = no_val;
   
   if ((error = nc_get_var_short(ncid, varid, zin)) != NC_NOERR) {
      fprintf(stderr,"\nError reading topography values.\n");
      return(NULL);
   }
   
   nc_close(ncid);
   
   /* adjust seafloor values to requested sign */
   
   error = xy2ij(&z_info, -55.0, 25.0, &i, &j);  /* coordinates in mid-ocean */
   z_is_neg = 1;
   if (zin[j*z_info.nx+i] > 0 )
       z_is_neg = 0;
       
   if (depths_neg) {
     if (!z_is_neg) {
        
	for (i = 0; i <nsqin; ++i) {
	  if (zin[i] != no_val)
	     zin[i] *= -1;
	}
     }
   }
   else {
      if (z_is_neg) {
	for (i = 0; i <nsqin; ++i) {
	   if (zin[i] != no_val)
	      zin[i] *= -1;
	}
      }
   }
   
/******************************/
/*    Output the data....    */

  /* Simplest case:  output entire dataset in same order */

   if ( (lon0to360 == z_info.lon0to360) && (z_info.ny == hptr->ny) && (z_info.nx == hptr->nx)) {
     
         *lonvec_addr = lonptr;
         *latvec_addr = latptr;
         return(zin);         /* done! */
   }
   
   
 /* Output a subset of topography */
    
   nsqout =  hptr->nx * hptr->ny;   
   seafloor = (short *) calloc( nsqout, sizeof(short));
      if (seafloor == NULL) {
         fprintf(stderr,"\nUnable to allocate memory for topography data.\n");
         return(NULL);
      }
      
   lonout = (double *) calloc(hptr->nx, sizeof(double));
   latout = (double *) calloc(hptr->ny, sizeof(double));


   if (!split_region) {
    
      /*  find start/end for subset range; 
	  transfer values element by element within box */
 	  
	  xy2ij(&z_info, lon1min, hptr->y_min, &start_col, &start_row);
	  xy2ij(&z_info, lon1max, hptr->y_max, &end_col, &end_row);
	  
	  error = (end_row - start_row + 1) != hptr->ny;
	  error += (end_col - start_col + 1) != hptr->nx;
	  
	  if (error) {
	     fprintf(stderr,"FATAL ERROR in hb_get_topo():  mismatch of output grid dimensions");
             return(NULL);
	  }
	  
	  indx = 0;
	  for (i = start_row; i <= end_row; ++i) {
	     for (j = start_col; j <= end_col; ++j) {
	        seafloor[indx++] = zin[i * z_info.nx + j];
	     }
	  }
  
          indx = 0;
	  for (i = start_row; i <= end_row; ++i) 
              latout[indx++] = hptr->y_min + hptr->node_offset * 0.5 * hptr->y_inc + hptr->y_inc * indx;
	      
	  
          indx = 0;
	  for (i = start_col; i <= end_col; ++i) 
              lonout[indx++] = hptr->x_min + hptr->node_offset * 0.5 * hptr->x_inc + hptr->x_inc * indx;
	  
	  free(latptr);
	  free(lonptr);
	  free(zin);
	  
         *lonvec_addr = lonout;
         *latvec_addr = latout;
	  return(seafloor);          /* done! */
    }
    
    
   /* Concatenate split regions for output 
      beginning with western region...*/
    
    error = xy2ij(&z_info, lon1min, hptr->y_min, &start_col, &start_row);   
    error += xy2ij(&z_info, lon1max, hptr->y_max, &end_col, &end_row); 
	  if (error) {
	     fprintf(stderr,"FATAL ERROR 1 in hb_get_topo():  mismatch in calculated input/output grid dimensions");
             return(NULL);
	  }
       
    strow_out = 0;
    stcol_out = 0;
    endrow_out = strow_out + (end_row - start_row);
    endcol_out = endcol_out + (end_col - start_col);
          
    indx = strow_out;
    for (i = start_row; i <= end_row; ++i) 
       latout[indx++] = hptr->y_min + hptr->node_offset * 0.5 * hptr->y_inc + hptr->y_inc * indx;
   
    
    indx = stcol_out;
    for (i = start_col; i <= end_col; ++i) 
       lonout[indx++] = hptr->x_min + hptr->node_offset * 0.5 * hptr->x_inc + hptr->x_inc * indx;
    
    
    i = strow_out - 1;
    for (irow = start_row; irow <= end_row; ++irow) {
        ++i;
        j = stcol_out - 1;
       for (icol = start_col; icol <= end_col; ++icol) {
          ++j;
	  indx = i * hptr->nx + j;
	  seafloor[indx] = zin[irow * z_info.nx + icol];
       }
    }

    /* add the eastern region.... */ 
   
    error = xy2ij(&z_info, lon2min, hptr->y_min, &start_col, &start_row);   
    error += xy2ij(&z_info, lon2max, hptr->y_max, &end_col, &end_row); 
	  if (error) {
	     fprintf(stderr,"FATAL ERROR 3 in hb_get_topo():  mismatch in calculated input/output grid dimensions");
             return(NULL);
	  }
       
    endrow_out = hptr->ny - 1;
    endcol_out = hptr->nx - 1;
    strow_out = endrow_out - (end_row - start_row);
    stcol_out = endcol_out - (end_col - start_col);
    
    indx = endcol_out;
    while (indx >= stcol_out) {
       lonout[indx] = hptr->x_max - (hptr->node_offset * 0.5 * hptr->x_inc) - (indx * hptr->x_inc);
       --indx;
    }
     
    indx = endrow_out;
    while (indx >= strow_out) {
        latout[indx]  = hptr->y_max - (hptr->node_offset * 0.5 * hptr->y_inc) - (indx * hptr->y_inc);
	--indx;
    }
    
    i = strow_out - 1;
    for (irow = start_row; irow <= end_row; ++irow) {
        ++i;
        j = stcol_out - 1;
       for (icol = start_col; icol <= end_col; ++icol) {
          ++j;
	  indx = i * hptr->nx + j;
	  seafloor[indx] = zin[irow * z_info.nx + icol];
       }
    }
   
    free(zin);
    free(latptr);
    free(lonptr);
    
    *lonvec_addr = lonout;
    *latvec_addr = latout;
    return (seafloor);
    
} /* end hb_get_topo() */

/*********************************************/
short find_nearest_topo_val(double lat, double lon, short *topo, struct GRID_INFO *gptr)
   /* Find bottom depth for lat/lon specified. Returns large negative
      value (NC_FILL_SHORT) if an error occurs. */
 
      
{
int i, j ;

   if (gptr->lon0to360) {   /* 0 -> 360 */
      if (lon < 0)
          lon += 360.0;
   }
      
   if (!gptr->lon0to360 ) {   /* mixed sign */
      if (lon < gptr->x_min )
        if (ABS(lon - gptr->x_min) < 0.001)
	    lon = gptr->x_min;
	else  
            lon += 360.0;
      
      if ( lon > gptr->x_max)
        if (ABS(lon - gptr->x_max) < 0.001)
            lon = gptr->x_max;
        else 
	    lon -= 360.0;
   }

   /* adjust lat for the poles */
   if (lat >= 89.9999)
       lat = 89.99;
   if  (lat <= -89.9999)
       lat = -89.99;
       
       
   if (xy2ij(gptr, lon, lat, &i, &j) < 0) {
       fprintf(stderr,"\nPosition exceeds geographic bounds of topo file Lat:%.2lf Lon: %.2lf \n", lat, lon);
       return(NC_FILL_SHORT);
   }
   
   if (i >= gptr->nx  || j >= gptr->ny || i<0 || j<0 )
          return(NC_FILL_SHORT);

   return(topo[j*gptr->nx+i]);
   

} /* end find_nearest_topo_val() */
