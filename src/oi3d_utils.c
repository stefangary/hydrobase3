/* oi3d_utils.c
................................................................................
                          *******  HydroBase3  *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             March 2011
			     
			     
................................................................................
.
.  A set of routines to facilitate the creation and reading of parameter files
.  (netcdf)  used in HydroBase 3D optimal interpolation routines.  

To create a new parameter file use:
        ncparms_init()
	ncparms_define()  with initialize = 1
	foreach lat/lon
	    ncparms_write()
	ncparms_close()
	
To update a parameter file use:
      	ncparms_update() 
	ncparms_getinfo()  
	foreach lat/lon
	    ncparms_write()
	ncparms_close()
	
To read a parameter file use:
        ncparms_open()
	ncparms_getinfo()
	foreach lat/lon
	    ncparms_read()
	ncparms_close()
	    
======================================================
.  int ncparms_init(char *filename)    
	Opens a new  file for writing.   

.  int ncparms_open(char *name)
           Open a nc file for reading.

.  int ncparms_update(char *name)
           Open a nc file for reading and writing.
	   
.  void ncparms_close(int ncid)
         Close a nc file.

   int ncparms_define(int ncid, struct NC_INFO *ncinfo, float *Z)
        Defines all variables and attributes for a parameters file.
	
   int ncparms_write(int ncid, float lat, float lon, struct NC_INFO *ncinfo, struct NC_PARMS *parms)

        Writes all parameters for a latitude,longitude to an already opened netcdf file. Returns 0  
         for a successful write or prints an error message and exits in case of an error.

    
   int ncparms_write_var(int ncid,char *vname, int *start, int *count, float *vptr)
        Writes contiguous block of float values to a predefined variable
   
   int ncparms_write_coord(int ncid, char *vname, float *vector, int n)
        Writes coordinate vector
.
.  int ncparms_getinfo(int ncid, struct NC_INFO *ncinfo, struct NC_PARMS *ncdata, int initialize)
         Fills in structure fields from variables/attributes in the netcdf parameters file.
         If initialize = 0, the fields in ncdata and ncinfo must be preallocated by calling function.
         if non-zero, memory for the arrays will be allocated here.
         
   int ncparms_getcoord(int ncid, char *dimname, char *vname, float *vptr, char *units)
          Reads coordinate variable from netcdf file.  Memory for the vector and units array 
          is already allocated.  Returns 0 if successful or -1 in case of an error
	  
   int ncparms_read(int ncid, float lat, float lon, struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)
         Reads in all parameters for a single lat/lon position (all depths) from an already
         open nc file.  The fields for ncinfo must have already been read in
         using  ncparms_getinfo(), and memeory for arrays in ncdata must be already allocated. 
	  
   int ncparms_getvar(int ncid,char *vname, int row, int col, int nz, float *vptr)
         Retrieves values of variable (vname) at specified lat/lon (all depths) 
         from netcdf file. Memory for the  array (at vptr)  is already allocated 
         Returns 0 if successful or -1 in case of an error

   int ncparms_rowcol( struct NC_INFO *ncinfo, float lat, float lon, int *rowptr, int *colptr)
         Returns row,col associated with lat,lon based on ncinfo  

   int ncparms_latlon( struct NC_INFO *ncinfo, int row, int col, float *latptr, float *lonptr)

   int ncparms_depth_indx(float theDepth, float *Z, int nz)
          Returns the standard level index corresponding to depth.  Or a
          negative number if some error occurs:
                 -1 :  depth not in Z 
                 -99:  depth exceeds greatest Z level.

   void ncparms_zero(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)
           initializes all fields in these structures 
	   
   void ncparms_alloc(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)
            allocates memory for structure fields

   void ncparms_prefill(struct NC_PARMS *ncdata, float fill_value)
            fills parameter arrays with specified value 
   
    void ncparms_free(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)
         Frees memory for structure fields 
====================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "hb_oi3d.h"


int ncparms_init(char *filename)
/* Creates a new netCDF file and returns the netCDF id for that file */
{
   int ncid, error;

   error = nc_create(filename, NC_NOCLOBBER|NC_64BIT_OFFSET, &ncid);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nUnable to create nc_file: %s", filename);
      fprintf(stderr, "\nPerhaps it exists already??");
      fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
      exit(1);
   }
   return (ncid);
   
}  /* end cdf_init() */

/***************************************************************************/

int ncparms_open(char *fname)

   /* opens an existing cdf file for reading and returns the file descriptor 
      to be used in future read operations. If file does not exist, a message
      is printed on stderr and -1 is returned.
      
   */
{
   int ncid;
   int i, error;


   error = nc_open(fname, NC_NOWRITE, &ncid); 
   
   if (error == NC_NOERR) { 
      fprintf(stderr,"\nOpened %s ...", fname);
      return(ncid);
   }
   else {
      fprintf(stderr,"\n Unable to open %s.\n", fname);
      fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
      return(-1);
   }


} /* end ncparms_open() */

/***************************************************************************/

int ncparms_update(char *fname)

   /* opens an existing nc file for reading/writing and returns the 
      file descriptor to be used in future read/write operations. 
      If file does not exist, a message is printed on stderr and -1 
      is returned.
   */
{
   int ncid;
   int i, error;

 error = nc_open(fname, NC_WRITE, &ncid);     
 if (error != NC_NOERR) {
    fprintf(stderr,"\n Unable to open %s.\n", fname);
    fprintf(stderr, "\n Error: \n%s", nc_strerror(error));
  }
  else {
     fprintf(stderr,"\nOpened %s ...", fname);
  }
  
   return (ncid);

} /* end ncparms_update() */

/***************************************************************************/

void ncparms_close(int ncid)
{
   nc_close(ncid);
   return;

} /* end ncparms_close() */

/*********************************************************************/ 
int ncparms_define(int ncid, struct NC_INFO *ncinfo, float *Z)
 /*  Defines all variables and attributes for a parameters file.
 */

{
   int lat_dim, lon_dim, depth_dim, depthseas_dim;
   int depth_id, lat_id, lon_id, depthseas_id;
   int varid, ival;
   int dims[3], n, error;
   char *descrip;
   float fval;

 /* define dimensions */

   error = nc_def_dim(ncid, "lat", (size_t) ncinfo->ny, &lat_dim);
   if (error != NC_NOERR) {
      fprintf(stderr, "\n Error defining lat_dim: \n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_def_dim(ncid, "lon", (size_t) ncinfo->nx, &lon_dim);
   error = nc_def_dim(ncid, "depth", (size_t) ncinfo->nz, &depth_dim);
   error = nc_def_dim(ncid, "depthseas", (size_t) ncinfo->nz_seas, &depthseas_dim);

/* define coordinate variables */

   dims[0] = lat_dim;
   error = nc_def_var(ncid, "latitude", NC_FLOAT, 1, dims, &lat_id);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nCan't define variable:\n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_put_att_text(ncid, lat_id, "units", 14, "degrees_north");
   error = nc_put_att_text(ncid, lat_id, "long_name", 9, "Latitude");
   error = nc_put_att_text(ncid, lat_id, "generic_name", 9, "latitude");

   dims[0] = lon_dim;
   error = nc_def_var(ncid, "longitude", NC_FLOAT, 1, dims, &lon_id);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nCan't define variable:\n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_put_att_text(ncid, lon_id, "units", 13, "degrees_east");
   error = nc_put_att_text(ncid, lon_id, "long_name", 10, "Longitude");
   error = nc_put_att_text(ncid, lon_id, "generic_name", 10, "longitude");

   dims[0] = depth_dim;
   error = nc_def_var(ncid, "depth", NC_FLOAT, 1, dims, &depth_id);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nCan't define variable:\n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_put_att_text(ncid, depth_id, "units", 7, "meters");
   error = nc_put_att_text(ncid, depth_id, "long_name", 9, "Depth (m)");
   error = nc_put_att_text(ncid, depth_id, "generic_name", 6, "depth");

   dims[0] = depthseas_dim;
   error = nc_def_var(ncid, "depthseas", NC_FLOAT, 1, dims, &depthseas_id);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nCan't define variable:\n%s", nc_strerror(error));
      return (-1);
   }
   error = nc_put_att_text(ncid, depthseas_id, "units", 7, "meters");
   error = nc_put_att_text(ncid, depthseas_id, "long_name", 27, "Depth of seasonal layer (m)");
   error = nc_put_att_text(ncid, depthseas_id, "generic_name", 6, "depth");

/* define other variables */

   dims[0] = lat_dim;
   dims[1] = lon_dim;

   error = nc_def_var(ncid, "dZdx", NC_FLOAT, 2, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 25, "Zonal bathymetry gradient");
   error = nc_put_att_text(ncid, varid, "units", 16, "meters/kilometer");

   error = nc_def_var(ncid, "dZdy", NC_FLOAT, 2, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 30, "Meridional bathymetry gradient");
   error = nc_put_att_text(ncid, varid, "units", 16, "meters/kilometer");

   dims[0] = lat_dim;
   dims[1] = lon_dim;
   dims[2] = depth_dim;

   error = nc_def_var(ncid, "Tparm0", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 23, "Parameter 0 temperature");
   error = nc_put_att_text(ncid, varid, "description", 43, "Nugget value (noise) for covariance function");

   error = nc_def_var(ncid, "Tparm1", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 23, "Parameter 1 temperature");
   error = nc_put_att_text(ncid, varid, "description", 52, "Sill - nugget value (signal) for covariance function");

   error = nc_def_var(ncid, "Tparm2", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 23, "Parameter 2 temperature");
   error = nc_put_att_text(ncid, varid, "description", 24, "Decorrelation scale (km)");

   error = nc_def_var(ncid, "Sparm0", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 0 salinity");
   error = nc_put_att_text(ncid, varid, "description", 43, "Nugget value (noise) for covariance function");

   error = nc_def_var(ncid, "Sparm1", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 1 salinity");
   error = nc_put_att_text(ncid, varid, "description", 52, "Sill - nugget value (signal) for covariance function");

   error = nc_def_var(ncid, "Sparm2", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 2 salinity");
   error = nc_put_att_text(ncid, varid, "description", 24, "Decorrelation scale (km)");

   error = nc_def_var(ncid, "Pparm0", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 0 pressure");
   error = nc_put_att_text(ncid, varid, "description", 43, "Nugget value (noise) for covariance function");

   error = nc_def_var(ncid, "Pparm1", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 1 pressure");
   error = nc_put_att_text(ncid, varid, "description", 52, "Sill - nugget value (signal) for covariance function");

   error = nc_def_var(ncid, "Pparm2", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 20, "Parameter 2 pressure");
   error = nc_put_att_text(ncid, varid, "description", 24, "Decorrelation scale (km)");

   error = nc_def_var(ncid, "Xlength", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 26, "X radius of search ellipse");
   error = nc_put_att_text(ncid, varid, "units", 10, "kilometers");

   error = nc_def_var(ncid, "Ylength", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 26, "Y radius of search ellipse");
   error = nc_put_att_text(ncid, varid, "units", 10, "kilometers");

   error = nc_def_var(ncid, "dTdx", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 26, "Zonal Temperature Gradient");
   error = nc_put_att_text(ncid, varid, "units", 19, "degrees C/kilometer");

   error = nc_def_var(ncid, "dTdy", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 31, "Meridional Temperature Gradient");
   error = nc_put_att_text(ncid, varid, "units", 19, "degrees C/kilometer");

   error = nc_def_var(ncid, "dPdx", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 24, "Zonal Pressure Gradient");
   error = nc_put_att_text(ncid, varid, "units", 18, "decibars/kilometer");

   error = nc_def_var(ncid, "dPdy", NC_FLOAT, 3, dims, &varid);
   error = nc_put_att_float(ncid, varid, "FillValue", NC_FLOAT, 1, &ncinfo->fill_value);
   error = nc_put_att_float(ncid, varid, "MaskValue", NC_FLOAT, 1, &ncinfo->mask_value);
   error = nc_put_att_text(ncid, varid, "long_name", 29, "Meridional Pressure Gradient");
   error = nc_put_att_text(ncid, varid, "units", 18, "decibars/kilometer");


/* assign global attributes */ 

   fval = ncinfo->lats[0] - 0.5 * ncinfo->yincr;
   error = nc_put_att_float(ncid, NC_GLOBAL, "latmin", NC_FLOAT, 1, &fval);

   n = ncinfo->ny - 1;
   fval = ncinfo->lats[n] + 0.5 * ncinfo->yincr;
   error = nc_put_att_float(ncid, NC_GLOBAL, "latmax", NC_FLOAT, 1, &fval);
   error = nc_put_att_float(ncid, NC_GLOBAL, "latincr", NC_FLOAT, 1, &ncinfo->yincr);
   
   fval = ncinfo->lons[0] - 0.5 * ncinfo->xincr;
   error = nc_put_att_float(ncid, NC_GLOBAL, "lonmin", NC_FLOAT, 1, &fval);

   n = ncinfo->nx - 1;
   fval = ncinfo->lons[n] + 0.5 * ncinfo->xincr;
   error = nc_put_att_float(ncid, NC_GLOBAL, "lonmax", NC_FLOAT, 1, &fval);
   error = nc_put_att_float(ncid, NC_GLOBAL, "lonincr", NC_FLOAT, 1, &ncinfo->xincr);

   ival = ncinfo->modelCode;
   error = nc_put_att_int(ncid, NC_GLOBAL, "Variogram_Model_Code", NC_INT, 1, &ival);
   
   descrip = (char *) calloc(300, sizeof(char));
   strcpy(descrip,"1=exponential, 2=gaussian, 3=gaussian");
   error = nc_put_att_text(ncid, NC_GLOBAL, "Variogram_Model_Values", strlen(descrip), descrip);
   free(descrip);
   
   fval = NCPARMS_FILL;
   error = nc_put_att_float(ncid, NC_GLOBAL, "FillValue", NC_FLOAT, 1, &fval);
   fval = NCPARMS_MASK;
   error = nc_put_att_float(ncid, NC_GLOBAL, "MaskValue", NC_FLOAT, 1, &fval);

   error = nc_put_att_text(ncid, NC_GLOBAL, "title", 45, "HydroBase3 covariance and gradient parameters");

   error = nc_enddef(ncid);       /* leave define mode */
   
   /* write out coordinate vectors */
   
   error = ncparms_write_coord(ncid, "latitude", ncinfo->lats, ncinfo->ny);
   error = ncparms_write_coord(ncid, "longitude", ncinfo->lons, ncinfo->nx);
   error = ncparms_write_coord(ncid, "depth", Z, ncinfo->nz);
   error = ncparms_write_coord(ncid, "depthseas", Z, ncinfo->nz_seas);
   
   return 0;

 } /* end ncparms_define */
 
/*********************************************************************/ 
int ncparms_write_coord(int ncid, char *vname, float *vector, int n)

{
int error, varid;
size_t start[1], count[1];
 
    start[0] = 0;
    count[0] = n;
    
    error = nc_inq_varid(ncid, vname, &varid);
    if (error != NC_NOERR) {
        fprintf(stderr,"\n %s variable undefined in the nc file", vname);
        fprintf(stderr, "\n %s", nc_strerror(error));
        exit(1);
    }
   error = nc_put_vara_float(ncid, varid, start, count, vector);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nError writing coordinate variable (%s) to nc parameter file.\n", vname);
      fprintf(stderr, "\n %s", nc_strerror(error));
      exit(1);
   }
   return(0);
    
}  /* end ncpars_write_coord */


/*********************************************************************/ 
	
int ncparms_write(int ncid, float lat, float lon, struct NC_INFO *ncinfo, struct NC_PARMS *parms)

 /*  Writes all parameters for a latitude,longitude to an already opened netcdf file. Returns 0  
   for a successful write or prints an error message and exits in case of an error.
 */
 
{
   size_t start[3], count[3];
   int row, col;
   int error;
   
   error = ncparms_rowcol(ncinfo, lat, lon, &row, &col);

   start[0] = row;
   start[1] = col;
   start[2] = 0;
   
   count[0] = 1;
   count[1] = 1;
   count[2] = ncinfo->nz;
   
   error = 0;
   error += ncparms_write_var(ncid, "Tparm0", start, count, parms->Tparm0);
   error += ncparms_write_var(ncid, "Tparm1", start, count, parms->Tparm1);
   error += ncparms_write_var(ncid, "Tparm2", start, count, parms->Tparm2);
   error += ncparms_write_var(ncid, "Sparm0", start, count, parms->Sparm0);
   error += ncparms_write_var(ncid, "Sparm1", start, count, parms->Sparm1);
   error += ncparms_write_var(ncid, "Sparm2", start, count, parms->Sparm2);
   error += ncparms_write_var(ncid, "Pparm0", start, count, parms->Pparm0);
   error += ncparms_write_var(ncid, "Pparm1", start, count, parms->Pparm1);
   error += ncparms_write_var(ncid, "Pparm2", start, count, parms->Pparm2);
   error += ncparms_write_var(ncid, "Xlength", start, count, parms->Lx);
   error += ncparms_write_var(ncid, "Ylength", start, count, parms->Ly);
   error += ncparms_write_var(ncid, "dTdx", start, count, parms->dTdx);
   error += ncparms_write_var(ncid,"dTdy", start, count, parms->dTdy);
   error += ncparms_write_var(ncid, "dPdx", start, count, parms->dPdx);
   error += ncparms_write_var(ncid, "dPdy", start, count, parms->dPdy);

   count[2] = 0;
   error += ncparms_write_var(ncid, "dZdx", start, count, &parms->dZdx);
   error += ncparms_write_var(ncid, "dZdy", start, count, &parms->dZdy);
      
   if (error) {
       fprintf(stderr,"Error writing to parameter file at lat: %.3f  lon: %.3f\n", lat, lon);
       exit(1);
   }
   
   return(0);

} /* end ncparms_write */
 
/*********************************************************************/ 
 int ncparms_write_var(int ncid,char *vname, size_t *start, size_t *count, float *vptr)
   
 {
 
  int error, varid;

     error = nc_inq_varid(ncid, vname, &varid);
    if (error != NC_NOERR) {
        fprintf(stderr,"\n %s  undefined in the nc file", vname);
        fprintf(stderr, "\n %s", nc_strerror(error));
        exit(1);
    }
   error = nc_put_vara_float(ncid, varid, start, count, vptr);
   if (error != NC_NOERR) {
      fprintf(stderr, "\nError writing %s to nc paramter file.\n", vname);
      fprintf(stderr, "\n %s", nc_strerror(error));
      exit(1);
   }
   
   return(0);

 } /* end ncparms_write_var */
 
/*********************************************************************/ 
 
int ncparms_getinfo(int ncid, struct NC_INFO *ncinfo, struct NC_PARMS *ncdata, int initialize )
  /*  Fills in structure fields from variables/attributes in the netcdf parameters file.
      If initialize = 0, the fields in ncdata and ncinfo are preallocated by calling function.
        if non-zero, memory for the arrays will be allocated here.
      Returns 0 for no error.
   */
     
{
  int error, varid;
  int ndims, nvars, ngatts, recdim;
  int lat_dim, lon_dim, depth_dim, depthseas_dim;
  float fval;
  int lval;
  size_t start[1], count[1], len;
  
  
   if ((error = nc_inq(ncid, &ndims, &nvars, &ngatts, &recdim)) != NC_NOERR){
      fprintf(stderr,"\nError reading netcdf file.\n");
      return (1);
   }
   
   /* get dimensions */
 
   if ((error = nc_inq_dimid(ncid, "lat", &lat_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo lat dimension.\n");
       return (2);
   }
   error = nc_inq_dimlen(ncid, lat_dim, &len );
   ncinfo->ny = (int) len;
 
 
   if ((error = nc_inq_dimid(ncid, "lon", &lon_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo lon dimension.\n");
       return (2);
   }
   error = nc_inq_dimlen(ncid, lon_dim, &len );
   ncinfo->nx = (int) len;
 
   if ((error = nc_inq_dimid(ncid, "depth", &depth_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo depth dimension.\n");
       return (2);
   }
   error = nc_inq_dimlen(ncid, depth_dim, &len );
   ncinfo->nz = (int) len;
   ncdata->nz = (int) len;
 
   if ((error = nc_inq_dimid(ncid, "depthseas", &depthseas_dim)) != NC_NOERR ) {
      fprintf(stderr,"\nNo seasonal depth dimension.\n");
       return (2);
   }
   error = nc_inq_dimlen(ncid, depthseas_dim, &len );
   ncinfo->nz_seas = (int) len;
   ncdata->nz_seas = (int) len;
   
   
   if (initialize) {
       ncparms_alloc(ncinfo, ncdata);
   }

   
   /* get coordinate variables */
   ncinfo->y_units = (char *) calloc(80, sizeof(char));
   ncinfo->x_units = (char *) calloc(80, sizeof(char));
   ncinfo->z_units = (char *) calloc(80, sizeof(char));
   error = ncparms_getcoord(ncid, "lat", "latitude", ncinfo->lats, ncinfo->y_units);
   error = ncparms_getcoord(ncid,"lon", "longitude", ncinfo->lons, ncinfo->x_units);
   error = ncparms_getcoord(ncid, "depth", "depth", ncdata->depths, ncinfo->z_units);
   error = nc_get_att_float(ncid, NC_GLOBAL, "lonincr", &fval);
   ncinfo->xincr = fval;
   error = nc_get_att_float(ncid, NC_GLOBAL, "latincr", &fval);
   ncinfo->yincr = fval;
   error = nc_get_att_float(ncid, NC_GLOBAL, "FillValue", &fval);
   ncinfo->fill_value = fval;
   error = nc_get_att_float(ncid, NC_GLOBAL, "MaskValue", &fval);
   ncinfo->mask_value = fval;
   error = nc_get_att_int(ncid, NC_GLOBAL, "Variogram_Model_Code", &lval);
   ncinfo->modelCode = lval;

   return(0);
}  /* end ncparms_getinfo() */

/*********************************************************************/ 
int ncparms_getcoord(int ncid, char *dimname, char *vname, float *vptr, char *units)
/* Reads coordinate variable from netcdf file.  Memory for the vector and units array 
   is already allocated 
   Returns 0 if successful or -1 in case of an error*/
{
   int vardim, varid, error;
   size_t  start[1], count[1];

   if ((error = nc_inq_dimid(ncid, dimname, &vardim)) != NC_NOERR) {
       fprintf(stderr,"\n Unable to retrieve %s dimension from cdf file\n", dimname);
       return (-1);
   }

   if (nc_inq_dimlen(ncid, vardim, &count[0]) != NC_NOERR)
       return (-1);

   if ( nc_inq_varid(ncid, vname, &varid) != NC_NOERR)
       return (-1);

   start[0] = 0;
   
   if ((error = nc_get_vara_float(ncid, varid, start, count, vptr)) != NC_NOERR){
       fprintf(stderr,"\n Error retrieving %s variable from nc file:", vname);
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return (-1);
   }
   error = nc_get_att_text(ncid, varid, "units", units);
   
   if (error != NC_NOERR) {
       fprintf(stderr,"\n Error retrieving %s units from nc file:", vname);
       fprintf(stderr,"\n%s\n", nc_strerror(error));
   }
   
   return (0);

} /* end ncparms_getcoord() */

/*********************************************************************/ 
int ncparms_read(int ncid, float lat, float lon, struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)
  /* Reads in all parameters for a single lat/lon position (all depths) from an already
     open nc file.  The fields for ncinfo must have already been read in
     using  ncparms_getinfo(), and memeory for arrays in ncdata must be already allocated. */
{
int row, col, error;

   error = ncparms_rowcol(ncinfo, lat, lon, &row, &col);
   if (error) {
      fprintf(stderr,"Lat [%.3f]/ Lon [%.3f] exceeds bounds of nc file \n", lat, lon);
      return(-1);
   }
   
   error = 0;
   error += ncparms_getvar(ncid, "Tparm0", row, col, ncinfo->nz, ncdata->Tparm0);
   error += ncparms_getvar(ncid, "Tparm1", row, col, ncinfo->nz, ncdata->Tparm1);
   error += ncparms_getvar(ncid, "Tparm2", row, col, ncinfo->nz, ncdata->Tparm2);
   error += ncparms_getvar(ncid, "Sparm0", row, col, ncinfo->nz, ncdata->Sparm0);
   error += ncparms_getvar(ncid, "Sparm1", row, col, ncinfo->nz, ncdata->Sparm1);
   error += ncparms_getvar(ncid, "Sparm2", row, col, ncinfo->nz, ncdata->Sparm2);
   error += ncparms_getvar(ncid, "Pparm0", row, col, ncinfo->nz, ncdata->Pparm0);
   error += ncparms_getvar(ncid, "Pparm1", row, col, ncinfo->nz, ncdata->Pparm1);
   error += ncparms_getvar(ncid, "Pparm2", row, col, ncinfo->nz, ncdata->Pparm2);
   error += ncparms_getvar(ncid, "Xlength", row, col, ncinfo->nz, ncdata->Lx);
   error += ncparms_getvar(ncid, "Ylength", row, col, ncinfo->nz, ncdata->Ly);
   error += ncparms_getvar(ncid, "dTdx", row, col, ncinfo->nz, ncdata->dTdx);
   error += ncparms_getvar(ncid, "dTdy", row, col, ncinfo->nz, ncdata->dTdy);
   error += ncparms_getvar(ncid, "dPdx", row, col, ncinfo->nz, ncdata->dPdx);
   error += ncparms_getvar(ncid, "dPdy", row, col, ncinfo->nz, ncdata->dPdy);
   error += ncparms_getvar(ncid, "dZdx", row, col, 1, &ncdata->dZdx);
   error += ncparms_getvar(ncid, "dZdy", row, col, 1, &ncdata->dZdy);
   
 
   if (error) 
     return(-1);
   
   return(0);
   
}   /* end ncparms_read() */

/*********************************************************************/ 
int ncparms_getvar(int ncid,char *vname, int row, int col, int nz, float *vptr)
/* Retrieves values of variable (vname) at specified lat/lon (all depths) 
   from netcdf file.  
   Memory for the  array (at vptr)  is already allocated 
   Returns 0 if successful or -1 in case of an error*/
{
   int varid, error;
   size_t  start[3], count[3];
   
   if ( nc_inq_varid(ncid, vname, &varid) != NC_NOERR)  {
      fprintf(stderr,"\nError reading %s \n", vname);
      fprintf(stderr,"\n%s\n", nc_strerror(error));
      return (-1);
   }
   
   start[0] = row;
   start[1] = col;
   start[2] = 0; 
   
   count[0] = 1;  /* retrieve a single lat/lon */
   count[1] = 1;
   count[2] = nz; /* retrieve all depths at this lat/lon */
   
   error = nc_get_vara_float(ncid, varid, start, count, vptr);
   if (error != NC_NOERR) {
     fprintf(stderr,"\nError reading values for property (%s)\n", vname);
     fprintf(stderr,"\n%s\n", nc_strerror(error));
     return(-1);
   }
   
   return(0);
   
}  /* end ncparms_getvar() */
/*********************************************************************/ 
/*********************************************************************/ 
int ncparms_rowcol( struct NC_INFO *ncinfo, float lat, float lon, int *rowptr, int *colptr)

   /* Returns row,col associated with lat,lon based on ncinfo.
      By definition, longitudes are on the range 0 to 360 
      and grid nodes are at center of squares */
   
{   

   if (lon < 0) 
      lon += 360.0;
      
     *rowptr = NINT(((lat + .0001) - ncinfo->lats[0]) / ncinfo->yincr);
     *colptr = NINT((lon + .0001 - ncinfo->lons[0]) / ncinfo->xincr);

   if ((*rowptr < 0) || (*rowptr >= ncinfo->ny))
      return(-1);
   if ((*colptr < 0) || (*colptr >= ncinfo->nx))
      return(-1);
      
   return 0;
     
} /* end ncparms_rowcol() */

/*********************************************************************/ 

int ncparms_latlon( struct NC_INFO *ncinfo, int row, int col, float *latptr, float *lonptr)
{
   if ((row >= ncinfo->ny) || (row < 0) || (col >= ncinfo->nx) || (col < 0))
      return (-1);
      
     
      *latptr = ncinfo->lats[0] + row * ncinfo->yincr;
      *lonptr = ncinfo->lons[0] + col * ncinfo->xincr;
  
  
  return (0);
}

/*********************************************************************/ 

int ncparms_depth_indx(float theDepth, float *Z, int nz)
   /* Returns the standard level index corresponding to depth.  Or a
      negative number if some error occurs:
                 -1 :  depth not in Z 
                 -99:  depth exceeds greatest Z level.
   */
{
int d, sd, i;

   d = (int) theDepth + 0.00001;
   i = -1;
   
   while (++i < nz) {

      if (d < (sd = (int) (Z[i] + 0.00001)))
          return (-1);

      if (d == sd)
          return(i);
   }

/* if this statement is reached, the depth exceeds greatest standard depth 
     defined */
   return (-99);
 
  return(0);
}

/*********************************************************************/ 
void ncparms_zero(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)

/* initializes all fields in these structures */

{
     ncinfo->nz = 0;
     ncinfo->nz_seas = 0;
     ncinfo->fill_value = 0;
     ncinfo->mask_value = 0;
     ncinfo->modelCode = 0;
     ncdata->nz = 0;
     ncdata->nz_seas = 0;
     ncdata->dZdx = 0;
     ncdata->dZdy = 0;
    
     ncinfo->lats = NULL;
     ncinfo->lons = NULL;
     ncinfo->x_units = NULL;
     ncinfo->y_units = NULL;
     ncinfo->z_units = NULL;

     ncdata->depths = NULL;
     ncdata->Lx = NULL;
     ncdata->Ly = NULL;
     ncdata->Tparm0 = NULL;
     ncdata->Tparm1 = NULL;
     ncdata->Tparm2 = NULL;
     ncdata->Sparm0 = NULL;
     ncdata->Sparm1 = NULL;
     ncdata->Sparm2 = NULL;
     ncdata->Pparm0 = NULL;
     ncdata->Pparm1 = NULL;
     ncdata->Pparm2 = NULL;
     ncdata->dTdx = NULL;
     ncdata->dTdy= NULL;
     ncdata->dPdx = NULL;
     ncdata->dPdy = NULL;
     
     
     return;
 }
/*********************************************************************/ 
void ncparms_alloc(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)

/* allocates memory for structure fields */

{
  
     ncinfo->lats = (float *)calloc(ncinfo->ny, sizeof(float));
     ncinfo->lons = (float *)calloc(ncinfo->nx, sizeof(float));
     ncdata->depths = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Lx = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Ly = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Tparm0 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Tparm1 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Tparm2 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Sparm0 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Sparm1 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Sparm2 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Pparm0 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Pparm1 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->Pparm2 = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->dTdx = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->dTdy= (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->dPdx = (float *)calloc(ncinfo->nz, sizeof(float));
     ncdata->dPdy = (float *)calloc(ncinfo->nz, sizeof(float));
     if (ncdata->dPdy == NULL) {
         fprintf(stderr,"Unable to allocate memory in ncparms_alloc()\n");
	 exit(1);
     }
          
     return;
 }
/*********************************************************************/ 
void ncparms_prefill(struct NC_PARMS *ncdata, float fill_value)
   /* fills parameter arrays with specified value */
{
   int ilev;
   
   for (ilev = 0; ilev < ncdata->nz; ++ilev) {
      
     ncdata->Lx[ilev] = fill_value;
     ncdata->Ly[ilev] = fill_value;
     ncdata->Tparm0[ilev] = fill_value;
     ncdata->Tparm1[ilev] = fill_value;
     ncdata->Tparm2[ilev] = fill_value;
     ncdata->Sparm0[ilev] = fill_value;
     ncdata->Sparm1[ilev] = fill_value;
     ncdata->Sparm2[ilev] = fill_value;
     ncdata->Pparm0[ilev] = fill_value;
     ncdata->Pparm1[ilev] = fill_value;
     ncdata->Pparm2[ilev] = fill_value;
     ncdata->dTdx[ilev] = fill_value;
     ncdata->dTdy[ilev] = fill_value;
     ncdata->dPdx[ilev] = fill_value;
     ncdata->dPdy[ilev] = fill_value;
   }

   return;
}
/*********************************************************************/ 
void ncparms_free(struct NC_INFO *ncinfo, struct NC_PARMS *ncdata)

/* Frees memory for structure fields */

{

     free(ncinfo->lats);
     free(ncinfo->lons);
     free(ncdata->depths);
     free(ncdata->Lx);
     free(ncdata->Ly);
     free(ncdata->Tparm0);
     free(ncdata->Tparm1);
     free(ncdata->Tparm2);
     free(ncdata->Sparm0);
     free(ncdata->Sparm1);
     free(ncdata->Sparm2);
     free(ncdata->Pparm0);
     free(ncdata->Pparm1);
     free(ncdata->Pparm2);
     free(ncdata->dTdx);
     free(ncdata->dTdy);
     free(ncdata->dPdx );
     free(ncdata->dPdy);
     
     ncinfo->lats = NULL;
     ncinfo->lons = NULL;
     ncdata->depths = NULL;
     ncdata->Lx = NULL;
     ncdata->Ly = NULL;
     ncdata->Tparm0 = NULL;
     ncdata->Tparm1 = NULL;
     ncdata->Tparm2 = NULL;
     ncdata->Sparm0 = NULL;
     ncdata->Sparm1 = NULL;
     ncdata->Sparm2 = NULL;
     ncdata->Pparm0 = NULL;
     ncdata->Pparm1 = NULL;
     ncdata->Pparm2 = NULL;
     ncdata->dTdx = NULL;
     ncdata->dTdy = NULL;
     ncdata->dPdx = NULL;
     ncdata->dPdy = NULL;
     
     return;
 }  /* end ncparms_free() */


