/*  hb_ncinfo.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             original 1993
			     Updated to ANSI Nov 1999
			     updated for HydroBase3M Sep 2009
................................................................................
.
.  Formerly hb_cdfinfo 
.  Reads a netCDF file created by the HydroBase 3 routines and prints a
.  summary of information about the file.
.
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include "hydro_cdf.h"


main (int argc, char **argv)
{
   int cdfid, i;
   int error, print_mess = 1;
   struct CDF_HDR cdf;
   void print_usage(char *);

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
   if (argv[1][0] == '-' && argv[1][1] == 'h') {
	   print_usage(argv[0]);
	   exit(0);
   }

   if ((cdfid = cdf_open("", argv[1], "", print_mess)) < 0)
      exit(1);

   if (error = read_cdf_hdr(cdfid, &cdf)) {
      fprintf(stderr,"\nError reading CDF file.\n");
      exit(1);
   }

   fprintf(stdout,"\n\nGrid Bounds...");
   fprintf(stdout,"\n    latitude: %f %f  increment: %f", cdf.ymin, cdf.ymax, cdf.yincr);
   fprintf(stdout,"\n    longitude: %f %f  increment: %f", cdf.xmin, cdf.xmax, cdf.xincr);
   fprintf(stdout,"\n    nrows: %d  ncols: %d ", cdf.ny, cdf.nx);
   fprintf(stdout,"\n    first gridnode centered at (lon, lat):");
   if (cdf.node_offset > 0)
       fprintf(stdout," (%f, %f)", (cdf.xmin + .5 * cdf.xincr), 
                               (cdf.ymax - .5 * cdf.yincr));
   else
       fprintf(stdout," (%f, %f)", cdf.xmin, cdf.ymax );

   fprintf(stdout,"\n\nYears :    ");
   for (i = 0; i < cdf.nt; ++i) {
      fprintf(stdout,"\n    %4d to %4d", cdf.tmin[i],cdf.tmax[i]);
   }

   fprintf(stdout,"\n\nProperties available...\n"); 
   for (i = 0; i < cdf.nprops; ++i)    
     fprintf(stdout,"  %2s ", cdf.prop_id[i]);  
   fprintf(stdout, "\n Counts for the above properties are ");
   if (!cdf.counts_included)  
      fprintf(stdout,"NOT ");
   fprintf(stdout,"included.");
   fprintf(stdout, "\n Error variances (sq root of) for the above properties are ");
   if (cdf.counts_included <= 1)  
      fprintf(stdout,"NOT ");
   fprintf(stdout,"included.");
   fprintf(stdout,"\n %d standard depths plus bottom depth at each lon/lat gridnode are included.", cdf.nz);

   fprintf(stdout,"\n\nCommand which generated this file ...");
   fprintf(stdout,"\n%s\n\n", cdf.command);

   exit(0);


} /* end main */


void print_usage(char *program)
{
   fprintf(stderr,"\n%s prints out various information about the \n", program);
   fprintf(stderr,"\ncontents of a HydroBase netcdf_file\n");
   fprintf(stderr,"\nUsage: %s nc_file_name ", program);
   fprintf(stderr, "\n       %s -h .... prints this message. \n", program);
   return;
} /* end print_usage() */
