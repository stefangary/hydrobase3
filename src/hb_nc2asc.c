/*  hb_nc2asc.c

................................................................................
                          *******  HydroBase 3  *******
                 ..............................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated to ANSI-C standards Feb 2000
			     Updated to accommodate mean and error variance 
................................................................................
.
.  USAGE:
.  hb_nc2asc filename [-S] > outfile
.
.  Reads a .nc file created by HydroBase3 and creates an ascii file
.  in HydroBase profile format with each gridnode in the nc file forming
.  a pseudo-station.
.
................................................................................
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hydro_cdf.h"
#include "hydrobase.h"


int lon0to360;

main (int argc, char **argv)
{
   int cdfid, n, i, j, index, col, row, tbin;
   int error, print_mess, update_depth, fg_flag ;
   int index_prop, no_values, include_stats;
   struct CDF_HDR cdf;
   struct HYDRO_DATA station;
   struct HYDRO_HDR hdr;
   short *cnt;
   char *varname;
   float *z, *x, *v;
   float HB_f_mask;
   double *p, *t, *s;
   double **propaddr, *prop1;
   void print_usage(char *);
   FILE *outfile;

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }

/*  set these default values */

   tbin = 0;
   print_mess = 1;
   error = 0;
   cdfid = -1;
   update_depth = FALSE;
   include_stats = FALSE;
   fg_flag = FALSE;
   HB_f_mask = HBMASK;
   lon0to360 = 0;
   outfile = stdout;

/* parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {

               case 'D':
                        update_depth = TRUE;
                        break;
               case 'F':
                        fg_flag = 1;
                        break;
               case 'L':        /* force sign of longitude  */
	           switch (argv[i][2]) {
	             case '+':
	                lon0to360 = 1;
		        break;
	             case '-':
	                lon0to360 = -1;
	             case '0':
	                lon0to360 = 0;
		     break;
	             default:
	                 fprintf(stderr,"\nUse -L+  or -L- to force sign of  longitude");
	                 fprintf(stderr,"\nUse -L0 (zero) to use mixed sign longitudes (implies crossing Greenwich ");
	                 error = 1;
	           }
	  
	           break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			
			}
                        break;
	  
	      case 'S':
	          include_stats = 1;
	          break;
	   
	      case 'h':  
	                print_usage(argv[0]);
	                exit(0);

              default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             fprintf(stderr,"\nFor usage instructions, type %s -h\n\n", argv[0]);
             exit(1);
          }

       }  /* end if */

       else {
           if ((cdfid = cdf_open("", argv[i], "", print_mess)) < 0)
               exit(1);
       }

   }  /* end for */

   if (cdfid < 0) {
       fprintf(stderr,"\nYou must specify an input file!\n");
       fprintf(stderr,"\nFor usage instructions, type %s\n\n", argv[0]);
       exit(1);
   }
/* initialize some variables ... */

   for (i = 0; i < MAXPROP; ++i) 
      station.observ[i] = (double *) NULL;
   
/* read cdf file ... */

   if (error = read_cdf_hdr(cdfid, &cdf)) {
       if (error == 2) {
         fprintf(stderr, "\n Error in read_cdf_hdr():");
         fprintf(stderr, " Not a HydroBase cdf file!\n");
         exit (1);
       }
       else {
         fprintf(stderr, "\n Unknown Error in read_cdf_hdr():");
         fprintf(stderr, "    error code = %d\n", error);
         exit (1);
       }
   }

/* set up header for creating HydroBase stations... */

   strncpy(hdr.country, "XX", 3);
   strncpy(hdr.ship, "XX", 3);
   hdr.cruise = 999;
   hdr.station = 999;
   hdr.year = 9999;
   hdr.month = 99;
   hdr.day = 99;
   hdr.origin = '0';
   hdr.instrument = 'u';
   hdr.qual[0] = '0';
   hdr.qual[1] = '0';
   hdr.qual[2] = '0';
   hdr.qual[3] = '0';
   hdr.nprops  = cdf.nprops + 1;   /* add depth as a property */
   
   if (include_stats && cdf.counts_included) {
      hdr.nprops += cdf.nprops;        /* add nobs for each property */
      if (cdf.counts_included > 1 )
         hdr.nprops += cdf.nprops;     /* add error variance (_err) variables */
	 
      include_stats = cdf.counts_included;	 
   }
   station.nprops = hdr.nprops;
   hdr.prop_id = (int *) calloc((size_t)hdr.nprops, sizeof(int));

/* get property IDs, insert depth at top of list ...
   allocate memory to store the station data ... */

   hdr.prop_id[0] = (int) DE;
   station.observ[(int)DE] = (double *) calloc((size_t)cdf.nz, sizeof(double));
   i = 0;
   for (j = 0; j < cdf.nprops; ++j ) {
      index =  get_prop_indx(cdf.prop_id[j]);
      hdr.prop_id[++i] = index;
      station.observ[index] = (double *) calloc((size_t)cdf.nz, sizeof(double));
      if (include_stats) {
        varname = (char *)calloc(5, sizeof(char));
	strncpy(varname, "N", 1); 
        strcat(varname, cdf.prop_id[j]);
	index =  get_prop_indx(varname);
	if (index >= 0) {
	   hdr.prop_id[++i]= index;
	   propaddr = get_prop_array(&station, index);
           *propaddr = (double *) calloc((size_t)cdf.nz, sizeof(double));
	   if (include_stats == 3) {
	      varname[0] = 'V';
	      index = get_prop_indx(varname);
	      if (index >= 0) {
	         hdr.prop_id[++i]= index;
		 propaddr = get_prop_array(&station, index);
                 *propaddr = (double *) calloc((size_t)cdf.nz, sizeof(double));
	      }
	   }
	}
	free(varname);
      }
   } /* end for j */
   
   x = (float *) calloc((size_t)cdf.nz, sizeof(float));
   v = (float *) calloc((size_t)cdf.nz, sizeof(float));
   z = (float *) calloc((size_t)cdf.nz, sizeof(float));
   cnt = (short *) calloc((size_t)cdf.nz, sizeof(short));


    /* get the standard depths used in the cdf file */
    
   n = read_cdf_depths(cdfid, z);

/* determine the best property to determine whether any data is available
   at a given depth in .cdf file */

   if (station.observ[(int)TE] != NULL)
       index_prop = (int)TE;
   else if (station.observ[(int)T90] != NULL)
       index_prop = (int)T90;
   else if (station.observ[(int)TH] != NULL)
       index_prop = (int)TH;
   else if (station.observ[(int)SA] != NULL)
       index_prop = (int)SA;
   else if (station.observ[(int)PR] != NULL)
       index_prop = (int)PR;
   else
       index_prop = hdr.prop_id[1];


/* visit each gridpoint in the cdf file and create a station profile */

   for (row = 0; row < cdf.ny; ++row) {
      for (col = 0; col < cdf.nx; ++col) {

/* first get the various properties, and convert from float to double ... */

         no_values = 1;
	 j = 1;
         for (i = 0; i < cdf.nprops; ++i) {
            error = read_cdf_prop(cdfid, cdf.prop_id[i], x, row, col,
                     tbin, 0, cdf.nz);
            if (cdf.counts_included) {
               error = read_cdf_prop_count(cdfid, cdf.prop_id[i], cnt, row, col,
                     tbin, 0, cdf.nz);
	       if (include_stats > 1) {
                  error = read_cdf_prop_err(cdfid, cdf.prop_id[i], v, row, col,
                     tbin, 0, cdf.nz);
	       }
	    }
            for (n = 0; n < cdf.nz; ++n) {
	       if (fg_flag && cdf.counts_included) {
	           if (cnt[n] == -1) {
		      x[n] = cdf.fill_value;
		      if (include_stats > 1) 
		         v[n] = cdf.fill_value;
		   }
               }
	       station.observ[hdr.prop_id[j]][n] = (double) x[n];
               no_values = no_values && (is_flagged(x[n], cdf.fill_value) || is_flagged(x[n], HB_f_mask));
	       if (include_stats) {
	          station.count[hdr.prop_id[j]][n] = (double) cnt[n];
		  if (include_stats == 3)
	            station.variance[hdr.prop_id[j]][n] = (double) v[n];
		     
	       }
            } /* end for */
	    
	    ++j;
	    if (include_stats) {
	       ++j;
	       if (include_stats > 1) 
	          ++j;
	    }   
         }

         if (no_values)   
            continue;

         /* compute depth from pressure, if option was specified */
   
         if (update_depth && station.observ[(int)PR] != NULL) {
            for (j = 0; j < cdf.nz-1; ++j) {
	       if (station.observ[(int)PR][j] > -1. && station.observ[(int)PR][j] < 10000) 
                  z[j] = (float) hb_depth(station.observ[(int)PR][j], (double)hdr.lat);
            }
	 }

/* get bottom depth associated with this gridpoint and insert it into
   the depth array ... */

         error = read_cdf_bottom_depth(cdfid, &z[cdf.nz-1], row, col, tbin);
         hdr.pdr = (int) z[cdf.nz-1];
	 
/*   convert depth from float to double.  If p, t, s are present, eliminate
    scans which do not have a value for the index property. */
         n = 0;
         for (j = 0; j < cdf.nz; ++j) {
            station.observ[(int)DE][j] = (double) z[j];

            if (! (is_flagged(station.observ[index_prop][j], cdf.fill_value) || is_flagged(station.observ[index_prop][j],HB_f_mask)) ) {
                 for (i = 0; i < hdr.nprops; ++i) {
                    index = hdr.prop_id[i];
		    prop1 = *get_prop_array(&station, index);
                    prop1[n] = prop1[j];
		    if ((is_flagged(prop1[n], cdf.fill_value)) || (is_flagged(prop1[n], HB_f_mask))) {
		       prop1[n] = HB_MISSING;
		    }
                 }
                 ++n;
            }

         }  /* end for j */

         station.nobs = hdr.nobs =  n;

         if (hdr.nobs == 0)
            continue;


/*  fill in header elements specific to this gridnode ... */

         error = get_lat_lon(&cdf, row, col, &hdr.lat, &hdr.lon);
         hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
	 
	 if (lon0to360 > 0) {  /* make lons uniformly positive */
	    if (hdr.lon < 0.0)
	        hdr.lon += 360.0;
	}
	else if (lon0to360 < 0) {  /* make lons uniformly negative */
	   if (hdr.lon > 0.0)
	      hdr.lon -= 360.0;	
	}

         error = write_hydro_station(outfile, &hdr, &station);

      }  /* end for col */
   }  /* end for row */

   cdf_close(cdfid);
   fprintf(stderr, "\n End of hb_nc2asc.\n");

   exit(0);


} /* end main */


void print_usage(char *program)
{
   fprintf(stderr,"\n%s converts a binary HydroBase gridded netcdf file into an ascii", program);
   fprintf(stderr,"\nHydroBase station file where each gridnode becomes ");
   fprintf(stderr,"\na profile.");
   fprintf(stderr,"\n\nUsage: %s nc_file_name [-D] [-F] [-L+|-|0] [-O<outfile>] [-S] [-h] \n", program);
   fprintf(stderr,"    OPTIONS:\n");
   fprintf(stderr,"-D update depths from standard index by converting pressure to depth at each standard level \n");
   fprintf(stderr,"-F  omit levels which are flagged as being a first-guess value  (i.e. count = -1) \n");
   fprintf(stderr,"-L  force sign of longitude ) \n");
   fprintf(stderr,"        -L+   force positive longitudes ) \n");
   fprintf(stderr,"        -L-   force negative longitudes ) \n");
   fprintf(stderr,"       Default:  use mixed longitudes as stored in cdf file ) \n");
   fprintf(stderr,"\n   [-O] : output filename, default is to stdout");
   fprintf(stderr,"-S  include statistics (# of obs and error variance) if available \n");
   fprintf(stderr,"-h help...... prints this message. \n");
   return;
} /* end print_usage() */
