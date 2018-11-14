/*  hb_getdist.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated Feb 2000 to ANSI standards
			     update for HydroBase3 Mar 2012
................................................................................
................................................................................
.  Determines cumulative distance for set of stations and writes 
   (station #, distance) pairs to the stdout device. Optional starting
   distance can be specified. Optional lat/lon columns can be added 
   to enable plotting lat or lon along distance axis.
.
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include "hydrobase.h"


#define    EXTENT    ""
#define    DIR       ""
#define    PRINT_MSG  1               /* 0 or 1 */

#define    PI     3.141592654
#define    RADperDEG 0.017453292             /* pi/180 */
#define    EarthRadius  3437.747    /* in nm */
#define    NMperKM  .539593
#define    KMperNM   1.853248652

#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))

/* prototypes for locally defined functions */

void print_usage(char *);
int vector(double, double, double, double, double *, double *);

main (int argc, char **argv)
{
   int      i, nfiles, curfile; 
   int      error, status;
   int      kilo, latout, lonout;
   char   *dir, *extent;
   double  hdg, lat, lon, phi, cumdist, dist;
   double prev_lat, prev_lon;
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA data;
   FILE *infile;



/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    nfiles = 0;
    curfile = 1;
    cumdist = 0.0;
    prev_lat = -999.0;
    infile = stdin;
    latout = lonout= 0;


/* parse command line arguments... */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent = &argv[i][2];
                        break;
	       case 'K':
		        kilo = 1;
			break;
	       case 'L':
		        latout = 1;
			lonout = 1;
			switch (argv[i][2]) {
			   case 'a':
			      lonout = 0;
			      break;
			   case 'o':
			      latout = 0;
			      break;
			   default:
			      break;
			}
			break;
               case 'P':                    /* set starting position */
                        error = (sscanf(&argv[i][2],"%lf/%lf", &prev_lat, &prev_lon) != 2);
                        break;
               case 'S':                    /* set starting distance */
                        error = (sscanf(&argv[i][2],"%lf", &cumdist) != 1);
                        break;
               case 'h':
                        print_usage(argv[0]);
			exit(0);
                        break;
               default :
                        print_usage(argv[0]);
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */

            if (error ) {
                print_usage(argv[0]);
                fprintf(stderr,"\nError parsing command line args.\n");
                fprintf(stderr,"     in particular: '%s'\n", argv[i]);
                exit(1);
            }
            
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */


/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   hdr.prop_id = (int *) NULL;

    
    
/* loop for each file */

   do {
   if ( !nfiles) {
      fprintf(stderr,"\n Expecting data from stdin....  ");
   }
   else {
      infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
      if (infile == NULL) 
       goto NEXTFILE;
   }
   
     /* loop for each station */

     
     while ((status = get_station(infile, &hdr, &data)) == 0) { 
          if (prev_lat > -91.0) {
	       lat = (double)hdr.lat;
	       lon = (double)hdr.lon;
	       if (lon * prev_lon < 0)  { /* change in sign,check for crossed dateline */
	          if (ABS(lon - prev_lon) > 180.0) 
		       prev_lon -= 360.0;
		  else
		       prev_lon += 360.0;	       
	       }
               hdg = vector(prev_lat, prev_lon, lat, lon, &phi, &dist);
               if (kilo)
                   dist *= KMperNM;
	       
               cumdist += dist;
	  }
	       
	       fprintf(stdout, "%5d %10.2lf ", hdr.station, cumdist);
	       if (latout)
	          fprintf(stdout, "%8.3f ", hdr.lat);
	       if (lonout)
	          fprintf(stdout, "%8.3f ", hdr.lon);
	       fprintf(stdout, "\n");
	       
	       prev_lat = (double)hdr.lat;
	       prev_lon = (double)hdr.lon;

     }  /* end while */ 

     report_status(status, stderr);
     fclose(infile);

NEXTFILE:
      ;
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n\n End of hb_getdist. \n");
   exit(0);

}  /* end of main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nDetermines cumulative distance for set of stations");
   fprintf(stderr,"\n    and writes [station #, distance] pairs to the stdout device.\n");
   fprintf(stderr,"\n   An optional starting offset can be specified.  \n");
   fprintf(stderr,"\nUsage:  %s list_of_filenames", program);

   fprintf(stderr," [-Ddirname] [-Eextent] [-K] [-L<a|o>] [-P] [-S<startdist>] ");
   fprintf(stderr,"\n If no infiles are specified, input is expected from stdin");
   fprintf(stderr,"\n    -D  : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -K  : output distance in km (default is nm)");
   fprintf(stderr,"\n    -L  : output lat (-La), lon (-Lo) or both (-L)");
   fprintf(stderr,"\n           to enable plotting position along distance axis");
   fprintf(stderr,"\n    -P  : specify optional starting lat/lon to compute cumulative distance (default is first station in file)");
   fprintf(stderr,"\n    -S  : specify starting distance.  (default = 0)");
   fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
int vector(double lat1, double lon1, double lat2, double lon2, double *phi_addr, double *dist_addr)
  /* Returns the direction from point1->point2 in degrees from north.
     The arctan of this direction in radians  (domain is -pi/2 to +pi/2)
     is returned at phi_addr. The distance in nautical miles is returned
     at dist_addr */
{
   double dx, dy;
   
   dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG*.5 *(lat1+lat2)) * EarthRadius ;
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   
   if (dy == 0) {  /* course is zonal */
     *dist_addr = dx;
     if (dx < 0) {
       *phi_addr = -PI/2;
       *dist_addr = -dx;
       return(270);
     }
     *phi_addr = PI /2;
     return (90);
   }
   
   if (dx == 0.) {   /* meridional */
     *dist_addr = dy;
     *phi_addr = 0.0;
     if (dy < 0) {
       *dist_addr = -dy;
       return(180);
     }
     return(0);
   }

   *phi_addr = atan(dx/dy);
   *dist_addr = sqrt(dx*dx + dy*dy);
   
   if (*phi_addr > 0) {   
      if (dx > 0)
         return (NINT(*phi_addr / RADperDEG));   /* 0 -> 90 */
      
      return (180 + NINT(*phi_addr / RADperDEG)); /* 180 -> 270 */
   }
   
   if (dx > 0) 
     return (180 + NINT(*phi_addr / RADperDEG));  /* 90 -> 180 */
     
   return (360 + NINT(*phi_addr / RADperDEG));  /* 270 -> 360 */

}  /* end vector() */
      
