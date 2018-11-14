/*   hb_gettrack.c
................................................................................
                              *  HydroBase3 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             July 2000
			     updated to HydroBase3 Jan 2012 
................................................................................
/*  gettrack takes a starting and ending position, computes the bearing.
    If a distance parameter is specified, station positions will be 
    generated along the trackline at that spacing.   
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hb_grids.h"


/* prototypes for locally defined functions */

void print_usage(char *);
int vector(double, double, double, double, double *, double *);
void output_pos(FILE *, double, double, int, int, int, int, int);
 
main (int argc, char **argv)
{
int i;
int stano, sflag,hflag, dflag;
int error, kilo, mflag;
int hdg;
double lat, lon, phi, dist_rem;
double start_lat, start_lon;
double end_lat, end_lon, co;
double prev_lat, prev_lon, dist;
double dx, dy;
char *st;
FILE  *outfile;
   

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
  
/*  set these default values */
 
    outfile = stdout;
    prev_lat = prev_lon = -999.0;
    error = 0;
    stano = 1;
    kilo = 0;
    dflag = mflag = sflag = hflag = 0;
     
for (i = 1; !error && i < argc; i++) {
   if (argv[i][0] == '-') {
     switch (argv[i][1]) {
        case 'A':
          st = &argv[i][2];
          if (*st == 'm') {
             ++st;
             error = sscanf(st,"%lf/%lf/%lf/%lf", &lat, &dy, &lon, &dx) == 4? 0:1;
             
             start_lat = ABS(lat) + dy/60;
             if (lat < 0)
                start_lat = -start_lat;
                
             start_lon = ABS(lon) + dx/60;
             if (lon < 0)
                start_lon = -start_lon;
             break;
          }
          error = sscanf(st,"%lf/%lf", &start_lat, &start_lon) == 2? 0:1;
          break;
        case 'B':
          st = &argv[i][2];
          if (*st == 'm') {
             ++st;
             error = sscanf(st,"%lf/%lf/%lf/%lf", &lat, &dy, &lon, &dx) == 4? 0:1;
             
             end_lat = ABS(lat) + dy/60;
             if (lat < 0)
                end_lat = -end_lat;
                
             end_lon = ABS(lon) + dx/60;
             if (lon < 0)
                end_lon = -end_lon;
             break;
          }
          error = sscanf(st,"%lf/%lf", &end_lat, &end_lon) == 2? 0:1;
          break;
        case 'D':
          dflag = 1;
          error = sscanf(&argv[i][2],"%lf", &dist) == 1 ? 0 : 1;
          break;
        case 'H':
          hflag = 1;
          break;
        case 'K':
          kilo = 1;
          break;
        case 'M':
          mflag = 1;
          break;
        case 'S':
          sflag = 1;
          if (argv[i][3] != '\0') {
             error = sscanf(&argv[i][2],"%d", &stano) == 1 ? 0 : 1;
          }
          break;
        case 'h':  
	  print_usage(argv[0]);
	  exit(0);

       default:
          error = 1;
          break;
     }  /* end switch */
   }
   else 
      error = 1;
      
} /* end for */
   
if (error) {
   fprintf(stderr,"\nError parsing command line");
   fprintf(stderr,"\n in particular: %s\n", argv[i]);
   exit(1);
}  

/*  start by getting bearing from point1 -> point2 
    hdg is degrees 0-360;  phi is in radians from north (-pi/2 < phi < pi/2) 
*/

   hdg = vector(start_lat, start_lon, end_lat, end_lon, &phi, &dist_rem);
   
/* write hdg and distance on screen */
   
   fprintf(stderr, "\n\n Bearing: %d  Distance: ", hdg);
   if (kilo) 
      fprintf(stderr," %10.2lf km\n", dist_rem/NMperKM);
   else 
      fprintf(stderr," %10.2lf nm\n", dist_rem);
   
   /* write out start position to outfile */
   output_pos(outfile, start_lat, start_lon, stano, hdg, sflag, hflag, mflag);
   
   if (dflag) {   /* determine station positions intermediate to end points */
   
      if (kilo)
        dist = dist * NMperKM;   /* convert to nautical miles */

      prev_lat = start_lat;
      prev_lon = start_lon;
      dx = ABS(dist* sin(phi));  
      dy = ABS(dist* cos(phi));
      if (hdg > 180)
         dx = -dx;
      if ((hdg > 90.) && (hdg < 270))
         dy = -dy;
         
      while (dist_rem > dist) {
        ++stano;
        co = cos(prev_lat * RADperDEG);
        lon = prev_lon + dx / (RADperDEG * EarthRadius * co);
        lat = prev_lat + dy / (RADperDEG * EarthRadius);
        
        /*write out new position */
        output_pos(outfile, lat, lon, stano, hdg, sflag, hflag, mflag);
        
        prev_lat = lat;
        prev_lon = lon;
        dist_rem -= dist;
      } /* end while */
     
   } /* end if dflag */
   
/* write out end position */
++stano; 
output_pos(outfile, end_lat, end_lon, stano, hdg, sflag, hflag, mflag);

fprintf(stderr,"\nEnd of gettrack.\n"); 
  
exit(0);
} /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nComputes course heading from pointA -> pointB ");       
   fprintf(stderr,"\n and positions at evenly spaced intervals along track.");
   fprintf(stderr,"\nOutputs  station # (optional), lat, lon, hdg "); 
   fprintf(stderr,"\n\nUsage:  %s -A[m]<lat/lon> -B[m]<lat/lon>  [-D<distance>] [-H] [-K] [-M] [-S<start_sta>] [-h]", program);
   fprintf(stderr,"\n-A :  lat/lon of starting point. If m is specified, lat and");
   fprintf(stderr,"\n         lon are each deg/mins");
   fprintf(stderr,"\n-B :  lat/lon of end point. If m is specified, lat and");
   fprintf(stderr,"\n         lon are each deg/mins");
   fprintf(stderr,"\n\n  OPTIONS:");
   fprintf(stderr,"\n-D :  distance between stations.");
   fprintf(stderr,"\n-H :  output the bearing from pointA to pointB.");
   fprintf(stderr,"\n-K :  distances are in km (default is nm).");
   fprintf(stderr,"\n-M :  output positions in deg/mins (default is decimal degrees).");
   fprintf(stderr,"\n-S :  output station # in 1st column sequentially starting with");
   fprintf(stderr,"\n         station specified.  ex:  -S1");
   fprintf(stderr,"\n-h help...... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
void output_pos(FILE *ofile,double  lat, double lon, int sta, int hdg, int sflag, int hflag, int mflag)
{
   int deg;
   double mins;

   if (sflag) 
     fprintf(ofile,"%5d ",sta);
  
   if (mflag) {
      deg =  (int)(fabs(lat) + .0001);
      mins = 60 * (fabs(lat) - (double) deg);
      fprintf(ofile,"%4d %4.1lf ",deg, mins);
      if (lat < 0)
           fprintf(ofile,"S ");
      else
         fprintf(ofile,"N ");
      
      deg =  (int)(fabs(lon) + .0001);
      mins = 60 * (fabs(lon) - (double) deg);
      fprintf(ofile,"%4d %4.1lf ",deg, mins);
      if (lon < 0)
           fprintf(ofile,"W ");
      else
         fprintf(ofile,"E ");
   }
   else 
      fprintf(ofile,"%8.2lf %8.2lf ", lat, lon);
  

   if (hflag) 
     fprintf(ofile,"%5d",hdg);
  
   fprintf(ofile,"\n");
   return;  
} /* end output_pos() */
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
      
/****************************************************************************/
 
