/*   hb_distdepth.c
................................................................................
                              *  HydroBase3 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             July 2000 
			     July 2010 updated for HB3
................................................................................
/*  hb_distdepth computes depth, distance between points, and distance along track for a list of station locations.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hb_paths.h"
#include "hb_memory.h"
#include "hb_grids.h"


/* global variables for topography grid*/
	
char *toponame;
int lon0to360, depth_is_neg;

/* prototypes for locally defined functions */

void print_usage(char *);
void output_pos(FILE *, double, double, int, int, short int, double, double, int, int);
 
main (int argc, char **argv)
{
int i;
int stano, sflag, lonfirst;
int error, kilo, moflag, miflag;
int hdg;
char lat_h, lon_h;
double lat_d, lat_m, lon_d, lon_m;
double lat, lon, phi, cumdist;
double start_lat, start_lon, temp;
double prev_lat, prev_lon, dist;
double dx, dy;
char *st;
char line[1000];
short int *ztopo, missing, depth;
double  *latvec, *lonvec;
struct GRID_INFO grid;
FILE  *outfile, *infile;
   

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
  
/*  set these default values */
 
    outfile = stdout;
    infile = stdin;
    prev_lat = prev_lon = -999.0;
    error = 0;
    stano= 1;
    kilo = 0;
    moflag = miflag = sflag = 0;
    lonfirst = 0;
    
/* define topography grid */
    grid.x_min = -180.0;
    grid.x_max = 180.0;
    grid.y_min = -90.0;
    grid.y_max = 90.0;
    lon0to360 = 0;
    depth_is_neg = 0; 
    toponame = BATHPATH;
     
for (i = 1; !error && i < argc; i++) {
   if (argv[i][0] == '-') {
     switch (argv[i][1]) {
        case 'I':
	   infile = fopen(&argv[i][2], "r");
	   if (infile == NULL) {
	      fprintf(stderr,"Unable to open %s for input.\n", &argv[i][2]);
	      exit(1);
	   }
          break;
        case 'K':
          kilo = 1;
          break;
	case 'L':
	   if (argv[i][2] == '+')
	       lon0to360 = 1;
	   else {
	      if (argv[i][2] == '-')
	         lon0to360 = 0;
	      else 
	         error = 1;
	   }
	   break;
        case 'M':
          switch (argv[i][2]) {
	    case 'i':
	       miflag = 1;
	       break;
	    case 'o':   /* fall through */
	    case '\0':
	       moflag = 1;
	       break;
	    default:
	       error = 1;
	  }
          break;
        case 'O':
          outfile = fopen(&argv[i][2],"w");
          if (outfile == NULL) {
                fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
                exit(1);
          }
          break;

        case 'S':
          sflag = 1;
          if (argv[i][2] != '\0') {
             error = sscanf(&argv[i][2],"%d", &stano) == 1 ? 0 : 1;
          }
          break;
        case ':':  
	  lonfirst = 1;
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
/* read in topography file ... */

    ztopo = hb_get_topo(toponame, &grid, &latvec, &lonvec, depth_is_neg, lon0to360, &missing);    

   if (ztopo == NULL)
      exit(1);
   
/*  Read first location and output it...*/

   fscanf(infile,"%[^\n]", line);
   fgetc(infile);
   
   if (miflag) {   /* degrees, minutes */
      if (lonfirst)  {
         if (sscanf(line, "%lf %lf %c %lf %lf %c", &lon_d, &lon_m, &lon_h, &lat_d, &lat_m, &lat_h)!= 6) {
	 
	    fprintf(stderr, "Expecting lon_deg lon_min hem lat_deg lat_min hem:\n%s\n", line);
	    exit(1);
	 }  
      }
      else  {
         if (sscanf(line, "%lf %lf %c %lf %lf %c", &lat_d, &lat_m, &lat_h, &lon_d, &lon_m, &lon_h) != 6) {
	 
	    fprintf(stderr, "Expecting lat_deg lat_min hem lon_deg lon_min hem:\n%s\n", line);
	    exit(1);
	 }
      } 
	 
      start_lon = lon_d + lon_m / 60.0;
	 if (lon_h == 'W' || lon_h == 'w')
	    start_lon = -start_lon;

      start_lat = lat_d + lat_m / 60.0;
	 if (lat_h == 'S' || lat_h == 's')
	    start_lat = -start_lat;
 
   }
   
   else {    /* decimal degrees */
   
      if (lonfirst)
         sscanf(line, "%lf %lf", &start_lon, &start_lat);
      else
         sscanf(line, "%lf %lf", &start_lat, &start_lon);
   }	 
      
   phi = 0;
   dist = 0.0;
   cumdist = 0.0;
   hdg = 0;
   
   depth =  find_nearest_topo_val(start_lat, start_lon, ztopo, &grid);

   if (sflag)
       fprintf(outfile, " Sta # "); 
   fprintf(outfile, "  Lat           Lon       Hdg  Depth(m)   Distance   Cum Dist ");
   if (kilo) 
      fprintf(outfile, "  (in km) \n\n");
   else
      fprintf(outfile, "  (in nm)  \n\n");
      
   
   output_pos(outfile, start_lat, start_lon, stano, hdg, depth, dist, cumdist, sflag, moflag);
  

/* loop for each position */

   prev_lat = start_lat;
   prev_lon = start_lon;
   
   while (fscanf(infile,"%[^\n]", line) == 1) {  
     
     fgetc(infile);
   
     if (miflag) {   /* degrees, minutes */
      if (lonfirst)  {
         if (sscanf(line, "%lf %lf %c %lf %lf %c", &lon_d, &lon_m, &lon_h, &lat_d, &lat_m, &lat_h)!= 6) {
	    fprintf(stderr, "Expecting lon_deg lon_min hem lat_deg lat_min hem:\n%s\n", line);
	    exit(1);
	 }  
      }
      else  {
         if (sscanf(line, "%lf %lf %c %lf %lf %c",&lat_d, &lat_m, &lat_h, &lon_d, &lon_m, &lon_h ) != 6) { 
	 
	    fprintf(stderr, "Expecting lat_deg lat_min hem lon_deg lon_min hem:\n%s\n", line);
	    exit(1);
	 }
      } 
	 
      start_lon = lon_d + lon_m / 60.0;
	 if (lon_h == 'W' || lon_h == 'w')
	    start_lon = -start_lon;

      start_lat = lat_d + lat_m / 60.0;
	 if (lat_h == 'S' || lat_h == 's')
	    start_lat = -start_lat;
 
      }
   
      else {    /* decimal degrees */

         if (sscanf(line, "%lf %lf", &start_lat, &start_lon) != 2) {
	    fprintf(stderr, "Error reading input positions\n");
	    exit(1);
    
         }
   
         if (lonfirst) {
          temp = start_lat;
	  start_lat = start_lon;
	  start_lon = temp;
         }
       
      }
       
       ++stano;
       dist = distance_c(prev_lat, prev_lon, start_lat, start_lon, kilo, &phi);
       hdg = NINT(phi);
       depth =  find_nearest_topo_val(start_lat, start_lon, ztopo, &grid);

   
       cumdist += dist;

       output_pos(outfile, start_lat, start_lon, stano, hdg, depth, dist, cumdist, sflag, moflag);
      
        prev_lat = start_lat;
        prev_lon = start_lon;
      
   }  /* end while fscanf */
   
   fprintf(stderr,"\nEnd of hb_distdepth.\n"); 
   exit(0);
} /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nDetermines distance between a list of locations and depth at each location");       
   fprintf(stderr,"\nOutputs  station # (optional), lat, lon, hdg, depth, distance, cumdist"); 
   fprintf(stderr,"\n\nUsage:  %s [-I{input_position_file>]   [-K] [-L<+/->]  [-M[o|i]] [-S<start_sta>] [-h]", program);
   fprintf(stderr,"\n\n  OPTIONS:");
   fprintf(stderr,"\n-I :  file containing input positions (lat/lon).");
   fprintf(stderr,"\n-K :  output distance in km (default is nm).");
   fprintf(stderr,"\n-L :  output West longitude either positive [-L+] or negative [-L-].");
   fprintf(stderr,"\n         default : W longitude is negative [-L-]");
   fprintf(stderr,"\n-M :   positions in deg/mins (default is decimal degrees).");
   fprintf(stderr,"\n       Use -M or -Mo for output positions");
   fprintf(stderr,"\n       Use -Mi for input positions");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
   fprintf(stderr,"\n-S :  output station # in 1st column sequentially starting with");
   fprintf(stderr,"\n         station specified.  ex:  -S1");
   fprintf(stderr,"\n-: :  input positions are lon lat.  default is lat lon");
   fprintf(stderr,"\n-h help...... prints this message. \n");
   fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
void output_pos(FILE *ofile, double  lat, double lon, int sta, int hdg, short int depth, double distance, double distsum, int sflag, int mflag)
{
   int deg;
   double mins;

   if (sflag) 
     fprintf(ofile,"%5d ", sta);
     
   if (lon0to360) {
      if (lon < 0.0)
         lon += 360.0;
      
   }
   else {
      if (lon > 180.)
          lon -= 360.0;
   }
  
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
   else {
      fprintf(ofile,"%12.3lf %12.3lf ", lat, lon);
   }
   
    fprintf(ofile,"%5d %8d  %10.2lf %12.2lf", hdg, depth, distance, distsum);
     
  
   fprintf(ofile,"\n");
   return;  
} /* end output_pos() */
/****************************************************************************/
