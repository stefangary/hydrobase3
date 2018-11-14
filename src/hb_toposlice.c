/* hb_toposlice.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2001
................................................................................
..................................................................  
.  Reads etopo1 bathymetry file supplied with HydroBase3
.  and writes out lon/lat/z triplets  corresponding to points along
.  pathway speci. Distance is an optional 4th column. For a pathway that is not 
.  zonal or meridional, the output will *not* be equally spaced, but instead 
.  correspond to points at the intersection of the pathway and the gridlines.
..................................................................  
*/ 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "hb_memory.h"
#include "hb_paths.h"
#include "hydrobase.h"
#include "hb_grids.h"

#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#define    ABS(x)       (((x) < 0) ? -(x) : (x))


struct POS_REC {
   double lat, lon, dist;
   int type;
};

   
/*  prototypes for locally defined functions */

void print_usage(char *);
int vector(double, double, double, double, double *, double *);
int compare_points(const void *, const void *);
double getlat(double, double, double, double, double *);
double getlon(double, double, double, double, double *);
double get_depth(struct POS_REC *, short **, double *, double *);
           
int main (int argc, char **argv)
{ 
   int i, j, n, npoints;
   int dflag, kflag, lflag, pflag;
   int mflag;
   int error, nalloc;
   int  output_neg_depth;   
   int hdg, zonal, meridional;
   int startpoint, quadrant;
   int xparm, topo_id;
   double zout, lat, lon, x;
   double phi, dist, startdist;
   double lat0, lon0, lonmin, lonmax;
   double startlat, endlat, startlon, endlon;
   double *latvec, *lonvec;
   double prevlat, prevlon;
   char *s, *toponame;
   FILE *outfile;
   FILE *posfile;
   struct POS_REC *point;
   struct GRID_INFO tinfo;
   short *topo, tmiss;

   if (argc < 1 ){
      print_usage(argv[0]);
   }
   
/* initialize these ... */

   pflag = 0;
   dflag = kflag = lflag = mflag = 0;
   output_neg_depth = 0; 
   error = 0;  
   outfile = stdout;
   posfile = stdin;
   toponame = BATHPATH_C;
   tinfo.x_min = -180.0;  /* get global data set */
   tinfo.x_max = 180.0;
   tinfo.y_min = -90.0;
   tinfo.y_max = 90.0;
   tinfo.lon0to360 = 0;
   
   
/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 
            case 'D':
	       dflag = 1;
	       ++s;
	       if (*s == 'k' || *s == 'K')
	          kflag = 1;
	       break;
           
	    case 'H':
               toponame = BATHPATH;
               break;

              
            case 'L':
	       lflag = 1;
	       posfile = fopen(&argv[i][2],"r");
	       if (posfile == NULL) {
	         fprintf(stderr,"\nUnable to open %s for input.\n", &argv[i][2]);
	         exit(1);
	       }
	       break;
	  
           case 'M':
	       mflag = 1;
               ++s;
               error += (sscanf(s,"%d", &xparm) != 1);
	       error += (xparm < 1) || (xparm > 3); 
	       break;
              
           case 'N' :
               output_neg_depth = 1;
               break;
	       
           case 'O':
               outfile = fopen(&argv[i][2],"w");
               if (outfile == NULL) {
                  fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                  exit(1);
               }
               break;
           case 'P':
	       pflag = 1;
	       break;
               
           case 'T':
               toponame = &argv[i][2];
               break;
               
            case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            default:  
               error = 1;
               
         } /* end switch */
             
      }
      else  {
        error = 1;
      }
      
      if (error ) { 
           fprintf(stderr,"\nError parsing command line args.\n");
           fprintf(stderr," in particular:  '%s'\n", argv[i]); 
           fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
           exit(1); 
      }
   } /* end for */

   topo = hb_get_topo(toponame, &tinfo, &latvec, &lonvec, output_neg_depth, tinfo.lon0to360, &tmiss);

    if (! lflag ) {
      fprintf(stderr,"Enter longitude latitude pairs ....\n  "); 
   } 

/*--------------------------------------------*/    
/* Read in first  lat/lon point */

   nalloc = 10000;

   point = (struct POS_REC *) get_memory((void *)NULL, (size_t)nalloc, sizeof(struct POS_REC));
   
   if (lflag)
       n = fscanf(posfile, "%lf%lf",  &startlon, &startlat);
   else 
       n = fscanf(stdin, "%lf%lf\n",  &startlon, &startlat);
       
   
   if (n != 2) {
      fprintf(stderr,"\nError reading startlon, startlat\n");
      exit(1);
   }
   
   point[0].lat = startlat;
   point[0].lon = startlon;
   point[0].dist = 0.0;
   point[0].type = 2;
   
   lonmin = 999.0;
   lonmax = -999.0;
   startdist = 0;
   
   if (startlon > lonmax) lonmax = startlon;
   if (startlon < lonmin) lonmin = startlon;

/*--------------------------------------------*/    
/* Get consecutive points, construct an array of gridline crossings between
   each pair of points */
          
   npoints = 1;
   startpoint = 0; 
   if (lflag)
      n = fscanf(posfile, "%lf%lf",  &endlon, &endlat);
   else
      n = fscanf(stdin, "%lf%lf\n",  &endlon, &endlat);
       
   while  (n == 2) {
      
      if (endlon > lonmax) lonmax = endlon;
      if (endlon < lonmin) lonmin = endlon;
     hdg = vector(startlat, startlon, endlat, endlon, &phi, &dist);
     dist += startdist;
     
     if (hdg >= 360) 
        hdg = 360 - hdg;
	
     meridional = zonal = 0;
     
     if (hdg == 180 || hdg == 0 )
        meridional = 1;
     if (hdg == 90 || hdg == 270)
        zonal = 1;
	
     if (hdg < 90)
        quadrant = 0;
     else if (hdg < 180)
        quadrant = 1;
     else if (hdg < 270)
        quadrant = 2;
     else 
        quadrant = 3;
    
     	
/*--------------------------------------------*/    
     if (!zonal) {
     
       /* get latitudinal gridline xings */
       
       lat0 = (double) NINT(startlat);
       if (quadrant == 0 || quadrant == 3) {
       
                /*heading north */
	 
	    if (lat0 > startlat)  /* find nearest whole degree */
	      lat0 -=  1;
	      
	      
	       /* find first gridline past startlat */  
	    while ((lat0 += tinfo.y_inc) <= startlat)  
	        ;  
		
		/* find all lat gridline xings until endlat */ 
	    while (lat0 < endlat) {
	       point[npoints].lat = lat0;
	       point[npoints].lon = getlon(lat0, startlon, startlat,  phi, &point[npoints].dist);
	       point[npoints].dist += startdist;
	       point[npoints].type =  1;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lat0 += tinfo.y_inc;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading south */
	 
	    if (lat0 < startlat)  /* find nearest whole degree */
	      lat0 +=  1;
	      
	      
	       /* find first gridline past startlat */  
	    while ((lat0 -= tinfo.y_inc) >= startlat)  
	        ;  
		
		/* find all lat gridline xings until endlat */ 
	    while (lat0 > endlat) {
	       point[npoints].lat = lat0;
	       point[npoints].lon = getlon(lat0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].dist += startdist;
	       point[npoints].type =  1;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lat0 -= tinfo.y_inc;
	    }
	   
       } /* end else */ 
     
     } /* end if !zonal */
     
     
/*--------------------------------------------*/    
     if (!meridional) { 
     
       /* get longitudinal gridline xings */
     
       lon0 = (double) NINT(startlon);
       if (quadrant == 0 || quadrant == 1) {
       
                /* heading east */
	 
	    if (lon0 > startlon)  /* find nearest whole degree */
	      lon0 -=  1;
	      
	      
	       /* find first gridline past startlon */  
	    while ((lon0 += tinfo.x_inc) <= startlon)  
	        ;  
		
		/* find all lon gridline xings until endlon */ 
	    while (lon0 < endlon) {
	       point[npoints].lon = lon0;
	       point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].dist += startdist;
	       point[npoints].type =  0;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lon0 += tinfo.x_inc;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading west */
	 
	    if (lon0 < startlon)  /* find nearest whole degree */
	      lon0 +=  1;
	      	      
	       /* find first gridline past startlon */  
	    while ((lon0 -= tinfo.x_inc) >= startlon)  
	        ;  
		
		/* find all lon gridline xings until endlon */ 
	    while (lon0 > endlon) {
	       point[npoints].lon = lon0;
	       point[npoints].lat = getlat(lon0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].dist += startdist;
	       point[npoints].type =  0;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lon0 -= tinfo.x_inc;
	    }
	   
       } /* end else */ 
     
     
     } /* end if !meridional */
     
/*--------------------------------------------*/    
     /* Add the endpoint to the array */
     
      point[npoints].lat = endlat;
      point[npoints].lon = endlon;
      point[npoints].dist = dist;
      point[npoints].type = 2;
      ++npoints;
      
      /* sort the array of points in order by increasing dist from startpoint */
            
      qsort((void *) &point[startpoint], npoints - startpoint, sizeof(struct POS_REC), compare_points);
      
      
      startlat = endlat;
      startlon = endlon;
      startpoint = npoints;
      startdist = dist;
      
      if (lflag)
         n = fscanf(posfile, "%lf%lf",  &endlon, &endlat);
      else
         n = fscanf(stdin, "%lf%lf\n",  &endlon, &endlat);
	 
   } /* end while sscanf(posfile) */
   
   if (lflag)
      fclose(posfile);
   
   if (n != EOF) {
     fprintf(stderr,"\n ERROR reading lon/lat pair.... continuing.\n");
   }
   
   if (npoints < 2) {
     fprintf(stderr,"\n You must provide at least 2 lon/lat end points.\n");
     fprintf(stderr,"\nUse -h for complete usage info. \n");
     exit(1);
   }
   
   point = (struct POS_REC *) get_memory((void *)point, (size_t)npoints, sizeof(struct POS_REC));

/*--------------------------------------------*/    
/*        end of defining pathway */
/*--------------------------------------------*/    

       
    fprintf(stderr,"\nWriting out values.....\n");

    prevlat = point[0].lat;
    prevlon = point[0].lon;
    if (xparm == 3) 
       dflag = 1;
       
    for (i = 0; i < npoints; ++i) {
        zout = find_nearest_topo_val(point[i].lat, point[i].lon, topo, &tinfo);
	   
	 if (pflag) 
	     zout = hb_p80( zout, point[i].lat);
	     
	if (dflag)  {
	    dist = point[i].dist;
	    if (kflag)
	       dist = dist * KMperNM;
	} 
	if (mflag) {
	  
	  x = point[i].lon;
	  if (xparm == 2)
	    x = point[i].lat;
	  if (xparm == 3)
	    x = dist;
	
	  if (i == 0) 
	     fprintf(outfile,">I\n%8.3lf  10000\n", x);
	     
	  fprintf(outfile,"%8.3lf %8.0lf\n", x, zout);
	   
	}
        else {
	   fprintf(outfile,"%8.3lf %8.3lf %8.0lf", point[i].lon, point[i].lat, zout);
	   if (dflag)  
	      fprintf(outfile," %10.2lf", dist);
	   fprintf(outfile,"\n");
	}
     }
     
     if (mflag) {
          --i;
	  x = point[i].lon;
	  if (xparm == 2)
	    x = point[i].lat;
	  if (xparm == 3)
	    x = dist;
	   fprintf(outfile,"%8.3lf 10000\n", x);
     }
     
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);

} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n%s outputs bathymetry along a pathway specified", program);
   fprintf(stderr,"\nas a list of lon/lat positions. The default bathymetry database");
   fprintf(stderr,"\nis %s ", BATHPATH_C);
   fprintf(stderr,"\nA higher resolution file <%s> is available with -H option.", BATHPATH);
   fprintf(stderr,"\nAn alternative file may be specified with the -T option. ");
   fprintf(stderr,"\nPathway must be a minimum of 2 points, but can be");  
   fprintf(stderr,"\nmultiple points. -L specifies filename containing these points.");  
   fprintf(stderr,"\nIf no posfile is specified, the points are expected from stdin.");
   fprintf(stderr,"\nLon/lat/z values are written to outfile or stdout.");
   fprintf(stderr,"\n -D causes the distance along the pathway to be written"); 
   fprintf(stderr,"\nout as a fourth column.");
   fprintf(stderr,"\nTo output a polygon mask file for plotting or masking topography,");
   fprintf(stderr,"\nspecify -M<xval> where xval is 1, 2, or 3 (lon, lat or distance)");
   fprintf(stderr,"\n\nUsage:  %s  [-T<global_topo_file>] [-D[k]] [-L<posfile>] [-M<xval>] [-N] [-O<outfile>] [-P]", program);
   fprintf(stderr,"\n\n   OPTIONS:"); 
   fprintf(stderr,"\n[-D] : include distance (nautical miles) in output. "); 
   fprintf(stderr,"\n       Append k for distance in km [-Dk]. "); 
   fprintf(stderr,"\n[-H] : Use high resolution topography file %s ", BATHPATH); 
   fprintf(stderr,"\n[-L] : position file contains lon/lat pairs (one per line) to define pathway"); 
   fprintf(stderr,"\n[-M] : outputs xy values for use as a mask file with hb_section."); 
   fprintf(stderr,"\n       specify xval = 1 for longitude"); 
   fprintf(stderr,"\n                    = 2 for latitude"); 
   fprintf(stderr,"\n                    = 3 for distance"); 
   fprintf(stderr,"\n       Info is added at beginning and end for polygon masking mode."); 
   
   fprintf(stderr,"\n[-N] : make seafloor values negative. ");
   fprintf(stderr,"\n          [default is seafloor positive] ");
   fprintf(stderr,"\n[-O] : output file. If not specified output goes to stdout"); 
   fprintf(stderr,"\n[-P] : output pressure instead of depth"); 
   fprintf(stderr,"\n[-T] : specify name of netcdf topography file. Default [%s]", BATHPATH_C);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");
   return;
   
} /*end print_usage() */
   
/****************************************************************************/
/***************************************************************************/
int compare_points(const void *pt1, const void *pt2)
   /* Routine for qsort (in <stdlib.h>) to sort structure by its .dist field.
      Returns -1, 0, or 1 for pt1 < pt2, pt1 == pt2, or pt1 > pt2 */
{
   struct POS_REC *p1, *p2;
   double dist1, dist2;
   
   p1 = (struct POS_REC *)pt1;
   p2 = (struct POS_REC *)pt2;
   dist1 = p1->dist;
   dist2 = p2->dist;
   
   if  (ABS(dist1 - dist2) <= 1) 
      return(0);
   
   if (dist1 > dist2)
      return(1);
    
    return(-1);
   
}
/**************************************************************************/
int vector(double lat1, double lon1, double lat2, double lon2, double *phi_addr, double *dist_addr)
  /* Returns the direction from point1->point2 in degrees from north.
     The arctan of this direction in radians  (domain is -pi/2 to +pi/2)
     is returned at phi_addr. The distance in nautical miles is returned
     at dist_addr */
{
   int hdg;
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
      if (dx > 0){
         hdg = (int) NINT(*phi_addr / RADperDEG);   /* 0 -> 90 */
	 return(hdg);
      }
      
      return ((180 + NINT(*phi_addr / RADperDEG))); /* 180 -> 270 */
   }
   
   if (dx > 0) 
     return (180 + NINT (*phi_addr / RADperDEG));  /* 90 -> 180 */
     
   return (360 + NINT(*phi_addr / RADperDEG));  /* 270 -> 360 */

}  /* end vector() */
      
/****************************************************************************/
double getlat(double lon2, double lon1, double lat1, double phi, double *dist_addr)
  /* Returns the latitude that intersects at lon2 the ray whose vertex
     is lon1,lat1 with angle phi in radians and domain -pi/2 to pi/2.
     The distance in naut miles between those points is also returned.  */
{
   double dx, dy, lat2;
   int i;
   
   lat2 = lat1;
   
   /* iterate a few times to converge on lat2 */
   
   for (i = 0; i < 10; ++i) {
     dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG * 0.5*(lat2+lat1) ) * EarthRadius;
     dy = dx / tan(phi);
     lat2 = dy /(RADperDEG * EarthRadius) + lat1;
   }

   *dist_addr = sqrt(dx*dx + dy*dy);
   return( dy /(RADperDEG * EarthRadius) + lat1);
   
}

/****************************************************************************/
double getlon(double lat2, double lon1, double lat1, double phi, double *dist_addr)
  /* Returns the longitude that intersects at lat2 the ray whose vertex
     is lon1,lat1 with angle phi in radians and domain -pi/2 to pi/2.
      The distance in naut miles between those points is also returned.    */
{
   double dx, dy;
   
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   dx =  dy * tan(phi);
   
   *dist_addr = sqrt(dx*dx + dy*dy);
   
   return(dx / (EarthRadius * cos(RADperDEG*0.5*(lat1+lat2)) * RADperDEG) + lon1);
  
}
