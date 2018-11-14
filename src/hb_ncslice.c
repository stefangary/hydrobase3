/*   hb_ncslice.c
................................................................................
                              *  HydroBase3 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             May 2001 
			     updated for HB3M Dec 2009
................................................................................
/*  For a list of positions (lon/lat pairs) and set of HydroBase gridded 
    netcdf files, produces a HydroBase station file containing "stations" 
    that lie along the line connecting the positions.  In essence the
    positions define a pathway along which a vertical slice is cut
    through the 3D gridded property fields.

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_memory.h"
#include "hb_gamma.h"
#include "hb_paths.h"
#include "hb_grids.h"


/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""


struct POS_REC {
   double lat, lon, dist;
   short seafloor;
   int type;
};


/* globally defined variables */

int prop_req[2*MAXPROP];        /* list of props requested for output */
int nrequested;                /* count of above */	
int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
int tindex;           /* temperature variable to use for derived properties */
int lon0to360;   

struct GAMMA_NC ginfo;   /* used for neutral density */

float xoffset, yoffset;
float xinc, yinc;
int lflag;
double *latvec, *lonvec;
double lonmin, lonmax;
  
/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *);
int vector(double, double, double, double, double *, double *);
struct POS_REC *get_pathway(FILE *, int *);
int compare_points(const void *, const void *);
double getlat(double, double, double, double, double *);
double getlon(double, double, double, double, double *);
int check_bounds(struct CDF_HDR *, float, float);
int load_properties(int, struct CDF_HDR *, struct HYDRO_DATA *, float, float);
int get_profile(struct POS_REC *, int, char **, char *, char *, struct HYDRO_HDR *, struct HYDRO_DATA *);
void compute_props(struct HYDRO_HDR *, struct HYDRO_DATA *);


main(int argc, char **argv)
{
  int i, j,  n, npoints;
  int nfiles;
  int error, nz;
  int pixel_grid, quadrant;
  int iflag, pflag;
  char *dir, *extent;
  char *st;
  FILE *posfile, *outfile;
  struct POS_REC *point;
  struct HYDRO_DATA data;
  struct HYDRO_HDR hdr;
  struct GRID_INFO *tinfo;
  short *topo, tmiss;
   
/* Set these default values. */

  error = 0;
  dir = DIR;
  extent = EXTENT;
  nfiles = 0;
  outfile = stdout;
  pixel_grid = 1;
  iflag = pflag = lflag = 0;
  posfile = stdin;
  window = 100;
  w_incr = 10;

/* initialize some variables ... */

   hdr.prop_id = NULL;
   tindex = (int)T90;
   
   nrequested = 0;
   for (i = 0; i < MAXPROP; ++i) { 
      data.observ[i] = (double *) NULL;
      data.variance[i] = (double *) NULL;
      data.count[i] = (double *)NULL;
      data.quality[i] = (double *)NULL;
      prop_req[i] = 0;
      prop_needed[i] = 0;
   }   
      
/* set up header for creating HydroBase stations... */

   strncpy(hdr.country, "XX", 3);
   strncpy(hdr.ship, "XX", 3);
   hdr.cruise = 999;
   hdr.station = 999;
   hdr.year = 9999;
   hdr.month = 99;
   hdr.day = 99;
   hdr.qual[0] = '0';
   hdr.qual[1] = '0';
   hdr.qual[2] = '0';
   hdr.qual[3] = '0';
   hdr.instrument = 'u';
   hdr.origin = '0';
   
/*----------------------------------------*/  
  
/* are there command line arguments? */

   if (argc < 1) {
      print_usage(argv[0]);
      exit(1);
   }
   
/*----------------------------------------*/  
/* parse command line arguments */

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
 
        case 'D':                   /* get input dir */
          dir = &argv[i][2];
          break;
	  
        case 'E':                    /* get file extent */
          extent = &argv[i][2];
          break;
	  
        case 'G':        /* set to gridnode registration */
	  pixel_grid = 0;
	  break;
	  
        case 'I':
          iflag = 1;
          error = (sscanf(&argv[i][2],"%f", &xinc) == 1) ? 0 : 1;
          yinc = xinc;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%f", &yinc) == 1) ? 0 : 1;
               break;
            }
          }
                        
          break;
	  
         case 'L':
	  lflag = 1;
	  posfile = fopen(&argv[i][2],"r");
	  if (posfile == NULL) {
	    fprintf(stderr,"\nUnable to open %s for input.\n", &argv[i][2]);
	    exit(1);
	  }
	  break;
	  
       case 'O':
          outfile = create_hydro_file(&argv[i][2], OVERWRITE);
	  if (outfile == NULL) {
	    fprintf(stderr,"\nUnable to open %s for output.", &argv[i][2]);
	    exit(1);
	  }
          break;
	  
        case 'P':
          pflag = 1;
          nrequested = parse_p_option(&argv[i][2]);
          break;

        case 'T':
	    switch (argv[i][2]) {
	       case '6':
		  tindex = (int)TE;
		 break;
	      case '9':   
		  tindex = (int) T90;
		  break;
	      default:
	            fprintf(stderr,"Specify -T68 if you want IPTS-68 used in computing derived properties\n");
		    ++error;
	    }
	break;
	  
         case 'W':
            error = (sscanf(&argv[i][2],"%d", &window) == 1) ? 0 : 1;
           st = &argv[i][2];
           while (*(st++) != '\0') {
               if (*st == '/') {
                   ++st;
                   error = (sscanf(st,"%d", &w_incr) == 1) ? 0 : 1;
                   break;
               }
           }
                        
           break;

	case 'h':  
	   print_usage(argv[0]);
	   exit(0);
	   
        default:
          error = TRUE;
      } /* end switch */
      
       
      if (error ) {
         fprintf(stderr,"\nError parsing command line args.\n");
         fprintf(stderr,"     in particular: '%s'\n", argv[i]);
         exit(1);
      }
    }  /* end if */
    
    else {
      ++nfiles;
    }  /* end else */
    
  }  /* end for */
  
  
/*--------------------------------------------*/    
/*  Check syntax of arguments */ 

   error = 0;
   
    if (!nfiles ) {
       fprintf(stderr,"\nYou must specify input cdf_files as first argument(s).");
       ++error;    
    }

    if (!iflag ) {
       fprintf(stderr,"\nYou must specify grid spacing with -I<xincr>[/yincr] ");
      ++error;
    }
    
    if (!lflag ) {
       fprintf(stderr,"\nExpecting list of points from stdin: \n");
    }
    
    if (!pflag ) {
       fprintf(stderr,"\nYou must specify a list of properties for output: -Ppr/te/sa/ox ");
      ++error;
    }
    
   if (error) {
     fprintf(stderr,"\nUse -h for complete usage info. \n");
     exit(1);
   }

   xoffset = yoffset = 0.0;  /* if gridlines on whole degrees */
   
   if (pixel_grid) {         /* if not...*/
     xoffset = 0.5 * xinc;
     yoffset = 0.5 * yinc;
   }
    
/*--------------------------------------------*/    

   point = get_pathway(posfile, &npoints);
   
/*--------------------------------------------*/    

  /* add seafloor depths */
   
   tinfo = (struct GRID_INFO *)calloc(1, sizeof(struct GRID_INFO));
   tinfo->x_min = -180;
   tinfo->x_max = 180;
   tinfo->y_min = -90;
   tinfo->y_max = 90;
   topo = hb_get_topo(BATHPATH, tinfo, &latvec, &lonvec, FALSE, FALSE, &tmiss);
   if (topo != NULL) {
   
      for (i = 0; i < npoints; ++i) {
          if (point[i].lon <= 180) 
             point[i].seafloor = find_nearest_topo_val(point[i].lat, point[i].lon, topo, tinfo);
	  else
             point[i].seafloor = find_nearest_topo_val(point[i].lat, (point[i].lon-360.), topo, tinfo);
   
      }
      
      free(latvec);
      free(lonvec);
      free(topo);
      free(tinfo);
      
   } /* end if topo */

   lon0to360 = 0;
   lon0to360 = (lonmin < 180. && lonmax > 180.);
   
   if (lon0to360) {
      fprintf(stderr,"Dateline is crossed:  converting longitudes 0 -> 360");
   }

/*  for each point along pathway, obtain a profile from cdf files */

   for (i = 0; i < npoints; ++i) {
   
      iflag = get_profile(&point[i], nfiles, argv, dir, extent, &hdr, &data);
      
      if (iflag ) {
            
         if (hdr.nobs > 0) {
             compute_props(&hdr, &data);
	     write_hydro_station(outfile, &hdr, &data);
	 }
      
         free((void *)hdr.prop_id);
	 hdr.prop_id = NULL;
         free_hydro_data(&data);
      }
   } /* end for i */
   
   fprintf(stderr,"\n\nEnd of %s\n", argv[0]);
   exit(0);
}  /* end main */


/************************************************************************/
void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s cdf_file(s) -I<delta_x[/delta_y]>  -P<list_of_properties> [-D<dirname>] [-E<file_extent>] [-O<outfile>] [-G] [-L<position_file>] [-T68] [-W<window>[/<w_incr>]] ", program);

   fprintf(stderr,"\n\n  List of filenames must be first argument ");
   fprintf(stderr,"\n -I  : xgrid[/ygrid] increment in degrees.  ");
   fprintf(stderr,"\n -P  : list of properties for output;");
   fprintf(stderr,"\n       ex:  -Ppr/th9/sa/ox/ht");
   fprintf(stderr,"\n            -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n   OPTIONS:  ");
   fprintf(stderr,"\n[-D] : directory for input cdf files  ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : file extent for input cdf files");  
   fprintf(stderr,"\n            ex: -E.nc ");
   fprintf(stderr,"\n[-G] : Gridnodes fall on whole degrees in netcdf files");
   fprintf(stderr,"\n       default is gridlines offset by 0.5 * gridinc (pixel registration)");
   fprintf(stderr,"\n[-L] : file containing list of lon/lat positions.");
   fprintf(stderr,"\n       If not specified, this list will be read from stdin ");
   fprintf(stderr,"\n[-O] : name of output file");
   fprintf(stderr,"\n[-T] : temperature scale (IPTS-68 or ITS-90) to use in computing derived properties");
   fprintf(stderr,"\n       Default is -T90.  Specify -T68 as option. ");
   fprintf(stderr,"\n[-W] : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n       gradient properties (bf, pv...)  ");
   fprintf(stderr,"\n       defaults: -W%d/%d", window, w_incr);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */

/*****************************************************************************/
int parse_p_option(char *st)
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index, nprops, i;
   char prop[6];
   double ref_val;

   nprops = 0;

   do {
      if (*st == '/')
         ++st;
 	 
      for(i = 0; i < 6; i++)
	    prop[i] = '\0';
      sscanf(st,"%[^'/']", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      st += strlen(prop);
      prop_req[nprops] = index;
      if ( index < MAXPROP)
          prop_needed[index] = 1;
      ++nprops;
      
      

     /* !**!  Special cases for properties */

      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == PE)) {
         if (sscanf(st, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }
         while (!(*st=='/' || *st=='\0' || *st==' '))
            ++st;

        switch ((enum property) index) {
           case S_:
              s_pref = ref_val;
              break;
           case PE:
              pe_pref = ref_val;
              break;
           case HT:
              ht_pref = ref_val;
              break;
           default:
              ;
        }  /* end switch */
      }

   /* end of special cases */

  } while (*st == '/');
  
  if (prop_needed[(int)GE] || prop_needed[(int)GN]) 
       gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
       
  prop_needed[(int)PR] = 1;
  prop_needed[(int)DE] = 1;
  prop_needed[(int)SA] = 1;
  prop_needed[tindex] = 1;

  return (nprops);
}  /* end parse_prop_list() */
/***************************************************************************/
struct POS_REC *get_pathway(FILE *posfile, int *n_ptr)
  /* Reads the posfile and generates and returns an array of positions
     The number of points is returned at n_ptr */
{
  int nalloc, n, npoints;
  int startpoint;
  int hdg, zonal, meridional, quadrant;
  double lat0, lon0;
  double startlat, endlat, startlon, endlon;
  double phi, dist, startdist;
  struct POS_REC *point;

/*--------------------------------------------*/    
/* Read in first  lat/lon point */

   nalloc = 10000;

   point = (struct POS_REC *) get_memory((void *)NULL, (size_t)nalloc, sizeof(struct POS_REC));
   
   n = fscanf(posfile, "%lf%lf",  &startlon, &startlat);
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
   startdist = 0.0;
   
   if (startlon > lonmax) lonmax = startlon;
   if (startlon < lonmin) lonmin = startlon;

/*--------------------------------------------*/    
/* Get consecutive points, construct an array of gridline crossings between
   each pair of points */
          
   npoints = 1;
   startpoint = 0;  
   while  ((n = fscanf(posfile, "%lf%lf",  &endlon, &endlat)) == 2) {
      
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
	 
	    lat0 -=  1;
	    lat0 += yoffset;      /* find nearest y-gridline */
	      
	       /* find first gridline past startlat */  
	    while ((lat0 += yinc) <= startlat)  
	        ;  
		
		/* find all lat gridline xings until endlat */ 
	    while (lat0 < endlat) {
	       point[npoints].lat = lat0;
	       point[npoints].lon = getlon(lat0, startlon, startlat, phi, &point[npoints].dist);
	       point[npoints].dist += startdist;
	       point[npoints].type =  1;
	       if (++npoints > nalloc) {
	          nalloc += 100;
		  point = (struct POS_REC *) get_memory((void *)point, (size_t) nalloc, sizeof(struct POS_REC));
	       }
	       
	       lat0 += yinc;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading south */
	 
	                          /* find nearest whole degree */
	    lat0 +=  1;
	    lat0 -= yoffset;      /* find nearest y-gridline */
	      
	       /* find first gridline past startlat */  
	    while ((lat0 -= yinc) >= startlat)  
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
	       
	       lat0 -= yinc;
	    }
	   
       } /* end else */ 
     
     } /* end if !zonal */
     
     
/*--------------------------------------------*/    
     if (!meridional) { 
     
       /* get longitudinal gridline xings */
     
       lon0 = (double) NINT(startlon);
       if (quadrant == 0 || quadrant == 1) {
       
                /* heading east */
	 
	    lon0 -=  1;
	    lon0 += xoffset;      /* find nearest x-gridline */
	      
	       /* find first gridline past startlon */  
	    while ((lon0 += xinc) <= startlon)  
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
	       
	       lon0 += xinc;
	    }
	    
	}  /* end if */
	
	
	else {                 /* heading west */
	 
	    lon0 +=  1;
	    lon0 -= xoffset;      /* find nearest x-gridline */
	      
	       /* find first gridline past startlon */  
	    while ((lon0 -= xinc) >= startlon)  
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
	       
	       lon0 -= xinc;
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
   
   point = get_memory((void *)point, (size_t)npoints, sizeof(struct POS_REC));
   
   *n_ptr = npoints;
   return(point);

} /* end of get_pathway() */

/***************************************************************************/
int compare_points(const void *pt1, const void *pt2)
   /* Routine for qsort to sort structure by its .dist field.
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
     dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG * 0.5*(lat2+lat1)) * EarthRadius;
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
   
   return(dx / (EarthRadius * cos(RADperDEG*.5 *(lat1+lat2)) * RADperDEG) + lon1);
  
}
/****************************************************************************/
int get_profile(struct POS_REC *pos, int nfiles, char **arglist, char *dir, char *extent, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* Returns 1 if successful, 0 if not */
{
   int i, j, n, cdfid, error, row, col;
   int hdg, ip, nz, np2;
   int found, dcl;
   struct CDF_HDR cdf;
   float lat, lon;
   double phi, dist, dist2;
   double flag, weight1, weight2;
   double **propaddr1, **propaddr2;
   double *prop1, *prop2;
   struct HYDRO_DATA *dptr2;
     
   int print_msg = 0;
   int curfile = 1;
   
   
   flag = (double) -8.9;
   
   if (lon0to360) {      /* ensure longitude range 0 -> 360 */
      if (pos->lon < 0)
         pos->lon += 360.0;
   }
   else {                  /* ensure longitude range -180 -> 180 */
      if (pos->lon > 180.)
         pos->lon -= 360.0;
   }
      
   hptr->lat = (float) pos->lat;
   hptr->lon = (float) pos->lon;
   hptr->ms10 = ms10(hptr->lat, hptr->lon, &hptr->ms1);
   hptr->pdr = (int) pos->seafloor;
   
   do {
   
      cdfid = cdf_open(dir, arglist[curfile], extent, print_msg);
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
      found = check_bounds(&cdf, hptr->lat, hptr->lon);   
	 
      if (!found) 
         cdf_close(cdfid);
      
   }  while (!found && curfile++ < nfiles);
   
   if (!found)   
      return (0);

 /*  utilize this point ...	  */
 
   hptr->nprops = load_properties(cdfid, &cdf, dptr, hptr->lat, hptr->lon);
   
   hptr->prop_id = (int *) calloc((size_t)hptr->nprops, sizeof(int));
   n = 0;
   for (i = 0; i < MAXPROP; ++i) {
      if (dptr->observ[i] != NULL) {
         hptr->prop_id[n++] = i;
         if (dptr->count[i] != NULL) 
            hptr->prop_id[n++] = 500+i;
         if (dptr->variance[i] != NULL) 
            hptr->prop_id[n++] = 700+i;
      }
  } 
   
   nz = cdf.nz; 
    
   /* Is this gridnode at the same position as lat/lon requested? */
   
   get_indices(&cdf, hptr->lat, hptr->lon, &row, &col);
   get_lat_lon(&cdf, row, col, &lat, &lon);
   hdg = vector((double)lat, (double)lon, pos->lat, pos->lon, &phi, &dist);

   if (dist > 2 ) {   /* pt does not coincide with gridnode  */
   
	/* find another gridpt and interpolate the two */
	
	switch (pos->type) {
	   case 0:         /* longitude stays the same */
	      if (hdg > 90 && hdg < 180 ) {
		   lat = lat - cdf.yincr;
		   break;
              }
	      lat = lat + cdf.yincr;
	      break;
	   
	   case 1:        /* latitude stays the same */
	      if (hdg > 180) {
		   lon = lon - cdf.xincr;
		   break;
              }
	      lon = lon + cdf.xincr;
	      break;
	   
	   default:   /* unlikely case where endpt doesn't fall 
	                 on a gridnode.  */
	      
	      if ((hdg > 45 && hdg < 135) || (hdg > 225 && hdg < 315)) {
	         if (hdg > 180) {
		   lon = lon - cdf.xincr;
		   break;
                 }
	         lon = lon + cdf.xincr;
	         break;
              }
	      
	      if (hdg > 90 && hdg < 180 ) {
		   lat = lat - cdf.yincr;
		   break;
              }
	      lat = lat + cdf.yincr;
	   
	} /* end switch */
	
	/* first check currently open cdf file for match */
	
        found = check_bounds(&cdf, lat, lon);
	
	curfile = 0;
	while (!found && ++curfile <= nfiles) {
           cdf_close(cdfid);
           cdfid = cdf_open(dir, arglist[curfile], extent, print_msg);
           if (error = read_cdf_hdr(cdfid, &cdf)) 
               exit (1);
           found = check_bounds(&cdf, lat, lon);
	}   
	
	if (found) {

          if (cdf.nz != nz) {
	     fprintf(stderr,"\nFATAL ERROR:  mismatch of standard levels in cdf files.\n");
	     exit(1);
	  }	
	  dptr2 = (struct HYDRO_DATA *) get_memory((void *)NULL, 1, sizeof(struct HYDRO_DATA));
	
          np2 = load_properties(cdfid, &cdf, dptr2, lat, lon);
	  
	  if (np2 != hptr->nprops){
	     fprintf(stderr,"\nFATAL ERROR:  mismatch of properties in cdf files.\n");
	     exit(1);
	  }	

          /* assign relative weights to each profile inversely proportional
	     to distance */
	  	   
          get_indices(&cdf, lat, lon, &row, &col);
          get_lat_lon(&cdf, row, col, &lat, &lon);
          hdg = vector((double)lat, (double)lon, pos->lat, pos->lon, &phi, &dist2);

          weight1 = dist2 / (dist + dist2);
	  weight2 = dist / (dist + dist2);
	  
	  /* find deepest common depth and stdlevel */

          	  
	  dcl = (int) dptr->observ[(int)DE][nz-1];
	  np2 = 0;   /* keep track of which profile */
	  if (dcl < 0) {
	     dcl = dptr2->observ[(int)DE][nz-1];
	     np2 = 1;
	  }
          else {
	     if (dptr2->observ[(int)DE][nz-1] > 0)
	        if (dcl > dptr2->observ[(int)DE][nz-1]) {
	           dcl = (int) dptr2->observ[(int)DE][nz-1];
		   np2 = 1;
		}
	  }

          if (dcl < 0)
	    np2 = 2;	  
	  i = 0;
	  while ((i < nz-1)  && (std_depth[i] < dcl)) {
	     ++i;
	  }
	  
	  dcl = i;
	   
	   /* interpolate each level except bottom to produce a single profile */
	   
          for (i = 0; i < hptr->nprops; ++i) {
	    ip = hptr->prop_id[i];
	    prop2 = *get_prop_array(dptr2, ip);
	    
	    if (prop2 != NULL ) {
	       prop1 = *get_prop_array(dptr, ip);
               for (j = 0; j < nz-1; ++j) {
	         if (prop1[j] <= flag) 
		     prop1[j] = prop2[j];
		 else if (prop2[j] <= flag)
		     ;
		 else 
		    prop1[j] = weight1 * prop1[j] + weight2 * prop2[j];
	    
	       }
	       
	       /* set levels below deepest common level to empty */
	       for (j = dcl; j < nz-1; ++j) {
	             prop1[j] = (double) HBEMPTY;
	       }
	       
	       if (np2 == 1)
		  prop1[nz-1] = prop2[nz-1];
	       
	       if (np2 == 2)
		  prop1[nz-1] = HBEMPTY;
	    }
	  } 
	  
          free_hydro_data(dptr2);
	  free((void *) dptr2);
	}
	 
   }  /* end if dist > 2 */
   
   cdf_close(cdfid);
	      
    /* Prepare station for output.
       Define an index property to check for missing data */
   i = 0;
   ip = hptr->prop_id[i];
    while ((ip == (int)DE) || (ip == (int)PR) || (ip > MAXPROP) ) {
       if (++i == hptr->nprops)
           continue;
       ip = hptr->prop_id[i];
    }
      
    /* Remove flagged levels from data arrays */
    
   n = 0;
   prop1 = *get_prop_array(dptr, ip);  /* the index property */
   for (j = 0; j < nz; ++j) { 
      if (prop1[j] > flag) { 
         for (i = 0; i < hptr->nprops; ++i) {
	   prop2 = *get_prop_array(dptr,hptr->prop_id[i]);
           prop2[n] = prop2[j];
         }
	 ++n;
      } 
   }
   
   hptr->nobs = n;
   hptr->pdr = dptr->observ[(int)DE][n-1] + 10;
   if (pos->seafloor > hptr->pdr)
      hptr->pdr = pos->seafloor;  
   dptr->nobs = n;
   dptr->nprops = hptr->nprops;
   
   return (1);
} /* end get_profile() */

/****************************************************************************/
int check_bounds(struct CDF_HDR *cdfptr, float lat, float lon)

      /* Returns 1 if this point falls into the cdf bounds
         or 0 if not */
      
{
   float xoff, yoff;
      

      xoff = yoff = 0.0;      
      if (!cdfptr->node_offset) {
         xoff = 0.5 * cdfptr->xincr;
         yoff = 0.5 * cdfptr->yincr;
      }

     /* within bounds? */


      if ((lat < cdfptr->ymin -yoff) || (lat >= cdfptr->ymax) )   /* lat NOT in bounds */
          return (0);    
	            
      if ((lon >= (cdfptr->xmin - xoff))
        && (lon < (cdfptr->xmax + xoff)) )
        return (1);                                          /* lat,lon both in bounds */
	
	/* make sure longitude ranges are comparable */
	
       if (lon >= 0.) {          /* positive longitude */
      
          if (cdfptr->xmax < 0) {  /*lon range in cdf file is negative, make it pos */
             cdfptr->xmax += 360;
             cdfptr->xmin += 360;
	  }
	  else {   /* cdf file crosses greenwich, make lon neg and compare */
	      lon -= 360;
	  }
	  
         if ((lon >= (cdfptr->xmin - xoff)) && (lon < (cdfptr->xmax + xoff)) )
               return (1);                                          /* lat,lon both in bounds */
	  
	  return (0);  
      }
      
      /* negative longitude */
         if (cdfptr->xmin >= 0) {  /* lon range in cdf file is positive, make it neg */
             cdfptr->xmin -= 360;
             cdfptr->xmax -= 360;
	 }    
          else {   /* make lon positive and compare */
	     lon += 360;
	  }
	  
         if ((lon >= (cdfptr->xmin - xoff)) && (lon < (cdfptr->xmax + xoff)) )
               return (1);                                          /* lat,lon both in bounds */
	  
	  return (0);  
        

}  /* end check_bounds() */
/****************************************************************************/

int load_properties(int cdfid, struct CDF_HDR *cdfptr, struct HYDRO_DATA *dptr, float lat, float lon)

/*  Loads properties from cdf_file into dptr for gridpt nearest lat/lon.
   Returns number of properties, including depth and statistical properties. */

{
   int *prop_avail, *var_avail, *count_avail, error;
   int i, j, n, nprops, index;
   float *x;
   short *count;
   int row, col;
   char *mne;
   
   
  prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int));
  count_avail = (int *) calloc((size_t)MAXPROP, sizeof(int));
  var_avail = (int *) calloc((size_t)MAXPROP, sizeof(int));

  for (i = 0; i < cdfptr->nprops; ++i) {
    index = get_prop_indx(cdfptr->prop_id[i]);
    prop_avail[index] = 1;
    
    if (cdfptr->counts_included)
         count_avail[index] = 1;
      
    if (cdfptr->counts_included > 1) 
         var_avail[index] = 1;
      
  }
   
  get_indices(cdfptr, lat, lon, &row, &col);
  
  n = cdfptr->nz;
  x = (float *) get_memory((void *)NULL, n, sizeof(float));
  count = (short *) get_memory((void *)NULL, n, sizeof(short));

  free_and_alloc(&dptr->observ[(int)DE], n);
  
  error = read_cdf_depths(cdfid, x);
  error = read_cdf_bottom_depth(cdfid, &x[n-1], row, col, 0);       
  for (j = 0; j < n; ++j) {
     dptr->observ[(int)DE][j] = (double) x[j];
     std_depth[j] = (double) x[j];
  }
  
  nprops = 1;   /* depth  is included */
  
  for (i=0; i < MAXPROP; ++i) {
    if (prop_avail[i] ) {
       free_and_alloc(&dptr->observ[i], n);
       mne = get_prop_mne(i);
       error = read_cdf_prop(cdfid, mne, x, row, col, 0, 0, n);
       for (j = 0; j < n; ++j) {
           if (is_flagged(x[j], cdfptr->fill_value) || is_flagged(x[j],cdfptr->mask_value) )
	        x[j] = (float) HB_MISSING;
           dptr->observ[i][j] = (double) x[j];
       } 
       ++nprops;      
       
      if (count_avail[i]) {
          free_and_alloc(&dptr->count[i], n);
          error = read_cdf_prop_count(cdfid, mne, count, row, col, 0, 0, n);
          for (j = 0; j < n; ++j) {
	      dptr->count[i][j] = (double) count[j];  
           if (count[j] == 0)	        
	      dptr->count[i][j] = (double) HB_MISSING;
          }
          ++nprops;      
       }  
           
       if (var_avail[i]) {
          free_and_alloc(&dptr->variance[i], n);
            error = read_cdf_prop_err(cdfid, mne, x, row, col, 0, 0, n);
            for (j = 0; j < n; ++j) {
               if (is_flagged(x[j], cdfptr->fill_value) || is_flagged(x[j], cdfptr->mask_value))
	           x[j] = (float) HB_MISSING;
               dptr->variance[i][j] = (double) x[j];
            } 
            ++nprops;      
	}
    }
  }
  
 free(prop_avail);
 free(var_avail);
 free(count_avail);
 free(x);
 free(count);
 return(nprops);
 
} /* end load_properties()*/

/**********************************************************************************/
void compute_props(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* For all properties in the global array prop_req, computes it 
     if necessary . */
{
   int i, j, k, n, ratio_done;
   int main_props_avail, tavail;
   double dlat;
   double *e;
   int *oflag, *cflag, *vflag;
   
/*  get both temperature scales */
   
   tavail = available(TE, hptr) || available(T90,hptr);
   
   if (!available((int)TE, hptr) && tavail) 
       create_t68_arrays(hptr, dptr);
       
   if (!available((int)T90, hptr) && tavail)
       create_t90_arrays(hptr, dptr);
  
/* determine if pr, de, te, and sa are available ... */

    main_props_avail = 1;
    if (!(available(PR, hptr) && available(DE, hptr) && tavail 
              && available(SA, hptr))) {
         fprintf(stderr,"\n>>>>> WARNING!!!  ");
         fprintf(stderr,"Station does not include pr, de, te, and sa.");
         main_props_avail = 0;
    }
   
   ratio_done = 0;
   
/* !**! Special cases for individual properties... */
   
   for (k = 0; k < nrequested; ++k) {
   
      i = prop_req[k];
      if (!available((enum property)i, hptr) && main_props_avail) {
      
           switch ((enum property) i) {
             case OX:
                 if (available(O2, hptr)) 
		    create_ox_arrays(hptr, dptr);
               break;
             case O2:
                 if (available(OX, hptr)) 
                    create_o2_arrays(hptr, dptr);
               break;
             case TH:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_theta(hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)TE],dptr->observ[(int)SA] );
               break;
             case TH9:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_theta(hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[(int)T90],dptr->observ[(int)SA] );
               break;
             case S0:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(0., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case S1:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(1000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case S2:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(2000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case S3:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(3000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case S4:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(4000., hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case S_:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sigma(s_pref, hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case HT:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_height(hptr->nobs, dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], ht_pref, dptr->observ[i]);
               break;
             case PE:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_energy(hptr->nobs, dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], pe_pref, dptr->observ[i]);
               break;
             case SV:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sp_vol( hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;

             case VA:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_svan( hptr->nobs, dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA]);
               break;
             case DR:
	     case AL:
	     case BE:
	         if (! ratio_done) {
                    free_and_alloc(&dptr->observ[(int)DR], hptr->nobs);
                    free_and_alloc(&dptr->observ[(int)AL], hptr->nobs);
                    free_and_alloc(&dptr->observ[(int)BE], hptr->nobs);
                    compute_ratio( hptr->nobs, dptr->observ[(int)DR], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], dptr->observ[(int)AL], dptr->observ[(int)BE]);
		    ratio_done = 1;
		  }
               break;

             case VS:
                 free_and_alloc(&dptr->observ[i], hptr->nobs);
                 compute_sound_vel( dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], hptr->nobs);
               break;
             case PV:
               free_and_alloc(&dptr->observ[i], hptr->nobs);
               buoy_freq(dptr->observ[i], dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], hptr->nobs, window, w_incr);
	       
	       dlat = (double) hptr->lat;
	       
	       e = (double *) calloc(hptr->nobs, sizeof(double));
	       
	       for (j = 0; j < hptr->nobs; ++j) {
	         e[j] = dptr->observ[i][j];
	         if (e[j] > -999.0)
		     e[j] *= e[j];
	       }
	       
               po_vort(dptr->observ[i], e, hptr->nobs, dlat);
	       free((void *)e);
               break;
	        
             case BF:
               free_and_alloc(&dptr->observ[i], hptr->nobs);
               buoy_freq(dptr->observ[i],  dptr->observ[(int)PR], dptr->observ[tindex],dptr->observ[(int)SA], hptr->nobs, window, w_incr);
               break;

             case GN:
	         if (!prop_req[(int)GE]) {
                   free_and_alloc(&dptr->observ[GN], hptr->nobs);
	           compute_gamma_n(&ginfo, hptr->nobs, dptr->observ[GN], 
		       dptr->observ[(int)PR], dptr->observ[tindex],                                   dptr->observ[(int)SA], (double) hptr->lon, 
		       (double) hptr->lat);
		 }
	         break;
	       
             case GE:
                 free_and_alloc(&dptr->observ[(int)GE], hptr->nobs);
                 free_and_alloc(&dptr->observ[(int)GN], hptr->nobs);
	         compute_gamma_nerr(&ginfo, hptr->nobs, dptr->observ[(int)GN],   		     dptr->observ[(int)GE],  dptr->observ[(int)PR],
 		     dptr->observ[tindex], dptr->observ[(int)SA], 
		     (double) hptr->lon, (double) hptr->lat);
	         break;
             default:
               break;
          } /* end switch */
      }  /* end if */
   } /* end for k */
   
  /* now determine which properties were not requested, free the memory */
  
   oflag = (int *)calloc(MAXPROP,sizeof(int));
   cflag = (int *)calloc(MAXPROP,sizeof(int));
   vflag = (int *)calloc(MAXPROP,sizeof(int));

   n = 0;
   for (k= 0; k < nrequested; ++k) {
      if (prop_req[k] < MAXPROP) { 
          i = prop_req[k];
          if (dptr->observ[i] != NULL) {
	     oflag[i] = 1;
	     ++n;
	  }
	}  
      else if (prop_req[k] < 600) {
          i = prop_req[k] - 500;
	  if (dptr->count[i] != NULL) {
             cflag[i] = 1;
	     ++n;
	   }
      }
      else if (prop_req[k] < 800) {
          i = prop_req[k] - 700;
	  if (dptr->variance[i] != NULL) {
             vflag[i] = 1;
	     ++n;
          }
      }
   }

  
   for (i = 0; i < MAXPROP; ++i) { 
      if (dptr->observ[i] != NULL  && !oflag[i]) { 
         free(dptr->observ[i]);
	 dptr->observ[i] = NULL;
      }
      if (dptr->count[i] != NULL  && !cflag[i]) { 
         free(dptr->count[i]);
	 dptr->count[i] = NULL;
      }
      if (dptr->variance[i] != NULL  && !vflag[i]) { 
         free(dptr->variance[i]);
	 dptr->variance[i] = NULL;
      }
   }
   
   /* update hptr */
   
   free(hptr->prop_id);
   hptr->prop_id = (int *) calloc(n, sizeof(int));
   
   hptr->nprops = n;
   dptr->nprops = n;
   n = 0;
   for (i= 0; i < MAXPROP; ++i) {
      if (oflag[i])
         hptr->prop_id[n++] = i;
      if (cflag[i])       
         hptr->prop_id[n++] = i+500;
      if (vflag[i])       
         hptr->prop_id[n++] = i+700;
   }
   
   
   free(oflag);
   free(cflag);
   free(vflag);
   return;
}  /* end compute_props() */




