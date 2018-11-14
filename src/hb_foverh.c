/* hb_foverh.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             May 2001
................................................................................
..................................................................  
.  Reads input file(s)
.  and writes out lon lat f/H triplets

..................................................................  
*/ 
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_memory.h"
#include "hb_paths.h"

#define MISSING  -9.0e9

struct GRID_INFO {
	int nx;			/* Number of columns */
	int ny;			/* Number of rows */
	int node_offset;	/* 0 for node grids, 1 for pixel grids */
	double x_min;		/* Minimum x coordinate */
	double x_max;		/* Maximum x coordinate */
	double y_min;		/* Minimum y coordinate */
	double y_max;		/* Maximum y coordinate */
	double x_inc;		/* x increment */
	double y_inc;		/* y increment */

};

/* globally defined variables */


double xmin = 0.0;        /* input grid bounds are entire globe */
double xmax = 360.0;
double ymin = -90.0;
double ymax = 90.0;
double xincr = 0.1;       /* default grid increment is 0.1 */
double yincr = 0.1;
int  depth_is_neg; 
int nrowsin, ncolsin;  

/*  prototypes for locally defined functions */

void print_usage(char *);
short **get_topo(FILE *);

           
int main (int argc, char **argv)
{ 
   int i, j, k, row, col;
   int i_flag, b_flag;
   int icase, skipit, msg_written;
   int nrowsout, ncolsout;
   int startcol, startcol2;
   int endcol, endcol2;
   int startrow, endrow;
   int row_inc, col_inc;
   int error, latfirst, neglon;
   int lon0to360, xgreenwich;   
   int dummy[4] = {0,0,0,0};
   short int **zin;
   double **zu;
   double zout, lat, lon, tmp;
   double fcor;
   char *s;
   FILE *outfile;
   FILE *topofile;
   FILE *ufile, *lfile, *hfile;
   struct GRID_INFO h;    /* for output grid */

   if (argc < 1 ){
      print_usage(argv[0]);
   }
   
/* initialize these ... */

   h.x_min = xmin;
   h.x_max = xmax;
   h.y_min = ymin;
   h.y_max = ymax;
   i_flag = 0;
   b_flag = 0;
   latfirst = 0;
   startcol = endcol = startcol2 = endcol2 = 0;
   startrow = endrow = 0;
   depth_is_neg = 1; 
   error = 0;  
   msg_written = 0;
   outfile = stdout;
   hfile = ufile = lfile = topofile = (FILE *)NULL;

   
/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 
            case 'B':                    /* get output grid bounds */
                s = &argv[i][2];
                if (*s == '/')
                      ++s;
                error = (sscanf(s,"%lf", &h.x_min) != 1);
                while (*(s++) != '/')
                    ;  
                error += (sscanf(s,"%lf", &h.x_max) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_min) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_max) != 1);

                if (h.x_min > h.x_max)  {
                  fprintf(stderr,"\nW bound cannot exceed E bound.\n");
                  error = 1;
                }
                                               
                if (h.y_min > h.y_max) { 
                  fprintf(stderr,"\nS bound cannot exceed N bound.\n");
                  error = 1;
                }
                b_flag = 1;          
                break;

           case 'H':
               hfile = fopen(&argv[i][2],"r");
               if (hfile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
	       
            case 'I':           /* get output grid increments */
               i_flag = 1;
               ++s;
               if (*s == '/')
                  ++s; 
               error += (sscanf(s,"%lf", &h.x_inc) != 1); 
               s = strchr(s,'/'); /* check for another delimiter*/ 
               h.y_inc = h.x_inc; 
               if (s != NULL) { 
                  ++s; /* move past delimiter */ 
                  error += (sscanf(s,"%lf", &h.y_inc) != 1); 
               } 
               break;
               
           case 'L':
               lfile = fopen(&argv[i][2],"r");
               if (lfile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
           case 'O':
               outfile = fopen(&argv[i][2],"w");
               if (outfile == NULL) {
                  fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                  exit(1);
               }
               break;
               
           case 'T':
               topofile = fopen(&argv[i][2],"r");
               if (topofile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
               
           case 'U':
               ufile = fopen(&argv[i][2],"r");
               if (ufile == NULL) {
                  fprintf(stderr,"\nError opening %s for input\n", &argv[i][2]);
                  exit(1);
               }
               break;
	       
             case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            case ':':
               latfirst = 1;
               break;
	       
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
   
   
 
   if (lfile && !ufile) {
        fprintf(stderr,"\nAn upper surface must be specified with -U when");
        fprintf(stderr,"\na lower surface is specified with the -L option. \n");
        fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
        exit(1); 
   }



/*  Determine type of input files */
 
    if ( !(hfile || ufile) ) {
       icase = 2;
       if (topofile == NULL) {
          topofile = fopen(BATHPATH,"r");
          if (topofile == NULL) {
             fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
             fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
             exit(1);
          }
       }
    }
    else if (hfile) {
       icase = 1;
    }
    else {
    
       icase = 4;
       
       if (!lfile) {
       
	 icase = 3;
         if (topofile == NULL) { 
           topofile = fopen(BATHPATH,"r");
           if (topofile == NULL) {
             fprintf(stderr,"\nUnable to open default topography file\n", BATHPATH);
             fprintf(stderr,"You can specify the name of a binary topography file for input with [-T]\n"); 
             exit(1);
           }
	 }
       }
       
    } /* end else */

/* Check for errors */
    
    if (icase == 2 && !b_flag)  {
          fprintf(stderr,"\nMust specify bounds [-Bw/e/s/n]  \n");
          fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
          exit(1);
    }
   
    if (icase == 4 && !(b_flag && i_flag)) {
          fprintf(stderr,"\nMust specify bounds and grid increment when using -U and -L together.\n");
          fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
          exit(1);
    }
    

/**************************/
/* CASE 1: H is specified */
/**************************/

  if (icase == 1) {
     while ((error = fscanf(hfile, "%lf %lf %lf", &lon, &lat, &zout)) == 3) {
       if (latfirst) {
          tmp = lat;
          lat = lon;
          lon = tmp;
       }
       if (zout < 0. && !msg_written) {
          fprintf(stderr,"\nNegative H values are skipped.");
	  msg_written = 1;
	  zout = 0.0;
       }
       fcor = hb_coriol(lat);
       if (lat < 0 && fcor >= 0) 
         fcor = -fcor;
       if (zout != 0.0)
          fprintf(outfile,"%8.3lf %8.3lf %.4e\n", lon, lat, fcor/zout);
        

    }  /* end while */
 
    if (error != EOF){
      fprintf(stderr,"\nError reading H file \n");
      exit(1);
    }
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);

 } /* end icase == 1 */

/**************************/
/* CASE 2: H is topography */
/**************************/


 if (icase == 2) {
 
    /* xincr, yincr are  grid increments for topography file.
       h.x_inc, h.y_inc are grid increments for output.  Output grid incr
       may be >= than the input incr (output grid will subsample input), but
       not smaller. */
    
     fprintf(stderr,"Increment for topofile: [%.2lf/%.2lf] \n", xincr, yincr);
   

      nrowsin = (int) (NINT((ymax - ymin)/ yincr)) + 1;
      ncolsin = (int) (NINT((xmax - xmin)/ xincr)) + 1; 
   
      /* read in topography file ... */

      zin = get_topo(topofile);
      fclose(topofile);
    
/*   output values according to gridbounds */
    
      startrow = (int) (NINT(( h.y_min - ymin) / yincr) + 0.0001);
      endrow = (int) (NINT(( h.y_max - ymin) / yincr) + 0.0001);
      
      row_inc = 1;
      if (h.y_inc > yincr) 
         row_inc = NINT (h.y_inc / yincr);
      else 
         h.y_inc = yincr;  

     /* find starting/ending columns.  Allow for crossing Greenwich meridian */
   
      lon0to360 = 1;   
      if (h.x_min < 0)
         lon0to360 = 0;
                             
      xgreenwich = 0;
      if (!lon0to360)
         if (h.x_max >= 0) 
             xgreenwich = 1;
      
      col_inc = 1;
      if (h.x_inc > xincr) 
        col_inc = NINT (h.x_inc / xincr);
      else
        h.x_inc = xincr;
   
      if (lon0to360) {   /* the simple case where lon bounds are all positive */
         startcol = (int)( NINT(h.x_min / xincr) + 0.0001);
         endcol = (int)( NINT(h.x_max  / xincr) + 0.0001);
      }
   
      if (!lon0to360) {  /* requested at least one negative longitude */
   
        startcol = (int) (NINT((h.x_min + 360.) / xincr) + 0.0001);
        endcol = (int) (NINT((h.x_max + 360.) / xincr) + 0.0001);
     
        if (xgreenwich) {
           endcol = ncolsin - 1;  /* include greenwich */
           startcol2 = 1;     /* don't repeat greenwich */
           endcol2 = (int) (NINT(h.x_max / xincr) + .0001);
        }
      }
   
    /* adjust output grid bounds to mirror input gridnodes */

      if (h.x_min >= 0) 
          h.x_min = startcol * xincr + xmin;
      else
          h.x_min = startcol * xincr + xmin - 360.0;
      
      if (xgreenwich)
           h.x_max = endcol2 * xincr + xmin;
      else if (h.x_max >= 0)
           h.x_max = endcol * xincr + xmin;
      else     
           h.x_max = endcol * xincr + xmin - 360.0;
       
      h.y_min = startrow * yincr + ymin;
      h.y_max = endrow * yincr + ymin;

      fprintf(stderr,"\nOutput grid bounds: %8.1lf/%8.1lf/%8.1lf/%8.1lf", h.x_min, h.x_max, h.y_min, h.y_max); 
      fprintf(stderr,"\n    grid increment: %6.3lf/%6.3lf\n", h.x_inc, h.y_inc); 
  
      h.node_offset = 0;
      
   /* allocate memory and initialize arrays...  */
   
      nrowsout = NINT((endrow - startrow) / row_inc) + 1;  
      ncolsout = (int) (NINT((h.x_max - h.x_min) / h.x_inc) + 1); 
   
      fprintf(stderr,"\n nrows: %4d   ncols: %4d\n", nrowsout, ncolsout); 
   
    /* subsample the input array and output values ... */
      
    fprintf(stderr,"\nWriting xyz file with seafloor as positive values.\n");
    
    for (row = startrow; row <= endrow; row+=row_inc) {
    
      for (col = startcol; col <= endcol; col+=col_inc) {
      
        zout = (double) zin[row][col];
        
        if (depth_is_neg)
           zout = -zout;
	   
       if (zout < 0. && !msg_written) {
          fprintf(stderr,"\nf/H for Land values are skipped.");
	  msg_written = 1;
	  zout = 0.0;
       }
	lat = ymin + row * yincr; 
	lon = xmin + col * xincr;
	if (!lon0to360 && lon > 180)
	    lon -=  360.;
	    
	fcor = hb_coriol(lat);
        if (lat < 0.0 && fcor > 0.0) 
	   fcor = -fcor;
	if (zout != 0.0 )
	   fprintf(outfile,"%8.3lf %8.3lf %.4e\n", lon, lat, fcor/zout);
      }
      
      if (startcol2 > 0) {
        for (col = startcol2; col <= endcol2; col+=col_inc) {
           zout = (double) zin[row][col];
           if (depth_is_neg)
              zout = -zout;
	   lat = ymin + row * yincr; 
	   lon = xmin + col * xincr;
	   if (!lon0to360 && lon > 180)
	       lon -= - 360.;
	   fcor = hb_coriol(lat);
           if (lat < 0.0 && fcor > 0.0) 
	      fcor = -fcor;
	   if (zout != 0.0 )
	      fprintf(outfile,"%8.3lf %8.3lf %.4e\n", lon, lat, fcor/zout);
        }
      }
    } 

    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);
    
} /* end icase == 2 */ 


   
/**************************************************************/
/* CASE 3:   H is distance between ufile surface and seafloor */
/**************************************************************/
 

 if (icase == 3) {

   /* read in topography file ... */

    fprintf(stderr,"Increment for topofile: [%.2lf/%.2lf] \n", xincr, yincr);
    nrowsin = (int) (NINT((ymax - ymin)/ yincr)) + 1;
    ncolsin = (int) (NINT((xmax - xmin)/ xincr)) + 1; 
   
    zin = get_topo(topofile);
    fclose(topofile);
    
   /* read values from ufile and subtract seafloor depth */
   
    while ((error = fscanf(ufile, "%lf %lf %lf", &lon, &lat, &zout)) == 3) {
      if (latfirst) {
          tmp = lat;
          lat = lon;
          lon = tmp;
       }

       /* find index to topo grid from lat/lon */
    
      row = NINT((lat - ymin ) / yincr);
      if (row >= nrowsin) 
          fprintf(stderr,"\nGridnode for latitude %.3lf is out of topography bounds.", lat);
    
      neglon = 0;
      if (lon < 0) {
        lon += 360.0;
	neglon = 1; 
      }
    
      col = NINT((lon - xmin) / xincr);
    
      if (col >= ncolsin) 
          fprintf(stderr,"\nGridnode for longitude %.3lf is out of topography bounds.", lon);


      fcor = hb_coriol(lat);
      if (lat < 0 && fcor > 0.0) 
           fcor = -fcor;

      tmp = (double) zin[row][col];
      
      if (depth_is_neg)
         tmp = -tmp;
     
      if (neglon) 
          lon = lon - 360. ; 
      if ((zout = (tmp-zout)) > 0.0)
         fprintf(outfile,"%8.3f %8.3f %.4e\n", lon, lat, fcor/zout);
        
    }  /* end while */
 
    if (error != EOF){
      fprintf(stderr,"\nError reading U file \n");
      exit(1);
    }
 
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);
    
  } /* end if icase == 3 */
  
/*********************************************************/
/* CASE 4:   H is distance between 2 specified surfaces  */
/*********************************************************/

  /* initialize arrays */
   
   nrowsin = (int) (NINT((h.y_max - h.y_min)/ h.y_inc)) + 1;
   ncolsin = (int) (NINT((h.x_max - h.x_min)/ h.x_inc)) + 1; 
   zu = (double **) malloc(nrowsin * sizeof(double *));
   
   for (row = 0; row < nrowsin; ++row) {
     zu[row] = (double *) malloc(ncolsin * sizeof(double));
     for (col = 0; col < ncolsin; ++col) 
      zu[row][col] = (double)MISSING;
   }
    
   /* read upper surface into zu arrays */ 
   
    while ((error = fscanf(ufile, "%lf %lf %lf", &lon, &lat, &zout)) == 3) {
      if (latfirst) {
          tmp = lat;
          lat = lon;
          lon = tmp;
       }

       col = NINT((lon - h.x_min) / h.x_inc);
       row = NINT((lat - h.y_min) / h.y_inc);
       
      skipit = 0;
      if ((col < 0) || (col >= ncolsin) || (row < 0) || (row >= nrowsin)) {
          fprintf(stderr,"\nPosition in lower surface out of gridbounds [%.3lf / %.3lf]", lon,lat);
          skipit = 1;
      }
      
      if (!skipit)
         zu[row][col] = zout;
	 
   } /* end while */
   
    
   if (error != EOF){
      fprintf(stderr,"\nError reading U file \n");
      exit(1);
   }
   
   /* read lower surface, subtract upper surface, and write out values */
    
    while ((error = fscanf(lfile, "%lf %lf %lf", &lon, &lat, &zout)) == 3) {
      if (latfirst) {
          tmp = lat;
          lat = lon;
          lon = tmp;
       }

       col = NINT((lon - h.x_min) / h.x_inc);
       row = NINT((lat - h.y_min) / h.y_inc);
       
      skipit = 0;
      if ((col < 0) || (col >= ncolsin) || (row < 0) || (row >= nrowsin)) {
          fprintf(stderr,"\nPosition in lower surface out of gridbounds  [%.3lf / %.3lf]", lon,lat);
          skipit = 1;
      }
      
      if (!skipit) {
        if ((tmp = zu[row][col]) > MISSING) {
           fcor = hb_coriol(lat);
           if (lat < 0 && fcor > 0) 
	       fcor = -fcor;
	   if (tmp < zout) tmp = zout;
	   if ((zout = (zout-tmp)) != 0.0)
	     fprintf(outfile,"%8.3f %8.3f %.4e\n", lon, lat, fcor/zout);
        }
      }
      	 
   } /* end while */
   
   if (error != EOF){
      fprintf(stderr,"\nError reading L file \n");
      exit(1);
   }
 
   fprintf(stderr,"\n End of %s.\n", argv[0]);
   exit(0);
      
} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n%s outputs {lon,lat,f/H} triplets, where H is specified in one of several ways:", program);
   fprintf(stderr,"\n - as an input file that directly specifies lon/lat/H triplets");
   fprintf(stderr,"\n - use bathymetric depth for H");
   fprintf(stderr,"\n - specify 2 xyz-files -- the depths of an upper and lower surface");
   fprintf(stderr,"\n      which will be subtracted to produce H (layer thickness)"); 
   fprintf(stderr,"\n - specify an upper surface (xyz-file) and use  ");
   fprintf(stderr,"\n      bathymetry as the lower surface.   ");
   fprintf(stderr,"\n\nUsage:  %s   [-B<w/e/s/n>] [-H<xyz_file>] [-I<xincr[/yincr]>] [-L<lower_xyz_file>] [-N] [-O<outfile>] [-T<topofile>] [-U<upper_xyz_file>] ", program);
   fprintf(stderr,"\n\n   OPTIONS:"); 
   fprintf(stderr,"\n[-B] : specify output grid bounds. [Default: 0/360/-90/90]");
   fprintf(stderr,"\n        Necessary when using -U and -L options to define H");        
   fprintf(stderr,"\n        OR when using bathymetery alone to define H.");        
   fprintf(stderr,"\n        Meaningless when using -H option or -U alone.");        
   fprintf(stderr,"\n[-H] : file of lon/lat/H triplets "); 
   fprintf(stderr,"\n[-I] : specify output x- and y-grid increments when output is gridded data. "); 
   fprintf(stderr,"\n        ex: -I.1/.5    default increment is .1 "); 
   fprintf(stderr,"\n[-L] : file defining lon/lat/depth of lower surface. "); 
   fprintf(stderr,"\n        Requires -U option.");    
   fprintf(stderr,"\n[-O] : output file. Default is stdout"); 
   fprintf(stderr,"\n[-T] : specify full pathname of binary topography file.  ");
   fprintf(stderr,"\n       Default is [%s]", BATHPATH);
   fprintf(stderr,"\n[-U] : file defining lon/lat/depth of upper surface.");
   fprintf(stderr,"\n        If no -L<lower_surface_file>,  seafloor depth");
   fprintf(stderr,"\n        will be used as lower surface."); 
   fprintf(stderr,"\n[-:] : lat is in column 1, lon column 2, for input files. "); 
   fprintf(stderr,"\n         Default order is lon,lat"); 
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");
   return;
   
} /*end print_usage() */
   
/****************************************************************************/
/****************************************************************************/

short **get_topo(FILE *fptr)
   /* allocates memory (nrowsin and ncolsin are global variables)
      and reads in topography values.  Returns a pointer to
      the start of the memory block. An error causes an error message to be
      printed and an exit.   */
{

   int row, n;
   short **z;
   
   /* Allocate memory */
   
   z = (short **) malloc(nrowsin * sizeof(short *));
   if (z == NULL) {
      fprintf(stderr,"\nError allocating memory.\n");
      exit(1);
   }
   
   for (row = 0; row < nrowsin; ++row) {
     z[row] = (short *) malloc(ncolsin * sizeof(short));
     if (z[row] == NULL) {
         fprintf(stderr,"\nError allocating memory.\n");
         exit(1);
      }
   }
   
   fprintf(stderr,"\nReading in topography values ");
   
   for (row = 0; row < nrowsin; ++row) {
     n = fread((short *)z[row], sizeof(short), ncolsin, fptr);
     if (n != ncolsin) {
         fprintf(stderr,"\nError reading the topofile at row %d\n", row);
         exit(1);
     }
     if ((row % 10) == 0)   /* signal progress */
        fprintf(stderr,".");
   }
   
   fprintf(stderr,"\nFinished reading topofile.\n");
   
   return(z);

} /* end get_topo() */
/****************************************************************************/
