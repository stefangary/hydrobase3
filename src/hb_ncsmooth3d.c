/*  hb_ncsmooth3d.c

................................................................................
                          *******  HydroBase 3 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Mar 2001
			     Updated for HB3 Feb 2012
...................................................
*
*  Smooths each point in HydroBase netcdf files 
*  using a gaussian filter along isopycnal
*  surfaces.  Masked and empty nodes are not changed.  X and Y radius 
*  in gridpoints contrain the size of the smoothing ellipse.  Weights
*  are computed as e^- PI *[dist/L]^2 where L is a user-specified lengthscale
*  default = 100 km.
*  
*  Output is a HydroBase3 nc gridded file.  

*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_grids.h"
#include "hb_memory.h"

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""


#define    PREFILL      1        /* turn off prefilling of cdf variables */


struct MIXLAYER {
    double density;
    double depth;
    int nobs;
};


/* globally defined variables */

int xrad, yrad;
double lengthscale; 
float HB_f_mask;        /* float values are used since the output */
float HB_f_empty;       /* cdf values are of this type */
int add_stats;          /* determines if count and error variables are output*/


/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *, int *);
int cdf_construct(struct GRID_INFO *, int, int *, char *, int, char **);
void get_x_profile(float, float, int, int, char **, char *, char *, double *, short *, float);
void get_sig_profiles(float, float, float *, int, char **, char *, char *, double *, double *, double *, short *, short *, short *, double *, double *, double *, double *, double *);

void define_mixed_layer(  double **,  double *, struct GRID_INFO *, int, int, float, struct MIXLAYER *);
void get_neighbors(double, double *, double,  float *, int , struct GRID_INFO *, double **, double **);
void x_to_stddepth(double *, short *, float *, double *, double *, float, int );

main(int argc, char **argv)
{

  int i, j, k, n, nprops, allprops,found;
  int bflag, iflag, oflag, pflag, mask_it;
  int nfiles, nz, nsqout, nwsq1, nwsq2;
  int row, col, sq, sq_out, error, tbin, datagap;
  int wsq1, wsq2, wsq_mid;
  int wrow, wcol, wcol1, wrow1;
  int dummy, is_te, is_sa, is_pr;
  int smooth_it, smooth_bottom;
  int infile, outfile, print_msg; 
  int startlev, endlev;
  int pixel_grid;           /* pixel or gridnode registration */
  int *prop_indx, *all_prop_indx, *prop_req, *prop_avail;
  float wlat, wlon, lat, lon, wlat2, wlon2;
  double xoffset, yoffset;
  double xradius, yradius;
  double alpha;
  double  isopyc_above;
  double pbot, curr_depth, pabove;
  double zmin, zmax, zval;
  double *zin, *weights;
  char *outfile_name;
  char *dir, *extent;
  char *st;
  struct CDF_HDR *cdfin;
  struct GRID_INFO hout, hwork1, hwork2;
  short *out_count;
  short **tcount1, **scount1, **pcount1;
  short **xcount1;
  float *bottom_de, *bdwork2, *bdwork1;
  float *xout, *xtmp;
  double *xsmooth;
  double *xtmp_s, *ptmp_s;
  double **xwork1, **pwork1;
  double **swork1, **twork1;
  double **sig0work1, **sig1work1;
  double **sig2work1, **sig3work1;
  double **sig4work1;
  double **sigptr, **xwork2;
  double **sig0work2, **sig1work2;
  double **sig2work2, **sig3work2;
  double **sig4work2;
  double **sout, **pout, **dsmooth;
  double **isopyc;
  struct MIXLAYER **mlptr;

/* Set these default values. */

  HB_f_mask = HBMASK;      /* HydroBase cdf file mask flag */
  HB_f_empty = HBEMPTY;    /*  and empty flag */
  smooth_bottom = 1;
  xrad = yrad = 1;

  error = 0;
  bflag = iflag = oflag = pflag =  0;
  pixel_grid = 1;          /* default pixel registration */

  dir = DIR;
  extent = EXTENT;
  nfiles = 0;
  print_msg = 1;
  tbin = 0;
  zmin = zmax = 0.0;
  lengthscale = 100.0;  /* km */
  add_stats = 1;   /* global variable determines whether variance is output */
  alpha = 0.0;
/*----------------------------------------*/  
  
/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
/*----------------------------------------*/  
/* parse command line arguments */

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &hout.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &hout.y_max) != 1);
                        	     
	  if (hout.x_min > hout.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (hout.y_min > hout.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
          break;
 
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
          error = (sscanf(&argv[i][2],"%lf", &hout.x_inc) == 1) ? 0 : 1;
          hout.y_inc = hout.x_inc;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%lf", &hout.y_inc) == 1) ? 0 : 1;
               break;
            }
          }
                        
          break;
	  
        case 'L':
          error = sscanf(&argv[i][2],"%lf", &lengthscale) != 1;
          break;
	  
        case 'O':
	  oflag = 1;
          outfile_name = &argv[i][2];
          break;
	  	      
        case 'P':
          pflag = 1;
	  prop_req = (int *) calloc(MAXPROP, sizeof(int));
          nprops = parse_p_option(&argv[i][2], prop_req);
          break;

        case 'S':
          error = sscanf(&argv[i][2],"%d", &xrad) != 1;
          yrad = xrad;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%d", &yrad) == 1) ? 0 : 1;
               break;
            }
          }
          break;
         break;
	 
        case 'W':
	   fprintf(stderr, "WARNING:  -W<weight_factor> not supported any longer.  Use combination of -L<lengthscale> and -S<xrad>/<yrad> to specify weights\n");
         break;

        case 'Z':
	  smooth_bottom = 0;
	  st = &argv[i][2];
	  if (*st == '\0') {
	     error = 1;
	     break;
	  }
	  if (*st == 'b' || *st == 'B') {
             smooth_bottom = 1;
	     ++st;
	  }
	  
	  if (*st == '\0') { /* this means smooth only the bottom */
	     smooth_bottom = -1;
	     break;
	  }
	  
	  error = (sscanf(st,"%lf/%lf", &zmin, &zmax) != 2);
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
/*  Check syntax of options */ 
   
   error = 0;
   
    if (!nfiles ) {
       fprintf(stderr,"\nYou must specify input cdf_files as first argument(s).");
       ++error;    
    }
    if (!bflag ) {
       fprintf(stderr,"\nYou must specify bounds with -B<w/e/s/n> ");
       ++error;    
    }
    if (!iflag ) {
       fprintf(stderr,"\nYou must specify grid spacing with -I<xincr>[/yincr] ");
      ++error;
    }
    if (!oflag ) {
       fprintf(stderr,"\nYou must specify -O<output_cdf_file> ");
      ++error;
    }
                                 
    if (!pflag ) {
       fprintf(stderr,"\nYou must specify a list of properties for output: -Ppr/te/sa/ox ");
      ++error;
    }
       
   if (zmin > zmax) {
     fprintf(stderr,"\nMinimum depth must exceed maximum depth in -Zrange\n");
      ++error;
   }
    
   if (error) {
    fprintf(stderr,"\nUse -h for complete usage info. \n");
    exit(1);
   }
 

/*-------------------------------------------*/   
/* Get array of standard depths from first cdf file in list */

   i = 1;
   infile = -1; 
   while (infile < 0 && i <= nfiles) {
     infile = cdf_open(dir, argv[i], extent, FALSE);
     ++i;
   }
   
   if (infile < 0) {
      fprintf(stderr,"Unable to open any cdf files.");
      exit(1);
    }
       
   cdfin = (struct CDF_HDR *) calloc(1, sizeof(struct CDF_HDR));
   
   if (error = read_cdf_hdr(infile, cdfin)) 
         exit (1);
   
   xtmp = (float *) calloc((size_t)cdfin->nz, sizeof(float));
   NSTDLEVS = read_cdf_depths(infile, xtmp);
   for (i = 0; i < NSTDLEVS; ++i) 
       std_depth[i] = (double) xtmp[i];
       
   std_depth_initialized = 1;
   free(xtmp);
   cdf_close(infile);
   nz = NSTDLEVS - 1;   /* index to bottom depth */   
   
/*-------------------------------------------*/ 
/* Determine depth limits */

  
     startlev = 0;          /* default: all levels including bottom */
     endlev = nz;
     
     if (!smooth_bottom)
        --endlev;
	
     if ( smooth_bottom < 0)
        startlev = nz;    /* smooth only the bottom */
     
     if (zmax > 0.0) {      /* limits were specified */
        i = 0;
	found = 0;
	while (i < NSTDLEVS && !found) {
	   if (std_depth[i] >= zmin) {
	      startlev = i;
	      found = 1;
	   }
	   else 
	      ++i; 
	}

	if (startlev >= NSTDLEVS) {
	   fprintf(stderr,"FATAL ERROR: Min depth in -Zrange exceeds maximum std depth in cdf_files\n");
	   exit(1);
	}
	
	i = startlev +1;
	found = 0;
	while (i < nz && !found) {
	   if (std_depth[i] > zmax) {
	      endlev = i-1;
	      found = 1;
	   }
	   else
	      ++i; 
	}
     }
   
/*-------------------------------------------*/   
/* Check for requested properties availability.    After this, prop_indx[] and
    nprops will hold the property info for the output cdf file. */

   prop_indx = (int *) calloc((size_t)nprops, sizeof(int)); 
   prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int)); 
   for (i = 0; i < cdfin->nprops; ++i)    
      prop_avail[get_prop_indx(cdfin->prop_id[i])] = 1;
      
   n = 0;   
   for (i = 0; i < nprops; ++i) {
      if ( prop_avail[prop_req[i]] ) {
        prop_indx[n++] = prop_req[i];
      }
      else {
        if (prop_req[i] == (int) PR || prop_req[i] == (int) TE || prop_req[i] == (int) SA)
          fprintf(stderr,"\n FATAL ERROR!! Property %.2s not available in cdf file (pr te sa are mandatory).\n", get_prop_mne(prop_req[i]));       
	else
          fprintf(stderr,"\n WARNING!! Property %.2s not available in cdf file and will not be output.\n", get_prop_mne(prop_req[i]));       
      }
   }
   nprops = n;
   free(prop_avail);    
   free(prop_req);
   
   /*  add variance variables if available */
   allprops = nprops;
   all_prop_indx = prop_indx;
   if (cdfin->counts_included > 1) {
      add_stats = 3;
      allprops = nprops * 2;
      all_prop_indx = (int *) calloc((size_t)allprops, sizeof(int)); 
      for (i = 0; i < nprops; ++i) {
         all_prop_indx[i] = prop_indx[i];
	 all_prop_indx[i+nprops] = prop_indx[i] + 700;
      }
   }
  
/*-------------------------------------------*/   
/*   set up CDF_HDR and output file */


   hout.nx = (int) NINT((hout.x_max - hout.x_min)/hout.x_inc);
   hout.ny = (int) NINT((hout.y_max - hout.y_min)/hout.y_inc);
   
   xoffset = 0.5*hout.x_inc;
   yoffset = 0.5*hout.y_inc;

   hout.node_offset = pixel_grid;   
   if (!pixel_grid) {
      xoffset = 0.0;
      yoffset = 0.0;
      ++hout.nx;
      ++hout.ny;
   }
   nsqout = hout.nx * hout.ny;

   outfile = cdf_construct(&hout, nprops, prop_indx, outfile_name, argc, argv);
   
/*-------------------------------------------*/   
/* Adjust order of properties so that sa is smoothed before te -- since 
   salinity is needed to convert potential temperature back to in situ
   temperature. 
*/
   prop_indx[1] = (int)SA;
   prop_indx[2] = (int)TE;
   all_prop_indx[1] = (int)SA;
   all_prop_indx[2] = (int)T90;
 
/*-------------------------------------------*/   
/*  Adjust xmin/ymin to be lower/left gridnode
    and xmax/ymax to upper/right gridnode depending upon the value of 
    node_offset.   
*/
    
   hout.x_min += xoffset;                         /* west gridnode */
   hout.x_max = hout.x_min + (hout.nx -1) * hout.x_inc;   /* east  */
   hout.y_min = hout.y_min + yoffset;            /* south gridnode */
   hout.y_max = hout.y_min + (hout.ny -1) * hout.y_inc;  /* north  */
   hout.node_offset = 0; /* now normalized to gridnode registration */   

/* xrad, yrad are # of gridnodes in smoothing ellipse*/

   xradius = xrad * hout.x_inc;   /* expansion in grid units */
   yradius = yrad * hout.y_inc;
     

   if (pixel_grid)
     fprintf (stderr, "\nUsing %s registration", "pixel");
   else
     fprintf (stderr, "\nUsing %s registration", "gridnode");
     
   fprintf (stderr, "\nGrid dimensions are nx = %d, ny = %d", hout.nx, hout.ny);
   fprintf (stderr, "\nOutput xmin/ymin gridnode: %8.3lf/%8.3lf ", hout.x_min, hout.y_min);
   
      
   fprintf (stderr, "\nSmoothing ellipse:  xradius (%5.2lf deg) /  yradius (%5.2lf deg) \nLengthscale: %5.2lf km", xradius, yradius, lengthscale);
/*-------------------------------------------*/   

/* Set up big work grid -- output grid + extended bounds --
   to store all profiles necessary */

   hwork1.x_min = hout.x_min - xradius;
   hwork1.x_max = hout.x_max + xradius;
   hwork1.y_min = hout.y_min - yradius;
   hwork1.y_max = hout.y_max + yradius;
   hwork1.x_inc = hout.x_inc;
   hwork1.y_inc = hout.y_inc;
   hwork1.node_offset = 0;     /* normalized to gridnode registration */
   
   hwork1.nx = hout.nx + 2 * xrad;
   hwork1.ny = hout.ny + 2 * yrad;
   nwsq1 =  hwork1.nx * hwork1.ny;  


   xwork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pwork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   twork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   swork1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig0work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig1work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig2work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig3work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   sig4work1 = (double **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(double *));
   pcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   tcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   scount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));
   xcount1 = (short **) get_memory((void *)NULL,(size_t)nwsq1, sizeof(short *));

   for (sq = 0; sq < nwsq1; ++sq) {
     pwork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     twork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     swork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig0work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig1work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig2work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig3work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     sig4work1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
     pcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     tcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
     scount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
   }

/*-------------------------------------------*/   
/* Set up small work grid to store isopycnal surface sent to laplacian filter
   function.   */

   hwork2.nx = xrad * 2 + 1;
   hwork2.ny = yrad * 2 + 1;
   hwork2.x_inc = hout.x_inc;
   hwork2.y_inc = hout.y_inc;
   hwork2.node_offset = 0;  /* normalized to gridnode registration */
   nwsq2 = hwork2.nx * hwork2.ny;
   wsq_mid  = yrad * hwork2.nx + xrad;  /* by definition the middle square of work2 */
   
   xwork2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig0work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig1work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig2work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig3work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   sig4work2 = (double **) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double *));
   weights = (double *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double));

 /*-------------------------------------------*/   
 
   bottom_de = (float *) get_memory ((void *)NULL, (size_t) nsqout, sizeof(float));
   bdwork1 = (float *) get_memory((void *)NULL,(size_t)nwsq1, sizeof(float));
   bdwork2 = (float *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(float));
   weights = (double *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(double));
   
   dsmooth = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   pout = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   sout = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   
   isopyc = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof (double *));
   mlptr = (struct MIXLAYER  **) get_memory((void *)NULL, (size_t)nsqout, sizeof (struct MIXLAYER *));
   
   for (sq = 0;  sq < nsqout; ++sq) {
      isopyc[sq] =  (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
      mlptr[sq] = (struct MIXLAYER *) get_memory((void *)NULL, 1, sizeof(struct MIXLAYER));
   }
 
/*-------------------------------------------*/ 
/* Visit each gridpoint in big workspace, read in profiles, compute sigmas. */

   fprintf(stderr,"\nReading in profiles...");
   
   for (row = 0; row < hwork1.ny; ++row ) {
      lat = hwork1.y_min + row * hwork1.y_inc;
      
      for (col = 0; col < hwork1.nx; ++col) {
         lon = hwork1.x_min + col * hwork1.x_inc;
  
	 sq = row * hwork1.nx + col;

	 get_sig_profiles(lat, lon, &bdwork1[sq], nfiles, argv, dir, extent, pwork1[sq], twork1[sq], swork1[sq], pcount1[sq], tcount1[sq], scount1[sq], sig0work1[sq], sig1work1[sq], sig2work1[sq], sig3work1[sq], sig4work1[sq]);

	  /* Convert in situ temperatures to potential temperatures for smoothing */

         for (j = 0; j < NSTDLEVS; ++j) {
           if (twork1[sq][j] > -9.0 && twork1[sq][j] < 100.)
	      twork1[sq][j] = hb_theta(swork1[sq][j], twork1[sq][j], pwork1[sq][j], 0.0);
	 }
	 
      } /* end for col */
   } /* end for row */
   
	 
   fprintf(stderr,"\nNow smoothing...");


/* Loop for each property */

   for (i = 0; i < allprops; ++i) {

       is_pr = (all_prop_indx[i] == (int)PR);
       is_sa = (all_prop_indx[i] == (int)SA);
       is_te = (all_prop_indx[i] == (int)T90 || all_prop_indx[i] == (int)TE);
       
       fprintf(stderr,"\n%s ", get_prop_mne(all_prop_indx[i]));
      
       switch ((enum property) all_prop_indx[i]) {
       
          case PR :              /* just set the pointers */
	     xwork1 = pwork1;
	     xcount1 = pcount1;
	     break;
	  case T90 :
	  case TE :
	     xwork1 = twork1;
	     xcount1 = tcount1;
	     break;
	  case SA :
	     xwork1 = swork1;
	     xcount1 = scount1;
	     break;
	  
	  default:
	  
	   /* Read property from input files */
	   
	   for ( sq = 0; sq < nwsq1; ++sq) {
	      wrow = NINT(sq / hwork1.nx);
	      wcol = sq - wrow * hwork1.nx;
	      wlat = hwork1.y_min + wrow * hwork1.y_inc;
	      wlon = hwork1.x_min + wcol * hwork1.x_inc;
              xwork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
              xcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
	      get_x_profile(wlat, wlon, all_prop_indx[i], nfiles, argv, dir, extent, xwork1[sq], xcount1[sq], bdwork1[sq]);
	   }
       }  /* end switch */
	           
        /* Visit each gridpoint in output file .  
           Determine isopycnal of each level
           and obtain surrounding values of property at that isopycnal.
           Apply the gaussian filter and produce output profile 
           that includes mask and empty flags where appropriate.   Use lat/lon 
           (not row/col) to translate between grids for output/input, and 
           work1 and work2 because each has its own dimension and order. */

        for (row = 0; row < hout.ny; ++row ) {
           lat = hout.y_max - row * hout.y_inc;
           fprintf(stderr,".");
      
           for (col = 0; col < hout.nx; ++col) {
              lon = hout.x_min + col * hout.x_inc;
	 
	      sq_out = row * hout.nx + col;
	 
	  /* reset these at each gridnode*/
	  
	      hwork2.x_min = lon - xrad * hwork2.x_inc;  
	      hwork2.x_max = lon + xrad * hwork2.x_inc; 
	      hwork2.y_min = lat - yrad * hwork2.y_inc;
	      hwork2.y_max = lat + yrad * hwork2.y_inc; 
	 
	 /* set pointers to define small work grids */
	 
	      for (wrow = 0; wrow < hwork2.ny; ++wrow) {
	         for (wcol = 0; wcol < hwork2.nx; ++wcol) {
	    
	           wsq2 = wrow * hwork2.nx + wcol; /* index to small work grids */
                   wlat2 = hwork2.y_min + wrow * hwork2.y_inc;
	           wlon2 = hwork2.x_min + wcol * hwork2.x_inc;
		   
	           wcol1 = NINT((wlon2 - hwork1.x_min) / hwork1.x_inc);
	           wrow1 = NINT((wlat2 - hwork1.y_min) / hwork1.y_inc);
	           wsq1 = wrow1 * hwork1.nx + wcol1;   /* index to big work grids */
	      
		   xwork2[wsq2] = xwork1[wsq1];
	           sig0work2[wsq2] = sig0work1[wsq1];
	           sig1work2[wsq2] = sig1work1[wsq1];
	           sig2work2[wsq2] = sig2work1[wsq1];
	           sig3work2[wsq2] = sig3work1[wsq1];
	           sig4work2[wsq2] = sig4work1[wsq1];
	    
	           bdwork2[wsq2] = bdwork1[wsq1];  /* store value, not ptr */
		   
		   if (wsq2 == wsq_mid)
		     out_count = xcount1[wsq1];
	      
	         } /* end for wcol */
	      }  /*end for wrow */
	      
	      bottom_de[sq_out] = bdwork2[wsq_mid];
	      
	      /* The first time each output square is visited,
	      construct array of isopycnals corresponding to std depth levels
	      and define a mixed-layer density and depth */
	      
	      if (is_pr) {

	          get_weights_c(weights, &hwork2, alpha,  lengthscale, lengthscale, xrad, yrad );

		  define_mixed_layer(sig0work2, weights, &hwork2, xrad, yrad,  bottom_de[sq_out], mlptr[sq_out]);

	          isopyc_above = 0.0;
		  n = 0;
		  if (mlptr[sq_out]->depth > -1.0 && mlptr[sq_out]->depth < 100000) {
		     while (std_depth[n] <= mlptr[sq_out]->depth) {
		        isopyc[sq_out][n] = mlptr[sq_out]->density;
		        ++n;
		     }
		     isopyc_above = mlptr[sq_out]->density + 0.02;  /* by definition, 
		                                      density at bottom of mixed layer */
		  }
	          for (j = n; j < NSTDLEVS; ++j) {
		  
	              curr_depth = std_depth[j];
		      if (j == nz) {
		         curr_depth = bottom_de[sq_out];
		       } 
		       
		      if (curr_depth <= 500 ) {
			      sigptr = sig0work2;
		       }
		       else if (curr_depth > 500 && curr_depth <= 1500) {
		              sigptr = sig1work2;
		       }
		       else if (curr_depth > 1500 && curr_depth <= 2500) {
		              sigptr = sig2work2;
			}
		        else if (curr_depth > 2500 && curr_depth <= 3500) {
		              sigptr = sig3work2;
			}
		        else {
		              sigptr = sig4work2;
			}
			 
		        isopyc[sq_out][j] =  sigptr[wsq_mid][j];
			
			if ( isopyc[sq_out][j] > 0 && isopyc[sq_out][j] < 100) {
			   /* check for density inversion */
			    if (isopyc[sq_out][j] < isopyc_above - 0.001)    
		                 isopyc[sq_out][j] = isopyc_above;
		            
		            isopyc_above = isopyc[sq_out][j];
			 }
	          } /* end for j */
		  
	     } /*end if is_pr */

	      mask_it = is_flagged(bottom_de[sq_out], HB_f_mask);
	      
	      pbot = HB_f_mask;
	      if ( ! mask_it)
	         pbot = hb_p80((double)bottom_de[sq_out], lat);

             xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));

             smooth_it = 0;
	     for (j = 0; j < NSTDLEVS; ++j) {
	     
	        if (mask_it) {
		   xout[j] = HB_f_mask;
		   out_count[j] = 0; 
		}
		else {
	           xout[j] =  (float) xwork2[wsq_mid][j];
		   if ( ! (is_flagged(xout[j], HB_f_empty) || is_flagged(xout[j], HB_f_mask) ) ) {
		      if (j >= startlev && j <= endlev)
		          ++smooth_it;
		      if (out_count[j] == 0)
		          out_count[j] = -1;  /* only when count variables not included in cdf file*/
		   }
		   else {  /* this level is flagged */
		      out_count[j] = 0;
		   }
		}
	     } /* end for j */
	     
	     if (smooth_it) {
	      if ( ! is_pr )
	           get_weights_c(weights, &hwork2, alpha,  lengthscale, lengthscale, xrad, yrad );
	       
               xsmooth = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	       for (j = 0; j < NSTDLEVS; ++j) {
	          xsmooth[j] = (double)xout[j];
               }
	       free((void *) xout);
	   
               /* smooth levels in mixed layer */
               j = 0;
	       while ( j < nz && std_depth[j] < mlptr[sq_out]->depth) {
	      
 	      	  if (! (is_flagged((float)xsmooth[j], HB_f_empty) ||  is_flagged((float)xsmooth[j], HB_f_mask)))  { 
		  
                    if (is_pr) { 
		       if (NINT(std_depth[j]) == 0)
		           xsmooth[j] = 0.0;
		       else 
	                   xsmooth[j] = hb_p80(std_depth[j], lat);  /* specifically set mixed layer pressures */
		       
	            } 
	            else { 
		 
	             if (j >= startlev && j <= endlev) {
		         curr_depth = std_depth[j];
                         zin = (double *) get_memory((void *)NULL, (size_t)nwsq2, sizeof(double));
                         get_neighbors(isopyc[sq_out][j], zin, curr_depth, bdwork2, nwsq2, &hwork2, xwork2, sig0work2);
			 zero_weights_2d(zin, weights,  HB_f_empty,HB_f_mask, &hwork2, xrad, yrad);
	                 xsmooth[j] = weighted_mean(zin, weights, HB_f_empty,HB_f_mask, nwsq2);
                         free((void *)zin);
		     }
		   }
	         }
		 ++j;
			
	       } /* end while stddepth < mixedlayer depth */
	      
	       k = j;
	           
               /* Loop for each depth level below mixed layer including bottom ...
                  determine corresponding isopycnal value, then interpolate to find
	          that isopycnal at all the neighboring gridnodes. Fill
	          in array zin[] to pass to the filter function.
	           		    
	          Smooth bottom properties but NOT bottom pressure.  */
	     
	       for (j = k; j < NSTDLEVS; ++j) {
	      
		  if ( ! ( is_flagged((float)xsmooth[j], HB_f_empty) || is_flagged((float)xsmooth[j], HB_f_mask) ) ) { 
		  
		    get_weights_c(weights, &hwork2, alpha,  lengthscale, lengthscale, xrad, yrad );

		 
	            if (j >= startlev && j <= endlev) {

			     curr_depth = std_depth[j];
		             if (j == nz) {
		                 curr_depth = bottom_de[sq_out];
		              } 
		             if (curr_depth <= 500)
			       sigptr = sig0work2;
		             else if (curr_depth > 500 && curr_depth <= 1500)
		               sigptr = sig1work2;
		             else if (curr_depth > 1500 && curr_depth <= 2500)
		               sigptr = sig2work2;
		             else if (curr_depth > 2500 && curr_depth  <= 3500)
		               sigptr = sig3work2;
		             else 
		               sigptr = sig4work2;
			       
                      
                      zin = (double *) get_memory((void *)NULL, (size_t)nwsq2, sizeof(double));
                      get_neighbors(isopyc[sq_out][j], zin, curr_depth, bdwork2, nwsq2, &hwork2, xwork2, sigptr);
		      zero_weights_2d(zin, weights,  HB_f_empty,HB_f_mask, &hwork2, xrad, yrad);
	              xsmooth[j] = weighted_mean(zin, weights, HB_f_empty,HB_f_mask, nwsq2);
		      		      
		      if (is_pr) {      /* don't allow smoothed pressures to be deeper than bottom pressure */
		         if ( xsmooth[j] > pbot )
			    xsmooth[j] = zin[wsq_mid];  /* revert to unsmoothed pr */
		      }
		      
		      free ((void *)zin);
		    }  /* end if j */ 
		    
		    else if (j == nz && smooth_bottom && !is_pr) {
	              sigptr = sig0work2;
		      curr_depth = bottom_de[sq_out];
		      if (curr_depth > 500 && curr_depth <= 1500)
		           sigptr = sig1work2;
		      else if (curr_depth > 1500 && curr_depth <= 2500)
		           sigptr = sig2work2;
		      else if (curr_depth > 2500 && curr_depth <= 3500)
		           sigptr = sig3work2;
		      else  if (curr_depth > 3500 )
		           sigptr = sig4work2;
		   
                      zin = (double *) get_memory((void *)NULL, (size_t)nwsq2, sizeof(double));
                      get_neighbors(isopyc[sq_out][j], zin, curr_depth, bdwork2, nwsq2, &hwork2, xwork2, sigptr);
		    
		      zero_weights_2d(zin, weights,  HB_f_empty,HB_f_mask, &hwork2, xrad, yrad);
	              xsmooth[j] = weighted_mean(zin, weights, HB_f_empty,HB_f_mask, nwsq2);
		      free ((void *)zin);
		    }
		    
		 } /* end if !is_flagged */
		 
	      } /* end for j */


              /* Interpolate smoothed profile back onto std depths */
	      	      
              xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
	      
	      if (is_pr) {
	        
	        /* Convert smoothed pressure to depth values (dsmooth)  
	          for interpolating other properties onto depths.
		  Convert  std_depths to pressure (pout)*/
	      
	         dsmooth[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	         pout[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));		 
		 
		 for (j = 0; j < nz; ++j) {
		     if (NINT(std_depth[j]) == 0) 
		       pout[sq_out][j] = 0.0;
		     else
		        pout[sq_out][j] = hb_p80(std_depth[j], lat);
		     
		     xout[j] = (float) xsmooth[j];
		     dsmooth[sq_out][j] = xsmooth[j];  /* in case its flagged */
		     if (xsmooth[j] > -1.0 && xsmooth[j] < 99999) {
		           if (NINT(xsmooth[j] == 0))  { /* trap roundoff errors */
			       xsmooth[j] = 0.0;   
		               xout[j] = 0.0;
			   }
			   dsmooth[sq_out][j] = hb_depth(xsmooth[j], lat);
	             }
		     
		     /* replace smoothed pr values with pressure of stddepth */
		     if (xout[j] >= 0.0 && xout[j] < 99999.)  {  
		        xout[j] = (float) ABS(pout[sq_out][j]); 
		     }
		     
		 } /* end for j */
		 
		 /* explicitly set bottom pressure and depth */
		 
		 pout[sq_out][nz] = pbot;
		 dsmooth[sq_out][nz] = bottom_de[sq_out] + 5.0; /* a bit deeper to help interpolation */
		 xout[nz] =  (float) pbot;
		 if ( ! (xsmooth[nz] >= 0.0  && xsmooth[nz] < 99999.0) )  {
		    xout[nz] = (float) xsmooth[nz];   /* either a mask or empty flag */
		    out_count[nz] = 0;
		 }
              } /* end if is_pr */
	      else {  /* all other props except pr */
	          if (dsmooth[sq_out] != NULL)
	             x_to_stddepth(xsmooth, out_count, xout, dsmooth[sq_out], std_depth, bottom_de[sq_out], NSTDLEVS);
	      }  /* end else */
	      
	      free((void *) xsmooth);
              
	      /*  special things to do for t and s */
	      
	      if (is_sa) {  /* save output sa values for later computing theta  */
	      
	         sout[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	         for (j = 0; j < NSTDLEVS; ++j) {
		     sout[sq_out][j] = (double) xout[j];
		 }
              }
	   
	   
	      if (is_te) {  /* adiabatically adjust potential temperature to pressure */
	         xtmp = xout;
	         xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
		 
	         for (j = 0; j < NSTDLEVS; ++j) {
	             xout[j] = xtmp[j];
	             if (xtmp[j] > -3.0 && xtmp[j] < 100.0  ) {
	                xout[j] = (float) hb_theta(sout[sq_out][j], (double)xtmp[j], 0.0, pout[sq_out][j]);
	                if (sout[sq_out][j] <= 0 || sout[sq_out][j] > 100. ) {
			    fprintf(stderr,"WARNING:  salinity (%.4e) has no value but temperature does (%.4e) \n",sout[sq_out][j], xtmp[j]);
			    xout[j] = HB_f_empty;
			    out_count[j] = 0;
			}
		    }
	         }
	         free((void *) sout[sq_out]);
	         free((void *) xtmp);
	      }  /* end if is_te */
	   
	   
           } /* end if smooth_it */
	   
	   
 	 if (all_prop_indx[i] < 700) {
	    write_prop_cdf(outfile, xout, get_prop_mne(all_prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	    write_prop_count_cdf(outfile, out_count, get_prop_mne(all_prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	 }
	 else {
	    write_prop_err_cdf(outfile, xout, get_prop_mne(all_prop_indx[i]-700), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	 }
	   
	 free((void *) xout);
	 
       }  /* end for col */
     } /* end for row */
     
     if (!(is_pr || is_te || is_sa)) {
        for (sq = 0; sq < nwsq1; ++sq) {
	   free((void *) xwork1[sq]);
	   free((void *) xcount1[sq]);
	}
       
     }
     
   } /* end for i*/
   
   write_bottom_depth_cdf(outfile, 0, 0, 0, hout.ny, hout.nx, 1, bottom_de);	 
   cdf_close(outfile); 
   
   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);
	
}  /* end main */
   
   
   
/************************************************************************/
void print_usage(char *program)
{

  fprintf(stderr,"\nUses a 2-D Gaussian filter with weights = e^-{ PI *[dist/L]^2} to smooth each");
  fprintf(stderr,"\npoint in HydroBase netcdf files along isopycnals.  ");
  fprintf(stderr,"\nMultiple cdf files can be input, but all");
  fprintf(stderr,"\nmust possess the same depth dimensions, and x- y-");
  fprintf(stderr,"\ngrid increments.  The output cdf file will have the");
  fprintf(stderr,"\nsame depth dimension, but the lat/lon dimensions specified");
  fprintf(stderr,"\nwith -B<w/e/s/n> and -I<x/yincr>  ");
  fprintf(stderr,"\nThe input domain is expanded relative to the output");
  fprintf(stderr,"\nso that data beyond the borders of the grid can be");
  fprintf(stderr,"\nincorporated into the smoothing; x and y radii of smoothing");
  fprintf(stderr,"\nellipse are specified with [-S])");
  fprintf(stderr,"\nSearches in each of 8 directions to find suitable values");
  fprintf(stderr,"\nbut if isopycnal runs into seafloor or outcrops, the search is");
  fprintf(stderr,"\nsuspended in that direction.");   
  fprintf(stderr,"\nThe properties  pr, te, sa  are mandatory ");
  fprintf(stderr,"\nin both input and output files -- regardless");
  fprintf(stderr,"\nof whether they are specified in the -P argument.");

  fprintf(stderr,"\n\nUSAGE:  %s cdf_file(s) -B<w/e/s/n> -I<x_inc>[/<y_inc>] -O<outfile> -P<properties> [-D<input_dir>] [-E<input_file_extent>] [-G] [-L<lengthscale>] [-S<xrad>/<yrad>] [-Z[b][<zmin/zmax>]] [-h] [-:]\n\n", program);
  fprintf(stderr," -B  sets the output grid bounds w/e/s/n.\n");
  fprintf(stderr," -I  sets the output grid spacing for x/y dimensions.\n");
  fprintf(stderr," -O  name of output cdf file.\n");
  fprintf(stderr," -P  list of properties to include in output file\n");
  fprintf(stderr,"        ex:  -Ppr/t90/sa/ox\n");
  fprintf(stderr,"       -P (by itself) produces a list of available properties.\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr,"-D directory for input cdf files.\n");
  fprintf(stderr,"-E file extent for input cdf files.\n");
  fprintf(stderr,"-G force gridnode registration for output(default is pixel registration).\n");
  fprintf(stderr,"-L lengthscale (L) for computing weights = e^-{pi * [dist/L]^2}  \n");
  fprintf(stderr,"      Default is [-L%.3lf]\n", lengthscale);
  fprintf(stderr,"-S x- and y-radii (integer gridnodes) of smoothing ellipse \n");
  fprintf(stderr,"      Default is [-S%1d/%1d]\n", xrad, yrad);
  fprintf(stderr,"-Z  zmin/zmax values limits smoothing to levels\n");
  fprintf(stderr,"   between those depths. \n"); 
  fprintf(stderr,"   'b' switches on option to smooth the bottom layer. \n");
  fprintf(stderr,"   ex: -Zb0/1500 smoothes depth levels 0->1500 meters plus\n");
  fprintf(stderr,"            the bottom layer. \n");
  fprintf(stderr,"       -Z1500/3000 smooths just 1500->3000 m leaving  \n");
  fprintf(stderr,"           others --including bottom level--unaltered. \n");
  fprintf(stderr,"       -Zb smooths just the bottom.  \n");
  fprintf(stderr,"   Default (no -Z argument): smooth all levels including bottom layer.\n");
  fprintf(stderr,"-h help...... prints this message. \n");
  return;
  
} /* end print_usage() */
/*****************************************************************************/
/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
  char prop[6];
  int n;

  if (*st == '\0') {
         print_prop_menu();
         exit(0);
  }
  
  
  /* Explicitly set the pr, t90, sa because they are mandatory */
  
  prop_indx[0] = (int)PR;
  prop_indx[1] = (int)T90;
  prop_indx[2] = (int)SA;
  
  n = 3;
  do {
     if (*st == '/')
         ++st;
	 
     prop[0]=prop[1]=prop[2]=prop[3]= prop[4] = prop[5]='\0';
     sscanf(st,"%[^'/']", prop);
     prop_indx[n] = get_prop_indx(prop);
     if (prop_indx[n] < 0)  {
       fprintf(stderr,"\n Unknown property '%s' specified in -P%s\n", prop, st);
       exit(1);
     }
     if (prop_indx[n] > 100)  {
       fprintf(stderr,"\n -P%s : No need to specify variance and count properties, they are automatically output.  Quality flags are not output.\n", prop);
       exit(1);
     }
     st += strlen(prop);

     /* !**!  Special cases for properties ...
     
       No need for pref specs for props like s_ , ht and pe.  No computation of 
       additional props will be done. Any property 
       requested must already exist in the cdf file.  */
     
         /* de, pr, t90, sa are automatically done so don't count them here */

     if ( !( (prop_indx[n] == (int)DE) || (prop_indx[n] == (int)PR) || (prop_indx[n] == (int)T90)|| (prop_indx[n] == (int)TE) || (prop_indx[n] == (int)SA)))  
            ++n;
     

   } while (*st == '/');
   return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
/*****************************************************************************/
int cdf_construct(struct GRID_INFO *hptr, int nprops, int *prop_indx,  char *filename, int nargs, char **arglist)

   /* Opens a cdf output file and writes an appropriate header,
      standard depths, and time bin info.  Returns the id associated
      with the open file */
{
   struct CDF_HDR cdfhdr;
   int i, error, cdfid;

   
   cdfhdr.xmax = (float) hptr->x_max;
   cdfhdr.xmin = (float) hptr->x_min;
   cdfhdr.ymax = (float) hptr->y_max;
   cdfhdr.ymin = (float) hptr->y_min;
   cdfhdr.xincr = (float) hptr->x_inc;
   cdfhdr.yincr = (float) hptr->y_inc;
   
   cdfhdr.nx = (int) NINT((hptr->x_max - hptr->x_min)/hptr->x_inc);
   cdfhdr.ny = (int) NINT((hptr->y_max - hptr->y_min)/hptr->y_inc);
   cdfhdr.node_offset = hptr->node_offset;
   if (! hptr->node_offset) {  
      ++cdfhdr.nx;     
      ++cdfhdr.ny;
   }
   cdfhdr.nz = NSTDLEVS;
   cdfhdr.nt = 1;
   cdfhdr.tmin = (int *) malloc(sizeof(int));
   cdfhdr.tmax = (int *) malloc(sizeof(int));
   cdfhdr.tmin[0] = 0;
   cdfhdr.tmax[0] = 9999;
   cdfhdr.counts_included = add_stats;
   cdfhdr.fill_value =  HB_f_empty;
   cdfhdr.mask_value =  HB_f_mask;
   strncpy(cdfhdr.x_units, "degrees", 8);
   strncpy(cdfhdr.y_units, "degrees", 8);
   strncpy(cdfhdr.z_units, "meters", 7);
   strncpy(cdfhdr.title,"HydroBase", 10);
   strcpy(cdfhdr.command, *arglist);
   for (i = 1; i < nargs; ++i) {
      strncat(cdfhdr.command, " ", 1);
      strcat(cdfhdr.command, arglist[i]);
   }
   cdfhdr.nprops = nprops;
   cdfhdr.prop_id = (char **) malloc(nprops * sizeof(char *));
   cdfhdr.prop_units = (char **) malloc(nprops * sizeof(char *));
   for (i = 0; i < nprops; ++i) {
      cdfhdr.prop_id[i] = (char *) calloc(5, sizeof(char));
      cdfhdr.prop_units[i] = (char *) calloc(80, sizeof(char));
      strcpy(cdfhdr.prop_id[i], get_prop_mne(prop_indx[i]));
      strcpy(cdfhdr.prop_units[i], get_prop_units(prop_indx[i]));
   }
   
/* Open output file and write out some info so we can free up some memory... */

   cdfid = cdf_init(filename);   
   error = cdf_define(cdfid, &cdfhdr, PREFILL, cdfhdr.counts_included);
   if (error)  exit(1);
   error = write_std_depths_cdf(cdfid, &cdfhdr);
   error = write_time_bins_cdf(cdfid, &cdfhdr);
   error = write_lat_vector(cdfid, &cdfhdr);
   error = write_lon_vector(cdfid, &cdfhdr);
   
   for (i = 0; i < nprops; ++i) {
      free((void *)cdfhdr.prop_id[i]);
      free((void *)cdfhdr.prop_units[i]);
   }
   free((void *) cdfhdr.prop_id);
   free((void *) cdfhdr.prop_units);
  
   
   
   return(cdfid); 
} /*end cdf_construct() */ 

	   
	   
/************************************************************************/

void get_x_profile(float lat, float lon, int pindex, int nfiles, char **argv, char *dir, char *ext, double *x, short *count, float bdepth )

 /*  Reads each cdf file in list and extracts depth profile of property
     indexed by pindex.  Returned array contains either a value, HB_f_empty
     of HB_f_mask.  If a count variable is present the info is 
     stored in the count array, otherwise the count array is assigned
     either 0 or 1.    */

{
   struct CDF_HDR cdf;
   int cdfid, curfile, print_msg = 0;
   int lon0to360;
   float *xin;
   short *nobs;
   char *mne;
   int error, i,iprop,   row, col;
   int tbin = 0, lev;
   
   

   lon0to360 = 1;
   if (lon < 0)
      lon0to360 = 0;

/* Loop for each input file */

   curfile = 1;
   
   do {
      cdfid = cdf_open(dir, argv[curfile], ext, print_msg);
      if (cdfid < 0)
         goto NEXTFILE;
	 
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
/* compare bounds of file to the lat/lon specified.  Try
   to match the convention for longitude -- usually W is neg, E is pos --
   but we cannot be certain it is always specified like this .*/
	 
      if (lon0to360) {    /* target lon is positive */
        if (cdf.xmax < 0) {  /* all input longitudes are neg, recast all as positive */
	   cdf.xmax += 360.0;
           if (cdf.xmin < 0)
	      cdf.xmin += 360.0;
        }
      }
      else {              /* target lon is negative */
           if (cdf.xmin > 0) {  /* all input longitudes are positive */
	      cdf.xmin -= 360.0;
              if (cdf.xmax > 0)
	         cdf.xmax -= 360.0;
           }
      } 

      if ((cdf.xmin > lon) || (cdf.xmax < lon) 
      || (cdf.ymin > lat) || (cdf.ymax < lat)) {
         cdf_close(cdfid);
	 goto NEXTFILE;
      }


      error = get_indices(&cdf, lat, lon, &row, &col);
      
      if (row < 0 || row >= cdf.ny || col < 0 || col >= cdf.nx) {
         cdf_close(cdfid);
         goto NEXTFILE;
      }
      
/* check that cdf file has same number of standard depths  */
     
      if (cdf.nz != NSTDLEVS) {
         fprintf(stderr, "\nFATAL ERROR:  Mismatch of standard depths in cdf files.\n");  
         exit(1);
      }


/* determine if property is available in this file */

     iprop = pindex;
      if (pindex >= 700 && cdf.counts_included > 1)
          iprop = pindex - 700;
      
      i = 0;
      while (i < cdf.nprops && !(iprop == get_prop_indx(cdf.prop_id[i])))
         ++i;
	 
      	 
      if (i >= cdf.nprops) {
      
            fprintf(stderr, "\nWARNING: Property [%s] not available in this file.", get_prop_mne(pindex));
            cdf_close(cdfid);
	    goto NEXTFILE;
      } 
	 
      
/* allocate space for property arrays; */

      xin = (float *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(float));
      nobs = (short *) get_memory((void *)NULL, (size_t)cdf.nz, sizeof(nobs));
      

 	   
    if (pindex >= 700 && cdf.counts_included > 1)   {
        mne = get_prop_mne(pindex-700);
        error = read_cdf_prop_err(cdfid, mne, xin, row, col, tbin, 0, cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read %s_%s at row,col =  %d,%d from cdf file.", mne,ERR_VAR_SUFFIX, row, col);
            exit(1);
        }
     }
     else {
         mne = get_prop_mne(pindex);
        error = read_cdf_prop(cdfid, mne, xin, row, col, tbin, 0, cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read %s at row,col =  %d,%d from cdf file.", mne, row, col);
            exit(1);
        }
     }
     
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,mne,nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read %.2s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
            exit(1);
        }

     }
	
     /* check for flagged levels.  in case of multiples entries,
        sum up the data at each level.  */ 
	   
     for (lev = 0; lev < cdf.nz; ++lev) {
     
       if (!is_flagged(xin[lev], cdf.fill_value) && !is_flagged(xin[lev], HB_f_mask)) {
       
	   if (!cdf.counts_included || nobs[lev] == 0) 
	      nobs[lev] = 1;
	   
	   x[lev] += (double) (nobs[lev] * xin[lev]);
	   count[lev] += nobs[lev];
	}	   
     } /* end for lev */

     free((void *) xin);
     free((void *) nobs);
     cdf_close(cdfid);
     
NEXTFILE:
      ;
      	 
   }  while (curfile++ < nfiles); 
   
   
   /* flag levels with zero count */

   for (lev = 0; lev < cdf.nz; ++lev) {
   
      if (count[lev]) {
         x[lev] /= count[lev];
      }
      else if ( is_flagged(bdepth, HB_f_mask) ) {
          x[lev] =  (double) HB_f_mask; 
      }
      else {
         x[lev] = HB_f_empty;
	 
         if (std_depth[lev] >= bdepth)
              x[lev] = (double) HB_f_mask; 
      }
   
   } /* end for lev */
   
      
   return;

}  /* end get_x_profile() */

/***************************************************/
void get_sig_profiles(float lat, float lon, float *bdepth_ptr, int nfiles, char **argv, char *dir, char * ext, double *ppro, double *tpro, double *spro, short *pcount, short *tcount, short *scount, double *s0pro, double *s1pro, double *s2pro, double *s3pro, double *s4pro)

/*  Search cdf files for lon/lat, read in pr, te, sa, and compute all sigmas.  
    Read in bottom depth into the pointer. Mark masked levels with HB_f_mask.  
    Missing values are flagged with HB_f_empty.

*/
{
   struct CDF_HDR cdf;
   int cdfid, curfile, print_msg = 0;
   int lon0to360;
   float deepest, bd, flag, maskflag;
   float *xin;
   short *nobs;
   int error, i, t90_avail, row, col;
   int tbin = 0, lev;
   int maxlev, bindx;

  /* initialize these */
  
   lon0to360 = 1;
   if (lon < 0)
      lon0to360 = 0;

   bindx = NSTDLEVS - 1;
   deepest = HB_f_empty;
   bd = HB_f_empty;
   flag = HB_f_empty + 10;
   maskflag = HB_f_mask - 10;
   
   xin = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
   nobs = (short *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(short));

/* Loop for each input file */

   curfile = 1;
   
   do {
      cdfid = cdf_open(dir, argv[curfile], ext, print_msg);
      if (cdfid < 0)
         goto NEXTFILE2;
	 
      if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);
	 
/* compare bounds of file to the lat/lon specified.  Try
   to match the convention for longitude -- usually W is neg, E is pos --
   but we cannot be certain it is always specified like this .*/
	 
      if (lon0to360) {    /* target lon is positive */
        if (cdf.xmax < 0) {
	   cdf.xmax += 360.0;
           if (cdf.xmin < 0)
	      cdf.xmin += 360.0;
         }
      }
      else {              /* target lon is negative */
           if (cdf.xmin > 0) {
	      cdf.xmin -= 360.0;
              if (cdf.xmax > 0)
	         cdf.xmax -= 360.0;
	   }
      } 

      if ((cdf.xmin > lon) || (cdf.xmax < lon) 
      || (cdf.ymin > lat) || (cdf.ymax < lat)) {
         cdf_close(cdfid);
	 goto NEXTFILE2;
      }


      error = get_indices(&cdf, lat, lon, &row, &col);
      
      if (row < 0 || row >= cdf.ny || col < 0 || col >= cdf.nx) {
         cdf_close(cdfid);
         goto NEXTFILE2;
      }
      
/* check that cdf file has same number of standard depths  */
     
      if (cdf.nz != NSTDLEVS) {
         fprintf(stderr, "\nFATAL ERROR:  Mismatch of standard depths in cdf files.\n");  
         exit(1);
      }

   read_cdf_bottom_depth(cdfid, &bd, row, col, tbin);
   if (is_flagged(bd, cdf.fill_value) ) 
       bd = HB_f_empty;
   
   maxlev = bindx;     /* initially don't include bottom obs in profile unless...*/

   if ((bd > flag) &&  (bd < maskflag) && (bd > deepest)) { 
      deepest = bd;       /* this profile is the deepest */   
      maxlev = NSTDLEVS;  /* include bottom obs in profile */
      ppro[bindx] = 0.0;   /* zero previous bottom obs to handle the case */
      pcount[bindx] = 0;   /* of multiple entries for same square */
      tpro[bindx] = 0.0;
      tcount[bindx] = 0;  
      spro[bindx] = 0.0;
      scount[bindx] = 0;  
   }
   
/******************/      
/* extract pr, t90, sa at all depth levels.  For now, set masked or missing values
       to zero  */
    
     t90_avail = 0;
     for (i = 0; i < cdf.nprops; ++i) {
        if (get_prop_indx(cdf.prop_id[i]) == (int)T90)
	t90_avail = 1;
     }
    if (t90_avail)
         error = read_cdf_prop(cdfid, "t90", xin, row, col, tbin, 0, cdf.nz);
     else  {
       error = read_cdf_prop(cdfid, "te", xin, row, col, tbin, 0, cdf.nz);
       for (lev = 0; lev < cdf.nz; ++lev) {   /* convert to T90 */
          if (! (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask)) )
	      xin[lev] *= 0.99976;
       }
     }
     
     if (cdf.counts_included) {
        if (t90_avail)
            error = read_cdf_prop_count(cdfid,"t90",nobs,row,col,tbin,0,cdf.nz);
	 else
	    error = read_cdf_prop_count(cdfid,"te",nobs,row,col,tbin,0,cdf.nz);
	   
	   if (error > 0)  {
               fprintf(stderr,"\nError attempting to read te_cnt from cdf file.");
               exit(1);
	    }
     }
     
     for (lev = 0; lev < maxlev; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
	  nobs[lev] = 0;
       }  
       else {
          if (cdf.counts_included)
	     tcount[lev] += nobs[lev];
	  else {
             ++tcount[lev];
	     nobs[lev] = 1;
	  }
          tpro[lev] += (double) (xin[lev] * nobs[lev]);
       }
     }
  
     error = read_cdf_prop(cdfid,"sa", xin, row, col, tbin, 0, cdf.nz);
     if (error) {
       fprintf(stderr, "\n No sa variable in cdf file:  %s\n", argv[curfile]);
       exit(1);
     }
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,"sa",nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read sa_cnt from cdf file.");
            exit(1);
        }

     }
     for (lev = 0; lev < maxlev; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
	  nobs[lev] = 0;
       }  
       else {
          if (cdf.counts_included)
	     scount[lev] += nobs[lev];
	  else {
             ++scount[lev];
	     nobs[lev] = 1;
	  }
          spro[lev] += (double) (xin[lev] * nobs[lev]);
       }

     }
     error = read_cdf_prop(cdfid,"pr", xin, row, col, tbin, 0, cdf.nz);
     if (error) {
       fprintf(stderr, "\n No pr variable in cdf file:  %s\n", argv[curfile]);
       exit(1);
     }
     if (cdf.counts_included) {
        error = read_cdf_prop_count(cdfid,"pr",nobs,row,col,tbin,0,cdf.nz);
        if (error > 0) {
            fprintf(stderr,"\nError attempting to read pr_cnt from cdf file.");
            exit(1);
        }

     }
     for (lev = 0; lev < maxlev; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
	  nobs[lev] = 0;
       }  
       else {
          if (cdf.counts_included)
	     pcount[lev] += nobs[lev];
	  else {
             ++pcount[lev];
	     nobs[lev] = 1;
	  }
          ppro[lev] += (double) (xin[lev] * nobs[lev]);
       }
    }
  
     
     cdf_close(cdfid);

NEXTFILE2:
      ;	 
   }  while (curfile++ < nfiles); 
   
/*****************************************************/ 
/* find average in case there were multiple entries */

   for (lev = 0; lev < NSTDLEVS; ++lev) {
   
     if (pcount[lev])
       ppro[lev] /= pcount[lev];
     else 
       ppro[lev] = HB_f_empty;
   
     if (tcount[lev])
       tpro[lev] /= tcount[lev];
     else 
       tpro[lev] = HB_f_empty;
   
     if (scount[lev])
       spro[lev] /= scount[lev];
     else 
       spro[lev] = HB_f_empty;
   
   }

   free((void *)xin);
   free((void *)nobs);
   
   
  /* gridpoint is masked */
  
   if ( is_flagged(bd, HB_f_mask) && deepest < 0.0) {
   
      *bdepth_ptr = HB_f_mask;
   
       for (lev = 0; lev <= bindx; ++lev) {
          s0pro[lev] = HB_f_mask;
          s1pro[lev] = HB_f_mask;
          s2pro[lev] = HB_f_mask;
          s3pro[lev] = HB_f_mask;
          s4pro[lev] = HB_f_mask;
          ppro[lev] = HB_f_mask;
          tpro[lev] = HB_f_mask;
          spro[lev] = HB_f_mask;
      }
      return;
   }
     
  /* Not masked......compute sigmas */
  
   compute_sigma(0.0, cdf.nz, s0pro, ppro, tpro, spro);
   compute_sigma(1000.0, cdf.nz, s1pro, ppro, tpro, spro);
   compute_sigma(2000.0, cdf.nz, s2pro, ppro, tpro, spro);
   compute_sigma(3000.0, cdf.nz, s3pro, ppro, tpro, spro);
   compute_sigma(4000.0, cdf.nz, s4pro, ppro, tpro, spro);

/* replace missing values with empty flag */

   for (lev = 0; lev < bindx; ++lev) {
        if (s0pro[lev] < 0) {
             s0pro[lev] = HB_f_empty;
             s1pro[lev] = HB_f_empty;
             s2pro[lev] = HB_f_empty;
             s3pro[lev] = HB_f_empty;
             s4pro[lev] = HB_f_empty;
             ppro[lev] = HB_f_empty;
             tpro[lev] = HB_f_empty;
             spro[lev] = HB_f_empty;
        }
   } /* end for lev */


    *bdepth_ptr = bd;  
    if (deepest > 0.0)
         *bdepth_ptr = deepest;
	 
          
   if (is_flagged(deepest, HB_f_empty)) {        /* is bottom missing? */
       s0pro[bindx] = HB_f_empty;
       s1pro[bindx] = HB_f_empty;
       s2pro[bindx] = HB_f_empty;
       s3pro[bindx] = HB_f_empty;
       s4pro[bindx] = HB_f_empty;
       ppro[bindx] = HB_f_empty;
       tpro[bindx] = HB_f_empty;
       spro[bindx] = HB_f_empty;
       return;      
   }

  /* flag values below the seafloor depth with HB_f_mask */
  
   for (lev = 0; lev < bindx; ++lev) {
      if (std_depth[lev] >= *bdepth_ptr) {
          s0pro[lev] = HB_f_mask;
          s1pro[lev] = HB_f_mask;
          s2pro[lev] = HB_f_mask;
          s3pro[lev] = HB_f_mask;
          s4pro[lev] = HB_f_mask;
          ppro[lev] = HB_f_mask;
          tpro[lev] = HB_f_mask;
          spro[lev] = HB_f_mask;
      }
   } /* end for lev */
    
 /* All profiles should now contain either a value, HB_f_empty or HB_f_mask */
         
   return;
    
} /* end get_sig_profiles() */

	   
/****************************************************************************/

void define_mixed_layer( double **sigwork, double *weights, struct GRID_INFO *hptr, int xrad, int yrad, float bdepth, struct MIXLAYER *mixlayptr)
     /* Determines density of mixed layer from average of surrounding sigwork
       profiles. Returns density and depth of mixed layer in the structure. 
    */     
{
   double delta, maxsigma, sigma, *xsurf , *mldepths;
   int j, nwsq, sq, sqmid, n;
   
   sqmid = yrad * hptr->nx + xrad;
   delta = 0.02;   /* arbitrary definition of mixed layer sigma range */
   mixlayptr->density = (double) HB_f_empty;
   mixlayptr->depth =  (double) HB_f_empty;
   mixlayptr->nobs = 0;
      
 /* If surface value of sigma profile being smoothed is flagged... flag the mixed layer and return */   
    j = 0;
   if (sigwork[sqmid][j] < 0.0 || sigwork[sqmid][j] > 100.0) {  
      mixlayptr->density = sigwork[sqmid][j];
       return;
   }

   /* Otherwise, get average mixed layer params from surrounding grid nodes */ 
    
   nwsq = hptr->ny * hptr->nx;
   xsurf = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   mldepths = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
    n = 0;
   for (sq = 0; sq < nwsq; ++sq) {
 
      j = 0;
      sigma = sigwork[sq][j];    
      xsurf[sq] = sigma;    
      mldepths[sq] = sigma;   /* in case it is flagged */
      
      /* find depth range over which density is within delta units of surface value*/  
      if (sigma > 0.0 && sigma < 100.0) {
        maxsigma = sigma + delta;

	++j;
        while ((j < NSTDLEVS) && (sigwork[sq][j] <= maxsigma) && (sigwork[sq][j] > 0.0))
	    ++j;
	    
	 mldepths[sq] = std_depth[--j]; 
	 ++n;
      }
      
   }  /* end for sq */
   
   zero_weights_2d(xsurf, weights, HB_f_empty, HB_f_mask, hptr, xrad, yrad);
   mixlayptr->density = weighted_mean(xsurf, weights, HB_f_empty, HB_f_mask, nwsq);
   mixlayptr->depth =  weighted_mean(mldepths, weights, HB_f_empty, HB_f_mask, nwsq);
   mixlayptr->nobs = n;
   
   free(xsurf);
   free(mldepths);
   return;
   
}  /* end define_mixed_layer() */
/************************************************************************/


void get_neighbors(double isopyc, double *xsurf, double depth_being_interpolated, float *bdpth, int nsq, struct GRID_INFO *hptr, double **xwork2, double **sigwork2)

/*
   Traverse each gridnode in work arrays. The bdpth array has the same
   dimension as the work grids and x grid -- it contains the seafloor
   depth at each gridnode. Interpolate to find prop value 
   associated with specified isopycnal and store it in xsurf[] array.
   If isopycnal runs into bottom or outcrops, mark xsurf[sq] with HB_f_mask.
*/

{
   double *d, *x, *sig, xval;
   double *sigptr, *xptr;
   double emptyflag;
   double diff, lastdiff, zval, best_zval;
   int n, nz, lev, sq, j;
   int no_bottom, datagap;

   if (depth_being_interpolated < 0) {
       fprintf(stderr,"FATAL ERROR:  depth_being_interpolated is negative in get_neighbors()\n");  
    } 

/* Allocate memory to hold profiles */

   d = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   x = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   sig = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   
   emptyflag = HB_f_empty + 10.0;
   nz = NSTDLEVS - 1;
   
   for (sq = 0; sq < nsq; ++sq) {
      xptr = xwork2[sq];
      sigptr = sigwork2[sq];
   
      
      /* load  work profiles into local profile array and weed out
         missing values for interpolation function */
	 
      n = 0;
      for (lev = 0; lev < nz; ++lev) {
        if (sigptr[lev] > 0.0 && sigptr[lev] < 100.0 && xptr[lev] > emptyflag) {
	  d[n] = std_depth[lev];
	  x[n] = xptr[lev];
	  sig[n] = sigptr[lev];
	  ++n;
	}
      }
      
      /* check seafloor for valid observation */
       
      no_bottom = 1;     
      if (sigptr[nz] > 0.0 && sigptr[nz] < 100.0 && xptr[nz] > emptyflag) {
          d[n] = bdpth[sq];  /* insert seafloor depth for this square */
	  x[n] = xptr[nz];
	  sig[n] = sigptr[nz];
	  ++n;
	  no_bottom = 0;
      }
      
      xsurf[sq] = HB_f_empty;
           
      if (n == 0) {   /* no values at this grid node */
      
         if (bdpth[sq] > 99999.0)  {  /* node is masked */
	    xsurf[sq] = (double) HB_f_mask;
	 }
      }
      
      else  {   /* check for all crossings and use closest to this stddepth */
 
           j = 0;
	   lastdiff = (double)  HB_f_mask;  /* initialize with very large number */ 
	   best_zval =  (double) HB_f_empty;  /* initialize with very negative number */ 
	   do {
               zval = hb_linterp(isopyc, &sig[j], &d[j], n-j);
	       if (zval > -9990.) {
	       
	           if ((diff = ABS(zval - depth_being_interpolated)) < lastdiff ) {
		            best_zval = zval;
			    lastdiff = diff;
		    }
		    
		    while ((d[j] <= zval) && (j < n))
	               ++j;
	       }
           } while (zval > -99990. && (j < n));
	   
	   if (best_zval >= 0 ) {   /* found an appropriate crossing */
	     
	       /* if this is at the top of a pycnostad, get as close as possible 
	       to depth being interpolated */
	       
	       if (lastdiff > 300) {
	          zval = best_zval;
	          diff =  depth_being_interpolated - best_zval;
		  if (diff > 0)  {   /* zval is located > 300 meters above the level being interpolated for*/
		     j = 0;
		    while ((d[j] <= zval) && (j < n))
	               ++j;
		     /* now j points to level just below best_zval */
		     
		     while ((j < n) && (d[j] < depth_being_interpolated) && ABS((sig[j] - isopyc) < 0.005) )
		         ++j;
			 
			 /*now j points to level closest to depth_being_interpolated that is 
	           within 0.005 sigma units of the isopycnal value being sought */
		       
		       if (j == n)  /* should never happen, but check anyway */
		             --j;
			 
		       best_zval = d[j];
		  } /* end if diff > 0 */
		  
		/*  else  best_zval is below depth being interpolated  and so can't get any closer  
		(assumes that density is monotonically increasing) */
		  
	       } /* end if lastdiff > 300 */
	       
	       zval = best_zval;

	       
               /* get xvalue and check for datagaps around this depth...*/
	    
	            xval = hb_linterp(zval, d, x, n);
                    j = 0;
	            while (d[j] <= zval && j < n)
	               ++j ;
	      
	            if  (j == n  && d[j] == zval)
	                 datagap = 0;
		     else if (j == 0)
		         datagap = 0;
	            else if ((d[j-1] == zval) || (d[j] == zval))
	                datagap = 0;
	            else if (zval < 1001)
	                datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
                    else
                        datagap = (d[j] - d[j-1]) > GAP_DEEP;
		
	            if (!datagap)
	                xsurf[sq] = xval;
	    }
	    else {  /* no crossing -- does it outcrop?*/
	    
	       if (!no_bottom && (isopyc > sig[n-1])) {
	             xsurf[sq] = HB_f_mask;   /* isopycnal runs into sea floor */
	       }
	    } 
      }   /* end else check for all crossings */
      
   }  /* end for sq */

   free(d);
   free(x);
   free(sig);
   return;

}  /* end get_neighbors() */
	   
	   
/*********************************************************/	   
/****************************************************/

void x_to_stddepth(double *xin, short *xcnt, float *xout, double *din, double *stdd, float btmd, int nlevs)
/* interpolates x,depth arrays (xin, din) onto depth levels in stdd.  
        xcnt :  < 0 for values that were fitted, > 0 for levels that did not need fitting,
	         = 0 for levels that are masked or missing. 
        btmd :  bottom depth. 
    Uses xcnt to determine which values need to be interpolated:  
    Returns array  xout ready to be written to netcdf file: depths below btmd are assigned HB_f_mask.
    Missing values are assigned HB_f_empty.
*/
{
   double *dtmp, *xtmp,  xval;
   double emptyflag, prevdepth;
   int n, nz, lev, j;
    
/*  case where there is no ocean */
   if ( is_flagged(btmd, HB_f_mask) ) {
      for (lev = 0; lev < nlevs; ++lev)  {
          if (xcnt[lev] != 0) {
	    fprintf(stderr," >> SOFTWARE BUG:  masked square has non-zero count for fitted values.\n");
	    xcnt[lev] = 0;
	  }
          xout[lev] =  (float) xin[lev];
      }	  
      return;
   }
   
/* Allocate memory to hold profiles */

   dtmp = (double *) get_memory((void *)NULL, (size_t) nlevs, sizeof(double));
   xtmp = (double *) get_memory((void *)NULL, (size_t) nlevs, sizeof(double));

   nz = nlevs - 1;    /* index to deepest level */
   emptyflag  = HB_f_empty + 10;
   
      /* load  work profiles into local profile array and weed out
         missing values for interpolation function */
      prevdepth = -10.0; 
      n = 0;
      for (lev = 0; lev < nz; ++lev) {
        if (xcnt[lev] != 0 ) {
	  dtmp[n] = din[lev];
	  if (dtmp[n] < prevdepth) {
/*	      fprintf(stderr," >> MESSAGE from x_to_stddepth(): depth array not monotonically increasing. prevdepth: %.1lf   thisdepth: %.1lf.\n",  prevdepth, dtmp[n]); */
	  }
	  else {
	     prevdepth = dtmp[n] -5.0;
	     xtmp[n] = xin[lev];
	     ++n;
	  }
	}
      }
      
      if (xcnt[nz] != 0) {     /* add bottom level */
         dtmp[n] =  din[nz];
	  if (dtmp[n] < prevdepth) {
/*	      fprintf(stderr," >> MESSAGE from x_to_stddepth(): bottom depth not monotonically increasing.  prevdepth: %.1lf   thisdepth: %.1lf\n", prevdepth, dtmp[n]); */
	  }
	 xtmp[n] = xin[nz];
	 ++n;
      }
           
      if (n < 2) {  /* not enough xvalues to interpolate, just copy input array to output array */
          for (lev = 0; lev <= nz; ++lev)  
	      xout[lev] =  (float) xin[lev];
	      
	  free(dtmp);
	  free(xtmp);
	  return;
      } 
      
      /* interpolate values onto stddepths. */
      
      for (lev = 0; lev < nz; ++lev) {
      
         if (xcnt[lev] == 0)  /* level is flagged */
	      xout[lev] = (float) xin[lev];
	 else {
	      xval = hb_linterp(stdd[lev], dtmp, xtmp, n);
	   
	      xout[lev] = (float) xval;
	      if (xval < -9999.0){
/*                 fprintf(stderr," >> SOFTWARE BUG in x_to_stddepth(): unable to interpolate a value back onto stddepth= %.1lf ! bottom depth = %.1f\n", stdd[lev], btmd); */
	      
	          xout[lev] = xin[lev];
	      }
	 }
       }
       
       xout[nz] = (float) xin[nz];   /* explicity set, don't interpolate, deepest value */
      free(dtmp);
      free(xtmp);
      return;

} /* end x_to_stddepth */

