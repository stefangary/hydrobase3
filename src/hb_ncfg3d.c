/*  hb_ncfg3d.c

................................................................................
                          *******  HydroBase 3 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Dec 2000
			     Updated for t90 Feb 2012
			     
...................................................
*
*  Interpolates for missing values in 3-dimensional grids of 
*  hydrographic properties using an interative Laplacian/spline fitting
*  algorithm.
*  Reads lat/lon/depth gridded fields of properties in HydroBase cdf format
*  and a binary topography file.   Points that are missing (determined from 
*  the topography data) are filled in by estimating an isopycnal value
*  value for the point (from surrounding data) and then fitting each
*  property on that isopycnal surface. Masked nodes are not
*  incorporated into the fit.
*  A x-/y- Search radius parameter determines the maximum acceptable distance an
*  interpolated value can be from an observed value. If a masked node is
*  encountered while searching along radius, the search is suspended in 
*  that direction.  If no observation is within that distance, the node is 
*  flagged as "empty". Depths below the seafloor depth are flagged as "masked".
*  
*  Output is a HydroBase netcdf file.  

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
#include "zgrid.h"
#include "hb_memory.h"
#include "hb_paths.h"

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""


#define    PREFILL      1        /* turn off prefilling of cdf variables */

/* global variables */

float HB_f_mask;          /* float values are used since the output */
float HB_f_empty;         /* cdf values are of this type */
double HB_zgrid_mask;  /* zgrid() expects this on input to mark masked nodes */
int verbose;
int add_stats;

/* needed by print_usage() */
    
double converge_limit;     /*  make this small */
int max_iterations;
double tension;
double xradius, yradius;


/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *, int *);
int cdf_construct(struct GRID_INFO *, int, int *, char *, int, char **);
void set_topo_info(struct GRID_INFO *, struct GRID_INFO *, int);   
void get_x_profile(float, float, int, int, char **, char *, char *, double *, short *, float);
int define_bottom(float *, short);
void get_sig_profiles(float, float, float *, short, int, char **, char *, char *, double *, double *, double *, short *, short *, short *, double *, double *, double *, double *, double *);
double find_isopycnal(int, int, int, int, double **);
double find_bottom_isopycnal(int, int, int, double **, float *, float);
void get_neighbors(int, double, double,  struct POINT *, int *, int *, double *, float *, int , struct GRID_INFO *, double **, double **);



main(int argc, char **argv)
{
  int i, j, n, nprops, allprops;
  int nfiles, nz, nsqout, nwsq1, nwsq2;
  int row, col, sq, sq_out, error, tbin;
  int wsq1, wsq2, wsq_mid;
  int wrow, wcol, wcol1, wrow1;
  int do_fit, gotprops, is_te, is_sa;
  int bflag, iflag, oflag, uflag, pflag;
  int infile, outfile, print_msg; 
  int xrad, yrad; 
  int pixel_grid;           /* pixel or gridnode registration */
  int n_empty, n_mask, n_set;
  int noldpts;
  int nalloc;
  int mask_it;
  int *prop_indx, *all_prop_indx, *prop_req, *prop_avail;
  float wlat, wlon, lat, lon, wlat2, wlon2;
  double toobig;          /* zgrid() representation of empty/masked nodes */
  double bigger;
  double xoffset, yoffset;
  double topo_xinc, topo_yinc;
  double isopyc_above, curr_depth;
  char modifier;
  char *st;
  char *outfile_name;
  char *dir, *extent;
  struct CDF_HDR *cdfin;
  struct POINT *old_data; 
  struct GRID_INFO hwork1, hout, hwork2;
  struct GRID_INFO topo_info;
  int *knxt, *imnew;
  short *seafloor, topoval;
  double *topolat, *topolon;
  double *new_grid;
  double *xin, *zpij;
  short *new_count;
  short **tcount1, **scount1, **pcount1;
  short **xcount1;
  float *bottom_de, *bdwork2, *bdwork1;
  float *xout, *xtmp;
  double **xwork1, **pwork1;
  double **swork1, **twork1;
  double **sig0work1, **sig1work1;
  double **sig2work1, **sig3work1;
  double **sig4work1;
  double **sigptr, **xwork2;
  double **sig0work2, **sig1work2;
  double **sig2work2, **sig3work2;
  double **sig4work2;
  double **isopyc, **sfit;
  float **pfit;

/* Set these default values. */

  toobig = 1.0E34;        /* flag is a  smaller than zgrid() flag */	
  bigger = 1.0E36;        /* larger than empty flag/smaller than mask flag */		      
  HB_f_mask = HBMASK;      /* HydroBase cdf file mask flag */
  HB_f_empty = HBEMPTY;    /*  and empty flag */
  HB_zgrid_mask = DBL_MAX; /* a machine dependent value */
    
  error = 0;
  bflag = iflag = oflag = uflag = pflag =  0;
  max_iterations = 10;
  converge_limit = 0.0000001;
  tension = 0.5;
  xradius = 2.0;          /* default search radius in degrees */
  yradius = 2.0;
  pixel_grid = 1;          /* default pixel registration */
  topo_xinc = topo_yinc = 0.1;
  dir = DIR;
  extent = EXTENT;
  nfiles = 0;
  print_msg = 1;
  tbin = 0;
  verbose = 0;
  add_stats = 1;   /* global variable determines whether variance is output */
  
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
                        	     
	  if (&hout.x_min > &hout.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (&hout.y_min > &hout.y_max) {
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
	  
        case 'O':
	  oflag = 1;
          outfile_name = &argv[i][2];
          break;
	  	      
        case 'P':
          pflag = 1;
	  prop_req = (int *) calloc(MAXPROP, sizeof(int));
          nprops = parse_p_option(&argv[i][2], prop_req);
          break;

        case 'Q':
          error = sscanf(&argv[i][2],"%d", &max_iterations) != 1;
          break;

        case 'S':
          error = sscanf(&argv[i][2],"%lf", &xradius) != 1;
          yradius = xradius;
          st = &argv[i][2];
          while (*(st++) != '\0') {
            if (*st == '/') {
               ++st;
               error = (sscanf(st,"%lf", &yradius) == 1) ? 0 : 1;
               break;
            }
          }
          break;

        case 'T':
	  tension = atof(&argv[i][2]);
          break;
	  
        case 'V':
	  verbose = TRUE;
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
   if (tension < 0. || tension > 1.0) {
	    fprintf(stderr,"ERROR:  specify tension parameter between 0->1\n");
	    ++error;
   }
   if (max_iterations < 1) {
      fprintf (stderr, "SYNTAX ERROR -Q option.  Must specify interations >= 1\n");
      error++;
   }
   if (xradius < 0. || yradius < 0.) {
      fprintf (stderr, "SYNTAX ERROR -S option.  Must specify a positive search radius.\n");
      error++;
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
/* Check for requested properties availability.    After this, prop_indx[] and
    nprops will hold the property info for the output cdf file. */
   
   prop_indx = (int *) calloc((size_t)nprops, sizeof(int)); 
   prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int)); 
   for (i = 0; i < cdfin->nprops; ++i)    
      prop_avail[get_prop_indx(cdfin->prop_id[i])] = 1;
   
   prop_avail[(int)T90] = prop_avail[(int)TE] || prop_avail[(int)T90];
	    
   n = 0;   
   for (i = 0; i < nprops; ++i) {
      if ( prop_avail[prop_req[i]] ) {
        prop_indx[n++] = prop_req[i];
      }
      else {
        if (prop_req[i] == (int) PR || prop_req[i] == (int) T90 || prop_req[i] == (int) SA)
          fprintf(stderr,"\n FATAL ERROR!! Property %s not available in cdf file (pr te sa are mandatory).\n", get_prop_mne(prop_req[i]));       
	else
          fprintf(stderr,"\n WARNING!! Property %s not available in cdf file and will not be output.\n", get_prop_mne(prop_req[i]));       
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
      allprops = n * 2;
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
/* Adjust order of properties so that sa is fitted before t90 -- since 
   salinity is needed to convert potential temperature back to in situ temperature. 
*/
   prop_indx[1] = (int)SA;
   prop_indx[2] = (int)T90;
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

/* save these values in gridnode units */

   xrad = NINT((xradius / hout.x_inc) + .0001 );  
   yrad = NINT((yradius / hout.y_inc) + .0001 );
           
  if (pixel_grid)
     fprintf (stderr, "\nUsing %s registration", "pixel");
  else
     fprintf (stderr, "\nUsing %s registration", "gridnode");
   fprintf (stderr, "\nGrid dimensions are nx = %d, ny = %d", hout.nx, hout.ny);
   fprintf (stderr, "\nOutput xmin/ymin gridnode: %8.3lf/%8.3lf ", hout.x_min, hout.y_min);
   
      
   fprintf (stderr, "\nNumber of gridnodes in search radii: %d/%d", xrad, yrad);
   xradius = xrad * hout.x_inc;
   yradius = yrad * hout.y_inc;
   fprintf (stderr, "\n Search radius must be a multiple of grid increment -- using: %.3lf/%.3lf ", xradius, yradius);
/*-------------------------------------------*/   
/* Set up big work grid -- output grid + extended bounds  */

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
   
   if ( ( i = (NINT((hwork1.x_max - hwork1.x_min) /  hwork1.x_inc + .0001) + 1) ) != hwork1.nx ) {
      fprintf (stderr, "\nMismatch between number of columns in work1 grid %d  %d", i, hwork1.nx);
      exit(1);
   }
   if ( ( i = (NINT((hwork1.y_max - hwork1.y_min) /  hwork1.y_inc + .0001) +1 )) != hwork1.ny ) {
      fprintf (stderr, "\nMismatch between number of rows in work1 grid  %d  %d", i, hwork1.nx);
      exit(1);
   }


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

   sfit = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   pfit = (float **) get_memory((void *)NULL, (size_t)nsqout, sizeof(float *));
   isopyc = (double **) get_memory((void *)NULL, (size_t)nsqout, sizeof(double *));
   for (sq = 0; sq < nsqout; ++sq)
      isopyc[sq] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
  

/*-------------------------------------------*/   
/* Set up small work grid to store isopycnal surface sent to zgrid */

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

 /*-------------------------------------------*/   
/*-------------------------------------------*/   
  /* Read in topography values */
  
    set_topo_info(&hwork1, &topo_info, 0);
    seafloor = hb_get_topo(BATHPATH_C, &topo_info, &topolat, &topolon, FALSE, topo_info.lon0to360, &topoval);
/*-------------------------------------------*/ 

  /* set up array to store bottom depths for output and work */
  
   bottom_de = (float *) get_memory((void *)NULL, (size_t) nsqout, sizeof(float));
   bdwork1 = (float *) get_memory((void *)NULL,(size_t)nwsq1, sizeof(float));
   bdwork2 = (float *) get_memory((void *)NULL,(size_t)nwsq2, sizeof(float));
/*-------------------------------------------*/   
/* Visit each gridpoint in big workspace, read in profiles, compute sigmas. */

   fprintf(stderr,"\nReading in profiles...");
   
   for (row = 0; row < hwork1.ny; ++row ) {
      lat = hwork1.y_min + row * hwork1.y_inc;
      
      for (col = 0; col < hwork1.nx; ++col) {
         lon = hwork1.x_min + col * hwork1.x_inc;
  
	 sq = row * hwork1.nx + col;
	 topoval = find_nearest_topo_val(lat, lon, seafloor, &topo_info);
	 get_sig_profiles(lat, lon, &bdwork1[sq], topoval, nfiles, argv, dir, extent, pwork1[sq], twork1[sq], swork1[sq], pcount1[sq], tcount1[sq], scount1[sq], sig0work1[sq], sig1work1[sq], sig2work1[sq], sig3work1[sq], sig4work1[sq]);
	   
	  /* Convert in situ temperatures to potential temperatures for smoothing */

         for (j = 0; j < NSTDLEVS; ++j) {
           if (twork1[sq][j] > -3.0 && twork1[sq][j] < 100.)
	      twork1[sq][j] = hb_theta(swork1[sq][j], twork1[sq][j], pwork1[sq][j], 0.0);
	 }
	 
      } /* end for col */
   } /* end for row */
   

   fprintf(stderr,"\nNow filling empty nodes ...  ");


/* Loop for each property except pressure */

   for (i = 1; i < allprops; ++i) {

       is_sa = (all_prop_indx[i] == (int)SA);
       is_te = (all_prop_indx[i] == (int)T90);
       
       fprintf(stderr,"\n%s ", get_prop_mne(all_prop_indx[i]));
      
       switch ((enum property) all_prop_indx[i]) {
       
 	  case T90 :          /* just set the pointers */
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
	      wrow = sq / hwork1.nx;
	      wcol = sq - wrow * hwork1.nx;
	      wlat = hwork1.y_min + wrow * hwork1.y_inc;
	      wlon = hwork1.x_min + wcol * hwork1.x_inc;
              xwork1[sq] = (double *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(double));
              xcount1[sq] = (short *) get_memory((void *)NULL,(size_t)NSTDLEVS, sizeof(short));
	      get_x_profile(wlat, wlon, all_prop_indx[i], nfiles, argv, dir, extent, xwork1[sq], xcount1[sq], bdwork1[sq]);
	   }
       }  /* end switch */
	           
        /* Visit each gridpoint in output file and determine whether 
	gaps need to be filled in by fitting.  Produce output profile 
	that includes mask and empty flags where appropriate.  
	Use lat/lon (not row/col) to translate between grids for 
	output/input, seafloor, and work because each has its own 
	dimension and order. */  


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
		   
	           wcol1 = NINT((wlon2 - hwork1.x_min + 0.0001) / hwork1.x_inc);
	           wrow1 = NINT((wlat2 - hwork1.y_min + 0.0001) / hwork1.y_inc);
	           wsq1 = wrow1 * hwork1.nx + wcol1;   /* index to big work grids */
	      
		   xwork2[wsq2] = xwork1[wsq1];
	           sig0work2[wsq2] = sig0work1[wsq1];
	           sig1work2[wsq2] = sig1work1[wsq1];
	           sig2work2[wsq2] = sig2work1[wsq1];
	           sig3work2[wsq2] = sig3work1[wsq1];
	           sig4work2[wsq2] = sig4work1[wsq1];
	    
	           bdwork2[wsq2] = bdwork1[wsq1];  /* store value, not ptr */
		   
		   if (wsq2 == wsq_mid)
		     new_count = xcount1[wsq1];
	      
	         } /* end for wcol */
	      }  /*end for wrow */
	      
	      bottom_de[sq_out] = bdwork2[wsq_mid];


             xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));

             do_fit = 0;
	     mask_it = is_flagged(bottom_de[sq_out], HB_f_mask);     
	     for (j = 0; j < NSTDLEVS; ++j) {
                if (mask_it) 
		  xout[j] = HB_f_mask;
		else {
	           xout[j] = (float) xwork2[wsq_mid][j];
		   if (is_flagged(xout[j], HB_f_empty) )
		      ++do_fit;
		}
	     }
	   

	     if (do_fit) {

            /* Loop for each depth level including bottom ...
             find corresponding isopycnal value from distance-weighted avg
	     of surrounding nodes. Then interpolate to find
	     that isopycnal at all the neighboring gridnodes.  Construct
	     an array of old_data to reflect the xproperty value on that
	     isopycnal at each of the neighboring gridnodes and send it to
	     zgrid().  Use the value returned at newgrid[wsq_mid] as the 
	     fitted value to be output. */

               isopyc_above = 0.0;		 
	       for (j = 0; j < NSTDLEVS; ++j) {
	      
		   if (j != nz) {
		      if (std_depth[j] <= 500)
			sigptr = sig0work2;
		      else if (std_depth[j] > 500 && std_depth[j] <= 1500)
		        sigptr = sig1work2;
		      else if (std_depth[j] > 1500 && std_depth[j] <= 2500)
		        sigptr = sig2work2;
		      else if (std_depth[j] > 2500 && std_depth[j] <= 3500)
		        sigptr = sig3work2;
		      else 
		        sigptr = sig4work2;
	           }
		   else {  /* special case for bottom */
		         if (bottom_de[sq_out] <= 500)
			   sigptr = sig0work2;
		         else if (bottom_de[sq_out] > 500 && bottom_de[sq_out] <= 1500)
		           sigptr = sig1work2;
		         else if (bottom_de[sq_out] > 1500 && bottom_de[sq_out] <= 2500)
		           sigptr = sig2work2;
		         else if (bottom_de[sq_out] > 2500 && bottom_de[sq_out] <= 3500)
		           sigptr = sig3work2;
		         else 
		           sigptr = sig4work2;
		   }
		   
		if (std_depth[j] >= bottom_de[sq_out]) {
		   xout[j] = HB_f_mask;
		}
				   
		else if (! is_flagged((float) xwork2[wsq_mid][j], HB_f_empty) ) {
		   xout[j] =  (float) xwork2[wsq_mid][j];
		   isopyc[sq_out][j] = sigptr[wsq_mid][j];
		   isopyc_above = isopyc[sq_out][j];
		}
		   
		else {   /* isopycnally fit this level */
					   
		   if (isopyc[sq_out][j] == 0.0) {
		   
		      if (j != nz) {
		         isopyc[sq_out][j] = find_isopycnal(j, xrad, yrad, hwork2.nx, sigptr);
		         
		      }
		      else  {         /* special case for bottom */
		         isopyc[sq_out][j] = find_bottom_isopycnal(xrad, yrad, hwork2.nx, sigptr, bdwork2, bottom_de[sq_out]);
		      }
		   }
		   
		   if (isopyc[sq_out][j] <= 0) {
		      xout[j] = HB_f_empty;
		      continue;
		   }
		   else {
		   
                       if (isopyc[sq_out][j] < isopyc_above - 0.001)    /* check for density inversion */
		            isopyc[sq_out][j] = isopyc_above;
		            
		        isopyc_above = isopyc[sq_out][j];
		   }	
		   	   
		   curr_depth = std_depth[j];
	           if (j == nz)
		           curr_depth = bottom_de[sq_out];
			   
	           /* set up array for input to zgrid() */
		   
                   xin = (double *) get_memory((void *)NULL, (size_t) nwsq2, sizeof(double));
		   nalloc = nwsq2 + 100;
		   old_data = (struct POINT *) get_memory((void *)NULL,(size_t)nalloc, sizeof(struct POINT));
		   noldpts = 0;
                   get_neighbors(prop_indx[i], isopyc[sq_out][j], curr_depth, old_data, &nalloc, &noldpts, xin, bdwork2, nwsq2, &hwork2, xwork2, sigptr);
		   
		   xout[j] = HB_f_empty;
		   
		   if (noldpts > 0) {
		   
                    /* set up other arrays needed for zgrid() */
		    
                      knxt = (int *) get_memory((void *)NULL, (noldpts + 1), sizeof(int));
                      zpij = (double *) get_memory((void *)NULL, (noldpts + 1), sizeof(double));
                      imnew = (int *) get_memory((void *)NULL, hwork2.ny, sizeof(int));
		
	              new_grid = zgrid(converge_limit, max_iterations, noldpts, xrad, yrad, tension, &hwork2, knxt, imnew, zpij, xin, old_data, &n_mask, &n_empty, &n_set);
		   
	              free((void *)xin);
	              free((void *)knxt );
	              free((void *)zpij);
	              free((void *)imnew);
		   
                     if (new_grid[wsq_mid] >= toobig ) {
	                new_grid[wsq_mid] =  (double)HB_f_empty;
	             }
	             else {
	                if (new_count[j] <= 0) 
	                    new_count[j] = 1;
	             }
		   
		     xout[j] = (float) new_grid[wsq_mid];
		  
		     free((void *) new_grid);
		  }
		  
	          free ((void *)old_data);
		   
	         } /* end else */
		 
		 
	      } /* end for j */
	      
	   } /* end if do_fit */
	   
	   
	   if (is_sa) {  
	      
	           /* save salinity profile with fitted values and compute 
		      pressure.  These will be used to convert potential
		      temperatures back to in situ values. */
		      
	         sfit[sq_out] = (double *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(double));
	         pfit[sq_out] = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
		      
	         for (j = 0; j < nz; ++j) {
		     sfit[sq_out][j] = (double) xout[j];
		     pfit[sq_out][j] =  xout[j];
		     
	             if (sfit[sq_out][j] > -9.0 && sfit[sq_out][j] < 100.0) {
		       if (NINT(std_depth[j]) == 0)
			    pfit[sq_out][j] = 0.0;
	               else
		            pfit[sq_out][j] =  (float) hb_p80(std_depth[j], (double)lat);
			
			
		     }
	         }
		 
		 sfit[sq_out][nz] = (double) xout[nz];
		 pfit[sq_out][nz] = xout[nz];
	         if (sfit[sq_out][nz] > -9.0 && sfit[sq_out][nz] < 100.0)
	            pfit[sq_out][nz] = (float) hb_p80((double)bottom_de[sq_out], (double)lat);
		 
	       write_prop_cdf(outfile, pfit[sq_out], "pr", row, col, tbin, 0, 1,1,1,NSTDLEVS);
	       write_prop_count_cdf(outfile, new_count, "pr", row, col, tbin, 0, 1,1,1,NSTDLEVS);
		 
	   }
	      
	   if (is_te) {  /* adiabatically adjust potential temperature to pressure */
	         xtmp = xout;
	         xout = (float *) get_memory((void *)NULL, (size_t)NSTDLEVS, sizeof(float));
		 
	         for (j = 0; j < NSTDLEVS; ++j) {
	    
	             if (xtmp[j] > -3.0 && xtmp[j] < 100.0) {
	                xout[j] = (float) hb_theta(sfit[sq_out][j], (double)xtmp[j], 0.0, (double)pfit[sq_out][j]);
			if ( ! (xout[j] > -3.0 && xout[j] < 100.0) ) {
			   xout[j] = HB_f_empty;
			}
		    }
		     else if (xtmp[j] < -9.0) {
	                xout[j] = HB_f_empty;
	             }
	  	     else {
	                xout[j] = HB_f_mask;
                     }
	       
	         }  /* end for j */
	    
	         free((void *) xtmp);
		 free((void *) sfit[sq_out]);
		 free((void *) pfit[sq_out]);
	    
	   }  /* end if is_te */
	   
 	   if (all_prop_indx[i] < 700) {
	         write_prop_cdf(outfile, xout, get_prop_mne(all_prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	         write_prop_count_cdf(outfile, new_count, get_prop_mne(all_prop_indx[i]), row, col, tbin, 0, 1,1,1,NSTDLEVS);
           }
	   else {
	   
	         write_prop_err_cdf(outfile, xout, get_prop_mne(all_prop_indx[i]-700), row, col, tbin, 0, 1,1,1,NSTDLEVS);
	   }

	   free((void *)xout);
	   
         }  /* end for col */
       } /* end for row */
     
       
       if (!(is_te || is_sa)) {
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

  fprintf(stderr,"\nUses an iterative Laplacian/spline algorithm to isopycnally interpolate missing");
  fprintf(stderr,"\ngridpoints in HydroBase cdf files. Location of missing values are determined ");
  fprintf(stderr,"\nfrom the topography file. A search radius");
  fprintf(stderr,"\nparameter determines the maximum acceptable distance an");
  fprintf(stderr,"\ninterpolated value can be from an observed value.  The");
  fprintf(stderr,"\ninput domain is expanded relative to the output domain ");
  fprintf(stderr,"\nso that data around borders of the grid can be incorporated.");
  fprintf(stderr,"\nAll input netcdf files must possess the same depth dimensions");
  fprintf(stderr,"\nand grid spacing.  The output netcdf file will have the");
  fprintf(stderr,"\nsame depth dimension, but the lat/lon dimensions specified");
  fprintf(stderr,"\nwith -B<w/e/s/n> and -I<x/yincr>");
  
  fprintf(stderr,"\n\nUSAGE:  %s nc_file(s) -B<w/e/s/n> -I<x_inc>[/<y_inc>] -O<outfile> -P<properties> [-D<input_dir>] [-E<input_file_extent>] [-G] [-Q<max_iterations>] [-Ssearch_radius] [-T<tension>[b|i]] [-V] [-h]\n\n", program);
  fprintf(stderr," -B  sets the output grid bounds w/e/s/n.\n");
  fprintf(stderr," -I  sets the output grid spacing for x/y dimensions.\n");
  fprintf(stderr," -O  name of output cdf file.\n");
  fprintf(stderr," -P  list of properties to include in output file\n");
  fprintf(stderr,"        ex:  -Ppr/t90/sa/o2\n");
  fprintf(stderr,"       -P (by itself) produces a list of available properties\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr,"-D directory for input cdf files.\n");
  fprintf(stderr,"-E file extent for input cdf files.\n");
  fprintf(stderr,"-G force gridnode registration for ouput (default is pixel registration).\n");
  fprintf(stderr,"-Q sets the max iterations to achieve convergence.  [%d]\n", max_iterations);
  fprintf(stderr,"-S search for a non-empty gridnode within <xradius/yradius> of an empty node\n");
  fprintf(stderr,"      in the initial grid.  Estimate a value only if a non-empty node is found.\n");
  fprintf(stderr,"      x- yradii are specified in xy-units (degrees) .\n");
  fprintf(stderr,"      <radius> = 0 interpolates no gridnodes but puts in masking info.\n"); 
  fprintf(stderr,"      Default is [%.1lf/%.1lf]\n", xradius, yradius);
  fprintf(stderr,"-T tension parameter -- range [0..1]\n");
  fprintf(stderr,"      A value of 1 produces a pure laplacian solution,\n");
  fprintf(stderr,"      while a 0 value gives a harmonic spline solution\n");
  fprintf(stderr,"      with a smoother field but the possibility of \n");
  fprintf(stderr,"      spurious peaks or valleys. Default: [-T%lf]\n", tension);
  fprintf(stderr,"-V verbose.  \n");
  fprintf(stderr,"-h help...... prints this message. \n");
  return;
  
} /* end print_usage() */
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

     if ( !( (prop_indx[n] == (int)DE) || (prop_indx[n] == (int)PR) || (prop_indx[n] == (int)TE)|| (prop_indx[n] == (int)T90) || (prop_indx[n] == (int)SA)))  
            ++n;
     

   } while (*st == '/');
   return (n);
}  /* end parse_p_option() */

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
      cdfhdr.prop_units[i] = (char *) calloc(80,sizeof(char));
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
void set_topo_info(struct GRID_INFO *ginfo, struct GRID_INFO *tinfo, int topobuf)   
   /* for now, if pole is crossed, limit the y bound to the pole */
{
   int error, xpole;
   double dummy;

    
    if (topobuf == 0) {
        tinfo->x_min = ginfo->x_min;
	tinfo->y_min = ginfo->y_min;
        tinfo->x_max = ginfo->x_max;
	tinfo->y_max = ginfo->y_max;
	
     }
     else {
        error = pointB(ginfo->y_min, ginfo->x_min, 270, (double) topobuf, TRUE, &dummy, &tinfo->x_min);
        error = pointB(ginfo->y_max, ginfo->x_max, 90, (double) topobuf, TRUE, &dummy, &tinfo->x_max);
        
	xpole = pointB(ginfo->y_max, ginfo->x_max, 0, (double) topobuf, TRUE, &tinfo->y_max, &dummy);
	if (xpole) 
	   tinfo->y_max = 90.0;
	   
        xpole = pointB(ginfo->y_min, ginfo->x_min, 180,(double) topobuf , TRUE, &tinfo->y_min, &dummy);
	if (xpole)
	   tinfo->y_min = -90.0;
     
     }
     
     tinfo->lon0to360 = tinfo->x_min < 0 ? 0: 1;
     tinfo->xgreenwich = (tinfo->x_min < 0 && tinfo->x_max >= 0) ? 1 : 0;

     if (tinfo->y_min  <= -90)
	    tinfo->y_min = -89.999;
     if (tinfo->y_max  >= 90)
	    tinfo->y_max = 89.999;
     return;
} /* end set_topo_info() */


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
   int error, i, iprop, row, col;
   int tbin = 0, lev, new_bd, bindex;
   float bd;
   
   

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
         if (cdf.xmax < 0) {
	   cdf.xmax += 360.0;
           if (cdf.xmin < 0)
	      cdf.xmin += 360.0;
	 }
     
      }
      else {              /* target lon is negative */
           if (cdf.xmin > 0)  {
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
	 
     
     read_cdf_bottom_depth(cdfid, &bd, row, col, tbin);
     new_bd = 0;
     if (bd < (bdepth-10.0))
          new_bd = 1;    /* this bottom depth is not the real bottom */
      
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
            fprintf(stderr,"\nError attempting to read %s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
            exit(1);
        }
     }
 	
     /* check for flagged levels.  in case of multiple entries,
        sum up the data at each level.  */ 
	
     bindex = cdf.nz -1;
     
     if (pindex >= 700) {
     
         for (lev = 0; lev < bindex; ++lev) {
             if ( ! is_flagged(xin[lev], cdf.fill_value) && ! is_flagged(xin[lev], HB_f_mask)) {
                 if (nobs[lev] > 3) {
	            x[lev] += (double) (nobs[lev] * xin[lev]);
	            count[lev] += nobs[lev];
	         }
	     }
          } /* end for lev */
	  
          if (!new_bd) {
            if ( ! is_flagged(xin[bindex], cdf.fill_value)  && ! is_flagged(xin[bindex], HB_f_mask)) {
              if (nobs[lev] > 3) {
	         x[bindex] += (double) (nobs[bindex] * xin[bindex]);
	         count[bindex] += nobs[bindex];
	      }
	    }
          }
     }
     else {
     	   
        for (lev = 0; lev < bindex; ++lev) {
          if ( ! is_flagged(xin[lev], cdf.fill_value) && ! is_flagged(xin[lev], HB_f_mask)) {
       
	   if (!cdf.counts_included || nobs[lev] == 0) 
	      nobs[lev] = 1;
	   
	   x[lev] += (double) (nobs[lev] * xin[lev]);
	   count[lev] += nobs[lev];
	   }
        } /* end for lev */

       /* add bottom value only if it is a true bottom value */     

        if (!new_bd) {
            if ( ! is_flagged(xin[bindex], cdf.fill_value)  && ! is_flagged(xin[bindex], HB_f_mask)) {
	      if (!cdf.counts_included || nobs[bindex] == 0) 
	         nobs[bindex] = 1;
	      
	      x[bindex] += (double) (nobs[bindex] * xin[bindex]);
	      count[bindex] += nobs[bindex];
	    }
        }
     }
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
      else {
         x[lev] = HB_f_empty;
        if ( !is_flagged(bdepth, HB_f_empty) && !is_flagged(bdepth, HB_f_mask)) {	 
           if (std_depth[lev] >= bdepth)
                 x[lev] = (double) HB_f_mask; 
	} 
      }
   
   } /* end for lev */
   
      
   return;

}  /* end get_x_profile() */

/***************************************************/
int define_bottom(float *bd_addr, short seafloor)
  /* Determines an appropriate bottom depth based on
      observed data  plus seafloor topography and sets *bd_addr.
      Returns 1 if the bottom depth was changed or
      0 if not changed.
  */
{

  int new_bd;
   
  /* Compare seafloor depth and bottom depth.  
     IF bottom_de is defined: 
         Seafloor negative?
	   do nothing;
         Seafloor positive:
	    Compare seafloor/bottom_de:  
	      Seafloor deeper:
	        set bottom_de to seafloor
	      bottom_de deeper:
	        do nothing 
     ELSE bottom_de is undefined:
         Seafloor negative?  set bottom_de to mask value.  
	 Seafloor positive:  set bottom_de to seafloor-10.
  */

     
  new_bd = 0; 
      
  if (!is_flagged(*bd_addr,HB_f_empty) && !is_flagged(*bd_addr, HB_f_mask)) {  
        
	/* bottom_depth is defined */
     
     if (seafloor > 0){           /* negative seafloor means land */
          if ( ((float)seafloor - *bd_addr) > 100) {
             *bd_addr = seafloor - 10;
	     new_bd = 1;
	  }
	  
     }
  }  
  else {      /* bottom_de not defined */
     
    if (seafloor < 1) 
      *bd_addr = HB_f_mask;
    else 
      *bd_addr =  seafloor;

  }
  return(new_bd);       
  
}  /* end define_bottom() */	

/************************************************************************/
void get_sig_profiles(float lat, float lon, float *bdepth_ptr, short seafloor, int nfiles, char **argv, char *dir, char * ext, double *ppro, double *tpro, double *spro, short *pcount, short *tcount, short *scount, double *s0pro, double *s1pro, double *s2pro, double *s3pro, double *s4pro)

/*  Search cdf files for lon/lat, read in pr, te, sa, and compute all sigmas.  
    Read in bottom depth into the pointer and compare to seafloor to define
    a realistic bottom depth.  Mark masked levels with HB_f_mask.  
    Missing values are flagged with HB_f_empty.  Returns 1 if a new bottom depth 
    was defined or 0 if bottom_depth remains same at this square.

*/
{
   struct CDF_HDR cdf;
   int cdfid, curfile, print_msg = 0;
   int lon0to360, new_bd;
   float deepest, bd;
   float *xin;
   short *nobs;
   int error, i,  row, col, t90_avail;
   int tbin = 0, lev;
   int maxlev, bindx;

  /* initialize these */
  
   lon0to360 = 1;
   if (lon < 0)
      lon0to360 = 0;

   bindx = NSTDLEVS - 1;
   bd = HB_f_empty;
   deepest = HB_f_empty;           
   
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
           if (cdf.xmin > 0)  {
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
       
   new_bd = define_bottom(&bd, seafloor);

   if (!is_flagged(bd, HB_f_mask) && (bd > deepest)) { 
      deepest = bd;        /* this profile is the deepest */   
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
   
     if (error) {
          fprintf(stderr, "\n No te or t90 variable in cdf file:  %s\n", argv[curfile]);
          exit(1);
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
     
     if (new_bd || (bd < (deepest - 50.)) ) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
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

     if (new_bd || (bd < (deepest - 50.)) ) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
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
     
     if (new_bd || (bd < (deepest - 50.)) ) {        /* cancel bottom observation */
       xin[bindx] = cdf.fill_value;
       nobs[bindx] = 0;
     }
     
     for (lev = 0; lev < NSTDLEVS; ++lev) {
     
       if (is_flagged(xin[lev], cdf.fill_value) || is_flagged(xin[lev], HB_f_mask) ) {
          xin[lev] = 0.0;
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
     
  /* gridpoint is on land, mask all arrays */
  
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
   
   
  /* Not on land ...compute sigmas */
  
   compute_sigma(0.0, cdf.nz, s0pro, ppro, tpro, spro);
   compute_sigma(1000.0, cdf.nz, s1pro, ppro, tpro, spro);
   compute_sigma(2000.0, cdf.nz, s2pro, ppro, tpro, spro);
   compute_sigma(3000.0, cdf.nz, s3pro, ppro, tpro, spro);
   compute_sigma(4000.0, cdf.nz, s4pro, ppro, tpro, spro);

/* replace missing values with empty flag */

   for (lev = 0; lev <= bindx; ++lev) {
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

    if (deepest > 0.0)
         *bdepth_ptr = deepest;   /* this value gets returned:  a true bottom depth  */
    else 
         *bdepth_ptr = bd;
	       

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

/*******************************************/
double find_isopycnal( int zindex, int xrad, int yrad, int ncols, double **sigptr)
/*
   Search in each direction (WNES) to find profile with data at this
   level.  Suspend search in that direction if run into a masked value
   or if xrad/yrad is reached.  Return average of four sigmas weighted 
   by distance from central gridnode.  Return HB_f_empty if no value could be
   computed.
*/

{
   int row, col, stop, sq, xradius, yradius, wflag;
   double sig, ssum, wght, wsum;

   
   
   ssum = 0.0;
   wsum = 0.0;
   wflag = 0;

   /* search west */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      --col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   
   /* search east */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      ++col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   
   /* search north */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  += wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
  
   /* search south */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   /* search northwest */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      --col;
      ++row;
      wght = wght * 0.5 ;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
    /* search northeast */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      ++row;
      wght = wght * 0.5 ;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   /* search southeast */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      --row;
      wght = wght * 0.5 ;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   /* search southwest */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   yradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius <= xrad) && (++yradius <= yrad) && !stop)  {
      --col;
      --row;
      wght = wght * 0.5 ;
      sq = row * ncols + col;
      sig = sigptr[sq][zindex];
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
   if (!wflag) 
      return((double) HB_f_empty);
      
   return (ssum / wsum);
   
}  /* end find_isopycnal() */

/****************************************************/
double find_bottom_isopycnal(int xrad, int yrad, int ncols, double **sigptr, float *bdwork, float bdepth)
/*
   Search in each direction (WNES) to find profile with data near this
   depth to determine an appriate bottom isopycnal value.  Suspend 
   search in that direction if run into a masked value
   or if xrad/yrad is reached.  Return average of four sigmas weighted 
   by distance from central gridnode.  Return HB_f_empty if no value could be
   computed.
*/

{
   int i, row, col, stop, sq, xradius, yradius, zindex, wflag;
   double sig, ssum, wght, wsum;
   float bd;

   
   
   ssum = 0.0;
   wsum = 0.0;
   wflag = 0;
   zindex = NSTDLEVS-1;

   /* search west */
   
   row = yrad;  /* by definition, row/col corresponding to gridnode being fit */
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      --col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 151.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {  /* either a sig value or a mask */
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
   
   /* search east */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++xradius <= xrad && !stop)  {
      ++col;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
   /* search north */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 201.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  += wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
  
   /* search south */
   
   row = yrad;
   col = xrad;
   yradius = 0;
   wght = 2.0; 
   stop = 0; 
   while (++yradius <= yrad && !stop)  {
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */

      if (ABS(bd - bdepth) > 201.) {
         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }

   /* search northeast */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
    /* search southeast */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      ++col;
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
   /* search southwest */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      --col;
      --row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
    /* search northwest */
   row = yrad;
   col = xrad;
   xradius = 0;
   wght = 1.5; 
   stop = 0; 
   while ((++xradius < xrad) && (++yradius < yrad) && !stop)  {
      --col;
      ++row;
      wght = wght * 0.5;
      sq = row * ncols + col;
      bd = bdwork[sq];
      
      sig = sigptr[sq][zindex];  /* first assume neighbor is close to same bottom depth */
      
      if (ABS(bd - bdepth) > 201.) {

         if (bd < bdepth) 
	    sig = (double) HB_f_mask;  /* neighbor is too shallow */
	    
	 else {
	    /* find 2 stddepths that bracket bdepth and interpolate sig value */
	    i = 0;
	    while (i < zindex && (bdepth >= std_depth[i]))
	      ++i;
	      
	    if (i == 0)
	      	sig = (double) HB_f_empty;  /* bdepth < first std_depth  */
		
            else {   /* i and i-1 bracket the value ... */

                if (sigptr[sq][i-1] < 0.0 && sigptr[sq][i] < 0.0) {
		   sig = (double) HB_f_empty;  /* neighbor has no value  */
		}
		else { 
		  if (sigptr[sq][i-1] < 0.0)  /* upper level empty */ 
		     sig = sigptr[sq][i];
		  else if (sigptr[sq][i] < 0.0 )  /* lower level empty */
		     sig = sigptr[sq][i-1];
		  else 	                          /* interpolate */
	             sig = sigptr[sq][i-1] + (sigptr[sq][i] - sigptr[sq][i-1]) * (bdepth - std_depth[i-1]) / (std_depth[i] - std_depth[i-1]);
		  
		}
	    }
	 }
      }
      
      if (sig > 0 ) {
        if (sig < 100.) {
            ssum += sig * wght;
	    wsum  +=  wght;
	    ++wflag;
	}
        stop = 1;
      }
   }
  
 
   if (!wflag) 
      return((double) HB_f_empty);
      
   return (ssum / wsum);
   
}  /* end find_bottom_isopycnal() */

/****************************************************/


void get_neighbors(int pindex, double isopyc, double depth_being_interpolated, struct POINT *old_data, int *nalloc_addr, int *noldpts_addr, double *xsurf, float *bdpth, int nsq, struct GRID_INFO *hptr, double **xwork2, double **sigwork2)

/*
   Traverse each gridnode in work arrays. The bdpth array has the same
   dimension as the work grids and x grid -- it contains the seafloor
   depth at each gridnode. Interpolate to find prop value 
   associated with specified isopycnal.  Add a record to old_data array.
   If isopycnal runs into bottom or outcrops, mark xsurf[sq] with HB_zgrid_mask.
   Adjust old_data array size -- update nalloc, noldpts.
*/

{
   double *d, *x, *sig, xval;
   double diff, lastdiff, zval, best_zval;
   double *sigptr, *xptr;
   double emptyflag, lat, lon;
   int n, nz, lev, sq, row, col;
   int no_bottom, j, datagap, sqmid;


/* Allocate memory to hold profiles */

   d = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   x = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   sig = (double *) get_memory((void *)NULL, (size_t) NSTDLEVS, sizeof(double));
   
   emptyflag = HB_f_empty + 10.0;
   nz = NSTDLEVS - 1;
   sqmid = nsq / 2;    /* by definition, the square to be fit -- just label it empty */
   if (depth_being_interpolated < 0) {
       fprintf(stderr,"FATAL ERROR:  depth_being_interpolated is negative in get_neighbors()\n");  
    } 
   for (sq = 0; sq < nsq; ++sq) {
      if (sq == sqmid)
          continue;
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
      
      if (n == 0) {   /* no values at this grid node */
      
         if (is_flagged(bdpth[sq], HB_f_mask))  {  /* node is masked */
	    xsurf[sq] = (double) HB_zgrid_mask;
	 }
      }
      
      else  {   /* check for all crossings and use closest to this stddepth */
           j = 0;
	   lastdiff = (double)  HB_f_mask;  /* initialize with very large number */ 
	   best_zval =  (double) HB_f_empty;  /* initialize with very negative number */ 
	   do {
               zval = hb_linterp(isopyc, &sig[j], &d[j], n-j);
	       if (zval > -9990.) {
	        
	           if ( (diff = ABS(zval - depth_being_interpolated)) < lastdiff ) {
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
	           within 0.01 sigma units of the isopycnal value being sought */
		       
		       if (j == n)  /* should never happen, but check anyway */
		             --j;
			 
		       best_zval = d[j];
		  }/* end if diff > 0 */
		  
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
		
		
	           if (!datagap) {
	               row = sq / hptr->nx;
	               col = sq % hptr->nx;
	               lat = (double) (hptr->y_min + row * hptr->y_inc);
	               lon = (double) (hptr->x_min + col * hptr->x_inc);
	               old_data[*noldpts_addr].XP = lon;
	               old_data[*noldpts_addr].YP = lat;
	               old_data[*noldpts_addr].ZP = xval;
	               if (++(*noldpts_addr) == *nalloc_addr) {
	                 *nalloc_addr += 100;
	                  old_data = (struct POINT *) get_memory((void *)old_data, (size_t) *nalloc_addr, sizeof(struct POINT));
	              }
	           }
		
	  }
	 else {  /* no crossing -- does it outcrop?*/
	    
	       if (!no_bottom && (isopyc > sig[n-1])) {
	          xsurf[sq] = HB_zgrid_mask;   /* isopycnal runs into sea floor */
	       }
	 }
      }
      
   }  /* end for sq */

   free(d);
   free(x);
   free(sig);

}  /* end get_neighbors() */
