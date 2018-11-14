/* hb_build_oi_parms.c
................................................................................
                          *******  HydroBase 3 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             April 2011
...................................................
*  Creates or updates a HydroBase3 parameters file used in constructing
*  climatology fields in hb_ncoi3d.
*       semi-variogram parameters (nugget,sill-nugget,range) for P,T,S 
*              at each  gridpoint, stddepth
*              and for each month in seasonal layer (0-160 meters)
*        local bathymetry gradient (dx,dy)
*        local temperature gradient (dx,dy)
*        local pressure gradient (dx,dy)
*
*  hb_build_oi_parms -Ooutput_file [-B<w/e/s/n>] [-I<xinc/yinc>] 
*     -F<fg_file_root> -G<bin3d_files_root> [-T<topofile>]  [-U] [-Zstddepth_file]
*
* 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_grids.h"
#include "hb_memory.h"
#include "hb_oi3d.h"
#include "hb_paths.h"
#include "hb_variograms.h"
#include "netcdf.h"

#define SEASONAL_LAYER_DEF 10
#define NC_EXTENT ".nc"
#define MAXDIST 500
#define MIN_DECORR_LENGTH 100 /* minimum length scale for grid search (approx > deformation radius) */
#define MAX_DECORR_LENGTH 500 /* */
#define TOPOBUF 700    /* extra distance around input region to extract topography values */
#define MODEL_CODE 1   /* exponential model */
#define BINSIZE  20  /* for lags, in km */
#define MINOBS 30   /* min obs required within each lag bin for valid covariance*/
#define NBINTHRESH  5  /* # of bins allowed to have less than MINOBS */
#define MAXLAG2CHK   300      /* max lag to apply minobs/nbinthresh criteria */

struct SN_RATIOS {
   int ndepths, nlats;
   float *depths, *latbands;
   float *T, *S, *P;
};

/* global variables */

int lon0to360, xgreenwich;
float testmask, testempty;
int print_msg = 1;
int nfilled, nempty, nmasked;

int *merid;  /* used by function is_in_range() */

char cmonths[NMONTHS][4] = { "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec","clim" };

/* prototypes for locally defined functions */	

void print_usage(char *);
void create_ncparms(struct NC_INFO *, struct NC_PARMS *, struct GRID_INFO *, double *);
void mask_land_nodes(int, struct NC_INFO *, struct NC_PARMS *, struct GRID_INFO *, short *);
int check_overlap(struct GRID_INFO *, struct GRID_INFO *);
void set_topo_info(struct GRID_INFO *, struct GRID_INFO *, int);   
int get_size(struct GRID_INFO *, struct GRID_INFO *, double);
int get_bounds(struct GRID_INFO *, int, double, double, double *, double *, double *, double *);
int set_work_grid(struct GRID_INFO *, struct GRID_INFO *,double *, double *, double *, double *, int);
void get_profiles(int, double *, double *, double *, double *, double *, double *,  double *, double *, double *);
int get_depths(FILE *, double **, double **, double **);
int get_sn_ratios(int,struct SN_RATIOS *);
float snratio(double,float,struct SN_RATIOS *, char);
void set_default_lengthscales(double **, double **, struct NC_PARMS *);
int *set_mask2d(double *, struct GRID_INFO *, double, double);
int *set_mask2d_dist( struct GRID_INFO *, double, double, double, double);
int getsurf2d(double, int, double **, double **,  double **, struct GRID_INFO *, int *, double *, double *);
float get_mean(struct GRID_INFO *, double *, int, float, float);
void mask_deeper_levs(int, struct NC_INFO *, struct NC_PARMS *);
void grid_search(int, struct VARIOGRAM *, double *, double *, float *, int);
void get_gradients(double *, struct GRID_INFO *, double, double, float, int, float *, float *);
void fill_missing_nodes(int, struct NC_INFO *, struct NC_PARMS *, struct GRID_INFO *, double);
double fill_it(int, double, double, double *, struct GRID_INFO *, struct GRID_INFO *, double, double, double, double, double);

main(int argc, char **argv)

{
   int ncid, sn_ncid, topofile, ncid_fg[NMONTHS], ncid_b3d[NMONTHS];
   int error, i, j, k, update, fill_in, newparms, nz, ngrids, iclim, ilev;
   int iz4000, ideepest, nbinthresh, minobs, maxlag2chk;
   int row, col, ngood, zflag, xpole;
   int strow_fg, erow_fg, stcol_fg, ecol_fg;
   int strow_b3d, erow_b3d, stcol_b3d, ecol_b3d;
   int sq, nsq, nsq_fg, nsq_b3d, isq, isq_fg;
   int theSq_fg, theSq_b3d;    /* gridsquare corresponding to theLon, theLat */
   int nwsq_fg, nwsq_b3d;       /*actual # of nodes in work grids */
   int bflag, oflag, iflag, plotflag, snflag;
   int inrange, imonth, imodel, mask_it;
   int tindex_fg, tindex_b3d, reflev;
   int iwindow;   /* for computing horizontal property gradients */
   short *seafloor, topoval;
   char *toponame, *st, *outfile_name;
   char *fg_dir, *fg_ext, *fg_root;
   char *b3d_dir, *b3d_ext, *b3d_root;
   char *buf;
   float *pin, *tin, *sin;
   float *lats_b3d, *lons_b3d;
   float params[3];
   float parm0, parm1, parm2;
   float flon, flat;
   double maxdist, *xdist, *ydist, *stddepth;
   double lat, lon, latmin, latmax, lonmin, lonmax;
   double theLat, theLon;
   double theFG_t, theFG_s, theFG_p;
   double maxlag, binsize, surfval, dcorr_min, dcorr_max;
   struct NC_INFO ncparms_info;
   struct NC_PARMS ncparms_data;
   struct CDF_HDR hdr_fg, hdr_b3d;
   struct GRID_INFO fg_info, b3d_info, grid_info;
   struct GRID_INFO fg_orig, b3d_orig;
   struct GRID_INFO topo_info;
   struct SN_RATIOS sn_data;
   double *topolat, *topolon;
   double **sigptr_fg, **sigptr_b3d;
   double **p_fg, **t_fg,**th_fg, **s_fg, **sig0_fg, **sig1_fg, **sig2_fg;
   double **sig3_fg, **sig4_fg;
   double **p_b3d, **t_b3d,**th_b3d, **s_b3d, **sig0_b3d, **sig1_b3d;
   double **sig2_b3d, **sig3_b3d, **sig4_b3d;
   double *psurf_fg, *tsurf_fg, *ssurf_fg;
   double *psurf_b3d, *tsurf_b3d, *ssurf_b3d;
   double ***dtsurf, ***dssurf, ***dpsurf;
   double *bpress_fg, *bpress_b3d;
   double range1[3], range2[3];  /* for get_parm_range() */
   int *nz_b3d, *nz_fg;  /* nzlevs in each location in work grids */
   int *surfmask, *surfmask2;
   FILE *z_file;
   FILE *plot1, *plot2, *plot3, *plot4, *plot5, *plot6, *plot7;
   struct VARIOGRAM *vgram_ptr;
/*----------------------------------------*/  
 
 /* set these default values */
    fg_dir = "";
    fg_ext = NC_EXTENT;
    fg_root = NULL;
    b3d_dir = "";
    b3d_ext = NC_EXTENT;
    b3d_root = NULL;
    toponame = BATHPATH_C;
    bflag = zflag = 0;
    snflag = 0;
    update = 0;
    newparms = 0;
    fill_in = 0;
    plotflag = 0;
    merid = NULL;   /* used by function is_in_range() */
    z_file = NULL;
    stddepth = NULL;
    xdist = ydist = NULL;
    grid_info.x_inc = 1.0;  /* info about input region */
    grid_info.y_inc = 1.0;
    grid_info.x_min = 0.0;
    grid_info.x_max = 360.0;
    grid_info.y_min = -90.0;
    grid_info.y_max = 90.0;
    grid_info.node_offset = 1;  /* pixel grid: nodes at center of grid cells */
    grid_info.lon0to360 = 1;   /* lons all positive */
    grid_info.xgreenwich = 0;  /* lons not mixed sign */
    
    lon0to360 = 1;     /* by definition for nc_parms file */
    xgreenwich = 0;
    nfilled = nempty = nmasked = 0;
    testmask = (float) NCPARMS_MASK * 0.01;
    testempty = (float) NCPARMS_FILL * 0.01;
    binsize = BINSIZE;
    imodel = MODEL_CODE;
    iclim = NMONTHS-1;   /* index for clim files */
    iwindow = 3;  
    maxdist = MAXDIST;
    dcorr_min = (double) MIN_DECORR_LENGTH;
    dcorr_max = (double) MAX_DECORR_LENGTH;
    nbinthresh = NBINTHRESH;
    minobs = MINOBS;
    ncparms_zero(&ncparms_info, &ncparms_data);
    
/*----------------------------------------*/  
 /* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
/*----------------------------------------*/  
/* parse command line arguments */

   for (i = 1; i < argc; i++) {
   
      error = argv[i][0] == '-'? 0 : 1;
      if (!error ) {
        switch (argv[i][1]) {
      
	  
         case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &grid_info.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &grid_info.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &grid_info.y_min) != 1);
          while (*(st++) != '/')
                           ;
          error += (sscanf(st,"%lf", &grid_info.y_max) != 1);
                        	     
	  if (&grid_info.x_min > &grid_info.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (&grid_info.y_min > &grid_info.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
	  
	  grid_info.lon0to360 = grid_info.x_min >= 0 ? 1: 0;
	  grid_info.xgreenwich = (grid_info.x_min < 0) && (grid_info.x_max >=0) ? 1:0;
          break;

         case 'C':
             newparms = TRUE;
          break;
     
          case 'D':   /* Maximum search radius */
             error = (sscanf(&argv[i][2],"%lf", &maxdist) == 1) ? 0 : 1;
          break;
	  
          case 'F':                   
             switch (argv[i][2])  {
	        case 'd':
		   fg_dir = &argv[i][3];
		   break;
		case 'e':
		   fg_ext = &argv[i][3];
		   break;
		case 'r':
		   fg_root = &argv[i][3];
		   break;
		default:
		    error = TRUE;
	            fprintf(stderr,"\nUnrecognized option: %3s", argv[i]);
                    fprintf(stderr,"      FG filenames have form: <dir>/<root>.<month>.<extent> \n");
                    fprintf(stderr,"      Default extent is '.nc' \n");
                    fprintf(stderr,"      Use -Fr to specify root or dir/root (required)\n");
                    fprintf(stderr,"      Use -Fd to specify directory (optional)\n");
                    fprintf(stderr,"      Use -Fe to specify file extent (optional)\n");
		    exit(1);
	     }
	     
          break;
  
         case 'G':       
             switch (argv[i][2])  {
	        case 'd':
		   b3d_dir = &argv[i][3];
		   break;
		case 'e':
		   b3d_ext = &argv[i][3];
		   break;
		case 'r':
		   b3d_root = &argv[i][3];
		   break;
		default:
		    error = TRUE;
	            fprintf(stderr,"\nUnrecognized option: %3s", argv[i]);
                    fprintf(stderr," bin3d filenames have form: <dir>/<root>.<month>.<extent> \n");
                    fprintf(stderr,"Default extent is '.nc' \n");
                    fprintf(stderr,"Use -Gr to specify root or dir/root (required)\n");
                    fprintf(stderr,"Use -Gd to specify directory (optional)\n");
                    fprintf(stderr,"Use -Ge to specify file extent (optional)\n");
		    exit(1);
	     }
	  break;
	  
          case 'H':   /* binsize for computing lags */
             error = (sscanf(&argv[i][2],"%lf", &binsize) == 1) ? 0 : 1;
          break;
	  
          case 'I':
             iflag = 1;
             error = (sscanf(&argv[i][2],"%lf", &grid_info.x_inc) == 1) ? 0 : 1;
             grid_info.y_inc = grid_info.x_inc;
             st = &argv[i][2];
             while (*(st++) != '\0') {
               if (*st == '/') {
                  ++st;
                error += (sscanf(st,"%lf", &grid_info.y_inc) == 1) ? 0 : 1;
                break;
               }
             }
          break;
          case 'L':   /* range of decorrelation scales to consider */
             error = (sscanf(&argv[i][2],"%lf", &dcorr_min) == 1) ? 0 : 1;
	     if ((st = strchr(&argv[i][2],'/')) != NULL) {
	       ++st;
               error = (sscanf(st,"%lf", &dcorr_max) == 1) ? 0 : 1;
	     }
          break;
	  
          case 'M':   /* model code: 1 is exponential, 2 gaussian, 3 power */
             error = (sscanf(&argv[i][2],"%d", &imodel) == 1) ? 0 : 1;
          break;
	  
         case 'O':
	  oflag = 1;
          outfile_name = &argv[i][2];
          break;
         case 'P':
	  plotflag = 1;
	  buf = (char *)calloc(1000,sizeof(char));
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".Tparms");
	  plot1 = fopen(buf,"w");
	  if (plot1 == NULL) {
	     fprintf(stderr,"Unable to open plotting files \n");
	     fprintf(stderr,"Does the directory exist? \n");
	     exit(1);
	  }
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".Sparms");
	  plot2 = fopen(buf,"w");
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".Pparms");
	  plot3 = fopen(buf,"w");
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".Tgrad");
	  plot4 = fopen(buf,"w");
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".VgramObs");
	  plot5 = fopen(buf,"w");
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".Bathgrad");
	  plot6 = fopen(buf,"w");
          strcpy(buf, &argv[i][2]);
	  strcat(buf, ".PgramObs");
	  plot7 = fopen(buf,"w");
	  fprintf(stderr,"Opened 6 plotting files \n");
          break;
         case 'T':
          toponame = &argv[i][2];
          break;
         case 'U':
          update = TRUE;
          if (argv[i][2] != '\0') {
	     switch (argv[i][2])  {
                case 'I':
		case 'i' :
		   fill_in = TRUE;
                   update = FALSE;
	           break;
		case 'N':
		case 'n':
 	          snflag = 1;
                 error = nc_open(&argv[i][3],NC_NOWRITE,&sn_ncid);
	         if (error != NC_NOERR) {
	           fprintf(stderr,"\nError opening %s\n",&argv[i][3]);
		   exit(1);
	         }
		break;
		default:
		   error = 1;
	     } /* end switch */
	  }  /* end if */
          break;
          case 'Z':
            z_file = fopen(&argv[i][2],"r");
            if (z_file == NULL) {
                fprintf(stderr,"\nError opening depth/lengthscales file: %s\n",&argv[i][2]);
                exit(1);
            }
	    zflag = 1;
          break;
	 case 'h':  
	   print_usage(argv[0]);
	   exit(0);
     
         default:
          error = TRUE;
 
        } /* end switch */
      } /* end if*/
       
      if (error ) {
         fprintf(stderr,"\nError parsing command line args.\n");
         fprintf(stderr,"     in particular: '%s'\n", argv[i]);
         exit(1);
      }
   } /* end for */    
   
   
/*--------------------------------------------*/    
/*  Check syntax of options */ 

    if (!oflag ) {
       fprintf(stderr,"\nYou must specify -O<output_file>\n ");
       exit(1);
    }

    if ( !fill_in && !newparms && !update) {
       fprintf(stderr,"\nYou must specify one of these options:  -C, -U, -Ui or -Un \n ");
       exit(1);
    }    
/*--------------------------------------------*/ 
       if ((z_file == NULL) && newparms) {
          fprintf(stderr,"\nYou must specify -Z<filename> :  a list of depths and x-,y-lengthscales (in km) for each level that parameters are to be computed. \n ");
          exit(1);
       }
/*--------------------------------------------*/    
/* Get list of depths and lengthscales to compute variograms  */

    if (z_file != NULL) {
       ncparms_info.nz = get_depths(z_file, &stddepth, &xdist, &ydist);
       fclose(z_file);
       i = ncparms_info.nz - 1;
       if (stddepth[i] < 0)
          --ncparms_info.nz;  /* remove non-standard bottom depth*/
	  
    
       while ( stddepth[i] < 4001 && i < ncparms_info.nz) 
          ++i;
    
       iz4000 = i;
    }

/*--------------------------------------------*/    
/* Create new parameter file  */    

    if (newparms ) {
    
       create_ncparms(&ncparms_info, &ncparms_data, &grid_info, stddepth);
       free(stddepth);
       ncparms_info.modelCode = imodel;
       ncid = ncparms_init(outfile_name);
       ncparms_define(ncid, &ncparms_info, ncparms_data.depths);
       
       fprintf(stderr,"Masking land areas...  \n");
       topo_info.x_min = grid_info.x_min;
       topo_info.x_max = grid_info.x_max;
       topo_info.y_min = grid_info.y_min;
       topo_info.y_max = grid_info.y_max;
       seafloor = hb_get_topo(toponame, &topo_info, &topolat, &topolon, FALSE, topo_info.lon0to360, &topoval);
       mask_land_nodes(ncid, &ncparms_info, &ncparms_data, &topo_info, seafloor);
       
       ncparms_close(ncid);
       
       fprintf(stderr,"Successfully created %s \n", outfile_name);
       fprintf(stderr," with bounds: %.1lf %.1lf %.1lf %1lf\n", grid_info.x_min, grid_info.x_max, grid_info.y_min, grid_info.y_max);
        fprintf(stderr,"First node centered at lat: %.1f  lon: %.1f \n", ncparms_info.lats[0], ncparms_info.lons[0]);
	
        exit(0);
    }
    
/*--------------------------------------------*/ 
 
    if (update || fill_in ) {

      ncid = ncparms_update(outfile_name);
      if (ncid < 0) 
          exit(1);  


       /*  initialize all the fields in ncparms structures */  
       error = ncparms_getinfo(ncid, &ncparms_info, &ncparms_data, 1);
       if (error) {
          fprintf(stderr,"FATAL ERROR:  cannot parse %s correctly\n", outfile_name);
	  exit(1);
       }
       
       error = 0;
       if (zflag) {
          for (i=0; i < ncparms_info.nz; ++i) {
             if (error = (NINT(stddepth[i]) != NINT(ncparms_data.depths[i])) ) {
	       fprintf(stderr,"Mismatch between depths in parameters file and %s\n", outfile_name);
	       fprintf(stderr,"%.1lf vs. %.1f \n", stddepth[i], ncparms_data.depths[i]);
	       exit(1);
	     }
          }
	  
          free(stddepth);
       }
       if (!zflag) {
          while ( ncparms_data.depths[i] < 4001 && i < ncparms_info.nz) 
             ++i;
    
          iz4000 = i;
       }
       
       if (xdist == NULL ) 
          set_default_lengthscales(&xdist, &ydist, &ncparms_data);
       
       
       /*compare bounds of work region to ncparms file bounds */
       inrange = 1;
       if (bflag) {
           i = ncparms_info.nx - 1;
           inrange = is_in_range((float) grid_info.x_min, (ncparms_info.lons[0]-0.5*ncparms_info.xincr), (ncparms_info.lons[i] + 0.5 *ncparms_info.xincr ), merid, lon0to360) && is_in_range((float)grid_info.x_max, (ncparms_info.lons[0]-0.5*ncparms_info.xincr), (ncparms_info.lons[i]+ 0.5 *ncparms_info.xincr ), merid, lon0to360);
	   
	   if (inrange) {
	       i = ncparms_info.ny -1;
	       inrange = ((float) grid_info.y_min >= ncparms_info.lats[0]-0.5*ncparms_info.yincr) && ((float) grid_info.y_max <= ncparms_info.lats[i]+0.5*ncparms_info.yincr);
	    }
	    
	    free(merid);
	    merid = NULL;
	 grid_info.y_inc = ncparms_info.yincr;
	 grid_info.x_inc = ncparms_info.xincr;
	 grid_info.node_offset = 1;
	 grid_info.lon0to360 = grid_info.x_min >= 0 ? 1 : 0;
	 grid_info.xgreenwich = (grid_info.x_min < 0) && (grid_info.x_max >= 0);
	 grid_info.nx = NINT((grid_info.x_max - grid_info.x_min) / grid_info.x_inc);
	 grid_info.ny = NINT((grid_info.y_max - grid_info.y_min) / grid_info.y_inc);
	 if ( !grid_info.node_offset) {
	   ++grid_info.nx;
	   ++grid_info.ny;
	 }
       }
       else {   /* set  work region to bounds of params file */
         grid_info.x_min = ncparms_info.lons[0] - 0.5 * ncparms_info.xincr;
         grid_info.x_max = ncparms_info.lons[ncparms_info.nx-1] + 0.5 * ncparms_info.xincr;
         grid_info.y_min = ncparms_info.lats[0] - 0.5 * ncparms_info.yincr;
         grid_info.y_max = ncparms_info.lons[ncparms_info.ny-1] + 0.5 * ncparms_info.yincr;
	 grid_info.y_inc = ncparms_info.yincr;
	 grid_info.x_inc = ncparms_info.xincr;
	 grid_info.node_offset = 1;
	 grid_info.nx = ncparms_info.nx;
	 grid_info.ny = ncparms_info.ny;
	 grid_info.lon0to360 = lon0to360;
	 grid_info.xgreenwich = xgreenwich;
       }
       
       if (! inrange) {  
         fprintf(stderr,"ERROR:%s does not include bounds specified here. \n", outfile_name);
         fprintf(stderr,"File dimensions cannot be changed. \n");
         fprintf(stderr,"Recommend creating a parameters file with global boundaries\n");
         exit(1);
       }
       
      lon0to360 = ncparms_info.lons[0] >= 0 ? 1: 0;
      xgreenwich = (ncparms_info.lons[0] < 0) && (ncparms_info.lons[ncparms_info.nx-1] >=0) ? 1:0;
     
    } /* end if update || fill_in */
/*--------------------------------------------*/ 
/*--------------------------------------------*/    
  /* At this point, info about the work region (grid_info) 
     and parameter file (ncparms_info) is set */
      
    nsq = grid_info.nx * grid_info.ny;
/*--------------------------------------------*/ 
/*--------------------------------------------*/ 
  
  
    if (fill_in) {
    
       fill_missing_nodes(ncid, &ncparms_info, &ncparms_data, &grid_info, maxdist);
       ncparms_close(ncid);
       fprintf(stderr,"\nNumber of nodes empty: %d", nempty);
       fprintf(stderr,"\nNumber of nodes filled: %d", nfilled);
       fprintf(stderr,"\nEnd of %s\n", argv[0]);
       exit(0);
    }

/*--------------------------------------------*/ 
/*--------------------------------------------*/ 
/*   Update region with new parameters  */ 

    if  ((fg_root == NULL) || (b3d_root == NULL)) {
         fprintf(stderr,"\nYou must specify files for first guess fields [-Fr<file_root>] and for bind3d fields [-Gr<file_root>] \n");
         exit(1);
    }
   
  /*--------------------------------------------*/    
/* Get info about FG grids and check for overlap with work region */

    buf = (char *)calloc(200, sizeof(char));
    strcpy(buf, fg_root);
    strncat(buf, ".clim", 5);
    
    ncid_fg[iclim] = cdf_open(fg_dir, buf, fg_ext, print_msg);
    error = read_cdf_hdr(ncid_fg[iclim], &hdr_fg);
    
    fg_orig.x_min = hdr_fg.xmin;
    fg_orig.x_max = hdr_fg.xmax;
    fg_orig.y_min = hdr_fg.ymin;
    fg_orig.y_max = hdr_fg.ymax;
    fg_orig.y_inc = hdr_fg.yincr;
    fg_orig.x_inc = hdr_fg.xincr;
    fg_orig.nx = hdr_fg.nx;
    fg_orig.ny = hdr_fg.ny;
    fg_orig.node_offset = hdr_fg.node_offset;
    fg_orig.lon0to360 = fg_orig.x_min >= 0 ? 1 : 0;
    fg_orig.xgreenwich = (fg_orig.x_min < 0) && (fg_orig.x_max >= 0) ? 1 : 0;
    
    
    if (!check_overlap(&fg_orig, &grid_info)) {
        fprintf(stderr,"ERROR: %s does not overlap the work region:\n", buf);
        fprintf(stderr,"  %.3f/%.3f/%.3f/%.3f\n", grid_info.x_min, grid_info.x_max, grid_info.y_min, grid_info.y_max);
        exit(1);
    }
    
    /* determine which Temperature property is available */
    tindex_fg = (int)T90;
    error = nc_inq_varid( ncid_fg[iclim],"t90", &j);
    if (error != NC_NOERR) {
       tindex_fg = (int)TE;
       if (nc_inq_varid (ncid_fg[iclim], "te", &j) != NC_NOERR) {
        fprintf(stderr,"FATAL ERROR: could not find te or t90 in %s\n", buf);
        exit(1);   
       }
    }
    
    free(buf);
/*--------------------------------------------*/    
/* Get info about bin3d grids and check for overlap with work region */
    buf = (char *)calloc(200, sizeof(char));
    strcpy(buf, b3d_root);
    strncat(buf, ".clim", 5);

    ncid_b3d[iclim] = cdf_open(b3d_dir, buf, b3d_ext, print_msg);
    error = read_cdf_hdr(ncid_b3d[iclim], &hdr_b3d);
    
    b3d_orig.x_min = hdr_b3d.xmin;
    b3d_orig.x_max = hdr_b3d.xmax;
    b3d_orig.y_min = hdr_b3d.ymin;
    b3d_orig.y_max = hdr_b3d.ymax;
    b3d_orig.y_inc = hdr_b3d.yincr;
    b3d_orig.x_inc = hdr_b3d.xincr;
    b3d_orig.nx = hdr_b3d.nx;
    b3d_orig.ny = hdr_b3d.ny;
    b3d_orig.node_offset = hdr_b3d.node_offset;
    b3d_orig.lon0to360 = b3d_orig.x_min >= 0 ? 1 : 0;
    b3d_orig.xgreenwich = (b3d_orig.x_min < 0) && (b3d_orig.x_max >= 0) ? 1 : 0;

    if (!check_overlap(&b3d_orig, &grid_info)) {
        fprintf(stderr,"ERROR: %s does not overlap the work region:\n", buf);
        fprintf(stderr,"  %.3f/%.3f/%.3f/%.3f\n", grid_info.x_min, grid_info.x_max, grid_info.y_min, grid_info.y_max);
        exit(1);
    }
    
    /* determine which Temperature property is available */
    tindex_b3d = (int)T90;
    error = nc_inq_varid( ncid_b3d[iclim],"t90", &j);
    if (error != NC_NOERR) {
       tindex_b3d = (int)TE;
       if (nc_inq_varid (ncid_b3d[iclim], "te", &j) != NC_NOERR) {
        fprintf(stderr,"FATAL ERROR: could not find te or t90 in %s\n", buf);
        exit(1);    
       }
    }
    free(buf);
/*--------------------------------------------*/    
/* Read in topography values for output region + search radius */

    set_topo_info(&grid_info, &topo_info, TOPOBUF);
    seafloor = hb_get_topo(toponame, &topo_info, &topolat, &topolon, FALSE, topo_info.lon0to360, &topoval);
    
/*--------------------------------------------*/    
/* Read in S/N ratios for latitude bands and depths*/

   if (snflag)  {
     error = get_sn_ratios(sn_ncid, &sn_data);
     nc_close(sn_ncid);
     if (error) {
       fprintf(stderr,"\n FATAL ERROR: in get_sn_ratios()\n");
       exit(1);
     }
   }
/*--------------------------------------------*/    
/* Determine size of work area for FG fields,
   and allocate memory */
  
   nsq_fg = get_size(&grid_info, &fg_orig, maxdist);

   p_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   t_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   s_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   th_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   sig0_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   sig1_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   sig2_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   sig3_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
   sig4_fg = (double **)get_memory(NULL,nsq_fg,sizeof(double *));
    
   for (i = 0; i < nsq_fg; ++i) {
        p_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
        t_fg[i] = (double *)get_memory(NULL,  hdr_fg.nz, sizeof(double));
        s_fg[i] = (double *)get_memory(NULL,  hdr_fg.nz, sizeof(double));
        th_fg[i] = (double *)get_memory(NULL, hdr_fg.nz , sizeof(double));
        sig0_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
        sig1_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
        sig2_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
        sig3_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
        sig4_fg[i] = (double *)get_memory(NULL, hdr_fg.nz, sizeof(double));
    }
/*--------------------------------------------*/   
/* Determine size of work area (max search radius) for bin3d fields,
   set up structures, allocate memory */

   nsq_b3d = get_size(&grid_info, &b3d_orig, maxdist);
     
   p_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   t_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   s_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   th_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   sig0_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   sig1_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   sig2_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   sig3_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
   sig4_b3d = (double **)get_memory(NULL,nsq_b3d,sizeof(double *));
     
   for (i = 0; i < nsq_b3d; ++i) {
        p_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
        t_b3d[i] = (double *)get_memory(NULL,  hdr_b3d.nz, sizeof(double));
        s_b3d[i] = (double *)get_memory(NULL,  hdr_b3d.nz, sizeof(double));
        th_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz , sizeof(double));
        sig0_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
        sig1_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
        sig2_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
        sig3_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
        sig4_b3d[i] = (double *)get_memory(NULL, hdr_b3d.nz, sizeof(double));
    }

/*--------------------------------------------*/ 
  
/*  Allocate memory for surface grids, residuals, and other work space */

    psurf_fg = (double *)get_memory(NULL,nsq_fg, sizeof(double));
    tsurf_fg = (double *)get_memory(NULL,nsq_fg, sizeof(double));
    ssurf_fg = (double *)get_memory(NULL,nsq_fg, sizeof(double));
    bpress_fg = (double *) get_memory(NULL,nsq_fg,sizeof(double ));
    nz_fg = (int *)get_memory(NULL,nsq_fg,sizeof(int ));
    
    psurf_b3d = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
    tsurf_b3d = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
    ssurf_b3d = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
    bpress_b3d = (double *) get_memory(NULL,nsq_b3d,sizeof(double ));
    nz_b3d = (int *)get_memory(NULL,nsq_b3d,sizeof(int ));
    
    /* this dimension is maximum of nlevs in seasonal layer vs
       deepest layer */
    ngrids = ncparms_info.nz - iz4000;
    if ( ncparms_info.nz_seas > ngrids)
        ngrids = ncparms_info.nz_seas;
	
	
    dpsurf = (double ***) get_memory(NULL,ngrids, sizeof(double));
    dssurf = (double ***) get_memory(NULL,ngrids, sizeof(double));
    dtsurf = (double ***) get_memory(NULL,ngrids, sizeof(double));
    
     for (ilev = 0; ilev < ngrids; ++ilev) {
       dpsurf[ilev] = (double **)get_memory(NULL,NMONTHS, sizeof(double *));
       dssurf[ilev] = (double **)get_memory(NULL,NMONTHS, sizeof(double *));
       dtsurf[ilev] = (double **)get_memory(NULL,NMONTHS, sizeof(double *));
     
       for (imonth = 0; imonth < NMONTHS; ++imonth) {
          dpsurf[ilev][imonth] = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
          dtsurf[ilev][imonth] = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
          dssurf[ilev][imonth] = (double *)get_memory(NULL,nsq_b3d, sizeof(double));
       }
    }
        
/*--------------------------------------------*/   
/*  Open each monthly FG and Bin3d file.  *.clim.nc files are already open */

    for (imonth = 0; imonth < iclim; ++imonth) {
        buf = (char *)calloc(200, sizeof(char));
        strcpy(buf, fg_root);
        strncat(buf, ".", 1);
        strncat(buf, cmonths[imonth], 4);
        ncid_fg[imonth] = cdf_open(fg_dir, buf, fg_ext, print_msg);
        free(buf);
    }
 
     for (imonth = 0; imonth < iclim; ++imonth) {
        buf = (char *)calloc(200, sizeof(char));
        strcpy(buf, b3d_root);
        strncat(buf, ".", 1);
        strncat(buf, cmonths[imonth], 4);
        ncid_b3d[imonth] = cdf_open(b3d_dir, buf, b3d_ext, print_msg);
        free(buf);
    }
    
    fprintf(stderr,"\nCompleted setup phase.... \n");
       
/**********************************************   
      END OF SET UP  
***********************************************/    

/* Loop for each gridpoint in the specified output region: grid_info 
   For each depth level, fill in all the fields of ncparms_data.  The seasonal
   layer is done month-by-month.  When all layers have been visited,
   write out the ncparms_data to the open parameters file (ncid)  */

    for (sq = 0; sq < nsq; ++sq) {
       fprintf(stderr,"\n%2d/%2d ", sq+1, nsq); /* show progress */
       
       ncparms_prefill(&ncparms_data, ncparms_info.fill_value);
       
       sq2rc(sq, &grid_info, &row, &col);
       ij2xy(&grid_info, col, row, &theLon, &theLat);
       
       xpole = get_bounds(&grid_info, sq, maxdist, maxdist, &latmin, &lonmin, &latmax, &lonmax);
       nwsq_fg = set_work_grid(&fg_orig, &fg_info,&latmin, &lonmin, &latmax, &lonmax, xpole);
       
       if (nwsq_fg <= 0)
          continue;
	  
       if (nwsq_fg > nsq_fg) {
	   fprintf(stderr,"Insufficient memory allocated to hold first guess profiles.\nNeed to rethink get_size\nExiting...\n");
	   exit(1);
       }
       
       error = xy2ij(&fg_info, theLon, theLat, &col, &row);
       if (error)
          continue;
	  
	  
	/* Do a quick check that central grid square corresponds to theLat,theLon */
	
       if (fg_info.nx % 2 == 0) { 
          ++fg_info.nx;
          if (col % 2 ) {  /* add column on east */
	    fg_info.x_max +=  fg_info.x_inc;
	  }
	  else {  /* add on west */
	    fg_info.x_min -= fg_info.x_inc;
	  }
       }

       if (fg_info.ny % 2 == 0) { 
          ++fg_info.ny;
          if (col % 2 ) {  /* add column on north */
	    fg_info.y_max +=  fg_info.y_inc;
	  }
	  else {  /* add on south */
	    fg_info.y_min -= fg_info.y_inc;
	  }
	  
       }
       
       error = xy2ij(&fg_info, theLon, theLat, &col, &row);
       theSq_fg = row * fg_info.nx + col;
             
       error = get_indices(&hdr_fg,(float)latmin, (float)lonmin, &strow_fg, &stcol_fg);
       if (error) {
          if (strow_fg >= hdr_fg.ny)
	    strow_fg = hdr_fg.ny-1;
          if (strow_fg < 0)
	    strow_fg = 0;
	  if (stcol_fg < 0)
	     stcol_fg = 0;
       }
       
       error = get_indices(&hdr_fg,(float)latmax, (float)lonmax, &erow_fg, &ecol_fg);
       if (error) {
          if (erow_fg >= hdr_fg.ny)
	    erow_fg = hdr_fg.ny-1;
          if (erow_fg < 0)
	    erow_fg = 0;
	  if (ecol_fg >= hdr_fg.nx)
	     ecol_fg = hdr_fg.nx -1;
       }

       xpole = get_bounds(&grid_info, sq, maxdist, maxdist, &latmin, &lonmin, &latmax, &lonmax);
       nwsq_b3d = set_work_grid(&b3d_orig, &b3d_info, &latmin, &lonmin, &latmax, &lonmax, xpole);
       if (nwsq_b3d <= 0) 
            continue;
       
       if ( nwsq_b3d > nsq_b3d) {
          fprintf(stderr,"Insufficient memory allocated to hold bin3d profiles.\nNeed to rethink get_size()\nExiting...\n");
	  exit(1);
        }
       
       error = xy2ij(&b3d_info, theLon, theLat, &col, &row);
       if (error)
           continue;
	   
  /* Do a quick check that central grid square corresponds to theLat,theLon */
	
       if (b3d_info.nx % 2 == 0) { 
          ++b3d_info.nx;
          if (col % 2 ) {  /* add column on east */
	    b3d_info.x_max +=  b3d_info.x_inc;
	  }
	  else {  /* add on west */
	    b3d_info.x_min -= b3d_info.x_inc;
	  }
       }

       if (b3d_info.ny % 2 == 0) { 
          ++b3d_info.ny;
          if (col % 2 ) {  /* add column on north */
	    b3d_info.y_max +=  b3d_info.y_inc;
	  }
	  else {  /* add on south */
	    b3d_info.y_min -= b3d_info.y_inc;
	  }
	  
       }
	   
       error = xy2ij(&b3d_info, theLon, theLat, &col, &row);
       theSq_b3d = row * b3d_info.nx + col;
            
       /*  lat/lon vectors needed for vgram_estimate() */
       
       lats_b3d = (float *)get_memory(NULL,b3d_info.ny, sizeof(float));
       lons_b3d = (float *)get_memory(NULL,b3d_info.nx, sizeof(float));
  
       for (i = 0; i < b3d_info.ny; ++i)
           lats_b3d[i] = (float) (b3d_info.y_min + i * b3d_info.y_inc);

       for (i = 0; i < b3d_info.nx; ++i)
           lons_b3d[i] = (float) (b3d_info.x_min + i * b3d_info.x_inc);
	   
       
       error = get_indices(&hdr_b3d,(float)latmin, (float)lonmin, &strow_b3d, &stcol_b3d);
       if (error) {
          if (strow_b3d < 0)
	    strow_b3d = 0;
          if (strow_b3d >= hdr_b3d.ny)
	    strow_b3d = hdr_b3d.ny-1;
	  if (stcol_b3d < 0)
	     stcol_b3d = 0;
       }
       
       error = get_indices(&hdr_b3d,(float)latmax, (float)lonmax, &erow_b3d, &ecol_b3d);
       if (error) {
          if (erow_b3d >= hdr_b3d.ny)
	    erow_b3d = hdr_b3d.ny-1;
          if (erow_b3d < 0 )
	    erow_b3d = 0;
	    
	  if (ecol_b3d >= hdr_b3d.nx)
	     ecol_b3d = hdr_b3d.nx -1;
       }
       
       /* Loop for each month, including clim */
       
       for (imonth = 0; imonth <= iclim; ++imonth) {
       
	  fprintf(stderr,"+");   /* show progress */

	  /*prefill the work arrays */
	  
          nz = hdr_fg.nz;   
	  for (isq = 0; isq < nwsq_fg; ++isq) {
	     for (i = 0; i < nz; ++i) {
		      p_fg[isq][i] = HBEMPTY;
		      t_fg[isq][i] = HBEMPTY;
		      s_fg[isq][i] = HBEMPTY;
	      }
	      nz_fg[isq] = 0;
              bpress_fg[isq] = 0;
	  }
	  
         /* Read in FG profiles for this month   */       

          pin = (float *)calloc(nz, sizeof(float));
          tin = (float *)calloc(nz, sizeof(float));
          sin = (float *)calloc(nz, sizeof(float));
         
          for (row = strow_fg; row >= erow_fg; --row) {
             for (col = stcol_fg; col <= ecol_fg; ++col) {
	        
	        error = read_cdf_prop(ncid_fg[imonth],"pr", pin, row, col, 0, 0, nz);
	        error = read_cdf_prop(ncid_fg[imonth],get_prop_mne(tindex_fg), tin, row, col, 0, 0, nz);
	        error = read_cdf_prop(ncid_fg[imonth],"sa", sin, row, col, 0, 0, nz);

                /* copy from float to double arrays and eliminate missing values  */
		
		get_lat_lon(&hdr_fg, row, col, &flat, &flon);
		lat = (double) flat;
		lon = (double) flon;
		xy2ij(&fg_info, lon, lat, &i, &j);
		isq = j * fg_info.nx + i;
		   		
		j = 0;
		for (i = 0; i < nz; ++i) {
		   if (pin[i] > testempty && pin[i] < testmask) {
		    
		      p_fg[isq][j] = (double) pin[i];
		      t_fg[isq][j] = (double) tin[i];
		      s_fg[isq][j++] = (double) sin[i];
		       
		    }
		}
		nz_fg[isq] = j;
		
		/* get bottom pressure value */
		
		topoval = find_nearest_topo_val(lat, lon, seafloor, &topo_info); 
		bpress_fg[isq] = hb_p80((double)topoval, lat);
		
	        get_profiles(nz_fg[isq], p_fg[isq], t_fg[isq], s_fg[isq],th_fg[isq], sig0_fg[isq], sig1_fg[isq], sig2_fg[isq], sig3_fg[isq], sig4_fg[isq]);
		
	     } /* end for col */
          }  /* end for row */
	  
	  
	  free(pin);
	  free(tin);
	  free(sin);
	  
         /* Read in bin3d profiles for this month  */

	  /*prefill the work arrays */
          nz = hdr_b3d.nz;
 	  for (isq = 0; isq < nwsq_b3d; ++isq) {
	     for (i = 0; i < nz; ++i) {
		      p_b3d[isq][i] = HBEMPTY;
		      t_b3d[isq][i] = HBEMPTY;
		      s_b3d[isq][i] = HBEMPTY;
	      }
	      nz_b3d[isq] = 0;
              bpress_b3d[isq] = 0;
	  }
	  
          pin = (float *)calloc(nz, sizeof(float));
          tin = (float *)calloc(nz, sizeof(float));
          sin = (float *)calloc(nz, sizeof(float));
	  
          for (row = strow_b3d; row >= erow_b3d; --row) {
             for (col = stcol_b3d; col <= ecol_b3d; ++col) {
	        error = read_cdf_prop(ncid_b3d[imonth],"pr", pin, row, col, 0, 0, nz);
	        error = read_cdf_prop(ncid_b3d[imonth],get_prop_mne(tindex_b3d), tin, row, col, 0, 0, nz);
	        error = read_cdf_prop(ncid_b3d[imonth],"sa", sin, row, col, 0, 0, nz);
		
		get_lat_lon(&hdr_b3d, row, col, &flat, &flon);
		lat = (double) flat;
		lon = (double) flon;
		xy2ij(&b3d_info, lon, lat, &i, &j);
		isq = j * b3d_info.nx + i;
		
               /* remove all missing and masked values */

		j = 0;
 		for (i = 0; i < nz; ++i) {
		   if (pin[i] > testempty && pin[i] < testmask) {
		    
		       p_b3d[isq][j] = (double) pin[i];
		       t_b3d[isq][j] = (double) tin[i];
		       s_b3d[isq][j++] = (double) sin[i];
		   }
		}
		nz_b3d[isq] = j;
		
		/* get bottom pressure value */
		
		topoval = find_nearest_topo_val(lat, lon, seafloor, &topo_info); 
		bpress_b3d[isq] = hb_p80((double)topoval, lat);
		
		get_profiles(nz_b3d[isq], p_b3d[isq], t_b3d[isq], s_b3d[isq],th_b3d[isq], sig0_b3d[isq], sig1_b3d[isq], sig2_b3d[isq], sig3_b3d[isq], sig4_b3d[isq]);
				
	     } /* end for col */
          }  /* end for row */
	  
	  free(pin);
	  free(tin);
	  free(sin);

  
/********************************************************************/
   /* Foreach surface in the seasonal layer evaluate and store residuals  */
/********************************************************************/
	  	  
	  reflev = -1;      /* Use pressure surfaces in the seasonal layer */
	  
	  for (ilev = 0; ilev < ncparms_info.nz_seas; ++ilev) {
	  
	     surfval = ncparms_data.depths[ilev];
	     
	     ngood = getsurf2d(surfval, reflev, p_fg, th_fg, p_fg, &fg_info, nz_fg, bpress_fg,tsurf_fg );

	     mask_it = tsurf_fg[theSq_fg] > testmask ? 1 : 0;
	     if (mask_it)  {
	        for (i = ilev; i < ncparms_info.nz_seas; ++i) {
		    dtsurf[i][imonth][theSq_b3d] = (double)HBMASK;
		    dssurf[i][imonth][theSq_b3d] = (double)HBMASK;
		}
		ilev = ncparms_info.nz_seas;
		continue;
	     } /* end if mask_it */
	     
	     theFG_t = tsurf_fg[theSq_fg]; 
	        if (theFG_t < testempty) {
	          theFG_t = get_mean(&fg_info, tsurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }
	   
	     ngood = getsurf2d(surfval, reflev, p_fg, s_fg,  p_fg, &fg_info, nz_fg, bpress_fg,ssurf_fg);
	     surfmask = set_mask2d(tsurf_fg, &fg_info,theLat,theLon);
	     surfmask2 = set_mask2d_dist(&fg_info, xdist[ilev], ydist[ilev], theLat, theLon);
	     for (j = 0; j < nwsq_fg; ++j) {
	         if (surfmask[j] > 0 || surfmask2[j] > 0) {
		     tsurf_fg[j] = (double)HBMASK;
		     ssurf_fg[j] = (double)HBMASK;
		 }
	     } 
	     free(surfmask);
	     free(surfmask2);
	     
	     theFG_s = ssurf_fg[theSq_fg]; 
	        if (theFG_s < testempty) {
	          theFG_s = get_mean(&fg_info, ssurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }

	     ngood = getsurf2d(surfval, reflev, p_b3d, th_b3d, p_b3d, &b3d_info,nz_b3d, bpress_b3d, tsurf_b3d);
		 
	     ngood = getsurf2d(surfval, reflev, p_b3d, s_b3d, p_b3d, &b3d_info, nz_b3d, bpress_b3d, ssurf_b3d );
	     surfmask = set_mask2d(tsurf_b3d, &b3d_info,theLat, theLon);
	     surfmask2 = set_mask2d_dist(&b3d_info, xdist[ilev], ydist[ilev], theLat, theLon);
	     for (isq = 0; isq < nwsq_b3d; ++isq) {
	         if (surfmask[isq] > 0 || surfmask2[isq] > 0) {
		     tsurf_b3d[isq] = (double)HBMASK;
		     ssurf_b3d[isq] = (double)HBMASK;
		 }
	     } 
	     free(surfmask);
	     free(surfmask2);

      /* Compute residuals */ 
	     
	     for (isq = 0; isq < nwsq_b3d; ++isq) {  
	        error = sq2rc(isq, &b3d_info, &row, &col);
	        error += ij2xy(&b3d_info, col, row, &lon, &lat);
	        error += xy2ij(&fg_info, lon, lat, &col, &row);
	        if (error) {
		  dtsurf[ilev][imonth][isq] = (double)HBEMPTY;
		  dssurf[ilev][imonth][isq] = (double)HBEMPTY;
		}
	        else {
		   isq_fg = row * fg_info.nx + col;
		   if (tsurf_b3d[isq] > testmask || tsurf_fg[isq_fg] > testmask) {
		    dtsurf[ilev][imonth][isq] = (double)HBMASK;
		    dssurf[ilev][imonth][isq] = (double)HBMASK;
		   }
		   else if (tsurf_b3d[isq] < testempty  ) {
		    dtsurf[ilev][imonth][isq] = (double)HBEMPTY;
		    dssurf[ilev][imonth][isq] = (double)HBEMPTY;
		   }
		   else {
		    dtsurf[ilev][imonth][isq] = tsurf_b3d[isq] - theFG_t;
		    dssurf[ilev][imonth][isq] = ssurf_b3d[isq] - theFG_s;
		   } 
	        }		
	     } /* end for isq */
	     
	     
	     if (imonth == iclim) {
                 /* Compute meridional and zonal property gradients */
	     
	        get_gradients(tsurf_fg, &fg_info, theLat, theLon, ncparms_info.fill_value, iwindow, &parm1, &parm2);
	        ncparms_data.dTdx[ilev] = parm1;
                ncparms_data.dTdy[ilev] = parm2;
		
	       if (plotflag) 
		   fprintf(plot4, "%5.1lf %5.1lf %5.0f %.3f %.3f \n", theLon, theLat, ncparms_data.depths[ilev], parm1, parm2);
	    
	       ncparms_data.dPdx[ilev] = ncparms_info.fill_value;
	       ncparms_data.dPdy[ilev] = ncparms_info.fill_value;
	       ncparms_data.Lx[ilev] = xdist[ilev];
	       ncparms_data.Ly[ilev] = ydist[ilev];
             }
	     
	  } /* end for ilev in the seasonal layer*/
        } /* end for each month */

    /* After residuals have been computed for all months and levels in seasonal layer, 
       compute variograms for each level (all months) and write to output files */ 

/*******************************/
       VGRAM_EST:
           ;
/*******************************/
	   
       for (ilev = 0; ilev < ncparms_info.nz_seas; ++ilev) {
       
           mask_it = dtsurf[ilev][iclim][theSq_b3d] > testmask ? 1 : 0;
	   if (mask_it) {
	       mask_deeper_levs(ilev, &ncparms_info, &ncparms_data);
	       goto NEXTSQUARE;
	   }    
       

         /* Estimate Temperature variogram and execute a grid search to find parameters that minimize errors */
	     
	     maxlag = xdist[ilev] > ydist[ilev] ? xdist[ilev] : ydist[ilev];
	     if (dcorr_max > maxlag)
	        maxlag = dcorr_max;
		
	     vgram_ptr = vgram_estimate(dtsurf[ilev], 12, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
	     
	     /* check for valid variogram */
	     
	     maxlag2chk = MAXLAG2CHK > maxlag ? (int)maxlag : MAXLAG2CHK;
	     if (vgram_check(vgram_ptr, nbinthresh, minobs, maxlag2chk) == 0 ) {
	        ++nempty;
	        continue;
	     }
	     
	     error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	     
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'T');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
	     
	     ncparms_data.Tparm0[ilev] = params[0];
	     ncparms_data.Tparm1[ilev] = params[1];
	     ncparms_data.Tparm2[ilev] = params[2];
	     ncparms_data.Lx[ilev] = (float) xdist[ilev];
	     ncparms_data.Ly[ilev] = (float) ydist[ilev];
            ++nfilled;
	     
	     if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
                   fprintf(plot1, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2]);	
                  fprintf(plot5, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f %5d", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2], vgram_ptr->nlags);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->Vh[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->h[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %5d", vgram_ptr->Nh[k]);
		   fprintf(plot5,"\n");
	     }        
	     free(vgram_ptr);
	     
    /* Estimate Salinity variogram */	     

	     vgram_ptr = vgram_estimate(dssurf[ilev],12, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
	     error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;

	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'S');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
	     
	     if (params[1] > testempty && params[2] > testempty) {
	        ncparms_data.Sparm0[ilev] = params[0];
	        ncparms_data.Sparm1[ilev] = params[1];
	        ncparms_data.Sparm2[ilev] = params[2];
	     
	        if (plotflag  && (params[1] < testmask) && (params[1] > testempty)) {
	           fprintf(plot2, "%5.1lf %5.1lf %5.0f %8.5f %8.5f %3.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
	        }
	     }
	     free(vgram_ptr);
	     
	  } /* end for ilev in the seasonal layer*/
	
/*--------------------------------------------*/    
/*--------------------------------------------*/    
   /* Loop for levels below seasonal layer and above 4000 */

	imonth = iclim;
	
	for (ilev = ncparms_info.nz_seas; ilev < iz4000; ++ilev) {
 	     fprintf(stderr,"-");   /* show progress */
          
	   reflev = 4000;
	   sigptr_fg = sig4_fg;
	   sigptr_b3d = sig4_b3d;
	   if (ncparms_data.depths[ilev] <= 500) {
	      reflev = 0;
	      sigptr_fg = sig0_fg;
	      sigptr_b3d = sig0_b3d;
	   }
	   else if (ncparms_data.depths[ilev] <= 1500) {
	      reflev = 1000;
	      sigptr_fg = sig1_fg;
	      sigptr_b3d = sig1_b3d;
	   }
	   else if (ncparms_data.depths[ilev] <= 2500) {
	      reflev = 2000;
	      sigptr_fg = sig2_fg;
	      sigptr_b3d = sig2_b3d;
	   }
	   else if (ncparms_data.depths[ilev] <= 3500) {
	      reflev = 3000;
	      sigptr_fg = sig3_fg;
	      sigptr_b3d = sig3_b3d;
	   }
	   
	   surfval = hb_linterp((double) ncparms_data.depths[ilev], p_fg[theSq_fg], sigptr_fg[theSq_fg], nz_fg[theSq_fg]);
	   
	   if (surfval < -9998) {
	   
	      mask_it = ncparms_data.depths[ilev] > bpress_fg[theSq_fg] ? 1 : 0;
	      if (mask_it)  {
	           mask_deeper_levs(ilev, &ncparms_info, &ncparms_data);
	           goto NEXTSQUARE;
	      }
	      if (ncparms_data.depths[ilev] < p_fg[theSq_fg][0]) {
	          ++nempty;
	      } 
	      continue; 
	   }
	       
	   ngood = getsurf2d(surfval,reflev,sigptr_fg, th_fg,p_fg, &fg_info, nz_fg, bpress_fg, tsurf_fg);

	   mask_it = tsurf_fg[theSq_fg] > testmask ? 1 : 0;
	   if (mask_it)  {
	           mask_deeper_levs(ilev, &ncparms_info, &ncparms_data);
	           goto NEXTSQUARE;
	   }
	   
	   
	   ngood = getsurf2d(surfval, reflev, sigptr_fg, s_fg,  p_fg, &fg_info, nz_fg, bpress_fg, ssurf_fg);
	   ngood = getsurf2d(surfval, reflev, sigptr_fg, p_fg,  p_fg, &fg_info, nz_fg, bpress_fg, psurf_fg);
	   surfmask = set_mask2d(tsurf_fg, &fg_info,theLat,theLon);
	   surfmask2 = set_mask2d_dist(&fg_info, xdist[ilev], ydist[ilev], theLat, theLon);
		
	   for (j = 0; j < nwsq_fg; ++j) {
	         if (surfmask[j] > 0 || surfmask2[j] > 0) {
		     tsurf_fg[j] = (double)HBMASK;
		     ssurf_fg[j] = (double)HBMASK;
		     psurf_fg[j] = (double)HBMASK;
		 }
	   } 
	   free(surfmask);
	   free(surfmask2);
	   
	   theFG_t = tsurf_fg[theSq_fg]; 
	        if (theFG_t < testempty) {
	          theFG_t = get_mean(&fg_info, tsurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }

	   theFG_s = ssurf_fg[theSq_fg]; 
	        if (theFG_s < testempty) {
	          theFG_s = get_mean(&fg_info, ssurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }
		
	   theFG_p = psurf_fg[theSq_fg]; 
	        if (theFG_p < testempty) {
	          theFG_p = get_mean(&fg_info, psurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }

		
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, th_b3d, p_b3d, &b3d_info,nz_b3d, bpress_b3d, tsurf_b3d);
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, s_b3d, p_b3d, &b3d_info, nz_b3d, bpress_b3d, ssurf_b3d );
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, p_b3d, p_b3d, &b3d_info, nz_b3d, bpress_b3d, psurf_b3d );

	   surfmask = set_mask2d(tsurf_b3d, &b3d_info,theLat, theLon);
	   surfmask2 = set_mask2d_dist(&b3d_info, xdist[ilev], ydist[ilev], theLat, theLon);

	   for (isq = 0; isq < nwsq_b3d; ++isq) {
	         if (surfmask[isq] > 0 || surfmask2[isq] > 0) {
		     tsurf_b3d[isq] = (double)HBMASK;
		     ssurf_b3d[isq] = (double)HBMASK;
		     psurf_b3d[isq] = (double)HBMASK;
		 }
	   } 
	   free(surfmask);
	   free(surfmask2);
	     
   /* Compute meridional and zonal property gradients */
	     
	     get_gradients(tsurf_fg, &fg_info, theLat, theLon, ncparms_info.fill_value, iwindow, &parm1, &parm2);
	        ncparms_data.dTdx[ilev] = parm1;
		ncparms_data.dTdy[ilev] = parm2;
	     if (plotflag ) {
		   fprintf(plot4, "%5.1lf %5.1lf %5.0f %.3f %.3f\n", theLon, theLat, ncparms_data.depths[ilev], parm1, parm2);
	     }
	     get_gradients(psurf_fg, &fg_info, theLat, theLon, ncparms_info.fill_value, iwindow, &parm1, &parm2);
	        ncparms_data.dPdx[ilev] = parm1;
	        ncparms_data.dPdy[ilev] = parm2;
	     
    /* Add other fields for this level */	     
	     ncparms_data.Lx[ilev] = (float) xdist[ilev];
	     ncparms_data.Ly[ilev] = (float) ydist[ilev];
	   
		
      /* Compute residuals */ 
	     
	   for (isq = 0; isq < nwsq_b3d; ++isq) {  
	        error = sq2rc(isq, &b3d_info, &row, &col);
	        error += ij2xy(&b3d_info, col, row, &lon, &lat);
	        error += xy2ij(&fg_info, lon, lat, &col, &row);
	     
	        if (error) {
		  dtsurf[0][imonth][isq] = (double)HBEMPTY;
		  dssurf[0][imonth][isq] = (double)HBEMPTY;
		  dpsurf[0][imonth][isq] = (double)HBEMPTY;
		}
	        else {
		   isq_fg = row * fg_info.nx + col;
		   if (tsurf_b3d[isq] > testmask || tsurf_fg[isq_fg] > testmask) {
		    dtsurf[0][imonth][isq] = (double)HBMASK;
		    dssurf[0][imonth][isq] = (double)HBMASK;
		    dpsurf[0][imonth][isq] = (double)HBMASK;
		   }
		   else if (tsurf_b3d[isq] < testempty ) {
		    dtsurf[0][imonth][isq] = (double)HBEMPTY;
		    dssurf[0][imonth][isq] = (double)HBEMPTY;
		    dpsurf[0][imonth][isq] = (double)HBEMPTY;
		   }
		   else {
		    dtsurf[0][imonth][isq] = tsurf_b3d[isq] - theFG_t;
		    dssurf[0][imonth][isq] = ssurf_b3d[isq] - theFG_s;
		    dpsurf[0][imonth][isq] = psurf_b3d[isq] - theFG_p;
		   } 
	        }
	   } /* end for isq */
	   
    /* Estimate Temperature variogram and execute a grid search to find parameters that minimize errors */
	     
	     maxlag = xdist[ilev] > ydist[ilev] ? xdist[ilev] : ydist[ilev];
	     if (dcorr_max > maxlag)
	        maxlag = dcorr_max;
		
	     
	     vgram_ptr = vgram_estimate(&dtsurf[0][imonth], 1, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
	     
	     /* check for valid variogram */
	     
	     maxlag2chk = MAXLAG2CHK > maxlag ? (int)maxlag : MAXLAG2CHK;
	     if (! vgram_check(vgram_ptr, nbinthresh, minobs, maxlag2chk) ) {
	        ++nempty;
	        continue;
	     }
	
	     error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'T');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
	     ncparms_data.Tparm0[ilev] = params[0];
	     ncparms_data.Tparm1[ilev] = params[1];
	     ncparms_data.Tparm2[ilev] = params[2];
	     if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
	        fprintf(plot1, "%5.1lf %5.1lf %5.0f %10.6f %10.6f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
		
		fprintf(plot5, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f %5d", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2], vgram_ptr->nlags);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->Vh[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->h[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %5d", vgram_ptr->Nh[k]);
		   fprintf(plot5,"\n");

	     }
	     ++nfilled;
	     free(vgram_ptr);
	     
    /* Estimate Salinity variogram */	     
	   
	     vgram_ptr = vgram_estimate(&dssurf[0][imonth], 1, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
	     error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'S');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);

	     ncparms_data.Sparm0[ilev] = params[0];
	     ncparms_data.Sparm1[ilev] = params[1];
	     ncparms_data.Sparm2[ilev] = params[2];
             if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
	       fprintf(plot2, "%5.1lf %5.1lf %5.0f %8.5f  %8.5f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
	     }
	     free(vgram_ptr);

    /* Estimate Pressure variogram */	     
	   
	     vgram_ptr = vgram_estimate(&dpsurf[0][imonth], 1, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
	     
	     error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'P');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
	     
	     ncparms_data.Pparm0[ilev] = params[0];
	     ncparms_data.Pparm1[ilev] = params[1];
	     ncparms_data.Pparm2[ilev] = params[2];
	     if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
	        fprintf(plot3, "%5.1lf %5.1lf %5.0f %5.0f  %5.0f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
		fprintf(plot7, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f %5d", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2], vgram_ptr->nlags);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %.6f", vgram_ptr->Vh[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %.6f", vgram_ptr->h[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %5d", vgram_ptr->Nh[k]);
		   fprintf(plot7,"\n");

	     }

	     free(vgram_ptr);
	     
        }  /* end for ilev below seasonal layer */
	
	
/*--------------------------------------------*/  
   /* Compute 1 variogram for all levels below 4000 */
/*--------------------------------------------*/  
      nz = 0;   /* count levels that have enough data points */
      ideepest = iz4000 -1;
      for (ilev = iz4000; ilev < ncparms_info.nz; ++ilev) {
      
  	     fprintf(stderr,"-");   /* show progress */
          
	   reflev = 4000;
	   sigptr_fg = sig4_fg;
	   sigptr_b3d = sig4_b3d;
	   
	   surfval = hb_linterp((double) ncparms_data.depths[ilev], p_fg[theSq_fg], sigptr_fg[theSq_fg], nz_fg[theSq_fg]);
	   
	   if (surfval < -9998) {
	   
	      mask_it = ncparms_data.depths[ilev] > bpress_fg[theSq_fg] ? 1 : 0;
	      if (mask_it)  {
	           mask_deeper_levs(ilev, &ncparms_info, &ncparms_data);
	           goto DEEP_VGRAM;
	      } /* end if mask_it */
	   }
	       
	   ngood = getsurf2d(surfval,reflev,sigptr_fg, th_fg,p_fg, &fg_info, nz_fg, bpress_fg, tsurf_fg);
	   ngood = getsurf2d(surfval, reflev, sigptr_fg, s_fg,  p_fg, &fg_info, nz_fg, bpress_fg, ssurf_fg);
	   ngood = getsurf2d(surfval, reflev, sigptr_fg, p_fg,  p_fg, &fg_info, nz_fg, bpress_fg, psurf_fg);

	   surfmask = set_mask2d(tsurf_fg, &fg_info,theLat,theLon);
	   surfmask2 = set_mask2d_dist(&fg_info, xdist[ilev], ydist[ilev], theLat, theLon);
		
	   for (j = 0; j < nwsq_fg; ++j) {
	         if (surfmask[j] > 0 || surfmask2[j] > 0) {
		     tsurf_fg[j] = (double)HBMASK;
		     ssurf_fg[j] = (double)HBMASK;
		     psurf_fg[j] = (double)HBMASK;
		 }
	   } 
	   free(surfmask);
	   free(surfmask2);
	   
      /* Compute meridional and zonal property gradients */
	     
	     get_gradients(tsurf_fg, &fg_info, theLat, theLon, ncparms_info.fill_value, iwindow, &parm1, &parm2);
	        ncparms_data.dTdx[ilev] = parm1;
		ncparms_data.dTdy[ilev] = parm2;
	     if (plotflag && (parm1 > testempty) && (parm1 < testmask)) {
		   fprintf(plot4, "%5.1lf %5.1lf %5.0f %f %f \n", theLon, theLat, ncparms_data.depths[ilev], parm1, parm2);
	     }
	     
    /* Add other fields for this level */	     
	     ncparms_data.Lx[ilev] = (float) xdist[ilev];
	     ncparms_data.Ly[ilev] = (float) ydist[ilev];
	     
	     
      /* Determine a First Guess value for this location */
	   
	     theFG_t = tsurf_fg[theSq_fg]; 
	        if (theFG_t < testempty) {
	          theFG_t = get_mean(&fg_info, tsurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }
		
	     theFG_s = ssurf_fg[theSq_fg]; 
	        if (theFG_s < testempty) {
	          theFG_s = get_mean(&fg_info, ssurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }
		
	     theFG_p = psurf_fg[theSq_fg]; 
	        if (theFG_p < testempty) {
	          theFG_p = get_mean(&fg_info, psurf_fg, nwsq_fg,  xdist[ilev], ydist[ilev]);
	        }
		
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, th_b3d, p_b3d, &b3d_info,nz_b3d, bpress_b3d, tsurf_b3d);
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, s_b3d, p_b3d, &b3d_info, nz_b3d, bpress_b3d, ssurf_b3d );
	   ngood = getsurf2d(surfval, reflev, sigptr_b3d, p_b3d, p_b3d, &b3d_info, nz_b3d, bpress_b3d, psurf_b3d );

	   surfmask = set_mask2d(tsurf_b3d, &b3d_info,theLat, theLon);
	   surfmask2 = set_mask2d_dist(&b3d_info, xdist[ilev], ydist[ilev], theLat, theLon);

	   for (isq = 0; isq < nwsq_b3d; ++isq) {
	         if (surfmask[isq] > 0 || surfmask2[isq] > 0) {
		     tsurf_b3d[isq] = (double)HBMASK;
		     ssurf_b3d[isq] = (double)HBMASK;
		     psurf_b3d[isq] = (double)HBMASK;
		 }
	   } 
	   free(surfmask);
	   free(surfmask2);
	   
	         /* Compute residuals */ 
	     
	   for (isq = 0; isq < nwsq_b3d; ++isq) {  
	        error = sq2rc(isq, &b3d_info, &row, &col);
	        error += ij2xy(&b3d_info, col, row, &lon, &lat);
	        error += xy2ij(&fg_info, lon, lat, &col, &row);
	     
	        if (error) {
		  dtsurf[0][nz][isq] = (double)HBEMPTY;
		  dssurf[0][nz][isq] = (double)HBEMPTY;
		  dpsurf[0][nz][isq] = (double)HBEMPTY;
		}
	        else {
		   isq_fg = row * fg_info.nx + col;
		   if (tsurf_b3d[isq] > testmask || tsurf_fg[isq_fg] > testmask) {
		    dtsurf[0][nz][isq] = (double)HBMASK;
		    dssurf[0][nz][isq] = (double)HBMASK;
		    dpsurf[0][nz][isq] = (double)HBMASK;
		   }
		   else if (tsurf_b3d[isq] < testempty ) {
		    dtsurf[0][nz][isq] = (double)HBEMPTY;
		    dssurf[0][nz][isq] = (double)HBEMPTY;
		    dpsurf[0][nz][isq] = (double)HBEMPTY;
		   }
		   else {
		    dtsurf[0][nz][isq] = tsurf_b3d[isq] - theFG_t;
		    dssurf[0][nz][isq] = ssurf_b3d[isq] - theFG_s;
		    dpsurf[0][nz][isq] = psurf_b3d[isq] - theFG_p;
		   } 
	        }
	   } /* end for isq */
	   
	   ideepest = ilev;
	   ++nz;
	   
      }  /* end for ilev */

/*--------------------------------------------*/  
      DEEP_VGRAM:
           ;
/*--------------------------------------------*/  
       /* After residuals have been computed for all deep levels  
          compute a single variogram and apply to all deep levels  */ 
	  
	   
      if (nz > 0) {   
         vgram_ptr = vgram_estimate(dtsurf[0], nz, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);

	     /* check for valid variogram */
	     
	     maxlag2chk = MAXLAG2CHK > maxlag ? (int)maxlag : MAXLAG2CHK;
	     if (! vgram_check(vgram_ptr, nbinthresh, minobs, maxlag2chk) ) {
	        nempty += nz;
	        goto NEXTSQUARE;
	     }
	

         error = get_parm_range(imodel, vgram_ptr, range1, range2);
	 
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'T');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);

         free(vgram_ptr);

         for (ilev = iz4000; ilev <= ideepest; ++ilev) {   
       	     ncparms_data.Tparm0[ilev] = params[0];
       	     ncparms_data.Tparm1[ilev] = params[1];
	     ncparms_data.Tparm2[ilev] = params[2];
	     ++nfilled;
	     if (plotflag) {
                 fprintf(plot1, "%5.1lf %5.1lf %5.0f %10.6f %10.6f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);	 
		 fprintf(plot5, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f %5d", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2], vgram_ptr->nlags);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->Vh[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %.6f", vgram_ptr->h[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot5, " %5d", vgram_ptr->Nh[k]);
		   fprintf(plot5,"\n");
       	     }
         }
           
         vgram_ptr = vgram_estimate(dssurf[0], nz, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
         error = get_parm_range(imodel, vgram_ptr, range1, range2);
	 
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
 	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'S');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
         free(vgram_ptr);

         for (ilev = iz4000; ilev < ideepest; ++ilev) {   
       	     ncparms_data.Sparm0[ilev] = params[0];
       	     ncparms_data.Sparm1[ilev] = params[1];
	     ncparms_data.Sparm2[ilev] = params[2];
	     if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
	        fprintf(plot2, "%5.1lf %5.1lf %5.0f %8.5f %8.5f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
	     }
         }
      
      
         vgram_ptr = vgram_estimate(dpsurf[0], nz, lats_b3d, lons_b3d, b3d_info.ny, b3d_info.nx, maxlag, binsize);
         error = get_parm_range(imodel, vgram_ptr, range1, range2);
	     range2[0] = dcorr_min;  /* adjust these to locally specified parameters */
	     range2[1] = dcorr_max;  
	     range2[2] = binsize;
	 
	     if (snflag)
	        params[0] = snratio(theLat, ncparms_data.depths[ilev], &sn_data, 'P');
	     
	     grid_search(imodel, vgram_ptr, range1, range2, params, snflag);
         free(vgram_ptr);

         for (ilev = iz4000; ilev < ideepest; ++ilev) {   
       	     ncparms_data.Pparm0[ilev] = params[0];
       	     ncparms_data.Pparm1[ilev] = params[1];
	     ncparms_data.Pparm2[ilev] = params[2];
	     if (plotflag && (params[1] < testmask) && (params[1] > testempty)) {
	        fprintf(plot3, "%5.1lf %5.1lf %5.0f %5.0f %5.0f %4.0f \n", theLon, theLat, ncparms_data.depths[ilev], params[0], params[1], params[2]);
		fprintf(plot7, "%5.1lf %5.1lf %5.0f  %10.6f  %10.6f %4.0f %5d", theLon, theLat, ncparms_data.depths[ilev], params[0] , params[1], params[2], vgram_ptr->nlags);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %.6f", vgram_ptr->Vh[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %.6f", vgram_ptr->h[k]);
		  for (k=0; k < vgram_ptr->nlags; ++k)
		    	fprintf(plot7, " %5d", vgram_ptr->Nh[k]);
		   fprintf(plot7,"\n");
	     }
         }
	 
      }  /* end if nz > 0 */

/*--------------------------------------------*/  
	NEXTSQUARE:
	   ;
/*--------------------------------------------*/  

    /*  Add bathymetry gradient : ncparms_data.dZdx */
    
        if (ilev > 0) {
	     get_gradients(bpress_fg, &fg_info, theLat, theLon, ncparms_info.fill_value, iwindow, &parm1, &parm2);
	     ncparms_data.dZdx = parm1;
	     ncparms_data.dZdy = parm2;
	   if (plotflag ) {
	     fprintf(plot6, "%5.1lf %5.1lf %5.0lf %5.0f %5.0f \n", theLon, theLat, bpress_fg[theSq_fg], ncparms_data.dZdx, ncparms_data.dZdy);
	   }
	}
	

	ncparms_write(ncid,(float)theLat, (float)theLon, &ncparms_info, &ncparms_data); 
	  
/*--------------------------------------------*/    
	
    } /* loop for each sq in region */ 

    ncparms_close(ncid);
    if (plotflag) {
       fclose(plot1);
       fclose(plot2);
       fclose(plot3);
       fclose(plot4);
       fclose(plot5);
       fclose(plot6);
       fclose(plot7);
    }
    
    fprintf(stdout,"\nNumber of nodes filled: %d\n", nfilled);
    fprintf(stdout,"Number of nodes masked: %d\n", nmasked);
    fprintf(stdout,"Number of nodes empty: %d\n", nempty);
    
    
    fprintf(stderr,"\nEnd of %s\n", argv[0]);
    
}  /* end main */

/************************************************************************/
void print_usage(char *program)
{
  fprintf(stderr,"\n\nUSAGE:  %s [-A |-C |-U] -B<w/e/s/n> [-D<maxdist>] -Fr<FG_rootname> [-Fd<FG_dirname>] [-Fe<FG_file_extent>] -Gr<bin3d_rootname> [-Gd<bin3d_dirname>] [-Ge<bin3d_file_extent>] [-H<binsize>] [-I<xinc/<yinc>] [-M<model_code>] [-N<filename>] -O<outfile> [-P<plotfile_rootname>] [-T<topofile>] -Z<stddepth_file> [-h]", program);
  
   fprintf(stderr,"\n\nCREATE a global parameters file for future updating: ");
   fprintf(stderr,"\n  %s -C -O<outfile> -Z<stddepth_file> [-I<xinc/<yinc>] ", program);
   fprintf(stderr,"\n\n UPDATE an existing parameters file with new covariances within the specified bounds\n");
   fprintf(stderr,"\n  %s -U -B<w/e/s/n> -Fr<FG_rootname>  -Gr<bin3d_rootname> -O<outfile>  -Z<stddepth_file> ", program);
   fprintf(stderr,"\n\n Use optional -Ui to interpolate missing nodes with parameters from nearest neighbors:");
   fprintf(stderr,"\n  %s -Ui -B<w/e/s/n> -O<outfile>  ", program);
   fprintf(stderr,"\n\n Use optional -Un<filename> to apply estimated signal-to-noise ratios in the specified netcdf file:");
   fprintf(stderr,"\n  %s -UnSNratios_atlantic.nc -B<w/e/s/n> -O<outfile> -Fr<FG_rootname>  -Gr<bin3d_rootname> -Z<stddepth_file> ", program);

   fprintf(stderr,"\n\n  OPTION FLAGS ", program);
   fprintf(stderr," -B  sets the boundaries of the region to work on: w/e/s/n.\n");
   fprintf(stderr," -C  CREATE a new parameters file.\n");
   fprintf(stderr," -D  sets the maximum search radius (in km). Default: %d\n", MAXDIST);
   fprintf(stderr," -F  specify names of HydroBase netcdf files to be used as first guess fields.\n");
   fprintf(stderr,"      One file for each of 12 months is required. \n");
   fprintf(stderr,"      filename has form: <dir>/<root>.<month>.<extent> \n");
   fprintf(stderr,"      Default extent is '.nc' \n");
   fprintf(stderr,"      Use -Fr to specify root or dir/root (required)\n");
   fprintf(stderr,"      Use -Fd to specify directory (optional)\n");
   fprintf(stderr,"      Use -Fe to specify file extent (optional)\n");
   fprintf(stderr," -G  names of HydroBase bin3d files.\n");
   fprintf(stderr,"      One file for each of 12 months is required. \n");
   fprintf(stderr,"      filename has form: <dir>/<root>.<month>.<extent> \n");
   fprintf(stderr,"      Default extent is '.nc' \n");
   fprintf(stderr,"      Use -Gr to specify root or dir/root (required)\n");
   fprintf(stderr,"      Use -Gd to specify directory (optional)\n");
   fprintf(stderr,"      Use -Ge to specify file extent (optional)\n");
   fprintf(stderr," -H  sets the binsize (in km) for computing lags, Default: %d\n", BINSIZE);
   fprintf(stderr," -I  sets the lat/lon grid ( in degrees)for paramters. Default: [1/1]\n");
   fprintf(stderr," -L  sets a range for decorrelation scales (in km) Default: %d/%d\n", MIN_DECORR_LENGTH, MAX_DECORR_LENGTH);
   fprintf(stderr," -M  model code:  1 = exponential, 2 = gaussian, 3 = power. Default: [1]\n");
   fprintf(stderr," -O  name of output parameters file.\n");
   fprintf(stderr," -P  Create 6 plotting files with specified rootname:\n");
   fprintf(stderr,"         *.Tparms   *.Sparms  *.Pparms\n");
   fprintf(stderr,"         *.Tgrad    *.Pgrad   *.Bathgrad\n");
   fprintf(stderr," -T  specify a topography file.\n");
   fprintf(stderr,"     default: %s\n", BATHPATH_C);
   fprintf(stderr," -U  UPDATE an existing parameters file with new covariance parameters.\n");
   fprintf(stderr,"     Append 'i' to interpolate parameters at missing nodes.\n");
   fprintf(stderr,"     Append 'n' followed by filename to specify netcdf file containing profiles of characteristic signal-to-noise ratios.\n");
   fprintf(stderr,"       Note that dimensions cannot be changed.  \n");
   fprintf(stderr," -Z  file containing list of standard depths,and x-, y-length scales to compute parameters\n");
   fprintf(stderr,"       Each line contains depth, xlength, ylength \n");
   fprintf(stderr,"       Within a parameter file, list of depths is set, but length scales can vary. \n");

} /* end print_usage() */

/*****************************************************************************/
void create_ncparms(struct NC_INFO *ncinfo_ptr, struct NC_PARMS *ncdata_ptr, struct GRID_INFO *ginfo_ptr, double *depths)

   /* sets field of ncinfo and ncdata before creating nc file. 
      Nodes are centered on grid cells:  
          lats[0] = ymin + 0.5 * yinc 
	  lons[0] = xmin + 0.5 * xinc
	  
	  By definition, longitudes are on the range 0 to 360.
      */
{
   float xmin, xmax, ymin, ymax, xinc, yinc;
   int i;
   
   xmin = (float) ginfo_ptr->x_min;
   xmax = (float) ginfo_ptr->x_max;
   ymin = (float) ginfo_ptr->y_min;
   ymax = (float) ginfo_ptr->y_max;
   xinc = (float) ginfo_ptr->x_inc;
   yinc = (float) ginfo_ptr->y_inc;
   
   if (xmin < 0)
      xmin += 360.0;
      
   if (xmax < 0)
      xmax += 360.0;
        
   ncinfo_ptr->xincr = xinc;
   ncinfo_ptr->yincr = yinc;	     
   ncinfo_ptr->ny = NINT( (ymax-ymin)/ yinc);
   ncinfo_ptr->nx = NINT( (xmax-xmin)/ xinc);
   /* nz is set in get_depths() */
   
   ncparms_alloc(ncinfo_ptr, ncdata_ptr);  /* allocates memory for structure fields */
   
   
   for (i= 0; i < ncinfo_ptr->nx; ++i)
       ncinfo_ptr->lons[i] = xmin + (i + 0.5) * xinc;  
       
   for (i= 0; i < ncinfo_ptr->ny; ++i)
       ncinfo_ptr->lats[i] = ymin + (i +0.5) * yinc;  
   
   ncinfo_ptr->fill_value = NCPARMS_FILL;
   ncinfo_ptr->mask_value = NCPARMS_MASK;
   
   
   for (i=0; i < ncinfo_ptr->nz; ++i) 
     ncdata_ptr->depths[i] = (float) depths[i];
   
   
   i = 0;
   while (ncdata_ptr->depths[i] < SEASONAL_LAYER_DEF)
      ++i;
      
   ncinfo_ptr->nz_seas = i;
   ncdata_ptr->nz_seas = i;
   ncdata_ptr->nz = ncinfo_ptr->nz;
   return;  

} /* end create_ncparms() */
/*****************************************************************************/
void mask_land_nodes(int ncid, struct NC_INFO *ncinfo_ptr, struct NC_PARMS *ncdata_ptr, struct GRID_INFO *topo_info, short *seafloor)

{
   int row, col, ilev, error;
   short bdepth;
   float flat, flon;
   
   ilev = 0;
   mask_deeper_levs(ilev, ncinfo_ptr, ncdata_ptr);   /* create a masked node */
   
   for (row = 0; row < ncinfo_ptr->ny;  ++row) {
      for (col = 0; col < ncinfo_ptr->nx; ++col) {
         error = ncparms_latlon(ncinfo_ptr, row, col, &flat, &flon);
	 
	  bdepth = find_nearest_topo_val((double)flat, (double)flon, seafloor, topo_info);   
          if (bdepth <= 0 )  /* seafloor depths are positive */
	       ncparms_write(ncid, flat, flon, ncinfo_ptr, ncdata_ptr);
      }
   }
   
   return;
} /* end mask_land_nodes */
/*****************************************************************************/

int check_overlap(struct GRID_INFO *area1, struct GRID_INFO *area2)
   /* determines whether area1 and area2 overlap.
     Returns 1 (TRUE) or 0 (FALSE) */
{

   double xmin1, xmin2, xmax1, xmax2;
   
   xmin1 = area1->x_min;
   xmin2 = area2->x_min;
   xmax1 = area1->x_max;
   xmax2 = area2->x_max;
   
   /* adjust for longitude conventions  */
   
   if (area1->lon0to360 != area2->lon0to360) {
     if (area1->lon0to360){
        if (xmin2 < 0)
	   xmin2 +=360;
	if (xmax2 < 0)
	   xmax2 += 360;
     }
   }
   
   if (xmin1 > xmax2 )  return(FALSE);
   if (xmax1 < xmin2 ) return(FALSE);
   if (area1->y_min > area2->y_max) return(FALSE);
   if (area1->y_max < area2->y_min) return(FALSE);
   return(TRUE);

} /* end check_overlap() */

/*****************************************************************************/
void set_topo_info(struct GRID_INFO *ginfo, struct GRID_INFO *tinfo, int topobuf)   
   /* for now, if pole is crossed, limit the y bound to the pole */
{
   int error, xpole;
   int row, col;
   double dummy;
   double lat, lon, lat_n, lat_s, lon_e, lon_w;

    
    if (topobuf == 0) {
        tinfo->x_min = ginfo->x_min;
	tinfo->y_min = ginfo->y_min;
        tinfo->x_max = ginfo->x_max;
	tinfo->y_max = ginfo->y_max;
	
	if (tinfo->x_max >= 360)
	   tinfo->x_max = 359.999;
        if (tinfo->y_max  >= 90)
	    tinfo->y_max = 89.999;
	   
	lon_w = tinfo->x_min;
	lon_e = tinfo->x_max;
	lat_n = tinfo->y_max;
	lat_s = tinfo->y_min;
	
     }
     else {
	
	xpole = pointB(ginfo->y_max, ginfo->x_max, 0, (double) topobuf, TRUE, &lat_n, &dummy);
	if (xpole) 
	   lat_n = 90.0;
	   
        xpole = pointB(ginfo->y_min, ginfo->x_min, 180,(double) topobuf , TRUE, &lat_s, &dummy);
	if (xpole)
	   lat_s = -90.0;
	   
	if (ABS(ginfo->y_max) > ABS(ginfo->y_min)) {
	   error = pointB(ginfo->y_max, ginfo->x_min, 270, (double) topobuf, TRUE, &dummy, &lon_w);
           error = pointB(ginfo->y_max, ginfo->x_max, 90, (double) topobuf, TRUE, &dummy, &lon_e);
	}
	else {
	   error = pointB(ginfo->y_min, ginfo->x_min, 270, (double) topobuf, TRUE, &dummy, &lon_w);
           error = pointB(ginfo->y_min, ginfo->x_max, 90, (double) topobuf, TRUE, &dummy, &lon_e);
        }
    
     }
     
    /* adjust bounds to correspond to actual grid nodes */
    
    xy2ij_nochk(ginfo, lon_e, lat_n, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_e, &lat_n);
    
    xy2ij_nochk(ginfo, lon_w, lat_s, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_w, &lat_s);
    
    
  
    tinfo->y_min = lat_s - ginfo->node_offset * 0.5 * ginfo->y_inc;
    tinfo->y_max = lat_n + ginfo->node_offset * 0.5 * ginfo->y_inc;
    tinfo->x_min = lon_w - ginfo->node_offset * 0.5 * ginfo->x_inc;
    tinfo->x_max = lon_e + ginfo->node_offset * 0.5 * ginfo->x_inc;
    tinfo->lon0to360 = tinfo->x_min < 0 ? 0: 1;
    tinfo->xgreenwich = (tinfo->x_min < 0 && tinfo->x_max >= 0) ? 1 : 0;

     if (tinfo->y_max  >= 90)
	    tinfo->y_max = 89.999;
     if (tinfo->y_min  < -90)
	    tinfo->y_max = -90;
    if (tinfo->lon0to360) {    
        if (tinfo->x_max >= 360)
	   tinfo->x_max = 359.999;
        if (tinfo->x_min < -360)
	   tinfo->x_max = -360;
    }
	  
     return;
} /* end set_topo_info() */

/*****************************************************************************/

int get_size(struct GRID_INFO *ginfo, struct GRID_INFO *g2info, double maxdist)

/* determines size of work array which is maxdist beyond bounds of ginfo,
  
   Returns the number of gridsquares in work array, or 0 if not successful,
   Does not alter either ginfo or g2info  */

{
   int error, ncols, nrows, col, row;
   double lat, lon, lat_n, lat_s, lon_e, lon_w;


/* determine new bounds */

   
   error = pointB(ginfo->y_max, ginfo->x_min, 0, maxdist, TRUE, &lat_n, &lon);
	   
   error = pointB(ginfo->y_min, ginfo->x_max, 180, maxdist, TRUE, &lat_s, &lon);

   if (ginfo->y_min >= 0)  { /* Northern Hemisphere */
       error = pointB(ginfo->y_max, ginfo->x_min, 270, maxdist, TRUE, &lat, &lon_w);
       error = pointB(ginfo->y_max, ginfo->x_max, 90, maxdist, TRUE, &lat, &lon_e);
    }
    else {  /* southern hemisphere */
    
       error = pointB(ginfo->y_min, ginfo->x_min, 270, maxdist, TRUE, &lat, &lon_w);
       error = pointB(ginfo->y_min, ginfo->x_max, 90, maxdist, TRUE, &lat, &lon_e);
    }
    
    /* adjust bounds to correspond to actual grid nodes */
    
    xy2ij_nochk(ginfo, lon_e, lat_n, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_e, &lat_n);
    
    xy2ij_nochk(ginfo, lon_w, lat_s, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_w, &lat_s);
    
    
    nrows = (lat_n - lat_s) / g2info->y_inc;
    ncols = (lon_e - lon_w) / g2info->x_inc;
    
    if (! ginfo->node_offset) {
       ++nrows;
       ++ncols;
    }
	  
    return(ncols*nrows);

} /* end get_size() */
/*****************************************************************************/
int get_bounds(struct GRID_INFO *ginfo_ptr, int sq, double xdist,double ydist, double *lat0, double *lon0, double *lat1, double *lon1)
/* Returns minlat,minlon, maxlat, maxlon of area around sq of length xdist,ydist. 
   Returns 1 if pole was crossed, 0 if not.
  */
{
  int row, col, xpole;
  double lon, lat, dummy;
  double minlon, minlat, maxlon, maxlat;
  
  col = sq % ginfo_ptr->nx;
  row = (sq - col) / ginfo_ptr->nx;

  ij2xy(ginfo_ptr, col, row, &lon, &lat);
  minlon = lon - 0.5 * ginfo_ptr->node_offset * ginfo_ptr->x_inc;
  minlat = lat - 0.5 * ginfo_ptr->node_offset * ginfo_ptr->y_inc;
  maxlon = lon + 0.5 * ginfo_ptr->node_offset * ginfo_ptr->x_inc;
  maxlat = lat + 0.5 * ginfo_ptr->node_offset * ginfo_ptr->y_inc;
  xpole = pointB(minlat, lon, 180, ydist, TRUE, lat0, &dummy);
  xpole += pointB(maxlat, lon, 0, ydist, TRUE, lat1, &dummy);
  pointB(lat, minlon, 270, xdist, TRUE, &dummy, lon0);
  pointB(lat, maxlon, 90, xdist, TRUE, &dummy, lon1);

    /* adjust bounds to correspond to actual grid nodes */
    
    xy2ij_nochk(ginfo_ptr, *lon0, *lat0, &col, &row);
    ij2xy_nochk(ginfo_ptr, col, row, lon0, lat0);
    xy2ij_nochk(ginfo_ptr, *lon1, *lat1, &col, &row);
    ij2xy_nochk(ginfo_ptr, col, row, lon1, lat1);

  return(xpole);
  
  
} /* get_bounds */
/*****************************************************************************/
int set_work_grid(struct GRID_INFO *gptr, struct GRID_INFO *g2ptr,double *lat0, double *lon0, double *lat1, double *lon1, int xpole)

/*  Sets the fields of g2ptr based on specified bounds. Returns the number of gridsquares in g2ptr.
   Next it checks and adjusts the specified bounds to remain within the boundaries of gptr. 
   Note that g2ptr bounds may be outside of gptr bounds, and can't be changed because the estimator 
   functions are set up by definition to work on the node in the center of the grid. */
{
   int row, col, error, ncols, nrows;
   double lat, lon, dy;

   g2ptr->node_offset = gptr->node_offset;  
   g2ptr->x_inc = gptr->x_inc;
   g2ptr->y_inc = gptr->y_inc;
   
   xy2ij_nochk(gptr, *lon0, *lat0, &col, &row);
   ij2xy_nochk(gptr, col, row,  &g2ptr->x_min, &g2ptr->y_min);
   xy2ij_nochk(gptr, *lon1, *lat1, &col, &row);
   ij2xy_nochk(gptr, col, row,  &g2ptr->x_max, &g2ptr->y_max);
    
   if (g2ptr->node_offset) {   /* for pixel grids, set bounds properly */
      
      g2ptr->x_min -=  0.5 * g2ptr->x_inc;
      g2ptr->x_max += 0.5 * g2ptr->x_inc;
      g2ptr->y_min -=  0.5 * g2ptr->y_inc;
      g2ptr->y_max +=  0.5 * g2ptr->y_inc;      
   }
   
   g2ptr->nx = NINT((g2ptr->x_max - g2ptr->x_min) / g2ptr->x_inc);
   g2ptr->ny = NINT((g2ptr->y_max - g2ptr->y_min) / g2ptr->y_inc);
      
   if (xpole) {
   
       if (g2ptr->y_min > 0) { /* northern hem */
        dy = 90 - g2ptr->y_max;
	g2ptr->y_max = 90 + dy;
       }
       else {
        dy =  g2ptr->y_min + 90 ;
	g2ptr->y_min = -90 - dy;
   
       }
       g2ptr->ny = NINT((g2ptr->y_max - g2ptr->y_min) / g2ptr->y_inc);
   }
   
   if (!g2ptr->node_offset) {  /* grid registered */ 
       ++g2ptr->nx;
       ++g2ptr->ny;
   }
   
   g2ptr->lon0to360 = g2ptr->x_min < 0 ? 0 : 1;
   g2ptr->xgreenwich = (g2ptr->x_min < 0) && (g2ptr->x_max >= 0);
   
   
   if (*lat0 < gptr->y_min)
       *lat0 = gptr->y_min + 0.5 * g2ptr->node_offset * g2ptr->y_inc;

   if (*lon0 < gptr->x_min)
       *lon0 = gptr->x_min + 0.5 * g2ptr->node_offset * g2ptr->x_inc;
       
   if (*lon1 > gptr->x_max)
       *lon1 = gptr->x_max - 0.5 * g2ptr->node_offset * g2ptr->x_inc;
       
   if (*lat1 > gptr->y_max)
       *lat1 = gptr->y_max - 0.5 * g2ptr->node_offset * g2ptr->y_inc;
   
   return(g2ptr->nx * g2ptr->ny);

} /* end set_work_grid() */

/*****************************************************************************/
void get_profiles(int nz, double *p, double *t, double *s, double *th, double *s0, double *s1,  double *s2, double *s3, double *s4)

   /* computes profiles accordingly */
{
  int i, j;
  if (nz > 0) {
     compute_theta(nz, th, p, t, s);
     compute_sigma(0., nz, s0, p, t, s);
     compute_sigma(1000., nz, s1, p, t, s);
     compute_sigma(2000., nz, s2, p, t, s);
     compute_sigma(3000., nz, s3, p, t, s);
     compute_sigma(4000., nz, s4, p, t, s);
   }
   
 return;  
		
} /* end get_profiles() */

/*****************************************************************************/
int get_depths(FILE *zfile, double **zptr, double **xlenptr, double **ylenptr)
   /* Reads an already opened file containing columns of 
           depth, xlength, ylength
      for each level to compute parameters -- x, y-lengths (search ellipse) are in km.
      Allocates memory for the arrays, and returns number of zlevels. */
{  
   int nz, i;
   char str[1000];
   double *z, *x, *y;

   nz = 0;
   
   
   while (fscanf(zfile, "%[^\n]", str) == 1) {
       ++nz;
       fgetc(zfile);
   }
   rewind(zfile);
   *zptr = (double *) calloc(nz, sizeof(double));
   *xlenptr = (double *) calloc(nz, sizeof(double));
   *ylenptr = (double *) calloc(nz, sizeof(double));
   
   z = *zptr;
   x = *xlenptr;
   y = *ylenptr;
   
   for (i = 0; i < nz; ++i) 
      fscanf(zfile,"%lf %lf %lf\n", &z[i], &x[i], &y[i]);
   
   
   return(nz);

} /* end get_depths() */
/*****************************************************************************/
int get_sn_ratios(int ncid,struct SN_RATIOS *sn)
  /* Reads signal-to-noise info from an open netcdf file and stores it in the structure.
     Returns 0 if successful or 1 if an error occurred. Error messages are written to
     stderr */
{
   int error, varid, n;
   int depth_dim, lat_dim;
   int dim[MAX_VAR_DIMS];
   size_t len;
   
   /* get dimensions */
   error = nc_inq_dimid(ncid,"DEPTH",&depth_dim);
   if (error != NC_NOERR) {
      fprintf(stderr,"\n Error reading signal-to-noise ratio file\n");
       fprintf(stderr,"\n%s\n", nc_strerror(error));
      return(1);
   }
   error = nc_inq_dimlen(ncid,depth_dim, &len);
   sn->ndepths = len;
   
   error = nc_inq_dimid(ncid,"LAT_BAND",&lat_dim);
   if (error != NC_NOERR) {
      fprintf(stderr,"\n Error reading signal-to-noise ratio file\n");
      fprintf(stderr,"\n%s\n", nc_strerror(error));
      return(1);
   }
   error = nc_inq_dimlen(ncid,lat_dim, &len);
   sn->nlats = len;
   
   /* allocate memory */
   
   n = sn->nlats * sn->ndepths;
   sn->depths = (float *)calloc(sn->ndepths, sizeof(float));
   sn->latbands = (float *)calloc(sn->nlats, sizeof(float));
   sn->T = (float *)calloc(n, sizeof(float));
   sn->S = (float *)calloc(n, sizeof(float));
   sn->P = (float *)calloc(n, sizeof(float));

   /* read values into arrays */
    
    error = nc_inq_varid(ncid, "depth", &varid);
    if (error) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    if ((error = nc_get_var_float(ncid,varid,sn->depths)) != NC_NOERR) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    
    error = nc_inq_varid(ncid, "latitude", &varid);
    if (error) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    if ((error = nc_get_var_float(ncid,varid,sn->latbands)) != NC_NOERR) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }

    error = nc_inq_varid(ncid, "Tparams", &varid);
    if (error) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    if ((error = nc_get_var_float(ncid,varid,sn->T)) != NC_NOERR) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }

    error = nc_inq_varid(ncid, "Sparams", &varid);
    if (error) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    if ((error = nc_get_var_float(ncid,varid,sn->S)) != NC_NOERR) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    
    error = nc_inq_varid(ncid, "Pparams", &varid);
    if (error) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    if ((error = nc_get_var_float(ncid,varid,sn->P)) != NC_NOERR) {
       fprintf(stderr,"\n%s\n", nc_strerror(error));
       return(1);
    }
    
    return(0);

}  /* end get sn_ratios() */
/*****************************************************************************/
float snratio(double theLat, float theDepth, struct SN_RATIOS *sn, char prop)
/* Returns the signal-to-noise ratio for the latitude and depth or exits if
  an error occurs. */
{
   int row, col, index;
   
   row = NINT(floor((theLat + 80.0001) / 10.0));
   col = 0;
   while (col < sn->ndepths && theDepth > sn->depths[col])
       ++col;
   if (col >= sn->ndepths)
      col = sn->ndepths -1;
      
   index = row * sn->ndepths + col;
   if (index < 0 || index >= (sn->ndepths * sn->nlats)) {
      fprintf(stderr,"\nFATAL index error in snratio()\n");
      exit(-1); 
   }
   switch (prop) {
     case 'p':
     case 'P':
       return(sn->P[index]);
     case 's':
     case 'S':
       return(sn->S[index]);
     case 't':
     case 'T':
       return(sn->T[index]);
     default:
       fprintf(stderr,"\nFATAL ERROR: Unknown property passed to snratio()\n");
       exit(-1); 
   }
}
/*****************************************************************************/
void set_default_lengthscales(double **xlen_ptr, double **ylen_ptr, struct NC_PARMS *pptr)
           /* allocates arrays and sets to default values */
{
  int i;
  double *xdist, *ydist;
  
	  
  *ylen_ptr = (double *) calloc(pptr->nz, sizeof(double));
  *xlen_ptr = (double *) calloc(pptr->nz, sizeof(double));
	  
  xdist = *xlen_ptr;
  ydist = *ylen_ptr;
	  
  for (i= 0; i < pptr->nz; ++i) {
     xdist[i] = 500;
     ydist[i] = 500;
     if (pptr->depths[i] < 151)  {
	  xdist[i] = ydist[i] = 200;
     }
     else if (pptr->depths[i] < 1001) {
	  xdist[i] = ydist[i] = 300;
     }
  }
	  
  return;

} /*end set_default_lengthscales() */
/*****************************************************************************/
/*****************************************************************************/
int *set_mask2d(double *x, struct GRID_INFO *hptr, double theLat, double theLon)

/*  Sets mask values to -1, 0, or 1 based on distribution of masked, 
    empty nodes in x (testmask and testempty are global variables).  
    
    Returns a pointer to mask or NULL if an error occurs.
*/
{
   int nsq, sq, row, col, xmid, ymid, error;
   int row2, col2, sq2, rad;
   int dx, dy, mask_rest;
   int *mask; 
   double phi_min, phi_max, dist_min;  
   double offset;
   double *phi, *dist;
   double *latv,*lonv;

   nsq = hptr->nx * hptr->ny;
   offset = hptr->node_offset * 0.5; /* 0 if gridnode, 0.5 for pixel grids*/

   latv = (double *) calloc(hptr->ny, sizeof(double));
   lonv = (double *) calloc(hptr->nx, sizeof(double));
   for (col = 0; col < hptr->nx; ++col) 
          lonv[col] = hptr->x_min + (col + offset) * hptr->x_inc;
   
   
   /* create and populate arrays to store distance, azimuth  info */
   
   phi = (double *)calloc(nsq, sizeof(double));
   dist = (double *)calloc(nsq, sizeof(double));
   mask = (int *)calloc(nsq, sizeof(int));
   
   error = xy2ij(hptr, theLon, theLat, &xmid, &ymid);
   if (error) 
      return(NULL);
      
   for (row = 0; row < hptr->ny; ++row) {
      latv[row] = hptr->y_min + (row + offset) * hptr->y_inc;
      for (col = 0; col < hptr->nx; ++col) {
         sq = row * hptr->nx + col;
	 /* get cartesian distance */
	 dist[sq] = distance_c(latv[row],lonv[col], theLat, theLon, TRUE, &phi[sq]);
	 /* but use gridpoint azimuth -180 to +180 */
	 phi[sq] = atan2(row-ymid, col-xmid) * DEGperRAD;  
         if (x[sq] < testempty)   
            mask[sq] = -1;
         else if (x[sq] > testmask) 
            mask[sq] = 1;
      }
   } 
      
 /* at this point, every element of mask is set to either -1, 0 or 1 */
 

/* search west for masked nodes */
   row = ymid;
   col = xmid;
   dx = 0;
   mask_rest = 0;
   while (++dx <= xmid) {
      --col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
   } /* end while */
  
/*search east */
   row = ymid;
   col = xmid;
   dx = 0;
   mask_rest = 0;
   while (++dx <= xmid) {
      ++col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
      
   } /* end while */

/*search north */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      ++row;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
   } /* end while */
   
/*search south */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      --row;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
   } /* end while */

/*search southwest */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      --row;
      --col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
      
   } /* end while */

/*search northwest */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      ++row;
      --col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
  } /* end while */

 /*search northeast */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      ++row;
      ++col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
   } /* end while */
  
 /*search southeast */
   row = ymid;
   col = xmid;
   dy = dx = 0;
   mask_rest = 0;
   while (++dy <= ymid && ++dx <= xmid) {
      --row;
      ++col;
      sq = row * hptr->nx + col;
      if (sq < nsq) {
         if (mask_rest)
            mask[sq] = 1;
         else if (mask[sq] == 1)  
	    mask_rest = 1;
      }
   } /* end while */

  /* Check for consecutive masked nodes along a row or a column (a topographic barrier),
     and, where found, mask all nodes at greater distance between the 2 azimuths connecting the 
    center node and those masked nodes */
     
 /* first do quadrants east and west of center node*/ 
      
   for (rad=1; rad <= xmid; ++rad) {
     col = xmid + rad;
     col2 = xmid - rad;
     row = ymid + rad;
     row2 = ymid - rad;
     if (row >= hptr->ny) {
        row = hptr->ny - 1;
	row2 = 0;
     }
         /* work top to bottom */
     while ( row > row2) {   
           /* eastern quadrant */
        sq = row * hptr->nx + col;
        sq2 = (row-1) * hptr->nx + col;
        if ( (mask[sq]==1) && (mask[sq2]==1) ) {
          phi_max = phi[sq];
	  phi_min = phi[sq2];
	  dist_min = dist[sq] > dist[sq2] ? dist[sq2]: dist[sq];
	  set_node(mask, phi, dist, nsq, phi_min, phi_max, dist_min, 1);
        }
	 /* western quadrant */
        sq = row * hptr->nx + col2;
        sq2 = (row-1) * hptr->nx + col2;
        if ( (mask[sq]==1) && (mask[sq2]==1) ) {
          phi_max = phi[sq2];
	  phi_min = phi[sq];
	  
	  if (phi_max > 0 && phi_min < 0)  /* special case for phi = +/-180 */
	      phi_max = -180.0;
	      
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_node(mask, phi, dist, nsq, phi_min, phi_max, dist_min, 1);
        }
	
	--row;
     } /* end while */
     
   } /* end for rad */
   
   /* quadrants north and south of center */
   
   for (rad=1; rad <= ymid; ++rad) {
     col = xmid - rad;
     col2 = xmid + rad;
     row = ymid + rad;
     row2 = ymid - rad;
     if (col < 0) {
	col = 0;
        col2 = hptr->nx - 1;
     }
     if (row2 < 0) 
         row2 = 0;
     if (row >= hptr->ny )
         row = hptr->ny -1;
	 
        /* work left to right */
     while ( col < col2) {
        /* north quadrant */
        sq = row * hptr->nx + col;
        sq2 = sq+1;
        if ( (mask[sq]==1) && (mask[sq2]==1) ) {
          phi_max = phi[sq];
	  phi_min = phi[sq2];
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_node(mask, phi, dist, nsq, phi_min, phi_max, dist_min, 1);
        }
        /* south quadrant */
        sq = row2 * hptr->nx + col;
        sq2 = sq + 1;
        if ( (mask[sq]==1) && (mask[sq2]==1) ) {
          phi_max = phi[sq2];
	  phi_min = phi[sq];
	  dist_min = dist[sq]> dist[sq2]? dist[sq2]: dist[sq];
	  set_node(mask, phi, dist, nsq, phi_min, phi_max, dist_min, 1);	  
        }
	++col;
     } /* end while */
     
   } /* end for rad */
   
   
   free(latv);
   free(lonv);
   free(phi);
   free(dist);
   
   return(mask);
   
} /* end set_mask2d */

/*****************************************************************************/
int *set_mask2d_dist(struct GRID_INFO *hptr, double xlen, double ylen, double theLat, double theLon)

/*  Returns array of values set to 0 or 1 for nodes that are beyond the 
    distance ellipse specified by xlen, ylen).  
    
    Returns a pointer to mask or NULL if an error occurs.
*/
{
   int *mask; 
   int nsq, sq, row, col, xmid, ymid, error;
   double offset, xlen2, ylen2; 
   double xdist, ydist, dist, phi;
   double latv, *lonv;
   
   xlen2 = xlen*xlen;
   ylen2 = ylen*ylen;
   
   nsq = hptr->nx * hptr->ny;
   offset = hptr->node_offset * 0.5; /* 0 if gridnode, 0.5 for pixel grids*/

   lonv = (double *) calloc(hptr->nx, sizeof(double));
   for (col = 0; col < hptr->nx; ++col) {
          lonv[col] = hptr->x_min + (col + offset) * hptr->x_inc;
   }
   mask = (int *)calloc(nsq, sizeof(int));
      
   for (row = 0; row < hptr->ny; ++row) {
      latv = hptr->y_min + (row + offset) * hptr->y_inc;
      for (col = 0; col < hptr->nx; ++col) {
         sq = row * hptr->nx + col;
	 dist = distance_c(latv,lonv[col], theLat, theLon, TRUE, &phi);
	 phi *= RADperDEG;
	 
	 xdist = dist * sin(phi);
         ydist = dist * cos(phi);
	 if (xdist*xdist/xlen2 + ydist*ydist/ylen2 > 1.0)
	     mask[sq] = 1;
      }
   }
   
   free(lonv); 
      
   return(mask);     
      
} /* end set_mask2d_dist() */
/*****************************************************************************/
int getsurf2d(double xval, int reflev, double **xprofiles, double **yprofiles, double **pprofiles, struct GRID_INFO *xinfo, int *nlevs, double *pbot, double *ysurf)
/* Extracts y values from yprofiles along surface corresponding to xval in xprofiles
   at each point in a lat/lon grid.  Returns the number of pts with actual values.  

     xval:  surface value to project y onto
     reflev:  for sigma surfaces (set to -1 for pressure surface)
     xprofiles, yprofiles, pprofiles:  dimension[xinfo->nx * xinfo->ny][xinfo->nz]
     xinfo:  information about the grid
     nlevs:   number of levels in each profile
     pbot:  pressure at seafloor
     ysurf: array of projected values:  dimension  [xinfo->nx * xinfo->ny]
*/

{
  int isq, ibotm, nsq, nvals;
  double yval, pval;
  
  nsq = xinfo->nx * xinfo->ny;
  nvals = 0;
  
  for (isq = 0; isq < nsq; ++isq) {
  
      ibotm = nlevs[isq] - 1;
      ysurf[isq] = HBEMPTY;
  
      
      if (nlevs[isq] <= 0 ) {
         if (pbot[isq] <= 1)
      	     ysurf[isq] = HBMASK;
      } 
      else if (xval > (xprofiles[isq][ibotm])) {   /* xval deeper than deepest val  */
         if (pprofiles[isq][ibotm] > (pbot[isq] - 50) )  /*  within 50 m of seafloor? */
	     ysurf[isq] = HBMASK;
      }
      else if (xval < xprofiles[isq][0]) {  /* xval not in profile */
          ;
      }
      else {
  
        yval = hb_linterp(xval,xprofiles[isq], yprofiles[isq], nlevs[isq]);
        if (yval > -9998) {
             ysurf[isq] = yval;
	     ++nvals;
	     	   
             if (reflev >= 0) {  /* for density surfaces */
	       pval = hb_linterp(xval, xprofiles[isq], pprofiles[isq],nlevs[isq]);
	       if (pval < (reflev-1500)) {
	          ysurf[isq] = HBEMPTY;
		  --nvals;
		}  
	       if ((pval > (reflev+1500)) && (reflev < 100)) {
	          ysurf[isq] = HBEMPTY;
		  --nvals;
	       }   
	     }
	}
      } /* end else */
  
  
  } /* end for */
  
  
  return(nvals);

} /* end getsurf2d() */
/****************************************************************************/
void mask_deeper_levs(int thislev, struct NC_INFO *info_ptr, struct NC_PARMS *parms_ptr)
   /* starting at thislev, sets all values at deeper levels to the mask value  */
{
  int ilev;

  
  for (ilev = thislev; ilev < info_ptr->nz; ++ilev) {
     parms_ptr->Lx[ilev] = info_ptr->mask_value;
     parms_ptr->Ly[ilev] = info_ptr->mask_value;
     parms_ptr->dTdx[ilev] = info_ptr->mask_value;
     parms_ptr->dTdy[ilev] = info_ptr->mask_value;
     parms_ptr->dPdx[ilev] = info_ptr->mask_value;
     parms_ptr->dPdy[ilev] = info_ptr->mask_value;
     parms_ptr->Tparm0[ilev] = info_ptr->mask_value;
     parms_ptr->Tparm1[ilev] = info_ptr->mask_value;
     parms_ptr->Tparm2[ilev] = info_ptr->mask_value;
     parms_ptr->Sparm0[ilev] = info_ptr->mask_value;
     parms_ptr->Sparm1[ilev] = info_ptr->mask_value;
     parms_ptr->Sparm2[ilev] = info_ptr->mask_value;
     parms_ptr->Pparm0[ilev] = info_ptr->mask_value;
     parms_ptr->Pparm1[ilev] = info_ptr->mask_value;
     parms_ptr->Pparm2[ilev] = info_ptr->mask_value;
     ++nmasked;
  }
  
  return;
} /* end mask_deeper_levs() */

/*****************************************************************************/
float get_mean(struct GRID_INFO *ginfo, double *xsurf, int nsq, float xdist, float ydist)
/* Estimates mean  at central node of gridded surface from surrounding values.  */
{
  int theRow, theCol;
  double *weights;
  double xbar;

  weights = (double *) calloc(nsq, sizeof(double));
  
  theCol = ginfo->nx / 2;
  theRow = ginfo->ny / 2; 
  
  get_weights_c(weights, ginfo, 0, (double)xdist, (double)ydist, theCol, theRow);
  
  xbar = weighted_mean(xsurf, weights, (double) HBEMPTY, (double) HBMASK, nsq);
  
  
  free(weights);
  return( (float) xbar);
  
} /* end get_mean() */

/****************************************************************************/
void grid_search(int imodel, struct VARIOGRAM *vptr, double *range1, double *range2, float *parms_out, int snflag)
   /*  Uses Weighted Least Squares to find parameter combination that 
       minimizes the difference between 
       observed and modeled V(h). Ranges are 3 element array giving min,max,incr for sill and decorr scale (p2) 
         Start by solving for sill (p0+p1) and p2 parameters, then break sill into p0 (nugget) and 
	 p1 (sill-nugget) parameters.
	  
	  semi-variogram = p0 + p1[1-exp(-h/p2)]
	  
	  parms_out[0] = p0 = nugget
	  parms_out[1] = p1 = sill - nugget
	  parms_out[2] = p2 = range
	  
	If snflag is set, on input parms_out[0] contains the ratio p1/p0. Use
	this to determine nugget and sill-nugget values
	  
	  */
{
  int np1, np2, i, j, k, ipwr;
  double *parms1, *parms2, *vh, *wght;
  double e_sum, e_min, w, w_sum, diff, parms_in[3];
  double sill, p0, p1, p2;
   
   /* generate vector of values corresponding to parms 1 & 2 based on ranges */
   	     
   np1 = NINT((range1[1] - range1[0]) / range1[2]);
   parms1 = (double *) calloc(np1, sizeof(double));
	     
   for (k = 0; k < np1; ++k)  
        parms1[k] = range1[0] + range1[2] * k;
	     
   np2 = NINT((range2[1] - range2[0]) / range2[2]);
   parms2 = (double *) calloc(np2, sizeof(double));
   
	     
   for (k = 0; k < np2; ++k)  
       parms2[k] = range2[0] + range2[2] * k;
       
   
  sill = NCPARMS_FILL;
  p0 = NCPARMS_FILL;
  p1 = NCPARMS_FILL;
  p2 = NCPARMS_FILL;
   
  /* assign weights based on number of pts in each bin. */
  
  wght = (double *) calloc(vptr->nlags, sizeof(double));
  for   (i = 0; i < vptr->nlags; ++i) {
    if (vptr->Nh[i] > 0)
        wght[i] = log10((double)vptr->Nh[i]);
  }  
      
  /* start with 2 parameter system:  sill (p0+p1) and range (p2) */ 
	     
  e_min = 9e34;
  parms_in[0] = 0;  /* set to zero nugget for now */
  
  for (k = 0; k < np2; ++k) {
    parms_in[2] = parms2[k];
    
    for (j = 0; j < np1; ++j) {
    
        parms_in[1] = parms1[j];
        
	e_sum = 0;
	w_sum = 0;
	vh = vgram_vals(imodel, vptr->h, vptr->nlags, parms_in);
        for (i = 0; i < vptr->nlags; ++i) {
	   
	   if (vptr->Nh[i] > 0 && vh[i] > 0) {
	      diff = vptr->Vh[i] - vh[i];
	      w = wght[i]/ (vh[i] * vh[i]);
	      e_sum += w * diff * diff;
	      w_sum +=  w ;
	   }   
	}  /* end for i */
	
	e_sum = e_sum / w_sum;
	if (e_sum < e_min) {
	   e_min = e_sum;
	   sill =  parms_in[1];
	   p2 =  parms_in[2];
	}
	
	free(vh);   
	   	 
    } /*end for j*/
  } /*end for k */
  

  if (snflag)  {
    /* compute nugget from signal-to-noise ratio, supplied at parms_out[0] */
    p0 = sill / (1 + parms_out[0]);
  }
  else {
  
 /* Estimate nugget value (noise) from average of observed covariances for 
    lags less than 60 km,  including a 0 at lag=0, but omitting values > sill */
 
     e_sum=0;
     i = 0;
     j = 1;
     while (i < vptr->nlags && vptr->h[i] < 60.0 ) {
       if (vptr->Vh[i] <= sill){
           e_sum += vptr->Vh[i];
	   ++j;
        }
	++i;
     }
     p0 = 0;
     if (j > 0) 
       p0 =  e_sum / j;
 
     if (p0 > sill)    /* not viable,assume a S/N ratio of 2:1 */
        p0 = sill * 0.33;
  }
  
   /* Now p0 and p1 are determined.
      Repeat grid search for optimal decorrelation scale, p2 */ 
 
    p1 = sill - p0;
    e_min = 9e34;
    parms_in[0] = p0;  /* set to nugget value  */
    parms_in[1] = p1;  /* sill - nugget */
  
  for (k = 0; k < np2; ++k) {
       
        parms_in[2] = parms2[k];
	e_sum = 0;
	w_sum = 0;
	vh = vgram_vals(imodel, vptr->h, vptr->nlags, parms_in);
        for (i = 0; i < vptr->nlags; ++i) {
	   
	   if (vptr->Nh[i] > 0 && vh[i] > 0) {
	      diff = vptr->Vh[i] - vh[i];
	      w = wght[i]/ (vh[i] * vh[i]);
	      e_sum += w * diff * diff;
	      w_sum +=  w ;
	   }   
	}  /* end for i */
	
	e_sum = e_sum / w_sum;
	if (e_sum < e_min) {
	   e_min = e_sum;
	   p2 =  parms_in[2];
	}
	free(vh);   
  } /*end for k */
  
  free(parms1);
  free(parms2);
  free(wght);
   
  parms_out[0] = (float) p0;
  parms_out[1] = (float) p1;
  parms_out[2] = (float) p2;
  return;
}  /* end grid_search() */

/*********************************************/
void get_gradients(double *prop, struct GRID_INFO *ginfo, double theLat, double theLon, float missing_val, int irad, float *dpdx, float *dpdy)
   /* Evaluates x and y gradients of property at location theLat, theLon over the radius specified.   
      Returns average gradient per 100 km computed from all pairs of points in the 
      window  or missing_val if gradient cannot be computed. 
      Testmask and testempty are global variables. */
   
{ 
   int error, ii, jj, kk, theSq, sq0, row0, col0, npts, stop;
   int  kilo, nwin;
   double sum, dist, phi;
   double *lat, *lon, *p;
   
   kilo = 1;
   nwin = irad * 2 + 1;
   lat = (double *) calloc(nwin, sizeof(double));
   lon = (double *) calloc(nwin, sizeof(double));
   p = (double *) calloc(nwin, sizeof(double));
   
   error = xy2ij(ginfo, theLon, theLat, &col0, &row0);
   if (error) {
      fprintf(stderr, "FATAL ERROR: lat/lon outside grid bounds in get_gradients()\n\n");
      exit(1);
   }
   sq0 = ginfo->nx * row0 + col0;
   
   /* load points in zonal direction */
   
   npts = 0;
   ii = 0;
   stop = 0;
   while (ii <= irad && !stop) {
       theSq = sq0 + ii;
       if (prop[theSq] > testmask ) {
           stop = 1;
       }   
       else {
           if (prop[theSq] > testempty) {
	      p[npts] = prop[theSq];
	      ij2xy(ginfo, col0+ii, row0, &lon[npts], &lat[npts]);
	      ++npts;
	   }
	   ++ii;
       }
   }
   
   if (!stop) {
       ii = 0 ;
       stop = 0;
       while (++ii <= irad && !stop) {
          theSq = sq0 - ii;
          if (prop[theSq] > testmask) {
	     stop = 1;
	  }
	  else {
	     if (prop[theSq] > testempty) {
		p[npts] = prop[theSq];
		ij2xy(ginfo, col0-ii, row0, &lon[npts], &lat[npts]);
		++npts;
	     }
	  }   
       }
   }
   
   *dpdx = missing_val;
   sum = 0;
   kk = 0;
   if (npts > 1) {
       for (jj = 0; jj < npts; ++jj) {
          for (ii = jj+1; ii < npts; ++ii) {
	     dist = distance_c( lat[jj], lon[jj], lat[ii], lon[ii], kilo, &phi);
	     sum += ( ABS(p[jj] - p[ii]) / dist);
	     ++kk;
	  }
       }
      *dpdx = (float) (sum * 100.0 / (double)kk );
   }
   
   /* load points in meridional direction */
   
    ii = 0 ;
    stop = 0;
    npts = 0;
    while (ii < irad && !stop) {
        theSq = sq0 + ii * ginfo->nx;
       if (prop[theSq] > testmask ) {
           stop = 1;
       }   
       else {
           if (prop[theSq] > testempty) {
	      p[npts] = prop[theSq];
	      ij2xy(ginfo, col0, row0+ii, &lon[npts], &lat[npts]);
	      ++npts;
	   }
	   ++ii;
       }	   
   }
   
   if (!stop) {
       ii = 0 ;
       stop = 0;
       while (++ii < irad && !stop) {
          theSq = sq0 - ii * ginfo->nx;
          if (prop[theSq] > testmask) {
	     stop = 1;
	  }
	  else {
	     if (prop[theSq] > testempty) {
		p[npts] = prop[theSq];
		ij2xy(ginfo, col0, row0-ii, &lon[npts], &lat[npts]);
		++npts;
	     }
	  }   
       }
   }
   
   *dpdy = missing_val;
   sum = 0;
   kk = 0;
   if (npts > 1) {
       for (jj = 0; jj < npts; ++jj) {
          for (ii = jj+1; ii < npts; ++ii) {
	     dist = distance_c(lat[jj],lon[jj], lat[ii], lon[ii], kilo, &phi);
	     sum += ( ABS(p[jj] - p[ii]) / dist);
	     ++kk;
	  }
       }
      *dpdy = (float) (sum  * 100.0/ (double)kk );
   }

   return;

} /* end get_gradients() */   
/*************************************************/
/*************************************************/
void fill_missing_nodes(int ncid, struct NC_INFO *ncinfo_ptr, struct NC_PARMS *parms_ptr, struct GRID_INFO *ginfo_ptr, double maxdist)
  /* On entry, the structures have all been allocated and filled according to the work to be done. 
     Here we visit each gridnode in the region specified by ginfo_ptr, read in the parameters and determine which nodes are flagged as EMPTY, and attempt to fill them using a near neighbor interpolation.
     The global counters, nempty and nfilled are updated. */
     
{
   int  row, col, nlevs, ilev, jlev, isq, iwsq, ismsq, nsqout, nwsq, nsmsq; 
   int ncrow, nccol;
   int xpole, wrow, wcol, error, masked;
   int *mask_ptr;
   size_t start[3], count[3];
   struct GRID_INFO parms_info, work_grid, small_grid;
   float theMean;
   double theLat, theLon;
   double dummy, latmin, lonmin, latmax, lonmax;
   double **Lx, **Ly, **dTdx, **dTdy;
   double **dPdx, **dPdy, **Tp0, **Tp1, **Tp2;
   double **Sp0,**Sp1, **Sp2, **Pp0,**Pp1, **Pp2;
   
   /* transfer info about the parameter file to a GRID_INFO structure */
   
   parms_info.nx = ncinfo_ptr->nx;
   parms_info.ny = ncinfo_ptr->ny;
   parms_info.x_min = ncinfo_ptr->lons[0];
   parms_info.x_max = ncinfo_ptr->lons[ncinfo_ptr->nx-1];
   parms_info.y_min = ncinfo_ptr->lats[0];
   parms_info.y_max = ncinfo_ptr->lats[ncinfo_ptr->ny-1];
   parms_info.x_inc = ncinfo_ptr->xincr;
   parms_info.y_inc = ncinfo_ptr->yincr;
   parms_info.node_offset = 1;
   parms_info.xgreenwich = 0;
   parms_info.lon0to360 = 1;
   
   count[0] = 1;    /* used by ncparms_write_var() */
   count[1] = 1;
   count[2] = 1;
   
   /* determine boundaries of work grid */
   
   nlevs = ncinfo_ptr->nz;
   nsqout = ginfo_ptr->nx * ginfo_ptr->ny;
   
   isq = 0;
   xpole = get_bounds(ginfo_ptr, isq, maxdist, maxdist, &latmin, &lonmin, &dummy, &dummy);
   if (xpole)
      latmin = -90;
   isq = nsqout - 1;
   xpole = get_bounds(ginfo_ptr, isq, maxdist, maxdist, &dummy, &dummy, &latmax, &lonmax);
   if (xpole)
      latmax = 90;
   xpole = 0;   
   nwsq = set_work_grid(ginfo_ptr, &work_grid, &latmin, &lonmin, &latmax, &lonmax, xpole);   
   
   /* allocate memory for work grids */
   
   Lx = (double **) calloc(nlevs, sizeof(double *));
   Ly = (double **) calloc(nlevs, sizeof(double *));
   dTdx = (double **) calloc(nlevs, sizeof(double *));
   dTdy = (double **) calloc(nlevs, sizeof(double *));
   dPdx = (double **) calloc(nlevs, sizeof(double *));
   dPdy = (double **) calloc(nlevs, sizeof(double *));
   Tp0 = (double **) calloc(nlevs, sizeof(double *));
   Tp1 = (double **) calloc(nlevs, sizeof(double *));
   Tp2 = (double **) calloc(nlevs, sizeof(double *));
   Sp0 = (double **) calloc(nlevs, sizeof(double *));
   Sp1 = (double **) calloc(nlevs, sizeof(double *));
   Sp2 = (double **) calloc(nlevs, sizeof(double *));
   Pp0 = (double **) calloc(nlevs, sizeof(double *));
   Pp1 = (double **) calloc(nlevs, sizeof(double *));
   Pp2 = (double **) calloc(nlevs, sizeof(double *));
   
   for (ilev = 0; ilev < nlevs; ++ilev) {
      Lx[ilev] = (double *) calloc(nwsq, sizeof(double));
      Ly[ilev] = (double *) calloc(nwsq, sizeof(double));
      dTdx[ilev] = (double *) calloc(nwsq, sizeof(double));
      dTdy[ilev] = (double *) calloc(nwsq, sizeof(double));
      dPdx[ilev] = (double *) calloc(nwsq, sizeof(double));
      dPdy[ilev] = (double *) calloc(nwsq, sizeof(double));
      Pp0[ilev] = (double *) calloc(nwsq, sizeof(double));
      Pp1[ilev] = (double *) calloc(nwsq, sizeof(double));
      Pp2[ilev] = (double *) calloc(nwsq, sizeof(double));
      Tp0[ilev] = (double *) calloc(nwsq, sizeof(double));
      Tp1[ilev] = (double *) calloc(nwsq, sizeof(double));
      Tp2[ilev] = (double *) calloc(nwsq, sizeof(double));
      Sp0[ilev] = (double *) calloc(nwsq, sizeof(double));
      Sp1[ilev] = (double *) calloc(nwsq, sizeof(double));
      Sp2[ilev] = (double *) calloc(nwsq, sizeof(double));
   }
   
   /* Read in parameters for each square in the work grid and 
      sort into the work arrays defined above */
      
   fprintf(stderr,"Reading parameters file...\n");
   for (row = 0; row < work_grid.ny; ++row) {
      for (col = 0; col < work_grid.nx; ++col) {
         isq = row * work_grid.nx + col;
         error = ij2xy(&work_grid, col, row, &theLon, &theLat);
         error = ncparms_read(ncid, (float)theLat, (float)theLon, ncinfo_ptr, parms_ptr);
	 if (!error) {
	    for (ilev = 0; ilev < nlevs; ++ilev) {
	       Lx[ilev][isq] = (double)parms_ptr->Lx[ilev];
	       Ly[ilev][isq] = (double)parms_ptr->Ly[ilev];
	       dTdx[ilev][isq] = (double)parms_ptr->dTdx[ilev];
	       dTdy[ilev][isq] = (double)parms_ptr->dTdy[ilev];
	       dPdx[ilev][isq] = (double)parms_ptr->dPdx[ilev];
	       dPdy[ilev][isq] = (double)parms_ptr->dPdy[ilev];
	       Pp0[ilev][isq] = (double)parms_ptr->Pparm0[ilev];
	       Pp1[ilev][isq] = (double)parms_ptr->Pparm1[ilev];
	       Pp2[ilev][isq] = (double)parms_ptr->Pparm2[ilev];
	       Tp0[ilev][isq] = (double)parms_ptr->Tparm0[ilev];
	       Tp1[ilev][isq] = (double)parms_ptr->Tparm1[ilev];
	       Tp2[ilev][isq] = (double)parms_ptr->Tparm2[ilev];
	       Sp0[ilev][isq] = (double)parms_ptr->Sparm0[ilev];
	       Sp1[ilev][isq] = (double)parms_ptr->Sparm1[ilev];
	       Sp2[ilev][isq] = (double)parms_ptr->Sparm2[ilev]; 
	    }
	 }
      }
   }
   
   /* Visit each square in the output grid and check for empty nodes */
   
   fprintf(stderr,"Now filling empty nodes...\n");
   
   nfilled = 0;
   nempty = 0;
   
   for (row = 0; row < ginfo_ptr->ny; ++row) {
      for (col = 0; col < ginfo_ptr->nx; ++col) {
      
         error = ij2xy(ginfo_ptr, col, row, &theLon, &theLat);
	 error = xy2ij(&work_grid, theLon, theLat, &wcol, &wrow);
	 iwsq = wrow * work_grid.nx + wcol;
	 xpole = get_bounds(&work_grid, iwsq, maxdist, maxdist, &latmin, &lonmin, &latmax, &lonmax); 
	 nsmsq = set_work_grid(&work_grid, &small_grid, &latmin, &lonmin, &latmax, &lonmax, xpole);
	 
	 ncparms_rowcol(ncinfo_ptr, theLat, theLon, &ncrow, &nccol);
	 start[0] = ncrow;
	 start[1] = nccol;
	 masked = 0;
	 ilev = -1;
	 while ( ++ilev < nlevs && !masked) {
	 
	    if (masked = (Tp1[ilev][iwsq] > testmask) ) {
	       for (jlev = ilev; jlev < nlevs; ++jlev) {
	           start[2] = jlev;
		   theMean = (float)NCPARMS_MASK;
		   error = ncparms_write_var(ncid,"Tparm1", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Tparm2", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Tparm0", start, count, &theMean);
	           Tp2[jlev][iwsq] = NCPARMS_MASK;
	           Tp0[jlev][iwsq] = NCPARMS_MASK;
	           Tp1[jlev][iwsq] = NCPARMS_MASK;
		   error = ncparms_write_var(ncid,"Sparm1", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Sparm2", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Sparm0", start, count, &theMean);
	           Sp2[jlev][iwsq] = NCPARMS_MASK;
	           Sp0[jlev][iwsq] = NCPARMS_MASK;
	           Sp1[jlev][iwsq] = NCPARMS_MASK;
		   error = ncparms_write_var(ncid,"Pparm1", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Pparm2", start, count, &theMean);
		   error = ncparms_write_var(ncid,"Pparm0", start, count, &theMean);
	           Pp2[jlev][iwsq] = NCPARMS_MASK;
	           Pp0[jlev][iwsq] = NCPARMS_MASK;
	           Pp1[jlev][iwsq] = NCPARMS_MASK;
		}
	        continue;
	    }
	       
	    start[2] = ilev;
	    if (Tp1[ilev][iwsq] < testempty) {
	         theMean = (float) fill_it(iwsq, theLat, theLon, Tp1[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		++nempty;
		
		if (theMean < testempty) {
		   if ((ilev > 0) && (Tp1[ilev-1][iwsq] > testempty) && (Tp1[ilev-1][iwsq] < testmask))
		       theMean = Tp1[ilev-1][iwsq];
		       
		}
		if (theMean > testempty) {
                   Tp1[ilev][iwsq] = theMean;
		   ++nfilled;
		   error = ncparms_write_var(ncid,"Tparm1", start, count, &theMean);
/*		   
fprintf(stderr,"Filled: %.0f m %8.3lf %8.3lf %.6f {%d,%d,%d}\n", parms_ptr->depths[ilev], theLat, theLon, theMean,start[0],start[1],start[2]); */
		   
		} 
		else {
		   fprintf(stderr,"Still empty: %.0f m %8.3lf %8.3lf \n", parms_ptr->depths[ilev], theLat, theLon);
		} 
	    }
	    
	    if (Tp2[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, Tp2[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Tp2[ilev-1][iwsq] > testempty && Tp2[ilev-1][iwsq] < testmask)
		       theMean = Tp2[ilev-1][iwsq];
		}
		if (theMean > testempty) {
		   Tp2[ilev][iwsq] = theMean;
		   error = ncparms_write_var(ncid,"Tparm2", start, count, &theMean);
		}  
	    }
	    if (Tp0[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, Tp0[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Tp0[ilev-1][iwsq] > testempty && Tp0[ilev-1][iwsq] < testmask)
		       theMean = Tp0[ilev-1][iwsq];
		}
		if (theMean > testempty) {
                   Tp0[ilev][iwsq] = theMean;		   
		   error = ncparms_write_var(ncid,"Tparm0", start, count, &theMean);
		}  
	    }
	    
	    if (ilev > ncinfo_ptr->nz_seas) {
	       
	       if (Pp0[ilev][iwsq] < testempty) {
	           theMean = (float) fill_it(iwsq, theLat, theLon, Pp0[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Pp0[ilev-1][iwsq] > testempty && Pp0[ilev-1][iwsq] < testmask)
		       theMean = Pp0[ilev-1][iwsq];
		}
		if (theMean > testempty) {
		      Pp0[ilev][iwsq] = theMean;
		      error = ncparms_write_var(ncid,"Pparm0", start, count, &theMean);
		}  
		   
	       }  
	       if (Pp1[ilev][iwsq] < testempty){
	           theMean = (float) fill_it(iwsq, theLat, theLon, Pp1[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
	       if (theMean < testempty) {
		   if (ilev > 0 && Pp1[ilev-1][iwsq] > testempty && Pp1[ilev-1][iwsq] < testmask)
		       theMean = Pp1[ilev-1][iwsq];
	       }
		   if (theMean > testempty) {
		      Pp1[ilev][iwsq] = theMean;
		      error = ncparms_write_var(ncid,"Pparm1", start, count, &theMean);
		   }  
	       }
	       
	       if (Pp2[ilev][iwsq] < testempty) {
	           theMean = (float) fill_it(iwsq, theLat, theLon, Pp2[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Pp2[ilev-1][iwsq] > testempty && Pp2[ilev-1][iwsq] < testmask)
		       theMean = Pp2[ilev-1][iwsq];
		}
		   if (theMean > testempty) {
		      Pp2[ilev][iwsq] = theMean;
		      error = ncparms_write_var(ncid,"Pparm2", start, count, &theMean);
		   }  
		   
	       }  
	    }
	    
	    
	    if (Sp0[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, Sp0[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Sp0[ilev-1][iwsq] > testempty && Sp0[ilev-1][iwsq] < testmask)
		       theMean = Sp0[ilev-1][iwsq];
		}
		   if (theMean > testempty) {
		      Sp0[ilev][iwsq]= theMean;
		      error = ncparms_write_var(ncid,"Sparm0", start, count, &theMean);
		   }  
	    }
	    
	    if (Sp1[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, Sp1[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Sp1[ilev-1][iwsq] > testempty && Sp1[ilev-1][iwsq] < testmask)
		       theMean = Sp1[ilev-1][iwsq];
		}
		   if (theMean > testempty) {
		      Sp1[ilev][iwsq] = theMean;
		      error = ncparms_write_var(ncid,"Sparm1", start, count, &theMean);
		   }  
	    }
	    
	    if (Sp2[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, Sp2[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean < testempty) {
		   if (ilev > 0 && Sp2[ilev-1][iwsq] > testempty && Sp2[ilev-1][iwsq] < testmask)
		       theMean = Sp2[ilev-1][iwsq];
		}
		   if (theMean > testempty) {
		      Sp2[ilev][iwsq] = theMean;
		      error = ncparms_write_var(ncid,"Sparm2", start, count, &theMean);
		   }  
	    }
	    
	    if (dTdx[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, dTdx[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean > testempty) {
		      error = ncparms_write_var(ncid,"dTdx", start, count, &theMean);
		}  		
	    }
	    if (dTdy[ilev][iwsq] < testempty) {
	        theMean = (float) fill_it(iwsq, theLat, theLon, dTdy[ilev], &work_grid, &small_grid, latmin, lonmin, latmax, lonmax, maxdist);
		if (theMean > testempty) {
		      error = ncparms_write_var(ncid,"dTdy", start, count, &theMean);
		}  		
	    }
	    
	    if (Lx[ilev][iwsq]< testempty) {
	       theMean = (float) maxdist;
	       error = ncparms_write_var(ncid,"Xlength", start, count, &theMean);
	    } 
	    if (Ly[ilev][iwsq]< testempty) {
	       theMean = (float) maxdist;
	       error = ncparms_write_var(ncid,"Ylength", start, count, &theMean);
	    } 
	    
	 
	 } /* end while */
	 
      } /* end for col */
   }  /* end for row */
   
/*******************************************/
   
   for (ilev = 0; ilev < nlevs; ++ilev) {
      free(Lx[ilev]);
      free(Ly[ilev]);
      free(dTdx[ilev]);
      free(dTdy[ilev]);
      free(dPdx[ilev]);
      free(dPdy[ilev]);
      free(Pp0[ilev]);
      free(Pp1[ilev]);
      free(Pp2[ilev]);
      free(Tp0[ilev]);
      free(Tp1[ilev]);
      free(Tp2[ilev]);
      free(Sp0[ilev]);
      free(Sp1[ilev]);
      free(Sp2[ilev]);
   }
   free(Lx );
   free(Ly );
   free(dTdx );
   free(dTdy );
   free(dPdx );
   free(dPdy );
   free(Tp0 );
   free(Tp1 );
   free(Tp2 );
   free(Sp0 );
   free(Sp1 );
   free(Sp2 );
   free(Pp0 );
   free(Pp1 );
   free(Pp2 );

   return;

} /* end fill_missing_nodes */
/*************************************************/
double fill_it(int theSq, double theLat, double theLon, double *grid1, struct GRID_INFO *g1info, struct GRID_INFO *g2info, double lat0, double lon0, double lat1, double lon1, double maxdist)
/*  extracts values from grid1 between lat0,lon0 and lat1,lon1 and puts them into grid2. 
    Performs a near-neighbor search of grid2 to estimate a value at theSq in grid1  */
{
   int nsq2, isq2, isq1, error, icount;
   int row, col, row2, col2, strow, endrow, stcol, endcol;
   double theMean, dummy, x;
   double *weights, *grid2;
   double *latv, *lonv;
   
   nsq2 = g2info->nx * g2info->ny;
   grid2 = (double *) calloc(nsq2, sizeof(double));
   weights = (double *) calloc(nsq2, sizeof(double));
   latv = (double *) calloc(g2info->ny, sizeof(double));
   lonv = (double *) calloc(g2info->nx, sizeof(double));
   
   error = xy2ij(g1info, lon0, lat0, &stcol, &strow);
   error += xy2ij(g1info, lon1, lat1, &endcol, &endrow);
   
   if (error) {
      fprintf(stderr,"FATAL ERROR: start/end position out of bounds in fill_it()..\n");
      fprintf(stderr,"lat0,lon0: %.3lf, %.3lf   lat1,lon1: %.3lf,%.3lf\n", lat0, lon0,lat1,lon1);
      exit(1);
   }
   
   /* set up lat, lon vectors */
   
   
   for (col = 0; col < g2info->nx; ++col)
         lonv[col] = g2info->x_min + col * g2info->x_inc;
   
   for (row = 0; row < g2info->ny; ++row)  
      latv[row] = g2info->y_min + row * g2info->y_inc;
   
   /* prefill grid2 with missing values */
   
   for (isq2 = 0; isq2 < nsq2; ++isq2) 
       grid2[isq2] =  (double) NCPARMS_FILL;
   
    /* transfer values from grid1 to grid2 */
  
   for (row = strow; row < endrow; ++row) {
      for (col = stcol; col < endcol; ++col) {
         isq1 = row * g1info->nx + col;
	 xy2ij(g2info, lonv[col], latv[row], &col2, &row2);
	 isq2 = row2 * g2info->nx + col2;
	 grid2[isq2] = grid1[isq1];
      }
   }
   
   error = xy2ij(g2info, theLon, theLat, &col2, &row2);
   get_weights_c(weights, g2info, 0, maxdist, maxdist, col2, row2);
   nn_interp2d(grid2, weights, (double) NCPARMS_FILL, (double) NCPARMS_MASK, col2, row2, g2info, &theMean, &dummy, &dummy, &icount);
   
   free(lonv);
   free(latv);
   free(grid2);
   free(weights);
   
   return(theMean);
   
}  /* end fill_it() */
