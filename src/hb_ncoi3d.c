/* hb_ncoi3d.c
................................................................................
                          *******  HydroBase 3 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Feb 2012


* Implements an optimal interpolation algorithm (aka ordinary kriging and 
  Gauss Markov mapping) based on estimated spatial covariance functions.
  Constructs 3D property fields,including error variance statistics, gridded 
  as a function of latitude,longitude and depth. The grid is user-specified.
  Observed properties are interpolated along isopycnal surfaces below the 
  (locally determined) mixed layer and the resulting property profile is
  interpolated back onto standard depths.
  
  This module requires a netcdf global parameters file (created by
   hb_build_oi_parms) containing variogram parameters and horizontal gradient
   information.  The observations must be pre-processed (sorted into geographic
   bins) with hb_bin3d.  
  
  The gridded properties (mean, sqrt(error variance), number of obs, 
  variogram parameters, x- and y-length scales for the search ellipse) are
  output to a netcdf file.  
  
  hb_ncoi3d -B<bounds> -I<xinc>/<yinc> -Gr<gridded_bin3d_files> [-Fr<first-guess_files>] -O<output_file> -P<proplist> [-V<variogram_file>] [-M<minobs>] [-S<max_search_radius>] [-Z<stddepth_file>]
  
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


#define MIX_LAYER_DEF 0.02    /* default sigma range for mixed layer */
#define SEASONAL_LAYER_DEF 200  /* depth of seasonal layer */
#define NC_EXTENT ".nc"
#define MAXDIST 600  /* in kilometers */
#define DIST_DEFAULT 200 /* if no length scale defined */
#define MINOBS 20   
#define PREFILL 0
#define TOPOBUF 50


/****** global variables ******/

float testmask, testempty;
int print_msg = 1;
int use_fg_in_pinch;
int nfilled, nempty, nmasked, nfg_used;
int tindex = (int)T90;  /* use T90 throughout*/
int is_pr, is_te, is_sa;

int *merid;  /* used by function is_in_range() */
int lon0to360, xgreenwich;  /* for ncparms file */
float mix_layer_delta;
float maxdist;
int nkrig, nempty, nmask, nother;
int minobs;
int nmonths, one_file;
char *parms_name;

char cmonths[NMONTHS][4] = { "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec","clim" };

/* These define the relationship between Tgrad and x-, y-lengthscales
   They are completely arbitrary and can be adjusted  */

float dtmax[8] = {0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 100};
float dist[8] =  {400, 350, 300, 250, 200, 150, 100, 50};

/******  define local data structures  **********/

struct MIXLAYER {
    float density;
    float depth;
    float theta;
    float salinity;
    float derr, terr, serr;
    int nobs;
};

/********************************************/
/* prototypes for locally defined functions */	

void print_usage(char *);
int parse_p_option(char *, int *);
int check_overlap(struct GRID_INFO *, struct GRID_INFO *);
void set_topo_info(struct GRID_INFO *, struct GRID_INFO *, int);   
int subgrid(struct GRID_INFO *, struct GRID_INFO *, double, double, double, double, float);
void cdf_construct( struct GRID_INFO *, int, int *,  char *, char *, char *, int *, int, char **);
void get_sig_profiles(int, struct CDF_HDR *, struct GRID_INFO *, int, float *, short *, struct GRID_INFO *, int, float **, float **, float **, float **, short **, float **, float **, float **, float **, float **, int *);
void get_monthly_profiles(int, struct CDF_HDR *, struct GRID_INFO *, int, float *, int, float **, float **, float **, float **, short **, float **, int *,struct MIXLAYER *);
void get_x_profiles(int, struct CDF_HDR *, struct GRID_INFO *, int, float *, float **, short **, int);
int define_bottom(float *, short );
void define_mixed_layer(double, double *, double *, double *, double *, short *, int, struct MIXLAYER *);
void define_length_scales(struct NC_PARMS *, float *, float *, double *, int, float);
void get_prop_parms(struct NC_PARMS *, int, double *, int, float, float *, float *, float *);
int get_bounds(struct GRID_INFO *, int, float, float, double *, double *, double *, double *);
int set_work_grid(struct GRID_INFO *, struct GRID_INFO *,double *, double *, double *, double *, int);
float find_isopycnal(float, int, struct GRID_INFO *, float **, float **, int *, float, float, float *);
float weighted_stats(struct GRID_INFO *, double *, int, float, float, float *);
int extract_surface(float, int, float, struct GRID_INFO *, float **, float **, short **, float **, int *, float *, double, double, int, float *, float *, struct GRID_INFO *, double **, short **, int *, int);
int getsurf2d_depth(float, struct GRID_INFO *, float **, float **,  short **, int *, float *,double, double, double, double, double *, short *, struct GRID_INFO *);
int getsurf2d(float, int, float, float **, float **, short **, float **, struct GRID_INFO *, int *, float *, double *, short *, struct GRID_INFO *, double, double, double, double);
float check_all_crossings(float, float *, float *, int , float);
short interp_short(float, float *, short *, int );
float interp_float(float, float *, float *, int );
int *set_mask2d(double *, struct GRID_INFO *, double, double);
int *set_mask2d_dist( struct GRID_INFO *, float, float, double, double);
float *compute_residuals(struct GRID_INFO *, double *,int, struct GRID_INFO *, double *, int, short * );
int krig(double, double, double *, short *,int, struct GRID_INFO *, int, float *, float *, float *,  short *);
int set_ml_props(struct MIXLAYER *, struct GRID_INFO *, double, double, int, float *, float *, struct MIXLAYER *);
void x_to_stddepth(struct MIXLAYER *, float *, float *, float *, short *, int, float *, float *, short *, float *, float *, short *, float, double *, int, double);

main(int argc, char **argv)

{
   int ncid_parms;
   int ncid_out[NMONTHS], ncid_fg[NMONTHS], ncid_b3d[NMONTHS];
   int error, n, i, j, nprops, nz, nstdlevs, iclim, ilev, iprop;
   int nsq_out, nsq_fg, nsq_b3d, isq, isq_fg;
   int nwsq_fg, nwsq_b3d;       /*actual # of nodes in work grids */
   int reflev, tbin,  nz_ml;
   int theSq_fg, theSq_b3d;    /* gridsquare corresponding to theLon, theLat */
   int row, col, theRow, theCol, row_out, col_out;
   int bflag, oflag, iflag, popt, do_krig;
   int inrange, imonth, imodel, mask_it;
   int adjust_lengths, fg_avail;
   int tindex_fg, tindex_b3d;
   int nz_seas_fg, nz_seas_b3d, nz_seas_out;
   int *prop_indx, *prop_req, *prop_avail;
   int *nz_fg, *nz_b3d, *nz_m_fg[NMONTHS], *nz_m_b3d[NMONTHS];
   short *seafloor, topoval;
   short **tcnt_b3d, **xcnt_b3d;
   short **tcnt_m_b3d[NMONTHS], **xcnt_m_b3d[NMONTHS];
   short *n_prof, *n_ml, *n_ncout, xcount;
   short *xcnt_surf_b3d;
   char *toponame,  *st;
   char *out_dir, *out_ext, *out_root;
   char *fg_dir, *fg_ext, *fg_root;
   char *b3d_dir, *b3d_ext, *b3d_root;
   char *buf;
   float xmean, xerr, theFG; 
   float xlength, ylength;
   float curr_depth, isopyc_above;
   float oi_params[3];
   float *bdepth_out, *bdepth_fg, *bdepth_b3d;
   float *parms0, *parms1, *parms2, *din;
   float *x_prof, *e_prof, *x_ml, *e_ml;
   float *x_ncout, *e_ncout;
   float *xlen_ncout, *ylen_ncout;
   float *resids, *tsave;
   float **xlen, **ylen;
   float **xlen_seas[NMONTHS], **ylen_seas[NMONTHS];
   float **sigptr_fg, **sigptr_b3d;
   float **p_fg, **d_fg, **t_fg, **s_fg, **x_fg;
   float **sig0_fg, **sig1_fg, **sig2_fg, **sig3_fg, **sig4_fg;
   float **up_fg[NMONTHS],**ud_fg[NMONTHS], **ut_fg[NMONTHS], **us_fg[NMONTHS], **ux_fg[NMONTHS];
   float **usig0_fg[NMONTHS];
   float **p_b3d, **d_b3d, **t_b3d, **s_b3d,**x_b3d;
   float **sig0_b3d, **sig1_b3d, **sig2_b3d, **sig3_b3d, **sig4_b3d;
   float **up_b3d[NMONTHS],**ud_b3d[NMONTHS], **ut_b3d[NMONTHS], **us_b3d[NMONTHS],**ux_b3d[NMONTHS];
   float **usig0_b3d[NMONTHS];
   float **isopyc;
   float **psave[NMONTHS], **ssave[NMONTHS], **isopyc_up[NMONTHS];
   double theLat, theLon;
   double pbml;
   double *topolat, *topolon;
   double *xsurf_fg, *xsurf_b3d;
   struct NC_INFO ncparms_info;
   struct NC_PARMS ncparms_data;
   struct CDF_HDR hdr_fg, hdr_b3d, hdr_out;    /* external grids (netcdf files) */
   struct GRID_INFO fgwork, b3dwork, outgrid;  /* internal grids */
   struct GRID_INFO fg_in, b3d_in;
   struct GRID_INFO fg_orig, b3d_orig;
   struct GRID_INFO topo_info;
   struct MIXLAYER *ml_out[NMONTHS],*ml_b3d[NMONTHS], *ml_fg[NMONTHS];
   struct VARIOGRAM vgram;
   FILE *z_file;
   
   
/* set these default values */
    fg_dir = "";
    fg_ext = NC_EXTENT;
    fg_root = NULL;
    fg_avail = 0;
    b3d_dir = "";
    b3d_ext = NC_EXTENT;
    b3d_root = NULL;
    out_dir = "";
    out_ext = NC_EXTENT;
    out_root = NULL;
    toponame = BATHPATH_C;
    bflag = iflag = oflag = popt = 0;
    one_file = 0;
    merid = NULL;   /* used by function is_in_range() */
    z_file = NULL;
    nfilled = nempty = nmasked = 0;
    testmask = HBMASK * 0.1;
    testempty = HBEMPTY * 0.1;
    maxdist = MAXDIST;
    minobs = MINOBS;
    nmonths = NMONTHS;
    parms_name = OI_PARMS_FILE;
    mix_layer_delta = MIX_LAYER_DEF;
    nkrig = nempty = nmask = nother = nfg_used = 0;
    use_fg_in_pinch = 0;
    tbin = 0;    /* index to only tbin in output files */
/*----------------------------------------*/  
 /* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
/*----------------------------------------*/  
   for (i = 1; i < argc; i++) {
   
      error = argv[i][0] == '-'? 0 : 1;
      if (!error ) {
        switch (argv[i][1]) {
         case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &outgrid.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &outgrid.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &outgrid.y_min) != 1);
          while (*(st++) != '/')
                           ;
          error += (sscanf(st,"%lf", &outgrid.y_max) != 1);
                        	     
	  if (&outgrid.x_min > &outgrid.x_max) {
	    fprintf(stderr,"\nWest bound must be numerically <= east bound");
	    exit(1);
	  }
	  
	  if (&outgrid.y_min > &outgrid.y_max) {
	    fprintf(stderr,"\nNorth bound must be numerically <= south bound");
	    exit(1);
	  }
	  
	  outgrid.lon0to360 = outgrid.x_min >= 0 ? 1: 0;
	  outgrid.xgreenwich = (outgrid.x_min < 0) && (outgrid.x_max >=0) ? 1:0;
          break;
	
          case 'F':    
	     fg_avail = 1;               
             switch (argv[i][2])  {
	        case 'd':
		   fg_dir = &argv[i][3];
		   break;
		case 'e':
		   fg_ext = &argv[i][3];
		   break;
		case 'o':
		   one_file = 1;
		   fg_root = &argv[i][3];
		   fg_ext = "";
		   break;
		case 'r':
		   fg_root = &argv[i][3];
		   break;
		case 's':
		   use_fg_in_pinch = 1;
		   break;
		
		default:
		    error = TRUE;
	            fprintf(stderr,"\nUnrecognized option: %3s\n", argv[i]);
                    fprintf(stderr," FG filenames have form: <dir>/<root>.<month>.<extent> \n");
                    fprintf(stderr,"      Default extent is '.nc' \n");
		    fprintf(stderr,"      <month> will be inserted into filenames automatically.\n");
                    fprintf(stderr,"      Use -Fr to specify dir/root (required)\n");
                    fprintf(stderr,"      Use -Fe to specify file extent (optional)\n");
                    fprintf(stderr,"  ex:  -Fr/data/gridded/Pacific.1deg.FG\n");
                    fprintf(stderr," specifies that 13 files (monthly + clim) called Pacific.1deg.FG.jan.nc, etc.... are to be found in the directory /data/gridded \n");
                    fprintf(stderr,"      Use -Fe to specify a file extent other than '.nc'\n");
                    fprintf(stderr,"      Use -Fo to specify a single input/output file\n");
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
		case 'o':
		   one_file = 1;
		   b3d_root = &argv[i][3];
		   b3d_ext = "";
		   break;
		case 'r':
		   b3d_root = &argv[i][3];
		   break;
		default:
		    error = TRUE;
	            fprintf(stderr,"\nUnrecognized option: %3s\n", argv[i]);
                    fprintf(stderr," bin3d filenames have form: <dir>/<root>.<month>.<extent> \n");
                    fprintf(stderr,"Default extent is '.nc' \n");
		    fprintf(stderr,"<month> will be inserted into filenames automatically.\n");
                    fprintf(stderr,"Use -Gr to specify dir/root (required)\n");
                    fprintf(stderr,"Use -Ge to specify file extent (optional)\n");
                    fprintf(stderr,"  ex:  -Gr/data/gridded/Pacific.quintdeg.bin3d\n");
                    fprintf(stderr," specifies that 13 files (monthly + clim) called Pacific.quintdeg.bin3d.jan.nc, etc.... are to be found in the directory /data/gridded \n");
                    fprintf(stderr,"Use -Go to specify a single input/output file\n");
		    exit(1);
	     }
	  break;

          case 'I':
             iflag = 1;
             error = (sscanf(&argv[i][2],"%lf", &outgrid.x_inc) == 1) ? 0 : 1;
             outgrid.y_inc = outgrid.x_inc;
             st = &argv[i][2];
             while (*(st++) != '\0') {
               if (*st == '/') {
                  ++st;
                error += (sscanf(st,"%lf", &outgrid.y_inc) == 1) ? 0 : 1;
                break;
               }
             }
          break;
	  
          case 'N':   /*minobs */
             error = (sscanf(&argv[i][2],"%d", &minobs) == 1) ? 0 : 1;
          break;
	  
          case 'O':
	     oflag = 1;
             switch (argv[i][2])  {
	        case 'd':
		   out_dir = &argv[i][3];
		   break;
		case 'e':
		   out_ext = &argv[i][3];
		   break;
		case 'o':
		   one_file = 1;
		   out_root = &argv[i][3];
		   out_ext = "";		   
		   break;
		case 'r':
		   out_root = &argv[i][3];
                   break;
		default:
		    error = TRUE;
	            fprintf(stderr,"\nUnrecognized option: %3s", argv[i]);
                    fprintf(stderr," output filenames have form: <dir>/<root>.<month>.<extent> \n");
                    fprintf(stderr,"Default extent is '.nc' \n");
                    fprintf(stderr,"Use -Or to specify dir/root (required)\n");
                    fprintf(stderr,"Use -Oe to specify file extent (optional)\n");
                    fprintf(stderr," Ex:  -Or/data/gridded/Pacific.quintdeg.oi3d\n");
                    fprintf(stderr," will create 13 files called Pacific.quintdeg.oi3d.jan.nc, Pacific.quintdeg.oi3d.feb.nc,... Pacific.quintdeg.oi3d.clim.nc, in the directory /data/gridded \n");
                    fprintf(stderr,"Use -Oo to specify a single output file\n");
		    exit(1);
	     }
	     break;
	  break;
	  
          case 'P':
             popt = 1;
	     prop_req = (int *) calloc(MAXPROP, sizeof(int));
             nprops = parse_p_option(&argv[i][2], prop_req);
          break;

          case 'S':
             error = (sscanf(&argv[i][2],"%f", &maxdist) == 1) ? 0 : 1;
          break;
	  
          case 'V':
	      parms_name = &argv[i][2];
	  break;
	
          case 'Z':
               if (argv[i][2] == '\0') {
                   nstdlevs = std_depth_init(z_file);
                   fprintf(stdout,"\nStandard depth levels: \n");
                   for (j = 0; j < nstdlevs; ++j) {
                       fprintf(stdout,"  %.1lf", std_depth[j]);
                    }
                    fprintf(stdout,"\n");
                    exit(0);
                 }
                 z_file = fopen(&argv[i][2],"r");
                 if (z_file == NULL) {
                    fprintf(stderr,"\nError opening standard depth file: %s\n",&argv[i][2]);
                    exit(1);
                 }
           break;

	   case 'h':  
	          print_usage(argv[0]);
	          exit(0);
		  
           default:
                        error = 1;

	
	
	} /* end switch */
      } /* end if */
   } /* end for */
    
/*--------------------------------------------*/    
/*  Check syntax of options */ 

    if (!oflag )  { 
       fprintf(stderr,"\nYou must specify output file(s) with -O\n ");
       exit(1);
    }
    if (!bflag || !iflag) {
       fprintf(stderr,"\nYou must specify grid bounds and increment with -B and -I\n ");
       exit(1);
    }  
    if  ((b3d_root == NULL)) {
         fprintf(stderr,"\nYou must specify files for input bin3d gridded fields [-Gr<file_root>] \n");
         exit(1);
    }
    if(fg_avail) {
      if (use_fg_in_pinch)
          fprintf(stderr,"\nFirst Guess values will be substituted where no value can be estimated from the bin3d data\n");
     
    }
    
    if (one_file)
       nmonths = 1;
       
    iclim = nmonths-1;   /* index for clim files */
   
 /*--------------------------------------------*/ 
/*  set up output grid for internal arrays and output files*/  
   
    outgrid.node_offset = 1;
    outgrid.lon0to360 = outgrid.x_min >= 0 ? 1 : 0;
    outgrid.xgreenwich = (outgrid.x_min < 0) && (outgrid.x_max >= 0);
    outgrid.nx = NINT((outgrid.x_max - outgrid.x_min) / outgrid.x_inc);
    outgrid.ny = NINT((outgrid.y_max - outgrid.y_min) / outgrid.y_inc);
    
    hdr_out.xmin = outgrid.x_min;   /* this is only needed to get index for output gridnode */
    hdr_out.xmax = outgrid.x_max;
    hdr_out.ymin = outgrid.y_min;
    hdr_out.ymax = outgrid.y_max;
    hdr_out.xincr = outgrid.x_inc;
    hdr_out.yincr = outgrid.y_inc;
    hdr_out.nx = outgrid.nx;
    hdr_out.ny = outgrid.ny;
    hdr_out.node_offset = outgrid.node_offset;

/*--------------------------------------------*/ 
/*  open parameters file and initialize all the fields in ncparms structures */  

    ncid_parms =  ncparms_open(parms_name);
    if (ncid_parms < 0) {
                fprintf(stdout,"\n");
                exit(0);
    }
    error = ncparms_getinfo(ncid_parms, &ncparms_info, &ncparms_data, 1);
    if (error) {
          fprintf(stderr,"FATAL ERROR:  cannot parse %s correctly\n", parms_name);
	  exit(1);
    }
    
    lon0to360 = ncparms_info.lons[0] >= 0 ? 1 : 0;
    xgreenwich = (ncparms_info.lons[0] < 0) && (ncparms_info.lons[ncparms_info.nx-1] >=0) ? 1:0;
    
/*compare bounds of work region to ncparms file bounds */

    i = ncparms_info.nx - 1;
    inrange = is_in_range((float) outgrid.x_min, (ncparms_info.lons[0]-0.5*ncparms_info.xincr), (ncparms_info.lons[i] + 0.5 *ncparms_info.xincr ), merid, lon0to360) && is_in_range((float)outgrid.x_max, (ncparms_info.lons[0]-0.5*ncparms_info.xincr), (ncparms_info.lons[i]+ 0.5 *ncparms_info.xincr ), merid, lon0to360);
   
    if (inrange) {
       i = ncparms_info.ny -1;
       inrange = ((float) outgrid.y_min >= ncparms_info.lats[0]-0.5*ncparms_info.yincr) && ((float) outgrid.y_max <= ncparms_info.lats[i]+0.5*ncparms_info.yincr);
    }

    free(merid); 
    merid = NULL;

    if (! inrange) {  
       fprintf(stderr,"ERROR:Parameters file [%s] does not include bounds specified here. \n", parms_name);
         exit(1);
    }
/*--------------------------------------------*/    
  /* At this point, info about the output region (outgrid) 
     and parameter file (ncparms_info) is set.   */
/*--------------------------------------------*/ 

   if (fg_avail) {
     /* Get info about FG grids and check for overlap with work region */

       buf = (char *)calloc(200, sizeof(char));
       strcpy(buf, fg_root);
       if (!one_file)
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
    
       if (!check_overlap(&fg_orig, &outgrid)) {
        fprintf(stderr,"ERROR: %s does not overlap the work region:\n", buf);
        fprintf(stderr,"  %.3f/%.3f/%.3f/%.3f\n", outgrid.x_min, outgrid.x_max, outgrid.y_min, outgrid.y_max);
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
       
    } /* end if fg_avail */
/*--------------------------------------------*/
/* Get info about bin3d grids and check for overlap with output region */

    buf = (char *)calloc(200, sizeof(char));
    strcpy(buf, b3d_root);
    if (!one_file)
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

    if (!check_overlap(&b3d_orig, &outgrid)) {
        fprintf(stderr,"ERROR: %s does not overlap the work region:\n", buf);
        fprintf(stderr,"  %.3f/%.3f/%.3f/%.3f\n", outgrid.x_min, outgrid.x_max, outgrid.y_min, outgrid.y_max);
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
/*--------------------------------------------*/   
 /* Check b3d files for requested properties availability.    After this, prop_indx[] and
    nprops will hold the property info for the output cdf files. */
   
   prop_indx = (int *) calloc((size_t)nprops, sizeof(int)); 
   prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int)); 
   
   for (i = 0; i < hdr_b3d.nprops; ++i)    
      prop_avail[get_prop_indx(hdr_b3d.prop_id[i])] = 1;
   
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
   
   
   if (fg_avail) {
   /* now check FG file for properties */
   
      prop_avail = (int *) calloc((size_t)MAXPROP, sizeof(int)); 
      for (i = 0; i < hdr_fg.nprops; ++i)    
         prop_avail[get_prop_indx(hdr_fg.prop_id[i])] = 1;
   
      prop_avail[(int)T90] = prop_avail[(int)TE] || prop_avail[(int)T90];
      n = 0;   
      for (i = 0; i < nprops; ++i) {
         if (!prop_avail[prop_indx[i]]) {
             if (prop_indx[i] == (int) PR || prop_indx[i] == (int) T90 || prop_indx[i] == (int) SA) {
                fprintf(stderr,"\n FATAL ERROR!! Property %s not available in cdf file (pr te sa are mandatory).\n", get_prop_mne(prop_indx[i]));  
	        exit(1);
          }     
          else {
             fprintf(stderr,"\n WARNING!! Property %s not available in First Guess netcdf file and will not be output.\n", get_prop_mne(prop_indx[i]));       
          }
         }
         else {
         prop_indx[n++] = prop_indx[i];
         }
      }
      free(prop_avail); 
      
   }
         
   nprops = n;
   free(prop_req);
      
/*--------------------------------------------*/   
/*  get list of standard depths for output and create output files  */

    nstdlevs = std_depth_init(z_file);   /* creates and fills global array std_depth */
    
    nz_seas_out = 0;
    while (nz_seas_out < nstdlevs && std_depth[nz_seas_out] <= SEASONAL_LAYER_DEF)
       ++nz_seas_out;
    
    cdf_construct(&outgrid, nprops, prop_indx, out_dir, out_root, out_ext, ncid_out, argc, argv);

/*--------------------------------------------*/   
/* Read in topography values for output region + search radius */

    set_topo_info(&outgrid, &topo_info, maxdist+TOPOBUF);
    seafloor = hb_get_topo(toponame, &topo_info, &topolat, &topolon, FALSE, topo_info.lon0to360, &topoval);
    free(topolat);
    free(topolon);

/*--------------------------------------------*/   
/*--------------------------------------------*/ 

   if (fg_avail) {

    /* Determine size of input area for FG fields,
   and allocate memory */
  
      nsq_fg = subgrid(&fg_orig, &fg_in, outgrid.y_min, outgrid.x_min, outgrid.y_max, outgrid.x_max, maxdist);
  
      bdepth_fg = (float *) calloc(nsq_fg, sizeof(float));
      nz_fg = (int *)  calloc(nsq_fg, sizeof(int));

      p_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      d_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      t_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      s_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      sig0_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      sig1_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      sig2_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      sig3_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
      sig4_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
    
      for (i = 0; i < nsq_fg; ++i) {
        d_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        p_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        t_fg[i] = (float *)get_memory(NULL,  hdr_fg.nz, sizeof(float));
        s_fg[i] = (float *)get_memory(NULL,  hdr_fg.nz, sizeof(float));
        sig0_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        sig1_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        sig2_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        sig3_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
        sig4_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));
       }
        /* determine number of upper ocean depth levels.  Add a buffer for interpolation purposes*/    
       din = (float *) calloc(hdr_fg.nz, sizeof(float));
    
       read_cdf_depths(ncid_fg[iclim], din);
       nz_seas_fg = 0;
       while ( nz_seas_fg < hdr_fg.nz && din[nz_seas_fg] <= (SEASONAL_LAYER_DEF+500)) 
         ++nz_seas_fg;
    
       free(din);
    
       for (imonth = 0; imonth < nmonths; ++imonth) {
          up_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
          ud_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
          ut_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
          us_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
          usig0_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
          nz_m_fg[imonth] = (int *)  calloc(nsq_fg, sizeof(int));
          ml_fg[imonth] = (struct MIXLAYER *)get_memory(NULL, nsq_fg, sizeof(struct MIXLAYER));
	
   
          for (i = 0; i < nsq_fg; ++i) {
             up_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
             ud_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
             ut_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
             us_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
             usig0_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
          }
       }
       
    } /* end if fg_avail */
/*--------------------------------------------*/   
/* Determine size of input area for bin3d fields,
   set up structures, allocate memory */

   nsq_b3d = subgrid(&b3d_orig, &b3d_in, outgrid.y_min, outgrid.x_min, outgrid.y_max, outgrid.x_max, maxdist);
   
   bdepth_b3d = (float *) calloc(nsq_b3d, sizeof(float));
   nz_b3d = (int *)  calloc(nsq_b3d, sizeof(int));
     
   p_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   d_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   t_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   s_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   sig0_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   sig1_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   sig2_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   sig3_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   sig4_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
   tcnt_b3d = (short **)get_memory(NULL,nsq_b3d,sizeof(short *));
     
   for (i = 0; i < nsq_b3d; ++i) {
        p_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
        d_b3d[i] = (float *)get_memory(NULL,  hdr_b3d.nz, sizeof(float));
        t_b3d[i] = (float *)get_memory(NULL,  hdr_b3d.nz, sizeof(float));
        s_b3d[i] = (float *)get_memory(NULL,  hdr_b3d.nz, sizeof(float));
        sig0_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
        sig1_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
        sig2_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
        sig3_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
        sig4_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz, sizeof(float));
	tcnt_b3d[i] = (short *)get_memory(NULL, hdr_b3d.nz, sizeof(short));
    }
    
     /* determine number of upper ocean depth levels  Add a  buffer for interpolation purposes*/    
    din = (float *) calloc(hdr_b3d.nz, sizeof(float));
    read_cdf_depths(ncid_b3d[iclim], din);
    
    nz_seas_b3d = 0;
    while ( nz_seas_b3d < hdr_b3d.nz && din[nz_seas_b3d] <= (SEASONAL_LAYER_DEF+500)) 
      ++nz_seas_b3d;
    
    free(din);
    
   
    for (imonth = 0; imonth < nmonths; ++imonth) {
        up_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
        ud_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
        ut_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
        us_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
        usig0_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
	tcnt_m_b3d[imonth] = (short **)get_memory(NULL, nsq_b3d, sizeof(short *));
	ml_b3d[imonth] = (struct MIXLAYER *)get_memory(NULL,nsq_b3d,sizeof(struct MIXLAYER));
        nz_m_b3d[imonth] = (int *) get_memory(NULL,nsq_b3d,sizeof(int));
       
	for (i = 0; i < nsq_b3d; ++i) {
            up_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d, sizeof(float));
            ud_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d, sizeof(float));
            ut_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d, sizeof(float));
            us_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d, sizeof(float));
            usig0_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d, sizeof(float));
            tcnt_m_b3d[imonth][i] = (short *)get_memory(NULL, nz_seas_b3d, sizeof(short));
        }
    }
/*--------------------------------------------*/ 
/*  Allocate memory for other work space */
/*--------------------------------------------*/    
    nsq_out = outgrid.nx * outgrid.ny;
    bdepth_out = (float *) calloc(nsq_out, sizeof(float));
    
    isopyc = (float **)calloc(nsq_out, sizeof(float *));
    xlen = (float **)calloc(nsq_out, sizeof(float *));
    ylen  = (float **)calloc(nsq_out, sizeof(float *));
    parms0 = (float *) calloc(nstdlevs, sizeof(float));
    parms1 = (float *) calloc(nstdlevs, sizeof(float));
    parms2 = (float *) calloc(nstdlevs, sizeof(float));
    
    for (i = 0; i < nsq_out; ++i) {
       isopyc[i] = (float *) calloc(nstdlevs, sizeof(float));
       xlen[i] = (float *) calloc(nstdlevs, sizeof(float));
       ylen[i] = (float *) calloc(nstdlevs, sizeof(float));
    }
    
    for (imonth = 0; imonth < nmonths; ++imonth) {
       isopyc_up[imonth] = (float **) calloc(nsq_out, sizeof(float *));
       xlen_seas[imonth] = (float **) calloc(nsq_out, sizeof(float *));
       ylen_seas[imonth] = (float **) calloc(nsq_out, sizeof(float *));
       psave[imonth] = (float **) calloc(nsq_out, sizeof(float *));
       ssave[imonth] = (float **) calloc(nsq_out, sizeof(float *));
       ml_out[imonth] = (struct MIXLAYER *)get_memory(NULL,nsq_out,sizeof(struct MIXLAYER));

       for (i = 0; i < nsq_out; ++i) {
          isopyc_up[imonth][i] = (float *) calloc(nz_seas_out, sizeof(float));
	  xlen_seas[imonth][i] = (float *) calloc(nz_seas_out, sizeof(float));
	  ylen_seas[imonth][i] = (float *) calloc(nz_seas_out, sizeof(float));
       }
       
    }
/*--------------------------------------------*/   
/*  Open each monthly FG and Bin3d file.  *.clim.nc files are already open */

    if (!one_file) {
    
      if (fg_avail) {
         for (imonth = 0; imonth < iclim; ++imonth) {
           buf = (char *)calloc(2000, sizeof(char));
           strcpy(buf, fg_root);
           strncat(buf, ".", 1);
           strcat(buf, cmonths[imonth]);
           ncid_fg[imonth] = cdf_open(fg_dir, buf, fg_ext, 0);
           free(buf);
	 }
      }
 
      for (imonth = 0; imonth < iclim; ++imonth) {
        buf = (char *)calloc(2000, sizeof(char));
        strcpy(buf, b3d_root);
        strncat(buf, ".", 1);
        strncat(buf, cmonths[imonth], 4);
        ncid_b3d[imonth] = cdf_open(b3d_dir, buf, b3d_ext, 0);
        free(buf);
      }
    }
    
    fprintf(stderr,"\nCompleted setup phase.... \n");
/*--------------------------------------------*/ 
/*--------------------------------------------*/ 
/*--------------------------------------------*/ 
/*--------------------------------------------*/ 


   if (fg_avail) {
      fprintf(stderr,"Reading in first-guess fields...");
   
      /* Get full depth profiles from .clim file */
      imonth = iclim;
      get_sig_profiles(ncid_fg[imonth], &hdr_fg, &fg_in, hdr_fg.nz, bdepth_fg, seafloor, &topo_info, tindex_fg, d_fg, p_fg, t_fg, s_fg, NULL, sig0_fg, sig1_fg, sig2_fg, sig3_fg, sig4_fg, nz_fg);

      /* Get upper ocean profiles from monthly files and close input files for now */
   
      for (imonth = 0; imonth < nmonths; ++imonth) {
       get_monthly_profiles(ncid_fg[imonth], &hdr_fg, &fg_in, nz_seas_fg, bdepth_fg, tindex_fg,  ud_fg[imonth], up_fg[imonth], ut_fg[imonth], us_fg[imonth], NULL, usig0_fg[imonth], nz_m_fg[imonth], ml_fg[imonth]);
       cdf_close(ncid_fg[imonth]);
      }
   }
   
   
   fprintf(stderr,"\nReading in bin3d profiles...");
   
   imonth = iclim;
   get_sig_profiles(ncid_b3d[imonth], &hdr_b3d, &b3d_in, hdr_b3d.nz, bdepth_b3d, seafloor, &topo_info, tindex_b3d, d_b3d, p_b3d, t_b3d, s_b3d, tcnt_b3d, sig0_b3d, sig1_b3d, sig2_b3d, sig3_b3d, sig4_b3d, nz_b3d);

   /* Get upper ocean profiles from monthly files */
   
   for (imonth = 0; imonth < nmonths; ++imonth) {
       get_monthly_profiles(ncid_b3d[imonth], &hdr_b3d, &b3d_in, nz_seas_b3d, bdepth_b3d, tindex_b3d, ud_b3d[imonth], up_b3d[imonth], ut_b3d[imonth], us_b3d[imonth], tcnt_m_b3d[imonth], usig0_b3d[imonth], nz_m_b3d[imonth], ml_b3d[imonth]);
       cdf_close(ncid_b3d[imonth]);
   }
   fprintf(stderr,"  Done. \n");

/*-----------------------------------------------------------------------------------*/    
/*-----------------------------------------------------------------------------------*/    
/*-----------------------------------------------------------------------------------*/    
/* Begin estimation portion of module.  

   Loop for each requested property.   The order of props is set to
  1) pr 2)sa 3) th9  Pressure must be first and the output pressure and salt profiles
  must be saved in memory in order to interpolate all props back onto stddepths
  and to convert th9 profiles back to t90 */
  
  fprintf(stderr,"Working on property fields (a symbol is printed for each grid square completed):");  
  for (iprop = 0; iprop < nprops; ++iprop) {
  
       is_pr = (prop_indx[iprop] == (int)PR);
       is_sa = (prop_indx[iprop] == (int)SA);
       is_te = (prop_indx[iprop] == (int)T90);
       fprintf(stderr,"\n%s ", get_prop_mne(prop_indx[iprop]));
       
       
       switch ((enum property) prop_indx[iprop]) {
  	  case PR :          /* just set the pointers */
             if (fg_avail)  x_fg = p_fg;
	     x_b3d = p_b3d;
	     xcnt_b3d = tcnt_b3d; /* tcnt applies to pressure  */
             for (imonth = 0; imonth < nmonths; ++imonth) {
	        if (fg_avail)  ux_fg[imonth] = up_fg[imonth];
		ux_b3d[imonth] = up_b3d[imonth];
		xcnt_m_b3d[imonth] = tcnt_m_b3d[imonth];
	     }
	     break;
  	  case T90 :          /* just set the pointers */
             if (fg_avail)  x_fg = t_fg;
	     x_b3d = t_b3d;
	     xcnt_b3d = tcnt_b3d;
             for (imonth = 0; imonth < nmonths; ++imonth) {
	        if (fg_avail)  ux_fg[imonth] = ut_fg[imonth];
		ux_b3d[imonth] = ut_b3d[imonth];
		xcnt_m_b3d[imonth] = tcnt_m_b3d[imonth];
	     }
	     break;
	  case SA :         /* just set the pointers */
             if (fg_avail)  x_fg = s_fg;
	     x_b3d = s_b3d;
	     xcnt_b3d = tcnt_b3d;  /* tcnt applies to salinity also */
             for (imonth = 0; imonth < nmonths; ++imonth) {
	        if (fg_avail)  ux_fg[imonth] = us_fg[imonth];
		ux_b3d[imonth] = us_b3d[imonth];
		xcnt_m_b3d[imonth] = tcnt_m_b3d[imonth];
	     }
	     break;
	  
	  default:
	     fprintf(stderr,"Computation of other properties not yet supported\n");
	     goto FINISH;
	     
	       
             xcnt_b3d = (short **) get_memory(NULL,nsq_b3d,sizeof(short *));
             x_b3d = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
             for (i = 0; i < nsq_b3d; ++i) {
                x_b3d[i] = (float *)get_memory(NULL, hdr_b3d.nz , sizeof(float));
	        xcnt_b3d[i] = (short *)get_memory(NULL, hdr_b3d.nz, sizeof(short));
             }
   
             for (imonth = 0; imonth < nmonths; ++imonth) {
                ux_b3d[imonth] = (float **)get_memory(NULL,nsq_b3d,sizeof(float *));
                xcnt_m_b3d[imonth] = (short **)get_memory(NULL,nsq_b3d,sizeof(short *));
	        for (i = 0; i < nsq_b3d; ++i) {
                   ux_b3d[imonth][i] = (float *)get_memory(NULL, nz_seas_b3d , sizeof(float));
                   xcnt_m_b3d[imonth][i] = (short *)get_memory(NULL, nz_seas_b3d, sizeof(short));
                 }
                 buf = (char *)calloc(2000, sizeof(char));
                 strcpy(buf, b3d_root);
                 if (!one_file)  {
                    strncat(buf, ".", 1);
                    strncat(buf, cmonths[imonth], 4);
                 }
		 ncid_b3d[imonth] = cdf_open(b3d_dir, buf, b3d_ext, 0);
                 free(buf);
            }
	     
	     get_x_profiles(ncid_b3d[iclim],&hdr_b3d, &b3d_in, hdr_b3d.nz, bdepth_b3d, x_b3d, xcnt_b3d, prop_indx[iprop]);
             for (imonth = 0; imonth < nmonths; ++imonth) {
	        get_x_profiles(ncid_b3d[imonth],&hdr_b3d, &b3d_in, nz_seas_b3d, bdepth_b3d, ux_b3d[imonth], xcnt_m_b3d[imonth], prop_indx[iprop]);
		cdf_close(ncid_b3d[imonth]);
             }	 
	     
	     if (fg_avail) {
                x_fg = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
                for (i = 0; i < nsq_fg; ++i) 
                  x_fg[i] = (float *)get_memory(NULL, hdr_fg.nz, sizeof(float));

                for (imonth = 0; imonth < nmonths; ++imonth) {
                   ux_fg[imonth] = (float **)get_memory(NULL,nsq_fg,sizeof(float *));
                   for (i = 0; i < nsq_fg; ++i) 
                      ux_fg[imonth][i] = (float *)get_memory(NULL, nz_seas_fg, sizeof(float));
                 
		   buf = (char *)calloc(2000, sizeof(char));
                   strcpy(buf, fg_root);
                   if (! one_file) {
                      strncat(buf, ".", 1);
                      strcat(buf, cmonths[imonth]);
		   }
                   ncid_fg[imonth] = cdf_open(fg_dir, buf, fg_ext, 0);
                   free(buf);
	        }
 	        get_x_profiles(ncid_fg[iclim],&hdr_fg, &fg_in, hdr_fg.nz, bdepth_fg, x_fg, NULL, prop_indx[iprop]);
	     
                for (imonth = 0; imonth < nmonths; ++imonth) {
	           get_x_profiles(ncid_fg[imonth],&hdr_fg, &fg_in, nz_seas_fg, bdepth_fg, ux_fg[imonth], NULL, prop_indx[iprop]);
		   cdf_close(ncid_fg[imonth]);
                }	
	    } /* end if fg_avail */  
      } /* end switch */
       
       /* Visit each square in the output region.  Produce an output profile
          that includes mask and empty flags by kriging at each level, then interpolate
	  the resulting profiles back onto stddepths */
	  
	for (isq = 0; isq < nsq_out; ++isq) {
	
	   if (isq % outgrid.nx == 0)  /* show progress */
	      fprintf(stderr,"*");
	   else
	      fprintf(stderr,"-");
	      
	   error = sq2rc(isq, &outgrid, &theRow, &theCol);
	   error = ij2xy(&outgrid, theCol, theRow, &theLon, &theLat);
	   
	   if (fg_avail) {
	      if (error = xy2ij(&fg_in, theLon, theLat, &col, &row)) {
	        fprintf(stderr,"FATAL ERROR translating output Lat, Lon to first-guess grid\n");
	        exit(1);
	      }
	      theSq_fg = row * fg_in.nx + col;
	   }
	   
	   if (error = xy2ij(&b3d_in, theLon, theLat, &col, &row)) {
	     fprintf(stderr,"FATAL ERROR translating output Lat, Lon to bin3d grid\n");
	     exit(1);
	   }
	   theSq_b3d = row * b3d_in.nx + col;
	   
	      
	   error = get_indices(&hdr_out,(float)theLat, (float)theLon, &row_out, &col_out);

	   
	   /* allocate memory for profile mean, error and count */
	   
	   x_prof = (float *) get_memory(NULL, nstdlevs, sizeof(float));
	   e_prof = (float *) get_memory(NULL, nstdlevs, sizeof(float));
	   n_prof = (short *) get_memory(NULL, nstdlevs, sizeof(short));
	   
	   /* first time through, define the bottom depth, initialize isopyc_above, alloc memory for
	      output length scale arrays */
	   
	   if (is_pr) {
	      topoval = find_nearest_topo_val(theLat, theLon, seafloor, &topo_info);
	      
	      bdepth_out[isq] = bdepth_b3d[theSq_b3d];
	      if (fg_avail) 
	          bdepth_out[isq] = bdepth_fg[theSq_fg] > bdepth_b3d[theSq_b3d] ? bdepth_fg[theSq_fg] : bdepth_b3d[theSq_b3d];
		  
	      define_bottom(&bdepth_out[isq], topoval);
	      isopyc_above = 0.0;
	      xlen_ncout = (float *) get_memory(NULL, nstdlevs, sizeof(float));
	      ylen_ncout = (float *) get_memory(NULL, nstdlevs, sizeof(float));
	   
	   }
	   
	   /* handle case where grid square is on land */
	   
	   if (bdepth_out[isq] > testmask) {
             x_ncout = (float *) calloc(nstdlevs, sizeof(float));
	     n_ncout = (short *) calloc(nstdlevs, sizeof(short));
	     
	     for (ilev = 0; ilev < nstdlevs; ++ilev) {
	        x_ncout[ilev] = HBMASK;
		xlen[isq][ilev] = HBMASK;
		ylen[isq][ilev] = HBMASK;
		if (is_pr) {
		   xlen_ncout[ilev] = HBMASK;
		   ylen_ncout[ilev] = HBMASK;
		}
	     }
	     
	     for (imonth = 0; imonth < nmonths; ++imonth) {
	        error = write_prop_cdf(ncid_out[imonth], x_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	        error = write_prop_err_cdf(ncid_out[imonth], x_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	        error = write_prop_count_cdf(ncid_out[imonth], n_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	     
	        if (is_pr) {
	           write_bottom_depth_cdf(ncid_out[imonth], row_out, col_out, tbin, 1, 1, 1, &bdepth_out[isq]);
		   write_xlength_cdf(ncid_out[imonth], row_out, col_out, tbin, 0, 1, 1, 1, nstdlevs, xlen_ncout);
		   write_ylength_cdf(ncid_out[imonth], row_out, col_out, tbin, 0, 1, 1, 1, nstdlevs, ylen_ncout);
	        }
	     }
	     
	     free(x_ncout);
	     free(n_ncout);
	     if (is_pr) {
	        free(xlen_ncout);
		free(ylen_ncout);
	     }
	     
             nmask += nstdlevs;
	     continue;
	   }
	   
       /* get oi_parameters for this location */
	   
	   error = ncparms_read(ncid_parms,(float)theLat, (float)theLon, &ncparms_info, &ncparms_data);
	   get_prop_parms(&ncparms_data, prop_indx[iprop],std_depth, nstdlevs, bdepth_out[isq], parms0, parms1, parms2);
	   if (is_pr) {
	      define_length_scales(&ncparms_data, xlen[isq], ylen[isq], std_depth, nstdlevs, bdepth_out[isq]);
	      for (ilev = 0; ilev < nz_seas_out; ++ilev) {
	         for (imonth = 0; imonth < nmonths; ++imonth) {
	            xlen_seas[imonth][isq][ilev] = xlen[isq][ilev];
	            ylen_seas[imonth][isq][ilev] = ylen[isq][ilev];
		 }
	      }
	   }
	   	   
       /********************************************************/
       /* Begin by working on levels below the seasonal layer */
	   
	   for (ilev = nz_seas_out; ilev < nstdlevs; ++ilev) {
	     
	     curr_depth = (float)std_depth[ilev];
	     if (ilev == nstdlevs-1)
	        curr_depth = bdepth_out[isq];
		
	     if (curr_depth > bdepth_out[isq]) {
	           isopyc[isq][ilev] = HBMASK;
		   x_prof[ilev] = HBMASK;
		   e_prof[ilev] = HBMASK;
		   n_prof[ilev] = 0;
		   ++nmask;
		   continue;
	     }
	     
	     /* Determine appropriate potential density */
	     
	       if (curr_depth >= 3500) {
	          if (fg_avail) sigptr_fg = sig4_fg;
		  sigptr_b3d = sig4_b3d;
		  reflev = 4000;
	       }
	       else if (curr_depth >= 2500 && curr_depth < 3500) {
	          if (fg_avail) sigptr_fg = sig3_fg;
		  sigptr_b3d = sig3_b3d;
		  reflev = 3000;
	       }
	       else if (curr_depth >= 1500 && curr_depth < 2500) {
	          if (fg_avail) sigptr_fg = sig2_fg;
		  sigptr_b3d = sig2_b3d;
		  reflev = 2000;
	       }
	       else if (curr_depth >= 500 && curr_depth < 1500) {
	          if (fg_avail) sigptr_fg = sig1_fg;
		  sigptr_b3d = sig1_b3d;
		  reflev = 1000;
	       }
	       else {
	          if (fg_avail) sigptr_fg = sig0_fg;
		  sigptr_b3d = sig0_b3d;
		  reflev = 0;
	       }
	       
	     /***************************************/
	       
	     if (is_pr) {     /* First time at this square....*/
	       
	     /* Define isopycnal for this level from bin3d values available */
                  		  
	         if ( xlen[isq][ilev] > testmask || xlen[isq][ilev] < testempty) 
		      xlen[isq][ilev] = DIST_DEFAULT;
	         if ( ylen[isq][ilev] > testmask || ylen[isq][ilev] < testempty) 
		      ylen[isq][ilev] = DIST_DEFAULT;
		      
		 isopyc[isq][ilev] = find_isopycnal(curr_depth, theSq_b3d, &b3d_in, d_b3d, sigptr_b3d, nz_b3d, xlen[isq][ilev], ylen[isq][ilev], bdepth_b3d);
		 
		      
		 /* if no values available, compute it from the first-guess fields */
		 if (isopyc[isq][ilev] <= 0  && fg_avail) {
		    isopyc[isq][ilev] = sigptr_fg[theSq_fg][ilev];
		    if (isopyc[isq][ilev] <= 0) 	
		       isopyc[isq][ilev] = find_isopycnal(curr_depth, theSq_fg, &fg_in, d_fg, sigptr_fg, nz_fg, xlen[isq][ilev], ylen[isq][ilev], bdepth_fg);
		 
		 }
		  
	       if (isopyc[isq][ilev] > 0 && isopyc[isq][ilev] < 100) {
	         /* check for density inversion */
		 
		  if (isopyc[isq][ilev] < isopyc_above )
		     isopyc[isq][ilev] = isopyc_above;
		     
		  isopyc_above =  isopyc[isq][ilev];  
	       }
	       
	     }  /* end if is_pr */
	     
	     /***************************************/
	     xlength = xlen[isq][ilev];
	     ylength = ylen[isq][ilev];
	     adjust_lengths = FALSE;
	     if (is_pr) adjust_lengths = TRUE;
	     
	     n_prof[ilev] = extract_surface(isopyc[isq][ilev], reflev, curr_depth, &b3d_in, sigptr_b3d, x_b3d, xcnt_b3d, d_b3d, nz_b3d, bdepth_b3d, theLat, theLon, theSq_b3d, &xlength, &ylength, &b3dwork, &xsurf_b3d,&xcnt_surf_b3d, &nwsq_b3d, adjust_lengths);
	     
	     xlen[isq][ilev] = xlength;
	     ylen[isq][ilev] = ylength;
	     theFG = 0.0;
	     
	     if (fg_avail) {
        	     n = extract_surface(isopyc[isq][ilev], reflev, curr_depth, &fg_in, sigptr_fg, x_fg, (short **)NULL, d_fg, nz_fg, bdepth_fg, theLat, theLon, theSq_fg, &xlength, &ylength, &fgwork, &xsurf_fg,(short **)NULL, &nwsq_fg, FALSE);
	     
	     
	        /* Save first-guess property to add back to the mean */
	     
	        theFG = xsurf_fg[nwsq_fg/2];
	        if (theFG < testempty) {
	          theFG = weighted_stats(&fgwork, xsurf_fg, nwsq_fg,  xlen[isq][ilev], ylen[isq][ilev], (float *)NULL);
		  xsurf_fg[nwsq_fg/2] = theFG;
	        }
	     
	        do_krig = theFG > testempty;
	     
	     
	        if (do_krig) {
	           resids = compute_residuals(&fgwork, xsurf_fg, nwsq_fg, &b3dwork, xsurf_b3d, nwsq_b3d, &n_prof[ilev]);
		   
		   if (resids != NULL) {  
		     for (i = 0; i < nwsq_b3d; ++i) 
	                 xsurf_b3d[i] = (double)resids[i];  
		   }   
		   
		   free(resids);
		   resids = NULL; 
		}
	     } /* end if fg_avail */
	     
	     do_krig = n_prof[ilev] >= minobs;
	      
	      if (parms1[ilev] < testempty || parms2[ilev] < testempty)
	          do_krig = 0;
	     
              if ( do_krig)  {
	           oi_params[0] = parms0[ilev];
	           oi_params[1] = parms1[ilev];
	           oi_params[2] = parms2[ilev];
	      
 	           error = krig(theLat, theLon, xsurf_b3d, xcnt_surf_b3d, nwsq_b3d, &b3dwork, ncparms_info.modelCode, oi_params, &xmean, &xerr, &xcount);
	           if (error) do_krig = 0;
		
		   ++nkrig;
	      }
		
	      if (!do_krig)  {
	         if (theFG < testempty)
		      theFG = 0.0;  
	         xmean = weighted_stats(&b3dwork, xsurf_b3d, nwsq_b3d, xlen[isq][ilev], ylen[isq][ilev], &xerr);
		 if (xerr > testempty)
		      xerr = sqrtf(xerr);
		 xcount = n_prof[ilev];
		 ++nother;
	      }
	     
	     if (xmean > testempty) {
	        x_prof[ilev] = (float)xmean + theFG;
	        e_prof[ilev] = (float)xerr;
	        n_prof[ilev] = xcount;
	     }
	     else {
	        if (fg_avail && use_fg_in_pinch) {
	           x_prof[ilev] = (float)theFG;
	           e_prof[ilev] = (float)HBEMPTY;
	           n_prof[ilev] = -1;   /* flags as FG value */
		   ++nfg_used;
		}
		else {
	           x_prof[ilev] = (float)HBEMPTY;
	           e_prof[ilev] = (float)HBEMPTY;
	           n_prof[ilev] = 0;
	           ++nempty;
		}
	     }
	     
	     if (is_pr) {
	        xlen_ncout[ilev] = xlen[isq][ilev];
		ylen_ncout[ilev] = ylen[isq][ilev];
	     }
	     
	     /***************************************/
	     free(xcnt_surf_b3d);
	     free(xsurf_b3d);
	     if (fg_avail)  free(xsurf_fg);
	     
	   }  /* end for ilev */
 /***********************************************************/
 /***********************************************************/
 /* work on each monthly upper ocean layer and append the levels beneath the seasonal
	  layer to produce a monthly profile for output */
	   
	   for (imonth = 0; imonth < nmonths; ++imonth) {
	      
          /************************************************/
          /* First determine the mixed layer properties  */
          /***********************************************/
/*	      fprintf(stderr,"%c",cmonths[imonth][0]);  */
	      ilev = 0;
	      curr_depth = (float)std_depth[ilev];
	      reflev = 0;  
	      sigptr_b3d = usig0_b3d[imonth];
	      if (fg_avail) 
	          sigptr_fg = usig0_fg[imonth];
	      
	      x_ml = NULL;
	      e_ml = NULL;
	      n_ml = NULL;
	      nz_ml = 0;
	      isopyc_above = 0;
	      
	      if (is_pr) {
		  ml_out[imonth][isq].theta = (float) HBEMPTY;
		  ml_out[imonth][isq].salinity = (float) HBEMPTY;
		  
		  /* determine density, theta, salinity, depth properties */
		  
		  xlength = xlen_seas[imonth][isq][ilev];
		  ylength = ylen_seas[imonth][isq][ilev];
	          n = set_ml_props(ml_b3d[imonth], &b3d_in, theLat, theLon, theSq_b3d, &xlength, &ylength, &ml_out[imonth][isq]);
		  
		  nz_ml = 0;
	          while (ml_out[imonth][isq].depth >= std_depth[nz_ml]) 
	              ++nz_ml;
		  
		  if (nz_ml) {    
		     x_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		     e_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		     n_ml = (short *) get_memory(NULL, nz_ml, sizeof(short));
		  
	      	     for (ilev=0; ilev < nz_ml; ++ilev) {
		        if (std_depth[ilev] == 0.0)  /* set this explicitly to avoid floating pt errors */
			   x_ml[ilev] = 0.0;
			else   
		           x_ml[ilev] = (float) hb_p80(std_depth[ilev], theLat); 
			   
		        e_ml[ilev] = 0; 
			if (ilev == nz_ml-1)  /* error associated with depth of mixed layer */
			   e_ml[ilev] = ml_out[imonth][isq].derr;
			   
		        n_ml[ilev] = (short) ml_out[imonth][isq].nobs;
			if( ilev < nz_seas_out) {
			   xlen_seas[imonth][isq][ilev] = xlength; 
			   ylen_seas[imonth][isq][ilev] = ylength; 
			}
			xlen_ncout[ilev] = xlength;
			ylen_ncout[ilev] = ylength;
	             }
		     
		     isopyc_above = ml_out[imonth][isq].density + mix_layer_delta;
	         }
	      }
	      else if (is_sa || is_te) {
	      
		  if (ml_out[imonth][isq].density > testempty) {
		    nz_ml = 0;
	            while (ml_out[imonth][isq].depth >= std_depth[nz_ml]) 
	              ++nz_ml;
		  
		    if (nz_ml) {    
		       x_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		       e_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		       n_ml = (short *) get_memory(NULL, nz_ml, sizeof(short));
		       
		       if (is_sa) {
	      	          for (ilev=0; ilev < nz_ml; ++ilev) {
		  
		             x_ml[ilev] = (float) ml_out[imonth][isq].salinity; 
		             e_ml[ilev] = (float) ml_out[imonth][isq].serr; 
		             n_ml[ilev] = ml_out[imonth][isq].nobs; 
	                  }
		       
		       }
		         	 
		       if (is_te) {	
			  for (ilev=0; ilev < nz_ml; ++ilev){
			     /* convert to in situ temperatures */
			     x_ml[ilev] = (float) hb_theta((double)ml_out[imonth][isq].salinity, (double)ml_out[imonth][isq].theta, 0.0, (double)psave[imonth][isq][ilev]);
		             e_ml[ilev] = (float) ml_out[imonth][isq].terr; 
		             n_ml[ilev] = ml_out[imonth][isq].nobs; 
			  }
		       }
		       
		    } /* end if nz_ml */
	      
	          } /* end if ml_out > testempty */
	      }
	      
	      else {  /* other properties besides pr, te, sa */
	           
		  if (ml_out[imonth][isq].density > testempty) {
		  
		     xlength = xlen_seas[imonth][isq][ilev];
		     ylength = ylen_seas[imonth][isq][ilev];
		     
	             n_prof[ilev] = extract_surface(ml_out[imonth][isq].density, reflev, curr_depth, &b3d_in, sigptr_b3d, ux_b3d[imonth], xcnt_m_b3d[imonth], ud_b3d[imonth], nz_m_b3d[imonth], bdepth_b3d, theLat, theLon, theSq_b3d, &xlength, &ylength, &b3dwork, &xsurf_b3d, &xcnt_surf_b3d, &nwsq_b3d, FALSE);
		     
		     do_krig = n_prof[ilev] >= minobs;
		     theFG = 0.0;
		     
                     if (fg_avail) {
	          
		        n = extract_surface(ml_out[imonth][isq].density, reflev, curr_depth, &fg_in, sigptr_fg, ux_fg[imonth], (short **)NULL, ud_fg[imonth], nz_m_fg[imonth], bdepth_fg, theLat, theLon, theSq_fg, &xlength, &ylength, &fgwork, &xsurf_fg,(short **)NULL, &nwsq_fg, FALSE);

	               theFG = xsurf_fg[nwsq_fg/2];
	               if (theFG < testempty && n > 0) {
	                  theFG = weighted_stats(&fgwork, xsurf_fg, nwsq_fg,  xlength, ylength, (float *)NULL);
		          xsurf_fg[nwsq_fg/2] = theFG;
	                }
		 
	               do_krig = theFG > testempty;
	               resids = (float *) NULL;

  
                       if (do_krig) {		
		         resids = compute_residuals(&fgwork, xsurf_fg, nwsq_fg, &b3dwork, xsurf_b3d, nwsq_b3d, &n_prof[ilev]);
			 
		         if (resids != NULL && n_prof[ilev] > 0 && n > 0) {  
		          for (i = 0; i < nwsq_b3d; ++i) 
	                     xsurf_b3d[i] = (double)resids[i];  
		         }
			 
			 if (resids != NULL)
	                     free(resids);
		      }    
		     
                    } /* end if fg_avail */

                    if (parms1[ilev]< testempty ||  parms2[ilev] < testempty)
		        do_krig = 0;
			
		    if ( do_krig)  {
	              oi_params[0] = parms0[ilev];
	              oi_params[1] = parms1[ilev];
	              oi_params[2] = parms2[ilev];
	      
 	           error = krig(theLat, theLon, xsurf_b3d, xcnt_surf_b3d, nwsq_b3d, &b3dwork, ncparms_info.modelCode, oi_params, &xmean, &xerr, &xcount);
		
	             if (error) do_krig = 0;
		     else    ++nkrig;
	            }
		 
	            if (!do_krig)  {
		    
	              if (theFG < testempty || n_prof[ilev] == 0) 
		          theFG = 0.0;  
			  
	              xmean = weighted_stats(&b3dwork, xsurf_b3d, nwsq_b3d, xlength, ylength, &xerr);
		      if (xerr > testempty)
		        xerr = sqrtf(xerr);
			
		      if (xmean > testempty) {
		         xcount = n_prof[ilev];
			 if (xcount == 0)
			    ++xcount;
		         ++nother;
		      }
		      else  {  /* couldn't compute any surface values this way */
		        ml_out[imonth][isq].depth = -1;
                        ml_out[imonth][isq].density  = 0;
		      }		      
	            }
		      
		    nz_ml = 0;
	            while (ml_out[imonth][isq].depth >= std_depth[nz_ml]) 
	              ++nz_ml;
		  
		    if (nz_ml) {    
		       x_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		       e_ml = (float *) get_memory(NULL, nz_ml, sizeof(float));
		       n_ml = (short *) get_memory(NULL, nz_ml, sizeof(short));
		  
	      	       for (ilev=0; ilev < nz_ml; ++ilev) {
		  
		          x_ml[ilev] = (float) xmean + theFG; 
		          e_ml[ilev] = (float) xerr; 
		          n_ml[ilev] = xcount; 
	               }
		    }
		    
	            free(xcnt_surf_b3d);
	            free(xsurf_b3d);
	            if (fg_avail) 
		         free(xsurf_fg);
		 } /* end if ml_out > testempty */    
	      }  /* end if is_pr - else */
	      
	  /************************************************/
          /* Next the levels below the mixed layer.....  */
	  /**********************************************/
	      for (ilev = nz_ml; ilev < nz_seas_out; ++ilev) {
	          
	         curr_depth = (float)std_depth[ilev];
	         reflev = 0;  
		 sigptr_fg = usig0_fg[imonth];
		 sigptr_b3d = usig0_b3d[imonth];
		 
	         if (curr_depth > bdepth_out[isq]) {
	           isopyc[isq][ilev] = HBMASK;
		   x_prof[ilev] = HBMASK;
		   e_prof[ilev] = HBMASK;
		   n_prof[ilev] = 0;
		   if (is_pr) {
		      xlen_ncout[ilev] = HBMASK;
		      ylen_ncout[ilev] = HBMASK;
		   }
		   ++nmask;
		   continue;
	         }
		 

	       /* Define isopycnal for this level  */
                 
		 if (is_pr || isopyc_up[imonth][isq][ilev] < 1.0) {
	            if ( xlen_seas[imonth][isq][ilev] > testmask ) 
		          xlen_seas[imonth][isq][ilev] = DIST_DEFAULT;
	            if ( ylen_seas[imonth][isq][ilev] > testmask) 
		          ylen_seas[imonth][isq][ilev] = DIST_DEFAULT;
	             
		    isopyc_up[imonth][isq][ilev] = find_isopycnal(curr_depth, theSq_b3d, &b3d_in, up_b3d[imonth], sigptr_b3d, nz_m_b3d[imonth], xlen[isq][ilev], ylen[isq][ilev], bdepth_b3d);
		    
		    /* if unable to identify an appropriate isopycnal in bin3d fields, use first guess */
		    if (isopyc_up [imonth][isq][ilev] <= 0  && fg_avail) {	
		       isopyc_up[imonth][isq][ilev] = sigptr_fg[theSq_fg][ilev];
		       if (isopyc_up [imonth][isq][ilev] <= 0)
		           isopyc_up[imonth][isq][ilev] = find_isopycnal(curr_depth, theSq_fg, &fg_in, up_fg[imonth], sigptr_fg, nz_m_fg[imonth], xlen[isq][ilev], ylen[isq][ilev], bdepth_fg);
		    } 
		    
		     		     
		    if (isopyc_up[imonth][isq][ilev] <= isopyc_above)
		        /* ensure isopycnal is denser than the level above */
		        isopyc_up[imonth][isq][ilev] = isopyc_above + 0.001;
			
		    
                    if (isopyc_up [imonth][isq][ilev] > 0) 
		       isopyc_above = isopyc_up[imonth][isq][ilev];
		 }
		 
		 adjust_lengths = FALSE;
		 if (is_pr) adjust_lengths = TRUE;
		 
		 xlength = xlen_seas[imonth][isq][ilev];
		 ylength = ylen_seas[imonth][isq][ilev];
		 
	         n_prof[ilev] = extract_surface(isopyc_up[imonth][isq][ilev], reflev, curr_depth, &b3d_in, sigptr_b3d, ux_b3d[imonth], xcnt_m_b3d[imonth], ud_b3d[imonth], nz_m_b3d[imonth], bdepth_b3d, theLat, theLon, theSq_b3d, &xlength, &ylength, &b3dwork, &xsurf_b3d, &xcnt_surf_b3d, &nwsq_b3d, adjust_lengths);
		 
		 xlen_seas[imonth][isq][ilev] = xlength;
		 ylen_seas[imonth][isq][ilev] = ylength;
		  
		 do_krig = n_prof[ilev] >= minobs;
		 theFG = 0.0;

	          
                 if (fg_avail) {
		    n = extract_surface(isopyc_up[imonth][isq][ilev], reflev, curr_depth, &fg_in, sigptr_fg, ux_fg[imonth], (short **)NULL, ud_fg[imonth], nz_m_fg[imonth], bdepth_fg, theLat, theLon, theSq_fg, &xlength, &ylength, &fgwork, &xsurf_fg,(short **)NULL, &nwsq_fg, FALSE);
                
	           theFG = xsurf_fg[nwsq_fg/2];
	           if (theFG < testempty) {
	             theFG = weighted_stats(&fgwork, xsurf_fg, nwsq_fg,  xlen_seas[imonth][isq][ilev], ylen_seas[imonth][isq][ilev], (float *)NULL);
		     xsurf_fg[nwsq_fg/2] = theFG;
	             if (theFG < testempty) 
		       theFG = 0.0;
	           }
		
	           do_krig = theFG > testempty;
	           resids = (float *) NULL;
	        
                   if (do_krig ) {		
		     resids = compute_residuals(&fgwork,xsurf_fg, nwsq_fg, &b3dwork, xsurf_b3d, nwsq_b3d, &n_prof[ilev]);
		     
		     
		      if (resids != NULL) {  
		        for (i = 0; i < nwsq_b3d; ++i) 
	                   xsurf_b3d[i] = (double)resids[i];  
		     
		      free(resids);
		     }    
	       
		     free(xsurf_fg);
 
	           } 
		
                 } /* end if fg_avail*/

                if (parms1[ilev] < testempty || parms2[ilev] < testempty) 
		   do_krig = 0;
		
		if ( do_krig)  {
	           oi_params[0] = parms0[ilev];
	           oi_params[1] = parms1[ilev];
	           oi_params[2] = parms2[ilev];
	      
 	           error = krig(theLat, theLon, xsurf_b3d, xcnt_surf_b3d, nwsq_b3d, &b3dwork, ncparms_info.modelCode, oi_params, &xmean, &xerr, &xcount);
		
		     
	           if (error) do_krig = 0;
		   else    ++nkrig;
	        }

	        if (!do_krig)  {
	            xmean = weighted_stats(&b3dwork, xsurf_b3d, nwsq_b3d, xlen_seas[imonth][isq][ilev], ylen_seas[imonth][isq][ilev], &xerr);
		    if (xerr > testempty)
		      xerr = sqrtf(xerr);
		    
		    xcount = n_prof[ilev];
		    ++nother;
	        }

		if (xmean > testempty) {
		
		     x_prof[ilev] = xmean + theFG;
		     e_prof[ilev] = xerr;
		     n_prof[ilev] = xcount;
		     
		     if (is_pr) {
		     
		        /* ensure pressure is deeper than level above */
			
		        if (ilev == nz_ml) {
			  pbml = hb_p80(ml_out[imonth][isq].depth, theLat);
			  if ( x_prof[ilev] < pbml)
			     x_prof[ilev] = pbml + 1.0;
			}
		         else {
			   if (x_prof[ilev] < x_prof[ilev-1])
			       x_prof[ilev] = x_prof[ilev-1] + 1.0;
			 }
			 
			 if (curr_depth == 0.0)   /* set this explicitly to counteract roundoff errors */
			     x_prof[ilev] = 0.0;
			     
		     }
		}
	        else {
	           if (fg_avail && use_fg_in_pinch) {
	              x_prof[ilev] = (float)theFG;
	              e_prof[ilev] = (float)HBEMPTY;
	              n_prof[ilev] = -1;   /* flags as FG value */
		      ++nfg_used;
		   }
		   else {
	              x_prof[ilev] = (float)HBEMPTY;
	              e_prof[ilev] = (float)HBEMPTY;
	              n_prof[ilev] = 0;
	              ++nempty;
		   }
	        }
		if (is_pr) {
		   xlen_ncout[ilev] = xlen_seas[imonth][isq][ilev];
		   ylen_ncout[ilev] = ylen_seas[imonth][isq][ilev];
                }
		  
		free(xcnt_surf_b3d);
	        free(xsurf_b3d);

	     }/* end for ilev */
	     
	     
	     /* save entire pr profile for each month to interpolate other properties later*/
	     
	     if (is_pr) {
	       psave[imonth][isq]= (float *)get_memory((void *)NULL, nstdlevs, sizeof(float));
	       for  (ilev = 0; ilev < nz_ml; ++ilev)
	          psave[imonth][isq][ilev] = x_ml[ilev];
		  
	       for (ilev = nz_ml; ilev < nstdlevs; ++ilev)
	          psave[imonth][isq][ilev] = x_prof[ilev];
	     }
	     
	     /* save entire sa profile for each month to convert TH9 to T90 later*/
	     if (is_sa) {
	       ssave[imonth][isq]= (float *)get_memory((void *)NULL, nstdlevs, sizeof(float));
	       for  (ilev = 0; ilev < nz_ml; ++ilev)
	          ssave[imonth][isq][ilev] = ml_out[imonth][isq].salinity;
	       for (ilev = nz_ml; ilev < nstdlevs; ++ilev)
	          ssave[imonth][isq][ilev] = x_prof[ilev];
	     }
	     
	     /* adiabatically adjust potential temperature to pressure at each level */
	     if (is_te) {
	     
	        tsave = x_prof;  /* save the original theta for subsequent months */
		x_prof = (float *)get_memory((void *)NULL, nstdlevs, sizeof(float));
	        for (ilev=nz_ml; ilev < nstdlevs; ++ilev) {
		   if (x_prof[ilev] > -3.0 && x_prof[ilev] < 100.0)
		      x_prof[ilev] = (float) hb_theta((double)ssave[imonth][isq][ilev], (double)tsave[ilev], 0.0, (double)psave[imonth][isq][ilev]);
		}
	        free(ssave[imonth][isq]);
	     }  
	     
             x_ncout = (float *) calloc(nstdlevs, sizeof(float));
	     e_ncout = (float *) calloc(nstdlevs, sizeof(float));
	     n_ncout = (short *) calloc(nstdlevs, sizeof(short));
	     
	     x_to_stddepth(&ml_out[imonth][isq], psave[imonth][isq], x_ml, e_ml, n_ml, nz_ml, x_prof, e_prof, n_prof, x_ncout, e_ncout, n_ncout, bdepth_out[isq], std_depth, nstdlevs, theLat); 

	      
	     error = write_prop_cdf(ncid_out[imonth], x_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	     error = write_prop_err_cdf(ncid_out[imonth], e_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	     error = write_prop_count_cdf(ncid_out[imonth], n_ncout, get_prop_mne(prop_indx[iprop]), row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
	     
	     if (is_pr) {
	        write_bottom_depth_cdf(ncid_out[imonth], row_out, col_out, tbin, 1, 1, 1, &bdepth_out[isq]);
		write_xlength_cdf(ncid_out[imonth], row_out, col_out, tbin, 0, 1, 1, 1, nstdlevs, xlen_ncout);
		write_ylength_cdf(ncid_out[imonth], row_out, col_out, tbin, 0, 1, 1, 1, nstdlevs, ylen_ncout);
	     }
		
	     
	     if (is_te) {
	        /* update the pressure count to be same as t90_cnt*/
	        error = write_prop_count_cdf(ncid_out[imonth], n_ncout, "pr", row_out, col_out, tbin, 0, 1,1,1,nstdlevs);
		
		free(x_prof);  /* restore potential temperature profile */
		x_prof = tsave;
	     }
	     
	     free(x_ncout);
	     free(e_ncout);
	     free(n_ncout);

	     if (x_ml != NULL) {
	        free(x_ml);
	        free(e_ml);
	        free(n_ml);
	     }
	     
	   } /* end for imonth */
	   
	   free(x_prof);
	   free(e_prof);
	   free(n_prof);
	   if (is_pr) {
	      free(xlen_ncout);
	      free(ylen_ncout);
	   }
	   
	}  /* end for isq */

        if (is_pr)  
	   free(seafloor);
	
	
	if (! ( is_pr || is_te || is_sa))  {
	
	   /* x-arrays were explicitly allocated -- free up memory */
	    if (fg_avail) {
	       for (isq = 0; isq < nsq_fg; ++isq) {
	          free(x_fg[isq]);
	       }
	       free(x_fg);
	    }
	    
	    for (isq = 0; isq < nsq_b3d; ++isq) {
               free(x_b3d[isq]);
               free(xcnt_b3d[isq]);
	    }
	    
            free(x_b3d);
            free(xcnt_b3d);
            for (imonth = 0; imonth < nmonths; ++imonth) {
	       if (fg_avail) {
	          for (isq = 0; isq < nsq_fg; ++isq) {
	             free(ux_fg[imonth][isq]);
	          }
	          free(ux_fg[imonth]);
	       }
		  
	       for (isq = 0; isq < nsq_b3d; ++isq) {
	          free(ux_b3d[imonth][isq]);
		  free(xcnt_m_b3d[imonth][isq]);
	       }
               free(ux_b3d[imonth]);
               free(xcnt_m_b3d[imonth]);
            }
	} /* end if !is_pr... */
	
  }  /* end for iprop */
   
   
FINISH:

   for (imonth = 0; imonth < nmonths; ++imonth) 
      nc_close(ncid_out[imonth]);
   
   fprintf(stderr, "\n\nNumber kriged: %d\n", nkrig);
   fprintf(stderr, "     averaged: %d\n", nother);
   fprintf(stderr, "      FG used: %d\n", nfg_used);
   fprintf(stderr, "        empty: %d\n", nempty);
   fprintf(stderr, "       masked: %d\n", nmask);
   fprintf(stderr,"\n End of %s\n", argv[0]);
 
   
   exit(0);
} /* end main */


/****************************************************************************/
void print_usage(char *program)
{
  fprintf(stderr,"\nUses an optimal interpolation algorithm (i.e.oridinary kriging) to estimate mean and error variance of measured properties based on estimated covariance functions. ");
  fprintf(stderr," ");
  fprintf(stderr,"\nSpatially varying covariance parameters (sill, length as a function of lag) ");
  fprintf(stderr,"are supplied in the HydroBase3 file lib/global_oi_parms.nc. X and Y length scales (search ellipse) ");
  fprintf(stderr,"are determined from the magnitude of the zonal and ");
  fprintf(stderr,"meridional temperature gradients along isopycnals estimated locally at each gridnode.  If a minimum number of ");
  fprintf(stderr,"observations is not available, the search ellipse is expanded until ");
  fprintf(stderr,"the minimum is achieved or a maximum distance is reached -- both can be optionally specified. ");
  fprintf(stderr,"Monthly first-guess fields (supplied in lib directory) are subtracted from the observed property values ");
  fprintf(stderr,"along an isopycnal surface to provide a field of residuals as a function of distance (lag) from the point being estimated. ");
  fprintf(stderr,"If the isopycnal runs into topography, observed points further away ");
  fprintf(stderr,"in that direction are not included in the estimate to prevent ");
  fprintf(stderr,"mixing of watermasses across topographic barriers. ");
  
  fprintf(stderr,"\n\nInput observations MUST be preprocessed with hb_bin3d. ");
  fprintf(stderr,"\nInput/output grids need NOT have the same dimensions. ");
  fprintf(stderr,"\nOutput grid dimensions are specified with -B<w/e/s/n>, ");
  fprintf(stderr,"\n-I<x/yincr> and -Z<stddepth_file>");
  
  fprintf(stderr,"\n\nUSAGE:  %s -B<w/e/s/n> [-F[d|e|r|o]<first_guess_name>]|s  -G[d|e|r|o]<bin3d_name> -I<xinc>[/<yinc>] [-N<minobs>]  -O[d|e|r|o]<outfile_name> -P<properties>  [-S<max_search_radius>] [-V<parameter_filename>] [-Z<stddepth_file>] [-h]\n\n", program);
  fprintf(stderr," -B  sets the output grid bounds w/e/s/n.\n");
  fprintf(stderr," -F  [OPTIONAL} specifies name of file(s) containing first-guess fields.\n");
  fprintf(stderr,"     Filenames have form <dir>/<root>.<month>.<ext>\n");
  fprintf(stderr,"     UNLESS only one input/output is being worked on.\n");
  fprintf(stderr,"     Use -Fo<filename> to specify a single first guess file.\n");
  fprintf(stderr,"     Use r,d or e to specify one or more of these components.\n");
  fprintf(stderr,"     -Fr<name> mandatory, gives the root name (or dir/root) for these files\n");
  fprintf(stderr,"     -Fd<dir> optionally specifies a directory where they are stored\n");
  fprintf(stderr,"     -Fe<ext> optionally specifies the file extent (suffix).  Default file extent is .nc\n");
  fprintf(stderr,"     The <month> component is automatically added to filenames.\n");
  fprintf(stderr,"     ex: -Fr/d1/HB3/gridded/Atlantic.1deg.FG  will open and read 13 files\n");
  fprintf(stderr,"         (monthly + clim) called Atlantic.1deg.FG.jan.nc, Atlantic.1deg.FG.feb.nc,  etc \n\n");
  fprintf(stderr,"     ex: -Fo/d1/HB3/gridded/Atlantic.quintdeg.FG.djf.nc will open a single first-guess file.\n");
  fprintf(stderr,"     Use -Fs to substitute First Guess values for gridnodes where properties cannot be estimated from bin3d data\n");
  fprintf(stderr,"     The default is no substitution.\n");
  fprintf(stderr," -G  specifies name of gridded observation files (produced by hb_bin3d). \n");
  fprintf(stderr,"     bin3d filenames have form <dir>/<root>.<month>.<ext>\n");
  fprintf(stderr,"     and are specified with r,d,e  as in -F (explained above)\n");
  fprintf(stderr,"     ex: -Gr/d1/HB3/gridded/Atlantic.quintdeg.bin3d  will open and read 13 files\n");
  fprintf(stderr,"          called Atlantic.quintdeg.bin3d.jan.nc, Atlantic.quintdeg.bin3d.feb.nc ... Atlantic.quintdeg.bin3d.clim.nc in the directory /d1/HB3/gridded\n\n");
  fprintf(stderr,"     ex:  -Go/d1/HB3/gridded/Atlantic.quintdeg.bin3d.djf.nc will open a single bin3d file.\n");
  fprintf(stderr," -I  sets the output grid spacing for x- & y-dimensions.\n");
  fprintf(stderr," -O  specifies output .nc files  as in -F and -G above.\n");
  fprintf(stderr,"     ex: -Or/d1/HB3/gridded/Atlantic.quintdeg.oi3d will create 13 output files \n");
  fprintf(stderr,"           (Atlantic.quintdeg.oi3d.jan.nc .. Atlantic.quintdeg.oi3d.clim.nc in the directory /d1/HB3/gridded\n");
  fprintf(stderr,"     ex: -Oo/d1/HB3/gridded/Atlantic.quintdeg.oi3d.djf.nc will create a single output file.\n");
  fprintf(stderr," -P  list of properties to include in output file\n");
  fprintf(stderr,"        ex:  -Ppr/t90/sa\n");
  fprintf(stderr,"       -P (by itself) produces a list of available properties\n");
  fprintf(stderr, "\n\tOPTIONS:\n");
  fprintf(stderr," -N  minimum number of observations within specified search ellipse\n");     fprintf(stderr,"          default is  =%d  \n", minobs);
  fprintf(stderr," -S  maximum search radius (distance in km) \n");
  fprintf(stderr,"        Default is [%.0lf]\n", maxdist);
  fprintf(stderr," -V  file containing covariance parameters %s\n", parms_name);
  fprintf(stderr," -Z  file containing list of standard depths for output\n");
  fprintf(stderr," -h help...... prints this message. \n\n");

} /* end print_usage() */
/****************************************************************************/
/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
  char prop[7];
  int n, index;
  double reflev;

  if (*st == '\0') {
          print_prop_menu();
         exit(0);
  }
  
  /* Explicitly set mandatory properties in this order */
  
  prop_indx[0] = (int)PR;
  prop_indx[1] = (int)SA;
  prop_indx[2] = (int)T90;
  
  n = 3;
  do {
     if (*st == '/')
         ++st;
     prop[0]=prop[1]=prop[2]=prop[3]= prop[4]= prop[5] = prop[6]='\0';
     sscanf(st,"%[^'/']", prop);
     index = get_prop_indx(prop);
     
     if (index < 0)  {
       fprintf(stderr,"\n Unknown property '%s' specified in -P%s\n", prop, st);
       exit(1);
     }
     if (index > 100)  {
       fprintf(stderr,"\n -P%s : No need to specify variance and count properties, they are automatically output.  Quality flags are not output.\n", prop);
       continue;
     }
     
     switch (index) {
         case (int)T90:
         case (int)TE:
         case (int)TH:   /* fall through */
         case (int)TH9:
	     tindex = index;
         case (int)PR:
         case (int)SA:
	    break;
	    
	 default:
	    prop_indx[n] = index;   
     }

     st += strlen(prop);
     
     /* !**!  Special cases for properties */
     
     if ((prop_indx[n] == (int)S_)  || (prop_indx[n] == (int)PE) || (prop_indx[n] == (int)HT) ) {
        fprintf(stderr,"\nWARNING density and derivative properties like dynamic height are not appropriate  to average isopycnally.");
        fprintf(stderr,"\nCompute from averaged pr,t90,sa output by this module.\n");
	
     }
     

     if ((prop_indx[n] == (int) GE) || (prop_indx[n] == (int) GN)) {
        fprintf(stderr,"\nWARNING neutral density (gamma-n) is not an appropriate property to average isopycnally.");
        fprintf(stderr,"\nCompute from averaged pr,t90,sa output by this module.\n");
	
     }
     
     /* Don't count these properties... */
     if ( (prop_indx[n] == (int) DE)  || (prop_indx[n] == (int) GE) || (prop_indx[n] == (int) GN) ||(prop_indx[n] == (int)S_)  || (prop_indx[n] == (int)PE) || (prop_indx[n] == (int)HT) || prop_indx[n] > 100 )
        --n;
	
     if (index != (int)PR && index != tindex && index != (int)SA)
         ++n;
	 
   } while (*st == '/');
   
   
   return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
void cdf_construct(struct GRID_INFO *ginfo, int nprops, int *prop_indx,  char *dir, char *root, char *ext, int *nc_ids, int nargs, char **arglist)

   /* Opens all HydroBase netcdf output files (number of files is specified by global
      variable : nmonths) writing appropriate header,standard depths, and time bin info to each.
      Note that std_depth_init() must have been called prior to this function. The grid parameters 
      are supplied in ginfo. 
      Returns the ids associated with the open files in nc_ids. */
{
   struct CDF_HDR cdfhdr;
   int i, error;
   char *filename;
  
   cdfhdr.xmax = (float) ginfo->x_max;
   cdfhdr.xmin = (float) ginfo->x_min;
   cdfhdr.ymax = (float) ginfo->y_max;
   cdfhdr.ymin = (float) ginfo->y_min;
   cdfhdr.xincr = (float) ginfo->x_inc;
   cdfhdr.yincr = (float) ginfo->y_inc;
   
   cdfhdr.nx = ginfo->nx;
   cdfhdr.ny = ginfo->ny;
   cdfhdr.node_offset = ginfo->node_offset;
   
   cdfhdr.nz = NSTDLEVS;
   cdfhdr.nt = 1;
   cdfhdr.tmin = (int *) malloc(sizeof(int));
   cdfhdr.tmax = (int *) malloc(sizeof(int));
   cdfhdr.tmin[0] = 0;
   cdfhdr.tmax[0] = 9999;
   cdfhdr.counts_included = 3;
   cdfhdr.fill_value =  HBEMPTY;
   cdfhdr.mask_value =  HBMASK;
   strncpy(cdfhdr.x_units, "degrees", 8);
   strncpy(cdfhdr.y_units, "degrees", 8);
   strncpy(cdfhdr.z_units, "meters", 7);
   strncpy(cdfhdr.title,"HydroBase3", 11);
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
   
/* Open output files and write out some info so we can free up some memory... */
   
   for (i = 0; i < nmonths; ++i) {
      if (one_file)
          nc_ids[i] = cdf_init(root);  
      else {
         filename = (char *) calloc(500, sizeof(char));
         strcpy(filename, dir);
         strcat(filename, root);
         strncat(filename, ".",1);
         strcat(filename,cmonths[i]);
         strcat(filename, ext);
         nc_ids[i] = cdf_init(filename);  
         free(filename);
      }
       
      error = cdf_define(nc_ids[i], &cdfhdr, PREFILL, cdfhdr.counts_included);
      if (error)  exit(1);
   
      error = write_std_depths_cdf(nc_ids[i], &cdfhdr);
      error = write_std_depths_cdf(nc_ids[i], &cdfhdr);
      error = write_time_bins_cdf(nc_ids[i], &cdfhdr);
      error = write_lat_vector(nc_ids[i], &cdfhdr);
      error = write_lon_vector(nc_ids[i], &cdfhdr);
      
   
   } /* end for i */
   
   for (i = 0; i < nprops; ++i) {
      free((void *)cdfhdr.prop_id[i]);
      free((void *)cdfhdr.prop_units[i]);
   }
   free((void *) cdfhdr.prop_id);
   free((void *) cdfhdr.prop_units);
   free(cdfhdr.tmin);
   free(cdfhdr.tmax);
  
   return; 
} /*end cdf_construct() */ 

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
   int error, xpole, row, col;
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
	   tinfo->y_max = 90.0;
	   
        xpole = pointB(ginfo->y_min, ginfo->x_min, 180,(double) topobuf , TRUE, &lat_s, &dummy);
	if (xpole)
	   tinfo->y_min = -90.0;
	   
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

    if (tinfo->y_min  < -90)
	    tinfo->y_min = -90;
    if (tinfo->y_max  >= 90)
	    tinfo->y_max = 89.999;
    if (tinfo->lon0to360) {    
        if (tinfo->x_max >= 360)
	   tinfo->x_max = 359.999;
        if (tinfo->x_min < -360)
	   tinfo->x_max = -360;
    }
     return;
} /* end set_topo_info() */

/*****************************************************************************/

int subgrid(struct GRID_INFO *ginfo, struct GRID_INFO *newgrid, double lat0, double lon0, double lat1, double lon1, float maxdist)

/* determines size of grid which is maxdist beyond specified bounds
   but within bounds of ginfo. Fills in the fields of newgrid accordingly.
   The newgrid does not cross the pole.  Returns the number of grid nodes
   or 0 if an error occurs  */

{
   int error, ncols, nrows, col, row;
   int xpole;
   double lat, lon, lat_n, lat_s, lon_e, lon_w;


/* determine new bounds */

   
   xpole = pointB(lat1, lon0, 0, (double)maxdist, TRUE, &lat_n, &lon);
   if (xpole)
       lat_n = 90.0;
	   
   xpole = pointB(lat0, lon1, 180, (double)maxdist, TRUE, &lat_s, &lon);
   if (xpole)
       lat_s = -90.0;

   if (ABS(lat1) > ABS(lat0))  {  /* use most poleward bound */
       error = pointB(lat1, lon0, 270, (double)maxdist, TRUE, &lat, &lon_w);
       error = pointB(lat1, lon1, 90, (double) maxdist, TRUE, &lat, &lon_e);
    }
    else {  
       error = pointB(lat0, lon0, 270, (double)maxdist, TRUE, &lat, &lon_w);
       error = pointB(lat0, lon1, 90, (double)maxdist, TRUE, &lat, &lon_e);
    }
    
    /* adjust bounds to correspond to actual grid nodes */
    
    xy2ij_nochk(ginfo, lon_e, lat_n, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_e, &lat_n);
    
    xy2ij_nochk(ginfo, lon_w, lat_s, &col, &row);
    ij2xy_nochk(ginfo, col, row, &lon_w, &lat_s);
    
    if (lat_s < ginfo->y_min)
       lat_s = ginfo->y_min + 0.5 * ginfo->node_offset * ginfo->y_inc;
    if (lat_n > ginfo->y_max)
       lat_n = ginfo->y_max - 0.5 * ginfo->node_offset * ginfo->y_inc;
    if (lon_w < ginfo->x_min)
        lon_w = ginfo->x_min + 0.5 * ginfo->node_offset * ginfo->x_inc;
    if (lon_e > ginfo->x_max)
        lon_e = ginfo->x_max - 0.5 * ginfo->node_offset * ginfo->x_inc;
       
    
    newgrid->y_min = lat_s - 0.5 * ginfo->node_offset * ginfo->y_inc;
    newgrid->y_max = lat_n + 0.5 * ginfo->node_offset * ginfo->y_inc;
    newgrid->x_min = lon_w - 0.5 * ginfo->node_offset * ginfo->x_inc;
    newgrid->x_max = lon_e + 0.5 * ginfo->node_offset * ginfo->x_inc;
    
    newgrid->y_inc = ginfo->y_inc;
    newgrid->x_inc = ginfo->x_inc;
    newgrid->node_offset = ginfo->node_offset;
    newgrid->xgreenwich = (newgrid->x_min < 0) && (newgrid->x_max >= 0);
    newgrid->lon0to360 = (newgrid->x_min> 0) && (newgrid->x_max > 0);
    
    newgrid->ny = NINT((newgrid->y_max - newgrid->y_min) / newgrid->y_inc);
    newgrid->nx = NINT((newgrid->x_max - newgrid->x_min) / newgrid->x_inc);
    
    if (! newgrid->node_offset) {
       ++newgrid->ny;
       ++newgrid->nx;
    }
	  
    return(newgrid->ny * newgrid->nx);

} /* end subgrid() */


/*****************************************************************************/
void get_sig_profiles(int cdfid, struct CDF_HDR *cdfhdr, struct GRID_INFO *ginfo, int maxlev, float *bdepth, short *topo, struct GRID_INFO *topo_info, int t_indx, float **d, float **p, float **t, float **s, short **tcnt, float **sig0, float **sig1, float **sig2, float **sig3, float **sig4, int *nzlevs)
   /* Reads in profiles down to maxlev for each square in grid (described by ginfo),  
      computes sigmas and theta, defines a realistic bottom depth, and creates continuous profiles, 
      setting nzlevs at each square to reflect the number of valid zlevels found. 
      Sets mask and empty flags appropriately: testempty and testmask are global variables 
      
      cdfid  :  open HydroBase netcdf file
      cdfhdr :  info about file
      ginfo:  info about output grid
      maxlev:  index of deepest level to extract 
      bdepth: value of bottom depth at each square
      topo:  topography values
      topo_info:  info about topography grid
      t_indx :  index of temperature variable in netcdf file
      d,p,t,s,sig0,sig1...  property profiles for each square
      nzlevs : number of valid zlevels at each square
      
       These arrays correspond to dimensions of ginfo:
         bdepth[nsquares]
	 p[nsquares][maxlev] (and other property, count arrays)
	 nzlevs[nsquares]
      
      */
{
   int isq, nsq, row, col, row_cdf, col_cdf;
   int error, new_bd, tbin;
   int nz, ilev;
   float bd, *pin, *tin, *sin, *din;
   short *tnobs, topoval;
   double theLat, theLon;
   double *t_d, *p_d, *s_d;
   double *sig0_d, *sig1_d, *sig2_d, *sig3_d, *sig4_d;

  /* initialize these */
   
   tbin = 0;
   nsq = ginfo->nx * ginfo->ny;
   
   din = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   pin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   tin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   sin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   tnobs = (short *) get_memory((void *)NULL, cdfhdr->nz, sizeof(short));
      
   t_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   p_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   s_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   sig0_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   sig1_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   sig2_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   sig3_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   sig4_d = (double *) calloc(cdfhdr->nz, sizeof(double));
   
   error = read_cdf_depths(cdfid, din);  

   /* visit each square in the ginfo area */
   
   for (isq = 0; isq < nsq; ++isq) {
   
      /*prefill the work arrays */
      
      for (ilev = 0; ilev < maxlev; ++ilev) {
      
	 d[isq][ilev] = HBEMPTY;
         p[isq][ilev] = HBEMPTY;
         t[isq][ilev] = HBEMPTY;
         s[isq][ilev] = HBEMPTY;
         sig0[isq][ilev] = HBEMPTY;
         sig1[isq][ilev] = HBEMPTY;
         sig2[isq][ilev] = HBEMPTY;
         sig3[isq][ilev] = HBEMPTY;
         sig4[isq][ilev] = HBEMPTY;
	 if (tcnt != NULL) 
	    tcnt[isq][ilev] = 0;
      }
      
      if (nzlevs != NULL)
           nzlevs[isq] = 0;

      bd = cdfhdr->fill_value;
      
      error = sq2rc(isq, ginfo, &row, &col);
      error = ij2xy(ginfo, col, row, &theLon, &theLat);
      error = get_indices(cdfhdr, (float)theLat, (float)theLon, &row_cdf, &col_cdf);
      read_cdf_bottom_depth(cdfid, &bd, row_cdf, col_cdf, tbin);
      
      topoval = find_nearest_topo_val(theLat, theLon, topo, topo_info);
      new_bd = define_bottom(&bd, topoval);
      
      bdepth[isq] = bd;
      
      if (bd > testmask) {   
      
        /* set to mask and go to next square */
	
         for (ilev = 0; ilev < maxlev; ++ilev) {
	    d[isq][ilev] = cdfhdr->mask_value;
	    p[isq][ilev] = cdfhdr->mask_value;
	    t[isq][ilev] = cdfhdr->mask_value;
	    s[isq][ilev] = cdfhdr->mask_value;
	    sig0[isq][ilev] = cdfhdr->mask_value;
	    sig1[isq][ilev] = cdfhdr->mask_value;
	    sig2[isq][ilev] = cdfhdr->mask_value;
	    sig3[isq][ilev] = cdfhdr->mask_value;
	    sig4[isq][ilev] = cdfhdr->mask_value;
	 }
         if (nzlevs != NULL)
             nzlevs[isq] = 0;
         continue;
      }
      
      error = read_cdf_prop(cdfid,get_prop_mne((int)PR), pin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop(cdfid,get_prop_mne(t_indx), tin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop(cdfid,get_prop_mne((int)SA), sin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop_count(cdfid,get_prop_mne(t_indx),tnobs, row_cdf, col_cdf, 0, 0, cdfhdr->nz);

      /* copy input arrays to double arrays and eliminate missing values  */
      
      nz = 0;
      for (ilev = 0; ilev < maxlev; ++ilev) {
        if (pin[ilev] < testmask && pin[ilev] > testempty) {
	   d[isq][nz] =  din[ilev];
	   p_d[nz] =  (double)pin[ilev];
           t_d[nz] =  (double)tin[ilev];
	   s_d[nz] = (double)sin[ilev];
	   if (tcnt != NULL)
	      tcnt[isq][nz] = tnobs[ilev];
	   ++nz;
	}
      } 
      
      if (d[isq][nz-1] < 0)
             d[isq][nz-1] =  bdepth[isq];
      
      if (nzlevs != NULL)
             nzlevs[isq] = nz;

      /* make sure temperature is T90 */

      switch (t_indx) {
	 case (int) TH:
	    for (ilev = 0; ilev < nz; ++ilev)
	        t_d[ilev] = hb_theta(s_d[ilev],t_d[ilev], 0.0, p_d[ilev]);
		
		/* fall through */
		
         case (int) TE:
	 
            t68_to_t90(t_d, t_d, nz);
	    break;
	 default:
	    break;   
      } /* end switch */
      
     compute_sigma(0., nz, sig0_d, p_d, t_d, s_d);
     compute_sigma(1000., nz, sig1_d, p_d, t_d, s_d);
     compute_sigma(2000., nz, sig2_d, p_d, t_d, s_d);
     compute_sigma(3000., nz, sig3_d, p_d, t_d, s_d);
     compute_sigma(4000., nz, sig4_d, p_d, t_d, s_d);
     
     /* convert temperature profiles to TH9 profiles */
     
     compute_theta(nz, t_d, p_d, t_d, s_d);
     
     /* copy work arrays to output arrays */
     
     for (ilev = 0; ilev < nz; ++ilev) {
	    p[isq][ilev] = p_d[ilev];
	    t[isq][ilev] = t_d[ilev];
	    s[isq][ilev] = s_d[ilev];
	    sig0[isq][ilev] = sig0_d[ilev];
	    sig1[isq][ilev] = sig1_d[ilev];
	    sig2[isq][ilev] = sig2_d[ilev];
	    sig3[isq][ilev] = sig3_d[ilev];
	    sig4[isq][ilev] = sig4_d[ilev];
     
     }
     
   } /* end for isq */
   
   free(pin);
   free(din);
   free(tin);
   free(sin);
   free(tnobs);
   free(t_d);
   free(s_d);
   free(p_d);
   free(sig0_d);
   free(sig1_d);
   free(sig2_d);
   free(sig3_d);
   free(sig4_d);
   return;
   
} /* end get_sig_profiles() */
/************************************************************************/
/*****************************************************************************/
void get_monthly_profiles(int cdfid, struct CDF_HDR *cdfhdr, struct GRID_INFO *ginfo, int maxlev, float *bdepth, int t_indx, float **d, float **p, float **t, float **s, short **tcnt, float **sig0, int *nzout, struct MIXLAYER *mixlay)
   /* Reads in profiles down to maxlev for each square in grid (described by ginfo),  
      computes sigma0 and theta, evaluates a mixed layer and creates continuous profiles 
      of the upper depths. 
      Sets mask and empty flags appropriately: testempty and testmask are global variables 
      
      cdfid  :  open HydroBase netcdf file
      cdfhdr :  info about file
      ginfo:  info about output grid
      maxlev:  index of deepest level to extract 
      bdepth: value of bottom depth at each square
      t_indx :  index of temperature variable in netcdf file
      d,p,t,s,sig0 prooperty profiles for each square
      nzout : number of levels with data at each square
      
      These arrays correspond to dimensions of ginfo:
         bdepth[nsquares]
	 p[nsquares][maxlev] (and other property, count arrays)
	 mixlay[nsquares]
      
      */
{
   int isq, nsq, row, col, row_cdf, col_cdf;
   int error,  tbin;
   int nz, ilev;
   float *pin, *tin, *sin, *din;
   short *tnobs;
   double theLat, theLon;
   double *pwork, *twork, *swork, *sig0work, *dwork;

  /* initialize these */
   
   tbin = 0;
   nsq = ginfo->nx * ginfo->ny;
   
   pin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   din = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   tin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   sin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   dwork = (double *) get_memory((void *)NULL, cdfhdr->nz, sizeof(double));
   pwork = (double *) get_memory((void *)NULL, cdfhdr->nz, sizeof(double));
   twork = (double *) get_memory((void *)NULL, cdfhdr->nz, sizeof(double));
   swork = (double *) get_memory((void *)NULL, cdfhdr->nz, sizeof(double));
   sig0work = (double *) get_memory((void *)NULL, cdfhdr->nz, sizeof(double));
   tnobs = (short *) get_memory((void *)NULL, cdfhdr->nz, sizeof(short));
   
   
   read_cdf_depths(cdfid, din);

   /* visit each square in the ginfo area */
   
   for (isq = 0; isq < nsq; ++isq) {
   
      /*prefill the property arrays */
      
      for (ilev = 0; ilev < cdfhdr->nz; ++ilev) {
         pwork[ilev] = HBEMPTY;
         dwork[ilev] = HBEMPTY;
         twork[ilev] = HBEMPTY;
         swork[ilev] = HBEMPTY;
         sig0work[ilev] = HBEMPTY;
      }
      
      for (ilev = 0; ilev < maxlev; ++ilev) {
          d[isq][ilev] = HBEMPTY;
         p[isq][ilev] = HBEMPTY;
         t[isq][ilev] = HBEMPTY;
         s[isq][ilev] = HBEMPTY;
         sig0[isq][ilev] = HBEMPTY;
	 if (tcnt != NULL) 
	    tcnt[isq][ilev] = 0;
      }
      
      mixlay[isq].density = (float) HBEMPTY;
      mixlay[isq].depth = (float) HBEMPTY;
      mixlay[isq].theta = (float) HBEMPTY;
      mixlay[isq].salinity = (float) HBEMPTY;
      mixlay[isq].nobs = 0;

      error = sq2rc(isq, ginfo, &row, &col);
      error = ij2xy(ginfo, col, row, &theLon, &theLat);
      error = get_indices(cdfhdr, (float)theLat, (float)theLon, &row_cdf, &col_cdf);
      
      if (bdepth[isq] > testmask) {   
      
        /* set to mask and go to next square */
	 
         for (ilev = 0; ilev < maxlev; ++ilev) {
	    d[isq][ilev] = cdfhdr->mask_value;
	    p[isq][ilev] = cdfhdr->mask_value;
	    t[isq][ilev] = cdfhdr->mask_value;
	    s[isq][ilev] = cdfhdr->mask_value;
	    sig0[isq][ilev] = cdfhdr->mask_value;
	 }
	 mixlay[isq].density = cdfhdr->mask_value;
	 mixlay[isq].depth = cdfhdr->mask_value;
	 mixlay[isq].theta = cdfhdr->mask_value;
	 mixlay[isq].salinity = cdfhdr->mask_value;
	 mixlay[isq].nobs = 0;
         continue;
      }
      
      error = read_cdf_prop(cdfid,get_prop_mne((int)PR), pin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop(cdfid,get_prop_mne(t_indx), tin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop(cdfid,get_prop_mne((int)SA), sin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop_count(cdfid,get_prop_mne(t_indx),tnobs, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      din[cdfhdr->nz-1] = bdepth[isq];
 
      /* copy from float to double arrays and eliminate missing values  */
      
      nz = 0;
      for (ilev = 0; ilev < cdfhdr->nz; ++ilev) {
        if (pin[ilev] < testmask && pin[ilev] > testempty) {
	   dwork[nz] = (double) din[ilev];
	   pwork[nz] = (double) pin[ilev];
           twork[nz] = (double) tin[ilev];
	   swork[nz] = (double) sin[ilev];
	   if (nz < maxlev) {       
	      if (tcnt != NULL)
	         tcnt[isq][nz] = tnobs[ilev];
	   }
	   ++nz;
	}
      } 
 
         /* make sure temperature is T90 */
      if (nz > 0) {
         switch (t_indx) {
	    case (int) TH:
	       for (ilev = 0; ilev < nz; ++ilev)
	          twork[ilev] = hb_theta(swork[ilev],twork[ilev], 0.0, pwork[ilev]);
		
		/* fall through */
		
            case (int) TE:
               t68_to_t90(twork, twork, nz);
	       break;
	    default:
	    break;   
         } /* end switch */
      
         compute_sigma(0., nz, sig0work, pwork, twork, swork);
         compute_theta(nz, twork, pwork, twork, swork);
	 
       } /* end if nz */
       
       define_mixed_layer(theLat, sig0work, pwork, twork, swork, tnobs, nz, &mixlay[isq]);
     
      /* copy upper part of work arrays to property arrays */
      
       nzout[isq] = 0;
       for (ilev = 0; ilev < maxlev; ++ilev) {
	   d[isq][nzout[isq]] = (float)dwork[ilev];
	   p[isq][nzout[isq]] = (float)pwork[ilev];
           t[isq][nzout[isq]] = (float)twork[ilev];
	   s[isq][nzout[isq]] = (float)swork[ilev];
	   sig0[isq][nzout[isq]] = (float)sig0work[ilev];
	   if (pwork[ilev] > testempty && pwork[ilev] < testmask)
	     ++nzout[isq];
       } 

     
   } /* end for isq */
   
   free(din);
   free(pin);
   free(tin);
   free(sin);
   free(pwork);
   free(dwork);
   free(twork);
   free(swork);
   free(sig0work);
   free(tnobs);
   return;
   
} /* end get_monthly_profiles() */
/*****************************************************************************/
void get_x_profiles(int cdfid, struct CDF_HDR *cdfhdr, struct GRID_INFO *ginfo, int maxlev, float *bdepth, float **x, short **xcount, int propindx)

/* Reads property profiles from netcdf file into the x and xcount arrays. 
   bdepth array has already been filled with topography values
   in ()  */
{
   int isq, nsq, row, col, row_cdf, col_cdf;
   int error, tbin;
   int nz, ilev;
   float *xin;
   short *xnobs;
   double theLat, theLon, bpress;
   
  /* initialize these */
   
   tbin = 0;
   nsq = ginfo->nx * ginfo->ny;
   
   xin = (float *) get_memory((void *)NULL, cdfhdr->nz, sizeof(float));
   xnobs = (short *) get_memory((void *)NULL, cdfhdr->nz, sizeof(short));
   
   /* visit each square in the ginfo area */
   
   for (isq = 0; isq < nsq; ++isq) {
       
      if (bdepth[isq] > testmask) { 
        
      /* set to mask and go to next square */
	
         for (ilev = 0; ilev < maxlev; ++ilev) {
	    x[isq][ilev] = cdfhdr->mask_value;
	    if (xcount != NULL) 
	        xcount[isq][ilev] = 0;
	 }
	 continue;
      }
      
       /* prefill the output profile */
      
      for (ilev = 0; ilev < maxlev; ++ilev) {
         x[isq][ilev] = HBEMPTY;
	 if (xcount != NULL) 
	    xcount[isq][ilev] = 0;
      }
      
      error = sq2rc(isq, ginfo, &row, &col);
      error = ij2xy(ginfo, col, row, &theLon, &theLat);
      error = get_indices(cdfhdr, (float)theLat, (float)theLon, &row_cdf, &col_cdf);
      error = read_cdf_prop(cdfid,get_prop_mne(propindx), xin, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      error = read_cdf_prop_count(cdfid,get_prop_mne(propindx),xnobs, row_cdf, col_cdf, 0, 0, cdfhdr->nz);
      
      /* eliminate missing values  */
      
      nz = 0;
      for (ilev = 0; ilev < maxlev; ++ilev) {
         if (xin[ilev] < testmask && xin[ilev] > testempty) {
	   x[isq][nz] = xin[ilev];
	   if (xcount != NULL)
	      xcount[isq][nz] = xnobs[ilev];
	      
	   ++nz;
         }
      }
   }  /* end for isq */
   
   free(xin);
   free(xnobs);
   return;

}  /* end get_x_profiles() */

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
      
  if (*bd_addr > testempty && *bd_addr < testmask) {  
        
	/* bottom_depth is defined */
     
     if (seafloor > 0){           /* negative seafloor means land */
          if ( ((float)seafloor - *bd_addr) > 10) {
             *bd_addr = (float)(seafloor );
	     new_bd = 1;
	  }
	  
     }
  }  
  else {      /* bottom_de not defined */
     
    if (seafloor < 1) 
      *bd_addr = HBMASK;
    else 
      *bd_addr =  (float) seafloor;

  }
  return(new_bd);       
  
}  /* end define_bottom() */	
/*****************************************************************************/
void define_mixed_layer(double theLat, double *sig0, double *pr, double *th, double *sa, short *nobs, int nz, struct MIXLAYER *mlptr)
   /* Determines the depth, density, temperature and salinity of the mixed layer 
     and sets the fields of the structure accordingly. 
     mix_layer_delta is a global variable */
 
{
   int i, ilev, blev;
   double maxsig, p_bml;
   double pr1, th1, sa1, tsum, ssum, pdiff, psum;
   
   if (nz == 0)
      return;
      
   ilev = 0;
   while (ilev < nz && sig0[ilev] < testempty)
        ++ilev;
   
   if (pr[ilev] > 20)
      return;
   
   	
   mlptr->density = (float)sig0[ilev];
   mlptr->depth = (float)pr[ilev];
   mlptr->theta = (float)th[ilev];
   mlptr->salinity = (float)sa[ilev];
   mlptr->nobs = nobs[ilev];
	
   maxsig = mlptr->density + mix_layer_delta;
	
   p_bml = hb_linterp(maxsig, sig0, pr, nz);  /* pressure at bottom of mixed layer */
   
   if (p_bml > -9999) {
       mlptr->depth = (float) hb_depth(p_bml, theLat);
       
       blev = ilev;
       while (blev < nz && pr[blev] < p_bml)
          ++blev;
       
       --blev;
       
       if (blev == ilev) {
       
          mlptr->theta = (float) (th[ilev]);
          mlptr->salinity = (float)(sa[ilev]);
          mlptr->nobs = nobs[ilev];
	  return;
       }
       
       pr1 = pr[ilev];
       th1 = th[ilev];
       sa1 = sa[ilev];
       tsum = ssum = psum = 0.0;
       
       for (i = ilev+1; i <= blev; ++i){
          pdiff = pr[i] - pr1;
          tsum += (th1 + th[i]) * 0.5 * pdiff;
          ssum += (sa1 + sa[i]) * 0.5 * pdiff;
	  psum += pdiff;
	  pr1 = pr[i];
	  th1 = th[i];
	  sa1 = sa[i];
	  
	  if (nobs[i] >  mlptr->nobs)
	         mlptr->nobs = nobs[i];
       }
       /* compute average t90 and sa, then recompute density for these values */
       mlptr->theta = (float) (tsum / psum);
       mlptr->salinity = (float)(ssum / psum);
       hb_svan((double)mlptr->salinity, (double)mlptr->theta, 0.0, &maxsig);
       mlptr->density = (float) maxsig;	
   } 
   return;
 
} /* end define_mixed_layer() */
/*****************************************************************************/     
 void define_length_scales(struct NC_PARMS *pptr, float *xdist, float *ydist, double *stdlev, int nlevs, float bdepth)
/*  Sets the x and y length scales at each standard depth level 
    based on temperature gradients in ncparms. dtmax and dist are globally defined 
    variables.  */

{
   int ilev, n, nranges, bindx;
   float *xlen, *ylen;
   
   xlen = (float *)malloc(pptr->nz * sizeof(float));
   ylen = (float *)malloc(pptr->nz * sizeof(float));
   
   nranges = 8; 
   
   /* First, translate gradient to distance scales for each level in the parameter file */
   
   for (ilev = 0; ilev < pptr->nz; ++ilev) {
          
       xlen[ilev] = pptr->dTdx[ilev];
       
       if (xlen[ilev] < testmask) {
          if (xlen[ilev] > testempty ) {
             n = 0;
             while ( n < nranges && pptr->dTdx[ilev] > dtmax[n])  
               ++n;
             xlen[ilev] = dist[n]; 
	  }
	  else {                     /* if dTdx is missing, use length scale above */
	     if (ilev > 0)
	        xlen[ilev] = xlen[ilev-1];
	     else
	        xlen[ilev] = (float)DIST_DEFAULT;
	  }
       }
       
       ylen[ilev] = pptr->dTdy[ilev];
       if (ylen[ilev] < testmask) {
       
          if (ylen[ilev] > testempty) {
              n = 0;
             while ( n < nranges && pptr->dTdy[ilev] > dtmax[n])  
               ++n;
             ylen[ilev] = dist[n]; 
	  }
	  else {                        /* if dTdy is missing, use length scale above */
	     if (ilev > 0)
	        ylen[ilev] = ylen[ilev-1];
	     else
	        ylen[ilev] = (float)DIST_DEFAULT;
	  }
       }
   }
   
   /* Now assign the appropriate distance scale to each stdlev according to depth range */
   
   bindx = nlevs-1;
   ilev = 0;
   for (n = 1; n < pptr->nz; ++n) {
      while (ilev < bindx && stdlev[ilev] < pptr->depths[n]) {
          if (stdlev[ilev] > bdepth) {
	     xdist[ilev] = HBMASK;
	     ydist[ilev] = HBMASK;
	  }
	  else {
             xdist[ilev] = xlen[n-1];
	     ydist[ilev] = ylen[n-1];
	  }
	  ++ilev;
       }   
   }
   
   /* Continue for stdlevs below deepest depth in parameter file */
   
   while (ilev < bindx) {           
          if (stdlev[ilev] > bdepth) {
	     xdist[ilev] = HBMASK;
	     ydist[ilev] = HBMASK;
	  }
	  else {
             xdist[ilev] = xlen[n-1];  /* assign last parameter */
	     ydist[ilev] = ylen[n-1];
	  }
	  ++ilev;
   }
   
   /* Lastly, find where the non-standard bottom value falls in the stdlev array
      and use the length scales of the last stdlev above it */
   
   ilev = 1;
   while ( ilev < bindx && bdepth >= stdlev[ilev]) 
      ++ilev;
      
   xdist[bindx] = xdist[ilev-1];   
   ydist[bindx] = ydist[ilev-1];
 
   free(xlen);
   free(ylen);  
   return;

} /* end define_length_scales() */
/*****************************************************************************/

void get_prop_parms(struct NC_PARMS *parms_ptr, int pindx, double *stdlev, int nlevs, float bdepth, float *parms0, float *parms1, float *parms2)

/* Assigns values to the parms1 and parms2 arrays appropriate to the property at each stdlev. */

{
   int ilev, n, bindx, istart;
   float *p0, *p1, *p2;
   
   switch (pindx) {
       case (int)PR:
           p0 = parms_ptr->Pparm0;
           p1 = parms_ptr->Pparm1;
	   p2 = parms_ptr->Pparm2;
	   break;
       case (int)T90:
       case (int)TH9:
       case (int)TE:
       case (int)TH:
           p0 = parms_ptr->Tparm0;
           p1 = parms_ptr->Tparm1;
	   p2 = parms_ptr->Tparm2;
	   break;
       case (int)SA:
           p0 = parms_ptr->Sparm0;
           p1 = parms_ptr->Sparm1;
	   p2 = parms_ptr->Sparm2;
	   break;
	   
        default:
	   p0 = (float *)NULL;
	   p1 = (float *)NULL;
	   p2 = (float *)NULL;
	   return;
	   
   } /* end switch */

   
   bindx = nlevs-1;
   ilev = 0;
   
   istart = 1;                         /* find first non-empty parameter value */
   while (p1[istart] < testempty)
      ++istart;
   
   for (n = istart; n < parms_ptr->nz; ++n) {
   
      while (ilev < bindx && stdlev[ilev] < parms_ptr->depths[n] ) {
         if (stdlev[ilev] > bdepth) {
	    parms0[ilev] = HBMASK;
	    parms1[ilev] = HBMASK;
	    parms2[ilev] = HBMASK;
	 }
	 else {
	    parms0[ilev] = p0[n-1];
	    parms1[ilev] = p1[n-1];
	    parms2[ilev] = p2[n-1];
	 }
         ++ilev;
      }
   }
   
   while (ilev < bindx ) {
       if (stdlev[ilev] > bdepth) {
	    parms0[ilev] = HBMASK;
	    parms1[ilev] = HBMASK;
	    parms2[ilev] = HBMASK;
       }
       else {
	    parms0[ilev] = parms0[ilev-1];
	    parms1[ilev] = parms1[ilev-1];
	    parms2[ilev] = parms2[ilev-1];
       }
       ++ilev;
   }
   
  /* Lastly, find where the non-standard bottom value falls in the stdlev array
      and use the parameters of the last stdlev above it */
   
   ilev = 1;
   while ( ilev < bindx && bdepth >= stdlev[ilev]) 
      ++ilev;

   parms0[bindx] =  parms0[ilev-1]; 
   parms1[bindx] =  parms1[ilev-1]; 
   parms2[bindx] =  parms2[ilev-1]; 

   return;

}  /* end get_prop_parms() */
/*****************************************************************************/
int get_bounds(struct GRID_INFO *ginfo_ptr, int sq, float xdist,float ydist, double *lat0, double *lon0, double *lat1, double *lon1)
/* Returns minlat,minlon, maxlat, maxlon of area of length xdist,ydist around gridnode sq. 
   Returns 1 if pole was crossed, 0 if not.
  */
{
  int row, col, xpole;
  double lon, lat, dummy;
  double minlon, minlat, maxlon, maxlat;
  
  
  if (xdist < 0) xdist = maxdist;
  if (ydist < 0) ydist = maxdist;
  
  
  col = sq % ginfo_ptr->nx;
  row = (sq - col) / ginfo_ptr->nx;

  ij2xy(ginfo_ptr, col, row, &lon, &lat);
  minlon = lon - 0.5 * ginfo_ptr->node_offset * ginfo_ptr->x_inc;
  minlat = lat - 0.5 * ginfo_ptr->node_offset * ginfo_ptr->y_inc;
  maxlon = lon + 0.5 * ginfo_ptr->node_offset * ginfo_ptr->x_inc;
  maxlat = lat + 0.5 * ginfo_ptr->node_offset * ginfo_ptr->y_inc;
  xpole = pointB(minlat, lon, 180, (double) ydist, TRUE, lat0, &dummy);
  xpole += pointB(maxlat, lon, 0, (double)ydist, TRUE, lat1, &dummy);
  pointB(lat, minlon, 270, (double) xdist, TRUE, &dummy, lon0);
  pointB(lat, maxlon, 90, (double) xdist, TRUE, &dummy, lon1);

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

float find_isopycnal(float theDepth, int theSq, struct GRID_INFO *ginfo, float **z, float **sig, int *nz, float xdist, float ydist, float *bdepth)
  /* use weighted  interpolation to find an appropriate 
    isopycnal value at this level 
     theDepth, theSq refer to first-guess fields
     ginfo :  grid for z, sig, nz (first guess fields)
     xdist, ydist: for search ellipse
     bdepth: 
    */
{
  int xpole, nwsq, wrow, wcol, nobs;
  double minlat, minlon, maxlat, maxlon;
  double dummy, xout;
  double *xwork, *weights;
  struct GRID_INFO work_info;
  
  if (xdist < 0 || xdist > maxdist)
     xdist = maxdist;
     
  if (ydist < 0 || ydist > maxdist)
     ydist = maxdist;
     
  
     
  xpole = get_bounds(ginfo, theSq, (double)xdist, (double)ydist, &minlat, &minlon, &maxlat, &maxlon);

  nwsq = set_work_grid(ginfo, &work_info, &minlat, &minlon, &maxlat, &maxlon, xpole); 
  
  if (nwsq < 0)  {
     fprintf(stderr,"\nset_work_grid() returned nwsq = %d \n", nwsq);
     fprintf(stderr," theDepth = %.1f  first_guess square = %d\n", theDepth, theSq);
     exit(1);
  }

  xwork = (double *) get_memory(NULL,nwsq, sizeof(double));
  
  nobs = getsurf2d_depth(theDepth, ginfo, z, sig, (short **)NULL, nz, bdepth, minlat, minlon, maxlat, maxlon, xwork, (short *)NULL, &work_info);

  weights = (double *) calloc(nwsq, sizeof(double));
  wcol = work_info.nx / 2;
  wrow = work_info.ny / 2; 
  
  get_weights_c(weights, &work_info, 0, (double)xdist, (double)ydist, wcol, wrow);
  
  xout =  (float) weighted_mean( xwork, weights, (double) HBEMPTY, (double) HBMASK, nwsq);

 /*********  nn_interp2d(xwork, weights, (double) HBEMPTY, (double) HBMASK, wcol, wrow, &work_info, &xout, &dummy, &dummy, &nobs);*************/

   
   free(xwork);
   free(weights);
   return(xout);

} /* end find_isopycnal() */
/*****************************************************************************/
float weighted_stats(struct GRID_INFO *ginfo, double *xsurf, int nsq, float xdist, float ydist, float *xvar_ptr)
/* Estimates mean and variance at central node of gridded surface from surrounding values.  */
{
  int theRow, theCol;
  double *weights;
  double xbar, xvar;

  weights = (double *) calloc(nsq, sizeof(double));
  
  theCol = ginfo->nx / 2;
  theRow = ginfo->ny / 2; 
  
  get_weights_c(weights, ginfo, 0, (double)xdist, (double)ydist, theCol, theRow);
  
  xbar = weighted_mean(xsurf, weights, (double) HBEMPTY, (double) HBMASK, nsq);
  
  if ( xvar_ptr != NULL) {
     xvar = weighted_variance(xsurf, weights, xbar,(double) HBEMPTY, (double) HBMASK, nsq);
     *xvar_ptr = (float) xvar;
  }
  
  free(weights);
  return( (float) xbar);
  
} /* end weighted_stats() */

/*****************************************************************************/
int extract_surface(float surfval, int reflev, float ztarget, struct GRID_INFO *ginfo, float **xprofiles, float **yprofiles, short **ycounts, float **zprofiles, int *nz_profiles, float *dbotm, double theLat, double theLon, int theSq, float *xdist, float *ydist, struct GRID_INFO *surf_info, double **ysurf_addr, short **ysurf_count_addr, int *nsq_ptr, int checknobs)

/* Allocates memory for the surface and count arrays. Determines if there are enough points
   (minobs) within the specified length scales. If not, it enlarges the search ellipse and
   tests again until maxdist is reached in the x and y directions.  Then the surface array
   masking is set. 
   Returns the number of nodes with actual values, the surface and count arrays, 
   the (potentially adjusted) xdist, ydist, and the  number of nodes in the surface array.

*/

{
   int i, xpole, nsq, ngood, do_it;
   int *surfmask, *surfmask2;
   short *ysurf_count;
   double minlat, minlon, maxlat, maxlon;
   double *ysurf;

   xpole = get_bounds(ginfo, theSq, maxdist, maxdist, &minlat, &minlon, &maxlat, &maxlon);
   nsq = set_work_grid(ginfo, surf_info, &minlat, &minlon, &maxlat, &maxlon, xpole); 
   ysurf = (double *)get_memory(NULL, nsq, sizeof(double));
     for (i = 0; i < nsq; ++i) 
        ysurf[i] = HBEMPTY;
	
   ysurf_count = NULL;
   if (ycounts != NULL)
       ysurf_count = (short *)get_memory(NULL, nsq, sizeof(short));


    if (reflev < 0) {
        getsurf2d_depth(surfval, ginfo, xprofiles, yprofiles, ycounts, nz_profiles, dbotm, minlat, minlon, maxlat, maxlon, ysurf, ysurf_count, surf_info);
    }
    else  {
        ngood = getsurf2d(surfval, reflev, ztarget, xprofiles, yprofiles, ycounts, zprofiles, ginfo, nz_profiles, dbotm, ysurf, ysurf_count, surf_info, minlat, minlon, maxlat, maxlon);
     }
     
     surfmask = set_mask2d(ysurf, surf_info, theLat, theLon);
     surfmask2 = set_mask2d_dist(surf_info, *xdist, *ydist, theLat, theLon);

     ngood = 0;
     for (i = 0; i < nsq; ++i) 
         if (!surfmask[i] && !surfmask2[i]) 
             ++ngood;
     
     
     do_it = 1;
	
     while (checknobs && ngood < minobs && do_it) {
	   do_it = 0;
	   if (*xdist +50 <= maxdist+1) {
	      *xdist += 50;
	      do_it = 1;
	   }
	   if (*ydist +50 <= maxdist+1) {
	      *ydist += 50;
	      do_it = 1;
	   }
	   
	   if (do_it) {
	      free(surfmask2);
	     surfmask2 = set_mask2d_dist(surf_info, *xdist, *ydist,theLat, theLon);
	  }
	  ngood = 0;   
          for (i = 0; i < nsq; ++i) 
             if (!surfmask[i] && !surfmask2[i]) 
                 ++ngood;
     } /* end while */
     
     
     for (i = 0; i < nsq; ++i) {
	 if (surfmask[i] > 0 || surfmask2[i] > 0)
		   ysurf[i] = (double) HBMASK;
     }
     free(surfmask);
     free(surfmask2);

     /*  set the pointers and return */
     
     *ysurf_addr = ysurf;
     if (ysurf_count != NULL)
        *ysurf_count_addr = ysurf_count;
     *nsq_ptr = nsq;
     return(ngood);
	
}  /* end extract_surface() */

/*****************************************************************************/
int getsurf2d_depth(float zval, struct GRID_INFO *ginfo, float **zprofiles, float **yprofiles, short **ycount, int *nz, float *bdepth, double lat0, double lon0, double lat1, double lon1, double *ysurf, short *ysurf_cnt, struct GRID_INFO *surf_info)

/*  extracts ysurf from yprofiles along zval in zprofiles (depth) within the
    specified bounds. If ycount and ysurf_cnt are supplied (!=NULL) count info will be returned
       ginfo refers to zprofiles, yprofiles, ycount, nz
       surf_info is a subset of ginfo and refers to ysurf
       
       Note that bounds may specify the entire surf_info area or a subset 
 */
{
  int isq, isq_surf, nsq, row, col;
  int strow, endrow, stcol, endcol;
  int strow_surf, stcol_surf, surf_col, surf_row;
  int error;
  float yval;
  
  error = xy2ij(ginfo, lon0, lat0, &stcol, &strow);
  error += xy2ij(ginfo, lon1, lat1, &endcol, &endrow);
  
  if (error) {
     fprintf(stderr,"FATAL ERROR: bounds outside of first_guess grid in getsurf2d_depth(): \n    %.3lf %.3lf %.3lf %.3lf\n", lat0,lon0, lat1,lon1); 
     exit(1);
   }  
   
   error = xy2ij(surf_info, lon0, lat0, &stcol_surf, &strow_surf);
  if (error) {
     fprintf(stderr,"FATAL ERROR: bounds outside of first_guess work grid in getsurf2d_depth(): \n    %.3lf %.3lf %.3lf %.3lf\n", lat0,lon0, lat1,lon1); 
     exit(1);
   } 
    /* Initialize surf vals */
   nsq = surf_info->nx * surf_info->ny;
   for (isq = 0; isq < nsq; ++isq) {
      ysurf[isq] = (double)HBEMPTY;
      if (ysurf_cnt != NULL)
         ysurf_cnt[isq] = 0;
   }
   
   surf_row = strow_surf - 1;   
   for (row = strow; row < endrow; ++row) {
        ++surf_row;
	surf_col = stcol_surf -1;
      for (col = stcol; col < endcol; ++col) {
      
          ++surf_col;
          isq = row * ginfo->nx + col;
  
	  if (nz[isq] < 1)
	     yval = (double) HBEMPTY;
	  else if  (zval > bdepth[isq] )
	     yval = (double) HBMASK;
          else 
	     yval = interp_float(zval, zprofiles[isq], yprofiles[isq], nz[isq]);
	  
	  if (yval > -9998) {
	     isq_surf = surf_row * surf_info->nx + surf_col;
	     if (isq_surf >= nsq)  {
	           fprintf(stderr,"INDEXING ERROR in getsurf2d_depth()\n\n");
		   exit(1);
	     }
	     ysurf[isq_surf] = (double) yval;
	     if (ysurf_cnt != NULL  && yval < (double) testmask) {
	         ysurf_cnt[isq_surf] = interp_short(zval, zprofiles[isq],ycount[isq], nz[isq]);
		 
	     }
	  }
      }
   }
   
   return(1);

} /* end getsurf2d_depth() */
/*****************************************************************************/
/*****************************************************************************/
int getsurf2d(float isopyc, int reflev, float ztarget, float **sigprofiles, float **yprofiles, short **yprof_nobs, float **zprofiles, struct GRID_INFO *ginfo, int *nz, float *dbotm, double *ysurf, short *ysurf_nobs, struct GRID_INFO *surf_info, double lat0, double lon0, double lat1, double lon1)

/* Extracts y values from yprofiles along surface corresponding to isopyc in sigprofiles
   at each point in a lat/lon grid.  Returns the number of squares with actual values.  

     isopyc:  surface value to project y onto
     reflev:  for sigma surfaces 
     ztarget:  depth of level being interpolated
     sigprofiles, yprofiles, yprof_nobs, zprofiles:  
     ginfo:  refers to profiles, nlevs, pbot grids
     nz:   number of depth levels in each profile
     dbotm:  depth at seafloor
     
     surf_info:  refers to ysurf, ynobs grids
     ysurf, ysurf_nobs: arrays of projected values
     lat0,lon0,lat1,lon1: limits of region common to both surf_info and ginfo
      
     Note that bounds may specify the entire surf_info area or a subset 
*/

{
  int isq, ibotm, error, ilev, row, col;
  int strow, endrow, stcol, endcol;
  int strow_surf, stcol_surf, surf_col, surf_row;
  int nsq_surf, isq_surf, nvals, datagap;
  short ynobs;
  float yval, best_zval;
  float *z;
  
  error = xy2ij(ginfo, lon0, lat0, &stcol, &strow);
  error += xy2ij(ginfo, lon1, lat1, &endcol, &endrow);
  if (error) {
     fprintf(stderr,"FATAL ERROR: bounds outside of input grid in getsurf2d(): \n   %.3lf %.3lf %.3lf %.3lf\n", lat0,lon0, lat1,lon1); 
     exit(1);
   }  

   error = xy2ij(surf_info, lon0, lat0, &stcol_surf, &strow_surf);
  if (error) {
     fprintf(stderr,"FATAL ERROR: bounds outside of input grid in getsurf2d(): \n    %.3lf %.3lf %.3lf %.3lf\n", lat0,lon0, lat1,lon1); 
     exit(1);
   }  
   
   nvals = 0;
  
   surf_row = strow_surf - 1;   
   for (row = strow; row < endrow; ++row) {
        ++surf_row;
	surf_col = stcol_surf -1;
      for (col = stcol; col < endcol; ++col) {
      
          ++surf_col;
          isq = row * ginfo->nx + col;
	  ibotm = nz[isq]-1;
	  yval = HBEMPTY;
	  ynobs = 0;
	  
	  if (nz[isq] < 1) {   /* no values */
	     if (ztarget > dbotm[isq]+100)
	         yval = HBMASK;
	  }
          else  {  
	     
	     best_zval = check_all_crossings(isopyc, sigprofiles[isq], zprofiles[isq], nz[isq], ztarget);
	     
	     if (best_zval < testempty) {
	     
	         if ((zprofiles[isq][ibotm] > (dbotm[isq]- 50)) && (ztarget > zprofiles[isq][ibotm]) && (isopyc > sigprofiles[isq][ibotm]))
	            yval = HBMASK;  /* isopycnal runs into topo */
		    
	     }
	     
	     else if ( best_zval < (reflev - 1500))  /* too shallow compared to reflev */
	            yval = HBEMPTY;
	     
	     else {
	         /* find yval, and check for vertical datagaps around this depth level */
		 
		 z = zprofiles[isq];   /* set this pointer */
		 
		 yval = interp_float(best_zval, z, yprofiles[isq], nz[isq]);
		 
		 ilev = 0;
		 while (z[ilev] <= best_zval && ilev < nz[isq])
		     ++ilev;
		 
		 if (ilev == ibotm || ilev == 0)
		    datagap = 0;
		    
		 else if ( (z[ilev-1] == best_zval) || (z[ilev] == best_zval))
		    datagap = 0;
		    
		 else if (best_zval < 1001)
		    datagap = (z[ilev] - z[ilev-1]) > GAP_SHALLOW;
		 else
		    datagap = (z[ilev] - z[ilev-1]) > GAP_DEEP;
		    
		 if (datagap)
		    yval = HBEMPTY;
		    
		 if (!datagap) {
		    ++nvals;
		    if (ysurf_nobs != NULL)
		         ynobs = interp_short(best_zval, z, yprof_nobs[isq], nz[isq]);
		 }
		    
	     } /* end else */
          } /* end else */
	  
	  isq_surf = surf_row * surf_info->nx + surf_col;
	  
	  ysurf[isq_surf] = (double)yval;
	  if (ysurf_nobs != NULL)
	      ysurf_nobs[isq_surf] = ynobs;
  
      }  /* for col */
   } /* for row */

  
  return(nvals);

} /* end getsurf2d() */
/***************************************************/
float check_all_crossings(float isopyc, float *sigprof, float *zprof, int nz, float ztarget)
/*  interpolates to find all crossings of isopyc in sigprof and returns the depth of isopyc 
    closest to the target depth.  sigprof and zprof must be monotonically increasing */
{
  int j;
  float zval, diff, lastdiff, best_zval;
  
  j =0;
  lastdiff = HBMASK;
  best_zval = HBEMPTY;
  
  do {
  
     zval =interp_float(isopyc, &sigprof[j], &zprof[j], nz-j);
     if (zval > -999.0) {
        if ( (diff = ABS(zval - ztarget)) < lastdiff ) {
	   best_zval = zval;
	   lastdiff = diff;
	}
	
	while ((zprof[++j] <= zval) && (j < nz))
	   ;
     
     }
     
  }  while (zval > -9990. && (j < nz));
  
  if (best_zval >= 0) {  /* found an appropriate crossing */
  
     /* if at the top of a pycnostad, get as close as possible to ztarget */
     
     if (lastdiff > 300) { 
         zval = best_zval;
         diff = ztarget - best_zval;
     
         if (diff > 0)  {   /* zval is located above the ztarget */
            j = 0;
	    while ((zprof[j] <= zval) && (j < nz))
	       ++j;
	   
	     /* now j points to level just below best_zval */
	 
	     while ( (j < nz) && (zprof[j] < ztarget) && (ABS(sigprof[j]-isopyc) < 0.005)) 
	        ++j;
	    
	    /* now j points to level closest to ztarget to within 0.01 sigma units
	      of isopycnal being sought */
	 
	     if (j >= nz)   /* should never happen, but check anyway */
	        j = nz - 1;
		
	     best_zval = zprof[j];
	     
         } /* end if diff > 0*/
	 
	 /* if diff <= 0, best_zval is below target depth, so we can't get any closer*/
	 	 	 
      } /* end if lastdiff */
  }
  
  return(best_zval);
  
}  /* check_all_crossings() */

/***************************************************/
short interp_short(float xval, float *x, short *nobs, int npts)

/*   linearly interpolates to find the position of xval in array x, and returns
.    the corresponding value in the array, nobs.  If xval does not appear in array
.    x, the value of "flag" is returned.  This routine assumes that both x and nobs
.    arrays are continuous (no missing values)*/
{
    int k, flag;
    float v1,v2;
    
    flag = -9999;
    
    if (npts <= 0)
       return (flag);
       
    if (npts == 1) {
       if (xval == x[0] )
            return(nobs[0]);
	    
	    
	return (flag);
    }
    
    for (k = 0; k < npts-1; k++) {
        v1 = xval - x[k];
        v2 = xval - x[k+1];
        
        if (v1 == 0)
             return nobs[k];
        if (v2 == 0)
             return nobs[k+1];

        if ((v1 < 0.) && (v2 < 0.))
            continue;
        else if ((v1 > 0.) && (v2 > 0.))
            continue;
        else if ((v1 == v2) && (v1 != 0.))
            continue;
        else if ((v1 == 0.) && (v2 ==  0.))
            return nobs[k];
        else {
            if (nobs[k] > nobs[k+1])
	       return(nobs[k]);
	       
	    return(nobs[k+1]);
	}   
    }

/*  if execution falls through to this point, xval was not found in x array */
    return flag;
}

/***************************************************/
float interp_float(float xval, float *x, float *y, int npts)

/*   linearly interpolates to find the position of xval in array x, and returns
.    the corresponding value in the array, nobs.  If xval does not appear in array
.    x, the value of "flag" is returned.  This routine assumes that both x and nobs
.    arrays are continuous (no missing values)*/
{
    int k, flag;
    float v1,v2;
    
    flag = -9999;
    
    if (npts <= 0)
       return (flag);
       
    if (npts == 1) {
       if (xval == x[0] )
            return(y[0]);
	    
	    
	return (flag);
    }
    
    for (k = 0; k < npts-1; k++) {
        v1 = xval - x[k];
        v2 = xval - x[k+1];
        
        if (v1 == 0)
             return y[k];
        if (v2 == 0)
             return y[k+1];

        if ((v1 < 0.) && (v2 < 0.))
            continue;
        else if ((v1 > 0.) && (v2 > 0.))
            continue;
        else if ((v1 == v2) && (v1 != 0.))
            continue;
        else if ((v1 == 0.) && (v2 ==  0.))
            return y[k];
        else {
            if (y[k] > y[k+1])
	       return(y[k]);
	       
	    return(y[k+1]);
	}   
    }

/*  if execution falls through to this point, xval was not found in x array */
    return flag;
}
/*****************************************************************************/
int *set_mask2d_dist(struct GRID_INFO *hptr, float xlen, float ylen, double theLat, double theLon)

/*  Returns array of values set to 0 or 1 for nodes that are beyond the 
    distance ellipse specified by xlen, ylen).  
    
    Returns a pointer to mask or NULL if an error occurs and returns the number of squares that 
    have valid data.
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
/*****************************************************************************/
int *set_mask2d(double *x, struct GRID_INFO *hptr, double theLat, double theLon)

/*  Sets mask values to -1, 0, or 1 based on distribution of masked, 
    empty nodes in x (testmask and testempty are global variables).  
    
    Returns a pointer to mask or NULL if an error occurs. 
*/
{
   int nsq, sq, row, col, xmid, ymid, error, npts;
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
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
      
   } /* end while */
  
/*search east */
   row = ymid;
   col = xmid;
   dx = 0;
   mask_rest = 0;
   while (++dx <= xmid) {
      ++col;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
      
   } /* end while */

/*search north */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      ++row;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
   } /* end while */
   
/*search south */
   row = ymid;
   col = xmid;
   dy = 0;
   mask_rest = 0;
   while (++dy <= ymid) {
      --row;
      sq = row * hptr->nx + col;
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
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
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
      
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
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
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
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
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
      if (mask_rest)
         mask[sq] = 1;
      else if (mask[sq] == 1)  
	 mask_rest = 1;
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
/*****************************************************************************/
float *compute_residuals(struct GRID_INFO *ginfo1, double *xsurf1,int nsq1, struct GRID_INFO *ginfo2, double *xsurf2, int nsq2, short *nobs_addr)
/*  Subtracts xsurf1 from corresponding geographic locations in xsurf2.
    Creates a new array and returns it. */
{
   int error, isq1, isq2, row, col;
   double lon, lat;
   float *diff;
   
   diff = (float *) get_memory(NULL, nsq2,sizeof(float));
   *nobs_addr = 0;
   for (isq2 = 0; isq2 < nsq2; ++isq2) {  
	error = sq2rc(isq2, ginfo2, &row, &col);
	error += ij2xy(ginfo2, col, row, &lon, &lat);
	error += xy2ij(ginfo1, lon, lat, &col, &row);
	if (error) {
	   diff[isq2] = (float) HBEMPTY;
	}
	else {
	   isq1 = row * ginfo1->nx + col;
	   if (xsurf2[isq2] > testmask || xsurf1[isq1] > testmask) {
	      diff[isq2] = (float) HBMASK;
	   }
	   else if (xsurf2[isq2] < testempty || xsurf1[isq1] < testempty) {
	      diff[isq2] = (float)HBEMPTY;
	   }
	   else {
	      diff[isq2] = (float) (xsurf2[isq2] - xsurf1[isq1]);
	     ++(*nobs_addr);
           } 
	}		
    } /* end for isq */
    
    return(diff);

} /*end compute_residuals() */
/*****************************************************************************/
/*****************************************************************************/
int krig(double theLat, double theLon, double *xsurf, short *counts,int nwsq, struct GRID_INFO *ginfo, int imodel, float *params, float *theMean, float *theError, short *theCount)

/*  Estimates mean, error variance at a geographic location (theLat, theLon) from weighted
    mean of all valid surrounding points. Weights are computed from the estimated covariance of
    all pairs of points (Cij), and the covariance of all points with the specified location (Cx0),
    subject to the further constraint that the weights sum to unity:
       
       for i = 1:npts
          weight[i] = inv(C) * Cx0

         SumOf weight[i] = 1.0;
	 
    This contraint is achieved by setting up the system to minimize the error variance plus
    an additional term involving a Lagrange parameter.
     
  Input :
    theLat, theLon:  location of point to be estimated
    xsurf, counts:  surrounding observations in a 1-D array representing 2-D grid
                     On entry, every element of xsurf is set to an observed value, HBMASK or HBEMPTY
    nwsq:  size of xsurf, counts
    ginfo: pertains to the above grids
    imodel:  variogram model code
    params:  3 element array of variogram parameters for estimating covariance 
    
  Output:
    theMean:   
    theError:  sqrt(error_variance)
    theCount:  number of observations contributing to theMean, theError
    
    NOTE:  if input counts array is NULL, theCount is the number of non-empty, non-masked
            elements in resids[];
	    otherwise theCount is the sum of all non-empty, non-masked elements of counts[];
   
           testmask and testempty are global variables  
*/
{
  int i, j, row, col, isq, m, npts;
  int error, kilo = 1;                      /* distance in kilometers */
  short *zcount;
  float  *zin, **C, *Cx0, *W;
  double vh, phi, sill;
  double *zlat, *zlon;
  double *latv, *lonv, *lags, dparams[3];

  /* construct vectors of latitude and longitude corresponding to 
     row/col dimensions of ginfo */
     
  latv = (double *) malloc(ginfo->ny * sizeof(double));
  lonv = (double *) malloc(ginfo->nx * sizeof(double));
  
  for (i = 0; i < ginfo->ny; ++i) 
      latv[i] = ginfo->y_min + i * ginfo->y_inc;
  
  for (i = 0; i < ginfo->nx; ++i) 
      lonv[i] = ginfo->x_min + i * ginfo->x_inc;
/****************************************************/   
/*  store valid points, counts and their lat/lon */
  
  zlat = (double *) calloc(nwsq, sizeof(double));
  zlon = (double *) calloc(nwsq, sizeof(double));
  zin = (float *) calloc(nwsq, sizeof(float));
  zcount = NULL;
  if (counts != NULL)
     zcount = (short *) calloc(nwsq, sizeof(short));
  
  npts = 0;
  for (row = 0; row < ginfo->ny; ++row ) {
     for (col = 0; col < ginfo->nx; ++col) {
        isq = row *ginfo->nx + col;
	if (xsurf[isq] > testempty && xsurf[isq] < testmask) {
	   zin[npts] = (float) xsurf[isq];
	   zlat[npts] = latv[row];
	   zlon[npts] = lonv[col];
           if (zcount != NULL)
	     zcount[npts] = counts[isq];
	   ++npts;
	}
     }
  }
  free(latv);
  free(lonv);
  
/****************************************************/   
/* compute lags for all pairs of input points */
  
  lags = (double *) calloc(npts*npts, sizeof(double));
  
  for (i = 0; i < npts; ++i) {  
     for (j = i; j < npts; ++j) {
        if (i != j)
           lags[i*npts+j] = distance_c(zlat[j],zlon[j],zlat[i],zlon[i], kilo, &phi);
     }
  }
/****************************************************/   
/* construct covariance matrix  */
  
  m = npts + 1;  /* dimension of arrays for kriging */
  
  C = (float **) calloc(m, sizeof(float *));  
  for (i = 0; i < m; ++i) {
     C[i] = (float *) calloc(m, sizeof(float));
  }
  
  /* fill last column and row with 1 (for computing Lagrange parameter later) */
  
  for (i = 0; i < npts; ++i) {
     C[i][npts] = 1;
     C[npts][i] = 1;
  }
  C[npts][npts] = 0;   /* bottom right corner of matrix gets zero */
  
  /* Fill remaining elements with estimated covariance */
  
  dparams[0] = params[0];
  dparams[1] = params[1];
  dparams[2] = params[2];
  sill = dparams[0] + dparams[1];
  for (i = 0; i < npts; ++i) {
     for (j = i; j < npts; ++j) {
          if (i == j)
	     C[i][j] = sill;  /* C(0) = sill */
	  else {
             error = vgram_vals2(imodel, &lags[i*npts+j], 1, dparams, &vh);
	     C[i][j] = C[j][i] = (float)(sill - vh) ;
	  }
     }
  }
  
  free(lags);
  
/** Compute the inverse of this matrix */  
 
  error = invert(C, m);
  if (error)  {
      for (i=0; i<m; ++i)
         free(C[i]);
      free(C);
      free(zlat);
      free(zlon);
      free(zin);
      if (zcount != NULL)
         free(zcount);
      return(-1);
  }
/****************************************************/   
/* Compute lags for all points relative to theLat, theLon
   and get estimated covariance for each lag */ 

   lags = (double *) calloc(npts, sizeof(double));
   Cx0 = (float *) calloc(m, sizeof(float));  /* dimension is npts+1 */
   
   for (i=0; i<npts; ++i) {
       lags[i] = distance_c(zlat[i],zlon[i], theLat, theLon, kilo, &phi);
       error = vgram_vals2(imodel, &lags[i], 1, dparams, &vh);
       Cx0[i] = (float)(sill - vh);
   }
   
   Cx0[npts] = 1.0;   /* last element in this vector is assigned 1 
                       for computing a Lagrange parameter */
		       
   free(lags);	      
/****************************************************/   
/* Calculate weights:  W = inv(C) * Cx0   */
   
   W = (float *) calloc(m, sizeof(float));
   for (i=0; i < m; ++i) {
      W[i] = 0.0;
      for (j=0; j < m; ++j) 
          W[i] += C[i][j] * Cx0[j];
   }
   
/* Lastly compute mean, error variance and number of observations */
   
   *theCount = 0;
   *theError = 0.0;
   *theMean = 0.0;
   for (i=0; i<npts; ++i) {
      *theMean += W[i] * zin[i];
      *theError += W[i] * Cx0[i];
      if (zcount != NULL)
        *theCount += zcount[i];
   }
   *theError = sill - *theError - W[npts];  /* C(0)- sumof(w*Cx0) - lagrange */
   if (*theError < 0) {
     fprintf(stderr,"\nFATAL ERROR:  sign of error variance is negative in krig().  Something is wrong.\n\n");
     exit(1);
   }
      
   *theError = sqrtf(*theError);
   if (zcount == NULL)
      *theCount = npts;

/*  Finish */
   
   free(zlat);
   free(zlon);
   free(zin);
   free(Cx0);
   free(W);
   for (i = 0; i < m; ++i)
      free(C[i]);
   free(C);
   if (zcount != NULL)
      free(zcount);
   
   return(0);

} /* end krig() */
/*****************************************************************************/
int set_ml_props(struct MIXLAYER *ml_in, struct GRID_INFO *ginfo, double theLat, double theLon, int theSq, float *xdist, float *ydist,struct MIXLAYER *ml_out)
   /* Determines the mixed layer depth and density for the 
      specified location as a weighted mean of mixed layers that fall within 
      the specified distance scales.  */
{
   
   int isq, xpole, gsq, nsq, ngood, do_it, npts;
   int col, row, error;
   int *surfmask, *surfmask2;
   int *nobs;
   double dlat,dlon;
   double minlat, minlon, maxlat, maxlon;
   double *den, *dep, *sal, *thet, sig;
   float var;
   struct GRID_INFO surf_info;


   ml_out->density = HBEMPTY;
   ml_out->depth = HBEMPTY;
   ml_out->theta = HBEMPTY;
   ml_out->salinity = HBEMPTY;
   ml_out->derr = HBEMPTY;
   ml_out->terr = HBEMPTY;
   ml_out->serr = HBEMPTY;
   ml_out->nobs = 0;
   
   xpole = get_bounds(ginfo, theSq, (double)maxdist, (double)maxdist, &minlat, &minlon, &maxlat, &maxlon);
   nsq = set_work_grid(ginfo, &surf_info, &minlat, &minlon, &maxlat, &maxlon, xpole); 

   den = (double *)calloc( nsq, sizeof(double));
   dep = (double *)calloc( nsq, sizeof(double));
   sal = (double *)calloc( nsq, sizeof(double));
   thet = (double *)calloc( nsq, sizeof(double));
   nobs = (int *) calloc( nsq, sizeof(int));
 
   for (isq = 0; isq < nsq; ++isq ) {
   
     /* translate surf_info square to ginfo square */
   
     error = sq2latlon(isq, &surf_info, &dlat, &dlon );
     error += xy2ij(ginfo, dlon, dlat, &col, &row);
     if (!error) {
        gsq = row * ginfo->nx + col;
        den[isq] = (double) ml_in[gsq].density;
        dep[isq] = (double) ml_in[gsq].depth;
        sal[isq] = (double) ml_in[gsq].salinity;
        thet[isq] = (double) ml_in[gsq].theta;
        nobs[isq] = ml_in[gsq].nobs;
     }
     else {
        den[isq] = (double) HBEMPTY;
        dep[isq] = (double) HBEMPTY;
        sal[isq] = (double) HBEMPTY;
        thet[isq] = (double) HBEMPTY;
        nobs[isq] = 0;
     }
   }
   
     surfmask = set_mask2d(den, &surf_info, theLat, theLon);
     surfmask2 = set_mask2d_dist(&surf_info, *xdist, *ydist, theLat, theLon);
     ngood = 0;
     for (isq = 0; isq < nsq; ++isq) 
         if (!surfmask[isq] && !surfmask2[isq]) 
             ++ngood;
        
     npts = 0;
     for (isq = 0; isq < nsq; ++isq) {
	 if (surfmask[isq] > 0 || surfmask2[isq] > 0) {
		   den[isq] = (double) HBMASK;
 		   dep[isq] = (double) HBMASK;
 		   sal[isq] = (double) HBMASK;
 		   thet[isq] = (double) HBMASK;
	  }
	  else {
	     npts += nobs[isq];
	  }
     }
     free(surfmask);
     free(surfmask2);
     
     if (npts > 0) {
        ml_out->salinity = weighted_stats(&surf_info, sal, nsq, *xdist, *ydist, &var);
        if (var > 0 && var < testmask)
	    ml_out->serr = sqrtf(var);
        ml_out->theta = weighted_stats(&surf_info, thet, nsq, *xdist, *ydist, &var);
        if (var > 0 && var < testmask)
            ml_out->terr = sqrtf(var);
        ml_out->depth = weighted_stats(&surf_info, dep, nsq, *xdist, *ydist, &var);
        if (var > 0 && var < testmask)
           ml_out->derr = sqrtf(var);
        ml_out->nobs = npts;
	hb_svan((double)ml_out->salinity, (double)ml_out->theta, 0.0, &sig);
	ml_out->density = sig;
     }
     
     free(den);
     free(sal);
     free(thet);
     free(dep);
     free(nobs);
     return(ngood);
   
}  /* end set_ml_props */
/*****************************************************************************/   

void x_to_stddepth(struct MIXLAYER *mlptr, float *pin_f, float *xml, float *eml, short *nml, int nzml, float *xin_f, float *ein_f, short *nin_f,float *xout, float *eout, short *nout, float btmd, double *zout, int nstdlevs, double theLat)
   /* Interpolates xin_f array onto standard depths, using pin_f (pressure at each lev) 
      converted to depths. Upon entry, xml, eml and nml contain the mixed layer properties,
      xin_f, ein_f and nin_f hold the profile values resulting from the optimal interpolation
      below the mixed layer.   These are copied into double arrays, then interpolated back
      onto standard depth levels. The deepest value is inserted into the last element of the 
      output array.  Returns the arrays (xout, eout, nout) ready to be written to the netcdf file */
{

   double xval, dprev, diff;
   double *din_d,*ein_d, *xin_d; 
   short *nin_d;
   int n, nzin, lev, j, closest, ibotm, datagap;

   if (btmd > testmask) {
       for (lev = 0; lev < nstdlevs; ++lev) {
          if (xin_f[lev] < testmask) 
	     fprintf(stderr,"\n SOFTWARE WARnin_dG: non-masked value detected at masked square\n");
          xout[lev] = HBMASK;
	  eout[lev] = HBMASK;
	  nout[lev] = 0;
       }
       return;
   }
   
   /* Allocate memory to hold temporary arrays */
   
   xin_d = (double *) calloc(nstdlevs, sizeof(double));
   din_d = (double *) calloc(nstdlevs, sizeof(double));
   ein_d = (double *) calloc(nstdlevs, sizeof(double));
   nin_d = (short *) calloc(nstdlevs, sizeof(short));
   
   ibotm = nstdlevs-1;  /* index to deepest level */
   
   /* write mixed layer values into upper levels of array */
 
   if (xml != NULL) {
      lev = 0;
      while (lev < nzml) {
	 xin_d[lev] = (double) xml[lev];
	 din_d[lev] = zout[lev]; 
	 ein_d[lev] = (double) eml[lev];
	 nin_d[lev] = nml[lev];
	 ++lev;
      }
   }

   /* start with first level in xin_d that is deeper than depth of mixed layer.
   /* copy input arrays to temporary arrays weeding out any missing values
      and adjusting to ensure that depth is monotonically increasing */   
   
   dprev = nzml > 0? din_d[nzml-1] : din_d[0];
   nzin = nzml;
   lev = nzml;
   while (lev < nstdlevs) {
      if (pin_f[lev] > testempty && pin_f[lev] < testmask) {
        din_d[nzin] = hb_depth((double)pin_f[lev], theLat);
	xin_d[nzin] = (double) xin_f[lev];   
	ein_d[nzin] = (double) ein_f[lev];
	nin_d[nzin] = nin_f[lev];
	if (din_d[nzin] < dprev) {
	   din_d[nzin] = dprev + 1.0;
	   if (is_pr) {
	      if (nzin == 0)
	         xin_d[nzin] = dprev + 1.0;
	      else
	         xin_d[nzin] = xin_d[nzin-1];
	   }
	}
	dprev = din_d[nzin];
	
        ++nzin;
      }
      ++lev;
   }
 
   /* Now we have transferred input profiles to continuous arrays 
      ready to interpolate back onto standard depth levels */   
      
   /* prefill the output arrays with empty values */
         
   for (lev = 0; lev < nstdlevs; ++lev) {
         xout[lev] = (float)HBEMPTY;
	 eout[lev] = (float)HBEMPTY;
	 nout[lev] = 0;
   }
   
   /* create appropriate output arrays */
   
   if (nzin == 0) {             /* no observations */
      if (btmd > testempty) {
      
         lev = 0;   
	 while (lev < ibotm && zout[lev] <= btmd )
            ++lev;
	    
	 for (j = lev; j < ibotm; ++j) {
	     xout[lev] = (float)HBMASK;
	     eout[lev] = (float)HBMASK;
	 }
      }
      
      free(din_d);
      free(ein_d);
      free(xin_d);
      free(nin_d);
      return;
   }
   
   if (nzin == 1) { /* not enough values to interpolate, just output what we've got */
      diff = 9999;
      if (btmd > testempty) {
      
         diff = ABS(din_d[0] - btmd);  
         if (diff <= 20) {                    /* bottom value */
	    xout[ibotm] = (float)xin_d[0];
	    eout[ibotm] = (float)ein_d[0];
	    nout[ibotm] = (float)nin_d[0];
	 }
	          
         lev = 0;   
         while (lev < ibotm && zout[lev] <= btmd )
               ++lev;
	 for (j = lev; j < ibotm; ++j) {
	       xout[lev] = (float)HBMASK;
	       eout[lev] = (float)HBMASK;
	 }
      }
      
      if (diff > 20 ) {                   /* insert the value at proper output depth */
         while (lev < ibotm && din_d[0] <= zout[lev])
	     ++lev;
	 
	 --lev;    
	 xout[lev] = (float) xin_d[0];
	 eout[lev] = (float) ein_d[0];
	 nout[lev] = nin_d[0];
      }
      
      free(din_d);
      free(xin_d);
      free(ein_d);
      free(nin_d);
      return;
   }  /* end if nz == 1 */
   
/***************************************/
/* interpolate profile onto stddepths */
   
   lev = 0;
   while (lev < ibotm && zout[lev] < btmd) {
   
      if (zout[lev] < din_d[0]) {  /* stddepth shallower than top of profile */
        if ((din_d[0] - zout[lev]) < 11) {  /* close enough */
	   xout[lev] = xin_d[0];
	   eout[lev] = ein_d[0];
	   nout[lev] = nin_d[0];
	   if (is_pr)
	     xout[lev] = (float)hb_p80(zout[lev], theLat);
	}
      
      }  /* end if zout */
      else {
         j = 1;
         while (j < nzin && din_d[j] <= zout[lev])
            ++j;
	 
         /* din_d[j-1] and din_d[j] now bracket this depth level */
	 
         if (j < nzin) {  /* stdlev is within depth range of profile */
         
            closest = ABS(din_d[j]-zout[lev]) > ABS(zout[lev]-din_d[j-1])? j-1 : j;
	    diff = din_d[j] - din_d[j-1];
            datagap = diff > GAP_DEEP;
            if (zout[lev] < 1200)
               datagap = diff > GAP_SHALLOW;
        
	    if (datagap) {
	      xval = -9999.0;
	      if ( ABS(zout[lev]-din_d[closest]) < 10  ) {  
	         /* depth close enough to stdlev */
	          xout[lev] = (float)xin_d[closest];
	          eout[lev] = (float)ein_d[closest];
	          nout[lev] = nin_d[closest];
	      }
	    }
	    else {     /* no datagap */
	  
	      xval = xin_d[j-1];  /* use if zout < estimated depth of 1st observation */
	  
              if (zout[lev] > din_d[j-1]){
	          xval = hb_linterp(zout[lev], &din_d[j-1], &xin_d[j-1], 2);
	      }
	      else {
	         if (is_pr)
	            xval = hb_p80(zout[lev], theLat);
	      }
	       
	      if (xval > -9998) { 
	        xout[lev] = (float)xval;
	        eout[lev] = (float)ein_d[closest];
		nout[lev] = nin_d[closest];
              }
	    }
         }
      }  /* end else */
      ++lev; 
      
   } /* end while */
   
   /* mask the elements deeper than btmd, except the last element */
      
   for (j = lev; j < ibotm; ++j) {
     xout[j] = HBMASK;
     eout[j] = HBMASK;
   }

   if (ABS(din_d[nzin-1] - btmd) <= 100) {
       xout[ibotm] = (float) xin_d[nzin-1];
       eout[ibotm] = (float) ein_d[nzin-1];
       nout[ibotm] = nin_d[nzin-1];
   }
   
   free(din_d);
   free(xin_d);
   free(ein_d);
   free(nin_d);
   return;

} /* end x_to_stddepth() */
