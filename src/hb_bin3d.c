/*  hb_bin3d.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             December 2008
			     updated Dec 2009 for HB3M
			     updated to use T90 instead of T68
                             
			     
................................................................................
..........................................................................
. Computes an average profile for each property specified at each lat/lon
. interval from randomly spaced stations and outputs a netCDF file of hydro 
. properties gridded as a function of lat, lon, and depth.  The depth
. dimension consists of a series of standard depths, but these values also
. can be optionally specified by the user with -Z.  The averaging is done 
. along isopycnal surfaces and interpolated back onto the depth levels.  
. Lists of potential density (sigma values) are specified with the -S option
. Appropriate lists, which are customized for each 10-deg square of the ocean,
. are included in the HydroBase3 distribution (in subdirectory
. /lists).   The goal is to specify sigma values which 
. span the entire water column with appropriate vertical sampling. 
.
.  For the seasonal layer (upper 200 m) profiles can be sorted by months and output 
.  in separate files.  -M<moonthlist> where monthlist is some subset of
    0/1/2/3/4/5/6/7/8/9/10/11/12.  0 = clim, the rest correspond to jan ... dec.

.   To output separate .nc files for each month use -O<root_name> to which .month.nc
.   will be appended :   
    ex: -M1/2 -O7305.1deg  will produce 
         7305.1deg.jan.nc
         7305.1deg.feb.nc
	 
    ex: -M0/1/2/3/4/5/6/7/8/9/10/11/12 -O7305.1deg  will produce a full suite of 12 monthly
         files plus 7305.1deg.clim.cdf (containing all months, each month equally weighted).
    Alternatively, use -M<monthlist> and -N<filename> to generate one nc output file 
    incorporating profiles for the months specified (all profiles equally weighted) 
    
    ex:  -M12/1/2/3 -N7305.1deg.djfm.cdf 
	 

. Output values:  mean, nobs, std deviation (sqrt(variance)  for each property specified with -P
.
.  hb_bin3d is intended to be the preparatory step for input to the optimal
.  interpolation module, hb_oi3d. The grid intervals need not be the same for 
.  both:  hb_bin3d can be done at shorter grid scales.   
..........................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hydro_cdf.h"
#include "hydrobase.h"
#include "hb_paths.h"

#define  RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  6371.0    /* in km */

#define    UI    int
#define    PREFILL      0        /* turn off prefilling of cdf variables */
#define    ADD_PROP_STATS  3    /* indicates variance and # of obs info is included */
#define    PRINT_MSG  1          /* turn on message printing to stderr */
#define    MAX_TIME_BINS  1      /* time dimension stores just one matrix  */
#define    MIX_LAYER_DEF 0.02    /* default sigma range for mixed layer */
#define    BOTTOM_LAYER_DEF 0.0005  /* default sigma range for bottom boundary layer */
#define    SEASONAL_LAYER_DEF 200  /* depth of seasonally heated layer (meters) */
#define    NMONTHS  13            /* for handling seasonal layer */

/* input_file pathnames */

#define    EXTENT   ""
#define    DIR      ""
#define    MAXEXT   10  /* maximum # of different file extents specified */

/***************************************************************/
/*  define the standard sigma levels on which data will be averaged;
     sigma series is divided by reference levels ...*/

#define  MAXSIGLEVS 800    /* max size of sigma series -- arbitrary number */
#define  NREFLEVS  5       /* # of ref levs defining sigma values */

double siglevs[MAXSIGLEVS];   /* stores the sigma series  */
int  nsiglevs;                /* number of elements in full sigma series */

int ref_prop_id[NREFLEVS] = {(int)S0,(int)S1,(int)S2,(int)S3,(int)S4};
int pref_id[NREFLEVS] = {0, 1, 2, 3, 4};   /* pref / 1000 */
int nsig[NREFLEVS];       /* number of elements of sigseries for each ref lev */
int zmin[NREFLEVS] = {-1,500,1500,2500,3500}; /* depth ranges for ref levs */
int zmax[NREFLEVS] = {500,1500,2500,3500,99999};
double *sigptr[NREFLEVS];   /* points to station arrays of referenced sigma values */
char *s_root, *s_dir;

/***************************************************************/

/* globally referenced variables for input station data */

struct HYDRO_DATA sta;
struct HYDRO_HDR hdr;
double *pr, *de, *te, *sa;
double s_pref;
double pe_pref;
double ht_pref;

int  use_month[NMONTHS];   /* list of months to create separate seasonal layers */
int report_gaps;           /* switch turns on reporting of vertical datagaps */
double mix_layer_def;      /* option to specify mixed layer definition */
double bottom_layer_def;      /* option to specify bottom boundary layer definition */
double bottom_seasonal_layer;  /* bottom of the seasonal layer */
float lengthscale;            /* L for weight function = e^[-(d/L)^2] */

/*  global variables used in computing properties. */

float latitude, longitude;
int window, w_incr;     /* used to specify pr window for gradient properties*/
int ratio_done;
int ite, ipr, isa;      /* ith order of requested properties */
int tindex = (int)T90;  /* use T90 throughout*/
/************ structures to represent grid node ****************/

struct surfrec {
          double  density, depth;
          double *prop;
          double *var;
	  double *wghtsum;
              UI  n;
  struct surfrec *next;
};

struct deepestrec {
          double depth, pressure, thickness;
	  double sig_0, sig_1, sig_2, sig_3, sig_4;
          double *prop;
	  double *var;
	     int *count;
          struct deepestrec *next;
};

/* The grid node structure holds all information necessary
 * at each grid node. */
struct gridnode {
         double **prop;  /* accumulating sums for each property avg */
	 double **dprop, **dpropsq;  /* for computing variance */
	 double **weightsum;
         double  *d;
	 double  *dweightsum;
	 double **month_prop[NMONTHS];
	 double **month_dprop[NMONTHS], **month_dpropsq[NMONTHS];
	 double **month_wght[NMONTHS];
	 double *month_d[NMONTHS];
	 double *month_dwght[NMONTHS];
             UI **month_count[NMONTHS];
             UI *month_dnobs[NMONTHS];
             UI **count;
             UI  *nobs;
 struct deepestrec *deepest;
 struct surfrec  *mix_layer[NMONTHS];
};
/***************************************************************/
/*  prototypes for functions defined within this module... */

struct gridnode **alloc_grid(int, int, int, int);

int get_sig_series(int, FILE *,  double *);
int get_time_bins(int, int, struct CDF_HDR *);
int parse_p_option(char *, int *);
int parse_monthlist(char *);
int do_std_depth(struct gridnode *, int, int, int *, int, int, int, float **, float **, short **, float *);
int mixed_layer_vals(int, int, struct surfrec *, double *, double *, double *, double *);
int define_deepest_level(int, struct deepestrec *);
int get_prop(int);
int monotonic(double *, UI  *, double **, double **, UI  **, int, int, double *, UI *, double **, double **, UI **, int , double *, double *, UI *, double **, double **, UI **);

char *print_month(int);
FILE *get_standard_sigfile(int);

void free_grid(int, int, int, struct gridnode **);
void parse_s_option(char *, FILE **);
void insert_data(double *,int, double *, double *, double, UI *, int);
void insert_surf_vals(struct gridnode *,int *,int, int);
void insert_seasonal_layer(double *, int, double *, double *, double, UI *, int);
void insert_deepest_vals(struct gridnode *,int *,int);
void compute_avg(double *,double *, int);
void sum_levels(double *, double *, int, double *, double *, UI *, double *, double *, UI *);
void print_usage(char *);
void define_sigma_transitions(struct gridnode *, int);
void delete_list(struct deepestrec *);
void delete_surfrec_list(struct surfrec *);
void do_variance_prep(double *, int , double *,  double, double *, double *, int);
void do_variance_prep_monthly(double *, int , double *,  double, double *, double *, int);
void compute_variance(double *, double *, double *, int *, int);
void compute_ml_stats(struct surfrec *, struct surfrec *, int );
void write_to_file(int,struct gridnode ***, struct CDF_HDR *, int *, int, int, int, int, int);
void do_pycnostads(int, struct gridnode ***, struct CDF_HDR *, int *, int , int , int , int, int );

struct surfrec *get_surfrec(struct surfrec *, double, int);
struct surfrec *create_surfrec(int);
struct surfrec *copy_surfrec(struct surfrec *, int );
struct surfrec *get_density_rec(int, struct surfrec *, int);
struct surfrec *define_avg_mixed_layer(struct surfrec *, int);
struct deepestrec *create_deepestrec(int);

double interpolate(double, double *, double *, int);
double get_weight(float, float, float, float, float); 


/***************************************************************/

int main (int argc, char **argv)
{
   short   bopt, iopt, popt, nopt, oopt, no_distance_weight;
   int     curfile = 1, nfiles = 0, nmonths, imonth; 
   int n_extents, curext; 
   int     print_mess, xdateline; 
   int     nlevs, nprops = 0, prop_indx[MAXPROP];
   int     i, j, error, test_outcrop, prop_avail, index;
   int     nstdlevs, tmin, tmax, n_filled;
   int     include_stats;
   int    *tmpi;
   FILE   *sigfile[NREFLEVS];
   FILE   *z_file;
   FILE   *t_file;
   FILE   *infile;
   int     row, col, nrows, ncols, tbin, ntbins;
   int     in_bounds, out_of_bounds=0;
   char   *dir, *st;
   char   *nc_filename[NMONTHS], *nc_fileroot;
   char  **extent_list;
   int     nc_file;
   float  latc, lonc;
   struct CDF_HDR  h;
   struct gridnode **grid[MAX_TIME_BINS];
   struct surfrec *m[NMONTHS], *newrecptr, *ml_clim;
   double weight;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;
    extent_list[0] = EXTENT;
    bopt = iopt = popt = nopt = oopt = 0;
    z_file = t_file = NULL;
    for (j = 0; j < NREFLEVS; ++j)
       sigfile[j] = NULL;
    include_stats = ADD_PROP_STATS;
    error = 0;
    print_mess = PRINT_MSG;
    report_gaps = 0;
    xdateline = 0;
    window = 100;
    w_incr = 10;
    mix_layer_def = MIX_LAYER_DEF;
    bottom_layer_def = BOTTOM_LAYER_DEF;
    bottom_seasonal_layer = SEASONAL_LAYER_DEF;
    tmin = 0;
    tmax = 9999;
    no_distance_weight = 0;
    lengthscale = 100.0;   /* L in km for distance weighting = e^-(d/L)^2 */
    s_dir = (char *)LISTPATH;
    s_root = NULL;   
    for (j=0; j < NMONTHS; ++j) {
       use_month[j] = 0;
       nc_filename[j] = NULL;
    }
    nmonths = 0;
    ipr = ite = isa = -1;
   
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &h.xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &h.ymax) != 1);
                        
                        if (h.xmin > 0 && h.xmax < 0)
                           h.xmax += 360.;
                        if (h.xmax > 180)
                           xdateline = 1;
                           
                        break;

               case 'D':                   /* get input directory */
                        dir = &argv[i][2];
                        break;
			
               case 'E':       /* multiple file extents are supported */
                        extent_list[n_extents] = &argv[i][2];
			if (++n_extents == MAXEXT) {
                             fprintf(stderr,"\nToo many different file extents. Max allowed is %d\n", MAXEXT);
			 exit(1);
			}
                        break;

               case 'I':
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &h.xincr) == 1) ? 0 : 1;
                        h.yincr = h.xincr;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%f", &h.yincr) == 1) ? 0 : 1;
                              break;
                            }
                        }
                        
                        break;
               case 'L':
                        error = (sscanf(&argv[i][2],"%f", &lengthscale) == 1) ? 0 : 1;
                        break;
                        
               case 'M':
                        nmonths = parse_monthlist(&argv[i][2]);
                        break;
                        
               case 'N':                                /* for single .nc output file */
			nc_fileroot = &argv[i][2];
                        nopt = 1;
                        break;

               case 'O':                              /* for multiple .nc output files */
			nc_fileroot = &argv[i][2];
                        oopt = 1;
                        break;

               case 'P':
                        popt = 1;
                        nprops = parse_p_option(&argv[i][2], prop_indx);
                        break;

               case 'S':
                        parse_s_option(&argv[i][2], sigfile);
                        break;


               case 'U':
                        error = (sscanf(&argv[i][2],"%lf", &bottom_seasonal_layer) == 1) ? 0 : 1;
		        fprintf(stderr,"\nBottom of seasonal layer is %lf",bottom_seasonal_layer);
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
                        
                case 'X':
                        error = (sscanf(&argv[i][2],"%lf", &mix_layer_def) == 1) ? 0 : 1;
                        break;
                case 'Y':
                        error = (sscanf(&argv[i][2],"%d/%d", &tmin,&tmax) == 2) ? 0 : 1;
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
                           fprintf(stderr,"\nError opening stddep file: %s\n",&argv[i][2]);
                           exit(1);
                        }
                        break;

	       case 'h':  
	          print_usage(argv[0]);
	          exit(0);
		  
               default:
                        error = 1;

          }    /* end switch */

          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             print_usage(argv[0]);
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!bopt || !iopt || !nfiles || !popt ) {
       fprintf(stderr,"\nYou must specify input file(s), bounds, properties and gridspacing");
       fprintf(stderr," \n Use  %s -h to print usage", argv[0]);   
       exit(1);
   }
   if (!(nopt || oopt) ) {
       fprintf(stderr,"\nYou must specify output file(s)with one of -N or -O\n");
       fprintf(stderr," \n Use  %s -h to print usage", argv[0]);   
       exit(1);
   }
   
   if (nopt && oopt) {
       fprintf(stderr,"\nChoose either -N or -O options, but not both.\n");
       fprintf(stderr," \n Use  %s -h to print usage", argv[0]);   
       exit(1);
   }
   
   if (nmonths == 0) {   /* no months were specified, default is use all months ... */
         nmonths = NMONTHS;
	 for (i = 0; i < nmonths; ++i)
	   use_month[i] = 1;
   }

   if (nopt) {
      use_month[0] = 1;  
      fprintf(stderr,"\nOne output file will be generated for ");
   } 
   else  {
      fprintf(stderr,"\nSeparate output files will be generated for "); 
   }
   for (i=0; i < NMONTHS; ++i) {
      if (use_month[i])
         fprintf(stderr, "%s ",print_month(i));
   }

   if (lengthscale == 0)
      no_distance_weight = 1;

/* initialize global array of std_depths */

     nstdlevs = std_depth_init(z_file);

/* define values for sigma series ... */

   nsiglevs = 0;
   for (i = 0; i < NREFLEVS; ++i ) {
      if (sigfile[i] == NULL) {
         if (s_root == NULL) {
           fprintf(stderr,"\nNo sigma-series file for ref level %1d000 was specified\n", i);
           print_usage(argv[0]);
           exit(1);
	 }
	 sigfile[i] = get_standard_sigfile(pref_id[i]);
      }
      nsig[i] = get_sig_series(ref_prop_id[i], sigfile[i], &siglevs[nsiglevs]);
      nsiglevs += nsig[i];
      if (nsiglevs >= MAXSIGLEVS) {
          fprintf(stderr,"\n # of sigma levels exceeds space allocated");
          fprintf(stderr,"\n Recompile program with larger MAXSIGLEVS.\n");
          exit(1);
      }
   }
   fprintf(stderr,"\n total number of sigma levels: %d\n", nsiglevs);
   
   fprintf(stderr," %s will be used. \n", get_prop_descrip(tindex));


/* compute dimensions of grid and allocate space for computation ...*/

   nrows = NINT( ceil((double)((h.ymax - h.ymin) / h.yincr)));
   ncols = NINT(ceil((double)((h.xmax - h.xmin) / h.xincr)));
   h.node_offset = 1;       /* 0 for node grid, 1 for pixel grid */

   ntbins = get_time_bins(tmin, tmax, &h);

/*  set these cdf header values */
   h.nx = ncols;
   h.ny = nrows;
   h.nz = nstdlevs;
   h.nt = ntbins;

   fprintf(stderr,"\n allocating memory for grid ... " );

   for (tbin = 0; tbin < ntbins; ++tbin) {
     fprintf(stderr,"\n tbin = %i",tbin);
     grid[tbin] = alloc_grid(ncols, nrows, nsiglevs, nprops);
   }
   fprintf(stderr,"\n");

   if (no_distance_weight)
      fprintf(stderr,"Distance weighting will not be applied\n");
   else 
      fprintf(stderr,"Lengthscale for distance weighting [= e^-(d/L)^2]  is %.1lf km\n", lengthscale);
   
/* initialize these... */

   for (i = 0; i < MAXPROP; ++i) {
      sta.observ[i] = (double *) NULL;
      sta.count[i] = (double *) NULL;   
      sta.variance[i] = (double *) NULL;
      sta.quality[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;

/* loop for each input_file */

   do {
   
     curext = 0;
    do {
        if ((infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_mess)) 
          == NULL) {
          goto NEXTFILE;
        }

     
     /* loop for each station */

    while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;
	  
       in_bounds = 1;
       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) < 0) {
           /* if it is on the border, accept it */
           if (row == -1)  {
	       row = 0;
	       if (col == -1)
	           col = 0;
		else if (col == h.nx)
		   --col;
	   }
	   if (row == h.ny) {
	       --row;
	       if (col == -1)
	           col = 0;
		else if (col == h.nx)
		   --col;
	   }
	   if (col == -1) {
	       col = 0;
	       if (row == -1)
	           row = 0;
	       else if (row == h.ny)
	           --row ;
	   }
	   if (col == h.nx) {
	       --col;
	       if (row == -1)
	           row = 0;
	       else if (row == h.ny)
	           --row ;
	   }
	   if (col < 0 || col >= h.nx || row < 0 || row >= h.ny) {
/*              fprintf(stderr,"Out of bounds: %.3f %.3f\n", hdr.lat, hdr.lon); */
              ++out_of_bounds;
	      in_bounds = 0;
	   }
       }
       
       if (in_bounds) {
       
	 /* Check for presence of T90. */
         if (sta.observ[(int)T90] == NULL ) {
	   /* If T68 is available, automatically convert to T90
	    * and delete the T68 data. */
	   if (available(TE,&hdr)) {
	     free_and_alloc(&sta.observ[(int)T90], hdr.nobs);
	     t68_to_t90(sta.observ[(int)TE], sta.observ[(int)T90], hdr.nobs);
	     free(sta.observ[(int)TE]);
	     sta.observ[(int)TE] = NULL;
		   
	     for (i=0; i< sta.nprops; ++i) {
	       if (hdr.prop_id[i] == (int)TE)
		 hdr.prop_id[i] = (int)T90;
	     }    
	   }	      
	 }
         ratio_done = 0;
         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)T90];
         sa = sta.observ[(int)SA];
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) {
           fprintf(stderr, "\nThe requisite parameters: pr, t90, sa, de are not available at this station.\n");
           write_hydro_hdr(stderr, &hdr);
           continue;
         }

         /* compute referenced sigmas for this station... */

         free_and_alloc(&sta.observ[(int)S0], hdr.nobs);
         compute_sigma(0., hdr.nobs, sta.observ[(int)S0], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S1], hdr.nobs);
         compute_sigma(1000., hdr.nobs, sta.observ[(int)S1], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S2], hdr.nobs);
         compute_sigma(2000., hdr.nobs, sta.observ[(int)S2], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S3], hdr.nobs);
         compute_sigma(3000., hdr.nobs, sta.observ[(int)S3], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S4], hdr.nobs);
         compute_sigma(4000., hdr.nobs, sta.observ[(int)S4], pr, te, sa);

	 
         for (i = 0; i < NREFLEVS; ++i) {
            sigptr[i] = sta.observ[ref_prop_id[i]];
         }

         /* compute weight as a function of distance from center of grid node */
	 if (no_distance_weight ) 
	     weight = 1.0;
	 else {
	    get_lat_lon(&h, row, col, &latc, &lonc);
	    weight = get_weight(hdr.lat, hdr.lon, latc, lonc, lengthscale);
	 }

         /* compute other appropriate properties for this station ... 
	    Use potential temperatures in binning process instead of in situ*/
	    
         free_and_alloc(&sta.observ[(int)TH9], hdr.nobs);
         compute_theta( hdr.nobs, sta.observ[(int)TH9], pr, te, sa);

         for (i = 0; i < nprops; ++i) {
            test_outcrop = 0;  /* only set to 1 for depth */
	    
	    index =  prop_indx[i]; 
	    if (index == (int)T90 || index == (int)TE || index == (int)TH)
	          index = (int) TH9;
		  
            prop_avail = get_prop(index);
	    
            for (tbin = 0; tbin < ntbins; ++tbin) {
              if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {

                if (prop_avail) {
                  insert_data(sta.observ[index], hdr.nobs, grid[tbin][row][col].prop[i], grid[tbin][row][col].weightsum[i], weight, grid[tbin][row][col].count[i], test_outcrop);
		  insert_seasonal_layer(sta.observ[index], hdr.nobs, grid[tbin][row][col].month_prop[hdr.month][i], grid[tbin][row][col].month_wght[hdr.month][i], weight, grid[tbin][row][col].month_count[hdr.month][i], test_outcrop);
                }

              }
            } /* end for tbin */

         } /* end for */

      /* insert depth values and deal with the sea surface/bottom ... */

         test_outcrop = (de[0] <= 11.) ? 1 : 0;
         for (tbin = 0; tbin < ntbins; ++tbin) {
             if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {
                insert_data( de, hdr.nobs, grid[tbin][row][col].d, grid[tbin][row][col].dweightsum, weight, grid[tbin][row][col].nobs, test_outcrop);
		insert_seasonal_layer(de, hdr.nobs, grid[tbin][row][col].month_d[hdr.month], grid[tbin][row][col].month_dwght[hdr.month], weight, grid[tbin][row][col].month_dnobs[hdr.month], test_outcrop);
                insert_surf_vals(&grid[tbin][row][col], prop_indx, nprops, hdr.month);
                insert_deepest_vals(&grid[tbin][row][col], prop_indx, nprops);
             }
         }
       } /* end if in_bounds */
       

       /* free up space... */

       free_hydro_data(&sta);

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }

     fclose(infile);
     
NEXTFILE:
        ;
     
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );             


/******** end of first input phase *********/

   
   fprintf(stderr,"  computing averages ...\n");

/* for each gridnode, compute means at all sigma levels ... */

   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
            for (i = 0; i < nprops; ++i) {
               compute_avg(grid[tbin][row][col].prop[i], 
                           grid[tbin][row][col].weightsum[i], nsiglevs);
	    
            }
	    
           compute_avg(grid[tbin][row][col].d, grid[tbin][row][col].dweightsum, 
                        nsiglevs);
	    for (imonth = 1; imonth < NMONTHS; ++imonth) {
                compute_avg(grid[tbin][row][col].month_d[imonth], grid[tbin][row][col].month_dwght[imonth], 
                        nsig[0]);
		for (i = 0; i < nprops; ++i) {
                   compute_avg(grid[tbin][row][col].month_prop[imonth][i], 
                           grid[tbin][row][col].month_wght[imonth][i], nsig[0]);
		
		}
	    
	    }
         }
      }
   }
   
 /* Second pass through input files to compile variance statistics.
    At this point, we have mean values computed for each property and sigma level. */


   fprintf(stderr,"\n Compiling variance statistics ...");
   curfile = 1;
   
    do {
    
      curext = 0;
      do {
   
        if ((infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_mess)) 
           == NULL) {
         goto NEXTFILE1;
         }

   
      
     /* loop for each station */

     while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;

       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) == 0) {

         if (sta.observ[(int)T90] == NULL ) {
		if (available(TE,&hdr)) {
		   free_and_alloc(&sta.observ[(int)T90], hdr.nobs);
		   t68_to_t90(sta.observ[(int)TE], sta.observ[(int)T90], hdr.nobs);
		   free(sta.observ[(int)TE]);
		   sta.observ[(int)TE] = NULL;
		   
		   for (i=0; i< sta.nprops; ++i) {
		      if (hdr.prop_id[i] == (int)TE)
		         hdr.prop_id[i] = (int)T90;
		   }    
		}	      
	 }
         pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)T90];
         sa = sta.observ[(int)SA];
	 
	 ratio_done = 0;
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) 
           continue;
         

         /* compute referenced sigmas for this station... */

         free_and_alloc(&sta.observ[(int)S0], hdr.nobs);
         compute_sigma(0., hdr.nobs, sta.observ[(int)S0], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S1], hdr.nobs);
         compute_sigma(1000., hdr.nobs, sta.observ[(int)S1], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S2], hdr.nobs);
         compute_sigma(2000., hdr.nobs, sta.observ[(int)S2], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S3], hdr.nobs);
         compute_sigma(3000., hdr.nobs, sta.observ[(int)S3], pr, te, sa);
         free_and_alloc(&sta.observ[(int)S4], hdr.nobs);
         compute_sigma(4000., hdr.nobs, sta.observ[(int)S4], pr, te, sa);

         for (i = 0; i < NREFLEVS; ++i) {
            sigptr[i] = sta.observ[ref_prop_id[i]];
         }

         /* compute weight as a function of distance from center of grid node */
	 if (no_distance_weight ) 
	     weight = 1.0;
	 else {
	    get_lat_lon(&h, row, col, &latc, &lonc);
	    weight = get_weight(hdr.lat, hdr.lon, latc, lonc, lengthscale);
	 }
         /* compute other appropriate properties for this station ...
	    Use potential temperatures in binning process instead of in situ*/
	    
         free_and_alloc(&sta.observ[(int)TH9], hdr.nobs);
         compute_theta( hdr.nobs, sta.observ[(int)TH9], pr, te, sa);

         for (i = 0; i < nprops; ++i) {
	    index = prop_indx[i];
	    if (index == (int) T90 || index == (int)TE || index == (int)TH )
	        index = (int) TH9;
            prop_avail = get_prop(index);
	    
            for (tbin = 0; tbin < ntbins; ++tbin) {
              if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {

                if (prop_avail) {
		    test_outcrop = 0;
		    do_variance_prep( sta.observ[index], hdr.nobs, grid[tbin][row][col].prop[i], weight, grid[tbin][row][col].dprop[i], grid[tbin][row][col].dpropsq[i], test_outcrop);
 		    do_variance_prep_monthly( sta.observ[index], hdr.nobs, grid[tbin][row][col].month_prop[hdr.month][i], weight, grid[tbin][row][col].month_dprop[hdr.month][i], grid[tbin][row][col].month_dpropsq[hdr.month][i], test_outcrop);
                 
                }

              }
            } /* end for tbin */

         } /* end for i */
	 
       } /* end if */

       /* free up space... */

       free_hydro_data(&sta);

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }
     fclose(infile);
NEXTFILE1:
         ;

    } while (++curext < n_extents);
   } while (curfile++ < nfiles );             
   
   fprintf(stderr,"  computing std deviations...\n");

/* for each gridnode, compute variances ... */

   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
            for (i = 0; i < nprops; ++i) {
               compute_variance(grid[tbin][row][col].dprop[i], grid[tbin][row][col].dpropsq[i],
                           grid[tbin][row][col].weightsum[i], grid[tbin][row][col].count[i], nsiglevs);
               for (j = 1; j < NMONTHS; ++j) {
	           compute_variance(grid[tbin][row][col].month_dprop[j][i], grid[tbin][row][col].month_dpropsq[j][i], grid[tbin][row][col].month_wght[j][i], grid[tbin][row][col].month_count[j][i], nsig[0]);
		}
			   
            }
            define_sigma_transitions(&grid[tbin][row][col], nprops);
         }
      }
   }

   fprintf(stderr,"\n  constructing nc files ...\n");

/* construct the netcdf header  */

   h.nprops = nprops;       
   h.fill_value = (float) HBEMPTY;
   h.mask_value = (float) HBMASK;
   h.counts_included = 3;
   strncpy(h.x_units, "degrees", 8);
   strncpy(h.y_units, "degrees", 8);
   strncpy(h.z_units, "meters", 7);
   strncpy(h.title,"HydroBase3", 10);
   strcpy(h.command, *argv);
   for (i = 1; i < argc; ++i) {
      if (strlen(h.command) < 900) {
        strncat(h.command, " ", 1);
        strcat(h.command, argv[i]);
      }
   }
   h.prop_id = (char **) malloc(h.nprops * sizeof(char *));
   h.prop_units = (char **) malloc(h.nprops * sizeof(char *));
   for (i = 0; i < nprops; ++i) {
      h.prop_id[i] = (char *) calloc(5, sizeof(char));
      h.prop_units[i] = (char *) malloc(50);
      strcpy(h.prop_id[i], get_prop_mne(prop_indx[i]));
      strcpy(h.prop_units[i], get_prop_units(prop_indx[i]));
   }
   
   if ( h.prop_units[h.nprops - 1] == NULL) {
      fprintf(stderr, "\nUnable to malloc memory for the netcdf header.\n");
      exit(1);
   }
   
   /* Visit each gridnode, evaluate the mixed layers for each month and the average mixed-layer
      for all months (clim). Delete the mixed layer lists 0-12 and replace with a single record 
      representing the average for each month. Copy sigma-0 averages into seasonal arrays for 
      clim. During this stage, all temperatures are still in TH9, and will be converted to requested
      output temperatures later. */
   
   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
	    ml_clim = NULL;
            for (imonth = 1; imonth < NMONTHS; ++imonth) { 
               m[imonth] = define_avg_mixed_layer(grid[tbin][row][col].mix_layer[imonth], nprops);
	       if (m[imonth] != NULL) {
	           newrecptr = copy_surfrec(m[imonth], nprops);
	           delete_surfrec_list(grid[tbin][row][col].mix_layer[imonth]);
	           grid[tbin][row][col].mix_layer[imonth] = newrecptr;
		   m[imonth] = newrecptr;
	       
	           if (use_month[0]) {      /* insert it into linked list at ml_clim*/
		       iopt = 0;
		       if (nopt && use_month[imonth])  iopt = 1;
		       if (oopt) iopt = 1;
		       if (iopt) {
		  /* each month's mixed layer will be weighted equally to compute
		     an average for all months (n, weightsum are incremented) */

	                   newrecptr = get_surfrec(ml_clim, m[imonth]->density, nprops);
		           newrecptr->depth += m[imonth]->depth;
		           ++newrecptr->n ;
		           for (i = 0; i < nprops; ++i) {
		                newrecptr->prop[i] += m[imonth]->prop[i];
			        newrecptr->var[i] += m[imonth]->var[i];
		               ++newrecptr->wghtsum[i];
		           }
			   if (ml_clim == NULL) 
			       ml_clim = newrecptr;
		       }
		   }  /*if use_month[0] */
	       } /* if m[imonth] */
           }  /* for imonth */
            if (use_month[0]) {
	       m[0] = define_avg_mixed_layer(ml_clim, nprops);
	       if (m[0] != NULL) {
	           newrecptr = copy_surfrec(m[0], nprops);
	           delete_surfrec_list(ml_clim);
	           grid[tbin][row][col].mix_layer[0] = newrecptr;	          
	       }
	       /* Now copy sigma-0 arrays into seasonal layer arrays */
	       imonth = 0;
	       for (j=0; j < nsig[0]; ++j) {
	           grid[tbin][row][col].month_d[imonth][j] = grid[tbin][row][col].d[j];
	           grid[tbin][row][col].month_dnobs[imonth][j] = grid[tbin][row][col].nobs[j];
		   for (i = 0; i < nprops; ++i) {
		      grid[tbin][row][col].month_prop[imonth][i][j] = grid[tbin][row][col].prop[i][j];
		      grid[tbin][row][col].month_dpropsq[imonth][i][j] = grid[tbin][row][col].dpropsq[i][j];
		      grid[tbin][row][col].month_count[imonth][i][j] = grid[tbin][row][col].count[i][j];
		   }
	       }
            } /* end if use_month[0] */
         }
      }
   }   /* end for tbin */
  
/*  Output stage:   Temperatures are converted to requested properties in write_to_file() */

   if (use_month[0]) {
       imonth = 0;
       nc_filename[imonth] = (char *) calloc(1000, sizeof(char));
       strcpy(nc_filename[imonth], nc_fileroot); /* for nopt, filename is stored at fileroot */
       if (oopt) 
         strcat(nc_filename[imonth], ".clim.nc");
	 
       nc_file = cdf_init(nc_filename[imonth]);
       write_to_file(nc_file, grid, &h,  prop_indx, nstdlevs, ntbins, nrows, ncols, imonth);   
       cdf_close(nc_file);
       
   }
   if (oopt) {
      for (imonth = 1; imonth < NMONTHS; ++imonth) { 
          nc_filename[imonth] = (char *) calloc(1000, sizeof(char));
          strcpy(nc_filename[imonth], nc_fileroot);
          strcat(nc_filename[imonth], ".");
	  strcat(nc_filename[imonth], print_month(imonth));
	  strcat(nc_filename[imonth], ".nc");
	  
	  nc_file = cdf_init(nc_filename[imonth]);
          write_to_file(nc_file, grid, &h, prop_indx, nstdlevs, ntbins, nrows, ncols, imonth);   
          cdf_close(nc_file);
       } /* end for imonth */   
   }
   
   for (tbin = 0; tbin < ntbins; ++tbin) 
      free_grid(ncols, nrows, nprops, grid[tbin]);
      

/*  Now that memory has been freed up, search for vertical data
    gaps in isopycnally averaged grid and determine which result from
    pycnostads.  
    
    Strategy:
    Re-read all the station files again, this time averaging data
    on stddepth surfaces.  Then re-open the netcdf file, visit each
    grid node and see if there are vertical data gaps which may
    have resulted from weak vertical gradients.  Replace empty levels
    with isobarically averaged datapoints if they exist. Rewrite the
    updated property values to the netcdf file. */
    
   fprintf(stderr,"\nChecking for pycnostads....");
   fprintf(stderr,"\n    allocating memory for grid ... " );

   for (tbin = 0; tbin < ntbins; ++tbin) 
     grid[tbin] = alloc_grid(ncols, nrows, nstdlevs, nprops);
   

/* loop for each input_file */

   /*fprintf(stderr,"\n    summing ...");*/
   curfile = 1;
   
   do {
     curext = 0;
     do {
        if ((infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_mess)) == NULL) {
          goto NEXTFILE2;
        }

      
     /* loop for each station */

     while ((error = get_station(infile, &hdr, &sta)) == 0) {

       /* is station within bounds ?   */
       
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;

       if (get_indices(&h, hdr.lat, hdr.lon, &row, &col) == 0) {
	 /* Check for presence of T90 and if only T68 is present,
	  * then auto convert to T90 and kill the T68. */
	 if (sta.observ[(int)T90] == NULL) {
	   if (available(TE,&hdr)) {
	     free_and_alloc(&sta.observ[(int)T90], hdr.nobs);
	     t68_to_t90(sta.observ[(int)TE], sta.observ[(int)T90], hdr.nobs);
		   
	     for (i=0; i< sta.nprops; ++i) {
	       if (hdr.prop_id[i] == (int)TE)
		 hdr.prop_id[i] = (int)T90;
	     }
	   }    
	 }
	 pr = sta.observ[(int)PR];   /* set frequently used pointers */
         de = sta.observ[(int)DE];
         te = sta.observ[(int)T90];
         sa = sta.observ[(int)SA];

	 
	 ratio_done = 0;
         if (pr == NULL || te == NULL || sa == NULL || de == NULL) 
           continue;
         
         /* compute other appropriate properties for this station ... */
	 
	 /* No need for potential temperature here, because summing is done on isobaric
	    surfaces and the netcdf file has already been converted to in situ temperatures */

         for (i = 0; i < nprops; ++i) {
	    index = prop_indx[i];
	    
	    prop_avail = get_prop(index);

            for (tbin = 0; tbin < ntbins; ++tbin) {
              if (hdr.year >= h.tmin[tbin] && hdr.year <= h.tmax[tbin]) {

                if (prop_avail) {
		  /*fprintf(stderr,"\n Got to sum_levels...");*/
                  sum_levels(sta.observ[index], de, hdr.nobs, grid[tbin][row][col].prop[i], grid[tbin][row][col].dpropsq[i], grid[tbin][row][col].count[i], grid[tbin][row][col].month_prop[hdr.month][i], grid[tbin][row][col].month_dpropsq[hdr.month][i], grid[tbin][row][col].month_count[hdr.month][i]);
                }

              }
            } /* end for tbin */

         } /* end for i */
	 
       } /* end if */

       /* free up space... */

       free_hydro_data(&sta);

     }  /*end while !eof */ 
     if (error > 0) {
         report_status(error, stderr);
         exit(1);
   }

     fclose(infile);
NEXTFILE2:
         ;
     
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );             

/******** end of third input phase *********/

   for (imonth = 0; imonth < NMONTHS; ++ imonth) {
      if  (nc_filename[imonth] != NULL) {
         nc_file = cdf_update("", nc_filename[imonth], "", print_mess);
         do_pycnostads(nc_file, grid, &h, prop_indx, nstdlevs, ntbins, nrows, ncols, imonth);
	 cdf_close(nc_file);   
      }
   }
     
   fprintf(stderr,"\n End of %s\n", argv[0]);
  
   exit(0);
} /* end main */


/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n    Computes an average profile for each property specified");
   fprintf(stderr,"\nat each lat/lon gridnode from randomly spaced stations.");
   fprintf(stderr,"\nInput profiles are weighted by distance from gridnode center");
   fprintf(stderr,"\nunless weighting is turned off with -L0.");
   fprintf(stderr,"\n\nOutputs one or more netCDF files of property mean, std deviation, nobs gridded ");
   fprintf(stderr,"\nas a function of lat, lon, and depth. "); 
   fprintf(stderr,"\n\n    For the upper ocean seasonal and mixed layers, averages are computed for individual months. ");
   fprintf(stderr,"  These monthly averages are combined to determine an overall mean for all 12 months ");
   fprintf(stderr,"\nwhich is output in the .clim.nc files.  Use -N OR -O to control how output files are named and structured.");
   fprintf(stderr,"\n\n    The depth dimension consists of a series of standard depths [default],");
   fprintf(stderr," which can be user specified with -Z option. Averaging is done along isopycnal surfaces");
   fprintf(stderr,"\nand interpolated back onto the depth levels.  ");
   fprintf(stderr,"\n\n    Lists of density surfaces, customized for each 10-deg WMO square, are included with the HydroBase distribution (in /lists)");
   fprintf(stderr,"\nTo utilize the customized density surfaces, the 10-deg WMO square must be specified with -Sr<ms10> ");
   fprintf(stderr,"\nAlternate density (sigma) files can be specified with -S<pref>/<filename>, one for each reference pressure: 0,1,2,3,4");
   fprintf(stderr,"\n");
   
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -B/west/east/south/north [-D<dirname>] [-E<file_extent>] -I<gridspacing> [-L<lengthscale>] [-M<month_list>]  [-N<cdf_output_file> OR -O<root_name_of_output_files>]  -P<list_of_properties> -S<d|r|<ref_id>/<file_info>  [-U<seasonal layer_def>] [-X<mixed_layer_def>] -Y[yearmin/yearmax] [-Z<std_depth_file>]  \n", program);
   fprintf(stderr,"\n\n -B   specifies grid bounds");
   fprintf(stderr,"\n -N   name of sole netCDF output file. USE EITHER -N or -O ");
   fprintf(stderr,"\n -O   root name for multiple netCDF output files to which ${month}.nc will be appended (see -M below)");

   fprintf(stderr,"\n -I   specifies grid increment in degrees;  ex: -I0.5");
   fprintf(stderr,"\n          OR specify separate x,y increments with a slash");
   fprintf(stderr,"\n          to delimit xincr/yincr;  ex: -I2.0/0.5\n");
   fprintf(stderr,"\n -P   list of properties to grid -- usually just measured properties;");
   fprintf(stderr,"\n          ex:  -Ppr/t90/sa/ox");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties\n");
   fprintf(stderr,"\n -S   Specify sigma surface list files (if not specified, defaults to standard files in %s)", s_dir);
   fprintf(stderr,"\n        ####.sig0list, ####.sig1list, ####.sig2list, ####.sig3list, ####.sig4list ");
   fprintf(stderr,"\n       where #### is the  4-digit WMO square ");
   fprintf(stderr,"\n   -Sd <dirname> alternate directory where sigma files are stored");
   fprintf(stderr,"\n   -Sr <rootname> where rootname = 10-deg WMO square corresponding to region being binned");
   fprintf(stderr,"\n   ex: -Sd/home/ruth/siglists -Sr7305  specifies the directory where 7305.sig0list 7305.sig1list, etc will be found ");
   fprintf(stderr,"\n\n OPTIONS:");
   fprintf(stderr,"\n -D  specifies directory for input data files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/\n ");

   fprintf(stderr,"\n -E : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n -L lengthscale (L, in km) for distance weighting =  e ^[-(dist/L)^2]");
   fprintf(stderr,"\n     -L0 turns off distance weighting (all points contribute equally to mean)" );
   fprintf(stderr,"\n          default is L =%.1f  km", lengthscale);
   fprintf(stderr,"\n -M  list of months to separately evaluate seasonal layer (upper 200 m): ");
   fprintf(stderr,"\n      0/1/2/3/4/5/6/7/8/9/10/11/12  0 = clim, the rest correspond to jan...dec.");

   fprintf(stderr,"\n      To output separate .nc files for each month use -O<root_name> ");
   fprintf(stderr,"\n      to which .month.nc will be appended : "); 
   fprintf(stderr,"\n          ex: -M1/2 -O7305.1deg  will produce ");
   fprintf(stderr,"\n               7305.1deg.jan.nc");
   fprintf(stderr,"\n               7305.1deg.feb.nc");
   fprintf(stderr,"\n          ex: -M0/1/2/3/4/5/6/7/8/9/10/11/12 -O7305.1deg will produce ");
   fprintf(stderr,"\n               a full suite of 12 monthly files plus");
   fprintf(stderr,"\n               7305.1deg.clim.nc (containing all months, each month equally weighted)");
   fprintf(stderr,"\n      Alternatively, use -M<monthlist> and -N<filename> to generate one output .nc file");
   fprintf(stderr,"\n      incorporating profiles for the months specified (months are equally weighted)"); 
   fprintf(stderr,"\n          ex:  -M12/1/2/3 -N7305.1deg.djfm.nc ");
   fprintf(stderr,"\n     The mixed layer may be deeper than the seasonal layer for a particular month, this condition is accommodated.");
   fprintf(stderr,"\n     Below the mixed and seasonal layers, the full set of profiles");
   fprintf(stderr,"\n     is used to compute means and variances.");
   fprintf(stderr,"\n -U  Set depth of the upper ocean seasonally heated layer (in meters).");
   fprintf(stderr,"\n          Separate averages for each month will be computed over this depth range");
   fprintf(stderr,"\n          default: %.1lf ", bottom_seasonal_layer);
   fprintf(stderr,"\n -X  specify the definition of the mixed layer as ");
   fprintf(stderr,"\n          sigma(sea surface) + this value.  default: %f ", MIX_LAYER_DEF);
   fprintf(stderr,"\n          ex:  -X.02");
   fprintf(stderr,"\n -Y  optional minyear/maxyear to constrain time interval of observations.");
   fprintf(stderr,"\n -Z  file containing list of standard depths.");
   fprintf(stderr,"\n          Values MUST be monotonically INCREASING.\n");
   fprintf(stderr,"\n -h  help ... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */
/****************************************************************************/

void parse_s_option(char *str, FILE **file_ptr)
   /* get the sigma values defining the isopycnal surfaces */
{
  int j;

  switch (*str) {
    case 'r':
           s_root = ++str;
	   return;
	   
    case 'd':
           s_dir = ++str;
	   return;
    case '0':
           j = 0;
           break;
    case '1':
           j = 1;
           break;
    case '2':
           j = 2;
           break;
        case '3':
           j = 3;
           break;
        case '4':
           j = 4;
           break;
        default:
           fprintf(stderr,"\n Error parsing -S option\n");
           fprintf(stderr,"USAGE: -S<ref_lev_id>/<file_name>");
           fprintf(stderr,"\n ex: -S4/sig4levs.list\n");
           exit(1);
           break;

    }  /* end switch */

    file_ptr[j] = fopen(str,"r");
    if (file_ptr[j] == NULL) {
      fprintf(stderr,"\nError opening %s \n", str);
      exit(1);
    }

   return;
}   /* end parse_s_option() */
/*****************************************************************************/
int parse_monthlist(char *st)
   /* Reads integer months separted by slashes and sets a flag in the global array
      use_month[].  Returns the number of months listed */  
{
   int imonth, count=0;

   if (*st == '\0')
      return(0);
      
   if (*st == '/') ++st;
      
   while (st != NULL) {
      if (*st != '\0') {
         if (sscanf(st,"%d", &imonth) != 1) {
            fprintf(stderr,"\nError parsing month list from %s\n", st);
            return(-1);
         }
        use_month[imonth] = 1;
         ++count;
         if ((st = strchr(st,'/')) != NULL)
 	    ++st;
      }   
   }
   return(count);
   
} /* end parse_monthlist */
/****************************************************************************/

char *print_month(int m_index)
/* Returns a string associated with month index */
{
      switch (m_index) {
         case 0:
	    return("clim");
	 case 1:
	    return("jan");
	 case 2:
	    return("feb");
	 case 3:
	    return("mar");
	 case 4:
	    return("apr");
         case 5:
	    return("may");
	 case 6:
	    return("jun");
	 case 7:
	    return("jul");
	 case 8:
	    return("aug");
	 case 9:
	    return("sep");
 	 case 10:
	    return("oct");
	 case 11:
	    return("nov");
	 case 12:
	    return("dec");
	 default:
	    return((char *)NULL);
     } /* end switch */

} /* end print_month() */
/*****************************************************************************/
int parse_p_option(char *st, int *prop_indx)
{
  char prop[6];
  int n;
  double reflev;

  if (*st == '\0') {
          print_prop_menu();
         exit(0);
  }

  /* n is counter for number of properties */
  n = 0;
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
     
     if (prop_indx[n] == (int)PR)
         ipr = n;
     if (prop_indx[n] == (int)T90 || prop_indx[n] == (int)TE || prop_indx[n] == (int)TH9 || prop_indx[n] == (int)TH) {
       /* Global variable ite is set to -1
	* when program is initialized, gets
	* set to sane value at first detection
	* of requested property.  Use this to
	* detect multiple temp requests, which
	* are not permitted.*/
       if (ite < 0) {
	 /* Based on whichever temperature that is
	  * requested, we now use that temperature.*/
	 tindex = prop_indx[n];
	 ite = n;
       }
       else {
	 fprintf(stderr,"\n Specify a single temperature property for output!");
         fprintf(stderr,"\n ie: t90(recommended), te, th9 or th.\n");
	 exit(1);
       } 
       if (tindex != (int)T90) {
	 fprintf(stderr,"\n T90 is the recommended temperature output variable.");
	 fprintf(stderr,"\n You have specified %s\n", prop);
       }
     }
     if (prop_indx[n] == (int)SA)
       isa = n;
     st += strlen(prop);
     
     /* !**!  Special cases for properties */
     if ((prop_indx[n] == (int)S_)  || (prop_indx[n] == (int)PE) || (prop_indx[n] == (int)HT) ) {
       fprintf(stderr,"\n WARNING: density and derivative properties");
       fprintf(stderr,"\n like dynamic height are not appropriate  to average isopycnally.");
       fprintf(stderr,"\nCompute from averaged pr,te,sa output by hb_bin3d.\n");
     }
     if ((prop_indx[n] == (int) GE) || (prop_indx[n] == (int) GN)) {
       fprintf(stderr,"\n WARNING: neutral density (gamma-n) is not an");
       fprintf(stderr,"\n appropriate property to average isopycnally.");
       fprintf(stderr,"\nCompute from averaged pr,te,sa output by hb_bin3d.\n");
     }
     /* Don't count these properties... */
     if ( (prop_indx[n] == (int) DE)  || (prop_indx[n] == (int) GE) || (prop_indx[n] == (int) GN) ||(prop_indx[n] == (int)S_)  || (prop_indx[n] == (int)PE) || (prop_indx[n] == (int)HT) || prop_indx[n] > 100 )
        --n;
      
     ++n;
  } while (*st == '/');
  fprintf(stderr,"\n Found %i properties requested.",n);
  return (n);
}  /* end parse_p_option() */

/*****************************************************************************/
struct gridnode **alloc_grid(int nx, int ny, int nz, int nprops)
{
   int i, j, k, n;
   struct gridnode **g;

   g = (struct gridnode **) malloc(ny * sizeof(struct gridnode *));
   if (g == NULL) {
      fprintf(stderr,"\nUnable to allocate memory for grid.\n");
      exit(1);
   }
   for (i = 0; i < ny; ++i) {
         g[i] = (struct gridnode *) malloc(nx * sizeof(struct gridnode));
         if (g[i] == NULL) {
           fprintf(stderr,"\nUnable to allocate memory for grid[%d]\n", i);
           exit(1);
         }
         for (j = 0; j < nx; ++j) {

            g[i][j].prop = (double **) malloc(nprops * sizeof(double *));
            if (g[i][j].prop == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop\n", i, j);
                exit(1);
            }
            g[i][j].weightsum = (double **) malloc(nprops * sizeof(double *));
            if (g[i][j].weightsum == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum\n", i, j);
                exit(1);
            }
            g[i][j].dprop = (double **) malloc(nprops * sizeof(double *));
            if (g[i][j].dprop == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].dprop\n", i, j);
                exit(1);
            }

             g[i][j].dpropsq = (double **) malloc(nprops * sizeof(double *));
            if (g[i][j].dpropsq == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].dpropsq\n", i, j);
                exit(1);
            }
           g[i][j].count = (UI **) malloc(nprops * sizeof(UI *));
            if (g[i][j].count == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count\n", i, j);
                exit(1);
            }

            g[i][j].d = (double *) calloc(nz, sizeof(double));
            if (g[i][j].d == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].d\n", i, j);
               exit(1);
            }
            g[i][j].dweightsum = (double *) calloc(nz, sizeof(double));
            if (g[i][j].dweightsum == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].dweightsum\n", i, j);
               exit(1);
            }
            g[i][j].nobs = (UI *) calloc(nz, sizeof(UI));
            if (g[i][j].nobs == NULL) {
               fprintf(stderr,"\nUnable to allocate memory for g[%d][%d].nobs\n", i, j);
               exit(1);
            }

            g[i][j].deepest = (struct deepestrec *)NULL;

	    
	    for (n = 0; n < NMONTHS; ++n) {  /* monthly arrays */
	      g[i][j].month_prop[n] = (double **)calloc(nprops, sizeof(double *));
	      g[i][j].month_dprop[n] = (double **)calloc(nprops, sizeof(double *));
	      g[i][j].month_dpropsq[n] = (double **)calloc(nprops, sizeof(double *));
	      g[i][j].month_wght[n] = (double **)calloc(nprops, sizeof(double *));
	      g[i][j].month_count[n] = (UI **)calloc(nprops, sizeof(UI *));
	      g[i][j].month_d[n] = (double *)calloc(nsig[0], sizeof(double));
	      g[i][j].month_dwght[n] = (double *)calloc(nsig[0], sizeof(double));
	      g[i][j].month_dnobs[n] = (UI *)calloc(nsig[0], sizeof(UI));
              g[i][j].mix_layer[n] = (struct surfrec *) NULL;
	    }

            for (k = 0; k < nprops; ++k) {

              g[i][j].prop[k] = (double *) calloc(nz, sizeof(double));
              if (g[i][j].prop[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].prop[%d]\n", i, j, k);
                exit(1);
              }
              g[i][j].weightsum[k] = (double *) calloc(nz, sizeof(double));
              if (g[i][j].weightsum[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].weightsum[%d]\n", i, j, k);
                exit(1);
              }

              g[i][j].dprop[k] = (double *) calloc(nz, sizeof(double));
              if (g[i][j].dprop[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].dprop[%d]\n", i, j, k);
                exit(1);
              }
              g[i][j].dpropsq[k] = (double *) calloc(nz, sizeof(double));
              if (g[i][j].dpropsq[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].dpropsq[%d]\n", i, j, k);
                exit(1);
              }
              g[i][j].count[k] = (UI *) calloc(nz, sizeof(UI));
              if (g[i][j].count[k] == NULL) {
                fprintf(stderr,"\nUnable to allocate memory for grid[%d][%d].count[%d]\n", i, j, k);
                exit(1);
              }
	      
	      for (n = 0; n < NMONTHS; ++n) {  /* monthly arrays */
	         g[i][j].month_prop[n][k] = (double *)calloc(nsig[0], sizeof(double));
	         g[i][j].month_dprop[n][k] = (double *)calloc(nsig[0], sizeof(double));
	         g[i][j].month_dpropsq[n][k] = (double *)calloc(nsig[0], sizeof(double));
	         g[i][j].month_wght[n][k] = (double *)calloc(nsig[0], sizeof(double));
	         g[i][j].month_count[n][k] = (UI *)calloc(nsig[0], sizeof(UI));
	      }

              
            } /* end for */
         } /* end for */
   }  /* end for */

   return g;
} /* end alloc_grid */
/*****************************************************************************/
void free_grid(int nx, int ny, int nprops, struct gridnode **gridptr)
  /* frees all the memory for ny * nx struct gridnodes */
{
   int row, col, k, n;
   struct gridnode *g;
   
   for (row = 0; row < ny; ++row) {
      for (col = 0; col < nx; ++col) {
         g = &gridptr[row][col];
         for (k = 0; k < nprops; ++k) {
           if (g->prop[k] != NULL)
	     free((void *)g->prop[k]);
           if (g->dprop[k] != NULL)
	     free((void *)g->dprop[k]);
           if (g->dpropsq[k] != NULL)
	     free((void *)g->dpropsq[k]);
           if (g->weightsum[k] != NULL)
	     free((void *)g->weightsum[k]);
           if (g->count[k] != NULL)
	     free((void *)g->count[k]);
	     
	   for (n = 0; n < NMONTHS; ++n) {  /* monthly arrays */
              if (g->month_prop[n][k] != NULL)
	         free((void *)g->month_prop[n][k]);
              if (g->month_dprop[n][k] != NULL)
	         free((void *)g->month_dprop[n][k]);
              if (g->month_dpropsq[n][k] != NULL)
	         free((void *)g->month_dpropsq[n][k]);
              if (g->month_wght[n][k] != NULL)
	         free((void *)g->month_wght[n][k]);
              if (g->month_count[n][k] != NULL)
	         free((void *)g->month_count[n][k]);
	   }
	   
         }
	 free((void *)g->prop);
	 free((void *)g->dprop);
	 free((void *)g->dpropsq);
	 free((void *)g->weightsum);
	 free((void *)g->count);
	 free((void *)g->nobs);
	 free((void *)g->d);
	 free((void *)g->dweightsum);
	 for (n = 0; n < NMONTHS; ++n) {  
	   free((void *)g->month_prop[n]);
	   free((void *)g->month_dprop[n]);
	   free((void *)g->month_dpropsq[n]);
	   free((void *)g->month_wght[n]);
	   free((void *)g->month_d[n]);
	   free((void *)g->month_dwght[n]);
	   free((void *)g->month_dnobs[n]);
	   free((void *)g->month_count[n]);
	   delete_surfrec_list(g->mix_layer[n]);
	 }
	 delete_list(g->deepest);
      } /* end for col */
      
      free((void *) gridptr[row]);
   } /* end for row */
   
   free((void *) gridptr);
   
   return;
}  /* end free_grid() */
/****************************************************************************/
FILE * get_standard_sigfile(int sindex)
   /* constructs name of standard index file, opensand returns a ptr to it.
      Global variables  s_dir, s_root are used */
{
  char *st;
  FILE *fptr;
  
	 /* construct name of standard siglist file */
	 st = (char *)calloc(500, sizeof(char));
	 strcpy(st, s_dir);
	 strcat(st, s_root);
	 strncat(st, ".sig", 5);
	 switch (sindex){
	    case 0:
	       strncat(st, "0", 2);
	       break;
	    case 1:
	       strncat(st, "1", 2);
	       break;
	    case 2:
	       strncat(st, "2", 2);
	       break;
	    case 3:
	       strncat(st, "3", 2);
	       break;
	    case 4:
	       strncat(st, "4", 2);
	       break;
	    } /* end switch */
	 strncat(st, "list", 5);
	 fptr = fopen(st,"r");
	 if (fptr == NULL) {
           fprintf(stderr,"\nError opening sigma-list file: %s\n",st );
           fprintf(stderr,"Common errors: spelling, pathname separators (slashes) \n");
	   exit(1);
	 }
	 free(st);
	 return(fptr);

} /* end get_standard_sigfile()*/
/*****************************************************************************/
int get_sig_series(int ref_id, FILE *fptr, double *siglist)
  
    /*  the file will be read for (min, max, incr) triplets
       from which a sigma series will be generated and inserted at siglist.
       The number of sigma values is returned. 
       
       arguments:
       ref_id:   defines the sigma variable: (int) enum property 
         fptr:   pointer to file containing list of sigma values OR nil 
      siglist:   starting addr of array to insert values 
   */
{
   double  min, max, incr;
   double  next;
   int i, count;
   int bigcount;


   bigcount = nsiglevs;

     count = 0;
     while( fscanf(fptr,"%lf %lf %lf", &min, &max, &incr) == 3) {
         if (++bigcount >= MAXSIGLEVS) {
           return (bigcount);
         }
         *siglist = min;
         ++count;
         while ( ( next = *siglist + incr) <= max) {
           ++count;
           if (++bigcount >= MAXSIGLEVS) { /* avoid a SEG FAULT */
              return (bigcount);
           }
           *(++siglist) = next; 
         }
         ++siglist;
     }

     fclose(fptr);
     return (count);

} /* end get_sig_series */
/*****************************************************************************/
int get_time_bins(int minyr, int maxyr, struct CDF_HDR *hptr)
    /* Sets time bin min/max values in the CDF header */
{

   hptr->tmin = (int *) malloc( sizeof(int));
   hptr->tmax = (int *) malloc( sizeof(int));

   hptr->tmin[0] = minyr;        
   hptr->tmax[0] = maxyr;

   return (1);

 } /* end get_time_bins */

/*****************************************************************************/
double get_weight(float lat1, float lon1, float lat2, float lon2, float L)
   /* Returns weight as a function of distance between points:
         weight = e ^[-(d/L)^2]
  */
{  
   double dx, dy, dist;
   
   dx = (lon2 - lon1)* RADperDEG * cos(RADperDEG*.5 *(lat1+lat2)) * EarthRadius ;
   dy = (lat2 - lat1) * RADperDEG * EarthRadius;
   dist = sqrt(dx*dx + dy*dy);
   dist = dist /L;
   dist = dist * dist;
   
   return (exp( -dist));
}	 
/*****************************************************************************/
int get_prop(int index)

{
   int prop_avail, main_props_avail, j;

    pr = sta.observ[(int)PR];   /* set frequently used pointers */
    de = sta.observ[(int)DE];
    te = sta.observ[(int)T90];
    sa = sta.observ[(int)SA];
	 
    prop_avail = available((enum property) index, &hdr);
    
    if (prop_avail)
         return(prop_avail);
       
       
    main_props_avail = !(pr== NULL || te == NULL || sa == NULL || de == NULL);
    prop_avail = main_props_avail;
    

    if (main_props_avail) {
	 
            switch ((enum property) index) {
	    
               case OX:  
                    if (available(O2, &hdr)) {
                      free_and_alloc(&sta.observ[(int)OX], hdr.nobs);
                      for (j=0; j < hdr.nobs; ++j) {
                        sta.observ[OX][j] = ox_kg2l(sta.observ[(int)O2][j], pr[j], te[j], sa[j]);
                      }
                    }
		    else
		       prop_avail = 0;
                  break;
               
               case O2:  
                    if (available(OX, &hdr)) {
                      free_and_alloc(&sta.observ[(int)O2], hdr.nobs);
                      for (j=0; j < hdr.nobs; ++j) {
                        sta.observ[O2][j] = ox_l2kg(sta.observ[(int)OX][j], pr[j], te[j], sa[j]);
                      }
                    }
		    else
		       prop_avail = 0;
		       
                  break;
               case S_:
                  free_and_alloc(&sta.observ[(int)S_], hdr.nobs);
                  compute_sigma(s_pref, hdr.nobs, sta.observ[(int)S_],pr,te,sa);
                  break;

               case S0:
                  free_and_alloc(&sta.observ[(int)S0], hdr.nobs);
                  compute_sigma(0.0, hdr.nobs, sta.observ[(int)S0],pr,te,sa);
                  break;

               case S1:
                  free_and_alloc(&sta.observ[(int)S1], hdr.nobs);
                  compute_sigma(1000.0, hdr.nobs, sta.observ[(int)S1],pr,te,sa);
                  break;

               case S2:
                  free_and_alloc(&sta.observ[(int)S2], hdr.nobs);
                  compute_sigma(2000.0, hdr.nobs, sta.observ[(int)S2],pr,te,sa);
                  break;

               case S3:
                  free_and_alloc(&sta.observ[(int)S3], hdr.nobs);
                  compute_sigma(3000.0, hdr.nobs, sta.observ[(int)S3],pr,te,sa);
                  break;

               case S4:
                  free_and_alloc(&sta.observ[(int)S4], hdr.nobs);
                  compute_sigma(4000.0, hdr.nobs, sta.observ[(int)S4],pr,te,sa);
                  break;

               case TE:
                  free_and_alloc(&sta.observ[(int)TE], hdr.nobs);
                  t90_to_t68(sta.observ[(int)T90],sta.observ[(int)TE], hdr.nobs);
                  break;
		  
               case TH:
                  if (sta.observ[(int)TE] == NULL) {
		     free_and_alloc(&sta.observ[(int)TE], hdr.nobs);
                     t90_to_t68(sta.observ[(int)T90],sta.observ[(int)TE], hdr.nobs);
		  }   
                  free_and_alloc(&sta.observ[(int)TH], hdr.nobs);
                  compute_theta(hdr.nobs, sta.observ[(int)TH], pr, te, sa);
                  break;

               case TH9:
                  if (sta.observ[(int)T90] == NULL) {
		     free_and_alloc(&sta.observ[(int)T90], hdr.nobs);
                     t68_to_t90(sta.observ[(int)T90],sta.observ[(int)TE], hdr.nobs);
		  }   
                     free_and_alloc(&sta.observ[(int)TH9], hdr.nobs);
                     compute_theta(hdr.nobs, sta.observ[(int)TH9], pr, sta.observ[(int)T90], sa);
                  break;

               case HT:
                  free_and_alloc(&sta.observ[(int)HT], hdr.nobs);
                  compute_height(hdr.nobs, pr, te, sa, ht_pref, sta.observ[(int)HT]);
                  break;

               case PE:
                  free_and_alloc(&sta.observ[(int)PE], hdr.nobs);
                  compute_energy(hdr.nobs, pr, te, sa, pe_pref, sta.observ[(int)PE]);
                  break;

               case DR:
	          if (!ratio_done) {
                     free_and_alloc(&sta.observ[(int)DR], hdr.nobs);
                     free_and_alloc(&sta.observ[(int)AL], hdr.nobs);
                     free_and_alloc(&sta.observ[(int)BE], hdr.nobs);
                     compute_ratio(hdr.nobs, sta.observ[(int)DR], pr, te, sa, sta.observ[(int)AL], sta.observ[(int)BE]);
		     ratio_done = 1;
		   }
                  break;

                case SV:
                  free_and_alloc(&sta.observ[(int)SV], hdr.nobs);
                  compute_sp_vol(hdr.nobs, sta.observ[(int)SV], pr, te, sa);
                  break;
              case VA:
                  free_and_alloc(&sta.observ[(int)VA], hdr.nobs);
                  compute_svan(hdr.nobs, sta.observ[(int)VA], pr, te, sa);
                  break;

               case VS:
                  free_and_alloc(&sta.observ[(int)VS], hdr.nobs);
                  compute_sound_vel( sta.observ[(int)VS], pr, te, sa, hdr.nobs);
                break;


               case BF:   
                  free_and_alloc(&sta.observ[(int)BF], hdr.nobs);
                    buoy_freq(sta.observ[(int)BF],pr,te,sa,hdr.nobs,window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) {
                    if (sta.observ[(int)BF][j] > -9998.)
                      sta.observ[(int)BF][j] *= 1.0e5;
                  }
                  break;
               case PV:
                  free_and_alloc(&sta.observ[(int)PV], hdr.nobs);
                  buoy_freq(sta.observ[(int)PV],pr,te,sa,hdr.nobs,window, w_incr);
                  for (j = 0; j < hdr.nobs; ++j) 
                    sta.observ[(int)PV][j] *= sta.observ[(int)PV][j];
                  po_vort(sta.observ[(int)PV],sta.observ[(int)PV], hdr.nobs, (double)hdr.lat);
                  break;

               default:
                    prop_avail = available((enum property) index, &hdr);
                   break;
            } /* end switch */
	 
  } /* end if */    
	    
return(prop_avail);

} /* end get_prop() */
/*****************************************************************************/
void insert_data(double * y, int ny, double *ysig, double *wghtsum, double wght, UI *count, int test_outcrop)

/* For the y property at a station with ny observation levels, interpolate
   to find the yvals at each level in the sigma series,  weight each
   yval, add it to its appropriate sum, sum the weights and increment the counter.  Checks for vertical
   data gaps using the globally defined de array and does not interpolate over
   gaps which exceed GAP_SHALLOW in the thermocline (upper 1000 m) or GAP_DEEP 
   elsewhere.  The sigma series is specified in the global array, siglevs.
   The observed sigma values at this station are subdivided by ref levels, 
   and stored at global addresses pointed to by sigptr. 
   
           y:    array of y observed values
          ny:    dimension of y 
        ysig:    sum of distance-weighted yvalues on each sigma surface
     wghtsum:    sum of weights
       wght:     weight for this profile
       count:    counts # of observations on each surface 
test_outcrop:    0 or 1: if set, tests for outcropping surfaces
*/
{
   int  i, j, datagap, n;
   int  refindx;
   double *xtmp[NREFLEVS], *ytmp, *dtmp, z;
   double reflev;
   
   for (i = 0; i < NREFLEVS; ++i )
      xtmp[i] = (double *) malloc((UI) (ny * sizeof(double)));
   ytmp = (double *) malloc((UI) (ny * sizeof(double)));
   dtmp = (double *) malloc((UI) (ny * sizeof(double)));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in insert_data()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         for (j = 0; j < NREFLEVS; ++j )
            xtmp[j][n] = sigptr[j][i];
         ytmp[n] = y[i];
         dtmp[n++] = de[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     for (i = 0; i < NREFLEVS; ++i)
        free((void *)xtmp[i]);
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nsiglevs; ++i) {

    /* determine which ref lev corresponds to this sigseries element ... */

      refindx = 0;      
      j = nsig[0];
      while (i > j) {
         if (++refindx == NREFLEVS-1)
             break;
         j += nsig[refindx];
      }

      reflev = refindx * 1000.;
      
      if ((test_outcrop) && (siglevs[i] < xtmp[refindx][0])) {   
        /* note:  this assumes x (density) increases with depth!! */           
             ysig[i] += -1 * wght;     /* use depth= -1 for outcropping surfaces */
	     wghtsum[i] += wght;
             ++count[i];
      }

        /* see if associated depth exists at this station ... */

      else if ((z = interpolate(siglevs[i], xtmp[refindx], dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
        /* exclude observations which are more than 1500 m above the reference pressure.
           This prevents bogus values from shelf currents from being added. */
           
           datagap = datagap || (z < reflev-1500);
  
            if (!datagap) {
               if ((z = interpolate(siglevs[i], xtmp[refindx], ytmp, n)) > -998.) {
                  ysig[i] += z * wght;
		  wghtsum[i] += wght;
                  ++count[i];
               }
            }
      }
   }

   for (i = 0; i < NREFLEVS; ++i)
      free((void *)xtmp[i]);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end insert_data() */

/*****************************************************************************/
void insert_seasonal_layer(double * y, int ny, double *ysig, double *wghtsum, double wght, UI *count, int test_outcrop)

/* For the y property at a station with ny observation levels, interpolate
   to find the yvals at each level in the sigma-0 part of the sigseries, weight each
   yval, add it to its appropriate sum, sum the weights and increment the counter.
   Checks for vertical data gaps using the globally defined de array and does not 
   interpolate over gaps which exceed 200 m.
   
   The observed sigma values at this station are subdivided by ref levels, 
   and stored at global addresses pointed to by sigptr. 
   
           y:    array of y observed values
          ny:    dimension of y 
        ysig:    sum of distance-weighted yvalues on each sigma-0 surface
     wghtsum:    sum of weights
       wght:     weight for this profile
       count:    counts # of observations on each surface 
test_outcrop:    0 or 1: if set, tests for outcropping surfaces should only be
                 set when y is depth.
*/
{
   int  i, j, datagap, n;
   double *xtmp, *ytmp, *dtmp, z;
   
   xtmp = (double *) malloc((UI) (ny * sizeof(double)));
   ytmp = (double *) malloc((UI) (ny * sizeof(double)));
   dtmp = (double *) malloc((UI) (ny * sizeof(double)));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in insert_seasonal_layer()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         xtmp[n] = sigptr[0][i];
         ytmp[n] = y[i];
         dtmp[n++] = de[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     free((void *)xtmp);
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nsig[0]; ++i) {

      
      if ((test_outcrop) && (siglevs[i] < xtmp[0])) {   
        /* note:  this assumes x (density) increases with depth!! */           
             ysig[i] += -1 * wght;     /* use depth= -1 for outcropping surfaces */
	     wghtsum[i] += wght;
             ++count[i];
      }

        /* see if associated depth exists at this station ... */

      else if ((z = interpolate(siglevs[i], xtmp, dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
            if (!datagap) {
               if ((z = interpolate(siglevs[i], xtmp, ytmp, n)) > -998.) {
                  ysig[i] += z * wght;
		  wghtsum[i] += wght;
                  ++count[i];
               }
            }
      }
   }

   free((void *)xtmp);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end insert_seasonal_layer() */

/*****************************************************************************/
void insert_deepest_vals(struct gridnode *g_ptr, int *prop_indx, int nprops)
   /* Maintains a linked list of info regarding the deepest observed 
      level. 
      All observations within the deepest 100 meters will be retained.
      The list will actively grow and contract as each new station is
      read. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
       
         g_ptr:   ptr to current gridnode  
      prop_indx:   contains index to properties being gridded 
         nprops:   number of properties being gridded 
   */
{
   int i, j, index,reflev, itop, ibotm, n;
   int density_greater, depth_greater, depth_diff;
   double key, botm_density, tbar, sbar, ref_pr;
   double *a_ptr;
   struct deepestrec  *dr_ptr, *r1, *r2;
 
   dr_ptr = g_ptr->deepest; 
   j = hdr.nobs - 1;
   if (sta.observ[(int)PR][j] < 0)
       return;
   
   /* determine appropriate sigma reference level */
   reflev = 0;
   if (sta.observ[(int)PR][j] > zmin[4])
          reflev = 4;
   else if (sta.observ[(int)PR][j] > zmin[3])
          reflev = 3;
   else if (sta.observ[(int)PR][j] > zmin[2])
          reflev = 2;
   else if (sta.observ[(int)PR][j] > zmin[1])
          reflev = 1;
   
   ref_pr = reflev *1000.0;
   
   /* Check linked list, if it exists: 
      if it is not denser or close in depth to deepest rec so far, skip it */
      
   if (dr_ptr != NULL) {
       depth_greater = sta.observ[(int)DE][j] > dr_ptr->depth ? 1 : 0;
       depth_diff =  sta.observ[(int)DE][j] - dr_ptr->depth;
       switch (reflev){
            case 0:
              density_greater = sta.observ[(int)S0][j] > dr_ptr->sig_0 ? 1 : 0;
              break;
            case 1:
              density_greater = sta.observ[(int)S1][j] > dr_ptr->sig_1 ? 1 : 0;
	      break;
            case 2:
               density_greater = sta.observ[(int)S2][j] > dr_ptr->sig_2 ? 1 : 0;
               break;
            case 3:
               density_greater = sta.observ[(int)S3][j] > dr_ptr->sig_3 ? 1 : 0;
               break;
            case 4:
               density_greater = sta.observ[(int)S4][j] > dr_ptr->sig_4 ? 1 : 0;
               break;
            default:
               fprintf(stderr,"FATAL ERROR in insert_deepest_vals():  reflev is %d", reflev);
	       exit(1); 
      } /* end switch */

      if (!density_greater && (depth_diff < -200.))	 
         return;
      	 
   } /* end if dr_ptr */
   
   /* Add a new record */
        
   dr_ptr = create_deepestrec(nprops);
      
   dr_ptr->depth = sta.observ[(int)DE][j];
   dr_ptr->pressure = sta.observ[(int)PR][j];
      
   a_ptr = sta.observ[(int)S0];
   if (reflev == 1)
         a_ptr = sta.observ[(int)S1];
   else if (reflev == 2)
         a_ptr = sta.observ[(int)S2];
   else if (reflev == 3)
         a_ptr = sta.observ[(int)S3];
   else if (reflev == 4)
         a_ptr = sta.observ[(int)S4];
	 
   /* find thickness of bottom boundary layer */	 
   ibotm = j;
   botm_density = a_ptr[j];
   key = botm_density - bottom_layer_def;
   --j;
   while ((j > 0) && (a_ptr[j] > key)) {
         if (a_ptr[j] > botm_density)  {   /* check for denser vals above bottom */
	    botm_density = a_ptr[j];
	    key = botm_density - bottom_layer_def;
	 }
	 --j;   
   }
   itop = j+1;
   dr_ptr->thickness = sta.observ[(int)DE][ibotm] - sta.observ[(int)DE][itop]; 
      
   /* find average t and s in layer, compute densities */   
   tbar = sbar = 0.0;
   n = 0;
   for (j = ibotm; j >= itop; --j) {
        tbar += sta.observ[(int)TH9][j];
	sbar += sta.observ[(int)SA][j];
	++n;
   }
   tbar /= n;
   sbar /= n;
   
   /* first revert to in situ temperature */
   tbar = hb_theta(sbar, tbar, 0, dr_ptr->pressure);
   
   compute_sigma(0, 1, &dr_ptr->sig_0, &dr_ptr->pressure, &tbar, &sbar);
   compute_sigma(1000, 1, &dr_ptr->sig_1, &dr_ptr->pressure, &tbar, &sbar);
   compute_sigma(2000, 1, &dr_ptr->sig_2, &dr_ptr->pressure, &tbar, &sbar);
   compute_sigma(3000, 1, &dr_ptr->sig_3, &dr_ptr->pressure, &tbar, &sbar);
   compute_sigma(4000, 1, &dr_ptr->sig_4, &dr_ptr->pressure, &tbar, &sbar);
      
    /* find average of each property over bottom layer */
      
   for (i=0; i< nprops; ++i) {
     dr_ptr->prop[i] = HBEMPTY;
     
     index = prop_indx[i];
     if (index == (int)TE || index == (int)T90 || index == (int)TH)        /* work with theta not in situ TE */
         index = (int) TH9;
	 
     if (sta.observ[prop_indx[i]] != NULL) {
        if (prop_indx[i] == (int)PR) {
	    dr_ptr->prop[i] = dr_ptr->pressure;
	    ++dr_ptr->count[i];
	}
	else {
           dr_ptr->prop[i] = 0.0;
           for (j = ibotm; j >= itop; --j) {
             if (sta.observ[index][j] > -8.9) {
               dr_ptr->prop[i] += sta.observ[index][j];
	       ++dr_ptr->count[i];
	     }
           }
	   if (dr_ptr->count[i] > 0 )
	      dr_ptr->prop[i] /= dr_ptr->count[i]; 
	   else
	      dr_ptr->prop[i] = HBEMPTY;
	}
     }
   } /* end for i */



/* Insert new record into the linked list */   
  
   if (g_ptr->deepest  == NULL) {    /* empty list */
       g_ptr->deepest = dr_ptr;
       dr_ptr->next = NULL;
       return;
   } 
   
     
/* If station is denser than previous densest:
    insert new rec at beginning of list
    check remainder of list and delete any recs more than 200 m shallower */
     
   if (density_greater ) {
      
      if ( !depth_greater) {   /* set depth to deeper of the two */
         dr_ptr->depth = g_ptr->deepest->depth;
         dr_ptr->pressure = g_ptr->deepest->pressure;
	 dr_ptr->prop[ipr] = g_ptr->deepest->pressure;
      }
      
      dr_ptr->next = g_ptr->deepest;
      
      key = dr_ptr->depth - 200.;
      r1 = dr_ptr;
      r2 = dr_ptr->next;
      while (r2 != NULL) {
         if (r2->depth < key) {
            delete_list(r2->next);
            free((void *) r2->prop);
	    free((void *) r2->var);
	    free((void *) r2->count);
            free((void *)r2);
            r2 = NULL;
            r1->next = NULL;
         }
         else {
            r1 = r2;
            r2 = r2->next;
         }
      }
      g_ptr->deepest = dr_ptr;
      return;
   }
   
   
   /* It's not denser, but if it is  deeper than the densest observation so far, 
       replace the deepest depth/pressure and then go on to insert a new record 
       into the linked list */
   
   if (depth_greater) {
     j = hdr.nobs - 1;
     r1 = g_ptr->deepest;
     r1->depth = sta.observ[(int)DE][j];
     r1->pressure = sta.observ[(int)PR][j];
     r1->prop[ipr] = r1->pressure;
   }
   
 /*  it's not the densest, if it is shallower than 200 m above the
     densest, we're done.....*/

   key = g_ptr->deepest->depth -200.;

   if (dr_ptr->depth  < key)  
      return;
      
      
  /*  It is within 200 m of densest observation, insert it into linked list*/
    
   r1 = g_ptr->deepest;
   r2 = r1->next;
   do {
   
     if (r2 == NULL) {   /* insert at end of list */
       r1->next = dr_ptr;
       dr_ptr->next = NULL;
       return;
     }
     
     if (dr_ptr->sig_2 > r2->sig_2) { /*insert before r2 */
       r1->next = dr_ptr;
       dr_ptr->next = r2;
       
       /* and check remainder of list against key */
       
       while (r2 != NULL) {
         if (r2->depth < key) {
            delete_list(r2->next);
            free((void *) r2->prop);
            free((void *) r2->var);
            free((void *) r2->count);
            free((void *)r2);
            r2 = NULL;
            r1->next = NULL;
         }
         else {
            r1 = r2;
            r2 = r2->next;
         }
       }
       return;
     }
     
     r1 = r2;
     r2 = r2->next;
   
   } while (1);    
} /* end insert_deepest_vals () */
/*****************************************************************************/
struct deepestrec * create_deepestrec(int nprops)
{
   struct deepestrec *r1;
   
   r1 = (struct deepestrec *) malloc(sizeof(struct deepestrec));
   if (r1 == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
   r1->count = (int *) calloc((size_t)nprops, sizeof(int));
   r1->prop = (double *) calloc((size_t)nprops, sizeof(double));
   r1->var = (double *) calloc((size_t)nprops, sizeof(double));
   if (r1->var == NULL) {
       fprintf(stderr,"\nOut of memory in create_deepestrec() \n");
       exit(1);
   }
  
   return(r1);
   
} /* end create_deepestrec() */
/*****************************************************************************/
void delete_list(struct deepestrec *rptr)

  /* Recursively traverses a linked list of records and frees
     up all the memory allocated to those records */
{
  /* end of list... */
   if (rptr == NULL)
      return;
      
  /* recursive part... */
   delete_list(rptr->next);
   free((void *) rptr->prop);
   free((void *) rptr->var);
   free((void *) rptr->count);
   free((void *) rptr);
   return;

}  /* end delete_list() */
/*****************************************************************************/
/*****************************************************************************/
void delete_surfrec_list(struct surfrec *rptr)

  /* Recursively traverses a linked list of records and frees
     up all the memory allocated to those records */
{
  /* end of list... */
   if (rptr == NULL)
      return;
      
  /* recursive part... */
   delete_surfrec_list(rptr->next);
   free((void *) rptr->prop);
    free((void *) rptr->var);
  free((void *) rptr->wghtsum);
   free((void *) rptr);
   return;

}  /* end delete_surfrec_list() */
/*****************************************************************************/
/*****************************************************************************/
void insert_surf_vals(struct gridnode *g, int *prop_indx, int nprops, int month)
   /* defines the depth of a mixed layer at the sea surface, determines 
      the average values of properties within it, and inserts the information
      into a linked list of records sorted by density. The station data are 
      accessed through the global variable: struct HYDRO_DATA sta 
              g:   ptr to gridnode  
      prop_indx:   contains index to properties being gridded 
         nprops:   number of properties being gridded 
	 month:    index to start of linked list
   */
{
   int i, j, n, weight, index;
   double x, v, vprev, dprev;
   double dens, depth;
   struct surfrec *m;

   if (de[0] > 10.0)
       return;

   /* round off surface density to nearest tenth */

   dens =  (double) (NINT(sta.observ[(int)S0][0] * 10.)) / 10.0;

  /* bottom of mixed layer is defined where 
    density = surface density + mix_layer_def */  
     
   depth = interpolate((sta.observ[(int)S0][0]+ mix_layer_def), sta.observ[(int)S0], de, hdr.nobs);

   if (depth < -999.)  /* all observations are within defined range of surface density */
       depth = de[hdr.nobs-1];

   m = get_surfrec(g->mix_layer[month], dens, nprops);

   if (g->mix_layer[month] == NULL ) {   /* empty list */           
      g->mix_layer[month] = m;
   }

   /* add depth info to this record */

   m->depth += depth;
   ++ m->n;

  /* compute average properties for mixed layer and add to appropriate fields.
     The average value is computed by summing the observed values weighted by
     the depth between observations.*/

   for (i = 0; i < nprops; ++i) {
   
      index = prop_indx[i];
      if ( index == (int) TE || index == (int)T90 || index ==(int)TH)
         index = (int) TH9;
	 
      if (sta.observ[index] != NULL) {
         j = 0;
         x = sta.observ[index][0];
         n = 1;               /* weight the observation at sea surface */
         if (x < -8.9 ) {     /* unless it is a missing value */
            x = 0.0;
            n = 0;
         }
         dprev = 0.0;
         vprev = sta.observ[index][0];
         while ( (j < hdr.nobs) && (de[j] <= depth)) {
            if ( (v = sta.observ[index][j]) > -8.9) {
               if (vprev < -8.9) 
                   vprev = v;
               weight = (int) (de[j] - dprev);
               x += (v + vprev) * .5 * weight;
               n += weight;
               dprev = de[j];
               vprev = v;
            }
            ++j;
         } 
         if (n > 0) {
            m->prop[i] += (x / (float) n);
	    m->wghtsum[i] = m->wghtsum[i] + 1.0;   /* this counts number of profiles in this mixed layer */
         }
      }
   }
   return;

}  /* end insert_surf_vals() */

/*****************************************************************************/
struct surfrec *get_surfrec(struct surfrec *rptr, double d, int nprops)

   /* Recursively searches a linked list of surfrecs to:
        a) find an existing record for the specified density;
     or b) create a new record and insert it into the ordered list.

      Returns a pointer to the appropriate record.
      
  arguments: 
       rptr:    pointer to start of list 
          d:    density to key on 
     nprops:    # of properties to allocate space for
   */
{
    struct surfrec *r1ptr;
    double *tempd, *tempd2, *tempd3;

    if (rptr == NULL) {         /* empty list */
       r1ptr = create_surfrec(nprops);
       r1ptr->density = d;
       return(r1ptr);
    }

    if (NINT(d * 10) == NINT(rptr->density * 10)) {  /* current rec */
        return (rptr);
    }

    if (d < (rptr->density - .00001)) {   /* insert before the current rec */

       r1ptr = create_surfrec(nprops);
       tempd = r1ptr->prop;
       tempd2 = r1ptr->wghtsum;
       tempd3 = r1ptr->var;

         /* copy all fields from rptr into r1ptr */
       r1ptr->density = rptr->density;
       r1ptr->depth = rptr->depth;
       r1ptr->prop = rptr->prop;
       r1ptr->var = rptr->var;
       r1ptr->wghtsum = rptr->wghtsum;
       r1ptr->n = rptr->n;
       r1ptr->next = rptr->next;

        /* initialize the fields of rptr and link it to r1ptr */
       rptr->density = d;
       rptr->depth = 0;
       rptr->prop = tempd;
       rptr->var = tempd3;
       rptr->wghtsum = tempd2;
       rptr->n = 0;
       rptr->next = r1ptr;
       
       return(rptr);
    }

    r1ptr = get_surfrec(rptr->next, d, nprops);  /* search rest of list */
    if (rptr->next == NULL)
          rptr->next = r1ptr;
    return (r1ptr);

}   /* end get_surfrec() */

/*****************************************************************************/
struct surfrec *create_surfrec(int nprops)

   /* Allocates memory and initializes the fields of a struct surfrec.
      Returns a pointer to the record */
{
   struct surfrec *r;
   int i;

   r = (struct surfrec *) malloc(sizeof(struct surfrec));
   if (r == NULL) {
      fprintf(stderr,"\nUnable to allocate memory in create_surfrec()\n");
      exit(1);
   }
   r->depth = 0;
   r->density = 0;
   r->n = 0;
   r->next = NULL;
   r->wghtsum = (double *) malloc(nprops * sizeof(double));
   r->prop = (double *) malloc(nprops * sizeof(double));
   r->var = (double *) malloc(nprops * sizeof(double));

   if (r->prop == NULL) {
      fprintf(stderr,"\nUnable to allocate memory in create_surfrec()\n");
      exit(1);
   }

   for (i = 0; i < nprops; ++i) {
      r->prop[i] = 0.0;
      r->var[i] = 0.0;
      r->wghtsum[i] = 0.0;
   }   

   return (r);
}   /* end create_surfrec() */
/*****************************************************************************/
struct surfrec *copy_surfrec(struct surfrec *r2cpy, int nprops)

/* copies the record into a new record and returns a pointer to it  */
{
    struct surfrec *newrecptr;
    int i;
    
   if (r2cpy == NULL)
       return (NULL);
       
   newrecptr = create_surfrec(nprops);
   newrecptr->depth = r2cpy->depth;
   newrecptr->density = r2cpy->density;
   newrecptr->n = r2cpy->n;
   newrecptr->next = NULL;   
   for (i = 0; i < nprops; ++i) {
      newrecptr->prop[i] = r2cpy->prop[i];
      newrecptr->var[i] = r2cpy->var[i];
      newrecptr->wghtsum[i] = r2cpy->wghtsum[i];
   }
       
   return(newrecptr);
}  /* end copy_surfrec() */

/*****************************************************************************/
void compute_avg(double *sum, double *weightsum, int nlevs)

   /*   sum:    array containing summed, weighted values 
  weightsum:    array containing sum of weights for each level
      nlevs:    # of elements in array 
   */
{
   int j;

   for (j = 0; j < nlevs; ++j) {
      sum[j] = ( weightsum[j] > 0) ?  sum[j] / weightsum[j] : (double) HBEMPTY;
   }
   return;

}  /* end compute_avg() */
/*****************************************************************************/
void do_variance_prep(double * y, int ny, double *ybar,  double wght, double *dprop, double *dpropsq, int test_outcrop)

/* For the profile of property values with ny observation levels: 
      interpolate to find the yvals at each sigma level in the sigma series,  
      compute difference (yval - ybar) at that level, 
      add it to its appropriate sum.  No need to sum weights or count obs as that was already done
      for computing the means. 
      
   Checks for vertical data gaps using the globally defined de array and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere.  The sigma series is specified in the global array, siglevs.
   The observed sigma values at this station are subdivided by ref levels, 
   and stored at global addresses pointed to by sigptr. 
   
           y:    array of y observed values
          ny:    dimension of y 
        ybar:    average yvalue on each sigma surface (previously computed).
       wght:     weight for this profile
       dprop:    weighted sum of property deviations from ybar
       dpropsq:  weighted sum of squared property deviations
test_outcrop:    0 or 1: if set, tests for outcropping surfaces
*/
{
   int  i, j, datagap, n;
   int  refindx;
   double *xtmp[NREFLEVS], *ytmp, *dtmp, z;
   double delta, reflev;
   
   for (i = 0; i < NREFLEVS; ++i )
      xtmp[i] = (double *) calloc(ny , sizeof(double));
      
   ytmp = (double *) calloc(ny, sizeof(double));
   dtmp = (double *) calloc(ny, sizeof(double));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in do_variance_prep()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         for (j = 0; j < NREFLEVS; ++j )
            xtmp[j][n] = sigptr[j][i];
         ytmp[n] = y[i];
         dtmp[n++] = de[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     for (i = 0; i < NREFLEVS; ++i)
        free((void *)xtmp[i]);
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nsiglevs; ++i) {

    /* determine which ref lev corresponds to this sigseries element ... */

      refindx = 0;      
      j = nsig[0];
      while (i > j) {
         if (++refindx == NREFLEVS-1)
             break;
         j += nsig[refindx];
      }

      reflev = refindx * 1000.;
      
      if ((test_outcrop) && (siglevs[i] < xtmp[refindx][0])) {   
        /* note:   this adds a zero to variance (may have to think about this) */           
        /*  although depth= -1 for outcropping surfaces, weightsum and count arrays are incremented*/
      }

        /* see if associated depth exists at this station ... */

      else if ((z = interpolate(siglevs[i], xtmp[refindx], dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
        /* exclude observations which are more than 1500 m above the reference pressure.
           This prevents alot of bogus values from being added. */
           
           datagap = datagap || (z < reflev-1500);
  
            if (!datagap) {
               if ((z = interpolate(siglevs[i], xtmp[refindx], ytmp, n)) > -998.) {
                  delta =  ABS(z - ybar[i]);
		  dprop[i] += delta * wght;
		  dpropsq[i] += delta * delta * wght;
               }
            }
      }
   }

   for (i = 0; i < NREFLEVS; ++i)
      free((void *)xtmp[i]);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end do_variance_prep() */

/*****************************************************************************/
void do_variance_prep_monthly(double * y, int ny, double *ybar,  double wght, double *dprop, double *dpropsq, int test_outcrop)

/* For the profile of property values with ny observation levels: 
      interpolate to find the yvals at each sigma level in the sigma series,  
      compute difference (yval - ybar) at that level, 
      add it to its appropriate sum.  No need to sum weights or count obs as that was already done
      for computing the means. 
      
   Checks for vertical data gaps using the globally defined de array and does not interpolate over
   gaps .  The sigma series is specified in the global array, siglevs --
   only the sigma-0 part is used.
   The observed sigma values at this station are subdivided by ref levels, 
   and stored at global addresses pointed to by sigptr. 
   
           y:    array of y observed values
          ny:    dimension of y 
        ybar:    average yvalue on each sigma surface (previously computed).
       wght:     weight for this profile
       dprop:    weighted sum of property deviations from ybar
       dpropsq:  weighted sum of squared property deviations
test_outcrop:    0 or 1: if set, tests for outcropping surfaces
*/
{
   int  i, j, datagap, n;
   double *xtmp, *ytmp, *dtmp, z;
   double delta;
   
   xtmp = (double *) calloc(ny, sizeof(double));
   ytmp = (double *) calloc(ny, sizeof(double));
   dtmp = (double *) calloc(ny, sizeof(double));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in do_variance_prep()\n");
       exit(1);
   }

/* Ensure continuous (no missing values) x,  y, and depth arrays */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         xtmp[n] = sigptr[0][i];
         ytmp[n] = y[i];
         dtmp[n++] = de[i];
      }
   }
   
   if (n <= 1) {               /* not enough values */
     free((void *)xtmp);
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   for (i = 0; i < nsig[0]; ++i) {

      
      if ((test_outcrop) && (siglevs[i] < xtmp[0])) {   
        /* note:   this adds a zero to variance (may have to think about this) */           
        /*  although depth= -1 for outcropping surfaces, weightsum and count arrays are incremented*/
      }

        /* see if associated depth exists at this station ... */

      else if ((z = interpolate(siglevs[i], xtmp, dtmp, n)) > -998.) {

            /* now check for datagaps around this depth...*/
            j = 0;
            while ((++j < n) &&( dtmp[j] < z)  ) 
                ;
	    if (j == n)
	        --j;
		
            if ((dtmp[j-1] == z) || (dtmp[j] == z) )
                datagap = 0;
            else if (z < 1001)
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
            else
                datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
                
  
            if (!datagap) {
               if ((z = interpolate(siglevs[i], xtmp, ytmp, n)) > -998.) {
                  delta =  ABS(z - ybar[i]);
		  dprop[i] += delta * wght;
		  dpropsq[i] += delta * delta * wght;
               }
            }
      }
   }

   free((void *)xtmp);
   free((void *)ytmp);
   free((void *)dtmp);
   return;
} /* end do_variance_prep_monthly() */
/*****************************************************************************/

void compute_variance(double *dprop, double *dpropsq, double *weightsum, int *count, int nlevs)

   /*   dprop:    array containing sum of (xi - xbar) * weight 
      dpropsq:    contains sum of (xi-xbar)^^2 * weight
  weightsum:    array containing sum of weights for each level
      count:    array containing nobs at each level
      nlevs:    # of elements in array 
      
      Returns variance in dpropsq array.  If < 4 observations, variance set to empty.
   */
   
{
   double av_dev_sq;
   int j;

   for (j = 0; j < nlevs; ++j) {
   
      if (count[j] < 4) {                   /* too few observations */
           dpropsq[j] = (double) HBEMPTY;
	   dprop[j] = (double) HBEMPTY;
      }
      else {
         dprop[j] /=  weightsum[j];  /* average deviation */
         av_dev_sq = dprop[j] * dprop[j];   /*  av_deviation squared =*/
         dpropsq[j] =  (dpropsq[j] / weightsum[j] - av_dev_sq) * count[j]/ (count[j]-1) ;
      }
   }
   return;

}  /* end compute_variance() */
/*****************************************************************************/
void define_sigma_transitions(struct gridnode *g, int nprops)

 /*  determines where the averaged depth values cross the depth ranges (zmin
    and zmax) for each sigma reference level and zeros the counters for the 
    accumulated sums of any level outside that range. 
          g:    ptr to a grid node 
     nprops:    number of properties stored at each gridnode 
  */
{
   int i, j, n, refstart, refstop, iref;
   
  /* trim seasonal layer to upper 200 m */
   
   /* For each month... */
   for (n = 0; n < NMONTHS; ++n) {

     /* For each reference level... */
     for (i = 0; i < nsig[0]; ++i) {

       /* If there are observations this month at this location... */
       if (g->month_dnobs[n][i] > 0) {
	   if (g->month_d[n][i] < 0 || g->month_d[n][i] >= bottom_seasonal_layer) {
	     g->month_dnobs[n][i] = 0;
	     for (j = 0; j < nprops;  ++j) {
	       g->month_count[n][j][i] = 0;
	     }
	   }
	 }
      }
   }
   /* trim sig0 levels above bottom_seasonal_layer */  
   refstart = 0;
   refstop = nsig[0];
   iref = 0;
   for (i = refstart; i < refstop; ++i) {
         if (g->nobs[i] > 0) {
            if (g->d[i] < bottom_seasonal_layer || g->d[i] >= zmax[iref]) {
                g->nobs[i] = 0;
                for (j = 0; j < nprops; ++j) {
                   g->count[j][i] = 0;
                }
            }
         }
   }

    /* trim remaining sig levels */  
        
    for (iref = 1; iref < NREFLEVS; ++iref) {
      refstart = refstop;
      refstop = refstart + nsig[iref];

      for (i = refstart; i < refstop; ++i) {
         if (g->nobs[i] > 0) {
            if (g->d[i] < zmin[iref] || g->d[i] >= zmax[iref]) {
                g->nobs[i] = 0;
                for (j = 0; j < nprops; ++j) {
                   g->count[j][i] = 0;
                }
            }
         }
      }

   }
   
   return;

}  /* end define_sigma_transitions() */
/*****************************************************************************/
struct surfrec *define_avg_mixed_layer(struct surfrec *listptr, int nprops)

   /* Traverses a linked list of surfrecs sorted by density and computes 
      the average density at the sea surface, then searches the linked list
      for that density or 2 density records bracketing the average density and
      returns a pointer to a surfrec containing the property values
      associated with the average density. The variance and number of obs
      are computed for each property and returned in the struct surfrec.
      
      listptr:   address of first record in linked list
       nprops:   number of properties being averaged 
   */ 
{
   struct surfrec *r1;
   double x;
   UI n;
   int key;
   
   if (listptr == NULL)                         /* empty list */
       return ((struct surfrec *) NULL);

   r1 = listptr;
   x = 0.0;
   n = 0;
   while (r1 != NULL) {
      x += r1->density ; 
      ++n;   
      r1 = r1->next;    
   }
   key = NINT((x / (double)n) * 10.);     /* average density * 10  */
   
   /* Now search list for the density value and compute statistics... */

   r1 = get_density_rec(key, listptr, nprops);
   compute_ml_stats(listptr, r1, nprops);
   
   return (r1);

}  /* end define_avg_mixed_layer() */
/***********************************************************************/
void compute_ml_stats(struct surfrec *listptr, struct surfrec *mlptr, int nprops)
   /* Traverse linked list of surfrecs, compute means, of each mixed layer, then variance and nobs for all
      observations, store variances in mlptr->var array and nobs in mlptr->wghtsum array. 

      listptr :  start of linked list
        mlptr:   record corresponding to the average mixed layer for this node
nprops:  # of properties stored
      
      */

{
   int i, *count;
   struct surfrec *r1;
   double *dprop, *dpropsq, delta;

   if (listptr == NULL || mlptr == NULL)        /* this should never happen */
        return;

   dprop = (double *) calloc(nprops, sizeof(double));	
   dpropsq = (double *) calloc(nprops, sizeof(double));	
   count = (int *) calloc(nprops, sizeof(int));	
   
   /* traverse list and compute averages in each mixed layer record */

   r1 = listptr;
   
   while (r1 != NULL) {
       compute_avg(r1->prop, r1->wghtsum, nprops);
       r1->depth /= (double)r1->n;
       r1 = r1->next;
   }
   /* traverse list again and sum the deviations from mlptr->prop */
   
   r1 = listptr;
   while (r1 != NULL) {
      for ( i = 0; i < nprops; ++i) {
          delta = ABS(r1->prop[i] - mlptr->prop[i]);
          dprop[i] += delta * r1->wghtsum[i];
	  dpropsq[i] += delta * delta * r1->wghtsum[i];
	  count[i] += r1->wghtsum[i];
      }
      r1 = r1->next;
   }
   
   /* now compute variance for each property */
   for ( i = 0; i < nprops; ++i) {
      mlptr->var[i] = (dpropsq[i] - dprop[i] * dprop[i]/count[i]) / (count[i] - 1);
      mlptr->wghtsum[i] = count[i];
   }
    
   free((void *)dprop);
   free((void *)dpropsq); 
   free((void *)count);
   return;  
}  /* end compute_ml_stats() */
/***********************************************************************/
struct surfrec *get_density_rec(int key, struct surfrec *rptr, int nprops)
    /* Recursively searches a linked list sorted by density for the record
       corresponding to key or for two records bracketing that density.
       Returns a pointer to that record. 
       
          key:  density  being searched for *10  
         rptr:  ptr to element of linked list 
       nprops:  # of properties stored in each surfrec 
 
    */
{
   struct surfrec *r1, *r2, *new;
   double x[2], y[2];
   UI n;
   int  key1, key2, i;

   if (rptr == NULL) {           /* YIKES! */
       fprintf(stderr,"\nError in get_density_rec(): ");
       fprintf(stderr,"\nEnd of linked list encountered before finding avg density!");
       fprintf(stderr,"\nError in program logic....exiting.\n");
       exit(1);
   }

   r1 = rptr;
   r2 = rptr->next;

   if (r2 == NULL)             /* end of list, so r1 must be the density */
      return (r1);

   key2 = NINT(r2->density * 10);

   if (key > key2)                /* recursive part */
       return( get_density_rec(key, r2, nprops));
   

   if (key == key2)             /* exact match! */
       return(r2);
   
     
   key1 = NINT(r1->density * 10);

   if (key < key1) {              /* YIKES! */
       fprintf(stderr,"\nError in get_density_rec(): ");
       fprintf(stderr,"\nAvg density is less than first density in list!");
       fprintf(stderr,"\nThis is a bug in the program logic....exiting.\n");
       exit(1);
   }
   if (key == key1)             /* exact match! */
       return(r1);
   

/* if we get this far, then the key must lie between r1 & r2... so 
      choose r1 or r2-- whichever is closest.   */


   if ((key - key1) < (key2 - key) )
       return(r1);
   
   else 
       return(r2);
   

}  /* end get_density_rec() */

/************************************************************/
void write_to_file(int nc_file,struct gridnode ***grid, struct CDF_HDR *hptr,int *prop_indx, int nstdlevs, int ntbins, int nrows, int ncols, int imonth)
   /* handles all the functions for writing data to the nc_file.
      On entry to this function, all temperature values are TH9 values. 
      TH9 is converted to requested properties in this function. */
{
   int error, i, nprops, nlevs, npts;
   int row, col, tbin;
   float  **data, **vari, *bottomdepth;  
   short  **count; 
   float *p, *t, *s;

   nprops = hptr->nprops;

   error = cdf_define(nc_file, hptr, PREFILL, ADD_PROP_STATS );
   error = write_std_depths_cdf(nc_file, hptr);
   error = write_time_bins_cdf(nc_file, hptr);
   error = write_lat_vector(nc_file, hptr);
   error = write_lon_vector(nc_file, hptr);


   data = (float **) calloc(nprops, sizeof(float *));
   vari = (float **) calloc(nprops, sizeof(float *));
   count = (short **) calloc(nprops, sizeof(short *));
   bottomdepth = (float *) calloc(ncols, sizeof(float));
   
  npts = nstdlevs * ncols;

   for (i = 0; i < nprops; ++i) {
      data[i] = (float *) calloc(npts, sizeof(float));
      vari[i] = (float *) calloc(npts, sizeof(float));
      count[i] = (short *) calloc(npts, sizeof(short));
   }

   if (data[nprops-1] == NULL) {
      fprintf(stderr, "\nUnable to allocate memory for data & count arrays.\n");
      exit(1);
   }

/* interpolate the sigma series back onto the standard depth levels and 
   output the property data to the netcdf file... */

   fprintf(stderr,"  writing data for %s ...\n", print_month(imonth));

   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
            get_lat_lon(hptr, row, col, &latitude, &longitude);
            nlevs = do_std_depth(&grid[tbin][row][col],imonth, nsiglevs, prop_indx, 
                  nprops, ncols, col, data, vari, count, &bottomdepth[col]);
         }
	 
	 /* At this point, all T variables are TH9 values. 
	    If necessary, convert to requested output temperature  */
	 
	 s = data[isa];
	 p = data[ipr];
	 t = data[ite];
	 
	 switch (tindex) {
	    case (int)TE:
	    case (int)TH:
	       for (i = 0; i < npts; ++i) {
	        if (t[i] > -3 && t[i] < 50) 
	          data[ite][i] = (float) hb_theta((double)s[i], (double)t[i], 0, (double)p[i]);
		  data[ite][i] *= 1.00024;  /* T68 */
		  
		  if (tindex == (int)TH)
		    data[ite][i] = (float) hb_theta((double)s[i], (double)t[i], (double)p[i], 0);
	       }
	       break;
	       
	    case (int)T90:
	       for (i = 0; i < npts; ++i) {
	        if (t[i] > -3 && t[i] < 50) 
	          data[ite][i] = (float) hb_theta((double)s[i], (double)t[i], 0, (double)p[i]);
	       }
	       break;
	       
	 } /* end switch */
	 	     
	 

         for (i = 0; i < nprops; ++i) {
            write_prop_cdf(nc_file, data[i], hptr->prop_id[i], row, 0, tbin, 0,
                           1, ncols, 1, nlevs);
            write_prop_count_cdf(nc_file, count[i], hptr->prop_id[i],
                           row, 0, tbin, 0, 1, ncols, 1, nlevs);
            write_prop_err_cdf(nc_file, vari[i], hptr->prop_id[i], row, 0, tbin, 0,
                           1, ncols, 1, nlevs);

         }
         write_bottom_depth_cdf(nc_file, row, 0, tbin, 1, ncols, 1,
                                bottomdepth);
      }
   }
   /* writing these a second time circumvented a bug in netcdf lib */
   error = write_time_bins_cdf(nc_file, hptr);
   error = write_std_depths_cdf(nc_file, hptr);
   error = write_lat_vector(nc_file, hptr);
   error = write_lon_vector(nc_file, hptr);
   
   free((void *) data);
   free((void *) count);
   free((void *) bottomdepth);
   free((void *) vari);
   return;
} /* end write_to_file() */

/****************************************************************************/
/***********************************************************************/
int do_std_depth(struct gridnode *g, int imonth, int npts, int *prop_id, int nprops, int ncols, int icol, float **dataptr, float **varptr, short **countptr, float *bottom)

/*  Takes the info for a gridnode and computes the value of each property 
    interpolated onto the std depth levels.  The arrays at dataptr, varptr and countptr
    are assumed to have dimensions data[nprops][nstdlevs*ncols].  This
    function deals with one gridpoint at a time, therefore icol must be
    specified to determine the starting position within the array for
    each property; i.e. starting address = data[iprop][nstdlevs * icol].  
    (The arrays have been setup to include an entire row of
    data to improve the efficiency of writing out to the netcdf file).
    THE CALLING PROGRAM IS RESPONSIBLE FOR ALLOCATING THE SPACE AT
     *dataptr, *varptr AND *countptr: 
   
    The # of observations at each standard level and for each property 
    are approximated and returned starting at the address pointed to by 
    countptr.

    The function returns the # of std levels (including the bottom) . 
    std_depth[MAXSTDLEVS] is globally defined and already initialized.

            g:   address of this gridnode  
       imonth:   month index for seasonal and mixed layers
         npts:   size of arrays for this gridnode  
      prop_id:   array of property identifiers  
       nprops:   # of properties in prop array  
        ncols:   # of cols in data and count arrays  
         icol:   column # associated with this gridpt  
      dataptr:   ptr to start of data arrays   
      varptr:   ptr to start of variance arrays; stddev = sqrt(var) is returned   
     countptr:   ptr to start of count arrays  
       bottom:   returned depth of deepest observation 

*/
{
   int i, ii, j, jj, k, datagap, room, is_denser;
   double *x, *y, *v, *w, *sigval, *sig, z, r, zv;
   double *d_tmp, **prop_tmp, **var_tmp;
   double tref, pref, sig_bml, sig_deepest, top_of_layer;
   float  *d, *var;
   short  *n;
   int size, start, nb;
   int deepestlev, reflev;
   UI *nobs_tmp, **count_tmp;
   struct deepestrec *btmptr, *b1;
   struct surfrec *m;
   extern double std_depth[];

   sigval = (double *) malloc(npts * sizeof(double));
   sig = (double *) malloc(npts * sizeof(double));
   d_tmp = (double *) calloc(npts, sizeof(double));
   nobs_tmp = (UI *) calloc(npts, sizeof(UI));
   prop_tmp = (double **) calloc(nprops, sizeof(double *));
   var_tmp = (double **) calloc(nprops, sizeof(double *));
   count_tmp = (UI **) calloc(nprops, sizeof(UI *));
   for (i = 0; i < nprops; ++i) {
      prop_tmp[i] = (double *) calloc(npts, sizeof(double));
      var_tmp[i] = (double *) calloc(npts, sizeof(double));
      count_tmp[i] = (UI *) calloc(npts, sizeof(UI));
     
   }      
   size = monotonic(g->d, g->nobs, g->prop, g->dpropsq, g->count, nprops, npts, g->month_d[imonth], g->month_dnobs[imonth], g->month_prop[imonth], g->month_dpropsq[imonth], g->month_count[imonth], nsig[0], sigval, d_tmp, nobs_tmp, prop_tmp, var_tmp, count_tmp);
   
/* find the deepest sigma_level for which there is data */

   deepestlev = -1;
   for (k = 0; k < size; ++k) {
      if (nobs_tmp[k] > 0 && d_tmp[k] >= 0)
           deepestlev = k;
   }
   
   if (size <= 1 || deepestlev < 0) {
      for (i = 0; i < nprops; ++i) {
           d = &dataptr[i][icol * NSTDLEVS];
           n = &countptr[i][icol * NSTDLEVS];
	   var = &varptr[i][icol *NSTDLEVS];
           for (j = 0; j < NSTDLEVS; ++j) {
               *(d++) = (float) HBEMPTY;
               *(n++) = 0;
	       *(var++) = (float) HBEMPTY;
           }
      }
      *bottom = (float) HBEMPTY;
      free((void *)sigval);
      free((void *)sig);
      free((void *)d_tmp);
      free((void *)nobs_tmp);
      for (i = 0; i < nprops; ++i) {
         free((void *)prop_tmp[i]);
         free((void *)var_tmp[i]);
         free((void *)count_tmp[i]);
      }
      free((void *)prop_tmp);
      free((void *)var_tmp);
      free((void *)count_tmp);
     
      return(NSTDLEVS);
   }

/* construct temporary arrays, x, y, v & w to store depth, y-property. v-variance and 
   w-# of obs for that property.  Enlarge the size of the arrays 
   to accommodate the mixed layer, seasonal layer and deepest observations.  */

   room = 100;
   x = (double *) malloc((size+room) * sizeof(double));
   y = (double *) malloc((size+room) * sizeof(double));
   v = (double *) malloc((size+room) * sizeof(double));
   w = (double *) malloc((size+room) * sizeof(double));
   if (w == NULL) {
     fprintf(stderr,"\nUnable to allocate memory in do_std_depth()\n");
     exit(1);
   }

   

/* Define the starting point in the sigma series which is heavier than
    the density at the bottom of the average mixed layer. Ensure that 
    the depth of this sigma level is deeper than the depth of the mixed 
    layer. */
    
   start = 0;
   m = g->mix_layer[imonth];
   if (m != NULL) {
   
      if (m->depth > d_tmp[deepestlev])
          m->depth = d_tmp[deepestlev];

      pref = 0.0;
      if (m->depth >= zmin[4])
         pref = 4000;
      else if (m->depth >= zmin[3])
         pref = 3000;
      else if (m->depth >= zmin[2])
         pref = 2000;
      else if (m->depth >= zmin[1]) 
         pref = 1000;
	
	/* revert to in situ temperature first */ 
      tref = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->depth);
      tref = hb_theta(m->prop[isa], tref, m->depth, pref);
      hb_svan(m->prop[isa], tref, pref, &sig_bml);  /* sigma at bottom of mixed layer */
       
      while ((sigval[start] < sig_bml) && (start < size))
           ++start;

      while ((m->depth > d_tmp[start]) && (start < size))
         ++start;
   }

/* Eliminate missing values from each y-property. Interpolate the property
   onto std depth surfaces.  Check for vertical datagaps and flag a standard
   depth as missing if the gap is too large.  Also interpolate to approximate 
   the # of observations and variance at each std level.  The interpolation routine returns 
   HBEMPTY when appropriate, so that each level is assigned a value
   whether or not any data is present.  The value of count at a level with 
   no data is assigned zero. */

   for (j = 0; j < nprops; ++j) {
      npts = 0;
      if (m != NULL) {
         if (m->wghtsum[j] > 0) {   
              npts = mixed_layer_vals(j, prop_id[j], m, x, y, v, w);
              hb_svan(m->prop[isa], m->prop[ite], 0.0, &sig[0]);
              for (i = 1; i < npts; ++i)
                  sig[i] = sig[0]; 
              sig[npts-1] = sig_bml;
         }  
      }

      if ((npts - start ) > room)  {
            fprintf(stderr, "\nIncrease the value of 'room' to at least %d", npts-start);
            fprintf(stderr, "\n in do_std_depth()\n");
            exit(1);
      }
	 
	 
      for (k = start; k < size; ++k) {
            if (count_tmp[j][k] > 0) {
	    
	          if(prop_tmp[j][k] < -3) 
	             fprintf(stderr,"bug in do_std_depth()\n" );
		  
                  y[npts] = prop_tmp[j][k];
	          v[npts] = var_tmp[j][k];
                  w[npts] = (double) count_tmp[j][k];
                  x[npts] = d_tmp[k];
                  sig[npts] = sigval[k];
                  ++npts;
            }  
       } 
	   
	   
       define_deepest_level(j, g->deepest);
       btmptr = g->deepest;
	 
       if (btmptr != NULL) {
 	 
	    /* Check deepest siglevel is shallower than 
	       top of bottom boundary layer.
	       Remove levels until this condition is met. */
	       
          *bottom = (float) btmptr->depth;
	  top_of_layer = btmptr->depth - btmptr->thickness;
	 
	 
	  k = npts-1;
	  while (k > start  &&  (x[k] > top_of_layer)) 
	         --k;
	   
	  npts = k+1;  /* index to next empty array element */
         
            /* add top and bottom of deepest layer to arrays....   */
         
	  reflev = 0;

	   if (top_of_layer >= zmin[4])
	       reflev = 4;
	   else if (top_of_layer >= zmin[3])
	       reflev = 3;
	   else if (top_of_layer >= zmin[2])
	       reflev = 2;
	   else if (top_of_layer >= zmin[1])
	       reflev = 1;
	
	   switch (reflev) {
	    case 0:
	       while ( (k > start) && sig[k] >= (btmptr->sig_0 - bottom_layer_def)  )
	          --k;
	       npts = k+1;
	       if (btmptr->thickness > 5)
	          sig[npts++] = btmptr->sig_0-bottom_layer_def;
	       sig[npts] = btmptr->sig_0;
	       break;
	    case 1:
	       while ( (k > start) && sig[k] >= (btmptr->sig_1 - bottom_layer_def)  )
	          --k;
	       npts = k+1;
	       if (btmptr->thickness > 5)
		       sig[npts++] = btmptr->sig_1-bottom_layer_def;
	       sig[npts] = btmptr->sig_1;
	       break;
	    case 2:
	       while ( (k > start) && sig[k] >= (btmptr->sig_2 - bottom_layer_def)  )
	          --k;
	       npts = k+1;
	       if (btmptr->thickness > 5)
		       sig[npts++] = btmptr->sig_2-bottom_layer_def;
	       sig[npts] = btmptr->sig_2;
	       break;
	    case 3:
	       while ( (k > start) && sig[k] >= (btmptr->sig_3 - bottom_layer_def)  )
	          --k;
	       npts = k+1;
	       if (btmptr->thickness > 5)
		       sig[npts++] = btmptr->sig_3-bottom_layer_def;
	       sig[npts] = btmptr->sig_3;
	       break;
	    case 4:
	       while ( (k > start) && sig[k] >= (btmptr->sig_4 - bottom_layer_def)  )
	          --k;
	       npts = k+1;
	       if (btmptr->thickness > 5)
		       sig[npts++] = btmptr->sig_4-bottom_layer_def;
	       sig[npts] = btmptr->sig_4;
	       break;
	   } /*end switch */
   
  /* add property value to arrays  */
	 
	   if (btmptr->thickness > 5) {
	     x[npts-1] = top_of_layer;
             y[npts-1] = btmptr->prop[j];
             v[npts-1] = btmptr->var[j];
	     w[npts-1] = (double) btmptr->count[j];
	   }
	   x[npts] = btmptr->depth;     
           y[npts] = btmptr->prop[j];
	   v[npts] = btmptr->var[j];
	   w[npts] = (double) btmptr->count[j];
           
           if (prop_id[j] == (int) PR) {       /* explicitly insert pressure if it is the y prop */
               y[npts] = btmptr->pressure;
	       if (btmptr->thickness > 5) 
		   y[npts-1] = btmptr->pressure - btmptr->thickness;
	   }
		
	   ++npts;   /* so far, npts has been an index var -- now make it a counter */
	   
	 } /* end if btmptr != NULL */
	 
         d = &dataptr[j][icol * NSTDLEVS];
         n = &countptr[j][icol * NSTDLEVS];
	 var = &varptr[j][icol * NSTDLEVS];
         
        for (i = 0; i < NSTDLEVS-1; ++i) {   /* interpolate all but bottom */
           if (npts <= 1) {
               z = (double) HBEMPTY;
	       zv = (double) HBEMPTY;
               k = 0;
           }
           else {
             z =  interpolate(std_depth[i], x, y, npts);
             k = 0;
	     zv = (double) HBEMPTY;
           }

           if (z > (HBEMPTY+10.0)) {                  /* check for vertical datagaps */
	   
       /*    if (z < -100.) 
	         fprintf(stderr,"interpolate returned %8.3g\n", z);
       */
		 
              jj = 0;
              while (x[++jj] < std_depth[i]) 
                   ;

              if (((int)x[jj-1] == (int)std_depth[i]) 
                  || ((int)x[jj] == (int)std_depth[i])) 
                  datagap = 0;

              else if (std_depth[i] < 1001) 
                  datagap = (x[jj] - x[jj-1]) > GAP_SHALLOW;

              else
                  datagap = (x[jj] - x[jj-1]) > GAP_DEEP;

              if (datagap && (ABS(sig[jj] - sig[jj-1]) >  0.04)) {      /* check for pycnostad */
                  z = (double) HBEMPTY;
                  if (j == 0 && report_gaps)
                    fprintf(stderr," datagap:  %.1lf  %.1lf\n", x[jj-1], x[jj]);
              }
              else {
                   r = interpolate(std_depth[i], x, w, npts);
                   if (r < 0)
                       k = 0;  /* this should never happen */
                   else {
                       k = (r - (int)r) > 0 ? (int) r + 1 : (int) r;
                   }
		   
		   /* weed out negative variances (from interpolation of missing values
		      and levels with fewer than 4 obs */
		      
                   zv = interpolate(std_depth[i], x, v, npts);
		   
		   if (zv < 0 || k < 4)
		      zv = (double) HBEMPTY; 
		   
              }
	     
           }
            if ((z > (HBEMPTY+10.0)) && (k == 0))
	         k = 1;      /* non-empty node gets at least 1 observation */      
           *(d++) = (float) z;
           *(n++) = (short) k;
	   if (zv > 0)
	     *(var++) = (float) sqrt(zv);
	   else
	     *(var++) = (float)HBEMPTY;
         }

/* add the deepest observation */
      if (btmptr != NULL) {
          *(d++) = (float) btmptr->prop[j];
          *(n++) = (short) btmptr->count[j];
           if ( btmptr->var[j] > 0)
             *(var++) = (float) sqrt(btmptr->var[j]);
           else
             *(var++) = (float) btmptr->var[j];
      }
   }  /* for j  */

/* clean up space no longer needed ... */

   free((void *)sigval);
   free((void *)sig);
   free ((void *)x);
   free ((void *)y);
   free ((void *)v);
   free ((void *)w);
   free((void *)d_tmp);
   free((void *)nobs_tmp);
   for (i = 0; i < nprops; ++i) {
         free((void *)prop_tmp[i]);
         free((void *)var_tmp[i]);
         free((void *)count_tmp[i]);
   }
   free((void *)prop_tmp);
   free((void *)var_tmp);
   free((void *)count_tmp);
   delete_surfrec_list(g->mix_layer[imonth]);
   g->mix_layer[imonth] = NULL;

   return(NSTDLEVS);

} /* end do_std_depth() */
/****************************************************************************/

int mixed_layer_vals(int iprop, int prop_id, struct surfrec *m, double *x, double *y, double *v, double *w)

   /* Creates arrays of depth, property, and counts at
      100 m intervals for the requested property.  Returns
      the number of points in each array. 
      
         iprop:          indicates the ith property of nprops 
       prop_id:       indicates the ith property of MAXPROPS 
             m:      info defining the mixed layer 
             x:       depth array 
             y:       property array 
	     v:       variance array
             w:       count array 
    ite, isa:       index to temperature, salt (or -1) 
 
      */
      
{
   int i, npts;
   double  *t, *s, pref;
   float deltap;

/* Define top and bottom only unless the mixed layer depth
   exceeds 200 m ... */

   npts = 2;       
   if (m->depth > 199.)
        npts = 1 + (int) m->depth / 100;

/*  Assign depth and count values ... */

   for (i = 0; i < npts-1; ++i) {
         w[i] = m->wghtsum[iprop];
         x[i] = i * 100.0;
	 v[i] = m->var[iprop];           
   }         
   w[npts-1] = m->wghtsum[iprop];
   v[npts-1] = m->var[iprop];
   x[npts-1] = m->depth;

/*  !**! Special cases for individual properties... */

/*    temperatures are theta, not in situ  */

   switch ((enum property) prop_id) {
       case PR :                   
                     for (i = 0; i < npts; ++i) {
                        y[i] = x[i];   
                     }         
                     break;
       case S0:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }

                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        pref = 0.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S1:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        pref = 1000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S2:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        pref = 2000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S3:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        pref = 3000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S4:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        pref = 4000.0;
                        compute_sigma(pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case S_:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        compute_sigma(s_pref, npts, y, x, t, s);
                        free(t);
                        free(s);
                        break;
       case VA:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        compute_svan(npts, y, x, t, s);
                        free(t);
                        free(s);
                       break;
       case SV:
                     if (ite < 0 || isa < 0) {
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
                     }
                        t = (double *) malloc(npts * sizeof(double));
                        s = (double *) malloc(npts * sizeof(double));
                        for (i = 0; i < npts; ++i) {
                           t[i] = hb_theta(m->prop[isa], m->prop[ite], 0.0, m->prop[ipr]);
                           s[i] = m->prop[isa];
                        }
                        compute_sp_vol(npts, y, x, t, s);
                        free(t);
                        free(s);
                       break;
       default:
                        for (i = 0; i < npts; ++i)
                          y[i] = m->prop[iprop];
                        break;
   }  /* end switch */

   return (npts);

}  /* end mixed_layer_vals() */
/****************************************************************************/
int define_deepest_level(int iprop, struct deepestrec *rptr)
   /* Traverse the linked list of near bottom property values to find 
    variance of property and number of obs.  Store values in rptr
    Determine thickness of the layer */
{
  int i, count;
  double delta, dprop, dpropsq,bottomdepth, bottomprop;
  struct deepestrec *r1;
  
  if (rptr == NULL)
     return 0;
  
  bottomdepth = rptr->depth;
  bottomprop = rptr->prop[iprop];
  dprop = dpropsq = 0.0;
  count = 0;     
  r1 = rptr;
  while (r1 != NULL) {
      delta = ABS(r1->depth - bottomdepth);
      if ((r1->prop[iprop] > -8.9) && (delta <= rptr->thickness)) {
         if (bottomprop < -8.9)
	     bottomprop = r1->prop[iprop];
         delta = ABS(r1->prop[iprop] - bottomprop);
	 dprop += delta;
	 dpropsq += delta * delta;
	 ++count;
      }
      r1 = r1->next;
  }
  rptr->prop[iprop]= bottomprop;
  rptr->var[iprop] = HBEMPTY;
  rptr->count[iprop] = count;
  if (count >= 4) {
     rptr->var[iprop] = (dpropsq - dprop * dprop/count) / (count-1);
  }
  return (count);

}  /* end define_deepest_level() */
/****************************************************************************/
int monotonic(double *din, UI  *nobs_in, double **xin, double **vin, UI  **count_in, int nprops, int npts_in, double *seas_d, UI *seas_dnobs, double **seas_prop, double **seas_var, UI **seas_count, int seas_npts, double *sigmas, double *dout, UI *nobs_out, double **xout, double **vout, UI **count_out )

/* sorts the inpurt arrays into increasing order by depth, removing levels with no obs and  ensuring there are no depth or density inversions. The seasonal layer (upper 200 m) is concatenated with the underlying sigma levels. Returns the # of depth levels in the output arrays. 
   
         din    starting addr of depth  
     nobs_in    starting addr of nobs array  
         xin    starting addr of other properties  
         vin    starting addr of variance arrays  
    count_in    starting addr of nobs arrays for other properties  
      nprops    # of rows (other properties) in xx  
     npts_in    # of cols (points) in each property array 
      seas_d    start of seasonal layer depths
  seas_dnobs    seasonal layer nobs at each level
   seas_prop    seasonal layer property arrays
    seas_var    seasonal layer variance arrays
  seas_count    seasonal layer nobs for each prop
   seas_npts    length of seasonal layer arrays
      sigmas    array to store output sigma values 
        dout    output depth array
    nobs_out    array of nobs at each depth level
        xout    arrays for other properties
        vout    arrays for variance of other props
   count_out    arrays for nobs of other props

   */
     
{
   int i, j, k, mono;
   int size;
   double sig_bsl;
   double  *diff;

   diff = (double *) calloc(npts_in,sizeof(double));
   if (diff == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in monotonic()\n");
     exit(1);
   }

   k = 0;
   diff[0] = 0.0;
   mono = 1;
   size = 0;
   sig_bsl = 0.0;

/* Start with isopycnally averaged levels in the seasonal layer */

   if (seas_npts > 0) {
      for (i = 0; i < seas_npts; ++i) {
         if (seas_dnobs[i] > 0){
             dout[k] = seas_d[i];
	     nobs_out[k] = seas_dnobs[i];
	     sigmas[k] = siglevs[i];
	     sig_bsl = sigmas[k];
             for (j = 0; j < nprops; ++j) {
               xout[j][k] = seas_prop[j][i];
               vout[j][k] = seas_var[j][i];
               count_out[j][k] = seas_count[j][i];
             }
	    if (k > 0) {
	       diff[k] = dout[k] - dout[k-1];
	       mono = mono && (diff[k] >= 0);
	    }
            ++k;
	 }    
      }
      
      /* now find first level in siglevs array that is deeper and denser
        than the bottom of seasonal layer */
      i = 0; 
      while ((i < npts_in) && (din[i] < bottom_seasonal_layer) )  
         ++i;
      while (i < npts_in && siglevs[i] < sig_bsl)
         ++i;
      size = i;          /* not really array size, just a temporary value */
   }
   
/* add levels deeper than seasonal layer */
            
   for (i = size; i < npts_in; ++i) {
      if (nobs_in[i] > 0 ) {
         dout[k] = din[i];
         nobs_out[k] = nobs_in[i];
         sigmas[k] = siglevs[i];
         for (j = 0; j < nprops; ++j) {
           xout[j][k] = xin[j][i];
           vout[j][k] = vin[j][i];
           count_out[j][k] = count_in[j][i];
         }
	 if (k > 0) {
	    diff[k] = dout[k] - dout[k-1];
	    mono = mono && (diff[k] >= 0);
	 }
         ++k;
      }
   }
   
   size = k;
   if ( mono || (size < 2)) {
       free((void *) diff);
       return(size);
   }
   
   /* Remove levels where depth inversions occur */

   while (!mono) {
     k = 0;
     for (i = 0; i < size; ++i) {
       if (diff[i] > -1 ) {
          dout[k] = dout[i];
          nobs_out[k] = nobs_out[i];
          sigmas[k] = sigmas[i];
          for (j = 0; j < nprops; ++j) {
             xout[j][k] = xout[j][i];
             vout[j][k] = vout[j][i];
             count_out[j][k] = count_out[j][i];
          }
	  ++k;
	}
	else if (i > 2) {
	     if ((dout[i] > dout[i-2]) && k > 0) {
	        --k;
                 dout[k] = dout[i];
                nobs_out[k] = nobs_out[i];
                sigmas[k] = sigmas[i];
                for (j = 0; j < nprops; ++j) {
                  xout[j][k] = xout[j][i];
                  vout[j][k] = vout[j][i];
                  count_out[j][k] = count_out[j][i];
                }
	        ++k;
	     }
	 
	     else {
	     
	       if (k >= 1) {
	          --k;
	          if ( i < (size-1)) {
	            diff[i+1] = dout[i+1] - dout[k];
	          }
	       }
	     }  
	 }
        
     } /* end for */
     
     size = k;
     mono = 1;
     diff[0] = 0.0;
     if (size > 2) {
        for (i = 1; i < size; ++i) {
	 diff[i] = dout[i] - dout[i-1];
	 mono = mono && (diff[i] >= -1);
        }
      }
   } /* end while !mono */
   
   free((void *) diff);

   return (size);

} /* end monotonic() */ 
/****************************************************************************/

double interpolate(double xval, double *x, double *y, int nypts)
/* Performs a linear interpolation to find the position of xval in array, x,
   and returns the corresponding value in the array, y.  If xval does not
   appear in array x, the value HBEMPTY is returned.  This routine assumes 
   that the x array is monotonic and continuous (no missing values); and it
   assumes the y array is continuous.    */
{
   int    k;
   double v1, v2;

   for (k = 0; k < nypts-1; ++k) {

      v1 = xval - x[k];
      v2 = xval - x[k+1];

      if (v1 == 0)             /* x[k] == xval */
          return (y[k]);
      if (v2 == 0)             /* x[k+1] == xval */
          return (y[k+1]);
      if (v1 < 0. && v2 < 0.)  /* xval not between x1 and x2 */  
          continue;
      if (v1 > 0. && v2 > 0.) 
          continue;

      return ( y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]) );
   }

   return ((double)HBEMPTY);

}   /* end interpolate() */

/****************************************************************************/

void sum_levels(double *y, double *d, int ny, double *ysum, double *ysq_sum, UI *count, double *seas_sum, double *seas_sqsum, UI *seas_count)
  /* Interpolates the y and d arrays with ny observation levels onto the stddepth levels
  (globally defined).  Checks for vertical data gaps and does not interpolate over
   gaps which exceed 200 m in the thermocline (upper 1000 m) or 1000 m 
   elsewhere. The interpolated values are added to the appropriate
  levels of the array defined by ysum -- or seas_sum for depths in the seasonal layer -- and the count array is incremented. */ 
{
   int  i, j, nz, datagap, n;
   double *ytmp, *dtmp, yint, flag;

   /*fprintf(stderr,"\n Started sum_levels...");*/

   ytmp = (double *) calloc((size_t)ny, sizeof(double));
   dtmp = (double *) calloc((size_t)ny, sizeof(double));
   if (dtmp == NULL) {
       fprintf(stderr,"\nUnable to allocate memory in sum_levels()\n");
       exit(1);
   }

   /*fprintf(stderr,"\n Allocated memory in sum_levels...");*/

   nz = NSTDLEVS - 1;   /* NSTDLEVS is globally defined and includes a bottom depth*/
   flag = HBEMPTY + 1.0;  
   
/* Ensure continuous (no missing values) for  y and d arrays 
 * by copying only valid data into ytmp and dtmp. Also, augment
 * n by one unit once we have copied the data. n stores the
 * number of valid points. */
   n = 0;
   for (i = 0; i < ny; ++i) {
      if ( y[i] > -8.9)  {
         ytmp[n] = y[i];
         dtmp[n++] = d[i];
      }
   }
   
   /* There are not enough values if there are at most 1 good observation
    * in this profile.  Cannot do any interpolation, so skip this one.*/
   if (n <= 1) {               /* not enough values */
     free((void *)ytmp);
     free((void *)dtmp);
     return;
   }

   /*fprintf(stderr,"\n nz= %i ",nz); 
     fprintf(stderr,"\n n= %i ",n);*/

   /* Loop over each standard level... */
   for (i = 0; i < nz; ++i) { 

     /* Interpolate the valid input data to this standard level.
      * The flag is HBEMPTY + 1 = -9 + 1 = -8, if value is not
      * present, then skip rest. */
     if ((yint = interpolate(std_depth[i], dtmp, ytmp, n)) > flag) {
      
	/* now check for datagaps around this depth...*/
	/*fprintf(stderr,"\n %i ",i);*/
	    
            j = 0;
            while (j < n && dtmp[j] < std_depth[i]) 
               ++j ;
	
	    /*fprintf(stderr,"Found j = %i",j);*/
	
	    datagap = 0;
            if ((dtmp[j] == std_depth[i]) )
                datagap = 0;
	    else if ( j > 0) {
	   
	       if  (dtmp[j-1] == std_depth[i]) 
                  datagap = 0;
               else if ( std_depth[i] < 1001)
                   datagap = (dtmp[j] - dtmp[j-1]) > GAP_SHALLOW;
               else
                   datagap = (dtmp[j] - dtmp[j-1]) > GAP_DEEP;
            }              
  
	    /*fprintf(stderr," datagap = %i",datagap);*/

            if (!datagap) {
	    
	       if (std_depth[i] <= bottom_seasonal_layer) {
		 /*fprintf(stderr," in seasonal layer ");*/
	          seas_sum[i] += yint;
		  seas_sqsum[i] += yint * yint;
		  ++seas_count[i];
	       }
	       else {
		 /*fprintf(stderr,"\n below seasonal layer ");
		 fprintf(stderr,"\n yint = %f",yint);
		 fprintf(stderr,"\n i = %i",i);
		 fprintf(stderr,"\n ysum[1] = %f",ysum[1]);
		 fprintf(stderr,"\n ysum[i] = %f",ysum[i]);
		 fprintf(stderr,"\n ysq_sum[i] = %f",ysq_sum[i]);
		 fprintf(stderr,"\n count[i] = %i",count[i]);
		 There is a strange seg fault error here
		 when -U is nonzero.  Ysum does not appear to
		 be allocated, will segfault with ysum[1] as
		 well as ysum[i] under some conditions.*/
		 ysum[i] += yint;
		 ysq_sum[i] += yint *yint;
		 ++count[i];

	       }
            }
       } /* end if */
   } /* end for i */
   
   /*fprintf(stderr,"\n End of work in sum_levels...");*/

   free((void *)ytmp);
   free((void *)dtmp);

   /*fprintf(stderr,"\n Done deallocating y and d in sum_levels...");*/

   return;

} /* end sum_levels() */
/************************************************************/
void do_pycnostads(int nc_file, struct gridnode ***grid, struct CDF_HDR *hptr,  int *prop_indx, int nstdlevs, int ntbins, int nrows, int ncols, int imonth)
   /* handles all the functions for writing data to the nc_file */
{
   int error, i, j, nprops, n_filled, nlevs, npts;
   int row, col, tbin;
   float  **data, **vari, *bottomdepth;  
   short  **count; 

   nprops = hptr->nprops;
   npts = nstdlevs * ncols;
   
   data = (float **) calloc(nprops, sizeof(float *));
   vari = (float **) calloc(nprops, sizeof(float *));
   count = (short **) calloc(nprops, sizeof(short *));
   bottomdepth = (float *) calloc (ncols, sizeof(float));

   for (i = 0; i < nprops; ++i) {
      data[i] = (float *) calloc(npts, sizeof(float));
      vari[i] = (float *) calloc(npts, sizeof(float));
      count[i] = (short *) calloc(npts, sizeof(short));
   }

   if (data[nprops-1] == NULL) {
      fprintf(stderr, "\nUnable to allocate memory for data & count arrays.\n");
      exit(1);
   }

   nlevs = nstdlevs - 1;

   n_filled = 0;

/* check each isopycnally averaged gridnode.
   If depth level is flagged as missing, check
   the isobarically averaged data to differentiate a vertical datagap
   from a pycnostad ... */

   for (tbin = 0; tbin < ntbins; ++tbin) {
      for (row = 0; row < nrows; ++row) {
         for (col = 0; col < ncols; ++col) {
	    read_cdf_bottom_depth(nc_file, bottomdepth, row, col, tbin);
            for (i = 0; i < nprops; ++i) {
               read_cdf_prop(nc_file, hptr->prop_id[i], data[i], row, col, tbin, 0, hptr->nz);
               read_cdf_prop_err(nc_file, hptr->prop_id[i], vari[i], row, col, tbin, 0, hptr->nz);
	       read_cdf_prop_count(nc_file, hptr->prop_id[i], count[i], row, col, tbin, 0, hptr->nz);
	       j = 0;
	       while (j < nlevs && std_depth[j] < *bottomdepth) {
		  if (is_flagged(data[i][j], (float)HBEMPTY)) {
		     if (std_depth[j] <= bottom_seasonal_layer) {
		        if (grid[tbin][row][col].month_count[imonth][i][j] > 0) {
		          count[i][j] = (short) grid[tbin][row][col].month_count[imonth][i][j];
		          data[i][j] = (float) (grid[tbin][row][col].month_prop[imonth][i][j] / (double) count[i][j]);
		          if (grid[tbin][row][col].month_count[imonth][i][j] > 3) {
		             vari[i][j] = (float) sqrt( ABS(grid[tbin][row][col].month_dpropsq[imonth][i][j] - (double) count[i][j] * (double)(data[i][j] * data[i][j])) / (double)(count[i][j]-1) );
		          }
		          if (i == 0) ++n_filled;
			}
		     }
		     
		     
		     else {
		        if (grid[tbin][row][col].count[i][j] > 0) {
		          count[i][j] = (short) grid[tbin][row][col].count[i][j];
		          data[i][j] = (float) (grid[tbin][row][col].prop[i][j] / (double) count[i][j]);
		          if (grid[tbin][row][col].count[i][j] > 3) {
		             vari[i][j] = (float) sqrt( ABS(grid[tbin][row][col].dpropsq[i][j] - (double) count[i][j] * (double)(data[i][j] * data[i][j])) / (double)(count[i][j]-1) );
		          }
		          if (i == 0) ++n_filled;
			}
		     }
		  } /* end if is_flagged */
		  ++j;
	       } /* end while */
               write_prop_cdf(nc_file, data[i], hptr->prop_id[i], row, col, tbin, 0,
                           1, 1, 1, nlevs);
               write_prop_err_cdf(nc_file, vari[i], hptr->prop_id[i], row, col, tbin, 0,
                           1, 1, 1, nlevs);
               write_prop_count_cdf(nc_file, count[i], hptr->prop_id[i],
                           row, col, tbin, 0, 1, 1, 1, nlevs);
	    } /* end for i */
	 } /* end for col */
      } /* end for row */
   } /* end for tbin*/
 
   fprintf(stderr,"\n  %d levels were identified and filled", n_filled);
   free((void *) data);
   free((void *) count);
   free((void *) bottomdepth);
   free((void *) vari);
   return;
} /* end do_pycnostads() */
