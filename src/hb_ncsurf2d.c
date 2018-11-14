/*  hb_ncsurf2d.c

_................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             
			     created for HB3M Jan 2010

................................................................................
____________________________________________________________________________
hb_ncsurf2dand extracts properties along one or more surfaces from HB3 gridded netcdf files. 

The surfaces may be defined by some value of any property supported
by HydroBase (depth, pressure, temperature, density, etc.)
The properties at each station are linearly interpolated onto the
surface; 
For each surface, an output file is generated containing :
     
     nrows, ncols, grid spacing, grid bounds         (1st line)
     list_of_properties                              (2nd line)
     lat lon   n1 p1 p1err  n2 p2 p2err ...  ( 1 line for each gridpt)
      where for each property:      n = # of obs
                                    p = mean 
                                 perr = sqroot of error variance    
                                                      
____________________________________________________________________________
___________________________________________________________________________
  USAGE:  

 hb_ncsurf2d filename(_roots) -B/west/east/south/north -I<gridspacing>  -S<surface_def_file> -P<property_list> [-D<dirname>] [-E<file_extent>] [-Z]

  list of filenames MUST be first argument or input is expected to come
  from stdin.

 -B : specifies grid bounds

 -I : specifies grid increment in degrees;  ex: -I0.5

 -S : file containing surface definitions. If this is not specified
      these values are expected to come from stdin .

 -P : list of properties to evaluate; ex: -Ppr/t90/th9/sa/ox
      -P (by itself) will print a list of the available properties;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

 -Z : use deepest occurance of surface (for cases where property defining
      surface is not monotonic).
*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hydro_cdf.h"
#include "hb_gamma.h"
#include "hb_paths.h"


/* input file pathnames */

#define    EXTENT   ""
#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    DIR      ""


/* set up data structures to store gridded data */

struct surface {
          double  value;    /* value of property on this surface */
          double  pref;
            int   data_ind;   /* index ID of property */
            char  density_flag ;  /* set if this is a density surface */
          double *prop[MAXPROP]; /* ptrs to grid of prop values */
          double *propsq[MAXPROP];
        unsigned *count[MAXPROP];
            FILE *fptr;    /* output file */
  struct surface *next;    
};

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int *prop_req;         /* list of props requested for output */
int nrequested;                /* count of above */	
int depth_req, pr_req;
int prop_needed[MAXPROP];  /* set of props requested || defining a
                                    surface || needed for computation */
short err_avail[MAXPROP];        /* error variance stored in .nc file */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double *err[MAXPROP];     /* store error variances */         
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
int use_deepest;        /* flag for non-monotonic properties */
float HB_f_mask;        /* flag for masked nodes in .cdf file */

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;            
int tindex;           /* specifies temperature scale for computed props*/   
  
     /* boundaries for grid */

float   xmin, xmax, ymin, ymax, delta_x, delta_y;  
int 	lon0to360;    /* = 1 for all positive longitudes, -1 for all negative, 0 for crossing Greenwich */
int *merid_range;  /* required by is_in_range() */
int     ncols, nrows;
int     tbin;

   /* prototypes for functions defined locally */
   
void print_usage(char *);
int parse_prop_list(char *);
struct surface *add_surf(FILE *, int);
void alloc_grids(struct surface *, unsigned int);
double project_prop(double *, int, double );
double project_deriv_prop(int, double, int, double);
void get_cdf_data(int, struct surface *);
void insert_counted_data(struct surface *, int, int, short **, double, short *);
void do_stats(struct surface *);

main (int argc, char **argv)
{
   short   bopt, iopt, popt, zopt;
   int     index;
   int     i, ncfile;
   int     curfile = 1, nfiles = 0; 
   int     n_extents, curext; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   unsigned gridsize;
   FILE    *def_file;
   char    *dir, *st;
   char   **extent_list;
   struct surface *list, *surf;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;
    extent_list[0] = EXTENT;
    dir = DIR;
    def_file = (FILE *) stdin;
    warnflag = 0;          /* set to 1 after warning is printed */
    tbin = 0;
    prompt = 1;            /* for surface definitions */
    bopt = iopt = popt  = 0;
    use_deepest = 0;
    error = 0;
    s_pref = -1;
    lon0to360 = 1;    /* all positive longitudes */
    window = 100;
    w_incr = 10;
    HB_f_mask = HBMASK;
    merid_range = (int *) NULL;
    list = (struct surface *) NULL;

/* initialize these ... */

   prop_req = (int *)calloc(MAXPROP, sizeof(int));
   nrequested = 0;
   depth_req = pr_req = 0;
   for (i = 0; i < MAXPROP; ++i) {
      prop_needed[i] = 0;
      err_avail[i] = 0;
      station.observ[i] = (double *) NULL;
      err[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
   }
   hdr.prop_id = (int *) NULL;
   tindex = (int)T90;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C':
                       fprintf(stderr,"\n WARNING: -C no longer an option\n");
		       break;
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;

               case 'E':       /* multiple file extents are supported */
                        extent_list[n_extents] = &argv[i][2];
			if (++n_extents == MAXEXT) {
                             fprintf(stderr,"\nToo many different file extents. Max allowed is %d\n", MAXEXT);
			 exit(1);
			}
                        break;
               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymax) != 1);
                        
                        if (xmin >= 0 && xmax < 0)
                            xmax += 360;
                        
			if (xmin < 0 ) {
			   lon0to360 = 0;
			   if (xmax <= 0)
			        lon0to360 = -1;
			}
			
                        break;

               case 'I':
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%f", &delta_y);
                        }
                        break;


               case 'S':
                        def_file = fopen(&argv[i][2],"r");
                        prompt = 0;              /* turn off prompt flag */
                        if (def_file == NULL) {
                           fprintf(stderr,"\nError opening %s.\n",&argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'P':
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        nrequested = parse_prop_list(&argv[i][2]);
                        if (nrequested <= 0)
                             error = 1;
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

               case 'Z':
                        use_deepest = 1;
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
             exit(1);
          }

       }  /* end if */

       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!bopt || !iopt ||  !popt) {
       fprintf(stderr,"\nYou must specify bounds, properties and gridspacing! \n");
       exit(1);
   }

   if (!nfiles) {
	  fprintf(stderr,"\nYou must specify a list of netcdf files for input. ");
           exit(1);
   }
   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin,ymax,xmin,xmax );
   fprintf(stderr,"\nUsing ITS-90 temperatures \n" );

   if (n_extents == 0) ++n_extents;

   /* compute dimensions of matrix formed by grid  */

   nrows = (int) (ceil((double)((ymax - ymin) / delta_y)) + .0001);
   ncols = (int) (ceil((double)((xmax - xmin) / delta_x)) + .0001);
   gridsize = nrows * ncols;

   fprintf(stderr,"\n Grid will be:  %d rows by %d cols\n",  nrows, ncols);


/* get info on each surface,  add surface to the linked list, 
    allocate space for computation, write heading to each output file ... */

   while ( (surf = add_surf(def_file,prompt))  != NULL) {
      surf->next = list;
      list = surf;
      alloc_grids(list, gridsize);
      fprintf(list->fptr,"%d %d %.2f %.2f %.3f %.3f %.3f %.3f %2d\n", 
             nrows, ncols, delta_x, delta_y, xmin, ymin, xmax, ymax, nrequested);
      for (i = 0; i < nrequested; ++i ) {
           fprintf(list->fptr, "%s ", get_prop_mne(prop_req[i]));
      }
      fprintf(list->fptr,"\n");
   }

   
   
   if (prop_needed[(int)GN] || prop_needed[(int)GE]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);


/* loop for each input file */

   do {

      curext = 0;
     do {

          ncfile = cdf_open(dir, argv[curfile], extent_list[curext], print_msg);
          if (ncfile >= 0 ){
     
             get_cdf_data(ncfile, list) ;
             cdf_close(ncfile);
	  }
     
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );


   fprintf(stderr,"\n    writing stats to outfile ...\n");

   do_stats(list);

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nExtracts properties from 3D gridded fields along horizontal surfaces.\n");
   fprintf(stderr,"\nUsage:  %s nc_file_list -B/west/east/south/north -I<delta_x[/delta_y]> -P<list_of_properties> [-S<surface_def_file>] [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>] [-Z]", program);

   fprintf(stderr,"\n\n  List of input *.nc filenames is first argument");
   fprintf(stderr,"\n -B  : specifies grid bounds");
   fprintf(stderr,"\n -I  : specifies grid increment in degrees.  ex: -I0.5");
   fprintf(stderr,"\n       If x and y increments differ, separate with");
   fprintf(stderr,"\n       a /       ex: -I0.5/1.0");
   fprintf(stderr,"\n -P  : list of properties to project onto surface;");
   fprintf(stderr,"\n       ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n            -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n -S  : file containing surface definitions.");
   fprintf(stderr,"\n       If this is not specified these values are expected");
   fprintf(stderr,"\n       to come from the standard input device.");
   fprintf(stderr,"\n\n    OPTIONS:");
   fprintf(stderr,"\n[-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n       ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.1deg.nc ");
   fprintf(stderr,"\n[-W] : Specifies window length and subdivision interval (both in db) for computing ");
   fprintf(stderr,"\n       gradient properties (bf, pv). ");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n[-Z] : use deepest occurrence of surface in each profile (default is use first occurrence of surface) ");
   fprintf(stderr,"\n[-h]  :  help -- prints this message.");
   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index, nprops, i;
   char prop[7];
   double ref_val;

   nprops = 0;

   do {
      if (*st == '/')
         ++st;
	 
      for(i = 0; i < 7; i++)
	    prop[i] = '\0';
      sscanf(st,"%[^'/']", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      
     if (index < MAXPROP) {   /* don't use count, variance, or quality arrays */
           prop_req[nprops] = index;
           prop_needed[index] = 1;
	   if (index == (int)PR)
	       pr_req = 1;
	   if (index == (int)DE)
	       depth_req = 1;
           ++nprops;
      }


     /* !**!  Special cases for properties */

      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == PE)) {
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }

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
      
      st += strlen(prop);

  } while (*st == '/');

  return (nprops);
}  /* end parse_prop_list() */
/*****************************************************************************/
struct surface *add_surf(FILE *infile,int  prompt)
/*  -  Allocates storage space for a record of type <struct surface>;
    -  queries the user for information about the surface (via stdin);
    -  sets the fields of the struct accordingly;
    -  opens the output file;
    -  returns a ptr to the struct  or NULL if no struct was added.
*/

{
   static n = 0;
   char id[6];
   int  index;
   struct surface *surf;
   char    fname[500];   


   if (!n && prompt) {
   fprintf(stderr,"\nDefine each projection surface ...\n");

   fprintf(stderr,"\n    -----  Surface property options  -------\n");
   print_prop_menu();
   fprintf(stderr,"\nbot:  (bottom observation)");
   fprintf(stderr,"\n\nend:  (end list)");
   fprintf(stderr,"\n    ----------------------------------------"); 

   }

   if (prompt)
         fprintf(stderr,"\nChoose a property type for surface #%d: ", ++n);
   if ( fscanf(infile, "%s", id) != 1) {
         fprintf(stderr,"\nError reading from surface definition file\n");
         exit(1);
   }
 
   if (strncmp(id, "end", 3) == 0) 
               return(NULL);          /* exit function */

   if (strncmp(id, "bot", 3) == 0) {
         surf = (struct surface *) calloc(1, sizeof(struct surface));
         surf->data_ind = -1;
         if (prompt) 
            fprintf(stderr,"Enter path/name of outfile for this surface: ");

         if (fscanf(infile, "%s", fname) != 1) {
            fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
            exit(1);
         }
   
         if ((surf->fptr = fopen(fname, "w")) == NULL) {
             fprintf(stderr,"\nError opening %s for output.", fname);
             exit(1);
         }

         return (surf);

   }

      index = get_prop_indx(id);
      if (index < 0 || index > MAXPROP) {       /* is choice appropriate? */
          fprintf(stderr,"\n%s is not an option!\n\n", id);
          exit(1);    
      }  /* end if */ 

 

       surf = (struct surface *) calloc(1, sizeof(struct surface));

       surf->data_ind = index;
       prop_needed[index] = 1;

  /*  and reference pressure, if appropriate... */

      if ( index == (int)S_  ){
         if (prompt) {
            fprintf(stderr,"enter reference Pr:  ");
         }
         if (fscanf(infile, "%lf", &surf->pref) != 1) {
           fprintf(stderr,"\nError reading pref for property %s from surface definition file\n", id);
           exit(1);
         }
      }
  /* get value of surface ... */

      if (prompt) {
         fprintf(stderr,"enter %s value", get_prop_descrip(index));
         fprintf(stderr,": ");
      }
      if (fscanf(infile, "%lf", &surf->value) != 1) {
        fprintf(stderr,"\nError reading surf.value from surface definition file\n");
        exit(1);
      }

       
/* !**!  Special cases for individual properties... */

       switch ((enum property) index)  {
          case GN: 
                 surf->pref = 0.;
                 surf->density_flag = 1;
                 break;
		 
          case S_:
                 surf->density_flag = 1;
                 if (s_pref >= 0) {   /* has this already been requested? */
                     if ( NINT( s_pref) != NINT(surf->pref)) {
                         fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                         fprintf(stderr,"You can only use one of those prefs");
                         fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                        exit(1);
                     }
                 }
                 s_pref = surf->pref;
                 break;
          case S0: 
                 surf->pref = 0.;
                 surf->density_flag = 1;
                 break;
          case S1: 
                 surf->pref = 1000.;
                 surf->density_flag = 1;
                 break;
          case S2: 
                 surf->pref = 2000.;
                 surf->density_flag = 1;
                 break;
          case S3: 
                 surf->pref = 3000;      
                 surf->density_flag = 1;
                 prop_needed[index] = 1;     
                 break;
          case S4: 
                 surf->pref = 4000.;
                 surf->density_flag = 1;
                 break;

          default:
                 surf->pref = 0.;       
                 surf->density_flag = 0;
                 break;   

       }  /* end switch */



/* open output file */

   if (prompt) 
      fprintf(stderr,"Enter path/name of outfile for this surface: ");

   if (fscanf(infile, "%s", fname) != 1) {
         fprintf(stderr,"\nError reading surf.filename from surface definition file\n");
         exit(1);
   }
   
   if ((surf->fptr = fopen(fname, "w")) == NULL) {
      fprintf(stderr,"\nError opening %s for output.", fname);
      exit(1);
   }

   return (surf);

}  /* end add_surf() */

/*****************************************************************************/
void alloc_grids(struct surface *s,unsigned int n)
/*  Allocates memory for each appropriate array in the struct * passed as 
.    an argument and initializes all values to 0.
.    The second argument contains the amount of memory to allocate.
.    The count array for PR is always initialized since it is the benchmark
.    which determines whether there are any stations at a gridpt.   */
{
   int i, indx;

   for (i = 0; i < MAXPROP; ++i ) {
      s->count[i] = NULL;
      s->prop[i] = NULL;
      s->propsq[i] = NULL;
   } /* end for */

   for (i = 0; i < nrequested; ++i) {
      indx = prop_req[i];
      if ((s->count[indx] = (unsigned *) calloc(n, sizeof(unsigned))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->count");
             exit(1);
      }
      if ((s->prop[indx] = (double *) calloc(n, sizeof(double))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->prop");
             exit(1);
      }
      if ((s->propsq[indx] = (double *) calloc(n, sizeof(double))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->propsq");
             exit(1);
      }
     } /* end for */
       
    if (s->count[(int)PR] == NULL)  
       s->count[(int)PR] = (unsigned *) calloc(n, sizeof(unsigned));

   return;

}  /* end alloc_grids() */

/****************************************************************************/
double project_prop(double *obs, int nobs, double p0)

   /*   Projects property values onto a specified surface...Returns the value
        of the property at the pressure level corresponding to this surface,
        or -9999 if the property cannot be projected.  The global pressure arrays 
        are used.  
        
            obs:   array of observations for property
           nobs:   number of observations
             p0:   interpolated pressure of surface
    */
{
   double val;
   double x[2], y[2];
   int    j;
   int   prevdepth, curdepth;


   prevdepth = curdepth = 0;
   j = 0;

   /* find first valid observation level. */

  while (obs[j] < -8.9) {
     if (++j == nobs)
          return (-9999.);
  }
  x[0] = station.observ[(int)PR][j];
  y[0] = obs[j];

 /* Now, check successive pairs of datapoints until the 
    surface is found or the last observation is encountered.   */

  while (++j < nobs) {
    if (obs[j]  > -8.9) {
        x[1] = station.observ[(int)PR][j]; 
        y[1] = obs[j];
        if ((val = hb_linterp(p0, x, y, 2)) > -9998.) {
              return (val);
        }
        x[0] = x[1];
        y[0] = y[1];
     }
  }  /* end while */

  return (-9999.);
} /* end project_prop() */
/****************************************************************************/
double project_deriv_prop(int index, double p0, int nobs, double lat)

   /* Computes the value of a derivative property at the specified pressure 
      level from the global station data arrays. 
      
       index:   identifies property to be projected 
          p0:   pressure at surface on which property will be projected 
        nobs:   number of p,t,s observation levels at this station 
         lat:   latitude of this station 

       */
{
   double y;

/*  !**!  Special cases for individual properties... */

   switch ((enum property)index) {
      case BF:
          y = buoyancy(p0, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], nobs, window, w_incr);
          if (y > -9998.)
             y = y * 1.0e5;
          break;
      case PV:
          y = buoyancy(p0, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], nobs, window, w_incr);
          if (y < -9998.0)
              break;
          y = potvort((y*y), lat);
          break;
      default:
          y = -9999.0;
          break;

   }  /* end switch */

   return (y);

}  /* end project_deriv_prop() */
/****************************************************************************/
               
void get_cdf_data(int cdfid, struct surface *listptr)
{
   struct CDF_HDR cdf;
   struct surface *surf;
   int startrow, startcol;
   int endrow, endcol;
   int cdf0to360, split, ratio_done, jj;
   float minlat, minlon;
   float maxlat, maxlon;
   float xoffset, yoffset;
   float lat, lon;
   double dlat, dlon;
   float *z, *x, *e;
   short **count;
   char *mne;
   int error, i, j, k, sq, row, col, npts;
   int index_prop;
   short prop_avail[MAXPROP];     /* set of props in .cdf file */
   short compute_flag[MAXPROP];   /* props NOT in .cdf which can be computed*/
   short tavail;


   if (error = read_cdf_hdr(cdfid, &cdf)) 
         exit (1);

   if (!cdf.counts_included) {
       fprintf(stderr,"FATAL ERROR: No observation counts stored in .nc file. \n"); 
       exit(1);
   }
   
/* compare bounds of file to range of surface grid to ensure overlap */


   cdf0to360 = 1;            /* cdf longitudes are all positive */
   if (cdf.xmin < 0) {
       cdf0to360 = 0;       /* cdf longitudes are mixed sign */
       if (cdf.xmax <= 0)
          cdf0to360 = -1;  /* cdf longitudes are all negative */
   }

   lon = cdf.xmin;
   while ((lon < cdf.xmax) && ! (i = is_in_range(lon, xmin, xmax, merid_range, lon0to360))) {
       lon += cdf.xincr;
   }

   if ((i == 0) || (cdf.ymin > ymax) || (cdf.ymax < ymin)) {
       fprintf(stderr, "\nNo data in cdf file is within grid bounds.\n");
       return;
   }

   /* adjust longitudes if necessary to match sign conventions */
   
   split = 0;
   if (cdf0to360 != lon0to360) {

      if (lon0to360 > 0) {   /* for surface grid lons all positive */
     
          if (cdf0to360 < 0) {  /*cdf lons all neg */
	       cdf.xmin += 360.0;
	       cdf.xmax += 360.0;
	       cdf0to360 = 1;
	  }
	  if (cdf0to360 == 0) 	{ /*cdf lons are mixed */
	       if (cdf.xmax < xmin) {
		   cdf.xmin += 360.0;
	           cdf.xmax += 360.0;
		   cdf0to360 = 1;
	       }
	       else {
		   if ( (xmax > (cdf.xmin+360)) && (xmax < (cdf.xmax+360)))
	              split = 1;
	       }
	  }       
      }
      else if (lon0to360 < 0) { /* for surface grid lons all negative */
          if (cdf0to360 > 0) {  /*cdf lons all pos */
	       cdf.xmin -= 360.0;
	       cdf.xmax -= 360.0;
	       cdf0to360 = -1;
	  }	       
	  if (cdf0to360 == 0) 	{ /*cdf lons are mixed */
	       if (cdf.xmin > xmax) {
		   cdf.xmin -= 360.0;
	           cdf.xmax -= 360.0;
	           cdf0to360 = -1;
	       }
	      else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
	   	 if ((xmin > (cdf.xmin-360)) && (xmin < (cdf.xmax-360))) 
	           split = 1;
	     }
	  }       
      }
      else {  /* lon0to360 == 0*/
          if (cdf0to360 > 0) {  /* cdf lons are are pos */
	     if (cdf.xmin > xmax) {
	        cdf.xmin -= 360;
		cdf.xmax -= 360;
		cdf0to360 = -1;
	     }
	     else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
	   	if ((xmin > (cdf.xmin-360)) && (xmin < (cdf.xmax-360))) 
	           split = 1;
	     }
	  }
	  else {  /* cdf lons are are neg */
	     if (cdf.xmax < xmin) {
	        cdf.xmin += 360;
		cdf.xmax += 360;
		cdf0to360 = 1;
	     }
	     else {
	       /* test for 2 split pieces of cdf file overlapping with grid */
		  if ( (xmax > (cdf.xmin+360)) && (xmax < (cdf.xmax+360)))
	              split = 1;
	     }
	  }
      }
   }
   
   xoffset = yoffset = 0.0;
   
   if (cdf.node_offset) {   /* for pixel gridded cdf file */
      xoffset = 0.5 * cdf.xincr;
      yoffset = 0.5 * cdf.yincr;
   }
   
for (jj = 0; jj <= split; jj++) {

   if (jj == 1) {
      if (cdf0to360 > 0) {
         cdf.xmin -= 360;
	 cdf.xmax -= 360;
	 cdf0to360 = -1;
      }
      else {
         cdf.xmin += 360;
	 cdf.xmax += 360;
	 cdf0to360 =  1;
      
      }
   }
   
   minlon = ((cdf.xmin + xoffset) < xmin) ? xmin : cdf.xmin + xoffset;
   maxlon = ((cdf.xmax - xoffset) > xmax) ? xmax : cdf.xmax - xoffset;
   minlat = ((cdf.ymin + yoffset) < ymin) ? ymin : cdf.ymin + yoffset;
   maxlat = ((cdf.ymax - yoffset) > ymax) ? ymax : cdf.ymax - yoffset;
   
   
   error = get_indices(&cdf, maxlat, minlon, &startrow, &startcol);
   error = get_indices(&cdf, minlat, maxlon, &endrow, &endcol);
   
   if (startrow < 0 ) startrow = 0;
   if (endrow >= cdf.ny) endrow = cdf.ny-1;
   if (startcol < 0 ) startcol = 0;
   if (endcol >= cdf.nx) endcol = cdf.nx-1;
   
/* determine which properties are available in the cdf file */

   for (i = 0; i < MAXPROP; ++i) {  /* initialize these flags */
      prop_avail[i] = 0;
      compute_flag[i] = 0;
      err_avail[i] = 0;
   }

   prop_avail[(int)DE] = 1;        /* depth is an index variable */
   prop_needed[(int)DE] = 1;
   prop_needed[(int)PR] = 1;      /* must always be in a cdf file */
   prop_needed[tindex] = 1;      
   prop_needed[(int)SA] = 1;      
   

   for (i = 0; i < cdf.nprops; ++i) {   
      prop_avail[get_prop_indx(cdf.prop_id[i])] = 1;
      if (cdf.counts_included > 1)
         err_avail[get_prop_indx(cdf.prop_id[i])] = 1;
   }
   
   tavail = prop_avail[(int)T90] || prop_avail[(int)TE];
   
/* determine which properties can and should be computed ...   */

/*  !**! Special cases for individual properties */

   for (i = 0; i < MAXPROP; ++i) {
      if (prop_needed[i] && !prop_avail[i]) {
         switch ((enum property) i) {
            case T90:
            case TH:
               compute_flag[i] =  tavail;
               prop_needed[(int)TE] = 1;
	       break;
            case TE:
            case TH9:
               compute_flag[i] =  tavail;
               prop_needed[(int)T90] = 1;
	       break;
            case OX: 
                compute_flag[i] = prop_avail[(int)O2] && prop_avail[(int)PR] &&
                                  tavail && prop_avail[(int)SA];
                if (compute_flag[i] == 0) {
                   fprintf(stderr,"No oxygens available in this file. \n"); 
                   prop_needed[i] = 0;
                }
                else {
                   prop_needed[(int)O2] = 1;
                   fprintf(stderr,"\n OX not available in cdf file, but will be computed from o2, pr, t90, sa values");
                }
                break; 
            case O2:
                compute_flag[i] = prop_avail[(int)OX] && prop_avail[(int)PR] &&
                                  tavail && prop_avail[(int)SA];
                if (compute_flag[i] == 0) {
                   fprintf(stderr,"No oxygens available in this file. \n"); 
                   prop_needed[i] = 0;
               }
                else {
                   prop_needed[(int)OX] = 1;
                   fprintf(stderr,"\n O2 not available in cdf file, but will be computed from ox, pr, t90, sa values");
               }
               break; 
            case S0: 
            case S1: 
            case S2:    /* fall through */
            case S3: 
            case S4: 
            case S_:
            case GN:
            case GE:
            case HT:
            case PE:
            case BF:
            case PV:
            case SV:
            case VS:
            case VA:
            case DR:
	    case AL:
	    case BE:
                compute_flag[i] = prop_avail[(int)PR] && tavail && prop_avail[(int)SA];

                if (compute_flag[i] == 0) {
                  fprintf(stderr,"\n\n FATAL ERROR!");
                  fprintf(stderr,"\n Unable to compute %s for this file.\n", 
                          get_prop_mne(i));
                  fprintf(stderr,"The .nc file must include pr, t90 (or te), sa \n"); 
                  exit (0);
                }

                fprintf(stderr,"\n %s not available in nc file, but will be computed from averaged p,t,s values.", get_prop_mne(i));
                break;

            default:
                fprintf(stderr,"\n\n FATAL ERROR!");
                fprintf(stderr,"\n  Property <%s> not available in this file.\n", get_prop_mne(i));
                exit(0);

         }  /* end switch */
      }
   }   /* end for */

/* get depth values */

   z = (float *) malloc ((size_t) cdf.nz * sizeof(float));

   if (read_cdf_depths(cdfid, z) < 0) {
      exit(1);
   }

/* allocate space for property, count & error arrays; initialize count arrays to 0 */

   for (i = 0; i < MAXPROP; ++i) {
      if (prop_needed[i]) {
         free_and_alloc(&station.observ[i], cdf.nz);
	 free_and_alloc(&err[i], cdf.nz);
      }
   }

   x = (float *) calloc ((size_t) cdf.nz, sizeof(float));  
   e = (float *) calloc ((size_t) cdf.nz, sizeof(float));
  
   count = (short **) calloc ((size_t) MAXPROP, sizeof(short *));
      if (count == NULL) {
          fprintf(stderr, "\nInsufficient memory for call to malloc()");
          exit(1);
      }
      for (i = 0; i < MAXPROP; ++i) {
         count[i] = NULL;
         if (prop_needed[i]) {
            count[i] = (short *) calloc((size_t) cdf.nz, sizeof(short));
            if (count[i] == NULL) {
              fprintf(stderr, "\nInsufficient memory for call to malloc()");
              exit(1);
            }
         }
      }
   
/* visit each appropriate gridpt in the cdf file  */

   for (row = startrow; row <= endrow; ++row) {
       for (col = startcol; col <= endcol; ++col) {

          
          /* convert lat/lon of this data from cdf file to the grid position
             for the surfaces. */

           error = get_lat_lon(&cdf, row, col, &lat, &lon);
           if ((lat > ymax) || (lat < ymin) || (lon > xmax) || (lon < xmin))
               continue;

           dlat = (double) lat;
	   dlon = (double) lon; 
	   
	   /* this assumes we are using pixel grid registration */
	     
           sq =  NINT((lat + .0001 - ymin) / delta_y - 0.5) * ncols + 
                 NINT((lon + .0001 - xmin) / delta_x - 0.5);

           /* get bottom depth for this square ...*/

           read_cdf_bottom_depth(cdfid, &z[cdf.nz-1], row, col, tbin);
           for (i = 0; i < MAXPROP; ++i) {

              /* extract available properties ... */

              if (prop_needed[i] && prop_avail[i]) {

                 mne = get_prop_mne(i);
                 error = read_cdf_prop(cdfid, mne, x, row, col, tbin, 0, cdf.nz);
                 if (error > 0) {
                    fprintf(stderr,"\nError attempting to read %s at row,col =  %d,%d from cdf file.", mne, row, col);
                    exit(1);
                 }
		 
		 if (i == (int)PR) {  /* check for roundoff errors in surface pressure */
		    if (x[0] < 1 && x[0] > -1)
		       x[0] = 0.0;
		 }
		 
                 if (i == (int)DE)   /* use pressure counts, errors for depth */
		    mne = get_prop_mne((int)PR);
		    
                 error = read_cdf_prop_count(cdfid, mne, count[i],
                             row, col, tbin, 0, cdf.nz);
                 if (error > 0) {
                       fprintf(stderr,"\nError attempting to read %s_cnt at row,col =  %d,%d from cdf file.", mne, row, col);
                       exit(1);
                 }
		 
		 if (err_avail[i]) {
		    error = read_cdf_prop_err(cdfid, mne, e, row, col, tbin, 0, cdf.nz);
                    if (error > 0) {
                       fprintf(stderr,"\nError attempting to read %s_%s at row,col =  %d,%d from cdf file.", mne, ERR_VAR_SUFFIX, row, col);
                       exit(1);
                    }
		 }

            /* place extracted data into appropriate column of station.observ */

                 for (j = 0; j < cdf.nz; ++j) {
                    station.observ[i][j] = (double) x[j];
		    if (err_avail[i]) {
		       err[i][j] = (double) e[j];
		    }
		 }
              } /* end if prop_needed */

           } /* end for i */

           /* choose the best property to determine whether any data is
              available at a given depth in .cdf file ... */

          if (prop_avail[(int)T90] )
               index_prop = (int)T90;
           else if (prop_avail[(int)TE])
               index_prop = (int)TE;
           else if (prop_avail[(int)SA] )
               index_prop = (int)SA;
           else if (prop_avail[(int)PR])
               index_prop = PR;
           else 
               index_prop = get_prop_indx(cdf.prop_id[0]);

           /* eliminate levels with no data */

           k = 0;
           for (i = 0; i < cdf.nz; ++i) {
              if (! (is_flagged((float)station.observ[index_prop][i], cdf.fill_value) 
	       || is_flagged((float)station.observ[index_prop][i], HB_f_mask)) ) {
                 for (j = 0; j < MAXPROP; ++j) {
                   if (prop_needed[j] && prop_avail[j] && (j != (int)DE)) {
                      station.observ[j][k] = station.observ[j][i];
                      count[j][k] = count[j][i];
		      if (err_avail[j]) {
		         err[j][k] = err[j][i];
			 if (is_flagged((float)err[j][k], cdf.fill_value))
			    err[j][k] = HB_MISSING;
		      }
		      if (is_flagged((float)station.observ[j][k], cdf.fill_value) || is_flagged((float)station.observ[j][k], HB_f_mask)) {
		        station.observ[j][k] = HB_MISSING;
			if (err_avail[j]) {
			   err[j][k] = HB_MISSING;
			}
		      }
                   }
                 }
                 station.observ[(int)DE][k] = (double) z[i];      
                 count[(int)DE][k] = count[index_prop][i];
                 ++k;
              }
           }
           if ((npts = k) <= 0)
                continue;

           /* next, compute properties that were flagged. Set count for
              each computed property.  Error variances do not apply  */
        
	/* First, take care of temperature */
	
	 i = (int)TE;
         if (prop_needed[i] && !prop_avail[i] && tavail) {
              t90_to_t68(station.observ[(int)T90],station.observ[(int)TE], npts);
	      compute_flag[i] = 0;
              for (j = 0; j < npts; ++j) 
                    count[i][j] = count[(int)T90][j];
			       
	      if (err_avail[(int)T90]) {
		 err_avail[i] = 1;
		for (j = 0; j < npts; ++j)    
		   err[i][j] = err[(int)T90][j];
	      }     
          }               

	 i = (int)T90;
         if (prop_needed[i] && !prop_avail[i] && tavail) {
              t68_to_t90(station.observ[(int)TE],station.observ[(int)T90], npts);
	      compute_flag[i] = 0;
              for (j = 0; j < npts; ++j) 
                    count[i][j] = count[(int)TE][j];
			       
	      if (err_avail[(int)TE]) {
		err_avail[i] = 1;
		for (j = 0; j < npts; ++j)    
		   err[i][j] = err[(int)TE][j];
	      }     
          }               

           ratio_done = 0;
      
/*  !**! Special cases for individual properties */
           for (i = 0; i < MAXPROP; ++i) {
              if (compute_flag[i]) {
                  switch ((enum property) i) {
                      case TH:
                        compute_theta(npts, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[(int)T90][j];
			       
			 if (err_avail[(int)T90]) {
			       err_avail[i] = 1;
                               for (j = 0; j < npts; ++j) 
			            err[i][j] = err[(int)T90][j];
			 }     
                         break;
                      case TH9:
                         compute_theta(npts, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[(int)T90][j];
			       
			 if (err_avail[(int)T90]) {
			      err_avail[i] = 1;
			      for (j = 0; j < npts; ++j)
			          err[i][j] = err[(int)T90][j];
			 }     
                         break;
                     case OX:
                        for (j=0; j < npts; ++j) {
                           station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[tindex][j],station.observ[(int)SA][j]);
                            count[i][j] = count[(int)O2][j];
                        }
                        break;
                      case O2:
                        for (j=0; j < npts; ++j) {
                           station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[tindex][j],station.observ[(int)SA][j]);
                           count[i][j] = count[(int)OX][j];
                        }   
                        break;
                   case S0: 
                         compute_sigma(0., npts, station.observ[(int)S0], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                         break;
                     case S1: 
                         compute_sigma(1000., npts, station.observ[(int)S1], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case S2:   
                         compute_sigma(2000., npts, station.observ[(int)S2], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case S3: 
                         compute_sigma(3000., npts, station.observ[(int)S3], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case S4: 
                         compute_sigma(4000., npts, station.observ[(int)S4], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case S_: 
                         compute_sigma(s_pref, npts, station.observ[(int)S_], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case HT:
                         compute_height(npts, station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], ht_pref, station.observ[(int)HT]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case PE:
                         compute_energy(npts, station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], pe_pref,  station.observ[(int)PE]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case SV: 
                         compute_sp_vol(npts, station.observ[(int)SV], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case VA: 
                         compute_svan( npts, station.observ[(int)VA], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA]);
                            for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                            
                         break;
                     case DR: 
		     case AL:
		     case BE:
		         if ( ! ratio_done ) {
                            compute_ratio( npts, station.observ[(int)DR], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
			 
                          for (j = 0; j < npts; ++j) {
			      if (prop_needed[(int)DR] )
                                   count[(int) DR][j] = count[index_prop][j];
			      if (prop_needed[(int)AL] )
                                   count[(int) AL][j] = count[index_prop][j];
			      if (prop_needed[(int)BE] )
                                   count[(int) BE][j] = count[index_prop][j];
                           }
                           
			   ratio_done = 1;
			 }
                         break;

                     case VS: 
                         compute_sound_vel( station.observ[(int)VS], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], npts);
                         for (j = 0; j < npts; ++j) 
                               count[i][j] = count[index_prop][j];
                         break;

                     case GN: 
		         if (!compute_flag[(int)GE]) { 
                            compute_gamma_n( &ginfo, npts, station.observ[(int)GN], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], dlon, dlat);
                             for (j = 0; j < npts; ++j) {
                               count[i][j] = count[index_prop][j];
                             }
			 }
                         break;

                     case GE:
	                compute_gamma_nerr(&ginfo, npts, station.observ[(int)GN],   
		          station.observ[(int)GE], station.observ[(int)PR],
		          station.observ[tindex], station.observ[(int)SA], dlon, dlat);
                          for (j = 0; j < npts; ++j) {
                             count[(int)GE][j] = count[index_prop][j];
                             count[(int)GN][j] = count[index_prop][j];
                          }
	                break;

                     default:
                         break;
              
                  }  /* end switch */

              }
           }

         /* traverse the linked list of surfaces & interpolate data onto each 
            surface */

           surf = listptr;
           while (surf != NULL) {
              insert_counted_data(surf, sq, npts, count, dlat, prop_avail);
              surf = surf->next;
           }  /* end while */

       }  /* end for col */
   }  /* end for row */

   free((void *) z);
   free((void *) e);
   free((void *) x);
   for (i = 0; i < MAXPROP; ++i) {
       if (station.observ[i] != NULL) {
          free((void *)station.observ[i]);
          station.observ[i] = NULL;
       }
       if (err[i] != NULL ) {
	  free((void *) err[i]);
	  err[i] = NULL;
       }
       if (count[i] != NULL) {
          free((void *)count[i]);
          count[i] = NULL;
       }
       
   }
   free((void *)count);
      
} /* end for jj */
   return;
}  /* end get_cdf_data() */
/****************************************************************************/

void insert_counted_data(struct surface *sptr, int sq, int nlevs, short **count, double lat, short *already_deriv)

   /*   Projects property values onto each surface and adds  values to the 
        appropriate grid element. The property and error variance arrays are globally defined,
        the address of the count arrays is passed as an argument. 
    */
{
   double  p0, val, last_val, err_val, r;
   double *x, *d;
   double sns[1], tns[1], pns[1], dummy;
   double dsns[1], dtns[1], dpns[1], glevels[1];
   int    nobs, datagap;                   
   int    j,k,index;

/*  bottom observations were requested... */

   if (sptr->data_ind == -1) {     
       for (j = 0; j < nrequested; ++j ) {

             index = prop_req[j];

             /* !**! Special cases for individual properties... */
             switch ((enum property) indx) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[index]) 
                       val = station.observ[index][nobs-1];
                    else
                       val = project_deriv_prop(index, station.observ[(int)PR][nobs-1], nobs, lat);
                    break;
              
               default:    /* all other properties */ 
                    if (station.observ[index] == NULL) {
                       val = -9999.;
                       break;
                    }
 
                    val = station.observ[index][nlevs-1];
		    if (err_avail[index]) {
		       err_val = err[index][nlevs-1];
		    }
		    
                    break;
             } /* end switch */

             if (val > -9998.) {
 	       nobs = ABS(count[index][nlevs-1]);
	       if (nobs == 0)     /* in case bottom has a zero count*/
	          nobs = 1;
               sptr->prop[index][sq] += nobs * val;
               sptr->count[index][sq] += nobs;
	       if (err_avail[index])
	           sptr->propsq[index][sq] += nobs * err_val;
               else
	           sptr->propsq[index][sq] += (nobs * val * val);
             }
          
        }  /* end for */
        return;
	
   } /* end if */
   
/* SURFACE is a neutral surface (gamma-n) */


   if (sptr->data_ind == (int)GN) {
   
      glevels[0] = sptr->value;
      
      neutral_surfaces(station.observ[(int)SA], station.observ[tindex], station.observ[(int)PR], station.observ[(int)GN], nlevs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
      if (pns[0] < 0.0) {   /* neutral surface not found */
      
         /* If it outcrops at sea surface, increment pressure and depth counters (effectively
	    add zero to the sums) and return */
	 
         if (sptr->value < station.observ[(int)GN][0] && station.observ[(int)PR][0] < 21) {
	    ++(sptr->count[(int)PR][sq]);
	    if (depth_req)
	       ++(sptr->count[(int)DE][sq]);
	 }
	 
	 return;
      }  /* end neutral surface not found */
      
      
      /* check for multiple xings of the neutral surface */
      
      if (dpns[0] != 0.0 && !warnflag) {
          fprintf(stderr,"\nWARNING: multiple crossings of gamma-n. The median will be used. Check <ns-multiples.dat> for additonal info.\n");
	  warnflag = 1;
      }
      
 /* find depth associated with pressure of neutral surface */
      
      p0 = pns[0];
      x = (double *) malloc(nlevs * sizeof(double));
      
    /* !**! Special cases for individual properties... */
      
      for (j = 0; j < nrequested; ++j) {
      
         index = prop_req[j];
	 
	 switch ((enum property) index)  {
	    
	       case DE: 
	          val = hb_depth(p0, lat);
		  break;
	       case T90:
	       case TE: 
	          val = tns[0];
		  break;
	       case SA: 
	          val = sns[0];
		  break;
	       case PR: 
	          val = pns[0];
		  break;
	       case GN: 
	          val = sptr->value;
		  break;
	       
	       case TH: 
	       case TH9: 
	          val = hb_theta(sns[0], tns[0], pns[0], 0.0);
		  break;
	    
	       case S0: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],0.0), 0.0, &val);
		  break;
		  
	       case S1: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],1000.0), 1000.0, &val);
		  break;
		  
	       case S2: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],2000.0), 2000.0, &val);
		  break;
		  
	       case S3: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],3000.0), 3000.0, &val);
		  break;
		  
	       case S4: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],4000.0), 4000.0, &val);
		  break;
		  
	       case S_: 
	          dummy = hb_svan(sns[0], hb_theta(sns[0],tns[0],pns[0],s_pref), s_pref, &val);
		  break;
		  
	       case SV: 
	          dummy = hb_svan(sns[0], tns[0], pns[0], &val);
		  val = 1.0e8 / (val + 1000.);
		  break;
		  
	       case VA: 
	          val = hb_svan(sns[0], tns[0], pns[0], &dummy);
		  break;
	       
	       case DR: 
	          val = hb_ratio(sns[0], tns[0], pns[0]);
		  break;
	       case AL: 
	          val = hb_alpha(sns[0], tns[0], pns[0]);
		  break;
	       case BE: 
	          val = hb_beta(sns[0], tns[0], pns[0]);
		  break;
		  
               case BF:    /* both derivative properties */
               case PV:
                   val = project_deriv_prop(index, pns[0], nlevs, lat);
		   for (k = 0; k < nlevs; ++k) {
		      count[index][k] = 1;
		   }
                   break;
		   
               default:   /* all other properties interpolate onto level of neutral surface */ 
                    if (station.observ[index] == NULL ) {
                       val = -9999.;
                       break;
                    }
                     val = project_prop(station.observ[index], nlevs, p0);
		    if (err_avail[index])
		      err_val = project_prop(err[index], nlevs, p0);
                    break;
  	    
	 }  /* end switch */
	 
         if (val > -9998.) {
	       /* determine number of obs which contributed to mean value */
	       
	       for  (k = 0; k < nlevs; ++k) {
	          x[k] = (double) ABS(count[index][k]);
	       }
	       d = station.observ[(int)PR];
	       r = hb_linterp(p0, d, x, nlevs);
	       nobs = (r - (int)r) > 0 ? (int)r + 1 : (int)r;
	       
               sptr->prop[index][sq] += (nobs * val);
               sptr->count[index][sq] += nobs;
	       if (err_avail[index] && nobs > 3)
	           sptr->propsq[index][sq] += nobs * err_val;
               else
	           sptr->propsq[index][sq] += (nobs * val * val);
         }
      } /* end for */
      
      if (!pr_req)                 /* even if pr is not being stored */
	++(sptr->count[(int)PR][sq]); /* the counter is always used */

      free((void *)x);
   
      return;
   }   


/* surface is not the bottom or gamma_n  ... */

   d = station.observ[(int)PR];     /* set this pointer for easier references */

 /* determine if the surface exists at this station
    and check for vertical datagaps.  An unacceptable datagap is defined
    in hydrobase.h:  GAP_SHALLOW
                     GAP_DEEP   */


   if ((p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, nlevs)) > -9998.) {


       /* after this loop, d[k-1] and d[k] will bracket the pressure value
          associated with the projection surface ... */

     k = 0;
     while (d[++k] < p0 ) 
            ;

     if (use_deepest && (k < (nlevs-1))) {   /* check for deeper surface */
        last_val = p0;
        while ( (p0 = hb_linterp(sptr->value, &station.observ[sptr->data_ind][k], &d[k], nlevs-k)) > -9998.) {
              last_val = p0;
              while ( d[++k] < p0 ) 
                  ;
        }
        p0 = last_val;
        k = 0;
        while (d[++k] < p0 )   /* find  d[i] associated with this deepest p0 */
            ;
     }
    
     if ((p0 - d[k-1] < 10) || (d[k] - p0 < 10) )
          datagap = 0;
     else if (p0 < 1001)
          datagap = (d[k] - d[k-1]) > GAP_SHALLOW;
     else
          datagap = (d[k] - d[k-1]) > GAP_DEEP;


      /* Exclude observations which are more than 1500 m above the 
         reference pressure, if one is specified.
         This prevents alot of bogus values from being added. */

      datagap = datagap || (p0 < sptr->pref-1500);

     if (!datagap) {

       --k;     /* k still points to index of surface in depth array */
       
        x = (double *) malloc(nlevs * sizeof(double));

       
        for (j = 0; j < nrequested; ++j ) {

          index = prop_req[j];
          switch ((enum property) index) {
                case BF:
                case PV:
                    if (already_deriv[index]) {
                       val = project_prop(station.observ[index], nlevs, p0);
		    }
                    else {
                       val = project_deriv_prop(index, p0, nlevs, lat);
		       for (k = 0; k < nlevs; ++k)
		          count[index][k] = 1;
		    }
                    break;
                case DE:
                    val = hb_depth(p0, lat);
                    break;
                case PR:
                    val = p0;
                    break;
                default:
                    val = project_prop(station.observ[index], nlevs, p0);
                    break;
          } /* end switch */

          if (val > -9998.) {
                 /* determine number of obs which contributed to mean value */
                  for (k = 0; k < nlevs; ++k) {
                      x[k] = (double) ABS(count[index][k]);
                  }
                  r = hb_linterp(p0, d, x, nlevs);
                  nobs = (r - (int) r) > 0 ? (int) r + 1: (int) r;
		  
                  sptr->prop[index][sq] += val * nobs;
                  sptr->count[index][sq]+= nobs;
		  
		  if (err_avail[index]) {
		      err_val = project_prop(err[index], nlevs, p0);
                      sptr->propsq[index][sq] += nobs * err_val;
		  }
		  else
                      sptr->propsq[index][sq] += val * val * nobs;

          }

        }  /* end for */
	
	if (!pr_req)
           ++(sptr->count[(int)PR][sq]);      /* the counter is always used */

        free((void *)x);
     } /* end if !datagap */
       
   }
   else {    /* surface not at this station */

     /* if this is a density surface, check if it outcrops at the sea surface.  
        Outcropping surfaces get a zero added to the pressure and depth arrays.
        First be sure that it isn't missing too many surface observations ... */

     if ((sptr->density_flag) && (d[0] < 21.)) {
       if (sptr->value < station.observ[sptr->data_ind][0]) {
           ++(sptr->count[(int) PR][sq]);    
           if (depth_req)
             ++(sptr->count[(int)DE][sq]);   
       }   
     }

   }

   return;

}  /* end insert_counted_data() */

/****************************************************************************/

void do_stats(struct surface *listptr)
{
   struct surface *surf;
   int    sq, row, col, j, index, n;
   int    print_it;
   float  lat, lon;
   float  mean, stddev;
   double var;

   surf = listptr;
   while (surf != NULL) {     /* traverse linked list of surfaces */
        sq = -1;      /* initialize, then simply visit each square*/

        for (row = 0; row < nrows; ++row) {
          for (col = 0; col < ncols; ++col) {
             lat = (row  + .5) * delta_y + ymin ;
             lon = (col + .5) * delta_x + xmin;
             ++sq;
	     
             print_it = 0;
             for (j = 0; j < nrequested; ++j ) { 
	       index = prop_req[j];             
               if (surf->count[index][sq] != 0) 
		     ++print_it;
	     }
	     
	     if (print_it) {
		  
               fprintf(surf->fptr,"%7.3f %8.3f", lat, lon);
               for (j = 0; j < nrequested; ++j ) {
                 index = prop_req[j];
                 mean = surf->prop[index][sq];
                 stddev = 0.;
                   if ( (n = ABS(surf->count[index][sq]) ) > 1) {    /*count can be negative if a first guess value */
                      mean = (float) surf->prop[index][sq] / (float) n;
		      
		      if (err_avail[index]) 
		         stddev = surf->propsq[index][sq] / (double) n;
		      else  { 
                         var = (surf->propsq[index][sq] - (double) (mean * mean * n)) / (double) (n-1);
                         stddev = (float) sqrt (ABS(var));
                      }
                   }
                 fprintf(surf->fptr,"  %8.4f %8.5f %4d", mean, stddev, n);
                 
               } /* end for */
               fprintf(surf->fptr,"\n");
             }

           }   /* end for */
        } /* end for */

       fclose(surf->fptr);
       surf = surf->next;

   }  /* end while */

   return;
}  /* end do_stats() */

