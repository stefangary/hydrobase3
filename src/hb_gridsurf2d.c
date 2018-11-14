/*  hb_gridsurf2d.c

_................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1993
                             Updated to ANSI standards Feb 2000
			     updated for HB3M Dec 2009

................................................................................
____________________________________________________________________________
hb_gridsurf2d projects hydrographic properties onto one or more surfaces
and computes mean and standard deviation values of those properties
for each node in a grid. The grid spacing and bounds are specified by 
the user. All points within a square are used to compute the mean and 
are weighted equally. Squares do not overlap.  The centers of each 
square form the x-y loci of elements of the matrix representing each
surface.

The surfaces may be defined by some value of any property supported
by HydroBase (depth, pressure, temperature, density, etc.)
The properties at each station are linearly interpolated onto the
surface; and the averaging, which produces the grid, is done on that
surface.  That is, if the surface is a 100 m depth surface, the
interpolated properties will be averaged on that depth surface.  
This produces results that are different from properties averaged on density surfaces and may introduce unrealistic watermass characteristics.
If the surface being gridded with gridsurf is of a type other than
density, we recommend that the data be gridded (averaged) first with
the modules hb_bin3d and hb_krig3d, and use that product as input to hb_ncsurf2d.

For each surface, an output file is generated containing :
     
     nrows, ncols, grid spacing, grid bounds         (1st line)
     list_of_properties                              (2nd line)
     lat lon   n1 p1 p1dev  n2 p2 p2dev ...  ( 1 line for each gridpt)
      where for each property:      n = # of obs
                                    p = mean 
                                 pdev = std dev    
                                                      
____________________________________________________________________________
___________________________________________________________________________
  USAGE:  

 hb_gridsurf2d filename(_roots)  -B/west/east/south/north -I<gridspacing>  -S<surface_def_file> -P<property_list> [-D<dirname>] [-E<file_extent>] [-Z]

  list of filenames MUST be first argument or input is expected to come
  from stdin.

 -B : specifies grid bounds

 -I : specifies grid increment in degrees;  ex: -I0.5

 -S : file containing surface definitions. If this is not specified
      these values are expected to come from stdin .

 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
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
int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
int use_deepest;        /* flag for non-monotonic properties */
float HB_f_mask;        /* flag for masked nodes in .cdf file */

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;    
int tindex;           /* specifies temperature variable for computed props*/        

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
void get_hydro_data(FILE *, struct surface *);
void insertdata(struct surface *, int, int, double, short *);
double project_prop(double *, int, double );
double project_deriv_prop(int, double, int, double);
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
   FILE    *def_file, *infile;
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
    merid_range = (int *) NULL;
    list = (struct surface *) NULL;

/* initialize these ... */

   prop_req = (int *)calloc(MAXPROP, sizeof(int));
   nrequested = 0;
   depth_req = pr_req = 0;
   for (i = 0; i < MAXPROP; ++i) {
      prop_needed[i] = 0;
      station.observ[i] = (double *) NULL;
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
               case 'N':
                       fprintf(stderr,"\n-C no longer an option. Use hb_ncsurf2d to extract surfaces from gridded *.nc files\n");
		       exit(0);
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

               case 'T':
	                switch (argv[i][2]) {
			   case '6':
			       tindex = (int)TE;
			       break;
			   case '9':   
			       tindex = (int) T90;
			       break;
			   default:
			       tindex = (int)TE;
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
       if (prompt) {
          fprintf(stderr,"\nYou must specify a surface definition file with -S");
	  fprintf(stderr,"\nwhen station files are input via stdin. ");
           exit(1);
       }
       
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = stdin;
   }
   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin,ymax,xmin,xmax );

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
           fprintf(list->fptr, "%2s ", get_prop_mne(prop_req[i]));
      }
      fprintf(list->fptr,"\n");
   }
   
   if (prop_needed[(int)GN] || prop_needed[(int)GE]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);


/* loop for each input file */

   do {

      curext = 0;
     do {
              
          if (nfiles > 0) {
            infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_msg);
            if (infile == NULL)
               goto NEXTFILE;
          }
     
            /* read each file completely */
         get_hydro_data(infile, list);

       if (nfiles > 0)
         fclose(infile);
NEXTFILE:
       ;
     
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );


/*  compute means and std deviations and write statistics to outfile ... */


   fprintf(stderr,"\n    writing stats to outfile ...\n");

   do_stats(list);

   fprintf(stderr,"\n\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nExtract properties along horizontal surfaces in gridded form\n ");
   fprintf(stderr,"\nUsage:  %s filename_root(s) -B/west/east/south/north -I<delta_x[/delta_y]> -P<list_of_properties> [-S<surface_def_file>] [-T<68|90>]  [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>] [-Z]", program);

   fprintf(stderr,"\n\n  Supply a list of filenames or station data is expected from stdin.");
   fprintf(stderr,"\n -B  : specifies grid bounds");
   fprintf(stderr,"\n -I  : specifies grid increment in degrees.  ex: -I0.5");
   fprintf(stderr,"\n       If x and y increments differ, separate with");
   fprintf(stderr,"\n       a /       ex: -I0.5/1.0");
   fprintf(stderr,"\n -P  : list of properties to project onto each surface;");
   fprintf(stderr,"\n       ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n            -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n    OPTIONS:");
   fprintf(stderr,"\n[-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n       ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : specifies input_file extent(s).  Maximum # of different extents is %d ", MAXEXT);  
   fprintf(stderr,"\n       (default is no extent)     ex: -E.dat ");
   fprintf(stderr,"\n[-S]  : file containing surface definitions.");
   fprintf(stderr,"\n       If this is not specified these values are expected");
   fprintf(stderr,"\n       to come from the standard input device.");
   fprintf(stderr,"\n[-T] : use IPTS-68  for computing derived variables (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
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
      
         if (index == (int)DE)
            depth_req = 1;
         if (index == (int)PR)
            pr_req = 1;

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
   int i,  k;

   for (i = 0; i < MAXPROP; ++i ) {
      s->count[i] = NULL;
      s->prop[i] = NULL;
      s->propsq[i] = NULL;
   } /* end for */

   for (k = 0; k < nrequested; ++k) {
      i = prop_req[k];
      if ((s->count[i] = (unsigned *) calloc(n, sizeof(unsigned))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->count");
             exit(1);
      }
      if ((s->prop[i] = (double *) calloc(n, sizeof(double))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->prop");
             exit(1);
      }
      if ((s->propsq[i] = (double *) calloc(n, sizeof(double))) == NULL) {
             fprintf(stderr,"\n\nUnable to allocate memory for s->propsq");
             exit(1);
      }

    } /* end for */
       
    if (s->count[(int)PR] == NULL)  
       s->count[(int)PR] = (unsigned *) calloc(n, sizeof(unsigned));

   return;

}  /* end alloc_grids() */

/****************************************************************************/
void get_hydro_data(FILE * file, struct surface *listptr)

   /*  Reads each station in a HydroBase file and adds property values
       to the appropriate surfaces.   This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j;
   int row, col, sq;
   int main_props_avail, ratio_done;
   double dlat;
   short derived[MAXPROP];
   struct surface  *surf;

/* initialize this so all derived properties get computed in insertdata() */

   for (i = 0; i < MAXPROP; ++i) {
       derived[i] = 0;
   }

/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

       /* is station within gridbounds? */
       if (hdr.lat <= ymax && hdr.lat >= ymin && is_in_range(hdr.lon, xmin, xmax, merid_range, lon0to360) ) {

          if (lon0to360 == 1) {
            if (hdr.lon < 0)
             hdr.lon += 360.0;
	   }  
	   else if (lon0to360 < 0) {
	      if (hdr.lon >= 0 )
	          hdr.lon -= 360.0;
	   }
	   else {   /* range crosses Greenwich */
	          if (hdr.lon < 0)   /* first make it positive */
	             hdr.lon += 360.;

	          if (hdr.lon > xmax ) /* make it neg if not in interval 0->xmax */
	             hdr.lon -= 360.;
	   }
          
   /* the extra smidgeon puts stations that fall on the border into
      the proper grid box */
      
       row = NINT( (hdr.lat + .0001 - ymin)/ delta_y - 0.5);
       col = NINT( (hdr.lon + .0001 - xmin)/ delta_x - 0.5);

  /* check again for  out of bounds stations   */

          if ((row < 0) || (row >= nrows) || (col < 0) || (col >= ncols))  {
             fprintf(stderr,"\n Software bug in get_hydro_data():  is_in_range() returned TRUE, but row,col out_of_bounds");
	     continue;
          }
         sq = row * ncols + col;
         
         
/* ensure that pr, de, te, and sa are available ... */

     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
	 
     if ( ! available(DE, &hdr) && available(PR, &hdr)) 
	 create_de_arrays( &hdr, &station);

	 
     switch (tindex) {
         case (int) T90:
             if (!available((int)T90, &hdr) && available(TE,&hdr) ) 
                create_t90_arrays( &hdr, &station);
	 break;
	 case (int)TE:
             if (!available(TE, &hdr) &&  available(T90,&hdr) ) 
                create_t68_arrays( &hdr, &station);
	 break;
     }
      
     main_props_avail = (available(PR, &hdr) && available((enum property)tindex, &hdr) && available(SA, &hdr));
     
     ratio_done = 0;

 /* compute appropriate properties at each level in station ... derivative 
    properties (pot.vorticity, buoyancy...) get computed later in 
    insertdata(). */

/* !**! Special cases for individual properties... */

    for (i = 0; i < MAXPROP; ++i) {
       if (prop_needed[i] && main_props_avail && !available((enum property)i, &hdr)) {
          switch ((enum property) i) {
             case OX:
               if (available(O2, &hdr)) 
                   create_ox_arrays( &hdr, &station);
               break;
             case O2:
               if (available(OX, &hdr)) 
                    create_o2_arrays( &hdr, &station);
               break;
	     case T90:
                create_t90_arrays( &hdr, &station);
	       break;
	     case TE:
                create_t68_arrays( &hdr, &station);
	        break;
             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;
             case TH9:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
               break;
             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], ht_pref, station.observ[i]);
               break;
             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], pe_pref, station.observ[i]);
               break;
             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case DR:
	     case AL:
	     case BE:
	       if (! ratio_done) {
                   free_and_alloc(&station.observ[(int)DR], hdr.nobs);
                   free_and_alloc(&station.observ[(int)AL], hdr.nobs);
                   free_and_alloc(&station.observ[(int)BE], hdr.nobs);
                   compute_ratio( hdr.nobs, station.observ[(int) DR], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
		   ratio_done = 1;
               }
	       break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
	       
             case VS:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs);
               break;


             case PV:   /* fall through */
             case BF:
                  prop_needed[(int)PR] = 1;
                  break;

             case GN:
	       if (!prop_needed[(int)GE]) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
	          compute_gamma_n(&ginfo, hdr.nobs, station.observ[i],
	           station.observ[(int)PR],station.observ[tindex],
		   station.observ[(int)SA], (double)hdr.lon, (double)hdr.lat);
	       }
	       break;
	       
             case GE:
               free_and_alloc(&station.observ[(int)GE], hdr.nobs);
               free_and_alloc(&station.observ[(int)GN], hdr.nobs);
	       compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[tindex], station.observ[(int)SA], 
		     (double)hdr.lon, (double)hdr.lat);
	       break;
            default:
               break;
          } /* end switch */
       } /* end if */

    }  /* end for */
      /* traverse the linked list of surfaces & interpolate data onto each 
         surface */

         dlat = (double) hdr.lat;
         surf = listptr;
         while (surf != NULL) {
	    
            insertdata(surf, sq, hdr.nobs, dlat, derived);
            surf = surf->next;
         }  /* end while */

       } /* end if is_in_range*/
       
      free_hydro_data(&station);

   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

/****************************************************************************/

void insertdata(struct surface *sptr, int sq, int nobs, double lat, short *already_deriv)

   /*   Projects property values onto a surface and adds each value to the 
        appropriate grid element.   The station data arrays are global. 
        
        arguments:
            sptr:   defines projection surface
              sq:   grid offset for this station 
            nobs:   number of observation levels in this station 
             lat:   latitude for this station 
   already_deriv:   flags derivative props which are already computed 
   */
{
   int j, index, datagap;
   double p0, val, last_val, *d;
   double sns[1], tns[1], pns[1], dummy;
   double dsns[1], dtns[1], dpns[1], glevels[1];
   double **propaddr, *prop;

/* return BOTTOM OBSERVATION ... */

   if (sptr->data_ind == -1) {     
   
        for (j = 0; j < nrequested; ++j ) {

             index = prop_req[j];

             /* !**! Special cases for individual properties... */
             switch ((enum property) index) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[index]) 
                       val = station.observ[index][nobs-1];
                    else
                       val = project_deriv_prop(index, station.observ[(int)PR][nobs-1], nobs, lat);
                    break;
              

               default:    /* all other properties */ 
	            propaddr = get_prop_array(&station, index);
	            if (*propaddr == NULL) {
                       val = -9999.;
                       break;
                    }
		    prop = *propaddr;
                    val = prop[nobs-1];
                    break;
             } /* end switch */


             if (val > -9998.) {
               sptr->prop[index][sq] += val;
               sptr->propsq[index][sq] += (val * val);
               ++(sptr->count[index][sq]);
             }
        }  /* end for */
	
	if (!pr_req)
	   ++sptr->count[(int)PR][sq];
        return;
   }
   
/* SURFACE is a neutral surface (gamma-n) */

   if (sptr->data_ind == (int)GN) {
   
      if (station.observ[(int)GN] == NULL) {
         fprintf(stderr,"\nWARNING: Unable to compute neutral surface.");
         fprintf(stderr,"   skipping...");
         return;
      }
   
      glevels[0] = sptr->value;
      
      neutral_surfaces(station.observ[(int)SA], station.observ[tindex], station.observ[(int)PR], station.observ[(int)GN], nobs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
      if (pns[0] < 0.0) {   /* neutral surface not found */
      
         /* If it outcrops at sea surface, increment pressure and depth counters
	    and return */
	 
         if (sptr->value < station.observ[(int)GN][0] && station.observ[(int)PR][0] < 21) {
	    ++(sptr->count[(int)PR][sq]);
	    if (depth_req)
	       ++(sptr->count[(int)DE][sq]);
	 }
	 
	 return;
      }  /* end neutral surface not found */
      
      
      /* check for multiple xings of the neutral surface */
      
      if (dpns[0] != 0.0 && !warnflag) {
          fprintf(stderr,"\nWARNING: multiple crossings of gamma-n. Median of the surfaces will be used. Check <ns-multiples.dat> for additonal info.\n");
	  warnflag = 1;
      }
      
      p0 = pns[0];
      
 /*   check for vertical datagaps.  An unacceptable datagap is defined
       as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP for all other surfaces   */
	       
      d = station.observ[(int)PR];  /* set this ptr for easier referencing */
      
        /* after this loop, d[j-1] and d[j] will bracket the pr value
          associated with the projection surface ... */
      j = 0;
      while (d[++j] < p0 ) 
            ;
	    
      if ((p0 - d[j-1] < 10) || (d[j] - p0 < 10) )
          datagap = 0;
      else if (p0 < 1001)
          datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
      else
          datagap = (d[j] - d[j-1]) > GAP_DEEP;

      if (datagap)  return;
	    
     
    /* !**! Special cases for individual properties... */
      
	 
      for (j = 0; j < nrequested; ++j) {
      
         index = prop_req[j];
	 switch ((enum property) index)  {
	       case DE: 
	          val = hb_depth(p0,lat);
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
	       case TH9: 
	       case TH: 
	          val = hb_theta(sns[0], tns[0], pns[0], 0.0);
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
	       
		  
               case BF:    /* both derivative properties */
               case PV:
                   val = project_deriv_prop(index, pns[0], nobs, lat);
                   break;
		   

               default:    /* all other properties interpolate onto level of neutral surface */ 
	            propaddr = get_prop_array(&station, index);
                    if (*propaddr == NULL ) {
                       val = -9999.;
                       break;
                    }
                    val = project_prop(*propaddr, nobs, p0);
                    break;
  	    
	    }  /* end switch */
	 
            if (val > -9998.) {
               sptr->prop[index][sq] += val;
               sptr->propsq[index][sq] += (val * val);
               ++(sptr->count[index][sq]);
            }
      } /* end for */

      if (!pr_req)              /* even if pr is not being stored */
        ++(sptr->count[(int)PR][sq]);      /* the counter is always used */

   
      return;
   } /* end if gamma-n */  

/* SURFACE IS NOT THE BOTTOM OR A NEUTRAL SURFACE... */

   if (station.observ[sptr->data_ind] == NULL) {
         fprintf(stderr,"\nWARNING: No pr,te,sa info to compute %2s surface", get_prop_mne(sptr->data_ind));
         return;
   }
   
   
   
   d = station.observ[(int)PR];  /* set this ptr for easier referencing */

 /* determine if the surface exists at this station
    and check for vertical datagaps.  An unacceptable datagap is defined
    as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP for all other surfaces   */

   
   if ((p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, nobs)) > -9998.) {
   
       /* after this loop, d[j-1] and d[j] will bracket the pr value
          associated with the projection surface ... */
     j = 0;
     while (d[++j] < p0 ) 
            ;

     if (use_deepest && (j < (nobs-1))) {   /* check for deeper surface */
        last_val = p0;
        while ( (p0 = hb_linterp(sptr->value, &station.observ[sptr->data_ind][j], &d[j], nobs-j)) > -9998.) {
              last_val = p0;
              while ( d[++j] < p0 ) 
                  ;
        }
        p0 = last_val;
        j = 0;
        while (d[++j] < p0 )   /* find  d[j] associated with this deepest p0 */
            ;
     }
    

     if ((p0 - d[j-1] < 10) || (d[j] - p0 < 10) )
          datagap = 0;
     else if (p0 < 1001)
          datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
     else
          datagap = (d[j] - d[j-1]) > GAP_DEEP;


      /* Exclude observations which are more than 1500 m above the 
         reference pressure, if one is specified.
         This prevents alot of bogus values from being added. */

      datagap = datagap || (p0 < sptr->pref-1500);

     if (!datagap) {
     
       --j;     /* j still points to index of surface in pressure array */
       
       /* note that use of j changes here! */

        for (j = 0; j < nrequested; ++j ) {  
            index = prop_req[j];

            /* !**! Special cases for individual properties... */

             switch ((enum property) index) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[index]) 
                       val = project_prop(station.observ[index], nobs, p0);
                    else
                       val = project_deriv_prop(index, p0, nobs, lat);
                    break;
 
               default:    /* all other properties */ 
 	            propaddr = get_prop_array(&station, index);
                    if (*propaddr == NULL ) {
                       val = -9999.;
                       break;
                    }
                   val = project_prop(*propaddr, nobs, p0);
                   break;
             } /* end switch */

             if (val > -9998.) {
               sptr->prop[index][sq] += val;
               sptr->propsq[index][sq] += (val * val);
               ++(sptr->count[index][sq]);
             }
           
          
        }  /* end for */
	
	if (!pr_req)                       /* even if pr is not being stored */
           ++(sptr->count[(int)PR][sq]);      /* the counter is always used */

      }   /* end if  !datagap */
       
   }  /* if p0 > -9998 */
   else {    /* surface not at this station */

     /* if this is a density surface, check if it outcrops at the sea surface. 
        Outcropping surfaces get a zero added to the pressure and depth
        arrays (operationally their counters just get incremented). 
        First be sure that the station isn't missing too many surface
        observations ... */

     if ((sptr->density_flag) && (d[0] < 21.)) {
       if (sptr->value < station.observ[sptr->data_ind][0]) {
           ++(sptr->count[(int) PR][sq]);    
           if (depth_req)
             ++(sptr->count[(int)DE][sq]); 
       }   
     }

   }

   return;

}  /* end insertdata() */

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
          y = buoyancy(p0, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], nobs, window, w_incr);
          if (y > -9998.)
             y = y * 1.0e5;
          break;
      case PV:
          y = buoyancy(p0, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], nobs, window, w_incr);
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

void do_stats(struct surface *listptr)
{
   struct surface *surf;
   int    sq, row, col, i, n, index;
   int    print_it;
   float  lat, lon;
   double  mean, stddev;
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
             for (i = 0; i < nrequested; ++i ) {
	       index = prop_req[i];
               if ((surf->count[index][sq] != 0)) 
		     ++print_it;
	     }
	     
	     if (print_it) {
		  
               fprintf(surf->fptr,"%7.3f %8.3f", lat, lon);
               for (i = 0; i < nrequested; ++i ) {
	           index = prop_req[i];
                   mean = surf->prop[index][sq];
                   stddev = 0.;
                   if ( (n = ABS(surf->count[index][sq]) ) > 1) {    /*count can be negative */
                      mean = (double) surf->prop[index][sq] / (double) n;
                      var = (surf->propsq[index][sq] - (double) (mean * mean * n)) / (double) (n-1);
                      stddev = sqrt (ABS(var));
                   }
                   fprintf(surf->fptr,"  %8.4lf %.5g %4d", mean, stddev, n);
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

