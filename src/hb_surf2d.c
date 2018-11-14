/*  hb_surf2d.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             1993
			     updated for HB3M Dec 2009
................................................................................
 
hb_surf2d reads a HydroBase station file, projects hydrographic properties 
onto one or more surfaces and outputs the values of those properties together 
with optional year, month, lat, and lon for each station in file.

The surfaces may be defined by some value of any property supported
by HydroBase (depth, pressure, temperature, density, etc.)
The properties at each station are linearly interpolated onto the
surface;

For each surface, an output file is generated containing :
     
        one or more
     {year month lat lon station_id}   p1 ... pn  ( 1 line for each station)
                                                      
____________________________________________________________________________
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



/* set up data structure to store info about each surface */

struct surface {
          double  value;    /* value of property on this surface */
          double  pref;
            int   data_ind;   /* index ID of property */
            char  density_flag ;  /* set if this is a density surface */
            FILE *fptr;    /* output file */
  struct surface *next;    
};
    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int *prop_req;         /* list props requested for output */
int nrequested;                /* count of above */				 int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
   

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double ht_pref;          /* ref lev for computing dynamic height */
double pe_pref;          /* ref lev for computing potential energy */
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;     
int tindex;              /* temperature scale for computing derived values */

     /* boundaries  */

float   xmin, xmax, ymin, ymax;     
int  xdateline;
int yflag, mflag, lflag, iflag, dopt;


   /* prototypes for locally defined functions */
   
struct surface *add_surf(FILE *, int);
void    print_usage(char *);
void    get_hydro_data(FILE *, struct surface *);
int     parse_prop_list(char *);
double project_deriv_prop(int, double, int, double);
double project_prop(double *, int, struct surface *);
double project_on_depth(double *, int, double);


main (int argc, char **argv)
{
   short   bopt, popt;
   int     index;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int n_extents, curext; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   FILE    *def_file, *infile;
   char    *dir, *st;
   char  **extent_list;
   struct surface *list = NULL, *surf;


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
    prompt = 1;            /* for surface definitions */
    bopt  = popt = dopt = 0;
    error = 0;
    xmin = -360;
    xmax = 360;
    ymin = -90;
    ymax = 90;
    s_pref = -1;
    xdateline = 0;
    window = 100;
    w_incr = 10;
    yflag = mflag = lflag = iflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */

/* initialize these ... */

   prop_req = (int *)calloc(MAXPROP, sizeof(int));
   nrequested = 0;
   for (i = 0; i < MAXPROP; ++i) {
       prop_needed[i] = 0;
       station.observ[i] = (double *)NULL; 
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
   }
   hdr.prop_id = (int *) NULL;
   tindex = (int)T90;



/* are there command line arguments? */

   if (argc <= 1) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
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
                        
                        if (xmin > 0 && xmax < 0)
                           xmax += 360;
                           
                        if (xmax > 180)
                           xdateline = 1;
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


               case 'S':
                        def_file = fopen(&argv[i][2],"r");
                        prompt = 0;              /* turn off prompt flag */
                        if (def_file == NULL) {
                           fprintf(stderr,"\nError opening %s.\n",&argv[i][2]);
                           exit(1);
                        }
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
               case 'M':
                        mflag = 1;
			if (argv[i][2] == 'd')
			    dopt = 1;
                        break;
               case 'Y':
                        yflag = 1;
                        break;
               case 'L':
                        lflag = 1;
                        break;
               case 'I':
                        iflag = 1;
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

   if (  !popt) {
       fprintf(stderr,"\nYou must specify properties with the -P flag.\n");
       exit(1);
   }

   if (  !(mflag || yflag || lflag || iflag)) {
       fprintf(stderr,"\nWARNING: You have not specified any output parameter from {month|year|position|station_id}.  Only properties will be output.\n");
   }
   
   if (!nfiles) {
       if (prompt) {
          fprintf(stderr,"\nYou must specify a surface definition file with -S");
	  fprintf(stderr,"\nwhen station files are input via stdin. ");
           exit(1);
       }
       
       fprintf(stderr,"\nExpecting input from stdin ... ");
   }

   fprintf(stderr,"\n\n  bounds: %.3f %.3f %.3f %.3f\n", ymin, ymax, xmin, xmax );

/* get info on each surface,  add surface to the linked list, 
    allocate space for computation, write heading to each output file ... */

   while ( (surf = add_surf(def_file, prompt))  != NULL) {
      surf->next = list;
      list = surf;
   }

   fclose(def_file);
   
    
   if (prop_needed[(int)GN] || prop_needed[(int)GE]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

  
   infile = stdin;
   if (n_extents == 0) ++n_extents;
   

/* loop for each input file */

   do {

     curext = 0;
     do {
         if (nfiles > 0) {
             infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_msg);
             if (infile == NULL)
                goto NEXTFILE;
         }

      get_hydro_data(infile, list);

      if (nfiles > 0) 
       fclose(infile);
       
NEXTFILE:
        ;
	
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\nEnd of %s\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
 fprintf(stderr,"\n%s projects hydrographic properties onto one or more surfaces", program);
 fprintf(stderr,"\nand outputs the values of those properties together with");
 fprintf(stderr,"\nsome combination of year, month, lat, and lon for each profile in a HydroBase station file.");
 fprintf(stderr,"\nSurfaces may be defined by some value of any property supported");
 fprintf(stderr,"\n(e.g. depth, pressure, temperature, density, etc.)");
 
   fprintf(stderr,"\n\nUsage:  %s filename_root(s) -P<list_of_properties> [-B/west/east/south/north] [-D<dirname>] [-E<file_extent>]  [-S<surface_def_file>][-T<68|90>] [-W<window>[/<w_incr>]] [-Y] [-M[d]] [-L] [-I]", program);

   fprintf(stderr,"\n\nList of filenames MUST be first argument.");
   fprintf(stderr,"\n  If no files are named, station data expected from stdin.");
   fprintf(stderr,"\n   -P  list of properties to project onto surface;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n   OPTIONS:");
   fprintf(stderr,"\n   -B  specifies grid bounds: w/e/s/n");
   fprintf(stderr,"\n   -D  specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   -E  specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n   -I  include station ID in output listing");
   fprintf(stderr,"\n   -L  include lat/lon in output listing");
   fprintf(stderr,"\n   -M  include month  in output listing.  Append 'd' to also output day. ex: -Md");
   fprintf(stderr,"\n   -Y  include year in output listing");
   fprintf(stderr,"\n   -S  file containing surface definitions.");
   fprintf(stderr,"\n       If this is not specified these values are expected");
   fprintf(stderr,"\n       to come from the standard input device.  This file");
   fprintf(stderr,"\n       must be specified if station data are input via stdin.");
   fprintf(stderr,"\n   -T  use IPTS-68  for computing derived variables (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
   fprintf(stderr,"\n   -W  Specifies pressure window (db) length and ");
   fprintf(stderr,"\n       subdivisions for computing gradient properties (bf and pv) ");
   fprintf(stderr,"\n          defaults: -W100/10 ");
   fprintf(stderr,"\n   -h  help (prints this message)");  

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
      prop_req[nprops] = index;
      if (index < MAXPROP)
           prop_needed[index] = 1;
      ++nprops;

      /* !**!  Special cases for properties */

      if (((enum property)index == S_ )|| ((enum property) index == HT) || ((enum property) index == PE)) {
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
struct surface *add_surf(FILE *infile, int prompt)
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
   char    fname[80];   


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
         surf = (struct surface *) malloc(sizeof(struct surface));
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
      if (index  < 0 || index > MAXPROP) {       /* is choice appropriate? */
          fprintf(stderr,"\n%s is not an option!\n\n", id);
          exit(1);    
      }  /* end if */ 


       surf = (struct surface *) malloc(sizeof(struct surface));

       surf->data_ind = index;
       prop_needed[index] = 1;

  /*  and reference pressure, if appropriate... */

      if ( (index== (int)S_ ) || (index == (int)HT) || (index == (int)PE)){
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
         fprintf(stderr,"enter %s value", get_prop_mne(index));
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

          case HT: 
                 ht_pref = surf->pref;       
                 surf->density_flag = 0;
                 break;   
	         
          case PE: 
                 pe_pref = surf->pref;       
                 surf->density_flag = 0;
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

void get_hydro_data(FILE *file, struct surface *listptr)
   /*  Reads each station in a HydroBase file and adds property values
       to the appropriate surfaces.   This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, main_props_avail, ratio_done, j;
   double dlat;
   short derived[MAXPROP];
   struct surface  *surf;
   void outputdata(struct surface *, short *);

/* initialize this to flag all derived properties */

   for (i = 0; i < MAXPROP; ++i) 
            derived[i] = 0;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check for and skip over out of bounds stations   */
  
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;

       if ((hdr.lat > ymax) || (hdr.lat < ymin) 
        || (hdr.lon < xmin) || (hdr.lon > xmax))  
             continue;

         
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
    properties (pot.vorticity, buoyancy...) get computed later . */

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
       }

    }

      /* traverse the linked list of surfaces & interpolate data onto each 
         surface */

         surf = listptr;
         while (surf != NULL) {
            outputdata(surf, derived);
            surf = surf->next;
         }  /* end while */


      free_hydro_data(&station);
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

/****************************************************************************/

void outputdata(struct surface *sptr, short *already_deriv)

   /*   Projects property values onto a surface and writes values to output 
        file.   The station data arrays are global. 
	
	arguments:
                 sptr: defines projection surface 
        already_deriv: flags derivative props which are already computed 
   */
{
   int j, index, datagap;
   char field_string[20];
   double p0, d0, val, dlat, *d;
   double sns[1], tns[1], pns[1], dummy;
   double dsns[1], dtns[1], dpns[1], glevels[1];
   double **propaddr, *prop;

   if (sptr->data_ind == -1) {     /* return bottom observations */
        if (yflag)
         fprintf(sptr->fptr,"%5d ", hdr.year);
        if (mflag)
          fprintf(sptr->fptr,"%2d ", hdr.month);
	if (dopt)
          fprintf(sptr->fptr,"%2d ", hdr.day);
	
        if (lflag)
           fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);
        if (iflag)
          fprintf(sptr->fptr,"%5d ", hdr.station);

        for (j = 0; j < nrequested; ++j ) {

             index = prop_req[j];
/* !**! Special cases for individual properties... */

             switch ((enum property) index) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[index]) 
                       val = station.observ[index][hdr.nobs-1];
                    else
                       dlat = (double) hdr.lat;
                       val = project_deriv_prop(index, station.observ[(int)PR][hdr.nobs-1], hdr.nobs, dlat);
                    break;
              

               default:    /* all other properties */ 
	            propaddr = get_prop_array(&station, index);
	            if (*propaddr == NULL) {
                       val = -9.9;
                       break;
                    }
		    prop = *propaddr;
                    val = prop[hdr.nobs-1];
                    break;
             } /* end switch */

              sprintf(field_string," %c%d.%dlf", '%', get_field_width(index), get_field_precis(index));
             fprintf(sptr->fptr, field_string, val);
        }  /* end for */
        fprintf(sptr->fptr, "\n");
        return;
   }
   
   
 /* SURFACE is a neutral surface (gamma-n) */

   if (sptr->data_ind == (int)GN) {
   
      if (station.observ[(int)GN] == NULL) {
         fprintf(stderr,"\nWARNING: No pr,te,sa info to compute neutral surface.");
         fprintf(stderr,"   skipping...");
         return;
      }
   
      glevels[0] = sptr->value;
      
      neutral_surfaces(station.observ[(int)SA], station.observ[tindex], station.observ[(int)PR], station.observ[(int)GN], hdr.nobs, glevels, 1, sns, tns, pns, dsns, dtns, dpns);
      
       if (pns[0] < 0.0) {   /* neutral surface not found */
	 return;
      }  
      
       /* check for multiple xings of the neutral surface */
      
      if (dpns[0] != 0.0 && !warnflag) {
          fprintf(stderr,"\nWARNING: multiple crossings of gamma-n. Median of the surfaces will be used. Check <ns-multiples.dat> for additonal info.\n");
	  warnflag = 1;
      }
      
      d0 = hb_linterp(pns[0], station.observ[(int)PR], station.observ[(int)DE], hdr.nobs);
      
 /*   check for vertical datagaps.  An unacceptable datagap is defined
       as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP for all other surfaces   */
	       
      d = station.observ[(int)DE];  /* set this ptr for easier referencing */
      
        /* after this loop, d[i-1] and d[i] will bracket the de value
          associated with the projection surface ... */
      j = 0;
      while (d[++j] < d0 ) 
            ;
	    
	    
      if ((d0 - d[j-1] < 10) || (d[j] - d0 < 10) )
          datagap = 0;
      else if (d0 < 1001)
          datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
      else
          datagap = (d[j] - d[j-1]) > GAP_DEEP;

      if (datagap)  return;
     
      if (yflag)
         fprintf(sptr->fptr,"%5d ", hdr.year);
      if (mflag)
          fprintf(sptr->fptr,"%2d ", hdr.month);
	if (dopt)
          fprintf(sptr->fptr,"%2d ", hdr.day);
      if (lflag)
           fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);
      if (iflag)
          fprintf(sptr->fptr,"%5d ", hdr.station);
     
    /* !**! Special cases for individual properties... */
      
      for (j = 0; j < nrequested; ++j) {
      
         index = prop_req[j];
	 
	    switch ((enum property) index)  {
	    
	       case DE: 
	          val = d0;
		  break;
	    
	       case TE: 
	       case T90: 
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
                   dlat = (double) hdr.lat;
                   val = project_deriv_prop(index, pns[0], hdr.nobs, dlat);
                   break;
		   

               default:    /* all other properties interpolate onto depth of neutral surface */ 
	            propaddr = get_prop_array(&station, index);
                    if (*propaddr == NULL ) {
                       val = -9.9;
                       break;
                    }
                    val = project_on_depth(*propaddr, hdr.nobs, d0);
                    break;
  	    
	    }  /* end switch */
	    
	    
            sprintf(field_string," %c%d.%dlf", '%', get_field_width(index), get_field_precis(index));
            fprintf(sptr->fptr, field_string, val);
      } /* end for */
      
      fprintf(sptr->fptr, "\n");
   
      return;
   }   


/* surface is not the bottom or a neutral surface ... */

   d = station.observ[(int)DE];  /* set this ptr for easier referencing */

 /* determine if the surface exists at this station
    and check for vertical datagaps.  An unacceptable datagap is defined
    as :       > GAP_SHALLOW for surfaces in the upper 1000 m
               > GAP_DEEP  for all other surfaces   */

   if ((d0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], d, hdr.nobs)) > -9998.) {

     if (prop_needed[(int)PR])    /* pr is required for derivative props */
        p0 = hb_linterp(sptr->value, station.observ[sptr->data_ind], station.observ[(int)PR], hdr.nobs);

       /* after this loop, d[j-1] and d[j] will bracket the de value
          associated with the projection surface ... */
     j = 0;
     while (d[++j] < d0 ) 
            ;

     if ((d[j-1] == d0) || (d[j] == d0) )
          datagap = 0;
     else if (d0 < 1001)
          datagap = (d[j] - d[j-1]) > GAP_SHALLOW;
     else
          datagap = (d[j] - d[j-1]) > GAP_DEEP;


      /* Exclude observations which are more than 1500 m above the 
         reference pressure, if one is specified.
         This prevents alot of bogus values from being added, but may cause 
         results to be different than expected if user is unaware of this. */

      datagap = datagap || (d0 < sptr->pref-1500);

     if (!datagap) {
        if (yflag)
         fprintf(sptr->fptr,"%5d ", hdr.year);
        if (mflag)
          fprintf(sptr->fptr,"%2d ", hdr.month);
	if (dopt)
          fprintf(sptr->fptr,"%2d ", hdr.day);
        if (lflag)
           fprintf(sptr->fptr,"%8.3f %8.3f ", hdr.lon, hdr.lat);
        if (iflag)
          fprintf(sptr->fptr,"%5d ", hdr.station);

        for (j = 0; j < nrequested; ++j ) {

            index = prop_req[j];

/* !**! Special cases for individual properties... */

             switch ((enum property) index) {

               case BF:    /* derivative properties */
               case PV:
                    if (already_deriv[index]) 
                       val = project_prop(station.observ[index], hdr.nobs, sptr);
                    else {
                       dlat = (double) hdr.lat;
                       val = project_deriv_prop(index, p0, hdr.nobs, dlat);
                    }
                    break;

               default:    /* all other properties */ 
	            propaddr = get_prop_array(&station, index);
                    if (*propaddr == NULL ) {
                       val = -9.9;
                       break;
                    }
                    val = project_prop(*propaddr, hdr.nobs, sptr);
                    break;
             } /* end switch */

             sprintf(field_string," %c%d.%dlf", '%', get_field_width(index), get_field_precis(index));
             fprintf(sptr->fptr, field_string, val);
        }  /* end for */
        fprintf(sptr->fptr, "\n");

      }   /* end if  !datagap */
       
   }

   return;

}  /* end outputdata() */

/****************************************************************************/
double project_prop(double *obs, int nobs, struct surface *sptr)

   /*   Projects property values onto a specified surface...Returns the value
        of the property at the pressure level corresponding to this surface,
        or -9.0 if the property cannot be projected.  The global arrays 
        pr & de are used.  
     arguments:
        obs:    array of observations for property 
        nobs:   number of observations 
        sptr:   specifies the surface on which to project the prop 
  */
{
   double val;
   double x[2], y[2];
   int    j;
   int   prevdepth, curdepth;

   prevdepth = curdepth = 0;
   j = 0;

   /* find first valid observation level. */

  while ((obs[j] < -8.9)||(station.observ[sptr->data_ind][j] < -8.9)) {
     if (++j == nobs)
          return (HB_MISSING);
  }
  x[0] = station.observ[sptr->data_ind][j];
  y[0] = obs[j];

 /* Now, check successive pairs of datapoints until the 
    surface is found or the last observation is encountered.   */

  while (++j < nobs) {
    if ((obs[j] > -8.9) && (station.observ[sptr->data_ind][j] > -8.9)) {
        x[1] = station.observ[sptr->data_ind][j]; 
        y[1] = obs[j];
        if ((val = hb_linterp(sptr->value, x, y, 2)) > -9998.) {
           return (val);
        }
        x[0] = x[1];
        y[0] = y[1];
     }
  }  /* end while */

  return (HB_MISSING);
} /* end project_prop() */
/****************************************************************************/
double project_deriv_prop(int index, double p0, int nobs, double lat)

   /* Computes the value of a derivative property at the specified pressure 
      level from the global station data arrays. 
      
    arguments:
      index:   identifies property to be projected 
      p0:      pressure at surface on which property will be projected 
      nobs:    number of p,t,s observation levels at this station 
      lat:     latitude of this station 
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
          y = HB_MISSING;
          break;

   }  /* end switch */

   return (y);

}  /* end project_deriv_prop() */

/****************************************************************************/
double project_on_depth(double *obs, int nobs, double d0)

   /*   Projects property values onto the depth surface (d0)...Returns the value
        of the property at the level corresponding to this surface,
        or HB_MISSING if the property cannot be projected.  The global arrays 
        station.observ are used.  
        
            obs:   array of observations for property
           nobs:   number of observations
             d0:   interpolated depth of surface
    */
{
   double val;
   double x[2], y[2];
   int    j;


   j = 0;

   /* find first valid observation level. */

  while (obs[j] < -8.9) {
     if (++j == nobs)
          return (HB_MISSING);
  }
  x[0] = station.observ[(int)DE][j];
  y[0] = obs[j];

 /* Now, check successive pairs of datapoints until the 
    surface is found or the last observation is encountered.   */

  while (++j < nobs) {
    if (obs[j]  > -8.9) {
        x[1] = station.observ[(int)DE][j]; 
        y[1] = obs[j];
        if ((val = hb_linterp(d0, x, y, 2)) > -9998.) {
              return (val);
        }
        x[0] = x[1];
        y[0] = y[1];
     }
  }  /* end while */

  return (HB_MISSING);
} /* end project_on_depth() */
