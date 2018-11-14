/*  hb_columns.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1995
                             updated to ANSI-C Feb 2000
			     formerly hb_4matlab
			     changed for HB3M Dec 2009
			     
................................................................................
____________________________________________________________________________
hb_columns computes properties at each pressure level in a station and 
outputs a columnar listing of each observed level containing all the properties
specified with the -P option plus user-specified info: { lat/lon year month station/cruise} ____________________________________________________________________________
  USAGE:  

 hb_columns filename(_roots) -P<property_list> [-L] [-Y] [-M] [-S] [-W<window>[/<w_incr>] [-D<dirname>] [-E<file_extent>] [-O<output_filename>]

  list of filenames must be first argument or input is expected from stdin
  
 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
      -P (by itself) will print a list of the available properties;
      
 -L : output lat/lon of observed level

 -Y : output year of observed level

 -M : output month of observed level

 -S : output station#/cruise# of observed level

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 
 -O : output file name

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
#define    DIR      ""
#define    MAXEXT   10  /* maximum # of different file extents specified */

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int prop_req[2*MAXPROP];         /* set of props requested for output */
int nrequested;                  /* count of above */
int prop_needed[MAXPROP];      /* set of props requested
                                     or needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
int tindex;              /* index to temperature variable */

int yopt, lopt, sopt, mopt, dopt, uopt;   /* option flags */
int station_id;          /* for use with -U */
FILE *outfile;

struct GAMMA_NC ginfo;


  /* prototypes for locally defined functions */
  
void print_usage(char *);
int parse_prop_list(char *);
void get_hydro_data(FILE *);

main (int argc, char **argv)
{
   short   popt;
   int     index;
   int n_extents, curext; 
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *dir, *st;
   char  **extent_list;
   FILE *infile;


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
    popt = yopt = lopt = sopt = mopt = dopt = uopt = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    outfile = stdout;
    station_id = 10000000;  /* Allow for 9 million profiles in a single file! */
    
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      prop_req[i] = 0;
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

               case 'L':
                        lopt = 1;
                        break;
                        
               case 'Y':
                        yopt = 1;
                        break;
                        
               case 'M':
                        mopt = 1;
			if (argv[i][2] == 'd')
			    dopt = 1;
                        break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			
			}
                        break;
                        
                        
               case 'S':
                        sopt = 1;
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
			       error = 1;
			}
			break;
	       case 'U':
                        uopt = 1;
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

   if ( !popt) {
       fprintf(stderr,"\nYou must specify a list of properties to output.\n");
       exit(1);
   }

   if (tindex == (int) T90)
      fprintf(stderr,"Using T90 temperatures\n");
   else
      fprintf(stderr,"Using T68 temperatures\n");

   if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = stdin;
   }

   if (! ( lopt || mopt || sopt || yopt || uopt)) {
     fprintf(stderr,"\nWARNING! no header info (year, month, position, cruise or station id) \nhas been specified for output.\n");
   }
   

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
            /* read each file completely */

 
     get_hydro_data(infile);

     if (nfiles > 0) 
         fclose(infile);

NEXTFILE:
            ;

     } while (++curext < n_extents);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\nEnd of hb_columns.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{  
 fprintf(stderr,"\n%s computes properties specified with -P option at each depth level in a station \n and outputs a column listing of each level plus any optional lat/lon/year/month/station/cruise of the observation.\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<list_of_properties>  [-D<dirname>] [-E<file_extent>] [-O<output_file>] [-L][-M[d]] [-S] [-T68] [-Y] [-W<window>/<w_incr>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -P  : list of properties to output;");
   fprintf(stderr,"\n          ex:  -Ppr/th/sa/ox/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n   [-L] : option to output lat/lon ");
   fprintf(stderr,"\n   [-Y] : option to output year ");
   fprintf(stderr,"\n   [-M] : option to output month.  Append d to optionally output day: -Md ");
   fprintf(stderr,"\n   [-S] : option to output station#/cruise# ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n   [-O] : output filename, default is output to stdout");
   fprintf(stderr,"\n   [-T] : Set temperature scale used in computing derived variables to either IPTS-68 (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
   fprintf(stderr,"\n   [-W] : Specifies pressure window length and subdivisions (db) for computing ");
   fprintf(stderr,"\n          gradient properties (bf, pv)  ");
   fprintf(stderr,"\n          defaults: -W100/10");
   fprintf(stderr,"\n   [-U] : output a unique cast number for each cast written.");
   fprintf(stderr,"\n          this is currently simply the count of each cast in");
   fprintf(stderr,"\n          the .hb file that is input to this routine.");
   fprintf(stderr,"\n   [-h] : help...... prints this message. \n");
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

  /* some special cases ... */

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
  
 
 /* !**!  Special cases for properties */
  
  if (prop_needed[(int)GE] && !prop_needed[(int)GN]) {
     prop_needed[(int)GN] = 1;
  }
  
  if (prop_needed[(int)GE] || prop_needed[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
  if (prop_needed[(int)TH])
     prop_needed[(int)TE] = 1;
     
  if (prop_needed[(int)TH9])
     prop_needed[(int)T90] = 1;
     
     prop_needed[tindex] = 1;
     
 /* end of special cases */
  
  return (nprops);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(FILE * file)

   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.    */
{
   int error, i, j, nbytes, pts_avail, ratio_done, len_str;
   int index;
   double dlon, dlat, deltap;
   double *scan, *prop;
   char *str, *st;

     str = calloc(100, sizeof(char));

/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check if basic properties pr, te, and sa are available ... */

     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
	 
     if (!available(T90, &hdr) && available(TE,&hdr) )  
         create_t90_arrays( &hdr, &station);

     if (prop_needed[(int)TE] && !available(TE,&hdr)) 
        create_t68_arrays( &hdr, &station);
	
    pts_avail = (available(PR, &hdr) && available((enum property)tindex, &hdr) 
                 && available(SA, &hdr));

    ratio_done = 0;
    
 /* compute appropriate properties at each level in station ... */

/* !**! Special cases for individual properties... */

    for (i = 0; i < nrequested; ++i) {
       index = prop_req[i];
       if (pts_avail && !available((enum property)index, &hdr)) {
          switch ((enum property) index) {
	    case DE:
                create_de_arrays( &hdr, &station);
	        break;
             case OX:
	        if (available(O2, &hdr)) 
                   create_ox_arrays( &hdr, &station);
                break;
             case O2:
               if (available(OX, &hdr)) 
                   create_o2_arrays( &hdr, &station);
               break;
             case TH:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;
             case TH9:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
               break;
             case S0:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S1:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S2:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S3:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S4:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case S_:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case HT:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], ht_pref, station.observ[index]);
               break;
             case PE:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], pe_pref, station.observ[index]);
               break;
             case SV:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;

             case VA:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;
             case DR:
             case AL:
             case BE:
	        if (! ratio_done ) {
                    free_and_alloc(&station.observ[(int)DR], hdr.nobs);
                    free_and_alloc(&station.observ[(int)AL], hdr.nobs);
                    free_and_alloc(&station.observ[(int)BE], hdr.nobs);
                    compute_ratio( hdr.nobs, station.observ[(int)DR], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
		     ratio_done = 1;
		 }
               break;

             case VS:
               free_and_alloc(&station.observ[index], hdr.nobs);
               compute_sound_vel( station.observ[index], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs);
               break;

             case PV: 
               free_and_alloc(&station.observ[index], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[index], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[index][j] > -99.) 
                   station.observ[index][j] *= station.observ[index][j];
               }
               po_vort(station.observ[index], station.observ[index], hdr.nobs, dlat);

               break;
               
             case BF:
               free_and_alloc(&station.observ[index], hdr.nobs);
               buoy_freq(station.observ[index], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[index][j] > -9998.)   /* convert the units */
                        station.observ[index][j] *= 1.0e5;
               }
               break;

             case GN:
                  free_and_alloc(&station.observ[index], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_n(&ginfo, hdr.nobs, station.observ[index], station.observ[(int)PR],station.observ[tindex],station.observ[(int)SA], dlon, dlat);
	       break;
	       
             case GE:
	       if (pts_avail) {
                  free_and_alloc(&station.observ[(int)GE], hdr.nobs);
                  free_and_alloc(&station.observ[(int)GN], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
		  compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[tindex], station.observ[(int)SA], dlon,
		     dlat);
               }
	       break;
             default:
	         /* all other properties (including count, variance, quality) must 
		    be available already to be output */
               break;
          } /* end switch */
       } /* end if */
     } /* end for */
     
     
 /* output the station */
    hdr.nprops = station.nprops = nrequested;
    free((void *)hdr.prop_id);
    hdr.prop_id = (int *) calloc((size_t)nrequested, sizeof(int));
    j = 0;
    for (i = 0; i < nrequested; ++i) {
       index = prop_req[i];
       hdr.prop_id[j++] = index;
    }

    /* calculate the station ID */
    station_id = station_id + 1;
    
    nbytes = 0;
    st = &str[0];
    if (lopt) { 
       sprintf(st,"%9.3f %9.3f ", hdr.lon, hdr.lat);
       nbytes += 20;
    }
    
    st = &str[nbytes];
    if (yopt) {
       sprintf(st,"%4d ", hdr.year);
       nbytes += 5;
    }
   
    st = &str[nbytes];
    if (mopt) { 
       sprintf(st,"%2d ", hdr.month);
       nbytes += 3;
    }
    st = &str[nbytes];
    if (dopt) { 
       sprintf(st,"%2d ", hdr.day);
       nbytes += 3;
    }
   
    st = &str[nbytes];
    if (sopt) {
       sprintf(st,"%5d %4d ", hdr.cruise, hdr.station);
       nbytes += 11;
    }

    st = &str[nbytes];
    if (uopt) {
       sprintf(st,"%8d ", station_id);
       nbytes += 9;
    }

    len_str = strlen(str);   
    scan = (double *) malloc(hdr.nprops * sizeof(double));
    for (i = 0; i < hdr.nobs; ++i) {
       for (j = 0; j < hdr.nprops; ++j) {
          prop = *get_prop_array(&station, hdr.prop_id[j]);
	  if (prop != NULL)
             scan[j] = prop[i];
	  else
	     scan[j] = (double)HB_MISSING;
       }
       if (len_str > 0)
           fprintf(outfile,"%s", str);
       write_hydro_scan(outfile, scan, hdr.nprops, hdr.prop_id);
    }
    free(scan);
    free_hydro_data(&station);
        
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

   
