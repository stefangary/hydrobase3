/*  hb_propcalc.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C 1999
			     updated to HydroBase 3: Dec 2008
			     updated for T90 May 2011
................................................................................
____________________________________________________________________________
  USAGE:  

 propcalc filename(_roots) -P<property_list> [-D<dirname>] [-E<file_extent>] [-W<window/w_incr>]

  list of filenames MUST be first argument!

 -P : list of properties to evaluate; ex: -Ppr/te/th/sa/ox
      -P (by itself) will print a list of the available properties;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 -W : specifies window width/increment for gradient properties (default is -W100/10)

____________________________________________________________________________
propcalc computes properties at each standard level in a station and 
outputs a (hydrobase format) ascii station file containing all the properties
specified with the -P option.                                                    ____________________________________________________________________________
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
int tindex;          /* index to temperature property (T90 vs. TE) */	

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
FILE *outfile;

struct GAMMA_NC ginfo;   /* used for neutral density */

main (int argc, char **argv)
{
   short   popt;
   int     index;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   int n_extents, curext; 
   char     *dir, *st;
   char  **extent_list;
   FILE *infile;
   void    print_usage(char *);
   void    get_hydro_data(FILE *);
   int     parse_prop_list(char *);


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
    
    infile = stdin;
    outfile = stdout;
    popt = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    tindex = (int)T90;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
   }
   hdr.prop_id = (int *) NULL;
   nrequested = 0;


/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                  /* allows multiple file extents */
                        extent_list[n_extents] = &argv[i][2];
			if (++n_extents == MAXEXT) {
                             fprintf(stderr,"\nToo many different file extents. Max allowed is %d\n", MAXEXT);
			 exit(1);
			}
                        break;

               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
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
			       error = 1;
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
       fprintf(stderr,"\nYou must specify input file(s) and properties.\n");
       exit(1);
   }
   
   if (n_extents == 0) ++n_extents;
   
   if (tindex == (int) T90)
      fprintf(stderr,"Using T90 temperatures\n");
   else
      fprintf(stderr,"Using T68 temperatures\n");


/* loop for each input file */

   do {
   
     curext = 0;
     do {
       if (nfiles) {
          infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_msg);
          if (infile == NULL)
             goto NEXTFILE;
        }
        else 
          fprintf(stderr,"\nStation input expected from stdin...\n");
     
            /* read each file completely */

 
         get_hydro_data(infile);

        if (nfiles)
         fclose(infile);
	 
NEXTFILE:
       ;

     } while (++curext < n_extents);
   } while (curfile++ < nfiles );

   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n\n***************************************************");
   fprintf(stderr,"\n hb_propcalc computes properties at each pressure level");
   fprintf(stderr,"\n in a station and outputs a HydroBase format"); 
   fprintf(stderr,"\n station file containing all the properties specified");
   fprintf(stderr,"\n with the -P option (where they exist or can be computed.");
   fprintf(stderr,"\n***************************************************");

                                                       
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<list_of_properties>   [-W<window>[/<w_incr>]] [-D<dirname>] [-E<file_extent>] [-T<68|90>]", program);

   fprintf(stderr,"\n\n  List of filename (roots) must be first arguments");
   fprintf(stderr,"\n-P : properties to list out;");
   fprintf(stderr,"\n          ex:  -Ppr/th9/sa/o2/ht");
   fprintf(stderr,"\n               -P (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : specifies directory input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-O : output filename, default is output to stdout");
   fprintf(stderr,"\n-T : Set temperature scale used in computing derived variables to either IPTS-68 (-T68) or ITS-90 (-T90) (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
   fprintf(stderr,"\n-W : Specifies vertical window length(db) and subdivision interval");
   fprintf(stderr,"\n           for computing gradient properties (bf and pv...)  ");
   fprintf(stderr,"\n          defaults: -W100/10 ");
                        

   fprintf(stderr,"\n-h help...... prints this message.");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)

   /*  Parses a list of property mnemonics and sets global flags accordingly.
       Returns the number of properties.  An error will cause an exit.
        
       char *st:    the list with -P stripped off  */
{
   int index, nprops, i, gamma_req;
   char prop[7];
   double ref_val;

   nprops = 0;
   gamma_req = 0;

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
      ++nprops;


      /* !**!  Special cases for properties */
      

      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == PE)) {
         
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }

        switch ((enum property) index) {
           case S_:
              if (s_pref >= 0) {   /* check for previous request? */
                  fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                  fprintf(stderr,"You can only use one of those prefs");
                  fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                  exit(1);
              }
              s_pref = ref_val;
              break;
           case PE:
              pe_pref = ref_val;
              break;
           case HT:
              ht_pref = ref_val;
              break;
           case GN:
              gamma_req = 1;
              break;
           default:
              ;
        }  /* end switch */
      }
      
   /* end of special cases */
      
      st += strlen(prop);


  } while (*st == '/');
  
 /* !**!  Special cases for properties */
  
  if (index == (int)GE && !gamma_req) {
     prop_req[nprops] = (int)GN;
     ++nprops;
     gamma_req = 1;
  }
  
  if (gamma_req) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
 /* end of special cases */

  return (nprops);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(FILE * file)
   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j, pts_avail, ratio_done, navail, index;
   int *tmp;
   double dlat, dlon;
   double **propaddr;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check if basic properties pr, te, and sa are available ... */
  
     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
     
     if (!available(T90, &hdr) && available(TE,&hdr) )  
         create_t90_arrays( &hdr, &station);

     if (!available(TE, &hdr) && available(T90,&hdr)) 
        create_t68_arrays( &hdr, &station);
      
      pts_avail = (available(PR, &hdr) && available((enum property)tindex,&hdr) 
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
     
 /* determine availability of requested properties... */
    navail = 0;
    for (i = 0; i < nrequested; ++i) {
       index = prop_req[i];
       if ( *(get_prop_array(&station,index)) != NULL )
          ++navail;
    }
    
 /* output the station */
    hdr.nprops = station.nprops = navail;
    free((void *)hdr.prop_id);
    hdr.prop_id = (int *) calloc((size_t)navail, sizeof(int));
    j = 0;
    for (i = 0; i < nrequested; ++i) {
       index = prop_req[i];
       propaddr = get_prop_array(&station, index);
       if (*propaddr != NULL)  {
         hdr.prop_id[j++] = index;
       }
    }
    write_hydro_station(outfile, &hdr, &station);
    free_hydro_data(&station);
    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

