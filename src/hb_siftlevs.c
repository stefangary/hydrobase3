/*  hb_siftlevs.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     updated to HB3M Dec 2009
................................................................................
  Evaluates each level in a HydroBase station file and separates scans 
  which meet the specified criteria from those that don't.  It is a tool to search for errant data.
____________________________________________________________________________
  USAGE:  

 hb_siftlevs filename(_roots) -P<property/min/max> [-O<output_file>]  [-R<reject_file>]  [-Zzmin/max] [-Ttmin/tmax] [-Ssigmin/sigmax] [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument!

 -P : property to evaluate with min/max values to search for; 
        ex: -Ps2/40/100
      -P (by itself) will print a list of the available properties;
      
 -O : specifies output file (default is to stdout) 
 -R : specifies file for data that do NOT meet criteria 
 -S : specifies min/max sigma-0 range
 -T : specifies theta range
 -Z : specifies min/max pressure range
 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
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

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int propindex;               /* for -P option */
int prop_req[MAXPROP];       /* set of  props requested for output */
int nrequested;      
int prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station, bstation, gstation;
struct HYDRO_HDR hdr, bhdr, ghdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/
double xmin, xmax;
double tmin, tmax;    /* theta limits over which to search */
double zmin, zmax;    /* pressure limits over which to search */
double sigmin, sigmax;    /* sigma limits over which to search */
int count;
int r_flag, t_flag, s_flag;
FILE *rejectfile, *outfile;

struct GAMMA_NC ginfo;   /* used for neutral density */

  /* prototypes for locally defined functions */
  
void print_usage(char *);
int parse_prop_list(char *);
void get_hydro_data(FILE *);

main (int argc, char **argv)
{
   short   popt;
   int     i;
   int     curfile,  nfiles; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   FILE *  infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    popt = 0;
    r_flag = s_flag = t_flag = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    count = 0;
    zmin = -10;
    zmax = 10000;
    tmin = -10;
    tmax = 40;
    sigmin = -10;
    sigmax = 100;
    outfile = stdout;
    curfile = 1;
    nfiles = 0;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
       prop_req[i] = 0;
       prop_needed[i] = 0;
      station.observ[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
      bstation.observ[i] = (double *) NULL;
      bstation.variance[i] = (double *) NULL;
      bstation.count[i] = (double *)NULL;
      bstation.quality[i] = (double *)NULL;
      gstation.observ[i] = (double *) NULL;
      gstation.variance[i] = (double *) NULL;
      gstation.count[i] = (double *)NULL;
      gstation.quality[i] = (double *)NULL;
   }
   hdr.prop_id = (int *) NULL;
   bhdr.prop_id = (int *) NULL;
   ghdr.prop_id = (int *) NULL;

   prop_needed[(int)S0] = 1;
   prop_needed[(int)TH9] = 1;
   prop_needed[(int)PR] = 1;


/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			
			}
                        break;
               case 'R':
                        r_flag = 1;
                        rejectfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                        if (rejectfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.  It may already exist?\n", &argv[i][2]);
                          exit(1);
                        }
                        break;
               case 'P':
                        popt = 1;
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        error = parse_prop_list(&argv[i][2]);
                        break;

                case 'S':
                        error = (sscanf(&argv[i][2],"%lf/%lf", &sigmin, &sigmax) == 2) ? 0 : 1;
 			prop_req[nrequested] = (int)S0;
			++nrequested;
                       s_flag = 1;
                        break;
                      
                case 'T':
                        error = (sscanf(&argv[i][2],"%lf/%lf", &tmin, &tmax) == 2) ? 0 : 1;
                        t_flag = 1;
			prop_req[nrequested] = (int)TH9;
			++nrequested;
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
                        error = (sscanf(&argv[i][2],"%lf/%lf", &zmin, &zmax) == 2) ? 0 : 1;
 			prop_req[nrequested] = (int)PR;
			++nrequested;
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

   if (  ! (popt && r_flag)) {
       fprintf(stderr,"\nYou must specify a property with -P and a reject file with -R.\n");
       exit(1);
   }

  if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
  
/* loop for each input file */

   do {

     infile = stdin;
     if (nfiles > 0) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile == NULL)
          goto NEXTFILE;
     }

     
            /* read each file completely */

 
         get_hydro_data(infile);
         fprintf(stderr,"\nNumber of stations that fit criteria: %d\n", count);


NEXTFILE:
        if (nfiles > 0 && infile != NULL) 
            fclose(infile);
  
   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_siftlevs.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s determines whether an observation level has a value of", program);
   fprintf(stderr,"\na property which falls between the specified min/max");
   fprintf(stderr,"\nranges.  If so, it prints the scans to stdout in HydroBase station format.  Scans that do not meet the criteria are written to the file specified by -R\n");
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<property/min/max>  [-O<outfile>] [-R<reject_file>] [-W<window>[/<w_incr>]][-Z<zmin>/<zmax>] [-S<sigmin>/<sigmax>] [-T<tmin>/<tmax>] [-D<dirname>] [-E<file_extent>] [-h]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n-P : property mnemonic/min/max;");
   fprintf(stderr,"\n     ex:  -Psa/0/30");
   fprintf(stderr,"\n          -P (by itself) produces a list of available    properties");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n     ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n     ex: -E.dat ");
   fprintf(stderr,"\n-O : output filename, default is to stdout");  
   fprintf(stderr,"\n-R : filename for rejected stations");
   fprintf(stderr,"\n-S : Optional sigma-0 limits min/max to search;");
   fprintf(stderr,"\n     ex:  -S27/27.8");
   fprintf(stderr,"\n-T : Optional theta limits min/max to search;");
   fprintf(stderr,"\n     ex:  -T3/4");
   fprintf(stderr,"\n-Z : Optional pressure limits min/max to search;");
   fprintf(stderr,"\n     ex:  -Z0/1000");
   fprintf(stderr,"\n-W : Specifies pressure window (db) for computing ");
   fprintf(stderr,"\n     gradient properties (bf, pv)  ");
   fprintf(stderr,"\n     defaults: -W100/10");
   fprintf(stderr,"\n-h : help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)
    /*  Parses the  property mnemonic, min, max values and sets global 
        variables accordingly.
        Returns 0.  An error will cause an exit. */
{
   int  error, i;
   char prop[7];
   double ref_val;
  
      if (*st == '/')
         ++st;
       
      for(i = 0; i < 7; i++)
	    prop[i] = '\0';
      sscanf(st,"%[^'/']", prop);

      propindex = get_prop_indx(prop);
      if (propindex < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      if (propindex >= MAXPROP) {
         fprintf(stderr,"\n Not an appropriate property: %s\n", prop);
         exit(1);
      }
      prop_req[nrequested] = propindex;
      ++nrequested;
      prop_needed[propindex] = 1;


      /* !**!  Special cases for properties */

      if (((enum property)propindex == S_ ) || ((enum property) propindex == PE) || ((enum property) propindex == HT)) {
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %s", prop);
             fprintf(stderr,"\n   ex: -P%s1500/34.35/34.5\n", prop);
             exit(1);
         }

         switch ((enum property) propindex) {
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
   if (*st == '/')
     ++st;
   error = (sscanf(st,"%lf/%lf", &xmin, &xmax) == 2) ? 0 : 1;
   if (error) {
      fprintf(stderr,"\nUnable to parse min/max values in -P option: %s\n", st);
      exit(1);
   }

  return (0);
}  /* end parse_prop_list() */
 

/*****************************************************************************/
void get_hydro_data(FILE * file)
   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i, j,found, found2, gotone, pts_avail;
   int nprops_out;
   int *new_propid, *props_out;
   double dlat, dlon;
   double **propaddr;
   double *prop, *prop2;


/* read each station in file ... */

    count = 0;
    while ((error = get_station(file, &hdr, &station)) == 0) {

     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
	 
     if (!available((int)T90, &hdr) && available(TE,&hdr) ) 
         create_t90_arrays( &hdr, &station);

     if (prop_needed[(int)TE] && !available(TE, &hdr) && available((int)T90, &hdr))
         create_t68_arrays( &hdr, &station);
      
  /* check if basic properties pr, t90, and sa are available ... */

    pts_avail = (available(PR, &hdr) && available(T90, &hdr) 
                 && available(SA, &hdr));


/* save original prop_id  expand hdr.prop_id to hold new properties */

      nprops_out = hdr.nprops;
      props_out = (int *) calloc((size_t)MAXPROP*2, sizeof(int));
      new_propid = (int *) calloc((size_t)MAXPROP*2, sizeof(int));
      for (i = 0; i < hdr.nprops; ++i)  {
         new_propid[i] = hdr.prop_id[i];
	 props_out[i] = hdr.prop_id[i];
      }
      free(hdr.prop_id);
      hdr.prop_id = new_propid;
       
 
 /* compute appropriate properties at each level in station ... */

/* !**! Special cases for individual properties... */

    for (i = 0; i < MAXPROP; ++i) {
       if ( prop_needed[i] && pts_avail && !available((enum property)i, &hdr)) {
          switch ((enum property) i) {
             case OX:
               if (available(O2, &hdr)) {
	          free_and_alloc(&station.observ[i], hdr.nobs);
	          for (j = 0; j < hdr.nobs; ++j)
		     station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j],station.observ[(int)T90][j],station.observ[(int)SA][j]); 
	          hdr.prop_id[hdr.nprops] = i;
	          station.nprops = ++hdr.nprops;
	       
	       }
               break;
             case O2:
               if (available(OX, &hdr)) {
	          free_and_alloc(&station.observ[i], hdr.nobs);
	          for (j = 0; j < hdr.nobs; ++j)
		     station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j],station.observ[(int)T90][j],station.observ[(int)SA][j]); 
	           hdr.prop_id[hdr.nprops] = i;
	           station.nprops = ++hdr.nprops;
                    
		}    
               break;
             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case TH9:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
	       
             case S0:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(0., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S1:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(1000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S2:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(2000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S3:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(3000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S4:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(4000., hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case S_:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sigma(s_pref, hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case HT:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA], ht_pref, station.observ[i]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA], pe_pref, station.observ[i]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA]);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case DR:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_ratio(station.observ[(int)PR][j], station.observ[(int)T90][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case AL:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_alpha(station.observ[(int)PR][j], station.observ[(int)T90][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case BE:
                   free_and_alloc(&station.observ[i], hdr.nobs);
		   for (j = 0; j < hdr.nobs; ++j ) {
                      station.observ[i][j] = hb_beta(station.observ[(int)PR][j], station.observ[(int)T90][j],station.observ[(int)SA][j]);
		   }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;
             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)T90], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.) 
                   station.observ[i][j] *= station.observ[i][j];
               }
               po_vort(station.observ[i], station.observ[i], hdr.nobs, dlat);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;

               break;
               
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[(int)TE], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -9998.)   /* convert the units */
                        station.observ[i][j] *= 1.0e5;
               }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
               break;

             case GN:
	       if (!prop_req[(int)GE]) {
                  free_and_alloc(&station.observ[i], hdr.nobs);
                  dlat = (double) hdr.lat;
                  dlon = (double) hdr.lon;
	          compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[(int)T90],station.observ[(int)SA], dlon, dlat);
	       }
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
	       break;
	       
             case GE:
               free_and_alloc(&station.observ[(int)GE], hdr.nobs);
               free_and_alloc(&station.observ[(int)GN], hdr.nobs);
               dlat = (double) hdr.lat;
               dlon = (double) hdr.lon;
	       compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[(int)T90], station.observ[(int)SA], dlon,
		     dlat);
		    hdr.prop_id[hdr.nprops] = i;
		    ++hdr.nprops;
		    hdr.prop_id[hdr.nprops] = (int)GN;
		    ++hdr.nprops;
	       break;
             default:
               break;
          } /* end switch */
       } /* end if */
     } /* end for */
     
    /* add requested properties to original list of props */
    
     for (i = 0; i < nrequested; ++i) {
       if (available((enum property) prop_req[i],&hdr)){
          found2 = 0;
	  j = -1;
	  while (++j < nprops_out && !found2) {
	    if (prop_req[i] == props_out[j])
	       found2 = 1;
	  }
	  
	  if (!found2) {
	      props_out[nprops_out++] = prop_req[i];
	  }
       }
     }
     gstation.nprops = bstation.nprops = nprops_out;

 /*  allocate memory for good and bad scans */
    for (i = 0; i < nprops_out; ++i) {
        propaddr = get_prop_array(&station, props_out[i]);
	
       if (*propaddr != NULL ) {
         free_and_alloc(get_prop_array(&bstation, props_out[i]), hdr.nobs);
         free_and_alloc(get_prop_array(&gstation, props_out[i]), hdr.nobs);
       }
    }
    
    ghdr.nobs = bhdr.nobs = 0;
        
 /* check for values in min/max range */
    found = 0;
    gotone = 0;
    for (i = 0; i < hdr.nobs; ++i) {
       if (available( (enum property) propindex, &hdr)) { 
          if ((station.observ[(int)PR][i] >= zmin) && (station.observ[(int)PR][i] <= zmax )) {
            if ((station.observ[(int)TH9][i] >= tmin) && (station.observ[(int)TH9][i] <= tmax )) {
              if ((station.observ[(int)S0][i] >= sigmin) && (station.observ[(int)S0][i] <= sigmax )) {
                if ((station.observ[propindex][i] >= xmin) && (station.observ[propindex][i] <= xmax)) {
	      
                   found = 1;
	           gotone = 1;
	      
	         /* store the scan in gstation */
	           for (j = 0; j < nprops_out; ++j) {
		     propaddr = get_prop_array(&gstation, props_out[j]);
	             if (*propaddr != NULL) {
		       prop2 = *propaddr;
		       prop = *get_prop_array(&station, props_out[j]);
		       prop2[ghdr.nobs] = prop[i];
		     }
	           }
	      
	           ++ghdr.nobs;
	        }
	      }
            }
          }
       }
       if ( ! found) {    /* store the scan in bstation */
       
	 for (j = 0; j < nprops_out; ++j) {
	     propaddr = get_prop_array(&station, props_out[j]);
	     if (*propaddr != NULL) {
		   prop = *propaddr;
		   prop2 = *get_prop_array(&bstation,props_out[j] );
		   prop2[bhdr.nobs] = prop[i];
	     }
	 }
         ++bhdr.nobs;
       }
       found = 0;
    } /* end for i */
    

   if (gotone)  
       ++count;   
    
       /* output the station */
       
     free(hdr.prop_id);
     hdr.prop_id = (int *)calloc(nprops_out, sizeof(int));
    if (ghdr.prop_id != NULL)
	      free(ghdr.prop_id);
    if (bhdr.prop_id != NULL)
              free(bhdr.prop_id);
    bhdr.prop_id = (int *) calloc(nprops_out, sizeof(int));
    ghdr.prop_id = (int *) calloc(nprops_out, sizeof(int));
	      
    for (i = 0; i < nprops_out; ++i) {
          hdr.prop_id[i] = props_out[i];
          bhdr.prop_id[i] = hdr.prop_id[i];
          ghdr.prop_id[i] = hdr.prop_id[i];
    }
	  
    hdr.nprops = station.nprops = nprops_out;
    bhdr.nprops = bstation.nprops = nprops_out;
    ghdr.nprops = gstation.nprops = nprops_out;
    
    strncpy(bhdr.country, hdr.country, 3);
    strncpy(ghdr.country, hdr.country, 3);   
    strncpy(bhdr.ship, hdr.ship, 3);
    strncpy(ghdr.ship, hdr.ship, 3); 
    strncpy(bhdr.qual, hdr.qual, NQUAL);
    strncpy(ghdr.qual, hdr.qual, NQUAL);   
    bhdr.origin = ghdr.origin = hdr.origin;  
    bhdr.instrument = ghdr.instrument = hdr.instrument;  
    bhdr.cruise = ghdr.cruise = hdr.cruise;  
    bhdr.station = ghdr.station = hdr.station;  
    bhdr.year = ghdr.year = hdr.year;  
    bhdr.month = ghdr.month = hdr.month;  
    bhdr.day = ghdr.day = hdr.day;  
    bhdr.lat = ghdr.lat = hdr.lat;  
    bhdr.lon = ghdr.lon = hdr.lon;  
    bhdr.pdr = ghdr.pdr = hdr.pdr;  
    bhdr.ms10 = ghdr.ms10 = hdr.ms10;  
    bhdr.ms1 = ghdr.ms1 = hdr.ms1;  
    
    bstation.nobs = bhdr.nobs;
    gstation.nobs = ghdr.nobs;
    if ((ghdr.nobs + bhdr.nobs) != hdr.nobs) {
         fprintf(stderr,"\nNumber of bad (%d) + good (%d) obs does not equal original number of obs (%d)\n", ghdr.nobs, bhdr.nobs, hdr.nobs);
	 exit(1);
    }
       
       /* output scans that meet criteria... */
    if (ghdr.nobs > 0)
       write_hydro_station(outfile, &ghdr, &gstation);
       
       /* output scans that don't meet criteria... */
       
    if (bhdr.nobs > 0 && r_flag) 
        write_hydro_station(rejectfile, &bhdr, &bstation);
      
      /* cleanup memory ... */
      
      free_hydro_data(&station);
      free_hydro_data(&bstation);
      free_hydro_data(&gstation);

   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

