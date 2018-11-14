/*  hb_layerav.c

................................................................................
                          *******  HydroBase3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1996
			     updated to ANSI standard Dec 2001
			     updated to HydroBase 3 Mar 2011
			     updated for T90 variables Nov 2011
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_layerav filename(_roots) -1mne/value -2mne/value -P<proplist> [-Wwindow/incr]  [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument

 -1 : property and value at upper level;
 -2 : property and value at lower level;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
hb_layerav computes average value of properties over a pressure interval
defined by two surfaces. Each observation contributing to the average
is weighted by the pressure over which it was observed.  Outputs 
year, month, lat, lon, prop1 .. propN .                                                  ____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int *prop_req;         /* list of props requested for output */
int nrequested;                /* count of above */	
int prop_needed[MAXPROP];   /* set of props requested || defining a
                                    surface || needed for computation */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
double x1, x2;
double pref1, pref2;
int prop1, prop2;
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */
int warnflag;
int yflag, mflag, lflag, iflag, thickflag, zflag;
double zmin, zmax;
int check_bottom_outcrop;
int tindex;           /* specifies temperature variable for computed props*/      
   /* prototypes for locally defined functions */

void    print_usage(char *);
void    get_hydro_data(FILE *);
int     parse_prop_list(char *);
double get_weighted_av(double, double, double *);
void compute_prop(int, double **, double *, double *, double *, int, double, double);

int main (int argc, char **argv)
{
   short   popt, opt1, opt2;
   int     i, n;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   char    *id;
   FILE *infile, *outfile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    popt = opt1 = opt2 = 0;
    error = 0;
    window = 50;
    w_incr = 10;
    zmin = 0.0;
    zmax = 10000.0;
    infile = stdin;
    outfile = stdout;
    yflag = mflag = lflag = iflag = thickflag = zflag = 0;
    warnflag = 0;           /* set to 1 after warning is printed */
    check_bottom_outcrop = 0;

/* initialize these ... */


   prop_req = (int *)calloc(MAXPROP, sizeof(int));
   nrequested = 0;
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
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;

               case '1':
                        opt1 = 1;
			id = (char *) calloc(6,sizeof(char));
			st = &argv[i][2];
                        error = (sscanf(st,"%[^'/']", id) == 1) ? 0 : 1;
			prop1 = get_prop_indx(id);
                        if (prop1 < 0 || prop1 > MAXPROP ) {
                           fprintf(stderr,"\n%s is not an appropriate property.", id);
                           ++error;
                        }
			st += strlen(id);
			if (*st == '/')
                           ++st;
                        
                        if (prop1 == (int)S_ || prop1 == (int)HT || prop1 == (int)PE ) {
                         if (sscanf(st,"%lf/%lf", &pref1, &x1) != 2)
                            ++error;
                        }
                        else {
                         if (sscanf(st,"%lf", &x1) != 1)
                            ++error;
                        }
			
			free(id);
                        break;

               case '2':
                        opt2 = 1;
			id = (char *)calloc(6, sizeof(char));
			st = &argv[i][2];
                        error = (sscanf(st,"%2s", id) == 1) ? 0 : 1;
			prop2 = get_prop_indx(id);
                        if ( prop2 < 0  || prop2 > MAXPROP) {
                           fprintf(stderr,"\n%s is not an appropriate property.", id);
                           ++error;
                        }
			
			st += strlen(id);
			if (*st == '/')
                           ++st;
                        
                        if (prop2 == (int)S_ || prop2 == (int)HT || prop2 == (int)PE ) {
                         if (sscanf(st,"%lf/%lf", &pref2, &x2) != 2)
                            ++error;
                        }
                        else {
                         if (sscanf(st,"%lf", &x2) != 1)
                            ++error;
                        }
			n = strlen(argv[i]);
			if (argv[i][n-1] == 'b')   /* check for optional b */
			    check_bottom_outcrop = 1;
                        
			free(id);
			break;
               case 'H':
                        thickflag = 1;
                        break;
               case 'M':
                        mflag = 1;
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
               case 'O':
                        outfile = fopen(&argv[i][2],"w");
                        if (outfile == NULL) {
                           fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
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
			       ++error;
			       fprintf(stderr,"\nOption to output thickness is now -H");
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
                        
               case 'Z':
	                zflag = 1;
                        error = (sscanf(&argv[i][2],"%lf", &zmin) == 1) ? 0 : 1;
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                            if (*st == '/') {
                              ++st;
                              error = (sscanf(st,"%lf", &zmax) == 1) ? 0 : 1;
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

   if (!popt ) {
       fprintf(stderr,"\nYou must specify output properties with -P.\n");
       exit(1);
   }
   if (!opt1 || !opt2 ) {
       fprintf(stderr,"\nYou must define the layer with -1 and -2.\n");
       exit(1);
   }
   if (  !(mflag || yflag || lflag || iflag)) {
       fprintf(stderr,"\nWARNING: You have not specified any output parameter from {month|year|position|station_id}.  Only properties will be output.\n");
   }
   
   if (zflag) {
       fprintf(stderr,"\nUsing pressure limits %.2lf - %.2lf db\n", zmin, zmax);
   }
   if (thickflag) {
       fprintf(stderr,"\nLayer thickness (in meters) will be output\n");
   }
   
   if ((prop1 == (int)GN) || (prop2 == (int)GN)) {
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
   }

/* loop for each input file */

   do {

    if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... \n");
    }
    else {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile == NULL)
          goto NEXTFILE;
    }

            /* read each file completely */

    get_hydro_data(infile);
    
    if (nfiles) fclose(infile);


NEXTFILE:
     ;

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_layerav.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nhb_layerav computes average value of properties over a pressure interval");
   fprintf(stderr,"\ndefined by two surfaces. Each observation contributing to the average");
   fprintf(stderr,"\nis weighted by the pressure over which it was observed.");      fprintf(stderr,"\nThe layer can be further limited to a pressure range using -Z<zmin/zmax>");     fprintf(stderr,"\nOutputs year, month, lat, lon, prop1 .. propN. ");
   
    fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -1<prop_id[pref]/value> -2<prop_id[pref]/value[b]> -P<list_of_properties> [-D<dirname>] [-E<file_extent>] [-H] [-I] [-L] [-M] [-T68] [-Y] [-W<window>[/<w_incr>] [-Z<zmin/zmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or input is expected to come from stdin...");
   fprintf(stderr,"\n    -1  : surface 1 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n    -2  : surface 2 property/value (if property = s_");
   fprintf(stderr,"\n          you need to specify pref after s_");
   fprintf(stderr,"\n          Append b to average down to the bottom if this");
   fprintf(stderr,"\n          surface is deeper than the bottom.  This will");
   fprintf(stderr,"\n          only occur if deepest observation is within 100m of the seafloor.");
   fprintf(stderr,"\n    -P  : list of properties: ex: -Pth/sa");
   fprintf(stderr,"\n          by itself, -P will list available properties");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current");
   fprintf(stderr,"\n          directory)  ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-H] : include thickness of layer in output listing");
   fprintf(stderr,"\n   [-I] : include station ID in output listing");
   fprintf(stderr,"\n   [-L] : include lat/lon in output listing");
   fprintf(stderr,"\n   [-M] : include month in output listing");
   fprintf(stderr,"\n   [-O] : name of output file  (default is stdout).");
   fprintf(stderr,"\n   [-T] : use IPTS-68  for computing derived variables (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
   fprintf(stderr,"\n   [-Y] : include year in output listing");
   fprintf(stderr,"\n   [-W] : pressure window in db for computing gradient");
   fprintf(stderr,"\n          properties (bvfreq, pv...) Default values: [%1d/%1d] ", window, w_incr);
   fprintf(stderr,"\n   [-Z] : specify pressure limits for layer.");
   fprintf(stderr,"\n          Exclude pressures that fall outside of these limits from weighted average.");
   fprintf(stderr,"\n   [-h] : help... prints this message.");
   fprintf(stderr,"\n\n");  
       
   return;
}


/*****************************************************************************/
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
     if (index < MAXPROP) {  /* don't use count, variance, or quality arrays */
           prop_req[nprops] = index;
           prop_needed[index] = 1;
      }
	    
     if (index == (int)TH)
         prop_needed[(int)TE] = 1;    
     if (index == (int)TH9)
         prop_needed[(int)T90] = 1;    

      ++nprops;

      /* !**!  Special cases for properties */

      if (((enum property)index == S_ ) || ((enum property) index == PE) || ((enum property) index == HT) ) {
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %.2s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }

         if (s_pref >= 0) {   /* has this already been requested? */
                     if ( NINT( s_pref) != NINT(ref_val)) {
                         fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                         fprintf(stderr,"You can only use one of those prefs");
                         fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                        exit(1);
                     }
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
     prop_req[nprops] = (int)GN;
     prop_needed[(int)GN] = 1;
     ++nprops;
  }
  
  if (prop_needed[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

 /* end of special cases */

  return (nprops);
  
}  /* end parse_prop_list() */
/*****************************************************************************/
void get_hydro_data(FILE * file)
   /*  Reads each station in a HydroBase file and computes property values
       to find value at each surface.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, index, j, nreq;
   int main_props_avail;
   double  p1, p2, avg, pref;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check if pr, de, te, and sa are available ... */
  
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
     } /* end switch */
     
     if (prop_needed[(int)TH9]) {
         free_and_alloc(&station.observ[(int)TH9], hdr.nobs);
         compute_theta(hdr.nobs, station.observ[(int)TH9], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
     }
     if (prop_needed[(int)TH]) {
         free_and_alloc(&station.observ[(int)TH], hdr.nobs);
         compute_theta(hdr.nobs, station.observ[(int)TH], station.observ[(int)PR], station.observ[(int)TH],station.observ[(int)SA] );
     }
  
       main_props_avail = 1;
       if (!(available(PR, &hdr) && available(DE, &hdr) && available((enum property)tindex, &hdr)  && available(SA, &hdr))) 
         main_props_avail = 0;
       

/* get pressure associated with 1st surface...*/  
  
      if (station.observ[prop1] == NULL && main_props_avail ) {

	 switch ((enum property) prop1) {
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
               free_and_alloc(&station.observ[prop1], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[prop1], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;
             case TH9:
               free_and_alloc(&station.observ[prop1], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[prop1], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
               break;
	       
              default:
               free_and_alloc(&station.observ[prop1], hdr.nobs);	       
	       compute_prop(prop1, &station.observ[prop1], station.observ[(int)PR], station.observ[(int)tindex], station.observ[(int)SA], hdr.nobs, pref1, (double) hdr.lat); 
	       
         }  /* end switch */
      } /* end if !available prop1 */
      
      p1 = -9.0;
      if ( station.observ[prop1] != NULL )
         p1 = hb_linterp(x1, station.observ[prop1], station.observ[(int)PR],
           hdr.nobs);
           
      if (p1 < -8.) {  /* didn't find surf 1 */
           
        switch ((enum property) prop1) {
         case S0:
         case S1:   /* if surface density obs > density we are */
         case S2:   /* seeking, set p1 to first pressure in array */
         case S3:
         case S4:
         case S_:
         case GN:
             if ((station.observ[prop1][0] > x1) 
             && (station.observ[(int)PR][0] < 150.)) {
                p1 = station.observ[(int)PR][0];
             }
             break;
          default:
             ;
        }  /* end switch */
      
      } /* end if */
      
      if (zflag && (p1 >= 0.)) {   /* check depth limits */
          if (p1 < zmin)
	     p1 = zmin;
	  if (p1 > zmax)
	     p1 = -9999.;
      } 
           
/* get pressure associated with 2nd surface...*/    
 
      if ((station.observ[prop2] == NULL) && (main_props_avail)) {

	 switch ((enum property) prop2) {
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
               free_and_alloc(&station.observ[prop2], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[prop2], station.observ[(int)PR], station.observ[(int)TE],station.observ[(int)SA] );
               break;
             case TH9:
               free_and_alloc(&station.observ[prop2], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[prop2], station.observ[(int)PR], station.observ[(int)T90],station.observ[(int)SA] );
               break;
              default:
               free_and_alloc(&station.observ[prop2], hdr.nobs);	       
	       compute_prop(prop2, &station.observ[prop2], station.observ[(int)PR], station.observ[(int)tindex], station.observ[(int)SA], hdr.nobs, pref2, (double) hdr.lat); 
         }  /* end switch */

      } /* end if */
      
      p2 = -9.0;
      if (station.observ[prop2] != NULL)     
         p2 = hb_linterp(x2, station.observ[prop2], station.observ[(int)PR], hdr.nobs);

      if (zflag && (p2 >= 0.)) {   /* check depth limits */
          if (p2 < zmin)
	     p2 = -9999.;
	  if (p2 > zmax)
	     p2 = zmax;
      } 

      if (check_bottom_outcrop) {
         if ((p1 > -8) && (p2 < -8.)) {  
           
           switch ((enum property) prop2) {
            case S0:
            case S1:   
            case S2:  
            case S3:
            case S4:
            case S_:
	    case PR:
	    case DE:
                if ((x2 > station.observ[prop2][hdr.nobs -1] )
                && (hdr.pdr != 0) 
                && (station.observ[(int)DE][hdr.nobs -1] > (hdr.pdr - 100))) {
                   p2 = station.observ[(int)PR][hdr.nobs -1];
                }
                break;
             default:
                ;
           }  
         }
      }

      if ((p1 > -1) && (p2 > -1) && (p2 >= p1)) {
        if (yflag)
           fprintf(stdout,"%5d ", hdr.year);
        if (mflag)
           fprintf(stdout,"%2d ", hdr.month);
        if (lflag)
           fprintf(stdout,"%8.3f %8.3f ", hdr.lon, hdr.lat);
        if (iflag)
          fprintf(stdout,"%5d ", hdr.station);
     
        if (thickflag) {
	  fprintf(stdout, "%10.3lf ",  (hb_depth(p2,(double)hdr.lat) - hb_depth(p1,(double)hdr.lat)) );
	}
       
        for (j = 0; j < nrequested; ++j) {
	     index = prop_req[j];
             if (station.observ[index] == NULL) {
                switch ((enum property) index) {
                  case S_:
                      pref = s_pref;
                      break;
                  case PE:
                      pref = pe_pref;
                      break;
                  case HT:
                      pref = ht_pref;
                      break;
                  default:
                      pref = 0;
                }  /* end switch */
                
                free_and_alloc(&station.observ[index], hdr.nobs);
                compute_prop(index, &station.observ[index], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], hdr.nobs, pref, (double) hdr.lat);
             }
             avg = get_weighted_av(p1, p2, station.observ[index]);
             fprintf(stdout, " %10.4lf", avg);
        } /* end for */
        fprintf(stdout, "\n");
      } /* end if */
    
  /* clean up...*/  
          
      free_hydro_data(&station);
      
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

/****************************************************************************/
double get_weighted_av(double p1, double p2, double *xptr)

/*  computes a weighted property average over the pressure interval
    specified by p1 and p2.  Returns the average or -999. for no value.  */
{
   int i, n, start, end;
   double x1, x2;
   double weight, w, sum;
   double *xtmp, *ptmp;
  
   if (xptr == NULL) 
      return (-999.0);
   
   xtmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   ptmp = (double *) calloc ((size_t)hdr.nobs, sizeof(double));
   n = 0;
   for (i = 0; i < hdr.nobs; ++i) {
      if (xptr[i] >= 0.0) {
         xtmp[n] = xptr[i];
         ptmp[n] = station.observ[(int)PR][i];
         ++n;
      }
   }
   if (n == 0){
      free(ptmp);
      free(xtmp);
      return (-999.);
   }   
      
   x1 = hb_linterp(p1, ptmp, xtmp, n);
   x2 = hb_linterp(p2, ptmp, xtmp, n);
   
   if ((x1 < -8.0) || (x2 < -8.0)) {
      free(ptmp);
      free(xtmp);
      return (-999.);
   }

   /* find index of first pressure greater than pressure at top of layer */   
   start = 0;
   while (p1 > ptmp[start])
     ++start;
   /*  find largest index of pressure array less than pressure at bottom of layer */
   end = n-1;
   while ((p2 <= ptmp[end]))
     --end;
   
   /* now do weighted average... */
   
   w = 0.0;
   sum = 0.0;
   for (i = start; i <= end; ++i) {
         weight = ptmp[i] - p1;
         w += weight;
         sum += (xtmp[i] + x1) * 0.5 * weight;
         p1 = ptmp[i];
         x1 = xtmp[i];
   }
   
   /* add last pressure interval ... */
   
   weight = p2 - p1;
   w += weight;
   sum += (x1 + x2) *.5 * weight;
   
   free(ptmp);
   free(xtmp);
   
   return (sum / w);
   
}  /* end get_weighted_av() */
/****************************************************************************/
void compute_prop(int i, double **xaddr, double *pptr, double *tptr, double *sptr, int n, double ref_val, double dlat)

/*
     int i;         index to enum property 
     double **xaddr address of array already allocated 
     double *pptr   pressure array 
     double *tptr   temperature array 
     double *sptr   salinity array 
     int n          number of obs 
     double ref_val ref pressure for s_, ht, or pe 
     double dlat    latitude 
 */
{

int j;
double *xptr, *xxptr;

   xptr = *xaddr;

/* !**! Special cases for individual properties... */
   switch ((enum property) i) {
       case S0:
               compute_sigma(0., n, xptr, pptr, tptr, sptr);
               break;
       case S1:
               compute_sigma(1000., n, xptr, pptr, tptr, sptr );
               break;
       case S2:
               compute_sigma(2000., n, xptr, pptr, tptr, sptr);
               break;
       case S3:
               compute_sigma(3000., n, xptr, pptr, tptr, sptr);
               break;
       case S4:
               compute_sigma(4000., n, xptr, pptr, tptr, sptr);
               break;
       case S_:
               compute_sigma(ref_val, n, xptr, pptr, tptr, sptr);
               break;
       case HT:
               compute_height(n, pptr, tptr, sptr, ref_val, xptr);
               break;
       case PE:
               compute_energy(n, pptr, tptr, sptr, ref_val, xptr);
               break;
       case SV:
               compute_sp_vol( n, xptr, pptr, tptr, sptr);
               break;

        case VA:
               compute_svan( n, xptr, pptr, tptr, sptr);
               break;
	       
        case DR:
               for ( j= 0; j< n; ++j) {
	          xptr[j] = hb_ratio(sptr[j], tptr[j], pptr[j]);
	       }
               break;
        case AL:
               for (j= 0; j< n; ++j) {
	          xptr[j] = hb_alpha(sptr[j], tptr[j], pptr[j]);
	       }
               break;
        case BE:
               for (j= 0; j< n; ++j) {
	          xptr[j] = hb_beta(sptr[j], tptr[j], pptr[j]);
	       }
               break;
       case VS:
               compute_sound_vel( xptr, pptr, tptr, sptr, n);
               break;

       case PV: 
               buoy_freq(xptr, pptr, tptr, sptr, n, window, w_incr);
               for (j = 0; j < n; ++j) {
                 if (xptr[j] > -99.) 
                   xptr[j] *= xptr[j];
               }
               po_vort(xptr, xptr, n, dlat);

               break;
               
       case BF:
               buoy_freq(xptr, pptr, tptr, sptr, n, window, w_incr);
               for (j = 0; j < n; ++j) {
                    if (xptr[j] > -9998.)   /* convert the units */
                        xptr[j] *= 1.0e5;
               }
               break;
       case GN:
               compute_gamma_n(&ginfo, n, xptr, pptr, tptr, sptr, (double) hdr.lon, dlat);
	       break;
       case GE:
               xxptr = (double *) calloc((size_t)n, sizeof(double));
               compute_gamma_nerr(&ginfo, n, xxptr, xptr, pptr, tptr, sptr, (double) hdr.lon, dlat);
	       free((void *)xxptr);
	       break;
       default: 
               free(xptr);
               xptr = NULL;
               break;
  } /* end switch */
          
  return;

}   
