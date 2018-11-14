/*  hb_getbottom.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             original 1995
			     Update to ANSI standard Feb 2000
			     update to HB3 Nov 2011
			     
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_getbottom filename(_roots) -X<prop> -Y<prop> -Z<prop> [-K] [-I] [-M] -S<sta_dist_file> [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument or input is expected to come 
  from stdin.

 -X : x-property: lat, lon or distance;
 -Y : y-property: 2-char HydroBase mnemonic;
 -Z : z-property: 2-char HydroBase mnemonic;
 
 -S : use file containing sta-dist-depth for x-positions of casts.
 -I : output station_id in first column
 -K : compute distance in km (nm is default)
 -M : add info to use output as a masking file
 -D : specifies directory for input files (default is current dir)  
      ex: -D/d1/hbase/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

____________________________________________________________________________
hb_getbottom finds the deepest observation level of Zproperty at each station and outputs an ascii listing of X,Y values.                                                    
____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

#define  NMperKM  .539593  /* nautical miles per km conversion */
#define    EXTENT   ""
#define    DIR      ""
#define   RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  3437.746873       /* in nm */

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
short prop_req[MAXPROP];         /* set of props requested for output */
short prop_needed[MAXPROP];      /* set of props requested || defining a
                                    surface || needed for computation */
int tindex;   /* specifies IPTS-68 or ITS-90 */	
int mflag;     /* add info to use output as a masking file */			 int firstline;
double lastx; 
FILE *outfile;  

int   xopt, yopt, zopt;
double *distance;
int   sopt, min, kilo, output_sta;
float prev_lat, prev_lon;
double cumdist;
int nstations;     /* flag for get_xprop() */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;          /* ref lev for computing a non-standard sigma level */
double pe_ref, ht_ref;   /* ref lev for pot energy and dyn height */
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */
       
/* prototypes for locally defined functions */

void    print_usage(char *);
void    get_hydro_data(FILE *);
int     parse_prop_list(char *);
   
   
main (int argc, char **argv)
{
   int     index, nsta;
   int     i, max, sta;
   int     curfile = 1, nfiles = 0;
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   float   depth;
   FILE    *stafile, *infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    xopt = yopt = zopt = -1;
    sopt = kilo = 0;
    output_sta = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    infile = stdin;
    outfile = stdout;
    stafile = NULL;
    cumdist = 0;
    mflag = 0;
    firstline = 1;

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
               case 'E':                    /* get file extent */
                        extent = &argv[i][2];
                        break;
               case 'I':
                        output_sta = 1;
                        break;
               case 'K':
                        kilo = 1;
                        break;
               case 'M':
                        mflag = 1;
                        break;
               case 'O':
                       outfile = fopen(&argv[i][2],"w");
                       if (outfile == NULL) {
                          fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                          exit(1);
                       }
                       break;
               case 'P':                    /* set starting position */
                        error = (sscanf(&argv[i][2],"%f/%f", &prev_lat, &prev_lon) != 2);
			nstations = 1;
                        break;
               case 'S':            /* get station/distance file */
                        sopt = 1;
                        stafile = fopen(&argv[i][2], "r");
                        if (stafile == NULL) {
                           fprintf(stderr,"\nUnable to open %s",  &argv[i][2]);
                           exit(1);
                        }
                        fprintf(stderr,"\nOpened %s",  &argv[i][2]);
                        break;
               case 'Z':
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        zopt = parse_prop_list(&argv[i][2]);
                        if (zopt < 0)
                             error = 1;
			
			     
                        break;
               case 'Y':
                        if ( argv[i][2] == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        yopt = parse_prop_list(&argv[i][2]);
                        if (yopt < 0)
                             error = 1;
			     
			if (yopt == (int)GE) {
                          fprintf(stderr,"-Yge (gamma error) is not a valid Y-axis property.");
                          fprintf(stderr,"Use -Ygn for neutral density");
			  exit(0);
			}
                        break;

               case 'X':
                        if ( argv[i][2] == '\0') {
                               fprintf(stderr,"\nX-Axis options:\n");
                               fprintf(stderr,"   la : latitude");
                               fprintf(stderr,"   lo : longitude");
                               fprintf(stderr,"   di : distance\n");
                               fprintf(stderr,"   yr : year\n");
                             exit(0);
                        }
                        st = &argv[i][2];
                        switch (*st) {
                           case 'l':
                              switch (*(++st)) {
                                 case 'a':
                                    xopt = 1;
                                    break;
                                 case 'o':
                                    xopt = 2;
                                    break;
                                 default:
                                    error = 1;
                              } /* end switch */
                              break;
                           case 'd':
                              switch (*(++st)) {
                                 case 'i':
                                    xopt = 3;
                                    break;
                                 default:
                                    error = 1;
                              } /* end switch */
                              break;
                            case 'y':
                              switch (*(++st)) {
                                 case 'r':
                                    xopt = 4;
                                    break;
                                 default:
                                    error = 1;
                              } /* end switch */
                              break;
                          default:
                              error = 1;
                           
                        } /* end switch */
                        if (xopt <= 0)
                             error = 1;
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

   if ( (xopt < 0) || (yopt < 0) || (zopt < 0)) {
       fprintf(stderr,"\nYou must specify X Y and Z properties.\n");
       exit(1);
   }
   
   if (prop_req[(int)GE] || prop_req[(int)GN]) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);

/* Read in station distance-depth file:
   first, count the stations to define the required memory space */
   
   if (sopt) {
      min = 9999999;
      max = 0;
      while (fscanf(stafile,"%d%*f%*f", &sta) == 1) {
         if (sta < min)
            min = sta;
         if (sta > max)
            max = sta;
      }
      nsta = max - min + 1;
      distance = malloc(nsta * sizeof(double));
      if (distance == NULL) {
       fprintf(stderr,"\nError allocating memory for distance array of size %d \n", nsta);
       exit(1);
      }
      
   /* now read the distance into the array.... */
   
      rewind(stafile);
      while (fscanf(stafile,"%d", &index) == 1) {
        index -= min;
        if (index >= nsta) {
         fprintf(stderr,"\nWe have a problem here...");
         fprintf(stderr,"\n    didn't count station index correctly?");
         exit(3);
        }
        fscanf(stafile,"%lf%f", &distance[index], &depth);
      }
      fclose(stafile);
   } /* end if */

/* loop for each input file */

   do {
   
       if (nfiles) {
          infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
          if (infile == NULL)
             goto NEXTFILE;
       }
 
       get_hydro_data(infile);

NEXTFILE:
         fclose(infile);

   } while (curfile++ < nfiles );
   
   if (mflag) 
   fprintf(outfile,"%lf 10000 \n", lastx);
   

   fprintf(stderr,"\n\nEnd of hb_getbottom.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n%s finds the deepest observation level of Zproperty", program); 
  fprintf(stderr,"\n at each station and outputs an ascii listing of X,Y values");                                    

   fprintf(stderr,"\n\nUsage:  %s filename_root(s)  -X<prop> -Y<prop> -Z<prop> [-K] [-I] [-S<distance_file>] [-D<dirname>] [-E<file_extent>] [-h] ", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or input is expected to come from stdin");
   fprintf(stderr,"\n-X  : x-property: [la]t,[lo]n, [di]stance or [yr]year;");
   fprintf(stderr,"\n      ex:  -Xlo");
   fprintf(stderr,"\n           -X (by itself) produces a list of available properties");
   fprintf(stderr,"\n-Y  : y-property:  2char mnemonic");
   fprintf(stderr,"\n      ex:  -Ypr");
   fprintf(stderr,"\n           -Y (by itself) produces a list of available properties");
   fprintf(stderr,"\n-Z  : Z-property:  2char mnemonic");
   fprintf(stderr,"\n      ex:  -Zpr");
   fprintf(stderr,"\n           -Z (by itself) produces a list of available properties");
   fprintf(stderr,"\n-K :  distances are in km (default is nm).");
   fprintf(stderr,"\n-I :  output the station_id in first column.");
   fprintf(stderr,"\n-M :  create a masking file for use with hb_gridsection.");
   fprintf(stderr,"\n-S :  name of file containing station/distance/depth for ");
   fprintf(stderr,"\n      assigning distance along track each cast.  If this option is not");
   fprintf(stderr,"\n      used, distance is computed from lat/lon.  ex: -Skn104.ctd.depth ");
   fprintf(stderr,"\n-D :  dirname for input files (default is current directory) ");
   fprintf(stderr,"\n      ex: -D../data/ ");
   fprintf(stderr,"\n-E :  input_file extent (default is no extent)");  
   fprintf(stderr,"\n      ex: -E.dat ");
   fprintf(stderr,"\n-h :  help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)

    /*  st : the list with -Y stripped off */
   
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index, i;
   char prop[7];
   double ref_val;

      for(i = 0; i < 7; i++)
	    prop[i] = '\0';

      if (*st == '/')
         ++st;
      sscanf(st,"%[^'/']", prop);
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      if (index > MAXPROP) {       
         fprintf(stderr,"\n Count and error properties are not valid options %s\n", prop);
         exit(1);
         
      }
      prop_req[index] = 1;
      prop_needed[index] = 1;
      if (index == (int)TH)
          prop_needed[(int)TE];

  /* !**! Special cases for individual properties... */

   /* s_ */
      if ((enum property)index == S_ ) {
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr,"\n Specify a ref pressure when requesting the property 's_'");
             fprintf(stderr,"\n   ex: -Ps_1500/th/sa\n");
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
          s_pref = ref_val;
      }

    /* ht */
      if ((enum property)index == HT ) {
         if (sscanf(st+2, "%lf", &ht_ref) != 1) {
             fprintf(stderr,"\n Specify a ref pressure when requesting the property 'ht'");
             fprintf(stderr,"\n   ex: -Pht1500/th/sa\n");
             exit(1);
         }
 
      }
      
     /* pe */
      if ((enum property)index == PE ) {
         if (sscanf(st+2, "%lf", &pe_ref) != 1) {
             fprintf(stderr,"\n Specify a ref pressure when requesting the property 'pe'");
             fprintf(stderr,"\n   ex: -Ppe2000/th/sa\n");
             exit(1);
         }

      }
      
  
	
   /* end of special cases */

      st += strlen(prop);

  return (index);
}  /* end parse_prop_list() */
 
/*****************************************************************************/
void get_prop(int i, int main_props_avail)

   /* arguments:
      i:         index corresponding to property 
      main_props_avail:  boolean: T if pr, te, sa are available 

      /* computes, if necessary and able, a property at each level in station ... */
{
   int j;
   double dlat, dlon;

/* !**! Special cases for individual properties... */

       if (prop_req[i] && !available((enum property)i, &hdr)) {
       
          if (!main_props_avail) {
              free_and_alloc(&station.observ[i], hdr.nobs);
	      for (j = 0; j < hdr.nobs; ++j) 
	        station.observ[i][j] = (double) HB_MISSING;
	      return;
	  }
	  
	  
          switch ((enum property) i) {
             case OX:
                 if (available(O2, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_kg2l(station.observ[(int)O2][j], station.observ[(int)PR][j], station.observ[tindex][j], station.observ[(int)SA][j]);
                    }
                 }
               break;
             case O2:
                 if (available(OX, &hdr)) {
                    free_and_alloc(&station.observ[i], hdr.nobs);
                    for (j=0; j < hdr.nobs; ++j) {
                       station.observ[i][j] = ox_l2kg(station.observ[(int)OX][j], station.observ[(int)PR][j], station.observ[tindex][j], station.observ[(int)SA][j]);
                    }
                 }
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
               compute_height(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], ht_ref, station.observ[i]);
               break;
             case PE:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_energy(hdr.nobs, station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], pe_ref, station.observ[i]);
               break;
             case SV:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sp_vol( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;

             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
               break;

             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.)  
                    station.observ[i][j] *= station.observ[i][j];
               }
               po_vort(station.observ[i],station.observ[i], hdr.nobs, dlat);

               break;
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -99.)   /* convert the units */
                        station.observ[i][j] *= 1.0e5;
               }
               break;
	       
	       
             case GN:
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               dlon = (double) hdr.lon;
	       compute_gamma_n(&ginfo, hdr.nobs, station.observ[i], station.observ[(int)PR],station.observ[tindex],station.observ[(int)SA], dlon, dlat);
	       break;
	       
             case GE:
               free_and_alloc(&station.observ[(int)GE], hdr.nobs);
               free_and_alloc(&station.observ[(int)GN], hdr.nobs);
               dlat = (double) hdr.lat;
               dlon = (double) hdr.lon;
	       compute_gamma_nerr(&ginfo, hdr.nobs, station.observ[(int)GN],   
		     station.observ[(int)GE], station.observ[(int)PR],
		     station.observ[tindex], station.observ[(int)SA], dlon,
		     dlat);
	       break;

             default:
               break;
          } /* end switch */
       } /* end if */
     
  return;
} /* end get_prop() */
/*****************************************************************************/
double get_x_prop(int option)
  /* option:  1 = lat, 2 = lon, 3 = distance, 4= year */
{
double dist, dx, dy;

   switch (option) {
      case 1:
         return ((double)hdr.lat);
      case 2:
         return ((double)hdr.lon);
      case 4:
         return ((double)hdr.year);
      case 3:
         if (sopt) {
           return(distance[hdr.station-min]);
         }
         
         if (nstations > 0) {
	   if (hdr.lon * prev_lon < 0 ) {  /* different signs, check if dateline is crossed */
	      if (ABS(hdr.lon - prev_lon) > 180) {
	           if (hdr.lon < 0)  /* make longitudes negative */
		       prev_lon -= 360.0;
		   else
		       prev_lon += 360.0;
	       }
	   }
           dy = (double) ABS(hdr.lat - prev_lat);
           dx = cos((double)(hdr.lat + prev_lat) *.5 * RADperDEG) * ABS(hdr.lon - prev_lon);
           dist = RADperDEG * EarthRadius * sqrt(dx * dx + dy * dy);  /* in nautical miles */
           if (kilo)
             dist /= NMperKM; 
           cumdist += dist;
         }
         prev_lat = hdr.lat;
         prev_lon = hdr.lon;
         ++nstations;
         return (cumdist);
      default:
        return (-1);
   } /* end switch */
}  /* end get_x_prop() */
/*****************************************************************************/
void get_hydro_data(FILE *file)

   /*  Reads each station in a HydroBase file and determines the deepest
       observation for the zproperty and outputs xprop,yrop pairs to outfile.  */
{
   int error, i, j, mainprops;
   double x, get_x_prop(int);
   void get_prop(int, int);


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {
    
     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
      
     if (!available((int)T90, &hdr) && available(TE,&hdr) ) 
         create_t90_arrays( &hdr, &station);

     if (prop_needed[(int)TE] && !available(TE, &hdr))
         create_t68_arrays( &hdr, &station);

  /* check that pr, te, and sa are available ... */
       mainprops = 1;
       if (!(available(PR, &hdr) && available((enum property)tindex, &hdr) && available(SA, &hdr))) {
         mainprops = 0;
       }
       
      get_prop(yopt, mainprops);
      get_prop(zopt, mainprops);
      x = get_x_prop(xopt);

    
 /* search for deepest level for the z property */
       if (station.observ[zopt] != NULL && station.observ[yopt] != NULL) {
          i = hdr.nobs-1;
          while ((i >= 0) && ((station.observ[zopt][i] < -8.) 
                 || (station.observ[yopt][i] < -8.))) {
             --i;
          }
    
          if (i >= 0) {
	  
	     if (mflag) {
	      lastx = x;
	      if (firstline) {
	        fprintf(outfile,">I\n %d 10000 \n", (int)x);
	        firstline = 0;
		}
	     }
	     
             if (output_sta)
               fprintf(outfile,"%5d ", hdr.station);
    
             if (xopt == 4) 
               fprintf(outfile,"%4d %10.4lf\n", (int)x, station.observ[yopt][i]);
             
             else 
               fprintf(outfile,"%10.3lf %10.4lf\n", x, station.observ[yopt][i]);
             
          }
       }    
       free_hydro_data(&station);
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

/****************************************************************************/
   
