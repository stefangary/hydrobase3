/*  hb_section.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1995
			     updated to ANSI Feb 2000
			     updated to HydroBase3 Jan 2010
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_section filename(_roots) -X<prop> -Y<prop> -Z<prop> [-K] [-P] [-W<window>[/<w_incr>] -S<sta_dist_file> [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument or input is expected to come 
  from stdin.

 -X : x-property: lat, lon, distance, or year;
 -Y : y-property: char HydroBase mnemonic;
 -Z : z-property: char HydroBase mnemonic;
 
 -S : use file containing sta-dist-depth for x-positions of casts.
 -K : compute distance in km (nm is default)
 -P:  specify lat/lon from which to compute distance to first station in file
 -D : specifies directory for input files (default is current dir)  
      ex: -D/d1/hbase/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum

 -W :  window(db) for gradient properties
____________________________________________________________________________
hb_section computes properties at each pressure level in a station and 
outputs an ascii listing of each observed level containing the properties
specified with the -X -Y and -Z options.                                                    
*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"

#define  NMperKM  .539593
#define  RADperDEG 0.017453292             /* pi/180 */
#define  EarthRadius  3437.746873       /* in nm */


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      
int gn_req;      /* gamma-n flag */
int tindex;      /* specifies temperature scale for derived variables  */

int   xopt, yopt, zopt;
int warning_given, lon0to360;
double *distance;
int   sopt, min, kilo, nstations;
float prev_lat, prev_lon;
double cumdist;
double FLAG;
FILE *outfile;

struct GAMMA_NC ginfo;   /* used for neutral density */

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double ht_pref;           /* ref lev for computing dynamic height */
double pe_pref;           /* ref lev for computing potential energy */
int window, w_incr;     /* used to specify pr window for gradient properties*/

   /* prototypes for locally defined functions  ... */
void    print_usage(char *);
void    get_hydro_data(FILE *);
int     parse_prop_list(char *);
double  get_x_prop(int);
double **get_prop(int, int);

main (int argc, char **argv)
{
   int     index, nsta;
   int     i, max, sta;
   int     curfile = 1, nfiles = 0;
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   FILE     *infile;
   float   depth;
   FILE    *stafile;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    FLAG = -8.9;
    xopt = yopt = zopt = -1;
    sopt = kilo = 0;
    error = 0;
    s_pref = -1;
    window = 100;
    w_incr = 10;
    infile = stdin;
    outfile = stdout;
    stafile = NULL;
    cumdist = 0;
    warning_given = 0;
    lon0to360 = 0;
    nstations = 0;
    gn_req = 0;
    tindex = (int) T90;
    
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *) NULL;
      station.quality[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;



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
               case 'F':
                        lon0to360 = 1;
                        break;
               case 'K':
                        kilo = 1;
                        break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr,"\nUnable to open %s",  &argv[i][2]);
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
                        break;

               case 'X':
                        if ( argv[i][2] == '\0') {
                               fprintf(stderr,"\nX-Axis options:\n");
                               fprintf(stderr,"   la : latitude");
                               fprintf(stderr,"   lo : longitude");
                               fprintf(stderr,"   di : distance\n");
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
                           default:
                              error = 1;
                           
                        } /* end switch */
                        if (xopt <= 0)
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

   if (gn_req) 
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
   
       if (nfiles > 0) {
          infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
          if (infile == NULL)
             goto NEXTFILE;
       }
 
       get_hydro_data(infile);

       if (nfiles > 0) {
         fclose(infile);
       }
NEXTFILE:
     ;

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_section.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nReads HydroBase station format profiles and outputs ");
   fprintf(stderr,"\nan ascii listing containing properties specified by");
   fprintf(stderr,"\n-X -Y and -Z options at each observed level.");

    fprintf(stderr,"\n\nUSAGE:  %s filename(_roots) -X<prop> -Y<prop> -Z<prop> [-D<dirname>] [-E<file_extent>] [-K] [-P] [-S<sta_dist_file>] [-T68|90] [-W<window>[/<w_incr>] ", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument or input is expected to come from stdin");
   fprintf(stderr,"\n-X : x-property: [la]t,[lo]n, or [di]stance;");
   fprintf(stderr,"\n          ex:  -Xlo");
   fprintf(stderr,"\n               -X (by itself) produces a list of available properties");
   fprintf(stderr,"\n-Y : y-property:  character property id");
   fprintf(stderr,"\n          ex:  -Ypr");
   fprintf(stderr,"\n               -Y (by itself) produces a list of available properties");
   fprintf(stderr,"\n-Z : Z-property:  character property id ");
   fprintf(stderr,"\n          ex:  -Zpr");
   fprintf(stderr,"\n               -Z (by itself) produces a list of available properties");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-F : force longitudes to be 0 - 360 for crossing the dateline (default is -180 to 180).");
   fprintf(stderr,"\n-K : distances are in km (default is nm).");
   fprintf(stderr,"\n-O : spedify name of output file. Default is stdout");
   fprintf(stderr,"\n-P : specify lat/lon to compute distance to first station in file (when X is distance).");

   fprintf(stderr,"\n-S : specifies name of file containing station/distance/depth for ");
   fprintf(stderr,"\n     assigning distance along track for each cast.  If this option is not");
   fprintf(stderr,"\n     used, distance is computed from lat/lon.  ex: -Skn104.ctd.depth ");
   fprintf(stderr,"\n-T : use IPTS-68  for computing derived variables (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
    fprintf(stderr,"\n-W : Specifies pressure window / pressure interval (in db)");
   fprintf(stderr,"\n     for computing vertical gradient properties (bf, pv)");
   fprintf(stderr,"\n     defaults: -W100/10 ");
   fprintf(stderr,"\n-h : help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)
    /*  Parses a list of property mnemonics and sets global flags accordingly.
        Returns the number of properties.  An error will cause an exit. */
{
   int index, i;
   char prop[6];
   double ref_val;

     for(i = 0; i < 6; i++)
	 prop[i] = '\0';

      if (*st == '/')
         ++st;
      
      sscanf(st,"%[^'/']", prop);
       
      index = get_prop_indx(prop);
      if (index < 0) {
         fprintf(stderr,"\n Unknown property requested: %s\n", prop);
         exit(1);
      }
      if (index == (int) GN || index == (int) GE)
            gn_req = 1;
      

  /* !**! Special cases for individual properties ... */

   /* properties requiring reference levels */
   
      if (((enum property)index == S_ ) || ((enum property) index == HT) || ((enum property) index == PE)) {
         if (sscanf(st+2, "%lf", &ref_val) != 1) {
             fprintf(stderr, "\n Specify a ref pressure for  %s", prop);
             fprintf(stderr,"\n   ex: -P%.2s1500/th/sa\n", prop);
             exit(1);
         }

        switch ((enum property) index) {
           case S_:
              if (s_pref >= 0) {   /* has this already been requested? */
                if ( NINT( s_pref) != NINT(ref_val)) {
                   fprintf(stderr,"Sorry. You have requested the property s_ with multiple reference levels.\n");
                   fprintf(stderr,"You can only use one of those prefs");
                   fprintf(stderr,"  You will probably have to do multiple runs for each different pref you want to associate with s_\n\n");
                   exit(1);
                }
               }
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


  return (index);
}  /* end parse_prop_list() */
 
/*****************************************************************************/
double ** get_prop(int i, int main_props_avail)

      /* computes, if necessary and able, a property at each level in station ... 
      
                   i:  index corresponding to property 
    main_props_avail:  > 0 (true) if pr, te, sa are available 
    returns the address of the array ptr.
    
    */
{
   int j;
   double dlat, dlon;
   

/* !**! Special cases for individual properties... */

       if ( !available((enum property)i, &hdr) && main_props_avail) {
          switch ((enum property) i) {
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
	     case TE:
                create_t68_arrays( &hdr, &station);
		break;
	     case T90:
                create_t90_arrays( &hdr, &station);
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
            case VS:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs);
               break;

             case DR:
	     case AL:
	     case BE:
                  free_and_alloc(&station.observ[(int)DR], hdr.nobs);
                  free_and_alloc(&station.observ[(int)AL], hdr.nobs);
                  free_and_alloc(&station.observ[(int)BE], hdr.nobs);
                  compute_ratio( hdr.nobs, station.observ[DR], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], station.observ[(int)AL], station.observ[(int)BE]);
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
	     
	        /* it must already exist */
		
               break;
          } /* end switch */
       } /* end if */
     
  return (get_prop_array(&station, i));
} /* end get_prop() */
/*****************************************************************************/
double get_x_prop(int option)

   /* option:  1 = lat, 2 = lon, 3 = distance */
{
double dist, dx, dy;

       if (lon0to360) {
           if (hdr.lon < 0)
              hdr.lon += 360.0;
       }

   switch (option) {
      case 1:
         return ((double)hdr.lat);
      case 2:
         return ((double)hdr.lon);
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
         if ((prev_lon * hdr.lon) < 0) {
            if ((ABS(hdr.lon - prev_lon) > 179.0) && !warning_given) {
               fprintf(stderr, "\nWARNING:  longitude changed sign.  If the section crosses the dateline, specify -F to force longitudes into range 0-360.");
               warning_given = 1;
            }
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

   /*  Reads each profile in a HydroBase file and computes property values
       at each standard level.    */
{
   int error, i, j, mainprops;
   double x;
   double  *yprop;
   double  *zprop;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* check that pr,  te, and sa are available ... */
  
     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
 	 
	 
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
 
       mainprops = 1;
       if (!(available(PR, &hdr) && available((enum property) tindex, &hdr) && available(SA, &hdr))) {
         mainprops = 0;
       }
       
      yprop = *get_prop(yopt, mainprops);
      zprop = *get_prop(zopt, mainprops);
      x = get_x_prop(xopt);

    
 /* output the station */

    if ((zprop != NULL) &&  (yprop != NULL) ) {
       for (i = 0; i < hdr.nobs; ++i) {
          if ((zprop[i] > FLAG) && (yprop[i] > FLAG) )
             fprintf(outfile,"%10.3lf %10.4lf %10.4lf\n", x, yprop[i], zprop[i]);
       }
    }
    
    free_hydro_data(&station);
    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

