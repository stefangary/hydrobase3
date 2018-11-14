/* hb_xyzprop.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
			     created for HydroBase 3 Dec 2009
			     updated for T90 variables Nov 2011
................................................................................
................................................................................
.   creates a list of  property triplets for each observation level in
.   a file of hydro stations.   Separates stations with a '>' to permit
.   points to be connected in GMT. 
.   -F option specifies a minimum value of X,Y,Z as a limit (default is -9 to avoid outputting 
.    missing values).
.
.  If no input, output files are specified, reads from STDIN; writes to STDOUT
. 
.                          
................................................................................
................................................................................
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "hb_paths.h"


/* input file pathnames */

#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    EXTENT   ""
#define    DIR      ""


   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;
double s_pref;           /* ref lev for computing a non-standard sigma level */
double pe_pref;           /* ref lev for computing potential energy */
double ht_pref;           /* ref lev for computing dynamic height */
int window, w_incr;     /* used to specify pr window for gradient properties*/

struct GAMMA_NC ginfo;   /* used for neutral density */

double flag;
char delimiter = '>';
int tindex;

/*  prototypes for locally defined functions */

void print_usage(char *);
int  parse_prop_list(char *);
double **get_property(int);

int main (int argc, char **argv)

{
   double  d;
   int    print_msg = 1;       /* print message in file open routines*/
   int  i, n, error;
   int nfiles, curfile;
   int n_extents, curext; 
   short zflag, xflag, yflag, mflag;
   int xindex, yindex, zindex;
   int scan_ok, any_scans, status;
   FILE *infile, *outfile;
   enum property X, Y, Z;
   char *s, str_format[30];
   char *dir;
   char  **extent_list;

   double **xaddr, **yaddr, **zaddr;
   double *xprop, *yprop, *zprop;

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/* initialize these */

   xflag = yflag = zflag = mflag = 0;
   flag = HB_MISSING + 0.1;
   
   extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
   n_extents = 0;
    extent_list[0] = EXTENT;
   dir = DIR;
   infile = stdin;
   outfile = stdout;
   
   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
   }   
   hdr.prop_id = (int *) NULL;
   nfiles = 0;
   curfile = 1;
   error = 0;
   s_pref = -1;
   window = 100;
   w_incr = 10;
   tindex = (int)T90;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':             /* multiple file extents supported */
                        extent_list[n_extents] = &argv[i][2];
			if (++n_extents == MAXEXT) {
                             fprintf(stderr,"\nToo many different file extents. Max allowed is %d\n", MAXEXT);
			 exit(1);
			}
                        break;

                case 'F':                    /* set minimum acceptable value*/
                        s = &argv[i][2];  
                        error += (sscanf(s,"%lf", &flag) != 1);
                        break;

               case 'M':                   /* specify multiple segment option */
                        mflag = 1;
                        break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
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
			       error = 1;
			}
			break;
               case 'X':                    /* get x-property */
                        xflag = 1;
                        xindex = parse_prop_list(&argv[i][2]);
                        if (xindex < 0) 
                             error = 1;
			 else
                             X = (enum property) xindex;
                        break;
               case 'Y':                    /* get y-property */
                        yflag = 1;
                        yindex = parse_prop_list(&argv[i][2]);
                        if (yindex < 0) 
                             error = 1;
                        else 
			     Y = (enum property) yindex;
                        break;
               case 'Z':
                        zflag = 1;
                        zindex = parse_prop_list(&argv[i][2]);
                        if (zindex < 0) 
                             error = 1;
			 else
                            Z = (enum property) zindex;
                        break;
               case 'h':                    /* help */
                        print_usage(argv[0]);
                        exit(0);
               default:
                       error = 1;

          }    /* end switch */
	  
          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             fprintf(stderr,"Use -h for a complete USAGE message.\n");
             exit(1);
          }
	  
       }  /* end if */
       
       else {
          ++nfiles;
       }

   }  /* end for */

   if (! (xflag && yflag && zflag) ) {
       fprintf(stderr,"\nYou must specify -X, -Y  and -Z properties ...\n");
       exit(1);
   }
   if (tindex == (int) T90)
      fprintf(stderr,"Using T90 temperatures\n");
   else
      fprintf(stderr,"Using T68 temperatures\n");
   
   sprintf(str_format," %c%d.%dlf  %c%d.%dlf  %c%d.%dlf\n", 
      '%',get_field_width(xindex), get_field_precis(xindex), 
      '%', get_field_width(yindex),get_field_precis(yindex), 
      '%', get_field_width(zindex),get_field_precis(zindex));
   
     if (n_extents == 0) ++n_extents;

      /* loop to read each file */
      
   do {
   
     curext = 0;
     do {
        if (nfiles) {
          infile = open_hydro_file(dir, argv[curfile], extent_list[curext], print_msg);
          if (infile == NULL)
             goto NEXTFILE;
        } 
        else {
             fprintf(stderr,"Expecting station input from stdin...\n" ); 
        }
     
      /* loop to read each station */
      
     while ((status = get_station(infile, &hdr, &station)) == 0) { 

       xaddr = get_property(xindex); 
       yaddr = get_property(yindex);
       zaddr = get_property(zindex);
        
       any_scans = 0;
       if (xaddr != NULL && yaddr!= NULL && zaddr != NULL) {
          xprop = *xaddr;
	  yprop = *yaddr;
          zprop = *zaddr;
	  
          for (n = 0; n < hdr.nobs; ++n) {
          
            scan_ok = (xprop[n] > flag) && ( yprop[n] > flag) &&  (zprop[n] > flag);
            
            if (scan_ok) {
               ++any_scans;
               fprintf(outfile,str_format, xprop[n], yprop[n], zprop[n]);
            }              

          }  /* end for */
          
          if (any_scans && mflag) {           /* station separator */
             fprintf(outfile,"%c \n", delimiter ); 
          }
                    
        } /* end if */ 
	
	free_hydro_data(&station);
      }  /*end while */
       
      report_status(status, stderr);
 
      if (nfiles) {
         fclose(infile);
      }
      
      NEXTFILE:
        ;

     } while (++curext < n_extents);
   
   } while (curfile++ < nfiles );
   

}  /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nCreates a listing of x,y,z property triplets. ");
   fprintf(stderr,"\n Use the -M option to separate profiles with a >.");
   fprintf(stderr,"\nUsage:  %s -X<x_property> -Y<y_property> -Z<z_property> [-D<dirname>] [-E<file_extent>] [-M<character>] [-O<outfile>] ", program);
   fprintf(stderr,"\n\n    -X  :  property mnemonic to specify X");
   fprintf(stderr,"\n    -Y  :  property mnemonic to specify Y");
   fprintf(stderr,"\n    -Z  :  property mnemonic to specify Z");
   fprintf(stderr,"\n   OPTIONS:");
   fprintf(stderr,"\n   [-D] : specifies directory for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n   [-F] : specify flag (minimum value) for output of all properties.  default [%.2lf]", flag);
   fprintf(stderr,"\n   [-M] : separate stations with > symbol (default) or a specified character  ex: -M%%");
   fprintf(stderr,"\n   [-O] : output filename, default is output to stdout");
   fprintf(stderr,"\n   [-T] : set temperature scale used in computing derived variables to either IPTS-68 (-T68) or ITS-90 (-T90) (default is ITS-90)");
   fprintf(stderr,"\n          ex: -T68  ");
   fprintf(stderr,"\n   [-h] : help -- prints this message");
   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
int parse_prop_list(char *st)

   /*  Parsesa property mnemonic and returns prop index.
         An error will cause an exit.
        
       char *st:    the list with -P stripped off  */
{
   int index, nprops, i;
   char prop[7];
   double ref_val;

   index = -1;
   
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
           default:
              ;
        }  /* end switch */
      }

   /* end of special cases */

  
 /* !**!  Special cases for properties */
  
  
  if (index ==(int)GE || index == (int)GN) 
     gamma_nc_init((char *) GAMMA_NC_PATH, &ginfo);
     
 /* end of special cases */
      

  return (index);
}  /* end parse_prop_list() */
 
/*****************************************************************************/
double **get_property(index)
   /*  checks for property, computes it if possible, and returns pointer to start of the property array  */
{
   int error, i, j,pts_avail, ratio_done;
   int *tmp;
   double dlat, dlon;
   double **propaddr;
   
   if (available((enum property)index, &hdr)) {
     return (get_prop_array(&station, index));
   }
   
   propaddr=(double **)NULL;

  /* Not available so attempt to compute it, 
     First check if basic properties pr, te, and sa are available ... */
  
     if ( ! available(PR, &hdr) && available(DE, &hdr)) 
         create_pr_arrays( &hdr, &station);
     
     if (!available((int)T90, &hdr) && available(TE,&hdr) ) 
         create_t90_arrays( &hdr, &station);
     
     
     if (!available(TE, &hdr) &&  available(T90,&hdr) ) 
         create_t68_arrays( &hdr, &station);

     pts_avail = (available(PR, &hdr) && available((enum property)tindex, &hdr) 
                 && available(SA, &hdr));

     if (!pts_avail)   /* no need to go further */
        return (propaddr);
     
     
  /* compute appropriate properties at each level in station ... */
 /* !**! Special cases for individual properties... */
 
     ratio_done = 0;
     i = index;   /* just an alternative variable */
     switch ((enum property)index) {
     
	    case DE:
                create_de_arrays( &hdr, &station);
             case OX:
               if (available(O2, &hdr)) 
                   create_ox_arrays( &hdr, &station);
               break;
             case O2:
               if (available(OX, &hdr)) 
                    create_o2_arrays( &hdr, &station);
               break;
             case TH:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_theta(hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA] );
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
             case VA:
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_svan( hdr.nobs, station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA]);
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
               free_and_alloc(&station.observ[i], hdr.nobs);
               compute_sound_vel( station.observ[i], station.observ[(int)PR], station.observ[tindex],station.observ[(int)SA], hdr.nobs);
               break;

             case PV: 
               free_and_alloc(&station.observ[i], hdr.nobs);
               dlat = (double) hdr.lat;
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                 if (station.observ[i][j] > -99.) 
                   station.observ[i][j] *= station.observ[i][j];
               }
               po_vort(station.observ[i], station.observ[i], hdr.nobs, dlat);

               break;
               
             case BF:
               free_and_alloc(&station.observ[i], hdr.nobs);
               buoy_freq(station.observ[i], station.observ[(int)PR], station.observ[tindex], station.observ[(int)SA], hdr.nobs, window, w_incr);
               for (j = 0; j < hdr.nobs; ++j) {
                    if (station.observ[i][j] > -9998.)   /* convert the units */
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
	       
             default:   /* nothing done */
               return (get_prop_array( &station, index));
      } /* end switch */
	  
   
     propaddr = get_prop_array( &station, index);
     if (propaddr == NULL || *propaddr == NULL)
         return((double **)NULL);

   /* update hdr.prop_id and nprops fields */
	     
     tmp = hdr.prop_id;
     hdr.prop_id = (int *)calloc(++hdr.nprops, sizeof(int));
     for (i=0; i < station.nprops; ++i)
	     hdr.prop_id[i] = tmp[i];
	  
     hdr.prop_id[i] = index;
     free(tmp);   
     ++station.nprops;
     return (propaddr);
  
} /* end get_property() */

