/*  hb_saltshift.c

................................................................................
                          *******  HydroBase3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C 1999
			     updated to HydroBase3 Feb 2012
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_saltshift filename(_roots) -S<offset> [-D<dirname>] [-E<file_extent>] [-L<smin/smax>] [-Z<depthmin/depthmax>]

  list of filenames MUST be first argument!

 -S : amount to shift salinity;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 -L : specifies salt range
 -Z : specifies depth range

____________________________________________________________________________
Shifts the salinity by the amount specified with -S.                                                    ____________________________________________________________________________
*/


#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


/* input file pathnames */

#define    EXTENT   ""
#define    DIR      ""

    /* flags to indicate which properties are to be output and which are
       needed for computation */
      

   /* global variables to store station ... */

struct HYDRO_DATA station;
struct HYDRO_HDR hdr;

int mflag;
double shift;
double zmin, zmax;
double smin, smax;

main (int argc, char **argv)
{
   short   sopt;
   int     index, nprops = 0;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   FILE *infile;
   void    print_usage(char *);
   void    get_hydro_data(FILE *);


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    sopt = mflag = 0;
    error = 0;
    shift = 0;
    zmin = 0;
    zmax = 10000;
    smin = 0;
    smax = 50;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
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

               case 'S':
                        sopt = 1;
			st = &argv[i][2];
			if (*st == 'm') { 
			     mflag = 1;
			     ++st;
			}
                        error = (sscanf(st,"%lf", &shift) == 1) ? 0 : 1;
                        break;
               case 'L':
                       error = (sscanf(&argv[i][2],"%lf/%lf", &smin, &smax) == 2) ? 0 : 1;
                        break;

               case 'Z':
                       error = (sscanf(&argv[i][2],"%lf/%lf", &zmin, &zmax) == 2) ? 0 : 1;
                        break;
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

   if (  !sopt) {
       fprintf(stderr,"\nYou must specify the -S[m]<salt_shift> argument.\n");
       exit(1);
   }
   
   fprintf(stderr,"Using depth limits %.1lf to %.1lf  and salt range %.3lf to %.3lf\n", zmin, zmax, smin, smax);
   infile = stdin;

/* loop for each input file */

   do {

     if (nfiles) {
       infile = open_hydro_file(dir, argv[curfile], extent, print_msg);
       if (infile == NULL)
          goto NEXTFILE;

     }
     else {
       fprintf(stderr, "\n Expecting input from stdin....");
     }
            /* read each file completely */

 
         get_hydro_data(infile);

NEXTFILE:

         fclose(infile);

   } while (curfile++ < nfiles );

   fprintf(stderr,"\n\nEnd of hb_saltshift.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n\n***************************************************");
   fprintf(stderr,"\n hb_saltshift shifts the salinities by a specified amount");
   fprintf(stderr,"\n***************************************************");

                                                       
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -S[m]<offset>  [-D<dirname>] [-E<file_extent>]  [-L<smin/smax>]  [-Z<depthmin/depthmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or station file is");
   fprintf(stderr,"\n  is expected to be piped from stdin.");
   fprintf(stderr,"\n    -S  : amount to shift salinity values OR -Sm to multiply salinities");
   fprintf(stderr,"\n          ex:  -S-.008");
    fprintf(stderr,"\n          ex:  -Sm10.0");
  fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-L] : specifies limits of salt range <smin/smax> to shift");  
   fprintf(stderr,"\n            ex: -L34.0/35.0 ");
   fprintf(stderr,"\n   [-Z] : specifies limits of depth range <dmin/dmax> to perform shift");  
   fprintf(stderr,"\n            ex: -Z0/1000 ");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(FILE * file)
   /*  Reads each station in a HydroBase file and computes property values
       at each standard level.  This module requires that the HydroBase
       file contains a minimum of pr, de, te, sa observations.  */
{
   int error, i;
   double dlat;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* ensure that pr, de, te, and sa are available ... */

       if (available(SA, &hdr)) {
       
          for (i = 0; i < hdr.nobs; ++i)  {
	     if (station.observ[(int)SA][i] >= smin && station.observ[(int)SA][i] <= smax  && station.observ[(int)DE][i] >= zmin && station.observ[(int)DE][i] <= zmax )  {
	     
	         if (mflag) 
                    station.observ[(int)SA][i] *= shift;
		 else
                    station.observ[(int)SA][i] += shift;
	     }
	  }
       }       
       else {
         fprintf(stderr,"\n>>>>> WARNING!!!  ");
         fprintf(stderr,"Station does not include salinity.");
       }
       

       write_hydro_station(stdout, &hdr, &station);

       for (i = 0; i < MAXPROP; ++i) {
          if (station.observ[i] != NULL) {
             free((void *)station.observ[i]);
             station.observ[i] = NULL;
          }
       }    
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

