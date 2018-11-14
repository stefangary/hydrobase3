/* hb_vrotate.c
................................................................................
                         *******  HydroBase3 *******
................................................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             July 2000 
			     updated to HydroBase3 Feb 2013
................................................................................
................................................................................
.  Rotates velocity vectors (vn and ve) through  by the specified angle (deg).
.  Ve is transformed into Vdo (downstream direction); Vn into Vcr (cross stream
.  direction)
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "hydrobase.h"

#define    DIR     ""
#define    EXTENT   ""
#define    PRINT_MSG 1
#define    PI   3.141592654

void print_usage(char *);

main (int argc, char **argv)
{
   int     curfile = 1, nfiles = 0; 
   int     i, m;
   int     status;
   int     error = 0;
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
   char   *dir, *extent;
   double rotation, phi, vdir, speed;
   double *vnorth, *veast, *vdown, *vcross;
   FILE   *infile;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   

/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    rotation = 0.0;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'A' :
                        error = sscanf(&argv[i][2],"%lf", &rotation) != 1;
                        break;
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent = &argv[i][2];
                        break;
               case 'h':
                        print_usage(argv[0]);
                        exit(0);
               default  :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

   if (!nfiles) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input file(s)!\n");
       exit(1);
   }
   if ( rotation == 0.0) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify the counterclockwise rotation angle from Veast -> Vdownstream direction in degrees [-180 to + 180].\n");
       exit(1);
   }

/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   h.prop_id = (int *) NULL;
   
    rotation *=  PI / 180.0;    /* convert to radians*/
      
 /* loop for each input file */

   do {

     infile = open_hydro_file(dir, argv[curfile], extent, PRINT_MSG);
     if (infile  == NULL) 
       goto NEXTFILE;
     

     /* loop for each station */

     while ((status = get_station(infile, &h, &data)) == 0) { 

        if (data.observ[(int)VN] != NULL && data.observ[(int)VE] != NULL) {
	    vdown = (double *) calloc(h.nobs, sizeof(double));
	    vcross = (double *) calloc(h.nobs, sizeof(double));
	    veast = data.observ[(int)VE];
	    vnorth = data.observ[(int)VN];
	   
	    for (i = 0; i < h.nobs; ++i) {
	       vdown[i] = vnorth[i];
	       vcross[i] = veast[i];
	       if (vnorth[i] > -99.0 && veast[i] > -99.0) {
	            speed = sqrt(veast[i] * veast[i] + vnorth[i] * vnorth[i]);
	            vdir = atan2(vnorth[i], veast[i]);
		    phi = vdir - rotation;
	       
	            vdown[i] = speed * cos(phi);
		    vcross[i] = speed * sin(phi);
	       }
	    }
	    
	   free(vnorth);
	   free(veast);
	   data.observ[(int)VN] = vcross;
	   data.observ[(int)VE] = vdown;
	
	}         

        error = write_hydro_station(stdout, &h, &data);
        if (error) {
          fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
          exit(1);
        }
	
        for (i = 0; i < MAXPROP; ++i) {
           if (data.observ[i] != NULL) {
              free((void *)data.observ[i]);
              data.observ[i] = NULL;
           }
        }    

     }  /* end while */

     report_status(status, stderr);
     fclose(infile);

NEXTFILE:
     ;
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n\nEnd of %s.\n\n", argv[0]);
   exit(0);

} /* end main() */




/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filename_root(s)", program);
   fprintf(stderr," -A<angle_of_rotation>  [-Ddirname] [-Eextent]  ");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n    -A  : countclockwise angle of rotation from Veast -> Vdownstream directions, in degrees [-180 to +180] ");
   fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -h  : help -- prints this message");
  fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
