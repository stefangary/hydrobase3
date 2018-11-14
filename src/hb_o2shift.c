/*  hb_o2shift.c

................................................................................
                          *******  HydroBase3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C 1999
			     updated to HydroBase3 2012
			     Supports variance,count,quality arrays 
................................................................................
____________________________________________________________________________
  USAGE:  

 hb_o2shift filename(_roots) -S<offset> [-D<dirname>] [-E<file_extent>] [-L<omin/omax>] [-X] [-Z<depthmin/depthmax>]

  list of filenames MUST be first argument!

 -S : amount to shift oxygen;

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.ctd
 -L : specifies oxygen range
 -X : oxygen shift specified in ml/l (default is micromoles / kg
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

int mflag, use_ml;
double shift;
double zmin, zmax;
double omin, omax;

main (int argc, char **argv)
{
   short   sopt;
   int     index, nprops = 0;
   int     i;
   int     curfile = 1, nfiles = 0; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    *extent, *dir, *st;
   FILE    *infile;
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
    omin = 0;
    omax = 500;
    infile = stdin;
    use_ml = 0;

/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      station.observ[i] = (double *) NULL;
      station.variance[i] = (double *) NULL;
      station.count[i] = (double *)NULL;
      station.quality[i] = (double *)NULL;
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
                       error = (sscanf(&argv[i][2],"%lf/%lf", &omin, &omax) == 2) ? 0 : 1;
                        break;

               case 'X':
	                use_ml = 1;
                        break;
	       
               case 'Z':
                       error = (sscanf(&argv[i][2],"%lf/%lf", &zmin, &zmax) == 2) ? 0 : 1;
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

   if (  !sopt) {
       fprintf(stderr,"\nYou must specify the -S[m]<o2_shift> argument.\n");
       exit(1);
   }
   
   fprintf(stderr,"Oxygen units are ");
   
   if (use_ml)
        fprintf(stderr,"ml/l \n");
   else
        fprintf(stderr,"micromoles/kg \n");
   	
   
   fprintf(stderr,"Using depth limits %.1lf to %.1lf  and oxygen range %.3lf to %.3lf\n", zmin, zmax, omin, omax);

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

   fprintf(stderr,"\n\nEnd of hb_o2shift.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n\n***************************************************");
   fprintf(stderr,"\n hb_o2shift shifts the oxygen variable by a specified amount");
   fprintf(stderr,"\n***************************************************");

                                                       
   fprintf(stderr,"\nUsage:  %s filename_root(s)  -S[m]<offset>  [-D<dirname>] [-E<file_extent>]  [-L<omin/omax>] [-X] [-Z<depthmin/depthmax>]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument or station file is");
   fprintf(stderr,"\n  is expected to be piped from stdin.");
   fprintf(stderr,"\n    -S  : amount to shift oxygen values OR -Sm to multiply oxygen");
   fprintf(stderr,"\n          ex:  -S-3");
    fprintf(stderr,"\n          ex:  -Sm10.0");
  fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-L] : specifies limits of oxygen range <omin/omax> to shift");  
   fprintf(stderr,"\n            ex: -L200/280 ");
   fprintf(stderr,"\n   [-X] : Use units of ml/l for oxygen shift (default is micromoles/kg");  
   fprintf(stderr,"\n   [-Z] : specifies limits of depth range <dmin/dmax> to perform shift");  
   fprintf(stderr,"\n            ex: -Z0/1000 ");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(FILE * file)
   /*  Reads each station in a HydroBase file and adjust oxygen values
       at each standard level.    */
{
   int error, i, tindex, oindex;
   double dlat;


/* read each station in file ... */

    while ((error = get_station(file, &hdr, &station)) == 0) {

  /* is oxygen available? */

       if (available(O2, &hdr) || available(OX,&hdr) ) {
       
          tindex = (int)T90;
	  
	  if (!available(T90, &hdr) && available(TE, &hdr) ) 
              create_t90_arrays( &hdr, &station);
	       
	 
          if (use_ml) {
	     
	     oindex = (int)OX;
	     
	     if (!available(OX, &hdr) ) {
	        /* convert  O2 -> OX*/
		
	        free_and_alloc(&station.observ[(int)OX], hdr.nobs);
	        
		for(i = 0; i < hdr.nobs; ++i)
		   station.observ[(int)OX][i] = ox_kg2l(station.observ[(int)O2][i], station.observ[(int)PR][i], station.observ[(int)T90][i], station.observ[(int)SA][i]);
	        
	        for (i = 0; i < hdr.nprops; ++i) {
		    if (hdr.prop_id[i] == (int)O2)
		       hdr.prop_id[i] = (int) OX;
		}
		free(station.observ[(int)O2]);
	     }
	  }
	  
	  if (!use_ml) {
	     oindex = (int)O2;
	     
	     if (!available(O2, &hdr) ) {
	      
	      /* convert OX  -> O2*/
	        free_and_alloc(&station.observ[(int)OX], hdr.nobs);
	        
		for(i = 0; i < hdr.nobs; ++i)
		   station.observ[(int)O2][i] = ox_l2kg(station.observ[(int)OX][i], station.observ[(int)PR][i], station.observ[(int)T90][i], station.observ[(int)SA][i]);
	        
	        for (i = 0; i < hdr.nprops; ++i) {
		    if (hdr.prop_id[i] == (int)OX)
		       hdr.prop_id[i] = (int) O2;
		}
		free(station.observ[(int)OX]);
	     }
	  
	  
	  }
       
          for (i = 0; i < hdr.nobs; ++i)  {
	     if (station.observ[oindex][i] >= omin && station.observ[oindex][i] <= omax  && station.observ[(int)DE][i] >= zmin && station.observ[(int)DE][i] <= zmax )  {
	     
	         if (mflag) 
                    station.observ[oindex][i] *= shift;
		 else
                    station.observ[oindex][i] += shift;
	     }
	  }
       }       

       write_hydro_station(stdout, &hdr, &station);
       free_hydro_data(&station);
       
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

