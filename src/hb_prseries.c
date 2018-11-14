/*  hb_prseries.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI Mar 2000
			     updated to HydroBase 3 Dec 2009
................................................................................
____________________________________________________________________________
  USAGE:  

hb_prseries filename(_roots) -P<prs_int> [-D<dirname>] [-E<file_extent>]

  list of filenames MUST be first argument or input expected from stdin

 -P : pressure interval for output series.

 -D : specifies directory for input files (default is current dir)  
      ex: -D/home/pete/data

 -E : specifies input_file extent (default is no extent) ex: -E.sum
 
 -O : specifies output file

____________________________________________________________________________
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
#define    MAXEXT   10  /* maximum # of different file extents specified */

   /* global variables to store station ... */

struct HYDRO_DATA data, newdata;
struct HYDRO_HDR hdr;
int prsint;
FILE *outfile;

  /* prototypes for locally defined functions */
  
void print_usage(char *);
void get_hydro_data(FILE * );


main (argc, argv)
int   argc;
char *argv[];
{
   int     index, nprops = 0;
   int      i, count;
   int     curfile = 1, nfiles = 0; 
   int n_extents, curext; 
   int     print_msg = 1;       /* print message in file open routines*/
   int     error, prompt;
   char    **extent_list, *dir, *st;
   FILE *infile;


/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
 
/*  set these default values */

    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;


    dir = DIR;
    extent_list[0] = EXTENT;
    prsint = 0;
    error = 0;
    outfile = stdout;
 
/* initialize these ... */

   for (i = 0; i < MAXPROP; ++i) {
      data.observ[i] = (double *) NULL;
      data.count[i] = (double *) NULL;
      data.variance[i] = (double *) NULL;
      data.quality[i] = (double *) NULL;
     newdata.observ[i] = (double *) NULL;
     newdata.count[i] = (double *) NULL;
     newdata.variance[i] = (double *) NULL;
     newdata.quality[i] = (double *) NULL;
   }
   hdr.prop_id = (int *) NULL;


/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent_list[n_extents] = &argv[i][2];
			++n_extents;
                        break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			}
                        break;
               case 'P':
                        
                       error = (sscanf(&argv[i][2],"%d", &prsint)==1) ? 0:1;
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

   if ( !prsint) {
       fprintf(stderr,"\nYou must specify input file(s) and a prs_int.\n");
       exit(1);
   }
   if (!nfiles) {
       fprintf(stderr,"\nExpecting input from stdin ... ");
       infile = stdin;
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

NEXTFILE:

     if (nfiles > 0) 
         fclose(infile);
	 
	 
     } while (++curext < n_extents);

   } while (curfile++ < nfiles );

   fprintf(stderr,"End of hb_prseries.\n");
   exit(0);

}   /* end main */

/****************************************************************************/

void print_usage(char *program)
{
 fprintf(stderr,"\n%s reads each station in a HydroBase file and outputs a new file with observations at the specified pressure interval\n", program);

   fprintf(stderr,"\nUsage:  %s filename_root(s)  -P<prsint> [-D<dirname>] [-E<file_extent>] [-O<outfile>] [-h]", program);

   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n  or station input is expected from stdin.");
   fprintf(stderr,"\n    -P  : pressure interval for output array");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-O] : output filename, default is stdout");
   fprintf(stderr,"-h help...... prints this message. \n");

   fprintf(stderr,"\n\n");  
   return;
}

/*****************************************************************************/
void get_hydro_data(FILE *file)

   /*  Reads each station in a HydroBase file and outputs observations at
      the specified pressure interval.    
   */
{
   int error, ii, j, npts, index;
   double  **newaddr;
   double *prop, *newprop, *p;
   double ptop, pbot; 
   
   
/* read each station in file ... */

    while ((error = get_station(file, &hdr, &data)) == 0) {

      if (available(PR, &hdr)) {
      
	    /* allocate space for new arrays */
            for (ii = 0; ii < hdr.nprops; ++ii) {
	       index = hdr.prop_id[ii];
	       newaddr = get_prop_array(&newdata, index);
	       *newaddr = (double *) calloc(6000, sizeof(double));
             }
	     
           /* generate new pressure series */
	   	     
	   p = data.observ[(int)PR];  /* pointer to old pressure array */
	   ptop = p[0];
	   pbot = p[hdr.nobs-1];
	   newdata.observ[(int)PR][0] = ptop;   
	   npts = 1;
	   while (ptop < pbot) {
	      ptop +=prsint;
	      newdata.observ[(int)PR][npts++] = ptop;
	   }
	   newdata.observ[(int)PR][npts++] = pbot;
	   	   
           for (ii = 1; ii < hdr.nprops; ++ii) {
	     if ((index = hdr.prop_id[ii]) != (int)PR) {
	       prop = *get_prop_array(&data, index);
	       newprop = *get_prop_array(&newdata, index);
	       newprop[0] = prop[0];
	       for (j = 1; j < npts-1; ++j) {
	          newprop[j] = hb_linterp(newdata.observ[(int)PR][j], p, prop, hdr.nobs);
	       }
	       newprop[npts-1] = prop[hdr.nobs-1];
	     }
	   
	   }  /* end for */
	   

	   
      /* output the station */
	   hdr.nobs = npts;
	   newdata.nobs = npts;
	   newdata.nprops = hdr.nprops;	
    
          if (hdr.nobs > 0 )
            write_hydro_station(outfile, &hdr, &newdata);
	    
     } /* end if available */

     free_hydro_data(&data);
     free_hydro_data(&newdata);
              
     free(hdr.prop_id); 
     hdr.prop_id = NULL;
       
   }  /* end while */

   if (error > 0)
      report_status(error, stderr);

   return;

} /* end get_hydro_data() */

