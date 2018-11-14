/*  hb_getpos.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated Feb 2000 to ANSI standards
			     updated to HB3 Dec 2009
................................................................................
................................................................................
.  Extracts the lat, lon from HydroBase data files and writes them out to 
.  the stdout device.
.
................................................................................
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


#define    EXTENT    ""
#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    DIR       ""
#define    PRINT_MSG  1               /* 0 or 1 */

main (int argc, char **argv)
{
   int     i, nfiles, curfile; 
   int n_extents, curext; 
   int     error, status;
   short   bopt, latfirst, output_stano;
   float  xmin, xmax, ymin, ymax;
   int     xdateline;
   int    lonpos, lonneg;
   char   *dir, *st;
   char  **extent_list;
   struct HYDRO_HDR hdr;
   struct HYDRO_DATA data;
   void print_usage(char *);
   FILE * infile, *outfile;



/*  set these default values */

    dir = DIR;
    bopt = 0;
    output_stano = 0;
    latfirst = 0;
    lonpos = lonneg = 0;
    error = 0;
    xdateline = 0;
    nfiles = 0;
    curfile = 1;
    outfile = stdout;
    
    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;
    extent_list[0] = EXTENT;

/* parse command line arguments... */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':       /* multiple file extents are supported */
                        extent_list[n_extents] = &argv[i][2];
			if (++n_extents == MAXEXT) {
                             fprintf(stderr,"\nToo many different file extents. Max allowed is %d\n", MAXEXT);
			 exit(1);
			}
                        break;
	       case 'F':
	             switch (argv[i][2]) {
		         case '+':
			    lonpos = 1;
			    break;
			 case '-':
			    lonneg = 1;
			    break;
			 default:
			    error = 1;
		     }
		     break;

               case 'B':                    /* get grid bounds */
                        bopt = 1;
                        st = &argv[i][2];
                           if (*st == '/')
                               ++st;
                        error = (sscanf(st,"%f", &xmin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &xmax) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymin) != 1);
                        while (*(st++) != '/')
                           ;  
                        error += (sscanf(st,"%f", &ymax) != 1);
                        if (xmin > 0 && xmax < 0)
                           xmax += 360;
                        if (xmax > 180)
                           xdateline = 1;
                        break;
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			
			}
                        break;
		case 'S':
		        output_stano = 1;
			break;
               case ':':
                        latfirst = 1;
                        break;
               case 'h':
                        print_usage(argv[0]);
			exit(0);
                        break;
               default :
                        print_usage(argv[0]);
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */

            if (error ) {
                print_usage(argv[0]);
                fprintf(stderr,"\nError parsing command line args.\n");
                fprintf(stderr,"     in particular: '%s'\n", argv[i]);
                exit(1);
            }
            
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */


/* initialize these... */

   for (i = 0; i < MAXPROP; ++i) {
      data.observ[i] = (double *) NULL;
      data.count[i] = (double *) NULL;
      data.variance[i] = (double *) NULL;
      data.quality[i] = (double *) NULL;
    
   }
   hdr.prop_id = (int *) NULL;

   if (n_extents == 0) ++n_extents;

/* loop for each file */

   do {
     curext = 0;
     do {
       if ( !nfiles) {
          infile = STDIN;
          fprintf(stderr,"\n Expecting data from stdin....  ");
       }
       else {
          infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
          if (infile  == NULL ) 
           goto NEXTFILE;
       }
   
     /* loop for each station */

     while ((status = get_station(infile, &hdr, &data)) == 0) { 
 
       if (xdateline && hdr.lon < 0)
          hdr.lon += 360;
  
       if (lonneg) {    /* force all longitudes negative */
          if (hdr.lon > 0)
	    hdr.lon -= 360.0;
       }
       
       if (lonpos) {
          if (hdr.lon < 0)
	    hdr.lon += 360.0;
       }
                
       if (bopt) {
          if ((hdr.lon <= xmax) && (hdr.lon >= xmin) && 
           (hdr.lat <= ymax) && (hdr.lat >= ymin) ) {
	     if (output_stano)
	        fprintf(outfile, "%4d ", hdr.station);
		
	     if (latfirst)  
               fprintf(outfile, "%8.3f %8.3f\n", hdr.lat, hdr.lon);
	     else
               fprintf(outfile, "%8.3f %8.3f\n", hdr.lon, hdr.lat);
	       
	       
	  }     
       }
       else {
	     if (output_stano)
	        fprintf(outfile, "%4d ", hdr.station);
	     if (latfirst)  
               fprintf(outfile, "%8.3f %8.3f\n", hdr.lat, hdr.lon);
	     else
               fprintf(outfile, "%8.3f %8.3f\n", hdr.lon, hdr.lat);
       }

     }  /* end while */ 

     report_status(status, stderr);
     if (nfiles)
        fclose(infile);

NEXTFILE:
      ;
      
     } while (++curext < n_extents);
      
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n End of hb_getpos. \n");
   exit(0);

}  /* end of main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nExtracts latitude/longitude from HydroBase data files");
   fprintf(stderr,"\n    and writes it to the stdout device -- or an output file if specified with -O .\n");
   fprintf(stderr,"\nUsage:  %s list_of_filenames", program);

   fprintf(stderr," [-Ddirname] [-Eextent] [-F<+|->] [-Bminlon/maxlon/minlat/maxlat] [-O<outfile>] [-S] [-:]");
   fprintf(stderr,"\n If no infiles are specified, input is expected from stdin");
   fprintf(stderr,"\n  [-B]  : specifies optional boundaries.");
   fprintf(stderr,"\n        ex: -B-90/0/0/65 ");
   fprintf(stderr,"\n  [-D]  : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n  [-E] : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n  [-F]  : Force longitudes to be positive or negative");
   fprintf(stderr,"\n        ex: -F+  OR  -F-  default is mixed sign ");
   fprintf(stderr,"\n  [-O] : output filename, default is output to stdout");
   fprintf(stderr,"\n  [-S]  : output station # in first column.");
   
   fprintf(stderr,"\n  [-:]  : outputs latitude/longitude. ");
   fprintf(stderr,"\n        default is lon/lat ");
   fprintf(stderr,"\n\n");  
   return;
}

