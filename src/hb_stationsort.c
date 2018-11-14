/* hb_stationsort.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             1993
                             Updated Mar 2000
			     Updated Dec 2009
................................................................................
............................................................................

Sorts HydroBase station files into individual station files.

____________________________________________________________________________
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"

#define   PRINT_MSG  1    /* 0 or 1 */
#define    DIR     ""
#define    EXTENT   ""
#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    SLASH   '/'

  /* global variables */
struct HYDRO_HDR hdr;
struct HYDRO_DATA data;
int nflag;

  /* prototypes for locally defined functions */
  
void    print_usage(char *);
FILE *openfile( char *, char *, int );


main (int argc, char **argv)
{
register int j;
int   mode;
int  nfiles, curfile;
int  i;
int  status;
int n_extents, curext; 
char *outpath, *rootout;
char  *dir;
char  **extent_list;
FILE *infile, *outfile;

/*  set these default values */

    dir = DIR;
    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;
    extent_list[0] = EXTENT;
    outpath = NULL;
    rootout = NULL;
    mode = NOCLOBBER;
    nflag = 0;
    curfile = 1;
    nfiles = 0;

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
               case 'O':
                        outpath = &argv[i][2];
                        break;
               case 'N':
                        nflag = 1;
                        rootout = &argv[i][2];
                        break;
               case 'T':
                        mode = OVERWRITE;
                        break;
               case 'h':
                        print_usage(argv[0]);
                        exit(0);
               default :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

      
/* initialize the following */

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
           infile = stdin;
           fprintf(stderr,"\n Expecting data from stdin....  ");
        }
        else {
           infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
           if (infile  == NULL) 
            goto NEXTFILE;
      }
      
      /* loop to read each station */

      i = 0;
      while ((status = get_station(infile, &hdr, &data)) == 0) {

          /* open an output file and write station into it */
         ++i; 
         outfile = openfile(outpath, rootout, mode);
	 if (outfile == NULL)
	   continue;
         write_hydro_station(outfile, &hdr, &data);
         fclose(outfile);

      }  /* end while */
   
      report_status(status, stderr);

      fclose(infile);
      
NEXTFILE:
      ;
      
     } while (++curext < n_extents);
      
   } while (curfile++ < nfiles);
   
   fprintf(stderr, "\n  %d stations written to separate files.", i);
   fprintf(stderr, "\n End of %s.\n\n", argv[0]);

 exit(0);
} /* end of main */


/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s ", program);

   fprintf(stderr,"list_of_input_files -[N]new_file_root_name [-Ooutpath] [-T] [-D<dir>] [-E<extent>]");
   
   fprintf(stderr,"\n\nInfile_list must be first arguments or else input is expected from stdin");  
   fprintf(stderr,"\n   [-N] : optionally specifies output file root ");  
   fprintf(stderr,"\n        to which will be appended a number reflecting ");
   fprintf(stderr,"\n        the ordinal position of that station. ");
   fprintf(stderr,"\n        ex: with -Npan1964  ");
   fprintf(stderr,"\n        stations will be named pan1964.0001, pan1964.0002, etc ");
   fprintf(stderr,"\n        Without this option (the default)files are");
   fprintf(stderr,"\n        named by countrycode, ship, cruise, and station");
   fprintf(stderr,"\n   [-D] : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent(s) (default is no extent)");  
   fprintf(stderr,"\n            Use separate -E arguments to specify multiple extents (up to 10 max)");
   fprintf(stderr,"\n            ex: -E.dat -E.ctd -E.flt");
   fprintf(stderr,"\n   [-O] : specifies directory of output files (default is ./) ");
   fprintf(stderr,"\n        ex: -O../natl/ ");
   fprintf(stderr,"\n   [-T] : truncate existing files (default is do not clobber an existing file.)");
   fprintf(stderr,"\n\n");  
   return;
   
} /* end print_usage() */


/****************************************************************************/
FILE *openfile( char *outpath, char *root, int mode)
/*
 Opens an existing file or creates a new file for output. If root is specified,  the files are named <outpath>/root.####, where #### is a 4-digit integer which is incremented each time this function is accessed.  Otherwise the filename is constructed from country/ship/cruise.station information. Returns a pointer to the file.
*/
{
   char fname[1000];
   static count = 0;
   int  i;
   FILE * fptr;    

   
   if (!nflag) {
     root = (char *) malloc(12 * sizeof(char));
     sprintf(root, "%2s%2s%d", hdr.country, hdr.ship, hdr.cruise);
     count = hdr.station;
   }
   else {
      ++count;
   }
   
   if (outpath == NULL) {
         sprintf(fname, "%s.%04d", root, count);  
   }
   else {
       i = strlen(outpath);
       if (outpath[i-1] != SLASH) 
         sprintf(fname, "%s%c%s.%04d", outpath, SLASH, root, count);  
       else
         sprintf(fname, "%s%s.%04d", outpath, root, count);  
   }

   if ((fptr = create_hydro_file(fname,mode)) == NULL)  {
       switch (mode) {
          case OVERWRITE:
                fprintf(stderr,"\nUnable to open %s.", fname);
                exit(1);
          default:
                fprintf(stderr,"\nUnable to open %s.", fname);
                fprintf(stderr,"  It may already exist...\n");
       }
   }
   if (!nflag) {
      free((void *) root);
   }
   return (fptr);

}  /* end of openfile */

