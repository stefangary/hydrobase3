/* hb_mysort.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             Updated March 2000
                             Updated June 2009
			     SFG modifying hb_mssort.c
................................................................................
............................................................................

Sorts file(s) of HydroBase stations by year and month.
Up to MAXF files can be simultaneously opened. Any stations
that do not fit into any of the open files are stored in a
file called msextra.dat -- which can be sorted
subsequently. A summary of numbers of stations in each file is written to
stdout.

Same as hb_mssort except that sorting is by year and month.
____________________________________________________________________________
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include "hydrobase.h"

#define SLASH    "/"
#define   PMODE    0666   /* read & write permission for output files */
#define   DIR     ""                 
#define   EXTENT   ""
#define   MAXEXT   10  /* maximum # of different file extents specified */
#define   PRINT_MSG  1    /* 0 or 1 */
#define   MAXF    1000   /* max number of open files at one time */

/* List of outfile ID's - an ID is assigned to each opened file. */
FILE *outfile[MAXF + 1];
char *myname[MAXF + 1];

/* Vector storing the year and month associated with each open file.*/
int year_list[MAXF + 1];
int month_list[MAXF + 1];

/* Vector storing the number of files written to each open file.*/
int count[MAXF + 1];

/* String vector storing the name of each month. */
char *month_name[12] = { "01", "02", "03",
			 "04", "05", "06",
			 "07", "08", "09",
			 "10", "11", "12" };

int nextfile;

struct HYDRO_HDR hdr;
struct HYDRO_DATA data;

  /* prototypes for locally defined functions */

void print_usage(char *);  
FILE *openfile(char *, char *, char *, int );
int sort10(char *, char *, int);
int sort5( char *, char *, int);
int sort2( char *, char *, int);
int sort1( char *, char *, int);


main (int argc, char **argv)
{
FILE *infile;
int j, error;
int  mode;
int  nfiles, i;
int n_extents, curext; 
int  curfile = 1,  status;
char *fname, name[80];
char *outpath, *extout, *dir;
char  **extent_list;

/* check for command line arguments... */

if (argc < 2) {
   print_usage(argv[0]);
   exit(1);
}

/*  set these default values */

    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));
    n_extents = 0;
    extent_list[0] = EXTENT;
    
    dir = DIR;
    outpath = DIR;
    extout = EXTENT;
    mode = NOCLOBBER;
    nfiles = 0;
    nextfile = 0;
    error = 0;

/* parse command line arguments... */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
          switch (argv[i][1]) {
               case 'D':
                        dir = &argv[i][2];
                        break;

               case 'E':
                        extent_list[n_extents] = &argv[i][2];
			++n_extents;
                        break;

               case 'O':
                        outpath = &argv[i][2];
                        break;
               case 'N':
                        extout = &argv[i][2];
                        break;
               case 'A':
                        mode = APPEND;
                        break;
               case 'T':
                        mode = OVERWRITE;
                        break;
               case 'h': 
                  print_usage(argv[0]); 
                  exit(0);
		  
               default :
                        error = 1;
                        
          }  /* end switch */
            
          if (error ) {
             fprintf(stderr,"\nError parsing command line args.\n");
             fprintf(stderr,"     in particular: '%s'\n", argv[i]);
             fprintf(stderr,"Use -h for a complete usage statement.\n");
             exit(1);
          }
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

   if (n_extents == 0) ++n_extents;
     
/* initialize the following */

 for (j=0; j <= MAXF; ++j){
       count[j] = 0;
       year_list[j] = -1;
       month_list[j] = -1;
       outfile[j] = NULL; 
}

for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
hdr.prop_id = (int *) NULL;


/*********************************************
 * Open the scratch file myextra.dat
 *********************************************/

strcpy(name, outpath);
j = strlen(name);              /* check for easily forgotten slash in path */
if (j && (name[j-1] != *SLASH)) 
   strncat(name,SLASH,2);
strncat(name, "myextra.dat", 12);

if ((outfile[MAXF] = create_hydro_file(name,mode)) == NULL) {
       switch (mode) {
          case APPEND:
                fprintf(stderr,"\nUnable to open %s in append mode.", name);
                exit(1);
          case OVERWRITE:
                fprintf(stderr,"\nUnable to open %s.", name);
                exit(1);
          default:
                fprintf(stderr,"\nUnable to open %s.", name);
                fprintf(stderr,"  It may already exist...\n");
                exit(1);
       } /* end switch */
}
      
/******************************************
 * Loop for each input file.
 ******************************************/

do {

  curext = 0;
  do {

   if (nfiles == 0) {
       infile = stdin;
       fprintf(stderr,"\nExpecting station input from stdin ....\n");
   }
   else {
   
      infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
      if (infile == NULL) 
         goto NEXTFILE;
   }

   /* Loop to read each station */
   while ((status = get_station(infile, &hdr, &data)) == 0) {
      /* Find which file to write to. */

      /* Search for an already open file with this year and month.
       * Note that nextfile is initialized at zero, so in the
       * first instance, this while loop will exit right away
       * and engage the open an output file lines below. If a file
       * is found with this exact match of year and month in
       * the year and month lists, then the loop exits storing
       * j as the index as the file to write to. */
      j = 0;
      while ((j < nextfile) && !((year_list[j] == hdr.year)&&(month_list[j] == hdr.month))) {
          ++j; 
      }   /* End of search for already open file. */

      /* We have not found an open file with this year and 
       * month, so open a new output file. */
      if ((j == nextfile) && (j < MAXF)) {
          year_list[j] = hdr.year;
          month_list[j] = hdr.month;
          myname[j] = (char *) calloc(100, sizeof(char));
	  sprintf(myname[j],"%s.%4d.%s", argv[curfile],hdr.year,month_name[hdr.month-1]);
          outfile[j] = openfile(outpath, myname[j], extout, mode);

	  /* Augment the number of open files. */
          ++nextfile;
      } 

      /* Write station to appropriate outfile,
       * file number j in the lists.*/
      write_hydro_station(outfile[j], &hdr, &data);
      ++count[j];

   }  /* end while */
   report_status(status, stderr);

   fclose(infile);
NEXTFILE:
    ;
    
  }   while (++curext < n_extents);

}  while (curfile++ < nfiles);


/*  write summary of distribution to stderr and close output files */

fprintf(stdout, "\n\n accounting of stations written to each file ...\n");
/* WORKING HERE */
for (j=0; j < nextfile; ++j) {
       fprintf(stdout, "   %s%s     %d \n", myname[j], extout, count[j]);
       fclose(outfile[j]);
}

fprintf(stderr,"  msextra.dat    %d \n", count[MAXF]);
fclose(outfile[MAXF]);                            /* close scratch file */

exit(0);
} /* end of main */


/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n%s sorts file(s) of HydroBase stations by year", program);
  fprintf(stderr,"\nand month.");
  fprintf(stderr,"\nUp to MAXF files can be simultaneously opened.");
  fprintf(stderr,"\nAny stations that do not fit into an open file");
  fprintf(stderr,"\nare stored in a file called msextra.dat -- which can be");
  fprintf(stderr,"\nsorted subsequently. An accounting of numbers of stations ");
  fprintf(stderr,"\nin each file is written to stdout.\n");
  fprintf(stderr,"\nUsage:  %s list_of_filename_root(s)", program);

   fprintf(stderr," [-Ddir] [-Eextent][-Ooutpath] [-Nnew_extent] [-A] [-T]");
   fprintf(stderr," > logfile\n");
  fprintf(stderr,"\n   [-D] : specifies directory of input files (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n   [-O] : specifies directory of output files (default is ./) ");
   fprintf(stderr,"\n        ex: -O../natl/ ");
   fprintf(stderr,"\n   [-N] : specifies output file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -N.dat ");
   fprintf(stderr,"\n   [-A] : append to existing files (default is do not alter an existing file.)");
   fprintf(stderr,"\n   [-T] : truncate existing files (default is do not alter an existing file.)");
   fprintf(stderr,"\n   [-h] : help... prints this message.)");
   fprintf(stderr,"\n\n");  
   return;
} /* end print_usage() */

/****************************************************************************/
FILE *openfile(char *outpath, char *n, char *extent, int mode)
/*
 Opens an existing file or creates a new file for output. The files are named
 <outpath><n>.<extent>.
 arguments:
          n:
    outpath:      null terminated strings 
     extent:
       mode:     APPEND, OVERWRITE or NOCLOBBER for output file 
 */
{
   char fname[200];
   int i;
   FILE * newfile;

/* move beyond any leading dot in file extent... */
   if (*extent == '.')
      ++extent;

   /* check pathname for easily forgotten slash */

   i = strlen(outpath);
   if (i  && (outpath[i-1] != *SLASH)) 
     sprintf(fname, "%s%c%s.%s", outpath, *SLASH, n, extent);  
   else
     sprintf(fname, "%s%s.%s", outpath, n, extent);  
      

   if ((newfile = create_hydro_file(fname,mode)) == NULL)  {
       switch (mode) {
          case APPEND:
                fprintf(stderr,"\nUnable to open %s in append mode.", fname);
                exit(1);
          case OVERWRITE:
                fprintf(stderr,"\nUnable to open %s.", fname);
                exit(1);
          default:
                fprintf(stderr,"\nUnable to open %s.", fname);
                fprintf(stderr,"  It may already exist...\n");
                exit(1);
       }
   }
   switch (mode) {
          case APPEND:
                fprintf(stderr,"\nOpening or appending %s.", fname);
                return(newfile);
          case OVERWRITE:
                fprintf(stderr,"\nOpening or overwriting %s.", fname);
                return(newfile);
          default:
                fprintf(stderr,"\nOpening %s.", fname);
                return(newfile);
   }

}  /* end of openfile */

