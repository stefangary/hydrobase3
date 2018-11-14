/* hb_random.c
................................................................................
                              *  HydroBase 3 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated to ANSI-C Feb 2000
			     Updated to HydroBase 3: Dec 2008
			     Supports multiple -E options Dec 2009
			     Supports variance,count,quality arrays 
................................................................................
................................................................................
.  Extracts randomized stations HydroBase files.
.
................................................................................
................................................................................
* SFG copied hb_extract and used that framework for hb_random.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>  /* Microseconds for random seed init instead of time.h */
#include "hydrobase.h"

#define    DIR     ""
#define    EXTENT   ""
#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    MAXCRIT  12   /* maximum # of criteria for any of the options */
#define    PRINT_MSG 1

const float DEFAULT_PROB = 0.5; /* Default filter probability */

/* prototypes for locally defined functions */

   void    print_usage(char *);

main (int argc, char **argv)
{

   int     i;
   int     status;
   int    *staOK;
   int    *staID;
   int     nout = 0, nskip = 0;
   int     error;

   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
  /* Microseconds for time, see:
   * http://stackoverflow.com/questions/7343833/srand-why-call-it-only-once/*/
   struct timeval t1;

   /* Variables set during defaults or command
    * line parsing. */
   FILE *infile, *rejectfile, *outfile;
   int     curfile = 1, nfiles = 0;
   int     cursta;
   int     curext, n_extents; 
   int     h_flag=0, f_flag=0, r_flag=0, p_flag=0;
   int     num_sta=0;
   int     num_sta_draw, draw_index, loop_count;
   float   probability;
   char   *dir, **extent_list;
   char   *s;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
/* allocate space for lists */

    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));

/*  set these default values */

    dir = DIR;
    extent_list[0] = EXTENT;
    n_extents = 0;
    outfile = stdout;
    probability = DEFAULT_PROB;

/*  parse the command line arguments */

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
	       case 'F':
		        f_flag = 1;
			/* Read in number of stations. */
                        s = &argv[i][2];
			error = sscanf(s,"%d", &num_sta);
			if ( error != EOF ) {
			  fprintf(stderr,"\nNumber stations set to %d",num_sta);
			}else{
			  fprintf(stderr,"\nERROR: reading -F%s\n",s);
			  exit(1);
			}
			break;
               case 'H':
                        h_flag = 1;
                        break;

               case 'h':
                        print_usage(argv[0]);
			exit(0);
                        break;

               case 'O':
                       outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
                          fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
                          exit(1);
                        }
                        break;
	       case 'P':
		        p_flag = 1;
			/* Read probability from command line. */
                        s = &argv[i][2];
			error = sscanf(s,"%f", &probability);
			if ( error != EOF ) {
			  fprintf(stderr,"\nProbability set to %f",probability);
			}else{
			  fprintf(stderr,"\nERROR: reading -P%s\n",s);
			  exit(1);
			}
			if ( probability > 0.5 ) {
			  fprintf(stderr,"\n ERROR probability > 0.5"); 
			  fprintf(stderr,"\n Due to limited random algorithm"); 
			  fprintf(stderr,"\n high probabilities will likely");
			  fprintf(stderr,"\n cause difficulty in searches.");
			  fprintf(stderr,"\n Instead, specify a lower prob.");
			  fprintf(stderr,"\n and use the rejected stations");
			  fprintf(stderr,"\n with the -R option since they");
			  fprintf(stderr,"\n are just as randomly rejected");
			  fprintf(stderr,"\n as the -O stations are randomly");
			  fprintf(stderr,"\n selected.");
			  exit(1);
			}

			break;
               case 'R':
                        r_flag = 1;
                        rejectfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                        if (rejectfile == NULL) {
                          fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
                          fprintf(stderr,"\nIt may exist already?\n");
                          exit(1);
                        }
                        break;
	        
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

   if ( !nfiles) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input file(s)!\n");
       exit(1);
   }
   
/* initialize some variables */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
     data.count[i] = (double *) NULL;
     data.variance[i] = (double *) NULL;
     data.quality[i] = (double *) NULL;
   }
   h.prop_id = (int *)NULL;

   if (n_extents == 0) n_extents = 1;

   /* Set random seed based on current system time. */
   gettimeofday(&t1,NULL);
   srand(t1.tv_usec*t1.tv_sec);

   /*===============================================
    * GENERAL PLAN:
    * 1) In order to ensure that the same
    *    number of stations are drawn for each
    *    iteration, the stations must all be counted.
    *    (Or, if in fast mode, the number of stations
    *     specified on the command line.)
    * 2) Then, the the probability is used to select
    *    the number of stations which must be randomly
    *    drawn.
    * 3) The stations are randomly selected.
    * 4) The selected stations are written to output.
    *================================================*/

   /*================================================
    * Count all available stations only
    * if we are not in fast mode.
    *================================================*/

   if ( !f_flag ) {
     do { /* loop for each input file */

       curext = 0;
       do { /* Loop for each extension */

	 infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
	 if (infile == NULL) 
	   goto NEXTFILE;
     
	 /* Loop for each station in each file */
	 while ((status = get_station(infile, &h, &data)) == 0) { 
	   /* Add up number of stations */
	   num_sta = num_sta + 1;
	 }  /* end while loop over each station */

	 report_status(status, stderr);
	 fclose(infile);

       NEXTFILE:
	 ;
       
       } while (++curext < n_extents);  /* End of looping over extensions.*/
     } while (curfile++ < nfiles);      /* End of looping over input files.*/
     fprintf(stderr,"\n Found %d input stations.",num_sta);
   }                                    /* End of Fast Mode check. */

   /*====================================================
    * Now that all stations are counted, or we
    * have been provided with the number of
    * stations by the user, we determine which
    * stations to randomly keep.
    *=====================================================*/

   /* Allocate a vector of integers that says yes or no
    * for selection for each station.*/
   staOK = (int *)malloc(sizeof(int)*num_sta);
   
   /* Initialize the vector to be all 0 -> not selected. */
   for (cursta = 0; cursta < num_sta; ++cursta) {
     staOK[cursta] = 0;
   }

   /* Allocate a vector of integers to store a station
    * ID for each station and populate, simply counting
    * from 0 up.  This will be used to index staOK.*/
   staID = (int *)malloc(sizeof(int)*num_sta);
   for (cursta = 0; cursta < num_sta; ++cursta) {
     staID[cursta] = cursta;
   }

   /* Determine the number of stations we wish to
    * randomly draw. Truncation is happening when
    * assigning the double values on the right to
    * the left, so add 0.5 to get rounding.*/
   num_sta_draw = probability*(float)num_sta + 0.5;

   if ( num_sta_draw > num_sta ) {
     fprintf(stderr,"\n ERROR: num_sta_draw > num_sta");
     fprintf(stderr,"\n num_sta_draw = %d",num_sta_draw);
     fprintf(stderr,"\n      num_sta = %d\n",num_sta);
     exit(1);
   }

   if ( num_sta_draw < 0 || num_sta < 0 ) {
     fprintf(stderr,"\n ERROR: num_sta_draw < 0 || num_sta < 0");
     fprintf(stderr,"\n num_sta_draw = %d",num_sta_draw);
     fprintf(stderr,"\n      num_sta = %d\n",num_sta);
     exit(1);
   }

   fprintf(stderr,"\n Will draw %d stations",num_sta_draw);

   /* Finally, pick which stations we want. By looping
    * over each station to be selected and drawing a
    * random number.*/
   cursta = 0;
   loop_count = 0;
   while ( cursta < num_sta_draw ) {
     
     /* Draw a random number, scale the random
      * number to the number of stations available
      * to be drawn.  Since rand() goes [0,RAND_MAX]
      * inclusive, do not round and allow for
      * truncation.  But, if rand() results in RAND_MAX,
      * then the result will be num_sta, which will cause
      * a segmentation fault in staOK (counting starting at
      * 0, going to num_sta-1 in vector).  So check for
      * this special case, and if it occurs, reset to
      * num_sta - 1.*/
     draw_index = (double)num_sta*(double)rand()/(double)RAND_MAX;
     if ( draw_index < 0 ) {
       fprintf(stderr,"\n WARNING: draw_index < 0");
       draw_index = 0;
     }
     if ( draw_index > num_sta ) {
       fprintf(stderr,"\n WARNING: draw_index > num_sta");
       draw_index = num_sta - 1;
     }
     if ( draw_index == num_sta ) {
       /* Normal operation, no cause for warning message. */
       draw_index = num_sta - 1;
     }

     /* Check whether this station has been drawn before
      * or not.  If it hasn't been drawn, then set the flag
      * to draw it and augment the loop counter.  Otherwise,
      * do nothing, which will force the loop to draw
      * another random number.*/
     if ( staOK[draw_index] == 0 ) {
       /* Not yet previously drawn, so set flag and
	* augment the counter.*/
       staOK[draw_index] = 1;
       cursta = cursta + 1;
     }

     /* Always count the number of loops as a check
      * to prevent infinite looping.*/
     loop_count = loop_count + 1;
     if ( loop_count > num_sta ) {
       fprintf(stderr,"\n ERROR: More random draws than input stations!\n");
       exit(1);
     }

   }

   /*=====================================================
    * Now that we know which stations to keep,
    * loop over stations and write out the ones
    * to keep.
    *=====================================================*/

   /* Reset the file and station counters */
   curfile = 1;
   cursta = 0;

   do { /* Loop over each input file basename... */

     curext = 0;
     do { /* Loop over each extension... */
     infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
     if (infile == NULL) 
       goto NEXTFILE2;

     /* Loop for each station and select the ones that are
      * staOK[cursta] == 1. */
     while ((status = get_station(infile, &h, &data)) == 0) { 

       if ( staOK[cursta] == 1 ) {
          if (h_flag) {
             error = write_hydro_hdr(outfile, &h);
          }
          else {
             error = write_hydro_station(outfile, &h, &data);
          }  
          if (error) {
             fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
             exit(1);
          }
          ++nout;
       }
       else {
          ++nskip;
          if (r_flag) {
              if (h_flag) {
                error = write_hydro_hdr(rejectfile, &h);
              }
              else {
                error = write_hydro_station(rejectfile, &h, &data); 
              }
          }
       }

       /* Move to the next station */
       cursta = cursta + 1;

     }  /* end while */

     report_status(status, stderr);
     fclose(infile);

NEXTFILE2:
     ;
     
    } while (++curext < n_extents);  
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n %d stations selected.",nout);
   fprintf(stderr,"\n %d stations skipped.",nskip);
   fprintf(stderr,"\n End of hb_random.\n");
   exit(0);

} /* end main() */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n %s will select random stations from",program);
  fprintf(stderr,"\n HydroBase station files.\n");

  fprintf(stderr,"\n Usage: %s filename_root(s)", program);
  fprintf(stderr,"\n [-Ddirname] [-Eextent] [-H] [-Ooutput_file] ");
  fprintf(stderr,"\n [-Rrejected_file] [-Fnum_input_stations]");
  fprintf(stderr,"\n [-Pprobability]\n");

  fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
  fprintf(stderr,"\n          ex: -D../data/ ");
  fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
  fprintf(stderr,"\n          ex: -E.dat ");
  fprintf(stderr,"\n    -H  : output header only");  
  fprintf(stderr,"\n    -O  : specifies output file for selected stations");
  fprintf(stderr,"\n          instead of standard output.");
  fprintf(stderr,"\n          ex: -Onpac_post2003.sta");
  fprintf(stderr,"\n    -R  : specifies file for rejected stations ");  
  fprintf(stderr,"\n          ex: -Rnpac_pre2003.sta");
  fprintf(stderr,"\n    -F  : fast mode works by giving %s the number",program);
  fprintf(stderr,"\n          of input stations so no counting is needed.");
  fprintf(stderr,"\n          Fast mode is useful for bootstapping analysis");
  fprintf(stderr,"\n          based on a single initial database.");
  fprintf(stderr,"\n    -P  : specify the probability [0,1] of accepting");
  fprintf(stderr,"\n          a station, default: -P%f",DEFAULT_PROB);
  fprintf(stderr,"\n    -h  : help ... prints this message.");
  fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/
