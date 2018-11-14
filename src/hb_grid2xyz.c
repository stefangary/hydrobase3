/*  hb_grid2xyz.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     updated to ANSI Dec 1999
................................................................................
_______________________________________________________________________________

   USAGE:  hb_grid2xyz <input_filename> -O<output_file_root> -P<property_list>
             [-B] [-N] [-S] 
_______________________________________________________________________________
.  Separates a HydroBase *.grid file (output by hb_gridsurf2d or hb_ncsurf2d) into individual
.  *.xyz files containing lon, lat, property triplets.   The number of observations
.  is included in the output if -N option is specified.  -S creates files
.  with error values.  -B creates blanking files (gridpts where property
.  has no observations)
.
.  The program outputs up to 4 files for each property specified: 
.     1) property info:   lon, lat, av_value [, n]
.     2) blank gridpts:   lon, lat
.     3) stddev info  :   lon, lat, stddev [, n]
.     4) error blanks:   lon, lat
.  The error information is only output if the -S option is specified.
.
.   Output files are named:     <root><property_id>.xyz
.                                  "       "       .blk
.                                  "       "       _err.xyz
.                                  "       "       _err.blk
.   where root is supplied by the -O option and property_id is the same as
.   the property_list_ids in the -P option.  
.   
.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"


/* data structures to represent matrix */

struct old_grid_pt {
       float mean;
      double var;
         int n;
 };



void  print_usage(char *);


int main (int argc, char **argv)
{
   int    row, col, n, i, j, sq;
   int    ncols, nrows;
   int    gridsize;
   int    error, index, nprops;
   int    prop_index[MAXPROP];     
   char   prop_request[MAXPROP], prop_avail[MAXPROP];                      
   int    nflag, stddev_flag;
   int    pflag, oflag, bopt;
   char   *id, field_descrip[12];
   double mean, stddev;
   float  lat, lon;
   float  xmin, ymin,xmax, ymax, xspacing, yspacing;
   FILE  *infile, *outfile[4];
   char  *fname, *root, *s;
   struct old_grid_pt *gridpt, *old_matrix[MAXPROP];
   struct new_grid_pt *newpt, *new_matrix[MAXPROP];

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/* initialize flags */

    for (i = 0; i < MAXPROP; ++i) {
        prop_request[i] = prop_avail[i] = 0;
    }
    id = (char *) calloc(6, sizeof(char));
    oflag = pflag = 0;
    nflag = stddev_flag = 0;
    bopt = 0;
    error = 0;
    
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'B':
                       bopt = 1;
                       break;
               case 'N':
                       nflag = 1;
                       break;
               case 'O':                    /* get output file root_name */
                        oflag = 1;
                        root = &argv[i][2];
                        break;

               case 'P':                    /* get list of properties */
                        pflag = 1;
                        s = &argv[i][2];
                        if (*s == '\0') {
                               print_prop_menu();
                               exit(0);
                        }
                        do {
                         if (*s == '/')
                               ++s;
                         sscanf(s,"%[^'/']", id);
                         index = get_prop_indx(id);
                         if (error = (index < 0) ? 1 : 0)
                            break;
                         prop_request[index] = '1';
                         s += strlen(id);
                        } while (*s == '/');
                        break;
               case 'E':
                       stddev_flag = 1;
                       break;

                case 'h':
                       print_usage(argv[0]);
                       exit(0);
               default:
                       error = 1;

          }    /* end switch */
	  
          if (error ) {
             print_usage(argv[0]);
             fprintf(stderr,"\nError parsing command line arg:\n");
             fprintf(stderr,"      '%s'\n", argv[i]);
             exit(1);
          }

       }  /* end if */

       else  {
          infile = fopen(argv[i],"r");
	  if (infile == NULL) {
	     fprintf(stderr,"Unable to open %s for reading.\n", argv[i]);
	     exit(1);
	  }
          fprintf(stderr, "\nOpened %s ... \n", argv[i]);
         
       }  /* end else */

   }  /* end for */

   if (!(oflag &&  pflag)) {
       fprintf(stderr,"\nYou must specify -O and -P arguments.\n");
       exit(1);
   }

/*   get row/column dimensions from infile */
 
   if (fscanf(infile,"%d%d%f%f%f%f%f%f%d", &nrows, &ncols, &xspacing, &yspacing, &xmin,
      &ymin, &xmax, &ymax, &nprops )  != 9) {
        fprintf(stderr,"\nError reading heading from infile\n\n");
        exit(1);
   }
   fprintf(stderr, "\nInput grid is %d rows by %d cols", nrows, ncols);

   gridsize = nrows * ncols;

   fprintf(stderr,"\n Properties available: ");
   for (i = 0; i < nprops; ++i) {
       fscanf(infile,"%s", id);
       prop_index[i] = get_prop_indx(id);
       prop_avail[prop_index[i]] = '1';
       fprintf(stderr," %s", id);
   }
   fprintf(stderr,"\n");

/* check that requested properties are available and allocate space for matrices */

   for (i = 0; i < MAXPROP; ++i ) {
       if (prop_request[i] && !prop_avail[i]) {
          fprintf(stderr,"\nRequested property: %s not available.", get_prop_mne(i));
          prop_request[i] = 0;
       }

       old_matrix[i] = NULL;
       if (prop_avail[i]) {
          old_matrix[i] = (struct old_grid_pt *) malloc(sizeof(struct old_grid_pt) * gridsize);
          if (old_matrix[i] == NULL) {
             fprintf(stderr,"\nUnable to allocate memory for old matrix.");
             exit(1);
          }
          for (sq= 0; sq < gridsize; ++sq) {
               gridpt = &old_matrix[i][sq];
               gridpt->mean = 0.;
               gridpt->var = 0.;
               gridpt->n = 0;
          }
       } 
   }

/*  read data from infile and store in matrices by row, col.  Multiply means
    by nobs to produce raw sums ; square the std dev before multiplying by
    (nobs - 1) and adding (sumx * sumx / n) to produce the raw sumofsquares */

   fprintf(stderr,"\nReading input file...\n");

   while ( fscanf(infile, "%f%f",  &lat, &lon) != EOF) {
     row = (int) (.0001 + (lat - ymin) / yspacing);
     col = (int) (.0001 + (lon - xmin) / xspacing);

     for (i = 0; i < nprops; ++i) {
        if (fscanf( infile,"%lf%lf%d", &mean, &stddev, &n) != 3) {
           fprintf(stderr,"\nError reading input file at ");
           fprintf(stderr,"values of (row, col) = ( %d , %d)\n", row, col);
        }

        j = prop_index[i];
        sq = row * ncols + col;       /* find index to old_matrix */
        gridpt = &old_matrix[j][sq];  /* get addr of gridpoint */
        gridpt->mean = mean;          /* store info at that addr */
        gridpt->var = stddev;
        gridpt->n = n;
     }

   }  /* end while */

/*  for each pt in the new grid, 
    Compute means and standard deviations from weighted sums   */

   fprintf(stderr,"\nWriting to output files ...\n");
   fname = (char *)malloc(500);
   free(id);

   for (i = 0; i < MAXPROP; ++i) {
     if (prop_request[i]) {

        id = get_prop_mne(i);    
        sprintf(field_descrip, " %c%d.%dlf", '%',get_field_width(i), get_field_precis(i)); 
        fprintf(stderr,"\n%s ", id);        /* print prop id on stderr */

        fname = strcpy(fname,root);
        fname = strcat(fname, id);
        fname = strncat(fname, ".xyz", 5);
        outfile[0] = fopen(fname, "w");       /* open output file for prop */
        
	if (bopt) {
           fname = strcpy(fname,root);
           fname = strcat(fname, id);
           fname = strncat(fname, ".blk", 4);
           outfile[1] = fopen(fname, "w");       /* open blanking file for prop */
	}   

        if (stddev_flag) {
           fname = strcpy(fname,root);
           fname = strcat(fname, id);
           fname = strncat(fname, "_err.xyz", 9);
           outfile[2] = fopen(fname, "w");    /* open stddev file for prop */

	   if (bopt) {
              fname = strcpy(fname,root);
              fname = strcat(fname, id);
              fname = strncat(fname, "_err.blk", 8);
              outfile[3] = fopen(fname, "w");    /* open stddev blanking file  */
	   }   
        }

        sq = -1;
        for (row=0; row < nrows; ++row) {
           fprintf(stderr,"*");
           for (col=0; col < ncols; ++col) {
              ++sq;
              lat = ymin + (row + .5) * yspacing;
              lon = xmin + (col + .5) * xspacing;

                  gridpt = &old_matrix[i][sq];

                 if ((n = gridpt->n) == 0 ) {       /* write to blanking files */
		    if (bopt) {
                       fprintf(outfile[1],"%.3f %.3f\n", lon, lat);
                       if (stddev_flag) {
                          fprintf(outfile[3],"%.3f %.3f\n", lon, lat);
                       }
		    }
                 }
                 else {                     /* write to property file */
                    fprintf(outfile[0],"%.3f %.3f ", lon, lat); 
                    fprintf(outfile[0], field_descrip, gridpt->mean);
                    if (nflag) {
                         fprintf(outfile[0]," %4d",  n);
                    }
                    fprintf(outfile[0],"\n");

                    if (stddev_flag) {         
                       if (gridpt->var >= 0) {         /* write to stddev file */
                          fprintf(outfile[2],"%.3f %.3f %.7lf", lon, lat, gridpt->var);
                          if (nflag) {
                             fprintf(outfile[2]," %4d",  n);
                           }
                          fprintf(outfile[2],"\n");
                       }
                       else  {         /* write to stddev blanking file */
		          if (bopt) 
                            fprintf(outfile[3],"%.3f %.3f\n", lon, lat);
                       }
                    }

                 } /* end else */
           } /* end for */

        } /* end for */

       fclose(outfile[0]);
       if (bopt)
          fclose(outfile[1]);
       if (stddev_flag) {         
         fclose(outfile[2]);
	 if (bopt)
            fclose(outfile[3]);
       }

     } /* end if */
   }  /* end for */

   fprintf(stderr,"\nEnd of %s.\n", argv[0]);
   exit(0);
} /* end main */




/****************************************************************************/

void print_usage(char *program)
{

  fprintf(stderr,"\nSeparates a HydroBase *.grid file ");
  fprintf(stderr,"\n(output by hb_gridsurf2d or hb_ncsurf2d) ");
  fprintf(stderr,"\ninto individual *.xyz files containing lon, lat, prop triplets.  ");
  fprintf(stderr,"\n\nThe program outputs up to 4 files for each property specified:"); 
  fprintf(stderr,"\n   1) property info:   lon, lat, prop [, n]");
  fprintf(stderr,"\n   2) blank gridpts:   lon, lat");
  fprintf(stderr,"\n   3) error info  :   lon, lat, error [, n]");
  fprintf(stderr,"\n   4) error blanks:   lon, lat");
  fprintf(stderr,"\nError information (stddev) is only output if the -E option is specified.");
  fprintf(stderr,"\nBlanking files are only output if the -B option is specified.");
  fprintf(stderr,"\nNumber of obs is included if the -N option is specified.");

  fprintf(stderr,"\n\n Output files are named:     <root><property_id>.xyz");
  fprintf(stderr,"\n                             <root><property_id>.blk");
  fprintf(stderr,"\n                             <root><property_id>_err.xyz");
  fprintf(stderr,"\n                             <root><property_id>_err.blk");
  fprintf(stderr,"\nwhere root is supplied by the -O option and property_id is the same as");
  fprintf(stderr,"\nthe property_list_ids in the -P option.");


  fprintf(stderr,"\nUsage:  %s input_filename -O<output_file_root> -P<property_list> [-B] [-E] [-N] [-h]", program);
   fprintf(stderr,"\n\n   specify input file_name (HydroBase *.grid file ");
   fprintf(stderr,"\n   output by hb_gridsurf) as first argument");
   fprintf(stderr,"\n    -O  : specifies root name for output files");
   fprintf(stderr,"\n          ex: -Opr200.1970_75.");
   fprintf(stderr,"\n    -P  : specifies properties to project onto surface");
   fprintf(stderr,"\n          Enter a list of your choices (separated by slashes)...");
   fprintf(stderr,"\n     ex: -Ppr/th9/sa/ox/de/ht");
   fprintf(stderr,"\n         -P (by itself) will print a menu of available properties");
   fprintf(stderr,"\n   [-B] : output list of blank gridpts.");
   fprintf(stderr,"\n   [-E] : output error file for each property.");
   fprintf(stderr,"\n   [-N] : output number of obs incorporated into each gridpt.");
   fprintf(stderr,"\n   [-h] : help -- prints this message");
   

   fprintf(stderr,"\n\n");  
   return;
}



