/*  hb_surfdiff2d.c  
................................................................................
                          *******  HydroBase3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             Updated to ANSI conformance July 2000
................................................................................

   Reads 2 files of gridded xyz values, and subtracts the values in the 
  second file from the first for gridpoints where both contain information
  and writes them out to the stdout device.
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

  /* boundaries for grid */

float   xmin, xmax, ymin, ymax, delta_x, delta_y;     
int     ncols, nrows;

int latfirst;

  /* prototypes for locally defined functions  */
     
void print_usage(char *);
int readprop(char *, float **);


int main(int argc, char **argv)
{
   FILE *infile, *outfile;
   short bopt,  iopt;
   char *st, *name[2];
   int   error, n, i, row, col;
   float lat, lon, pos;
   float **x, x1; 

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/* set these default values... */

   bopt  = iopt  = 0;
   latfirst = 0;
   error = 0;
   outfile = stdout;
   n = 0;


/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
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
                        break;

               case 'I':
                        iopt = 1;
                        error = (sscanf(&argv[i][2],"%f", &delta_x) == 1) ? 0 : 1;
                        delta_y = delta_x;
                        st = strchr(&argv[i][2],'/');
                        if (st != NULL) {
                          sscanf(++st,"%f", &delta_y);
                        }
                        break;
			
               case 'O':
                        outfile = fopen(&argv[i][2], "w");
			if (outfile == NULL) {
                          fprintf(stderr, "\nUnable to open %s for writing.\n", &argv[i][2]);
                          exit(1);
			
			}
                        break;

               case ':':
                        latfirst = 1;
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
       else {
          if (n < 2) 
             name[n++] = argv[i];
          
          else {
             fprintf(stderr,"\nToo many input files specified!\n");
             fprintf(stderr,"Ignoring %s\n", argv[i]);
          }
       }

   }  /* end for */

   if (!bopt || !iopt || (n < 2) ) {
       fprintf(stderr,"\nYou must specify input files, bounds, and gridspacing!\n");
       exit(1);
   }
   /* compute dimensions of matrix formed by grid  */

   nrows = (int) (ceil((double)((ymax - ymin) / delta_y)) + .0001);
   ncols = (int) (ceil((double)((xmax - xmin) / delta_x)) + .0001);
      

/*   allocate space for gridded values and initialize... */

   x = (float **) malloc(nrows * sizeof(float *));
   for (i = 0; i < nrows; ++i ) {
      x[i] = (float *) malloc(ncols * sizeof(float));
   }
   for (row = 0; row < nrows; ++row) {
      for (col = 0; col < ncols; ++col) {
          x[row][col] = -99999.0;
      }
   }

/*   read in values from second file ... */
 
   n = readprop(name[1], x);

/* now read first file and output difference ... */

   infile = fopen(name[0], "r");
   if (infile == NULL) {
       fprintf(stderr, "\nUnable to open %s\n", name[0]);
       exit(1);
   }
   fprintf(stderr," Opened %s ...\n", name[0]);

   while (fscanf(infile,"%f%f%f", &lon, &lat, &x1) != EOF) {
     if (latfirst) {
        pos = lon;
	lon = lat;
	lat = pos;
     }
     row = (int) (.0001 + (lat - ymin) / delta_y);
     col = (int) (.0001 + (lon - xmin) / delta_x);
     if ((row >= 0) && (row < nrows) && (col >= 0) && (col < ncols)) {
        if (x[row][col] > -99998.)
           fprintf(outfile,"%8.3f %8.3f %8.3f\n", lon, lat, (x1-x[row][col]));
     }
   }
   
   close(infile);
   fflush(outfile);
   fprintf(stderr, "Done.\n");
   exit(0);
}  /* end main */

/****************************************************************************/

void print_usage(char *program)
{
  fprintf(stderr,"\n %s reads 2 files of gridded xyz values, and subtracts the values ", program);
  fprintf(stderr,"\nin the second file from the first for gridpoints where both contain");
  fprintf(stderr,"\ninformation and writes them out to the stdout device.");
   fprintf(stderr,"\n\nUsage:  %s file1 file2 -Bwest/east/south/north -Ideltax/deltay  [-O<outfile>] [-:] > output_file", program);
   fprintf(stderr,"\n  -B   specifies grid bounds. Ex: -B-80/0/0/65");
   fprintf(stderr,"\n  -I   specifies grid increments. Ex: -I1/.5");
   fprintf(stderr,"\n [-O]  specify name of output file");
   fprintf(stderr,"\n [-:]  order of first 2 columns is lat/lon (default is lon/lat)");
   fprintf(stderr,"\n\n");  
   return;
}
/***********************************************************************/
int readprop(char *filename, float **x)
{
   FILE *infile;
   float  y;
   float lont ,lon, lat;
   int row, col;


   if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "\n Unable to open %s for input\n\n", filename);
        exit(1);
   }
   fprintf(stderr," Opened %s ...\n", filename);

   while (fscanf(infile,"%f%f%f", &lon, &lat, &y) != EOF) {
     if (latfirst) {
        lont = lon;
	lon = lat;
	lat = lont; 
     }
     row = (int) (.0001 + (lat - ymin) / delta_y);
     col = (int) (.0001 + (lon - xmin) / delta_x);


     /* check that it is within the bounds */

     if ((row >= 0) && (row < nrows) && (col >= 0) && (col < ncols)) 
            x[row][col] = y;
   }

   close (infile);

   return (0);

}  /* end readprop() */




