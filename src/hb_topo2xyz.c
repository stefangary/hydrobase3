/* hb_topo2xyz.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Institution
                             April 2010
................................................................................
..................................................................  
.  Reads etopo1 gridded bathymetry 
.  and writes out ascii lon/lat/z triplets

..................................................................  
*/ 
#include <stdio.h> 
#include <stdlib.h> 
#include "hb_grids.h"
#include "hb_paths.h"


/* globally defined variables */


int  depth_is_neg; 
char *toponame;

/*  prototypes for locally defined functions */

void print_usage(char *);

           
int main (int argc, char **argv)
{ 
   int i, j, k, row, col;
   int error;
   int lon0to360;   
   short int *topo, missing, z;
   double *latvec, *lonvec;
   char *s;
   FILE *outfile;
   struct GRID_INFO h;

   if (argc < 1 ){
      print_usage(argv[0]);
   }
   
/* initialize these ... */

   h.x_min = -180.;       /* ETOPO1 ranges for globe */
   h.x_max = 180.0;
   h.y_min = -90.0;
   h.y_max = 90.0;
   lon0to360 = 0;
   depth_is_neg = 0; 
   error = 0;  
   outfile = stdout;
   toponame = BATHPATH_C;

   
/* parse the command line arguments */

   for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') {
         s = &argv[i][1]; 
         switch (*s) { 
            case 'B':                    /* get output grid bounds */
                s = &argv[i][2];
                if (*s == '/')
                      ++s;
                error = (sscanf(s,"%lf", &h.x_min) != 1);
                while (*(s++) != '/')
                    ;  
                error += (sscanf(s,"%lf", &h.x_max) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_min) != 1);
                while (*(s++) != '/')
                   ;  
                error += (sscanf(s,"%lf", &h.y_max) != 1);

                if (h.x_min > h.x_max)  {
                  fprintf(stderr,"\nW bound cannot exceed E bound.\n");
                  error = 1;
                }
                                               
                if (h.y_min > h.y_max) { 
                  fprintf(stderr,"\nS bound cannot exceed N bound.\n");
                  error = 1;
                }
                          
                break;

            case 'H' :
               toponame = BATHPATH;
	       break;
              
           case 'N' :
               depth_is_neg = 1;
               break;
	       
           case 'O':
               outfile = fopen(&argv[i][2],"w");
               if (outfile == NULL) {
                  fprintf(stderr,"\nError opening %s for output\n", &argv[i][2]);
                  exit(1);
               }
               break;
	       
           case 'T':
               toponame = &argv[i][2];
                break;
               
             case 'h': 
               print_usage(argv[0]);
               exit(0);
               
            default:  
               error = 1;
               
         } /* end switch */
             
      }
      else  {
        error = 1;
      }
      
      if (error ) { 
           fprintf(stderr,"\nError parsing command line args.\n");
           fprintf(stderr," in particular:  '%s'\n", argv[i]); 
           fprintf(stderr,"\nType %s -h for help\n", argv[0]); 
           exit(1); 
      }
   } /* end for */

   lon0to360 = h.x_min < 0 ? 0 : 1;
   
   topo = hb_get_topo(toponame, &h, &latvec, &lonvec, depth_is_neg, lon0to360, &missing);
   if (topo == NULL)
      exit(1);
 
   fprintf(stderr,"\nWriting xyz file.... \n");

   for (i = 0; i < h.ny; ++i) {
       for (j = 0; j <	h.nx; ++j) {  
           if ((z = topo[i*h.nx+j]) != missing )
 	      fprintf(outfile,"%8.3lf %8.3lf %7d\n", lonvec[j], latvec[i], z);
        }
    }
    
    fprintf(stderr,"\n End of %s.\n", argv[0]);
    exit(0);

} /* end main */


/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n%s reads gridded bathymetry values from ", program);
   fprintf(stderr,"\ntopography file  %s", BATHPATH_C);
   fprintf(stderr,"\nHigh resolution topography (1-minute) is also available (-H)");
  fprintf(stderr,"\nAn alternate file can be specified with -T ");
   fprintf(stderr,"\nif it has a compatible netcdf format. ");
   fprintf(stderr,"\nOptional bounds can be specified for the output area."); 
   fprintf(stderr,"\nLon/lat/z values are written to outfile or stdout.");
   fprintf(stderr,"\n\nUsage:  %s  [-B<w/e/s/n>] [H] [-N] [-O<outfile>][-T<global_topo_file>] ", program);
   fprintf(stderr,"\n\n   OPTIONS:"); 
   fprintf(stderr,"\n[-B] : specifies output bounds, default: -180/180/-90/90");
   fprintf(stderr,"\n[-H] : Use high resolution topography file %s ", BATHPATH); 
    fprintf(stderr,"\n[-N] : Output negative seafloor values. default: seafloor positive");
   fprintf(stderr,"\n[-O] : output filename. If not specified output goes to stdout"); 
   fprintf(stderr,"\n[-T] : specify name of netcdf topography file. Default[%s]", toponame);
   fprintf(stderr,"\n[-h] : help -- prints this message.");
   fprintf(stderr,"\n\n");
   return;
   
} /*end print_usage() */
   
