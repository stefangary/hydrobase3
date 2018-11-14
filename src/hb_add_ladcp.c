/* hb_add_ladcp.c
................................................................................
                         *******  HydroBase3 *******
................................................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Inst
                             Nov 2006 
			     updated to HydroBase3 Feb 2013
................................................................................
................................................................................
.  Reads in files of LADCP data and adds properties ve and vn to a HydroBase station file
.  The number of header lines in the LADCP file can be specified with -H
.  (default is 7 with station # in first line. 
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "hydrobase.h"

#define    PRINT_MSG 1

/* global variables store default values */

int nhdrlines, sta_id_line, field_no;
  /* prototypes for functions declared locally */
  
void print_usage(char *);

main (int argc, char **argv)
{
   int     i, j, npts;
   int     status, error, cflag;
   int     station_id;
   int  *tmpi;
   double de_in[6000], ve_in[6000], vn_in[6000];
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
   char   *st, buffer[1000];
   FILE   *ladcp_file;
   FILE   *infile, *outfile;


/*  set these default values */

    infile = stdin;
    outfile = stdout;
    ladcp_file = NULL;
    nhdrlines = 7;
    sta_id_line = 1;
    field_no = 3;
    error = 0;
    cflag = 0;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'C' :            /* LADCP velocities are in cm/sec */
                        cflag = 1;
                        break;
               case 'H' :
                       error = (sscanf(&argv[i][2], "%d", &nhdrlines) == 1) ? 0 : 1;
		       if (!error) {
		          if ( (st = strchr(&argv[i][2], '/') ) != NULL) {
		             ++st;
			      error = (sscanf(st, "%d", &sta_id_line) == 1) ? 0 : 1;
			   }
			}   
                        break;
               case 'I' :
                        infile = open_hydro_file("", &argv[i][2], "", 1);
			if (infile == NULL) {
                           fprintf(stderr,"\nUnable to open %s \n", &argv[i][2]);
			   exit(1);
			}   
                        break;
               case 'L' :
                        ladcp_file = fopen(&argv[i][2], "r");
                        break;
               case 'O' :
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
			if (outfile == NULL) {
                           fprintf(stderr,"\nUnable to open %s \n", &argv[i][2]);
			   exit(1);
			}   
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
       else {
          ++error;
       } 
       if (error) { 
            fprintf(stderr,"\nError parsing command line");
            fprintf(stderr,"\n in particular: %s\n", argv[i]);
            exit(1);
       }
   }  /* end for */

   if (ladcp_file == NULL) {
            fprintf(stderr,"\nNo LADCP file specified.");
	    exit(1);
   }
   
   if (infile == STDIN) {
            fprintf(stderr,"\nExpecting HydroBase file from STDIN ");
   }
   
 /* read ladcp header lines*/

   fprintf(stderr,"LADCP file has %d header lines...\n", nhdrlines);
   if (sta_id_line > nhdrlines) {
	     fprintf(stderr, "Station ID and header lines mismatch:  nhdrlines: %2d, station id on line: %d\n", nhdrlines, sta_id_line );
	     exit(1);
   }
   
     --sta_id_line;
     for (i = 0; i < nhdrlines; ++i) {
         if (fscanf(ladcp_file, "%[^\n]", buffer) != 1) {
	     fprintf(stderr, "Error reading ladcp file at line %1d\n", i);
	     exit(1);
	 }
	 
	 error = getc(ladcp_file);
	 if (i == sta_id_line) {

           st = &buffer[0];
           while ( ! (*st >= '0' && *st <= '9') ) {
              if (*st == '\0')  {
	           fprintf(stderr, "Unable to find station id on line %1d\n", i);
	           exit(1);
	       }
	       ++st;
           }
           sscanf(st, "%d", &station_id);
	   fprintf(stderr,"LADCP station id = %2d\n", station_id); 
        }
     }
     
     /* read data lines and store de, ve, vn */
     npts = 0;
     while ( fscanf(ladcp_file, "%[^\n]", buffer) == 1)  {
	 error = getc(ladcp_file);
         if ( sscanf(buffer, "%lf %lf %lf", &de_in[npts], &ve_in[npts], &vn_in[npts]) != 3) {
	           fprintf(stderr, "Error reading LADCP data at line %d\n", i);
	           fprintf(stderr, "%s\n", buffer);
	           exit(1);
	 }
	 if (cflag) {                /* convert units from cm/sec to m/sec */
	    ve_in[npts] /= 100.0;
	    vn_in[npts] /= 100.0;
	 }
	 ++npts;
	 ++i;
     }  /* end while */
     
     
     fprintf(stderr,"Read in %d data levels from LADCP file\n", npts);
     fclose(ladcp_file);   
     
     if (de_in[0] <= 10.)   /* adjust shallowest depth to 0 if it is in upper 10 meters */
        de_in[0] = 0.0;  

/* initialize these... */

   for (i = 0; i < MAXPROP; ++i)
      data.observ[i] = (double *) NULL;
   h.prop_id = (int *) NULL;
   
     while ((status = get_station(infile, &h, &data)) == 0) { 
     
       if (station_id == h.station) {
           free_and_alloc(&data.observ[(int)VE], h.nobs);
           free_and_alloc(&data.observ[(int)VN], h.nobs);
	   for (j = 0; j < h.nobs; ++j) {
	      data.observ[(int)VE][j] = hb_linterp(data.observ[(int)DE][j], de_in, ve_in, npts);
	      data.observ[(int)VN][j] = hb_linterp(data.observ[(int)DE][j], de_in, vn_in, npts);
	      if (data.observ[(int)VN][j] < -9998.0)   /* make flagged values more obvious */
	          data.observ[(int)VN][j] = -99.99;
	      if (data.observ[(int)VE][j] < -9998.0)
	          data.observ[(int)VE][j] = -99.99;
	   }
	   
	   h.nprops +=2;
	   tmpi = h.prop_id;
	   h.prop_id = (int *) calloc(h.nprops, sizeof(int));
	   for (j = 0; j < data.nprops; ++j) 
	       h.prop_id[j] = tmpi[j];
	   h.prop_id[j] = (int) VE;
	   h.prop_id[++j] = (int) VN;    
	   data.nprops += 2;
           free(tmpi);
       }
     
      error = write_hydro_station(outfile, &h, &data);
           if (error) {
             fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
             exit(1);
      }
      
       
       
     }  /* end while */

     report_status(status, stderr);
     fclose(infile);
     fclose(outfile);

   fprintf(stderr,"\n\nEnd of %s.\n\n", argv[0]);
   exit(0);

} /* end main() */

/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s updates a Hydrobase station file with LADCP velocity (east and north)", program);
   fprintf(stderr,"\n The LADCP file contains a single profile.  Only one HydroBase profile, therefore, ");
   fprintf(stderr,"\nwill be updated, but any other stations will simply be copied to the output file.");
   
   fprintf(stderr,"\n\nUsage:  %s  ", program);
   fprintf(stderr,"  -L<ladcp_file>  [-C] [-H<#of_hdr_lines>/<line_#_of staid>]  [-I<HydroBase input file> [-O<outfile>]");
   fprintf(stderr,"\n");
   fprintf(stderr,"\n L : name of LADCP file" );
   fprintf(stderr,"\tOPTIONS:");
   fprintf(stderr,"\n[-C] : LADCP velocity units are cm/sec (default is m/sec)");
   fprintf(stderr,"\n[-H] : defines # of header lines in LADCP file and line containing station id");
   fprintf(stderr,"\n         Defaults:  %1d/%1d", nhdrlines, sta_id_line);
   fprintf(stderr,"\n[-I] : input name of HydroBase station file to be updated ");
   fprintf(stderr,"\n[-O] : name of output file. Default is STDOUT");
   fprintf(stderr,"\n[-h] : help -- prints this message");
   fprintf(stderr,"\n\n");  
   return;
}
