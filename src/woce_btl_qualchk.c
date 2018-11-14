/*  woce_btl_qualchk.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Institution
................................................................................
.   Reads Woce hydrographic data files
.   Creates xy-files for checking salinity and oxygen versus pressure 
     and theta vars.
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"

#define   MISSING   -9.0   /* missing value flag */
#define   BUFSIZE   1600   /* buffer for a line of data */
#define   DIR    ""
#define   EXTENT ""



  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
char *buffer;
char *expocode;
char *outfile_root;
short *prop_avail; 
int *column;
int *qual_offset;
int cast;
int iso2;
int qual_word_col;
int col_sta, col_cast, col_type, col_date;
int col_code, col_lat, col_lon, col_pdr, col_nobs;
int dcol_sta, dcol_cast;
int sq2, sq3, sq4, sq9, sq6;
int oq2, oq3, oq4, oq9, oq6;
int t68flag;
FILE *qual2_sapr;
FILE *qual3_sapr;
FILE *qual4_sapr;
FILE *qual9_sapr;
FILE *qual2_sath;
FILE *qual3_sath;
FILE *qual4_sath;
FILE *qual9_sath;
FILE *qual2_deepsath;
FILE *qual3_deepsath;
FILE *qual4_deepsath;
FILE *qual9_deepsath;
FILE *qual2_o2pr;
FILE *qual3_o2pr;
FILE *qual4_o2pr;
FILE *qual9_o2pr;
FILE *qual2_o2th;
FILE *qual3_o2th;
FILE *qual4_o2th;
FILE *qual9_o2th;
FILE *qual2_deepo2th;
FILE *qual3_deepo2th;
FILE *qual4_deepo2th;
FILE *qual9_deepo2th;


/* prototypes for locally defined functions */

void print_usage(char *);
void open_output_files(int);
void close_output_files();
FILE *openfile(char *, char *, char *);
void parse_headers(FILE *, FILE *);
int find_next_station(FILE *) ;
int readdata(FILE *);

int main ( int argc, char **argv)
{
   int error, nobs, nprops, npout;
   int  i, j, curfile = 1, nfiles = 0;
   short sflag, qflag;
   char *dir, *extent, *st;
   int  npts, nsalt, nox;
   FILE *infile, *sumfile, *openfile();
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    sflag = 0;
    sq2 = sq3 = sq4 = sq6= sq9 = 0;
    oq2 = oq3 = oq4 = oq6= oq9 = 0;
    iso2 = 1;
	
    expocode = (char *) calloc((size_t)15, sizeof(char));
    prop_avail = (short *) calloc((size_t)MAXPROP, sizeof(short));
    column = (int *) calloc((size_t)MAXPROP, sizeof(int));
    qual_offset = (int *) calloc((size_t)MAXPROP, sizeof(int));

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
               case 'R':                    /* get rootname for 24 standard output files  */
	                outfile_root = &argv[i][2];
                        break;

               case 'S':                    /* get summary file  */
                        sflag = 1;
                        sumfile = fopen(&argv[i][2], "r");
                        if (sumfile == NULL) {
                           fprintf(stderr,"\nUnable to open .sum file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;
               case 'X':                    /* oxygen in ml/l  */
	            iso2 = 0;
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

   if (! (nfiles && sflag )) {
       fprintf(stderr,"\nYou must specify an input file and a summary file.\n");
       exit(1);
   }
      

  infile = openfile(dir, argv[curfile], extent);
  if (infile == NULL) 
        exit(1);

   /* parse the header info from datafile and sumfile*/
      
  parse_headers(infile, sumfile);

  /* loop for each station */
  
  fprintf(stderr,"\nSearching for station ");
  
    
  while ((nobs = find_next_station(sumfile)) >= 0) {
  
      open_output_files(hdr.station);
       nobs = readdata(infile);
       if (nobs == EOF) {
          fprintf(stderr, "\n Station %d not found in data file.\n ", hdr.station);
       }
       close_output_files();
     
   }  /* end while */

   nsalt = sq2 + sq3 + sq4 + sq6 + sq9;
   nox = oq2 + oq3 + oq4 + oq6 + oq9;
   fprintf(stderr,"\n\nTally of Quality Codes read in:");
   fprintf(stderr,"\n salinity '2' : %3d  [%.1lf%%] ", sq2, (double) sq2 / (double)nsalt * 100.0);
   fprintf(stderr,"\n salinity '3' : %3d  [%.1lf%%] ", sq3,  (double) sq3 /(double) nsalt * 100.0);
   fprintf(stderr,"\n salinity '4' : %3d  [%.1lf%%] ", sq4,  (double) sq4 / (double)nsalt * 100.0);
   fprintf(stderr,"\n salinity '6' : %3d  [%.1lf%%] ", sq6,  (double) sq6 / (double)nsalt * 100.0);
   fprintf(stderr,"\n salinity '9' : %3d  [%.1lf%%] ", sq9,  (double) sq9 / (double)nsalt * 100.0);
   fprintf(stderr,"\n Total salt values: %7d  ",  nsalt );
   fprintf(stderr,"\n oxygen '2' : %3d  [%.1lf%%] ", oq2, (double) oq2 / (double)nox * 100.0);
   fprintf(stderr,"\n oxygen '3' : %3d  [%.1lf%%] ", oq3, (double) oq3 / (double)nox * 100.0);
   fprintf(stderr,"\n oxygen '4' : %3d  [%.1lf%%] ", oq4, (double) oq4 / (double)nox * 100.0);
   fprintf(stderr,"\n oxygen '6' : %3d  [%.1lf%%] ", oq6, (double) oq6 / (double)nox * 100.0);
   fprintf(stderr,"\n oxygen '9' : %3d  [%.1lf%%] ", oq9, (double) oq9 / (double)nox * 100.0);
   fprintf(stderr,"\n Total O2 values: %7d \n\n ",  nox );
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -Ssumfile [-D<dirname>] [-E<file_extent>] [-Ooutfile] [-T] [-X] [-h]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument");
   fprintf(stderr,"\n -S  : summary info file. ex:  316N151.sum ");
   fprintf(stderr,"\n -R : root name for output files");
   fprintf(stderr,"\n         suffixes for each quality code and property-property plot will be appended:  ");
   fprintf(stderr,"\n   e.g.  _qual2.sapr    _qual2.sath  _qual2.deepsath  _qual2.o2pr");
   
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n[-D] : dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-X] : oxygens units are ml/l (default is umole/kg");  
   fprintf(stderr,"\n-h : help...... prints this message.");

   fprintf(stderr,"\n\n");  
   return;
}
/*****************************************************************************/
 void open_output_files(int sta_id)
 {
   char fname[10000];
   
    fprintf(stderr,"\nOpening output files....");  
    
    
    
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.sapr");
   qual2_sapr = fopen(fname,"w" );
   if (qual2_sapr == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.sapr");
   qual3_sapr = fopen(fname,"w" );
   if (qual3_sapr == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
     
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.sapr");
   qual4_sapr = fopen(fname,"w" );
   if (qual4_sapr == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.sapr");
   qual9_sapr = fopen(fname,"w" );
   if (qual9_sapr == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }

   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.sath");
   qual2_sath = fopen(fname,"w" );
   if (qual2_sath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.sath");
   qual3_sath = fopen(fname,"w" );
   if (qual3_sath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
     
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.sath");
   qual4_sath = fopen(fname,"w" );
   if (qual4_sath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.sath");
   qual9_sath = fopen(fname,"w" );
   if (qual9_sath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }

   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.deepsath");
   qual2_deepsath = fopen(fname,"w" );
   if (qual2_deepsath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.deepsath");
   qual3_deepsath = fopen(fname,"w" );
   if (qual3_deepsath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
     
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.deepsath");
   qual4_deepsath = fopen(fname,"w" );
   if (qual4_deepsath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.deepsath");
   qual9_deepsath = fopen(fname,"w" );
   if (qual9_deepsath == NULL) {
      fprintf(stderr,"\nUnable to open %s\n", fname);
      exit(1);
   }
   if (iso2) {
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.o2pr");
      qual2_o2pr = fopen(fname,"w" );
      if (qual2_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.o2pr");
      qual3_o2pr = fopen(fname,"w" );
      if (qual3_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.o2pr");
      qual4_o2pr = fopen(fname,"w" );
      if (qual4_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.o2pr");
      qual9_o2pr = fopen(fname,"w" );
      if (qual9_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }

      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.o2th");
      qual2_o2th = fopen(fname,"w" );
      if (qual2_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.o2th");
      qual3_o2th = fopen(fname,"w" );
      if (qual3_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.o2th");
      qual4_o2th = fopen(fname,"w" );
      if (qual4_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.o2th");
      qual9_o2th = fopen(fname,"w" );
      if (qual9_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }

      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.deepo2th");
      qual2_deepo2th = fopen(fname,"w" );
      if (qual2_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.deepo2th");
      qual3_deepo2th = fopen(fname,"w" );
      if (qual3_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.deepo2th");
      qual4_deepo2th = fopen(fname,"w" );
      if (qual4_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.deepo2th");
      qual9_deepo2th = fopen(fname,"w" );
      if (qual9_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
    }
    
    
    else {
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.oxpr");
      qual2_o2pr = fopen(fname,"w" );
      if (qual2_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.oxpr");
      qual3_o2pr = fopen(fname,"w" );
      if (qual3_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.oxpr");
      qual4_o2pr = fopen(fname,"w" );
      if (qual4_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.oxpr");
      qual9_o2pr = fopen(fname,"w" );
      if (qual9_o2pr == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }

      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.oxth");
      qual2_o2th = fopen(fname,"w" );
      if (qual2_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.oxth");
      qual3_o2th = fopen(fname,"w" );
      if (qual3_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.oxth");
      qual4_o2th = fopen(fname,"w" );
      if (qual4_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.oxth");
      qual9_o2th = fopen(fname,"w" );
      if (qual9_o2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }

      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual2.deepoxth");
      qual2_deepo2th = fopen(fname,"w" );
      if (qual2_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual3.deepoxth");
      qual3_deepo2th = fopen(fname,"w" );
      if (qual3_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
     
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual4.deepoxth");
      qual4_deepo2th = fopen(fname,"w" );
      if (qual4_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
      sprintf(fname,"%s.%04d%s", outfile_root, sta_id, ".qual9.deepoxth");
      qual9_deepo2th = fopen(fname,"w" );
      if (qual9_deepo2th == NULL) {
         fprintf(stderr,"\nUnable to open %s\n", fname);
         exit(1);
      }
    
    } 
  return;
 }  
/*****************************************************************************/
void close_output_files()
{
fclose(qual2_sapr);
fclose(qual3_sapr);
fclose(qual4_sapr);
fclose(qual9_sapr);
fclose(qual2_sath);
fclose(qual3_sath);
fclose(qual4_sath);
fclose(qual9_sath);
fclose(qual2_deepsath);
fclose(qual3_deepsath);
fclose(qual4_deepsath);
fclose(qual9_deepsath);
fclose(qual2_o2pr);
fclose(qual3_o2pr);
fclose(qual4_o2pr);
fclose(qual9_o2pr);
fclose(qual2_o2th);
fclose(qual3_o2th);
fclose(qual4_o2th);
fclose(qual9_o2th);
fclose(qual2_deepo2th);
fclose(qual3_deepo2th);
fclose(qual4_deepo2th);
fclose(qual9_deepo2th);
return;
}
/*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)
{
   char st[80];
   int i;
   FILE *infile;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   infile = fopen(st,"r");
   if (infile == NULL)
         fprintf(stderr,"\n Unable to open %s \n", st);
   else
         fprintf(stderr,"\n   opened %s ...\n", st);
   
   return(infile);
   
}  /* end openfile() */
/*****************************************************************************/
void parse_headers(FILE *fptr, FILE *sumfile)

   /* reads the header info from an already opened .sea data file and sum file. 
      The info is used to fill in appropriate values of global variables.  */
{
   int n, i, error;
   char *line, *line2, *line3, *xx, x[9];

   line = (char *) calloc(1600, sizeof(char));
   line2 = (char *) calloc(1600, sizeof(char));
   line3 = (char *) calloc(1600, sizeof(char));
   
/* Read expocode from bottle data file... */
/* line 1 */
   xx = NULL;
   while (xx == NULL) {

     if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 1st header line. \n");
        exit(1);   
      } 
      xx = strstr(line, "EXPOCODE"); 
   } /* end while */
     
   if (sscanf(&xx[9],"%s", expocode) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read expocode \n");
        exit(1);
   }
   if (sscanf(&expocode[2],"%2s", hdr.ship) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read shipcode \n");
        exit(1);
   }
   if (sscanf(&expocode[0],"%2s", &hdr.country) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read country # \n");
        exit(1);
   } 
   sscanf(&expocode[4],"%[^/]", x);
   if (sscanf(x,"%d", &hdr.cruise) != 1) {
        hdr.cruise = 9999;
        fprintf(stderr,"\n cruise # not digital: %s ", x);
        fprintf(stderr,"\n assigning it the value 9999 and continuing \n");
   } 
	if (hdr.cruise > 99999) {
	   hdr.cruise = hdr.cruise % 10000;
	   fprintf(stderr,"\n cruise number too large, using %d", hdr.cruise);
	}
   
   error = getc(fptr);

/* PROPERTIES/PARAMETERS */

   i = 2;    
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #%d. \n", i);
        exit(1);   
   }         
   error = (getc(fptr) != LF);
      
/* UNITS.. */

   ++i;
   
   if (fscanf(fptr,"%[^\n]", line2) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #%d. \n", i);
        exit(1);   
   }         
   error = (getc(fptr) != LF);
   
/* the next line indicates which properties have quality bytes */
   
   if (fscanf(fptr,"%[^\n]", line3) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #%d. \n", i);
        exit(1);   
   }         
   error = getc(fptr);
   
/* bottle file is now positioned at start of data records. Read in the first
   line of data to a global buffer... */
   
   buffer = (char *) malloc((size_t)BUFSIZE * sizeof(char));
   if (fscanf(fptr,"%[^\n]", buffer) != 1) {
        fprintf(stderr,"\n Error attempting to read first data line. \n");
        exit(1);   
   }
   error = getc(fptr);
   
/* substitute ordinal number for asterisks where quality byte is indicated */

   xx = line3;
   n = 0;
   while (xx = strstr(xx, "*******") ) {
      sprintf(x,"%7d", n);
      for (i = 0; i < 7; ++i) {
         xx[i] = x[i];
      }
      ++n;
   }
   
/* assimilate lines 2,3 & 4 to determine column positions and quality byte pos */   
   
   xx = strstr(line, "QUALT1");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find QUALT1 in datafile header \n");
           exit(1);
   }
   qual_word_col = (int) (xx - line - (n - 6));
   
   xx = strstr(line, "STNNBR");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find STNNBR in datafile header \n");
           exit(1);
   }
   dcol_sta = (int) (xx - line);
   
   xx = strstr(line, "CASTNO");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CASTNO in datafile header \n");
           exit(1);
   }
   dcol_cast = (int) (xx - line);
   
    xx = strstr(line, "CTDPRS");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CTDPRS in datafile header \n");
           exit(1);
   }
   prop_avail[(int)PR] = 1;
   column[(int)PR] = (int) (xx - line - 1);
   qual_offset[(int)PR] = 0;
  
   xx = strstr(line, "CTDTMP");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CTDTMP in datafile header \n");
           exit(1);
   }
   prop_avail[(int)T90] = 1;
   column[(int)T90] = (int) (xx - line - 1);
   qual_offset[(int)T90] = 0;
   
   xx = strstr(line, "THETA");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find THETA in datafile header \n");
           exit(1);
   }
   prop_avail[(int)TH9] = 1;
   column[(int)TH9] = (int) (xx - line - 2);
   qual_offset[(int)TH9] = 0;
  
   xx = strstr(line, "CTDSAL");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CTDSAL in datafile header \n");
           exit(1);
   }
   prop_avail[(int)S2] = 1;     /* use S2 for CTD Salinity in place of bad bottle salts */
   column[(int)S2] = (int) (xx - line - 1);
   sscanf(&line3[column[(int)S2]],"%d", &qual_offset[(int)S2]);
   
   
   xx = strstr(line, "SALNTY");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find SALNTY in datafile header \n");
           exit(1);
   }
   prop_avail[(int)SA] = 1;
   column[(int)SA] = (int) (xx - line - 1);
   sscanf(&line3[column[(int)SA]],"%d", &qual_offset[(int)SA]);
   
   xx = strstr(line, "OXYGEN");
   if (xx != NULL) {
      prop_avail[(int)O2] = 1;
      column[(int)O2] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)O2]],"%d", &qual_offset[(int)O2]);
      
       xx = strstr(line, "CTDOXY");
       if (xx != NULL) {
          prop_avail[(int)VA] = 1;     /* use VA for CTD Oxygen in place of bad bottle oxygens */
          column[(int)VA] = (int) (xx - line - 1);
          sscanf(&line3[column[(int)VA]],"%d", &qual_offset[(int)VA]);
       }
       
        /* check units */
      line2[column[(int)O2] + 7] = '\0';  /* temporary EOString marker */  
      xx = strstr(&line2[column[(int)O2]], "KG");
      if (xx == NULL) {
        xx = strstr(&line2[column[(int)O2]], "ML");
        if (xx == NULL) 
          fprintf(stderr, "\nCan't tell if OXYGEN units are umol/kg or ml/l.");
        
	else 
            fprintf(stderr, "\nSEA file says oxygen units are ml/l .");
	
     }
     else 
         fprintf(stderr, "\nSEA file says oxygen units are umol/kg .");
       
   }
   
   
/*************    SUMMARY FILE PART *************************** */   
/* Now get column info from .sum file:      */
   
   for (i = 1; i <= 3; ++i) {
      if (fscanf(sumfile,"%[^\n]", line3) != 1) {
        fprintf(stderr,"\nError reading summary file header line #%d. \n", i);
        exit(1);   
      }         
      error = getc(sumfile);
   }
    
   xx = strstr(line3, "STNNBR");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find STNNBR in summary header \n");
           exit(1);
   }
   col_sta = (int) (xx - line3);
   
   xx = strstr(line3, "CAST");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CAST in summary header \n");
           exit(1);
   }
   col_cast = (int) (xx - line3);  /*find offset of this string in line 3 */   
   
   
   xx = strstr(line3, "TYPE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CAST TYPE in summary header \n");
           exit(1);
   }
   col_type = (int) (xx - line3);  /*find offset of this string in line 3 */
   
   
   xx = strstr(line3, "DATE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find DATE in summary header \n");
           exit(1);
   }
   col_date = (int) (xx - line3);  /*find offset of this string in line 3 */
   
   xx = strstr(line3, " CODE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CODE in summary header \n");
           exit(1);
   }
   col_code = (int) (xx - line3 + 1);  /*find offset of this string in line 3 */
   
   xx = strstr(line3, "LATITUDE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find LATITUDE in summary header \n");
           exit(1);
   }
   col_lat = (int) (xx - line3);  /*find offset of this string in line 3 */
   
   xx = strstr(line3, "LONGI");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find LONGITUDE in summary header \n");
           exit(1);
   }
   col_lon = (int) (xx - line3 - 1);  /*find offset of this string in line 3 */
   
   xx = strstr(line3, "CDEPTH");
   col_pdr = 0;
   if (xx != NULL) {
      col_pdr = (int) (xx - line3);  /*find offset of this string in line 3 */
   }
   
   xx = strstr(line3, "BOTTLES");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find NO. OF BOTTLES in summary header \n");
           col_nobs = 0;
   }
   else
      col_nobs = (int) (xx - line3);  /*find offset of this string in line 3 */
   
   /* Move past dashed line */
   
   error = fscanf(sumfile,"%[^\n]", line);      
   error = getc(sumfile);
  
   free(line);
   free(line2);
   free(line3); 
   return;
}  /* end parse_headers() */

/*****************************************************************************/
  
  
int find_next_station(FILE *sumfile) 

/* Finds next bottle station in the summary file and returns the number of
  bottles.  Returns 0 at end of summary file.*/
{ 
   char *line, *line2, found, ct[4], code[3], hem;
   char *st;
   int i, n, error, sta, c;
   long int offset;
   float deg, min;
   
   line = (char *) calloc((size_t)1600, sizeof(char));
   found = 0;
   n = 0;
        
   while (!found) {
   
      if (fscanf(sumfile,"%[^\n]", line) != 1) {
        return (EOF);                                    /* end of file */
      }
      error = getc(sumfile);                           /* move past LF */
   
      if (sscanf(&line[col_sta],"%d", &sta) != 1) {      /* station # */
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(1);   
      }
      if (sscanf(&line[col_type],"%s", &ct[0]) != 1) {
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(1);   
      }
      if (sscanf(&line[col_cast],"%d", &c) != 1) {
           c = 0;   
      }
      if (sscanf(&line[col_code],"%s", &code[0]) != 1) {
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(1);   
      }
      
 /*  Is this a bottle station? */ 
      
      if ( (strncmp(ct,"ROS",3) == 0) || (strncmp(ct,"LVS",3) == 0) ||
           (strncmp(ct,"BOT",3) == 0) ) {
           
        /* Is it a new station? */  
              
        if ((sta != hdr.station) || (c > cast)) {
        
           /* We prefer the bottom info, but we'll settle for beginning
              of cast if that's all that is available.... */
	      
           hdr.pdr = 0;  
              
           if (strncmp(code,"BE",2) == 0) {
              if (col_pdr > 0) {    /* try to get pdr depth  */ 
		    st = &line[col_pdr];
		    i = 0;
		    while ((++i < 6) && (*st == ' '))  /* check for an entry */
		        ++st;
                    if (sscanf(&line[col_pdr],"%d", &hdr.pdr) != 1) 
                    hdr.pdr = 0;
	      }
              line2 = line;   /* save this, but check next line for BO code */
              offset = ftell(sumfile);
              line = (char *) calloc((size_t)1600, sizeof(char));
              
              if (fscanf(sumfile,"%[^\n]", line) != 1) {   /* end of file */
                 line = line2;                
                 code[1] = 'O';               /* set code to BO and continue */
              }
              else {
                 sscanf(&line[col_code],"%s", &code[0]); 
                 if  (strncmp(code,"BO",2) == 0) {  /* found bottom of cast */
                    free(line2);
                    error = getc(sumfile);          /* move past LF */
                 }
                 else {          /* bottom not available,revert to saved line */
                    free(line);
                    line = line2;
                    code[0] = 'B';     /* set code to BO and continue to next section */
                    code[1] = 'O';       
                    fseek(sumfile, offset, SEEK_SET); /* reset file ptr */
                 }
              }
           }  /* end if code == BE */
           
           if ( (strncmp(code,"BO",2) == 0) || (strncmp(code,"MR",2) == 0)){
           
              found = 1;
              hdr.station = sta;
              cast = c;
	      if (hdr.pdr == 0) {  /* haven't found a bottom depth yet */
                 if (col_pdr > 0) {    /* try to get pdr depth  */ 
		    st = &line[col_pdr];
		    i = 0;
		    while ((++i < 6) && (*st == ' '))  /* check for an entry */
		        ++st;
                    if (sscanf(&line[col_pdr],"%d", &hdr.pdr) != 1) 
                    hdr.pdr = 0;
		 }
	      }
              
              if (sscanf(&line[col_lat],"%f %f %c", &deg, &min, &hem) != 3) {
                 fprintf(stderr,"\n Error attempting to read latitude \n");
                 exit(1);   
              }
              hdr.lat = deg + min / 60.;
              if (hem == 'S')
                 hdr.lat = -hdr.lat;
                 
              if (sscanf(&line[col_lon],"%f %f %c", &deg, &min, &hem) != 3) {
                 fprintf(stderr,"\n Error attempting to read longitude \n");
                 exit(1);   
              }
              hdr.lon = deg + min / 60.;
              if (hem == 'W')
                 hdr.lon = -hdr.lon;
              
              
              if (sscanf(&line[col_date],"%2d%2d%2d",&hdr.month,&hdr.day,
               &hdr.year) != 3) {
                fprintf(stderr,"\n error in sscanf attempt to read date \n");
                exit(1);
              } 
              if (hdr.year < 1000)
                 hdr.year += 1900;

              if (hdr.year < 1950)  /* correct for year 2000+ data */
	         hdr.year += 100;   
 
              n = 50;  /* a nominal value that should be large enough */
	      
             
            } /* end if code == BO || MR */
        } /* end if */
      } /* end if */   
        
   } /* end while */
          
   free(line);
   return(n);

}  /* end find_next_station() */
/*****************************************************************************/
int readdata(FILE *fptr)

/* Reads bottle data and quality codes for an entire station and returns the 
   number of observations which have a minimum of pr, te, sa measurements. 
   Returns EOF if end of file and no station match is found. Returns 0 if no
   valid obs or if it cannot find the station and cast.  It is assumed that
   the station/cast info is in ascending order in the data file!  */
 
{
   int i, n, error, scanOK, sta, ca, indx, flag;
   char *line, quality[100];
   char squal, tqual, oqual;
   double salt, pressure, oxygen, theta;
   
   fprintf(stderr," %d/%d ", hdr.station, cast);
   
   line = buffer;
   error = sscanf(&line[dcol_sta],"%d", &sta);
   error = sscanf(&line[dcol_cast],"%d", &ca);
   
   flag = 0;
   n = 0;
   
   while (sta < hdr.station) {
     buffer = (char *) malloc((size_t)BUFSIZE * sizeof(char));
     if (buffer == NULL) {
       fprintf(stderr, "\n Unable to malloc memory for buffer.");
       exit(10);
     }
     if (fscanf(fptr, "%[^\n]", buffer) != 1 ) {
       fprintf(stderr, "\n Unable to find station #%d in data file.", hdr.station);
       fprintf(stderr, "\n EOF was reached.\n");
       return(EOF);
     }
     error = getc(fptr) != LF;
     free(line);
     line = buffer;
     error = sscanf(&line[dcol_sta],"%d", &sta);
     error = sscanf(&line[dcol_cast],"%d", &ca);
   }
   
   if (sta > hdr.station) {
      fprintf(stderr, "\n Unable to find station #%d cast #%d in data file.", hdr.station, cast);
      return(0);
   }
   
   while ((ca < cast) && (sta == hdr.station)) {
     buffer = (char *) malloc((size_t)BUFSIZE * sizeof(char));
     if (buffer == NULL) {
       fprintf(stderr, "\n Unable to malloc memory for buffer.");
       exit(10);
     }
     if ( fscanf(fptr, "%[^\n]", buffer) != 1) {
       fprintf(stderr, "\n Unable to find cast #%d of station #%d in data file.",
            cast, hdr.station);
       fprintf(stderr, "\n EOF was reached.\n");
       return(EOF);
     }
     error = getc(fptr) != LF;
     free(line);
     line = buffer;
     error = sscanf(&line[dcol_sta],"%d", &sta);
     error = sscanf(&line[dcol_cast],"%d", &ca);
   }
 
   /* sta matches hdr.station */         

   while ((sta == hdr.station) && (ca == cast)) {      
      flag = 1;
      error = sscanf(&line[qual_word_col], "%s", quality);
      
      squal = quality[qual_offset[(int)SA]];
      tqual = quality[qual_offset[(int)T90]];

      error = sscanf(&line[column[(int)PR]], "%lf", &pressure);
      error = sscanf(&line[column[(int)TH9]], "%lf", &theta);
      error = sscanf(&line[column[(int)SA]], "%lf", &salt);
      switch (squal) {
          case '1': 
          case '2': 
	     ++sq2;
	     fprintf(qual2_sapr,"%8.4lf %8.1lf\n", salt, pressure);
	     fprintf(qual2_sath,"%8.4lf %8.4lf\n", salt, theta);
	     if (pressure >= 1000) 
	        fprintf(qual2_deepsath,"%8.4lf %8.4lf\n", salt, theta);
	     break;
	  case '3':
	     ++sq3;
	     fprintf(qual3_sapr,"%8.4lf %8.1lf\n", salt, pressure);
	     fprintf(qual3_sath,"%8.4lf %8.4lf\n", salt, theta);
	     if (pressure >= 1000) 
	        fprintf(qual3_deepsath,"%8.4lf %8.4lf\n", salt, theta);
	  
	     break;
	  case '4':
	     ++sq4;
	     fprintf(qual4_sapr,"%8.4lf %8.1lf\n", salt, pressure);
	     fprintf(qual4_sath,"%8.4lf %8.4lf\n", salt, theta);
	     if (pressure >= 1000) 
	        fprintf(qual4_deepsath,"%8.4lf %8.4lf\n", salt, theta);
			  
	     break;
	  case '6':
	     ++sq6;
	      break;
	  case '9':
	     /* get CTDSAL */
	     ++sq9;
	     if (prop_avail[(int)S2] )  {
                  error = sscanf(&line[column[(int)S2]], "%lf", &salt);
	          fprintf(qual9_sapr,"%8.4lf %8.1lf\n", salt, pressure);
	          fprintf(qual9_sath,"%8.4lf %8.4lf\n", salt, theta);
	          if (pressure >= 1000) 
	             fprintf(qual9_deepsath,"%8.4lf %8.4lf\n", salt, theta);
	     }	  
	     break;
      }
      
      if (prop_avail[(int)O2]) {
          error = sscanf(&line[column[(int)O2]], "%lf", &oxygen);
          oqual = quality[qual_offset[(int)O2]];
          switch (oqual) {
             case '1':
             case '2':
	        ++oq2;
	        fprintf(qual2_o2pr,"%8.4lf %8.1lf\n", oxygen, pressure);
	        fprintf(qual2_o2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        if (pressure >= 1000) 
	           fprintf(qual2_deepo2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        break;
             case '3':
	        ++oq3;
	        fprintf(qual3_o2pr,"%8.4lf %8.1lf\n", oxygen, pressure);
	        fprintf(qual3_o2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        if (pressure >= 1000) 
	           fprintf(qual3_deepo2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        break;
             case '4':
	        ++oq4;
	        fprintf(qual4_o2pr,"%8.4lf %8.1lf\n", oxygen, pressure);
	        fprintf(qual4_o2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        if (pressure >= 1000) 
	           fprintf(qual4_deepo2th,"%8.4lf %8.4lf\n", oxygen, theta);
	        break;
	     case '6':
	         ++oq6;
	          break;
             case '9':
	        /* get CTDOXY */
	        ++oq9;
		if (prop_avail[(int)VA] ) {
                      error = sscanf(&line[column[(int)VA]], "%lf", &oxygen);
	             fprintf(qual9_o2pr,"%8.4lf %8.1lf\n", oxygen, pressure);
	             fprintf(qual9_o2th,"%8.4lf %8.4lf\n", oxygen, theta);
	             if (pressure >= 1000) 
	                fprintf(qual9_deepo2th,"%8.4lf %8.4lf\n", oxygen, theta);
		}
	        break;
          }
      }
      
     
      /* read next line into buffer and get station number */
            
      buffer = (char *) malloc((size_t)BUFSIZE * sizeof(char));
     if (buffer == NULL) {
       fprintf(stderr, "\n Unable to malloc memory for buffer.");
       exit(10);
     }
      if ( fscanf(fptr, "%[^\n]", buffer) != 1) {
         /* must be end of file */
         return(n);
      }
       error = getc(fptr) != LF;
       free(line);
       line = buffer;
       error = sscanf(&line[dcol_sta],"%d", &sta);
       error = sscanf(&line[dcol_cast],"%d", &ca);
      
   }  /* end while */
   
   if (n == 0) {
     if (flag == 0) {
       fprintf(stderr,"\nWARNING:Station #%d cast #%d was not found.\n", 
         hdr.station, cast);
     }
   }
   return(n);
   
}  /* end readdata() */
