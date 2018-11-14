/*  woce_btl2hb.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Institution
			     Updated for T90 variables Nov 2011
................................................................................
.   Reads Woce hydrographic data files
.   extracts :    header info
.                 p,d,t,s,ox,n2,n3,si,p4,f1,f2,he,tu
.   
.   USAGE: woce_btl2hb infile_list -Ooutfile -Ssum_file  [-D<dir>] [-E<extent>]
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
struct HYDRO_DATA data;
char *buffer;
char *expocode;
short *prop_avail; 
int *column;
int *qual_offset;
int cast;
int qual_word_col;
int replace_ox, replace_sa;
int col_sta, col_cast, col_type, col_date;
int col_code, col_lat, col_lon, col_pdr, col_nobs;
int dcol_sta, dcol_cast;
int is_o2;                /* 1 if input Oxygen is in ml/l */
int oxsub;                /* 1 to substitute ctd oxygen for bad bottle ox */
int oxunits_kg;         /* 1 to output oxygen in micromoles/kg */
int t68flag;               /* 1 if input temperature is T68 not T90 */ 
int tindex;

/* prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);
void parse_headers(FILE *, FILE *);
int find_next_station(FILE *) ;
int readdata(FILE *);

int main ( int argc, char **argv)
{
   int error, nobs, nprops, npout;
   int  i, j, curfile = 1, nfiles = 0;
   short sflag;
   int n_ox_available; 
   char *dir, *extent, *st;
   double dlat;
   int  npts;
   FILE *infile, *sumfile, *openfile();
   FILE * outfile;
   
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
    is_o2 = 1;  /* assume input ox units are micromole/kg */
    oxsub = 1;  /* substitute ctd ox for bad bottle values */
    oxunits_kg = 1;  /* output units in micromoles/kg */
    replace_sa = replace_ox = 0;
    expocode = (char *) calloc((size_t)15, sizeof(char));
    outfile = stdout;
    prop_avail = (short *) calloc((size_t)MAXPROP, sizeof(short));
    column = (int *) calloc((size_t)MAXPROP, sizeof(int));
    qual_offset = (int *) calloc((size_t)MAXPROP, sizeof(int));
    tindex = (int)T90;
    
    for (i = 0; i < MAXPROP; ++i) {  /* initialize */
	data.observ[i] = NULL;
    }


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
               case 'O':                    /* get output file  */
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;
               case 'R':
	       
	             st = &argv[i][2];
		     switch (*st) {
		        case '\0':
			   replace_ox = 1;
			   replace_sa = 1;
			   break;
			case 'O':
			case 'o':
			   replace_ox = 1;
			   break;
			case 'S':
			case 's':
			   replace_sa = 1;
			   break;
			default:
			   ++error;   
		     }
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
               case 'T':           /* input temperature is T68 not T90 */
                        t68flag = 1;
                        break;
               case 'X':           /* output oxygen units in ml/l  */
                        oxunits_kg = 0;                       
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
   if (replace_ox) {
      fprintf(stderr,"\nWill substitute CTD oxygen values for all bottle oxygen values ");
   
   }
   if (replace_sa) {
      fprintf(stderr,"\nWill substitute CTD salinity values for all bottle  values ");
   
   }
   if (!replace_ox) {
      if (oxsub)
         fprintf(stderr,"\nWill substitute bad bottle oxygen values with CTD O2  ");
      else
         fprintf(stderr,"\nThis will only extract oxygen bottle values with a quality code of 2");
   }
   if (oxunits_kg)
      fprintf(stderr,"\nO2 units will be output in micromoles/kg  ");
   else
      fprintf(stderr,"\nO2 units will be output in ml/l");
   
   if (t68flag)
      fprintf(stderr,"\nInput temperatures are T68 ");

  infile = openfile(dir, argv[curfile], extent);
  if (infile == NULL) 
        exit(1);

  /* set these header parameters */
  
  hdr.qual[0] = '1';
  hdr.qual[1] = '1';
  hdr.qual[2] = '1';
  hdr.qual[3] = '1';
  hdr.origin = '4';
  hdr.instrument = 'b';

      
  /* parse the header info from datafile and sumfile*/
      
  parse_headers(infile, sumfile);
  npout = 0;
  for (i = 0; i < MAXPROP; ++i) {
      
       if ( i != (int)S2 && (i != (int)VA) && prop_avail[i] )  /* S2 stores CTDsalt not sigma-2, VA stores CTDoxygen */
          ++npout;
  }
  ++npout;  /* add depth to output list */
  
  /* initialize these... */
  
  hdr.station = cast = 0; 
  hdr.nprops = npout;
       
  /* loop for each station */
  
  fprintf(stderr,"\nSearching for station ");
    
  while ((nobs = find_next_station(sumfile)) >= 0) {
  
    /* allocate space and read in each property */
    
    hdr.prop_id = (int *) calloc((size_t)npout, sizeof(int));
    j = 0; 
    prop_avail[(int)DE] = 1;       /* set this so space for depth is allocated */
    for (i = 0; i < MAXPROP; ++i) {
        if (prop_avail[i]) {
           data.observ[i] = (double *) calloc((size_t)nobs, sizeof(double));
	   if ( i != (int) S2 && i != (int)VA) {  /* skip CTDsalt and CTDoxy */
              hdr.prop_id[j] = i;
	      if (hdr.prop_id[j] == (int)OX && oxunits_kg) {
	         hdr.prop_id[j] = (int)O2;
                 data.observ[(int)O2] = (double *) calloc((size_t)nobs, sizeof(double));
	      }
	      
	      if (hdr.prop_id[j] == (int)O2 && !oxunits_kg) {
	         hdr.prop_id[j] = (int)OX;
                 data.observ[(int)OX] = (double *) calloc((size_t)nobs, sizeof(double));
	      }
	      
              ++j;
	   }
        }
    }
    prop_avail[(int)DE] = 0;      /* restore value */
    
    nobs = readdata(infile);
    if (nobs == EOF) {
       fprintf(stderr, "\n Station %d not found in data file.\n ", hdr.station);
    }
    
    if (nobs > 0) {
        dlat = (double) hdr.lat;
        j = 0;
        for (i = 0; i < nobs; ++i) {
          data.observ[(int)DE][i] = hb_depth(data.observ[(int)PR][i], dlat);
        }
     
     /* adjust some of the header info ...   */
     
        hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
        if (hdr.station > 9999)
           hdr.station = hdr.station % 10000;
     
        hdr.nobs = data.nobs = nobs;
        data.nprops = hdr.nprops; 
	hdr.pdr = data.observ[(int)DE][hdr.nobs-1] + 10.0;
        
        write_hydro_station(outfile, &hdr, &data);
        
    }  /* end if */
     
     for (i = 0; i < MAXPROP; ++i) {
        if ( data.observ[i] != NULL) 
           free(data.observ[i]);
     }
     free(hdr.prop_id);   
     
   }  /* end while */

   fprintf(stderr,"\nEnd of conversion.\n");
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s infile -Ssumfile [-D<dirname>] [-E<file_extent>] [-Ooutfile][-Rs] [-Ro] [-T] [-X] [-h]", program);
   fprintf(stderr,"\n -S  : summary info file. ex:  316N151.sum ");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n[-D] : dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n[-E] : input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n[-O] : name of output file [default is output to screen] ");
   fprintf(stderr,"\n[-Ro] : Replace bottle oxygens with CTD oxygens ");
   fprintf(stderr,"\n[-Rs] : Replace bottle salts with CTD salts ");
   fprintf(stderr,"\n[-T] : input temperatures are T68 (default is T90) ");
   fprintf(stderr,"\n[-X] : output oxygn units in ml/l [default is micromoles/kg] ");
   fprintf(stderr,"\n-h : help...... prints this message.");

   fprintf(stderr,"\n\n");  
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
   int n, i, error, col;
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
   
   if (sscanf(&expocode[0],"%2s", hdr.country) != 1) {
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
   
   col = (int) (xx - line - 1);
   
   /* check for T90 in units */
   line2[col +7] = '\0';   /* temporary EOString marker */
   xx = strstr(&line2[col], "90");
   if (xx == NULL) {
      xx = strstr(&line2[col], "68");
      if (xx == NULL) {
        fprintf(stderr, "\nSEA file doesn't say if temperature is T90 or T68.");
	if (t68flag) {
	  fprintf(stderr,"\n...will assume it is T68 ");
	  tindex = (int)TE;
	}
	else {
	  fprintf(stderr,"\n...will assume it is T90 ");
	  tindex = (int) T90;
	}
      }
      else {
        t68flag = 1;
	tindex = (int) TE;
      }
   }
   else {
     if (t68flag) {
       fprintf(stderr, "\nSEA file says input temperature is T90 ");
       t68flag = 0;
       tindex = (int) T90;
     }
   }
   prop_avail[tindex] = 1;
   column[tindex] = col;
   qual_offset[tindex] = 0;
  
   
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
        if (xx == NULL) {
	  is_o2 = 1;
          fprintf(stderr, "\nCan't tell if input OXYGEN units are umol/kg or ml/l.");
          fprintf(stderr, "\nWill assume umole/kg.  Check to make certain and update headers records if necessary.");
        }
        else {
                 fprintf(stderr, "\nInput oxygen units are ml/l.");
		 is_o2 = 0;
	}
     }
     else {
         is_o2 = 1;
         fprintf(stderr, "\nInput oxygen units are umol/kg.");
      } 
   }

   xx = strstr(line, "SILCAT");
   if (xx != NULL) {
      prop_avail[(int)SI] = 1;
      column[(int)SI] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)SI]],"%d", &qual_offset[(int)SI]);
   }
   
   xx = strstr(line, "NITRIT");
   if (xx != NULL) {
      prop_avail[(int)N2] = 1;
      column[(int)N2] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)N2]],"%d", &qual_offset[(int)N2]);
   }
   
   xx = strstr(line, "NITRAT");
   if (xx != NULL) {
      prop_avail[(int)N3] = 1;
      column[(int)N3] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)N3]],"%d", &qual_offset[(int)N3]);
   }
   
   xx = strstr(line, "PHSPHT");
   if (xx != NULL) {
      prop_avail[(int)P4] = 1;
      column[(int)P4] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)P4]],"%d", &qual_offset[(int)P4]);
   }
   
   xx = strstr(line, "CFC-11");
   if (xx != NULL) {
      prop_avail[(int)F1] = 1;
      column[(int)F1] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)F1]],"%d", &qual_offset[(int)F1]);
   }
   
   xx = strstr(line, "CFC-12");
   if (xx != NULL) {
      prop_avail[(int)F2] = 1;
      column[(int)F2] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)F2]],"%d", &qual_offset[(int)F2]);
   }
   xx = strstr(line, "CFC113");
   if (xx != NULL) {
      prop_avail[(int)F3] = 1;
      column[(int)F3] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)F3]],"%d", &qual_offset[(int)F3]);
   }
   
   
   xx = strstr(line, "HELIUM");
   if (xx != NULL) {
      prop_avail[(int)HE] = 1;
      column[(int)HE] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)HE]],"%d", &qual_offset[(int)HE]);
   }
   
   xx = strstr(line, "TRITUM");
   if (xx != NULL) {
      prop_avail[(int)TU] = 1;
      column[(int)TU] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)TU]],"%d", &qual_offset[(int)TU]);
   }

   xx = strstr(line, "  CCL4");
   if (xx != NULL) {
      prop_avail[(int)CL4] = 1;
      column[(int)CL4] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)CL4]],"%d", &qual_offset[(int)CL4]);
   }
   
   xx = strstr(line, "   SF6");
   if (xx != NULL) {
      prop_avail[(int)SF6] = 1;
      column[(int)SF6] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)SF6]],"%d", &qual_offset[(int)SF6]);
      
   }
   else {
      xx = strstr(line, "  SF-6");
      if (xx != NULL) {
         prop_avail[(int)SF6] = 1;
         column[(int)SF6] = (int) (xx - line - 1);
         sscanf(&line3[column[(int)SF6]],"%d", &qual_offset[(int)SF6]);
      }
   }
   
   xx = strstr(line, "TCARBN");
   if (xx != NULL) {
      prop_avail[(int)CO2] = 1;
      column[(int)CO2] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)CO2]],"%d", &qual_offset[(int)CO2]);
   }
   
   xx = strstr(line, "ALKALI");
   if (xx != NULL) {
      prop_avail[(int)ALK] = 1;
      column[(int)ALK] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)ALK]],"%d", &qual_offset[(int)ALK]);
   }
   
   xx = strstr(line, "   DOC");
   if (xx != NULL) {
      prop_avail[(int)DOC] = 1;
      column[(int)DOC] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)DOC]],"%d", &qual_offset[(int)DOC]);
   }
   
   xx = strstr(line, "   TDN");
   if (xx != NULL) {
      prop_avail[(int)TDN] = 1;
      column[(int)TDN] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)TDN]],"%d", &qual_offset[(int)TDN]);
   }
   
   xx = strstr(line, "DELC13");
   if (xx != NULL) {
      prop_avail[(int)C13] = 1;
      column[(int)C13] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)C13]],"%d", &qual_offset[(int)C13]);
   }
   
    xx = strstr(line, "DELC14");
   if (xx != NULL) {
      prop_avail[(int)C14] = 1;
      column[(int)C14] = (int) (xx - line - 1);
      sscanf(&line3[column[(int)C14]],"%d", &qual_offset[(int)C14]);
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
   char squal, tqual;
   
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
      tqual = quality[qual_offset[tindex]];

      scanOK = 1;
      if (replace_sa) {
      
      }
      switch (squal) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '9':
            scanOK = 0;  /* check CTDSAL quality */
	    if ( prop_avail[(int)S2] &&  quality[qual_offset[(int)S2]] == '2' ) {
		  scanOK = 1;
	    }
            break;
         default:
            break;
      } /* end switch */
              
      switch (tqual) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '9':
            scanOK = 0;
            break;
         default:
            break;
      } /* end switch */
      
      if (scanOK) {     
         for (indx = 0; indx < MAXPROP; ++indx) {
	   if (indx == (int)S2) {
	       ;   /* skip CTDSAL */
           }
	   else if (prop_avail[indx]) {
              error = sscanf(&line[column[indx]], "%lf", &data.observ[indx][n]);
              if (indx != (int)PR) {
	      
                switch (quality[qual_offset[indx]]) {
                  case '3': 
                  case '4': /* these are codes for unacceptable data */
                  case '9':
                    data.observ[indx][n] = MISSING;
		    
		    if (indx == (int)SA)  {   /* Substitute CTD salt */
                         error = sscanf(&line[column[(int)S2]], "%lf", &data.observ[indx][n]);
		    }
		    if (indx == (int)O2  && oxsub)  {   /* Substitute CTD oxygen */
		        if (prop_avail[(int)VA] && quality[qual_offset[(int)VA]] == '2' )
                            error = sscanf(&line[column[(int)VA]], "%lf", &data.observ[indx][n]);
		    }
                    break;
                  default:
                    break;
                } /* end switch */
		
		if ((indx == (int)SA) && replace_sa)  {   /* Substitute CTD salt */
                         error = sscanf(&line[column[(int)S2]], "%lf", &data.observ[indx][n]);
		 }
		if  ((indx == (int)O2) || (indx == (int)OX ))  {   
		   if (replace_ox) { /* Substitute CTD oxygen */
		        if (prop_avail[(int)VA] && quality[qual_offset[(int)VA]] == '2' )
                            error = sscanf(&line[column[(int)VA]], "%lf", &data.observ[indx][n]);
	           }
		 }  
		
              }
           }
         }  /* end for */ 
	 
	 /* convert T90 to T68 if appropriate */
	 
         
         /* adjust oxygen units, if necessary */
         
         if (prop_avail[(int)O2]) {
            if ( ! oxunits_kg)  {   /* output ox in ml/l */
	       if ( is_o2)  /* input ox is in micromoles/kg */
                  data.observ[(int)OX][n] = ox_kg2l(data.observ[(int)O2][n], data.observ[(int)PR][n],data.observ[tindex][n], data.observ[(int)SA][n]);
		  
	       if (!is_o2) {
	          data.observ[(int)OX][n] = data.observ[(int)O2][n];
	       }
	    }
	    else {    /* output ox in umol/kg */
	        if (!is_o2) {  /* input ox is in ml/l */
                  data.observ[(int)O2][n] = ox_l2kg(data.observ[(int)O2][n], data.observ[(int)PR][n],data.observ[tindex][n], data.observ[(int)SA][n]);
		}
		
	    }
         }
         
         if ((data.observ[tindex][n] > -8.) && (data.observ[(int)SA][n] > -8.)&&(data.observ[(int)PR][n] > -8.))     
             ++n;
      } /* end if scanOK */

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
     else {
       fprintf(stderr,"\nWARNING:Station #%d cast #%d had no valid observations.\n",
         hdr.station, cast);
     }
   }
   return(n);
   
}  /* end readdata() */
