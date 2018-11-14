/*  woce_ctd2hb.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                    author:  Ruth Curry
                             Woods Hole Oceanographic Instituls tion
			     Updated for t90 variables
................................................................................
................................................................................
.   Reads Woce ctd data files
.   extracts :    header info
.                 p,t,s,o2 
.   applies an optional filter to salinity and/or o2
.   and outputs p,d,t,s and o2 at the specified prs_int to the output
.   file.
.
.   USAGE: woce_ctd2hb infile_list -Ooutfile -Ssum_file  -P<prs_int>
          [-F<filtwidth>[/<ox_filtwidth>]] [-D<dir>] [-E<extent>]
................................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"


#define   MISSING   -9.9   /* missing value flag */
#define   FILTWIDTH  0     /* # of pts in gaussian filter */
#define   OFILTWIDTH  0     /* # of pts in gaussian filter */
#define   NPOUT     5
#define   DELTAP    2     /* pressure interval defined by WOCE format */
#define   DELTAPOUT 10    /* default pressure interval for HydroBase */
#define   BUFSIZE   10000  /* should be large enough to hold all scans */

#define   DIR    ""
#define   EXTENT ""


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
char *expocode;
int quality_byte_pos;
int p_pos, t_pos, s_pos, o_pos;
int pq_pos, tq_pos, sq_pos, oq_pos;
int is_o2;                 /* 1 if input ox units are umol/kg */
int oxflag;                /* 1 if output ox units are ml/l */
int t68flag;               /* 1 if input temperature is T68 not T90 */ 
int nrecords, nobs_read_in;
int tindex;



/*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);
int parse_header(FILE *, FILE *);
int readdata(FILE *, int *);
void gauss(double *, int, int);   


main (int argc, char **argv)
{
   int error, nobs, nprops, npout;
   int  i, j, curfile = 1, nfiles = 0;
   int filtwidth, ofiltwidth; 
   short sflag, oflag; 
   int n_ox_available; 
   char *dir, *extent, *st;
   double dlat;
   int prsint, npts, deltap;
   FILE *infile, *sumfile;
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
    sflag = oflag  = 0;
    npout = NPOUT;
    prsint = DELTAPOUT;
    filtwidth = FILTWIDTH;
    ofiltwidth = OFILTWIDTH;
    expocode = (char *) calloc((size_t)15, sizeof(char));
    oxflag = t68flag = 0;
    tindex = (int)T90;
    is_o2 = 1;
    outfile = stdout;
 
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
               case 'F':                    /* get filter width */
                        error = (sscanf(&argv[i][2],"%d", &filtwidth) != 1);
                        st = &argv[i][2];
                        while (*(st++) != '\0') {
                         if (*st == '/') {
                           ++st;
                           error = (sscanf(st,"%d", &ofiltwidth) == 1) ? 0 : 1;
                           break;
                         }
                        }
                        break;

                case 'O':                    /* get output file  */
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'P':                    /* get pressure interval */
                        error = (sscanf(&argv[i][2],"%d", &prsint) != 1);
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
                case 'X':           /* output ox in in ml/l  */
                        oxflag = 1;
                        break;
                case 'h':                    /* help  */
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

   if (! ( nfiles && sflag )) {
       fprintf(stderr,"\nYou must specify input file(s) and a sumfile.\n");
       exit(1);
   }
   
   fprintf(stderr,"\nUsing pressure interval = %d db.", prsint);
   
   if (oxflag)
      fprintf(stderr,"\nOutput oxygens units will be ml/l");
   else
      fprintf(stderr,"\nOutput oxygens units will be micromoles/kg ");
   
      
  /* set these header parameters */
  
  hdr.qual[0] = '0';
  hdr.qual[1] = '0';
  hdr.qual[2] = '0';
  hdr.qual[3] = '0';
  hdr.origin = '4';
  hdr.instrument = 'c';

/* initialize these */

  for (i= 0; i < MAXPROP; ++i)
     data.observ[i] = NULL;
     
/* loop for each input file */

   do {
   
     infile = openfile(dir, argv[curfile], extent);
     if (infile == NULL) goto NEXTFILE;
      
      /* parse the header info */
      
     nobs = parse_header(infile, sumfile);
       
     
     /*  allocate space */
     
     p = (double *) calloc((size_t)BUFSIZE, sizeof(double));
     d = (double *) calloc((size_t)BUFSIZE, sizeof(double));
     t = (double *) calloc((size_t)BUFSIZE, sizeof(double));
     s = (double *) calloc((size_t)BUFSIZE, sizeof(double));
     o = (double *) calloc((size_t)BUFSIZE, sizeof(double));
        
    /* read and filter each property */    
        
     n_ox_available = readdata(infile, &nobs);
     
     if (nobs == 0)
        goto NEXTFILE;
	
     nobs_read_in = nobs;    /* save this for later comparison */
     npout = NPOUT;
     if (n_ox_available == 0)
        npout = 4;
	
     if (filtwidth > 2) {
        gauss(s, nobs, filtwidth);
     }
     
     if (ofiltwidth > 2) {
       if (n_ox_available == nobs)
          gauss(o, nobs, ofiltwidth);
   
     }
        

     /* decimate the arrays according to specified prsint */
     
     deltap = (int) (p[2] - p[1]);
     npts = prsint / deltap;  /*  get # of pts per interval*/
     if (npts == 0) {
        fprintf(stderr,"Bad number of points per interval: %d \n", npts);
        exit(1);
     }
    
     dlat = (double) hdr.lat;
     j = 0;
     for (i = 0; i < nobs; i += npts) {
       p[j] = p[i];
       d[j] = hb_depth(p[i], dlat);
       t[j] = t[i];
       s[j] = s[i];
       o[j] = o[i];
       ++j;
     }
     
     /* add bottom observation  ... */
     
     if ((i - npts) != (--nobs)) {
        p[j] = p[nobs];
        d[j] = hb_depth(p[nobs], dlat);
        t[j] = t[nobs];
        s[j] = s[nobs];
        o[j] = o[nobs];
        ++j;
     }
     
     /* adjust some of the header info ...   */
     
     if (hdr.station > 9999)
           hdr.station = hdr.station % 10000;
           
     hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
     hdr.prop_id = (int *) calloc((size_t)npout, sizeof(int));
     hdr.prop_id[0] = (int)PR;
     hdr.prop_id[1] = (int)DE;
     hdr.prop_id[2] = tindex;
     hdr.prop_id[3] = (int)SA;
     if (npout > 4) {
        hdr.prop_id[4] = (int)O2;
        if (oxflag)
           hdr.prop_id[4] = (int)OX;
     }
     hdr.nprops = npout;
     
     hdr.nobs = data.nobs = j;
     data.nprops = hdr.nprops; 
     hdr.pdr = d[hdr.nobs-1];
     
     data.observ[(int)PR] = p;   /* set these pointers */
     data.observ[(int)DE] = d;
     data.observ[tindex] = t;
     data.observ[(int)SA] = s;
     if (npout > 4) {
      if (oxflag)
         data.observ[(int)OX] = o;
      else
         data.observ[(int)O2] = o;
    }
        
     if (hdr.nobs > 0 )
        write_hydro_station(outfile, &hdr, &data);
     
     free(hdr.prop_id);   
     free(p);
     free(d);
     free(t);
     free(s);
     free(o);
     for (i= 0; i < MAXPROP; ++i)
        data.observ[i] = NULL;
        
        
     if (nobs_read_in > nrecords)    
         fprintf(stderr,"\n%s: NRECORDS(%d) in ctd file differed from actual nobs (%d)\n", argv[curfile], nrecords, nobs_read_in); 

NEXTFILE:
         fclose(infile);

 
   } while (curfile++ < nfiles );
   

   fprintf(stderr,"\nEnd of %s\n", argv[0]);
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s converts WOCE ctd format files to HydroBase2 station format", program);
   fprintf(stderr,"\n By default, the 2db pressure series is decimated to 10 db, but the -P option permits the user to specify the output pressure interval.  A gaussian filter is by default applied to the temperature, salinity, and oxygen profiles.  To omit the filtering, specify filterwidth of 1/1 with the -F option.  You may also specify alternative filterwidths.  The necessary arguments are a list of ctd filenames, the name of the sumfile, and an output filename. \n");
   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile -Ssumfile [-Pprs_int]  [-F<filtwidth>[/<ox_filtwidth>]] [-D<dirname>] [-E<file_extent>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -O  : specifies output file ");
   fprintf(stderr,"\n    -S  : specifies summary info file ");
   fprintf(stderr,"\n   [-P] : specifies pressure interval of output series (default is 10) ");
   fprintf(stderr,"\n   [-F] : specifies filterwidth (npts) or t,s and ox (default is 5/11) ");
   fprintf(stderr,"\n   [-D] : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n   [-E] : specifies input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n   [-T] : input temperatures are T68 (default is T90) ");
   fprintf(stderr,"\n   [-X] : output oxygen units in ml/l (default is umoles/kg) ");
   fprintf(stderr,"\n   [-h] : help -- prints this message");  

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
int parse_header(FILE *fptr, FILE *sumfile)
   /* reads the header info from an already opened ctd data file and returns
      the number of observation records.  The info is used to fill in
      appropriate values in the struct HYDRO_HDR hdr, a global variable.  */
{
   int n, c, nobs, i, error, cast;
   int pos_s, pos_c, pos_t, pos_d;
   char *line, *line3, *xx, x[9], ct[4], e[15], hem;
   char csave;
   float deg, min;

   line = (char *) calloc(256, sizeof(char));
   line3 = (char *) calloc(256, sizeof(char));
   
/* Read from ctd data file... */


/* find line with EXPOCODE keyword */

   xx = NULL;
   while (xx == NULL) {
     fscanf(fptr,"%[^\n]", line);
     xx = strstr(line, "EXPOCODE");
   } 
   
   if (sscanf(&xx[9],"%s", expocode) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read expocode \n");
        exit(1);
   }
      
   error = getc(fptr);
   if (error != (int)'\n')
        fprintf(stderr,"\n expected a LF(0x0a) but got %x in parse_header()\n",
                error);
   
/* line 2 */

   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read 2nd header line. \n");
        exit(1);   
   }    
   if (sscanf(&line[7],"%d", &hdr.station) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read station# \n");
        exit(1);
   } 
   
   xx = strstr(line,"CASTNO");
   if (sscanf(&xx[6],"%d", &cast) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read cast# \n");
        exit(1);
   } 
   xx = strstr(line,"RECORDS=");
   if (xx == NULL) {
      xx = strstr(line,"Records=");
   }
   if (sscanf(&xx[8],"%d", &nobs) != 1) {
        fprintf(stderr,"\n error in sscanf attempt to read # of records.\n");
        exit(1);
   } 
   error = getc(fptr);
   nrecords = nobs;   /* save this for later comparison */
    
    
   /* line 3 */
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #3. \n");
        exit(1);   
   }         
   error = (getc(fptr) != LF);
      
   /* line 4  -- determine positions of parameters */
   
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #4. \n");
        exit(1);   
   }         
   error = (getc(fptr) != LF);
   
   xx = strstr(line,"QUALT");
   if (xx == NULL) {
     fprintf(stderr,"\nUnable to find QUALT (quality byte position). \n");
     exit(1);
   }
   quality_byte_pos = xx - line - 1;
   
   xx = strstr(line,"CTDPRS");
   if (xx == NULL) {
     fprintf(stderr,"\nUnable to find CTDPRS in header line 4. \n");
     exit(1);
   }
   p_pos = (xx -line) - 1;
   pq_pos = (xx - line) / 8;
   
   
   xx = strstr(line,"CTDTMP");
   if (xx == NULL) {
     fprintf(stderr,"\nUnable to find CTDTMP in header line 4. \n");
     exit(1);
   }
   t_pos = xx - line - 1;
   tq_pos = (xx - line) / 8;
   
   
   xx = strstr(line,"CTDSAL");
   if (xx == NULL) {
     fprintf(stderr,"\nUnable to find CTDSAL in header line 4. \n");
     exit(1);
   }
   s_pos = xx - line - 1;
   sq_pos = (xx - line) / 8;
   
   
   xx = strstr(line,"CTDOXY");
   if (xx == NULL) {
     fprintf(stderr,"\nUnable to find CTDOXY in header line 4. \n");
     o_pos = -1;
   }
   else {
      o_pos = xx - line - 1;
      oq_pos = (xx - line) / 8;
   }
     
   /* line 5 -- check units */
   
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #5. \n");
        exit(1);   
   }         
   error = (getc(fptr) != LF);

   csave = line[t_pos + 7];
   line[t_pos + 7] = '\0';           /* temporary EOString marker */  
   xx = strstr(&line[t_pos], "90");
   if (xx == NULL) {
      xx = strstr(&line[t_pos], "68");
      if (xx == NULL) {
        fprintf(stderr, "\nHeader doesn't say if temperature is T90 or T68.");
	if (t68flag) {
	  fprintf(stderr,"\n...will assume it is T68");
	  tindex = (int)TE;
	}
	else {
	  fprintf(stderr,"\n...will assume it is T90 ");
	  tindex = (int)T90;
	}
      }
   }
   else {
     tindex = (int) T90;
     if (t68flag) {
       fprintf(stderr, "\nHeader says input temperature is T90 ");
       t68flag = 0;
     }
   }  
   line[t_pos + 7] = csave;   /* restore saved character */  
      

   if (o_pos >= 0) {
     csave = line[o_pos + 7];
     line[o_pos + 7] = '\0';           /* temporary EOString marker */  
     xx = strstr(&line[o_pos], "KG");
     if (xx == NULL) {
        xx = strstr(&line[o_pos], "ML");
        if (xx == NULL) {
          fprintf(stderr, "\nCan't tell if CTDOXY units are umol/kg or ml/l.");
	  is_o2 = 1;
	  if (oxflag)
	    fprintf(stderr,"\n...will assume it is umol/kg and apply conversion factor");
	  else
	    fprintf(stderr,"\n...will assume it is umol/kg and needs no conversion factor");
        }
	else {
	  is_o2 = 0;
	}
     }
     else {
        is_o2 = 1;
        if (oxflag) {
         fprintf(stderr, "\nHeader says oxygen units are umol/kg -- will convert to ml/l for output.");
       }
     }  
     line[o_pos + 7] = csave;   /* replace EOStr */ 
   } 
      
/* line 6 -- skip */
   
   if (fscanf(fptr,"%[^\n]", line) != 1) {
        fprintf(stderr,"\n Error attempting to read header line #5. \n");
        exit(1);   
   }         
   error = (getc(fptr) != LF);
   
/* ctd data file is now positioned at start of data records. */
   
/* Now search for station in .sum file:
   first, get to and save header line #3 for its column info...       */
   
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
   pos_s = (int) (xx - line3);
   
   xx = strstr(line3, "TYPE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CAST TYPE in summary header \n");
           exit(1);
   }
   pos_t = (int) (xx - line3);  /*find offset of this string in line 3 */
   
   xx = strstr(line3, "CAST");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find CAST in summary header \n");
           exit(1);
   }
   pos_c = (int) (xx - line3);  /*find offset of this string in line 3 */   
   
   xx = strstr(line3, "DATE");
   if (xx == NULL) {
           fprintf(stderr,"\n Cannot find DATE in summary header \n");
           exit(1);
   }
   pos_d = (int) (xx - line3);  /*find offset of this string in line 3 */   



   /* Move past dashed line */
   
   error = fscanf(sumfile,"%[^\n]", line);      
   error = getc(sumfile);
  
   
/* Check each line until an appropriate match is found.  Use the
      starting lat/lon */
       
   while (fscanf(sumfile,"%[^\n]", line) == 1) {
   
      if (sscanf(line,"%s", e) != 1) {  /* expocode */
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(2);
      }
      if (sscanf(&line[pos_s],"%d", &n) != 1) { /* station # */
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(1);   
      }
      if (sscanf(&line[pos_t],"%s", ct) != 1) {
        fprintf(stderr,"\n Error attempting to parse line in summary file.\n%s\n",
         line);
        exit(1);   
      }
      if (sscanf(&line[pos_c],"%d", &c) != 1) {
           c = 0;   
      }
      
  /* if station is present, we should find one entry that meets all these criteria: */
      
      if (  (n == hdr.station) 
        && (c == cast) 
        && ((strncmp(ct,"ROS",3) == 0) || (strncmp(ct,"CTD",3) == 0)) ) {
        
        if (sscanf(&e[2],"%2s", hdr.ship) != 1) {
          fprintf(stderr,"\n error in sscanf attempt to read shipcode \n");
          exit(1);
        }
        if (sscanf(&e[0],"%2s", hdr.country) != 1) {
           fprintf(stderr,"\n error in sscanf attempt to read country # \n");
           exit(1);
        } 
        sscanf(&e[4],"%[^/]", x);
          if (sscanf(x,"%d", &hdr.cruise) != 1) {
             hdr.cruise = 9999;
            fprintf(stderr,"\n cruise # not digital: %s ", x);
            fprintf(stderr,"\n assigning it the value 9999 and continuing \n");
        } 
	if (hdr.cruise > 99999) {
	   hdr.cruise = hdr.cruise % 10000;
	   fprintf(stderr,"\n cruise number too large, using %d", hdr.cruise);
	}
   
       if (sscanf(&line[pos_d],"%2d%2d%2d", &hdr.month, &hdr.day, &hdr.year) != 3) {
           fprintf(stderr,"\n Error attempting to read date in sum file \n");
           exit(1);   
        }
        if (hdr.year < 1900)  /* adjust 2-digit year */
           hdr.year += 1900;
	   
        if (hdr.year < 1950)  /* correct for year 2000+ data */
	   hdr.year += 100;   
   
        xx = strstr(line3, "LATITUDE");
        if (xx == NULL) {
           fprintf(stderr,"\n Cannot find LATITUDE in summary header \n");
           exit(1);
        }
        i = (int) (xx - line3);  /*find offset of this string in line 3 */
        if (sscanf(&line[i],"%f %f %c", &deg, &min, &hem) != 3) {
           fprintf(stderr,"\n Error attempting to read latitude \n");
           exit(1);   
        }
        hdr.lat = deg + min / 60.;
        if (hem == 'S')
           hdr.lat = -hdr.lat;
 
        xx = strstr(line3, "LONGITUDE");
        if (xx == NULL) {
           fprintf(stderr,"\n Cannot find LONGITUDE in summary header \n");
           exit(1);
        }
        i = (int) (xx - line3);  /*find offset of this string in line 3 */
        if (sscanf(&line[i],"%f %f %c", &deg, &min, &hem) != 3) {
           fprintf(stderr,"\n Error attempting to read longitude \n");
           exit(1);   
        }
        hdr.lon = deg + min / 60.;
        if (hem == 'W')
           hdr.lon = -hdr.lon;
      
        hdr.pdr = 0;     
        xx = strstr(line3, "CDEPTH");
        if (xx != NULL) {
         i = (int) (xx - line3);  /*find offset of this string in line 3 */
         if (sscanf(&line[i],"%d", &hdr.pdr) != 1) 
            hdr.pdr = 0;
        }
       
          
        rewind(sumfile);
        free(line);
        return(nobs);
      }  /* end if */

         /* continue search */      
      error = getc(sumfile);
      
   }  /* end while */
   
   fprintf(stderr,"\n Unable to find station #%d in summary file \n",
           hdr.station);
   exit(2);
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr, int *nobs_addr)
/* Reads ctd p,t,s,ox,and quality codes for an entire station and returns the 
   number of oxygen observations. */
 
{
   int i, n, error, nobs;
   int ox_avail;
   char line[120], *xx;
   char pqual, tqual, squal, oqual;
   
   n = 0;
   i = 0;
   nobs = *nobs_addr;
   ox_avail = 0;
   
   while ( fscanf(fptr, "%[^\n]", line) != EOF) {
      ++i;
      if ((error = getc(fptr)) != LF){
        fprintf(stderr,"\nWarning:  unexpected char (%x) read instead of LF (0x10).\n", error);
      }
      
      xx = &line[quality_byte_pos];
      while (*xx == ' ') {
        ++quality_byte_pos;
	++xx;
      }
      
      pqual = *(xx+pq_pos);
      tqual = *(xx+tq_pos);
      squal = *(xx+sq_pos);
      oqual = '9';
      if (o_pos >= 0) 
         oqual = *(xx+oq_pos);
      
              
      switch (pqual) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '5': /* fall through */
         case '9':
            continue;
	 case '1':
	 case '2':
	 case '6':    /* data okay */
	 case '7':    /* fall through */
	 case '8':    
	    break;
         default:
	    fprintf(stderr,"\nUnknown pr quality code [%c]", pqual);
            exit(1);
      } /* end switch */
     
      switch (squal) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '5': /* fall through */
         case '9':
            continue;
	 case '1':
	 case '2':    /* data okay */
	 case '6':    /* fall through */
	 case '7':    
	 case '8':    
	    break;
         default:
	    fprintf(stderr,"\nUnknown sa quality code [%c]", squal);
            exit(1);
      } /* end switch */
      
      switch (tqual) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '5': /* fall through */
         case '9':
            continue;
	 case '1':
	 case '2':
	 case '6':    /* data okay */
	 case '7':    /* fall through */
	 case '8':    
	    break;
         default:
	    fprintf(stderr,"\nUnknown te quality code [%c]", tqual);
            exit(1);
      } /* end switch */
      
      error = sscanf(&line[p_pos], "%lf", &p[n]) != 1;
      if (error){
         fprintf(stderr,"\nError parsing pressure at data line #%d.\n", i);
         exit(1);
      }
      error = sscanf(&line[t_pos], "%lf", &t[n]) != 1;
      if (error){
         fprintf(stderr,"\nError parsing temperature at data line #%d.\n", i);
         exit(1);
      }
     
      error = sscanf(&line[s_pos], "%lf", &s[n]) != 1;
      if (error){
         fprintf(stderr,"\nError parsing salinity at data line #%d.\n", i);
         exit(1);
      }
      
      switch (oqual) {
         case '3': 
         case '4': /* these are codes for unacceptable data */
         case '5': /* fall through */
         case '9':
            o[n] = MISSING;
            break;
	    
         default:     /* 1, 2, or 6 */
           error = sscanf(&line[o_pos], "%lf", &o[n]) != 1;
           if (error){
              fprintf(stderr,"\nError parsing oxygen at data line #%d.\n", i);
              exit(1);
           }
            ++ox_avail;
            if (oxflag && is_o2)
               o[n] = ox_kg2l(o[n], p[n],t[n],s[n]);   /* convert from micromoles/kg to ml/l */
	    if (!oxflag && !is_o2)
               o[n] = ox_l2kg(o[n], p[n],t[n],s[n]);   /* convert from micromoles/kg to ml/l */
	    
            
      } /* end switch */
      
      ++n;
   } /* end while */
   
   *nobs_addr = n;
   return(ox_avail);
   
}  /* end readdata() */

/*******************************************/
void gauss(double *x, int nx, int width)

   /*  Uses a weighted arithmetic average to smooth an array of values.
       The filter width is the total number of pts which contribute to
       any one smoothed value and should be an odd number to incorporate an
       equal number of points above and below each point */
{
   double *tmp, *wght;
   double pi = 3.1415926;
   double alpha, sum;
   int i, j, k, halfwidth, window, end;
      
   tmp = (double *) calloc((size_t)nx, sizeof(double));
   
   halfwidth = width >> 1;      /* divide by 2 */
   width = (halfwidth << 1) + 1;   /* ensure that width is an odd number */
   wght = (double *) calloc((size_t) width, sizeof(double));
   
   alpha = (pi / (double) halfwidth);
   alpha = alpha * alpha / 4.5;
   sum = 0.0;
   
   /* first , calculate weights for bottom half of filter window... */
   
   i = 0;
   while (i < halfwidth) {
       wght[i] = halfwidth - i;
       wght[i] = exp(-alpha * wght[i] * wght[i]);
       sum += wght[i];
       ++i;
   }
   wght[halfwidth] = 1.0;    /* weight at central element */
   sum = sum + sum + 1.0;    /* sum of weights for entire filter */
   
   /* normalize the weights */
   
   for (i = 0; i <= halfwidth; ++i) 
      wght[i] /= sum;
   
   /* now set the weights in the top half of the filter to mirror the bottom */
   
   j = halfwidth - 1;
   for (i = halfwidth+1; i < width; ++i ) {
      wght[i] = wght[j--];
   }
   
   /* copy the data array into temporary space */
   
   for (i = 0; i < nx; ++i) {
      tmp[i] = x[i];
   }
   
   /* apply the filter -- this leaves the values at the top and bottom
      of the array (distance of halfwidth) unchanged */
   
   end = nx - halfwidth;
   for (i = halfwidth; i < end; ++i) {
      x[i] = 0.0;
      window = i + halfwidth;
      k = 0;
      for (j = i-halfwidth; j <= window; ++j) 
         x[i] += tmp[j] * wght[k++];
   }
   free((void *)tmp);
   free((void *)wght);
   return;
   
}  /* end gauss() */
/******************************************************************/
