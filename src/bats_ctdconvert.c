/*  bbsr2_ctdconvert.c
................................................................................
                         *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated 2009 to handle new ctd format from BBSR
................................................................................
................................................................................
.   Reads Bermuda Bio Station ctd format file(s)
.   extracts :    header info
.                 p,t,s,o2,  
.   
.   and outputs p,d,t,s and o at the specified prs_int to the output
.   file.
.
.   USAGE: ctdconvert infile_list -Ooutfile -Eerrorfile -P<prs_int> [-D<dir>]
................................................................................
*/



#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"


#define   LF        0x0a   /* ascii code for linefeed char */
#define   MISSING   -9.0   /* missing value flag */

#define   DIR    ""


  /* globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;
double *p, *d, *t, *s, *o;
char pr_avail, te_avail, sa_avail, o2_avail;
char line_buffer[200];
int deltap;          /* pressure interval of ctd cast */

   /* prototypes for locally define functions */

void print_usage(char *);
int readdata(FILE *, int);
int parse_header(FILE *, FILE *);
FILE *openfile(char *, char *, char *);

int main (int argc, char **argv)
{
  int error, nobs, nprops, cast_id;
  FILE *outfile, *errorfile;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, pflag, eflag; 
   char *dir, *extent;
   double dlat;
   int prsint, npts;
   FILE *infofile, *datafile;
   char *infofile_extension = "_info.txt";
   char *datafile_extension = "_ctd.txt";
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    error = 0;
    pflag = oflag = eflag = 0;
    deltap = 2.0;
 
 
/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        dir = &argv[i][2];
                        break;

               case 'O':                    /* get output file  */
                        oflag = 1;
                        outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'P':                    /* get pressure interval */
                        pflag = 1;
                        error = (sscanf(&argv[i][2],"%d", &prsint) != 1);
                        break;

	       case 'E':                    /* get error output file */
	                eflag = 1;
		        errorfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (errorfile == NULL) {
                           fprintf(stderr,"\nError opening error file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
			}
			break;

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

   if (! (oflag && nfiles && pflag && eflag)) {
       fprintf(stderr,"\nYou must specify input file(s), output file, pressure interval and error file.\n");
       exit(1);
   }

/* loop for each input file */

   do {
     infofile = openfile(dir, argv[curfile], infofile_extension);
     if (infofile == NULL) goto NEXTFILE;
     datafile = openfile(dir, argv[curfile], datafile_extension);
     if (datafile == NULL) {
       fprintf(stderr,"\nData file %s not found.\n", argv[curfile]);
       goto NEXTFILE;
     }
     error = fscanf(datafile, "%[^\n]", line_buffer);

     while (parse_header(infofile, datafile) == 1) {
      /* parse the header info */
         
       nprops = 3;   /* assume pr, de, te are available. */
       if (sa_avail)   
         ++nprops;
       if (o2_avail)   
         ++nprops;
       hdr.nprops = nprops;
            
       /* adjust some of the header info ...   */
     
       dlat = (double) hdr.lat;
       strncpy(hdr.country,"32",3);
       strncpy(hdr.ship,"0G",3);
       for (i = 0; i < NQUAL; ++i)
	 hdr.qual[i] = '0';
       hdr.origin = '5';
       hdr.instrument = 'c';
       hdr.pdr = 0;
       hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
       hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
       i = 0;
       hdr.prop_id[i] = (int)PR;
       hdr.prop_id[++i] = (int)DE;
       hdr.prop_id[++i] = (int)T90;
       if (sa_avail)
	 hdr.prop_id[++i] = (int)SA;
       if (o2_avail) 
         hdr.prop_id[++i] = (int)O2;
     
     /*  allocate space */
     
       nobs = 6000;
     
       p = (double *) calloc(nobs, sizeof(double));
       d = (double *) calloc(nobs, sizeof(double));
       t = (double *) calloc(nobs, sizeof(double));
       if (sa_avail)
	 s = (double *) calloc(nobs, sizeof(double));
       if (o2_avail) 
         o = (double *) calloc(nobs, sizeof(double));
        
    /* read and filter each property */    
        
       nobs = readdata(datafile, hdr.station);

     /* decimate the arrays according to specified prsint */
     
       npts = prsint / deltap;  /*  get # of pts per interval*/
       if (prsint % deltap)
        fprintf(stderr,"WARNING: prsint (%d) requested is not an even multiple of the pressure sorted ctd file: output will be a %d db series\n", prsint, (npts*deltap));
     
       if (npts == 0) {
         fprintf(stderr,"Bad number of points per interval: %d \n", npts);
         exit(1);
       }
    
       j = 0;
       for (i = 0; i < nobs; i += npts) {
         p[j] = p[i];
         d[j] = d[i];
	 t[j] = t[i];
	 if (sa_avail)
	   s[j] = s[i];
         if (o2_avail) 
	   o[j] = o[i];
         ++j;
       }
     
     /* add bottom observation  ... */
     
       if ((i - npts) != (--nobs)) {
         p[j] = p[nobs];
         d[j] = d[nobs];
	 t[j] = t[nobs];
	 if (sa_avail)
	   s[j] = s[nobs];
         if (o2_avail) 
	   o[j] = o[nobs];
         ++j;
       }
     
       hdr.nobs = data.nobs = j;
       data.nprops = hdr.nprops; 
     
       data.observ[(int)PR] = p;   /* set these pointers */
       data.observ[(int)DE] = d;
       data.observ[(int)T90] = t;
       if(sa_avail)
	 data.observ[(int)SA] = s;
       if(o2_avail)
	 data.observ[(int)O2] = o;
        
       if (pr_avail && te_avail && sa_avail)
         write_hydro_station(outfile, &hdr, &data);
       else
	 write_hydro_station(errorfile, &hdr, &data);
     
       free(hdr.prop_id);   
       free(p);
       free(d);
       free(t);
       if (sa_avail)
	 free(s);
       if (o2_avail)
         free(o);
       }

     error = fclose(infofile);
     error = fclose(datafile);


NEXTFILE: 
     ;

   } while (curfile++ < nfiles );

   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s station_list  -E<error_file> -Ooutfile -Pprs_int [-D<dirname>] ", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file for incomplete profiles");  
   fprintf(stderr,"\n    -O  : specifies output filename");
   fprintf(stderr,"\n    -P  : specifies pressure interval for output file");  
   fprintf(stderr,"\n            ex: -P10 ");

   fprintf(stderr,"\n\n");  
   return;
}
   
/*****************************************************************************/
FILE *openfile(char *dir, char *root, char *extent)
{
   char st[500];
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
int parse_header(FILE *fptr, FILE *fptr2)

   /* reads the header info from an already opened ctd data file 
     The info is used to fill in appropriate values in
      the struct HYDRO_HDR hdr, a global variable.  */
{
  int error, i;
  char line[1000], st1[20], st2[20], st3[20], st4[20], st5[20], st6[20], st7[20], st8[20], st9[20], st10[20], st11[20], st12[20];

  double data_pr, data_te, data_sa, data_o2;
  pr_avail = 1;
  te_avail = 1;
  sa_avail = 1;
  pr_avail = 1;
  o2_avail = 1;
  error = fscanf(fptr, "%[^\n]", line);
  error = getc(fptr);
   
   error = sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s%s", st1, st2, st3, st4, st5, st6, st7, st8, st9, st10, st11, st12);
   if (error ==12) {

   /*parsing strings*/

     error = sscanf(st1, "%5d%2d", &hdr.cruise, &hdr.station);
     error = sscanf(st2, "%4d%2d%2d", &hdr.year, &hdr.month, &hdr.day);
     error = sscanf(st9, "%f", &hdr.lat);
     error = sscanf(st11, "%f", &hdr.lon);
     if (hdr.lon >0 )
       hdr.lon *= -1; /* longitude shoud be negative */

     error = sscanf(line_buffer, "%s%s%s%s%s%s%s%s%s%s%s", st1, st2, st3, st4, st5, st6, st7, st8, st9, st10, st11);
     error = sscanf(st5, "%lf", &data_pr);
     error = sscanf(st7, "%lf", &data_te); 
     error = sscanf(st8, "%lf", &data_sa);
     error = sscanf(st9, "%lf", &data_o2);

     if(data_pr < -8) {
       pr_avail = 0;
       fprintf(stderr, "\nWARNING: no PR in cast %d.\n", hdr.station);
     }

     if(data_te < -8) {
       te_avail = 0;
       fprintf(stderr, "\nWARNING: no T90 in cast %d.\n", hdr.station);
     }

     if(data_sa < -8) {
       sa_avail = 0;
       fprintf(stderr, "\nWARNING: no SA in cast %d.\n", hdr.station);
     }

     if(data_o2 < -8) {
       o2_avail = 0;
       fprintf(stderr, "\nWARNING: no O2 in cast %d.\n", hdr.station);
     }

     return 1;
   }
   else
     return -1;
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr, int cast_id)
{
  int nobs, i, error, propcount, nread, cruise, scan_id, new_cast = 0;
   double buffer[5], unimportant[3];

   nobs = 0;  
   propcount = 1;
   if (pr_avail)
      ++propcount;
   if (te_avail) 
      ++propcount;
   if (sa_avail) 
      ++propcount;
   if (o2_avail) 
      ++propcount;

   /* read in from line buffer first time */

   nread = sscanf(line_buffer,"%d%lf%lf%lf%lf%lf%lf%lf%lf",&scan_id, &unimportant[0], &unimportant[1], &unimportant[2], &buffer[0],&buffer[1],&buffer[2],&buffer[3], &buffer[4]);
      


   /*loop through the rest of the station*/

   do{
      error = getc(fptr);
      
      nread = sscanf(line_buffer,"%d%lf%lf%lf%lf%lf%lf%lf%lf",&scan_id, &unimportant[0], &unimportant[1], &unimportant[2], &buffer[0],&buffer[1],&buffer[2],&buffer[3], &buffer[4]);
      
      if(scan_id % 100 == cast_id) {

	i = 0;
        p[nobs] = buffer[0];

        d[nobs] = buffer[1];

        if (te_avail) {
          t[nobs] = buffer[2];
        } 

        if (sa_avail){ 
          s[nobs] = buffer[3];
        }  
        if (o2_avail) {
          o[nobs] = buffer[4];
        }
      
        if (nread < propcount + 4){
          fprintf(stderr,"\nOnly able to parse %d properties at data line: \n%s\n", nread, line_buffer);
          exit(1);
        }
	++nobs;
      }

      else 
	new_cast = 1;

   }   while ((fscanf(fptr, "%[^\n]", line_buffer) != EOF) && !new_cast);
   return (nobs);
}  /* end readdata() */
