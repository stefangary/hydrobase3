/*  bats_btlconvert.c
................................................................................
                         *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Matthew McGowen
                             Woods Hole Oceanographic Inst
                             2009
................................................................................
................................................................................
.   Reads Bermuda Bio Station btl format file(s)
.   extracts :    header info
.                 d,t,s,ox,n2,p4,si 
.   determines whether nutrient information is available in a btl file
.   and outputs p,d,t,s,o and available nutrient data to the output
.   file.
.
.   USAGE: bbsr_btl infile -Ooutfile -Eerrorfile [-D<dir>]
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
char line_buffer[300];
double *p, *d, *t, *s, *o, *n, *phos, *sil;

   /* prototypes for locally define functions */

void print_usage(char *);
int readdata(FILE *);
int parse_header(FILE *);
FILE *openfile(char *, char *);

int main (int argc, char **argv)
{
   int error, nobs, nprops, cast_id, bad_d, bad_t, bad_s, bad_o, bad_n, bad_phos, bad_sil, is_error;
   FILE *outfile, *errorfile;
   int  i, j, curfile = 1, nfiles = 0; 
   short oflag, eflag; 
   char *dir, st[20];
   double dlat;
   int prsint, npts;
   FILE *infile;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    error = 0;
    oflag = eflag = 0;
 
 
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

               case 'E':                    /* get error file  */
                        eflag = 1;
                        errorfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
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

   if (! (eflag && oflag && nfiles)) {
       fprintf(stderr,"\nYou must specify input file(s), output file, and error file\n");
       exit(1);
   }

   /* get the file pointer to the first line of data */
   do {
     infile = openfile(dir, argv[curfile]);
     if (infile == NULL) goto NEXTFILE;
     do {
       error = fscanf(infile, "%[^\n]", line_buffer);
       error = getc(infile);
       error = sscanf(line_buffer, "%s", st);
     } while (strcmp(st, "/Variables") != 0);

     for(i = 0; i < 2; ++i) {
       error = fscanf(infile, "%[^\n]", line_buffer);
       error = getc(infile);
     }

     while (parse_header(infile) == 1) {
       is_error = 0;
      /* parse the header info */
         
       nprops = 8;   /* assume all props are available. */
            
       /* adjust some of the header info ...   */
     
       dlat = (double) hdr.lat;
       strncpy(hdr.country,"32",3);
       strncpy(hdr.ship,"0G",3);
       for (i = 0; i < NQUAL; ++i)
	 hdr.qual[i] = '0';
       hdr.origin = '5';
       hdr.instrument = 'c';
       hdr.pdr = 3015;
       hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
       hdr.prop_id = (int *) malloc(nprops * sizeof(int));
       i = 0;

     
     /*  allocate space */
     
       nobs = 600;
     
       p = (double *) calloc(nobs, sizeof(double));
       d = (double *) calloc(nobs, sizeof(double));
       t = (double *) calloc(nobs, sizeof(double));
       s = (double *) calloc(nobs, sizeof(double));
       o = (double *) calloc(nobs, sizeof(double));
       n = (double *) calloc(nobs, sizeof(double));
       phos = (double *) calloc(nobs, sizeof(double));
       sil = (double *) calloc(nobs, sizeof(double));

        
    /* read and filter each property */    
        
       nobs = readdata(infile);
     
       if (npts == 0) {
         fprintf(stderr,"Bad number of points per interval: %d \n", npts);
         exit(1);
       }
       bad_d = bad_t = bad_s = bad_o = bad_n = bad_phos = bad_sil = 0;
       for (i = 0; i < nobs; i++) {
	 if(d[i]< -8)
	   ++bad_d;
	 if(t[i]< -8)
	   ++bad_t;
	 if(s[i]< -8)
	   ++bad_s;
	 if(o[i]< -8)
	   ++bad_o;
	 if(n[i]< -8)
	   ++bad_n;
	 if(phos[i]< -8)
	   ++bad_phos;
	 if(sil[i]< -8)
	   ++bad_sil;
       }
     
       hdr.nobs = data.nobs = nobs;
       hdr.nprops = 0;
       i = 0;
       if(bad_d < nobs) {
         hdr.prop_id[i] = (int)PR;
         hdr.prop_id[++i] = (int)DE;
         data.observ[(int)PR] = p;   /* set these pointers */
         data.observ[(int)DE] = d;
	 hdr.nprops += 2;
       }
       else 
	 is_error = 1;
       if(bad_t < nobs) {
         hdr.prop_id[++i] = (int)TE;
         data.observ[(int)TE] = t;
	 hdr.nprops++;
       }
       else
	 is_error = 1;
       if(bad_s < nobs) {
	 hdr.prop_id[++i] = (int)SA;
	 data.observ[(int)SA] = s;
	 hdr.nprops++;
       }
       else
	 is_error = 1;
       if(bad_o < nobs) {
         hdr.prop_id[++i] = (int)O2;
	 data.observ[(int)O2] = o;
	 hdr.nprops++;
       }
       if(bad_n < nobs) {
         hdr.prop_id[++i] = (int)N2;
	 data.observ[(int)N2] = n;
 	 hdr.nprops++;
      }
       if(bad_phos < nobs) {
         hdr.prop_id[++i] = (int)P4;
   	 data.observ[(int)P4] = phos;
	 hdr.nprops++;
       }
       if(bad_sil < nobs) {
         hdr.prop_id[++i] = (int)SI;
	 data.observ[(int)SI] = sil;    
	 hdr.nprops++;
       } 
       data.nprops = hdr.nprops; 
       if(!is_error)
         write_hydro_station(outfile, &hdr, &data);
       else
         write_hydro_station(errorfile, &hdr, &data);

       free(hdr.prop_id);   
       free(p);
       free(d);
       free(t);
       free(s);
       free(o);
       free(n);
       free(phos);
       free(sil);
     }

     error = fclose(infile);
     if(error) {
       fprintf(stderr, "\nError closing infofile\n");
       exit(1);
     }
NEXTFILE: 
     fprintf(stderr,"\nCompleted reading %s\n", argv[curfile]);

   } while (curfile++ < nfiles );

   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nUsage:  %s filelist -Ooutfile [-D<dirname>]", program);
   fprintf(stderr,"\n\n  List of filenames MUST be first argument!");
   fprintf(stderr,"\n    -D  : specifies dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n    -O  : specifies output filename");

   fprintf(stderr,"\n\n");  
   return;
}
   
/*****************************************************************************/
FILE *openfile(char *dir, char *root)
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
   infile = fopen(st,"r");
   if (infile == NULL)
         fprintf(stderr,"\n Unable to open %s \n", st);
   else
         fprintf(stderr,"\n   opened %s ...\n", st);
   
   return(infile);
   
}  /* end openfile() */
/*****************************************************************************/
int parse_header(FILE *fptr)

   /* reads the header info from an already opened btl data file 
     The info is used to fill in appropriate values in
      the struct HYDRO_HDR hdr, a global variable.  */
{
  int error, niskin, nscan;
  char line[1000], st1[20], st2[20], st3[20], st4[20], st5[20], st6[20], st7[20], st8[20], st9[20], st10[20], st11[20], st12[20], st13[20], st14[20], st15[20], st16[20], s17[20], st18[20], st19[20], st20[20], st21[20], st22[20], st23[20], st24[20], st25[20];
  double de_dat, te_dat, sa_dat, o2_dat, n2_dat, p4_dat, si_dat;

  nscan = sscanf(line_buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s", st1, st2, st3, st4, st5, st6, st7, st8, st9, st10, st11, st12, st13, st14, st15, st16, s17, st18, st19, st20, st21, st22, st23, st24, st25); 
   
  if (nscan >= 14 && !feof(fptr)) {

   /*parsing strings*/
     fprintf(stderr, "\nReading cast id# %s\n", st1);
     error = sscanf(st1, "%5d%2d%2d", &hdr.cruise, &hdr.station, &niskin);
     error = sscanf(st2, "%4d%2d%2d", &hdr.year, &hdr.month, &hdr.day);
     error = sscanf(st5, "%f", &hdr.lat);
     error = sscanf(st6, "%f", &hdr.lon);
     if (hdr.lon >0 )
       hdr.lon *= -1; /* longitude shoud be negative */

     return 1;
   }
   else
     return -1;
}  /* end parse_header() */
/*****************************************************************************/
int readdata(FILE *fptr)
{
  int nobs, i, error, propcount, nread, scan_cruise, scan_station, niskin, new_cast = 0;
   double de_dat, te_dat, sa_dat, o2_dat, n2_dat, p4_dat, si_dat;
   char st1[20], st2[20], st3[20], st4[20], st5[20], st6[20], st7[20], st8[20], st9[20], st10[20], st11[20], st12[20], st13[20], st14[20], st15[20], st16[20], s17[20], st18[20], st19[20], st20[20], st21[20], st22[20], st23[20], st24[20], st25[20];

   nobs = 0;  
   /*loop through the rest of the station*/
    do{
      nread = sscanf(line_buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s", st1, st2, st3, st4, st5, st6, st7, st8, st9, st10, st11, st12, st13, st14, st15, st16, s17, st18, st19, st20, st21, st22, st23, st24, st25);
      sscanf(st1, "%5d%2d%2d",&scan_cruise, &scan_station, &niskin);
      if((scan_cruise != hdr.cruise) || (scan_station != hdr.station))
	return nobs;
      error = sscanf(st7, "%lf", &de_dat); 
      error = sscanf(st8, "%lf", &te_dat);
      error = sscanf(st10, "%lf", &sa_dat);
      if (sa_dat < -8)
      error = sscanf(st9, "%lf", &sa_dat);
      error = sscanf(st12, "%lf", &o2_dat);
      d[nobs] = de_dat;
      p[nobs] = hb_p80(d[nobs], hdr.lat);
      t[nobs] = te_dat;
      s[nobs] = sa_dat;
      o[nobs] = o2_dat;
      if (t[nobs] > -8) {
        t[nobs] *= 1.00024; /* all temperatures are ITS-90 */
      }
      if(nread == 25) {
        error = sscanf(st18, "%lf", &n2_dat);
        error = sscanf(st19, "%lf", &p4_dat);
        error = sscanf(st20, "%lf", &si_dat);
	n[nobs] = n2_dat;
	phos[nobs] = p4_dat;
	sil[nobs] = si_dat;
      }
      else {
        n[nobs] = MISSING;
        phos[nobs] = MISSING;
        sil[nobs] = MISSING;
      }
      ++nobs;
    } while ((fscanf(fptr, "%[^\n]", line_buffer) != -1) && (getc(fptr) != EOF));
   return (nobs);
}  /* end readdata() */
