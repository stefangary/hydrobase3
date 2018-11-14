/*  wod_convert.c
................................................................................
.   Reads World Ocean Database  observed level hydrographic data files
   UPDATED to distinguish IPTS-68 and ITS-90 temperature scales
.   extracts :    header info
.                 p,d,t,s,ox,n2,n3,si,p4
.   
.   USAGE: wod_convert infile_list -B<badfile> -O<outfile> -T<instr_type> [-D<dir>] [-E<extent>]
...............................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_paths.h"

#define   BLANK     0x20   /* ascii code for blank */
#define   MISSING   -9.9   /* missing value flag */
#define   DIR    ""
#define   EXTENT ""



/* definitions and global vars from NODC's wodC.c program */

#define maxtax 25
#define MAXSEC 100
#define maxsec 100
#define MAXBIO 50
#define maxbio 50
#define MAXPARM 100
#define maxparm 100
#define maxpsec 25 * maxparm

FILE *fp;   
char cc[2];
int icruise=0, ostation=0, year=0, month=0, day=0;
int hour, longitude, latitude;
int levels, isoor, nparm, ip2[MAXPARM], iperror[MAXPARM];
int htotfig[3], hsigfig[3], hrightfig[3];
int origcfig, origsfig;
char origc[30], origs[30];
int ipip[MAXPARM], ipi[MAXPARM], npi;
int nsec;
int stotfig[MAXSEC], ssigfig[MAXSEC], srightfig[MAXSEC];
int seccode[MAXSEC], secval[MAXSEC];
int npsec;
int pstotfig[maxpsec],pssigfig[maxpsec],psrightfig[maxpsec];
int psecparm[maxpsec],pseccode[maxpsec],psecval[maxpsec];
int nbio;
int btotfig[MAXBIO], bsigfig[MAXBIO], brightfig[MAXBIO];
int biocode[MAXBIO], bioval[MAXBIO];
int ntsets;
int *ntloc,*ntcode,*ntval,*nterr,*ntoerr,*nttotfig,*ntsigfig,*ntrightfig;
int *depth,*zerr,*zoerr,*ztotfig,*zsigfig,*zrightfig;
int *dataval,*derr,*doerr,*dtotfig,*dsigfig,*drightfig;
int isize, zsize;
int ntsetsmax=0, isizemax=0, zsizemax=0;
float tenp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
int sdepth[] = { 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250,
                  300, 400, 500, 600, 700, 800, 900, 1000, 1100,
                  1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000,
                  3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000,
                  7500, 8000, 8500, 9000 };
int maxscode;
char **nodc_ship;   /* shipcode table indexed by OCL shipcode */

  /*  HydroBase globally defined variables */
  
struct HYDRO_HDR  hdr, bad_hdr;
struct HYDRO_DATA data, bad_data;

  /*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);
int checksta();
int check_profile_flag(int);
int check_obs_flag(int);
int oclread();
void spacer(int);
int extractc(int, int *, char * );
int extracti(int, int *, int *, int *,  int *, int);            


main (int  argc, char **argv)
{
   int error, ropt, topt;
   int staread, staout, stabad;
   int  i, j, curfile = 1, nfiles = 0;
   char buf[200];   /* to read in shipcode data */ 
   char *dir, *extent, *st ;
   FILE *shipfile;
   FILE *outfile, *badfile;
   
/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }
   
 
/*  set these default values */

    dir = DIR;
    extent = EXTENT;
    error = 0;
    staout = staread = stabad = 0;
    outfile = stdout;
    topt = ropt = 0;
    shipfile = NULL;
    maxscode = 9999;
    
    for (i = 0; i < MAXPROP; ++i) {
       data.observ[i] = NULL;
       bad_data.observ[i] = NULL;
   }
    hdr.prop_id = NULL;
    bad_hdr.prop_id = NULL;

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
                        fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
                        break;
               case 'R':                    /* open  file for rejected profiles */
                        badfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                        if (badfile == NULL) {
                           fprintf(stderr,"\nError opening output file: %s\n", 
                                  &argv[i][2]);
                            fprintf(stderr,"Does file already exist?\n"); 
                          exit(1);
                        }
                        fprintf(stderr,"\nRejected data will be written to: %s", &argv[i][2]);

			ropt = 1;
                        break;
               case 'S':                    /* get shipcode file  */
                        shipfile = fopen(&argv[i][2], "r");
                        if (shipfile == NULL) {
                           fprintf(stderr,"\nError opening ship file: %s\n", 
                                  &argv[i][2]);
                           exit(1);
                        }
                        break;

               case 'T':                    /* get datatype  */
                      switch (argv[i][2]) {
                        case 'c':  /* fall through */
                        case 'b':
                        case 'f':
                        case 's':
                        case 'u':
                        case 'm':
                             hdr.instrument = argv[i][2];
                             bad_hdr.instrument = argv[i][2];
                            break;
                        default:
                            fprintf(stderr,"\nError in datatype specification: %s\n", argv[i]);
                            fprintf(stderr,"Legitimate choices: \n");
                            fprintf(stderr,"   [c]td  \n");
                            fprintf(stderr,"   [b]ottle \n");
                            fprintf(stderr,"   [f]loat  \n");
                            fprintf(stderr,"   [s]easoar  \n");
                            fprintf(stderr,"   [m]oored profiler  \n");
                            fprintf(stderr,"   [u]nknown  \n\n");
                            exit(1);
                        }
			topt = 1;
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
   
   if (! topt) {
       fprintf(stderr,"\nYou must specify a data type: c,b,f,s,m,or u\n");
       exit(1);
   }

   if (! ropt) {
       fprintf(stderr,"\nBad stations will not be saved in a separate file\n");
   }
   
   if (! nfiles) {
       fprintf(stderr,"\nExpecting data to come from STDIN...");
       fp = stdin;
   }
 
 
 /* read in file shipcode.asc which is derived from the OCL supplied code
    file SHIPCODE.TXT */
    
    if (shipfile == NULL) {
        
       shipfile = fopen(SHIPCODE_PATH, "r");
       if (shipfile == NULL) {
          fprintf(stderr,"\nUnable to open default shipcode file: %s \n",SHIPCODE_PATH );
          exit(1);
       }
    }

/* initialize shipcode table */
    nodc_ship = (char **) calloc (maxscode, sizeof(char *));
    if (nodc_ship == NULL) {
          fprintf(stderr,"Unable to allocate memory for shipcode table\n");
          exit(1);
       }
    for (i = 0; i < maxscode; ++i) {
       nodc_ship[i] = (char *) calloc(4, sizeof(char));
       if (nodc_ship[i] == NULL) {
          fprintf(stderr,"Unable to allocate memory for shipcode table\n");
          exit(1);
       }
    }
        
    fprintf(stderr,"\nReading shipcode file.... ");
    while ((error = fscanf(shipfile,"%[^\n]", buf)) == 1) {
           getc(shipfile);  /* move past LF */
        
          if (sscanf(buf,"%d", &i) == 1) {  /* first col is index to shipcode table*/
             if (i >= 0  && i < maxscode) {
                 st = &buf[5];
	         while (*st == ' ')
	             ++st;
                 for (j = 0; j < 4; ++j) {  /*  store four chars at index in table */
                    nodc_ship[i][j] = *st;
                    ++st;
                 }
             }
          }
    }
      
   
   /* initialize dynamic arrays */
   
  spacer(1);

/********************************************************

 ENTER STANDARD LEVEL DEPTHS, IN CASE THIS IS STANDARD LEVEL DATA

********************************************************/

  for ( j = 0; j < 40; j++ ) 
       *(depth+j)= *(sdepth+j);

  /* loop for each file ...*/   
  
  do {
  
     if (nfiles) {
  	fp = openfile(dir, argv[curfile], extent);
  	if (fp == NULL) 
        	goto NEXTFILE;
     }
     
     /* loop for each station */
     
      
     while ( oclread() >= 0 && !feof(fp) ) {
      
           ++staread;
        
	
	   if (check_sta() > 0) {        
              write_hydro_station(outfile, &hdr, &data);
              ++staout;
	   }
	   
	   if (ropt && (bad_hdr.nobs > 0)) {
	      write_hydro_station(badfile, &bad_hdr, &bad_data);
              ++stabad;
	   }
        
       /* clean up ... */
       
            if (hdr.prop_id != NULL) {
               free((void *)hdr.prop_id);
               hdr.prop_id = NULL;
            } 
                 
            if (bad_hdr.prop_id != NULL) {
               free((void *)bad_hdr.prop_id);
               bad_hdr.prop_id = NULL;
            } 
            
            for (i=0; i< MAXPROP; ++i) {
               if (data.observ[i] != NULL) {
                  free((void *) data.observ[i]);
                  data.observ[i] = NULL;
               }
               if (bad_data.observ[i] != NULL) {
                  free((void *) bad_data.observ[i]);
                  bad_data.observ[i] = NULL;
               }
            }  
       }  /* end while */
       
NEXTFILE:
     if (nfiles) 
       fclose(fp);
       
    }while (curfile++ < nfiles);       

   fprintf(stderr,"\nEnd of conversion.");
   fprintf(stderr,"\n  %d stations read in", staread);
   fprintf(stderr,"\n  %d stations accepted", staout);
   fprintf(stderr,"\n  %d stations rejected ", staread-staout);
   fprintf(stderr,"\n  %d stations contained bad scans \n\n", stabad - (staread-staout));
   exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nConverts .NBD .BD .CTD and .NCT files (observed levels)from WOD98 to HydroBase\n");
   fprintf(stderr,"\nUsage:  %s filelist  [-Ooutfile] [-D<dirname>] [-E<file_extent>] [-S<shipcode_file>] [-h]", program);
   fprintf(stderr,"\n\n  List of filenames must be first argument");
   fprintf(stderr,"\n  Otherwise, input is expected from stdin.");
   fprintf(stderr,"\n-T : specifies data_type code [b,c,f,s,m, or u]");
   fprintf(stderr,"\n\n   OPTIONS: ");
   fprintf(stderr,"\n-D : dirname for input files (default is current directory) ");
   fprintf(stderr,"\n            ex: -D../data/ ");
   fprintf(stderr,"\n-E : input_file extent (default is no extent)");  
   fprintf(stderr,"\n            ex: -E.dat ");
   fprintf(stderr,"\n-O : output file (overwrites an existing file)");
   fprintf(stderr,"\n      defaults to stdout " );
   fprintf(stderr,"\n-R : specifies file for rejected stations (noclobber)");
   fprintf(stderr,"\n-S : path/filename for OCL/NODC shipcode table.");
   fprintf(stderr,"\n      defaults to %s ", SHIPCODE_PATH );
   fprintf(stderr,"\n-h : help...... prints this message. ");

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
         fprintf(stderr,"\n   opened %s ... ", st);
   
   return(infile);
   
}  /* end openfile() */

/*****************************************************************************/
int check_sta()
   /* Returns 1 if the station has acceptable data, 0 if not */
{
   int *prop_avail, *prop_OK, *dindex;
   int **data_flagged;
   int offset, i, j,ii, nobs, delta;
   int found, index, tindex;
   float prsint;
   
   /* allocate some space */
   prop_avail = (int *) calloc(MAXPROP, sizeof(int));
   prop_OK = (int *) calloc(MAXPROP, sizeof(int));
   dindex = (int *) calloc(MAXPROP, sizeof(int));
   data_flagged = (int **) calloc(MAXPROP, sizeof(int *));
   
   /* store metadata in HydroBase station header ... */
      
   hdr.lat = (float) latitude / tenp[*(hrightfig+1)];
   bad_hdr.lat = hdr.lat;
   hdr.lon = longitude/ tenp[ *(hrightfig+2)];
   bad_hdr.lon = hdr.lon;
   hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
   bad_hdr.ms10 = hdr.ms10;
   bad_hdr.ms1 = hdr.ms1;
   hdr.origin = '1';
   bad_hdr.origin = '1';
   strncpy(hdr.country, cc, 2);
   hdr.country[2] = '\0'; 
   hdr.year = year;
   bad_hdr.year = year;
   hdr.month = month;
   bad_hdr.month = month;
   hdr.day = day;
   bad_hdr.day = day;
   
   hdr.cruise =  icruise - ((icruise/100000) * 100000);
   hdr.station = ostation - ((ostation/10000) * 10000);
   
   /* use originator's cruise and station #s if available */
   
   *(origc+origcfig) = '\0';  /* set end-of-string delimiters */
   *(origs+origsfig) = '\0';
   if (origcfig > 0) {      
      if ((i = atoi(origc)) > 0)
         hdr.cruise = i - ((i/100000) * 100000);
   }
   
   if (origsfig > 0) {
      if ((i = atoi(origs)) > 0)
         hdr.station = i - ((i/10000) * 10000);
   }
   
   bad_hdr.cruise = hdr.cruise;
   bad_hdr.station = hdr.station;
   
   /*search for a platform code and PDR depth in secondary header arrays  */
   
   found = 0;
   i = 0;
   while (!found && (i < nsec)) {
     /* secondary header code for OCL Platform is 3
      * so if seccode[i] == 3, then found is assigned
      * a value of true, otherwise, stays at 0. */
      found = seccode[i] == 3;
      ++i;
   }
   
   strncpy(hdr.ship, "XX", 3);
   if (found) {
      --i;   
      index = secval[i] / tenp[srightfig[i]];
      if (index >= maxscode) {
         fprintf(stderr,"\nError parsing index to shipcodes.asc table: %d\n",index);
         exit(1);
      }
      if (index > 0 && (nodc_ship[index][2] != '\0')) {
         hdr.country[0] = nodc_ship[index][0];
         hdr.country[1] = nodc_ship[index][1];
         hdr.country[2] = '\0'; 
	 
         hdr.ship[0] = nodc_ship[index][2]; 
         hdr.ship[1] = nodc_ship[index][3]; 
         hdr.ship[2] = '\0';
      }
   }
   
   strncpy(bad_hdr.ship, hdr.ship, 3);
   strncpy(bad_hdr.country, hdr.country, 3);
   
   found = 0;
   i = 0;
   while (!found && (i < nsec)) {
      found = seccode[i] == 10;    /* secondary header code for seafloor depth */
      ++i;
   }

   hdr.pdr = bad_hdr.pdr = 0;   
   if (found) {
      --i;   
      hdr.pdr = secval[i] / tenp[srightfig[i]]; 
      bad_hdr.pdr = hdr.pdr; 
   }
   
   
    
   /* Determine temperature scale */
   
     tindex = (int)TE;
     found = 0;
     ii = 0;
     while (!found  && (ii < npsec)) {
        if ( (psecparm[ii] == 1) && (pseccode[ii] == 3)) {  /* these are OCL codes for Temperature and code table 3, respectively */
	   found = 1;
	   if (psecval[ii] == 103)  tindex = (int)T90;    /* OCL code for T90 == 103, for T68 == 102 */
	}
	++ii;
     }
   
  /* next check to see what properties are available */
     
   for (i = 0; i < nparm; ++i) {
      switch (ip2[i]) {
	      case 1:
	         prop_avail[tindex] = 1;	         
	         prop_OK[tindex] = check_profile_flag(iperror[i]);
		 dindex[tindex] = i;
	         break;
	      case 2:
                 prop_avail[(int)SA] =1;
	         prop_OK[(int)SA] = check_profile_flag(iperror[i]);
		 dindex[(int)SA] = i;
	         break;
	      case 3:
	         prop_avail[(int)OX] = 1;
	         prop_OK[(int)OX]  = check_profile_flag(iperror[i]);
		 dindex[(int)OX] = i;
	         break;
	      case 4:
	         prop_avail[(int)P4] = 1;
	         prop_OK[(int)P4] = check_profile_flag(iperror[i]);
		 dindex[(int)P4] = i;
	         break;
	      case 6:
	         prop_avail[(int)SI] = 1;
	         prop_OK[(int)SI] = check_profile_flag(iperror[i]);
		 dindex[(int)SI] = i;
	         break;
	      case 7:
	         prop_avail[(int)N2] = 1;
	         prop_OK[(int)N2] = check_profile_flag(iperror[i]);
		 dindex[(int)N2] = i;
	         break;
	      case 8:
	         prop_avail[(int)N3] = 1;
	         prop_OK[(int)N3] = check_profile_flag(iperror[i]);
		 dindex[(int)N3] = i;
	         break;
	      case 25:
	         prop_avail[(int)PR] = 1;
	         prop_OK[(int)PR] = check_profile_flag(iperror[i]);
		 dindex[(int)PR] = i;
	         break;
	      default:
	            ;
	   
      } /* end switch */
   } /* end for */
   
   /* load property data into HydroBase station structures */
   
   hdr.nprops = 0;	
   bad_hdr.nprops = 0;	
   
   for (i = 0; i < MAXPROP; ++i) {
      if (prop_avail[i] ) {
	 data.observ[i] = (double *) malloc(levels * sizeof(double));
	 bad_data.observ[i] = (double *) malloc(levels * sizeof(double));
	 data_flagged[i] = (int *) calloc(levels, sizeof(int));
	 
	 for (j=0; j < levels; ++j) {
	    offset = dindex[i]*levels;
	    data.observ[i][j] = (double) (*(dataval+offset+j) /tenp[*(drightfig+offset+j)]);
	    data_flagged[i][j] = check_obs_flag(*(derr+offset+j)) ? 0 : 1;
	    if (data.observ[i][j] < -9.) 
	       data_flagged[i][j] = 1;
 	 }
 	 ++bad_hdr.nprops;
         if ( prop_OK[i])
             ++hdr.nprops;
      }
   }
   
	   
   /* load depths */
   	   
   prop_avail[(int)DE] = 1;  /* by definition */
   prop_OK[(int)DE] = 1;  
   ++hdr.nprops;
   ++bad_hdr.nprops;
   data.observ[(int)DE] = (double *)malloc (levels * sizeof(double));
   bad_data.observ[(int)DE] = (double *)malloc (levels * sizeof(double));
   data_flagged[(int)DE] = (int *) calloc(levels, sizeof(int));
   for (j=0; j < levels; ++j) {
      data.observ[(int)DE][j] = (double) (depth[j] / tenp[zrightfig[j]]);
      data_flagged[(int)DE][j] = check_obs_flag(zerr[j]) ? 0 : 1;
     if (data.observ[(int)DE][j] < -9.) 
         data_flagged[(int)DE][j] = 1;
   }  
    
    /* if pressure wasn't a stored parameter, compute it... */
   
   if (!prop_avail[(int)PR]) {
      data.observ[(int)PR] = (double *)malloc (levels * sizeof(double));
      bad_data.observ[(int)PR] = (double *)malloc (levels * sizeof(double));
      data_flagged[(int)PR] = (int *) calloc(levels, sizeof(int));
      for (j=0; j<levels; ++j) {
         data.observ[(int)PR][j] = hb_p80(data.observ[(int)DE][j], (double)hdr.lat);
         data_flagged[(int)PR][j] = data_flagged[(int)DE][j];
      }
      prop_avail[(int)PR] = 1;
      prop_OK[(int)PR] = 1;  
      ++hdr.nprops;
      ++bad_hdr.nprops;
   } 
    
   /* check for requisite parameters... */
      
   if ( ! (prop_OK[tindex] && prop_OK[(int)SA])) {
   
           /* just write station in bad file */
           
      for (i=0; i< MAXPROP; ++i) {
         for (j=0; j < levels; ++j) {
            if (prop_avail[i])
               bad_data.observ[i][j] = data.observ[i][j];
         }
      }
          
      bad_hdr.prop_id = (int *) malloc(bad_hdr.nprops * sizeof(int));
      j = 0;
      for (i=0; i< MAXPROP; ++i) {
         if (prop_avail[i]) {
            bad_hdr.prop_id[j++] = i;
            free((void *)data_flagged[i]);
         }
      }
      bad_hdr.nobs = bad_data.nobs = levels; 
      bad_data.nprops = bad_hdr.nprops;
      free((void *)data_flagged);
      free((void *)prop_avail);
      free((void *)prop_OK);
      free((void *)dindex);
      return(0);        
   }
      
   
   /* Tally number of good/bad levels. A bad level = p,t,or s flagged 
      in which case the level is moved to bad_data structure.  */
   
   hdr.nobs = 0;  
   bad_hdr.nobs = 0; 
   for (j=0; j<levels; ++j) {
   
      if (data_flagged[(int)DE][j] || data_flagged[(int)PR][j] || data_flagged[tindex][j] || data_flagged[(int)SA][j] || data.observ[(int)PR][j] < 0 ) {
         for (i=0; i< MAXPROP; ++i) {
            if (prop_avail[i])
               bad_data.observ[i][bad_hdr.nobs] = data.observ[i][j];
         }
         ++bad_hdr.nobs;
      } 
      else  {
         for (i=0; i< MAXPROP; ++i) {
            if (prop_OK[i]) {
               data.observ[i][hdr.nobs] = data.observ[i][j];
               if (data_flagged[i][j])
                  data.observ[i][hdr.nobs] = MISSING;
            }
         }
         ++hdr.nobs;
      }
      
   }  /*end for */ 
   data.nobs = hdr.nobs ;  
   bad_data.nobs = bad_hdr.nobs; 
   
   
   hdr.prop_id = (int *) malloc(hdr.nprops * sizeof(int));
   bad_hdr.prop_id = (int *) malloc(bad_hdr.nprops * sizeof(int));
   j = 0;
   i = 0;
   for (index = 0; index < MAXPROP; ++index) {
      if (prop_OK[index]) {
         hdr.prop_id[j] = index;
         if (++j > hdr.nprops) {
           fprintf(stderr,"\nError counting hdr.nprops in check_sta()!\n");
           exit(1);
         }
      }
      if (prop_avail[index]) {
         bad_hdr.prop_id[i] = index;
         free((void *)data_flagged[index]);
         if (++i > bad_hdr.nprops) {
           fprintf(stderr,"\nError counting bad_hdr.nprops in check_sta()!\n");
           exit(1);
         }
      }
   }
   
   
   /*check interval of pressure series .. decimate to 10 db if delta-p is
   smaller */
   
   
   if (hdr.nobs > 30) {
   
      prsint = 0;      /* determine avg prs int */
      for (i = 15; i < 30; ++i) {
         prsint += (float)(data.observ[(int)PR][i+1] - data.observ[(int)PR][i]);
      }
      prsint = prsint / 15;
      
      if (prsint < 10) {
         delta = (int)(10 / prsint);
         for (i = 0; i < hdr.nprops; ++i) {
            nobs = 0;
            for (j = 0; j < hdr.nobs; j += delta) {
               data.observ[hdr.prop_id[i]][nobs] = data.observ[hdr.prop_id[i]][j];
               ++nobs;
            }
            if ((j - delta) < (hdr.nobs-1)) {
               data.observ[hdr.prop_id[i]][nobs] = data.observ[hdr.prop_id[i]][hdr.nobs-1];
               ++nobs;
            }
         }
         
         hdr.nobs = data.nobs = nobs;
      }
   }
      
   free((void *)prop_avail);
   free((void *)prop_OK);
   free((void *)dindex);
   free((void *)data_flagged);
   
   return (hdr.nobs) ;

}  /* end check_sta() */
/*****************************************************************************/
int check_profile_flag(int flag)
{
   switch (flag) {
   
      case 0:
      case 1:  /* fall through */ 
      case 2:
           return 1;
      default:
           return 0;
   } 

} /* end check_profile_flag() */
/*****************************************************************************/
int check_obs_flag(int flag)
{
   switch (flag) {
   
      case 0:   /* fall through */
      case 2:
      case 3:
           return 1;
      default:
           return 0;
   } 

} /* end check_obs_flag() */

/*****************************************************************************/
/*****************************************************************************/
int oclread()
  /* This code was directly imported from the wesite distributing WOD2009 */
{

 int i,j;

 char wodform;
 int totfig=0, sigfig=0, rightfig=0;
 int nbytet,ntypec,nbytec,nbytes,nbyteb;
 int ninfc,ntoff, doff;
 int missing=-9999;
 int npinfs=0,npinfe=0,npinf;
 int iend=0;

/**********************************************************

 READ IN WOD FORMAT CODE: 'A' FOR WOD05 FORMAT

***********************************************************/

 totfig= 1;
 if ( (iend = extractc(0,&totfig,&wodform)) == -1 ) return iend;


/**********************************************************

 READ IN NUMBER OF BYTES IN ASCII STATION

***********************************************************/

 if ( ( iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytet,missing))
          == -1 ) return iend;

/**********************************************************

 READ IN OCL STATION NUMBER

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ostation,missing))
       == -1 ) return iend;

/*********************************************************

 READ IN NODC COUNTRY CODE

**********************************************************/

 totfig= 2;
 if ( (iend = extractc(0,&totfig,cc)) == -1 ) return iend;

/**********************************************************

 READ IN CRUISE NUMBER

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&icruise,
             missing)) == -1) return iend;

/**********************************************************

 READ IN YEAR, MONTH, DAY, TIME

***********************************************************/

 totfig=4;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&year,missing))
       == -1) return iend;
 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&month,missing))
       == -1) return iend;
 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&day,missing))
       == -1) return iend;

 if ( (iend = extracti(1,htotfig,hsigfig,hrightfig,
              &hour, 9999)) == -1) return iend;

/**********************************************************

 READ IN LATITUDE AND LONGITUDE

***********************************************************/

 if ( (iend = extracti(1,(htotfig+1),(hsigfig+1),(hrightfig+1),
              &latitude, -9999)) == -1) return iend;

 if ( (iend = extracti(1,(htotfig+2),(hsigfig+2),(hrightfig+2),
              &longitude, -99999)) == -1) return iend;

/**********************************************************

 READ IN NUMBER OF LEVELS

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&levels,
             missing)) == -1) return iend;

/**********************************************************

 READ IN OBSERVED (0) OR STANDARD (1) LEVELS

***********************************************************/

 totfig=1;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&isoor,missing))
       == -1) return iend;

/**********************************************************

 READ IN NUMBER OF PARAMETERS AT THIS STATION (NOT INCLUDING
 BIOLOGY)

***********************************************************/

 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&nparm,missing))
       == -1) return iend;

/**********************************************************

 READ IN EACH PARAMETER CODE AND PROFILE ERROR CODE

***********************************************************/

 for ( i =0; i< nparm; i++ ) {

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ip2+i),
             missing)) == -1) return iend;
  totfig=1;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(iperror+i),
             missing)) == -1) return iend;

/*******************************************************************

 READ NUMBER OF PARAMETER SPECIFIC SECOND HEADER VARIABLES

********************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&npinf,missing))
       == -1 ) return iend;
  npinfe += npinf;

/*******************************************************************

 READ IN EACH PARAMETER SPECIFIC SECOND HEADER VARIABLE

********************************************************************/

  for ( j = npinfs; j < npinfe; j++ ) {

   *(psecparm+j) = *(ip2+i);
   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(pseccode+j),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(pstotfig+j),(pssigfig+j),(psrightfig+j),
                        (psecval+j), missing)) == -1) return iend;

  }

  npinfs += npinf;

 }

 npsec = npinfe;

/***************************************************************

 READ IN NUMBER OF BYTES IN CHARACTER AND PRIMARY INVESTIGATOR FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytec,
             missing)) == -1) return iend;

/**********************************************************

 READ IN NUMBER OF INFORMATION TYPES (MAX 3: FOR CRUISE CODE,
 STATION CODE AND PI INFORMATION)

***********************************************************/

 origcfig= 0;
 origsfig= 0;
 npi= 0;

 if ( nbytec > 0 ) {

  totfig=1;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ninfc,missing))
       == -1) return iend;

/**********************************************************

 READ IN TYPE OF INFORMATION : 1= CRUISE CODE,
 2=STATION CODE AND 3=PI INFORMATION)

***********************************************************/

  for ( i= 0; i < ninfc; i++) {

   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ntypec,missing))
       == -1) return iend;

/***********************************************************

 READ IN ORIGINATORS CRUISE CODE

************************************************************/

   if ( ntypec == 1) {

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origcfig,
     missing)) == -1) return iend;
    if ( (iend = extractc(0,&origcfig,origc)) == -1 ) return iend;

   }

/***********************************************************

 READ IN ORIGINATORS STATION CODE

************************************************************/

   else if ( ntypec == 2 ) {

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origsfig,
     missing)) == -1) return iend;
    if ( (iend = extractc(0,&origsfig,origs)) == -1 ) return iend;

   }

/***********************************************************

 READ IN PRIMARY INVESTIGATOR INFORMATION

************************************************************/

   else if ( ntypec == 3 ) {

/**********************************************************

 READ IN NUMBER OF PRIMARY INVESTIGATORS AT THIS STATION

***********************************************************/

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&npi,missing))
         == -1) return iend;

/**********************************************************

 READ IN EACH PARAMETER CODE AND PRIMARY INVESTIGATOR CODE
 FOR THAT PARAMETER 

***********************************************************/

    for ( j =0; j< npi; j++ ) {

     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipip+j),
              missing)) == -1) return iend;
     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipi+j),
              missing)) == -1) return iend;

    }

   }

  }

 }

/***************************************************************

 READ IN NUMBER OF BYTES IN SECONDARY HEADER FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytes,
             missing)) == -1) return iend;

 nsec = 0;
 if ( nbytes > 0 ) {

/**************************************************************

 READ IN NUMBER OF SECONDARY HEADER PARAMETERS

***************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nsec,
              missing)) == -1) return iend;

/**********************************************************

 READ IN EACH SECONDARY HEADER PARAMETER CODE AND VALUE

***********************************************************/

  for ( i =0; i< nsec; i++ ) {

   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(seccode+i),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(stotfig+i),(ssigfig+i),(srightfig+i),
                        (secval+i), missing)) == -1) return iend;

  }
 
 }

/***************************************************************

 READ IN NUMBER OF BYTES IN BIOLOGY HEADER AND TAXA SET FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbyteb,
             missing)) == -1) return iend;

 nbio= 0;
 ntsets= 0;
 if ( nbyteb > 0 ) {

/**************************************************************

 READ IN NUMBER OF BIOLOGY HEADER PARAMETERS

***************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbio,
              missing)) == -1) return iend;

/**********************************************************

 READ IN EACH BIOLOGY HEADER PARAMETER CODE AND VALUE

***********************************************************/

  for ( i =0; i< nbio; i++ ) {

   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(biocode+i),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(btotfig+i),(bsigfig+i),(brightfig+i),
                        (bioval+i), missing)) == -1) return iend;

  }
  
/***************************************************************

 READ IN NUMBER OF TAXA SETS

****************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ntsets,
             missing)) == -1) return iend;

/***************************************************************

 IF NUMBER OF TAXA SETS IS GREATER THAN THE PREVIOUS MAXIMUM,
 REALLOCATE SPACE ACCORDINGLY

****************************************************************/

  if ( ntsets > ntsetsmax ) {

   ntsetsmax= ntsets;
   spacer(2);

  }

  if ( ntsets > 0 ) {

/**************************************************************

 READ IN NUMBER OF ENTRIES FOR THIS TAXA SET, CALCULATE OFFSET
 FOR READING IN DATA

***************************************************************/
  
   for ( j = 0; j < ntsets; j++) {

    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntloc+j),
               missing)) == -1) return iend;

    ntoff= maxtax * j;

/**********************************************************

 READ IN EACH TAXA SET PARAMETER CODE AND VALUE,
 ERROR FLAG AND ORIGINATORS FLAG

***********************************************************/

    for ( i =0; i< *(ntloc+j); i++ ) {

     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntcode+ntoff+i),
              missing)) == -1) return iend;
     if ( (iend = extracti(1,(nttotfig+ntoff+i),(ntsigfig+ntoff+i),
                 (ntrightfig+ntoff+i), (ntval+ntoff+i), missing))
                  == -1) return iend;
     totfig=1;
     if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(nterr+ntoff+i),
              missing)) == -1) return iend;
     totfig=1;
     if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(ntoerr+ntoff+i),
              missing)) == -1) return iend;

    }

   }
  
  }

 }

/***************************************************************

 IF NUMBER OF DEPTHS IS GREATER THAN THE PREVIOUS MAXIMUM,
  REALLOCATE SPACE ACCORDINGLY

****************************************************************/

 zsize= levels;
 if ( isoor == 0 && zsize > zsizemax ) {

  zsizemax=zsize;
  spacer(3);

 }

/***************************************************************

 IF NUMBER OF DEPTHS MULTIPLIED BY NUMBER OF PARAMETERS IS
 GREATER THAN THE PREVIOUS MAXIMUM, REALLOCATE SPACE ACCORDINGLY

****************************************************************/

 isize= nparm * levels;
 if ( isize > isizemax ) {

  isizemax=isize;
  spacer(4);

 }

/**********************************************************

 READ IN EACH DEPTH VALUE, ERROR FLAG, AND ORIGINATORS FLAG

***********************************************************/

 for ( j = 0; j < levels; j++ ) {

  if ( isoor == 0 ) {

   if ( (iend = extracti(1,(ztotfig+j),(zsigfig+j),
               (zrightfig+j), (depth+j), missing))
               == -1) return iend;
   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zerr+j),
            missing)) == -1) return iend;
   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zoerr+j),
            missing)) == -1) return iend;

  }

/**********************************************************

 READ IN EACH DATA VALUE

***********************************************************/

  for ( i =0; i< nparm; i++ ) {

   doff= i * levels;


   if ( (iend = extracti(1,(dtotfig+doff+j),(dsigfig+doff+j),
               (drightfig+doff+j), (dataval+doff+j), missing))
                == -1) return iend;

   if ( *(dtotfig+doff+j) > 0 ) {

    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(derr+doff+j),
             missing)) == -1) return iend;
    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(doerr+doff+j),
             missing)) == -1) return iend;
   }

   else {

    *(derr+doff+j)=0;
    *(doerr+doff+j)=0;
    *(drightfig+doff+j)=2;

   }


  }

 }

/***********************************************************

 READ TO END OF STATION

************************************************************/

 while ( ( i = fgetc(fp)) != '\n' && !feof(fp) );
 return iend;

}

/************************************************************

 SPACER.C SETS UP ORIGINAL SPACING FOR ALL DYNAMIC ARRAYS

*************************************************************/

void spacer(int intime)
 
/* int intime       SET TO ONE TO INITIALIZE ALL DYNAMIC ARRAYS,
                    SET TO TWO TO REDIMENSION TAXA ARRAYS,
                    SET TO THREE TO REDIMENSION DEPTH, 
                    SET TO FOUR TO REDIMENSION MEASURED PARAMETER ARRAYS
                 */

      

{

 if ( intime == 1 ) { 

/***************************************************************

 ALLOCATE SPACE FOR TAXA SET DATA

****************************************************************/

  ntsets=1;
  ntsetsmax=1;
  if ( (ntcode =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntloc =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nttotfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntsigfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntrightfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntval =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nterr =calloc( ntsets * maxtax, sizeof(int)) )
        == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntoerr =calloc( ntsets * maxtax, sizeof(int)) )
        == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
 
/***************************************************************

 ALLOCATE SPACE FOR DEPTH AND MEASURED PARAMETERS

****************************************************************/

  zsize= 40;
  zsizemax= 40;
  if ( (ztotfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zsigfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zrightfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (depth =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zerr =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zoerr =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
 
  isize=40;
  isizemax= 40;
  if ( (dtotfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (dsigfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (drightfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (dataval =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (derr =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (doerr =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
 
 }

/***************************************************************

 REALLOCATE SPACE FOR TAXA DATA

****************************************************************/

 else if ( intime == 2 ) {

   ntsetsmax = ntsets;
   if ( (ntcode =realloc( ntcode, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntloc =realloc( ntloc, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nttotfig =realloc( nttotfig, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntsigfig =realloc( ntsigfig, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntrightfig =realloc( ntrightfig, ntsets * maxtax * sizeof(int)))
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntval =realloc( ntval, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nterr =realloc( nterr, ntsets * maxtax * sizeof(int)) )
         == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntoerr =realloc( ntoerr, ntsets * maxtax * sizeof(int)) )
         == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);

 }

/***********************************************************

 REALLOCATE SPACE FOR DEPTH

************************************************************/

 else if ( intime == 3 ) {

  if ( (ztotfig =realloc( ztotfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #1 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zsigfig =realloc( zsigfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #2 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zrightfig =realloc( zrightfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #3 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (depth =realloc( depth, zsize * sizeof(int)) )
       == NULL )
    printf( " #4 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zerr =realloc( zerr, zsize * sizeof(int)) )
       == NULL )
    printf( " #5 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zoerr =realloc( zerr, zsize * sizeof(int)) )
       == NULL )
    printf( " #6 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");

 }

/************************************************************

 REALLOCATE SPACE FOR MEASURED PARAMETERS

*************************************************************/

 else if ( intime == 4 ) {

  if ( (dtotfig =realloc( dtotfig, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (dsigfig =realloc( dsigfig, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (drightfig =realloc( drightfig, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (dataval =realloc( dataval, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (derr =realloc( derr, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (doerr =realloc( derr, isize * sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
 
 }

}
/**********************************************************

                 FUNCTION EXTRACTI

 EXTRACT1 EXTRACTS AN INTEGER VALUE FROM FILE.  THERE ARE THREE
 TYPES OF EXTRACTIONS.

 THE FIRST IS ACCOMPLISHED BY EXTRACTING THE NUMBER OF BYTES IN
 THE INTEGER VALUE AND THEN USING THIS INFORMATION
 TO EXTRACT THE INTEGER VALUE ITSELF.  THIS IS FOR INTEGERS USED
 INTERNALLY TO READ THE DATA.

 THE SECOND TYPE EXTRACTS NUMBER OF SIGNIFICANT FIGURES, NUMBER
 OF TOTAL FIGURES, AND NUMBER OF FIGURES TO THE RIGHT OF THE
 DECIMAL BEFORE USING NUMBER OF TOTAL FIGURES TO EXTRACT THE ACTUAL
 VALUE.  IF THE NUMBER OF SIGNIGICANT FIGURES IS NOT AN INTEGER,
 BUT A NEGATIVE SIGN (-), THE DATA VALUE IS A MISSING VALUE.

 THE THIRD TYPE EXTRACTS A VALUE USING AN INPUT NUMBER OF
 FIGURES (TOTFIG)

***********************************************************/

extracti(

 int type,                    /* TYPE OF EXTRACTION:
                                0=INTERNAL INTEGER VALUE
                                1=OUTPUT DATA VALUE 
                                2= OUTPUT DATA VALUE FROM
                                   GIVEN NUMBER OF FIGURES */

 int *totfig,                  /* NUMBER OF FIGURES IN THE
                                 ASCII REPRESENTATION OF THE VALUE 
                                 BEING EXTRACTED */

 int *sigfig,                  /* NUMBER OF SIGNIFICANT FIGURES IN
                                 THE VALUE BEING EXTRACTED */

 int *rightfig,                /* NUMBER OF PLACES TO THE RIGHT OF
                                 THE DECIMAL IN THE VALUE BEING
                                 EXTRACTED */

 int *value,                   /* ACTUAL VALUE BEING EXTRACTED, IN
                                 INTEGER FORM */

 int missing                  /* MISSING VALUE MARKER */

           )

{

 int sign,j,i;

/********************************************************

 SKIP IF THIS IS END OF LINE CHARACTER

*********************************************************/

 if ( type != 2 ) i=nocrfgetc();

/********************************************************

 IF THIS IS THE END OF FILE, SET EOF TO ONE TO NOTIFY
 MAIN PROGRAM

*********************************************************/

 if ( feof(fp) ) return -1;

 else {

/*******************************************************

 IF THIS IS A MISSING VALUE (I='-'), SET VALUE ACCORDINGLY.

********************************************************/

  if ( i == '-' ) {

   *value = missing;
   *totfig = 0;
   *sigfig = 0;
   *rightfig = 0;
   return 0;

  }

  else {

/*********************************************************

 ELSE READ IN AND/OR SET NUMBER OF SIGNIFICANT FIGURES,
 NUMBER OF TOTAL FIGURES, AND NUMBER OF FIGURES RIGHT
 OF THE DECIMAL IF A TYPE 1, AND NUMBER OF TOTAL FIGURES
 IF TYPE 0.

**********************************************************/

   if ( type == 1 ) {

    *sigfig = i - '0';

    i = nocrfgetc();
    *totfig= i - '0';

    i = nocrfgetc();
    *rightfig= i - '0';

   }

   else if ( type == 0 ) *totfig= i - '0'; 

/**********************************************************

 READ IN VALUE, INCLUDING SIGN

***********************************************************/

   *value= 0;
   for ( j = 1; j <= *totfig; j++ ) {

    i = nocrfgetc();

    if ( j > 1 ) *value= 10 * *value + ( i - '0' );
    else {
     sign = (i == '-') ? -1 : 1;
     if ( sign == 1 && i != ' ') *value = (i - '0');
    }

   }

   *value *= sign;

  }

 }

return 0;

}
 
/*************************************************************

                 FUNCTION EXTRACTC

 EXTRACTC EXTRACTS CHARACTER DATA FROM GIVEN FILE

**************************************************************/

extractc(

 int type,              /* SET TO ZERO IF TOTFIG SUPPLIED,
                           SET TO ONE IF TOTFIG READ IN */

 int *totfig,            /* NUMBER OF BYTES IN CHARACTER DATA */

 char *cdata            /* ARRAY OF CHARACTER DATA */

        )

{

 int i,j;

 if ( type !=0 ) i = nocrfgetc();
 
 if ( feof(fp) ) return -1;

 else {

  if ( type !=0 ) *totfig = i - '0';

  for ( j = 0; j < *totfig; j++ ) {

   i = nocrfgetc();
   *(cdata+j) = i;

  }

 }

 if ( feof(fp) ) return -1;
 else return 0;

}

/***********************************************

          FUNCTION NOCRFGETC

 NOCRFGETC (NO CARRIAGE RETURN FGETC) READS OCL FORMAT
 SKIPPING END OF LINE CHARACTERS IN A WAY WHICH WILL 
 WORK IN PC AND UNIX ENVIRONMENT.  THE PC END OF
 LINE CHARACTER (^M), WHICH IS PRESENT IN FILES ON
 THE CDs, IS FOLLOWED BY THE \n CHARACTER ON SOME
 UNIX PLATFORMS

 FP - FILE IDENTIFIER (GLOBAL VARIABLE)

 RETURNS INPUT CHARACTER OR -1 FOR END OF FILE

************************************************/

nocrfgetc()

{

 int i;

 while ( !feof(fp) && isprint ( (i=fgetc(fp)) ) == 0 );

 if ( feof(fp) ) i = -1;

 return i;

} 
