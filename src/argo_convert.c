/*  argo_convert.c
................................................................................
.   Reads ARGO NetCDF profile format files
.   extracts :    header info
.                 p,t90,s
.   
.   USAGE: argo_convert infile_list  [-B<badfile>] [-D<dir>] [-E<extent>] 
                        [-O<outfile>] [-Q<argo_quality_file>] [-V]
...............................................................................
*/

#define     LENNAME     1000
#define     M_LEVELS   5000
#define     MPROPS  4  /* pr de t90 sa */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "netcdf.h"
#include "hydrobase.h"


  /*  HydroBase globally defined variables */
  
struct HYDRO_HDR  hdr;
struct HYDRO_DATA data;

int    fdcdf;
int    lverbose;
int    npropsout;
int    staread, staout, stabad, sta_tonly, sta_err;
int    topt, qopt, sta_noqc, ropt;
FILE *outfile, *tfile, *no_qcfile, *badfile;

  /*  prototypes for locally defined functions */
void   print_usage (char *);
int    open_ncfile (char *, char *, char *);
int    apnc_read (void);
int    get_param_ad (char *, int, int, double *, int *);
int    check_qual (double *, double *, char *, char *, int);
int    check_sta (char *, int);

main (int  argc, char **argv)
{

   int    i,j, error, nfiles, curfile, status, oflag;
   char   *dir, *extent;
   FILE   *fp;

/* check for command line arguments */

   if (argc < 2) {
     print_usage(argv[0]);
     exit(1);
   }

/*  set these default values */

    curfile = 1;
    dir = "";
    error = 0;
    extent = "";
    lverbose = 0;
    nfiles = 0;
    topt = qopt = ropt = 0;
    oflag = 0;
    staout = staread = stabad = sta_tonly = sta_noqc= 0;
    sta_err = 0;

     /* Set the maximum property value, as define for the system 
        in MAXPROP in the include file, and the maximum number for this
        application.
      */
    npropsout = MPROPS;
    for (i = 0; i < MAXPROP; ++i) {
       data.observ[i] = NULL;
       }
 

/*  parse the command line arguments */
   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'B':                    /* open  file for output */
              fprintf (stderr,"\nThis option no longer supported\n");
              fprintf (stderr,"Use -R for rejected profiles.\n");
	      error = 1;
            break;

          case 'D':                   /* get input dir */
            dir = &argv[i][2];
            break;

          case 'E':                    /* get file extent */
            extent = &argv[i][2];
            break;

          case 'O':                    /* get output file  */
	    oflag = 1;
            outfile = create_hydro_file (&argv[i][2], OVERWRITE);
            if (outfile == NULL) {
              fprintf (stderr,"\nError opening output file: %s\n", &argv[i][2]);
              exit(1);
              }
            fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
            break;

          case 'Q':                    /* save profiles with no QC  */
             no_qcfile = create_hydro_file (&argv[i][2], OVERWRITE);
             if (no_qcfile == NULL) {
                fprintf (stderr,"\nError opening output file: %s\n", &argv[i][2]);
                exit(1);
              }
            fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
	    qopt = 1;
            break;
           case 'R':                    /* rejected profiles  */
             badfile = create_hydro_file (&argv[i][2], OVERWRITE);
             if (badfile == NULL) {
                fprintf (stderr,"\nError opening output file: %s\n", &argv[i][2]);
                exit(1);
              }
            fprintf(stderr,"\nOpening or overwriting: %s", &argv[i][2]);
	    ropt = 1;
            break;
         case 'T':                    /* open  file for temp_only output */
            tfile = create_hydro_file(&argv[i][2], OVERWRITE);
            if (tfile == NULL) {
              fprintf (stderr,"\nError opening output file: %s\n", 
                       &argv[i][2]);
              exit(1);
              }
            fprintf(stderr,"\nTemperature-only profiles will be written to: %s", &argv[i][2]);
            topt = 1;
            break;

          case 'V':                    /* set verbose/debugging.  */
            lverbose = 1;
            if (strlen(argv[i]) > 2)
              sscanf (&argv[i][2], "%d", &lverbose);
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

        }  /* end if option flag. */

       else  {
        ++nfiles;
       }
     }  /* end for each argument. */

   
   if (! oflag) {
       fprintf(stderr,"\nYou MUST specify an output file with -O\n");
       exit(1);
   }
  
   if (! topt) {
       fprintf(stderr,"\nTemperature-only profiles will not be saved in a separate file\n");
   }
   if ( !nfiles) {
       fprintf(stderr,"\nYou must specify an ARGO nc file for input.\n");
       exit(1);
   } 
   
   /* set these hard-wired fields */  

    hdr.origin = '3';        
    hdr.instrument = 'f';     

    hdr.prop_id = (int *) calloc(MPROPS, sizeof(int));
    
 /*  loop for each input file */
  do {
  
      error = open_ncfile (dir, argv[curfile], extent);
      if (error != NC_NOERR)  {
        goto NEXTFILE;
      }

     /* loop for each  profile. */
     
       status = apnc_read ();
        	          
    if (status != EOF ) {
       fprintf(stderr, "WARNING:  argo_read() returned value other than EOF\n");
    }
    
    error = nc_close (fdcdf);
       
NEXTFILE:
    ;
  } while   (curfile++ < nfiles);       

fprintf(stderr,"\nEnd of conversion.");
fprintf(stdout,"\n  %d files processed", nfiles);
fprintf(stdout,"\n  %d stations read in", staread);
fprintf(stdout,"\n  %d stations accepted", staout);
fprintf(stdout,"\n  %d stations rejected ", stabad);
fprintf(stdout,"\n  %d stations were not QC'd ", sta_noqc);
fprintf(stdout,"\n  %d stations had only pr,te profiles ", sta_tonly);
fprintf(stdout,"\n  %d stations had time/position errors \n\n", sta_err);
exit(0);

} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"%s converts QC'd observed profile data from ARGO \n", program);
   fprintf(stderr," into HydroBase3 station format. \n");
   fprintf(stderr,"\nUSAGE:  %s infile_list -O<outfile> ", program);
   fprintf(stderr,"\n   [-D<dir>] [-E<extent>] ");
   fprintf(stderr,"\n  [-T<filename>] [-V<detail>] ");
   fprintf(stderr,"\n infile_list: List of filenames must be first arguments.");
   fprintf(stderr,"\n      -O:  specify output file for HydroBase station files ");
   fprintf(stderr,"\n OPTIONS:     ");
   fprintf(stderr,"\n    [-D]:  directory for input files ");
   fprintf(stderr,"\n    [-E]:  extent for input files ");
   fprintf(stderr,"\n               --or stations are output to STDOUT");
   fprintf(stderr,"\n    [-Q]:  file to store profiles with no QC");
   fprintf(stderr,"\n    [-R]:  file to store rejected profiles ");
   fprintf(stderr,"\n    [-T]:  file to store temperature-only profiles");
   fprintf(stderr,"\n    [-V]:  verbose flag [1 2]  default for -V: 1");
   fprintf(stderr,"\n    [-h]:  help");
   fprintf(stderr," \n");
   return;
}
   
/*****************************************************************************/
   
/*****************************************************************************/
int open_ncfile (char *dir, char *root, char *extent)
{
   char st[80];
   int i, iop;
   
   strcpy(st, dir);
   if ((i=strlen(dir)) != 0) {
      if (dir[i-1] != '/')
         strncat(st,"/",1);
   }
   strcat(st, root);
   strcat(st, extent);
   iop = nc_open (st, NC_NOWRITE, &fdcdf);
   if (iop != NC_NOERR)  {
     fprintf(stderr,"\n Unable to open %s \n", st);
     }
    else
     fprintf(stderr,"   opened %s ... \n", st);
   
   return (iop);
   
}  /* end open_ncfile() */
/*****************************************************************************/
int apnc_read ()

   /*  Reads one entire  file. 
       with all its properties and flags into global Hydrobase structures 
      Returns status EOF for successful read, 
     An error reading or parsing causes an error message and an exit from program.
 */
{
   char   theVar[30], ctem[M_LEVELS], ctem_fill[M_LEVELS];
   char   data_state[5];
   int    pr_avail, te_avail, sa_avail, nlevs;
   int    status, staOK;
   int    i, k_prof, out_code;
   size_t    ib[3], ic[3];
   int    item, item_fill;
   int    idate[3], itime[3];
   int    dimid, varid;
   size_t    n_prof, n_levels, n_param;
   double *xprof;
   double dtem, dtem_fill;
   double juldate, juldate_ref;

/* Begin by reading  parameters that apply to entire file */   

     /* Get the number of profiles in this file group. */
     
     strcpy (theVar, "N_PROF\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension id for %s\n.", theVar);
       return (9);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_prof);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (n_prof < 1)  {
       fprintf(stderr,"Error checking apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %d\n",  theVar, (int)n_prof);


     /* Get the number of pressure levels.
        This is the maximum number of the set and used for the array 
        size so each individual station will have to be tested.
      */
     strcpy (theVar, "N_LEVELS\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension id for %s\n.", theVar);
       return (9);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_levels);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (n_levels < 1)  {
       fprintf(stderr,"Error checking apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %d\n",  theVar, (int)n_levels);
     if (n_levels >= M_LEVELS)  { 
       fprintf (stderr, " FATAL ERROR: Maximum levels in array exceeded \n");
       fprintf (stderr, "      Set M_LEVELS to larger value in program.\n");
       exit (1);
       }


     /* Get the number of parameters in this group. 
        This is again a maximum dimension number.
      */
     strcpy (theVar, "N_PARAM\0");
     status = nc_inq_dimid (fdcdf, theVar, &dimid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension id for %s\n.", theVar);
       return (9);
       } 
     status = nc_inq_dimlen (fdcdf, dimid, &n_param);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (n_param <= 1)  {
       fprintf(stderr,"Error checking apnc dimension value for %s\n.", theVar);
       return (9);
       } 
     if (lverbose > 0)  fprintf (stderr, "  %s  =  %d\n",  theVar, (int)n_param);


     /* We will need the reference date for later.
        As originally explained the files are grouped by day but
        they could be ordered any way in general so don't rely on the
        file name for the date.
      */
     strcpy (theVar, "REFERENCE_DATE_TIME\0");
     status = nc_inq_varid (fdcdf, theVar, &varid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       return (9);
       } 
     status = nc_get_att_text (fdcdf, varid, "_FillValue", ctem_fill);
     ctem_fill[1] = '\0';
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc fill attribute for %s.\n", theVar);
       return (9);
       } 
     status = nc_get_var_text (fdcdf, varid, ctem);
     ctem[14] = '\0';
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       return (9);
       } 
     if (strcmp (ctem, ctem_fill) == 0)  {
       fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       fprintf(stderr,"[%s]\n", ctem);
       } 

     sscanf (ctem, "%4d%2d%2d%2d%2d%2d",
             &idate[0], &idate[1], &idate[2],
             &itime[0], &itime[1], &itime[2]);
     gdate2dday (idate, itime, &juldate_ref);
     if (lverbose > 0)  
       fprintf (stderr, "  %s  = %lf    %s  \n",  theVar, juldate_ref, ctem);
    
 /********** End file initialization. ***********/

  /* Loop for each profile: check for quality and adjusted values. */

   for (k_prof=0; k_prof < n_prof; k_prof++)  {
   
     if (lverbose > 0)  fprintf (stderr, "\n\n  Start of Station.\n");

     ++staread;
     
     /* Get the locations. */
     ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
     ic[0] = 1;      ic[1] = 0; ic[2] = 0;
     
     strcpy (theVar, "LATITUDE\0");
     status = nc_inq_varid (fdcdf, theVar, &varid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     if (dtem > 90 || dtem < -90)  {
       fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     hdr.lat = (float) dtem;
     if (lverbose > 0)
       fprintf (stderr, "  %s: %lf   %lf\n", theVar, dtem, dtem_fill);

     strcpy (theVar, "LONGITUDE\0");
     status = nc_inq_varid (fdcdf, theVar, &varid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     hdr.lon = (float) dtem;
     if (lverbose > 0)  
       fprintf (stderr, "  %s: %lf   %lf\n", theVar, dtem, dtem_fill);

     /* get DATE */
     strcpy (theVar, "JULD_LOCATION\0");
     status = nc_inq_varid (fdcdf, theVar, &varid);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc fill attribute for %s.\n", theVar);
       ++sta_err;
       continue;
       } 
     status = nc_get_var1_double (fdcdf, varid, ib, &dtem);
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
     if (dtem == dtem_fill || dtem == 0.0)  {
       fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       fprintf(stderr,"[%lf]\n", dtem);
       /* try to use date in last header, or make a false date */
       if (hdr.year < 1990 || hdr.year > 2020) {  
           hdr.year = 9999;
	   hdr.month = 99;
	   hdr.day = 99;
       }
      } 
     else {
         item = dtem + juldate_ref;
         jdate (item, &hdr.day, &hdr.month, &hdr.year);
     }
   if (lverbose > 0)
     fprintf (stderr, "  %s: %d   %d  %d\n", theVar, hdr.year, hdr.month, hdr.day);


     /* There really is no country in the ARGO profile format so 
        we leave that the default "XX\0".
        The ARGO PLATFORM_NUMBER is commonly the 7 character WMO identifier 
        but some of the earlier profiles seem to have the 5 character 
        PTT identifier. Some PTT's were reused later but that is
        probably beyond the scope of this application. As the WMO identifer
        was specially composed of the 2 digit ocean area code and 5 digit
        sequence number we will try to split it up the best we can,
        hoping there is only numerics in the platform_number (as the
        data management team assured us would always be the case). This
        will map into the HB 2 character ship and HB interger cruise. The
        CYCLE_NUMBER will map into the HB integer station.
        The PLATFORM_NUMBER is defined to be eight characters so 
        force the [8] = '\0'.
      */
      
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 8; ic[2] = 0;
   
   strcpy (theVar, "PLATFORM_NUMBER\0");
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   status = nc_get_att_text (fdcdf, varid, "_FillValue", ctem_fill);
   ctem_fill[1] = '\0';
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc fill attribute for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   strcpy (ctem, "        \0");
   status = nc_get_vara_text (fdcdf, varid, ib, ic, ctem);
   ctem[8] = '\0';
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   if (strcmp (ctem, ctem_fill) == 0)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
          /* Round-about check for the length  of the string. */
   item = strspn (ctem, "0123456789");
   if (lverbose > 0)  fprintf (stderr, "  PLATFORM: %s   %d\n", ctem, item);
   ctem[item] = '\0';
   strcpy (hdr.country, "XX\0");
   strcpy (hdr.ship, "XX\0");
   if (item <= 5)  {
     ctem[7] = '0';
     }
    else if (item == 6)  {
     hdr.ship[0] = ' ';
     hdr.ship[2] = ctem[0];
     strcpy (ctem, ctem+1);
     }
    else  {
     strncpy (hdr.ship, ctem, 2);
     strcpy (ctem, ctem+2);
     }
   status = sscanf (ctem, "%d", &hdr.cruise);
   if (hdr.cruise <= 0) {
     hdr.cruise = 99;
     }
   if (lverbose > 0)
     fprintf (stderr, "  SHIP/CRUISE: %s   %d\n", hdr.ship, hdr.cruise);

     /* Map CYCLE_NUMBER to station. */
   strcpy (theVar, "CYCLE_NUMBER\0");
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 0; ic[2] = 0;
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   status = nc_get_att_int (fdcdf, varid, "_FillValue", &item_fill);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc fill attribute for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   status = nc_get_var1_int (fdcdf, varid, ib, &item);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   if (item == item_fill)  {
     fprintf(stderr,"Invalid/missing variable value for %s.\n", theVar);
       ++sta_err;
       continue;
     } 
   hdr.station = item;
   if (hdr.station > 9999) {
       ++sta_err;
       continue;
     }
   if (lverbose > 0) 
     fprintf (stderr, "  %s: %d   %d\n", theVar, item, item_fill);
     
     
     /* Data State Indicator (Argo Ref Table #6. */
   strcpy (theVar, "DATA_STATE_INDICATOR\0");
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 4; ic[2] = 0;
   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"WARNING: Error getting apnc variable id for %s.\n", theVar);
   } 
   status = nc_get_vara_text (fdcdf, varid, ib, ic, data_state);
   if (status != NC_NOERR) {
       fprintf(stderr,"WARNING: Error getting apnc variable value for %s.\n", theVar);
       strncpy(data_state,"    ",4);
   } 
   

     /* Get the  profile data. */
    
    hdr.nprops = MPROPS;
    hdr.prop_id[0] = (int)PR; 
    hdr.prop_id[1] = (int)DE;
    hdr.prop_id[2] = (int)T90;
    hdr.prop_id[3] = (int)SA;
     

   /* ib = start vector, ic = count vector */
   
   ib[0] = k_prof; ib[1] = 0; ib[2] = 0;
   ic[0] = 1;      ic[1] = 1; ic[2] = 4;
   pr_avail = te_avail = sa_avail = 0;
   staOK = 1;
        
   status = nc_inq_varid (fdcdf,"STATION_PARAMETERS" , &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
       ++sta_err;
       continue;
   } 

   while (ib[1] < n_param)  {
     status = nc_get_vara_text (fdcdf, varid, ib, ic, ctem);
     ctem[4] = '\0';
     if (status != NC_NOERR) {
       fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
       ++sta_err;
       continue;
       } 
     if (strncmp (ctem, "PRES", 4) == 0)  {
       pr_avail = 1;
     }
     if (strncmp (ctem, "TEMP", 4) == 0) {
        te_avail = 1;
     }
     if (strncmp (ctem, "PSAL", 4) == 0) {
        sa_avail = 1;

     }
     ib[1]++;
   } /* end while */
      
 /* PRESSURE */
   if (pr_avail)  {
     xprof = (double *) calloc(M_LEVELS, sizeof(double));
     status = get_param_ad ("PRES", k_prof, n_levels, xprof, &nlevs);
     if (status != NC_NOERR)  exit (1);

     if (nlevs > 0)
        data.observ[(int)PR] = xprof;  
     else {
        staOK = 0;
	free(xprof);
     }   

   }  /* End if pressure. */

   if (te_avail && staOK)  {
     xprof = (double *) calloc(M_LEVELS, sizeof(double));
     status = get_param_ad ("TEMP", k_prof, n_levels, xprof, &nlevs);
     if (status != NC_NOERR)  exit (1);
     
     if (nlevs > 0)
        data.observ[(int)T90] = xprof;  
     else {
        staOK = 0;
	free(xprof);
	free(data.observ[(int)PR]);
	data.observ[(int)PR] = NULL;
     }   
     
   }  /* End if temperature. */
   

   if ( staOK)  {
     xprof = (double *) calloc(M_LEVELS, sizeof(double));
     if (sa_avail) {
        status = get_param_ad ("PSAL", k_prof, n_levels, xprof, &nlevs);
        if (status != NC_NOERR)  exit (1);
     }
     else {
       for (i = 0; i < n_levels; ++i) 
          xprof[i] = HB_MISSING;
     }
     data.observ[(int)SA] = xprof;  

   } 

   if (!staOK) {
       ++stabad;
       if (ropt)
	  write_hydro_station(badfile, &hdr, &data);
    }
   
   if (staOK) {
   
      data.observ[(int)DE] = (double *) calloc(M_LEVELS, sizeof(double));
       
      out_code = check_sta(data_state, n_levels);

      switch (out_code) {
         case 0:
             if (ropt)
	          write_hydro_station(badfile, &hdr, &data);
	     ++stabad;
         case 1:
	    if (topt) 
	       write_hydro_station(tfile, &hdr, &data);
	    ++sta_tonly;
	    break;
         case 2:
	    write_hydro_station(outfile, &hdr, &data);
	    ++staout;
	    break;
         case 3:
	    if (qopt)
	        write_hydro_station(no_qcfile, &hdr, &data);
	    ++sta_noqc;
	    break;
      }
      
      free(data.observ[(int)PR] );
      free(data.observ[(int)DE] );
      free(data.observ[(int)T90] );
      free(data.observ[(int)SA] );
      data.observ[(int)PR] = NULL;
      data.observ[(int)DE] = NULL;
      data.observ[(int)T90] = NULL;
      data.observ[(int)SA] = NULL;
   } /* end if sta_OK */
  

 } /* end for k_prof */
 
   return(EOF);
}  /* end apnc_read() */
/*****************************************************************************/
int  get_param_ad (char *theVar, int k_prof, int count, double *theArray, int *ngood_addr)

/* Read an entire profile for specified parameter, its associated adjusted array and quality flags.  
   Returns best values in theArray, with HB_MISSING at levels of unacceptable quality.  Counts the
   number of good values and returns as argument.
   Returns NC_NOERR for successful read, or error code and error message if not. */
{
int     i;
size_t ib[3], ic[3];
int     status;
int     varid;
char    *qc, *qc_adj;
char    *varname;
double  *xadj;
double  dtem_fill;

   if (lverbose > 0)  
     fprintf (stderr, "  LOAD %s: %d  %d\n", theVar, k_prof, count);

   status = nc_inq_varid (fdcdf, theVar, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable id for %s.\n", theVar);
     return (status);
   } 
   status = nc_get_att_double (fdcdf, varid, "_FillValue", &dtem_fill);
   if (status != NC_NOERR) {
      fprintf(stderr,"Error getting apnc fill attribute for %s.\n", theVar);
     return (status);
   } 
   dtem_fill = ABS(dtem_fill) - 1 ;

   ib[0] = k_prof;   ib[1] = 0;        ib[2] = 0;
   ic[0] = 1;        ic[1] = count;    ic[2] = 0;

   status = nc_get_vara_double (fdcdf, varid, ib, ic, theArray);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable value for %s.\n", theVar);
     fprintf (stderr, "  ERROR STATUS: %d\n", status);
     return (status);
   } 
 
  /* Replace filled values with HB_MISSING flag */ 
  for (i = 0; i < count; i++)  {
     if (theArray[i] > dtem_fill)  
       theArray[i] = HB_MISSING;
     
  }  /* End for */

/* get qc flags */  
  varname = (char *) calloc(50, sizeof(char));
  strcpy(varname, theVar);
  strncat(varname,"_QC", 4);
   status = nc_inq_varid (fdcdf, varname, &varid);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting apnc variable id for %s.\n", varname);
     return (status);
   } 
   qc = (char *) calloc(count, sizeof(char));
   status = nc_get_vara_text (fdcdf, varid, ib, ic, qc);
   if (status != NC_NOERR) {
     fprintf(stderr,"Error getting array variable  for %s.\n", varname);
     fprintf (stderr, "  ERROR STATUS: %d\n", status);
     return (status);
   } 
/* see if adjusted array exists */  
   xadj = (double *) NULL;
   qc_adj = (char *) NULL;
   varname = (char *) calloc(50, sizeof(char));
   strcpy(varname, theVar);
   strncat(varname,"_ADJUSTED", 9);
   status = nc_inq_varid (fdcdf, varname, &varid);
   if (status == NC_NOERR) {
      xadj = (double *) calloc(count, sizeof(double));
      qc_adj = (char *) calloc(count, sizeof(char));
      status = nc_get_vara_double (fdcdf, varid, ib, ic, xadj);
      strncat(varname,"_QC", 4);
      status = nc_inq_varid (fdcdf, varname, &varid);
      if (status != NC_NOERR) {
         free(xadj);
         free(qc_adj);
	 xadj = NULL;
	 qc_adj = NULL;
      }
      else {
         status = nc_get_vara_text (fdcdf, varid, ib, ic, qc_adj);
        /* Replace filled values with HB_MISSING flag */ 
         for (i=0; i < count; i++)  {
           if (xadj[i] > dtem_fill)  
              xadj[i] = HB_MISSING;
         }  /* End for */
	 
      }
   }
   
   *ngood_addr = check_qual(theArray, xadj, qc, qc_adj, count);
   if (xadj != NULL)
       free(xadj);
   if (qc_adj != NULL)
       free(qc_adj);
           
   free(qc);
   free(varname);
   return (NC_NOERR);
}  /* End of get_param_ad. */

/*****************************************************************************/
int check_qual (double *prop,  double *prop_adj, char *qc, char *qc_adj, int count)
   /* Returns best choice of property values (either good values
   or missing) based on quality flags and existence of adjusted data.
    */
{
int i, ngood;

/* Replace original prop value with  adjusted data if flagged as good */
   for (i = 0; i < count; ++i) {
       if (qc_adj[i] == '1') {
	   prop[i] = prop_adj[i];
	   qc[i] = qc_adj[i];
       } 
   }
   
/* Now replace levels flagged as bad */

   ngood = 0;
   for (i = 0; i < count; ++i) {
       switch (qc[i]) {
          case '3':
	  case '4':  /* fall through */
	  case '5':
	  case '8':
	  case '9':
	     prop[i] = HB_MISSING;
	     break;
	  case ' ':
	     break;
	  default:  
             ++ngood;
       } /* end switch */
   }
   return (ngood);
} /* end check_qual */
/*****************************************************************************/
int check_sta (char *dsi, int nlevs)

 /* dsi = 4-char DATA_STATE_INDICATOR from Argo ref table #6
    nlevs = max number of obs, needs to be adjusted for 
    
    Sets global variables hdr, data for output.  
     Returns  0 if profile entirely bad
              1 for pr,te only
	      2 for pr, te, sa
	      3 for no quality control done
    */
{
   int    i, j, code, npropsout;
   

   j = 0;
   for (i = 0; i < nlevs; ++i) {
   
      if ( (data.observ[(int)PR][i] >= 0.0) &&(data.observ[(int)T90][i] >= -3.0) && (data.observ[(int)SA][i] >= 0.0)) {
         data.observ[(int)PR][j] = data.observ[(int)PR][i];
         data.observ[(int)T90][j] = data.observ[(int)T90][i];
         data.observ[(int)SA][j] = data.observ[(int)SA][i];
         data.observ[(int)DE][j] = hb_depth(data.observ[(int)PR][i], (double)hdr.lat);
	 ++j;
      }
   } /*end for*/
   
   if (j > 0) {
       hdr.nobs = j;
       code = 2;
       npropsout = 4;
   } 
   else {
      j = 0;
      for (i = 0; i < nlevs; ++i) {
   
         if ( (data.observ[(int)PR][i] >= 0.0) &&(data.observ[(int)T90][i] >= -3.0) ) {
            data.observ[(int)PR][j] = data.observ[(int)PR][i];
            data.observ[(int)T90][j] = data.observ[(int)T90][i];
            data.observ[(int)SA][j] = data.observ[(int)SA][i];
            data.observ[(int)DE][j] = hb_depth(data.observ[(int)PR][i], (double)hdr.lat);
	    ++j;
         }
      } /* end for i */
      
      if (j > 0) {
           code = 1;
	   hdr.nobs = j;
	   npropsout = 3;
      }
      else {
	   return(0);
      }
   }  
             
   if(dsi[1]  == 'A' || dsi[0] == '0') {  /* no QC done */
       if (code == 2) 
          code = 3;
   }    
       
   hdr.nprops = npropsout;
   data.nprops = npropsout;
   data.nobs = hdr.nobs;
 
   hdr.qual[0] = '0';
   hdr.qual[1] = '0';
   hdr.qual[2] = '0';
   hdr.qual[3] = '1';
   if (code == 3)
       hdr.qual[3] = '0';
   
	  
   hdr.ms10 = ms10(hdr.lat, hdr.lon, &hdr.ms1);
                   
   return (code);

} /* end check_sta */


