/*  woa01_convert.c
................................................................................
.   Reads World Ocean Atlas gridded standard level hydrographic data files
.   and creates HydroBase *.cdf gridded files.
.   
.   USAGE: woa01_convert -O<outdir> -D<indir> 
...............................................................................
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hydrobase.h"
#include "hydro_cdf.h"


#define NO_PREFILL 0
#define NSD  33      /* No of standard depths in WOA98 */
#define NPROPS 7     /* No of properties not including de */


float stdd[NSD] = { 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250,
                  300, 400, 500, 600, 700, 800, 900, 1000, 1100,
                  1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000,
                  3500, 4000, 4500, 5000, 5500 };
		  

float HB_f_mask = (float)HBMASK;       /* float values are used since the output */
float HB_f_empty = (float)HBEMPTY;     /* cdf values are of this type */

struct CDF_HDR cdf;

 /*  prototypes for locally defined functions */

void print_usage(char *);
FILE *openfile(char *, char *, char *);

int main (int argc, char **argv)
{
   int error, dummy, any_data;
   int i, j, jjj, ip, ii, jj;
   int nsq, indx, msq10;
   int row, col, zlev;
   int cdf_file;
   int include_counts;
   int prop_indx[NPROPS];
   FILE *infile;
   char *indir, *outdir, *monthcode;
   char *sptr, buffer[100], smallbuf[9];
   char filename[200];
   float latmin, lonmin;
   float *xptr, *xout;
   float *bottom;
   float ***p;
   float ***t;
   float ***s;
   float ***o;
   float ***ni;
   float ***si;
   float ***p4;
   float ****data;
   double dlat;
   
 
/*  set these default values */

    error = 0;
    indir = "";
    outdir = "";
    monthcode = NULL;
   
    
/*  parse the command line arguments */

   if (argc > 1) {
      for (i = 1; i < argc; i++) {
         if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':                   /* get input dir */
                        indir = &argv[i][2];
                        break;
               case 'O':                    /* get output dir  */
                        outdir = &argv[i][2];
                        break;
               case 'M':                    /* get month :  00 = clim  01 = jan, 02=feb, etc  */
                        monthcode = &argv[i][2];
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
      }  /* end for */
   } /* end if argc > 1*/

   if (monthcode == NULL) {
        fprintf(stderr,"\nSpecify 2-char month code.  ex:  00 = clim, 01=jan, 02=feb....\n");
        exit(1);
   }   
/* allocate memory for global property arrays */

   fprintf(stderr,"Allocating memory...");

   p = (float ***) calloc((size_t)180, sizeof(float **));
   t = (float ***) calloc((size_t)180, sizeof(float **));
   s = (float ***) calloc((size_t)180, sizeof(float **));
   o = (float ***) calloc((size_t)180, sizeof(float **));
   p4 = (float ***) calloc((size_t)180, sizeof(float **));
   si = (float ***) calloc((size_t)180, sizeof(float **));
   ni = (float ***) calloc((size_t)180, sizeof(float **));

/*   pcount = (int ***) calloc((size_t)180, sizeof(int **));
   tcount = (int ***) calloc((size_t)180, sizeof(int **));
   scount = (int ***) calloc((size_t)180, sizeof(int **));
   ocount = (int ***) calloc((size_t)180, sizeof(int **));
   p4count = (int ***) calloc((size_t)180, sizeof(int **));
   sicount = (int ***) calloc((size_t)180, sizeof(int **));
   nicount = (int ***) calloc((size_t)180, sizeof(int **));  */
   
   for (i = 0; i < 180; ++i)  {
      p[i] = (float **) calloc((size_t)360, sizeof(float *));
      t[i] = (float **) calloc((size_t)360, sizeof(float *));
      s[i] = (float **) calloc((size_t)360, sizeof(float *));
      o[i] = (float **) calloc((size_t)360, sizeof(float *));
      p4[i] = (float **) calloc((size_t)360, sizeof(float *));
      si[i] = (float **) calloc((size_t)360, sizeof(float *));
      ni[i] = (float **) calloc((size_t)360, sizeof(float *));
      
      for (j = 0; j < 360; ++j) {
        p[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        t[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        s[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        o[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        p4[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        si[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
        ni[i][j] = (float *) calloc((size_t)NSD, sizeof(float ));
     
      }
   }
  
/* hardwire the properties that will be stored */

   prop_indx[0] = (int)PR;
   prop_indx[1] = (int)TE;
   prop_indx[2] = (int)SA;
   prop_indx[3] = (int)OX;
   prop_indx[4] = (int)P4;
   prop_indx[5] = (int)SI;
   prop_indx[6] = (int)N3;

   data = (float ****)calloc((size_t)NPROPS, sizeof(float ***));
   data[0] = p;   
   data[1] = t;
   data[2] = s;
   data[3] = o;
   data[4] = p4;
   data[5] = si;
   data[6] = ni;
  
/*  open and read (plow through) temperature analyzed fields */

   sprintf(filename, "t%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) 
      exit(1);
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         dlat = (double) row - 89.5;
         for (col = 0; col < 360; ++col) {
	    xptr = &t[row][col][zlev];
	    if (fscanf(infile, "%f", xptr ) != 1){
	      fprintf(stderr,"\nError reading t00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
	      exit(1);
	    }
	    
	    /* set pressure value if t is defined */
	    
	    if (*xptr < -99.0) {
	       *xptr =  HB_f_mask;
	       p[row][col][zlev] = HB_f_mask;
	    }
	    else
	       p[row][col][zlev] = (float) hb_p80((double)stdd[zlev],dlat);
	 }
      }
   } 
   
   fclose(infile);
   
   
 /*   salinity */
   
   sprintf(filename, "s%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) 
      exit(1);
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         for (col = 0; col < 360; ++col) {
	    xptr = &s[row][col][zlev];
	    if (fscanf(infile, "%f", xptr) != 1){
	      fprintf(stderr,"\nError reading s00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
	      exit(1);
	    }
	    
	    /* if no salinity, mark p,t,s masked */
	    
	    if (*xptr < -99.0) {
	       *xptr = HB_f_mask;
	       t[row][col][zlev] = HB_f_mask;
	       p[row][col][zlev] = HB_f_mask;
	    }
	 }
      }
   } 
   
   fclose(infile);
   
   
 /*   oxygen */
   
   sprintf(filename, "o%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) {
      exit(1);
   }
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         for (col = 0; col < 360; ++col) {
	    xptr = &o[row][col][zlev];
	    if (fscanf(infile, "%f", xptr) != 1){
	      fprintf(stderr,"\nError reading o00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
	      exit(1);
	    }
	    if (*xptr < -99.0)
	       *xptr = HB_f_mask;
	 }
      }
   } 
   
   fclose(infile);
   
   
 /*   nitrate */
   
   sprintf(filename, "n%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) {
      exit(1);
   }
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         col = 0;
         for (jjj = 0; jjj < 36; ++jjj) {
	    if (fscanf(infile, "%[^\n]", buffer) != 1) {
	      fprintf(stderr,"\nError reading n00an1 at row[%d],line[%d],zlev[%d]\n", row, jjj, zlev);
	      exit(1);
	    }
	    fgetc(infile);  /* move past newline */
	    
	    sptr = buffer;
	    for (i = 0; i < 10; ++i) { 
	       strncpy(smallbuf,sptr, 8);
	       smallbuf[8] = 0;
	       xptr = &ni[row][col][zlev];
	       
	       if (sscanf(smallbuf, "%f", xptr) != 1){
	         fprintf(stderr,"\nError parsing n00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
		 fprintf(stderr,"buffer: %s \n smallbuf: %s\n", buffer, smallbuf);
	         exit(1);
	       }
	       if (*xptr < -99.0)
	          *xptr =  HB_f_mask;
		  
	       sptr += 8;
	       ++col;
	    }
	 }
      }
   } 
   
   fclose(infile);
   
   
 /*   silicate */
   
   sprintf(filename, "i%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) {
      exit(1);
   }
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         col = 0;
         for (jjj = 0; jjj < 36; ++jjj) {
	    if (fscanf(infile, "%[^\n]", buffer) != 1) {
	      fprintf(stderr,"\nError reading i00an1 at row[%d],line[%d],zlev[%d]\n", row, jjj, zlev);
	      exit(1);
	    }
	    fgetc(infile);  /* move past newline */
	    
	    sptr = buffer;
	    for (i = 0; i < 10; ++i) { 
	       
	       strncpy(smallbuf,sptr,8);
	       smallbuf[8] = 0;
	       xptr = &si[row][col][zlev];
	       
	       if (sscanf(smallbuf, "%f", xptr) != 1){
	         fprintf(stderr,"\nError parsing i00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
	         exit(1);
	       }
	       if (*xptr < -99.0)
	          *xptr =  HB_f_mask;
		  
	       sptr += 8;
	       ++col;
	    }
	 }
      }
   } 
   
   fclose(infile);
   
   
 /*   phosphate */
   
   sprintf(filename, "p%c%can1", monthcode[0], monthcode[1]);   
   infile = openfile(indir,filename, "");
   if (infile == NULL) {
      exit(1);
   }
  
   for (zlev = 0; zlev < NSD; ++zlev) {
      for (row = 0; row < 180; ++row) {
         col = 0;
         for (jjj = 0; jjj < 36; ++jjj) {
	    if (fscanf(infile, "%[^\n]", buffer) != 1) {
	      fprintf(stderr,"\nError reading p00an1 at row[%d],line[%d],zlev[%d]\n", row, jjj, zlev);
	      exit(1);
	    }
	    fgetc(infile);  /* move past newline */
	    
	    sptr = buffer;
	    for (i = 0; i < 10; ++i) { 
	       
	       xptr = &p4[row][col][zlev];
	       
	       strncpy(smallbuf,sptr,8);
	       smallbuf[8] = 0;
	       if (sscanf(smallbuf, "%f", xptr) != 1){
	         fprintf(stderr,"\nError parsing p00an1 at row[%d],col[%d],zlev[%d]\n", row, col, zlev);
	         exit(1);
	       }
	       if (*xptr < -99.0)
	          *xptr =  HB_f_mask;
		  
	       sptr += 8;
	       ++col;
	    }
	 }
      }
   } 
   
   fclose(infile);
   
/******************* End of Input Phase ************/

/* Construct a cdf header */

  cdf.nx = 10;
  cdf.ny = 10;
  cdf.nz = NSD+1;     /* includes a variable bottom depth */
  cdf.nprops = NPROPS;
  cdf.counts_included = 0;
  cdf.node_offset = 1;
  cdf.fill_value = HB_f_empty;
  cdf.tmin = (int *) calloc(1, sizeof(int)); 
  cdf.tmin[0] = 0;  
  cdf.tmax = (int *) calloc(1, sizeof(int));   
  cdf.tmax[0] = 9999;
  cdf.xincr = 1.0;
  cdf.yincr = 1.0;
  strncpy(cdf.x_units, "degrees", 8);
  strncpy(cdf.y_units, "degrees", 8);
  strncpy(cdf.z_units, "meters", 7);
  strncpy(cdf.title,"World Ocean Atlas", 18);
  strcpy(cdf.command, *argv);
  cdf.prop_id = (char **) malloc(cdf.nprops * sizeof(char *));
  cdf.prop_units = (char **) malloc(cdf.nprops * sizeof(char *));
  for (i = 0; i < cdf.nprops; ++i) {
      cdf.prop_id[i] = (char *) malloc(3);
      cdf.prop_units[i] = (char *) malloc(50);
      strncpy(cdf.prop_id[i], get_prop_mne(prop_indx[i]),3);
      strcpy(cdf.prop_units[i], get_prop_units(prop_indx[i]));
   }
  
  
/* initialize standard depths */

  NSTDLEVS = cdf.nz;
  for (zlev = 0; zlev < NSD; ++zlev)  /* loop for WOA std levs */
     std_depth[zlev] = (double) stdd[zlev];
     
  std_depth[zlev] = -99.0;   /* add the variable bottom depth */  
  std_depth_initialized = 1;

/* Since everything is standard depths, mark the bottom depth empty */

  nsq = cdf.nx * cdf.ny;
  bottom = (float *) calloc((size_t)nsq, sizeof(float));
  
  for (i = 0; i < nsq; ++i)
     bottom[i] = HB_f_empty;
  
/* allocate memory to output a latitude band at a time */

  xout = (float *) calloc((size_t)(cdf.nz*cdf.nx), sizeof(float));

/*  loop for each output cdf file */

   latmin = -90.0;
   for (jj = 0; jj < 18; ++jj) {
   
      lonmin = 0.0;

      for (ii = 0; ii < 36; ++ii) {
      
        msq10 = ms10(latmin, lonmin, &dummy);
	if (*outdir != '\0') 
	   sprintf(filename, "%s/%4d.woa_%2s.nc", outdir, msq10, monthcode);
	else
	   sprintf(filename, "%4d.woa_%2s.nc", msq10, monthcode);
	   
        cdf_file = cdf_init(filename);
	fprintf(stderr,"\nOpened %s for output", filename); 
	   
        cdf.xmin = lonmin;
        cdf.xmax = lonmin + 10.0;
        cdf.ymin = latmin;
        cdf.ymax = latmin + 10.0;
	   
        error = cdf_define(cdf_file, &cdf, NO_PREFILL, cdf.counts_included);
        error = write_std_depths_cdf(cdf_file, &cdf);
        error = write_time_bins_cdf(cdf_file, &cdf);
	error = write_bottom_depth_cdf(cdf_file, 0, 0, 0, cdf.ny, cdf.nx, 1, bottom);
	
	
	/* file is now ready for writing data 
	
	  ip = index to property
	  [col][row] are indices to input (global) arrays   
	  i, j are indices to 10-deg files
	  indx is for the xout array.
	  first point in ncfile is at {lonmin,latmax} 
	  so start at top row and work down.        */
	
	/* do each property  */
	
	any_data = 0;
        for (ip = 0; ip < cdf.nprops; ++ip) {
	   row = jj * 10 + 10; 
	             
	   for (j = 0; j < cdf.ny; ++j) {    
	      --row;
	      col = ii * 10;
	      indx = 0;
	      
	      for (i = 0; i < cdf.nx; ++i) {
	      
	         for (zlev = 0; zlev < NSD; ++zlev) {  /* standard levels */
		    xout[indx] = data[ip][row][col][zlev];
		    if (ip == 2 && !any_data){
		      if (xout[indx] < 40 && xout[indx] > 0)
		        any_data = 1;
		    }
		    ++indx;
		 } 
		 
		 xout[indx++] = HB_f_empty;   /* mark bottom level empty */
	         ++col;
	      }  /* end for i */
	      
	      write_prop_cdf(cdf_file, xout, cdf.prop_id[ip], j, 0, 0, 0, 1, cdf.nx, 1, cdf.nz);
	      
	   }  /* end for j */
	}  /* end for ip */
	
        cdf_close(cdf_file);
	
	lonmin += 10.0;
        if (lonmin >= 179.9990)
           lonmin -= 360;
	
      } /* end for ii */
      
      latmin += 10.0;
   
   } /* end for jj */
    
  fprintf(stderr,"\n\nThat's it.... we're done.\n");
  exit(0);
        
} /* end main() */ 
  
/*****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\nConverts World Ocean Atlas 1-deg analyzed fields");
   fprintf(stderr,"\nto HydroBase3 gridded *.nc files.  Output files are named by");
   fprintf(stderr,"\nto 4-digit WMOsq with extent .woa.cdf");
   fprintf(stderr,"\nUsage:  %s [-Ooutdir] [-D<dirname>] [-E<file_extent>] [-h]", program);
   fprintf(stderr,"\n-D : directory for input files t00an1, s00an1, o00an1, etc... ");
   fprintf(stderr,"\n     (default is current directory) ex: -D../data/ ");
   fprintf(stderr,"\n-O : directory for output files.");
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
