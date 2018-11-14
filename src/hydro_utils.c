/*  hydro_utils.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated to ANSI C Nov 1999
			     Updated for HB3 Dec 2009
			     Updated with T90 variables
................................................................................
.
.  A set of routines to facilitate reading, writing files in HydroBase
.  format.
................................................................................
................................................................................
 FILE *open_hydro_file(char *dir, char *root, char *extent, int print_msg):
       Opens an existing hydro file for reading.
 int  get_station(FILE *fptr, struct HYDRO_HDR *h_addr, struct HYDRO_DATA *d_ptr):
       Reads a complete dptr->
 int  read_hydro_hdr(FILE *fptr, struct HYDRO_HDR *haddr):    
       Gets station header info.
 void report_status(int status, FILE *file)
       Reports any errors reported by read_hydro_hdr().
 int  get_data_scan(FILE *fptr, double *scan, int n, int *prop_order) :  
       Reads the next scan in an open hydro file
 int  get_separator(FILE *fptr) :  
       Advances file beyond the station separator.
 FILE *create_hydro_file(char *fname, int mode):  
       Opens a HydroBase file for writing with specified mode.
 int  write_hydro_station(FILE *fptr, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr):
       Writes an entire station to an open hydro file .  
 int  write_hydro_hdr(FILE *fptr, struct HYDRO_HDR *hptr) : 
       Writes a header to an open hydro file 
 int  write_hydro_scan(FILE *fptr, double *scan, int n, int *prop_order) : 
       Writes an individual data scan 
 int  write_separator(FILE *fptr) : 
       Writes a separator line.
 int  ms10(float lat, float lon, int *ms1_ptr) : 
       Computes the 10-degree & 1-degree Marsden Square # for a lat,lon pair.
 int  available(enum property p, struct HYDRO_HDR *hptr) :
       Returns a 1 if the specified property is available for a dptr-> 
 int  is_in_range(float lon, float xmin, float xmax, int *merid, int lon0to360)
       Returns a 1 if specified lon is in the range xmin -> xmax 
 void free_hydro_data(struct HYDRO_DATA *staptr)
       Frees all arrays and fields in the structure.
 void free_and_alloc(double **dptr, int n)
       Frees up space pointed to by dptr and allocates space for n elements.
 void list_origin(FILE *fptr);
       Lists origination code table
 void list_instrument(FILE *fptr);
       Lists instrument code table
 void create_t68_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates t68 data (observ,count,variance, quality) from t90 arrays
 void create_t90_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates t90 data (observ,count,variance, quality) from t68 arrays
 void create_pr_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates pr arrays from de arrays
 void create_de_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates de arrays from pr arrays
 void create_o2_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates o2 arrays from ox arrays
 void create_ox_arrays(struct HYDRO_HDR *, struct HYDRO_DATA *);
       Creates ox arrays from o2 arrays
................................................................................
................................................................................
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"

#define NORIG 11
char orig_code[NORIG] = {   '0',  /* miscellaneous */
                            '1',  /* NODC / WOD01 */
                            '2',  /* ICES */
                            '3',  /* WOD05 */
                            '4',  /* WHPO */
                            '5',  /* BBSR */
                            '6',  /* WHOI */
                            '7',  /* SIO */
                            '8',  /* WGHC */
                            '9',   /* individual PI */
                            ' '   /* unknown */
 };
                         
char *orig_name[NORIG] = {  "miscellaneous",
                            "NODC/WOD2001" ,
                            "ICES",
                            "WOD2005",
                            "WHPO",
                            "BIOS",
                            "WHOI",
                            "SIO",
                            "WOCE Global Hydrographic Climatology (Gouretski & Koltermann, 2004)",
                            "individual PI",
                            "unknown"
}; 

#define NINSTR 7
char instr_code[NINSTR] = {'b',  /* bottle */
                           'c',  /* ctd */
                           'f',  /* profiling float */
                           's',  /* seasoar */
                           'm',  /* moored profiler */
                           'u',  /* unknown */
                           ' '   /* unspecified */
};
                   
char *instr_name[NINSTR] = {"bottle",
                            "ctd",
                            "float",
                            "seasoar",
                            "moored_profiler",
                            "unknown",
                            "unspecified"
};

#define NQUAL_FLAGS 7

int quality_code_flag[NQUAL_FLAGS] = {0, 1, 2, 3, 4, 6, 9};

char * quality_code_name[NQUAL_FLAGS] = 
          {"no quality assessment",
	   "good measurement",
	   "probably good measurement",
	   "probably bad data point",
	   "bad data point",
	   "interpolated data point",
	   "missing data"
};

      
/**********************************************************************/

FILE * open_hydro_file(char *dir, char *root, char *extent, int print_mess)

  /* opens an existing hydro file for READING only. Returns a file
   descriptor if successful or -1 if an error occurs. 
    
     char *dir       name of directory   or "" 
     char *root      root of filename    or full path to file 
     char *extent    extent for file     or "" 
     int  print_mess 0 to suppress printing messages on stderr device 
  */
{ 
   char fname[5000];
   int   i;
   FILE *fptr;

    strcpy(fname, dir);
     if ((i = strlen(dir)) != 0) {
        if (dir[i-1] != '/')
           strncat(fname,"/",1);
     }
     strcat(fname, root);
     strcat(fname, extent);

     fptr = fopen(fname, "r");

     if (print_mess) {
        if (fptr == NULL) {
          fprintf (stderr,"\nunable to open %s for input\n\n", fname);
        }
        else {
          fprintf(stderr,"\nOpened %s ...\n", fname);
        }
     }

     return (fptr);
   
}  /* end open_hydro_file() */

/**********************************************************************/

int get_station(FILE *fptr, struct HYDRO_HDR *h_addr, struct HYDRO_DATA *d_ptr)

  /* retrieves an entire station and places the data in struct HYDRO_DATA. 
     It advances beyond the station separator at the end of the current dptr-> 
     Returns 0 for a successful read, -1 for end_of_file,  
     and an error code if an error occurs.

      error code =  1 :  error parsing header.
                    2 :  error parsing datascan.
                    3 :  error reading station separator.
                    4 :  unexpected end of file occurred.
      arguments:              
      fptr:    already opened file handle 
      h_addr:  station header info (already allocated) 
      d_ptr:   pointer to station data (already allocated)
  */
{
   int i, j, iprop;
   char line[NBSEP+1];
   double *scan;
   double **propaddr, *prop;
 
   if ((i = read_hydro_hdr(fptr, h_addr)) < 0) 
      return (-1);
   else if (i > 0)
      return (1);
   else
      ;

/* assign values at d_ptr and allocate memory for data scans */

   d_ptr->nobs = h_addr->nobs;
   d_ptr->nprops = h_addr->nprops;
   for (i = 0; i < d_ptr->nprops; ++i) {
      iprop = h_addr->prop_id[i];
      propaddr = get_prop_array(d_ptr, iprop); 
      if (*propaddr != NULL) {
           free((void *) *propaddr);
           *propaddr = NULL;
      }
      *propaddr = (double *) calloc((size_t)h_addr->nobs, sizeof(double));
      if (*propaddr == NULL) {
         fprintf(stderr,"\nUnable to allocate memory in get_station()\n");
         exit(1);
      }
   }

    
   scan = (double *) calloc((size_t)h_addr->nprops, sizeof(double));

/* read each scan from file */
   
   for (i = 0; i < h_addr->nobs; ++i) {
      if (get_data_scan(fptr, scan, h_addr->nprops, h_addr->prop_id) > 0) {
           free((void *)scan);
           return (2);
      }
      for (j = 0; j < d_ptr->nprops; ++j) {
         iprop = h_addr->prop_id[j]; 
         propaddr = get_prop_array(d_ptr, iprop);
	 prop = *propaddr;
         prop[i] = scan[j];
      }
   }  /* end for */

   free((void *)scan);

   if (fscanf(fptr,"%[^\n]", line) != 1) {   /* move past station separator */
       return (3);
   }
   getc(fptr);  /* move past linefeed */

   return (0);
}  /* end get_station() */

/**********************************************************************/

int read_hydro_hdr(FILE *fptr, struct HYDRO_HDR *haddr)

  /* reads a header from an already opened hydro file, parses the
    information and stores it at haddr.  Returns 0 for a successful
    read, -1 for end-of-file, or a code > 0 to indicate an error.

          error codes returned :
                1 : error parsing header.
                2 : error parsing property codes
  */
{
    char *s, line[1000];
    int i, nbytes;

    if ( fscanf(fptr, "%[^\n]", line) != 1)
         return (-1);
	 
    getc(fptr);  /* move past linefeed */
	 
    if (sscanf(&line[0],"%s", &haddr->country[0]) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[3],"%s", &haddr->ship[0]) != 1){
         fprintf(stderr,"%s",line);
          return (1);
     }
    if (sscanf(&line[6],"%d", &haddr->cruise) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[12],"%d", &haddr->station) != 1){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (sscanf(&line[17],"%d %d %d",&haddr->year,&haddr->month,&haddr->day) != 3){
         fprintf(stderr,"%s",line);
         return (1);
     }
    if (haddr->year < 1000)
       haddr->year += 1900;
       
    sscanf(&line[68],"%d", &haddr->ms10);  /* offset of ms10 is diagnostic of hdr */
    
    if (haddr->ms10 >= 1000) {               /*  hydrobase3 hdr */
    
      if (sscanf(&line[28],"%c%c", &haddr->origin, &haddr->instrument) != 2)
         return (1);
      if (sscanf(&line[31],"%f %f %d", &haddr->lat, &haddr->lon, &haddr->pdr) != 3)
         return (1);
      if (sscanf(&line[59],"%d %d %d %d", &haddr->nobs, &haddr->nprops, &haddr->ms10, &haddr->ms1) != 4)
         return (1);
         
      for (i = 0; i < 4; ++i)   {    /* load quality codes this way in case of blanks */
         haddr->qual[i] = line[54+i];  
         if (haddr->qual[i] == ' ')
           haddr->qual[i] = '0';
      } 
     
    }
    else {                                  /* old hydrobase hdr */
    
       haddr->origin = '0';
       haddr->instrument = 'u';
       if (sscanf(&line[30],"%f %f %d %d %d %d %d", &haddr->lat, &haddr->lon, &haddr->pdr, &haddr->nobs, &haddr->nprops, &haddr->ms10, &haddr->ms1) != 7)
         return (1);
     }     
   
    
    if ( fscanf(fptr, "%[^\n]", line) != 1)
         return (2);
    getc(fptr);  /* move past linefeed */

/* allocate space to cross-reference properties for this station */

    if (haddr->prop_id != NULL)
         free((void *)haddr->prop_id);

    haddr->prop_id = (int *) calloc((size_t)haddr->nprops, sizeof(int));

/* now get the properties */
    s = line;
    for (i = 0; i < haddr->nprops; ++i) {
        if ((haddr->prop_id[i] = get_prop_indx(s)) < 0) {
           fprintf(stderr,"\n Unknown property listed in file: %.2s\n", s);
           return (2);
        }
	while (*s != ' ')  /* move past this property id */
	    ++s; 
	if (i < (haddr->nprops -1))  /* advance to next property */
	    ++s;
    }
    return(0);

} /* end read_hydro_hdr() */

/**********************************************************************/

void report_status(int status, FILE *file)
  /* Translates and reports errors during read operations in get_station() 
  */
{
  switch (status) {
       case -1:
            break;
       case 0:
            break;
       case 1:
            fprintf(file,"\nError parsing header in get_station(). \n");
            break;
       case 2:
            fprintf(file,"\nError parsing datascan in get_station(). \n");
            break;
       case 3:
            fprintf(file,"\nError reading station separator in get_station(). \n");
            break;
       case 4:
            fprintf(file,"\nUnexpected eof encountered in get_station(). \n");
            break;
       default:
            fprintf(file,"\n Unknown error code returned by get_station() : %d\n", status);

  } /* end switch */
  return;
}  /* end report_status() */

/**********************************************************************/

int get_data_scan(FILE *fptr, double *scan, int n, int *prop_order)

  /* Reads and parses one data scan.  Returns 0 for a successful
     read  or an error code if an error occurs.

      error code =  2 :  error parsing datascan.
                    4 :  unexpected end of file occurred.
                    
     arguments:
       fptr:         already opened file 
       scan:         array to store data values 
       n:            number of properties per scan 
       prop_order:   index to enumerated property list 
  */
{
   char *s, *line;
   int  i;

/* determine number of bytes in a scan allowing for the LF at the end... */

   line = (char *) calloc(5000, sizeof(char));
   if (fscanf(fptr, "%[^\n]", line) != 1 ) {
          free((void *)line);
          return(4);  
   }

   getc(fptr);  /* move past linefeed */
   
/* parse the line into property fields... */
   s = line;
   for (i = 0; i < n; ++i) {             
      if (sscanf(s,"%lf", &scan[i] ) != 1) {
          free((void *)line);
          return(2);
      }
      s += get_field_width(prop_order[i]) + 1;
   }

   free((void *)line);
   return (0);

} /* get_data_scan() */
/**********************************************************************/

int get_separator(FILE *fptr)

   /* Advances a file beyond the station separator.  Returns 0 for
      a successful read or 1 if an error occurs.  */
{
   char line[NBSEP];

   if (fscanf(fptr, "%[^\n]", line) != 1 ) {
          return(1);  
   }
   getc(fptr);  /* move past linefeed */

   return (0);

}  /* end get_separator() */

/**********************************************************************/
FILE * create_hydro_file(char *fname, int mode)

  /* opens a file for writing in specified mode. An error will occur
    if noclobber is specified and the file already exists. Returns a file  
    pointer for a successful open,  or NULL for an error.
    
    arguments:
       fname: string containing name of file to be opened 
       mode:  0=overwrite if file exists, 1=noclobber, 2=append 
   */
{
  FILE *fptr;
  int exists;

   switch (mode) {
       case APPEND:
               fptr = fopen(fname, "a");
                  return(fptr);
               break;
               
       case OVERWRITE:
               fptr = fopen(fname, "w");
                  return(fptr);
               break;

       case NOCLOBBER:  /* fall through */
       default:
           fptr = fopen(fname,"r");
          
	  if (fptr  != NULL) {  /* check if file exists already*/
             return (NULL);     /* return an error if it does */
          }
	  
          fptr = fopen(fname, "w");
          return(fptr);
          break;


   } /* end switch */
   
}  /* end create_hydro_file() */
/**********************************************************************/
int write_hydro_station(FILE *fptr, struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)

   /* writes a complete station to an open file.  Returns 0 for
      a successful write or a number > 0 if an error occurs.
   */
{
   int i, j;
   double *scan, *prop, **propaddr;

   if (write_hydro_hdr(fptr, hptr) > 0)
      return (1);

   scan = (double *) calloc((size_t) hptr->nprops, sizeof(double));

   for (i = 0; i < hptr->nobs; ++i) {
      for (j = 0; j < hptr->nprops; ++j) {
         propaddr = get_prop_array(dptr,hptr->prop_id[j]);
	 
	 if (*propaddr == NULL) {   /* error */
	     free((void *)scan);
	     return(1);
	 }
	 prop = *propaddr;
         scan[j] = prop[i];
      }
      if (write_hydro_scan(fptr, scan, hptr->nprops, hptr->prop_id) > 0) {
            free((void *)scan);
            return(2);
      }
   }
   free((void *)scan);

   write_separator(fptr);

   return(0);

} /* write_hydro_station() */
 
/**********************************************************************/

int write_hydro_hdr(FILE *fptr, struct HYDRO_HDR *hptr)

   /* writes a station header. Returns 0 for a successful operation or
      number > 0 if an error occurs. 
               error code : 1 :  error constructing header line.
               error code : 2 :  error writing line to file. */
{

   char *line, *s;
   int i;

/* check that year is upgraded to 4 digits, not 2 */

   if (hptr->year < 1000)
     hptr->year += 1900;
      
/* make sure quality bytes are not blank */

   for (i = 0; i < NQUAL; ++i) {
      if (hptr->qual[i] < '0' || hptr->qual[i] > '9')
         hptr->qual[i] = '0';
   }
      
/* make sure origin and instrument are not blank */

   if (hptr->origin == '\0') 
      hptr->origin = ' ';
      
   if (hptr->instrument == '\0') 
      hptr->instrument = ' ';

/* compose first header line and write to file ...*/

    line = (char *)calloc(1000, sizeof(char));
    sprintf(line,"%.2s %.2s %5d %4d %4d %2d %2d %c%c %7.3f %8.3f %5d %c%c%c%c %4d %3d %4d %2d", hptr->country, hptr->ship, hptr->cruise, hptr->station, hptr->year, hptr->month, hptr->day, hptr->origin, hptr->instrument, hptr->lat, hptr->lon, hptr->pdr, hptr->qual[0], hptr->qual[1], hptr->qual[2], hptr->qual[3], hptr->nobs, hptr->nprops, hptr->ms10, hptr->ms1); 

   if ( line[NBHEAD-1] != '\0') 
     return (1);
     
   fprintf(fptr,"%s", line); 

   fprintf(fptr,"\n");
   free(line);
   
/*  write line of property descriptors to file ... */

   for (i = 0; i < hptr->nprops - 1; ++i) {
      fprintf(fptr,"%s ", get_prop_mne(hptr->prop_id[i]));
   }
   fprintf(fptr,"%s", get_prop_mne(hptr->prop_id[i]));


   fprintf(fptr,"\n");

   return(0);
} /* end write_hydro_hdr() */

/**********************************************************************/

int write_hydro_scan(FILE *fptr, double *scan, int n, int *prop_id)

   /* writes a single data scan to an open hydro file.  Returns 0 for
      a successful write or a number > 0 if an error occurs.
           error code =   2 :  error writing to output file.
    arguments:
         fptr:    file opened for writing 
         scan:   array of property values in appropriate order 
         n:      number of properties 
        prop_id: index to enum property for each value in scan 

   */
{

   char *line, *s, str[10];
   int nbytes, i, jj, kk, k;

/* determine number of bytes needed for scan */

   nbytes = 0;
   for (i = 0; i < n; ++i) {
      nbytes += get_field_width(prop_id[i]) + 1;
   }
   ++nbytes;         /* allow for LF */
   line = (char *) calloc((size_t)nbytes, sizeof(char));
   s = line;

/* now write value in format appropriate to each property */

   for (i = 0; i < n; ++i ) {
      if (prop_id[i] < 0 || prop_id[i]%100 >= MAXPROP) {
          fprintf(stderr,"FATAL ERROR: prop_id index passed to write_hydro_scan() is out of range!!");
          exit(1);
      }
      sprintf(str," %c%d.%dlf", '%', get_field_width(prop_id[i]), 
              get_field_precis(prop_id[i]));
              
      /* some  values are too large to fit in the field width: */
      
      kk = get_field_width(prop_id[i])- get_field_precis(prop_id[i]) - 1;
           
      jj = 1;
      for (k = 2; k < kk; ++k) {
           jj *= 10;
      }
      if (scan[i] > (99.999 * jj) || scan[i] < (-9.999 * jj)) {
        scan[i] = HB_MISSING;
      }
        
      sprintf(s, str, scan[i]);
      s += get_field_width(prop_id[i]) + 1;
   } /* end for */

   line[nbytes-1] = '\0';
   fprintf(fptr, "%s", line);
   fprintf(fptr, "\n");

   free((void *)line);
   return (0);

}  /* end write_hydro_scan() */

/**********************************************************************/
int write_separator(FILE *fptr)

   /* writes a separator line to a hydro file.  Returns 0 for a successful
      write or a code > 0 if an error occurs.
   */
{
   char line[NBSEP];

   line[0] = '*';
   line[1] = '*';
   line[2] = '\0';
   fprintf(fptr, "%s", line);
   fprintf(fptr, "\n");
   return (0);     

}  /* end write_separator() */

/**********************************************************************/
int ms10(float lat, float lon, int *ms1_ptr)

   /*  Computes the 10-degree and 1-degree Marsden Square designations for a 
       particular lat/lon combination. */
{
   int  quadrant[2][2] = { 7000, 1000,
                           5000, 3000 };
   int  k, kk;
   int  ilat, ilon;


   k  = (lat < 0.) ? 1 : 0;   /* determine earth quadrant */
   kk = (lon < 0.) ? 0 : 1;
   if (lon >= 180.) {
      lon = lon - 360;
      kk = 0;
   }
   
   ilat = (int) (lat + .00001);  /* this smidgeon handles borderline */
   ilon = (int) (lon + .00001);  /*  cases correctly in all hemispheres */

   if (ilat < 0) ilat = -ilat;
   if (ilon < 0) ilon = -ilon;

   if (ilat >= 90) ilat = 89;    /* special case at the poles */

   *ms1_ptr = (ilat % 10) * 10 + (ilon % 10);

   return (quadrant[k][kk] + (ilat / 10) * 100 + (ilon / 10));

}  /* end ms10() */

/**********************************************************************/
int available(enum property prop, struct HYDRO_HDR *hptr)

   /*  Determines if the specified property is available at a particular
       dptr->  Returns 0 if not, 1 if available.  */
{
   int i, found;

   found = 0;
   i = -1;
   while (++i < hptr->nprops) {
         if (hptr->prop_id[i] == (int) prop)
             return (1);
   }
   return (0);

}  /* end available() */
/****************************************************************************/
int  is_in_range(float lon, float xmin, float xmax, int *merid, int lon0to360)

    /*  Returns 1 if lon falls within the longitide range xmin -> xmax
       or 0 if it does not
      for any combination of negative/positive longitude values.
      lon0to360 must be set upon entry into this function:
             = 1  if longitude range is all positive values
	     = 0   if mixed and crosses Greenwich meridian
	     = -1  if all negative values
       merid is an array used to quickly check whether a lon is in the interval
          by assigning a 0 or 1 values to each of its 360 elements. If merid is NULL
	  upon entry into the function, the array is created and its values are set
	  based upon the values of xmin,xmax and lon0to360.  Subsequent entries use
	  this array
	  
   */	   
{
   int i, index1, index2;
   
   if (merid == NULL) {
     merid = (int *) calloc(360, sizeof(int));
     
     if (lon0to360 > 0) {              /*lon range is all positive */
        index1 = NINT(xmin);
	index2 = NINT(xmax -1);
	for (i = index1; i <= index2; ++i) 
	    merid[i] = 1;
	
     }
     else if (lon0to360 < 0) {  /*lon range is all negative */
        index1 = NINT(xmin + 360.0);
	index2 = NINT(xmax + 360.0 - 1);
	for (i = index1; i <= index2; ++i) 
	    merid[i] = 1;
     
     }
     else {   /*  range crosses Greenwich */
          index1 = NINT(xmin + 360.0);
	  for (i = index1; i < 360; ++i) 
	     merid[i] = 1;
	  index2 = NINT(xmax - 1);
	  for (i = 0; i <= index2; ++i)
	     merid[i] = 1;
     }
   }  /* end if merid == NULL */
   
   i = NINT(lon);
   if (lon < 0)
       i += 360;
       
   if (i >= 360) 
       i = 359; 
         
   if (merid[i] ) 
        return (1);
   return (0);

}  /* end is_in_range() */
/****************************************************************************/

void free_and_alloc(double **dptr, int n)

   /* Frees up the space pointed to by dptr.  It MUST have been allocated
      using malloc()!!  Then mallocs space for n double values and sets
      dptr pointing to it.  If insufficient memory, a message is written to
      stderr and an exit() is made. */
{

   if (*dptr != NULL) {
      free((char *)*dptr);
      *dptr = NULL;
   }

   *dptr = (double *) malloc(n * sizeof(double));

   if (*dptr == NULL) {
      fprintf(stderr, "\nInsufficient memory for call to malloc()\n");
      exit(1);
   }

   return;
}  /* end free_and_alloc() */
/****************************************************************************/

void free_hydro_data(struct HYDRO_DATA *staptr)

   /* Frees all arrays and fields in struct HYDRO_DATA.  */
{
   int i;
   
    for (i = 0; i < MAXPROP; ++i) {
       if (staptr->observ[i] != NULL) {
          free((void *)staptr->observ[i]);
          staptr->observ[i] = NULL;
       }
       if (staptr->count[i] != NULL) {
          free((void *)staptr->count[i]);
          staptr->count[i] = NULL;
       }
       if (staptr->variance[i] != NULL) {
          free((void *)staptr->variance[i]);
          staptr->variance[i] = NULL;
       }
        if (staptr->quality[i] != NULL) {
          free((void *)staptr->quality[i]);
          staptr->quality[i] = NULL;
       }
   }    

   return;
}  /* end free_hydro_data() */
/****************************************************************************/
void list_origin(FILE *fptr)
  /* List out the origination code to file */
{
 int i;
 
 for (i = 0; i < NORIG; ++i ){
   fprintf(fptr,"  %c  %s\n", orig_code[i], orig_name[i]);
 }
 return;
}
/****************************************************************************/
void list_instrument(FILE *fptr)
  /* List out the instrument code to file */
{
 int i;
 
 for (i = 0; i < NINSTR; ++i ){
   fprintf(fptr,"  %c  %s\n", instr_code[i], instr_name[i]);
 }
 return;
}
/****************************************************************************/
void list_quality_flags(FILE *fptr)
  /* List out the quality code to file */
{
 int i;
 
 for (i = 0; i < NQUAL_FLAGS; ++i ){
   fprintf(fptr,"  %c  %s\n", quality_code_flag[i], quality_code_name[i]);
 }
 return;
}
/****************************************************************************/
void create_t68_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, .count, .variance .quality
     arrays for TE from T90 */
{
int i;
int *tmp;

  free_and_alloc(&dptr->observ[(int)TE], hptr->nobs);
  t90_to_t68(dptr->observ[(int)T90],dptr->observ[(int)TE],hptr->nobs);
  
  dptr->nprops = hptr->nprops;   /* make sure this is accurate */
  tmp = hptr->prop_id;
  hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
  for (i=0; i < dptr->nprops; ++i)
      hptr->prop_id[i] = tmp[i];
   hptr->prop_id[i] = (int) TE;
   free(tmp);   
   ++dptr->nprops;
	  
 if (dptr->count[(int)T90] != NULL) {
	    free_and_alloc(&dptr->count[(int)TE], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->count[(int)TE][i] = dptr->count[(int)T90][i];

	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Nte");
	    free(tmp);   
	    ++dptr->nprops;
 }
 if (dptr->variance[(int)T90] != NULL) {
	    free_and_alloc(&dptr->variance[(int)TE], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->variance[(int)TE][i] = dptr->variance[(int)T90][i];

	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Vte");
	    free(tmp);   
	    ++dptr->nprops;
 }
 if (dptr->quality[(int)T90] != NULL) {
	    free_and_alloc(&dptr->quality[(int)TE], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)TE][i] = dptr->quality[(int)T90][i];
	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Qte");
	    free(tmp);   
	    ++dptr->nprops;
 }

} /* end create_t68_arrays */
/*************************************/
void create_t90_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, count, variance quality
     arrays for T90 from T68 */
{
int i;
int *tmp;
         free_and_alloc(&dptr->observ[(int)T90], hptr->nobs);
	  t68_to_t90(dptr->observ[(int)TE],dptr->observ[(int)T90],hptr->nobs);
       
         dptr->nprops = hptr->nprops;   /* make sure this is accurate */
	 tmp = hptr->prop_id;
	 hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	 for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
	 hptr->prop_id[i] = (int) T90;
	 free(tmp);   
	 ++dptr->nprops;
	 
	 if (dptr->count[(int)TE] != NULL) {
	    free_and_alloc(&dptr->count[(int)T90], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->count[(int)T90][i] = dptr->count[(int)TE][i];

	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Nt90");
	    free(tmp);   
	    ++dptr->nprops;
	 }
	 if (dptr->variance[(int)TE] != NULL) {
	    free_and_alloc(&dptr->variance[(int)T90], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->variance[(int)T90][i] = dptr->variance[(int)TE][i];

	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Vt90");
	    free(tmp);   
	    ++dptr->nprops;
	 }
	 if (dptr->quality[(int)TE] != NULL) {
	    free_and_alloc(&dptr->quality[(int)T90], hptr->nobs);
	    for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)T90][i] = dptr->quality[(int)TE][i];

	    tmp = hptr->prop_id;
	    hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
	    for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
	    hptr->prop_id[i] = get_prop_indx("Qt90");
	    free(tmp);   
	    ++dptr->nprops;
	 }
} /* end create_t90_arrays */
/***************************************/
void create_pr_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, count, variance and quality
     arrays for pr from de */
{
int i;
int *tmp;

   if (dptr->observ[(int)DE] == NULL)
      return;
      
      
   free_and_alloc(&dptr->observ[(int)PR], hptr->nobs);
   for (i=0; i < hptr->nobs; ++i) 
       dptr->observ[(int)PR][i] = hb_p80(dptr->observ[(int)DE][i], (double)hptr->lat);
       
   dptr->nprops = hptr->nprops;   /* make sure this is accurate */
   tmp = hptr->prop_id;
   hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
   hptr->prop_id[0] = (int) PR;
   for (i=1; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
   free(tmp);   
   ++dptr->nprops;
	 
   if (dptr->count[(int)DE] != NULL) {
      free_and_alloc(&dptr->count[(int)PR], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
          dptr->count[(int)PR][i] = dptr->count[(int)DE][i];

      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("Npr");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->variance[(int)DE] != NULL) {
      free_and_alloc(&dptr->variance[(int)PR], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
         dptr->variance[(int)PR][i] = dptr->variance[(int)DE][i];

      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("Vpr");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->quality[(int)DE] != NULL) {
       free_and_alloc(&dptr->quality[(int)PR], hptr->nobs);
       for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)PR][i] = dptr->quality[(int)DE][i];

       tmp = hptr->prop_id;
       hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
       for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
       hptr->prop_id[i] = get_prop_indx("Qpr");
       free(tmp);   
       ++dptr->nprops;
   }
   return;
} /* end create_pr_arrays */
/*********************************************/
void create_de_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, count, variance and quality
     arrays for pr from de */
{
int i;
int *tmp;

   if (dptr->observ[(int)PR] == NULL)
      return;
      
      
   free_and_alloc(&dptr->observ[(int)DE], hptr->nobs);
   for (i=0; i < hptr->nobs; ++i) 
       dptr->observ[(int)DE][i] = hb_depth(dptr->observ[(int)PR][i], (double)hptr->lat);
       
   dptr->nprops = hptr->nprops;   /* make sure this is accurate */
   tmp = hptr->prop_id;
   hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
   for (i=0; i < dptr->nprops; ++i)
	hptr->prop_id[i] = tmp[i];
   hptr->prop_id[i] = (int) DE;
   free(tmp);   
   ++dptr->nprops;
	 
   if (dptr->count[(int)PR] != NULL) {
      free_and_alloc(&dptr->count[(int)DE], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
          dptr->count[(int)DE][i] = dptr->count[(int)PR][i];

      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("Nde");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->variance[(int)PR] != NULL) {
      free_and_alloc(&dptr->variance[(int)DE], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
         dptr->variance[(int)DE][i] = dptr->variance[(int)PR][i];

      dptr->nprops = hptr->nprops;   /* make sure this is accurate */
      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("Vde");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->quality[(int)PR] != NULL) {
       free_and_alloc(&dptr->quality[(int)DE], hptr->nobs);
       for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)DE][i] = dptr->quality[(int)PR][i];

       tmp = hptr->prop_id;
       hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
       for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
       hptr->prop_id[i] = get_prop_indx("Qde");
       free(tmp);   
       ++dptr->nprops;
   }
   return;
} /* end create_de_arrays */
/******************************************/
void create_o2_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, count, variance and quality
     arrays for O2 from OX */
{
int i, tindex;
int *tmp;

   if (dptr->observ[(int)OX] == NULL)
      return;
   if (dptr->observ[(int)PR] == NULL)
      return;
   if (dptr->observ[(int)SA] == NULL)
      return;
   tindex = (int)T90; 
   if (dptr->observ[(int)T90] == NULL)
      tindex = (int)TE;
   if (dptr->observ[tindex] == NULL)
      return;
      
   free_and_alloc(&dptr->observ[(int)O2], hptr->nobs);
   for (i=0; i < hptr->nobs; ++i) 
       dptr->observ[(int)O2][i] = ox_l2kg(dptr->observ[(int)OX][i], dptr->observ[(int)PR][i],dptr->observ[tindex][i], dptr->observ[(int)SA][i]);
       
   dptr->nprops = hptr->nprops;   /* make sure this is accurate */
   tmp = hptr->prop_id;
   hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
   for (i=0; i < dptr->nprops; ++i)
	hptr->prop_id[i] = tmp[i];
   hptr->prop_id[i] = (int) O2;
   free(tmp);   
   ++dptr->nprops;
	 
   if (dptr->count[(int)OX] != NULL) {
      free_and_alloc(&dptr->count[(int)O2], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
          dptr->count[(int)O2][i] = dptr->count[(int)OX][i];

      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("No2");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->quality[(int)OX] != NULL) {
       free_and_alloc(&dptr->quality[(int)O2], hptr->nobs);
       for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)O2][i] = dptr->quality[(int)OX][i];

       tmp = hptr->prop_id;
       hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
       for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
       hptr->prop_id[i] = get_prop_indx("Qo2");
       free(tmp);   
       ++dptr->nprops;
   }
   return;
} /* end create_o2_arrays */

/******************************************/
void create_ox_arrays(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr)
  /* generates dptr->observ, count,  and quality
     arrays for OX from O2 */
{
int i, tindex;
int *tmp;

   if (dptr->observ[(int)O2] == NULL)
      return;
   if (dptr->observ[(int)PR] == NULL)
      return;
   if (dptr->observ[(int)SA] == NULL)
      return;
   tindex = (int)T90; 
   if (dptr->observ[(int)T90] == NULL)
      tindex = (int)TE;
   if (dptr->observ[tindex] == NULL)
      return;
      
   free_and_alloc(&dptr->observ[(int)OX], hptr->nobs);
   for (i=0; i < hptr->nobs; ++i) 
       dptr->observ[(int)OX][i] = ox_kg2l(dptr->observ[(int)O2][i], dptr->observ[(int)PR][i],dptr->observ[tindex][i], dptr->observ[(int)SA][i]);
       
   dptr->nprops = hptr->nprops;   /* make sure this is accurate */
   tmp = hptr->prop_id;
   hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
   for (i=0; i < dptr->nprops; ++i)
	hptr->prop_id[i] = tmp[i];
   hptr->prop_id[i] = (int) OX;
   free(tmp);   
   ++dptr->nprops;
	 
   if (dptr->count[(int)O2] != NULL) {
      free_and_alloc(&dptr->count[(int)OX], hptr->nobs);
      for (i=0; i < hptr->nobs; ++i)
          dptr->count[(int)OX][i] = dptr->count[(int)O2][i];

      tmp = hptr->prop_id;
      hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
      for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
      hptr->prop_id[i] = get_prop_indx("Nox");
      free(tmp);   
      ++dptr->nprops;
   }
   if (dptr->quality[(int)O2] != NULL) {
       free_and_alloc(&dptr->quality[(int)OX], hptr->nobs);
       for (i=0; i < hptr->nobs; ++i)
	      dptr->quality[(int)OX][i] = dptr->quality[(int)O2][i];

       tmp = hptr->prop_id;
       hptr->prop_id = (int *)calloc(++hptr->nprops, sizeof(int));
       for (i=0; i < dptr->nprops; ++i)
	        hptr->prop_id[i] = tmp[i];
		
       hptr->prop_id[i] = get_prop_indx("Qox");
       free(tmp);   
       ++dptr->nprops;
   }
   return;
} /* end create_ox_arrays */
