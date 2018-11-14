/*   hydro_cdf.h
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated to conform to ANSI standards Nov 1999
			     Revamped for version 3 Dec 2009
................................................................................
.
.  Contains the definition for netcdf files used with HydroBase 3
.       Main features added:
              variables to store an error variance (or std deviation)
.           
.
*/

#ifndef HYDRO_CDF
#define HYDRO_CDF 1

#define COUNT_VAR_SUFFIX  "_cnt"   /* nomenclature for count variables */
#define ERR_VAR_SUFFIX  "_err"   /* nomenclature for error variables  */
#define MAXSTDLEVS  1000       /* max # of standard depth levels */

#define HBMASK  0.9e35          /* marks a masked node in cdf file */
#define HBEMPTY -0.9e35        /* marks an empty node in cdf file */


/*  Global variables */

double std_depth[1000];  /* initialize this with std_depth_init()  */
int  NSTDLEVS ;                /* actual number of standard depth levels to be
                                  set in std_depth_init() */

int  std_depth_initialized;    /* set to 1 after initializing */

    /****************************************************************/

struct CDF_HDR {
    int nx;             /* Number of columns */
    int ny;             /* Number of rows */
    int nz;             /* Number of depth levels */
    int nt;             /* Number of time bins */
    int nprops;         /* Number of properties */
    int counts_included; /* 0 = no stats (counts + errs)
                            1 = counts only
			    3 = counts + errs) */
    int node_offset;    /* 0 for node grids, 1 for pixel grids */
    float fill_value;   /* to denote missing data */
    float mask_value;   /* to denote masked nodes (no ocean) */
    float xmin;         /* Minimum x coordinate */
    float xmax;         /* Maximum x coordinate */
    float ymin;         /* Minimum y coordinate */
    float ymax;         /* Maximum y coordinate */
    float xincr;        /* x increment */
    float yincr;        /* y increment */
    int *tmin;          /* min year for each time bin */
    int *tmax;          /* max year for each time bin */
    char x_units[80];   /* units in x-direction */
    char y_units[80];   /* units in y-direction */
    char z_units[80];   /* units in z-direction */
    char t_units[80];   /* units in time dimension */
    char **prop_id;     /* char mnemonic of various properties */
    char **prop_units;  /* units of various properties */
    char title[80];     /* name of data set */
    char command[5000];  /* name of generating command */
};

   /*  Description of node_offset:

       Given: xmin = ymin = 0
              xmax = ymax = 10 
              xincr = yincr = 1

        For a node grid we have:
                (1) nx = (xmax - xmin) / xincr + 1 = 11
                    ny = (ymax - ymin) / yincr + 1 = 11
                (2) node # 0 is at (x,y) = (xmin, ymax) = (0,10) and 
                    represents the value for a box of size (1,1)
                    centered on the node.

        For a pixel grid we have:
                (1) nx = (xmax - xmin) / xincr = 10
                    ny = (ymax - ymin) / yincr = 10
                (2) node # 0 is at (x,y) = (xmin + 0.5*xincr, ymax - 
                    0.5*yincr) = (0.5, 9.5) and represents the value  
                    for a box with dimensions (1,1) centered on the node.
*/


/*  prototypes of functions defined in hydro_cdf.c */

extern int cdf_init(char *);
extern int cdf_open(char *, char *, char *, int);
extern int cdf_update(char *, char *, char *, int);
extern void cdf_close(int);
extern int cdf_define(int, struct CDF_HDR *, int, int);
extern int write_prop_cdf(int, float *, char *, int, int, int, int, int, int, int, int);
extern int write_prop_err_cdf(int, float *, char *, int, int, int, int, int, int, int, int);
extern int write_prop_count_cdf(int, short *, char *, int, int,int, int, int, int, int, int);
extern int write_lat_vector(int, struct CDF_HDR *);
extern int write_lon_vector(int, struct CDF_HDR *);
extern int write_std_depths_cdf(int, struct CDF_HDR *);
extern int write_bottom_depth_cdf(int, int, int, int, int, int, int, float *);
extern int write_time_bins_cdf(int, struct CDF_HDR *);
extern int read_cdf_hdr(int, struct CDF_HDR *);
extern int read_lat_vector(int, float *);
extern int read_lon_vector(int, float *);
extern int read_cdf_depths(int, float *);
extern int read_cdf_bottom_depth(int, float *, int, int, int);
extern int read_cdf_prop(int, char *, float *, int, int, int, int, int);
extern int read_cdf_prop_err(int, char *, float *, int, int, int, int, int);
extern int read_cdf_prop_count(int, char *, short *, int, int, int, int, int);
extern int get_indices(struct CDF_HDR *, float, float, int *, int *);
extern int get_lat_lon(struct CDF_HDR *, int, int, float *, float *);
extern int std_depth_init(FILE *);
extern int d2stdlev(double);
extern double stdlev2depth(int);
extern int is_flagged(float, float);


#endif /*ifndef HYDRO_CDF */
