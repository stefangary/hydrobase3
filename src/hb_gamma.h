/*  hb_gamma.h
................................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Mar 2001
...................................................
*/


/* macros */

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))     
#endif

#ifndef ABS
#define    ABS(x)       (((x) < 0) ? -(x) : (x))
#endif

#ifndef NINT
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#endif

#define MAX_SCV 100           /* arbitrary work array sizes */
#define MAX_INTVLS 50  

/* dimensions of netcdf arrays */

#define NX 90  /* # of longitudes */
#define NY 45  /* # of latitudes */
#define NZ 33  /* # of std pressure levels */
#define NDX 4  /* size of x-grid spacing in gamma.nc  */
#define NDY 4  /* size of y-grid spacing in gamma.nc  */

/* structure to store info about the netCDF file gamma.nc */

struct GAMMA_NC {
   int id_gnc;   /* netcdf file descriptor */
   int id_s;     /* netcdf variable ids */
   int id_t;
   int id_gamma;
   int id_a;
   int iocean[NY][NX]; 
   int n[NY][NX]; 
   double along_d[NX];
   double alat_d[NY];
   double p0[NZ];
};


/*********************************************
   Prototypes for functions in gamma_subs.c 
************************************************/

   /* This function must be called before any others 
      to initialize information from
      the netcdf file gamma.nc */
  
extern void gamma_nc_init(char *, struct GAMMA_NC *);

   /* The next 2 functions are the primary entry points 
      to label each level in a cast of p,t,s with 
      neutral density:  with or without error bars */
   
extern void compute_gamma_n(struct GAMMA_NC *, int, double *,  double *, double *, double *, double, double);

extern void compute_gamma_nerr(struct GAMMA_NC *, int, double *, double *, double *, double *, double *, double, double);

   /* This is the entry point to compute 
      properties on a neutral surface */
   
void neutral_surfaces(double *, double *, double *, double *, int, double *, int, double *, double *, double *, double *, double *, double *);


   /* The remaining functions are called internally 
      by the above functions... */

extern void depth_ns(double *, double *, double *, int, double, double, double, double *, double *, double *);

extern void depth_scv(double *, double *, double *, int, double, double, double, double *, double *, double *, int *);

extern void e_solve(double *, double *, double *, double *, int, int, double, double, double, double *, double *, double *, int *);

extern void eosall(double, double, double, double *, double *, double *, double *, double *, double *);

extern void gamma_errors(double *, double *, double *, double *, double *, int, double, double, double, double, double, double, double, double, int, double, double *, double *, double *);

extern int gamma_n(struct GAMMA_NC *, double *, double *, double *, int, double, double, double *, double *, double *);

extern double gamma_qdr(double, double, double, double, double, double);

extern void goor(double *, double *, double *, double *, int, double, double, double, double *, double *, double *, double *);

extern void goor_solve(double, double, double, double, double, double, double, double, double, double, double, double *, double *);

extern int indx(double *, int, double);

extern int ocean_test(double, double, int, double, double, int, double);

extern void read_nc(struct GAMMA_NC *, double, double, double ***, double ***, double *, double ***, double ***, int **, double *, double *, int **);

extern void scv_solve(double *, double *, double *, double *, int, int, double, double, double, double *, double *, double *, int *);

extern void sig_vals(double, double, double, double, double, double, double *, double *);

extern void stp_interp(double *, double *, double *, int n, double *, double *, double *);
