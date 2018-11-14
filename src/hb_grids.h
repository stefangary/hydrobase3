/* hb_grids.h
................................................................................
                          *******  HydroBase3  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2010
................................................................................
*/


#ifndef HB_GRID_PARMS
#define HB_GRID_PARMS 1


#ifndef TRUE
#define TRUE  1
#endif

#ifndef  FALSE
#define  FALSE 0
#endif

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))     /* min and max value macros */
#endif

#ifndef TOOBIG
#define TOOBIG 0.9e+35     /* value larger than any oceanographic property */
#endif                     /* but small enough to fit in a float variable */

#define MEM_CHUNK  2000
#define SSIZE       500    /* arbitrary work array size */

#ifndef NINT
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#endif

#ifndef ABS
#define    ABS(x)       (((x) < 0) ? -(x) : (x))
#endif

#define PI     3.141592654
#define RADperDEG 0.017453292             /* pi/180 */
#define DEGperRAD 57.29577951    /* 180/pi */
#define EarthRadius  3438    /* in nm */
#define EarthRadiusKm  6371   /* in km */
#define NMperKM  .539593
#define KMperNM   1.853248652
#define NCROSS     50              /* used by inside() */


typedef int BOOLEAN;                    /* used for logical variables */

struct GRID_INFO {
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	int node_offset;		/* 0 for node grids, 1 for pixel grids */
	int lon0to360;                  /* 1 if lons all positive */
	int xgreenwich;                 /* 1 if lons are mixed sign */
	double x_min;			/* Minimum x coordinate */
	double x_max;			/* Maximum x coordinate */
	double y_min;			/* Minimum y coordinate */
	double y_max;			/* Maximum y coordinate */
	double x_inc;			/* x increment */
	double y_inc;			/* y increment */

};


/* in grid_subs.c */
extern int ij2xy(struct GRID_INFO *, int, int, double *, double *);
extern int xy2ij(struct GRID_INFO *, double, double, int *, int *);
extern void ij2xy_nochk(struct GRID_INFO *, int, int, double *, double *);
extern void xy2ij_nochk(struct GRID_INFO *, double, double, int *, int *);
extern int sq2rc(int, struct GRID_INFO *, int *, int *);
extern int sq2latlon(int, struct GRID_INFO *, double *, double *);
extern int inside (double, double, double *, double *, int);
extern int get_mask(FILE *, struct GRID_INFO *, char *, BOOLEAN);
extern void nn_interp2d(double *, double *,  double, double, int , int , struct GRID_INFO *, double *, double *, double *, int *);
extern double weighted_mean(double *, double *, double, double,int);
double weighted_variance(double *, double *, double, double, double, int);
extern void zero_weights_2d(double *x, double *wghts,double empty_val, double mask_val, struct GRID_INFO *hptr, int xmid, int ymid);
extern void set_to_zero(int *, double *, double *, int, double, double, double);

extern double distance_c(double, double, double, double, int, double *);
extern int pointB(double, double, double, double, int, double *, double *);
extern int pointB_hav(double, double, double, double, int, double *, double *);
extern void get_weights_c(double *, struct GRID_INFO *, double, double, double, int, int);
extern void get_weights_g( double *, struct GRID_INFO *, double, int, int);


/*   in topo_subs.c */

extern short *hb_get_topo(char *, struct GRID_INFO *, double **, double **, int, int, short *);
extern short find_nearest_topo_val(double, double, short *, struct GRID_INFO *);

#endif  /* ifndef HB_GRID_PARMS */


/*
-----------------------------------------------------------------------------------------
 	Notes on node_offset:

	Assume x_min = y_min = 0 and x_max = y_max = 10 and x_inc = y_inc = 1.
	For a normal node grid we have:
		(1) nx = (x_max - x_min) / x_inc + 1 = 11
		    ny = (y_max - y_min) / y_inc + 1 = 1
		(2) node # 0 is at (x,y) = (x_min, y_max) = (0,10) and represents the surface
		    value in a box with dimensions (1,1) centered on the node.
	For a pixel grid we have:
		(1) nx = (x_max - x_min) / x_inc = 10
		    ny = (y_max - y_min) / y_inc = 10
		(2) node # 0 is at (x,y) = (x_min + 0.5*x_inc, y_max - 0.5*y_inc) = (0.5, 9.5)
		    and represents the surface value in a box with dimensions (1,1)
		    centered on the node.
-------------------------------------------------------------------------------------------*/
