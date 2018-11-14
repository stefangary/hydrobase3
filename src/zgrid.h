/* hb_zgrid.h
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
*/


#ifndef HB_FIT_PARMS
#define HB_FIT_PARMS 1


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

#define ZGRID_EMPTY 1.0e+35
#define MEM_CHUNK  2000
#define SSIZE       500    /* arbitrary work array size */


 
#ifndef NINT
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#endif

#ifndef ABS
#define    ABS(x)       (((x) < 0) ? -(x) : (x))
#endif


struct POINT {               /* Stores input data */
  double XP, YP, ZP;
} ;


extern double * zgrid (double, int, int, int, int, double,  struct GRID_INFO *, int *, int *, double *, double *, struct POINT *, int *, int *, int *);


#endif  /* ifndef HB_FIT_PARMS */


