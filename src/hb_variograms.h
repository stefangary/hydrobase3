/*  hb_variograms.h
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             2011
			     
................................................................................
.
.  Contains the definition for functions and structures used
.   in constructing semi-variograms.
.           
.
*/
#ifndef HB_variogram
#define HB_variogram 1


#define N_VGRAM_MODELS 3   /* for computing semi-variograms */

/*  
    model = 1 : exponential
          = 2 : gaussian
	  = 3 : power
 */

struct VARIOGRAM {
   int nlags;   /* number of distance bins */
   double *Vh;  /* variogram values for each lag */
   double *h;   /* average distance for each bin */
   int *Nh;     /* Number of points in each bin */
};


/* prototypes for functions in vgram_utils.c */

extern double *vgram_vals(int, double *, int, double *);
extern int vgram_vals2(int, double *, int, double *, double *);
extern struct VARIOGRAM *vgram_estimate(double **, int, float *, float *, int, int, double, double);
extern int get_parm_range(int, struct VARIOGRAM *, double *, double *);
extern int vgram_check(struct VARIOGRAM *, int, int, int);

#endif /*ifndef HB_variogram*/
