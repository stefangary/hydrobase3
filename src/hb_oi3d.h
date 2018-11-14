/*  hb_oi3d.h
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             2011
			     
................................................................................
.
.  Contains the definition for functions and structures  used
.   in optimal interpolation modules.
.   
.          
.
*/

#ifndef HB_oi3d
#define HB_oi3d 1

 
#define NCPARMS_FILL -9.0e034
#define NCPARMS_MASK  9.0e034
#define NMONTHS 13  /* 12 months + clim */

#ifndef NINT
#define   NINT(x)       (((x) < 0) ? ((int)((x) - 0.5)) : ((int)((x) + 0.5)))
#endif

/* By definition, these grids are pixel registered and
   all longitudes are on the range 0 to 360 */

struct NC_INFO {
    int nx;             /* Number of columns */
    int ny;             /* Number of rows */
    int nz;             /* Number of depth levels */
    int nz_seas;        /* Number of levels in seasonal layer */
    int modelCode;      /* for semi-variograms */
    float *lats;         /* Vector of latitude coordinates */
    float *lons;         /* Vector of longitude coordinates */
    float xincr;        /* x increment */
    float yincr;        /* y increment */
    float fill_value;   /* to denote missing data */
    float mask_value;   /* to denote seafloor mask */
    char *x_units;   /* units in x-direction */
    char *y_units;   /* units in y-direction */
    char *z_units;   /* units in z-direction */  
};

struct NC_PARMS {
     int nz;                   /* Number of depth levels */
     int nz_seas;              /* and for seasonal layer */
     float dZdx, dZdy;         /* bathymetry gradient */  
     float *depths;            /* depth at each level */ 
     float *Lx, *Ly;          /* x, y length scales */
     float *dTdx, *dTdy;      /* temperature gradient */
     float *dPdx, *dPdy;      /* pressure gradient */  
     float *Tparm0, *Tparm1, *Tparm2;  /* fitted parameters */
     float *Sparm0, *Sparm1, *Sparm2;
     float *Pparm0, *Pparm1, *Pparm2;
};


/*  prototypes for functions in oi3d_utils.c */


extern int ncparms_init(char *);
extern int ncparms_open(char *);
extern int ncparms_update(char *);
extern void ncparms_close(int);

extern int ncparms_define(int, struct NC_INFO *, float *);
extern int ncparms_write(int, float, float, struct NC_INFO *, struct NC_PARMS *);
extern int ncparms_write_coord(int, char *, float *, int );
extern int ncparms_write_var(int, char *, size_t *, size_t *, float *);

extern int ncparms_getinfo(int, struct NC_INFO *, struct NC_PARMS *, int);
extern int ncparms_read(int, float, float, struct NC_INFO *, struct NC_PARMS *);
extern int ncparms_getcoord(int, char *, char *, float *, char *);
extern int ncparms_getvar(int, char *, int, int, int, float *);

extern int ncparms_rowcol( struct NC_INFO *, float, float, int *, int *);
extern int ncparms_latlon( struct NC_INFO *, int, int, float *, float *);
extern int ncparms_depth_indx(float, float *, int);

extern void ncparms_zero(struct NC_INFO *, struct NC_PARMS *);
extern void ncparms_alloc(struct NC_INFO *, struct NC_PARMS *);
extern void ncparms_free(struct NC_INFO *, struct NC_PARMS *);
extern void ncparms_prefill(struct NC_PARMS *, float);

#endif /*ifndef hb_oi3d */
