/*   phyprops.c
...........................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             
............................................................................
. 
.    Algorithms for oceanographic computations 
             by N Fofonoff & R Millard
.    Converted from FORTRAN by Julie Allen, 9/99
............................................................................
. 
.   	double hb_sal78(double cond,double t,double p,int iflag) 
.    		salt from cond ratio or vice versa
.
.   	double hb_svan(double s,double tref,double pref,double *sigma_ptr) 
.    		specific volume anomaly
.
.   	double hb_depth(double p,double lat) 
.    		depth from pressure
.
.   	double hb_tf(double s,double p)
.    		freezing point of seawater
.
.  	double hb_cpsw(double s, double t, double p0)
.    		specific heat of seawater
.
.  	double hb_atg(double s,double t,double p)
.    		adiabatic temp gradient
.
.  	double hb_theta(double s,double t0,double p0,double prref)
.    		potential temperature referenced to pref
.
.  	double hb_svel(double s,double t,double p0)
.    		speed of sound in seawater
.
.  	double hb_bvfrq(double *s,double *t,double *p,int n,double *pav,double *e_ptr)
.    		brunt-vaisala frequency from sp.vol.anom.
.
.  	double hb_bvgrdts(double *p0_ptr,double *s,double *t,double *p,int n,double *e_ptr)
.    		brunt-vaisala frequency at a specified 
. 		pressure from t, s gradients
.
.  	double hb_grady(double *y,double *p,int nobs,double *pav_ptr,double *ybar_ptr)
.    		least squares slope of y vs. p
.
.  	double hb_p80(double dpth,double xlat)
.    		pressure from depth using eos80
.
.  	double hb_gravity(double xlat)
.    		grav acceleration at latitude (cm/sec**2)
.
.  	double hb_coriol(double xlat)
.    		coriolus parameter from latitude
.
. 	double hb_linterp(double xval, double *x, double *y,int npts)
.    		linearly interpolates to find the position of xval in array x
. 
. ........................................................................
. 
.  Used functions (rather than macros) to replace the Fortran internal 
.	functions sal, dsal and rt35 
.     sal(xr,xt) : 
. 	practical salinity scale (1978 definition w/ temp correction)
.                   xt = t - 15.0
.                   xr = sqrt(rt) 
*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "hydrobase.h"

/* define constants used locally */
#define R3500 1028.1063
#define R4 4.8314e-4		/* referred to as C in Millero and Poisson 1981 */
#define DR350 28.106331
#define RPS2CPH 572.9578 	/* 3600/2*pi changes rad/sec to cph */
#define GP 9.80655		/* Gravitational Constant in m/sec^2 */
#define PI 3.141592654
#define FLAG -99999.


/* Local functions, not labeled as hb_ */

double sal(double xr, double xt)
{
    return ((((2.7081*xr-7.0261)*xr+14.0941)*xr+25.3851)*
	xr-0.1692)*xr+0.0080+(xt/(1.0+0.0162*xt))*
	(((((-0.0144*xr+0.0636)*xr-0.0375)*xr-0.0066)*xr-0.0056)*xr+0.0005);
}

/* dsal(xr,xt) : returns derivative of sal(xr,xt) with respect to xr */
double dsal(double xr, double xt)
{
    return ((((13.5405*xr-28.1044)*xr+42.2823)*xr+50.7702)*
	xr-0.1692)+(xt/(1.0+0.0162*xt))*
	((((-0.0720*xr+0.2544)*xr-0.1125)*xr-0.0132)*xr-0.0056);
}

/* rt35(xt)  : c(35,t,0) / c(35,15,0) variation with temperature */
double rt35(double xt)
{
    return ((((1.0031e-9*xt-6.9698e-7)*xt+1.104259e-4)*
	xt+2.00564e-2)*xt+0.6766097);
}

/* polynomials of rp: c(s,t,p)/c(s,t,0) variation with pressure
.   c(xp) polynomial corresponds to a1-a3 constants: lewis 1980
.   a(xt) polynomial corresponds to b3 and b4 constants: lewis 1980 */
double c(double xp)
{
      return (((3.989e-15*xp-6.370e-10)*xp+2.070e-5)*xp);
}

double b(double xt) 
{
    return ((4.464e-4*xt+3.426e-2)*xt + 1.0);
}

double a(double xt)
{
     return (-3.107e-3*xt + 0.4215);
}

/*********************************************************/
/*********************************************************/
double hb_sal78(double x,double t,double p,int iflag)
{
    double xr,xt,xp,si,dels;
    double salt,set,rtt,at,bt,cp,r;
    double rt,dt;
    double sal78;
    int n;
/* 
.  function to convert conductivity ratio to salinity (iflag = 0)
.                   OR salinity to conductivity ratio (iflag = 1)
.      references:   also located in unesco report # 37 1981
.         practical salinity scale 1978: e.l. lewis ieee ocean eng. jan. 1980
. ...........................................................................
.    args:
.        x      :  conductivity ratio (iflag = 0) OR salinity (iflag = 1)
.        p      :  pressure ( decibars)
.        t      :  temperature  (deg celsius)
.        iflag  :  =0 (to convert cond ratio to salinity)
.                  =1 (to convert salinity to cond ratio)
.    returns salinity (iflag = 0)
.         or conductivity ratio (iflag = 1)
. ...........................................................................
.   check for values too small...
.       returns zero for cond ratio < 0.0005     (iflag = 0)
.                 or for salinity < 0.02         (iflag = 1)
*/

      sal78 = 0.0;
      if((iflag == 0) && (x <= 5e-4)) return sal78;
      if((iflag == .1) && (x <= 0.02)) return sal78;

      dt = t - 15.0;

/* convert conductivity ratio to salinity */
      if(iflag == 0) { 
         r = x;
         rt = r/(rt35(t) * (1.0 + c(p)/(b(t) + a(t)*r)));
         rt = sqrt(fabs(rt));
         sal78 = sal(rt,dt);
         return sal78;
      }

/* salinity to conductivity ratio
.  invert salinity to conductivity by the
.   newton-raphson iterative method.  */
      salt = x;

/*  set initial values ... */
      rt = sqrt(salt/35.0);
      si = sal(rt,dt);
      n = 0;
      dels = 1.;

/* iteratively refine salinity inversion ...  */
     while ((dels > 1.0e-4) && (n < 10))
     {
         rt = rt + (salt - si)/dsal(rt,dt);
         si = sal(rt,dt);
         n = n + 1;
         dels = fabs(si - salt);
    }

/* compute conductivity ratio ...  */
      rtt = rt35(t) * rt * rt;
      at = a(t);
      bt = b(t);
      cp = c(p);
      cp = rtt * (cp + bt);
      bt = bt - rtt * at;

/* solve quadratic equation for r: r = rt35 * rt * (1+c/ar+b) */
      r = sqrt(fabs(bt*bt + 4.0*at*cp)) - bt;

/* return conductivity */
      sal78 = 0.5 * r / at;
      return sal78;
}

/*********************************************************/
/*********************************************************/


double hb_svan(double s,double t,double p0,double *sigma_ptr) 
/*  ...............................................................
.  Specific Volume Anomaly (steric anomaly) based on 1980 equation
.  of state for seawater and 1978 practical salinity scale.
.  References:
.         Millero, et al (1980) deep-sea res.,27a,255-264
.         Millero and Poisson 1981,deep-sea res.,28a pp 625-629.
.    (both  references are also found in unesco report 38 -- 1981)
.  
.   The type of density anomaly (sigma) computed depends on the type of
.   temperature and pressure supplied:
. 
.    in situ density:  t is in situ temperature
.                      p0 is in situ pressure
. 
.    potential density: p0 is reference pressure
.                       t is temperature referenced already to p0
. 
.    sigma-t: : p0 = 0
.               t is in situ temperature
.  ................................................................
.  units:
.        p0  :    pressure  [or ref pressure] (decibars)
.         t  :    temperature  [or pot temp]  (deg C)
.         s  :    salinity (ipss-78)
.      svan  :    spec. vol. anom.  (m**3/kg *1.0e-8)
.      sigma :    density anomaly   (kg/m**3)
.  ................................................................
.  check value: svan=981.3021 e-8 m**3/kg.  for s = 40 (ipss-78) ,
.  t = 40 deg c, p0= 10000 decibars.
.  check value: sigma = 59.82037  kg/m**3 for s = 40 (ipss-78) ,
.  t = 40 deg c, p0= 10000 decibars. 
. ............................................................... */
{
    double p,sig,sr,r1,r2,r3;
    double a,b,c,d,e,a1,b1,aw,bw,k,k0,kw,k35;
    double svan,sva,v350p,dk,gam,pk,dr35p,dvan;

/*..............................................................
.  convert pressure to bars and take square root of salinity 
. ..............................................................   */
    p = p0 / 10.;
    sr = sqrt(fabs(s));

/*..............................................................
.  pure water density at atmospheric pressure
.    bigg p.h.,(1967) br. j. applied physics 8 pp 521-537 
. ..............................................................   */
    r1 = ((((6.536332e-9* t-1.120083e-6) * t +1.001685e-4) 
    	* t -9.095290e-3) * t +6.793952e-2) * t -28.263737;
    
/*.............................................................. 
.  seawater density at atmospheric press.
.     coefficients involving salinity:
.       r2 = A   in notation of Millero and Poisson 1981
.       r3 = B 
. ..............................................................   */
    r2 = (((5.3875e-9 * t -8.2467e-7) * t +7.6438e-5)
	* t -4.0899e-3) * t +8.24493e-1;
    r3 = (-1.6546e-6 * t +1.0227e-4) * t -5.72466e-3;

/*..............................................................
.  international one-atmosphere equation of state of seawater 
. ..............................................................   */
    sig = (R4*s + r3*sr + r2)*s + r1;

/*.............................................................. 
.   specific volume at atmospheric pressure 
. ..............................................................   */
    v350p = 1.0 / R3500;
    sva = -sig * v350p / (R3500+sig);
    *sigma_ptr = sig + DR350;

/*.............................................................. 
.   scale specific vol. anomaly to normally reported units 
. ..............................................................   */
    svan = sva * 1.0e+8;
    if (p == 0.) return svan;

/*..............................................................
.   high pressure equation of state for seawater
.   Millero, et al , 1980 dsr 27a, pp 255-264
.   constant notation follows article
.  
.   compute compression terms ...  
. ..............................................................   */
    e = (9.1697e-10 * t +2.0816e-8) * t -9.9348e-7;
    bw = (5.2787e-8 * t -6.12293e-6) * t +3.47718e-5;
    b = bw + e * s;

    d = 1.91075e-4;
    c = (-1.6078e-6 * t -1.0981e-5) * t +2.2838e-3;
    aw = ((-5.77905e-7 * t +1.16092e-4) * t +1.43713e-3) * t -0.1194975;
    a = (d * sr + c) * s + aw;

    b1 = (-5.3009e-4 * t +1.6483e-2) * t +7.944e-2;
    a1 = ((-6.1670e-5 * t +1.09987e-2) * t -0.603459) * t +54.6746;
    kw = (((-5.155288e-5 * t +1.360477e-2) * t -2.327105) *
	t +148.4206) * t -1930.06;
    k0 = (b1*sr + a1) * s + kw;

/*.............................................................. 
.  evaluate pressure polynomial 
. 
.   k equals the secant bulk modulus of seawater
.   dk = k(s,t,p) - k(35,0,p)
.   k35 = k(35,0,p) 
. ..............................................................   */
    dk = (b * p + a) * p + k0;
    k35 = (5.03217e-5*p+3.359406)*p+21582.27;
    gam = p / k35;
    pk = 1.0 - gam;
    sva = sva*pk + (v350p+sva)*p*dk/(k35*(k35+dk));

/*.............................................................. 
.   scale specific vol. anamoly to normally reported units...  
. ..............................................................   */
    svan = sva*1.0e+8;
    v350p = v350p * pk;

/*.............................................................. 
.  compute density anamoly with respect to 1000.0 kg/m**3
.  1) DR350: density anomaly at 35 (ipss-78), 0 deg. c and 0 decibars
.  2) dr35p: density anomaly 35 (ipss-78), 0 deg. c ,  pres. variation
.  3) dvan : density anomaly variations involving specfic vol. anamoly 
. ..............................................................   */

    dr35p = gam / v350p;
    dvan = sva / (v350p * (v350p+sva));
    *sigma_ptr = DR350 + dr35p - dvan;

    return svan;
}
/*********************************************************/
/*********************************************************/
double hb_depth(double p,double lat) 
/*********************************************************/
/* depth in meters from pressure in decibars using
.    saunders and fofonoff's method.
.        deep-sea res., 1976,23,109-111.
.    formula refitted for 1980 equation of state
. 
.        pressure        p        decibars
.        latitude        lat      degrees
.        depth           depth    meters
. 
.  checkvalue: depth = 9712.653 m for p=10000 decibars, latitude=30 deg
.      above for standard ocean: t=0 deg. celsius ; s=35 (ipss-78)
. ..................................................................   */
{
    double x,gr,depth;

    x = sin(lat/57.29578);
    x = x * x;

/* gr = gravity variation with latitude: anon (1970) bulletin geodesique */

    gr = 9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x) + 1.092e-6*p;
    depth = (((-1.82e-15*p+2.279e-10)*p-2.2512e-5)*p+9.72659)*p;
    return (depth / gr);
}

/*********************************************************/
/*********************************************************/
double hb_tf(double s,double p)
/*********************************************************/
/*   function to compute the freezing point of seawater
. 
.    reference: unesco tech. papers in the marine science no. 28. 1978
.    eighth report jpots
.    annex 6 freezing point of seawater f.j. millero pp.29-35.
. 
.   units:
.          pressure      p          decibars
.          salinity      s          pss-78
.          temperature   tf         degrees celsius
.          freezing pt.
. ..................................................................
.   checkvalue: tf= -2.588567 deg. c for s=40.0, p=500. decibars */
{
    return (-.0575+1.710523e-3*sqrt(fabs(s))-2.154996e-4*s)*s-7.53e-4*p;    
}

/*********************************************************/
/*********************************************************/
double hb_cpsw(double s, double t, double p0)
/*********************************************************/
/*.................................................................
.        pressure        p0       decibars
.        temperature     t        deg celsius (ipts-68)
.        salinity        s        (ipss-78)
.        specific heat   cpsw     j/(kg deg c)
. ..................................................................
.  ref: millero et al,1973,jgr,78,4499-4507
.        millero et al, unesco report no. 38 1981 pp. 99-188.
.  pressure variation from least squares polynomial
.  developed by fofonoff 1980.
.  check value: cpsw = 3849.500 j/(kg deg. c) for s = 40 (ipss-78),
.  t = 40 deg c, p0= 10000 decibars
.    scale pressure to bars */
{
    double a,b,c,p,sr,cp0,cp1,cp2;

    p=p0/10.;

/* sqrt salinity for fractional terms */
    sr = sqrt(fabs(s));

/* specific heat cp0 for p=0 (millero et al ,unesco 1981) */
    a = (-1.38385e-3*t+0.1072763)*t-7.643575;
    b = (5.148e-5*t-4.07718e-3)*t+0.1770383;
    c = (((2.093236e-5*t-2.654387e-3)*t+0.1412855)*t-3.720283)*t+4217.4;
    cp0 = (b*sr + a)*s + c;

/* cp1 pressure and temperature terms for s = 0 */
    a = (((1.7168e-8*t+2.0357e-6)*t-3.13885e-4)*t+1.45747e-2)*t-0.49592;
    b = (((2.2956e-11*t-4.0027e-9)*t+2.87533e-7)*t-1.08645e-5)*t+2.4931e-4;
    c = ((6.136e-13*t-6.5637e-11)*t+2.6380e-9)*t-5.422e-8;
    cp1 = ((c*p+b)*p+a)*p;

/* cp2 pressure and temperature terms for s > 0 */
    a = (((-2.9179e-10*t+2.5941e-8)*t+9.802e-7)*t-1.28315e-4)*t+4.9247e-3;
    b = (3.122e-8*t-1.517e-6)*t-1.2331e-4;
    a = (a+b*sr)*s;
    b = ((1.8448e-11*t-2.3905e-9)*t+1.17054e-7)*t-2.9558e-6;
    b = (b+9.971e-8*sr)*s;
    c = (3.513e-13*t-1.7682e-11)*t+5.540e-10;
    c = (c-1.4300e-12*t*sr)*s;
    cp2 = ((c*p+b)*p+a)*p;

/* specific heat return */
    return (cp0 + cp1 + cp2);
}

/*********************************************************/
/*********************************************************/
double hb_atg(double s,double t,double p)
/*********************************************************/
/* adiabatic temperature gradient deg c per decibar
.  ref: bryden,h.,1973,deep-sea res.,20,401-408
.  units:
.        pressure        p        decibars
.        temperature     t        deg celsius (ipts-68)
.        salinity        s        (ipss-78)
.        adiabatic       atg      deg. c/decibar
.  checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
.  t=40 deg c,p0=10000 decibars */
{
    double ds,atg;

    ds = s - 35.0;
    atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
        +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
        +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
        +(-4.2393e-8*t+1.8932e-6)*ds
        +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5;
    return atg;
}

/*********************************************************/
/*********************************************************/
double hb_theta(double s,double t0,double p0,double pr)
/*********************************************************/
/* to compute local potential temperature at pr
.  using bryden 1973 polynomial for adiabatic lapse rate
.  and runge-kutta 4-th order integration algorithm.
.  ref: bryden,h.,1973,deep-sea res.,20,401-408
.  fofonoff,n.,1977,deep-sea res.,24,489-491
. 
.        p0      (in situ) pressure     [db]
.        t0      (in situ) temperature  [deg C]
.        s        salinity              [ipss-78]
.        pr       reference pressure    [db]
.   returns:
.        theta    potential temperature [deg C]
. 
.  checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
.  p0=10000 decibars,pr=0 decibars */
{
    double p,t,h,xk,q;
    p = p0;
    t = t0;

    h = pr - p;
    xk = h * hb_atg(s,t,p);
    t = t + 0.5*xk;
    q = xk;
    p = p + 0.5*h;
    xk = h * hb_atg(s,t,p);
    t = t + 0.29289322*(xk-q);
    q = 0.58578644*xk + 0.121320344*q;
    xk = h * hb_atg(s,t,p);
    t = t + 1.707106781*(xk-q);
    q = 3.414213562*xk - 4.121320344*q;
    p = p + 0.5*h;
    xk = h * hb_atg(s,t,p);
    return (t + (xk-2.0*q)/6.0);
}
/*********************************************************/
/*********************************************************/
double hb_svel(double s,double t,double p0)
/*********************************************************/

/* sound velocity in seawater:  chen and millero 1977,jasa,62,1129-1135
. 
.        pressure        p0       decibars
.        temperature     t        deg celsius (ipts-68)
.        salinity        s        (ipss-78)
.        sound speed     svel     meters/second
. 
.  checkvalue: svel=1731.995 m/s, s=40 (ipss-78),t=40 deg c,p=10000 dbar */
{
    double p,sr;
    double a,a0,a1,a2,a3,b,b0,b1,c,c0,c1,c2,c3,d;

/*   scale pressure to bars */
    p=p0/10.;

    sr = sqrt(fabs(s));

/* s**2 term */
    d = 1.727e-3 - 7.9836e-6*p;

/* s**3/2 term */
    b1 = 7.3637e-5 +1.7945e-7*t;
    b0 = -1.922e-2 -4.42e-5*t;
    b = b0 + b1*p;

/* s**1 term */
    a3 = (-3.389e-13*t+6.649e-12)*t+1.100e-10;
    a2 = ((7.988e-12*t-1.6002e-10)*t+9.1041e-9)*t-3.9064e-7;
    a1 = (((-2.0122e-10*t+1.0507e-8)*t-6.4885e-8)*t-1.2580e-5)*t
          +9.4742e-5;
    a0 = (((-3.21e-8*t+2.006e-6)*t+7.164e-5)*t-1.262e-2)*t
          +1.389;
    a = ((a3*p+a2)*p+a1)*p+a0;

/* s**0 term */
    c3 = (-2.3643e-12*t+3.8504e-10)*t-9.7729e-9;
    c2 = (((1.0405e-12*t-2.5335e-10)*t+2.5974e-8)*t-1.7107e-6)*t
          +3.1260e-5;
    c1 = (((-6.1185e-10*t+1.3621e-7)*t-8.1788e-6)*t+6.8982e-4)*t
          +0.153563;
    c0 = ((((3.1464e-9*t-1.47800e-6)*t+3.3420e-4)*t-5.80852e-2)*t
          +5.03711)*t+1402.388;
    c = ((c3*p+c2)*p+c1)*p+c0;

/* sound speed return */
    return (c + (a+b*sr+d*s)*s);
}

/*********************************************************/
/*********************************************************/
double hb_grady(double *y,double *p,int nobs,double *pav_ptr,double *ybar_ptr)
/*********************************************************/
/*  Returns least squares slope 'grady' of y versus p.
.   The gradient is representive of the interval centered at pav;
.   ybar is the arithmetic average of y values over the entire interval.  */
{
    double grady,a0,cxx,cx,cxy,cy;
    int k;

    grady = 0.0;
    a0 = 0.0;
    cxx = 0.0;
    cx = 0.0;
    cxy = 0.0;
    cy = 0.0;
    *ybar_ptr = y[0];
    *pav_ptr  = p[0];
    
    if(nobs > 1) {
        for (k=0; k < nobs; k++) {
            cx =cx+p[k];
        }

        *pav_ptr = cx / nobs;

        for (k=0; k < nobs; k++) {
            cxy = cxy + y[k]*(p[k] - *pav_ptr);
            cy =cy+y[k];
            cxx=cxx+(pow(p[k]- *pav_ptr,2));
        }

        if(cxx == 0.0) {
            return grady;
        }

        a0 = cxy / cxx;
        *ybar_ptr = cy / nobs;
    }
    grady = a0;
    return grady;
}
/*********************************************************/
/*********************************************************/
double hb_bvfrq(double *s,double *t,double *p,int n,double *pav_ptr,double *e_ptr)
/*********************************************************/

/*  brunt-vaisala frequency  (uses eos80)
.   after formulation of breck owen's & n.p. fofonoff
. 
.        p      pressure     [dbars]
.        t      temperature  [deg C]
.        s      salinity     [ipss-78]
.        n      # of pts
.        pav    mean press over entire interval  [dbars] (returned)
.        e      n**2         [radians/sec]   (returned)
.        bvfrq  bouyancy freq  [cph]  (returned)
. ...............................................................
.  checkvalue: bvfrq=14.57836 cph e=6.4739928e-4 rad/sec.
.             s(1)=35.0, t(1)=5.0, p(1)=1000.0
.             s(2)=35.0, t(2)=4.0, p(2)=1002.0
.   >>> note result centered at pav=1001.0 dbars <<<
. ...............................................................   */
{
    double gr,gr2,tgrd,v350p,dvdp,vbar, bvf;
    double *data, mult; 
    double tav, sig, theta;
    int k;
    
    data = (double *) calloc(n, sizeof(double));

/* convert gravity to m/sec */
      gr = 9.80655;
      gr2 = gr * gr *1.e-4;

/* db to pascal conversion = 10^-4
.  get center pressure for interval */

    tgrd = hb_grady(t,p,n,pav_ptr,&tav);
    for (k=0; k < n; k++) {
	theta = hb_theta(s[k],t[k],p[k],*pav_ptr);
      	data[k] = hb_svan(s[k],theta,*pav_ptr,&sig) * 1.0e-8;
    }

/*  get v(35,0,pav) */
    v350p = (1./(sig + 1000.))- data[n-1];

/*  compute potential density anomaly gradient */
    dvdp= hb_grady(data,p,n,pav_ptr,&vbar);
    vbar += v350p;

    *e_ptr = -gr2 * dvdp/ (vbar * vbar);
    mult = *e_ptr < 0 ? -1.0 : 1.0;
    
    bvf = RPS2CPH * mult * sqrt(fabs(*e_ptr));
    return(bvf);
}
/*********************************************************/
/*********************************************************/
double hb_linterp(double xval, double *x, double *y, int npts)
/*********************************************************/

/*   linearly interpolates to find the position of xval in array x, and returns
.    the corresponding value in the array, y.  If xval does not appear in array
.    x, the value of "flag" is returned.  This routine assumes that both x and y
.    arrays are continuous (no missing values)*/
{
    int k;
    double v1,v2;
    
    if (npts <= 0)
       return (FLAG);
       
    if (npts == 1) {
       if (xval == x[0] )
            return(y[0]);
	    
	    
	return (FLAG);
    }
    
    for (k=0; k<npts-1; k++) {
        v1 = xval - x[k];
        v2 = xval - x[k+1];
        
        if (v1 == 0)
             return y[k];
        if (v2 == 0)
             return y[k+1];

        if ((v1 < 0.) && (v2 < 0.))
            continue;
        else if ((v1 > 0.) && (v2 > 0.))
            continue;
        else if ((v1 == v2) && (v1 != 0.))
            continue;
        else if ((v1 == 0.) && (v2 ==  0.))
            return y[k];
        else
            return (y[k] + (y[k+1] - y[k]) * v1 / (x[k+1] - x[k]));
    }

/*  if execution falls through to this point, xval was not found in x array;
.      return a flag value ....  */
    return FLAG;
}

/*********************************************************/
/*********************************************************/
double hb_bvgrdts(double *p0_ptr,double *s,double *t,double *p,int n,double *e_ptr)
/*********************************************************/
/* arguments:
.p0_ptr:    pressure level at which bvgrdts is evaluated (if p0 < 0,
.           (p(0)+p(n))/2 is used and returned in this argument ).
.     p:    observed pressures  ( decibars)
.     t:    observed temperatures   (deg celsius (ipts-68))
.     s:    observed salinities  (ipss-78)
.     n:    number of observations   (size of p,t,s arrays)
. e_ptr:    stability parameter  (radians/second^2) returned
.    bvgrdts:    bouyancy frequency (cph) returned
. ....................................
.  checkvalue: bvfrq=14.5xx cph e=6.47xxxe-4 rad/sec^2.
.             s(0)=35.0, t(0)=5.0, p(0)=1000.0
.             s(1)=35.0, t(1)=4.0, p(1)=1002.0
.      (note result centered at pav=1001.0 dbars)
. ....................................
.   Brunt-Vaisala frequency calculation with eos80 gradients of temp and salt.
.   The derivative quantities are computed in function eos8d() and returned in
.   the array drv.
. 
.   see The Oceans - Sverdrup, Johnson & Fleming p. 417-418
.   also Hesselberg & Sverdrup (1915)  ref. page 430 The Oceans
. ......................................................................... */
{
    double **drv;
    double bvgrdts;
    double gp2, tgrd, sgrd, pav, tav, sav;
    double z, gt, gs, v, v2;
    int i; 

/* allocate memory for 2-D array */
    
    drv = (double **)calloc(3, sizeof(double *));
    for (i=0; i < 3; ++i) {
       drv[i] = (double *) calloc(8, sizeof(double));
    }

/* 10^-4:pressure db to pascals (newton/m^2) */
    gp2 = GP * GP * 1.e-4;

/* compute gradients and values of t & s at pav ... */
    tgrd = hb_grady(t,p,n,&pav,&tav);
    sgrd = hb_grady(s,p,n,&pav,&sav);

    if (*p0_ptr >= 0) {
        pav = *p0_ptr;
        tav = hb_linterp(*p0_ptr, p, t, n);
        sav = hb_linterp(*p0_ptr, p, s, n);
    }

    z = hb_eos80d(sav, tav, pav, drv);

/* derivatives of specific vol. with respect to temp and salt... */
    gt = -drv[1][0];
    gs = -drv[2][7];

    v = (1 / (1000.+ drv[1][3]));
    v2 = v * v;
    *e_ptr = (gp2/v2*(gt*(tgrd - hb_atg(sav,tav,pav))+ gs * sgrd));
    bvgrdts = RPS2CPH * sqrt(fabs(*e_ptr));
    
    if (*e_ptr < 0.)  /* apply sign of *e_ptr to final value */
       bvgrdts = -bvgrdts;
       
    *p0_ptr = pav;

/* free up memory and return */
    
    for (i=0; i < 3; ++i) {
       free((void *)drv[i]);
    }
    free((void *) drv);
    
    return bvgrdts;
}

/*********************************************************/
/*********************************************************/
double hb_p80(double dpth,double xlat)
/*********************************************************/
/* pressure from depth from saunder's formula with eos80.
.  reference: saunders,peter m., practical conversion of pressure
.             to depth., j.p.o. , april 1981.
.  r millard
.  march 9, 1983
.  check value: p80 = 7500.004 dbars
.               for lat = 30 deg, depth = 7321.45 meters */
{
    double plat,d,c1;

    plat = fabs(xlat*PI/180.);
    d = sin(plat);
    c1 = 5.92e-3 + pow(d,2) * 5.25e-3;
    return ((1-c1)-sqrt((pow((1-c1),2))-(8.84e-6*dpth)))/4.42e-6;
}

/*********************************************************/
/*********************************************************/
double hb_gravity(double xlat)
/*********************************************************/
/*  Acceleration due to gravity in cm/sec^2 as a function of latitude */
{
    double plat, g;

    plat = fabs(xlat*PI/180.);
    g = 978.0318 *(1.0 + 5.3024e-3 * pow(sin(plat),2)
              - 5.9e-6 * pow((sin(2.*plat)),2));
    return(g);	      
}

/*********************************************************/
/*********************************************************/
double hb_coriol(double xlat)
/*********************************************************/
/*  Coriolis parameter as a function of latitude in sec^-1 */
{
    double plat;

    plat = fabs(xlat * PI / 180.);
    return (14.5842e-5 * sin(plat));
}
