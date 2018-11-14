/*  gamma_subs.c

........................................................................
                          *******  HydroBase 2 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
                             Mar 2001
...................................................
   Functions for computing neutral density (gamma_n).  Adapted from
   Fortran code by David Jackett and Trevor McDougall.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hydrobase.h"
#include "hb_gamma.h"
#include "netcdf.h"

/**************************************************/
void compute_gamma_n(struct GAMMA_NC *ginfo, int nobs, double *gamma, double *p, double *t, double *s, double dlon, double dlat)

   /* Label each observation level with neutral density (gamma-n) but
      no error bars.  Note that gamma_init_nc() must have been called
      previously to initialize *ginfo. */

{
   int    i, error;
   double *dg_lo, *dg_hi;
   
   /* Check that gamma.nc has been opened and read */
   
   if (ginfo ->id_s == ginfo->id_t) {
     fprintf(stderr,"\nFATAL ERROR:  gamma info not initialized from gamma.nc\n");
     exit(1);
   }
   
   dg_lo = (double *) calloc((size_t) nobs, sizeof(double));
   dg_hi = (double *) calloc((size_t) nobs, sizeof(double));
   if (dg_hi == NULL) {
     fprintf(stderr,"\nFATAL ERROR:  Unable to allocate memory in compute_gamma_n()\n");
     exit(1);
   }
   error = gamma_n(ginfo, s, t, p, nobs, dlon, dlat, gamma, dg_lo, dg_hi);

   for (i = 0; i < nobs; ++i) {
      if (gamma[i] < 0 || error < 0)
         gamma[i] = (double) HB_MISSING;
   
   }   
   free((void *) dg_lo);
   free((void *) dg_hi);
   
   return;

} /* end compute_gamma_n() */

/**************************************************/
void compute_gamma_nerr(struct GAMMA_NC *ginfo, int nobs, double *gamma, double *gerror, double *p, double *t, double *s, double dlon, double dlat)

   /* Label each observation level with neutral density (gamma-n) and
      computes associated error bar. Note that gamma_init_nc() must have been
      called previously to initialize *ginfo.*/ 
{
   int    i, error;
   double *dg_hi;
   
   /* Check that gamma.nc has been opened and read */
   
   if (ginfo ->id_s == ginfo->id_t) {
     fprintf(stderr,"\nFATAL ERROR:  info from gamma.nc not initialized\n");
     exit(1);
   }
   
   dg_hi = (double *) calloc((size_t) nobs, sizeof(double));
   if (dg_hi == NULL) {
     fprintf(stderr,"\nFATAL ERROR:  Unable to allocate memory in compute_gamma_nerr()\n");
     exit(1);
   }
   error = gamma_n(ginfo, s, t, p, nobs, dlon, dlat, gamma, gerror, dg_hi);
   
   for (i = 0; i < nobs; ++i) {
      if (gamma[i] < 0 || error < 0) {
         gamma[i] = (double) HB_MISSING;
         gerror[i] = (double) HB_MISSING;
      }
   }
      
   free((void *) dg_hi);
   
   return;

} /* end compute_gamma_nerr() */

/**************************************************/
void depth_ns(double *s, double *t, double *p, int n, double s0, double t0, double p0, double *sns_addr, double *tns_addr, double *pns_addr)

  /* Finds the position at which the neutral surface through a
       specified bottle intersects a neighbouring cast.
       adapted from code by D.Jackett and T.McDougall:
         
	INPUT :		s(n)		array of cast salinities
			t(n)		array of cast in situ temperatures
			p(n)		array of cast pressures
			n		length of cast
			s0		the bottle salinity
			t0		the bottle in situ temperature
			p0		the bottle pressure

	returns :	*sns_addr	salinity of the neutral surface
					intersection with the cast
			*tns_addr	in situ temperature of the intersection
			*pns_addr	pressure of the intersection

	UNITS :		salinities	psu (IPSS-78)
			temperatures	degrees C (IPTS-68)
			pressures	db
  */
{
   int i, k, ncr, iter, niter, last;
   int success;
   double sigl, sigu;
   double pc0, tc0, sc0, ec0, ecz0;
   double pc_0, ec_0, ecz_0;
   double p1, p2, ez1, ez2, pc1;
   double eps, r; 
   double *e;
   
   e = (double *) calloc((size_t) n, sizeof(double));

/* compute sigma difference between bottle and each level of cast */

   for (k = 0; k < n; ++k) {
      sig_vals(s0, t0, p0, s[k], t[k], p[k], &sigl, &sigu);
      e[k] = sigu - sigl;
   }

/* special case for cast length == 1 */

   if ( n == 1) {
   
      *sns_addr = s[0];
      *tns_addr = t[0];
      *pns_addr = p[0];
      free((void *) e);
      return;
   }
         
/* find the bottle pairs containing a crossing */

   last = n-1;
   ncr = 0;
   for (k = 1; k <= last; ++k) { 
   
     i = k-1; 
     if (e[i] == 0.0) {   /* an exact crossing at k-1 */
         ++ncr;
	 *sns_addr = s[i];
	 *tns_addr = t[i];
	 *pns_addr = p[i];
     } /* end if */
     
     else {
     
       if ((e[k]*e[i]) < 0.0) {
          ++ncr;
	  
          /* some Newton-Raphson iterations to find the crossing */
	  
	  pc0 = p[i] - e[i] *(p[k]-p[i])/(e[k]-e[i]);
	  iter = 0;
	  success = 0;
	  
	  do {
	     ++iter;
	     stp_interp(&s[i], &t[i], &p[i], 2, &sc0, &tc0, &pc0);
             sig_vals(s0, t0, p0, sc0, tc0, pc0, &sigl, &sigu);
             ec0 = sigu - sigl;
		 
	     p1 = (p[i] + pc0) * 0.5;
             ez1 = (e[i] - ec0)/(pc0 - p[i]);
	     p2 = (pc0 + p[k]) * 0.5;
	     ez2 = (ec0 - e[k])/(p[k] - pc0);
	     r = (pc0 - p1) / (p2 - p1);
	     ecz_0 = ez1 + r *(ez2-ez1);
	      
             ecz0 = ecz_0;
	     if (iter > 1) {
		  ecz0 = -(ec0 - ec_0) / (pc0 - pc_0);
		  if (ecz0 == 0.0) 
		     ecz0 = ecz_0;
	     }
		
	     pc1 = pc0 + ec0 / ecz0;
		
     /*  strategy when the iteration jumps out of the interval */
     
             if (pc1 <= p[i] || pc1 >= p[k]) {
	        e_solve(s, t, p, e, n, k, s0, t0, p0, sns_addr, tns_addr, pns_addr, &niter);
                if (*pns_addr < p[i] || *pns_addr > p[k]) {
		    fprintf(stderr,"\nFATAL ERROR (1) in depth-ns()\n");
		    exit(1);
		}
		  
		success = 1;
		  
             } /* end if */
		
	     else {
  
      /*  otherwise, test the accuracy of the iterate ... */

     		eps = ABS(pc1-pc0);
		
		if ((ABS(ec0) <= 5.0E-5) && (eps <= 5.0E-3)) {
		   
		   *sns_addr = sc0;
		   *tns_addr = tc0;
		   *pns_addr = pc0;
		   success = 1;
		}
		else if (iter > 10) {
		   e_solve(s, t, p, e, n, k, s0, t0, p0, sns_addr, tns_addr, pns_addr, &niter);
		   success = 1;
		} 
		else {
		   pc_0 = pc0;
		   ec_0 = ec0;
		   pc0 = pc1;
		}

             } /* end else */
		
	  } while (!success);
	  
       }  /* end if */
     
     } /* end else */
     
   } /* end for k */

  /* check last bottle for exact crossing */  
   
  if (e[last] == 0.0) {
	 ++ncr;
	*sns_addr = s[last];
	*tns_addr = t[last];
	*pns_addr = p[last];
   }
   
   free((void *) e);
    
   if (!ncr) {          /*  no crossings */
	*sns_addr = -99.0;
	*tns_addr = -99.0;
	*pns_addr = -99.0;
        return;
   }
    
   if (ncr > 1) {     /*  multiple crossings  */
	*sns_addr = -99.2;
	*tns_addr = -99.2;
	*pns_addr = -99.2;
        return;
   }
    
   return;
   
} /* end depth_ns() */

/**************************************************/

void depth_scv(double *s, double *t, double *p, int n, double s0, double t0, double p0, double *sscv, double *tscv, double *pscv, int *nscv_addr)

/*	Find the position at which the scv surface through a
	specified bottle intersects a neighbouring cast
        adapted from code by D.Jackett and T.McDougall:

	INPUT :		s[n]		array of cast salinities
			t[n]		array of cast in situ temperatures
			p[n]		array of cast pressures
			n		length of cast
			s0		the bottle salinity
			t0		the bottle in situ temperature
			p0		the bottle pressure
			sscv[MAX_SCV]   salinity array already allocated
			tscv[MAX_SCV]	temperature  "           "
			pscv[MAX_SCV]	pressure   "             "

	returns :	sscv[MAX_SCV]	salinities of the scv surface
					intersections with the cast
			tscv[MAX_SCV]	temperatures of the intersections
			pscv[MAX_SCV]	pressures of the intersections
			*nscv_addr	number of intersections

	UNITS :		salinities	psu (IPSS-78)
			temperatures	degrees C (IPTS-68)
			pressures	db

 */
{
   int i, k, iter, niter, last;
   int success;
   double sigl, sigu;
   double pc0, tc0, sc0, ec0, ecz0;
   double pc_0, ec_0, ecz_0;
   double p1, p2, ez1, ez2, pc1;
   double eps, r; 
   double dummy, tref;
   double *e;
   
   
   e = (double *) calloc((size_t)n, sizeof(double));

/* compute sigma difference between bottle and each level of cast */

   for (k = 0; k < n; ++k) {
      tref = hb_theta(s0,t0, p0, p[k]);
      dummy = hb_svan(s0, tref, p[k], &sigl);
      dummy = hb_svan(s[k], t[k], p[k], &sigu);
      e[k] = sigu - sigl;
   }
   
/* find the bottle pairs containing a crossing */
      
   last = n-1;
   *nscv_addr = 0;
   
   for (k = 1; k <= last; ++k) { 
     i = k-1; 
     if (e[i] == 0.0) {   /* an exact crossing at k-1 */
	 sscv[*nscv_addr] = s[i];
	 tscv[*nscv_addr] = t[i];
	 pscv[*nscv_addr] = p[i];
	 ++(*nscv_addr);
     } /* end if */
     
     else {
     
       if ((e[k]*e[i]) < 0.0) {
	  
          /* some Newton-Raphson iterations to find the crossing */

	  pc0 = p[i] - e[i] *(p[k]-p[i])/(e[k]-e[i]);
	  iter = 0;
	  success = 0;
	  
	  do {
	     ++iter;
	     stp_interp(&s[i], &t[i], &p[i], 2, &sc0, &tc0, &pc0);
             tref = hb_theta(s0,t0, p0, pc0);
             dummy = hb_svan(s0, tref, pc0, &sigl);
             dummy = hb_svan(sc0, tc0, pc0, &sigu);
             ec0 = sigu - sigl;
	  
	     p1 = (p[i] + pc0) * 0.5;
             ez1 = (e[i] - ec0)/(pc0 - p[i]);
	     p2 = (pc0 + p[k]) * 0.5;
	     ez2 = (ec0 - e[k])/(p[k] - pc0);
	     r = (pc0 - p1) / (p2 - p1);
	     ecz_0 = ez1 + r *(ez2-ez1);
	      
             ecz0 = ecz_0;
	     if (iter > 1) {
		  ecz0 = -(ec0 - ec_0) / (pc0 - pc_0);
		  if (ecz0 == 0.0) 
		     ecz0 = ecz_0;
	     }
		
	     pc1 = pc0 + ec0 / ecz0;
	     
     /*  strategy when the iteration jumps out of the interval */
     
             if (pc1 <= p[i] || pc1 >= p[k]) {
	       scv_solve(s, t, p, e, n, k, s0, t0, p0, &sscv[*nscv_addr], &tscv[*nscv_addr], &pscv[*nscv_addr], &niter);
	       
	       if (pscv[*nscv_addr] < p[i] || pscv[*nscv_addr] > p[k]) {
		    fprintf(stderr,"\nFATAL ERROR (1) in depth-scv()");
		    fprintf(stderr,"\nscv_solve() returned bad pressure. \n");
		    exit(1);
	       }
		success = 1;
	        ++(*nscv_addr);
	       
             } /* end if */
	     
	     else {
  
      /*  otherwise, test the accuracy of the iterate ... */

     		eps = ABS(pc1-pc0);
		
		if ((ABS(ec0) <= 5.0E-5) && (eps <= 5.0E-3)) {
		   
		   sscv[*nscv_addr] = sc0;
		   tscv[*nscv_addr] = tc0;
		   pscv[*nscv_addr] = pc0;
		   success = 1;
	           ++(*nscv_addr);
		}
		else if (iter > 10) {
		   scv_solve(s, t, p, e, n, k, s0, t0, p0, &sscv[*nscv_addr], &tscv[*nscv_addr], &pscv[*nscv_addr], &niter);
		   success = 1;
	           ++(*nscv_addr);
		} 
		else {
		   pc_0 = pc0;
		   ec_0 = ec0;
		   pc0 = pc1;
		}

             } /* end else */
	     
	  } while (!success);
	  
       }  /* end if */
     
     } /* end else */
     
     
     /* check array size to prevent overflow */
     
     if ( *nscv_addr >= MAX_SCV) {
       fprintf(stderr,"\n FATAL ERROR::depth_scv():  SCV arays not big enough.");
       fprintf(stderr,"\n Increase size of MAX_SCV in hb_gamma.h\n");
       exit(2);
     }
     
   } /* end for k */

  /* check last bottle for crossing */ 
   
   if (e[last] == 0.0) {
	sscv[*nscv_addr] = s[last];
	tscv[*nscv_addr] = t[last];
	pscv[*nscv_addr] = p[last];
	 ++(*nscv_addr);
   }
   
   free((void *)e);
   
   
   if (*nscv_addr == 0) {
      sscv[0] = -99.0;
      tscv[0] = -99.0;
      pscv[0] = -99.0;
   }
   
   return;

}  /* end depth_scv() */


/***********************************************************/

void derthe(double s, double t, double p0, double *dthedt, double *dtheds, double *dthedp)

  /* Uses the Bryden (1973) polynomial for potential temperature as
     a function of (s,t,p) to obtain the partial derivatives of theta
     WRT t,s,p.  Pressure is in dbars.  Adapted from code by D.Jackett 
     and T.McDougall.
  */
{

   double A0 = -0.36504E-4;
   double A1 = -0.83198E-5;
   double A2 = +0.54065E-7;
   double A3 = -0.40274E-9;
   double B0 = -0.17439E-5;
   double B1 = +0.29778E-7;
   double D0 = +0.41057E-10;
   double C0 = -0.89309E-8;
   double C1 = +0.31628E-9;
   double C2 = -0.21987E-11;
   double E0 = +0.16056E-12;
   double E1 = -0.50484E-14;
   
   double ds, p, pp, ppp;
   double part, tt, ttt;

   ds = s - 35.0;
   p = p0;
   pp = p * p;
   ppp = pp * p;
   tt = t * t;
   ttt = tt * t;
   part = 1.0 + p * (A1 + 2.0 * A2 * t + 3.0 * A3 * tt + ds * B1);
   *dthedt = part + pp *(C1 + 2.0*C2*t) + ppp * E1;
   *dtheds = p * (B0 + B1 * t) + pp * D0;
   part = A0 + A1 * t + A2 * tt + A3 * ttt + ds * (B0 + B1*t);
   *dthedp = part + 2.0 * p * (ds *D0 + C0 + C1*t + C2*tt) + 3.0*pp * (E0+ E1*t);

   return;
} /*end derthe() */

/***********************************************************/

void e_solve(double *s, double *t, double *p, double *e, int n, int k, double s0, double t0, double p0, double *sns_addr, double *tns_addr, double *pns_addr, int *iter)

/*  adapted from code by D.Jackett and T.McDougall:

 	DESCRIPTION :	Find the zero of the e function using a 
 			bisection method
 
 	PRECISION :	Double
 
 	INPUT :		s(n)		array of cast salinities
 			t(n)		array of cast in situ temperatures
 			p(n)		array of cast pressures
 			e(n)		array of cast e values
 			n		length of cast
 			k		interval (k-1,k) contains the zero
 			s0		the bottle salinity
 			t0		the bottle in situ temperature
 			p0		the bottle pressure
 
 	Returns :	*sns_addr	salinity of the neutral surface
 					intersection with the cast
 			*tns_addr	in situ temperature of the intersection
 			*pns_addr	pressure of the intersection
 
 
 	UNITS :		salinities	psu (IPSS-78)
 			temperatures	degrees C (IPTS-68)
 			pressures	db
 
 
 	AUTHOR :	David Jackett
 
 	CREATED :	June 1993
 
 	REVISION :	1.1		30/6/93
 

*/

{
   int i;
   double pl, pu, pm;
   double el, eu, em;
   double sm, tm;
   double sigl, sigu; 
   
   i = k-1;
   pl = p[i];
   el = e[i];
   pu = p[k];
   eu = e[k];
   
   *iter = 0;
   
   do {
     ++(*iter);
     pm = (pl + pu) * 0.5;
     
     stp_interp(&s[i], &t[i], &p[i], 2, &sm, &tm, &pm);
     sig_vals(s0, t0, p0, sm, tm, pm, &sigl, &sigu);
     em = sigu - sigl;
     
     if (el * em < 0.0) {
        pu = pm;
	eu = em;
     }
     else if (em *eu < 0.0) {
        pl = pm;
	el = em;
     }
     else {
        if (em == 0.0) {
	  *sns_addr = sm;
	  *tns_addr = tm;
	  *pns_addr = pm;
	  return;
	}
     }
     
     if ((ABS(em) <= 5.0E-5) && (ABS(pu-pl) <= 5.0E-3)) {
        *sns_addr = sm;
	*tns_addr = tm;
	*pns_addr = pm;
	return;
     }
     
   } while (*iter < 20);


   fprintf(stderr, "\nWARNING from e_solve()");   
   fprintf(stderr, "\niteration #%d  em: %lf  dp: %lf - %lf = %lf", *iter, (double) ABS(em), pl, pu, (double)ABS(pu-pl));   
   *sns_addr = -99.0;
   *tns_addr = -99.0;
   *pns_addr = -99.0;
   return;

}  /* end e_solve() */

/**********************************************************/

void eosall(double s, double t, double p0, double *theta_ptr, double *sigma_ptr, double *alpha_ptr, double *beta_ptr, double *gamma_ptr, double *soundv_ptr)

/*  Adapted from fortran code: 

C  WRITTEN 8 JULY 1985 BY TREVOR J. McDOUGALL
C EOSALL STANDS FOR "EQUATION OF STATE ALL"
C THIS SUBROUTINE USES THE FOLLOWING FUNCTIONS WRITEN BY BOB MILLARD,
C  - THETA(S,T,P0,PR) ; SVAN(S,T,P0,SIGTHE) ; EOS8D(S,T,P0,DRV)
C THE NEW EXPANSION COEFFICIENT , *alpha_ptr  , DUE TO HEAT , AND THE 
C CORRESPONDING SALINE CONTRACTION COEFFICIENT ARE DEFINED IN 
C TERMS OF THE TWO CONSERVATIVE PROPERTIES OF SEA WATER, 
C NAMELY POTENTIAL TEMPERATURE (REFERED TO ANY REFERENCE LEVEL)
C AND SALINITY. THESE COEFFICIENTS ARE DEFINED IN GILL(1982)
C AND HE LABELS THEM WITH DASHES, SEE HIS SECTION 3.7.4
 
  Check values of these new "expansion coefficients"
  at s=40.0,t=40.0,theta =36.8907, p0=10000.0 dbars:
  
 *alpha_ptr = 4395.6E-7 ; (alpha_old = 4181.1E-7)
 *beta_ptr = 6646.9E-7 ;  (beta_old = 6653.1E-7)
 *gamma_ptr = 31.4E-7   ; (gamma_old = 32.7E-7)
 *soundv_ptr = 1734.8 m/s.

 in-situ density is (drv[0][2] + 1000.0)

 */
{
   int i;
   double **drv;
   double dthedt, dtheds, dthedp;
   double pref, dummy;
   double alpha_old, beta_old, gamma_old;

    drv = (double **)calloc(3, sizeof(double *));
    for (i=0; i < 3; ++i) {
       drv[i] = (double *) calloc(8, sizeof(double));
    }
      
/* Reference pressure (pref) is kept general but will be equal to
   zero for all perceived applications of neutral surfaces. */

   pref = 0.0;
   *theta_ptr = hb_theta(s, t, p0, pref);
   dummy = hb_eos80d(s, t, p0, drv);

   alpha_old = - drv[1][2] / (drv[0][2] + 1000.0);
   beta_old = drv[0][7] / (drv[0][2] + 1000.0);
   gamma_old = drv[0][6] / (drv[0][2] + 1000.0);
   
/* calculate specific volume anomaly and sigma-theta */

   dummy = hb_svan(s, *theta_ptr, pref, sigma_ptr);
   derthe(s, t, p0, &dthedt, &dtheds, &dthedp);
   
   *alpha_ptr = alpha_old / dthedt;
   *beta_ptr = beta_old + *alpha_ptr * dtheds;
   *gamma_ptr = gamma_old + *alpha_ptr * dthedp;
   *soundv_ptr= sqrt((double) (ABS(1.0E+4 / (*gamma_ptr *(drv[0][2] + 1000.0)))));

/* free up memory */
    
   for (i=0; i < 3; ++i) 
       free((void *)drv[i]);
       
   free((void *) drv);
    
   return;
}  /* end eosall() */
/**********************************************************/

void gamma_errors(double *s, double *t, double *p, double *gamma, double *a, int n, double along, double alat, double s0, double t0, double p0, double sns, double tns, double pns, int kns, double gamma_ns, double *pth_error_ptr, double *scv_l_error_ptr, double *scv_h_error_ptr)

/*  Adapted from fortran code by D.Jackett and T.McDougall:

 	DESCRIPTION :	Find the p-theta and the scv errors associated 
 			with the basic neutral surface calculation
 
 	PRECISION :	Double
 
 	INPUT :		s(n)		array of Levitus cast salinities
 			t(n)		array of cast in situ temperatures
 			p(n)		array of cast pressures
 			gamma(n)	array of cast neutral densities
 			a(n)		array of cast quadratic coefficients
 			n		length of cast
 			along		longitude of Levitus cast
 			alat		latitude of Levitus cast
 			s0		bottle salinity
 			t0		bottle temperature
 			p0		bottle pressure
 			sns		salinity of neutral surface on cast
 			tns		temperature of neutral surface on cast
 			pns		pressure of neutral surface on cast
 			kns		index of neutral surface on cast
 			gamma_ns	gamma value of neutral surface on cast
 
 	returns :	*pth_error_ptr	 p-theta gamma error bar
 			*scv_l_error_ptr lower scv gamma error bar
 			*scv_h_error_ptr upper scv gamma error bar
 
 	UNITS :		salinity	psu (IPSS-78)
 			temperature	degrees C (IPTS-68)
 			pressure	db
 			gamma		kg m-3
 
 
 	AUTHOR :	David Jackett
 
 	CREATED :	March 1995
 
 	REVISION :	1.1		9/3/95

*/

{
  int kns1, kscv, kscv1, nscv;
  double sscv_m[MAX_SCV], tscv_m[MAX_SCV], pscv_m[MAX_SCV];
  double pscv, pscv_mid, gamma_scv;
  double pr0, Tb, gamma_limit, test_limit;
  double th0, thns, dummy, rho_ns, sig_ns;
  double  dp, dth, sig_l, sig_h, b, test, drldp;
  
  pr0 = 0.0;
  Tb = 2.7e-8;
  gamma_limit = 26.845;
  test_limit = 0.1;
  
  /* p - theta error */
  
  th0 = hb_theta(s0, t0, p0, pr0);
  thns = hb_theta(sns, tns, pns, pr0);
  dummy = hb_svan(sns, tns, pns, &sig_ns);
  rho_ns = 1000. + sig_ns;

  kns1 = kns + 1;  
  sig_vals(s[kns], t[kns], p[kns], s[kns1], t[kns1], p[kns1], &sig_l, &sig_h);
  
  b = (gamma[kns1] - gamma[kns]) / (sig_h - sig_l);
  
  dp = pns - p0;
  dth = thns - th0;
  *pth_error_ptr = rho_ns * b * Tb * ABS(dp*dth) / 6.0;
  
  /* scv error */
  
  *scv_l_error_ptr = 0.0;
  *scv_h_error_ptr = 0.0;
  
  if (alat <= -60. || gamma[0] >= gamma_limit) {
    drldp = (sig_h - sig_l) / (rho_ns * (p[kns1] - p[kns]));
    test = Tb * dth / drldp;
  
  
/*  approximation  */

     if (ABS(test) <= test_limit) {
  
        if (dp * dth >= 0.) 
	   *scv_h_error_ptr = (3 * *pth_error_ptr) / (1.0 - test);
        else
           *scv_l_error_ptr = (3 * *pth_error_ptr) / (1.0 - test);
     
     }
     else {      /* explicit scv solution, when necessary */
  
        depth_scv(s, t, p ,n, s0, t0, p0, sscv_m, tscv_m, pscv_m, &nscv);

        if (nscv > 0) {
     
	   pscv = pscv_m[0];
	   if (nscv > 1) {
	     pscv_mid = pscv_m[nscv/2];
	  
	     pscv = pscv_m[0];
	     if (p0 > pscv_mid)
	       pscv = pscv_m[nscv-1];
	    
	   } /* end if nscv > 1 */
	
	   kscv = indx(p, n, pscv);
	   kscv1 = kscv+1;
	   gamma_scv = gamma_qdr(p[kscv], gamma[kscv], a[kscv], p[kscv1], gamma[kscv1], pscv);
	
	   if (pscv <= pns) 
	      *scv_l_error_ptr = gamma_ns - gamma_scv;
	   else
	      *scv_h_error_ptr = gamma_scv - gamma_ns;
	
        } /* end if nscv > 0 */
     } /* end if ABS(test) */
  }  /* end if alat <= -60. || gamma[0] >= gamma_limit */
  
  /*  check for negative gamma errors */
  
  if (*pth_error_ptr < 0. || *scv_l_error_ptr < 0. || *scv_h_error_ptr < 0.) {
     fprintf(stderr,"\nFATAL ERROR in gamma_errors():  negative scv \n");
     exit(2);
  }

  return;

} /* end gamma_errors() */

/*******************************************************/


int gamma_n(struct GAMMA_NC *gfile, double *s, double *t, double *p, int n, double along, double alat, double *gamma, double *dg_lo, double *dg_hi)

/*

   NOTE:  in translating from the fortran code to C, the ordering of the
   multi-dimensioned array has been inverted to reflect the way
   the data are stored in the gamma.nc file.  The function read_nc() returns
   4 neighboring casts of {s, t, gamma, and a} each dimensioned:
            [lat][lon][pressure]
	    
   Thus the index variables in this code are reversed [j0][i0][nz] relative to the
   original fortran (nz,i0,j0) to reflect this.
   

 	DESCRIPTION :	Label a cast of hydrographic data at a specified 
 			location with neutral density
 
 
 	INPUT :	        gfile           ptr to a struct GAMMA_NC 
	        	s[n]		array of cast salinities
 			t[n]		array of cast in situ temperatures
 			p[n]		array of cast pressures
 			n		length of cast (n=1 for single bottle)
 			along	        longitude of cast (0-360)
 			alat		latitude of cast (-80,64)
 
 	returns :       1 for success, -1 for error
	        	*gamma		array of cast gamma values
 			*dg_lo		array of gamma lower error estimates
 			*dg_hi		array of gamma upper error estimates
 
 			NOTE:		-99.0 denotes algorithm failed
 					-99.1 denotes input data is outside
 					      the valid range of the present
 					      equation of state
 
 	UNITS :		salinity	psu (IPSS-78)
 			temperature	degrees C (IPTS-68)
 			pressure	db
 			gamma		kg m-3
 
 
*/
{
   int nz, ndx, ndy, nij;
   int i, j, k, ij, j0, i0, i_min, j_min;
   int ioce, ialtered;
   int itest, kns;
   int **n0, **iocean0;
   double along0[2], alat0[2];
   double dx, dy, rx, ry, rw;
   double dist2_min, dist2, dgw_max;
   double dgamma_0, dgamma_1, dgamma_2_l, dgamma_2_h, dgamma_3;
   double gw, gn, g1_err, g2_l_err, g2_h_err;
   double pns, sns, tns;
   double pr0, thk, wsum, wt;
   double *p0;
   double gwij[4], wtij[4];
   double ***s0, ***t0, ***gamma0, ***a0;


/* initialize */
   
   nz = NZ;
   ndx = NDX;
   ndy = NDY;
   
   pr0 = 0.0;
   p0 = gfile->p0;
   dgamma_0 = 0.0005;
   dgw_max = 0.3;


/* check for errors */

   ialtered = 0;
   
   if (along < 0.) {
      along += 360.0;
      ialtered = 1;
   }

   if (along == 360.0) {
      along = 0.0;
      ialtered = 2;
   }
   
   if (along < 0. || along > 360. || alat < -80. || alat > 64.)
      return (-1);
   
/*  allocate memory */

   n0 = (int **) calloc(2, sizeof(int *));   
   iocean0 = (int **) calloc(2, sizeof(int *));   
   s0 = (double ***)calloc(2, sizeof(double **));
   t0 = (double ***)calloc(2, sizeof(double **));
   a0 = (double ***)calloc(2, sizeof(double **));
   gamma0 = (double ***)calloc(2, sizeof(double **));
   
   for (i = 0; i < 2; ++i) {
      n0[i] = (int *) calloc(2, sizeof(int));
      iocean0[i] = (int *) calloc(2, sizeof(int));
      s0[i] = (double **)calloc(2, sizeof(double *));
      t0[i] = (double **)calloc(2, sizeof(double *));
      a0[i] = (double **)calloc(2, sizeof(double *));
      gamma0[i] = (double **)calloc(2, sizeof(double *));
      
      for (j = 0; j < 2; ++j) {
        s0[i][j] = (double *) calloc(NZ, sizeof(double));
        t0[i][j] = (double *) calloc(NZ, sizeof(double));
        a0[i][j] = (double *) calloc(NZ, sizeof(double));
        gamma0[i][j] = (double *) calloc(NZ, sizeof(double));
      }
   }
   

   for (k = 0; k < n; ++k) {
   
      gamma[k] = 0.0;
      dg_lo[k] = 0.0;
      dg_hi[k] = 0.0;
   
      if (s[k] < 0. || s[k] > 42. || t[k] < -2.5 || t[k] > 40. 
       || p[k] < 0. || p[k] > 10000.) {
      
         gamma[k] = -99.1;
	 dg_lo[k] = -99.1;
	 dg_hi[k] = -99.1;
      }
   }
   
/*   read records from the netCDF data file  */

   read_nc(gfile, along, alat, s0, t0, p0, gamma0, a0, n0, along0, alat0, iocean0);
   
/*   find the closest cast  */

    dist2_min = 1.e10;
    
    i_min = -1;
    
   for (j0 = 0; j0 < 2; ++j0) {
     for (i0 = 0; i0 < 2; ++i0) {
     
        if (n0[j0][i0] != 0) {
	  dist2 = (along0[i0] - along) * (along0[i0] - along) + (alat0[j0] - alat) * (alat0[j0] - alat);
	  
	  if (dist2 < dist2_min) {
	     i_min = i0;
	     j_min = j0;
	     dist2_min = dist2;
	  }
	}
     }
   }

   if (i_min < 0)  /* no gamma_n info for this location in gamma.nc */
      return(-1);
      
      
   ioce = iocean0[j_min][i_min];

/*   label the cast  */

   dx = ABS(fmod(along,(double)ndx));
   dy = ABS(fmod((alat + 80.),(double)ndy));
   rx = dx / (double) ndx;
   ry = dy / (double) ndy;

   for (k = 0; k < n; ++k) {
      if (gamma[k] > -99.0) {
      
         thk = hb_theta(s[k], t[k], p[k], pr0);
	 dgamma_1 = 0.0;
	 dgamma_2_l = 0.0;
	 dgamma_2_h = 0.0;
	 wsum = 0.0;
	 nij = 0;
	 
	/* average the gammas over the box */
	
	 for (j0 = 0; j0 < 2; ++j0) {
	    for (i0 = 0; i0 < 2; ++i0) {
	       if (n0[j0][i0]) {
	       
	          if (j0 == 0) {
		    wt = rx * (1.0 - ry);
		    if (i0 == 0)
		       wt = (1.0 - rx) * (1.0 - ry);
		  }
		  else  {
		    wt = rx * ry;
		    if (i0 == 0)
		       wt = (1.0 - rx) *  ry;
		  }
		  
		  wt += 1.0e-6;
		  
		  itest = ocean_test(along, alat, ioce, along0[i0], alat0[j0], iocean0[j0][i0], p[k]);
		  
		  if (!itest) 
		     wt = 0.0;
			
		  depth_ns(s0[j0][i0], t0[j0][i0], p0, n0[j0][i0], s[k], t[k], p[k], &sns, &tns, &pns);
		  
		  if (pns > -98.9) {
		  
		    kns = indx(p0, n0[j0][i0], pns);
		    gw = gamma_qdr(p0[kns], gamma0[j0][i0][kns], a0[j0][i0][kns], p0[kns+1], gamma0[j0][i0][kns+1], pns);
		    
		    gamma_errors(s0[j0][i0], t0[j0][i0], p0, gamma0[j0][i0], a0[j0][i0], n0[j0][i0], along0[i0], alat0[j0], s[k], t[k], p[k], sns, tns, pns, kns, gw, &g1_err, &g2_l_err, &g2_h_err);
		  }
		  
		  else if (pns < -99.0) {
		  
		     gw = 0.0;
	             g1_err = 0.0;
	             g2_l_err = 0.0;
	             g2_h_err = 0.0;
		  }
		  
		  else {
		     goor(s0[j0][i0], t0[j0][i0], p0, gamma0[j0][i0], n0[j0][i0], s[k], t[k], p[k], &gw, &g1_err, &g2_l_err, &g2_h_err);

                     /*  adjust weight for gamma extrapolation */
		  
		     gn = gamma0[j0][i0][n0[j0][i0]-1];  /* deepest gamma0 */
                     if (gw > gn) {
		       rw = (dgw_max <= (gw - gn)) ? 1.0 : (gw - gn)/ dgw_max;
		       wt = (1.0 - rw) * wt;
		     }

		  }
		  
		  if (gw > 0.0) {
		     gamma[k] += wt *gw;
		     dgamma_1 += wt * g1_err;
		     dgamma_2_l = dgamma_2_l >= g2_l_err ? dgamma_2_l: g2_l_err;
		     dgamma_2_h = dgamma_2_h >= g2_h_err? dgamma_2_h: g2_h_err;
		     wsum += wt;
		     wtij[nij] = wt;
		     gwij[nij] = gw;
		     ++nij;
		  }
	       
	       } /* end if n0 */
	    }  /* end for i0 */
	 }   /* end for j0 */
	 
	 
	 /*  the average */
	 
	 if (wsum != 0.) {
	    gamma[k] /= wsum;
	    dgamma_1 /= wsum;
	    
	    /* the gamma errors */
	    
	    dgamma_3 = 0.0;
	    for (ij = 0; ij < nij; ++ij) 
	       dgamma_3 += wtij[ij] * ABS(gwij[ij] - gamma[k]);
	       
            dgamma_3 /= wsum;
	    
	    dg_lo[k] = dgamma_0;
	    if (dg_lo[k] < dgamma_1) dg_lo[k] = dgamma_1;
	    if (dg_lo[k] < dgamma_2_l) dg_lo[k] = dgamma_2_l;
	    if (dg_lo[k] < dgamma_3) dg_lo[k] = dgamma_3;
	    
	    dg_hi[k] = dgamma_0;
	    if (dg_hi[k] < dgamma_1) dg_hi[k] = dgamma_1;
	    if (dg_hi[k] < dgamma_2_h) dg_hi[k] = dgamma_2_h;
	    if (dg_hi[k] < dgamma_3) dg_hi[k] = dgamma_3;
	    
	 }
	 
	 else {
	    gamma[k] = -99.0;
	    dg_lo[k] = -99.0;
	    dg_hi[k] = -99.0;
	 }
	 
      
      }  /* end if gamma[k] > -99 */
   }  /* end for k */
   
   if (ialtered == 1)
      along -= 360.0;
      
   if (ialtered == 2)
      along = 360.0;
      
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2;  ++j) {
         free((void *) gamma0[i][j]);
	 free((void *) s0[i][j]);
	 free((void *) t0[i][j]);
	 free((void *) a0[i][j]);
      }
   }
   for (i = 0; i < 2; ++i) {
      free((void *) gamma0[i]);
      free((void *) s0[i]);
      free((void *) t0[i]);
      free((void *) a0[i]);
      free((void *) n0[i]);
      free((void *) iocean0[i]);
   }
   
   free((void *) gamma0);
   free((void *) s0);
   free((void *) t0);
   free((void *) a0);
   free((void *) n0);
   free((void *) iocean0);

   return(1);
   
}  /* end gamma_n() */
/*******************************************************/


void gamma_nc_init(char *filename, struct GAMMA_NC *gfile)

/*
 	DESCRIPTION :	Opens file and initializes the variable id's
	                and record variables in struct GAMMA_NC.
  
 	INPUT :		filename        full pathname to file gamma.nc
	                gfile		pointer to a structure

        An error causes an exit with an appropriate message 
	written to stderr.
*/

{
   int row, col, k, nx, ny, nz;
   int id, ierr;
   size_t len;
   int ix[NX * NY]; 
   float x[NX];    /* netcdf file has float values, NX is largest dimension */
   
/* initialize from definitions in hb_gamma.h */

   nx = (int) NX;
   ny = (int) NY;
   nz = (int) NZ;
   
/* open file */
   
   ierr = nc_open(filename, NC_NOWRITE, &gfile->id_gnc);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: Unable to open %s\n", filename);
      exit(1);
   }


/* check dimensions */
   
   ierr = nc_inq_dimid(gfile->id_gnc, "lat",&id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No lat dimension in %s\n", filename);
      exit(1);
   }
   ierr = nc_inq_dimlen(gfile->id_gnc, id,  &len);
   if (len != ny) {
      fprintf(stderr,"\nFATAL ERROR: lat dimension [%d] differs from NY [%d]\n", (int)len, ny);
      exit(1);
   }
   
   ierr = nc_inq_dimid(gfile->id_gnc, "lon", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No lon dimension in %s\n", filename);
      exit(1);
   }
   ierr = nc_inq_dimlen(gfile->id_gnc, id,  &len);
   if (len != nx) {
      fprintf(stderr,"\nFATAL ERROR: lon dimension [%d] differs from NX [%d]\n", (int)len, nx);
      exit(1);
   }
   
   ierr = nc_inq_dimid(gfile->id_gnc, "pressure",&id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No pressure dimension in %s\n", filename);
      exit(1);
   }
   ierr = nc_inq_dimlen(gfile->id_gnc, id,  &len);
   if (len != nz) {
      fprintf(stderr,"\nFATAL ERROR: pressure dimension [%d] differs from NZ [%d]\n", (int)len, nz);
      exit(1);
   }
   
      
   /* get longitude values and store in structure */

   ierr = nc_inq_varid(gfile->id_gnc, "lon", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [lon] in %s\n", filename);
      exit(1);
   }
   
   ierr = nc_get_var_float(gfile->id_gnc, id, x);
   for (k = 0; k < nx; ++k) {
      gfile->along_d[k] = (double) x[k];
   }
   
   /* get latitude values  */
   
   ierr = nc_inq_varid(gfile->id_gnc, "lat", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [lat] in %s\n", filename);
      exit(1);
   }
   
   ierr = nc_get_var_float(gfile->id_gnc, id, x);
   for (k = 0; k < ny; ++k) {
      gfile->alat_d[k] = (double) x[k];
   }
   
   /* get pressure values  */
   
   ierr = nc_inq_varid(gfile->id_gnc, "pressure", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [pressure] in %s\n", filename);
      exit(1);
   }
   
   ierr = nc_get_var_float(gfile->id_gnc, id, x);
   for (k = 0; k < nz; ++k) {
      gfile->p0[k] = (double) x[k];
   }
   
   /* get iocean values  */
   
   ierr = nc_inq_varid(gfile->id_gnc, "iocean", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [iocean] in %s\n", filename);
      exit(1);
   }
   ierr = nc_get_var_int(gfile->id_gnc, id, ix);
   k = 0;
   for (row = 0; row < ny; ++row) {
      for (col = 0; col < nx; ++col) {
         gfile->iocean[row][col] = ix[k++];
      }
   }
   
   /* get n values  */
   
   ierr = nc_inq_varid(gfile->id_gnc, "n", &id);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [n: #_of_bottles] in %s\n", filename);
      exit(1);
   }
   ierr = nc_get_var_int(gfile->id_gnc, id, ix);
   k = 0;
   for (row = 0; row < ny; ++row) {
      for (col = 0; col < nx; ++col) {
         gfile->n[row][col] = ix[k++];
      }
   }
   
/* get remaining variable ids */

   ierr = nc_inq_varid(gfile->id_gnc, "s", &gfile->id_s);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [s] in %s\n", filename);
      exit(1);
   }

   ierr = nc_inq_varid(gfile->id_gnc, "t", &gfile->id_t);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [t] in %s\n", filename);
      exit(1);
   }
   
   ierr = nc_inq_varid(gfile->id_gnc, "a", &gfile->id_a);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [a] in %s\n", filename);
      exit(1);
   }
   
   ierr = nc_inq_varid(gfile->id_gnc, "gamma", &gfile->id_gamma);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"\nFATAL ERROR: No variable [gamma] in %s\n", filename);
      exit(1);
   }
   
   
   return;
} /* end gamma_nc_init() */


/**************************************************************/

double gamma_qdr(double pl, double gl, double a, double pu, double gu, double p)
/*

 	DESCRIPTION :	Evaluate the quadratic gamma profile at a pressure
 			between two bottles
 
 	PRECISION :	Double
 
 	INPUT :		pl, pu		bottle pressures
 			gl, gu		bottle gamma values
 			a		quadratic coefficient
 			p		pressure for gamma value
 
 	returns :	gamma value at p
 
 	UNITS :		pressure	db
 			gamma		kg m-3
 			a		kg m-3
 
 
 	AUTHOR :	David Jackett
*/
{
   double p1, p2;

  p1 = (p - pu) / (pu - pl);
  p2 = (p - pl) / (pu - pl);
  
  return((a * p1 + (gu-gl)) * p2 + gl);

} /* end gamma_qdr() */
/*******************************************************/

void goor_solve(double sl, double tl, double el, double su, double tu, double eu, double p, double s0, double t0, double p0, double sigb, double *sns_ptr, double *tns_ptr)

/*
	DESCRIPTION :	Find the intersection of a potential density surface 
			between two bottles using a bisection method

	PRECISION :	Double

	INPUT :		sl, su		bottle salinities
			tl, tu		bottle in situ temperatures
			el, eu		bottle e values
			p		bottle pressures (the same)
			s0		emanating bottle salinity
			t0		emanating bottle in situ temperature
			p0		emanating bottle pressure

	returns :	*sns_ptr	salinity of the neutral surface
					intersection with the bottle pair
			*tns_ptr	in situ temperature of the intersection


	UNITS :		salinities	psu (IPSS-78)
			temperatures	degrees C (IPTS-68)
			pressures	db

*/
{
  int    iter;
  double rm, thm, sm, tm, em;
  double rl, ru, pmid;
  double thl, thu;
  double sd, sigma;
  
  rl = 0.0;
  ru = 1.0;
  pmid = (p + p0) * 0.5;
  
  thl = hb_theta(sl, tl, p, pmid);
  thu = hb_theta(su, tu, p, pmid);
  
  iter = 0;
  
  do {
     rm = (rl + ru) * 0.5;
     sm = sl + rm * (su - sl);
     thm = thl + rm * (thu - thl);
     tm = hb_theta(sm, thm, pmid, p);
     sd = hb_svan(sm, thm, pmid, &sigma);
     em = sigma - sigb;
     
     if (em == 0.0) {
        *sns_ptr = sm;
	*tns_ptr = tm;
	return;
     }
     
     if (el*em < 0.0) {
        ru = rm;
	eu = em;
     }
     else {
       if (em*eu < 0.0) {
         rl = rm;
	 el = em;
       }
     }
      
     if (ABS(em) <= 5.0e-5 && ABS(ru-rl) <= 5.0e-3) {
        *sns_ptr = sm;
	*tns_ptr = tm;
	return;
     }
  } while (++iter <= 20);
  
  *sns_ptr = sm;
  *tns_ptr = tm;
  fprintf(stderr,"\n WARNING from goor_solve(): max iterations exceeded\n");
  return;

}  /* end goor_solve() */


/********************************************************************/

void goor(double *s, double *t, double *p, double *gamma, int n, double sb, double tb, double pb, double *gammab_ptr, double *g1_err_ptr, double *g2_l_err_ptr, double *g2_h_err_ptr)

/*
 	DESCRIPTION :	Extend a cast of hydrographic data so that 
 			a bottle outside the gamma range of the cast 
 			can be labelled with the neutral density variable
 
 	PRECISION :	Double
 
 	INPUT :		s(n)		array of cast salinities
 			t(n)		array of cast in situ temperatures
 			p(n)		array of cast pressures
 			gamma(n)	array of cast gammas
 			n		length of cast
 			sb		bottle salinity
 			tb		bottle temperature
 			pb		bottle pressure
 
 	returns :	*gammab_ptr	bottle gamma value
 			*g1_err_ptr	bottle Type i error estimate
 			*g2_l_err_ptr	bottle Type ii lower error estimate
 			*g2_h_err_ptr	bottle Type ii upper error estimate
 
 	UNITS :		salinity	psu (IPSS-78)
 			temperature	degrees C (IPTS-68)
 			pressure	db
 			gamma		kg m-3
 

*/

{

   int n_sth, i, j;
   double delt_b, delt_t, slope, pr0, Tbp;
   double pmid, bmid, sigma, sd, sigb, tref;   
   double s_new, t_new, e_new;
   double s_old, t_old, e_old;
   double sns, tns, pns, sigu, sigl;
   double thb, b, thns, rho_ns, sig_ns;
   double dth, dp, g2err;

      
   delt_b = -0.1;
   delt_t =  0.1;
   slope = -0.14;
   pr0 = 0.0;
   Tbp = 2.7e-8;
   
   /* determine if its bottom data  */

   i = n-1;   
   pmid = (p[i] + pb) * 0.5;
   tref = hb_theta(s[i], t[i],p[i],pmid);
   sd = hb_svan(s[i], tref, pmid, &sigma);
   tref = hb_theta(sb, tb, pb, pmid);
   sd = hb_svan(sb, tref, pmid, &sigb);
   
   if (sigb > sigma) { 
      
       /* extend the cast data until it is denser */
       
      n_sth = 0;
      s_new = s[i];
      t_new = t[i];
      e_new = sigma - sigb;
      
      while (sigma < sigb) {
 	 s_old = s_new;
	 t_old = t_new;
         e_old = e_new;
         ++n_sth;
	 s_new = s[i] + n_sth * delt_b * slope;
	 t_new = t[i] + n_sth * delt_b;
	 
	 sd = hb_svan(s_new, hb_theta(s_new, t_new, p[i], pmid), pmid, &sigma);
	 e_new = sigma - sigb;
	 
      } /* end while */
      
      
      if (sigma == sigb) {
        sns = s_new;
	tns = t_new;
      }
      else
        goor_solve(s_old, t_old, e_old, s_new, t_new, e_new, p[i], sb, tb, pb, sigb, &sns, &tns);
	
	
      /*  compute the new gamma value */
      
      j = i-1;
      sig_vals(s[j], t[j], p[j], s[i], t[i], p[i], &sigl, &sigu);
      bmid = (gamma[i] - gamma[j]) / (sigu - sigl);
      
      sd = hb_svan(s[i], t[i], p[i], &sigl);
      sd = hb_svan(sns, tns,p[i], &sigu);
      
      *gammab_ptr = gamma[i] + bmid * (sigu - sigl);
      pns = p[i];
   
   } /* end if (sigb > sigma) */
   
   else {   /* determine if the extension is at the top */
      pmid = (p[0] + pb) * 0.5;
      sd = hb_svan(s[0], hb_theta(s[0], t[0], p[0], pmid), pmid, &sigma);
      sd = hb_svan(sb, hb_theta(sb, tb, pb, pmid), pmid, &sigb);
      
      if (sigb < sigma) {
      
      /* extend the cast until it is lighter */
      
         n_sth = 0;
	 s_new = s[0];
	 t_new = t[0];
	 e_new = sigma - sigb;
	 
	 while (sigma > sigb) {
	   s_old = s_new;
	   t_old = t_new;
	   e_old = e_new;
	   ++n_sth;
	   s_new = s[0];
	   t_new = t[0] + n_sth * delt_t;
	   
	   sd = hb_svan(s_new, hb_theta(s_new, t_new, p[0], pmid), pmid, &sigma);
	   e_new = sigma - sigb;
	   
	 }  /* end while */
	 
	 if (sigma == sigb) {
	    sns = s_new;
	    tns = t_new;
	 }
	 else{
	    goor_solve(s_new, t_new, e_new, s_old, t_old, e_old, p[0], sb, tb, pb, sigb, &sns, &tns);
	 } 
	 
	 /* compute the new gamma value */  
	    
         sig_vals(s[0], t[0], p[0], s[1], t[1], p[1], &sigl, &sigu);
	 bmid = (gamma[1] - gamma[0]) / (sigu - sigl);
	 
	 sd = hb_svan(sns, tns, p[0], &sigl);
	 sd = hb_svan(s[0], t[0], p[0], &sigu);
	 
	 *gammab_ptr = gamma[0] - bmid * (sigu - sigl);
	 
	 pns = p[0];
      
      } /* end if sigb < sigma */
      
      else {
      
        fprintf(stderr,"\nFATAL ERROR in goor(): gamma out of range\n");
	exit(1);
      }
   
   } /* end else */
   
   
   /*  compute an error estimate */
   
   thb = hb_theta(sb, tb, pb, pr0);
   thns = hb_theta(sns, tns, pns, pr0);
   
   sd = hb_svan(sns, tns, pns, &sig_ns);
   rho_ns = 1000. + sig_ns;
   
   b = bmid;
   
   dp = pns - pb;
   dth = thns - thb;
   
   *g1_err_ptr = rho_ns * b * Tbp * ABS(dp *dth) / 6.0;
   
   g2err = rho_ns * b * Tbp * dp * dth * 0.5;
   
   if (g2err <= 0.0) {
      *g2_l_err_ptr = -g2err;
      *g2_h_err_ptr = 0.0;
      return;
   }
   
   *g2_l_err_ptr = 0.0;
   *g2_h_err_ptr = g2err;
   return;
   
}  /* end goor() */
/********************************************************************/

int indx(double *x, int n, double z)
/*
	DESCRIPTION :	Find the index of a real number in a
			monotonically increasing real array

	PRECISION :	Double

	INPUT :		x		array of increasing values
			n		length of array
			z		real number

	returns :	index k - if x(k) <= z < x(k+1), or
			n-1     - if z = x(n)

*/

{
   int k, kl, ku, km;
   
   if (x[0] < z && z < x[n-1]) {
   
      kl = 0;
      ku = n-1;
      
      while ((ku - kl) > 1)  {
         km = (ku + kl) / 2;
	 if (z > x[km])
	   kl = km;
	 else
	   ku = km;
      }
      k = kl;
      
      if (z == x[k+1])
         ++k;
      return(k);
   }
   
   if (z == x[0]) 
      return(0);
      
   if (z == x[n-1])
      return(n-2);
      
   /* if we get here something is wrong */
   
   fprintf(stderr,"\nFATAL ERROR in indx():  out of range\n");
   exit(1);
   
}  /* end indx() */

/*******************************************************/
void neutral_surfaces(double *s, double *t, double *p, double *gamma, int n, double *glevels, int ng, double *sns, double *tns, double *pns, double *dsns, double *dtns, double *dpns)

/*
	DESCRIPTION :	For a cast of hydrographic data which has been 
			labelled with the neutral density variable gamma,
			find the salinities, temperatures and pressures
			on ng specified neutral density surfaces.  Cast
			must be continuous -- no missing values allowed.


	INPUT :		s[n]		array of cast salinities
			t[n]		array of cast in situ temperatures
			p[n]		array of cast pressures
			gamma[n]	array of cast gamma values
			n		length of cast
			glevels[ng]	array of neutral density values
			ng		number of neutral density surfaces

	Returns :	sns[ng]		salinity on the neutral density surfaces
			tns[ng]		in situ temperature on the surfaces
			pns[ng]		pressure on the surfaces
			dsns[ng]	surface salinity errors
			dtns[ng]	surface temperature errors
			dpns[ng]	surface pressure errors

			NOTE:		sns, tns and pns values of -99.0
					denotes under or outcropping

					non-zero dsns, dtns and dpns values
					indicates multiply defined surfaces,
					and file 'ns-multiples.dat' contains
					information on the multiple solutions

	UNITS :		salinity	psu (IPSS-78)
			temperature	degrees C (IPTS-68)
			pressure	db
			gamma		kg m-3

*/

{
   int intvl[MAX_INTVLS];
   int i, ii, ig, k, k1, nintvls;
   int nlast, n2 = 2;
   int halfcastlen, middle;
   int mult;   
   
   double alfa_l, alfa_u, alfa_mid, beta_l, beta_u, beta_mid;
   double thl, thu, delth, dels, pl, pu, delp, delp2;
   double gmin, gmax;
   double a, b, c, q;
   double pr0 = 0.0;
   double ptol = 1.0e-3;
   double rhomid, bmid, sigmid, smid, tmid, pmid;
   double sns_top, tns_top, pns_top;
   double sns_middle, tns_middle, pns_middle;
   double sdum, thdum, sthdum, gdum;
   double pns1, pns2;
   double bden, rg, plast;
   double *swork, *twork, *pwork, *gwork;
   
   FILE *mfile; 
   
   mfile = (FILE *) NULL; 
   
   swork = (double *) calloc((size_t)n, sizeof(double));
   twork = (double *) calloc((size_t)n, sizeof(double));
   pwork = (double *) calloc((size_t)n, sizeof(double));
   gwork = (double *) calloc((size_t)n, sizeof(double));
   
/* check for missing gammas and monotonic pressure series */

   plast = -10.0;
   i = 0;
   for (k = 0; k < n; ++k) {
      if (gamma[k] > 0.0) {
         if (p[k] > plast) {
	    plast = p[k];
            swork[i] = s[k];
            twork[i] = t[k];
            pwork[i] = p[k];
            gwork[i] = gamma[k];
	    ++i;
	 }
      }
   }
   
   if (i <= 1) {
      for (ii = 0; ii < ng; ++ii) {
         sns[ii] = -99.0;
         tns[ii] = -99.0;
         pns[ii] = -99.0;
         dsns[ii] = 0.0;
         dtns[ii] = 0.0;
         dpns[ii] = 0.0;
      }
   }
   
   if (i == 0) {        /* no gamma-n information */
      free((void *)swork);
      free((void *)twork);
      free((void *)pwork);
      free((void *)gwork);
      return;
   }
      
   n = i;
   nlast = n - 1;
   halfcastlen = n * 0.5;
   
   /* loop for each neutral surface */

   for (ig = 0; ig < ng; ++ig) {
      
      nintvls = 0;
      
      for (k = 0; k < nlast; ++k) {
      
        gmin = gwork[k];
	gmax = gwork[k];
	
	if (n > 1) {
          k1 = k + 1;
          gmin = MIN(gwork[k], gwork[k1]); 
          gmax = MAX(gwork[k], gwork[k1]);
	}
	
	if (glevels[ig] >= gmin && glevels[ig] <= gmax) {
	  if (nintvls == MAX_INTVLS) {
	    fprintf(stderr,"\nFATAL ERROR:  too many crossings in neutral_surfaces()\n");
	    exit(1);
	  }
	  
	  intvl[nintvls] = k;
	  ++nintvls;
	}
      }  /* end for k */
      
      
      sns[ig] = -99.0;
      tns[ig] = -99.0;
      pns[ig] = -99.0;
      dsns[ig] = 0.0;
      dtns[ig] = 0.0;
      dpns[ig] = 0.0;


    /* unusual case of an exact crossing with castlength of 1 */
    
      if ((n == 1) && (nintvls == 1)) {  
         sns[ig] = swork[0];
	 tns[ig] = twork[0];
	 pns[ig] = pwork[0];
      }
                  
    /* more usual case:  neutral surface exists, castlength > 1 */
      else if (n > 1 && nintvls > 0) {
  
        /* if more than 1 interval, choose the median  */
	
	middle = nintvls  * 0.5;
	
	if (nintvls == 2)
	  middle = 1;
		
	  
	/* loop over all intersections */
	
	for (i = 0; i < nintvls; ++i) {
	   k = intvl[i];
	   k1 = k + 1;
	   
           /* coefficients of a quadratic for gamma  */
	   
	   eosall(swork[k], twork[k], pwork[k], &thdum, &sthdum, &alfa_l, &beta_l, &gdum, &sdum);
	   eosall(swork[k1], twork[k1], pwork[k1], &thdum, &sthdum, &alfa_u, &beta_u, &gdum, &sdum);
	   
	   alfa_mid = (alfa_l + alfa_u) * 0.5;
	   beta_mid = (beta_l + beta_u) * 0.5;
	   pmid = (pwork[k] + pwork[k1]) * 0.5;
	   
	   stp_interp(&swork[k], &twork[k], &pwork[k], n2, &smid, &tmid, &pmid);
	   
	   sdum = hb_svan(smid, tmid, pmid, &sigmid);
	   rhomid = 1000.0 + sigmid;
	   
	   thl = hb_theta(swork[k], twork[k], pwork[k], pr0);
	   thu = hb_theta(swork[k1], twork[k1], pwork[k1], pr0);
	   delth = thu - thl;
	   
	   dels = swork[k1] - swork[k];
	   
	   pl = pwork[k];
	   pu = pwork[k1];
	   delp = pu - pl;
	   delp2 = delp * delp;
	   
	   bden = rhomid * (beta_mid * dels - alfa_mid * delth);
	   
	   if (ABS(bden) < 1.0e-6)
	      bden = 1.0e-6;
	      
	   bmid = (gwork[k1] - gwork[k]) / bden;
	   
           /*	coefficients */
	   
	   a = dels * (beta_u-beta_l) - delth * (alfa_u-alfa_l);
	   a = (a * bmid * rhomid) / (2 * delp2);
	   
	   b = dels * (pu*beta_l - pl*beta_u) - delth * (pu*alfa_l - pl*alfa_u);
	   b = (b * bmid * rhomid) / delp2;
	   
	   c = dels * (beta_l * (pl - 2.* pu) + beta_u * pl) - delth * (alfa_l * (pl - 2.0 * pu) + alfa_u * pl); 
	   
	   c = gwork[k] + (bmid * rhomid * pl * c)/ (2*delp2);
	   c -=  glevels[ig];
   
	   /* solve the quadratic */
	   
	   if (a != 0.0 && bden > 1.0e-6) {
	   
	       mult = (b >= 0.0) ? 1 : -1;
	   
	       q = -(b + mult * sqrt(b*b - 4*a*c)) * 0.5;

	       pns1 = q / a;
	       pns2 = c / q;
	       
	       if (pns1 >= (pwork[k]-ptol) && pns1 <= (pwork[k1]+ptol)) 
	         pns[ig] = MIN(pwork[k1], (MAX(pns1, pwork[k])));
	       else if (pns2 >= (pwork[k]-ptol) && pns2 <= (pwork[k1]+ptol)) 
	         pns[ig] = MIN(pwork[k1], (MAX(pns2, pwork[k])));
	       else {
	         pns[ig] = -99.0;
		 fprintf(stderr,"Error 3 (quadratic sol'n) in neutral_surfaces()\n");
	       }
	   }
	   else {
	      rg = (glevels[ig] - gwork[k]) / (gwork[k1] - gwork[k]);
	      pns[ig] = pwork[k] + rg *(pwork[k1] - pwork[k]);
	   
	   } 
	   
	   if (pns[ig] >= 0.0)
	      stp_interp(swork, twork, pwork, n, &sns[ig], &tns[ig], &pns[ig]); 
	   
   /* case of multiple intersections */
   
           if (nintvls <= 1)  {
	      dsns[ig] = 0.0;
	      dtns[ig] = 0.0;
	      dpns[ig] = 0.0;
	   }
	   else {
	   
	     /* write multiples to file */
	     
	     if (mfile == NULL)
	          mfile = fopen("ns-multiples.dat", "a");
	     if (mfile == NULL) {
	        fprintf(stderr,"\nWARNING:  unable to open ns-multiples.dat for writing\n");
	     }
	     else {
	        if (i == 0) 
	           fprintf(mfile, "\n level: [%8.4lf]   # of intersections: [%d]", glevels[ig], nintvls);
	        fprintf(mfile,"\n%8.4lf %8.4lf %8.1lf ", sns[ig], tns[ig], pns[ig]); 
	     }
	     
	     /* find median values and errors */
	     
	     
	     if (i == 0) {
	  	  sns_top = sns[ig];
	  	  tns_top = tns[ig];
	  	  pns_top = pns[ig];
	     
	     }
	     
	     if (i == middle) {
	  	  sns_middle = sns[ig];
	  	  tns_middle = tns[ig];
	  	  pns_middle = pns[ig];
	     }
	     
	     if (i == (nintvls-1)) {
	     
		  if((pns_middle-pns_top) > (pns[ig]-pns_middle)) {
		    dsns[ig] = sns_middle - sns_top;
		    dtns[ig] = tns_middle - tns_top;
		    dpns[ig] = pns_middle - pns_top;
		  }
		  else {
		    dsns[ig] = sns[ig] - sns_middle;
		    dtns[ig] = tns[ig] - tns_middle;
		    dpns[ig] = pns[ig] - pns_middle;
		  }
		  
	  	  sns[ig] = sns_middle;
		  tns[ig] = tns_middle;
		  pns[ig] = pns_middle;
	      }
	      
	   } /* end else */
	   
	} /* end for i */

      }  /* end if - else if */
   
   }  /* end for ig */
   
   if (mfile != NULL) {
      fclose(mfile);
      mfile = (FILE *) NULL;
   }
   
   free((void *)swork);
   free((void *)twork);
   free((void *)pwork);
   free((void *)gwork);
   
   return;
}  /* end neutral_surfaces() */
/*******************************************************/

int ocean_test(double x1, double y1, int io1, double x2, double y2, int io2, double z)

/*
 	DESCRIPTION :	Test whether two locations are connected by ocean
 
 	PRECISION :	Double
 
 	INPUT :		x1		longitude of first location
 			y1		latitude of first location
 			io1		ocean of first location
 			x2		longitude of second location
 			y2		latitude of second location
 			io2		ocean of second location
 			z		depth of connection
 
 	returns :	0 or 1		success of connection (0 = no connect)
 

*/

{
   int itest, isj1, isj2;
   double y, x_js[3], y_js[3];
   double em1, em2, c1, c2;
   
   x_js[0] =  129.87;
   x_js[1] =  140.37;
   x_js[2] =  142.83;

   y_js[0] =  32.75;
   y_js[1] =  37.38;
   y_js[2] =  53.58;
   
   
   y = (y1 + y2) * 0.5;

/*	same ocean talks  */

   if (io1 == io2)
      return(1);
      

/*	Atlantic talks  */
      
   if (((io1 == 5) || (io1 == 6)) && ((io2 == 5) || (io2 == 6)))
      return(1);
      
/*  exclude Antarctic tip  */

   if( ((io1 * io2) == 12) && (y < -60.)) 
     return(0);
     

      
   if (y <= -20.) {
   
/*	land of South America doesn't talk  */

      if (y >= -48. && (io1 * io2) == 12)
          return(0);
   
/*	everything else south of -20 talks  */

      return(1);
   }
   

/* test multiple conditions ... */

   itest = 0;
   
     
/*      Pacific talks  */

   if (((io1 ==1) || (io1 == 2)) && ((io2 == 1) || (io2 == 2)))
      itest = 1;
   
   
/*	Indian talks  */

   if (((io1 == 3) || (io1 == 4)) && ((io2 == 3) || (io2 == 4)))
      itest = 1;
      
      
/*	Indonesian throughflow  */
      
   if (((io1 * io2) == 8) && (z <= 1200.) && (x1 >= 124.) && (x1 <= 132.) && (x2 >= 124.) && (x2 <= 132)) 
      itest = 1; 


/* 	exclude Japan Sea from talking  */

   if ( (x_js[0] <= x1 && x1 <= x_js[2] && y_js[0] <= y1 && y1 <= y_js[2])
      || (x_js[0] <= x2 && x2 <= x_js[2] && y_js[0] <= y2 && y2 <= y_js[2]) )  {
      
      em1 = (y_js[1] - y_js[0]) / (x_js[1] - x_js[0]);
      c1 = y_js[0] - em1 * x_js[0];
      
      em2 = (y_js[2] - y_js[1]) / (x_js[2] - x_js[1]);
      c2 = y_js[1] - em2 * x_js[1];

      isj1 = 0;      
      if (((y1 - em1 * x1 - c1) >= 0) && ((y1-em2*x1-c2) >= 0.))
         isj1 = 1;
	 
      isj2 = 0;
      if (((y2-em1*x2-c1) >= 0.) && ((y2-em2*x2-c2) >= 0.)) 
         isj2 = 1;
	 
   
      if (isj1 == isj2)
         return (1);
	
      return(0);
   }

   return(itest);
   
}  /* end ocean_test() */

/*******************************************************/

void read_nc(struct GAMMA_NC* gfile, double along, double alat, double ***s0, double ***t0, double *p0, double ***gamma0, double ***a0, int **n0, double *along0, double *alat0, int **iocean0)

/*
 	DESCRIPTION :	Read variables from the netcdf labelled data file. 
	The file must be opened by the HydroBase function gamma_nc_init() 
	-- which initializes struct GAMMA_NC *gfile with various info -- 
	and closed by the calling routine using nc_close() from the netcdf 
	C library functions version 3.5 or later.
  
 	INPUT :		gfile           ptr to various info (structure)
	                along		longitude of record
 			alat		latitude of record


        Memory for these arrays is allocated by the calling routine --
	not here.  The NZ dimension is defined in hb_gamma.h.  The order of
	dimensions in the gamma.nc file is lat, lon, pressure;  therefore
	the dimension of each of these C-arrays is [ny][nx][nz] with nz
	varying fastest.
	
	
 	Returns :	s0[2][2][NZ]	arrays of cast salinities
 			t0[2][2][NZ]	arrays of cast in situ temperatures
 			gamma0[2][2][NZ]	array of cast gamma values
 			a0[2][2][NZ]	arrays of cast a values
 			n0[2][2]	length of casts
 			along0[2]	array of cast longitudes
 			alat0[2]	array of cast latitudes
 			iocean0[2][2]	array of cast oceans
 
 	UNITS :		salinity	psu (IPSS-78)
 			temperature	degrees C (IPTS-68)
 			pressure	db
 			gamma		kg m-3
 
*/

{
   int k, nx, ny, nz, ndx, ndy;
   int start[3], count[3], ierr;
   int row0, col0, row1, col1;  /* row0 corresponds to j0 
                                   col0 corresponds to i0  */
   int i0, i1, j0, j1;				   
   
   float s0_s[NZ], t0_s[NZ];    /* netcdf file has float values */
   float a0_s[NZ], gamma0_s[NZ];
   

/* initialize from definitions in hb_gamma.h */

   nx = (int) NX;
   ny = (int) NY;
   nz = (int) NZ;
   ndx = (int) NDX;
   ndy = (int) NDY;


/* find indices for corners of box containing alat,along */ 
  
   col0 = (int) (along / ndx);           /* left  */
   row0 = (int) ((alat + 88.) / ndy);    /* lower  */
   row1 = row0 + 1;                      /* upper  */
   col1 = col0 + 1;                      /* right */
   if (col1 == nx)
      col1 = 0;

         
   alat0[0] = gfile->alat_d[row0];
   alat0[1] = gfile->alat_d[row1];
   along0[0] = gfile->along_d[col0]; 
   along0[1] = gfile->along_d[col1];
   
   i0 = j0 = 0;
   i1 = j1 = 1;
   
/* get lower/left corner (row0,col0) (j0,i0) */
      
   start[0] = row0;
   start[1] = col0;
   start[2] = 0;
   
   count[0] = 1;
   count[1] = 1;
   count[2] = nz;
   
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_s, (const size_t *) start, (const size_t *)count, s0_s);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"ERROR #[%d] reading gamma.nc variable\n", ierr);
      exit(1);
   }
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_t, (const size_t *)start, (const size_t *)count, t0_s);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"ERROR #[%d] reading gamma.nc variable\n", ierr);
      exit(1);
   }
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_a, (const size_t *)start, (const size_t *)count, a0_s);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"ERROR #[%d] reading gamma.nc variable\n", ierr);
      exit(1);
   }
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_gamma, (const size_t *)start, (const size_t *)count, gamma0_s);
   if (ierr != NC_NOERR) {
      fprintf(stderr,"ERROR #[%d] reading gamma.nc variable\n", ierr);
      exit(1);
   }

   /* store as double values in returned arrays */
      
   for (k = 0; k < nz; ++k) {
      s0[j0][i0][k] = (double) s0_s[k];
      t0[j0][i0][k] = (double) t0_s[k];
      a0[j0][i0][k] = (double) a0_s[k];
      gamma0[j0][i0][k] = (double) gamma0_s[k];
   }
   
/* get lower/right corner (row0,col1) (j0,i1) */
      
   start[0] = row0;
   start[1] = col1;
   start[2] = 0;
      
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_s, (const size_t *)start, (const size_t *)count, s0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_t, (const size_t *)start, (const size_t *)count, t0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_a, (const size_t *)start, (const size_t *)count, a0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_gamma, (const size_t *)start, (const size_t *)count, gamma0_s);

   /* store as double values in returned arrays */
      
   for (k = 0; k < nz; ++k) {
      s0[j0][i1][k] = (double) s0_s[k];
      t0[j0][i1][k] = (double) t0_s[k];
      a0[j0][i1][k] = (double) a0_s[k];
      gamma0[j0][i1][k] = (double) gamma0_s[k];
   }
   
/* get upper/left corner (row1,col0) (j1,i0) */
      
   start[0] = row1;
   start[1] = col0;
   start[2] = 0;
   
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_s, (const size_t *)start, (const size_t *)count, s0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_t, (const size_t *)start, (const size_t *)count, t0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_a, (const size_t *)start, (const size_t *)count, a0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_gamma, (const size_t *)start, (const size_t *)count, gamma0_s);

   /* store as double values in returned arrays */
      
   for (k = 0; k < nz; ++k) {
      s0[j1][i0][k] = (double) s0_s[k];
      t0[j1][i0][k] = (double) t0_s[k];
      a0[j1][i0][k] = (double) a0_s[k];
      gamma0[j1][i0][k] = (double) gamma0_s[k];
   }
   
   
/* get upper/right corner (row1,col1) (j1,i1) */
      
   start[0] = row1;
   start[1] = col1;
   start[2] = 0;
   
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_s, (const size_t *)start, (const size_t *)count, s0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_t, (const size_t *)start, (const size_t *)count, t0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_a, (const size_t *)start, (const size_t *)count, a0_s);
   ierr = nc_get_vara_float(gfile->id_gnc, gfile->id_gamma, (const size_t *)start, (const size_t *)count, gamma0_s);

   /* store as double values in returned arrays */
      
   for (k = 0; k < nz; ++k) {
      s0[j1][i1][k] = (double) s0_s[k];
      t0[j1][i1][k] = (double) t0_s[k];
      a0[j1][i1][k] = (double) a0_s[k];
      gamma0[j1][i1][k] = (double) gamma0_s[k];
   }
   
   iocean0[j0][i0] = gfile->iocean[row0][col0];
   iocean0[j0][i1] = gfile->iocean[row0][col1];
   iocean0[j1][i0] = gfile->iocean[row1][col0];
   iocean0[j1][i1] = gfile->iocean[row1][col1];
   
   n0[j0][i0] = gfile->n[row0][col0];
   n0[j0][i1] = gfile->n[row0][col1];
   n0[j1][i0] = gfile->n[row1][col0];
   n0[j1][i1] = gfile->n[row1][col1];
   
   return;
   
} /* end read_nc() */

/*****************************************************/

void scv_solve(double *s, double *t, double *p, double *e, int n, int k, double s0, double t0, double p0, double *sscv_ptr, double *tscv_ptr, double *pscv_ptr, int *iter)

/*
 	DESCRIPTION :	Find the zero of the v function using a 
 			bisection method
 
 	PRECISION :	Double
 
 	INPUT :		s(n)		array of cast salinities
 			t(n)		array of cast in situ temperatures
 			p(n)		array of cast pressures
 			e(n)		array of cast e values
 			n		length of cast
 			k		interval (k-1,k) contains the zero
 			s0		the bottle salinity
 			t0		the bottle in situ temperature
 			p0		the bottle pressure
 
 	returns :	*sscv_ptr	salinity of the scv surface
 					intersection with the cast
 			*tscv_ptr	in situ temperature of the intersection
 			*pscv_ptr	pressure of the intersection
                        *iter           number of iterations
 
 	UNITS :		salinities	psu (IPSS-78)
 			temperatures	degrees C (IPTS-68)
 			pressures	db
*/

{

   int i;
   double pl, pu, pm;
   double el, eu, em;
   double sm, tm;
   double sigl, sigu, sd; 

   i = k-1;
   pl = p[i];
   el = e[i];
   pu = p[k];
   eu = e[k];
   
   *iter = 0;

   do {
   
     (*iter)++;
     
     pm = (pl + pu) * 0.5;
     
     stp_interp(&s[i], &t[i], &p[i], 2, &sm, &tm, &pm);
     sd = hb_svan(s0, hb_theta(s0,t0,p0,pm), pm, &sigl);
     sd = hb_svan(sm, tm, pm, &sigu);
     em = sigu - sigl;

     
     if (el * em < 0.0) {
        pu = pm;
	eu = em;
     }
     else if (em *eu < 0.0) {
        pl = pm;
	el = em;
     }
     else {
        if (em == 0.0) {
	  *sscv_ptr = sm;
	  *tscv_ptr = tm;
	  *pscv_ptr = pm;
	  return;
	}
     }
     
     if ((ABS(em) <= 5.0e-5) && (ABS(pu-pl) <= 5.0e-3)) {
        *sscv_ptr = sm;
	*tscv_ptr = tm;
	*pscv_ptr = pm;
	return;
     }
    
    } while (*iter < 20);

   fprintf(stderr, "\nWARNING from scv_solve()");   
   fprintf(stderr, "\niteration #%d  em: %lf  dp: %lf - %lf = %lf", *iter, (double) ABS(em), pl, pu, (double)ABS(pu-pl));   
   *sscv_ptr = -99.0;
   *tscv_ptr = -99.0;
   *pscv_ptr = -99.0;
   return;
     
} /* end scv_solve() */
/********************************************************************/

void sig_vals(double s1, double t1, double p1, double s2, double t2, double p2, double *sig1_ptr, double *sig2_ptr)
/*
 	DESCRIPTION :	Computes the sigma values of two neighbouring 
 			bottles w.r.t. the mid pressure
 
 	PRECISION :	Double
 
 	INPUT :		s1,s2		bottle salinities
 			t1,t2		bottle in situ temperatures
 			p1,p2		bottle pressures
 
 	Returns :	*sig1_ptr	bottle potential density values
	                *sig2_ptr
 
 	UNITS :		salinity	psu (IPSS-78)
 			temperature	degrees C (IPTS-68)
 			pressure	db
 			density		kg m-3
 
*/

{
  double pmid, sd;
  
  pmid = (p1 + p2) * 0.5;
  
  sd = hb_svan(s1, hb_theta(s1, t1, p1, pmid), pmid, sig1_ptr);
  sd = hb_svan(s2, hb_theta(s2, t2, p2, pmid), pmid, sig2_ptr);
  return;

}  /* end sig_vals() */
/********************************************************************/

void stp_interp (double *s, double *t, double *p, int n, double *s0_addr, double *t0_addr, double *p0_addr)

/*
 	DESCRIPTION :	Interpolate salinity and in situ temperature
 			on a cast by linearly interpolating salinity
 			and potential temperature
 
 	PRECISION :	Double
 
 	INPUT :		s(n)		array of cast salinities
 			t(n)		array of cast in situ temperatures
 			p(n)		array of cast pressures
 			n		length of cast
 			*p0_addr	pressure for which salinity and
 					in situ temperature are required
 
 	returns :	*s0_addr	interpolated value of salinity
 			*t0_addr	interpolated value of situ temperature
 
 	UNITS :		salinities	psu (IPSS-78)
 			temperatures	degrees C (IPTS-68)
 			pressures	db
 

*/
{
  int k, k1;
  double r, pr0, th0, thk;
  
  k = indx(p, n, *p0_addr);
  k1 = k+1;
  pr0 = 0.0;
  
  r = (*p0_addr - p[k]) / (p[k1] - p[k]);
  *s0_addr = s[k] + r * (s[k1] - s[k]);
  thk = hb_theta(s[k], t[k], p[k], pr0);
  th0 = thk + r * (hb_theta(s[k1],t[k1], p[k1], pr0) - thk);
  *t0_addr = hb_theta(*s0_addr, th0, pr0, *p0_addr);
  
  return;

} /* end stp_interp() */



/********************************************************************/




