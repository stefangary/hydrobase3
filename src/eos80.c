/* eos80.c

...........................................................................
                          *******  HydroBase 3  *******
               ...................................................
                          
                    author:  Ruth Gorski Curry
                             Woods Hole Oceanographic Inst
               ...................................................
	******* EOS80 derivatives temperature and salinity *******
                          by N Fofonoff & R Millard
...........................................................................
		Converted from EOS80.f by Julie Allen
		Woods Hole Oceanographic Institution 
		July 1999
____________________________________________________________________________

   double hb_eos80d(double S, double T, double P0, double **DRV)
   double hb_alpha(double S, double T, double P)
   double hb_beta(double S, double T, double P)
   double hb_gamma(double S, double T, double P)
   double hb_ratio(double S, double T, double P)
____________________________________________________________________________

	DRV Matrix Format (transposed in converting Fortran to C
      	0    	1     	2 -- rows
cols	
  0   	V,	VT,	VTT    spvol wrt temperature
  1   	V0,	VOT,	V0TT   For S,T,0 wrt temperature
  2   	RO,	ROT,	ROTT   For S,T,P  Density deriv wrt temp 
  3   	K0,	K0T,	K0TT   For S,T,0 Sec bulk mod  wrt temp
  4   	A,	AT,	ATT
  5   	B,	BT,	BTT    Bulk mod press coeffs 
  6 	DRDP,	K,	DVDP   derivatives wrt pressure  
  7   	R0S,      ,	VS     derivatives wrt salinity 

	Check value: for S = 40 (IPSS-78) , T = 40 DEG C, P0= 10000 Decibars.
        	 DR/DP                DR/DT                 DR/DS
       		DRV(0,6)              DRV(1,2)             DRV(0,7)

	Finite difference with 34d order correction done in double precision
       		3.46969238E-3       -.43311722           .705110777

 	Explicit differentiation single precision formulation EOS80
       		3.4696929E-3        -.4331173            .7051107

____________________________________________________________________________

*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>

#define EOS80 1
#include "hydrobase.h"

/***************************************************/
double hb_eos80d(double S, double T, double P0, double **DRV)
/***************************************************/

/*  S 	- salinity (IPSS-78) 
    T 	- temperature deg celcius (IPTS-68) 
    P0 	- pressure (decibars) 
    DRV - DRV matrix (3 rows * 8 cols) SPACE ALREADY ALLOCATED 
*/
{
    double P,SIG,SR,R1,R2,R3,R4;
    double A,B,C,D,E,A1,B1,AW,BW,K,K0,KW,K35;
    double R3500,SAL,V350P,SVA,SIGMA,DR350,V0;
    double RHOT,RHO2,V0T,V0S,RHOS;
    double RHOTT,VOTT,SVAN,DBDS,BT,BTT;
    double DKDP,DVDP,VT,VTT,R0TT,VS,V2,V,KTT,K0T,KT;
    double VP,DVAN,DR35P,PK,GAM,DK,K0TT,KS,K0S,ATT,AT,DADS;
    double V0TT,RHO1,DRDT,R4S;
    double eos8d;
    int i,j;

/* Define constants; 
   R4 is referred to as C in Millero and Poisson, 1981 */
   
    R3500=1028.1063;
    DR350=28.106331;
    R4=4.8314e-4;

/* Convert pressure to bars and take sqrt of salinity */
    P = P0 /10.0;
    R3500=1028.1063;
    SAL=S;
    SR=sqrt(fabs(S));

/* Pure water density at atmospheric pressure
	Bigg, P.H., 1967, Br. J. Applied Physics 8, pp. 521-537. */
    R1=((((6.536332E-9*(T)-1.120083E-6)*(T)+
	1.001685E-4)*(T)-9.095290E-3)*(T)+6.793952E-2)*(T)-28.263737;

/* Seawater density atm pressure - coefficients involving salinity
	R2=A in notatino of Millero and Poisson, 1981 */
    R2 = (((5.3875E-9*(T)-8.2467E-7)*(T)+7.6438E-5)*(T)-4.0899E-3)*(T)+8.24493E-1;

/* R3=B in notation of Millero and Poisson, 1981 */
    R3 = (-1.6546E-6*(T)+1.0227E-4)*(T)-5.72466E-3;

/* International one-atmosphere equation of state of seawater */
    SIG = (R4*(S) + R3*SR + R2)*(S) + R1;

/* Specific volume at atmospheric pressure */
    V350P = 1.0/R3500;
    SVA = -SIG*V350P/(R3500+SIG);
    SIGMA=SIG+DR350;
    DRV[0][2] = SIGMA;
    V0 = 1.0/(1000.0 + SIGMA);
    DRV[0][1] = V0;

/* Compute derivative wrt SALT of RHO */
    R4S=9.6628E-4;
    RHOS=R4S*SAL+1.5*R3*SR+R2;

/* Compute derivative w/respect to temperature of RHO */
    R1 =(((3.268166E-8*(T)-4.480332E-6)*(T)+3.005055E-4)*(T)-1.819058E-2)*(T)+6.793952E-2;
    R2 = ((2.155E-8*(T)-2.47401E-6)*(T)+1.52876E-4)*(T)-4.0899E-3;
    R3 = -3.3092E-6*(T)+1.0227E-4;
    RHOT = (R3*SR + R2)*SAL + R1;
    DRDT=RHOT;
    DRV[1][2] = RHOT;
    RHO1 = 1000.0 + SIGMA;
    RHO2 = RHO1*RHO1;
    V0T = -RHOT/(RHO2);

/* Specific volume derivative wrt S */
    V0S=-RHOS/RHO2;
    DRV[0][7]=RHOS;
    DRV[1][1] = V0T;

/* Compute second derivative of RHO */
    R1 = ((1.3072664E-7*(T)-1.3440996E-5)*(T)+6.01011E-4)*(T)-1.819058E-2;
    R2 = (6.465E-8*(T)-4.94802E-6)*(T)+1.52876E-4;
    R3 = -3.3092E-6;
    RHOTT = (R3*SR + R2)*SAL + R1;
    DRV[2][2] = RHOTT;
    V0TT = (2.0*RHOT*RHOT/RHO1 - RHOTT)/(RHO2);
    DRV[2][1] = V0TT;

/* Scale specific anomaly to normally reported units */
    SVAN=SVA*1.0E+8;
    eos8d=SVAN;

/* New high pressure equation of state for seawater
	Millero, et al, 1980 DSR 27A, pp. 255-264.
	Constant notation follows article */ 

/* Compute compression terms */
    E = (9.1697E-10*(T)+2.0816E-8)*(T)-9.9348E-7;
    BW = (5.2787E-8*(T)-6.12293E-6)*(T)+3.47718E-5;
    B = BW + E*(S);

/* Derivative of B wrt SALT */
    DBDS=E;

/* Correct B for anomaly bias change */
    DRV[0][5] = B + 5.03217E-5;

/* Derivative of B */
    BW = 1.05574E-7*(T)-6.12293E-6;
    E = 1.83394E-9*(T) +2.0816E-8;
    BT = BW + E*SAL;
    DRV[1][5] = BT;

/* Coefficients of A
	second derivative of B */
    E = 1.83394E-9;
    BW = 1.05574E-7;
    BTT = BW + E*SAL;
    DRV[2][5] = BTT;
    D = 1.91075E-4;
    C = (-1.6078E-6*(T)-1.0981E-5)*(T)+2.2838E-3;
    AW = ((-5.77905E-7*(T)+1.16092E-4)*(T)+1.43713E-3)*(T)-0.1194975;
    A = (D*SR + C)*(S) + AW;

/* Correct A for anomaly bias change */
    DRV[0][4] = A + 3.3594055;

/* Derivative of A wrt SALT */
    DADS=2.866125E-4*SR+C;

/* Derivative of A */
    C = -3.2156E-6*(T) -1.0981E-5;
    AW = (-1.733715E-6*(T)+2.32184E-4)*(T)+1.43713E-3;
    AT = C*SAL + AW;
    DRV[1][4] = AT;

/* Second derivative of A */
    C = -3.2156E-6;
    AW = -3.46743E-6*(T) + 2.32184E-4;
    ATT = C*SAL + AW;
    DRV[2][4] = ATT;

/* Coefficient K0 */
    B1 = (-5.3009E-4*(T)+1.6483E-2)*(T)+7.944E-2;
    A1 = ((-6.1670E-5*(T)+1.09987E-2)*(T)-0.603459)*(T)+54.6746;
    KW = (((-5.155288E-5*(T)+1.360477E-2)*(T)-2.327105)*(T)+148.4206)*(T)-1930.06;
    K0 = (B1*SR + A1)*(S) + KW;

/* Add bias to output K0 value */
    DRV[0][3] = K0+21582.27;

/* Derivative of K0 wrt SALT */
    K0S=1.5*B1*SR+A1;

/* Derivative of K wrt SALT */
    KS=(DBDS*P+DADS)*P+K0S;

/* Derivative of KO */
    B1 = -1.06018E-3*(T)+1.6483E-2;
    A1 = (-1.8501E-4*(T)+2.19974E-2)*(T)-0.603459;
    KW = ((-2.0621152E-4*(T)+4.081431E-2)*(T)-4.65421)*(T)+148.4206;
    K0T = (B1*SR+A1)*SAL + KW;
    DRV[1][3] = K0T;

/* Second derivative of K0 */
    B1 = -1.06018E-3;
    A1 = -3.7002E-4*(T) + 2.19974E-2;
    KW = (-6.1863456E-4*(T)+8.162862E-2)*(T)-4.65421;
    K0TT = (B1*SR + A1)*SAL + KW;
    DRV[2][3] = K0TT;

/* Evaluate pressure polynomial 
	K equals the secant bulk modulus of seawater
	DK=K(S,T,P)-K(35,0,P)
	K35=K(35,0,P) */
    DK = (B*P + A)*P + K0;
    K35  = (5.03217E-5*P+3.359406)*P+21582.27;
    GAM=P/K35;
    PK = 1.0 - GAM;
    SVA = SVA*PK + (V350P+SVA)*P*DK/(K35*(K35+DK));

/* Scale specific volume anomaly to normally reported units */
    SVAN=SVA*1.0E+8;
    eos8d=SVAN;
    V350P = V350P*PK;

/* Compute density anomaly wrt 1000.0 kg/m**3
	1. DR350: density anomaly at 35 (IPSS-78), O deg C and 0 decibars
	2. DR35P: density anomaly at 35 (IPSS-78), O deg C, pressure variation
	3. DVAN:  density anomaly variations involving specific volume anomaly

   Check value: SIGMA = 59.82037 kg/m**3 FOR S = 40 (IPSS-78),
	T = 40 deg C, P0= 10000 decibars */
    DR35P=GAM/V350P;
    DVAN=SVA/(V350P*(V350P+SVA));
    SIGMA=DR350+DR35P-DVAN;
    DRV[0][2]=SIGMA;
    K=K35+DK;
    VP=1.0-P/K;
    KT = (BT*P + AT)*P + K0T;
    KTT = (BTT*P + ATT)*P + K0TT;
    V = 1.0/(SIGMA+1000.0e0);
    DRV[0][0] = V;
    V2=V*V;

/* Derivative specific volume wrt SALT */
    VS=V0S*VP+V0*P*KS/(K*K);
    RHOS=-VS/V2;
    DRV[2][7]=VS;
    DRV[0][7]=RHOS;
    VT = V0T*VP + V0*P*KT/(K*K);
    VTT = V0TT*VP+P*(2.0*V0T*KT+KTT*V0-2.0*KT*KT*V0/K)/(K*K);
    R0TT=(2.0*VT*VT/V-VTT)/V2;
    DRV[2][2]=R0TT;
    DRV[1][0] = VT;
    DRV[2][0] = VTT;
    RHOT=-VT/V2;
    DRDT=RHOT;
    DRV[1][2]=RHOT;

/* Pressure derivative DVDP 
  	Set A and B to unbiased values */
    A=DRV[0][4];
    B=DRV[0][5];
    DKDP = 2.0*B*P + A;

/* Correct DVDP to per decibar by multiplying *.1 */
    DVDP = -.1*V0*(1.0 - P*DKDP/K)/K;
    DRV[0][6] = -DVDP/V2;
    DRV[1][6] = K;
    DRV[2][6] = DVDP;

    return eos8d;
}

/********************************************************/
/********************************************************/

/* Calling functions to make above routine more user friendly */


/********************************************************/
double hb_alpha(double S,double T,double P)
/********************************************************/
{
    double **DRV;	
    double SV, alpha;
    int i;

    DRV = (double **) calloc(3, sizeof(double *));
    for (i=0; i< 3; ++i) {
       DRV[i] = (double *) calloc(8, sizeof(double));
    }
    SV = hb_eos80d(S,T,P,DRV);      
   alpha = -1.0E7 * DRV[1][2] / (DRV[0][2] +1000.);
 
   for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
   free(DRV);
   return (alpha);
}


/********************************************************/
double hb_beta(double S, double T, double P)
/********************************************************/
/* Salinity contraction coefficient
	hb_beta returns density derivative with salinity scaled by density */
{
    double **DRV;  
    double SV, beta;
    int i;

    DRV = (double **) calloc(3, sizeof(double *));
    for (i=0; i< 3; ++i) {
       DRV[i] = (double *) calloc(8, sizeof(double));
    }
    SV = hb_eos80d(S,T,P,DRV);
    beta = 1.0E7*DRV[0][7]/(DRV[0][2]+1000.);
    
   for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
   free(DRV);
    return beta;
}


/********************************************************/
double hb_gamma(double S, double T, double P)
/********************************************************/
/* Isothermal compressibility coefficient
   hb_gamma returns density derivative with pressure scaled by density */
{
    double **DRV;   
    double SV, gamma;
    int i;

    DRV = (double **) calloc(3, sizeof(double *));
    for (i=0; i< 3; ++i) {
       DRV[i] = (double *) calloc(8, sizeof(double));
    }

    SV = hb_eos80d(S,T,P,DRV);
    gamma =  1.0E7*DRV[0][6]/(DRV[0][2]+1000.);
    for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
    free(DRV);
    return gamma;
}


/********************************************************/
double hb_ratio(double S, double T, double P)
/********************************************************/
/* Density Ratio:  - drho/dT / drho/dS = -dS/dT
   hb_ratio returns ratio of salinity to temperature change per unit density */
{
    double **DRV;   
    double SV, ratio;
    int i;

    DRV = (double **) calloc(3, sizeof(double *));
    for (i=0; i< 3; ++i) {
       DRV[i] = (double *) calloc(8, sizeof(double));
    }

    SV = hb_eos80d(S,T,P,DRV);
    ratio = (-DRV[1][2] / DRV[0][7]);
    for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
    free(DRV);
    return ratio;
}
