

/* prop_subs.c
................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
        author:  Ruth Curry
                 Woods Hole Oceanographic Inst
                 1993
                 Updated to ANSI standards:  Nov 1999
		 Updated 2008 to accommodate std dev & nobs variables
		 Updated 2009 to accommodate Quality flags
		 updated 2011 for t90 variables
................................................................................
.
.  Functions which facilitate use of properties listed in "hydrobase.h"
.
.  void print_prop_menu() : 
.           prints a list of available properties on the stderr device.
.
.  int get_prop_indx(char *s) :
.           translates character string into an integer
.           representing an index to the list of properties.
.
.  double **get_prop_array(struct HYDRO_DATA *data, int i) :
.           Returns address of array associated with property index.          
.
.  char *get_prop_mne(int i)  :
.           returns the mnemonic associated with property index i.
.
.  char *get_prop_descrip(int i) :
.           Returns a char description of property[i].
.
.  char *get_prop_units(int i) :
.           Returns the units associated with property[i].
.
.  int get_field_width(int i) :
.           Returns the field width associated with property[i].
.
.  int get_field_precis(int i) :
.           Returns the decimal precision associated with property[i].
.
.  void compute_sigma(double pref, int nobs, double *sigma, double *p, 
.                     double *t, double *s):
.           Computes an array of sigma values for each level in p,t,s arrays.
.
.  void compute_ratio( int nobs, double *ratio, double *p, double *t, double *s)
.  	 Computes density ratio = -ds/dt from  arrays of p, t, s at each observation
.    	level.     
.	
.  void compute_svan(int nobs, double *sva, double *p, 
.                      double *t, double *s).
.           Computes an array of specific volume anomaly values for each level 
.           in p,t,s arrays.
.
.  void compute_sp_vol( int nobs, double *sv, double *p,  double *t, double *s)
.           Computes an array of spec. volume values (in situ) for each level  
.            in p,t,s arrays.
.
.  void compute_height(int nobs, double *p, double *t, double *s, double pref, double *ht) :
.           Computes dynamic height relative to pref for each level in the
.           p, t, s arrays.
.
.  void compute_energy(int nobs, double *p, double *t, double *s, double pref, double *chi)
.           Computes potential energy relative to reference level pref for each
.           level in the p,t,s arrays.
.
.  void compute_theta(int nobs, double *th, double *p, double *t, double *s);
.           Computes potential temp rel to the sea surface for each level 
.           in the p, t, s arrays.
.
. void buoy_freq(double *bf, double *p, double *t, double *s, 
.                 int  nobs,  int window, int w_incr);
.           Computes brunt-vaisala frequency in radians/sec at each pressure 
.           level, for a set of p,t,s observations. Window and w_incr together 
.           specify the pressure interval and how that interval is divided up
.           for computing the specific vol anomaly gradient near each observed
.           pressure. 
.
. void po_vort(double *pv, double *e, int nobs, double lat)
.           Computes potential vorticity for each value of the stability 
.           parameter (n^2 in (radians/sec)^2) and latitude.
.
.  void compute_sound_vel( double *sv, double *p,  double *t, double *s, int nobs)
.           Computes an array of sound velocity values (in m/sec) for each level  
.            in p,t,s arrays.
.
. double buoyancy(double p0, double *p, double *t, double *s, 
.                 int  nobs,  int window, int w_incr);
.           Computes buoyancy one value at a time.
.           Computes brunt-vaisala frequency in radians/sec at pressure, 
.           p0 given a set of p,t,s observations. The window specifies 
.           the pressure interval surrounding p0 in which observations 
.           will be considered; w_incr is the increment by which the window
.           is subdivided into smaller pieces to approximate the gradient.
.
. double potvort(double e, double lat)
.           Computes PV one value at a time.
.           Returns potential vorticity for the specified value of stability 
.           parameter,e, (= N^2 in (radians/sec)^2) and latitude.
.
. double ox_kg2l(double oxkg, double pr, double te, double sa)
.           Returns oxygen in ml/l from ox in umole/kg.
.
. double ox_l2kg(double oxl, double pr, double te, double sa)
.     Returns oxygen in umole/kg from ox in ml/l.
.    
.  void t68_to_t90(double *t68, double *t90, int nobs)
.    Converts array of temperatures from IPTS-68 to ITS-90 
.
.  void t90_to_t68(double *t90, double *t68, int nobs)
.    Converts array of temperatures from ITS-90 to IPTS-68
.
*/
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#define PROP_SUBS 1
#include "hydrobase.h"

/*************************************************************************/
/*************************************************************************/

/* character string mnemonics which correspond to properties */

char *prop_mne[MAXPROP] =  { "pr",   
                             "de",  
                             "te",
                             "th",   
                             "t90",
                             "th9",   
                             "sa", 
                             "ox", 
                             "o2", 
                             "n2",
                             "n3",
                             "p4",
                             "si", 
                             "ht",  
                             "pe",  
                             "s0",  
                             "s1",  
                             "s2", 
                             "s3", 
                             "s4",
                             "s_",
                             "bf",  
                             "pv",
                             "sv",
                             "va",   
                             "f1",
                             "f2",
                             "f3",
                             "he",
                             "tu",
                             "gn",
                             "ge",
                             "vn",
                             "ve",
			     "vs",
			     "dr",
			     "al",
			     "be",
			     "cl4",
			     "sf6",
			     "co2",
			     "alk",
			     "doc",
			     "tdn",
			     "c13",
			     "c14"
};

char *prop_descrip[MAXPROP] = { "pressure",
                                "depth",
                                "temperature IPTS-68",
                                "potential temperature IPTS-68: pref= 0.",
                                "temperature ITS-90",
                                "potential temperature ITS-90: pref= 0.",
                                "salinity",
                                "oxygen",
                                "oxygen",
                                "nitrite",
                                "nitrate",
                                "phosphate",
                                "silicate",
                                "dynamic height",
                                "potential energy anomaly",
                                "potential density: pref = 0.",
                                "potential density: pref = 1000.",
                                "potential density: pref = 2000.",
                                "potential density: pref = 3000.",
                                "potential density: pref = 4000.",
                                "potential density: pref = ?.",
                                "buoyancy frequency",
                                "potential vorticity",
                                "specific volume",
                                "specific volume anomaly",
                                "cfc-11",
                                "cfc-12",
                                "cfc-113",
                                "helium",
                                "tritium",
                                "neutral density (gamma-n)",
                                "gamma-n errorbar",
				"velocity north",
				"velocity east",
				"sound velocity",
				"density ratio",
				"thermal expansion",
				"haline contraction",
				"carbon tetrachloride",
				"sulfur hexafluoride",
				"total carbon dioxide",
				"total alkalinity",
				"dissolved organic carbon",
				"total dissolved nitrogen",
				"radiocarbon-14",
				"radiocarbon-13"
								
} ;

char *prop_units[MAXPROP] = { "dbars",
                              "meters",
                              "degrees C",
                              "degrees C",
                              "degrees C",
                              "degrees C",
                              "psu",
                              "ml/liter",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "dyn. meters (= *10 m**2/s**2)",
                              "10**6 g/sec**2 (ergs/cm**2)",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "kg/m**3",
                              "* 1.e-5 radians/sec",
                              "* 1.e-12 m^-1 sec^-1",
                              "* 1.e-8 m**3/kg",
                              "* 1.e-8 m**3/kg",
                              "picomole/kg",
                              "picomole/kg",
                              "picomole/kg",
                              "micromole/kg",
                              "tritium units",
                              "kg/m**3",
                              "kg/m**3",
                              "m/sec",
                              "m/sec",
                              "m/sec",
			      "",
			      "10**7 alpha",
			      "10**7 beta",
                              "picomole/kg",
			      "femtomoles/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "micromole/kg",
                              "/mille",
			      "/mille"			      
};

int field_width[MAXPROP] = {     
                                 8, /* pr */
                                 8, /* de */
                                 8, /* te */
                                 8, /* th */
                                 8, /* t90 */
                                 8, /* th9 */
                                 8, /* sa */
                                 8, /* ox */
                                 8, /* o2 */
                                 8, /* n2 */
                                 8, /* n3 */
                                 8, /* p4 */
                                 8, /* si */
                                 8, /* ht */
                                10, /* pe */
                                 8, /* s0 */
                                 8, /* s1 */
                                 8, /* s2 */
                                 8, /* s3 */
                                 8, /* s4 */
                                 8, /* s_ */
                                10, /* bf */
                                10, /* pv */
                                10, /* sv */
                                 8, /* va */
                                 8, /* f1 */
                                 8, /* f2 */
                                 8, /* f3 */
                                 8, /* he */
                                 8, /* tu */
                                 8, /* gn */
                                 8, /* ge */
                                 9, /* vn */
                                 9, /* ve */
                                 9, /* vs */
				 8, /* dr */
				 9, /* al*/
				 9, /* be */
				 8, /* CL4 */
				 8, /* sf6 */
				 8, /* co2 */
				 8, /* alk */
				 8, /* doc */
				 8, /* tdn */
				 8, /* c13 */
				 8  /* c14 */
};

int field_precis[MAXPROP] =  {   
                                 1, /* pr */
                                 1, /* de */
                                 4, /* te */
                                 4, /* th */
                                 4, /* t90 */
                                 4, /* th9 */
                                 4, /* sa */
                                 3, /* ox */
                                 1, /* o2 */
                                 3, /* n2 */
                                 3, /* n3 */
                                 3, /* p4 */
                                 3, /* si */
                                 3, /* ht */
                                 1, /* pe */
                                 4, /* s0 */
                                 4, /* s1 */
                                 4, /* s2 */
                                 4, /* s3 */
                                 4, /* s4 */
                                 4, /* s_ */
                                 4, /* bf */
                                 3, /* pv */
                                 3, /* sv */
                                 4, /* va */
                                 3, /* f1 */
                                 3, /* f2 */
                                 3, /* f3 */
                                 3, /* he */
                                 3, /* tu */
                                 4, /* gn */
                                 5, /* ge */
				 4, /* vn */
				 4, /* ve */
				 3, /* vs */
				 4, /* dr */
				 1, /* al */
				 1, /* be */
				 3, /* CL4 */
				 3, /* sf6 */
				 1, /* co2 */
				 1, /* alk */
				 2, /* doc */
				 2, /* tdn */
				 1, /* c14 */
				 2 /* c13 */				 
};
char *prop_mne_N[MAXPROP] =  { "Npr",   
                             "Nde",  
                             "Nte",
                             "Nth",   
                             "Nt90",
                             "Nth9",   
                             "Nsa", 
                             "Nox", 
                             "No2", 
                             "Nn2",
                             "Nn3",
                             "Np4",
                             "Nsi", 
                             "Nht",  
                             "Npe",  
                             "Ns0",  
                             "Ns1",  
                             "Ns2", 
                             "Ns3", 
                             "Ns4",
                             "Ns_",
                             "Nbf",  
                             "Npv",
                             "Nsv",
                             "Nva",   
                             "Nf1",
                             "Nf2",
                             "Nf3",
                             "Nhe",
                             "Ntu",
                             "Ngn",
                             "Nge",
                             "Nvn",
                             "Nve",
			     "Nvs",
			     "Ndr",
			     "Nal",
			     "Nbe",
			     "Ncl4",
			     "Nsf6",
			     "Nco2",
			     "Nalk",
			     "Ndoc",
			     "Ntdn",
			     "Nc13",
			     "Nc14"
};
char *prop_descrip_N[MAXPROP] = { "Number of pressure obs",
                                "Number of depth obs",
                                "Number of in situ temperature obs IPTS-68",
                                "Number of potential temperature obs IPTS-68: pref= 0.",
                                "Number of in situ temperature obs ITS-90",
                                "Number of potential temperature obs ITS-90: pref= 0.",
                                "Number of salinity obs",
                                "Number of oxygen obs",
                                "Number of oxygen obs",
                                "Number of nitrite obs",
                                "Number of nitrate obs",
                                "Number of phosphate obs",
                                "Number of silicate obs",
                                "Number of dynamic height obs",
                                "Number of potential energy anomaly obs",
                                "Number of potential density obs: pref = 0.",
                                "Number of potential density obs: pref = 1000.",
                                "Number of potential density obs: pref = 2000.",
                                "Number of potential density obs: pref = 3000.",
                                "Number of potential density obs: pref = 4000.",
                                "Number of potential density obs: pref = ?.",
                                "Number of buoyancy frequency obs",
                                "Number of potential vorticity obs",
                                "Number of specific volume obs",
                                "Number of specific volume anomaly obs",
                                "Number of cfc-11 obs",
                                "Number of cfc-12 obs",
                                "Number of cfc-113 obs",
                                "Number of helium obs",
                                "Number of tritium obs",
                                "Number of neutral density (gamma-n) obs",
                                "Number of gamma-n errorbar obs",
				"Number of velocity north obs",
				"Number of velocity east obs",
				"Number of sound velocity obs",
				"Number of density ratio obs",
				"Number of thermal expansion obs",
				"Number of haline contraction obs",
				"Number of carbon tetrachloride obs",
				"Number of sulfur hexafluoride obs",
				"Number of total carbon dioxide obs",
				"Number of total alkalinity obs",
				"Number of dissolved organic carbon obs",
				"Number of total dissolved nitrogen obs",
				"Number of radiocarbon-13 obs",
				"Number of radiocarbon-14 obs",
								
} ;
char *prop_mne_V[MAXPROP] =  { "Vpr",   
                             "Vde",  
                             "Vte",
                             "Vth",   
                             "Vt90",
                             "Vth9",   
                             "Vsa", 
                             "Vox", 
                             "Vo2", 
                             "Vn2",
                             "Vn3",
                             "Vp4",
                             "Vsi", 
                             "Vht",  
                             "Vpe",  
                             "Vs0",  
                             "Vs1",  
                             "Vs2", 
                             "Vs3", 
                             "Vs4",
                             "Vs_",
                             "Vbf",  
                             "Vpv",
                             "Vsv",
                             "Vva",   
                             "Vf1",
                             "Vf2",
                             "Vf3",
                             "Vhe",
                             "Vtu",
                             "Vgn",
                             "Vge",
                             "Vvn",
                             "Vve",
			     "Vvs",
			     "Vdr",
			     "Val",
			     "Vbe",
			     "Vcl4",
			     "Vsf6",
			     "Vco2",
			     "Valk",
			     "Vdoc",
			     "Vtdn",
			     "Vc13",
			     "Vc14"
};

char *prop_descrip_V[MAXPROP] = { "variance of pressure",
                                "variance of depth",
                                "variance of in situ temperature IPTS-68",
                                "variance of potential temperature IPTS-68: pref= 0.",
                                "variance of in situ temperature ITS-90",
                                "variance of potential temperature ITS-90",
                                "variance of salinity",
                                "variance of oxygen",
                                "variance of oxygen",
                                "variance of nitrite",
                                "variance of nitrate",
                                "variance of phosphate",
                                "variance of silicate",
                                "variance of dynamic height",
                                "variance of potential energy anomaly",
                                "variance of potential density: pref = 0.",
                                "variance of potential density: pref = 1000.",
                                "variance of potential density: pref = 2000.",
                                "variance of potential density: pref = 3000.",
                                "variance of potential density: pref = 4000.",
                                "variance of potential density: pref = ?.",
                                "variance of buoyancy frequency",
                                "variance of potential vorticity",
                                "variance of specific volume",
                                "variance of specific volume anomaly",
                                "variance of cfc-11",
                                "variance of cfc-12",
                                "variance of cfc-113",
                                "variance of helium",
                                "variance of tritium",
                                "variance of neutral density (gamma-n)",
                                "variance of gamma-n errorbar",
				"variance of velocity north",
				"variance of velocity east",
				"variance of sound velocity",
				"variance of density ratio",
				"variance of thermal expansion",
				"variance of haline contraction",							"variance of carbon tetrachloride",		
				"variance of sulfur hexafluoride",		
				"variance of total carbon dioxide",
				"variance of total alkalinity",	
				"variance of dissolved organic carbon",		
				"variance of total dissolved nitrogen",		
				"variance of radiocarbon-13",			
				"variance of radiocarbon-14"				} ;

char *prop_mne_Q[MAXPROP] =  { "Qpr",   
                             "Qde",  
                             "Qte",
                             "Qth",   
                             "Qt90",
                             "Qth9",   
                             "Qsa", 
                             "Qox", 
                             "Qo2", 
                             "Qn2",
                             "Qn3",
                             "Qp4",
                             "Qsi", 
                             "Qht",  
                             "Qpe",  
                             "Qs0",  
                             "Qs1",  
                             "Qs2", 
                             "Qs3", 
                             "Qs4",
                             "Qs_",
                             "Qbf",  
                             "Qpv",
                             "Qsv",
                             "Qva",   
                             "Qf1",
                             "Qf2",
                             "Qf3",
                             "Qhe",
                             "Qtu",
                             "Qgn",
                             "Qge",
                             "Qvn",
                             "Qve",
			     "Qvs",
			     "Qdr",
			     "Qal",
			     "Qbe",
			     "Qcl4",
			     "Qsf6",
			     "Qco2",
			     "Qalk",
			     "Qdoc",
			     "Qtdn",
			     "Qc13",
			     "Qc14"
};


char *prop_descrip_Q[MAXPROP] = { "pressure quality flag",
                                "depth quality flag",
                                "in situ temperature IPTS-68 quality flag",
                                "potential temperature IPTS-68: pref= 0. quality flag",
                                "in situ temperature ITS-90 quality flag",
                                "potential temperature ITS-90 quality flag",
                                "salinity quality flag",
                                "oxygen quality flag",
                                "oxygen quality flag",
                                "nitrite quality flag",
                                "nitrate quality flag",
                                "phosphate quality flag",
                                "silicate quality flag",
                                "dynamic height quality flag",
                                "potential energy anomaly quality flag",
                                "potential density: pref = 0. quality flag",
                                "potential density: pref = 1000. quality flag",
                                "potential density: pref = 2000. quality flag",
                                "potential density: pref = 3000. quality flag",
                                "potential density: pref = 4000. quality flag",
                                "potential density: pref = ?. quality flag",
                                "buoyancy frequency quality flag",
                                "potential vorticity quality flag",
                                "specific volume quality flag",
                                "specific volume anomaly quality flag",
                                "cfc-11 quality flag",
                                "cfc-12 quality flag",
                                "cfc-113 quality flag",
                                "helium quality flag",
                                "tritium quality flag",
                                "neutral density (gamma-n) quality flag",
                                "gamma-n errorbar quality flag",
				"velocity north quality flag",
				"velocity east quality flag",
				"sound velocity quality flag",
				"density ratio quality flag",
				"thermal expansion quality flag",
				"haline contraction quality flag",
				"carbon tetrachloride quality flag",
				"sulfur hexafluoride quality flag",
				"total carbon dioxide quality flag",
				"total alkalinity quality flag",
				"dissolved organic carbon quality flag",
				"total dissolved nitrogen quality flag",
				"radiocarbon-13 quality flag",
				"radiocarbon-14 quality flag"				} ;

/******************************************************/

void print_prop_menu()
{
   int i;

   for (i=0; i < MAXPROP; ++i) {
     /* skip properties dealing with quality and/or statistics (they begin with
     capital letters */
     if ( prop_mne[i][0] >= 'a' && prop_mne[i][0] <= 'z') {
        fprintf(stderr,"\n%s : %s [%s]", prop_mne[i], prop_descrip[i], prop_units[i]);
     }
   }
   fprintf(stderr,"\n\nCapital letter in front of observed property mnemonic :");
   fprintf(stderr,"\n  'N' for number of obs,");
   fprintf(stderr,"\n  'V' for standard deviation,");
   fprintf(stderr,"\n  'Q' for quality flag\n");
   
   return;
} /* end print_prop_menu() */

/******************************************************/

int get_prop_indx(char *str)
{
   int error = -1;
   int i;
   char *s;

   s = str;
   switch (*s) {
      case 'a':
            switch (*(++s)) {
               case 'l':
	           if (*(s+1) == 'k')
		      return((int)ALK);
                   return ((int) AL);  
               default:      
                   return (error);
             } 
            break;
      
      case 'b':
            switch (*(++s)) {
               case 'e':
                   return ((int) BE);  
               case 'f':
                   return ((int) BF);  
               default:      
                   return (error);
             } 
            break;
      case 'c':
            switch (*(++s)) {
               case 'l':
                   return ((int) CL4);
               case 'o':
                   return ((int) CO2);
		     
               case '1':
	           if (*(s+1) == '3')
                      return ((int) C13); 
		   return(C14); 
               default:      
                   return (error);
             } 
            break;
      case 'd':
            switch (*(++s)) {
               case 'e':
                   return ((int) DE);  
               case 'o':
                   return ((int) DOC);  
               case 'r':
                   return ((int) DR);  
               default:      
                   return (error);
             } 
            break;
      case 'f':
            switch (*(++s)) {
               case '1':
                   return ((int) F1);
               case '2' :
                   return ((int) F2);  
               case '3' :
                   return ((int) F3);  
               default:      
                   return (error);
             } 
            break;
      case 'g':
            switch (*(++s)) {
               case 'e' :
                   return ((int) GE);  
               case 'n' :
                   return ((int) GN);  
               default:      
                   return (error);
             } 
            break;
      case 'h':
            switch (*(++s)) {
               case 'e':
                   return ((int) HE);  
               case 't':
                   return ((int) HT);  
               default:      
                   return (error);
             } 
            break;
      case 'N':
            if ((i= get_prop_indx(++s)) < 0)
	       return (error);
	    return(i+ 500);
            break;
	    
      case 'n':
            switch (*(++s)) {
              case '2':
                 return ((int) N2);
              case '3':
                 return ((int) N3);
             default:
                 return error;
            }
            break;
      case 'o':
            switch (*(++s)) {
               case 'x':
                   return ((int) OX);  
               case '2':
                   return ((int) O2);  
              default:      
                   return (error);
             } 
            break;
      case 'p':
            switch (*(++s)) {
               case 'e':
                   return ((int) PE);
               case 'r':
                   return ((int) PR);
               case 'v':
                   return ((int) PV);  
               case '4':
                   return ((int) P4);  
               default:      
                   return (error);
             } 
            break;
      case 'Q':
            if ((i= get_prop_indx(++s)) < 0)
	       return (error);
	    return(i+ 900);
            break;
      case 's':
            switch (*(++s)) {
              case 'a':
                 return ((int) SA);
              case 'i':
                 return ((int) SI);
              case 'f':
                 return ((int) SF6);
              case 'v':
                 return ((int) SV);
              case '_':
                 return ((int) S_);
              case '0':
                 return ((int) S0);
              case '1':
                 return ((int) S1);
              case '2':
                 return ((int) S2);
              case '3':
                 return ((int) S3);
              case '4':
                 return ((int) S4);
              default:
                 return error;
            } /* end switch */
            break;
      case 't':
            switch (*(++s)) {
	      case '9': 
	         return ((int) T90);
              case 'd':
                 return ((int) TDN);
              case 'e':
                 return ((int) TE);
              case 'h':
	         if (*(s+1) == '9')
		    return((int)TH9);
		    
                 return ((int) TH);
              case 'u':
                 return ((int) TU);
     
              default:
                 return error;
            } /* end switch */
            break;
      case 'V':
            if ((i= get_prop_indx(++s)) < 0)
	       return (error);
	    return(i+ 700);
            break;
      case 'v':
            switch (*(++s)) {
              case 'a':
                 return ((int) VA);
               case 'e' :
                   return ((int) VE);  
               case 'n' :
                   return ((int) VN);  
               case 's' :
                   return ((int) VS);  
              default:
                 return error;
            } /* end switch */
            break;

      default:
            return error;
   } /* end switch */
} /* end get_prop_indx() */

/*****************************************************************/
double **get_prop_array(struct HYDRO_DATA *dataptr, int i)
/* Returns pointer to property array associated with index i */
{
   
   if (i < 100) 
      return(&dataptr->observ[i]);
      
   
   if (i < 600)
      return(&dataptr->count[i-500]);
      
   if (i < 800)
       return(&dataptr->variance[i-700]);
   
   if (i >= 900 && i < 1000)
       return(&dataptr->quality[i-900]);
   
   fprintf(stderr,"unknown property ID passed to get_prop_array()\n");
   return(NULL);
   
   
}  /* end get_prop_array() */
/*****************************************************************/
char *get_prop_mne(int i)
/*  returns the property mnemonic associated with index i */
{
   
   if (i < 100)
      return(prop_mne[i]);
      
   if (i >= 500 && i < 600)
      return(prop_mne_N[i-500]);
      
   if (i >= 700 && i < 800)
      return(prop_mne_V[i-700]);
      
   if (i >= 900)
      return(prop_mne_Q[i-900]);
      
   fprintf(stderr,"Unknown property index\n");
   return(NULL);
      
}
/*****************************************************************/
char *get_prop_descrip(int i)
   /* returns the ith property description */
{
   if (i < 100)
      return prop_descrip[i];
      
   if (i >= 500 && i < 600)
      return(prop_descrip_N[i-500]);
      
   if (i >= 700 && i < 800)
      return(prop_descrip_V[i-700]);
      
   if (i >= 900)
      return(prop_descrip_Q[i-900]);
      
   fprintf(stderr,"Unknown property index\n");
   return(NULL);
      
}
/*****************************************************************/
int get_field_width(int i)
   /* returns the ith property field width */
{
   if (i < 100)                /* observ */
      return (field_width[i]);
      
   if (i >= 500 && i < 600)   /* count */
      return (6);
      
   if (i >= 700 && i < 800)   /* variance */
      return (field_width[i-700]);
   
   if (i >= 900)              /* quality  */
      return (2);
    
   fprintf(stderr,"Unknown property index\n");
   return(-1);
}
/*****************************************************************/
int get_field_precis(int i)
   /* returns the decimal precision used to describe property i */
{
   if (i < 100)                  /* observ */
      return field_precis[i];
   if (i >= 700 && i < 800)     /* variance */
      return field_precis[4];
      
    return(0);                  /* count & quality  */
}
/*****************************************************************/
char *get_prop_units(int i)
    /* returns the units for property[i] */
{
   if (i < 100)
      return prop_units[i];
   if (i >= 700 && i < 800)
      return prop_units[i-700];
   return("");   
}
/*****************************************************************/
void compute_sigma(double pref, int nobs, double *sigma, double *p, double *t, double *s)
  /* Computes sigma values from  arrays of p, t, s at each observation
    level.     */
{
   int    j;
   double tref, sv;

      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    sigma[j] = HB_MISSING;
	 else {
            tref = hb_theta(s[j], t[j], p[j], pref);
            sv = hb_svan(s[j], tref, pref, &sigma[j]);
	 }
      }

      return;
}  /* end compute_sigma() */
/*****************************************************************/
void compute_ratio( int nobs, double *ratio, double *p, double *t, double *s, double *alpha, double *beta)
  /* Computes density ratio = -ds/dt, thermal expansion (alpha) and haline contraction (beta) coefficients  from  arrays of p, t, s at each observation level. Set ratio, alpha, or beta to NULL if either of these quantities are not returned   */
{
   int    i, j;
   double **DRV;
   double SV;
   
   DRV = (double **) calloc(3, sizeof(double *));
   for (i= 0; i < 3; ++i)
       DRV[i] = (double *) calloc(8, sizeof(double));
       
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0) {
	    if (alpha != NULL)
	         alpha[j] = HB_MISSING;
	    if (beta != NULL)	 
	         beta[j] = HB_MISSING;
	    if (ratio != NULL)	 
	        ratio[j] = HB_MISSING;
	 }
	 else {
	 
	    SV = hb_eos80d(s[j], t[j], p[j], DRV);
	    
	    if (alpha != NULL)
	         alpha[j] = -1.0e7 * DRV[1][2] / (DRV[0][2] + 1000.0);
	    if (beta != NULL)	 
	         beta[j] = 1.0e7 * DRV[0][7] / (DRV[0][2] + 1000.0);
	    if (ratio != NULL)	 
	        ratio[j] = -DRV[1][2] / DRV[0][7]  ;
	 }
      }

   for (i = 0; i < 3; ++i)       
       free( DRV[i]); 
   free(DRV);
   return;
}  /* end compute_ratio() */

/****************************************************************************/
void compute_sp_vol(int nobs, double *spvol, double *p, double *t, double *s)
  /*  Computes in situ Specific Volume from p, t, s at each observation
    level.     */
{
   int    j;
   double  sigma, sv;

   /* Units are  * 1e-8 m**3/kg 
      Convert the sigma value (kg/m**3) to a rho value by adding 1000.
      Specific volume is 1/rho, but is multiplied here by 10^8 to 
      convert to same units as sp.vol. anomaly. */
 
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    spvol[j] = HB_MISSING;
	 else {
           sv = hb_svan(s[j], t[j], p[j], &sigma);
           spvol[j] = 1.0e8 / (sigma + 1000.);
	 }
      }

      return;

}  /* end compute_sp_vol() */
/****************************************************************************/
void compute_svan(int nobs, double *sva, double *p, double *t, double *s)
  /* Computes specific volume anomaly ( = spvol(p,t,s) - spvol(p,0,35)
     from  arrays of p, t, s at each observation level.  
     Units are in 1.0e-8 m**3/kg 
               or 1.0e-5 cm**3/g  */
{
   int    j;
   double sigma;

      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    sva[j] = HB_MISSING;
	 else
            sva[j] = hb_svan(s[j], t[j], p[j], &sigma);
      }

      return;
}  /* end compute_svan() */
/****************************************************************************/
void compute_height(int nobs, double *p, double *t, double *s, double pref, double *h)

/* computes dynamic height relative to pref at each level in pressure
   array.  The units are :  
             dynamic height  : dyn meters = 1/10 * m**2/s**2 .  
           spec vol anomaly  :  1e-8 * m**3/kg
                   pressure  :  dbars = 1e4 N/m**2 = 1e4 kg/m s**2

   If a vertical datagap is encountered, 
    no height is computed beneath that level
 */
{
  int j, start, i, datagap;
  double sref, tref;
  double sva1, sva0, sig, last_h, last_p;
  

/* first find the reference pressure level ... */

   if ((tref = hb_linterp(pref, p, t, nobs)) < -8.) {
      for (i = 0; i < nobs; ++i) {
         h[i] = -99999.9;
      }
      return;
   }

   sref = hb_linterp(pref, p, s, nobs);

   start = 0;
   while (p[start] < pref)
      ++start;
      
   if (start > 0)
      --start;  /* start now points to the first pr level above the ref pr */

   /* check for datagap around the ref level */
      
   datagap = (p[start+1] - p[start]) > 600 ;
   if (datagap) {
     for (i = 0; i < nobs; ++i) {
        h[i] = -99999.9;   /* can't compute height for this station */
     }
     return;
   }
   

/* determine height between the ref lev and the first pr level above...*/   
   
   sva0 = hb_svan(sref, tref, pref, &sig);
   h[start] = 0.0;
   if (start >= 0) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      h[start] = (sva0 + sva1) *.5e-5 * (pref - p[start]);
      sva0 = sva1;
   }
   last_h = h[start];
   last_p = p[start];

/* now integrate upward through the station, check for missing values in p,t,s
   and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         h[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = ( last_p - p[j]) > 250;
      else
        datagap = ( last_p - p[j]) > 610;

      if (datagap) {
        for (i= j; i >= 0; --i) {
          h[i] = -99999.9;
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        sva1 = hb_svan(s[j], t[j], p[j], &sig);

/* the 1e-5 term corrects the units : 10^-8  * 10^4 * 10^-1 
                                      (svan)   Pa/db   dyn meters   */
        h[j] = last_h + ((sva0 + sva1) * 0.5e-5 * (last_p - p[j]));
        last_h = h[j];
        last_p = p[j];
        sva0 = sva1;
      }
   } /* end for */

   
 /* find ht between ref level and the first observation beneath that level... */
  
   ++start;   /* start now points to first level at or below ref pr */
   sva0 = hb_svan(sref, tref, pref, &sig);
   h[start] = 0.0;
   if (start < nobs) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      h[start] = (sva0 + sva1) * 0.5e-5 * (pref - p[start]) ;
      sva0 = sva1;
   }
   last_h = h[start];
   last_p = p[start];
  
/* now integrate downward through the station, check for missing values 
    in p,t,s and for vertical datagaps ... */
    
   for (j = start+1; j < nobs; ++j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         h[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = (  p[j] - last_p) > 250;
      else
        datagap = (  p[j] - last_p) > 610;

      if (datagap) {
        for (i= j; i < nobs; ++i) {
          h[i] = -99999.9;
        }
        return;  /* don't bother integrating downward any farther */
      }
      
      sva1 = hb_svan(s[j], t[j], p[j], &sig);

      h[j] = last_h + ((sva0 + sva1) *.5e-5 * (last_p - p[j]) );
      last_h = h[j];
      last_p = p[j];
      sva0 = sva1;
   }

   return;
      
} /* end compute_height() */

/****************************************************************************/
void compute_energy(int nobs, double *p, double *t, double *s, double pref, double *chi)

/* computes potential energy anomaly relative to the specified ref level at 
   each level in pressure array.  
        chi = 1/g * Integral of (Pr * sp_vol_anom * delta-P)
   The units are :  
           potential energy  : * 1e6 ergs/cm**2 (or g/sec**2) = 1e3 J/m**2
           spec vol anomaly  : * 1e-5 cm**3/g
                   pressure  :  dbars = 1e5 dyn/cm**2 = 1e5 g/cm/s**2
                          g  : 980 cm/s**2
                          
          
      value of constant comes from :
              correcting spvolanom units: 1e-5
              getting mean spvolanom:    .5
              correcting delta-pr units: 1e5
              getting mean pr:           .5
              correcting pr units:       1e5
              1/g:                       .00102040816
              adjust decimal place:      1e-6

    If a vertical datagap is encountered, the integration is stopped and
    no potential energy is computed above (beneath) that level.
 */
{
   int j, start, i, datagap;
   double sref, tref;
   double sva1, sva0, sig, last_chi, last_p;
   double C = .255102e-4;  /*constant defined above */

/* first find the reference pressure level ... */

   if ((tref = hb_linterp(pref, p, t, nobs)) < -8.) {
      for (i = 0; i < nobs; ++i) {
         chi[i] = -99999.9;
      }
      return;
   }
   
   sref = hb_linterp(pref, p, s, nobs);
   
   start = 0;
   while (p[start] < pref)
      ++start;
      
   --start;  /* start now points to the first pr level above the ref pr */

   /* check for datagap around the ref level */
      
   datagap = (p[start+1] - p[start]) > 600 ;
   if (datagap) {
     for (i = 0; i < nobs; ++i) {
        chi[i] = -99999.9;   /* can't compute energy for this station */
     }
     return;
   }
   
/* determine chi between the ref lev and the first pr level above...*/   
   
   sva0 = hb_svan(sref, tref, pref, &sig);
   chi[start] = 0.0;
   if (start >= 0) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      chi[start] = (sva0 + sva1) * (pref - p[start]) * (pref + p[start]) * C;
      sva0 = sva1;
   }
   last_chi = chi[start];
   last_p = p[start];

/* now integrate upward through the station, check for missing values in p,t,s
   and for vertical datagaps ... */

   for (j = start-1; j >= 0; --j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         chi[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = ( last_p - p[j]) > 250;
      else
        datagap = ( last_p - p[j]) > 600;

      if (datagap) {
        for (i= j; i >= 0; --i) {
          chi[i] = -99999.9;
        }
        j = 0;  /* don't bother integrating upward any farther */
      }
      else {
        sva1 = hb_svan(s[j], t[j], p[j], &sig);

            /* C is a constant defined above */
        chi[j] = last_chi + ((sva0 + sva1) * (last_p - p[j])*(p[j]+ last_p) * C);
        last_chi = chi[j];
        last_p = p[j];
        sva0 = sva1;
      }
   }
   
 /* find chi between ref level and the first observation beneath that level... */
  
   ++start;   /* start now points to first level at or below ref pr */
   sva0 = hb_svan(sref, tref, pref, &sig);
   chi[start] = 0.0;
   if (start < nobs) {
      sva1 = hb_svan(s[start], t[start], p[start], &sig);
      chi[start] = (sva0 + sva1) * (pref - p[start]) * (pref + p[start]) * C;
      sva0 = sva1;
   }
   last_chi = chi[start];
   last_p = p[start];
  
/* now integrate downward through the station, check for missing values 
    in p,t,s and for vertical datagaps ... */
    
   for (j = start+1; j < nobs; ++j) {
      if (s[j] < -8. || t[j] < -8. || p[j] < -8) {
         chi[j] = -99999.9;
         continue;
      }

      if (p[j] < 1001.)
        datagap = (  p[j] - last_p) > 250;
      else
        datagap = (  p[j] - last_p) > 600;

      if (datagap) {
        for (i= j; i < nobs; ++i) {
          chi[i] = -99999.9;
        }
        return;  /* don't bother integrating downward any farther */
      }
      
      sva1 = hb_svan(s[j], t[j], p[j], &sig);

            /* C is a constant defined above */
      chi[j] = last_chi + ((sva0 + sva1) * (last_p - p[j])*(p[j]+ last_p) * C);
      last_chi = chi[j];
      last_p = p[j];
      sva0 = sva1;
   }
   return;
      
} /* end compute_energy() */

/***************************************************************************/
void compute_theta(int n, double *th, double *p, double *t, double *s)
/* int n;       # of levels in arrays 
   double *th;  pointer to array of theta values 
   double *p;   pointer to array of pressure 
   double *t;   pointer to array of in situ temperature 
   double *s;   pointer to array of salinity 
*/
{
   int i;
   double pref = 0.0;

   for (i = 0; i < n; ++i) {
       if (s[i] < 0.0 || p[i] < 0.0)
         th[i] = HB_MISSING;
	 
       else if (s[i] > 100000 || p[i] > 100000)  /* simulate a land mask */  
         th[i] = s[i];
       else
         th[i] = hb_theta(s[i], t[i], p[i], pref);
   }
   return;
}
/***************************************************************************/
void buoy_freq(double *bf, double *p, double *t, double *s,int nobs,int window, int w_incr)
/*  double *bf;         array to hold returned buoyancy values
    double *p, *t, *s;  observed pressure, temperature, salinity values 
    int nobs;           number of observations in p,t,s arrays
    int window;         size of  pressure window (db) into which
                        observations will be incorporated for estimating t,s
                        gradients
    int w_incr;         subdivide the window into increments (db) for
                          approximating the gradients 
*/

/* Computes buoyancy frequency in radians/sec using gradients of 
   specific vol anomaly near each observed pressure. A pressure series of 
   p,t,s points is generated at w_incr increments from the shallowest observed
   pressure to the deepest. A sp.vol.anom. gradient is then estimated over 
   the specified pressure window and buoyancy is computed at each
   increment of the pressure series.  These values of buoyancy are interpolated
   to obtain a value at each pressure level of the original observations.
*/
{
   double *p_win, *t_win, *s_win, *b_win, p0;
   double  e, bflast;
   double cph2rps = 0.001745329;   /* (2*pi)/3600 */
   int mid, npts_win, npts;
   int i, j, k, n;

/* determine maximum size of arrays to hold pressure series for a window */

   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   
   e = (p[nobs-1] - p[0]) / (double) w_incr;
   npts = (int) (e + .00001);
   ++npts;
   ++npts; /* add one more  for bottom point */
   
   p_win = (double *) malloc(npts * sizeof(double));
   t_win = (double *) malloc(npts * sizeof(double));
   s_win = (double *) malloc(npts * sizeof(double));
   b_win = (double *) malloc(npts * sizeof(double));
   if (b_win == NULL) {
     fprintf(stderr, "\nUnable to allocate memory in buoy_freq()\n");
     exit(2);
   }

/* set up a pressure series of t and s at the specified increment... */

   p_win[0] = p[0];
   t_win[0] = t[0];
   s_win[0] = s[0];
   --npts;  /* to avoid subtracting everytime test is performed */
   for (i = 1; i < npts; ++i) {
      p_win[i] = p_win[i-1] + w_incr;
      t_win[i] = hb_linterp(p_win[i], p, t, nobs);
      s_win[i] = hb_linterp(p_win[i], p, s, nobs);
   }
        
   p_win[npts] = p[nobs-1];
   t_win[npts] = t[nobs-1];
   s_win[npts] = s[nobs-1];
   
   ++npts;  /* back to correct npts */

/* Compute bvfreq at each point in the pressure series sliding the pressure
   window incrementally through.  Adjust the window at the top 
   of the profile to use all the points that are available.  This will decrease
   to 2 points at the top. Do not do the same at the bottom because it
   leads to unusually high pv values  */

   mid = window / w_incr;  /* index of midpoint of window */
   npts_win = mid * 2 + 1;
   
   n = 2;
   j = mid;
   if (mid > npts)
     j = npts;
   for (i = 0; i < j; ++i) {
     n = i * 2 + 1;
     if (i == 0)
        n = 2;
     if (n >= npts)
        n = npts;
     b_win[i] = hb_bvfrq(&s_win[0], &t_win[0], &p_win[0], n, &p0, &e);
   }
   n = npts_win;
   j = npts - mid;
   k = 0;
   for (i = mid; i < j; ++i) {
        b_win[i] = hb_bvfrq(&s_win[k], &t_win[k], &p_win[k], n, &p0, &e);
        ++k;
   }
   
     /* set bottom of window to last computed b_win value */
   e = b_win[i-1];   
   if (j > mid) {
     for (i = j; i < npts; ++i) 
        b_win[i] = e;
   }
   
   for (i = 0; i < nobs; ++i) {
      bf[i] = hb_linterp(p[i], p_win, b_win, npts);
      if (bf[i] > -999.)
         bf[i] *= cph2rps;   /* correct the units */
   }
 
    
   free((void *)b_win);
   free((void *)p_win);
   free((void *)t_win);
   free((void *)s_win);
   return;

} /* end buoy_freq() */

/****************************************************************************/
void po_vort( double *pv, double *e, int nobs, double lat)
/* Returns potential vorticity for the specified value of stability 
   parameter (n-squared) and latitude. 
    pv;    pointer to array of returned povo values  
    e;     pointer to array of n-squared (buoyancy freq)values in (rad/sec)^2  
    nobs;  dimension of pv and e arrays 
    lat;   latitude in degrees 
*/
{
   int i;
   for (i = 0; i < nobs; ++i) {
      
      if (e[i] < -9990.)       /* check for missing value flag ... */
        pv[i] =  -999.;
      else {
        pv[i] = e[i] * 1.0e14  * hb_coriol(lat) / hb_gravity(lat);
     }
   }
   return ;

} /* end potvort() */


/***************************************************************************/
void compute_sound_vel(double *svel, double *p, double *t, double *s, int nobs)
  /*  Computes sound velocity from p, t, s at each observation
    level.     */
{
   int    j;
 
      for (j = 0; j < nobs; ++j) {
         if (s[j] < 0 || p[j] < 0)
	    svel[j] = HB_MISSING;
	 else 
           svel[j] = hb_svel(s[j], t[j], p[j]);
      }

      return;

}  /* end compute_sp_vol() */
/****************************************************************************/
double buoyancy(double p0, double *p, double *t, double *s, int nobs, int window, int w_incr)

/* Computes a single value of buoyancy frequency (radians/sec) using gradients of 
   t and s near the pressure designated by p0.  The t,s gradients are 
   estimated by linearly interpolating t and s as a function of p from 
   the observations which which fall within a pressure range of window
   on either side of p0.  A pressure series of p,t,s points is generated in
   this manner at w_incr intervals from p0.  If the p0 specified does 
   not exist in the pressure series represented by p, the value -9999 is
   returned.
   
   double p0;          pressure level at which buoyancy is computed 
   double *p, *t, *s;  observed pressure, temperature, salinity values 
   int nobs;           number of observations in p,t,s arrays 
   int window;         size of pressure interval (db) on either side of p0 in which
                      observations will be incorporated into estimating t,s
                      gradients
   int w_incr;         subdivide the window into increments (db) for
                      approximating the gradients 
*/
{
   double *p_win, *t_win, *s_win;
   double bv, e;
   double cph2rps = 0.001745329;   /* (2*pi)/3600 */
   int mid, npts_win;
   int i, n, end, start;

/* determine maximum size of arrays to hold pressure series for a window */

   window = window / 2;
   if (w_incr > window)
       w_incr = window;
   mid = window / w_incr;  /* index of midpoint of array */
   npts_win = mid * 2 + 1;
   
   p_win = (double *) malloc(npts_win * sizeof(double));
   t_win = (double *) malloc(npts_win * sizeof(double));
   s_win = (double *) malloc(npts_win * sizeof(double));

/* determine whether specified pressure level exists at this station */

   if (hb_linterp(p0, p, t, nobs) < -9990.) {
           return (-9999.);
   }

/* set up the window of pressure around the specifed pressure level...
   Account for cases where window juts above top observed pressure or 
   below bottom observed pressure ...  
   start represents the index at the top of window, end is
   the index at the bottom.  Keep the top and bottom as equidistant as
   possible from the midpoint because bvf will be computed at the mid-pressure
   of the window.*/

   p_win[mid] = p0;
   
   start = 0;
   i = mid;
   while (--i >= start) {
     if ((p_win[i] = p_win[i+1] - w_incr) < p[0] )
       start = i+1;
   }
   
   end = mid + mid - start;
   if (end == mid)
      ++end;
      
   i = mid;
   n = nobs-1;
   while (++i <= end) {
      if ((p_win[i] = p_win[i-1] + w_incr) > p[n])
        end = i-1;
   }
   
   if (start < mid)   
      start = mid - (end - mid);
   
   if (end == mid)
      start = mid - 1;
      
   npts_win = (end - start) + 1;

   for (i = start; i <= end; ++i) {
      t_win[i] = hb_linterp(p_win[i], p, t, nobs);
      s_win[i] = hb_linterp(p_win[i], p, s, nobs);
   }
   
   bv = hb_bvfrq(&s_win[start], &t_win[start], &p_win[start], npts_win, &p0, &e);

   free((void *)p_win);
   free((void *)t_win);
   free((void *)s_win);
   return (bv * cph2rps);

} /* end buoyancy() */
/****************************************************************************/
double potvort(double e, double lat)

/* Returns a single value of potential vorticity for the specified value of 
stability  parameter, e, (= buoyancy frequency squared) and latitude.
   
    e :    n-squared (buoyancy freq) in (rad/sec)^2 
   lat:   latitude in degrees
*/
{

   if (e < -9990.)       /* check for missing value flag ... */
       return (-99.);

    return (e * 1.0e14  * hb_coriol(lat) / hb_gravity(lat) );

} /* end potvort() */

/****************************************************************************/
double ox_kg2l(double ox, double pr, double te, double sa)
/*  Returns oxygen in units of ml/l 
   double ox:   oxygen value in micromoles/kg 
   double pr, te, sa:    in situ properties 
*/
{
    double x, pt;
    double sig, zero = 0.0;
    
    if (ox <= 0) 
       return ox;
    
    pt = hb_theta(sa, te, pr, zero);
    x = hb_svan(sa, pt, zero, &sig);
    return (ox * 0.022392 * (1 + sig /1000.));

}
/****************************************************************************/
double ox_l2kg(double ox, double pr, double te, double sa)
/*  Returns oxygen in units of umole/kg 
   double ox:   oxygen value in ml/l 
   double pr, te, sa:    in situ properties 
*/
{
    double x, pt;
    double sig, zero = 0.0;
    
    if (ox <= 0) 
       return ox;
    pt = hb_theta(sa, te, pr, zero);
    x = hb_svan(sa, pt, zero, &sig);
    return (ox / (0.022392 * (1 + sig /1000.)));

}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
void t68_to_t90(double *t68, double *t90, int nobs)

   /*  Converts array of temperatures from IPTS-68 to ITS-90  */
{
  int i;
  
  for (i = 0; i < nobs; ++i)
  
     t90[i] = t68[i] * 0.99976;
     
   return; 
}

/****************************************************************************/
void t90_to_t68(double *t90, double *t68, int nobs)

   /*  Converts array of temperatures from ITS-90 to IPTS-68  */
{
  int i;
  
  for (i = 0; i < nobs; ++i)
  
     t68[i] = t90[i] * 1.00024;
     
   return; 
}
