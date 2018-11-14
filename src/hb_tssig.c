/*   hb_tssig.c
................................................................................
                              *  HydroBase2 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
                             updated July 2000 to ANSI
................................................................................
................................................................................
.  usage:  tssig pref {smin smax sincr tmin tmax tincr} > outputfile.name
................................................................................
.   The min, max, increment
.   of theta and salt which define the grid are optional command line arguments.
.   If these values are not specified, these default values  are used ...
.          salt: 30, 40, .1       theta:  -2, 28, .25
.
.  Reference Pressure, for which all sigma values will be computed,is the first 
.   command line argument.
.
.  The result of this program is a matrix of sigma values with z(1,1) 
.   corresponding to the point (saltmin, thetamin)  or the lower left
.   hand corner of a theta-salt plot.  Sigma is computed from specific
.   volume anomaly based on 1980 equation of state (Millero, et al., 1980).
................................................................................

*/

#include <stdio.h>
#include <stdlib.h>
#include "hydrobase.h"



main(int argc, char **argv)
{
   double th, tmin = -2.0, tmax = 28.0, tincr = .25;
   double s, smin = 32.0, smax = 38.0, sincr = .1;
   double sig, pref, tref, zero = 0.0;

   if (argc < 2 || argv[1][1] == 'h') {
     fprintf(stderr, "\nUSAGE:  %s pref smin smax sincr th_min th_max th_incr", argv[0]);
     fprintf(stderr, " > outfile\n");
     exit(1);
   }

   pref = atof(argv[1]);
   if (argc > 2) {
      smin = atof(argv[2]);
      smax = atof(argv[3]);
      sincr = atof(argv[4]);
      tmin = atof(argv[5]);
      tmax = atof(argv[6]);
      tincr = atof(argv[7]);
   }

   fprintf(stderr, "\n  Pref: %6.1lf\n  salt: %6.3lf,%6.3lf,%6.3lf\n theta: %6.2lf,%6.2lf,%6.2lf\n", pref, smin, smax, sincr, tmin, tmax, tincr);

/*  compute sigma at each point in the grid and write to stdout in a
    flat format ... */

    for (s = smin; s <= (smax + .01); s += sincr) {
      for (th = tmin; th <= tmax; th += tincr) {
          tref = hb_theta(s, th, zero, pref);
          hb_svan(s, tref, pref, &sig);
          fprintf(stdout, " %6.3lf %6.3lf %6.3lf\n", s, th, sig);
      }
    }

   exit(0);
}   /* end main */
