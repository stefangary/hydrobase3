/*   hb_sigma.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
			     March 2000
................................................................................
................................................................................
.  Computes potential density given a pressure, temp, salinty 
................................................................................

*/

#include <stdio.h>
#include <stdlib.h>
#include "hydrobase.h"


main (int argc, char **argv)
{
   double tref, sigma, p, t, s, pref;
   double zero = 0.0;


/* are there command line arguments? */

   if (argc < 5) {
      fprintf(stderr, "\nComputes potential density rel to pref from");
      fprintf(stderr, "\nspecified values of pressure, temp, and salt");
       fprintf(stderr, "\n\nUsage:  %s  p t s  pref\n\n", argv[0]);
     exit(1);
   }

    p = strtod(argv[1],(char **)NULL);
    t = strtod(argv[2],(char **)NULL);
    s = strtod(argv[3],(char **)NULL);
    pref = strtod(argv[4],(char **)NULL);

    tref = hb_theta(s, t, zero, pref);
/*    fprintf(stdout," \n pot temp = %.3lf  tref = %.3lf ", t, tref); */
    hb_svan(s, tref, pref, &sigma);
    fprintf(stdout, " %lf \n",  sigma);

}  /* end main() */
