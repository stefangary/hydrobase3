/*   hb_theta.c

................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
			     March 2000
			     upated to HydroBase 3 Nov 2011
................................................................................
................................................................................
.  Computes potential temperature given a pressure, in situ temp, salinty 
................................................................................

*/

#include <stdio.h>
#include <stdlib.h>
#include "hydrobase.h"


main (int argc, char **argv)
{
   double tref, p, t, s, pref;


/* are there command line arguments? */

   if (argc < 4) {
      fprintf(stderr, "\nComputes potential temperature (at pref) from specified");
      fprintf(stderr, "\nvalues of pressure, in situ temp, and salt");
       fprintf(stderr, "\n\nUsage:  %s  p  t  s pref  \n\n", argv[0]);
     exit(1);
   }

    p = strtod(argv[1],(char **)NULL);
    t = strtod(argv[2],(char **)NULL);
    s = strtod(argv[3],(char **)NULL);
   pref = strtod(argv[4],(char **)NULL);

    tref = hb_theta(s, t, p, pref);
    fprintf(stdout, " %.3lf \n",  tref);

}  /* end main() */
