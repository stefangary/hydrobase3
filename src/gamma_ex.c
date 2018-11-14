/*  gamma_ex.c

  test the HydroBase C version of David Jackett's neutral density code
 
*/

#include <stdio.h>
#include <stdlib.h>
#include "hb_gamma.h"
#include "hb_paths.h"

#define MAXZ 6000
#define NLEVS 3

int main(int argc, char **argv) 

{
   int i, nobs, error;
   double lat, lon;
   double s[MAXZ], t[MAXZ], p[MAXZ], gamma[MAXZ];
   double dgl[MAXZ], dgh[MAXZ];
   double glevels[NLEVS] = {26.8, 27.9, 28.1};
   double sns[NLEVS], tns[NLEVS], pns[NLEVS];
   double dsns[NLEVS], dtns[NLEVS], dpns[NLEVS];
   
   struct GAMMA_NC gamma_info;

   char *gamma_nc_path;   
   FILE *datfile;
   
   
   gamma_nc_path = (char *) GAMMA_NC_PATH;
   
   /* an example of labelling data */
   
   datfile = fopen("/d5/ruth/HB2/Jackett/gamma/example.dat", "r");
   if (datfile == NULL) {
      fprintf(stderr, "\nUnable to open /d5/ruth/HB2/Jackett/gamma/example.dat\n");
      exit(1);
   }
   
   fscanf(datfile, "%lf%lf%d", &lon, &lat, &nobs);
   for (i = 0; i < nobs; ++i) {
      fscanf(datfile,"%lf%lf%lf", &s[i], &t[i], &p[i]);
   }

   /* label */
      
   gamma_nc_init(gamma_nc_path, &gamma_info);
   
   error = gamma_n(&gamma_info, s, t, p, nobs, lon, lat, gamma, dgl, dgh);
   if (error < 0) {
      fprintf(stderr,"\nERROR returned by gamma_n()\n");
      exit(1);
   }
   
   fprintf(stdout,"\nLocation:  %8.3lf %8.3lf\n", lon, lat);
   fprintf(stdout,"\nLabels:");
   
   for (i = 0; i < nobs; ++i) {
      fprintf(stdout, "\n %8.1lf %11.6lf %10.6lf %10.6lf", p[i], gamma[i], dgl[i], dgh[i]);
   }
   
   
   /* fit some surfaces */
   
   neutral_surfaces(s, t, p, gamma, nobs, glevels, (int)NLEVS, sns, tns, pns, dsns, dtns, dpns);
   
   fprintf(stdout,"\n\nSurfaces:");
      
   for (i = 0; i < NLEVS; ++i) {
      fprintf(stdout, "\n%8.2lf %12.6lf %12.6lf %14.6lf %6.1lf", glevels[i], sns[i], tns[i], pns[i], dpns[i]);
   }
   
   fprintf(stderr,"\n%s successfully completed.\n", argv[0]);
   exit(0);
} /* end main */





