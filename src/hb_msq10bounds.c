/* hb_msq10bounds.c
   
   Accepts a 10-deg MSq, an (optional) 5-deg suband (optional) increment as arguments and
   returns bounds in the proper format: w/e/s/n .  Specifying
   an increment will enlarge the bounds by half the increment
   -- useful to generate node-grids as opposed to pixel-grids.
   Use -L+ to force longitudes in range 0->360 (default is -180->+180)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_usage(char *);

main(int argc, char **argv)
{

   int msq10, i, hem, rem, error, no_offset, all_pos_lon;
   int sflag;
   float incr;
   float west, east, south, north;
   float lat, lon;
   char *s, *subarea;
   
   
   if (argc < 2) { 
      print_usage(argv[0]); 
      exit(1); 
   }
  
  incr = 0.0;
  no_offset = 1;
  error = 0;
  all_pos_lon = 0;
  sflag = 0;
  subarea = (char *) calloc(3, sizeof(char));
  
/* parse the  command line arguments */
   
    for (i = 1; i < argc; i++) { 
      if (argv[i][0] == '-') { 
         s = &argv[i][1]; 
         switch (*s) { 
            case 'I': 
               ++s;
               incr = (float) atof(s);
               incr *= .5;
               no_offset = 0; 
               break;
            case 'L': 
	       all_pos_lon = 1;
	       break;
            case 'S': 
	       sflag = 1;
	       ++s;
	       subarea[0]= *s;
	       subarea[1] = *(++s);
	       break;
            case 'h': 
               print_usage(argv[0]); 
               exit(0);
	       
            default:  
               error = 1;
         } /* end switch */
             
      } /* end if */
      else {

         if (! (msq10 = atoi(argv[i]))) {
          fprintf(stderr,"\nUnable to parse MSQ10: %s \n", argv[i]);
          exit(1);
         }
      }

      if (error ) { 
        fprintf(stderr,"\nError parsing command line args.\n");
        fprintf(stderr," in particular:  '%s'\n", argv[i]); 
        exit(1); 
      }
   } /* end for */
 
   hem = msq10 / 1000; 
   rem = msq10 - (hem * 1000);
   lat = (float)(rem / 100 * 10);
   lon = (float) (rem % 100 * 10);
   
      switch (hem) {   
         case 7: 
           if (all_pos_lon) {
	      west = 350-lon- incr;
	      east = 360-lon + incr;
	      south = lat - incr;
	      north = lat+10 + incr;
	   }
	   else {
	      west = -(lon +10 + incr);
	      east = -lon + incr;
	      south = lat - incr;
	      north = lat+10 + incr;
           }

	   if (sflag)  {
              switch (subarea[0]) {
	  
	        case '0':
		   west += 5;
		   north -= 5;
	        break;
	        case '1':
		   east -= 5;
		   north -= 5;
	        break;
	        case '2':
		   west += 5;
		   south += 5;
	        break;
	        case '3':
		   east -= 5;
		   south += 5;
	        break;
	        default:
                  fprintf(stderr, "\nCould not recognize subarea %c\n", subarea[0]);
                  exit(3);  
	        } /* end switch */
	   }	
         break;
      
        case 5: 
           if (all_pos_lon) {
	      west = 350-lon- incr;
	      east = 360-lon + incr;
	      south =-(lat+10 + incr);
	      north = -lat+ incr;
	   }
	   else {
	      west = -(lon +10 + incr);
	      east = -lon + incr;
	      south = -(lat+10 + incr);
	      north = -lat+ incr;
           }
	   if (sflag)  {
              switch (subarea[0]) {
	  
	        case '0':
		   south += 5;
		   west += 5;
	        break;
	        case '1':
		   south += 5;
		   east -=  5;
	        break;
	        case '2':
		   north -= 5;
		   west += 5;
	        break;
	        case '3':
		   north -= 5;
		   east -= 5;
	        break;
	        default:
                  fprintf(stderr, "\nCould not recognize subarea %c\n", subarea[0]);
                  exit(3);  
	        } /* end switch */
	   }	
        break;
      
        case 3: 
           west = lon - incr;
           east = lon + 10 + incr;
           south = -lat -10 - incr;
           north = -lat + incr;
	   if (sflag)  {
              switch (subarea[0]) {
	  
	        case '0':
		   south +=  5;
		   east -=  5; 
	        break;
	        case '1':
		   south += 5;
		   west += 5;
	        break;
	        case '2':
		   east -= 5; 
		   north -= 5;
	        break;
	        case '3':
		   west += 5;
		   north -= 5;
	        break;
	        default:
                  fprintf(stderr, "\nCould not recognize subarea %c\n", subarea[0]);
                  exit(3);  
	        } /* end switch */
	   }	
        break;
      
        case 1:
	   west = lon - incr;
	   east = lon+10 + incr;
	   south = lat - incr;
	   north = lat+10+incr;
	   if (sflag)  {
              switch (subarea[0]) {
	  
	        case '0':
		   north -= 5;
		   east -= 5;
	        break;
	        case '1':
		   north -= 5;
		   west += 5;
	        break;
	        case '2':
		   south += 5;
		   east -= 5;
	        break;
	        case '3':
		   south += 5;
		   west += 5;
	        break;
	        default:
                  fprintf(stderr, "\nCould not recognize subarea %c\n", subarea[0]);
                  exit(3);  
	        } /* end switch */
	   }	
        break;
	default:
           fprintf(stderr, "\nCould not recognize the hemisphere in MSQ #%d\n", msq10);
           exit(2);  
	
      } /* end switch */
      
      if (no_offset) {
         fprintf(stdout,"%.0f/%.0f/%.0f/%.0f \n", west , east, south, north);
         exit(0);
      }
      
      fprintf(stdout,"%.1f/%.1f/%.1f/%.1f \n", west , east, south, north);
      exit(0);
    
} /* end main */
/****************************************************************************/

void print_usage(char *program) 
{ 
   fprintf(stderr,"\n Accepts a 10-deg MSq, an optional 5-deg subarea, and (optional) increment as arguments and");
   fprintf(stderr,"\nreturns bounds in the proper format: w/e/s/n  .  Specifying");
   fprintf(stderr,"\nan increment will enlarge the bounds by half the increment");
   fprintf(stderr,"\n-- useful to generate node-grids as opposed to pixel-grids.");

   fprintf(stderr,"\n\nUsage:  %s msq_10 [-I<incr>] [-L] [-S<subarea>] ", program);
   fprintf(stderr,"\n\n msq_10 is the 4-digit WMO (Marsden Square)."); 
   fprintf(stderr,"\n \n  OPTIONS:"); 
   fprintf(stderr,"\n[-I]  : specifies grid increment."); 
   fprintf(stderr,"\n        Bounds will be offset by half the increment."); 
   fprintf(stderr,"\n[-L]  : Output longitudes positive [0 - 360]"); 
   fprintf(stderr,"\n[-S]  : specify 5-deg subarea [0, 1, 2, 3]."); 
   fprintf(stderr,"\n        bounds will reflect 5-deg square."); 
   fprintf(stderr,"\n[-h]  : help -- prints this message."); 
  fprintf(stderr,"\n\n");
   return;
}
   
   
