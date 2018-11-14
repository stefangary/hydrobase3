/* mem_subs.c
................................................................................
                          *******  HydroBase2  *******
               ...................................................
                          
                    author:  Ruth  Curry
                             Woods Hole Oceanographic Institution
                             2000
................................................................................
*
*   Handles allocation and reallocation of dynamic memory
*
................................................................................
 * Modules in this file:
 *
 *  void *get_memory(void *prev_addr, size_t nelem, size_t size)  
 * ................................................................................
*/ 
#include <stdio.h>
#include <stdlib.h>
#include "hb_memory.h"

void *get_memory (void *prev_addr, size_t nelem, size_t size)

/*       Allocates memory if prev_addr is NULL or reallocates chunks of memory
       if prev_addr points to some address  */
{
	void *tmp;

	if (nelem == 0) return((void *)NULL); /* Take care of n = 0 */
	
	if (prev_addr) {
          if ((tmp = realloc ((void *) prev_addr, (size_t)(nelem * size))) == NULL) {
	    fprintf(stderr, "\nError:  unable to reallocate  memory, n = %d\n", (int)nelem);
            exit (1);
	  }
	}
	else {
	   if ((tmp = calloc ((size_t) nelem, (unsigned) size)) == (void *)NULL) {
	      fprintf (stderr, "Error: Unable to allocate memory, n = %d\n", (int)nelem);
            exit (1);
	      
	   }
	}
	return (tmp);
}
/*********************************************************************/
