/* hb_memory.h 

................................................................................
                          *******  HydroBase 2 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             2001  
................................................................................
*/

/****************************************************************************/

extern void *get_memory(void *prev_addr, size_t nelem, size_t size);  
/*     Allocates memory if prev_addr is NULL or reallocates chunks of memory
       if prev_addr points to some address
*/
