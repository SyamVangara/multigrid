/*


*/


#ifndef _ARRAY_H_
#define _ARRAY_H_ 

#include <stdlib.h>
// Function prototypes

// Allocate contiguous memory for 2d array
extern int malloc2d(double ***a, int n, int m);

// Allocate contigous memory for 2d array with variable row length
extern int malloc2dY(double ***a, int n, int *m);

// Free the allocated memory of 2d arrays
extern void free2dArray(double ***a);

#endif
