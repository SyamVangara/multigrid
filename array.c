#include "array.h"

#define SUM_ARRAY(result, a, size ) {int i; result=0; for (i=0;i<size;i++)\
	{result += a[i];}}


int malloc2d(double ***a, int n, int m) {
	// Allocate contiguous memory for 2d array "A" of (nxm) 
	//
	// Pointer to the pointer of array - "&A" (= a)
	// Number of rows (n), columns - m; 
	// Location of function call - fle_name, line_num
	//
	// Array "A" can be accessed using A[i][j]                           

	*a = (double **)malloc(n*sizeof(double *));
	**a = (double *)malloc(n*m*sizeof(double));
	if (*a==NULL||**a==NULL) return 1;
	
	// Assign a pointer to each row
	for(int i=1;i<n;i++){
		*(*a+i) = *(*a+i-1) + m;
	}

	return 0;
}

int malloc2dY(double ***a, int n, int *m) {
	// Allocate contiguous memory for 2d array "A" with variable row length (nxm(i))
	//
	// Pointer to the pointer of array - "&A" (= a)
	// Number of rows (n), columns - m(i) in i^th row 
	// Location of function call - fle_name, line_num
	//
	// Array "A" can be accessed using A[i][j]                           
	
	int aTotal;

	// aTotal - total number of elements in A
	aTotal=0; 
	for (int i=0;i<n;i++) aTotal += m[i];

	//SUM_ARRAY(aTotal, m, n); // aTotal - total number of elements in A
	
	*a = (double **)malloc(n*sizeof(double *));
	**a = (double *)malloc(aTotal*sizeof(double));
	if (*a==NULL||**a==NULL) return 1;
	
	// Assign a pointer to each row
	for(int i=1;i<n;i++){
		*(*a+i) = *(*a+i-1) + *(m+i-1);
	}

	return 0;
}

void free2dArray(double ***a) {
	//free the allocated memory
	free(**a);
	free(*a);
}

