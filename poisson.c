#include "header.h"
#include <time.h>

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}

#define PI 3.14159265358979323846
#define DIMENSION 2
#define FUNC(i,j) (-2*PI*PI*sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))
#define SOL(i,j) (sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))

void GetFuncValues2d(double **coord, int *n, double **f);
void GetError(double **coord, int *n, double **u, double *error);
void UpdateBC(double **coord, double **u, int *n);
void OpA(double *A, double *h);
int ipow(int base, int exp);
int JacobiMalloc(double ***f, double ***u, int *n);
int MultigridMalloc(double ***f, double ***u, int *n, int levels);

int main() {
	
	double weight=(2.0/3.0);
	int    n[DIMENSION], ierr=0, levels, numIter;
	double **coord, h[DIMENSION], bounds[DIMENSION*2];
	double **f, **u, *rnorm, error[3], As[5];// ,**r;
	FILE   *solData, *resData, *errData;
	
	freopen("poisson.in", "r", stdin);
	freopen("poisson.out", "w", stdout);
	freopen("poisson.err", "w", stderr);
	// Inputs
	//printf("Enter the no .of points in each dimension = ");
	scanf("%d",n); // unTotal is used temporarily
	//printf("Enter the no .of iterations = ");
	scanf("%d",&numIter);
	//printf("Enter the no .of Multigrid levels = ");
	scanf("%d",&levels);
	
	clock_t begin = clock();

	for (int i=1;i<DIMENSION;i++) {
		n[i]  = n[0];      // No. of points in each dimension
	}
	for (int i=0;i<DIMENSION;i++) {
		bounds[i*2] = 0.0;    // Lower bound in each dimension
		bounds[i*2+1] = 1.0;  // Upper bound in each dimension
	}

	
	// Memory allocation of RHS, solution and residual
	//ierr = JacobiMalloc(&f,&u,n); CHKERR_PRNT("malloc failed");
	ierr = MultigridMalloc(&f,&u,n,levels); CHKERR_PRNT("malloc failed");
	rnorm = (double *)malloc((numIter+1)*sizeof(double));if (rnorm==NULL) ERROR_MSG("malloc failed");

	clock_t memT = clock();
	// Meshing
	ierr = UniformMesh(&coord,n,bounds,h,DIMENSION); CHKERR_PRNT("meshing failed");
	
	clock_t meshT = clock();
	
	// f values
	GetFuncValues2d(coord,n,f);
	
	// Stencil operator coefficients
	OpA(As,h);
	
	// Update 'u' with boundary conditions
	UpdateBC(coord,u,n);
	
	clock_t initT = clock();
	// Solver
	//Jacobi(u,f,As,weight,rnorm,numIter,n); // Weighted Jacobi
	
	// Multigrid V-cycle
	Multigrid(u,f,As,weight,rnorm,levels,n,numIter);

	clock_t solverT = clock();
	
	// Error computation
	GetError(coord,n,u,error);
	
	// Output
	solData = fopen("uData.dat","w");
	resData = fopen("rData.dat","w");
	errData = fopen("eData.dat","w");

	for(int i=0;i<3;i++){
		printf("\nerror[%d] = %.16e\n",i,error[i]);
		fprintf(errData,"%.16e\n",error[i]);
	}

	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			fprintf(solData,"%.16e ",u[i][j]);
		}
		fprintf(solData,"\n");
	}

	for (int i=0;i<numIter+1;i++) {
		fprintf(resData,"%.16e ",rnorm[i]);
	}
	fprintf(resData,"\n");

	clock_t ppT = clock();
	
	printf("Total time:             %lf\n",(double)(ppT-begin)/CLOCKS_PER_SEC);
	printf("Memory allocation time: %lf\n",(double)(memT-begin)/CLOCKS_PER_SEC);
	printf("Meshing time:           %lf\n",(double)(meshT-memT)/CLOCKS_PER_SEC);
	printf("Initialization time:    %lf\n",(double)(initT-meshT)/CLOCKS_PER_SEC);
	printf("Solver time:            %lf\n",(double)(solverT-initT)/CLOCKS_PER_SEC);
	printf("Post processing time:   %lf\n",(double)(ppT-solverT)/CLOCKS_PER_SEC);
	
	fclose(solData);
	fclose(resData);
	fclose(errData);
	free2dArray(&coord);
	free2dArray(&f);
	free2dArray(&u);
	free(rnorm);
	
	return 0;
}


void GetFuncValues2d(double **coord, int *n, double **f) {

	// f(x,y) = -2*PI^2*sin(Pi*x)*sin(pi*y)	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			f[i][j] = FUNC(i,j);
		}
	}

}

void GetError(double **coord, int *n, double **u, double *error) {
	
	// u(x,y) = sin(Pi*x)*sin(pi*y)	
	double diff;
	error[0] = 0.0;
	error[1] = 0.0;
	error[2] = 0.0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			diff = fabs(u[i][j]-SOL(i,j));
			error[0] = fmax(diff,error[0]);
			error[1] = error[1] + diff;
			error[2] = error[2] + diff*diff;
		}
	}
	error[2] = sqrt(error[2]);
}

void UpdateBC(double **coord, double **u, int *n) {

	int iend;
	
	for (int j=0;j<n[0];j++) {
		u[0][j] = SOL(0,j);
	}
	
	iend = n[0]-1;
	for (int i=1;i<n[1]-1;i++) {
		u[i][0] = SOL(i,0);
		u[i][iend] = SOL(i,iend);
	}

	iend = n[1]-1;
	for (int j=0;j<n[0];j++) {
		u[iend][j] = SOL(iend,j);
	}
	
}

void OpA(double *A, double *h) {
	
	double hy2, hx2;
	
	hx2 = h[0]*h[0];
	hy2 = h[1]*h[1];
	A[0] = 1/hy2;
	A[1] = 1/hx2;
	A[2] = -2*((1/hx2)+(1/hy2));
	A[3] = 1/hx2;
	A[4] = 1/hy2;

}

int ipow(int base, int exp) {

	int result = 1;
	while (exp) {
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

int JacobiMalloc(double ***f, double ***u, int *n) {
	
	int ierr = 0;

	ierr = malloc2d(f,n[1],n[0]); CHKERR_RETURN("malloc failed");
	ierr = malloc2d(u,n[1],n[0]); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
		}
	}

	return ierr;
}

int MultigridMalloc(double ***f, double ***u, int *n, int levels) {
	
	int TotalRows, n1, n0, *m, k, ierr = 0;

	TotalRows = (2*(n[1]-1)*(ipow(2,levels)-1))/(ipow(2,levels))+levels;
	m = (int *)malloc(TotalRows*sizeof(int)); if (m==NULL) ERROR_RETURN("malloc failed"); 
	k = 0;
	for (int i=0;i<levels;i++) {
		n1 = (n[1]+ipow(2,i)-1)/(ipow(2,i));
		n0 = (n[0]+ipow(2,i)-1)/(ipow(2,i));
		for (int j=0;j<n1;j++) {
			m[k] = n0;
			k = k+1;
		}
	}
	
	ierr = malloc2dY(f,TotalRows,m); CHKERR_RETURN("malloc failed");
	ierr = malloc2dY(u,TotalRows,m); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<TotalRows;i++) {
		for (int j=0;j<m[i];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
		}
	}
	free(m);
	return ierr;
}


