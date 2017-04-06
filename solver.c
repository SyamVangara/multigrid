#include "solver.h"

void UpdateRHS(double *A, double **u, double **r, int *n) {
	
	int iend;
	
	for (int j=1;j<n[0]-1;j++) {
		r[1][j] = r[1][j] - (u[0][j])*A[0];
		//f[0][j-1] = f[0][j-1] - (SOL(0,j))*A[0];
	}
	
	iend = n[0]-3;
	for (int i=1;i<n[1]-1;i++) {
		r[i][1]    = r[i][1]    - (u[i][0])*A[1];
		r[i][iend] = r[i][iend] - (u[i][iend+1])*A[1];
		//f[i-1][0]    = f[i-1][0]    - (SOL(i,0))*A[1];
		//f[i-1][iend] = f[i-1][iend] - (SOL(i,iend+2))*A[1];
	}

	iend = n[1]-3;
	for (int j=1;j<n[0]-1;j++) {
		r[iend][j] = r[iend][j] - (u[iend+1][j])*A[0];
		//f[iend][j-1] = f[iend][j-1] - (SOL(iend+2,j))*A[0];
	}
	

}

double Residual(double **u, double **f, double **r, double *As, int *n) {

	double resLinf;
//	*rnorm = 0.0;
	//*(rnorm+1) = 0.0;
	//*(rnorm+2) = 0.0;
	resLinf = 0.0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			r[i][j] = f[i][j] - (As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[2]*u[i][j]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			//res = fabs(r[i][j]);
/*
			*rnorm = fmax(res,*rnorm);
			*(rnorm+1) = *(rnorm+1) + res;
			*(rnorm+2) = *(rnorm+2) + res*res;
*/
			//rnorm = rnorm + res*res;
			resLinf = fmax(resLinf,fabs(r[i][j]));
		}
	}
	//rnorm = sqrt(rnorm);
	return resLinf;
}

void JacobiStep(double **u, double **f, double *As, double w, int *n) {
	
	double temp;
	
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			
			temp = f[i][j] - (As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			u[i][j] = (1-w)*u[i][j] + (w/As[2])*temp;
			
			//u[i][j] = u[i][j] + (w/As[2])*r[i][j];
		}
	}
}

void Jacobi(double **u, double **f, double **r, double *As, double w, double *rnorm, int v,int *n) {
	
	int i=0;
	rnorm[0] = Residual(u,f,r,As,n);
	while (i<v && (1.0+0.5*rnorm[i])!=1.0) {
		i = i+1;
		JacobiStep(u,f,As,w,n);
		rnorm[i] = Residual(u,f,r,As,n);
		//GetResidual(*u,*f,As,shift,*r,nt);
		//res = norm(*r,nt);
	}
	printf("residual = %.16e\n",rnorm[i]);

}

void ResidualRestriction(double **f, double **r, int *n) {
	
	for (int i=2;i<n[1]-1;i=i+2) {
		for (int j=2;j<n[0]-1;j=j+2) {
			/*	
			f[n[1]+i/2][j/2] = f[i][j]-(As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[2]*u[i][j]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			*/
			f[n[1]+i/2][j/2] = r[i][j];
		}
	}
}

void ErrorCorrection(double **u, int *n, int flag) {
	
	double Iop[3][3];
	int    im, jm;

	for (int lj=0;lj<3;lj++) {
 		Iop[0][lj]= 0.5 - 0.25*fabs(1-lj);
 		Iop[1][lj]= 1.0 - 0.5*fabs(1-lj);
 		Iop[2][lj]= 0.5 - 0.25*fabs(1-lj);
	}
	for (int i=2;i<n[1]-1;i=i+2) {
		for (int j=2;j<n[0]-1;j=j+2) {
			im = n[1]+i/2;
			jm = j/2;
			for (int li=0;li<3;li++) {
				for (int lj=0;lj<3;lj++) {
			 		u[i+li-1][j+lj-1] = u[i+li-1][j+lj-1]+Iop[li][lj]*u[im][jm];
				}
			}
			if (flag==0) u[im][jm] = 0.0;
		}
	}
}

void SweepAndRestrict(double **u, double **f, double **r, double *As, double w, int v,int *n) {
	
	double tempRes;
	for (int i=0;i<v;i++) {
		//JacobiStep(u,f,As,w,n);
		//tempRes = Residual(u,f,r,As,n);
		JacobiStep(u,f,As,w,n);
	}
	tempRes = Residual(u,f,r,As,n);
	ResidualRestriction(f,r,n);
}

void CorrectAndSweep(double **u, double **f, double *As, double w, int v,int *n) {
	
	//double tempRes;
	ErrorCorrection(u,n,0);
	for (int i=0;i<v;i++) {
		//tempRes = Residual(u,f,r,As,n);
		JacobiStep(u,f,As,w,n);
	}

}

void Vcycle(double **u, double **f, double **r, double *As, double w, int *v,int levels,int *n) {
	
	double AsH[levels][5];//, res;
	int    nH[levels][2], nid[levels];
	
	for (int j=0;j<5;j++) {
		AsH[0][j] = As[j];
	}
	
	nH[0][0] = n[0];
	nH[0][1] = n[1];
	nid[0] = 0;
	for (int i=1;i<levels;i++) {
		for (int j=0;j<5;j++) {
			AsH[i][j] = 0.25*AsH[i-1][j];
		}
		nH[i][0] = (nH[i-1][0]+1)/2;
		nH[i][1] = (nH[i-1][1]+1)/2;
		nid[i] = nid[i-1] + nH[i-1][1];
	}
	
	for (int i=0;i<levels-1;i++) {
		SweepAndRestrict((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	
	for (int i=0;i<v[1];i++) {
		//res = Residual((u+nid[levels-1]),(f+nid[levels-1]),(r+nid[levels-1]),AsH[levels-1],nH[levels-1]);
		JacobiStep((u+nid[levels-1]),(f+nid[levels-1]),AsH[levels-1],w,nH[levels-1]);
	}
	for (int i=levels-2;i>=0;i=i-1) {
		CorrectAndSweep((u+nid[i]),(f+nid[i]),AsH[i],w,v[0],nH[i]);
	}
}

void Multigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int *n,int m) {

	int v[2];

	//printf("Enter the number of fine grid sweeps = ");
	scanf("%d",v);
	//printf("Enter the number of coarse grid sweeps = ");
	scanf("%d",v+1);
	
	rnorm[0] = Residual(u,f,r,As,n);
	for (int i=1;i<m+1;i++) {
		Vcycle(u,f,r,As,w,v,levels,n);
		rnorm[i] = Residual(u,f,r,As,n); 
	}
	printf("residual = %.16e\n",rnorm[m]);

}

void Pcycle(double **u, double **f, double **r, double *As, double w,int levels,int *n) {
	
	double AsH[levels][5], res;
	int    nH[levels][2], nid[levels], flag;
	
	for (int j=0;j<5;j++) {
		AsH[0][j] = As[j];
	}
	
	nH[0][0] = n[0];
	nH[0][1] = n[1];
	nid[0] = 0;
	for (int i=1;i<levels;i++) {
		for (int j=0;j<5;j++) {
			AsH[i][j] = 0.25*AsH[i-1][j];
		}
		nH[i][0] = (nH[i-1][0]+1)/2;
		nH[i][1] = (nH[i-1][1]+1)/2;
		nid[i] = nid[i-1] + nH[i-1][1];
	}
	
	res = Residual(u,f,r,As,n);
	ResidualRestriction(f,r,n);
	for (int i=1;i<levels-1;i++) {
		res = Residual((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],nH[i]);
		ResidualRestriction((f+nid[i]),(r+nid[i]),nH[i]);
		//SweepAndRestrict((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	flag = 0;
	for (int i=levels-2;i>=0;i=i-1) {
		ErrorCorrection(u+nid[i],nH[i],flag);
		//CorrectAndSweep((u+nid[i]),(f+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	for (int i=0;i<levels;i++) {
		//res = Residual((u+nid[levels-1]),(f+nid[levels-1]),(r+nid[levels-1]),AsH[levels-1],nH[levels-1]);
		JacobiStep((u+nid[i]),(f+nid[i]),AsH[i],w,nH[i]);
	}
}
void PMultigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int*n, int m) {
	
	int i, flag;
	double AI[5];	
/*
	ResidualRestriction(f,f,n); // Building f-tilda
	OpAI(As,AI);

	AsyncRres(u,f,r,As,n);
	AsyncCorrection(u,r,AI,n,flag);
	AsyncRestriction(u,r,As,n);
	rnorm[0] = AsyncResNorm(u,r,As,n);
*/
	i = 0;
	rnorm[i] = Residual(u,f,r,As,n);
	while (i<m && (1.0+0.5*rnorm[i])!=1.0) {
		i = i+1;
		Pcycle(u,f,r,As,w,levels,n);
		rnorm[i] = Residual(u,f,r,As,n);
		//GetResidual(*u,*f,As,shift,*r,nt);
		//res = norm(*r,nt);
	}
	printf("residual = %.16e\n",rnorm[i]);
	
}

double L2norm(double *a, int n) {
	
	double result;
	result = a[0]*a[0];
	for (int i=1;i<n;i++) {
		result = result + a[i]*a[i];
	}
	return sqrt(result);
}

double L1Norm(double *a, int n) {
	
	double result;
	result = fabs(a[0]);
	for (int i=1;i<n;i++) {
		result = result + fabs(a[i]);
	}
	return result;
}

double LiNorm(double *a, int n) {
	
	double result;
	result = fabs(a[0]);
	for (int i=1;i<n;i++) {
		result = fmax(result,fabs(a[i]));
	}
	return result;
}

void Initialization(double **u, int *n) {
	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			u[i][j] = 0.0;
		}
	}
}

