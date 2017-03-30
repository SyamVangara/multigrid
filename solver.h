



#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <math.h>

void UpdateRHS(double *A, double **u, double **r, int *n);
double Residual(double **u, double **f, double **r, double *As, int *n);
void JacobiStep(double **u, double **f, double *As, double w, int *n);
void Jacobi(double **u, double **f, double **r, double *As, double w, double *rnorm, int v,int *n);
void ResidualRestriction(double **f, double **r, int *n);
void ErrorCorrection(double **u, int *n, int flag);
void Vcycle(double **u, double **f, double **r, double *As, double w, int *v,int levels,int *n);
void Multigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int *n,int m);
void AsyncMultigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int*n, int m);
double L2norm(double *a, int n);
double L1Norm(double *a, int n);
double LiNorm(double *a, int n);
double norm(double *a, int n);
void Initialization(double **u, int *n);

#endif

