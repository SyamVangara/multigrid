



#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdio.h>
#include <math.h>

void UpdateRHS(double *A, double **u, double **r, int *n);
double ResidualNorm(double **u, double **f, double *As, int *n);
void JacobiStep(double **u, double **f, double *As, double w, int *n);
void Jacobi(double **u, double **f, double *As, double w, double *rnorm, int v,int *n);
void ResidualRestriction(double **u, double **f, double *As, int *n);
void ErrorCorrection(double **u, int *n);
void Vcycle(double **u, double **f, double *As, double w, int *v,int levels,int *n);
void Multigrid(double **u, double **f, double *As, double w, double *rnorm, int levels, int *n,int m);
double norm(double *a, int n);
void Initialization(double **u, int *n);

#endif

