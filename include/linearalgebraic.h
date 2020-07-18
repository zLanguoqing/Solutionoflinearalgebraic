#ifndef _LINEARALGEBRAIC_H_
#define _LINEARALGEBRAIC_H_
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define MIN (1e-30)
int Gauss(double a[], double b[], int n);
int GaussJordan(double a[], double b[], int n, int m);
int ChaseTridiagonal(double b[], int n, int m, double d[]);
int Band(double b[],double d[],int n,int l,int il,int m);
int LDLT(double a[],int n,int m,double c[]);
int Cholesky(double a[],int n,int m,double d[]);
  int GaussSeidel(double a[],double b[],int n,double x[],double eps);
#endif