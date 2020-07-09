#ifndef _LINEARALGEBRAIC_H_
#define _LINEARALGEBRAIC_H_
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define MIN (1e-30)
int Gauss(double a[], double b[], int n);
int GaussJordan(double a[], double b[], int n, int m);
int TriMatrix(double b[], int n, int m, double d[]);
#endif