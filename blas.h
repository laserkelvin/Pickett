#ifndef BLAS_H
int idamax (int  n, double *dx, int  incx);
double dasum (int  n, double *dx, int  incx);
int daxpy (int  n, double  da, double *dx, int  incx, double *dy, int  incy);
int dcopy (int  n, double *dx, int  incx, double *dy, int  incy);
double ddot (int  n, double *dx, int  incx, double *dy, int  incy);
double dnrm2(int  n, double *dx, int  incx);
int dscal (int  n, double  da, double *dx, int  incx);
int dswap (int  n, double *dx, int  incx, double *dy, int  incy);
int drot (int  n, double *dx, int  incx, double *dy, int  incy, double  c, double  s);
int isamax (int  n, float *dx, int  incx);
float sasum (int  n, float *dx, int  incx);
int saxpy (int  n, float  da, float *dx, int  incx, float *dy, int  incy);
int scopy (int  n, float *dx, int  incx, float *dy, int  incy);
float sdot (int  n, float *dx, int  incx, float *dy, int  incy);
float snrm2(int  n, float *dx, int  incx);
int sscal (int  n, float  da, float *dx, int  incx);
int sswap (int  n, float *dx, int  incx, float *dy, int  incy);
int srot (int  n, float *dx, int  incx, float *dy, int  incy, float  c, float  s);
#define BLAS_H
#endif

