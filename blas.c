#include "blas.h"
/*   Herbert M. Pickett, 22 March 1999 */
/*   Adapted from subset of BLAS routines in LINPAK */
int idamax(n, dx, incx)
int n;
double *dx;
int incx;
{/*  FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
  double dmax, dtmp;
  int i, iret;
  if (n <= 1) return n;
  incx *= sizeof(double);
  dmax = (*dx); 
  if (dmax < 0.) dmax = -dmax;
  iret = 1;
  for (i = 2; i <= n; ++i) {
    dx = (double *) ((char *) dx + incx); 
    dtmp = (*dx); 
    if (dtmp < 0.) dtmp = -dtmp;
    if (dtmp > dmax) {
      dmax = dtmp; iret = i;
    }
  }
  return iret;
} /* idamax */

double dasum(n, dx, incx)
int n;
double *dx;
int incx;
{/* TAKES THE SUM OF THE ABSOLUTE VALUES. */
  double dtmp, dret;
  if (n > 0) {
    incx *= sizeof(double);
    dret = (*dx);
    if (dret < 0.) dret = -dret;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dtmp = (*dx); 
      if (dtmp < 0.) dtmp = -dtmp;
      dret += dtmp;
    }
  } else {
    dret = 0;
  }
  return dret;
} /* dasum */

int daxpy(n, da, dx, incx, dy, incy)
int n, incx, incy;
double da, *dx, *dy;
{/* CONSTANT TIMES A VECTOR PLUS A VECTOR. LOOP UNROLLED */
  double dtmp0, dtmp, *dx0, *dy0;
  if (n > 0 && da != 0.){
    incx *= sizeof(double); incy *= sizeof(double);
    *dy += da * (*dx);
    if((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      *dy += da * (*dx);
      --n;
    }
    n = n >> 1;
    while(n > 0) {
      dy0 = (double *) ((char *) dy + incy);
      dy = (double *) ((char *) dy0 + incy);
      dtmp0 = (*dy0); dtmp = (*dy); 
      dx0 = (double *) ((char *) dx + incx); 
      dx = (double *) ((char *) dx0 + incx);
      *dy0 = dtmp0 + da * (*dx0); *dy = dtmp + da * (*dx);
      --n;
    }
  }
  return 0;
} /* daxpy */

int dcopy(n, dx, incx, dy, incy)
int n, incx, incy;
double *dx, *dy;
{/* COPIES A VECTOR, X, TO A VECTOR, Y. */
  double da;
  if (n > 0) {
    incy *= sizeof(double);
    if (incx == 0) {
      da = dx[0];
      *dy = da;
      while ((--n) > 0) {
        dy = (double *) ((char *) dy + incy);
        *dy = da;
      }
    } else {
      incx *= sizeof(double);
      *dy = (*dx);
      while ((--n) > 0) {
        dx = (double *) ((char *) dx + incx);
        dy = (double *) ((char *) dy + incy);
        *dy = (*dx);
      }
    }
  }
  return 0;
} /* dcopy */

double ddot(n, dx, incx, dy, incy)
int n, incx, incy;
double *dx, *dy;
{/*  FORMS THE DOT PRODUCT OF TWO VECTORS. LOOP UNROLLED */
  double dret, *dx0, *dy0;
  if (n <= 0)  return 0.;
  if (incx == incy && dx == dy) { /* CODE FOR SINGLE VECTOR PRODUCT */
    incx *= sizeof(double);
    dret = (*dx) * (*dx);
    if((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx); 
      dret += (*dx) * (*dx);
      --n; 
    }
    n = n >> 1;
    while(n > 0) {
      dx0 = (double *) ((char *) dx + incx); 
      dx = (double *) ((char *) dx0 + incx);
      dret += (*dx0) * (*dx0) + (*dx) * (*dx);
      --n; 
    }
  } else { /* CODE FOR TWO DIFFERENT VECTORS */
    incx *= sizeof(double); incy *= sizeof(double);
    dret = (*dx) * (*dy);
    if((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dret += (*dx) * (*dy);
      --n; 
    }
    n = n >> 1;
    while(n > 0) {
      dx0 = (double *) ((char *) dx + incx); 
      dx = (double *) ((char *) dx0 + incx);
      dy0 = (double *) ((char *) dy + incy);
      dy = (double *) ((char *) dy0 + incy);
      dret += (*dx0) * (*dy0) + (*dx) * (*dy);
      --n;
    }
  }
  return dret;
} /* ddot */

double dnrm2(n, dx, incx)
int n, incx;
double *dx;
{
  /* calculates the length of vector x */
  double dsum, xabs, xmax, r, sqrt();
  int i, nsum;
  incx *= sizeof(double);
  nsum = 0; dsum = 1.; xmax = 0.;
  for (i = 0; i < n; ++i) {
    xabs = (*dx);
    if (xabs < 0.) xabs = -xabs;
    if (xabs > xmax) {
      if (nsum != 0) { 
        r = xmax / xabs;
        dsum = 1. + (dsum * r) * r;
      }
      xmax = xabs; ++nsum;
    } else if (xabs != 0.){
      r = xabs / xmax;
      dsum += r * r;
      ++nsum;
    }
    dx = (double *)((char *)dx + incx);
  }
  if (nsum > 1) xmax *= sqrt(dsum);
  return xmax;
} /* dnrm2 */

int dscal(n, da, dx, incx)
int n, incx;
double da, *dx;
{/* SCALES A VECTOR BY A CONSTANT. LOOP UNROLLED */
  double dtmp0, dtmp, *dx0;
  if (n > 0 && da != 1.) {
    incx *= sizeof(double);
    *dx *= da;
    if((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx); 
      *dx *= da;
      --n; 
    }
    n = n >> 1;
    while(n > 0) {
      dx0 = (double *) ((char *) dx + incx); 
      dx = (double *) ((char *) dx0 + incx);
      dtmp0 = (*dx0); dtmp = (*dx);
      *dx0 = da * dtmp0; *dx = da * dtmp;
      --n;       
    }
  }
  return 0;
} /* dscal */

int dswap(n, dx, incx, dy, incy)
int n, incx, incy;
double *dx, *dy;
{/* INTERCHANGES TWO VECTORS. */
  double dtmp;
  if (n > 0) {
    incx *= sizeof(double); incy *= sizeof(double);
    dtmp = (*dx); *dx = (*dy); *dy = dtmp;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmp = (*dx); *dx = (*dy); *dy = dtmp;
    }
  }
  return 0;
} /* dswap */

int drot(n, dx, incx, dy, incy, c, s)
int n, incx, incy;
double *dx, *dy, c, s;
{/* APPLIES A PLANE ROTATION. */
  double dtmpx, dtmpy;
  if (n > 0) {
    incx *= sizeof(double); incy *= sizeof(double);
    dtmpx = (*dx); dtmpy = (*dy);
    *dx = dtmpx * c + dtmpy * s;
    *dy = dtmpy * c - dtmpx * s;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmpx = (*dx); dtmpy = (*dy);
      *dx = dtmpx * c + dtmpy * s;
      *dy = dtmpy * c - dtmpx * s;
    }
  }
  return 0;
} /* drot */
