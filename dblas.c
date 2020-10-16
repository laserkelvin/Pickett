#include <stdlib.h>
#include <math.h>
#include "cblas.h"
/*   Herbert M. Pickett, 22 March 1999 */
/*   Adapted from subset of BLAS routines in LINPAK */
CBLAS_INDEX cblas_idamax(const int n, const double *dx, const int incX)
{  /*  FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
  double dmax, dtmp;
  int i, iret, incx;
  if (n <= 1)
    return 0;
  incx = incX * (int) sizeof(double);
  dmax = (*dx);
  if (dmax < 0.)
    dmax = -dmax;
  iret = 0;
  for (i = 1; i < n; ++i) {
    dx = (double *) ((char *) dx + incx);
    dtmp = (*dx);
    if (dtmp < 0.)
      dtmp = -dtmp;
    if (dtmp > dmax) {
      dmax = dtmp;
      iret = i;
    }
  }
  return (size_t) iret;
}                               /* idamax */

double cblas_dasum(const int N, const double *dx, const int incX)
{ /* TAKES THE SUM OF THE ABSOLUTE VALUES. */
  double dtmp, dret;
  int n, incx;
  n = N;
  if (n > 0) {
    incx = incX * (int) sizeof(double);
    dret = (*dx);
    if (dret < 0.)
      dret = -dret;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dtmp = (*dx);
      if (dtmp < 0.)
        dtmp = -dtmp;
      dret += dtmp;
    }
  } else {
    dret = 0;
  }
  return dret;
}    /* dasum */

void cblas_daxpy(const int N, const double da, const double *dx,
                 const int incX, double *dy, const int incY)
{  /* CONSTANT TIMES A VECTOR PLUS A VECTOR. LOOP UNROLLED */
  double dtmp0, dtmp, *dx0, *dy0;
  int n, incx, incy;
  n = N;
  if (n > 0 && da != 0.) {
    incx = incX * (int) sizeof(double);
    incy = incY * (int) sizeof(double);
    *dy += da * (*dx);
    if ((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      *dy += da * (*dx);
      --n;
    }
    n = n >> 1;
    while (n > 0) {
      dy0 = (double *) ((char *) dy + incy);
      dy = (double *) ((char *) dy0 + incy);
      dtmp0 = (*dy0);
      dtmp = (*dy);
      dx0 = (double *) ((char *) dx + incx);
      dx = (double *) ((char *) dx0 + incx);
      *dy0 = dtmp0 + da * (*dx0);
      *dy = dtmp + da * (*dx);
      --n;
    }
  }
}  /* daxpy */

void cblas_dcopy(const int N, const double *dx, const int incX,
                 double *dy, const int incY)
{  /* COPIES A VECTOR, X, TO A VECTOR, Y. */
  double da;
  int n, incx, incy;
  n = N;
  if (n > 0) {
    *dy = (*dx);
    incy = incY * (int) sizeof(double);
    if (incX == 0) {
      da = (*dy);
      while ((--n) > 0) {
        dy = (double *) ((char *) dy + incy);
        *dy = da;
      }
    } else {
      incx = incX * (int) sizeof(double);
      while ((--n) > 0) {
        dx = (double *) ((char *) dx + incx);
        dy = (double *) ((char *) dy + incy);
        *dy = (*dx);
      }
    }
  }
} /* dcopy */

double cblas_ddot(const int N, const double *dx, const int incX,
                  const double *dy, const int incY)
{   /*  FORMS THE DOT PRODUCT OF TWO VECTORS. LOOP UNROLLED */
  double dret, *dx0, *dy0;
  int n, incx, incy;
  n = N;
  if (n <= 0)
    return 0.;
  if (incX == incY && dx == dy) {       /* CODE FOR SINGLE VECTOR PRODUCT */
    incx = incX * (int) sizeof(double);
    dret = (*dx) * (*dx);
    if ((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      dret += (*dx) * (*dx);
      --n;
    }
    n = n >> 1;
    while (n > 0) {
      dx0 = (double *) ((char *) dx + incx);
      dx = (double *) ((char *) dx0 + incx);
      dret += (*dx0) * (*dx0) + (*dx) * (*dx);
      --n;
    }
  } else {                      /* CODE FOR TWO DIFFERENT VECTORS */
    incx = incX * (int) sizeof(double);
    incy = incY * (int) sizeof(double);
    dret = (*dx) * (*dy);
    if ((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dret += (*dx) * (*dy);
      --n;
    }
    n = n >> 1;
    while (n > 0) {
      dx0 = (double *) ((char *) dx + incx);
      dx = (double *) ((char *) dx0 + incx);
      dy0 = (double *) ((char *) dy + incy);
      dy = (double *) ((char *) dy0 + incy);
      dret += (*dx0) * (*dy0) + (*dx) * (*dy);
      --n;
    }
  }
  return dret;
}                               /* ddot */

double cblas_dnrm2(const int n, const double *dx, const int incX)
{
  /* calculates the length of vector x */
  double dsum, xabs, xmax, r;
  int i, nsum, incx;
  incx = incX * (int) sizeof(double);
  nsum = 0;
  dsum = 1.;
  xmax = 0.;
  for (i = 0; i < n; ++i) {
    xabs = (*dx);
    if (xabs < 0.)
      xabs = -xabs;
    if (xabs > xmax) {
      if (nsum != 0) {
        r = xmax / xabs;
        dsum = 1. + (dsum * r) * r;
      }
      xmax = xabs;
      ++nsum;
    } else if (xabs != 0.) {
      r = xabs / xmax;
      dsum += r * r;
      ++nsum;
    }
    dx = (double *) ((char *) dx + incx);
  }
  if (nsum > 1)
    xmax *= sqrt(dsum);
  return xmax;
}                               /* dnrm2 */

void cblas_drot(const int N, double *dx, const int incX,
                double *dy, const int incY, const double c, const double s)
{                               /* APPLIES A PLANE ROTATION. */
  double dtmpx, dtmpy;
  int n, incx, incy;
  n = N;
  if (n > 0) {
    incx = incX * (int) sizeof(double);
    incy = incY * (int) sizeof(double);
    dtmpx = (*dx);
    dtmpy = (*dy);
    *dx = dtmpx * c + dtmpy * s;
    *dy = dtmpy * c - dtmpx * s;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmpx = (*dx);
      dtmpy = (*dy);
      *dx = dtmpx * c + dtmpy * s;
      *dy = dtmpy * c - dtmpx * s;
    }
  }
}                               /* drot */

void cblas_drotg(double *sa, double *sb, double *c, double *s)
{
  double asa, asb, r, t;
  asa = (*sa);
  if (asa < 0.)
    asa = -asa;
  asb = (*sb);
  if (asb < 0.)
    asb = -asb;
  if (asa > asb) {              /* c > 0.707 */
    t = (*sb) / (*sa);
    r = sqrt(1. + t * t);
    *sa *= r;
    *c = 1. / r;
    *s = *sb = (*c) * t;
  } else if (asb != 0.) {       /* s > 0.707 */
    t = (*sa) / (*sb);
    r = sqrt(1. + t * t);
    *s = 1. / r;
    *c = (*s) * t;
    *sa = (*sb) * r;
    *sb = 1.;
    if ((*c) != 0.)
      *sb /= (*c);
  } else {
    *sa = *s = 0.;
    *sb = *c = 1.;
  }
}
void cblas_dscal(const int N, const double da, double *dx, const int incX)
{  /* SCALES A VECTOR BY A CONSTANT. LOOP UNROLLED */
  double dtmp0, dtmp, *dx0;
  int n, incx;
  n = N;
  if (n > 0 && da != 1.) {
    incx = incX * (int) sizeof(double);
    *dx *= da;
    if ((n & 1) == 0) {
      dx = (double *) ((char *) dx + incx);
      *dx *= da;
      --n;
    }
    n = n >> 1;
    while (n > 0) {
      dx0 = (double *) ((char *) dx + incx);
      dx = (double *) ((char *) dx0 + incx);
      dtmp0 = (*dx0);
      dtmp = (*dx);
      *dx0 = da * dtmp0;
      *dx = da * dtmp;
      --n;
    }
  }
}                               /* dscal */

void cblas_dswap(const int N, double *dx, const int incX,
                 double *dy, const int incY)
{  /* INTERCHANGES TWO VECTORS. */
  double dtmp;
  int n, incx, incy;
  n = N;
  if (n > 0) {
    incx = incX * (int) sizeof(double);
    incy = incY * (int) sizeof(double);
    dtmp = (*dx);
    *dx = (*dy);
    *dy = dtmp;
    while ((--n) > 0) {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmp = (*dx);
      *dx = (*dy);
      *dy = dtmp;
    }
  }
}  /* dswap */
