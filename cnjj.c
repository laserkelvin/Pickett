#include <math.h>
#include "cnjj.h"

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */

#define NFAC 1024
static double fac[NFAC];

int cjjini()
{ /**** set up log factorials ****/
  int nfac;
  double v;
  v = 0.;
  fac[0] = 0.;
  fac[1] = 0.;
  for (nfac = 2; nfac < NFAC; ++nfac) {
    v += log((double) nfac);
    fac[nfac] = v;
  }
  return 0;
} /* cjjini */

double c3jj(j1, j2, j3, m1, m2, m3)
const int j1, j2, j3, m1, m2, m3;
{ /*    3J COEFFICIENT CALCULATION */
  double trm, val;
  int jsum, k, k1, k2, jm1, jm2, jm3, lmin, la, lb, lc;
  if ((m1 + m2 + m3) != 0)
    return 0.;
  jsum = deljj(j1, j2, j3, &val);
  if (jsum < 0)
    return 0.;
  la = jsum - j3;
  if (m1 == 0 && m3 == 0) {     /*   CHECK FOR M1=M2=M3=0 CASE */
    if ((jsum & 1) != 0)
      return 0.;
    k = jsum >> 1;
    la = la >> 1;
    lb = (jsum - j2) >> 1;
    lc = (jsum - j1) >> 1;
    val = exp(0.5 * val + fac[k] - fac[la] - fac[lb] - fac[lc]);
    if ((k & 1) != 0)
      val = -val;
    return val;
  }

  k = j1 + m1;
  if (k < 0 || (k & 1) != 0)
    return 0.;
  jm1 = k >> 1;
  val += fac[jm1];
  jm1 = j1 - jm1;
  if (jm1 < 0)
    return 0.;
  val += fac[jm1];
  lb = jm1;
  k1 = jsum - j2 - jm1;         /* k1 = j3 - j2 + m1 */
  k = j2 - m2;
  if (k < 0 || (k & 1) != 0)
    return 0.;
  jm2 = k >> 1;
  val += fac[jm2];
  jm2 = j2 - jm2;
  if (jm2 < 0)
    return 0.;
  val += fac[jm2];
  lc = jm2;
  k2 = jsum - j1 - jm2;         /* k2 = j3 - j1 - m2 */
  jm3 = (j3 - m3) >> 1;
  if (jm3 < 0)
    return 0.;
  val += fac[jm3];
  jm3 = j3 - jm3;
  if (jm3 < 0)
    return 0.;
  val += fac[jm3];
  jsum += jm3 - j2;

  if (k2 < 0 && k2 < k1) {      /* summation starts at -k2 */
    k2 = -k2;
    k1 += k2;
    jsum += k2;                 /* implicit swap of k and k2 */
    la -= k2;
    lb -= k2;
    lc -= k2;
  } else if (k1 < 0) {          /* summation starts at -k1 */
    k1 = -k1;
    k2 += k1;
    jsum += k1;                 /* implicit swap of k and k1 */
    la -= k1;
    lb -= k1;
    lc -= k1;
  }
  val = exp(0.5 * val - fac[la] - fac[lb] - fac[lc] - fac[k1] - fac[k2]);
  if ((jsum & 1) != 0)
    val = -val;
  trm = val;
  k = 0;
  lmin = la;      /* lmin = min(la,lb,lc) */
  if (lmin > lb)
    lmin = lb;
  if (lmin > lc)
    lmin = lc;
  while (lmin > 0) {
    trm = -trm * (double) (la--) * (double) (lb--) * (double) (lc--) /
        ((double) (++k) * (double) (++k1) * (double) (++k2));
    val += trm;
    lmin--;
  }
  return val;
} /* c3jj */

double c6jj(j1, j2, j3, jj1, jj2, jj3)
const int j1, j2, j3, jj1, jj2, jj3;
{ /* 6J COEFFICIENT CALCULATION    USING INTEGERS = 2 J */
  double trm, val;
  int k, k1, k2, k3, k4, lmin, la, lb, lc;
  k1 = deljj(j1, j2, j3, &val);
  if (k1 < 0)
    return 0.;
  k2 = deljj(j1, jj2, jj3, &trm);
  if (k2 < 0)
    return 0.;
  la = k2 - j1;
  val += trm;
  k3 = deljj(jj1, j2, jj3, &trm);
  if (k3 < 0)
    return 0.;
  lb = k3 - j2;
  val += trm;
  k4 = deljj(jj1, jj2, j3, &trm);
  if (k4 < 0)
    return 0.;
  lc = k4 - j3;
  val += trm;
  k = k4;                       /* k = max(k1,k2,k3,k4) */
  if (k3 > k) {                 /* move max to k4 */
    k = k3;
    k3 = k4;
    k4 = k;
  }
  if (k2 > k) {                 /* move max to k4 */
    k = k2;
    k2 = k4;
    k4 = k;
  }
  if (k1 >= k) {                /* leave max in k1 */
    k = k1;
  } else {                      /* move max to k1 */
    k4 = k1;
    k1 -= k;
    la += k1;
    lb += k1;
    lc += k1;
  }
  k1 = 0;
  k2 = k - k2;
  k3 = k - k3;
  k4 = k - k4;
  val = exp(0.5 * val + fac[++k] - fac[la] - fac[lb] - fac[lc]
            - fac[k2] - fac[k3] - fac[k4]);
  if ((k & 1) == 0)
    val = -val;
  trm = val;
  lmin = la;   /* lmin = min(la,lb,lc) */
  if (lmin > lb)
    lmin = lb;
  if (lmin > lc)
    lmin = lc;
  while (lmin > 0) {
    trm = -trm * (double) (++k) * (double) (la--) * (double) (lb--) *
                 (double) (lc--) / ((double) (++k1) * (double) (++k2) *
                                   (double) (++k3) * (double) (++k4));
    val += trm;
    lmin--;
  }
  return val;
}  /* c6jj */

int deljj(jxa, jxb, jxc, val)
const int jxa, jxb, jxc;
double *val;
{
  int ia, ib, ic, isum;

/* JXA,JXB,JXC=J+J */

/* ON RETURN JRET=J1+J2+J3 , UNLESS TRIANGLE RELATIONS NOT OBEYED */
/*                                (IN WHICH CASE JRET < 0) */

/*     JRET IS ASSUMED TO BE LESS THAN 1023 */
  isum = jxa + jxb + jxc;
  *val = 0.;
  if ((isum & 1) != 0)
    return -1;
  isum = isum >> 1;
  ia = isum - jxa;
  if (ia < 0)
    return ia;
  ib = isum - jxb;
  if (ib < 0)
    return ib;
  ic = isum - jxc;
  if (ic < 0)
    return ic;
  *val = fac[ia] + fac[ib] + fac[ic] - fac[isum + 1];
  return isum;
} /* deljj */

double c9jj(jj)
int jj[];
{ /* 9J COEFFICIENT CALCULATION */
  double val;
  int *jjp;
  int jmin, jsum, j, lperm, kmin, kmax, kk;

/* FIND SMALLEST VALUE OF J */

  jsum = jj[0];
  jmin = jsum;
  kmin = 0;
  for (kk = 1; kk < 9; ++kk) {
    j = jj[kk];
    jsum += j;
    if (j < jmin) {
      jmin = j;
      kmin = kk;
    }
  }

  lperm = 0;
  if (kmin >= 3) {/* EXCHANGE COLUMN SO JMIN IS IN FIRST COLUMN */
    if (kmin < 6) {
      kmin -= 3;
      jjp = &jj[3];
    } else {
      kmin -= 6;
      jjp = &jj[6];
    }
    for (kk = 0; kk < 3; ++kk) {
      j = jj[kk];
      jj[kk] = jjp[kk];
      jjp[kk] = j;
    }
    lperm = jsum;
  }
  if (kmin > 0) { /* EXCHANGE ROW SO JMIN IS IN FIRST ROW AND COLUMN */
    jjp = &jj[kmin];
    for (kk = 0; kk < 9; kk += 3) {
      j = jj[kk];
      jj[kk] = jjp[kk];
      jjp[kk] = j;
    }
    lperm += jsum;
  }
  if (jmin == 0) { /* SPECIAL FOR  J=0 */
    if (jj[3] != jj[6])
      return 0.;
    if (jj[1] != jj[2])
      return 0.;
    kmin = jj[1] + 1;
    kmax = jj[3] + 1;
    val = (double) kmax;
    if (kmin != kmax)
      val = sqrt(kmin * val);
    val = c6jj(jj[8], jj[7], jj[3], jj[4], jj[5], jj[1]) / val;
    lperm += jj[1] + jj[3] + jj[5] + jj[7];
  } else {
    kmin = jj[8] - jj[0];
    kmax = jj[8] + jj[0];

    kk = jj[1] - jj[5];
    if (kk < 0)
      kk = -kk;
    if (kk > kmin)
      kmin = kk;
    kk = jj[1] + jj[5];
    if (kk < kmax)
      kmax = kk;
    kk = jj[3] - jj[7];
    if (kk < 0)
      kk = -kk;
    if (kk > kmin)
      kmin = kk;
    kk = jj[3] + jj[7];
    if (kk < kmax)
      kmax = kk;
    val = 0.;
    for (kk = kmin; kk <= kmax; kk += 2) {
      val += (double) (kk + 1) *
          c6jj(jj[0], jj[1], jj[2], jj[5], jj[8], kk) *
          c6jj(jj[3], jj[4], jj[5], jj[1], kk, jj[7]) *
          c6jj(jj[6], jj[7], jj[8], kk, jj[0], jj[3]);
    }
    lperm += kmin << 1;
  }
  if ((lperm & 2) != 0)
    val = -val;
  return val;
} /* c9jj */
