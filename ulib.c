/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */ 
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: fixed potential underflow condition in ordblk */
/*   22 August 2000: removed use of binary file in getvar, putvar */
/*   22 August 2000: improved code in p to recognize + - as delimiters */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "calpgm.h"

static double zero = 0.;
static double tiny = 1.e-36;
static double machep = 1e-15;

int ordblk(ndm, n, iqnsep, t, e, isblk, p, ipasgn)
const int ndm, n;
double *t, *e, *p;
short *iqnsep, *isblk, *ipasgn;
{
  static double one = 1.;
  double pbgn, pcmp, etmp, pmax, ptmp, *tbase, *tp;
  long ndml;
  int iend, k, iswap, iq, is, ns, iz, eigmax, eig, ibk, knt, iperm;
  int iend0, eiglast, nm;
  short isep, nsep, ksep;

  /*   subroutine to arrange eigenvectors and eigenvalues so the */
  /*               block projection matrix is dominant diagonal  */
  /*               and energies are ordered within block         */
  /*   on input -                                                */

  /*     t is the eigenvector matrix and e is the eigenvalue     */
  /*     ndm is the dimensioned length of t and                  */
  /*     n is the actual length                                  */
  /*     isblk contains the starting index of a block            */
  /*             (note: isblk(nblk+1)=nn+1   )                   */
  /*     where nblk is the number of blocks                      */
  /*     iord <0 meads reverse sort , iord = 0 means no sort     */
  /*   returns -                                                 */
  /*     sorted t and e                                          */
  /*     p is the projection of the state on its block           */
  /*     ipasgn is a scratch vector which gives the assignment   */
  /*                positive value means provisional assignment  */
  /*                negative value means final assignment        */
  if (n <= 1) {
    p[0] = 1.;
    return 0;
  }
  iperm = 0;
  ndml = ndm;
  nm = n - 1;
  isep = iqnsep[0]; nsep = 1; iqnsep[0] = nsep; 
  for (eig = 1; eig < n; ++eig) {
    ksep = iqnsep[eig];
    if (ksep != isep) {
      isep = ksep;  
      ++nsep;
    }
    iqnsep[eig] = nsep;
  }
  if ((int) isblk[1] >= n && nsep == 1) { 
    /* skip projection finding for only 1 block */
    dcopy(n, &one, 0, p, 1);
  } else {  /* find largest projection for each eigenvector */
    if (nsep > 1) { /* isolate subsets defined by hdiag */
      for (eig = 0; eig < n; ++eig){
        ipasgn[eig] = 0; p[eig] = -1.; 
      }
      iend = 0; eiglast = 0;
      for (isep = 1; isep <= nsep; ++isep) {
        iend0 = iend;
        for (iend = iend0; iend < n; ++iend) { /* find number in subset */
          if (iqnsep[iend] != isep) break;
        }
        /* calculate projections */
        for (eig = eiglast; eig < n; ++eig) {
          if (ipasgn[eig] != 0) continue;
          iq = -1; pcmp = -1.;
          k = iend0; tp = &t[eig * ndml]; 
          for (ksep = isep; ksep <= nsep; ++ksep) {
            ptmp = 0.;
            while (k < n && iqnsep[k] == ksep) {
              ptmp += tp[k] * tp[k]; 
              ++k;
            }
            if (ptmp > pcmp) {
              pcmp = ptmp; iq = (int)ksep;
              if (pcmp > 0.5) break;
            }
          }
          ipasgn[eig] = (short)iq; p[eig] = pcmp;
        }
        /* look for largest projection */
        k = iend0; knt = 0;
        for (ksep = isep; ksep <= nsep; ++ksep) {
          iq = 0;
          for( ; k < n && iqnsep[k] == ksep; ++k) {
            if (ipasgn[k] < 0 || iq < 0) continue;
            iq = -1; pcmp = -1.; 
            for (eig = eiglast; eig < n; ++eig) {
              if (ipasgn[eig] == ksep && p[eig] > pcmp) {
                pcmp = p[eig]; iq = eig; 
              }
            }
            if (iq < 0) continue; /* subset not filled when iq < 0 */
            if (iq != k) {  /* swap columns */
              etswap(ndm, n, k, iq, t, e, p);
              ipasgn[iq] = ipasgn[k];   
            }
            ipasgn[k] = -ksep; 
            while (eiglast < n && ipasgn[eiglast] < 0) ++eiglast;
          }
          if (iq < 0) continue;
          /* remove excess members of subset */
          for (eig = eiglast; eig < n; ++eig) { 
            if (ipasgn[eig] == ksep) {
              ipasgn[eig] = 0; ++knt;
            }
          }
        } /* ksep loop */
        if (knt == 0) break;
      } /* isep loop */
    } /* end of sep sorting */
    pmax = 0.;
    eigmax = 0;
    iq = 0;
    tbase = t;
    for (eig = 0; eig < n; ++eig) { /* first pass at finding projection */
      pbgn = 1.;
      pcmp = -1;
      for (ibk = 1; (is = isblk[ibk]) < n; ++ibk) {
        ns = isblk[ibk + 1] - is;
        tp = tbase + is;
        ptmp = ddot(ns, tp, 1, tp, 1);
        /* accumulate projection for first block specially */
        pbgn -= ptmp;
        if (ptmp > pcmp) {
          pcmp = ptmp + tiny;
          iq = ibk;
          /* check for automatic largest projection */
          if (pcmp >= 0.5) break;
        }
      }
      if (pbgn > pcmp) {
        pcmp = pbgn;
        iq = 0;
      }
      p[eig] = pcmp;
      ipasgn[eig] = (short) iq;
      if (pcmp >= pmax) {
        pmax = pcmp;
        eigmax = eig;
      }
      tbase += ndm;
    }
    knt = nm; iend = nm;
    while (knt > 0) {  /* search for unassigned states */
      ibk = ipasgn[eigmax];
      if (ibk < 0) break;
      is = isblk[ibk];
      iz = isblk[ibk + 1] - 1;
      if (is <= eigmax && eigmax <= iz) {
        /* current position is within block */
        ipasgn[eigmax] = -1; --knt; eiglast = -1;
      } else {
        /* find next available position in block */
        isep = iqnsep[eigmax]; eiglast = eigmax;
        for (iswap = iz; iswap >= is; --iswap) {
          if (ipasgn[iswap] >= 0 && iqnsep[iswap] == isep) {
            ++iperm; --knt; eiglast = -1;
            ipasgn[eigmax] = ipasgn[iswap]; ipasgn[iswap] = -1;
            etswap(ndm, n, eigmax, iswap, t, e, p);
            break;
          }
        }
      }
      /* find next unassigned state with largest projection */
      pmax = 0.;
      tbase = t;
      iend0 = -1;
      for (eig = 0; eig <= iend; ++eig, tbase += ndm) {
        if (ipasgn[eig] < 0) continue;
        iend0 = eig;
        if (eig == eiglast) {       /* recalculate projections */
          isep = iqnsep[eig];
          pcmp = -1; iq = -1;
          for (ibk = 0; (is = isblk[ibk]) < n; ++ibk) {
            iz = isblk[ibk + 1];
            while (iqnsep[is] != isep && is < iz) 
              ++is;
            ns = iz - is;
            if (ns > 0 && ipasgn[is] >= 0) {
              tp = tbase + is; 
              ptmp = ddot(ns, tp, 1, tp, 1);
              if (ptmp > pcmp) {
                pcmp = ptmp + tiny;
                iq = ibk;
              }                
            }
          }
          p[eig] = pcmp;
          ipasgn[eig] = (short) iq;
        }
        if (p[eig] >= pmax ) {
          pmax = p[eig];
          eigmax = eig;
        }
      }
      /*  update iend to point to last unassigned state */
      iend = iend0;
    }
  }
  /* energy order blocks by minimum search and swap */
  tbase = t;
  ibk = 1;
  iq = 0;
  eig = 0;
  iz = isblk[1];
  for (is = 1; is < n; ++is, tbase += ndm) {
    if (is == iz) {
      iz = isblk[++ibk];
    } else {                    /*  find min e */
      etmp = e[eig];
      iswap = eig;
      isep = iqnsep[eig];
      for (k = is; k < iz; ++k) {
        if (iqnsep[k] != isep) continue;
        ++iq;
        if (e[k] < etmp) {
          etmp = e[k];
          iswap = k;
        }
      }
      if (iswap > eig) {        /* perform swap with min */
        ++iperm;
        etswap(ndm, n, eig, iswap, t, e, p);
      }
    }
    eig = is;
  }
  /*  set phase for positive diagonal */
  tbase = t;
  for (eig = 0; eig < n; ++eig, tbase += ndm) {
    if (tbase[eig] < 0.) {
      ++iperm;
      for (k = 0; k < n; ++k)
        tbase[k] = -tbase[k];
    }
  }
  iperm &= 1;
  if (iq != 0) iperm += 2;
  return iperm;
}                               /* ordblk */
void etswap(ndm, nsize, ix1,ix2,t,e,q)
const int ndm, nsize, ix1, ix2;
double *t, *e, *q;
{
   double tmp;
   long ndml;
   if (ix1 == ix2) return;
   if (q != NULL) {
     tmp = q[ix1]; q[ix1] = q[ix2]; q[ix2] = tmp;
   }
   if (e != NULL) {
     tmp = e[ix1]; e[ix1] = e[ix2]; e[ix2] = tmp;
   }
   ndml = ndm;
   dswap(nsize, &t[ix1 * ndml], 1, &t[ix2 * ndml], 1);
}
int hdiag(ndm, n, z, d, e, iqsep)
const int ndm, n;
double *z, *d, *e;
short *iqsep;
{
  double f, g, t, bb, *zim, *zk, *zkk, *ebgn, *dbgn;
  int isgn, k, im, nk, nv, nm, ndiag, nz, nvz, nret;

  /*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,  */
  /*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. */
  /*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */
  /*     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A            */
  /*     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING             */
  /*     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                          */
  /*     THE SUBROUTINE THEN CALLS TRIAG WHICH DIAGONALIZES              */
  /*     ON INPUT-                                                       */
  /*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL       */
  /*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM        */
  /*          DIMENSION STATEMENT,                                       */
  /*        NX IS THE ORDER OF THE MATRIX,                               */
  /*        Z CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE        */
  /*          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.             */
  /*     ON CALL TO TRIAG -                                              */
  /*        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,  */
  /*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL       */
  /*          MATRIX IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY,     */
  /*        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX              */
  /*          PRODUCED IN THE REDUCTION.                                 */
  /*     ON OUTPUT-                                                      */
  /*        D CONTAINS THE EIGENVALUES IN ARBITRARY ORDER                */
  /*        Z CONTAINS THE EIGENVECTORS                                  */
  /*        RETURN VALUE IS SET TO                                       */
  /*          ZERO       FOR NORMAL RETURN,                              */
  /*          J          IF THE J-TH EIGENVALUE HAS NOT BEEN             */
  /*                     DETERMINED AFTER 30 ITERATIONS.                 */
  /*---------------------------------------------------------------------*/
  if (n <= 0 || ndm < n)
    return 0;
  ndiag = ndm + 1;
  isgn = 1;
  nm = n - 1;
  zim = &z[nm * (long) ndm];
  iqsep[nm] = (short) nm;
  e[nm] = (double) nm;
  for (im = nm; im >= 2; --im) {/*** FOR I=N STEP -1 UNTIL 3 DO -- ****/
    zim -= ndm;
    nv = im - 1;
    dcopy(im, &z[im], ndm, e, 1);
    k = (int) idamax(im, e, 1);
    iqsep[nv] = (short) k;
    f = e[k];
    if (k < nv) {
      if (k > 0)
        dswap(k, &z[k], ndm, &z[nv], ndm);
      zk = &z[k * (long) ndm];
      nz = nv  - k - 1;
      if (nz > 0)
        dswap(nz, &zk[k + 1], 1, &zk[nv + ndm], ndm);
      t = zk[k]; zk[k] = zim[nv]; zim[nv] = t;
      e[k] = zk[im] = zim[im]; e[nv] = zim[im] = f;
      isgn = -isgn;
    }
    g = dnrm2(im, e, 1);
    if (g > (fabs(f) + tiny)) {
      if (f > 0.)
        g = -g;
      t = f - g;
      bb = t / g;
      t = 1. / t;
      dscal(nv, t, e, 1);
      for (nz = 0; nz < nv; ++nz)
        if (e[nz] != zero)  break;
    } else {
      nz = nv;
      bb = 0.;
    }
    e[nv] = (double) nz;
    if (nz < nv) {
      dcopy(nv, e, 1, zim, 1);
      zim[im] = g;
      /********* FORM D = Z * E ***************/
      dcopy(im, &z[nv], ndm, d, 1);
      zkk = &z[nz];
      ebgn = &e[nz];
      nvz = nv - nz;
      for (k = 0; k < nz; ++k) {
        d[k] += ddot(nvz, zkk, 1, ebgn, 1);
        zkk += ndm;
      }
      for (k = nz; k < nv; ++k) {
        nk = nv - k;
        t = e[k];
        daxpy(nk, t, &zkk[1], 1, &d[k + 1], 1);
        d[k] += ddot(nk, zkk, 1, &e[k], 1);
        zkk += ndiag;
      }
      dscal(im, bb, d, 1);
      dbgn = &d[nz];
      /*********** FORM D * E ***************/
      t = 0.5 * bb * (d[nv] + ddot(nvz, ebgn, 1, dbgn, 1));
      daxpy(nvz, t, ebgn, 1, dbgn, 1);
      d[nv] += t;
      /*********** FORM REDUCED A ***********/
      zkk = &z[nz];
      for (k = 0; k < nz; ++k) {
        t = d[k];
        daxpy(nvz, t, ebgn, 1, zkk, 1);
        zkk[nvz] += t;
        zkk += ndm;
      }
      for (k = nz; k < nv; ++k) {
        nk = nv - k;
        t = (*dbgn);
        daxpy(nk, t, ebgn, 1, zkk, 1);
        zkk[nk] += t;
        ++nk;
        t = (*ebgn);
        daxpy(nk, t, dbgn, 1, zkk, 1);
        zkk += ndiag;
        ++ebgn;
        ++dbgn;
      }
      zkk[0] += 2. * d[nv];
      d[nv] = bb;
      isgn = -isgn;
    }
  }

  d[0] = z[0];
  z[0] = 1.;
  iqsep[0] = 0;
  if (n == 1)
    return 1;
  if (isgn < 0) { /* CHANGE SIGN TO MAKE HOUSEHOLDER REFLECTION INTO ROTATION */
    z[0] = -z[0];
    z[1] = -z[1];
  }
  /*********** ACCUMULATION OF TRANSFORMATION MATRICES ***********/
  nret = 0; iqsep[0] = 0;
  zim = z;
  for (nv = 1; nv < n; ++nv) {
    nz = (int) e[nv];
    f = zim[nv];
    zim += ndm;
    bb = d[nv];
    d[nv] = zim[nv];
    if (f == zero || fabs(f) < machep * (fabs(d[nv-1]) + fabs(d[nv])) ) {
      ++nret; f = zero;
    }
    e[nv] = f;
    if (nz < nv) {
      zk = &z[nz];
      zkk = &zim[nz];
      nvz = nv - nz;
      for (k = 0; k < nv; ++k) {
        g = bb * ddot(nvz, zkk, 1, zk, 1);
        daxpy(nvz, g, zkk, 1, zk, 1);
        zk[nvz] = g;
        zk += ndm;
      }
      dscal(nvz, bb, zkk, 1);
      zim[nv] = bb + 1.;
    } else {
      dcopy(nv, &zero, 0, zim, 1);
      dcopy(nv, &zero, 0, &z[nv], ndm);
      zim[nv] = 1.;
    }
    nz = iqsep[nv];
    if (nz < nv)
      dswap(nv + 1, &z[nz], ndm, &z[nv], ndm);
    iqsep[nv] = (short) nret;
  }
  ++nret;
  nv = triag(ndm, n, n, z, d, e);
  if (nv != 0) nret = -nv;
  return nret;
}                               /* hdiag */

int triag(ndm, n, nz, z, d, e)
const int ndm, n, nz;
double *z, *d, *e;
{
  /* Local variables */
  int i, k, m, nm, m0, itr;
  double b, c, f, g, p, r, s, *zindx1, *zindx2, *zindx0;

  /********************************************************************
   THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE imtql2,
        NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON.
   THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
                    OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.

   ON INPUT-
       NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
                  ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
                  DIMENSION STATEMENT,
       N  IS THE ORDER OF THE MATRIX,
       NZ IS THE ACTUAL COLUMN LENGTH OF Z,
       D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
       E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
                  E(0) IS ARBITRARY,
       Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
                  REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
                  OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
                  THE IDENTITY MATRIX.

   ON OUTPUT-
       D CONTAINS THE EIGENVALUES IN ARBITRARY ORDER.  IF AN
                  ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
                  UNORDERED FOR INDICES > IERR
       E HAS BEEN DESTROYED,
       Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
                  TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
       Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
                  EIGENVALUES,
   RETURNS
       ZERO       FOR NORMAL RETURN,
       J          IF THE J-TH EIGENVALUE HAS NOT BEEN
                       DETERMINED AFTER 30 ITERATIONS.
   MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
                  THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
  ****************************************************************************/

  nm = n - 1;
  m0 = 0;
  zindx0 = z;
  for (k = nm; k > 0; --k) {
    for (itr = 30; itr > 0; --itr) {
      g = fabs(d[k]);
      for (m = k; m > 0; --m) { /* LOOK FOR SMALL SUB-DIAGONAL ELEMENT  */
        b = e[m];
        if (b == zero)
          break;
        f = g;
        g = fabs(d[m - 1]);
        if (fabs(b) < machep * (f + g))
          break;
      }
      if (m == k)
        break;
      if (m != m0) {
        m0 = m;
        zindx0 = z;
        if (m != 0)
          zindx0 += m * (long) ndm;
      }
      zindx2 = zindx0;
      /*********** FORM SHIFT ***********/
      c = d[k];
      s = e[k];
      g = 0.5 * (d[k - 1] - c) / s;
      r = sqrt(g * g + 1);
      if (g < 0)
        r = -r;
      p = d[m];
      g = p - c + s / (g + r);
      /*********** QL TRANSFORMATION ***********/
      for (i = m; i < k; ++i) {
        b = f = e[i + 1];
        if (i > m) {
          b *= c;
          f *= s;
        }
        e[i] = r = sqrt(f * f + g * g);
        s = f / r;
        c = g / r;
        r = d[i + 1];
        f = c * r - s * b;
        g = c * b - s * p;
        r += p;
        p = c * f - s * g;
        g = s * f + c * g;
        d[i] = r - p;
        zindx1 = zindx2;
        zindx2 += ndm;
        drot(nz, zindx1, 1, zindx2, 1, c, s);
      }
      d[k] = p;
      e[k] = g;
      e[m] = zero;
    }
    if (itr == 0)
      break;                    /*  NO CONVERGENCE TO AN EIGENVALUE */
  }
  return k;
}                               /* triag */

int getpar(lu, luout, nfit, npar, idpar, par, erpar, plbl, plblen)
FILE *lu, *luout;
int *nfit, *npar, plblen;
bcd_t *idpar;
double *erpar, *par;
char *plbl;
{ /*  SUBROUTINE TO READ PARAMETERS */
#define NCARD 130
  char *ic, *tlbl, *sbcd;
  char card[NCARD];
  double vec[3], parbase, ermin, erdef;
  size_t lblen;
  int i, n, kfit, kk, ibcd, ndbcd, nd;

  ermin = 1.0e-37;
  erdef = 1.0e+37;
  lblen = (size_t) plblen;
  if (plblen == 0) {            /*  plblen=0 for CALCAT */
    erdef = ermin;
    lblen = strlen(plbl);
  }
  n = *npar;
  kfit = 0;
  for (i = 0; i < 30; ++i)
    fputc(' ', luout);
  fputs("PARAMETERS - A.PRIORI ERROR \n", luout);
  tlbl = plbl;
  parbase = 1.; ibcd = 0; ndbcd = (int) idpar[0]; nd = ndbcd + ndbcd;
  sbcd = (char *) mallocq ((size_t) nd);
  for (i = 0; i < n; ++i) {         /* loop over input lines */
    if (fgetstr(card, NCARD, lu) <= 0)
      break;                    /* end file */
    kk = getbcd(card, &idpar[ibcd], ndbcd);
    if (kk <= 0) 
      break;
    if (i == 0) 
      idpar[0] = (bcd_t) ndbcd;
    vec[0] = 0.;
    vec[1] = erdef;
    vec[2] = 2.;
    if (pcard(&card[kk], vec, 3, NULL) == 0)
      break;                    /* nothing in line */
    if (vec[2] < 1.01)
      break;                    /* reading into covariances */
    ic = strchr(card, '/');
    if (ic == NULL) {
      *tlbl = '\0';
    } else {
      memcpy(tlbl, &ic[1], lblen);
      tlbl[lblen] = '\0';
    }
    par[i] = vec[0];
    putbcd(sbcd, nd, &idpar[ibcd]);
    if (NEGBCD(idpar[ibcd]) == 0) 
      ++kfit;
    fprintf(luout,"%6d %6d %s %21.13E ", i+1, kfit, sbcd, vec[0]);
    if (NEGBCD(idpar[ibcd]) == 0) {
      if (vec[1] < ermin)
        vec[1] = ermin;
      fprintf(luout,"%15.6E %s\n", vec[1], tlbl);
      parbase = vec[0];
      if (fabs(parbase) < ermin) 
        parbase = ermin;
    } else {
      vec[0] = vec[0] / parbase;
      fprintf(luout,"%15.6f %s\n", vec[0], tlbl);
      vec[1] = ermin; par[i] = vec[0];
    }
    erpar[i] = vec[1];
    if (plblen != 0)
      tlbl += plblen;
    ibcd += ndbcd;
  }
  free(sbcd);
  *npar = i;
  *nfit = kfit;
  fprintf(luout, "%d parameters read, %d independent parameters\n", i,
          kfit);
  return (i - n);
}                               /* getpar */

int getvar(lu, nfit, var, idpar, erpar, iflg)
FILE *lu;
bcd_t *idpar;
double *var, *erpar;
const int nfit, iflg;
{                               /*  SUBROUTINE TO READ VARIANCE INTO PACKED UPPER TRIANGULAR */
#define NCARD 130
  static double dval[8] = { 2., 2., 2., 2., 2., 2., 2., 2. };
  double *pvar, *pcol, val;
  int i, j, k, n, nv, ibcd, ndbcd;
  char card[NCARD];

  if (var == NULL)
    return -1;
  n = 0;
  if (iflg == 0) {              /* number of requested parameters == actual */
    n = nfit;
    k = nv = 0;
    pvar = var;
    for (j = 0; j < n; ++j) {
      for (i = 0; i <= j; ++i) {
        if (k >= nv) {
          if (fgetstr(card, NCARD, lu) <= 0) {
            nv = 0;             /* end file */
          } else {
            nv = pcard(card, dval, 8, NULL);    /* get covariances */
          }
          k = 0;
        }
        if (k < nv && fabs(*pvar = dval[k]) < 1.01) {
          ++k;
          ++pvar;
        } else {
          j = n = 0;            /* break out of both loops */
        }
      }
    }
  }
  ibcd = 0; ndbcd = (int) idpar[0];
  if (n == 0) {                 /* MAKE DEFAULT VAR */
    pvar = var;
    k = 0;
    val = 0.;
    *pvar = erpar[0]; 
    for (j = 1; j < nfit; ++j) {
      ++k; ibcd += ndbcd;
      ++pvar;
      dcopy(j, &zero, 0, pvar, 1);
      pvar += j;
      while (NEGBCD(idpar[ibcd]) != 0) {
        ++k; ibcd += ndbcd;
      }
      *pvar = erpar[k];
    }
  } else {                      /* scale by erpar */
    pcol = var;
    k = 0;
    for (i = 0; i < nfit; ++i) {
      pcol += i;
      pvar = pcol;
      while (NEGBCD(idpar[ibcd]) != 0) {
        ++k; ibcd += ndbcd;
      }
      val = erpar[k];
      ++k; ibcd += ndbcd;
      for (j = i; j < nfit; ++j) {
        pvar += j;
        *pvar *= val;
      }
    }
  }
  return n;
}                               /* getvar */

int putvar(lu, nfit, var, erpar)
FILE *lu;
const int nfit;
double *var, *erpar;
{                               /*  SUBROUTINE TO WRITE VARIANCE FROM PACKED UPPER TRIANGULAR */
  double val, *pvar, *pcol;
  int i, j, ndiag, knt;
  ndiag = nfit + 1;
  pcol = var;
  for (i = 0; i < nfit; ++i) {
    pcol += i;
    pvar = pcol;
    val = 1. / (*erpar);
    ++erpar;
    for (j = i; j < nfit; ++j) {
      pvar += j;
      *pvar *= val;
    }
  }
  pvar = var;
  knt = 0;
  for (j = 0; j < nfit; ++j) {
    for (i = 0; i <= j; ++i) {
      fprintf(lu, "%10.7f", *pvar);
      knt = (knt + 1) & 7;
      if (knt == 0)
        fputc('\n', lu);
      ++pvar;
    }
  }
  if (knt != 0)
    fputc('\n', lu);
  return 0;
}                               /* putvar */

double calerr(nfit, var, derv)
const int nfit;
const double *var, *derv;
{
  double err, sum;
  int n;

  err = 0.;
  if (var != NULL) {
    for (n = 1; n <= nfit; ++n) {
      sum = ddot(n, derv, 1, var, 1);
      err += sum * sum;
      var += n;
    }
    err = sqrt(err);
  }
  return err;
}                               /* calerr */

int deflin(iqnfmt, idqn)
int iqnfmt;
short *idqn;
{
  short km, i, k, nqn;

  nqn = (short) (iqnfmt % 10);
  if (nqn == 0) nqn = 10;
  for (k = 0; k < 20; ++k) {
    idqn[k] = -1;
  }
  i = (short) ((iqnfmt / 100) % 5);
  km = 0;
  for (k = i; k < nqn; ++k) {
    idqn[k] = km;
    idqn[k + nqn] = (short) (km + nqn);
    km = k;
  }
  return nqn;
}                               /* deflin */

int getlin(lu, nqn, idqn, iqn, xfreq, xerr, xwt, card, ncard)
FILE *lu;
const int nqn, ncard;
short *idqn, *iqn;
double *xfreq, *xerr, *xwt;
char *card;
{                               /*  read in .LIN card from LU */
  static double terr = 1.e-07;
  static double txwt = 1.e-30;
  static double big = 1000.;
  static int fmt[20] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                         3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  double val[20];
  int j, jd, k, nn, nc;
  char ctmp;

  if (nqn <= 6) {
    nn = 12; nc = 36;
  } else {
    nn = nqn + nqn; nc = nn * 3;
  }
  if (ncard <= nc) 
    return -2;
  if (fgetstr(card, ncard, lu) < nc)
    return -1;
  dcopy(nn, &big, 0, val, 1);
  ctmp = card[nc]; card[nc] = '\0';
  pcard(card, val, nn, fmt);
  card[nc] = ctmp;
  for (k = 0; k < nn; ++k) {    /* FIXUP QUANTA */
    j = (int) val[k];
    if (j > 999) {
      j = 0;
      jd = idqn[k];
      if (jd >= 0)
        j = iqn[jd];
    }
    iqn[k] = (short) j;
  }
  val[0] = 0.;
  val[1] = 0.01;
  val[2] = 1.;
  if (pcard(card + nc, val, 3, NULL) == 0)
    return -1;
  *xfreq = val[0];
  *xerr = val[1];
  if (fabs(*xerr) < terr)
    *xerr = terr;
  *xwt = val[2];
  if (*xwt < txwt)
    *xwt = txwt;
  return nc;
}                               /* getlin */

int getbcd(const char *line, bcd_t *ivbcd, int nbcd)
{
  int k, kk, kbgn, kend, na;
  char c;
  na = nbcd;
  if (na < 0)
    na = -na;
  if (na > 127) 
    na = 0; 
  ivbcd[0] = (bcd_t) na;
  kbgn = -1;
  do { 
    /* ignore leading whitespace */
    c = line[++kbgn];
    if (c >= '0' && c <= '9') break;
    if (c ==',') 
      c = '\0';
  } while (c != '\0' && kbgn < 1024);
  kend = kbgn;
  if (c != '\0') {
    if (kbgn > 0 && line[kbgn - 1] == '-')
      ivbcd[0] |= 0x80;
    do {
      /* find first trailing whitespace */
      c = line[++kend];
      if (c < '0' || c > '9') break;
    } while (kend < 1024);
  }
  /* fill bcd vector */
  kk = kend;
  for (k = 1; k < na; ++k) {
    if (--kk < kbgn) break; 
    ivbcd[k] = (bcd_t)(line[kk] - '0');
    if (--kk >= kbgn)  
      ivbcd [k] |= (bcd_t) (line[kk] - '0') << 4;
  }
  for (; k < na; ++k) {
    /* fill with 0 */
    ivbcd[k] = (bcd_t) 0;
  }
  if (line[kend] == '.') {
    do { 
      /* ignore trailing decimal point and fraction */
      c = line[++kend];
      if (c < '0' || c > '9') break;
    } while (kend < 1024);
  }
  if (nbcd > 0) {
    kk = 0;
    do { 
      /* ignore trailing whitespace */
      c = line[kend];
      if (c == '\0' || c == '-'|| c == '.' || c == '/') break;
      if (c >= '0' && c <= '9') break;
      if (c == ',') { /* ignore first comma */
        if (kk != 0) break;
        ++kk;
      }
    } while (++kend <= 1024);
  }
  return kend;
} /* getbcd */
int putbcd(char *line, int nlen, const bcd_t *ivbcd)
{
  int k, first, i, kk, n;
  k = (int)(ivbcd[0] & 0x7f); 
  n = k + k;
  if (nlen < n)
    n = nlen;
  --n;
  line[0] = ' ';
  first = 0;
  for (kk = 1; kk < n; ++kk) {
    if ((kk & 1) != 0) {
      --k;
      i = (int)(ivbcd[k] >> 4) & 0x0f;
    } else {
      i = (int) ivbcd[k] & 0x0f;
    }
    if (i == 0) {
      if (first == 0) 
        line[kk] = ' ';
      else
        line[kk] = '0';
    } else {
      if (first == 0) 
        first = kk;
      line[kk] = (char) i + '0';
    }     
  }
  if (first == 0) {
    if (n < 2)
      n = 2;
    first = n - 1;
    line[first] = '0';
  } 
  if (NEGBCD(ivbcd[0]) != 0)
    line[first - 1] = '-';
  line[n] = '\0';
  return n;
} /* putbcd */

int bcd2i(bcd_t btmp)
{
  if (btmp > (bcd_t) 9)
    btmp -= (bcd_t) 6 * (btmp >> 4);
  return (int) btmp;
} /* bcd2i */
bcd_t i2bcd(int ival)
{
  if (ival > 9) {
    if (ival >= 99)
      ival = 0x99;
    else
      ival += 6 * (ival / 10);
  }
  return (bcd_t)ival;
} /* i2bcd */
