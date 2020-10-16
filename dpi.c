/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */ 

/*   PACKAGE FOR DOUBLET PI  */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "calpgm.h"
/* Common Declarations */

static double zero = 0.;
static int nvib, iwhole, isdgn, nqn;

int hamx(iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump)
const int iblk, nsize, npar;
const BOOL ifdump;
const bcd_t *idpar;
const double *par;
double *egy, *t, *dedp, *pmix;
{
#define NDPAR 27
#define NZPAR 11
  /*@owned@*/ static double *wp[NDPAR + 1];
  /*@owned@*/ static int *ixpar;
  /*@owned@*/ static short *jq;
  /*@owned@*/ static short *isblk;
  /*@owned@*/ static short *iqnsep;
  static double wz[NZPAR];
  static int ioff, nfit, init;
  double *w, *ww, *pt, *pd;
  double g, gg, q, qq, x, z, xjq, rx, ry, sqj, sqj0, egy0;
  double aa, ab, ba, bb, ai, rt, rt2, fg, dj1, dj2;
  size_t nl;
  int idgn, nsq, ipar, kpar, ierr, itmp, f, i, j, js, pty;
  int ii, jj, kk, kq, aji, ajf, aii, iv, ibcd, ndbcd;
  bcd_t ibtmp;
  BOOL ofdiag;
  
  ndbcd = (int) idpar[0];
  if (iblk == 0) {
    wz[0] = 1.;
    for (i = 0; i <= NDPAR; ++i) {
      if (init != 0 && wp[i] != NULL)
        free(wp[i]);
      wp[i] = NULL;
    }
    if (init != 0 && isblk != NULL) {
      free(isblk);
      isblk = NULL;
    }
    if (init != 0 && iqnsep != NULL) {
      free(iqnsep);
      iqnsep = NULL;
    }
    if (init != 0 && jq != NULL) {
      free(jq);
      jq = NULL;
    }
    if (init != 0 && ixpar != NULL) {
      free(ixpar);
      ixpar = NULL;
    }
    init = 1;
    if (npar <= 0)
      return 0;
    nl = (size_t) npar * sizeof(int);
    ixpar = (int *) mallocq(nl);
    j = nsize << 1;
    if (j < MAXQN) j = MAXQN;
    nl = (size_t) j * sizeof(short);
    jq = (short *) mallocq(nl);
    iqnsep = (short *) mallocq(nl);
    nl += sizeof(short);
    isblk = (short *) mallocq(nl);
    nl = (size_t) nsize;
    nl *= nl * sizeof(double);
    ioff = 0;
    nfit = npar; 
    for (j = 0, ibcd = 0; j < npar; ++j, ibcd += ndbcd) {
      if (NEGBCD(idpar[ibcd]) != 0) {
        --nfit;
      }
      if (nvib == 1) {
        iv = 0;
        ibtmp = idpar[ibcd + 1];
      } else {
        iv = bcd2i(idpar[ibcd + 1]);
        ibtmp = idpar[ibcd + 2];
      }
      kpar = bcd2i(ibtmp);
      if (kpar > NDPAR) {
        ixpar[j] = 100;
        continue;
      }
      ixpar[j] = (kpar << 7) + iv;
      if (kpar != 0 && wp[kpar] == NULL) 
        wp[kpar] = (double *) mallocq(nl);
      if (iv == 99)
        iv = 0;
      if (kpar == 1 && iv == 0 && par[j] < 0.)
        ioff = -1;
    }
    for (kk = 15; kk <= 19; ++kk) { /* fill in required elements */
	    jj = kk - 4;
      if (kk == 19) --jj;
      if (wp[kk] != NULL && wp[jj] == NULL) 
        wp[jj] = (double *) mallocq(nl);
    }
    for (kk = 23; kk <= 26; ++kk) { /* fill in required elements */
	    jj = kk - 3;
      if (kk == 26) --jj;
      if (wp[kk] != NULL && wp[jj] == NULL) 
        wp[jj] = (double *) mallocq(nl);
    }
    return nsize;
  }
  getqn(iblk, 1, MAXQN, jq, &idgn);
  f = idgn - 1;
  aii = isdgn - 1;
  js = jq[0];
  pty = jq[2];
  iv = jq[3];
  if (nvib == 1)
    iv = 0;
  jj = js;
  if (pty < 0)
    ++jj;
  for (i = 0; i < nsize; i += 2) {      /*  SET UP J QUANTUM NUMBER */
    jq[i] = (short) js;
    jq[i + 1] = jq[i];
    if ((jj & 1) != 0)
      jq[i] = (short) (-js);
    ++js;
    ++jj;
  }
  nsq = nsize * nsize;
  for (i = 0; i <= NDPAR; ++i) {        /* ZERO DERIVATIVE W */
    if ((w = wp[i]) != NULL)
      dcopy(nsq, &zero, 0, w, 1);
  }
  for (i = 0; i < nsize; i += 2) {      /* SET UP FINE STRUCTURE TERMS */
    ii = i + 1;
    xjq = (double) (jq[i] + 1);
    rx = (double) (jq[i] * jq[i]);
    ry = rx - 2.;
    sqj = rx - 0.25;
    x = rx - 1.;
    z = sqrt(x);
    for (kk = 0; kk <= 11; ++kk) {
      if (kk == 11) kk = 27;
      w = wp[kk];
      if (w == NULL)
        continue;
      ofdiag = FALSE;
      aa = bb = ab = 0.;
      switch (kk) {
      case 0:                   /* Energy */
        break;
      case 1:                   /* A */
        aa = (double) ioff;
        bb = aa + 1;
        break;
      case 2:                   /* AJ */
        aa = -rx;
        bb = ry;
        break;
      case 3:                   /* AH */
        aa = rx * rx + x;
        bb = -ry * ry - x;
        break;
      case 4:                   /* B+0.5*Q */
        aa = rx;
        bb = ry;
        ab = z;
        ofdiag = TRUE;
        break;
      case 27:
        ab = z;
        ofdiag = TRUE;
        break;
      case 5:                   /* D */
        aa = -rx * rx - x;
        bb = -ry * ry - x;
        ab = -2 * z * x;
        ofdiag = TRUE;
        break;
      case 6:                   /* H */
        aa = (rx * rx + x) * rx + 2. * x * x;
        bb = (ry * ry + x) * ry + 2. * x * x;
        ab = ((3. * x + 1.) * x - 1.) * z;
        ofdiag = TRUE;
        break;
      case 7:                   /* P */
        aa = 0.5 * xjq;
        ab = 0.25 * z;
        ofdiag = TRUE;
        break;
      case 8:                   /* Q */
        aa = jq[i] + 0.5;
        bb = 0.5;
        ab = 0.5 * z * jq[i];
        ofdiag = TRUE;
        break;
      case 9:                   /* PD */
        aa = 0.5 * xjq * sqj;
        ab = 0.5 * z * sqj;
        ofdiag = TRUE;
        break;
      case 10:                  /* QD */
        aa = 0.5 * xjq * xjq * sqj;
        bb = 0.5 * x * sqj;
        ab = 0.5 * z * xjq * sqj;
        ofdiag = TRUE;
      }
      if (ofdiag)
        w[ii + i * nsize] = ab;
      if (kk < NZPAR) {
        if (i == 0) wz[kk] = 0.5 * (aa + bb);
        aa -= wz[kk]; bb -= wz[kk];
      }
      w[i + i * nsize] = aa;
      w[ii + ii * nsize] = bb;
    }
  }
  /* CHI PARAMETERS (MAGNETIC INTERACTIONS) */
  rt2 = sqrt(2.);
  ai = 0.5 * aii;
  fg = ai * (ai + 1.) * (ai + ai + 1.);
  ba = 0.;
  for (i = 0; i < nsize; i += 2) {
    ii = i + 1;
    aji = (jq[ii] << 1) - 1;
    dj1 = (double) (aji + 1);
    sqj0 = (double)(jq[i] * jq[i] - 1);
    for (j = i; j < nsize; j += 2) {
      jj = j + 1;
      ajf = (jq[jj] << 1) - 1;
      dj2 = (double) (ajf + 1);
      if ((ajf - aji) > 2) break;
      sqj = 0.5 * (sqj0 + (double)(jq[i] * jq[i] - 1));
      g = sqrt(dj1 * dj2 * fg) * c6jj(aji, 2, ajf, aii, f, aii);
      itmp = aji + f + aii;
      if ((itmp & 2) != 0) g = -g;
      gg = g;
      if ((ajf & 2) != 0) gg = -gg;
      if (pty < 0) g = -g;
      for (kq = 11; kq <= 14; ++kq) {
        w = wp[kq];
        if (w == NULL) continue;
        ww = wp[kq + 4];
        switch (kq) {
        case 11: /* and 15 */
          aa = 2. * gg * c3jj(aji, 2, ajf, -1, 0, 1);
          w[j + i * nsize] = aa;
          if (ww != NULL) 
            ww[j + i * nsize] = sqj * aa;
          break;
        case 12: /* and 16 */
          aa = rt2 * g * c3jj(aji, 2, ajf, -1, 2, -1);
          w[j + i * nsize] = aa;
          if (ww != NULL) 
            ww[j + i * nsize] = sqj * aa;
          break;
        case 13: /* and 17 */
          bb = -gg * c3jj(aji, 2, ajf, -3, 0, 3) / 1.5;
          w[jj + ii * nsize] = bb;
          if (ww != NULL) 
            ww[jj + ii * nsize] = sqj * bb;
          break;
        case 14: /* and 18 and 19 */
          ab = rt2 * gg * c3jj(ajf, 2, aji, -3, 2, 1);
          w[jj + i * nsize] = ab;
          if (ww != NULL) 
            ww[jj + i * nsize] = sqj * ab;
          if (j != i) {
            ba = rt2 * gg * c3jj(aji, 2, ajf, -3, 2, 1);
            w[j + ii * nsize] = ba; 
            if (ww != NULL) 
              ww[j + ii * nsize] = sqj * ba;
          }
          ww = wp[19];
          if (ww == NULL) break;
          z = -0.5 * (double)(jq[i] + jq[j]);
          ww[jj + i * nsize] = z * ab;
          if (j != i) 
            ww[j + ii * nsize] = z * ba;
          break;
        }
      }
    }
  }
  /* ETA PARAMETERS (QUADRUPOLE INTERACTIONS) */
  if (aii >= 2) {
    fg = fg * (ai + 1.5) / (ai * ai * (ai - 0.5));
    rt = pty / sqrt(24.);
    for (i = 0; i < nsize; i += 2) {
      ii = i + 1;
      aji = (jq[ii] << 1) - 1;
      dj1 = (double) (aji + 1);
      sqj0 = (double)(jq[i] * jq[i] - 1);
      for (j = i; j < nsize; j += 2) {
        jj = j + 1;
        ajf = (jq[jj] << 1) - 1;
        dj2 = (double) (ajf + 1);
        if ((ajf - aji) > 4) break;
        sqj = 0.5 * (sqj0 + (double)(jq[i] * jq[i] - 1));
        q = sqrt(dj1 * dj2 * fg) * c6jj(aji, 4, ajf, aii, f, aii);
        itmp = aji + f + aii;
        if ((itmp & 2) != 0) q = -q;
        qq = q;
        if ((ajf & 2) != 0) qq = -qq;
        q *= rt;
        qq *= 0.25;
        for (kq = 20; kq <= 22; ++kq) {
          w = wp[kq];
          if (w == NULL) continue;
          ww = wp[kq + 3];
          switch (kq) {
          case 20: /* and 23 */
            aa = qq * c3jj(aji, 4, ajf, -1, 0, 1);
            w[j + i * nsize] = aa;
            if (ww != NULL) 
              ww[j + i * nsize] = sqj * aa;
            break;
          case 21: /* and 24 */
            bb = -qq * c3jj(aji, 4, ajf, -3, 0, 3);
            w[jj + ii * nsize] = bb; 
            if (ww != NULL) 
               ww[jj + ii * nsize] = sqj * bb;
            break;
          case 22: /* and 25 and 26 */
            ab = -q * c3jj(aji, 4, ajf, -1, 4, -3);
            w[jj + i * nsize] = ab;
            if (ww != NULL) 
              ww[jj + i * nsize] = sqj * ab;
            if (i != j) {
              ba = q * c3jj(aji, 4, ajf, -3, 4, -1);
              w[j + ii * nsize] = ba;
              if (ww != NULL) 
                ww[j + ii * nsize] = sqj * ba;
            }
            ww = wp[26];
            if (ww == NULL) break;
            z = -0.5 * (double)(jq[i] + jq[j]);
            ww[jj + i * nsize] = z * ab;
            if (j != i) 
              ww[j + ii * nsize] = z * ba;
          } /* case */
        } /* kq loop */
      } /* j loop */
    } /* i loop */
  }
  /*  BRING H TOGETHER */
  dcopy(nsq, &zero, 0, t, 1);
  egy0 = zero; kk = 0;
  for (ipar = 0, ibcd = 0; ipar < npar; ++ipar, ibcd += ndbcd) {
    if (NEGBCD(idpar[ibcd]) == 0) 
      kk = ipar;
    kpar = ixpar[ipar];
    ii = kpar & 127;
    kpar = kpar >> 7;
    if (ii == iv || ii == 99) {
      aa = par[ipar];
      if (kk != ipar)
        aa *= par[kk];
      if (kpar < NZPAR) {
        egy0 += aa * wz[kpar];
        if (kpar == 0)
          continue;
      }
      w = wp[kpar];
      daxpy(nsq, aa, w, 1, t, 1);
    }
  }
  if (ifdump) {                 /* special to dump Hamiltonian */
    for (i = 0; i < nsize; ++i) {
      t[i + i * nsize] += egy0;
      egy[i] = t[i + i * nsize];
      pmix[i] = 0.;
      dcopy(nfit, &zero, 0, &dedp[i], nsize);
    }
    return 0;
  }

  /* DIAGONALIZE AND ORDER EIGENVECTORS */

  for (i = 0; i < nsize; ++i) {
    isblk[i] = (short) i;
  }
  isblk[nsize] = (short) nsize;
  ierr = hdiag(nsize, nsize, t, egy, pmix, iqnsep);
  if (ierr < 0) {
    printf("diagonalization failure %d\n", ierr);
    exit(EXIT_FAILURE);
  }
  ordblk(nsize, nsize, iqnsep, t, egy, isblk, pmix, jq);
  for (i = 0; i < nsize; ++i)
    egy[i] += egy0;
  /*  TRANSFORM W TO DEDP */
  kk = 0;
  ab = 1.;
  pd = dedp;
  for (ipar = 0, ibcd = 0; ipar < npar; ++ipar, ibcd += ndbcd) {
    if (NEGBCD(idpar[ibcd]) == 0) {
      if (ipar != 0)
        pd += nsize;
      dcopy(nsize, &zero, 0, pd, 1);
      kk = ipar;
    } else {
      ab = par[ipar];
    }
    kpar = ixpar[ipar];
    ii = kpar & 127;
    kpar = kpar >> 7;
    if (ii == iv || ii == 99) {
      if (kpar < NZPAR) {          /* corrrect for constant part */
        bb = wz[kpar];
        if (ipar != kk)
          bb *= ab;
        for (i = 0; i < nsize; ++i)
          pd[i] += bb;
        if (kpar == 0)
          continue;
      }
      pt = t;
      w = wp[kpar];
      for (i = 0; i < nsize; ++i) {
        bb = pt[0] * ddot(nsize, w, 1, pt, 1);
        for (j = 1; j < nsize; ++j) {
          aa = ddot(j, &w[j], nsize, pt, 1) +
            ddot(nsize - j, &w[j + j * nsize], 1, &pt[j], 1);
          bb += aa * pt[j];
        }
        if (ipar != kk)
          bb *= ab;
        pd[i] += bb;
        pt += nsize;
      }
    }
  }
  return 0;
}                               /* hamx */


int intens(iblk, isiz, jblk, jsiz, ndip, idip, dip, s)
const int iblk, isiz, jblk, jsiz, ndip;
const bcd_t *idip;
const double *dip;
double *s;
{
  double fac, aji, ajj, omeg, val;
  int idgn, jdgn, i, j, k, fi, fj, ix, jx, idif, jji, jjj, ii;
  int iv, jv, kdip, ibcd;
  bcd_t ivb, jvb;

  /* START OF INTENSITY CALCULATION */
  if (((iblk + jblk) & 1) == 0)
    return 0;
  iv = (iblk - 1) >> 1;
  fi = iv / nvib;
  iv -= fi * nvib;
  jv = (jblk - 1) >> 1;
  fj = jv / nvib;
  jv -= fj * nvib;
  idif = fi - fj;
  if (idif > 1 || idif < -1)
    return 0;
  if (iv < jv) {
    i = iv; iv = jv; jv = i;
  }
  ivb = i2bcd(iv); jvb = i2bcd(jv);
  for (kdip = 0, ibcd = 1; kdip < ndip; ++kdip, ibcd += NDECDIP) {
    if (idip[ibcd] == ivb && idip[ibcd + 1] == jvb) break;
  }
  if (kdip >= ndip)
    return 0;
  idgn = fi + fi + iwhole;
  fi = idgn - 1;
  jdgn = fj + fj + iwhole;
  fj = jdgn - 1;
  fac = dip[kdip] * sqrt((double) (idgn * jdgn));
  if (fi != fj) fac = -fac;
  ix = idgn - isdgn;
  if (ix < 0)
    ix = -ix;
  jx = jdgn - isdgn;
  if (jx < 0)
    jx = -jx;
  jx += 2 - ix;
  ii = isdgn - 1;
  for (i = 0; i < isiz; ++i) {
    dcopy(jsiz, &zero, 0, &s[i], isiz);
    k = i & 1;
    omeg = k + 0.5;
    jji = ix + i - k;
    aji = 0.5 * jji;
    if (k == 0)
      fac = -fac;
    j = i - jx;
    if (j >= jsiz)
      continue;
    if (j >= 0) {
      jjj = jji - 2;
      val = fac * sqrt((aji - omeg) * (aji + omeg) / aji)
           * c6jj(jji, 2, jjj, fj, ii, fi);
      s[i + j * isiz] = val;
    }
    j += 2;
    if (j >= jsiz)
      continue;
    if (j >= 0 && jji > k) {
      val = fac * omeg * sqrt((aji + aji + 1.) / (aji * aji + aji))
           * c6jj(jji, 2, jji, fj, ii, fi);
      s[i + j * isiz] = val;
    }
    j += 2;
    if (j >= jsiz)
      continue;
    if (j >= 0) {
      jjj = jji + 2;
      ajj = aji + 1.;
      val = -fac * sqrt((ajj - omeg) * (ajj + omeg) / ajj)
           * c6jj(jji, 2, jjj, fj, ii, fi);
      s[i + j * isiz] = val;
    }
  }
  return idgn;
}                               /* intens */


int setint(lu, ifdiag, nsav, ndip, idip, isimag)
FILE *lu;
BOOL *ifdiag;
int *nsav, *isimag;
const int ndip;
bcd_t *idip;
{
  bcd_t iv,jv;
  int kdip, ibcd;
  *ifdiag = TRUE;
  *nsav = 1;
  if (nvib == 1 && ndip > 0) {
    idip[1] = (bcd_t) 0; idip[2] = (bcd_t) 0;
  }
  for (kdip = 0, ibcd = 1; kdip < ndip; ++kdip, ibcd += NDECDIP) {
    isimag[kdip] = 0;
    iv = idip[ibcd]; jv = idip[ibcd + 1];
    if (iv < jv) {
      idip[ibcd] = jv; idip[ibcd + 1] = iv;
    }
  }
  return 0;
}                               /* setint */

int getqn(iblk, indx, maxqn, iqn, idgn)
const int iblk, indx, maxqn;
short *iqn;
int *idgn;
{
  int kf, idif, tdgn, iv, ix;

  /* ..PACKAGE FOR DOUBLET PI */
  iv = (iblk - 1) >> 1;
  kf = iv / nvib;
  iv -= kf * nvib;
  tdgn = kf + kf + iwhole;
  if (indx == 0) {
    if (tdgn > isdgn)
      tdgn = isdgn;
    *idgn = tdgn + tdgn;
    return nqn;
  }
  ix = indx + 1;
  idif = tdgn - isdgn;
  if (idif < 0)
    idif = -idif;
  iqn[0] = (short) ((idif >> 1) + (ix >> 1));
  iqn[1] = (short) ((ix & 1) + 1);
  iqn[2] = (short) (((iblk & 1) != 0) ? -1 : 1);
  iqn[3] = (short) kf;
  if (nvib > 1) {
    iqn[4] = iqn[3];
    iqn[3] = (short) iv;
  }
  if (iqn[0] < iqn[1])
    tdgn = 0;
  *idgn = tdgn;
  return 1;
}                               /* getqn */


int setopt(luin, nfmt, itd, ndbcd, namfil)
FILE *luin;
int *nfmt, *itd, *ndbcd;
char *namfil;
{                               /* SET OPTIONS */
  char card[82];
  double val[2];

  cjjini();
  if (fgetstr(card, 82, luin) <= 0)
    return -1;
  val[1] = val[0] = 1.;
  pcard(card, val, 2, NULL);
  isdgn = (int) val[0];
  if (isdgn < 1)
    isdgn = 1;
  nvib = (int) val[1];
  if (nvib < 1 || nvib > 99)
    nvib = 1;
  nqn = 4;
  if (nvib > 1)
    ++nqn;
  strcpy(namfil, "dpi.nam");
    iwhole = (isdgn + 1) & 1;
  *itd = 2;
  *ndbcd = 3;
  *nfmt = 1;
  return 1;
}                   /* setopt */

int setfmt(iqnfmt, nfmt)
int *iqnfmt, nfmt;
{
  *iqnfmt = 804;
  if (iwhole == 0)
    *iqnfmt += 10;
  if (isdgn == 1)
    *iqnfmt = 803;
  if (nvib > 1)
    *iqnfmt += 1101;
  return nqn;
}

int setblk(lu, npar, idpar, par, nblkpf, negy)
FILE *lu;
const int npar;
int *nblkpf, *negy;
bcd_t *idpar;
const double *par;
{
  int nsiz;
  double ai;
  if (npar > 0) {
    *nblkpf = nvib + nvib;
    ai = 0.5 * (isdgn - 1);
    fprintf(lu, "doublet PI lines for I = %5.1f\n", ai);
    if (nvib > 1)
      fprintf(lu, "  and %d vibrations\n", nvib);
    nsiz = isdgn + isdgn;
    if (*negy > nsiz)
      *negy = nsiz;
  } else {
    nsiz = 0;
  }
  hamx(0, nsiz, npar, idpar, par, &ai, &ai, &ai, &ai, 0);
  nsiz = 1;
  if (nvib > 1)
    nsiz = 100;
  return nsiz;
}                               /* setblk */
