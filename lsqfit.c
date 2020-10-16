#include <math.h>
#include "lsqfit.h"
#include "cblas.h"
#define idamax cblas_idamax
#define dasum  cblas_dasum
#define daxpy  cblas_daxpy
#define dcopy  cblas_dcopy
#define ddot   cblas_ddot
#define dnrm2  cblas_dnrm2
#define drot   cblas_drot
#define drotg  cblas_drotg
#define dscal  cblas_dscal
#define dswap  cblas_dswap

static double zero = 0.;
static double tiny = 1.e-34;
static double tmach = 1.e-10;
int lsqfit(dk, ndm, nr, nvec, marqp, varv, 
           dkold, ediag, enorm, iperm)
int ndm, nvec, nr, *iperm;
double *dk, *dkold, *marqp, *varv, *ediag, *enorm;
{ /* function to do a least squares fit */
  /* nr is the number of rows (number of fitted parameters) */
  /* nvec is the number of solution vectors */
  /* dk contains the lower triangular matrix */
  /*    in the first nr rows  followed by and nvec rows of vectors */
  /* on entry, dk contains the lower triagular derivative matrix  */
  /* on exit,  dk contains the pseudo inverse and nvec solutions */
  /* lsqfit returns -1 if fit is diverging and 0 otherwise */
  /* marqp[0] contains the marquart-levenberg parameter */
  /* marqp[1] contains the trust region (if marqp[1] < 0., find maximum  */
  /*                                     trust region for marqp[0] ) */
  /* marqp[2] contains ratio of new to old trust region */
  /* on input to lsqfit, varv[0] contains the input variance */
  /* lsqfit sets varv[1] to varv[0] for the last good fit */
  /* lsqfit sets varv[2] to the predicted constrained variance change */
  /* lsqfit sets varv[3] to the predicted unconstrained variance change */
  /* lsqfit sets varv[4] to varv[0] - varv[1] */
  /* on output ediag contains the diagonal elements of the transformed */
  /*                 derivative matrix before constraints are applied */
  /* dkold is used to save non-dero parts of dk with dimension of at least */
  /*          nr * (nr + 1) / 2 + nvec * nr */
  static int nrest;
  double *pdk, *pdkold;
  double scale, val, slope, curv, dif;
  int k, nt, kk, nrank;
  nt = nr + nvec;
  if (marqp[1] > 0.) { /* update trust region */
    dif = varv[0] - varv[1];
    slope = varv[2];
    curv = dif + slope;
    if (nrest != 0) {
      scale = marqp[2];
      val = scale * varv[4];
      if (val > dif) {
        curv = scale * (val - dif);
        slope = scale * val - dif;
      }
    }
    if (curv < tiny) {
      scale = 0.382;
    } else {
      scale = 0.5 * slope / curv;
    }
    val = 0.0002 * varv[1];
    if (dif > val) { /* restore parameters */
      if (scale < 0.382) {
        scale = 0.382;
      } else if (scale > 0.5) {
        scale = 0.5;
      }
      varv[2] *= scale; varv[4] = dif; 
      marqp[1] *= scale;
      marqp[2] = scale;
      pdkold = dkold; pdk = dk;
      pdkold[nr] *= scale;
      dcopy(nt, pdkold, 1, pdk, 1);
      kk = ndm + 1;
      for (k = 1; k < nr; ++k) {
        pdk += kk; pdkold += nt; --nt;
        pdkold[nr - k] *= scale;
        dcopy(nt, pdkold, 1, pdk, 1);
      } 
      dscal(nr, scale, &dk[nr], ndm);
      ++nrest;
      return -1;
    }
    if (scale > 1. && marqp[0] > 0.){
      scale = 4.47; /* sqrt(20) */
      marqp[1] *= scale; marqp[2] = scale;
    } 
  }
  varv[1] = varv[0]; nrest = 0;
  /*  normalize parameters */
  for (k = 0; k < nr; ++k) {
    pdk = &dk[k]; kk = k + 1;
    val = dnrm2(kk, pdk, ndm);
    if (val < tiny) {
      ediag[k] = 0.; enorm[k] = 1.;
    } else {
      ediag[k] = 1.; enorm[k] = scale = 1. / val;
      dscal(kk, scale, pdk, ndm);
    }
  }
  /*  convert to lower triangle */
  nrank = dqrfac(dk, ndm, nr, nr, nvec, ediag, iperm);
  /* nrank is the rank of the derivative matrix */
  pdk = &dk[nr]; 
  varv[3] = ddot(nrank, pdk, ndm, pdk, ndm);
  /* solve for marquardt parameter */
  val = marqp[1]; marqp[2] = 1.;
  nt = lmpar(nr, nrank, nvec, dk, ndm, ediag, marqp, dkold);
  if (nt == 0) {
    marqp[1] = -1.; val = -1.;
    nrank = lmpar(nr, nrank, nvec, dk, ndm, ediag, marqp, dkold);
  }
  if (val > tiny) 
    marqp[2] = marqp[1] / val;
  varv[2] = 2. * ddot(nrank, pdk, ndm, pdk, ndm) + tiny;
  if (nr > nrank) { /* clear rank-deficient part */
    pdk = &dk[nrank * (ndm + 1)];
    nt = nr + nvec - nrank;
    for (k = nrank; k < nr; ++k) {
      dcopy(nt, &zero, 0, pdk, 1);
      pdk += ndm;
    }
  }
  dqrsolv(dk, ndm, nrank, nr, nvec, iperm);  /*  invert lower triangle */
  /* save fit matrix */
  pdkold = dkold; pdk = dk;
  scale = enorm[0]; ediag[0] = fabs(*pdkold);
  kk = ndm + 1; nt = nr + nvec;
  dscal(nt, scale, pdk, 1);
  dcopy(nt, pdk, 1, pdkold, 1);
  for (k = 1; k < nr; ++k) {
    pdk += kk; pdkold += nt; --nt;
    scale = enorm[k]; ediag[k] = fabs(*pdkold);
    dscal(nt, scale, pdk, 1);
    dcopy(nt, pdk, 1, pdkold, 1);
  }
  for (k = nr - 1; k >= 0; --k) {
    kk = iperm[k];
    if (kk != k) {
      scale = ediag[k]; ediag[k] = ediag[kk]; ediag[kk] = scale;
    }
  }
  return 0;
} /* lsqfit */

int jelim(t, vec, ndm, n, ns)
double *t, *vec;
int ndm, n, ns;
{ /*   do orthogonal transformations to rotate vector VEC into matrix T */
  /*   NDM  = dimensioned column length of T, (NDM >= N + NS) */
  /*   N    = number of columns of T */
  /*   NS   = number of solution vectors, (NS > 0) */
  double c, s, r, rr, as, ac, q, tv, *td;
  int nn, k;

  q = 1.; nn = n + ns;
  for (k = 0; k < n; ++k) {
    if (k != 0) 
      t += ndm;
    s = (*vec);
    as = fabs(s);
    c = (*t);
    ac = fabs(c);
    if (as > ac) {
      if (q < 1.) {
        dscal(nn, q, vec, 1);
        q = 1.;
      }
      dswap(nn, t, 1, vec, 1);
      as = ac; s = (*vec); c = (*t);
    }
    td = t; ++t; ++vec; --nn;
    if (as < tiny)
      continue;
    tv = s / c;
    s = -tv;
    daxpy(nn, s, t, 1, vec, 1);
    tv *= q;
    rr = 1. + tv * tv;
    tv *= (q / rr);
    daxpy(nn, tv, vec, 1, t, 1);
    if (rr > 1.) {
      r = sqrt(rr);
      q /= r;                   /* rr <= 2 */
      dscal(nn + 1, r, td, 1);
      if (q < tiny) {           /* requires n > 240 */
        dscal(nn, q, vec, 1);
        q = 1.;
      }
    }
  }
  return 0;
} /* jelim */

int dqrfac(t, ndm, m, n, ns, wk, iperm)
double *t, *wk;
int ndm, m, n, ns, *iperm;
{
  double diag, temp, tcmp, hcmp, sqdiag;
  double *pcol, *pt, *ptt, *pdiag, *pbgn;
  int j, k, nt, nwk, nwkm, swap, nrank;
  /* DQRFAC COMPUTES LOWER TRIANGULAR FACTORIZATION WITH PIVOTING
     T       REAL(NDM,MAX(M,N))

     NDM     INTEGER
     NDM IS THE LEADING DIMENSION OF THE ARRAY T (NDM >= N + NS).

     M       INTEGER
     M IS THE ROW LENGTH OF THE SYSTEM.

     N       INTEGER
     N IS THE COLUMN LENGTH OF THE SYSTEM.

     NS      INTEGER
     NS IS THE NUMBER OF SOLUTION VECTORS

     WK IS A REAL WORKING VECTOR OF LENGTH N 
        ON ENTRY WK IS SET TO SUM OF SQUARES OF ROWS OF T

     IPERM IS A WORKING VECTOR OF LENGTH N

     RETURN VALUE IS RANK OF SYSTEM

     ON RETURN

     T       FACTOR OF ORIGINAL MATRIX
   */
  nrank = n;
  if (m < n)
    nrank = m;
  nt = n + ns;
  hcmp = tmach; tcmp = 1.;
  pcol = pbgn = t;
  for (k = 0; k < nrank; ++k) {
    nwk = m - k;
    if (k > 0) {
      pt = pbgn + k;
      for (j = k; j < n; ++j) {
        wk[j] -= (*pt) * (*pt);
        ++pt;
      }
    }   
    tcmp = wk[k]; swap = k;
    for (j = k + 1; j < n; ++j) { /* get largest square */
      temp = wk[j];
      if (temp > tcmp) {
        tcmp = temp; swap = j;
      }
    }
    if (tcmp <= hcmp) {         /* recompute squares */
      tcmp = 0.; swap = k; 
      pt = pcol + k;
      for (j = k; j < n; ++j) {
        wk[j] = temp = ddot(nwk, pt, ndm, pt, ndm);
        if (temp > tcmp) {
          tcmp = temp; swap = j;
        }
        ++pt;
      }
      hcmp = tcmp * tmach;
    }
    iperm[k] = swap;
    if (swap > k) {             /* swap rows */
      wk[swap] = wk[k]; wk[k] = tcmp;
      dswap(m, &t[k], ndm, &t[swap], ndm);
    }
    pdiag = pcol + k;
    diag = (*pdiag);
    sqdiag = diag * diag;
    nwkm = nwk - 1; pt = pdiag + ndm; 
    tcmp = ddot(nwkm, pt, ndm, pt, ndm);
    if (tcmp > tiny) {        /* need to make lower triangular */
      sqdiag += tcmp;
      tcmp = sqrt(sqdiag);
      if (diag < 0.)
        tcmp = -tcmp;
      temp = diag + tcmp;
      *pdiag = temp;
      diag = -tcmp;
      tcmp = -1 / (temp * tcmp);
      ptt = pdiag;
      for (j = k + 1; j < nt; ++j) {
        ++ptt;
        temp = tcmp * ddot(nwk, pdiag, ndm, ptt, ndm);
        daxpy(nwk, temp, pdiag, ndm, ptt, ndm);
      }
      *pdiag = diag;
    }
    if (sqdiag < tiny) {
      nrank = k;
      break;
    }
    temp = 0.; /* zero rest of row */
    dcopy(nwkm, &temp, 0, pt, ndm);
    pbgn = pcol;
    pcol += ndm;
  }
  temp = 0.;
  nwk = m - nrank;
  for (k = nrank; k < n; ++k) { /* zero rest of matrix */
    wk[k] = 0.;
    iperm[k] = k;
    dcopy(nwk, &temp, 0, pcol + k, ndm);
  }
  return nrank;
} /* dqrfac */

int dqrsolv(t, ndm, nrank, n, ns, iperm)
double *t;
int ndm, nrank, n, ns, *iperm;
{
  double diag, temp, tcmp, *pcol, *pt, *pdiag, *pbgn;
  long ndml;
  int j, k, nwk, nt, swap;
  /* DQRSOLV COMPUTES THE INVERSE OF A PERMUTED  MATRIX WITH ROW SWAPPING
     ON ENTRY
     T       REAL(NDM,MAX(M,N))

     NDM     INTEGER
     NDM IS THE LEADING DIMENSION OF THE ARRAY T (NDM >= N + NS.

     NRANK   INTEGER
     NRANK IS THE RANK OF THE SYSTEM (NRANK <= N).

     N       INTEGER
     N IS THE COLUMN LENGTH OF THE SYSTEM.

     NS      INTEGER
     NS IS THE NUMBER OF SOLUTION VECTORS

     WK IS A REAL WORKING VECTOR OF LENGTH N

     IPERM IS A WORKING VECTOR OF LENGTH N

     ON RETURN

     T       INVERSE OF ORIGINAL MATRIX PLUS SOLUTIONS
   */
  ndml = ndm;
  nt = n + ns;
  pcol = &t[(nrank - 1) * ndml];
  nwk = nt - nrank;
  for (k = nrank - 1; k >= 0; --k) {    /* solve lower triangular */
    pdiag = pcol + k;
    diag = 1 / (*pdiag);
    *pdiag = diag;
    ++pdiag;
    pbgn = &t[k];
    dscal(nwk, diag, pdiag, 1);
    for (j = 0; j < k; ++j) {
      temp = -(*pbgn);
      *pbgn = temp * diag;
      daxpy(nwk, temp, pdiag, 1, pbgn + 1, 1);
      pbgn += ndm;
    }
    ++nwk;
    swap = iperm[k];
    if (swap > k) {             /* swap old columns */
      dswap(nt, pcol, 1, &t[swap * ndml], 1);
    }
    if (k == 0) break;
    pcol -= ndm;
  }
  nt = n - 1;
  pcol = &t[nt * ndml];
  for (k = nt; k > 0; --k) {
    tcmp = ddot(k, pcol, 1, pcol, 1);
    if (tcmp > tiny) {   /* make lower triangular */
      pdiag = pcol + k;
      temp = *pdiag;
      tcmp = sqrt(tcmp + temp * temp);
      if (temp < 0.)
        tcmp = -tcmp;
      temp += tcmp;
      *pdiag = temp;
      diag = -tcmp;
      tcmp = -1 / (temp * tcmp + tiny);
      pt = t;
      nwk = k + 1;
      for (j = 0; j < k; ++j) {
        temp = tcmp * ddot(nwk, pcol, 1, pt, 1);
        daxpy(nwk, temp, pcol, 1, pt, 1);
        pt += ndm;
      }
      *pdiag = diag;
    }
    pcol -= ndm;
  }
  return nrank;
}     /* dqrsolv */

int lmpar(nr, nc, ns, r, ldr, x, par, wk)
int nr, nc, ns, ldr;
double *r, *x, *par, *wk;
/* this subroutine is a paraphrase of netlib MINPACK routine lmpar.f */
/* diag is assumed to be unity and no further swapping attempted     */
/* it returns an optimum Marquardt Levenberg parameter in par[0]     */
/* the trust region is in par[1]                                     */
/* use upper triangle for storage, dimensions of r should be > n     */
{
  static double rerr = 0.1;
  static double tmin = 1.e-7;
  static double tmax = 1000.;
  double *pdiag, *pcol, *prow, *pdel;
  double dxnorm, gnorm, delta, val, sum, fp, fpx, parlo, parhi, parmin;
  int ndiag, nm, k, itr, initl, nt;

  if (ns < 0) return 0;
  ndiag = ldr + 1;
  pdel = &r[nr];
  initl = 1;
  delta = par[1];
  if (delta <= tmin || ns == 0) {
    initl = -1;
    if (par[0] <= 0.) {
      initl = 2; par[0] = 0.;
    }
  }
  prow = r;
  for (k = 0; k < nc; ++k) {     /* calculate gradient */
    wk[k] = ddot(k + 1, prow, ldr, pdel, ldr);
    ++prow;
  }
  gnorm = dnrm2(nc, wk, 1);
  /* save matrix */
  pdiag = r; pcol = wk; nt = nr + ns;
  dcopy(nt, pdiag, 1, pcol, 1);
  for (k = 1; k < nc; ++k) { 
    pdiag += ndiag; pcol += nt; --nt;
    dcopy(nt, pdiag, 1, pcol, 1);
  }
  parmin = tmin * fabs(*pdiag);
  if (parmin < tiny) parmin = tiny;
  parhi = parlo = fp = 1.;
  for (itr = 10; itr > 0; --itr) {
    if (initl <= 0 && par[0] < parmin) {
      initl = 2; par[0] = 0.;
    }
    if (initl > 0) {            /* Marquardt-Levenberg par = 0 */
      fpx = 0.;
    } else {
      if (par[0] > tmax) {
        par[0] = tmax; initl = -1;
      }
      pdiag = r;
      fpx = sqrt(par[0]); nt = nr + ns;
      for (k = 0; k < nc; ++k) {
        /* include Marquardt-Levenberg par */
        val = (*pdiag);
        *pdiag = sum = sqrt(par[0] + val * val);
        --nt;
        if (nt == 0) break;
        prow = pdiag + 1;
        dcopy(nt, prow, 1, x, 1);
        val /= sum;
        dscal(nt, val, prow, 1);
        val = fpx / sum;
        dscal(nt, val, x, 1);
        pdiag += ndiag;
        jelim(pdiag, x, ldr, nt - ns, ns);
      }
    }
    dcopy(nc, pdel, ldr, x, 1);
    nm = nc - 1;
    /* calculate solution */
    pdiag = &r[nm * ndiag];
    x[nm] /= (*pdiag);
    for (k = nm - 1; k >= 0; --k) { 
      pdiag -= ndiag;
      x[k] = (x[k] - ddot(nm - k, pdiag + 1, 1, &x[k + 1], 1)) / (*pdiag);
    }
    dxnorm = dnrm2(nc, x, 1);
    if (initl < 0 || initl > 1) { /* initialize delta */
      delta = dxnorm;
      if (delta == zero)
        delta = 1.;
      par[1] = delta;
      break;
    }
    fpx = fp;
    fp = dxnorm / delta - 1.;
    if (initl == 0) {
      if (fabs(fp) < rerr)
        break;
      if (parlo == zero && fp <= fpx && fpx < 0.)
        break;
    } else if (fp < rerr) {
      par[0] = 0.;
      break;
    }
    val = 1. / (dxnorm + tiny);
    dscal(nc, val, x, 1);
    pdiag = prow = r;
    x[0] = val = x[0] / (*pdiag);
    sum = val * val;
    for (k = 1; k < nc; ++k) {
      pdiag += ndiag;
      ++prow;
      x[k] = val = (x[k] - ddot(k, prow, ldr, x, 1)) / (*pdiag);
      sum += val * val;
    }
    val = fp / (sum + tiny);
    if (initl > 0) {
      initl = 0;
      parlo = val;
      parhi = gnorm / delta;
      if (nr > nc)
        parlo = 0.;
      if (parhi == zero) {
        val = rerr;
        if (val > delta)
          val = delta;
        parhi = tiny / val;
      }
      if (par[0] > parhi)
        par[0] = parhi;
    } else {
      /* restore matrix */
      pdiag = r; pcol = wk; nt = nr + ns;
      dcopy(nt, pcol, 1, pdiag, 1);
      for (k = 1; k < nc; ++k) { 
        pdiag += ndiag; pcol += nt; --nt;
        dcopy(nt, pcol, 1, pdiag, 1);
      }
      if (fp > 0.) {
        if (par[0] > parlo)
          parlo = par[0];
      } else {
        if (par[0] < parhi)
          parhi = par[0];
      }
      par[0] += val;
    }
    if (par[0] < parlo)
      par[0] = parlo;
  }
  if (itr == 0) return 0;
  if (nr > nc)
    dcopy(nr - nc, &zero, 0, &x[nc], 1);
  return nc;
} /* lmpar */
