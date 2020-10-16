#include <math.h>
#include <stdlib.h>

int ftran(nftx, xpec, nxpec, fc, flg)
int nftx;
double *xpec, *fc;
int nxpec, flg;
{
  /* dimension of fc >= nft+1 if nft even */
  static double cnorm, cnorm2, *c, *s;
  static int nft, nhalf, nfc;
  double *pfc, *pfs;
  double cval, sval, hval, tmpm, tmp;
  int j, k, jj, jk, nx;
  if (nft != nftx) {            /* initialize cosines */
    if (nft)
      free(c);
    nft = nftx;
    c = (double *) malloc(2 * nft * sizeof(double));
    s = c + nft;
    nhalf = (nft + 1) >> 1;
    nfc = nft >> 1;
    c[0] = 1.;
    c[nhalf] = -1.;
    s[nhalf] = s[0] = 0.;
    cnorm = 1. / (double) nft;
    cnorm2 = 2 * cnorm;
    tmpm = cnorm2 * acos(-1.);
    for (k = 1; k < nhalf; ++k) {
      tmp = k * tmpm;
      c[k] = cos(tmp);
      c[nft - k] = c[k];
      s[k] = sin(tmp);
      s[nft - k] = -s[k];
    }
  }
  pfc = fc;
  switch (flg) {
  case 0:                       /* analyze */
    nx = nxpec;
    cval = xpec[0];
    sval = 0;
    for (k = 1; k < nhalf; ++k) {
      sval += xpec[nx];
      nx += nxpec;
      cval += xpec[nx];
      nx += nxpec;
    }
    hval = 0.;
    if ((nft & 1) == 0) {
      sval += xpec[nx];
      nx += nxpec;
      pfc[nft - 1] = (cval - sval) * cnorm;
      pfc[nft] = 0.;
      hval = xpec[nx >> 1];
    }
    *pfc = (cval + sval) * cnorm;
    for (k = 1; k < nhalf; ++k) {
      jj = 0;
      jk = 0;
      hval = -hval;
      cval = xpec[0] + hval;
      sval = 0.;
      for (j = 1; j < nhalf; ++j) {
        jj += nxpec;
        tmp = xpec[jj];
        tmpm = xpec[nx - jj];
        jk += k;
        if (jk >= nft)
          jk -= nft;
        cval += (tmp + tmpm) * c[jk];
        sval += (tmp - tmpm) * s[jk];
      }
      ++pfc;
      *pfc = cval * cnorm2;
      ++pfc;
      *pfc = sval * cnorm2;
    }
    return 0;
  case -1:                      /* synthesize from fc */
    cval = (*pfc);
    ++pfc;
    for (k = 1; k <= nfc; ++k) {
      cval += (*pfc);
      pfc += 2;
    }
    break;
  case 1:                       /* find derivitive from fc */
    *pfc = cval = 0.;
    for (k = 1; k <= nfc; ++k) {
      ++pfc;
      pfs = pfc;
      ++pfc;
      tmp = k * (*pfc);
      *pfc = -k * (*pfs);
      *pfs = tmp;
      cval += tmp;
    }
    break;
  case 2:                       /* find integral from fc */
    *pfc = cval = 0.;
    for (k = 1; k <= nfc; ++k) {
      ++pfc;
      pfs = pfc;
      ++pfc;
      tmp = -(*pfc) / k;
      *pfc = (*pfs) / k;
      *pfs = tmp;
      cval += tmp;
    }
    break;
  default:
    cval = 0.;
  }
  nx = nxpec * nft;
  xpec[0] = cval;
  jj = 0;
  for (j = 1; j <= nfc; ++j) {
    jj += nxpec;
    pfc = fc;
    cval = (*pfc);
    sval = 0.;
    jk = 0;
    for (k = 1; k <= nfc; ++k) {
      jk += j;
      if (jk >= nft)
        jk -= nft;
      ++pfc;
      cval += (*pfc) * c[jk];
      ++pfc;
      sval += (*pfc) * s[jk];
    }
    xpec[jj] = cval + sval;
    xpec[nx - jj] = cval - sval;
  }
  return 0;
}                               /* ftran */
