/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 3 March 2004 */
/*   28 Sept.  2005: fix phase bug in getzitot for reduced matrix elements */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "calpgm.h"
#include "spinit.h"
#define ODD(i)  (((int)(i) & 1) != 0)
#define ODD2(i) (((int)(i) & 2) != 0)
#define MOD(i,n) ((int)(i) % (int)(n))

/*@null@*/ static ITMIX *itmix_head = NULL;
/*@null@*/ static SITOT *ithead = NULL;
static int fc_neqi;
static double zrfc[MAXITOT], zifc[MAXITOT];
static double zero = 0.;

/* spin coupling pairs, npair = (neqi - 1)
       value >= npair is spin_n with n = (value - npair) = 0,.., npair
       value < npair is intermediate spin, Itot is last intermediate spin
*/  
static int itpairs[] = {1,2,                      /* neqi = 2, npair = 1 */
                        3,4, 0,2,                 /* neqi = 3, npair = 2 */
                        4,6, 3,5, 0,1,            /* neqi = 4, npair = 3 */
                        5,8, 6,7, 0,1, 2,4,       /* neqi = 5, npair = 4 */
                        6,10,7,9, 5,8, 0,1, 3,2}; /* neqi = 6, npair = 5 */
/*                      ---  ---  ---- ---  ---                          */
/*                       0    1    2    3    4                           */
static int *itpair = itpairs;

double ctrans(int i1, int i2, int ii, short *qn, int nqn, 
                      /*@null@*/ /*@out@*/ short *swk);

void getzitot(double *zfac, int ltot, int ii, const int *lscomv, 
              const int *iscomv, const int* jscomv, int alpha, int neqi)
{ /* get reduced matrix elements of nuclear spins */
  /* called from tensor subroutine */
  SITOT *pitot;
  ITMIX *pret;
  EITMIX *ipnow, *jpnow, *itbase;
  short *qnj, *qni;
  double val;
  int i, itot, isym, ni;
  int j, jtot, jsym, nj;
  int nm1, nm2, nn, k, isgn, lv1, lv2, itdif, ioff;
  if (ii <= 0 || neqi < 2) return;
  nm1 = neqi - 1; nm2 = nm1 - 1;
  lv1 = lscomv[0]; lv2 = lscomv[1];
  if (lv1 == 0) {
    if (lv2 == 0) /* no dependence on Itot */
      return;
    lv1 = lv2; lv2 = 0;
  }
  /* search for desired calculation */
  for (pitot = ithead; pitot != NULL; pitot = pitot->next) {
    if (pitot->lv1  != lv1)  continue;
    if (pitot->lv2  != lv2)  continue;
    if (pitot->ltot != ltot) continue;
    if (pitot->ii   != ii)   continue;
    if (pitot->neqi != neqi) continue;
    break;
  }
  if (pitot == NULL) {
    *zfac = 0.;
    return;
  }
  itot = (int)((unsigned int) iscomv[nm2] >> 1); isym = iscomv[0];
  jtot = (int)((unsigned int) jscomv[nm2] >> 1); jsym = jscomv[0];
  itdif = itot - jtot; isgn = 0;  
  if (itdif < 0) {
    itdif = -itdif; isgn = itdif;
    i = itot; itot = jtot; jtot = i;
    i = isym; isym = jsym; jsym = i;
  }
  itbase = NULL; nn = -1;
  pret = pitot->pret; 
  if (pret != NULL) {
    itbase = pret->mix; nn = pret->nitot;
  }
  ioff = -1;
  if (itot < nn && itdif <= (ltot >> 1)) 
    ioff = pitot->ioff[itdif * nn + itot];
  if (itbase == NULL || ioff < 0) {
    *zfac = 0.;
    return;
  }
  ipnow = itbase + itot; ni = ipnow->n;
  jpnow = itbase + jtot; nj = jpnow->n;
  if (ni == 1 && nj == 1) { /* only A,A */
    val = pitot->val[ioff];
    if (lscomv[0] == 0) /* second of two equivalent spins */ 
      isgn += itdif;
  } else{  
    qni = ipnow->qsym; qnj = jpnow->qsym;
    i = j = -1;
    for (k = 0; k < ni; ++k) {
      if ((int)qni[k] == isym) {
        i = k; break;
      }
    }
    for (k = 0; k < nj; ++k) {
      if ((int)qnj[k] == jsym) {
        j = k; break;
      }
    }
    if (i < 0 || j < 0) {
      *zfac = 0.; 
      return; /* no match */
    }
    if (alpha >= neqi) { /* pseudo-diagonal */
      if (isgn == 0) 
        i += ipnow->offset[i];
      else
        j += jpnow->offset[j];
    }
    val = pitot->val[i + j * ni + ioff];
  }
  if (ODD(isgn)) val = -val;
  *zfac *= val;
} /* getzitot */

void setzitot(int lv1, int lv2, int lvt, int ii, int neqi)
{ /* calculate reduced matrix elements of nuclear spins */
  /* used by getzitot */
  /* lv1 = tensor order of first spin */
  /* lv2 = tensor order of second spin */
  /* lvt = tensor order of product */
  /* ii = I * 2 */
  SITOT *pitot;
  ITMIX *pret;
  EITMIX *ipnow, *jpnow, *itbase;
  double *twk, *vec, *vec0, *vec_bgn, *twkodd, *eigvec;
  short *qnj, *qni;
  double val, val1, val2, tfac, rthalf;
  int *iscomv, *jscomv, *lscomv;
  int i, ni, nieven, ia, j, nj, njeven, ja, kee, koo, ioff;
  int jdif, itdif, nn, ltot, kblk, itot, itmax, jtot;
  int sii1, sii2, nm1, nm2, nt, nodd, k, ifhalf;
  if (ii <= 0 || neqi < 2) return;
  lv1 = lv1 + lv1; lv2 = lv2 + lv2; ltot = lvt + lvt;
  if (lv2 > lv1) {
    i = lv1; lv1 = lv2; lv2 = i;
  }
  /* search if calculation already done */
  for (pitot = ithead; pitot != NULL; pitot = pitot->next) {
    if (pitot->lv1  != lv1)  continue;
    if (pitot->lv2  != lv2)  continue;
    if (pitot->ltot != ltot) continue;
    if (pitot->ii   != ii)   continue;
    if (pitot->neqi != neqi)   continue;
    break;
  }
  if (pitot != NULL) return;
  /* check for algorithm limits */
  if (lv1 <= 0 || lvt < 0) return;
  if (lv2 == 0) { /* single spin */
    if (ltot != lv1) return;
  } else {/* paired spins */
    if (lv1 != lv2 || lvt > lv2 || ODD(lvt)) return;
  }
  /* search for spin eigenvalue sets */
  for (pret = itmix_head; pret != NULL; pret = pret->next) { 
    if (pret->ii == ii && pret->neqi == neqi) break;
  }
  if (pret == NULL) return;
  /* allocate container object for calculation */
  pitot = (SITOT *)mallocq(sizeof(SITOT));
  pitot->next = ithead; ithead = pitot;
  pitot->ii = ii; pitot->ltot = ltot; 
  pitot->lv1 = lv1; pitot->lv2 = lv2;
  pitot->neqi = neqi;
  /* find dimensions: */
  /* nn = total number of matrix elements, kblk = number of matrices */
  itmax = pret->nitot; kblk = itmax * (lvt + 1);
  pitot->ioff = (int *)mallocq((size_t)kblk * sizeof(int));
  pitot->ioff[0] = -1;
  ifhalf = (ii * neqi) & 1;
  nn = 0; itbase = pret->mix; 
  for (itdif = 0; itdif <= lvt; ++itdif){
    for(itot = 0; itot < itmax; ++itot){
      k = itdif * itmax + itot;
      pitot->ioff[k] = -1;
      jtot = itot - itdif;
      if (jtot < 0 || itot + jtot + ifhalf < lvt) continue;
      ipnow = itbase + itot; jpnow = itbase + jtot;
      ni = ipnow->n; nj = jpnow->n;
      if (ni == 0 || nj == 0) continue;
      pitot->ioff[k] = nn;
      nn += ni * nj;
    }
  }
  ipnow = NULL; jpnow = NULL;
  pitot->pret = pret;
  pitot->val = (double *) mallocq((size_t)nn * sizeof(double));
  pitot->val[0] = 0.;
  dcopy(nn, &zero, 0, pitot->val, 1);
  /* set up scratch vector */
  twk = (double *)mallocq(pret->nditot * sizeof(double));
  twk[0] = 0.;
  nm1 = neqi - 1; nt = neqi + nm1; nm2 = nm1 - 1;  
  /* allocate space for spin tensor quanta */
  iscomv = (int *)mallocq((size_t)(3 * nt) * sizeof(int));
  jscomv = &iscomv[nt]; lscomv = &jscomv[nt];
  for (k = nm1; k < nt; ++k) { 
    iscomv[k] = ii; jscomv[k] = ii;
  }
  rthalf = sqrt(0.5);
  tfac = 1.;
  if (neqi >= 3) {
    tfac = (double) neqi;
    if (lv2 > 0) tfac *= 0.5 * (neqi - 1);
  }
  if (lv2 > 0 && ltot == 4) {/* paired spins */ 
    tfac *= -sqrt(0.6); /* correct for order of tensor product */
  }
  kblk = 0; ioff = 0; ni = 0; nj = 0;
  for (itdif = 0; itdif <= lvt; ++itdif){
    for(itot = itdif; itot < itmax; ++itot){
      /* loop over itot and jtot */
      ioff = pitot->ioff[itdif * itmax + itot];
      if (ioff < 0) continue;
      jtot = itot - itdif;
      ipnow = itbase + itot; jpnow = itbase + jtot;
      ni = ipnow->n; nj = jpnow->n;
      vec0 = &pitot->val[ioff];
      nieven = ipnow->neven; njeven = jpnow->neven;
      /* start create */
      sii1 = 0; sii2 = -1;
      if (lv2 > 0) { /* 2 spins */
        if (neqi >= 3) sii1 = 1;
        sii2 = nm1;
      } 
      for (k = nm1; k < nt; ++k) { 
        lscomv[k] = 0;
      }
      lscomv[nm1 + sii1] = lv1;
      if (sii2 > 0) 
        lscomv[nm1 + sii2] = lv2;
      for (k = 0; k < nm1; ++k) {
        i = itpair[k + k]; j = itpair[k + k + 1];
        if (lscomv[i] == 0) {
          lscomv[k] = lscomv[j];
        } else if (lscomv[j] == 0) {
          lscomv[k] = lscomv[i];
        } else {
          lscomv[k] = ltot;
        }
      }
      vec = vec0;
      for (j = 0; j < nj; ++j) {
        if (j > 0) 
          vec += ni;
        qnj = &jpnow->qn[nm1 * j];
        for (k = 0; k < nm1; ++k) 
          jscomv[k] = qnj[k];
        for (i = 0; i < ni; ++i) {
          qni = &ipnow->qn[nm1 * i];
	  for (k = 0; k < nm2; ++k) { /* check intermediate spins */
	    iscomv[k] = qni[k];
	    jdif = qni[k] - qnj[k];
	    if (jdif < 0)
	      jdif = -jdif;
	    if (jdif > lscomv[k] || (int)(qni[k] + qnj[k]) < lscomv[k]) 
	      break;
	  }
	  if (k == nm2) {
	    iscomv[k] = qni[k]; val = tfac;
	    tensor(&val, iscomv, jscomv, lscomv, itpair, nm1, -1); 
	    vec[i] += val;
	  }
	}
      }
      if (ni > 1) { 
	twkodd = &twk[nieven]; nodd = ni - nieven;
	vec_bgn = vec0; 
	for (j = 0; j < nj; ++j) { /* left transformation */
	  if (j != 0) vec_bgn += nj;
	  vec = vec_bgn;
	  dcopy(ni, vec, 1, twk, 1); 
	  eigvec = ipnow->eigvec;
	  if (nieven > 1) {
	    *vec = ddot(nieven, twk, 1, eigvec, 1);
	    for (i = 1; i < nieven; ++i) {
	      vec += 1; eigvec += nieven;
	      *vec = ddot(nieven, twk, 1, eigvec, 1);
	    }
	  }
	  if (nodd > 0) {
	    if (nieven > 0) {
	      vec += 1; eigvec += nieven;
	    }
	    *vec = ddot(nodd, twkodd, 1, eigvec, 1);
	    for (i = 1; i < nodd; ++i) {
	      vec += 1; eigvec += nodd;
	      *vec = ddot(nodd, twkodd, 1, eigvec, 1);
	    }
	  }
	}
      }    
      if (nj > 1) { 
	twkodd = &twk[njeven]; nodd = nj - njeven;
	vec_bgn = vec0; 
	for (i = 0; i < ni; ++i) { /* right transformation */
	  if (i != 0) vec_bgn += 1;
	  vec = vec_bgn;
	  dcopy(nj, vec, ni, twk, 1);
	  eigvec = jpnow->eigvec;
	  if (njeven > 1) {
	    *vec = ddot(njeven, twk, 1, eigvec, 1);
	    for (j = 1; j < njeven; ++j) {
	      vec += ni; eigvec += njeven;
	      *vec = ddot(njeven, twk, 1, eigvec, 1);
	    }
	  }
	  if (nodd > 0) {
	    if (njeven > 0) {
	      vec += ni; eigvec += njeven;
	    }
	    *vec = ddot(nodd, twkodd, 1, eigvec, 1);
	    for (j = 1; j < nodd; ++j) {
	      vec += ni; eigvec += nodd;
	      *vec = ddot(nodd, twkodd, 1, eigvec, 1);
	    }
	  }
	}
      }
      twkodd = NULL;
      /* convert to complex basis */
      vec = vec0;
      for (i = 0; i < ni; ++i) {
	ia = ipnow->offset[i]; 
	if (ia < 0) continue;
	for (j = 0; j < nj; ++j) {
	  ja = jpnow->offset[j]; 
	  if (ja < 0) continue;
	  if (ia == 0) {
	    if (ja == 0) continue;
	    kee = i + j * ni; koo = kee + ja * ni;
	    if (i < nieven) { /* A1,E  */
	      val = rthalf * vec[kee]; vec[kee] = val; vec[koo] =  val;
	    } else { /* A2,E (implicit i) */
	      val = rthalf * vec[koo]; vec[kee] = val; vec[koo] = -val;  
	    }
	  } else if (ja == 0) { 
	    kee = i + j * ni; koo = kee + ia;
	    if (j < njeven) {/* E,A1 */
	      val = rthalf * vec[kee]; vec[kee] =  val; vec[koo] = val;
	    } else {/* E,A2 (implicit i) */
              val = rthalf * vec[koo]; vec[kee] = -val; vec[koo] = val;  
            }
          } else {/* E,E (ia > 0, ja > 0) */
            kee = i + j * ni; k = kee + ja * ni; koo = k + ia; 
            val1 = vec[kee]; val2 = vec[koo]; 
            val = 0.5 * (val1 + val2); 
            vec[kee] = val; vec[koo] = val;  
            val = 0.5 * (val1 - val2); 
            vec[kee + ia] = val; vec[k] = val;
          }
        }
      }
      vec = NULL;
    } /* ITOT loop */
  } /* ITOT - JTOT loop */
  lscomv = NULL; jscomv = NULL; free(iscomv); 
  free(twk);
} /* setzitot */

ITMIX *get_itmix(ii, neqi)
const int ii, neqi;
{
  ITMIX *pret;
  EITMIX *pitmix;
  double *eigval, *twk, *eigvec, *eiga, *eigb;
  short  *swk, *qn, *qn0, *ix, *ixx, *iqsep;
  double val, sum, sumsq, cfac;
  size_t nl;
  int itot, nt, ns, nm1, nm2, nitot, nodd, neven, i, k, nx, ntot;
  int ia, ib, ioff, kk, j, nq, neigv, ibgn, iend, ktot, kcmp;

  for (pret = itmix_head; pret != NULL; pret = pret->next) {
    if (pret->ii == ii && pret->neqi == neqi) 
      return pret;
  }
  if (neqi < 2 || ii <= 0) return NULL;
  /* allocate container and its elements */
  pret = (ITMIX *)mallocq(sizeof(ITMIX));
  pret->next = itmix_head; itmix_head = pret; pret->ii = ii;
  pret->neqi = neqi;
  iend = neqi * ii; ibgn = iend & 1; 
  nq = iend - ibgn + 2; ntot = nq >> 1; pret->nitot = ntot; 
  pret->mix = (EITMIX *)mallocq((size_t) ntot * sizeof(EITMIX));
  if (fc_neqi != neqi) {
    zrfc[0] = 1.; zifc[0] = 0.;
    sum = 2. * acos(-1.) / neqi;
    zrfc[0] = 1.; zifc[0] = 0.;
    for (i = 1; i < neqi; ++i) {
      zrfc[i] = cos(i * sum); zifc[i] = sin(i * sum);
    }
    fc_neqi = neqi;
  } 
  nm1 = neqi - 1; nm2 = nm1 - 1; nl = (size_t) nm1 * sizeof(short); 
  itpair = &itpairs[nm2 * nm2 + nm2];
  swk = (short *)mallocq((size_t)(neqi << 2) * sizeof(short));
  ixx = swk; ix = &ixx[nm1]; ix[0] = 0;
  for (k = 0; k < neqi; ++k) {
    ix[k + nm1] = (short)ii;
  }
  /* find all intermediate spin combinations */
  k = 0; ix[0] = (short)0; ixx[0] = (short)(ii + ii); ntot = 0;
  do {
    while(k < nm2) {
      ++k;
      ia = itpair[k + k]; ib = itpair[k + k + 1];
      ia = ix[ia]; ib = ix[ib];
      ixx[k] = (short) (ia + ib);
      ia -= ib;
      if (ia < 0) 
        ia = -ia;
      ix[k] = (short) ia;
    }
    ntot += (int)(ixx[k] - ix[k] + 2) >> 1; 
    ix[k] = ixx[k];
    while (k > 0 && ix[k] >= ixx[k])
      --k;
    ix[k] += 2;
  } while (ix[0] <= ixx[0]);
  pret->qnv = (short *)mallocq((size_t)(ntot * (nm1 + 2)) * sizeof(short));
  pret->qnv[0] = (short) 0; qn = pret->qnv;  qn0 = qn;
  pitmix = itmix_head->mix; pitmix->neven = 0; pitmix->n = 0; pitmix->qn = qn; 
  pitmix->offset = qn; pitmix->qsym = qn;
  ns = 1; neigv = 0; neven = 0; nodd = ntot; itot = ibgn; kcmp = 0;
  for (ktot = 0; ktot < nq; ++ktot){
    k = 0; ix[0] = (short)0; ixx[0] = (short)(ii + ii); nt = 0; nx = 0; 
    if (nodd > 0) {
      do {
        while(k < nm2) {
          ++k;
          ia = itpair[k + k]; ib = itpair[k + k + 1];
          ia = ix[ia]; ib = ix[ib];
          ixx[k] = (short) (ia + ib);
          ia -= ib;
          if (ia < 0) 
            ia = -ia;
          ix[k] = (short) ia;
        }
        i = ix[k];
        if(i == itot) {
          /* check perpendicular C2 symmetry */
          if (neqi < 5) 
            i = ii + ii - ix[0];
          else
            i = ix[0] + ix[1];
          if ((i & 2) == kcmp) {
            memcpy(qn, ix, nl); 
            qn += nm1; ++nt;
            if (nt >= nodd) break; 
          } else {
            ++nx;
          }
          ix[k] = ixx[k];
        } else if (i > itot) {
          ix[k] = ixx[k];
        } else if (itot >= 2) {
          ix[k] = (short) (itot - 2);
        }
        while (k > 0 && ix[k] >= ixx[k])
          --k;
        ix[k] += 2;
      } while (ix[0] <= ixx[0]);
    }
    if (nt > 0) {
      if (ns < nt) 
        ns = nt;
      neigv += nt * nt;
    }
    if (kcmp == 0) {
      neven = nt; pitmix->neven = neven; kcmp = 2; nodd = nx;
    } else {
      nodd = nt; nitot = nodd + neven; pitmix->n = nitot;
      pitmix->qn = qn0; 
      pitmix->offset = qn; qn += nitot;  
      pitmix->qsym = qn; qn += nitot;
      for (k = 0; k < neven; ++k) {
        pitmix->qsym[k] = (short) 0; pitmix->offset[k] = (short) 0;
      }
      for (k = neven; k < nitot; ++k) {
        pitmix->qsym[k] = (short) 2; pitmix->offset[k] = (short) 0;
      }
      ++pitmix; qn0 = qn; itot += 2; kcmp = 0; 
      ntot -= nitot; nodd = ntot;
    }
  }
  nl = (size_t)(ns + ns); itmix_head->nditot = nl; 
  qn = NULL; qn0 = NULL; ixx = NULL; ix = NULL; 
  itmix_head->eigvecv = (double *)mallocq((size_t)neigv * sizeof(double));
  eigvec = itmix_head->eigvecv; pitmix = itmix_head->mix;
  *eigvec = 1.; pitmix->eigvec = eigvec; 
  if (neqi == 2) {
    for (k = 1; k < neigv; ++k) {
      ++eigvec; ++pitmix;
      *eigvec = 1.; pitmix->eigvec = eigvec; 
    }
    free (swk);
    return itmix_head;
  }
  eigval = (double *)mallocq(nl * sizeof(double));
  iqsep  = (short *) mallocq(nl * sizeof(short));
  eigval[0] = 0.;
  twk = &eigval[ns];
  ns = neqi >> 1; nitot = 0;
  for (ktot = 0; ktot < nq; ++ktot) { 
    if (ODD(ktot)) {
      ioff = pitmix->neven; nitot = pitmix->n - ioff;
      eigvec += ioff * ioff; qn = &pitmix->qn[ioff * nm1];
    } else {
      if (ktot > 0) {
        eigvec += nitot * nitot; ++pitmix;
        pitmix->eigvec = eigvec;
      }
      ioff = 0; nitot = pitmix->neven; qn = pitmix->qn;
    }
    if (nitot == 0) continue;
    for (i = 0; i < nitot; ++i) {
      for (k = 0; k <= i; ++k) { 
        val = ctrans(i, k, ii, qn, neqi, swk); 
        eigvec[i + k * nitot] = val;
        if (k == i) continue;
        val -= ctrans(k, i, ii, qn, neqi, swk);
        if (fabs(val) > 1.e-10) {
          puts(" center Itot matrix not symmetric");
          exit(EXIT_FAILURE);
        }
      }
    }
    if (nitot == 1) {
      ia = 0; eigval[0] = eigvec[0]; eigvec[0] = 1.;
    } else {
      ia = hdiag(nitot, nitot, eigvec, eigval, twk, iqsep);
    }
    nt = neqi << 2;
    for (k = 0; k < nitot; ++k) {
      val = eigval[k] + 0.1;
      for (i = 0; i < ns; ++i) {
        if (val > zrfc[i]) break;
      }
      kk = k + ioff;
      pitmix->qsym[kk] += (short)(i << 2);
      ia = pitmix->qsym[kk];
      for (i = kk - 1; i >= ioff; --i) {
        ib = pitmix->qsym[i];
        if (MOD(ib, nt) == ia) {
          pitmix->qsym[kk] = (short)(ib + nt);
          break;
        }
      }
      if (ioff == 0) continue;
      ia = ia >> 2;
      if (ia == 0 || (ia + ia) == neqi) continue;
      ib = pitmix->qsym[kk] - 2;
      for (i = 0; i < ioff; ++i) {
        if((int)pitmix->qsym[i] == ib) {
          pitmix->offset[i]  = (short)(kk - i); 
          pitmix->offset[kk] = (short)(i - kk);
          break;
        }
      }
      if (i == ioff) continue;
      qn0 = pitmix->qn; cfac = 1. / zifc[ia];
      eiga = &pitmix->eigvec[i * ioff]; 
      eigb = &eigvec[k * nitot]; sumsq = 0.;
      for (j = 0; j < nitot; ++j) {
        sum = 0.;
        for (i = 0; i < ioff; ++i) {
          val = ctrans(i, j + ioff, ii, qn0, neqi, swk); 
          sum += val * eiga[i];
          val += ctrans(j + ioff, i, ii, qn0, neqi, swk);
          if (fabs(val) > 1.e-10) {
            puts(" non-center Itot matrix is not antisymmetric");
            exit(EXIT_FAILURE);
          }
        }
        sum *= cfac;
        eigb[j] = sum; sumsq += sum * sum;
      }
      // printf("Itot %3d, %10.6f\n", (int)qn[nm2], sumsq);
    }
    if (ODD(ktot)) { /* change qsym */
      for (k = 0; k < nitot; ++k) {
        kk = k + ioff;
        if (pitmix->offset[kk] < 0) {
          ia = MOD(pitmix->qsym[kk], nt) >> 2;
          pitmix->qsym[kk] += (short)(((neqi - ia - ia) << 2) - 2);
        }
      }
    }
  }
  twk = NULL;
  free(eigval); free(swk); free(iqsep);
  return itmix_head;
} /* get_itmix */

double ctrans(i1, i2, ii, qn, nqsym, swk)
int i1, i2, ii, nqsym;
short *qn, *swk;
{
  short *qn1, *qn2, *iperm, *mv, *m;
  double val, ele, fac;
  int i, nm1, nm2, k, mt, isgn, kk, ia, ib, imin, ma, mb, itot;
  nm1 = nqsym - 1; nm2 = nm1 - 1;
  qn1 = &qn[i1 * nm1];
  if (i1 == i2) 
    qn2 = qn1;
  else  
    qn2 = &qn[i2 * nm1];
  val = 0.; itot = qn1[nm2];
  if (nqsym == 3) {
    val = (double)(qn1[0] + 1);
    if (i1 != i2) 
      val = sqrt(val * (qn2[0] + 1));
    if (ODD2(qn2[0])) val = -val;
    val *= c6jj(ii, ii, (int)qn1[0], ii, itot, (int)qn2[0]);
  } else if (nqsym == 4)  {
    if (qn1[0] == qn2[1] && qn1[1] == qn2[0]) {
      val = 1.; isgn = ii + ii + qn1[1] - itot;
      if (ODD2(isgn)) val = -val;
    }
  } else { /* nqsym = 5,6 */
    if (swk == NULL) return val;
    imin = -ii;
    iperm = swk; mv = &swk[nqsym + nqsym]; m = &mv[nm1];
    for (i = 0; i < nm1; ++i) {
      iperm[i] = (short)i;
      iperm[i + nm1] = (short)(i + nqsym);
    }
    iperm[nm1 + nm1] = (short) nm1;
    k = -1; 
    do { /* loop over m */
      while(k < nm2) {
        ++k;
        m[k] = (short) imin;
      }
      mt = itot;
      for (i = 0; i < nm1; ++i)
        mt -= m[i];
      m[nm1] = (short) mt;
      while (mt >= imin && mt <= ii) { 
        for (kk = 0; kk < nm1; ++kk) {
          ia = itpair[kk + kk]; ib = itpair[kk + kk + 1];
          ma = mv[ia]; mb = mv[ib]; mt = ma + mb; 
          i = qn1[kk]; 
          if (mt < -i || mt > i) break;
          mv[kk] = (short) mt;
          if (ma == 0 && mb == 0) {
            ia = (ia >= nm1)? ii: qn1[ia]; i += ia;
            ib = (ib >= nm1)? ii: qn1[ib]; i += ib;
            if (ODD2(i)) break;
          }
        }
        if (kk < nm1) break;
        for (kk = 0; kk < nm1; ++kk) {
          ia = itpair[kk + kk]; ib = itpair[kk + kk + 1];
          ia = iperm[ia]; ib = iperm[ib];
          ma = mv[ia]; mb = mv[ib]; mt = ma + mb; 
          i = qn2[kk];
          if (mt < -i || mt > i) break;
          mv[kk] = (short) mt;
          if (ma == 0 && mb == 0) {
            ia = (ia >= nm1)? ii: qn2[ia]; i += ia;
            ib = (ib >= nm1)? ii: qn2[ib]; i += ib;
            if (ODD2(i)) break;
          }
        }
        if (kk < nm1) break;
        ele = 1.; fac = 1.; isgn = 0; 
        for (kk = 0; kk < nm1; ++kk) {
          ia = itpair[kk + kk]; ib = itpair[kk + kk + 1];
          ma = mv[ia]; mb = mv[ib]; mt = ma + mb; 
          ia = (ia >= nm1)? ii: qn1[ia];
          ib = (ib >= nm1)? ii: qn1[ib];
          i = qn1[kk]; mv[kk] = (short) mt;
          ele *= c3jj(ia, ib, i, ma, mb, -mt);
          fac *= (i + 1); isgn += ia - ib + mt;
        }
        for (kk = 0; kk < nm1; ++kk) {
          ia = itpair[kk + kk]; ib = itpair[kk + kk + 1];
          ia = iperm[ia]; ib = iperm[ib];
          ma = mv[ia]; mb = mv[ib]; mt = ma + mb; 
          ia = (ia >= nm1)? ii: qn2[ia];
          ib = (ib >= nm1)? ii: qn2[ib];
          i = qn2[kk]; mv[kk] = (short) mt;
          ele *= c3jj(ia, ib, i, ma, mb, -mt);         
          fac *= (i + 1); isgn += ia - ib + mt ;
        }
        if (ODD2(isgn)) ele = - ele;
        val += ele * sqrt(fac);
        break;
      }
      while (k > 0 && (int) m[k] >= ii)
        --k;
      m[k] += 2;
    } while (k > 0 || (int) m[0] <= ii);
    iperm = NULL; m = NULL; mv = NULL;
  }
  return val;
} /* ctrans */
