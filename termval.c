#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "calpgm.h"

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
/*   Revised version in c, 22 March 1999 */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */ 

/**************************************************************************/
#define NDQN 10
typedef struct egys {
  /*@notnull@*/ /*@owned@*/ /*@reldef@*/ struct egys *next;
  /*@notnull@*/ /*@dependent@*/ struct egys *prev;
  double egyx, errx;
  int nfreq, iorg;
  short iqn[NDQN];
} EGYS;
typedef struct line { /* lines in arbitrary order */
  /*@null@*/ /*@owned@*/ /*@reldef@*/ struct line *next;
  /*@notnull@*/ /*@dependent@*/ struct egys *pupper, *plower;
  double freqx, errfq;
  int seq;
} LINE;

static EGYS head_egy;
static LINE head_line;

static EGYS *lposn(short *ikey, const int nkey, int *stat);
static void emove(EGYS *pold, EGYS *pend);

int main(argc, argv)
int argc;
char *argv[];
{
#define NFILE 3
#define NDLINE 128
  static const char *cext[NFILE] = {"lin","trm","lgy"};
  enum efile {elin, etrm, elgy};
  static short idqn[NDQN + NDQN], iqnv[NDQN + NDQN];
  static char *cfil[NFILE+1];
  static char tmpstr[NDLINE];
  FILE *lutrm, *lulin, *luegy;
  EGYS *pegy, *pegyup, *pegylo, *pegyorg;
  LINE *last_line, *bgn_line, *crit_line, *pline, *tmpline;
  int *tagv;
  short *pqn;
  double xerr, clight, xfrq, xfrq0, big, xwt, cdif, egyorg, egy; 
  double sum, xerr1, xsq, xsq0, dif, dtmp, tiny, xerr0;
  size_t ndline;
  int k, i, nloop, icmp, nqn, norg, nline, oldorg, nfreq, kloop;
  int ndtag, ktag, critag;

  clight = 29979.2458;
  tiny = 1.e-15;
  big = 9999.9999;
  filget(argc, argv, NFILE, cfil, cext);
  rqexit(-1);
  lulin = fopenq(cfil[elin], "r");
  lutrm = fopenq(cfil[etrm], "w");
  luegy = fopenq(cfil[elgy], "w");
  fputs("enter number of quanta per state: ", stdout);
  fflush(stdout);
  fgetstr(tmpstr, 82, stdin);
  xwt = 3;
  pcard(tmpstr, &xwt, 1, NULL);
  nqn = (int) xwt;
  if (nqn > NDQN)
    nqn = NDQN;
  if (nqn <= 0)
    nqn = 1;
  k = 400 + nqn;
  if (nqn < 4) k = nqn * 101;
  nqn = deflin(k, idqn);
  /* initialize dummy energy */
  head_egy.next = &head_egy;
  head_egy.egyx = head_egy.errx = 0.;
  head_egy.iorg = head_egy.nfreq = 0;
  pqn = head_egy.iqn;
  pqn[0] = -1;
  for (k = 1; k < NDQN; ++k)
    pqn[k] = 0;
  /* initialize dummy line */
  head_line.next = NULL; 
  head_line.plower = head_line.pupper = &head_egy;
  head_line.freqx = 0.; head_line.errfq = tiny;
  head_line.seq = 0;
  ndline = sizeof(LINE);
  norg = nline = 0;
  /*  read and sort experimental lines */
  while (getlin(lulin, nqn, idqn, iqnv, &xfrq, &xerr, &xwt,
                tmpstr, NDLINE) > 0) {
    if (xerr > 0.) {
      xfrq /= clight;
      xerr /= clight;
    } else {
      xerr = -xerr;
    }
    if (xerr < tiny) xerr = tiny;
    if (iqnv[0] < 0) continue;
    icmp = 0;
    pegyup = lposn(iqnv, nqn, &icmp);
    icmp = icmp + icmp;
    pegylo = lposn(iqnv + nqn, nqn, &icmp);
    if (pegyup == pegylo) continue;
    pline = (LINE *) mallocq(ndline);
    pline->next = head_line.next; head_line.next = pline;
    pline->plower = pegylo; pline->pupper = pegyup;
    pline->freqx = xfrq; pline->errfq = xerr;
    pline->seq = (++nline);
    pegylo->nfreq += 1; pegyup->nfreq += 1;
    /* set up index for coupled states */
    switch (icmp) {
    case 0:                     /* both energies old */
      i = pegyup->iorg; k = pegylo->iorg;
      if (i == k) break;
      if (i > k) {
        icmp = i; i = k; k = icmp;
      }                         /* make i < k */
      --norg;
      for (pegy = head_egy.next; pegy != &head_egy; pegy = pegy->next) {
        if (pegy->iorg >= k) {
          if (pegy->iorg == k)
            pegy->iorg = i;
          else
            pegy->iorg -= 1;
        }
      }
      break;
    case 1:                     /* lower energy new */
      pegylo->iorg = pegyup->iorg;
      break;
    case 2:                     /* upper energy new */
      pegyup->iorg = pegylo->iorg;
      break;
    case 3:                     /* both energies new */
      norg++;
      pegylo->iorg = pegyup->iorg = norg;
      break;
    }
    fprintf(lutrm, " %5d old: ", nline);
    pqn = pegyup->iqn;
    for (i = 0; i < nqn; ++i)
      fprintf(lutrm, "%3d", (int) pqn[i]);
    pqn = pegylo->iqn;
    for (i = 0; i < nqn; ++i)
      fprintf(lutrm, "%3d", (int) pqn[i]);
    fprintf(lutrm, "%15.6f %10.6f\n", pline->freqx, pline->errfq);
  } /* end loop over experimental lines */
  pegylo = &head_egy; ndtag = 2;
  for (pegy = head_egy.next; pegy != &head_egy; pegy = pegy->next) {
    /* setup reverse references */
    pegy->prev = pegylo;
    pegylo = pegy;
    if (ndtag < pegy->nfreq) ndtag = pegy->nfreq;
  }
  head_egy.prev = pegyup = pegylo;
  while (pegyup != &head_egy) { /* sort energies according to iorg */ 
    k = 0; 
    for (pegy = pegyup; pegy != &head_egy; pegy = pegy->prev) {
      if (k < pegy->iorg) {
        k = pegy->iorg; pegylo = pegy;
        if (k == norg) break;
      }
    }
    if (k == 0) break;
    if (k < norg)
      norg = k;
    if (pegylo != pegyup) { /* move best egy to end */
      emove(pegylo, pegyup);
    } else {
      pegyup = pegyup->prev;
    }
  }
  tagv = (int *) mallocq((size_t) ndtag * sizeof(int));
  bgn_line = head_line.next; head_line.next = NULL;
  for (pegyup = head_egy.prev; pegyup != &head_egy; pegyup = pegyup->prev) {
    nfreq = pegyup->nfreq;
    if (nfreq <= 0) { /* dangling energy */
      fprintf(lutrm, "   new egy (%4d): ", pegyup->iorg);
      pqn = pegyup->iqn;
      for (i = 0; i < nqn; ++i)
        fprintf(lutrm, "%3d", (int) pqn[i]);
      fputc('\n', lutrm);
      continue;
    }
    /* find line with energy pegyup */
    last_line = NULL; k = 0;
    for (pline = bgn_line; pline != NULL; 
         last_line = pline, pline = tmpline) {
      tmpline = pline->next;
      /* find relevant frequencies */
      if (pline->plower == pegyup) {
        pline->plower = pline->pupper;
        pline->pupper = pegyup;
        pline->freqx = -pline->freqx;
      } else if (pline->pupper != pegyup) {
        continue;
      } 
      /* pline->pupper == pegyup */
      if (last_line == NULL) {
        bgn_line = tmpline;
      } else {
        last_line->next = tmpline;
      }
      pline->next = head_line.next; head_line.next = pline;
      pline = last_line; 
      if (++k >= nfreq) break;
    }
    nloop = 0; 
    assert (nfreq == k);
    if (nfreq > 1) {
      if (ndtag < nfreq) {
        free (tagv);
        ndtag = nfreq;
        tagv = (int *) mallocq((size_t) ndtag * sizeof(int));
      }
      for (k = 0; k < nfreq; ++k) 
        tagv[k] = 0;
      nloop = 0; ktag = critag = 1; xfrq = 0.; 
      xfrq0 = 0.; xerr0 = 0.;
      tmpline = head_line.next; 
      xsq = xsq0 = tmpline->errfq;
      pegy = tmpline->plower; 
      pline = crit_line = &head_line;
      for (k = 0; k < nfreq; ++k) {
        last_line = pline; pline = pline->next;
        if (pline == NULL) break;
        assert(pline->pupper == pegyup);
        if (tagv[k] != 0) continue;
        if (pline->plower != pegy) {
          pegy = pline->plower; ++ktag;
        }
        kloop = 0; xwt = 0.; sum = 0.; xerr = 0;
        tmpline = &head_line; 
        for (i = 0; i < nfreq; ++i) {
          tmpline = tmpline->next;
          if (tagv[i] != 0) continue;
          if (tmpline == NULL) break;
          if (tmpline->plower != pegy) continue;
          if (kloop == 0) {
            xfrq0 = tmpline->freqx; 
            xerr0 = tmpline->errfq; 
            xwt = 1.;
          } else {
            xerr1 = xerr0 / tmpline->errfq; 
            dtmp = xerr1 * (tmpline->freqx - xfrq0);
            xerr += dtmp * dtmp; sum += dtmp * xerr1;
            xwt += xerr1 * xerr1;
          }
          tagv[i] = ktag;
          ++kloop;
        }
        if (kloop == 0) {
          if (pline->errfq < xsq0) {
            xsq0 = pline->errfq;
            crit_line = last_line;
          }  
        } else if (kloop >= nloop){
          dtmp = 0.;
          if (kloop > 1) {
            dtmp = sum; sum = sum / xwt;
            xerr  = xerr0 * xerr0 * (xerr - sum * dtmp);
          }
          if (kloop > nloop || xerr < xsq) {
            nloop = kloop;;
            xfrq = dtmp + xfrq0; xsq = xerr; critag = ktag;
          }
        }
      }
      if (nloop == 0) {
        tmpline = crit_line;
      } else {
        tmpline = &head_line;
        dif = 0.; kloop = 0;
        pline = tmpline;
        for (k = 0; k < nfreq; ++k) {
          last_line = pline; pline = pline->next;
          if (tagv[k] != critag) continue;
          if (pline == NULL) break;
          dtmp = fabs(pline->freqx - xfrq);
          if (kloop == 0 || dif > dtmp) {
            ++kloop; tmpline = last_line; dif = dtmp;
            if (kloop == nloop) break;
          }
        }
      }  
      if (tmpline != &head_line) {
        last_line = tmpline; crit_line = tmpline->next;
        last_line->next = crit_line->next;
        crit_line->next = head_line.next; head_line.next = crit_line;
      }
    }
    crit_line = head_line.next;
    assert (crit_line->pupper == pegyup);
    pegy = crit_line->plower; 
    pegy->nfreq -= 1;
    while (nfreq > 1) { 
      --nfreq;
      pline = crit_line->next;
      crit_line->next = pline->next;
      /* do combination differences */
      xfrq = pline->freqx - crit_line->freqx;
      xerr = pline->errfq; xwt = crit_line->errfq;
      xerr = sqrt(xerr * xerr + xwt * xwt);
      pegylo = pline->plower;
      if (pegylo == pegy) {     /* loop */
        pegy->nfreq -= 1;
        fprintf(lutrm, "      loop (%6d, %6d): ", crit_line->seq,
                pline->seq);
        free(pline);
      } else {                  /* combination difference */
        pegy->nfreq += 1;
        ++nline;
        fprintf(lutrm, "%6d new (%6d, %6d): ", nline, crit_line->seq,
                pline->seq);
        pqn = pegy->iqn;
        for (i = 0; i < nqn; ++i)
          fprintf(lutrm, "%3d", (int) pqn[i]);
        pline->pupper = pegy;
        pline->freqx  = xfrq;
        pline->errfq  = xerr;
        pline->seq    = nline;
        pline->next = bgn_line; bgn_line = pline;
      }
      pqn = pegylo->iqn;
      for (i = 0; i < nqn; ++i)
        fprintf(lutrm, "%3d", (int) pqn[i]);
      fprintf(lutrm, "%15.6f %10.6f\n", xfrq, xerr);
    }
  } /* end backwards loop over energy */
  bgn_line = pline = last_line = tmpline = NULL;
  pegyup = NULL;
  while ((crit_line = head_line.next) != NULL) {
    /* back substitute energies from critical lines */
    pegyup = crit_line->pupper;
    pegylo = crit_line->plower;
    pegyup->egyx = pegylo->egyx + crit_line->freqx;
    xerr = crit_line->errfq; xwt = pegylo->errx;
    pegyup->errx = sqrt(xerr * xerr + xwt * xwt);
    head_line.next = crit_line->next;
    free(crit_line);
  }
  oldorg = -1; pegy = pegyorg = &head_egy; egyorg = 0.;
  do {
    /* print energies */
    norg = pegy->iorg;
    fprintf(lutrm, "%4d: ", norg);
    pqn = pegy->iqn;
    if (oldorg != norg) { /* find lowest energy in stack */
      pegyorg = pegy; oldorg = norg; egyorg = pegy->egyx;
      for (pegyup = pegy->next; pegyup->iorg == norg; pegyup = pegyup->next){
        if (pegyup == &head_egy) break;
        if (egyorg > pegyup->egyx) {
          egyorg = pegyup->egyx; pegyorg = pegyup;
        }
      }
    }
    egy = pegy->egyx - egyorg;
    for (i = 0; i < nqn; ++i)
      fprintf(lutrm, "%3d", (int) pqn[i]);
    fprintf(lutrm, "%15.6f %10.6f\n", egy, pegy->errx);
    if (pegy != pegyorg) {
      for (i = 0; i < nqn; ++i)
        fprintf(luegy, "%3d", (int) pqn[i]);
      pqn = pegyorg->iqn;
      for (i = 0; i < nqn; ++i)
        fprintf(luegy, "%3d", (int) pqn[i]);
      for (i = nqn + nqn; i < 12; ++i) 
        fputs("   ", luegy);
      if (pegy->errx < 1.e-6) pegy->errx = 1.e-6;
      fprintf(luegy, "%15.6f %10.6f\n", egy, -pegy->errx);
    }
  } while ((pegy = pegy->next) != &head_egy);
  pegyorg = NULL;
  rewind(lulin);
  pegyup = NULL; pegylo = NULL;
  nline = 0;
  while (getlin(lulin, nqn, idqn, iqnv, &xfrq, &xerr, &xwt,
                tmpstr, NDLINE) > 0) {
    k = -2;
    pegyup = lposn(iqnv, nqn, &k);
    pegylo = lposn(iqnv + nqn, nqn, &k);
    if (k != -2) continue;
    cdif = pegylo->egyx - pegyup->egyx;
    if (xerr > 0.)
      cdif *= clight;
    cdif += xfrq;
    xwt = fabs(cdif / xerr);
    fprintf(lutrm, " %5d new: ", ++nline);
    pqn = pegyup->iqn;
    for (i = 0; i < nqn; ++i)
      fprintf(lutrm, "%3d", (int) pqn[i]);
    pqn = pegylo->iqn;
    for (i = 0; i < nqn; ++i)
      fprintf(lutrm, "%3d", (int) pqn[i]);
    fprintf(lutrm, "%16.6f %10.6f %10.6f %8.2f\n", xfrq, xerr, cdif, xwt);
  }                             /* end loop over experimental lines */
  pegyup = NULL; pegylo = NULL;
  fclose(lulin);
  fclose(lutrm);
  free(tagv);
  return 0;
}                               /* main */

EGYS *lposn(ikey, nkey, iret)
short *ikey;
const int nkey;
int *iret;
{
  EGYS *now, *plast;
  short *iqnx;
  int k, icmp, iord;
  plast = now = &head_egy;
  if (ikey[0] < 0) return now;
  iord = (*iret);
  while ((now = now->next) != &head_egy) {
    iqnx = now->iqn;
    icmp = 0;
    for (k = 0; k < nkey; ++k) {
      icmp = iqnx[k] - ikey[k];
      if (icmp != 0) break;
    }
    if (icmp == 0) return now;
    if (iord >= 0 && icmp >= 0) break;
    plast = now;
  }
  /* allocate new value */
  now = (EGYS *) mallocq(sizeof(EGYS));
  now->next = plast->next; plast->next = now;
  now->prev = plast;
  now->egyx = now->errx = 0.;
  now->iorg = now->nfreq = 0;
  now->iqn[0] = ikey[0];
  iqnx = now->iqn;
  for (k = 1; k < nkey; ++k)
    iqnx[k] = ikey[k];
  while (k < NDQN)
    iqnx[k++] = 0;
  *iret += 1;
  return now;
} /* lposn */

void emove(pold, pend)
EGYS *pold, *pend;
{
  /* remove pold from present position and put it after pend */
  EGYS *plast, *pnext;
  if (pold == pend)  return;
  plast = pold->prev; 
  pnext = pold->next;
  pnext->prev = plast; 
  plast->next = pnext;
  /* make new link */
  pnext = pend->next; pnext->prev = pold;
  pold->prev = pend;
  pend->next = pold; 
  pold->next = pnext;
}
