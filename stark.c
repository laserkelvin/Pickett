/*   Copyright (C) 1989, California Institute of Technology            */
/*   All rights reserved.  U. S. Government Sponsorship under          */
/*   NASA Contract NAS7-918 is acknowledged.                           */

/*   Herbert M. Pickett, 19 December 1989                              */

/*     get command line file names (.str = calc input, .stk = output)  */
/*     open .str and read first line for NQN                           */
/*     DO UNTIL end of file on .str                                    */
/*        IF Fup or Flow > Fmax THEN                                   */
/*            printout stark coefficients for Fmax-2                   */
/*        ENDIF                                                        */
/*        find lower level or create new entry                         */
/*        find upper level or create new entry                         */
/*        IF degeneracy THEN                                           */
/*            print info                                               */
/*        ELSE                                                         */
/*            calc stark effect and add to upper and lower             */
/*        ENDIF                                                        */
/*        read in line from .str                                       */
/*     END DO                                                          */
/*   STATE = -1 (initialize first line)                                */
/*   STATE =  0 (read lines)                                           */
/*   STATE =  1 (EOF detected, process last line)                      */
/*   STATE =  2 (compete printout and exit)                            */
/***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calpgm.h"
#define NINDX 16
#define NDSQ  ((NINDX * NINDX + NINDX) >> 1)
typedef struct ab_str {
  double a[NDSQ], b[NDSQ], c[NDSQ];
  /*@owned@ */ /*@null@*/ struct ab_str *next;
  char iqn[13];
} AB_STR;
typedef /*@owned@*/ /*@null@*/ AB_STR *PAB_STR;

static PAB_STR head[3];
static double zero = 0;

static /*@dependent@*/ AB_STR *lposn(char *iqnu, int idx, int nn);
static int starks(double *a, /*@out@*/ double *b,/*@out@*/ double *c,
                  int iqfx, int iqdel, int ifhalf);
static int ival(char *str);

int main(argc, argv)
int argc;
char *argv[];
{
#define NFILE 2
  static const char *ext[NFILE] = { "str", "stk" };
  enum efile {estr, estk};
  FILE *lustk, *lustr;
  AB_STR *now;
  double *paupp, *pbupp, *pcupp, *palow, *pblow, *pclow;
  char *fname[NFILE+1];
  double strv[NINDX], str, frqc, strtot, big, cnv, eps, ta, tb, tc, frqx, strx;
  size_t nl;
  int iflg, iqfl, indx, iqfu, iqfx, i, j, k, n, iqfmx, kk;
  int nindx, ifhalf, nt, nqn, istate, iqnfmt, nnindx;
  int koff[NINDX + 1], indxv[NINDX];
  char iqnu[13], iqnl[13], iqnux[13], iqnlx[13], card[82];

  maxmem(&nl);
  strv[0] = 0.;
  indxv[0] = 0;
  big = 100.;
  eps = 1e-12;
  cnv = 0.503403;
  head[0] = head[1] = head[2] = NULL;
  kk = 0;
  for (i = 0; i <= NINDX; ++i) {
    kk += i;
    koff[i] = kk;
  }
  istate = -1;
  filget(argc, argv, NFILE, fname, ext);
  lustr = fopenq(fname[estr], "r");
  lustk = fopenq(fname[estk], "w");
  if (fgetstr(card, 82, lustr) < 60) {
    puts("bad STR file");
    exit(EXIT_FAILURE);
  }
  /* parse "(2E15.5,I5,1X,A,A,I5)" */
  memcpy(iqnux, card + 36, 12);
  memcpy(iqnlx, card + 48, 12);
  indx = atoi(card + 60);
  if (indx <= 0)
    indx = 1;
  card[35] = '\0';
  iqnfmt = atoi(card + 30);
  card[30] = '\0';
  strx = atof(card + 15);
  card[15] = '\0';
  frqx = atof(card);
  iqnu[12] = iqnl[12] = iqnux[12] = iqnlx[12] = '\0';
  nindx = nnindx = nt = iqfu = iqfl = iqfx = 0;
  strtot = 0.;
  ifhalf = iqnfmt / 10;
  nqn = iqnfmt - ifhalf * 10;
  if ((ifhalf / 10) % 5 == nqn)
    nqn = 1;
  nqn = nqn + nqn - 2;
  ifhalf &= 1;
  iqfmx = 1;
  do {                          /*  begin master loop */
    memcpy(iqnu, iqnux, 12);
    memcpy(iqnl, iqnlx, 12);
    frqc = frqx;
    str = strx;
    if (istate == 0) {
      if (fgetstr(card, 82, lustr) <= 60) {
        istate = 1;             /* end file */
      } else {
        /* parse "(2E15.5,6X,A,A,I5)" */
        memcpy(iqnux, card + 36, 12);
        memcpy(iqnlx, card + 48, 12);
        indx = atoi(card + 60);
        if (indx <= 0)
          indx = 1;
        card[30] = '\0';
        strx = atof(card + 15);
        card[15] = '\0';
        frqx = atof(card);
      }
    }
    if (istate <= 0) {
      if (strncmp(iqnux, iqnu, 12) == 0 && strncmp(iqnlx, iqnl, 12) == 0) {
        /*  same quanta, different dipole OR first line */
        if (indx > nindx) {
          if (indx > NINDX) {
            printf("NINDX, NDT = %d, %d\n", indx, NINDX);
            exit(EXIT_FAILURE);
          }
          nindx = indx;
          k = nnindx;
          nnindx = koff[indx];
          n = nnindx - k;
          for (j = 0; j <= 1; ++j) {    /*  zero out extra A and B values */
            now = head[j];
            while (now != NULL) {
              dcopy(n, &zero, 0, now->a + k, 1);
              dcopy(n, &zero, 0, now->b + k, 1);
              dcopy(n, &zero, 0, now->c + k, 1);
              now = now->next;
            }
          }
        }
        strv[nt] = strx;
        indxv[nt] = indx - 1;
        ++nt;
        strtot += strx;
        istate = 0;
        continue;
      }
    }
    if (istate < 2) {
      iqfu = ival(iqnu + nqn);
      iqfl = ival(iqnl + nqn);
      iqfx = (iqfu > iqfl) ? iqfu : iqfl;
      if (iqfx < iqfmx) {
        printf("F = %d is out of order\n", iqfx);
        fflush(stdout);
        istate = 2;
        iqfx = iqfmx + 1;
      }
    }
    while (iqfx > iqfmx) {      /*   purge data for old F */
      while ((now = head[0]) != NULL) { /*  more data for F */
        k = 0;
        palow = now->a;
        pblow = now->b;
        pclow = now->c;
        for (i = 1; i <= nindx; ++i) {
          for (j = 1; j <= i; ++j) {
            ta = palow[k];
            tb = pblow[k];
            tc = pclow[k];
            ++k;
            if (fabs(ta) + fabs(tb) > eps) {
              fprintf(lustk, "%15.6E %15.6E %15.6E %s %5d %5d\n",
                      ta, tb, tc, now->iqn, j, i);
            }
          }
        }
        head[0] = now->next;
        now->next = head[2];
        head[2] = now;
      }
      ++iqfmx;
      head[0] = head[1];
      head[1] = NULL;
      if (head[0] == NULL)
        iqfmx = iqfx;
      if (iqfx == iqfmx) {
        printf("Starting F= %d\n", iqfx);
        fflush(stdout);
      }
    }
    if (istate == 2)
      exit(0);
    /* do contributions for last transition */
    strtot = fabs(strtot) * cnv;
    if (frqc < big * strtot) {  /* degenerate */
      for (i = 0; i < nt; ++i) {
        ta = strv[i] * cnv;
        ta *= ta;
        if (ta > eps) {
          starks(&ta, &tb, &tc, iqfx, iqfu - iqfl, ifhalf);
          str = sqrt((ta > tb) ? ta : tb);
          fprintf(lustk, "%15.4f %15.6E %s %s %5d\n", frqc, str, iqnu,
                  iqnl, i + 1);
        }
      }
    } else { /* sum contributions to a and b and c*/
      iflg = iqfu - iqfmx;
      now = lposn(iqnu, iflg, nnindx);
      paupp = now->a;
      pbupp = now->b;
      pcupp = now->c;
      iflg = iqfl - iqfmx;
      now = lposn(iqnl, iflg, nnindx);
      palow = now->a;
      pblow = now->b;
      pclow = now->c;
      iflg = iqfu - iqfl;
      for (i = 0; i < nt; ++i) {
        kk = koff[indxv[i]];
        for (j = 0; j <= i; ++j) {
          ta = strv[i] * cnv;
          if (i == j) {
            ta *= ta;
          } else {
            ta *= 2. * strv[j] * cnv;
          }
          starks(&ta, &tb, &tc, iqfx, iflg, ifhalf);
          ta /= frqc;
          tb /= frqc;
          tc /= frqc;
          k = indxv[j] + kk;
          paupp[k] += ta;
          palow[k] -= ta;
          pbupp[k] += tb;
          pblow[k] -= tb;
          pcupp[k] += tc;
          pclow[k] -= tc;
        }
      }
    }
    strv[0] = strx;
    indxv[0] = indx - 1;
    nt = 1;
    strtot = strx;
    if (istate == 1) {
      istate = 2;
      iqfx = iqfmx + 1;
    }
  } while (istate >= 0);
  head[0] = NULL;
  head[1] = NULL;
  head[2] = NULL;
  return 0;
}                               /* main */

AB_STR *lposn(iqnu, idx, nn)
char *iqnu;
int idx, nn;
{
  AB_STR *now, *last;
  int i;

  i = idx;
  if (i != 0)
    i = 1;
  last = NULL;
  now = head[i];
  while (now != NULL) {
    if (strncmp(now->iqn, iqnu, 12) == 0)
      return now;
    last = now;
    now = now->next;
  }
  /*  take from top of free and put at top of new */
  if (head[2] == NULL) {
    now = (AB_STR *) mallocq(sizeof(AB_STR));
  } else {
    now = head[2];
    head[2] = now->next;
  }
  if (last == NULL) {
    head[i] = now;
  } else {
    last->next = now;
  }
  now->next = NULL;
  dcopy(nn, &zero, 0, now->a, 1);
  dcopy(nn, &zero, 0, now->b, 1);
  dcopy(nn, &zero, 0, now->c, 1);
  memcpy(now->iqn, iqnu, 13);
  return now;
}                               /* lposn */

int starks(a, b, c, iqfx, iqdel, ifhalf)
double *a, *b, *c;
int iqfx, iqdel, ifhalf;
{
  int idgn;
  double f, fsq, aa, den;

  idgn = iqfx + iqfx - ifhalf; // idgn = 2F
  f = 0.5 * idgn; // f = F
  fsq = f * f;    // fsq = (F * F)
  aa = (*a);
  if (iqdel == 0) {
    den = (f + fsq) * (idgn + 1);
    aa = aa / den;
    *c = (fsq + f) * aa;
    *b = aa;
    *a = 0.;
  } else {
    den = f * (4. * fsq - 1);
    aa = aa / den;
    *c = (fsq - f) * aa;
    *b = -aa;
    *a = fsq * aa;
  }
  return 0;
}                               /* starks */

int ival(str)
char *str;
{
  int iret, i;

  i = (int) str[0] & 0xff;
  iret = i - (int) '0';
  if (iret > 9)
    iret = i - ((int) 'A' - 10);
  if (iret < 0)
    iret = 0;
  i = ((int) str[1] & 0xff) - (int) '0';
  if (i >= 0)
    iret = 10 * iret + i;
  return iret;
}                               /* ival */
