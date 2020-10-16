#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "calpgm.h"


/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
/*   Revised version in c, 22 March 1999 */
/*   fixed bug in that freed structure before data needed */
/*   funny comments with @ are modifiers for splint syntax checker */
/*   8/12/04 fixed format for big frequencies */
/**************************************************************************/
#define NDQN  20
#define NHASH 100
#define KFRQ  13
#define PERR  KFRQ
#define KERR  8
#define PSTR  PERR+KERR
#define KSTR  8
#define PEGY  PSTR+KSTR+2
#define KEGY  10
#define PDGN  PEGY+KEGY
#define KDGN  3
#define PTAG  PDGN+KDGN
#define KTAG  7
#define PQNF  PTAG+KTAG
#define KQNF  4
#define PQNU  PQNF+KQNF

typedef struct expt {
  /*@owned@*/ /*@null@*/ struct expt *next;
  double frqx;
  float errx;
  short iqn[NDQN], nbad;
} EXPT;
typedef /*@owned@*/ /*@null@*/ EXPT *PEXPT;
typedef /*@dependent@*/ /*@null@*/ EXPT *PDEXPT;

void prline(FILE *lu, int nqn, short *pqn, double frq, double err, double st,
            char *label);
int lposn(const short *ikey, const int nkey, const int *k12,
          PEXPT *head, /*@out@*/ PDEXPT *ppexpt);
void dolink(PEXPT *head, PDEXPT pexpt,
            const double xfrq, const double xerr, const short *iqn);
void unlink(PEXPT *head, /*@notnull@*/ EXPT *pexpt, 
            /*@out@*/ int *jmin,/*@out@*/ int *jmax);
BOOL mrglin(char *line2, char *line);

int main(argc, argv)
int argc;
char *argv[];
{
#define NFILE 6
#define NDLINE 130
  static PEXPT jhash[NHASH + 1];
  static const char *cext[NFILE] = {"cat","lin","mrg","new","log","opt"};
  enum efile {ecat, elin, emrg, enew, elog, eopt};
  static const double clight = 29979.2458;
  static const double big = 9999.9999;
  static int fmtlen[] = { KFRQ, KERR, KSTR };
  static const int blank = (int)' ';
  static int fmtqnf = KQNF;
  static short iqn1[NDQN], iqn2[NDQN], iqmap[NDQN], iqmaplin[NDQN];
  static char *cfil[NFILE+1];
  static char line1[NDLINE], line2[NDLINE], tmpstr[NDLINE];
  static char qnstr[NDQN+NDQN+1], tagstr[8];
  short *pqn, tqn;
  double cfrq, cfrq0, errc, xerr, st, st0, xfrq, frqdif, xwt, ratio, dtmp[3];
  int j, k, nbad, nout, nline, iqfmt, jmax, nqn2, nline2, catqn, pqnlo, jx;
  int minb, k1, k2, k12[2], nqn, nmatch, nmerge, icmp, nqnlin, nqnlin2;
  int jmax0, jmin, jmin0, nmap, kcmp, nxhash;
  BOOL match, noeof;
  FILE *lucat, *lunew, *lulog, *lumrg, *lulin, *luopt;
  EXPT *pexpt, *pexptp, *plast, *pnext;
  nmerge = nmatch = nout = nqn = nqn2 = 0;
  filget(argc, argv, NFILE, cfil, cext);
  rqexit(-1);
  for (k = 0; k <= NHASH; ++k)
    jhash[k] = NULL;
  lulin = fopenq(cfil[elin], "r");
  lulog = fopenq(cfil[elog], "w");
  lunew = fopenq(cfil[enew], "w");
  lucat = fopenq(cfil[ecat], "r");
  /* read qn format */
  nline2 = fgetstr(line2, NDLINE, lucat);
  iqfmt = 1;
  if (nline2 > PQNU) {
    memcpy(tagstr, line2 + PTAG, KTAG);
    tagstr[KTAG] = '\0';
    for (k = 1; k < KTAG; k++) {
      if (tagstr[k] == '-')
        break;
      if (tagstr[k] != ' ') {
        tagstr[k - 1] = '-';
        break;
      }
    }
    pcard(line2 + PQNF, dtmp, 1, &fmtqnf);
    iqfmt = (int) dtmp[0];
  }
  nqn = deflin(iqfmt, iqn1);
  nqn2 = nqn + nqn;
  catqn = 12;
  if (nqn > 6) catqn = nqn2;
  pqnlo = PQNU + catqn;
  minb  = pqnlo + nqn2;
  if (nline2 < minb) {
    puts("Bad CAT File");
    exit(EXIT_FAILURE);
  }
  /* read options */
  ratio = 0.; nqnlin = nqn; 
  luopt = fopen(cfil[eopt],"r");
  if (luopt == NULL) {
    printf("enter MAX (obs-calc)/(exp. error), no. quanta for lin (%2d) "
           , nqn);
    fflush(stdout);
    luopt = stdin;
  }
  if (fgetstr(tmpstr, NDLINE, luopt) > 0) {
    dtmp[0] = 100.; dtmp[1] = (double)nqn;
    pcard(tmpstr, dtmp, 2, NULL);
    ratio = dtmp[0]; nqnlin = (int) dtmp[1];
    if (nqnlin > 10) nqnlin = 10;
  }
  nmap = 0;
  if (nqnlin != nqn) {
    if (nqnlin < 0) {
      nqnlin = -nqnlin;
      if (nqnlin > nqn)
        nqnlin = nqn;
      --nqnlin;
      iqmap[0] = (short)nqnlin; iqmaplin[0] = (short)nqnlin;
      nmap = 1; nqnlin = nqn;
    } else {
      nmap = (nqnlin < nqn)? nqnlin: nqn;
      nmap = nmap - 2;
      for (k = 0; k < nmap; ++k) {
        iqmap[k] = (short)k; iqmaplin[k] = (short)k;
      }
      iqmap[nmap] = (short)(nqn - 1); 
      iqmaplin[nmap] = (short)(nqnlin - 1);
      ++nmap;
    }
    for (k = 0; k < nmap; ++k) {
      iqmap[k + nmap] = (short)(iqmap[k] + nqn); 
      iqmaplin[k + nmap] = (short)(iqmaplin[k] + nqnlin);
    }
    nmap = nmap << 1;
  }     
  nqnlin2 = nqnlin << 1;
  if (luopt != stdin) fclose(luopt);
  if (nqn <= 6) {
    memset(line1 + pqnlo, blank, 13);
  }
  k1 = 1; k2 = k1 + nqn;
  k12[0] = -1; k12[1] = k2;
  /*  read and sort experimental lines */
  nline = 0; jmax = 0; jmin = NHASH;
  for (;;) {
    if (getlin(lulin, nqnlin, iqn1, iqn2, &xfrq, &xerr, &xwt,
               tmpstr, NDLINE) < 0) {
      fclose(lulin);
      break;
    }
    if (iqn2[nqnlin] < 0)
      continue;
    errc = xerr;
    if (xerr < 0.) {
      xfrq *= clight;
      errc *= clight;
      xerr = -errc;
    }
    if (xerr > 999.999)
      continue;
    if (xfrq < 0.) { 
      /* swap quanta so frequency is positive */
      for (k = 0; k < nqnlin; ++k) {
        tqn = iqn2[k];
        iqn2[k] = iqn2[k + nqnlin];
        iqn2[k + nqnlin] = tqn;
      }
      xfrq = -xfrq;
    }
    /* set up HASH table */
    j = (int)iqn2[0]; 
    if (j > jmax)
      jmax = j;
    if (j < jmin)
      jmin = j;
    icmp = lposn(iqn2, nqnlin2, k12, jhash, &pexpt);
    if (icmp == 0) {
      /* lines have identical quanta */
      xwt = pexpt->errx;
      if (xwt < 0.)
        xwt = -xwt;
      if (xerr < xwt) {
        pexpt->frqx = xfrq;
        pexpt->errx = (float) xerr;
      }
      continue;
    } else {
      dolink(jhash, pexpt, xfrq, xerr, iqn2);
    }
    ++nline;
  }
  if (nline == 0) ratio = 0.;
  /* open .mrg output file */
  lumrg = fopenq(cfil[emrg], "w");
  printf(" %d  experimental lines read\n", nline);
  fflush(stdout);
  memcpy(qnstr, line2 + PQNU, (size_t) nqn2);
  memcpy(qnstr + nqn2, line2 + pqnlo, (size_t) nqn2);
  readqn(qnstr, iqn2, nqn2);
  noeof = TRUE;
  while (noeof) {
    memcpy(line1, line2, (size_t) minb);
    memcpy(iqn1, iqn2, sizeof(iqn2));
    if (rqexit(0) != 0)
      noeof = FALSE;
    nline2 = fgetstr(line2, NDLINE, lucat);
    if (nline2 >= minb) {
      memcpy(qnstr, line2 + PQNU, (size_t) nqn2);
      memcpy(qnstr + nqn2, line2 + pqnlo, (size_t) nqn2);
      readqn(qnstr, iqn2, nqn2);
    } else {                    /* set NOEOF = false for end of file */
      noeof = FALSE;
    }
    if (strncmp(line1 + PTAG, tagstr, KTAG) == 0) {
      /* write through experimental lines */
      ++nout;
      fputs(line1, lumrg); fputc('\n',lumrg);
      continue;
    }
    match = noeof;
    /* ...test for equality of frequency and quanta except K */
    if (match)
      match = (strncmp(line1, line2, KFRQ) == 0);
    if (match) {
      icmp = 0; kcmp = k1;
      for (k = 0; k < nqn2; ++k) {
        icmp = iqn1[k] - iqn2[k];
        if (k == kcmp) {
          if (icmp != 0) 
            icmp = iqn1[k] + iqn2[k];
          kcmp = k2;
        }
        if (icmp != 0)
          break;
      }
      match = (icmp == 0);
    }
    if (match)
      match = mrglin(line1, line2);
    if (match) {                /*  K's are equivalent */
      if (iqn2[k1] < 0) {
        line2[PQNU+2] = line1[PQNU+2];
        iqn2[k1] = iqn1[k1];
      }
      if (iqn2[k2] < 0) {
        line2[pqnlo+2] = line1[pqnlo+2];
        iqn2[k2] = iqn1[k2];
      }
      k12[0] = k1;
      ++nmerge; 
    } else {                    /* ready to ouput line */
      j = (int)iqn1[0];         /* test for quantum match */
      if (j >= jmin && j <= jmax && nqn == nqnlin && 
          lposn(iqn1, nqn2, k12, jhash, &pexpt) == 0) {
        pcard(line1, dtmp, 3, fmtlen);
        cfrq = dtmp[0];
        errc = dtmp[1];
        st   = dtmp[2];
        xfrq = pexpt->frqx;
        xerr = pexpt->errx;
        if (xfrq < 99999999.) {
          sprintf(tmpstr,"%13.4f%8.4f", xfrq, xerr);
        } else {
          sprintf(tmpstr,"%13.3f%8.3f", xfrq, xerr);
        }
        memcpy(tmpstr + PSTR, line1 + PSTR, PTAG - PSTR);
        memcpy(tmpstr + PTAG, tagstr, KTAG);
        strcpy(tmpstr + PQNF, line1 + PQNF);
        frqdif = xfrq - cfrq;
        if (frqdif > big) {
          frqdif = big;
        } else if (frqdif < -big) {
          frqdif = -big;
        }
        fprintf(lulog, "%.25s %14.4f%10.4f%10.4f%10.4f\n", line1 + PQNU,
                xfrq, frqdif, xerr, errc);
        if (fabs(frqdif) < xerr * ratio) {
          st = pow(10., st);
          prline(lunew, nqn2, iqn1, xfrq, xerr, st, "q");
          strcpy(line1, tmpstr);
          ++nmatch;
          unlink(jhash, pexpt, &jmin, &jmax);
          if (k12[0] >= 0 && lposn(iqn1, nqn2, k12, jhash, &pexpt) == 0) {
            ++nmatch;
            unlink(jhash, pexpt, &jmin, &jmax);
          }
        } else {
          fputs(" Above Experiment not used\n", lulog);
        }
      }
      ++nout; k12[0] = -1;
      fputs(line1, lumrg); fputc('\n',lumrg);
    }
  }                             /* end loop over calculated lines */
  printf(       "%d  lines written %d experimental lines \n", nout, nmatch);
  fprintf(lulog,"%d  lines written %d experimental lines \n", nout, nmatch);
  if (nmerge > 0) {
    printf(       " %d lines merged\n", nmerge);
    fprintf(lulog," %d lines merged\n", nmerge);
  }
  fclose(lucat);
  fclose(lumrg);
  nbad = nline - nmatch;
  nxhash = jmax;
  if (nxhash > NHASH)
    nxhash = NHASH;
  if (nbad > 0) {     /* find bad lines */
    printf(        " %d  lines to be matched by frequency\n", nbad);
    fprintf(lulog, " %d  lines to be matched by frequency\n", nbad);
  }
  nmatch = 0;
  lumrg = fopenq(cfil[emrg], "r");
  cfrq0 = -999999.; st0 = -500.; 
  if (ratio < 0.001)
    nbad = 0;
  while (nbad > 0) {
    if (jmin > jmax) break;
    if (fgetstr(line1, NDLINE, lumrg) < minb) break;
    if (strncmp(line1 + PTAG, tagstr, KTAG) == 0) 
      continue;  /* line already assigned */
    memcpy(qnstr, line1 + PQNU, (size_t) nqn2);
    memcpy(qnstr + nqn2, line1 + pqnlo, (size_t) nqn2);
    readqn(qnstr, iqn1, nqn2);
    jmin0 = jmin; jmax0 = nxhash;
    if (nmap >  0) {
      if (iqmap[0] == 0) {
        jmin0 = iqn1[0];
        if (jmin0 < jmin) continue;
        if (jmin0 > jmax) continue;
        if (jmin0 > NHASH)
          jmin0 = NHASH;
        jmax0 = jmin0;
      }
      for (k = 0; k < nmap; ++k) 
        iqn2[k] = iqn1[iqmap[k]];
    } 
    pcard(line1, dtmp, 3, fmtlen);
    //    cfrq0 = cfrq; st0 = st; 
    cfrq = dtmp[0];
    errc = dtmp[1];
    st = dtmp[2];
    //    if (fabs(cfrq - cfrq0) < 0.00015 && fabs(st - st0) < 0.00015) continue;
    st = pow(10., st);
    for (jx = jmin0; jx <= jmax0; ++jx) {
      for (pexpt = jhash[jx]; pexpt != NULL; pexpt = pexpt->next) {
        xfrq = pexpt->frqx;
        frqdif = cfrq - xfrq;
        xerr = pexpt->errx;
        if (fabs(frqdif) < ratio * xerr) {
          pqn = pexpt->iqn;
          if (nmap > 0) {
            icmp = 0;
            for (k = 0; k < nmap; ++k) {
              icmp = (int) (iqn2[k] - pqn[iqmaplin[k]]);
              if (icmp != 0) break;
            }
            if (icmp != 0) continue;
          }
          if (nqn2 != nqnlin2) {
            prline(lunew, nqnlin2,  pqn, xfrq, xerr, st, "f");
          } else {
            prline(lunew, nqnlin2, iqn1, xfrq, xerr, st, "f");
            pexptp = pexpt;
            for(;;) {
              plast = pexptp; pexptp = pexptp->next;
              if (pexptp == NULL) break;
              if (fabs(xfrq - pexptp->frqx) < 0.00015) {
                pnext = pexptp->next; 
                free(plast->next);
                plast->next = pnext; pnext = NULL; 
                pexptp = plast;
              }
            }
          }
          pexpt->nbad += (short) 1;
          ++nmatch;
        }
      } 
    }
  }
  if (nbad > 0) {
    printf(       "%d  lines found\n", nmatch);
    fprintf(lulog,"%d  lines found\n", nmatch);
    st = 0.;
    for (jx = jmin; jx <= nxhash; ++jx) {
      pexpt = jhash[jx];
      while (pexpt != NULL) {
        pexptp = pexpt->next;
        if (pexpt->nbad == (short) 0) {
          pqn = pexpt->iqn;
          xfrq = pexpt->frqx;
          xerr = pexpt->errx;
          prline(lunew, nqnlin2, pqn, xfrq, xerr, st, "u");
          pqn = NULL;
        }
        unlink(jhash, pexpt, &jmin, &jmax);
        pexpt = pexptp;
      }
    }
    pexptp = NULL;
  }
  fclose(lumrg);
  fclose(lunew);
  fclose(lulog);
  return 0;
}                               /* main */

void prline(FILE *lu, int nqn, short *pqn, double frq, double err, double st,
            char *label)
{
  int k;
  for (k = 0; k < nqn; ++k) {
    fprintf(lu, "%3d", (int) pqn[k]);
  }
  for (k = nqn; k < 12; ++k) {
    fputs("   ", lu);
  }
  fprintf(lu, "%15.4f %10.4f %10.3E %s\n", frq, err, st, label);
}

int lposn(ikey, nkey, k12, head, ppexpt)
const short *ikey;
const int nkey, *k12;
PEXPT *head;
PDEXPT *ppexpt;
{
  EXPT *now, *last;
  short *iqnx;
  int k, icmp, kcmp, jx;
  last = NULL;
  jx = ikey[0];
  if (jx > NHASH)
    jx = NHASH;
  now = head[jx];
  icmp = -1; 
  while (now != NULL) {
    iqnx = now->iqn; 
    icmp = iqnx[0] - ikey[0];
    if (icmp == 0) {
      kcmp = k12[0];
      for (k = 1; k < nkey; ++k) {
        icmp = iqnx[k] - ikey[k];
        if (k == kcmp) {
          if (icmp != 0)
            icmp = -iqnx[k] - ikey[k];
          kcmp = k12[1];
        }
        if (icmp != 0) break;
      }
    }
    if (icmp >= 0)
      break;
    last = now;
    now = now->next;
  }
  if (icmp != 0){ 
    *ppexpt = last;
  } else {
    *ppexpt = now;
  }
  return icmp;
}                               /* lposn */

void dolink(head, pexpt, xfrq, xerr, iqn)
PEXPT *head;
PDEXPT pexpt;
const double xfrq, xerr;
const short *iqn;
{
  /* insert new link after pexpt */
  EXPT *now;
  int jx;
  now = (EXPT *) mallocq(sizeof(EXPT));
  now->frqx = xfrq;
  now->errx = (float) xerr;
  memmove(now->iqn, iqn, (size_t) NDQN * sizeof(short));
  now->nbad = 0;
  if (pexpt == NULL) {
    jx = (int)iqn[0];
    if (jx > NHASH)
      jx = NHASH;
    now->next = head[jx];
    head[jx] = now;
  } else {
    now->next = pexpt->next;
    pexpt->next = now;
  }
}

void unlink(head, pexpt, jmin, jmax)
PEXPT *head;
EXPT *pexpt;
int *jmin, *jmax;
{
  EXPT *now, *last;
  int jx, jj;
  jx = (int)(pexpt->iqn[0]);
  if (jx > NHASH)
    jx = NHASH;
  now = head[jx];
  if (now == NULL) return;
  if (now == pexpt) { /* pexpt == head */
    head[jx] = now->next;
    if (head[jx] == NULL) {
      if (*jmax == jx) {
        jj = (*jmax) - 1;
        while(*jmin < jj && head[jj] == NULL) --jj;
        *jmax = jj;
      } else if (*jmin == jx) {
        jj = (*jmin) + 1;
        while(*jmax > jj && head[jj] == NULL) ++jj;
        *jmin = jj;
      }
    }
    free(pexpt);
    return;
  }
  do { /* find pexpt */
    last = now; now = now->next;
    if (now == pexpt) {
      now = pexpt->next;
      if (now == NULL && jx == NHASH) {
        *jmax = (int)(last->iqn[0]);
      }
      last->next = now;
      free(pexpt);
      return;
    } 
  } while (now != NULL);
}

BOOL mrglin(line2, line)
char *line2, *line;
{
  static double cnv = 2.3025851; /* 1 / log10(e) */
  static int fmtstr = KSTR;
  static char tmpstr[20];
  short ig1, ig2;
  int igup;
  double dif, str1, str2;

  pcard(line + PSTR, &str1, 1, &fmtstr);
  readqn(line + PDGN, &ig1, 0);
  igup = ig1;
  if (strncmp(line + PSTR, line2 + PSTR, KSTR) == 0) {
    str1 += 0.30103;
    igup += ig1;
  } else {
    pcard(line2 + PSTR, &str2, 1, &fmtstr);
    readqn(line2 + PDGN, &ig2, 0);
    if (str2 < str1) {
      dif = str2 - str1;
    } else {
      dif = str1 - str2;
      str1 = str2;
    }
    if (dif > -5.) {
      dif = exp(dif * cnv) + 1.;
      str1 += log10(dif);
    }
    igup += ig2;
  }
  if (igup >= 3600)
    return FALSE;
  sprintf(tmpstr,"%8.4f", str1);
  memcpy(line + PSTR, tmpstr, KSTR);
  gupfmt(igup, line + PDGN);
  return TRUE;
}                               /* mrglin */



