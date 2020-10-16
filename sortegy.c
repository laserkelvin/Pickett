#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "calpgm.h"


/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
/**************************************************************************/
#define NDQN  10

typedef struct egy {
  /*@owned@*/ /*@null@*/ struct egy *next;
  /*@dependent@*/ /*@null@*/ struct egy *cross;
  double egy;
  int block;
  short iqn[NDQN];
} EGY;
typedef /*@owned@*/ /*@null@*/ EGY *PEGY;
typedef /*@dependent@*/ /*@null@*/ EGY *PDEGY;


int qncmp(const int nqn, const short *qn1, const short *qn2);
PDEGY qnfind(/*@partial@*/ /*@null@*/EGY *head, const short *iqn,
             const int nqn);
void qncopy(const int nqn, const short *qn1, /*@out@*/ short *qn2);
void egymove(/*@partial@*/ PDEGY plast,/*@partial@*/ PDEGY porg);

int main(argc, argv)
int argc;
char *argv[];
{
#define NFILE 6
#define NDLINE 128
  static const char *cext[NFILE] = {"egy","lgy","inp","lin","opt","log"};
  enum efile {eegy, elgy, einp, elin, eopt, elog};
  static const double clight = 29979.2458;
  static const double tiny = 1.e-6;
  static int iblkv[2], nqnv[2];
  static short iqnq[20], iqn2[20];
  static char *cfil[NFILE+1];
  static char tmpstr[NDLINE];
  static EGY jhead1, jhead2;
  FILE *luegy, *lulgy, *luinp, *lulin, *luopt, *lulog;
  EGY  *pnow, *plast2, *pbest, *pnext, *plast;
  PDEGY pegyv[2], pegy1, pegy2;
  short *pqn, *pqn1;
  double egyc, xfrq, xerr, xwt, dif, dif2, dtmp[6+NDQN+NDQN];
  int j, k, kk, nline, iqfmt, nqn2, nqn, nblk, iblk, iopt, knt;
  nqn = 4;
  filget(argc, argv, NFILE, cfil, cext);
  rqexit(-1);
  if (jhead1.next != NULL) free(jhead1.next);
  jhead1.next = NULL; jhead1.cross = NULL; jhead1.egy = 0.; jhead1.block = -1;
  if (jhead2.next != NULL) free(jhead2.next);
  jhead2.next = NULL; jhead2.cross = NULL; jhead2.egy = 0.; jhead2.block = -1; 
  jhead1.iqn[0] = -1; jhead2.iqn[0] = -1;
  /* read options */
  luopt = fopen(cfil[eopt],"r");
  if (luopt == NULL) {
    printf("number of quanta per state, option ");
    fflush(stdout);
    luopt = stdin;
  }
  nqn = 4; iopt = 0;
  if (fgetstr(tmpstr, NDLINE, luopt) > 0) {
    dtmp[0] = 4.; dtmp[1] = 0.;
    pcard(tmpstr, dtmp, 2, NULL);
    nqn = (int) dtmp[0];
    iopt = (int) dtmp[1];
    if (nqn > NDQN) nqn = NDQN;
  }
  nqnv[0] = nqn; nqnv[1] = 0;
  if (luopt != stdin) fclose(luopt);
  iqfmt = nqn + 400;
  nqn = deflin(iqfmt, iqnq);
  nqn2 = nqn + nqn;
  nline = 0;
  luegy = fopenq(cfil[eegy], "r");
  /* read egy file */
  while (fgetstr(tmpstr, NDLINE, luegy) > 0) {
    dtmp[0] = -1.;
    if (pcard(tmpstr, dtmp, nqn + 6, NULL) < 1) break;
    iblk = (int)dtmp[0];
    egyc = dtmp[2];
    for (k = 0; k < nqn; ++k){
      iqn2[k] = (short) dtmp[k + 6];
    }
    if (*iqn2 < 0) continue;
    /* look for redundant energy */
    pegy1 = qnfind(jhead1.next, iqn2, nqn);
    if (pegy1 != NULL) continue;
    if (iopt == 1) {
      iblk = iqn2[0] + iqn2[0];
      if (((int) iqn2[3] % 3) != 0) ++iblk;
    }
    plast = &jhead1;
    for (pnext = jhead1.next; pnext != NULL; pnext = pnext->next) {
      j = pnext->block;
      if (j > iblk) break;
      if (j == iblk && pnext->egy > egyc) break;
      plast = pnext;
    }
    pnow =(EGY *) mallocq(sizeof(EGY));
    pnow->next = pnow;
    egymove(pnow, plast);
    pnow->cross = NULL; pnow->egy = egyc; 
    pnow->block = iblk; 
    qncopy(nqn, iqn2, pnow->iqn);
    pnow = NULL;
    ++nline;
  }
  fclose(luegy);
  nblk = 0;
  if (jhead1.next != NULL) {
    iblk = jhead1.block;
    for (pegy1 = jhead1.next; pegy1 != NULL; pegy1 = pegy1->next) {
      kk = iblk; iblk = pegy1->block;
      if (iblk != kk) ++nblk;
      pegy1->block = nblk - 1;
    }
  }
  printf(" %4d calculated energies read\n", nline); 
  if (nline <= 0) exit (0);
  /*  read and sort experimental lines from lgy file */
  lulgy = fopenq(cfil[elgy], "r");
  nline = 0; 
  for (;;) {
    if (getlin(lulgy, nqn, iqnq, iqn2, &xfrq, &xerr, &xwt,
               tmpstr, NDLINE) < 0) {
      break;
    }
    if (xerr >= 0.)
      xfrq /= clight;
    knt = 0;
    for (kk = 1; kk >= 0; --kk) { 
      iblkv[kk] = -1; pegyv[kk] = NULL;
      pqn = iqn2 + nqnv[kk];
      if (*pqn < 0) continue;
      pegy1 = qnfind(jhead1.next, pqn, nqn);
      if (pegy1 == NULL) break;
      iblkv[kk] = pegy1->block;
      pegy2 = qnfind(jhead2.next, pqn, nqn); pegyv[kk] = pegy2;
      if (pegy2 != NULL) knt += kk + 1;
    }
    if (kk >= 0 || knt == 3 || iblkv[1] < 0) 
      continue;
    if (iblkv[0] < 0 && knt == 2) { // new absolute energy
      pegy2 = pegyv[1];
      if (pegy2 == NULL) continue;
      xfrq -= pegy2->egy;
      while(pegy2 != NULL) {
        pegy2->egy += xfrq;
        plast2 = pegy2->cross; pegy2->cross = NULL; pegy2 = plast2;
      }
      continue;
    }
    for (kk = 1; kk >= 0; --kk) { 
      if (pegyv[kk] == NULL) {
        iblk = iblkv[kk];
        pqn = iqn2 + nqnv[kk];
        plast = &jhead2;
        for (pnext = jhead2.next; pnext != NULL; pnext = pnext->next) {
          if (pnext->block >= iblk) break;
          plast = pnext;
        }
        ++nline;
        pnow = (EGY *) mallocq(sizeof(EGY));
        pnow->next = pnow;
        egymove(pnow, plast);
        pnow->cross = NULL; pnow->egy = xfrq;
        pnow->block = iblkv[kk];
        qncopy(nqn, pqn, pnow->iqn);
        pegyv[kk] = pnow; 
        pnow = NULL; plast = NULL; pnext = NULL;
        if (iblkv[0] < 0) break;
      }
    }
    if (iblkv[0] < 0) continue;
    plast = pegyv[0]; 
    if (plast == NULL) continue;
    pegy2 = pegyv[1];
    if (pegy2 == NULL) continue;
    if (knt == 0) { // isolated pair
      plast->cross = pegy2; pegy2->cross = plast;
      plast->egy = 0.;
    } else if (knt == 1) { // common lower state (plast)
      pnext = plast->cross; pegy2->cross = pnext; 
      if (pnext != NULL) 
        plast->cross = pegy2;
      pegy2->egy = plast->egy + xfrq; 
    } else { // common upper state (pegy2)
      pnext = pegy2->cross; plast->cross = pnext; 
      if (pnext != NULL)
        pegy2->cross = plast;
      plast->egy = pegy2->egy - xfrq;
    } 
  } /* end loop over input */
  fclose(lulgy);
  printf(" %4d experimental energies read\n", nline);
  if (nline == 0) exit(0);
  fflush(stdout);
  /* fill in origins with calc */
  for (pegy2 = jhead2.next; pegy2 != NULL; pegy2 = pegy2->next) { 
    if (pegy2->cross == NULL) continue;
    plast = pegy2; pbest = plast; xfrq = pbest->egy;
    do {
      plast = plast->cross;
      if (plast == NULL) break;
      if (plast->egy < xfrq) {
        xfrq = plast->egy; pbest = plast;
      }
    } while (plast != pegy2);
    if (pbest == NULL) continue;
    pegy1 = qnfind(jhead1.next, pbest->iqn, nqn);
    if (pegy1 == NULL) pegy1 = &jhead1;
    xfrq = pegy1->egy - xfrq;
    plast = pegy2;
    while (plast != NULL) {
      plast->egy += xfrq;
      plast2 = plast->cross; plast->cross = NULL; plast = plast2;
    } 
  }
  /* sort experimental frequencies */
  plast = &jhead2;
  for (pegy2 = jhead2.next; pegy2 != NULL; pegy2 = pegy2->next) { 
    iblk = pegy2->block;
    plast2 = pegy2; pbest = plast; dif = pegy2->egy; 
    for (pnext = pegy2->next; pnext != NULL && pnext->block == iblk; 
         pnext = pnext->next) {
      if (dif > pnext->egy) {
        dif = pnext->egy; pbest = plast2;
      }
      plast2 = pnext;
    } 
    if (plast != pbest) {
      egymove(pbest, plast);
      pegy2 = plast;
    }
    plast = pegy2;
  }
  /* find gaps */
  pegy1 = &jhead1;
  for (pegy2 = jhead2.next; pegy2 != NULL; pegy2 = pegy2->next) {
    iblk = pegy2->block;
    while (pegy1 != NULL && pegy1->block < iblk) 
      pegy1 = pegy1->next;
    while (pegy1 != NULL && pegy1->cross != NULL) 
      pegy1 = pegy1->next;
    if (pegy1 == NULL) break;
    if (pegy1->block != iblk) continue;
    pbest = pegy1;
    egyc = pegy2->egy;
    dif = fabs(egyc - pegy1->egy); 
    for (pnext = pegy1->next; pnext != NULL && pnext->block == iblk;
	 pnext = pnext->next) {
      if (pnext->cross != NULL) continue;
      xwt = fabs(egyc - pnext->egy);
      if (xwt < dif) {
        dif = xwt; pbest = pnext;
      }
    }
    pbest->cross = pegy2;
  }
  /* swap near degeneracies */
  pegy1 = &jhead1;
  while (pegy1 != NULL) {
    plast = pegy1; pegy1 = pegy1->next;
    if (pegy1 == NULL) break;
    if (pegy1->block != plast->block) continue;
    pegy2 = pegy1->cross; plast2 = plast->cross;
    if (pegy2 == NULL || plast2 == NULL) continue;
    dif  = pegy1->egy - plast->egy;
    if (dif > 0.1) continue;
    dif2 = pegy2->egy - plast2->egy;
    if (fabs(dif2) > tiny && fabs(dif) > tiny) {
      if (dif2 > 0.) continue;
    } else {
      if (qncmp(nqn, pegy1->iqn,  pegy2->iqn) == 0 &&
          qncmp(nqn, plast->iqn, plast2->iqn) == 0) continue;
    }
    pegy1->cross = plast2; plast->cross = pegy2;
  }
  pegy1 = NULL;
  /* print log */
  lulog = fopenq(cfil[elog], "w");
  for (pegy1 = jhead1.next; pegy1 != NULL; pegy1 = pegy1->next) {
    fprintf(lulog, "%3d %15.6f", pegy1->block, pegy1->egy);
    pqn1 = pegy1->iqn;
    for (k = 0; k < nqn; ++k) 
      fprintf(lulog,"%3d",(int)pqn1[k]);
    pegy2 = pegy1->cross;
    if (pegy2 != NULL) {
      pqn = pegy2->iqn;
      /* make back link */
      if (qncmp(nqn, pqn, pqn1) != 0){
        pegy2->cross = pegy1; 
      }
      dif = pegy1->egy - pegy2->egy;
      fprintf(lulog, " %15.6f %15.6f", pegy2->egy, dif);
      for (k = 0; k < nqn; ++k) 
        fprintf(lulog,"%3d",(int)pqn[k]);
      if (pegy2->cross != 0)
        fputc('*' ,lulog);
    } 
    fputc('\n' ,lulog);
  }    
  fclose(lulog);
  luinp = fopen(cfil[einp], "r");
  if (luinp != NULL) { 
    lulin = fopenq(cfil[elin], "w");
    /* match inp file */
    for (;;) {
      if (getlin(luinp, nqn, iqnq, iqn2, &xfrq, &xerr, &xwt,
                 tmpstr, NDLINE) < 0) {
        break;
      }
      // find block for these quanta
      for (kk = 1; kk >= 0; --kk){
        pqn = iqn2 + nqnv[kk];
        pegy2 = qnfind(jhead2.next, pqn, nqn);
        if (pegy2 != NULL) {
          pegy1 = pegy2->cross;
          if (pegy1 != NULL) {
            qncopy(nqn, pegy1->iqn, pqn); 
          }
        }
      }
      for (j = 0; j < nqn2; ++j) 
        fprintf(lulin,"%3d",(int)iqn2[j]);
      for (; j < 12; ++j) 
        fputs("   ",lulin);
      if (fabs(xerr) < 1.e-6) {
        xerr = (xerr > 0.)? 1.e-6: -1.e-6;
      }
      fprintf(lulin,"%15.6f %10.6f %10.3E\n", xfrq, xerr, xwt);
    }
    fclose(luinp);
    fclose(lulin);
  }
  while (jhead1.next != 0) {
    pnow = jhead1.next;
    jhead1.next = pnow->next;
    free(pnow);
  }
  while (jhead2.next != 0) {
    pnow = jhead2.next;
    jhead2.next = pnow->next;
    free(pnow);
  }
  return 0;
}                               /* main */

int qncmp(const int nqn, const short *qn1, const short *qn2)
{
  int j, iret;
  iret = -1;
  for (j = 0; j < nqn; ++j) {
    iret = (int) (qn1[j] - qn2[j]);
    if (iret != 0) break;
  }
  return iret;
}
PDEGY qnfind(EGY *head, const short *iqn, const int nqn)
{
  PDEGY pegy;
  if (*iqn < 0) return NULL;
  for (pegy = head; pegy != NULL; pegy = pegy->next) {
    if (qncmp(nqn, pegy->iqn, iqn) == 0) break;
  }
  return pegy;
}

void qncopy(const int nqn, const short *qn1, short *qn2)
{
  int j;
  qn2[0] = qn1[0];
  for (j = 1; j < nqn; ++j) {
    qn2[j] = qn1[j];
  }
}

void egymove(PDEGY plast, PDEGY pdest)
{ 
  PEGY psrc;
  if (plast != NULL) {
    psrc = plast->next;
    if (psrc != NULL) {
      plast->next = psrc->next;
      if (pdest != NULL) {
        psrc->next = pdest->next;
        pdest->next = psrc;
      }
    }
  }
  psrc = NULL;
}
