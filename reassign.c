#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "calpgm.h"

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
#define NDQN  10
typedef struct rule {
  /*@owned@*/ /*@null@*/ struct rule *next;
  short iqn1[NDQN], iqn2[NDQN], jbgn, jend, knt;
} RULE;
typedef /*@owned@*/ /*@null@*/ RULE *PRULE;
typedef /*@dependent@*/ /*@null@*/ RULE *DRULE;

int fixqn(const int nqn, short *iqn, int *key, const DRULE head);

int main(argc, argv)
int argc;
char *argv[];
{
#define NDLINE 128
  PRULE head;
  static int key[NDQN];
  static double dvec[NDQN+NDQN+2];
  static short iqn[NDQN+NDQN], idqn[NDQN+NDQN];
  static char line[NDLINE];
  PRULE now;
  FILE *luin, *luout;
  double xfrq, xerr, xwt;
  short *pqn1, *pqn2;
  int k, nqn, nqn2, nn, nc, jmax, nrule, k1, k2;
  if (argc < 4) {
    puts("usage: reassign keyfile infile outfile");
    exit(0);
  }
  luin = fopen(argv[1], "r");
  if (luin == NULL) exit(0);
  /* read first line of keyfile */
  if (fgetstr(line, NDLINE, luin) <= 0) exit(0);
  nn = NDQN + 2;
  for (k = 0; k < nn; ++k) 
    dvec[k] = 0.; 
  if (pcard(line, dvec, nn, NULL) < 2) exit(0);
  nqn = (int)dvec[0]; jmax = (int) dvec[1];
  if (nqn > NDQN || nqn <= 0) exit(0);
  for (k = 0; k < nqn; ++k) {
    key[k] = (int)dvec[k+2];
    idqn[k] = -1; idqn[k + nqn] = -1;
  }
  nqn2 = nqn + nqn; nn = nqn2 + 2;
  head = NULL; now = NULL; nrule = 0;
  while (fgetstr(line, NDLINE, luin) > 0) {
    /* read assignment code from keyfile */
    for (k = 0; k < nn; ++k) 
      dvec[k] = 0.; 
    dvec[1] = (double)jmax;
    if (pcard(line, dvec, nn, NULL) == 0) continue;
    now = (PRULE)mallocq(sizeof(RULE)); ++nrule;
    now->next = head; head = now;
    now->jbgn = (short)dvec[0]; now->jend = (short)dvec[1];
    now->knt = (short) nrule;
    now->iqn1[0] = 0; now->iqn2[0] = 0;
    pqn1 = now->iqn1; pqn2 = now->iqn2;
    for (k = 0; k < nqn; ++k) {
      pqn1[k] = (short) dvec[k + 2];
      pqn2[k] = (short) dvec[k + 2 + nqn];
      if (key[k] == 0) continue;
      pqn1[k] -= now->jbgn;
      pqn2[k] -= now->jbgn;
    }
  }
  fclose(luin);
  pqn1 = NULL; pqn2 = &iqn[nqn]; nn = 0;
  luin = fopen(argv[2],"r"); luout = fopen(argv[3],"w");
  if (luin == NULL) exit(0);
  if (luout == NULL) exit(0);
  for(;;) {
    nc = getlin(luin, nqn, idqn, iqn, &xfrq, &xerr, &xwt, line, NDLINE);
    if (nc < 0) {
      fclose(luin); fclose(luout);
      break;
    }
    ++nn;
    k1 = fixqn(nqn, iqn,  key, head);
    k2 = fixqn(nqn, pqn2, key, head);
    for (k = 0; k < nqn2; ++k) {
      fprintf(luout, "%3d", (int) iqn[k]);
    }
    for (k = nqn2; k < 12; ++k) {
      fputs("   ", luout);
    }
    fprintf(luout, "%s\n", line + nc);
  }
  fclose (luout); fclose(luin);
  while (head != NULL) {
    now = head->next;
    free (head);
    head = now;
  }
  head = NULL;
  return 0;
}

int fixqn(const int nqn, short *iqn, int* key, const DRULE head)
{
  RULE *now;
  short *pqn;
  int k;
  short j;
  j = 0;
  for (now = head; now != NULL; now = now->next) {
    pqn = now->iqn1; j = -1;
    for(k = 0; k < nqn; ++k){
      if (key[k] == 0) {
        if (pqn[k] != iqn[k]) break;
      } else {
        if (j < 0) { 
          j = iqn[k] - pqn[k];
          if (j < now->jbgn || j > now->jend) break;
        }
        if ((j + pqn[k]) != iqn[k]) break;
      }
    }
    if (k == nqn) break; 
  }
  if (now != NULL) {
    pqn = now->iqn2; /* match */
    for (k = 0; k < nqn; ++k){
      if (key[k] == 0) {
        iqn[k] = pqn[k];
      } else {
        iqn[k] = j + pqn[k];
      }
    }
    k = now->knt;
    return k;
  }
  return 0;
}

