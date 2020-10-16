/* calbak.f -- translated by f2c */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calpgm.h"

int main(int argc, char *argv[])
{
#define NDLINE 130
  static double dval[8];
  static int fmt[8] = { 13, 8, 8, 2, 10, 3, 7, 4 };
  static char line[NDLINE];
  FILE *finp, *fout;
  double errx, fq, str;
  int nqn, nqnx, i, iqfmt, id, pqnlo, pqnhi;
  short iqn[20];

  /* PROGRAM TO TAKE CATALOG FORMAT AND CREATE FILE OF EXP. LINES FOR CALFIT */
  if (argc > 1) {
    finp = fopen(argv[1],"r");
  } else {
    puts("ENTER MERGED CATALOG NAME");
    fgetstr(line, NDLINE, stdin);
    finp = fopen(line, "r");
  }
  if (finp == NULL)
    exit(1);
  if (argc > 2) {
    fout = fopen(argv[2],"w");
  } else {
    puts("ENTER OUTPUT LIN FILE NAME");
    fgetstr(line, NDLINE, stdin);
    fout = fopen(line, "w");
  }
  if (fout == NULL)
    exit(1);
  puts("reading catalog file");
  nqn = nqnx = 0; pqnhi = 55; pqnlo = pqnhi + 12;
  while (fgetstr(line, NDLINE, finp) > 55) {
    pcard(line, dval, 8, fmt);
    fq = dval[0];
    errx = dval[1];
    str = dval[2];
    id = (int) dval[6];
    iqfmt = (int) dval[7];
    if (id >= 0)
      continue;
    if (nqn == 0) {
      nqn = iqfmt % 10;
      if (nqn == 0) nqn = 10;
      nqnx = nqn << 1;
      if (nqnx > 12) pqnlo = pqnhi + nqnx;
    }
    readqn(line + pqnhi, iqn, nqn);
    readqn(line + pqnlo, &iqn[nqn], nqn);
    for (i = 0; i < nqnx; ++i)
      fprintf(fout, "%3d", (int) iqn[i]);
    for (i = nqnx; i < 12; ++i)
      fputs("   ", fout);
    str = pow(10., str);
    if (errx > 0.0001) {
      fprintf(fout, "%15.4f %10.4f %10.3E\n", fq, errx, str);
    } else {
      fprintf(fout, "%15.6f %9.6f %10.3E\n", fq, errx, str);
    }
  }
  return 0;
}                               /* main */
