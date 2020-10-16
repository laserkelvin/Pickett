#include <stdio.h>
#include <stdlib.h>
#include "catread.h"

main()
{
/*      PROGRAM TO LIST CATALOG BETWEEN FREQUENCY LIMITS */
  char fqhi[16], buf[82], *molnam, fqlow[16], resp, *xtra, line[82];
  double c, freq, err, dfqlow, dfqhi;
  FILE *flist;
  float q[7];
  int moltag, cnvt, ierr, iver, iresp, nline;
  xtra = buf + 21;
  c = 29979.2458;
  printf(" Enter list file\n");
  gets(buf);
  flist = fopen(buf, "w");
  printf(" MHz units, or  wavenumbers ? (m,w) [M] ");
  gets(line);
  cnvt = (line[0] == 'w');
  if (cnvt) {
    printf(" ENTER FREQUENCY LIMITS IN WAVENUMBERS \n");
    dfqlow = 0.;
    dfqhi = 10000.;
    gets(buf);
    iresp = sscanf(buf, "%lf %lf", &dfqlow, &dfqhi);
    dfqlow *= c;
    dfqhi *= c;
  } else {
    printf(" ENTER FREQUENCY LIMITS IN GHz \n");
    dfqlow = 0.;
    dfqhi = 10000.;
    gets(buf);
    iresp = sscanf(buf, "%lf %lf", &dfqlow, &dfqhi);
    dfqlow *= 1000.;
    dfqhi *= 1000.;
  }
  sprintf(fqlow, "%13.4f", dfqlow);
  sprintf(fqhi, "%13.4f", dfqhi);
  moltag = 0;
  while (nxtdir(&moltag)) {
    molnam = catdir(moltag, &nline, q, &iver);
    nline = catfrq(moltag, fqlow, buf);
    printf("%7d %s %d\n", moltag, molnam, nline);
    if (nline > 0)
      puts(buf);
    if (nline > 0 && strcmp(buf, fqhi) < 0) {
      printf("%s is next. OK? (y/n/a) [y] ", molnam);
      gets(line);
      resp = line[0];
      if (resp == 'a')
        exit(0);
      if (resp == 'y') {
        fprintf(flist, "              %6d %s\n", moltag, molnam);
        do {
          if (cnvt) {
            sscanf(buf, "%13.4lf%8.4lf", &freq, &err);
            err /= c;
            freq /= c;
            fprintf(flist, "%13.5f%8.5f%s\n", freq, err, xtra);
          } else {
            fputs(buf, flist);
            fputc('\n', flist);
          }
          nline++;
          ierr = catrd(moltag, nline, buf);
        } while (ierr == 0 && strcmp(buf, fqhi) < 0);
      }
    }
  }
  printf("LAST SPECIES READ\n");
}
