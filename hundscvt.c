/* hundscvt.f -- translated by f2c (version 19970219). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calpgm.h"

int main(argc, argv)
int argc;
char *argv[];
{
#define NFNAM 36
  static int icv[4] = {1,2,4,0};
  static char infil[NFNAM], outfil[NFNAM], work[82], fmt[5];

  FILE *luin, *luout;
  char *finp, *fout, *pstr;
  int iqfmt, iq, nqn, ipj, i, k, wasb;
  char omega;

  if (argc >= 2){
    finp = argv[1];
  }else{
    finp = infil;
    puts("ENTER CATALOG INPUT FILE NAME");
    fgetstr(finp, NFNAM, stdin);
  }  
  if (argc >= 3){
    fout = argv[2];
  }else{
    fout = outfil;
    puts("ENTER CATALOG OUTPUT FILE NAME");
    fgetstr(fout, NFNAM, stdin);
  }
  luin = fopen(finp,"r");
  if (luin == NULL || fgetstr(work, 82, luin) < 55) {
    puts("bad input file");
    exit (1);
  }
  rewind(luin);
  luout = fopen(fout,"w");
  if (luout == NULL) {
    puts("bad output file");
    exit (1);
  }
  memcpy(fmt,&work[51],4); fmt[4] = 0;
  if (fmt[2] == ' ') fmt[2] = 0;
  iqfmt = atoi(fmt); nqn = iqfmt % 10; iq = iqfmt / 100;
  ipj = 4;
  if (fmt[0] == '1'){
    --nqn; ipj = 6;
  }
  if (nqn < 3) iq = 0;
  if (iq == 2 || iq == 13){ /* case b to a */
    wasb = 1; fmt[1] += 6; fmt[2] -= icv[nqn-3];
  }else if (iq == 8 || iq == 19){ /* case a to b */
    wasb = 0; fmt[1] -= 6; fmt[2] += icv[nqn-3];
  }else{
    puts("bad quantum format");
    exit(1);
  }
  while (fgetstr(work, 82, luin) > 70) {
    memcpy(&work[51], fmt, 4);
    for (i = 55; i < 70; i += 12) { /* check both states */ 
      pstr = &work[i];
      if (wasb) {
        if (strncmp(pstr, &pstr[ipj], 2)){ /* Omega = 3/2 */
          pstr[0] = pstr[ipj]; pstr[1] = pstr[ipj+1]; omega = '2';
        }else{
          omega = '1';
        }
        for (k = ipj - 1; k >= 2; --k) pstr[k+2] = pstr[k];
        pstr[2] = ' '; pstr[3] = omega;
      }else{
        omega = pstr[3];
        for (k = 2; k < ipj; ++k) pstr[k] = pstr[k+2];
        pstr[ipj] = pstr[0]; pstr[ipj+1] = pstr[1];
        if (omega == '2') { /* N = J - 1/2 */
          if (pstr[1] == '0') {
            if (pstr[0] == 'A') 
              pstr[0] = '9';
            else if (pstr[0] == '1')
              pstr[0] = ' ';
            else
              pstr[0] -= 1;
            pstr[1] = '9';
          }else{
            pstr[1] -= 1;
          }
        }  
      }
    }
    fputs(work, luout); fputc('\n',luout);
  }
  fclose(luout); fclose(luin);
  return 0;
} /* main */
