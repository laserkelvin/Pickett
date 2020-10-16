/* iambak.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calpgm.h"
#include "readopt.h"
#define NLABL 10
#define MNEG   1
#define MDIP   2
#define MUSED  4
#define MWN    8
 
typedef struct spar {
  double parval, parerr, parerrv;
  int  idval,mv,nv,iflg,v21;
  char labl[NLABL + 1];
} SIPAR;

static double zero = 0.;

int main(argc, argv)
int argc;
char *argv[];
{
#define NFNAM 36
#define NDBUF 130
#define NDBCD 11
#define NDPAR1 7
#define ND2TOP 7
#define NDVAL  8  /* NDVAL >= ND2TOP  and NDVAL >= NDPARI*/
  static char infil1[NFNAM], infil2[NFNAM], buf[NDBUF];
  static int ival[NDPAR1-1], ifval[10], irho[10], icval[10];
  static SIPAR par0;
  /*@out@*/ SIPAR *par;
  SIPAR *parp, *parc;
  FILE *luin1,*luin2, *luout;
  double *ptmp;
  short *kmap;
  double dval[NDVAL], pref, rho, del, tmp, err, clight;
  char *fin1, *fin2, *fout, *pstr;
  long idtmp, nxpar, ixpar;
  int npar, ipar, kpar, ndec, nvib, iv1, iv2, mv, nv, iend, ifc;
  int k, nval, ivv, ndecfac, idrho, nxv, nvib2, nsym, kk;
  bcd_t bcdpar[NDBCD];

  /* usage: calbak old.iam old.par new.iam */
  if (argc >= 2){
    fin1 = argv[1];
  }else{
    fin1 = infil1;
    puts("ENTER IAM INPUT FILE NAME");
    fgetstr(fin1, NFNAM, stdin);
  }  
  if (argc >= 3){
    fin2 = argv[2];
  }else{
    fin2 = infil2;
    puts("ENTER PAR INPUT FILE NAME");
    fgetstr(fin2, NFNAM, stdin);
  }
  if (argc >= 4){
    fout = argv[3];
  }else{
    fout = fin1;
    /* puts("ENTER IAM OUTPUT FILE NAME"); */
    /* fgetstr(fout, NFNAM, stdin); */
  }
  luin1 = fopen(fin1,"r");
  if (luin1 == NULL || fgetstr(buf, NDBUF, luin1) < 0) {
    puts("bad input file #1");
    exit (EXIT_FAILURE);
  }
  for (npar = 0; npar < 1000; ++npar) {
    if (fgetstr(buf, NDBUF, luin1) < 0) break;
  }
  rewind(luin1);
  clight = 29979.2458;
  bcdpar[0] = (bcd_t) NDBCD;
  par0.parval = 0.; par0.parerr = 1.e+36; par0.parerrv = 0.;
  /* read old iam file */
  if (fgetstr(buf, NDBUF, luin1) < 0) exit(EXIT_FAILURE);
  dcopy(NDPAR1, &zero, 0, dval, 1);
  nval = pcard(buf, dval, NDPAR1, NULL) - 1;
  rho = dval[0]; 
  for (k = 0; k < nval; ++k) 
    ival[k] = (int) dval[k + 1];
  nvib2 = ival[4]; nsym = ival[0];
  if (nval < 5) 
    nvib2 = 0;
  npar -= nvib2;
  par = (SIPAR *)mallocq((size_t) (npar + 1) * sizeof(SIPAR));
  memcpy(par, &par0, sizeof(SIPAR));
  kmap = (short *)mallocq((size_t)(nvib2 * ND2TOP + 1) * sizeof(short));
  kmap[0] = 0; 
  ivv = 0;
  for (k = 0; k < nvib2; ++k) {
    if (fgetstr(buf, NDBUF, luin1) < 0) break; 
    dcopy(ND2TOP, &zero, 0, dval, 1);
    if (pcard(buf,dval,ND2TOP,NULL) < 0) continue;
    for (kk = 0; kk < ND2TOP; ++kk) { 
      kmap[ivv] = (short) dval[kk];
      ++ivv;
    }
  }
  for (k = 0; k < 10; ++k) {
    ifval[k] = -1; icval[k] = -1; irho[k] = - 1;
  }
  idrho = 9100; irho[9] = 0;
  if (ival[2] > 0 && ival[3] < 0) { /*  K sampled FC */
    idrho = -1; irho[9] = -1; par->iflg = MUSED;
  }
  par->parval = rho; par->parerr = par->parerrv = 1.e-6;
  par->idval = idrho; par->v21 = 99;
  strcpy(par->labl,"rho");
  nxv = 0;
  pref = 1.; kpar = 1;
  for (ipar = 1; ipar <= npar; ++ipar) {
    if (fgetstr(buf, NDBUF, luin1) < 0) break;
    parp = &par[kpar]; 
    parp->iflg = 0;
    dval[1] = 0.; dval[0] = 0.; dval[2] = 0.; dval[3] = 1.; dval[4] = 1.e+30;
    if (pcard(buf,dval,5,NULL) < 4) continue;
    idtmp = (long) dval[0]; 
    parp->nv = (int) dval[1];
    mv = (int) dval[2]; parp->mv = mv;
    parp->parval = dval[3];
    parp->parerr = dval[4];
    parp->parerrv = 0.;
    if (dval[4] < 0) {
      parp->parval *= clight; parp->parerr *= -clight;
      parp->iflg |= MWN;
    }
    if (idtmp < 0) {
      idtmp = -idtmp; parp->iflg |= MNEG;
      parp->parval /= pref;
    } else {
      pref = parp->parval;
    }
    if (mv >= 100 || (idtmp < 9L && mv < 100)) {
      parp->iflg |= (MDIP | MUSED);
      parp->idval = (int) idtmp; parp->v21 = 0;
    } else {
      ivv = (int)(idtmp % 100); idtmp /= 100;
      parp->idval = (int) idtmp; parp->v21 = ivv;
      iv2 = ivv / 10; iv1 = ivv - iv2 * 10;
      if (iv1 < 9 && iv1 > nxv) 
        nxv = iv1;
      if (iv2 < 9 && iv2 > nxv) 
        nxv = iv2;
      if (parp->nv == 0 && iv1 == iv2) {
        if (mv == 2 && parp->idval == 0 && ifval[iv1] < 0)
          ifval[iv1] = kpar;
        if (mv == 1 && parp->idval == 2000 && icval[iv1] < 0)
          icval[iv1] = kpar;
        if (mv == 0 && parp->idval == idrho && irho[iv1] < 0) {
          if (ivv == 99) continue;
          irho[iv1] = kpar;
        }
      }
    }
    pstr = strchr(buf, '/');
    if (pstr == NULL) {
      parp->labl[0] = '\0';
    } else {
      memcpy(parp->labl, pstr + 1, NLABL);
      parp->labl[NLABL] = '\0';
    }
    ++kpar;
  }
  npar = kpar; fclose(luin1);
  luin2 = fopen(fin2,"r");
  if (luin2 == NULL || fgetstr(buf, NDBUF, luin2) < 0) {
    puts("bad input file #2");
    exit (EXIT_FAILURE);
  }
  /* read old par file */
  if (npar <= 0 || fgetstr(buf, NDBUF, luin2) < 0) exit(EXIT_FAILURE);
  pcard(buf,dval,1,NULL);
  nxpar = (long) dval[0];
  nvib = 0;
  for (;;) {/* read option cards */
    if (fgetstr(buf, NDBUF, luin2) < 0) break;
    dval [0] = 0.; dval[7] = 1.;
    if (readopt(buf, dval, 8) <= 0) break;
    if (nvib == 0) nvib = (int)fabs(dval[0]);
    if (dval[7] > -0.5) break;
  }
  if (nvib >= 99)  exit(EXIT_FAILURE);
  ndecfac = 10;
  if (nvib > 9)
    ndecfac = 100;
  ndec = ndecfac - 1; ndecfac *= ndecfac; 
  for (ixpar = 0; ixpar < nxpar; ++ixpar) {
    if (fgetstr(buf, NDBUF, luin2) < 0) break;
    iend = getbcd(buf, bcdpar, NDBCD);
    if (NEGBCD(bcdpar[0]) != 0) continue;
    iv1 = bcd2i(bcdpar[1]); ifc = 5;
    if (nvib <= 9) {
      iv2 = iv1 / 10; iv1 = iv1 - iv2 * 10; 
    } else {
      iv2 = bcd2i(bcdpar[2]); ++ifc;
    } 
    mv = iv2;
    idtmp = (int)(bcdpar[ifc] & 0x0f);
    idtmp = bcd2i(bcdpar[ifc - 1]) + idtmp * 100;
    idtmp = bcd2i(bcdpar[ifc - 2]) + idtmp * 100;
    idtmp = bcd2i(bcdpar[ifc - 3]) + idtmp * 100;
    if (iv1 < ndec) continue;
    nv = (int)(bcdpar[ifc] >> 4) & 0x0f;
    nv += (int)(bcdpar[ifc + 1] & 0x0f) * 10;
    ipar = (nv - 1) / 10;
    nv -= (ipar >> 1) * 10;
    if ((ipar & 1) != 0) 
      nv = 10 - nv;
    kpar = (int) idtmp;
    for (ipar = 1; ipar < npar; ++ipar) {
      parp = &par[ipar];
      if ((parp->iflg & (MDIP | MUSED)) != 0) continue; 
      if (parp->idval != kpar) continue;
      if (parp->mv != mv) continue;
      if (parp->nv != nv) continue;
      parp->iflg |= MUSED;
      dval[0] = parp->parval; dval[1] = 0.;
      pcard(&buf[iend], dval, 2, NULL);
      del = dval[0] - parp->parval; parp->parval = dval[0]; parp->parerrv = dval[1];
      k = 0;
      if (idtmp == 2000 && nv == 0 && mv == 1) k = 1;
      if (k == 0 && idtmp == 0 && nv > 0 && mv == 0) k = 1;
      if (k == 0 && idtmp == 0 && nv == 0 && mv == 2) k = 1;
      if (k != 0) {  
        printf("%2d %3d %2d, new=%12.3f, change=%12.3f ", parp->v21,
               nv, mv, parp->parval, del);
        err = parp->parerrv;
        if (fabs(err) < 1.e+08) 
          printf("err= %12.3f", err);
        printf(" %s\n", parp->labl);
      }
      break;
    }
  }
  parp = NULL; fclose(luin2);
  /* resolve defaults */
  for (k = 0; k <= nxv; ++k) {
    ivv = irho[k];
    if (ivv  < 0) ivv  = irho[9]; 
    if (ivv >= 0) par[ivv].iflg &= (MNEG | MDIP); /* clear MUSED */
    irho[k] = ivv;
    if (icval[k] < 0) icval[k] = icval[9];
    if (ifval[k] < 0) ifval[k] = ifval[9];
    if (icval[k] < 0 || ifval[k] < 0) irho[k] = -1;
  }
  /* update rho */
  for (k = 0; k <= nxv; ++k) {
    ivv = irho[k];
    if (ivv < 0) continue;
    parp = &par[ivv];
    if ((parp->iflg & MUSED) != 0) continue;
    parc = &par[ifval[k]];
    del = parc->parval;
    if (del > 1.) {
      err = parc->parerrv;
      parp->iflg |= MUSED;
      parc = &par[icval[k]];
      ptmp = &parc->parval;
      tmp = 0.5 / del; err *= (*ptmp) / del;
      del = parc->parerrv; 
      err = tmp * sqrt(err * err + del * del);
      if (fabs(err) > 1.e+08) err = 1.e+08;
      del = -tmp * (*ptmp); *ptmp = 1.e-06;
      parc = &par[irho[k]]; 
      ptmp = &parc->parval;
      tmp = (*ptmp) + del; *ptmp = tmp;
      if (ivv == 0) 
        rho = tmp;
      printf("new rho(%2d)= %12.9f, change= %12.9f (%12.9f)\n", k, tmp, del, err);
    }
  }
  /* write to output */
  luout = fopen(fout,"w");
  if (luout == NULL) {
    puts("bad output file");
    exit (EXIT_FAILURE);
  }
  fprintf(luout, "%12.9f ", rho);
  for (k = 0; k < nval; ++k) 
    fprintf(luout," %3d", ival[k]);
  fputc('\n',luout);
  ivv = 0;
  for (k = 0; k < nvib2; ++k) {
    for (kk = 0; kk < ND2TOP; ++kk) {
      fprintf(luout, "%5d", (int) kmap[ivv]);
      ++ivv;
    }
    fputc('\n',luout);
  }
  for (ipar = 1; ipar < npar; ++ipar) {
    parp = &par[ipar];
    idtmp = parp->idval;
    if (idtmp < 0) continue;
    if ((parp->iflg & MDIP) == 0) {
      idtmp = 100 * idtmp + parp->v21;
      if ((parp->iflg & MNEG) != 0) {
        parp->parval *= pref; idtmp = -idtmp;
      } else  {
        pref = parp->parval;
      }
       
    } else if ((parp->iflg & MNEG) != 0) {
      idtmp = -idtmp;
    }
    err = parp->parerr;
    if ((parp->iflg & MWN) != 0) {
      parp->parval /= clight;
      parp->parerr /= -clight;
    }
    fprintf(luout, "%12ld  %3d %2d %23.15e", 
            idtmp, parp->nv, parp->mv, parp->parval);
    if (parp->parerr < 1.e+30) 
      fprintf(luout, "%12.3e ",parp->parerr);
    fprintf(luout, " /%s\n", parp->labl);
  }
  parp = NULL;
  free(par);
  fclose(luout);
  exit(EXIT_SUCCESS);
} /* main */
