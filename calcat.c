/***********************************************************************/
/*   Copyright (C) 1989, California Institute of Technology            */
/*   All rights reserved.  U. S. Government Sponsorship under          */
/*   NASA Contract NAS7-918 is acknowledged.                           */

/*   Herbert M. Pickett, 20 Mar 1989                                   */
/*   Revised version in c, 22 March 1999                               */

/*   HMP, code added for binary read of variance                       */
/*   HMP, constained derivatives calculated within HAMX                */
/*   HMP, IDPAR is now REAL*8                                          */
/*   HMP, identical dipoles are summed on input                        */
/*   HMP, str output now will put out component dipoles                */
/*   HMP, large uncertainties do not inhibit output to .CAT file       */
/*   HMP, sort cat file at end                                         */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint              */ 
/*   19 July  2004: bug fix for multiple dipoles                       */

/* THIS IS A GENERALIZED INTENSITY AND FREQUENCY CALCULATOR            */

/* THE HAMILTONIAN IS ASSUMED TO BE ORGANIZED INTO BLOCKS              */
/*     SUCH THAT ADJACENT SETS OF BLOCKS ARE CONNECTED BY TRANSITIONS  */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calpgm.h"
#define MAXQNX  13
typedef struct {
  /*@owned@*/ /*@null@*/ double *eigblk;
  /*@owned@*/ /*@null@*/ double *egyblk;
  unsigned int nsizblk;
  int ixblk;
} SBLK;

/************** CALCAT interfaces ***********************************/
int qnfmt(short *iqu, int nqn, /*@out@*/char *sqn);
int simt(int isiz, int jsiz, double *s, double *t, double *u,
                double *wk);
/*@dependent@*/ SBLK *ibufof(const int iblk,
                                    const unsigned int ndel, 
                                    /*@out@*/ SBLK *blk);
SBLK *sblk_alloc(const int nstruct, const unsigned mxdm);

/********************************************************************/
#define PR_DELAY 6   /* seconds delay between informational messages */

int main(int argc, char* argv[])
{
#define NFILE 6
#define NTEMP 8
#define NCARD 130
#define NDVEC 10
 /* Initialized data */
  static double temp[8] =
    { 300., 225., 150., 75., 37.5, 18.75, 9.375, 300. };
  static double qsum[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };
  static const char *ext[NFILE] =
    { "int", "var", "out", "cat", "str", "egy" };
  enum efile {eint, evar, eout, ecat, estr, eegy};
  static char sqn[4*MAXQN + 2];
  static char warn[] = " WARNING: THERE WAS NO DIAGONALIZATION\n";
  static char headq[] = "TEMPERATURE - Q(SPIN-ROT.) - log Q(SPIN-ROT.)\n";
  /* Local variables */
  FILE *luint, *luvar, *luout, *lucat, *lustr, *luegy;
  /*@only@*/ SBLK *blk;
  /*@dependent@*/ SBLK *pblk;
  double **s;
  double *par, *derv, *dip, *egy, *sij, *var, *dedp, *teig;
  double *teigp, *egyp, *pmix, *dvec;
  double dgnf, qlog, elow, tmql, zero, qrot, strq, strr, diff,
    fqmax;
  double fqmin, starg, scomp, thrsh, strlg, strmn, thrshf, thrsh1, te,
    bigerr;
  double egymin, stcomp, cmc, dgn, fac, tmc, egx, err, frq, tmq, str, tmp;
  double telim;
  size_t nl, nlsq;
  long itag, nline, ioff;
  /*@only@*/ bcd_t *idip, *idpar;
  /*@only@*/ int *nvdip;
  /*@only@*/ int *isimag;
  /*@only@*/ int *iqnfmtv;
  int iflg, iblk, jblk, lblk, idgn, kbgn, jdgn, ndip, nbsav, nfit, igup, ixdp;
  int nsav, jxdp, itmp, jnxt, i, j, k, inblk, nbkpj, npdip, ntemp, ij, jj, jx;
  int iqnfmt, idf, itd, nqn, npar, isiz, jsiz, nsize_p, maxqn, ktsp, ibcd, ndbcd;
  int jmax, nfmt, iposv, iv, newfmt, globfmt, isneg, catqn, maxv;
  BOOL prir, prder, prfrq, preig, pregy, prstr, diag, first, ifdump;
  unsigned int ndel, maxdm;
  short iqni[MAXQN + MAXQNX];
  char *fname[NFILE+1], titl[NCARD], sgup[4];

  zero = 1.5e-38;
  bigerr = 999.9999;
  cmc = 29979.2458;
  tmc = 1.43878;
  nsize_p = maxmem(&nl);
  pregy = FALSE;
  prstr = FALSE;
  filget(argc, argv, NFILE, fname, ext);
  luout = fopenq(fname[eout], "w");
  luegy = lustr = luout;
  luint = fopenq(fname[eint], "r");
  nline = 0;
  egymin = 0.;
  /* find out how many lines are in the file */
  for (i = 0; i >= 0; i++) {
    if (fgetstr(titl, NCARD, luint) <= 0)
      break;
  }
  ndip = i - 2; /* all but first 2 lines are dipoles */
  rewind(luint);
  /* read first two lines of .int file */
  first = (fgetstr(titl, NCARD, luint) <= 0);
  if (!first) {
    chtime(titl, 82); 
    fputs(titl, stdout);
    fputs(titl, luout);
    first = (fgetstr(titl, NCARD, luint) <= 0);
  }
  dvec = (double *) mallocq((size_t)NDVEC * sizeof(double));
  dvec[0] = 0;
  dvec[1] = 999;
  dvec[2] = 1000;
  dvec[3] = 0;
  dvec[4] = 0;
  dvec[5] = -100;
  dvec[6] = -100;
  dvec[7] = 9999.99;
  dvec[8] = 300.;
  dvec[9] = -1;
  if (!first && pcard(titl, dvec, NDVEC, NULL) == 0)
    first = TRUE;
  iflg = (int) dvec[0];
  itag = (int) dvec[1];
  qrot = dvec[2];
  inblk = (int) dvec[3];
  lblk = (int) dvec[4];
  thrsh = dvec[5];
  thrsh1 = dvec[6];
  fqmax = dvec[7];
  tmq = dvec[8];
  maxv = (int) dvec[9];
  ntemp = NTEMP;
  for (k = NTEMP - 2; k >= 0; --k) {
    if (fabs(tmq - temp[k]) < 0.01) {
      --ntemp;
      break;
    }
  }
  if (ntemp == NTEMP) { /* sort temperatures in decreasing order */
    temp[NTEMP - 1] = tmq;
    for (k = NTEMP - 2; k >= 0; --k) {
      if (temp[k] > tmq)
        break;
      temp[k + 1] = temp[k];
      temp[k] = tmq;
    }
  }
  telim = log(1.e-18);
  if (first) {
    puts(" ERROR IN READING  FIRST CARDS OF .INT");
    fclose(luint);
    fclose(luout);
    exit(EXIT_FAILURE);
  }
  if (qrot < 1)
    qrot = 1;
  fprintf(luout, "ID=%6ld QSPINROT= %14.4f MIN, MAX QN= %3d %3d\n", itag,
          qrot, inblk, lblk);
  fprintf(luout, "MIN.LOG.STR= %9.3f MIN.LOG.STR(FRQ/300GHZ)**2= %9.3f",
          thrsh, thrsh1);
  fprintf(luout, " MAX FREQ %10.1f GHZ, TEMP= %8.2f\n", fqmax, tmq);
  fac = (4.16231e-5) / qrot;
  prir = (iflg >= 1000);
  /*     initialize intensities for infrared or microwave cases */
  if (prir) {
    fac *= cmc;
    cmc = 1.;
    iflg %= 1000;
    fqmin = 5e-7;
    bigerr *= 0.01;
  } else {
    fqmax *= 1000;
    fqmin = 5e-5;
  }
  if (fqmax < fqmin)
    fqmax = fqmin;
  /*   set up flags for different output options */
  prfrq = (iflg >= 100);
  iflg %= 100;
  if (iflg >= 10) {
    lustr = fopenq(fname[estr], "w");
    prstr = TRUE;
  }
  stcomp = 1. / zero;
  npdip = 1;
  if (prstr) {
    stcomp = 1e-7;
    if (iflg >= 20)
      npdip = -1;
    iflg %= 10;
  }
  ifdump = (iflg >= 5);
  prder = (iflg == 2 || iflg == 4);
  preig = (iflg >= 3);
  if (iflg != 0) {
    luegy = fopenq(fname[eegy], "w");
    pregy = TRUE;
  }
  nl = (size_t) ndip * sizeof(double);
  dip = (double *) mallocq(nl);
  dip[0] = 0.;
  nl = (size_t) ndip * sizeof(int);
  nvdip = (int *) mallocq(nl);
  nvdip[0] = 0;
  isimag = (int *) mallocq(nl);
  isimag[0] = -1;
  nl = (size_t) (ndip * NDECDIP);
  idip = (bcd_t *) mallocq(nl);
  idip[0] = (bcd_t) NDECDIP;
  k = -1; ibcd = 0;
  /*  read dipole moments */
  for (j = 0; j < ndip; ++j) {
    nvdip[j] = 0; isimag[j] = -1;
    if (fgetstr(titl, NCARD, luint) <= 0)
      break;
    jj = getbcd(titl, &idip[ibcd], NDECDIP);
    if (jj <= 0) 
      break;
    dvec[0] = 0.;
    if (pcard(&titl[jj], dvec, 1, NULL) <= 0)
      break;
    dip[j] = dvec[0];
    if (NEGBCD(idip[ibcd]) == 0 || j == 0) {
      nvdip[++k] = 1;
    } else {
      nvdip[k] += 1;
    }
    ibcd += NDECDIP;
    fprintf(luout, "%6d %s\n", k + 1, titl);
  }
  ndip = j;
  if (npdip > 0) {
    nvdip[0] = ndip;
  } else {
    npdip = k + 1;
    if (npdip > NDVEC) {
      free(dvec);
      dvec = (double *) mallocq((size_t)npdip * sizeof(double));
      dvec[0] = 0.;
    }
  }
  fclose(luint);                /* close .INT file */
  luvar = fopenq(fname[evar], "r");
  /* read var file */
  npar = 0; catqn = MAXCAT;
  if (fgetstr(titl, NCARD , luvar) > 0) {
    fputs(".VAR FILE TITLE :", luout);
    fputs(titl, luout);
    fputc('\n', luout);
    if (fgetstr(titl, NCARD, luvar) > 0) {
      dvec[0] = 0; dvec[1] = 1;
      pcard(titl, dvec, 2, NULL);
      npar = (int) dvec[0];
      if (dvec[1] < 0.) catqn = MAXQN;
    }
  }
  /* read option line(s) */
  nfmt = catqn;
  if (npar > 0 && setopt(luvar, &nfmt, &itd, &ndbcd, titl) <= 0)
    npar = 0;
  if (npar <= 0) {
    puts("Error reading .VAR file");
    exit(EXIT_FAILURE);
  }
  itmp = 10;
  iqnfmtv = (int *)mallocq((size_t) (nfmt << 1) * sizeof(int));
  iqnfmtv[0] = 0;
  nqn = setfmt(iqnfmtv, nfmt);
  if (nqn > MAXCAT) {
    catqn = nqn;
  } else {
    sqn[24] = ' ';
  }
  catqn = catqn << 1;
  memset(sqn, (int)' ', (size_t)(catqn << 1));
  iqnfmt = iqnfmtv[0];
  iposv = 0;
  if (nfmt > 1 || ((iqnfmt / 1000) & 1) != 0) 
    iposv = ((iqnfmt / 100) % 5) - 1;
  if (iposv == 0 || maxv < 0) 
    maxv = 1000;
  globfmt = 0;
  for (i = 0; i < nfmt; ++i) {
    j = iqnfmtv[i] / 1000;
    if (j != 4 && j != 5)  j = 0;
    iqnfmtv[i + nfmt] =  j;
    globfmt += j;
  }
  newfmt = iqnfmtv[nfmt];
  /* read parameters and variance */
  nl = (size_t) npar * sizeof(double);
  par = (double *) mallocq(nl);
  derv = (double *) mallocq(nl);
  nl = (size_t) (npar * ndbcd);
  idpar = (bcd_t *) mallocq(nl);
  idpar[0] = (bcd_t) ndbcd;
  strcpy(titl, "0123456789");
  itmp = getpar(luvar, luout, &nfit, &npar, idpar, par, derv, titl, 0);
  ndel = (unsigned) (nfit + 1);
  if (nfit > nsize_p || nfit <= 0) {
    puts(" memory allocation error for var matrix");
    exit(EXIT_FAILURE);
  }
  if ((nfit & 1) == 0) {
    nlsq = (size_t) ((unsigned) nfit >> 1) * sizeof(double);
    nlsq *= (size_t) ndel;
  } else {
    nlsq = (size_t) ((unsigned) ndel >> 1) * sizeof(double);
    nlsq *= (size_t) nfit;
  }
  var = (double *) mallocq(nlsq);
  getvar(luvar, nfit, var, idpar, derv, itmp);
  fclose(luvar);                /* close .VAR file */
  nbkpj = lblk;
  i = nsize_p;
  fprintf(luout, "IQNFMT = %d\n", iqnfmt);
  setblk(luout, npar, idpar, par, &nbkpj, &i);
  maxdm = (unsigned) i;
  inblk = nbkpj * inblk + 1;
  lblk = nbkpj * (lblk + 1);
  nl = (size_t) maxdm * sizeof(double);
  pmix = (double *) mallocq(nl);
  pmix[0] = 1.;
  nl *= (size_t) maxdm;
  s = (double **) mallocq((size_t) npdip * sizeof(double *));
  s[0] = (double *) mallocq(nl);
  for (k = 1; k < npdip; ++k) {
    s[k] = (double *) mallocq(nl);
  }
  /* set up DIAG,NSAV */
  setint(luout, &diag, &nsav, ndip, idip, isimag);
  if (prstr) {
    k = 0; i = nvdip[0]; ij = 0;
    /* check dipole moments */
    for (j = 1; j < ndip; ++j) {
      if (j == i) { /* new set */
        isimag[++k] = isimag[j];
        i += nvdip[k];
      } else if (isimag[j] != isimag[k]) {
        if (isimag[k] < 0) {
          isimag[k] = isimag[j];
        } else if (isimag[j] >= 0) {   
          ++ij;
        }
      }
      ibcd += NDECDIP;
    }
    if (ij > 0)
      fputs(" mixed magnetic and electric dipoles in str file\n",luout);
  }
  nsav = (nsav + 1) * nbkpj;
  if (diag && (ifdump || maxdm == 1))
    diag = 0;
  if (!ifdump && preig)
    preig = diag;
  /* set block structure */
  blk = sblk_alloc(nsav + 1, maxdm);
  nl = (size_t) maxdm;
  nl = ((size_t) ndel + nl) * nl;
  maxmem(&nlsq);
  nlsq = nlsq / nl;
  nbsav = nsav;
  if (nlsq < (size_t) nbsav)
    nbsav = (int) nlsq;
  nlsq = nl * sizeof(double);
#ifdef NDHEAPC
#if    NDHEAPC
  nl = (size_t) NDHEAPC;
  nl = nl / nlsq;
  if (nl < (size_t) nbsav)
    nbsav = (int) nl;
#endif
#endif
  if (nbsav < 2)
    nbsav = 2;
  nl = (size_t) maxdm *sizeof(double);
  /* initialize block structure */
  nlsq = nl * (size_t) maxdm;
  nl *= (size_t) ndel;
  pblk = blk;
  for (i = 0; i <= nbsav; ++i) {
    pblk->egyblk = (double *) mallocq(nl);
    if (diag) {                 /* allocate eigenvector space */
      pblk->eigblk = (double *) mallocq(nlsq);
    }
    ++pblk;
  }
  pblk = NULL;
  rqexit(1);
  /*  set up intensity constants */
  tmc = -tmc;
  tmq = tmc / (tmq * cmc);
  tmql = tmq * 0.43429448;      /*  tmql= log10(exp(1))*tmq */
  thrsh1 -= 10.9542425;
  starg = thrsh - log10(fqmax * fac);
  scomp = -38;
  strmn = (starg > scomp) ? starg : scomp;
  strmn = pow(10., strmn);
  maxqn = getqn(inblk, 0, 0, iqni, &isiz);
  if (maxqn > MAXQNX)
    maxqn = MAXQNX;
  ktsp = -1;
  lucat = fopenq(fname[ecat], "w");
  /**********************************************************************/

  /* START MAJOR LOOP OVER BLOCKS */

  for (iblk = inblk; iblk <= lblk; ++iblk) {
    first = prfrq;
    j = (iblk - 1) / nbkpj;
    if ((iblk - j * nbkpj) == 1) {
      if (rqexit(0) != 0)
        break;
      pblk = blk;
      jblk = iblk + nbkpj - nsav;
      for (i = 0; i <= nsav; ++i) {     /* free blocks not needed */
        if (pblk->ixblk < jblk)
          pblk->ixblk = 0;
        ++pblk;
      }
      if (caldelay(PR_DELAY) != 0) {
        printf(" STARTING QUANTUM %3d\n", j);
        fflush(stdout);
      }
    }
    /*  move first block to unused position */
    pblk = ibufof(-1, ndel, blk);
    /* get energy and derivatives for new block */
    getqn(iblk, 0, 0, iqni, &isiz);
    if (isiz <= 0)
      continue;
    if ((unsigned) isiz > maxdm)
      break;
    pblk->ixblk = iblk;
    pblk->nsizblk = (unsigned) isiz;
    teig = NULL;
    if (diag)
      teig = pblk->eigblk;
    if (teig == NULL)
      teig = s[0];
    egy = pblk->egyblk;
    if (egy == NULL)
      break;
    dedp = egy + isiz;
    hamx(iblk, isiz, npar, idpar, par, egy, teig, dedp, pmix, ifdump);
    /* print out energies and compute partition function */
    if (prfrq) {
      fprintf(luout, " ENERGIES FOR BLOCK NUMBER %3d, INDEX-DEGEN-ENERGY-",
              iblk);
      fputs("-EST.ERROR-MIXING-   QUANTUM NUMBERS\n", luout);
      first = TRUE;
    }
    kbgn = 0;
    teigp = teig;
    for (i = 0; i < isiz; ++i, teigp += isiz) {
      if (ifdump)
        kbgn = i;
      getqn(iblk, i + 1, nqn, iqni, &idgn);
      if (idgn <= 0) 
        continue;
      iv = iqni[iposv];
      if (iv > maxv)
        continue;
      if (i == 0 && nfmt > 1 && iv < nfmt) {
        iqnfmt = iqnfmtv[iv]; newfmt = iqnfmtv[iv + nfmt];
      }
      if (globfmt != 0) {
        ktsp = -1;
        if (newfmt != 0) 
          ktsp = iqni[nqn - 2];
        getqn(iblk, i + 1, maxqn, iqni, &idgn);
      }
      if (diag) {   /* check if eigenvector is zero */
        if (fabs(teigp[i]) < 0.01) {
          k = (int) idamax(isiz, teigp, 1);
          if (fabs(teigp[k]) < 0.01)
            idgn = 0;
        }
      }
      /* calculate uncertainty */
      dcopy(nfit, &dedp[i], isiz, derv, 1);
      err = calerr(nfit, var, derv) / cmc;
      if (err > bigerr)
        err = bigerr;
      egx = egy[i] / cmc;
      if (egymin > egx) {     /*  adjust energy zero */
        te = tmc * (egymin - egx);
        egymin = egx;
        for (k = 0; k < ntemp; ++k) {
          tmp = te / temp[k];
          if (tmp < telim)
            tmp = telim;
          qsum[k] *= exp(tmp);
        }
      }
      /*  accumulate partition function values */
      te = tmc * (egx - egymin);
      dgn = (double) idgn;
      for (k = 0; k < ntemp; ++k) {
        tmp = te / temp[k];
        if (tmp < telim)
          break;
        qsum[k] += dgn * exp(tmp);
      }
      if (prfrq) {
        fprintf(luout, "%6d %4d %17.6f %17.6f %10.6f",
                iblk, i + 1, egx, err, pmix[i]);
        if (globfmt != 0) 
          fprintf(luout, " (%3d)", ktsp);
        for (k = 0; k < maxqn; ++k) {
          fprintf(luout, "%3d", (int) iqni[k]);
        }
        fputc('\n', luout);
      }
      if (pregy) {
        fprintf(luegy, "%6d %4d %17.6f %17.6f %10.6f %4d:",
                iblk, i + 1, egx, err, pmix[i], idgn);
        if (globfmt != 0) 
          fprintf(luegy, "(%3d)", ktsp);
        for (k = 0; k < maxqn; ++k) {
          fprintf(luegy, "%3d", (int) iqni[k]);
        }
        fputc('\n', luegy);
      }
      if (prder) {
        j = 0;
        for (k = 0; k < nfit; ++k) {
          if (j == 0)
            fputs("        ", luegy);
          fprintf(luegy, " %5d %13.5E", k + 1, derv[k]);
          if ((++j) == 6) {
            fputc('\n', luegy);
            j = 0;
          }
        }
        if (j > 0)
          fputc('\n', luegy);
      }
      if (preig) {
        j = 0;
        for (k = kbgn; k < isiz; ++k) {
          if (j == 0)
            fputs("        ", luegy);
          fprintf(luegy, " %5d %13.5E", k + 1, teigp[k]);
          if ((++j) == 6) {
            fputc('\n', luegy);
            j = 0;
          }
        }
        if (j > 0)
          fputc('\n', luegy);
      }
    }
    jblk = 0;
    do {                        /* loop over previous blocks for intensity */
      jnxt = iblk + 1;
      pblk = blk;
      jj = 0;
      for (j = 0; j <= nsav; ++j) {     /*  find oldest block .GT. jblk */
        jx = pblk->ixblk;
        if (jx > jblk && jx < jnxt) {
          jnxt = jx;
          jj = j;
        }
        ++pblk;
      }
      jblk = jnxt;
      pblk = ibufof(jj, ndel, blk);
      jsiz = (int) pblk->nsizblk;
      teigp = pblk->eigblk;
      if (teigp == NULL)
        teigp = s[0];
      idgn = 0;
      ij = 0;
      for (i = 0; i < npdip; ++i) {     /*  get intensity */
        sij = s[i];
        idf = intens(iblk, isiz, jblk, jsiz, nvdip[i],
                     &idip[ij * NDECDIP], &dip[ij], sij);
        if (idf != 0) {
          if (diag)
            simt(isiz, jsiz, sij, teig, teigp, pmix);
          idgn = idf;
        }
        sij = NULL;
        ij += nvdip[i];
      }
      if (idgn <= 0)
        continue;
      dgnf = (double) idgn;
      egyp = pblk->egyblk;
      if (egyp == NULL)
        continue;
      idf = (int) maxdm;
      jmax = jsiz;
      for (i = 0; i < isiz; ++i) {
        getqn(iblk, i + 1, nqn, iqni, &idgn);
        if (idgn <= 0 || iqni[iposv] > (short)maxv)
          continue;
        dgn = idgn / dgnf;
        if (iblk == jblk)
          jmax = i + 1;
        for (j = 0; j < jmax; ++j) {
          if (i - j + idf < 0)
            break;
          getqn(jblk, j + 1, nqn, &iqni[MAXQN], &jdgn);
          if (jdgn <= 0 || iqni[iposv + MAXQN] > (short) maxv)
            continue;
          frq = egy[i] - egyp[j];
          ixdp = i;
          jxdp = j;
          for (k = 0; k < nfit; ++k) {
            ixdp += isiz;
            jxdp += jsiz;
            derv[k] = egy[ixdp] - egyp[jxdp];
          }
          isneg = 0;
          if (frq < 0.) {
            isneg = 1;
            frq = -frq;
            qnfmt(&iqni[MAXQN], nqn, sqn);
            qnfmt(iqni, nqn, sqn + catqn);
            igup = jdgn;
            elow = egy[i];
          } else {
            qnfmt(iqni, nqn, sqn);
            qnfmt(&iqni[MAXQN], nqn, sqn + catqn);
            igup = idgn;
            elow = egy[i] - frq;
          }
          ioff = isiz;
          ioff = i + ioff * j;
          strr = 0.; strq = 0.;
          for (ij = 0; ij < npdip; ++ij) {
            sij = s[ij];
            str = sij[ioff];
            dvec[ij] = str; 
            strr += str;
            strq += str * str;
          }
          if (strq > stcomp) {
            for (ij = 0; ij < npdip; ++ij) {
              str = dvec[ij];
              if (isneg != 0 && isimag[ij] > 0) str = -str;
              fprintf(lustr, "%15.4f%15.6E %4d %s %4d\n", frq, str,
                      iqnfmt, sqn, ij + 1);
            }
          }
          strr *= strr;
          if (frq > fqmax)
            continue;
          if (frq < fqmin)
            continue;
          if (strr < strmn)
            continue;
          str = dgn * strr * fac * frq * (1. - exp(tmq * frq));
          strlg = log10(str + zero) + tmql * elow;
          if (strlg < thrsh)
            continue;
          thrshf = thrsh1 + log10(frq + zero) * 2.;
          diff = thrsh - thrshf;
          if (fabs(diff) < 4.) {
            diff = pow(10., diff) + 1.;
            thrshf = log10(diff) + thrshf;
          }
          if (strlg < thrshf)
            continue;
          /* calculate errors */
          err = calerr(nfit, var, derv);
          if (err > bigerr)
            err = bigerr;
          elow /= cmc;
          if (first) {
            fputs(" FREQUENCY-EST.ERROR.-LINE.STR. DIP**2-LGSTR.-ITD,",
                  luout);
            fputs("-GUP-I.D.-QNFORM-QUANTUM NUMBERS\n", luout);
            first = FALSE;
          }
          if (prfrq) {
            fprintf(luout, "%13.4f %8.4f %12.5E %8.4f %2d",
                    frq, err, strr, strlg, itd);
            fprintf(luout, "%10.4f %3d %7ld %4d %s\n",
                    elow, igup, itag, iqnfmt, sqn);
          }
          ++nline;
          gupfmt(igup, sgup); sgup[3] = '\0';
          if (prir) {
            fprintf(lucat, "%13.6f%8.6f%8.4f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, itd, elow, sgup, itag, iqnfmt);
          } else if (frq < 99999999.) {
            fprintf(lucat, "%13.4f%8.4f%8.4f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, itd, elow, sgup, itag, iqnfmt);
          } else {
            fprintf(lucat, "%13.3f%8.3f%8.3f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, itd, elow, sgup, itag, iqnfmt);
          }
          fputs(sqn,lucat); fputc('\n',lucat);
        }
      }
    } while (jblk != iblk);
  }
  teig = teigp = NULL;
  /****************************************************************/

  /* print partition function */

  if (ifdump) {
    fputs(warn, stdout);
    fputs(warn, luout);
  }
  printf(       "INITIAL Q = %14.4f, NEW Q IS RELATIVE TO MIN.EGY.= %14.4f\n",
                qrot, egymin);
  fprintf(luout,"INITIAL Q = %14.4f, NEW Q IS RELATIVE TO MIN.EGY.= %14.4f\n",
          qrot, egymin);
  printf(" NUMBER OF LINES = %6ld\n", nline);
  fprintf(luout," NUMBER OF LINES = %6ld\n", nline);
  fputs(headq, stdout);
  fputs(headq, luout);
  for (i = 0; i < ntemp; ++i) {
    qlog = -100;
    if (qsum[i] > zero)
      qlog = log10(qsum[i]);
    printf(       " %10.3f %14.4f %9.4f\n", temp[i], qsum[i], qlog);
    fprintf(luout," %10.3f %14.4f %9.4f\n", temp[i], qsum[i], qlog);
  }
  if (pregy)
    fclose(luegy);
  if (prstr)
    fclose(lustr);
  setblk(luout, 0, idpar, par, &nbkpj, &i);     /* release storage */
  fclose(luout);
  fclose(lucat);
  /* free up memory */
  pblk = blk;
  for (i = 0; i <= nsav; ++i) {
    if (pblk->eigblk != NULL)
      free(pblk->eigblk);
    if (pblk->egyblk != NULL)
      free(pblk->egyblk);
    ++pblk;
  }
  pblk = NULL;
  free(blk);
  sij = NULL;
  for (i = 0; i < npdip; ++i) {
    free(s[i]); s[i] = NULL;
  }
  free(s);
  free(pmix);
  if (var != NULL)
    free(var);
  free(idpar);
  free(derv);
  free(par);
  free(iqnfmtv);
  free(idip); 
  free(isimag); 
  free(nvdip);
  free(dip);
  free(dvec);
  /* sort cat file */
  sortn(fname[ecat], fname[ecat], FALSE);
  return 0;
}                               /* main */

int qnfmt(iqu, nqn, sqn)
short *iqu;
int nqn;
char *sqn;
{ /* subroutine to do quick conversion of quantum integers to characters */
  static int czero = (int) '0';
  int i, k, ix;
  char tqn1, tqn2;

  for (k = 0; k < nqn; ++k) {
    i = iqu[k];
    if (i < 0) {
      i = -i;
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char) (i + czero);
      if (ix == 0) {
        tqn2 = '-';
      } else if (ix < 27) {
        tqn2 = (char) (ix + ((int) 'a' - 1));
      } else {
        tqn1 = tqn2 = '*';
      }
    } else {
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char) (i + czero);
      if (ix == 0) {
        tqn2 = ' ';
      } else if (ix < 10) {
        tqn2 = (char) (ix + czero);
      } else if (ix < 36) {
        tqn2 = (char) (ix + ((int) 'A' - 10));
      } else {
        tqn1 = tqn2 = '*';
      }
    }
    ix = k + k;
    sqn[ix] = tqn2;
    sqn[ix + 1] = tqn1;
  }
  return 0;
}                               /* qnfmt */

int simt(isiz, jsiz, s, t, u, wk)
int isiz, jsiz;
double *s, *t, *u, *wk;
{
  int i, j;
  double *sij, *si, *tc, *uc;
  /*    subroutine to do similarity transform */
  /*      S' = transpose(T) * S * U */
  if (isiz > 1) { /* do transform on the left */
    sij = s;
    for (j = 0; j < jsiz; ++j) {
      dcopy(isiz, sij, 1, wk, 1);
      tc = t;
      for (i = 0; i < isiz; ++i) {
        *sij = ddot(isiz, tc, 1, wk, 1);
        tc += isiz;
        ++sij;
      }
    }
  }
  if (jsiz > 1 && u != NULL) {  /* do transform on the right */
    si = s;
    for (i = 0; i < isiz; ++i) {
      dcopy(jsiz, si, isiz, wk, 1);
      sij = si;
      uc = u;
      for (j = 0; j < jsiz; ++j) {
        *sij = ddot(jsiz, wk, 1, uc, 1);
        uc += jsiz;
        sij += isiz;
      }
      ++si;
    }
  }
  return 0;
}                               /* simt */

SBLK *ibufof(ipos, ndel, blk)
const int ipos;
const unsigned int ndel;
SBLK *blk;
{
  /* Initialized data */
  static FILE *scratch;
  static long maxrec, lsizb;
  static int mempos, orgpos, nbsav;
  static unsigned int maxdm;

  /* Local variables */
  /*@dependent@*/ SBLK *pmblk, *pdblk;
  /*@dependent@*/ double *pvec;
  long ldisk, len, n, lret;
  int k;

  if (ipos >= 0) {
    pmblk = &blk[ipos];
    if (pmblk->egyblk != NULL)
      return pmblk;
  }
  k = mempos;
  mempos = 0;
  pmblk = blk;
  while (pmblk->ixblk > 0) {    /* find first empty block */
    ++mempos;
    ++pmblk;
  }
  if (pmblk->egyblk == NULL) {  /* write memory contents to disk */
    if (maxrec == 0) {
      scratch = tmpfile();
      if (scratch == NULL) {
        puts("scratch file open error");
        exit(EXIT_FAILURE);
      }
      nbsav = mempos;
      maxdm = pmblk->nsizblk;
      ldisk = (long) ndel;
      if (blk->eigblk != NULL)
        ldisk += (long) maxdm;
      lsizb = (long) maxdm *ldisk * (int) sizeof(double);
    }
    ldisk = (long) (mempos - nbsav);
    if (ldisk != 0)
      ldisk *= lsizb;
    fseek(scratch, ldisk, SEEK_SET);
    mempos = k + 1;
    if (mempos == orgpos)
      ++mempos;
    if (mempos >= nbsav)
      mempos = (orgpos != 0) ? 0 : 1;
    pdblk = pmblk;
    pmblk = &blk[mempos];
    len = (long) ndel;
    n = (long) pmblk->nsizblk;
    len *= n;
    pvec = pmblk->egyblk;
    lret = 0;
    if (pvec != NULL) {
      lret = (long) fwrite(pvec, sizeof(double), (size_t) len, scratch);
      if (lret == len && pmblk->eigblk != NULL) {
        if (ldisk == maxrec)
          n = (long) maxdm;
        len = n * n;
        lret = (long) fwrite(pmblk->eigblk, sizeof(double), (size_t) len,
                             scratch);
      }
      if (ldisk == maxrec && lret == len) {
        /*  fill in gaps with junk FIRST time */
        n = (long) maxdm - (long) pmblk->nsizblk;
        len = (long) ndel;
        if (n > 0) {
          len *= n;
          lret =
              (long) fwrite(pvec, sizeof(double), (size_t) len, scratch);
        }
        maxrec += lsizb;
      }
    }
    if (lret != len) {
      puts("heap write error ");
      exit(EXIT_FAILURE);
    }
    pdblk->nsizblk = pmblk->nsizblk;
    pdblk->ixblk = pmblk->ixblk;
    pmblk->ixblk = 0;
  }
  if (ipos < 0) {
    orgpos = mempos;
  } else {                      /* read from disk */
    ldisk = (long) (ipos - nbsav);
    if (ldisk != 0)
      ldisk *= lsizb;
    fseek(scratch, ldisk, SEEK_SET);
    pdblk = &blk[ipos];
    len = (long) ndel;
    n = (long) pdblk->nsizblk;
    len *= n;
    pvec = pmblk->egyblk;
    lret = 0;
    if (pvec != NULL) {
      lret = (long) fread(pvec, sizeof(double), (size_t) len, scratch);
      if (lret == len && (pvec = pmblk->eigblk) != NULL) {
        len = n * n;
        lret = (long) fread(pvec, sizeof(double), (size_t) len, scratch);
      }
    }
    if (lret != len) {
      puts("heap read error ");
      exit(EXIT_FAILURE);
    }
    pmblk->nsizblk = pdblk->nsizblk;
    pmblk->ixblk = pdblk->ixblk;
    pdblk->ixblk = 0;
  }
  return pmblk;
}                               /* ibufof */

SBLK *sblk_alloc(nstruct, mxdem)
const int nstruct;
const unsigned mxdem;
{
  SBLK *pret, *pblk;
  int k;
  pret = (SBLK *) mallocq((size_t) nstruct * sizeof(SBLK));
  pblk = pret;
  k = nstruct;
  do {
    pblk->eigblk = NULL;
    pblk->egyblk = NULL;
    pblk->ixblk = 0;
    pblk->nsizblk = mxdem;
    ++pblk;
  } while (--k > 0);
  return pret;
}

