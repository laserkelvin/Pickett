/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: read option cards with fgetstr */
/*   30 Dec.  1999: include changes for dlsq */
/*   10 Oct.  2001: change fit diverging code */
/*   21 Sept. 2002: fix NRJ criterion */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */ 

/**************************************************************************/
/*                                                                        */
/*   THIS IS A GENERALIZED LINE FITTING PROGRAM                           */
/*   IT FITS LINES TO PARAMETERS IN A MODEL HAMILTONIAN                   */
/*   BLENDED LINES ARE TREATED SPECIALLY:                                 */
/*     IF THE EXPTL.FREQ. IS THE SAME TO 1 HZ THEN ONLY THE INVERSE ERROR */
/*             AVERAGED FREQUENCIES ARE USED IN THE FIT FOR THE BLEND     */
/*                                                                        */
/**************************************************************************/
/* "trust region" Marquardt fitting is described in John. E. Dennis and   */
/* Robert B. Schnabel, Numerical Methods for Unsconstrained Optimization  */
/* and Non-linear Equations, Prentice-Hall, 1983.                         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "calpgm.h"
#include "lsqfit.h"
#define NDCARD 130
#define PR_DELAY 6    /* seconds delay between informational messages */
/************** CALFIT interfaces ***********************************/
int qnfmt2(int nqn, short *qnum, /*@out@*/ char *aqnum);
int parer(double par, double errx, double dif, 
                 /*@out@*/ char *ptmp);
int linein(FILE * luin, int *nline, int iqnfmt);
int lineix(FILE * lu, int flg, int nline, int nblkpf, int iqnfmt);
/********************************************************************/
static char card[NDCARD];

int main(int argc, char *argv[])
{
#define NFILE 5
#define LBLEN 10
#define FMT_xbgnMW "%5d: %s%14.5f%14.5f %10.5f %10.5f %9.5f"
#define FMT_xblnMW "%14.5f %10.5f %6.4f\n"
#define FMT_xbgnIR "%5d: %s%14.7f%14.7f %10.7f %10.7f %9.7f"
#define FMT_xblnIR "%14.7f %10.7f %6.4f\n"
  static double zero = 0.;
  static double tiny = 1.5e-38;
  static const char *ext[NFILE] = { "par", "lin", "fit", "bak", "var" };
  enum efile {epar, elin, efit, ebak, evar};
  /*@dependent@*/ SXLINE *linebuf;
  /*@owned@*/ bcd_t *idpar;
  /*@owned@*/ double *par, *erp, *oldpar, *erpar, *dpar, *delbgn;
  /*@owned@*/ double *fitbgn, *var, *oldfit, *fit, *teig, *pmix;
  /*@owned@*/ int *iperm;
  /*@owned@*/ char *parlbl;
  FILE *lupar, *lulin, *lufit, *lubak, *luvar;
  double *egy, *egyder;
  double *pvar, *pfit, *pfitb, *pfitd;
  char *fname[NFILE+1], *tlbl, *tlblnxt;
  double dvec[8], fqfac[4], varv[5], adif, cerr, afrq, rerr, xsqbest, cwid;
  double parfac, ex, xsqir, xsqmw, xerr, xfrq, xsqt, scale, avgir, xwid;
  double avgmw, fqfacq, xerrmx, dif, marqp[3], frq, val, xwt, marqlast, parfac0;
  double bigd, bigf, bigdIR, bigfIR, bigdMW, bigfMW;
  size_t nl, nlsq;
  int ifac, iblk, ndfit, ndiag, lblk, iflg, line, icnt, nfit, indx, nfir;
  int lstf, nitr, i, k, iblnd, lblnd, nline, initl, inpcor, n, marqflg;
  int nsize, nxpar, nxfit, lnext, noptn, ibase, nf, nblkpf, itd, ibcd, nfmt;
  int limlin, iqnfmt, maxdm, maxf, nrj, nqn, itr, npar, nsize_p, ndbcd;
  int supblnd, catqn, ndfree, ndfree0;
  short qnum[2*MAXQN];
  char ch, pare[64], aqnum[6*MAXQN+2], namfil[NDCARD];

  bigdIR = 9.9999999;
  bigdMW = 100. * bigdIR;
  bigfIR = 99999.9999999;
  bigfMW = 100. * bigfIR;
  fqfac[0] = 1;
  fqfac[1] = -1;
  marqflg = 0;
  nsize_p = maxmem(&nl);

  /*     open read and write files */

  filget(argc, argv, NFILE, fname, ext);
  filbak(fname[epar], fname[ebak]);
  lubak = fopenq(fname[ebak], "r");
  lufit = fopenq(fname[efit], "w");

  /*     read in run parameters */

  if (fgetstr(card, NDCARD, lubak) <= 0) {
    puts(" Unable to read title of .par file");
    exit(EXIT_FAILURE);
  }
  chtime(card, 82);
  fputs(card, lufit);
  puts(card);
  catqn = MAXCAT;
  dvec[0] = 100.;
  dvec[1] = 32767.;
  dvec[2] = 1.;
  dvec[3] = 0.;
  dvec[4] = 0.;
  dvec[5] = 1e6;
  dvec[6] = 1.;
  dvec[7] = 1.;
  n = fgetstr(card, NDCARD, lubak);
  if (n != 0)
    n = pcard(card, dvec, 8, NULL);
  if (n == 0) {
    puts(" Unable to read second line of .par file");
    exit(EXIT_FAILURE);
  }
  npar = (int) dvec[0];
  if (dvec[1] < 0.) {
    catqn = MAXQN; dvec[1] = -dvec[1];
  }
  nl = (size_t) dvec[1];
  nitr = (int) dvec[2];
  nxpar = (int) dvec[3];
  marqp[0] = dvec[4];
  xerrmx = dvec[5];
  parfac = dvec[6];
  parfac0 = parfac;
  fqfacq = dvec[7];
  fqfac[2] = fqfacq / 29979.2458;
  fqfac[3] = -fqfac[2];
  marqp[1] = -1;
  if (marqp[0] < 0.)
    marqp[0] = 0.;
  limlin = nline = (int) nl;
  if (fabs(parfac - 1.) > 1e-10) {
    fprintf(lufit, "PARAMETER ERRORS SCALED BY %15.6f", fabs(parfac));
    if (parfac < 0.) fputs(" times the standard error",lufit);
    fputc('\n', lufit);
  }
  if ((long) nline != (long) nl) {
    puts(" number of lines too big for this computer ");
    exit(EXIT_FAILURE);
  }
  /*  read in option card(s) */
  iqnfmt = 0; nfmt = catqn;
  noptn = setopt(lubak, &nfmt, &itd, &ndbcd, namfil);
  if (noptn <= 0) {
    puts("Error reading option lines");
    exit(EXIT_FAILURE);
  }
  nqn = setfmt(&iqnfmt, 1);
  if (nqn > MAXCAT) catqn = nqn;
  nqn = nqn + nqn;
  printf(       "LINES REQUESTED=%5d NUMBER OF PARAMETERS=%3d", limlin, npar);
  fprintf(lufit,"LINES REQUESTED=%5d NUMBER OF PARAMETERS=%3d", limlin, npar);
  printf(       " NUMBER OF ITERATIONS=%3d\n  MARQUARDT PARAMETER =", nitr);
  fprintf(lufit," NUMBER OF ITERATIONS=%3d\n  MARQUARDT PARAMETER =", nitr);
  printf(       "%11.4E max (OBS-CALC)/ERROR =%10.4E\n", marqp[0], xerrmx);
  fprintf(lufit,"%11.4E max (OBS-CALC)/ERROR =%10.4E\n", marqp[0], xerrmx);
  fflush(stdout);
  nl = (size_t) (LBLEN * npar + 1);
  parlbl = (char *) mallocq(nl);
  nl = (size_t) npar * sizeof(double);
  par = (double *) mallocq(nl);
  erp = (double *) mallocq(nl);
  nl = (size_t) (npar * ndbcd + ndbcd + 3);
  idpar = (bcd_t *) mallocq(nl);
  idpar[0] = (bcd_t) ndbcd;
  /* read in parameters  ( and return true if variance ) */
  inpcor =
      getpar(lubak, lufit, &nfit, &npar, idpar, par, erp, parlbl, LBLEN);
  nl = (size_t) npar * sizeof(double);
  oldpar = (double *) mallocq(nl);
  erpar = (double *) mallocq(nl);
  ndfit = nfit + 1;
  ndiag = ndfit + 1;
  if (ndfit > nsize_p) {
    printf("number of independent parameters is too big: %d %d\n",
           nfit, nsize_p - 1);
    exit(EXIT_FAILURE);
  }
  nl = (size_t) nfit * sizeof(int);
  iperm = (int *) mallocq(nl);
  nl = (size_t) ndfit * sizeof(double);
  dpar = (double *) mallocq(nl);
  delbgn = (double *) mallocq(nl);
  if ((nfit & 1) == 0) {
    nlsq = (size_t) ((unsigned) nfit >> 1) * sizeof(double);
    nlsq *= (size_t) ndfit;
  } else {
    nlsq = (size_t) ((unsigned) ndfit >> 1) * sizeof(double);
    nlsq *= (size_t) nfit;
  }
  fitbgn = (double *) mallocq(nlsq);
  var = (double *) mallocq(nlsq);
  nl = (size_t) nfit * sizeof(double);
  nlsq += nl;
  oldfit = (double *) mallocq(nlsq);
  nl *= (size_t) ndfit;
  fit = (double *) mallocq(nl);
  inpcor = getvar(lubak, nfit, var, idpar, erp, inpcor);
  /*  fixup supplied values */
  if (nxpar > npar || nxpar < 0)
    nxpar = 0;
  if (nxpar > 0) {
    fprintf(lufit, "NUMBER OF PARAMETERS EXCLUDED IF INDEX < 0 =%5d\n",
            nxpar);
  }
  nxpar = npar - nxpar;
  nxfit = nfit;
  for (i = npar - 1; i >= nxpar; --i) {
    ibcd = i * ndbcd;
    if (NEGBCD(idpar[ibcd]) == 0)
      --nxfit;
  }
  if (xerrmx < tiny)
    xerrmx = 1e6;
  if (fabs(fqfacq - 1.) >= 1e-10) {
    fprintf(lufit, " IR Frequencies Scaled by %12.10f\n", fqfacq);
  }
  ndfree0 = 0;
  if (inpcor != 0) {            /* initialize fit from supplied variance */
    pvar = var;
    pfit = fit;
    for (n = 1; n <= nfit; ++n) {      /* copy from var */
      dcopy(n, pvar, 1, pfit, ndfit);
      pvar += n; ++pfit;
    }
    pvar = NULL;
    fit[0] = var[0]; dpar[0] = 1.;
    pfit = fit;
    for (i = 0; i < nfit; ++i) { /* scale */
      if (i != 0) {
        pfit += ndfit;
        dcopy(i, &zero, 0, pfit, 1);
      }
      n = nfit - i;
      pfitd = pfit + i;
      val = dnrm2(n, pfitd, 1);
      if (val < tiny) {
        dpar[i] = 1.; scale = 0.;
      } else {
        dpar[i] = scale = 1. / val;
      }
      dscal(n, scale, pfitd, 1);
    }
    pfitd = NULL;
    pfit = fit;
    for (i = 0; i < nfit; ++i) {  
      n = i + 1;
      erpar[i] = ddot(n, pfit, ndfit, pfit, ndfit);
      ++pfit;
    }
    /*  invert lower triangle */
    n = dqrfac(fit, ndfit, nfit, nfit, 0, erpar, iperm);
    if (nfit > n) {
      fputs("supplied variance is singular\n", lufit);
    }
    dqrsolv(fit, ndfit, n, nfit, 0, iperm);
    pfit = fit;
    for (n = 1; i <= nfit; ++n) {       /* unscale */
      scale = dpar[n - 1];
      dscal(n, scale, pfit, ndfit);
      ++pfit;
    }
    fitbgn[0] = fit[0];
    pfit = fit;
    pfitb = fitbgn;
    for (i = 0; i < nfit; ++i) {        /* store in fitbgn */
      if (i != 0)
        pfit += ndiag;
      n = nfit - i;
      dcopy(n, pfit, ndfit, pfitb, 1);
      pfitb += n;
    }
    pfit = pfitb = NULL;
  } else {   /* initialize fit when there is no supplied variance */
    pvar = var;
    pfitb = fitbgn;
    fitbgn[0] = 1. / var[0];
    var[0] = 0.;
    for (n = 1; n < nfit; ++n) {
      ++pfitb;
      dcopy(n, &zero, 0, pfitb, 1);
      pvar += (n + 1);
      pfitb += n;
      *pfitb = 1. / (*pvar);
      if (*pfitb < 1.e-15) 
        --ndfree0;
      *pvar = 0.;
    }
    pfitb = pvar = NULL;
  }
  oldfit[0] = fitbgn[0];
  oldpar[0] = par[0];
  dcopy(nfit, &zero, 0, delbgn, 1);
  linebuf = lbufof(-nfit, nline);
  /* read lines */
  lulin = fopenq(fname[elin], "r");
  maxf = linein(lulin, &nline, iqnfmt);
  if (limlin > nline)
    limlin = nline;
  fclose(lulin);
  fflush(stdout);
  if (nline <= 0) {
    puts("no lines read");
    exit(EXIT_FAILURE);
  }  
  /* initialize block structure */
  nblkpf = maxf;
  maxdm = nsize_p;
  k = setblk(lufit, npar, idpar, par, &nblkpf, &maxdm);
  /* set parameter labels */
  getlbl(npar, idpar, parlbl, namfil, k, LBLEN);
  /* convert lines */
  line = lineix(lufit, nitr, nline, nblkpf, iqnfmt);
  if (line > 0)
    printf("%d bad lines\n", line);
  fflush(stdout);
  /* allocate memory for energy, derivatives, and eigenvector */
  if (maxdm > nsize_p) {
    printf(" Hamiltonian dinension is too big: %d %d\n", maxdm, nsize_p);
    exit(EXIT_FAILURE);
  }
  nlsq = (size_t) maxdm;
  nl = nlsq * sizeof(double);
  nlsq *= nl;
  teig = (double *) mallocq(nlsq);
  nl *= (size_t) (nfit + 2);
  pmix = (double *) mallocq(nl);
  pmix[0] = 1.;
  egy = pmix + maxdm;
  egyder = egy + maxdm;
  rqexit(1);
  xsqbest = zero;
  marqlast = marqp[0];
  /**********************************************************************/
  /* START ITERATION */
  itr = 0;
  dpar[0] = 0.; 
  if (nitr < 0)
    nitr = -nitr;
  nsize = 0;
  /* find energies and derivatives */
  do {
    k = 0; ibcd = nxpar * ndbcd;
    for (i = nxpar; i < npar; ++i, ibcd += ndbcd) {
      if (NEGBCD(idpar[ibcd]) == 0)
        dpar[k++] = par[i];
    }
    lstf = lblk = lnext = 0;
    getdbk(&lnext, &iblk, &indx, &initl, &ifac);
    do {
      line = getdbk(&lnext, &iblk, &indx, &initl, &ifac);
      if (iblk != lblk) {       /*  get size of block */
        if (rqexit(0) != 0)
          break;                /*  check operator interrupt */
        getqn(iblk, 0, 0, qnum, &nsize);
        if (nsize == 0)
          continue;
        lblk = iblk;
        if (nsize > maxdm) {
          printf
              ("WARNING .. SIZE OF BLOCK %d2 IS %d AND EXCEEDS DIMENSIONS\n",
               iblk, nsize);
          exit(EXIT_FAILURE);
        }
        k = (iblk - 1) / nblkpf;
        if (lstf != k && caldelay(PR_DELAY) != 0) {
          printf("Starting Quantum %3d\n", k);
          fflush(stdout);
          lstf = k;
        }
        /*  get energies and derivatives */
        egy = pmix + nsize;
        egyder = egy + nsize;
        hamx(iblk, nsize, npar, idpar, par, egy, teig, egyder, pmix,
             FALSE);
      }
      /* save energies and derivatives for lines in this block */
      dnuadd(nfit, nxfit, initl, indx, ifac, egy, egyder, nsize, line,
             dpar, fqfac);
    } while (lnext != 0);       /* repeat until no more energies */
    if (lnext != 0)
      break;
    k = (lblk - 1) / nblkpf;
    printf("Finished Quantum %3d\n", k);
    fflush(stdout);
    /*   initialize least squares matrix */
    xsqir = xsqmw = avgmw = avgir = 0.;
    pfitb = fitbgn;
    pfit = fit;
    pfitd = fit + nfit;
    dcopy(nfit, &zero, 0, pfitd, ndfit);
    for (n = 1; n <= nfit; ++n) {
      dcopy(n, pfitb, 1, pfit, ndfit);
      val = delbgn[n - 1];
      daxpy(n, val, pfitb, 1, pfitd, ndfit);
      pfitb += n;
      ++pfit;
    }
    xsqt = ddot(nfit, pfitd, ndfit, pfitd, ndfit);
    pfitb = pfitd = pfit = NULL;
    nf = nrj = nfir = 0; 
    icnt = -1;
    for (i = 0; i < 40; ++i)
      fputc(' ', lufit);
    fputs("EXP.FREQ.  -  CALC.FREQ. -   DIFF.  - EXP.ERR.- ", lufit);
    fputs("EST.ERR.-AVG. CALC.FREQ. -  DIFF. - WT.\n", lufit);
    line = 1; ex = 1.;
    do {                        /*   form least squares matrix, LOOP over lines */
      if (icnt <= 0 && caldelay(PR_DELAY) != 0) {
        printf("Fitting Line %d\n", line);
        fflush(stdout);
        icnt += 50;
      }
      lblnd = line;
      scale = rerr = 0;
      iflg = 0; supblnd = 0; xwid = 1.;
      do {                      /*    UNTIL all elements of blend are found */
        i = frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum);
        if ((iblnd & 1) != 0) {
          supblnd = 1; xwid = xerr;
          val = 1. / (scale * xwid);
          val = sqrt(1. - val * val); scale *= val;
          dscal(nfit + 1, val, dpar, 1);
        } else if (i != 0) {
          ex = fabs(xwt / xerr);
          scale += ex;
          /*  accumulate line contributions */
          val = dnuget(iflg, nfit, ex, lblnd, dpar);
          rerr = dpar[nfit];
          ++iflg;
        }
        ++lblnd;
        --icnt;
      } while (iblnd < 0);
      if (iflg == 0) {
        ++line;
        continue;
      }
      /* calculate errors */
      cerr = calerr(nfit, var, dpar) / scale;
      adif = dpar[nfit] / scale;
      afrq = xfrq - adif;
      if (fabs(rerr) < xerrmx) {
        xsqt += rerr * rerr;
        if (xerr < 0.) {
          avgir += adif;
          xsqir += adif * adif;
          ++nfir;
        } else {
          avgmw += adif;
          xsqmw += adif * adif;
        }
        /*  rotate line into FIT matrix */
        jelim(fit, dpar, ndfit, nfit, 1);
        ++nf;
      } else {
        fputs(" ***** NEXT LINE NOT USED IN FIT\n", lufit);
        ++nrj; supblnd = 0;
      }
      if (xerr < 0.) {
        bigf = bigfIR;
        bigd = bigdIR;
      } else {
        bigf = bigfMW;
        bigd = bigdMW;
      }
      if (fabs(afrq) > bigf) {
        afrq = (afrq > 0.) ? bigf : -bigf;
      }
      if (fabs(adif) > bigd) {
        adif = (adif > 0.) ? bigd : -bigd;
      }
      if (iflg == 1) {
        qnfmt2(nqn, qnum, aqnum);
        if (xerr < 0.) {
          fprintf(lufit, FMT_xbgnIR, line, aqnum, xfrq, afrq, adif, xerr,
                  cerr);
        } else {
          fprintf(lufit, FMT_xbgnMW, line, aqnum, xfrq, afrq, adif, xerr,
                  cerr);
        }
        fputc('\n', lufit);
        lblnd = line + 1;
      } else {
        lblnd = line; cwid = 0.;
        do {                    /* UNTIL all elements of blend are printed */
          i = frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum);
          if ((iblnd & 1) != 0) {
            xfrq = 0.; frq = sqrt(cwid);
            if (supblnd != 0) 
              xsqt += cwid / (xwid * xwid);
            qnfmt2(0, qnum, aqnum); 
            if (xerr < 0.) {
              fprintf(lufit, FMT_xbgnIR, lblnd, aqnum, xfrq, frq, frq, xwid,
                      xfrq);
            } else {
              fprintf(lufit, FMT_xbgnMW, lblnd, aqnum, xfrq, frq, frq, xwid,
                      xfrq);
            }
            fputc('\n', lufit);
          } else if (i != 0) {
            if (supblnd != 0) {
              ex = sqrt(xwt) / xwid;
              frq = dnuget(0, nfit, ex, lblnd, dpar);
              jelim(fit, dpar, ndfit, nfit, 1);
              ++nf;
            } else {
              frq = dnuget(-1, nfit, ex, lblnd, dpar);
            }
            dif = xfrq - frq;
            if (fabs(frq) > bigf)
              frq = (frq > 0.) ? bigf : -bigf;
            if (fabs(dif) > bigd)
              dif = (dif > 0.) ? bigd : -bigd;
            qnfmt2(nqn, qnum, aqnum);
            if (xerr < 0.) {
              fprintf(lufit, FMT_xbgnIR, lblnd, aqnum, xfrq, frq, dif, xerr,
                      cerr);
              fprintf(lufit, FMT_xblnIR, afrq, adif, xwt);
            } else {
              fprintf(lufit, FMT_xbgnMW, lblnd, aqnum, xfrq, frq, dif, xerr,
                      cerr);
              fprintf(lufit, FMT_xblnMW, afrq, adif, xwt);
            }
            dif = frq - afrq;
            cwid += xwt * dif * dif;
          }
          ++lblnd;
        } while (iblnd < 0);
      }
      line = lblnd;
    } while (line <= nline);    /* end loop over lines */
    if (nrj > 0) {
      printf(       "%5d Lines rejected from fit\n", nrj);
      fprintf(lufit,"%5d Lines rejected from fit\n", nrj);
    }
    if (nf < 1)
      nf = 1;
    /* zero upper triangle of fit matrix */
    pfit = fit;
    for (k = 1; k < nfit; ++k) {
      pfit += ndfit;
      dcopy(k, &zero, 0, pfit, 1);
    }
    varv[0] = xsqt + nrj * xerrmx * xerrmx; 
    marqlast = marqp[0];
    marqflg = lsqfit(fit, ndfit, nfit, 1, marqp, varv, 
                     oldfit, erpar, dpar, iperm);
    if (marqflg != 0) {
      dcopy(npar, oldpar, 1, par, 1);
      strcpy(card, "Fit Diverging: restore parameters\n");
      fputs(card, lufit);
      fputs(card, stdout);
    } else {
      xsqbest = (xsqt - varv[3]) / nf;
      if (xsqbest > 0.)
        xsqbest = sqrt(xsqbest);
      /* print normalized diagonals */
      fputs("NORMALIZED DIAGONAL:\n", lufit);
      icnt = 0;
      for (k = 0; k < nfit; ++k) {
        fprintf(lufit, "%5d %13.5E", k + 1, erpar[k]);
        if ((++icnt) == 6) {
          icnt = 0;
          fputc('\n', lufit);
        }
      }
      if (icnt > 0)
        fputc('\n', lufit);
    }
    printf(       "MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n",
              marqp[0], marqp[2]);
    fprintf(lufit,"MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n",
              marqp[0], marqp[2]);
    if (parfac0 < 0. && xsqt >= 0.) {
      ndfree = nf + ndfree0;
      if (ndfree <= 0) 
        ndfree = 1;
      parfac = -parfac0 * sqrt(xsqt / ndfree);
      fputs("WARNING: parameter errors multiplied by ",lufit);
      fprintf(lufit, "%15.5f for %8d degrees of freedom\n",
              parfac, ndfree);
    }
    xsqt = xsqt / nf;
    if (xsqt > 0.)
      xsqt = sqrt(xsqt);
    if (xsqbest > xsqt) 
      xsqbest = xsqt;
    /*   get estimated errors  and print parameters */
    for (i = 0; i < 32; ++i)
      fputc(' ', lufit);
    fputs("NEW PARAMETER (EST. ERROR) -- CHANGE THIS ITERATION\n", lufit);
    pfitd = fit;
    tlblnxt = parlbl;
    k = ibase = 0;
    for (i = 0, ibcd = 0; i < npar; ++i, ibcd += ndbcd) {
      tlbl = tlblnxt;
      tlblnxt += LBLEN;
      oldpar[i] = par[i];
      if (NEGBCD(idpar[ibcd]) == 0) {
        if (k != 0)
          pfitd += ndiag;
        n = nfit - k;
        dif = pfitd[n];
        par[i] += dif;
        delbgn[k] -= dif;
        erpar[i] = dpar[k] = parfac * dnrm2(n, pfitd, 1);
        dscal(n, parfac, pfitd, 1);
        parer(par[i], erpar[i], dif, pare);
        putbcd(card, NDCARD, &idpar[ibcd]);
        ch = (*tlblnxt); *tlblnxt = '\0';
        fprintf(lufit, "%4d %s %10.10s %s\n", ++k, card, tlbl,
                pare);
        *tlblnxt = ch;
        ibase = i;
      }
    }
    pfitd = NULL;
    pfit = fit;
    pvar = var;
    for (n = 1; n <= nfit; ++n) {       /* copy to var */
      dcopy(n, pfit, ndfit, pvar, 1);
      ++pfit;
      pvar += n;
    }
    pfit = NULL;
    if (nf > nfir) {
      scale = 1. / (double) (nf - nfir);
      avgmw *= scale;
      xsqmw = sqrt(xsqmw * scale);
    }
    if (nfir > 0) {
      scale = 1. / (double) nfir;
      avgir *= scale;
      xsqir = sqrt(xsqir * scale);
    }
    printf(" MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n",
                  avgmw, avgir);
    fprintf(lufit," MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n",
                  avgmw, avgir);
    printf(       " MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n",
                  xsqmw, xsqir);
    fprintf(lufit," MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n",
                  xsqmw, xsqir);
    ++itr;
    printf(       " END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n",
                  itr, xsqt, xsqbest);
    fprintf(lufit," END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n",
                  itr, xsqt, xsqbest);
    fflush(stdout);
  } while (itr < nitr && 0.999999 * xsqt > xsqbest);

  /*  end of iteration */

  /************************************************************************/

  linebuf = lbufof(-1, 0);      /* release storage */
  setblk(lufit, 0, idpar, par, &nblkpf, &maxdm);        /* release storage */
  if (itr == 0) {
    puts(" output files not updated");
    exit(0);
  }
  rewind(lubak);
  fgetstr(card, NDCARD, lubak);
  chtime(card, 82);
  fputs(card, stdout);
  puts("FIT COMPLETE");
  dcopy(npar, oldpar, 1, par, 1);

  /*  compute correlation matrix */

  prcorr(lufit, nfit, fit, ndfit, dpar);
  fputs(card, lufit);
  fclose(lufit);
  /* save results in output files */
  lupar = fopenq(fname[epar], "w");
  luvar = fopenq(fname[evar], "w");
  fputs(card, lupar);
  fputs(card, luvar);
  k = fgetstr(card, NDCARD, lubak); /* ignore second line and replace */
  if (catqn > MAXCAT) limlin = -limlin;
  nxpar = npar - nxpar;
  fprintf(lupar,"%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
           npar, limlin, nitr, nxpar, marqlast, xerrmx, parfac0, fqfacq);
  fprintf(luvar,"%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
           npar, limlin, nitr, nxpar, marqlast, xerrmx, parfac0, fqfacq);
  for (icnt = 0; icnt < noptn; ++icnt) {        /*  save option lines */
    k = fgetstr(card, NDCARD, lubak);
    card[k] = '\n';
    card[k + 1] = '\0';
    fputs(card, lupar);
    fputs(card, luvar);
  }
  card[LBLEN + 1] = '\0';
  tlbl = parlbl; ibase = 0;
  for (i = 0, ibcd = 0; i < npar; ++i, ibcd += ndbcd) {
    if (k > 0)
      k = fgetstr(card, NDCARD, lubak);     /*  skip parameter lines */
    if (tlbl[0] == '\0') {
      card[0] = '\0';
    } else {
      card[0] = '/';
      memcpy(card + 1, tlbl, LBLEN);
      card[LBLEN + 1] = '\0';
    }
    if (NEGBCD(idpar[ibcd]) == 0) {
      ibase = i;
    } else {
      scale = fabs(par[i]);
      erpar [i] = scale * erpar[ibase];
      par[i] *= par[ibase];
    }
    putbcd(pare, 64, &idpar[ibcd]);
    fprintf(lupar, "%s %23.15E %14.8E %s\n", pare, par[i], erp[i], card);
    fprintf(luvar, "%s %23.15E %14.8E %s\n", pare, par[i], erpar[i], card);
    tlbl += LBLEN;
  }
  while (k > 0) {               /* copy correlation lines or extra parameter lines */
    k = fgetstr(card, NDCARD, lubak);
    if (k > 0) {
      card[k] = '\n';
      card[k + 1] = '\0';
      fputs(card, lupar);
    }
  }
  fclose(lubak);
  fclose(lupar);
  putvar(luvar, nfit, var, dpar);
  fclose(luvar);

  free(pmix);
  free(teig);
  free(fit);
  free(oldfit);
  free(var);
  free(fitbgn);
  free(delbgn);
  free(dpar);
  free(iperm);
  free(erpar);
  free(oldpar);
  free(erp);
  free(idpar);
  free(par);
  free(parlbl);
  return 0;
}                               /* MAIN */

int qnfmt2(int nqn, short *qnum, char *aqnum)
{
  /* Local variables */
  int i;

  /*  formats quantum numbers for output */
  /*     NQN   = number of quanta */
  /*     QNUM  = vector of quanta */
  /*     AQNUM = string of quanta in character form */
  for (i = 0; i < nqn; ++i) {
    sprintf(aqnum, "%3d", (int) qnum[i]);
    aqnum += 3;
  }
  for (i = nqn; i < 12; ++i) {
    aqnum[2] = aqnum[1] = aqnum[0] = ' ';
    aqnum += 3;
  }
  aqnum[0] = '\0';
  return 0;
}                               /* qnfmt2 */

int parer(par, errx, dif, ptmp)
double par, errx, dif;
char *ptmp;
{
  static int czero = (int) '0';
  char *pfmt;
  double adif, apar, aten, aerr;
  char chexp[6], fmt[34];
  int msd, id, ie, efield, ip, lsd, k;


  /*      sets up special format for parameters and errors */
  /*     PAR  = parameter value */
  /*     ERRX = parameter error */
  /*     DIF  = parameter change */
  /*     PTMP = output string for printing */

  apar = par;
  aerr = errx;
  adif = dif;
  efield = 0;
  aten = 1.;
  /*     compute exponent fields */
  ie = (int) (log10(fabs(aerr) + 1.e-37) - 102.5) + 100;
  id = (int) (log10(fabs(adif) + 1.e-37) - 100.0) + 100;
  ip = (int) (log10(fabs(apar) + 1.e-37) - 100.0) + 100;
  lsd = -ie;
  if (lsd < 0)
    lsd = 0;
  msd = (ip > id) ? ip : id;
  /*  check for too many digits */
  k = 14 - ip;
  if (k < lsd)
    lsd = k;
  k = 10 - id;
  if (k < lsd)
    lsd = k;
  if (msd <= -2) {              /* number too small without exponent */
    k = (1 - msd) / 3;
    efield = -3 * k;
    while ((--k) >= 0)
      aten *= 1000;
  } else if (lsd < 0) {         /* number too big without exponent */
    k = (1 + msd) / 3;
    if (k > 0)
      efield = 3 * k;
    while ((--k) >= 0)
      aten *= 0.001;
  }
  if (efield != 0) {            /* E format */
    lsd += efield;
    memcpy(chexp, "0fE+00", 6);
    if (efield < 0) {
      chexp[3] = '-';
      efield = -efield;
    }
    msd = efield / 10;
    if (msd > 0) {
      efield -= msd * 10;
      chexp[4] = (char) (msd + czero);
    }
    chexp[5] = (char) (efield + czero);
    apar *= aten;
    aerr *= aten;
    adif *= aten;
  } else {                      /* F format */
    memcpy(chexp, "0f    ", 6);
  }
  if (lsd > 9)
    lsd = 9;
  if (lsd > 0)
    chexp[0] = (char) (lsd + czero);
  while ((lsd--) > 0)
    aerr *= 10.;
  ie = (int) (aerr + 0.5);
  pfmt = fmt;
  memcpy(pfmt, "%16.", 4);
  pfmt += 4;
  memcpy(pfmt, chexp, 2);
  pfmt += 2;
  memcpy(pfmt, "(%3d)", 5);
  pfmt += 5;
  memcpy(pfmt, chexp + 2, 4);
  pfmt += 4;
  memcpy(pfmt, " %12.", 5);
  pfmt += 5;
  memcpy(pfmt, chexp, 6);
  pfmt += 6;
  *pfmt = '\0';
  sprintf(ptmp, fmt, apar, ie, adif);
  return 0;
}                               /* parer */


int linein(luin, nline, iqnfmt)
FILE *luin;
int *nline;
int iqnfmt;
{
  /* Local variables */
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xfrqx, xerrx;
  int nqn, nqnu, nqnl, kqnu, kqnl, i, iqf, ipace, mxline, mxqn, isblnd, icmp;
  short nbln, nqnt[20], *iqnum;

  /*   get lines from input  and stores them */

  /*     LUIN= unit for finding lines */
  /*     NLINE = number of lines */
  /*     IQNFMT= qunatum number format for line input */
  /*     RETURN: largest quantum number */
  /*******************************************************************/

  mxline = *nline;
  mxqn = 1;
  nbln = 1;
  nqn = deflin(iqnfmt, nqnt);
  nqnu = nqn - 1;
  if (nqnt[nqnu] < 0)
    nqnu = 0;
  kqnu = nqnt[nqnu];
  nqnl = nqnu + nqn;
  kqnl = nqnt[nqnl];
  ipace = 100;
  xfrqx = xerrx = 0.;
  icmp = 0;
  for (i = 1; i <= mxline; ++i) {       /*  loop for reading lines */
    xline = lbufof(1, i);
    iqnum = xline->qn;
    if (getlin(luin, nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn, 
               card, NDCARD) < 0) {
      *nline = i - 1;
      return mxqn;
    }
    iqf = iqnum[nqnu];
    if (iqf == -1) {
      if (kqnu >= 0) {
        iqf = -iqnum[kqnu];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnu] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (iqf > mxqn)
      mxqn = iqf;
    iqf = iqnum[nqnl];
    if (iqf == -1) {
      if (kqnl > 0) {
        iqf = -iqnum[kqnl];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnl] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (mxqn < iqf)
      mxqn = iqf;
    xline->xfrq = xfrqn;
    xline->xerr = (float) xerrn;
    xline->xwt = (float) fabs(xwtn);
    xline->linku = 0;
    xline->linkl = 0;
    isblnd = 0;
    if (icmp != 0 && fabs(xfrqn - xfrqx) < fabs(xfrqn) * 1.e-14 + 1.e-8) { 
      /* frq match */
      if (fabs(xerrn - xerrx) < 1e-7) {
        isblnd = 1;
      } else if ((xerrn / xerrx) > 2.0 && nbln > 2) {
        isblnd = 1; ++nbln; icmp = 0;
        xline->xwt = (float)0.;
        iqnum[0] = (short)-1;
        iqnum[nqn] = iqnum[0];
      }
    }
    if (isblnd != 0) {
      xline->bln = nbln;
      xline = lbufof(1, i - 1);
      xline->bln = -2;
      nbln += 2;
    } else {
      xline->bln = 0;
      nbln = 2; icmp = 1;
    }
    if (ipace <= i) {
      ipace += 100;
      printf("Reading Line %d\n", i);
      fflush(stdout);
    }
    xerrx = xerrn;
    xfrqx = xfrqn;
  }
  return mxqn;
}                               /* linein */

int lineix(lu, flg, nline, nblkpf, iqnfmt)
FILE *lu;
int flg, nline, nblkpf, iqnfmt;
{  /*   get lines from input and store them */
  /*     LU = unit for printout of lines ( if > 0 ) */
  /*     NLINE = number of lines */
  /*     NBLKPF= number of blocks per F */
  /*     IQNFMT= qunatum number format for line input */
  /******************************************************************/
  static int nsort = 2048;
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xnorm;
  int nblk, ipos, i, j, ipace, nread, iblkl, iblku, ncat;
  int linkx, indxl, linky, indxu, orgblk, nqn, nqn2, nbad;
  /*@owned@*/ int *prvblk;
  short *iqnum;
  char aqnum[6*MAXQN+2];

  nbad = 0;
  nblk = 0;
  prvblk = (int *) mallocq((size_t) (nsort + 1) * sizeof(int));
  prvblk[0] = 0;
  for (i = 1; i <= nsort; ++i) {
    prvblk[i] = 0;
  }
  nqn = iqnfmt % 10;
  if (nqn == 0) nqn = 10;
  nqn2 = nqn + nqn; ncat = nqn2;
  if (ncat < 12) ncat = 12;
  i = (iqnfmt / 100) % 5;
  if (i >= nqn) {
    ipos = 1;
  } else {
    ipos = nqn;
  }
  if (flg < 0) {
    fputs(" LINE,BLKU,INDXU,BLKL,INDXL,QUANTUM NUMBERS", lu);
    for (i = 0; i < 19; ++i)
      fputc(' ', lu);
    fputs("ENERGY    EXP. ERROR    WEIGHTS\n", lu);
  }
  xnorm = 0.;
  ipace = 50;
  /*       loop for converting lines */
  for (nread = 1; nread <= nline; ++nread) {
    xline = lbufof(1, nread);
    xfrqn = xline->xfrq;
    xerrn = xline->xerr;
    xwtn = xline->xwt;
    /* find blocks and index for upper and lower states */
    iqnum = xline->qn;
    getblk(&iblku, &indxu, iqnum, nblkpf, ipos, nqn);
    getblk(&iblkl, &indxl, &iqnum[nqn], nblkpf, ipos, nqn);
    if (iblkl == 0 && iqnum[nqn] >= 0)
      iblku = 0;
    xline->ibu = iblku;
    xline->inu = (short) indxu;
    xline->ibl = iblkl;
    xline->inl = (short) indxl;
    if (iblku == 0 && (xline->bln & 1) == 0) {
      /*  print out bad line and try for next */
      ++nbad;
      xline->xwt = 0.;
      xwtn = 0.;
      qnfmt2(nqn2, iqnum, aqnum);
      printf(    "Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
      fprintf(lu,"Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
    } else {
      /*  set up links for calculating in order of block */
      if (iblku <= iblkl) {
        lnlink(prvblk, nsort, iblku, nread);
        lnlink(prvblk, nsort, iblkl, -nread);
        if (nblk < iblkl)
          nblk = iblkl;
      } else {
        lnlink(prvblk, nsort, iblkl, -nread);
        lnlink(prvblk, nsort, iblku, nread);
        if (nblk < iblku)
          nblk = iblku;
      }
      if (flg < 0) {
        iqnum = xline->qn;
        fprintf(lu," %4d%4d%4d%4d%4d:", nread, iblku, indxu, iblkl, indxl);
        for (i = 0; i < ncat; ++i) {
          j = iqnum[i];
          fprintf(lu, "%3d", j);
        }
        fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
        j = xline->bln;
        if (j != 0) {
          fprintf(lu, "   Line Blended with %3d\n", nread - (j >> 1));
        } else {
          fputc('\n', lu);
        }
      }
    }
    /* let the user know something is happening */
    if (ipace <= nread || nread == nline) {
      ipace += 50;
      printf("Converting Line %d\n", nread);
      fflush(stdout);
    }
    j = xline->bln;
    if (j != 0) {
      xnorm += xwtn;      
      if (j > 0) {     /* normalize weights */
        xnorm = 1. / xnorm; 
        for (j = nread - (j >> 1); j <= nread; ++j) {
          xline = lbufof(1, j);
          xline->xwt = (float) (xline->xwt * xnorm);
        }
        xnorm = 0.;
      }
    } else {
      xline->xwt = 1.;
    }
  }                             /* end loop for converting lines */
  orgblk = nsort;
  while (nblk > orgblk) {       /* finish up links */
    --orgblk;
    linkx = 0;
    for (j = 0; j < nsort; ++j) {
      if (prvblk[j] != 0) {
        linkx = prvblk[j];
        prvblk[j] = 0;
      }
    }
    prvblk[0] = linkx;
    prvblk[nsort] = 0;
    linkx = lnlink(prvblk, nsort, 1, 0);
    while (linkx != 0) {
      linky = linkx;
      getdbk(&linkx, &iblkl, &j, &j, &j);
      iblkl -= orgblk;
      lnlink(prvblk, nsort, iblkl, linky);
    }
    orgblk += nsort;
  }
  free(prvblk);
  return nbad;
}                               /* lineix */
