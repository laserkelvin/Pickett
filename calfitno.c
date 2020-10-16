/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: read option cards with fgetstr */
/*   30 Dec.  1999: include changes for dlsq */

/**************************************************************************/
/*      HMP, change code for computing constraints                        */
/*      HMP, change PRCOR                                                 */
/*      HMP, change for generating binary file of variance                */
/*      HMP, constained derivatives calculated within HAMX                */
/*      HMP, IDPAR is now REAL*8                                          */
/*      HMP, fixed bug in LINEIN for negative F values in lower state     */
/*      HMP, changed to Marquart-Levenberg iteration control              */
/*      HMP, fixed bug when bad lines and blends                          */
/*                                                                        */
/*   THIS IS A GENERALIZED LINE FITTING PROGRAM                           */
/*   IT FITS LINES TO PARAMETERS IN A MODEL HAMILTONIAN                   */
/*   BLENDED LINES ARE TREATED SPECIALLY:                                 */
/*     IF THE EXPTL.FREQ. IS THE SAME TO 1 HZ THEN ONLY THE INVERSE ERROR */
/*             AVERAGED FREQUENCIES ARE USED IN THE FIT FOR THE BLEND     */
/**************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "calpgm.h"
int jelim(double *t, double *vec, int ndm, int nz);
int qnfmt(int nqn, short *qnum, char *aqnum);
int parer(double par, double errx, double dif, char *ptmp);
int linein(FILE * luin, int *nline, int iqnfmt);
int lineix(FILE * lu, int flg, int nline, int nblkpf, int iqnfmt);

int main(argc, argv)
int argc;
char *argv[];
{
#define LPAR 0
#define LLIN 1
#define LFIT 2
#define LBAK 3
#define LVAR 4
#define LBIN 5
#define NFILE 6
#define LBLEN 10
  /* Initialized data */

  static double zero = 0.;
  static char *ext[NFILE] = { "par", "lin", "fit", "bak", "var", "unf" };
  static char fmt_xbgn[] = "%5d: %s %13.5f %13.5f %10.5f %9.5f %9.5f";

  /* Local variables */
  FILE *lupar, *lulin, *lufit, *lubak, *luvar, *lubin;
  double *fit, *par, *dpar, *erp, *erpar, *oldpar, *idpar;
  double *var, *fitbgn, *delbgn, *egy, *egyder, *teig, *oldfit;
  double *pvar, *pfit, *pmix, *pfitb, *pfitd;
  int *iperm;
  char *fname[NFILE], *parlbl, *tlbl, *tlblnxt;
  double dvec[8], fqfac[4], adif, bigd, bigf, cerr, afrq, rerr, marqmin;
  double parfac, ex, xsqir, xsqmw, xerr, xfrq, xsqt, scale, avgir, xsqold;
  double avgmw, fqfacq, xerrmx, dif, marqp, frq, xsqtold, sum, xsq, xwt;
  unsigned long nl;
  int ifac, iblk, ndfit, ndiag, lblk, iflg, line, icnt, nfit, indx, nfir;
  int lstf, nitr, i, k, iblnd, lblnd, nline, initl, inpcor, n, marqflg;
  int nsize, nxpar, nxfit, lnext, noptn, ibase, nf, nblkpf, itd, nvdiag;
  int limlin, iqnfmt, maxdm, maxf, nrj, nrjold, nqn, itr, npar, nsize_p;
  short qnum[12];
  char ch, card[82], pare[64], aqnum[38], namfil[82];

  bigf = 9999999.9;
  bigd = 999.99999;
  fqfac[0] = 1;
  fqfac[1] = -1;
  marqmin = 0.;
  nsize_p = maxmem(&nl);

  /*     open read and write files */

  filget(argc, argv, NFILE, fname, ext);
  filbak(fname[LPAR], fname[LBAK]);
  lubak = fopenq(fname[LBAK], "r");
  lufit = fopenq(fname[LFIT], "w");

  /*     read in run parameters */

  if (fgetstr(card, 82, lubak) <= 0) {
    puts(" Unable to read title of .par file");
    exit(1);
  }
  chtime(card, 82);
  fputs(card, lufit);
  puts(card);
  dvec[0] = 100.;
  dvec[1] = 32767.;
  dvec[2] = 1.;
  dvec[3] = 0.;
  dvec[4] = 1e-14;
  dvec[5] = 1e6;
  dvec[6] = 1.;
  dvec[7] = 1.;
  n = fgetstr(card, 82, lubak);
  if (n)
    n = pcard(card, dvec, 8);
  if (n == 0) {
    puts(" Unable to read second line of .par file");
    exit(1);
  }
  npar = (int) dvec[0];
  nl = (long) dvec[1];
  nitr = (int) dvec[2];
  nxpar = (int) dvec[3];
  marqp = dvec[4];
  xerrmx = dvec[5];
  parfac = dvec[6];
  fqfacq = dvec[7];
  fqfac[2] = fqfacq / 29979.2458;
  fqfac[3] = -fqfac[2];
  limlin = nline = (int) nl;
  parfac = fabs(parfac);
  if (fabs(parfac - 1.) > 1e-10) {
    printf("PARAMETER ERRORS SCALED BY %f", parfac);
  }
  if ((long) nline != (long) nl) {
    puts(" number of lines too big for this computer ");
    exit(1);
  }
  /*  read in option card(s) */
  iqnfmt = 0;
  noptn = setopt(lubak, &iqnfmt, &itd, namfil);
  if (noptn <= 0) {
    puts("Error reading option lines");
    exit(1);
  }
  nqn = (iqnfmt % 10) << 1;
  sprintf(card, "LINES REQUESTED=%5d NUMBER OF PARAMETERS=%3d",
	  limlin, npar);
  fputs(card, lufit);
  fputs(card, stdout);
  sprintf(card, " NUMBER OF ITERATIONS=%3d\n  MARQUARDT PARAMETER=", nitr);
  fputs(card, lufit);
  fputs(card, stdout);
  sprintf(card, "%11.4E max (OBS-CALC)/ERROR =%10.4E\n", marqp, xerrmx);
  fputs(card, lufit);
  fputs(card, stdout);
  /* read in parameters  ( and return true if variance ) */
  nl = LBLEN * npar + 1;
  parlbl = (char *) mallocq(nl);
  nl = npar;
  nl *= (unsigned) sizeof(double);
  par = (double *) mallocq(nl);
  idpar = (double *) mallocq(nl);
  erp = (double *) mallocq(nl);
  inpcor =
      getpar(lubak, lufit, &nfit, &npar, idpar, par, erp, parlbl, LBLEN);
  nl = npar;
  nl *= (unsigned) sizeof(double);
  oldpar = (double *) mallocq(nl);
  nl += sizeof(double);
  erpar = (double *) mallocq(nl);
  ndfit = nvdiag = nfit + 1;
  ndiag = ndfit + 1;
  if (ndfit > nsize_p) {
    printf("number of independent parameters is too big: %d %d\n",
	   nfit, nsize_p - 1);
    exit(1);
  }
  nl = nfit * sizeof(int);
  iperm = (int *) mallocq(nl);
  nl = nfit * sizeof(double);
  delbgn = (double *) mallocq(nl);
  nl += sizeof(double);
  dpar = (double *) mallocq(nl);
  nl *= (unsigned) nfit;
  fitbgn = (double *) mallocq(nl);
  var = fitbgn + nfit;
  lubin = NULL;
  inpcor = getvar(lubak, lubin, nfit, var, idpar, erp, inpcor);
  nl = ndfit * sizeof(double);
  nl *= (unsigned) ndfit;
  fit = (double *) mallocq(nl);
  nl = nfit;
  nl = (nfit * (unsigned) ndfit) >> 1;
  nl += nfit;
  nl *= sizeof(double);
  oldfit = (double *) mallocq(nl);
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
    if (idpar[i] >= 0.)
      --nxfit;
  }
  if (xerrmx < 1.5e-38)
    xerrmx = 1e6;
  if (fabs(fqfacq - 1.) >= 1e-10) {
    fprintf(lufit, " IR Frequencies Scaled by %12.10f\n", fqfacq);
  }
  if (inpcor) {			/* initialize fit from supplied variance */
    pvar = var;
    pfit = fit;
    for (i = 0; i < nfit; ++i) {	/* copy from var and scale */
      n = nfit - i;
      pfitd = pfit + i;
      dcopy(i, &zero, 0, pfit, 1);
      dcopy(n, pvar, nfit, pfitd, 1);
      dpar[i] = scale = 1 / (enorm(n, pfitd, 1) + 1.5e-38);
      dscal(n, scale, pfitd, 1);
      pvar += nvdiag;
      pfit += ndfit;
    }
    /*  invert lower triangle */
    if (nfit > dlsq(fit, ndfit, nfit, nfit, 0, erpar, iperm)) {
      puts("supplied variance is singular");
      exit(1);
    }
    pfit = fit;
    pfitb = fitbgn;
    for (n = 1; n <= nfit; ++n) {
      dcopy(n, pfit, ndfit, pfitb, nfit);
      scale = dpar[n - 1];
      dscal(n, scale, pfitb, nfit);
      ++pfit;
      ++pfitb;
    }
  } else {			/* initialize fit when there is no supplied variance */
    dcopy(nfit, &zero, 0, fitbgn, 1);
    pvar = var;
    pfitb = fitbgn;
    for (i = 0; i < nfit; ++i) {
      pvar = pfitb + nfit;
      *pfitb = 1. / (*pvar);
      *pvar = 0.;
      pfitb += nvdiag;
    }
  }
  dcopy(nfit, &zero, 0, delbgn, 1);
  if (lbufof(-nfit, nline) == NULL) {
    puts(" memory allocation error in lbufof");
    exit(1);
  }

  /* read lines */
  lulin = fopenq(fname[LLIN], "r");
  maxf = linein(lulin, &nline, iqnfmt);
  if (nline < limlin)
    limlin = nline;
  fclose(lulin);
  /* initialize block structure */
  nblkpf = maxf;
  maxdm = nsize_p;
  k = setblk(lufit, npar, idpar, par, &nblkpf, &maxdm);
  /* set parameter labels */
  getlbl(npar, idpar, parlbl, namfil, k, LBLEN);
  /* convert lines */
  line = lineix(lufit, nitr, nline, nblkpf, iqnfmt);
  if (line)
    printf("%d bad lines\n", line);
  /* allocate memory for energy, derivatives, and eigenvector */
  nl = maxdm;
  if (maxdm > nsize_p) {
    printf(" Hamiltonian dinension is too big: %d %d\n", maxdm, nsize_p);
    exit(1);
  }
  nl *= (unsigned) (nl * sizeof(double));
  teig = (double *) mallocq(nl);
  nl = maxdm;
  nl *= (unsigned) (ndiag * sizeof(double));
  pmix = (double *) mallocq(nl);
  egy = pmix + maxdm;
  egyder = egy + maxdm;
  rqexit(-1);
  xsqtold = xsqold = zero;
  nrjold = 0;
  /**********************************************************************/
  /* START ITERATION */
  itr = 0;
  if (nitr < 0)
    nitr = -nitr;

  /* find energies and derivatives */
  do {
    k = 0;
    for (i = nxpar; i < npar; ++i) {
      if (idpar[i] >= 0.)
	dpar[k++] = par[i];
    }
    lstf = lblk = lnext = 0;
    getdbk(&lnext, &iblk, &indx, &initl, &ifac);
    do {
      line = getdbk(&lnext, &iblk, &indx, &initl, &ifac);
      if (iblk != lblk) {	/*  get size of block */
	if (rqexit(0))
	  break;		/*  check operator interrupt */
	getqn(iblk, 0, qnum, &nsize);
	if (nsize == 0)
	  continue;
	lblk = iblk;
	if (nsize > maxdm) {
	  printf
	      ("WARNING .. SIZE OF BLOCK %d2 IS %d AND EXCEEDS DIMENSIONS\n",
	       iblk, nsize);
	  exit(1);
	}
	k = (iblk - 1) / nblkpf;
	if (lstf != k) {
	  printf("Starting Quantum %3d\n", k);
	  lstf = k;
	}
	/*  get energies and derivatives */
	egy = pmix + nsize;
	egyder = egy + nsize;
	hamx(iblk, nsize, npar, idpar, par, egy, teig, egyder, pmix, 0);
      }
      /* save energies and derivatives for lines in this block */
      dnuadd(nfit, nxfit, initl, indx, ifac, egy, egyder, nsize, line,
	     dpar, fqfac);
    } while (lnext);		/* repeat until no more energies */
    if (lnext)
      break;
    /*   initialize least squares matrix */
    xsqt = xsqir = xsqmw = avgmw = avgir = 0.;
    pfitb = fitbgn;
    pfitd = fit;
    for (k = 0; k < nfit; ++k) {
      n = nfit - k;
      dcopy(n, pfitb, 1, pfitd, 1);
      sum = ddot(n, &delbgn[k], 1, pfitb, 1);
      pfitd[n] = sum;
      xsqt += sum * sum;
      pfitb += nvdiag;
      pfitd += ndiag;
    }
    nf = nrj = nfir = 0;
    icnt = -1;
    for (i = 0; i < 40; ++i)
      fputc(' ', lufit);
    fputs("EXP.FREQ.  -  CALC.FREQ. -   DIFF.  - EXP.ERR.- ", lufit);
    fputs("EST.ERR.-AVG. CALC.FREQ. -  DIFF. - WT.\n", lufit);
    line = 1;
    do {			/*   form least squares matrix, LOOP over lines */
      if (icnt <= 0) {
	printf("Fitting Line %d\n", line);
	icnt += 50;
      }
      lblnd = line;
      scale = 0.;
      iflg = 1;
      do {			/*    UNTIL all elements of blend are found */
	--icnt;
	if (frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum)) {
	  ex = fabs(xwt / xerr);
	  scale += ex;
	  /*  accumulate line contributions */
	  dnuget(iflg, nfit, ex, lblnd, dpar);
	  iflg = -1;
	  if (iblnd == 0)
	    ++lblnd;
	}
      } while (iblnd == 0);
      /* calculate errors */
      if (iflg > 0) {
	++line;
	continue;
      }
      cerr = calerr(nfit, var, dpar) / scale;
      afrq = dpar[nfit] / scale;
      adif = xfrq - afrq;
      rerr = adif * scale;
      dpar[nfit] = rerr;
      if (fabs(rerr) < xerrmx) {
	++nf;
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
	jelim(fit, dpar, ndfit, nfit);
      } else {
	fputs(" ***** NEXT LINE NOT USED IN FIT\n", lufit);
	++nrj;
      }
      if (fabs(afrq) > bigf) {
	afrq = (afrq > 0.) ? bigf : -bigf;
      }
      if (fabs(adif) > bigd) {
	adif = (adif > 0.) ? bigd : -bigd;
      }
      if (iblnd > 0) {
	qnfmt(nqn, qnum, aqnum);
	fprintf(lufit, fmt_xbgn, line, aqnum, xfrq, afrq, adif, xerr,
		cerr);
	fputc('\n', lufit);
	lblnd = line + 1;
      } else {
	lblnd = line;
	do {			/* UNTIL all elements of blend are printed */
	  if (frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum)) {
	    frq = dnuget(0, nfit, xfrq, lblnd, dpar);
	    dif = xfrq - frq;
	    if (fabs(frq) > bigf)
	      frq = (frq > 0.) ? bigf : -bigf;
	    if (fabs(dif) > bigd)
	      dif = (dif > 0.) ? bigd : -bigd;
	    qnfmt(nqn, qnum, aqnum);
	    fprintf(lufit, fmt_xbgn, lblnd, aqnum, xfrq, frq, dif, xerr,
		    cerr);
	    fprintf(lufit, "%15.5f %9.5f %6.4f\n", afrq, adif, xwt);
	  }
	  ++lblnd;
	} while (iblnd == 0);
      }
      line = lblnd;
    } while (line <= nline);	/* end loop over lines */
    if (nrj > 0) {
      sprintf(card, "%5d Lines rejected from fit\n", nrj);
      fputs(card, lufit);
      fputs(card, stdout);
    }
    if (nf > 1)
      xsqt = xsqt / nf;
    marqflg = 0;
    if (itr > 0) {		/* update Marquardt parameter */
      if (nrj <= nrjold && xsqt <= xsqtold) {
	marqp *= 0.1;
	marqflg = -1;
      } else if (nrj >= nrjold && xsqt > xsqtold * 1.21) {
	puts("Fit Diverging: restore parameters");
	marqp *= 10.;
	if (marqp < marqmin)
	  marqp = marqmin;
	pfitb = oldfit;
	pfitd = fit;
	for (k = 0; k < nfit; ++k) {
	  n = ndfit - k;
	  dcopy(n, pfitb, 1, pfitd, 1);
	  pfitb += n;
	  pfitd += ndiag;
	}
	dcopy(npar, oldpar, 1, par, 1);
	marqflg = 1;
      }
    }
    /* save fit matrix  and zero upper triangle */
    pfit = fit;
    pfitb = oldfit;
    for (k = 0; k < nfit; ++k) {
      n = ndfit - k;
      dcopy(n, pfit + k, 1, pfitb, 1);
      dcopy(k, &zero, 0, pfit, 1);
      pfit += ndfit;
      pfitb += n;
    }
    pfit = fit + nfit;
    xsq = ddot(nfit, pfit, ndfit, pfit, ndfit);
    if (nf > 1)
      xsq = xsq / nf;

    /*  normalize parameters */
    marqmin = 1.;
    pfit = pfitd = fit;
    i = 0;
    for (k = 0; k < nfit; ++k) {
      n = k + 1;
      dpar[k] = scale = 1 / enorm(n, pfit, ndfit);
      dscal(n, scale, pfit, ndfit);
      sum = fabs(*pfitd);
      if (sum < marqmin)
	marqmin = sum;
      ++pfit;
      pfitd += ndiag;
    }
    if (marqp > 0.) {
      ex = sqrt(marqp);
      for (k = nfit - 1; k >= 0; --k) {
	dcopy(ndfit, &zero, 0, erpar, 1);
	erpar[k] = ex;
	jelim(fit, erpar, ndfit, nfit);
      }
    }
    sprintf(card, "MARQUARDT PARAMETER = %g\n", marqp);
    fputs(card, lufit);
    fputs(card, stdout);

    dlsq(fit, ndfit, nfit, nfit, 1, erpar, iperm);	/*  invert lower triangle */

    /* print normalized diagonals */
    fputs("NORMALIZED DIAGONAL:\n", lufit);
    pfit = fit;
    pvar = var;
    icnt = 0;
    for (k = 0; k < nfit; ++k) {
      fprintf(lufit, "%5d %13.5E", k + 1, fabs(erpar[k]));
      if ((++icnt) == 6) {
	icnt = 0;
	fputc('\n', lufit);
      }
    }
    if (icnt)
      fputc('\n', lufit);

    /*   get estimated errors  and print parameters */
    for (i = 0; i < 32; ++i)
      fputc(' ', lufit);
    fputs("NEW PARAMETER (EST. ERROR) -- CHANGE THIS ITERATION\n", lufit);
    pvar = var;
    pfitd = fit;
    tlblnxt = parlbl;
    k = ibase = 0;
    for (i = 0; i < npar; ++i) {
      oldpar[i] = par[i];
      tlbl = tlblnxt;
      tlblnxt += LBLEN;
      if (idpar[i] >= 0.) {
	n = nfit - k;
	dif = pfitd[n] * dpar[k];
	par[i] += dif;
	delbgn[k] -= dif;
	scale = dpar[k] * parfac;
	erpar[i] = fabs(scale) * enorm(n, pfitd, 1);
	parer(par[i], erpar[i], dif, pare);
	ch = (*tlblnxt);
	*tlblnxt = 0;
	fprintf(lufit, "%4d%16.0f %10.10s %s\n", ++k, idpar[i], tlbl,
		pare);
	*tlblnxt = ch;
	dscal(n, scale, pfitd, 1);
	dcopy(n, pfitd, 1, pvar, nfit);
	pvar += nvdiag;
	pfitd += ndiag;
	ibase = i;
      } else {
	scale = oldpar[i] / oldpar[ibase];
	par[i] = par[ibase] * scale;
	erpar[i] = erpar[ibase] * fabs(scale);
      }
    }
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
    sprintf(card, " MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n",
	    avgmw, avgir);
    fputs(card, lufit);
    fputs(card, stdout);
    sprintf(card, " MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n",
	    xsqmw, xsqir);
    fputs(card, lufit);
    fputs(card, stdout);
    if (marqflg > 0) {
      xsq = sqrt(xsqtold - xsqold);
    } else {
      xsqtold = xsqt;
      xsqold = xsq;
      nrjold = nrj;
      xsq = sqrt(xsqt - xsq);
    }
    xsqt = sqrt(xsqt);
    ++itr;
    sprintf(card,
	    " END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n",
	    itr, xsqt, xsq);
    fputs(card, lufit);
    fputs(card, stdout);
  } while (itr < nitr && 0.999999 * xsqt > xsq);

  /*  end of iteration */

  /************************************************************************/

  lbufof(-1, 0);
  if (itr == 0) {
    puts(" output files not updated");
    exit(1);
  }
  rewind(lubak);
  fgetstr(card, 82, lubak);
  chtime(card, 82);
  fputs(card, stdout);
  puts("FIT COMPLETE");

  /*  compute correlation matrix */

  prcorr(lufit, nfit, fit, ndfit, erpar, idpar);
  fputs(card, lufit);
  fclose(lufit);
  /* save results in output files */
  lupar = fopenq(fname[LPAR], "w");
  luvar = fopenq(fname[LVAR], "w");
  fputs(card, lupar);
  fputs(card, luvar);
  fgets(card, 82, lubak);	/* ignore second line and replace */
  nxpar = npar - nxpar;
  sprintf(card, "%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
	  npar, limlin, nitr, nxpar, marqp, xerrmx, parfac, fqfacq);
  fputs(card, lupar);
  fputs(card, luvar);
  for (icnt = 0; icnt < noptn; ++icnt) {	/*  save option lines */
    k = fgetstr(card, 82, lubak);
    card[k] = '\n';
    card[k + 1] = 0;
    fputs(card, lupar);
    fputs(card, luvar);
  }
  card[LBLEN + 1] = 0;
  tlbl = parlbl;
  for (i = 0; i < npar; ++i) {
    if (inpcor > 0)
      inpcor = fgetstr(card, 82, lubak);	/*  skip parameter lines */
    if (tlbl[0] == 0) {
      card[0] = 0;
    } else {
      card[0] = '/';
      memcpy(card + 1, tlbl, LBLEN);
    }
    fprintf(lupar, "%16.0f %23.15E %14.8E %s\n", idpar[i], par[i], erp[i],
	    card);
    fprintf(luvar, "%16.0f %23.15E %14.8E %s\n", idpar[i], par[i],
	    erpar[i], card);
    tlbl += LBLEN;
  }
  while (inpcor > 0) {		/* copy correlation lines */
    inpcor = fgetstr(card, 82, lubak);
    if (inpcor > 0) {
      card[inpcor] = '\n';
      card[inpcor + 1] = 0;
      fputs(card, lupar);
    }
  }
  fclose(lubak);
  fclose(lupar);
  lubin = fopenq(fname[LBIN], "wb");
  putvar(luvar, lubin, nfit, idpar, var, erpar);
  fclose(luvar);
  fclose(lubin);

  return 0;
}				/* MAIN */

int jelim(t, vec, ndm, n)
double *t, *vec;
int ndm, n;
{				/*   do orthogonal transformations to rotate vector VEC into matrix T */
  /*   NDM  = dimensioned column length of T */
  /*   IBGN = first element of VEC */
  /*   N    = length of VEC */

  static double tiny = 1.e-37;
  double c, s, sa, ca, r, *td;
  int nn;

  for (nn = n; nn > 0; --nn, t += ndm) {
    s = (*vec);
    c = (*t);
    td = t;
    ++t;
    ++vec;
    sa = fabs(s);
    ca = fabs(c);
    if (sa > ca) {
      r = c;
      c = s;
      s = r;
      *td = c;
      sa = ca;
      dswap(nn, t, 1, vec, 1);
    }
    if (sa < tiny)
      continue;
    s = s / c;
    r = 1. + s * s;
    if (r > 1.) {
      r = sqrt(r);
      *td = r * c;
      c = 1. / r;
      s *= c;
      drot(nn, t, 1, vec, 1, c, s);
    } else {
      s = -s;
      daxpy(nn, s, t, 1, vec, 1);
    }
  }
  return 0;
}				/* jelim */


int qnfmt(nqn, qnum, aqnum)
int nqn;
short *qnum;
char *aqnum;
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
  aqnum[0] = 0;
  return 0;
}				/* qnfmt */

int parer(par, errx, dif, ptmp)
double par, errx, dif;
char *ptmp;
{
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
  if (msd <= -2) {		/* number too small without exponent */
    k = (1 - msd) / 3;
    efield = -3 * k;
    while ((--k) >= 0)
      aten *= 1000;
  } else if (lsd < 0) {		/* number too big without exponent */
    k = (1 + msd) / 3;
    if (k > 0)
      efield = 3 * k;
    while ((--k) >= 0)
      aten *= 0.001;
  }
  if (efield) {			/* E format */
    lsd += efield;
    memcpy(chexp, "0fE+00", 6);
    if (efield < 0) {
      chexp[3] = '-';
      efield = -efield;
    }
    msd = efield / 10;
    if (msd > 0) {
      efield -= msd * 10;
      chexp[4] = (char) (msd + CTOI('0'));
    }
    chexp[5] = (char) (efield + CTOI('0'));
    apar *= aten;
    aerr *= aten;
    adif *= aten;
  } else {			/* F format */
    memcpy(chexp, "0f    ", 6);
  }
  if (lsd > 9)
    lsd = 9;
  if (lsd > 0)
    chexp[0] = (char) (lsd + CTOI('0'));
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
  *pfmt = 0;
  sprintf(ptmp, fmt, apar, ie, adif);
  return 0;
}				/* parer */


int linein(luin, nline, iqnfmt)
FILE *luin;
int *nline;
int iqnfmt;
{
  /* Local variables */
  SXLINE *xline;
  double xfrqn, xerrn, xwtn;
  int nqn, nqnu, nqnl, kqnu, kqnl, i, iqf, ipace, mxline, mxqn;
  short nqnt[12], *iqnum;

  /*   get lines from input  and stores them */

  /*     LUIN= unit for finding lines */
  /*     NLINE = number of lines */
  /*     IQNFMT= qunatum number format for line input */
  /*     RETURN: largest quantum number */
  /*******************************************************************/

  mxline = *nline;
  mxqn = 1;
  nqn = deflin(iqnfmt, nqnt);
  nqnu = nqn - 1;
  if (nqnt[nqnu] < 0)
    nqnu = 0;
  kqnu = nqnt[nqnu];
  nqnl = nqnu + nqn;
  kqnl = nqnt[nqnl];
  ipace = 100;

  for (i = 1; i <= mxline; ++i) {	/*       loop for reading lines */
    xline = lbufof(1, i);
    iqnum = xline->qn;
    if (getlin(luin, nqnt, iqnum, &xfrqn, &xerrn, &xwtn)) {
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
    xline->xwt = (float) xwtn;
    xline->linku = 0;
    xline->linkl = 0;
    xline->bln = 1;
    if (ipace <= i) {
      ipace += 100;
      printf("Reading Line %d\n", i);
    }
  }
  return mxqn;
}				/* linein */

int lineix(lu, flg, nline, nblkpf, iqnfmt)
FILE *lu;
int flg, nline, nblkpf, iqnfmt;
{				/*   get lines from input  and stores them */
  /*     LU = unit for printout of lines ( if > 0 ) */
  /*     NLINE = number of lines */
  /*     NBLKPF= number of blocks per F */
  /*     IQNFMT= qunatum number format for line input */
  /******************************************************************/
  static int nsort = 2048;
  static char fmt[] = " %4d%4d%4d%4d%4d:";
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xfrqx, sum;
  float xnorm;
  int nblk, ipos, i, j, ipace, more, blend, nread, nblnd, iblkl, iblku;
  int linkx, indxl, linky, blnbgn, indxu, orgblk, nqn, nqn2, nbad;
  int *prvblk;
  short *iqnum;
  char aqnum[38], card[90];

  nbad = 0;
  nblk = 0;
  i = (nsort + 1) * sizeof(int);
  prvblk = (int *) mallocq((unsigned long) i);
  for (i = 0; i <= nsort; ++i) {
    prvblk[i] = 0;
  }
  nqn = iqnfmt % 10;
  nqn2 = nqn << 1;
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
  sum = 0.;
  blnbgn = 0;
  xfrqx = 1e38;
  nblnd = 0;
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
    xline->ibu = (short) iblku;
    xline->inu = (short) indxu;
    xline->ibl = (short) iblkl;
    xline->inl = (short) indxl;
    if (iblku == 0) {
      /*     .. print out bad line and try for next */
      ++nbad;
      xline->xwt = 0.;
      xwtn = 0.;
      qnfmt(nqn2, iqnum, aqnum);
      sprintf(card, "Bad Line(%d): %s %14.5f %8.5f %8.5f\n",
	      nread, aqnum, xfrqn, xerrn, xwtn);
      fputs(card, lu);
      fputs(card, stdout);
      continue;
    }
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
    /* let the user know something is happening */
    more = (nread < nline);
    if (ipace <= nread || more == 0) {
      ipace += 50;
      printf("Converting Line %d\n", nread);
    }
    if (flg < 0) {
      xline = lbufof(0, nread);
      iqnum = xline->qn;
      fprintf(lu, fmt, nread, iblku, indxu, iblkl, indxl);
      for (i = 0; i < 12; ++i) {
	j = iqnum[i];
	fprintf(lu, "%3d", j);
      }
      fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
    }
    blend = (fabs(xfrqn - xfrqx) < 1e-6);
    if (blend) {
      ++nblnd;
      if (flg < 0) {
	fprintf(lu, "   Line Blended with %3d\n", blnbgn);
      }
      sum += xwtn;
      blend = more;
    } else {
      if (flg < 0)
	fputc('\n', lu);
    }
    if (blend == 0) {
      if (nblnd) {
	xnorm = (float) (1. / sum);
	while (nblnd >= 0) {
	  /*   go back and normalize weights for blend */
	  /*   IBLND=1 (no blend) =-1 (last of blend) =0 (blend) */
	  xline = lbufof(1, blnbgn);
	  xline->bln = 0;
	  xline->xwt *= xnorm;
	  ++blnbgn;
	  --nblnd;
	}
	xline->bln = -1;
      }
      nblnd = 0;
      blnbgn = nread;
      sum = xwtn;
    }
    xfrqx = xfrqn;
  }				/* end loop for converting lines */
  orgblk = nsort;
  while (nblk > orgblk) {	/* finish up links */
    --orgblk;
    linkx = 0;
    for (j = 0; j < nsort; ++j) {
      if (prvblk[j]) {
	linkx = prvblk[j];
	prvblk[j] = 0;
      }
    }
    prvblk[0] = linkx;
    prvblk[nsort] = 0;
    linkx = lnlink(prvblk, nsort, 1, 0);
    while (linkx) {
      linky = linkx;
      getdbk(&linkx, &iblkl, &j, &j, &j);
      iblkl -= orgblk;
      lnlink(prvblk, nsort, iblkl, linky);
    }
    orgblk += nsort;
  }
  free(prvblk);
  return nbad;
}				/* lineix */
