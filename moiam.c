/*  moiam .. program to compute internal rotation structural parameters */
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "calpgm.h"

#define NFT  21
#define CINC 12
#define C0 '\0'
static double zero = 0.;
static double rad, pi2;
static char card[82];

int ftran(int nft, double *xpec, int pec, double *fc, int flg);
static int rdinp(FILE * luin, int ndatm, int *nnatm, int *nx, double *wtx,
                 double *dist, double *theta, double *phi, char *istau, 
                 char **xyztag);
static int coords( /*@null@*/ FILE * lu, int noatm, double *coord, double *wt,
                  int *nx, double *dist, double *theta, double *phi);
static int gcalc(int natm, double *wt, double *coord, double *dcoord,
                 /*@null@*/ double *gval, double *cval);
static double determ(double *a, int ndd, int nr, int nc);
static int iamsolv(double *cval, double *fc);
static int pamsolv(int natm, double *x, double *gval, double *cval);
static int micvt(double *tval, double *gval, double vmuz, double vmuy,
                 int flg);
static double xinterp(int n, double *tauval, double *xtrans);
/*----------------------------------------------------------------------------
  input file: name.inp
    first line => nstruct, angram, mu_top, mu_z, mu_y
               nstruct= number of structures vs tau
               nstruct= -2 means two structures, non-cyclic vibration 
               angram = angle to rotate from IAS to RIAS (degrees)
               mu_top = dipole moment of top perpendicular to z
               mu_z   = dipole moment parallel to z
               mu_y   = dipole moment of frame perpendicular to z
                       (where z is defined by the first two atoms and
                        y is perpendicular to z in the plane of the
                        third atom)
    second line => values of tau for structures
    for each of nstruct structures
                set of structure parameter lines
                blank line
    end for
    NSYM,MAXM,MAXV,NFT (parameters for first line of name.iam)
    n,Vn
      if NSYM <0 then ityp=17,18 parameters included
      if NSYM =0 then PAM found, rather than IAM, ityp=17,18,19 included

  structure parameter lines:
           cwt, n,na,nb,nc,dist,theta,phi,istau
     OR    'X', wt n,na,nb,nc,dist,theta,phi,istau
         cwt(up to 4 characters) => symbol for mass (C,O,C13,..)
         n= atom number (if atom number is in sequence, distances are
                calculated from dist, theta, phi, otherwise
                   dist, theta, phi calculated from structure)
         na,nb,nc atoms defining refernce axes
         abs(dist) is length from atom n to atom na
             if (dist >0) then
               theta is the n-na-nb angle and
               phi is the dihedral angle from the na-nb-nc plane
             if (dist <0) then
                   angles are defined from the nb-na-nc angle bisector
               theta is the n-na-bisector angle
               phi is the angle between the bisector and the projection
                   of n-na into the na-nb-nc plane
          istau ne 0 if phi is linear in tau else if istau = 0 phi is assumed
          to br fourier series in tau
   NOTE: for best results the first three atoms should define the
           plane of symmetry of the trans state

  output file: name.iam
           first line => rho, NSYM,MAXM,MAXV
             other lines => ityp, nf, nm, value
               nf = coefficient in Fourier series in tau
                         (negative for sin, positive for cosine)
               ityp, nm=   99,0 for V                      delta K=0
                           99,2 for F                      delta K=0
                         1099,0 for A - (B+C)/2            delta K=0
                          199,0 for (B+C)/2                delta K=0
                        40099,0 for (B-C)/4                delta K=2
                       210099,0 for D_bc                   delta K=2
                       410099,0 for D_ac                   delta K=1
                       610099,0 for D_ab                   delta K=1
                            1,0 for mu_a (top)             delta K=0
                            2,0 for mu_b (top)             delta K=1
                            3,0 for mu_c (top)             delta K=1
                       200099,1 for Pa                     delta K=0
                       400099,1 for Pb (RAM)               delta K=1
                       600099,1 for Pc (PAM)               delta K=1
   NOTE: frame z is the axis of the first two atoms, and the yz plane
           is the plane of the first three atoms
----------------------------------------------------------------------------*/
int main(argc, argv)
int argc;
char *argv[];
{
#define NFILE 4
#define NDATM 50
#define ND3ATM (3 * NDATM)
#define NDCON 20
  static struct spar {
    int idv;
    char *labl;
  } sparv[] = {
    {
    99, "Vp"}, {
    99, "F"}, {
    1099, "A-(B+C)/2"}, {
    199, "(B+C)/2"}, {
    40099, "(B-C)/4"}, {
    210099, "D_bc"}, {
    410099, "D_ac"}, {
    610099, "D_ab"}, {
    1, "mu_a(top)"}, {
    2, "mu_b(top)"}, {
    3, "mu_c(top)"}, {
    1, "mu_a(z)"}, {
    2, "mu_b(z)"}, {
    3, "mu_c(z)"}, {
    1, "mu_a(y)"}, {
    2, "mu_b(y)"}, {
    3, "mu_c(y)"}, {
    200099, "P_a"}, {
    400099, "P_b"}, {
    600099, "P_c"}};
  static const char *cext[NFILE] = { "inp", "iam", "moi", "xyz" };
  enum efile {einp, eiam, emoi, exyz};
  static char *cfil[NFILE+1];
  static char fmtv[3] = "V0";
  static double gder[2 * NFT], cval[CINC * NFT], gval[NDCON * NFT],
      tauval[NFT];
  static double xcoord[ND3ATM * NFT], dcoord[ND3ATM * NFT], fc[2 * NFT],
      dval[5];
  static double coord[ND3ATM], gtmp[NDCON], vpot[10];
  static double wt[NDATM], dist[NDATM], theta[NDATM], phi[NDATM];
  static char istau[NDATM];
  static char *xyztag[NDATM];
  static int nx[4 * NDATM];
  FILE *luin, *luout, *lulist, *luxyz;
  double *pdbl;
  double ctau, dtau, stau, vmuy, vmuz, dtaud, angram, angorg;
  double vmutop, rho, tmp, cmp, xlintau;
  int iatm, k, mval, itau, nsym, n3atm, i, nnatm, noatm;
  int ia, ii, nval, nstruct, idpar, maxm, maxv, nft, nftc, nxyz;
  BOOL more, cyclic;

  rad = asin(1.) / 90.;
  pi2 = rad * 360.;
  rho = 0.;
  filget(argc, argv, NFILE, cfil, cext);
  luin = fopenq(cfil[einp], "r");
  luout = fopenq(cfil[eiam], "w");
  lulist = fopenq(cfil[emoi], "w");
  luxyz  = fopenq(cfil[exyz], "w");
  //  for (i = 0; i < NDATM; ++i) xyztag[i] = xyztag0;
  /* read: nstruct, angram, mutop, muz, muy */
  dval[0] = 3;
  dval[1] = 0.;
  dval[2] = dval[3] = dval[4] = 1.;
  if (fgetstr(card, 82, luin) <= 0)
    exit(EXIT_FAILURE);
  pcard(card, dval, 5, NULL);
  nstruct = (int) dval[0];
  cyclic = TRUE; nftc = NFT;
  if (nstruct == -2) {
    cyclic = FALSE; nftc = 2; nstruct = 2;
    fputs("non-cyclic vibration\n", lulist);
  }
  angram = dval[1];
  vmutop = dval[2];
  vmuz = dval[3];
  vmuy = dval[4];
  if (nstruct <= 0 || nstruct > NFT) {
    puts("bad nstruct");
    exit(EXIT_FAILURE);
  }
  n3atm = 0;
  noatm = nstruct * ND3ATM;
  dcopy(noatm, &zero, 0, xcoord, 1);
  /* read: tau values */
  i = 0;
  do {
    if (fgetstr(card, 82, luin) <= 0)
      exit(EXIT_FAILURE);
    ii = pcard(card, &tauval[i], nstruct - i, NULL);
    i += ii;
  } while (i < nstruct);
  /* solve for coordinate FT */
  nnatm = NDATM;
  if (cyclic) {
    xlintau = xinterp(nstruct, tauval, dcoord);
  } else {
    tauval[0] = 0.; tauval[1] = 180.; xlintau = 0.;
  }
  pdbl = dcoord;
  for (i = 0; i < nstruct; ++i) { /* read in parameters */
    fprintf(lulist,"coordinates for tau = %8.1f\n", tauval[i]);
    iatm = rdinp(luin, NDATM, &noatm, nx, wt, dist, theta, phi,
                 istau, xyztag);
    if (iatm == 0) {
      printf("coordinate problem for tau = %8.1f\n", tauval[i]);
      exit(EXIT_FAILURE);
    }
    if (iatm < nnatm)
      nnatm = iatm;
    /*   get coordinates */
    coords(lulist, noatm, coord, wt, nx, dist, theta, phi);
    n3atm = 3 * nnatm;
    if (cyclic) {
      for (iatm = 3; iatm < n3atm; ++iatm) {
        ii = iatm / 3; k = iatm % 3;
        if (ii <= k) continue;
        if (k == 0) {
          tmp = dist[ii];
        } else if (k == 1) {
          tmp = theta[ii];
        } else {
          tmp = phi[ii];
          if (istau[ii] != C0) {
            if (i == 0) xcoord[iatm] = tmp;
            tmp -= xcoord[iatm];
            tmp -= 360. * floor(tmp / 360. + 0.5);
            xcoord[iatm + ND3ATM] += xlintau * tmp * (tauval[i] - tauval[0]);
            continue;
          }
        }
        daxpy(nstruct, tmp, pdbl, nstruct, &xcoord[iatm], ND3ATM);
      }
     ++pdbl;
    } else if (i == 0) {
      dcopy(n3atm, coord, 1, xcoord, 1);
    } else {
      for (iatm = 0; iatm < n3atm; ++iatm) 
        dcoord[iatm] = xcoord[iatm] - coord[iatm];
      dcopy(n3atm,  coord, 1, &xcoord[ND3ATM], 1);
      dcopy(n3atm, dcoord, 1, &dcoord[ND3ATM], 1);
    }
  }
  nxyz = 0;
  for (i = 0; i < nnatm; ++i) 
    if (wt[i] >= 0.5) ++nxyz; 
  /* read: nsym, maxm, maxv, nft */
  nsym = 1;
  maxm = 30;
  maxv = 3;
  nft = NFT;
  dval[0] = (double) nsym;
  dval[1] = (double) maxm;
  dval[2] = (double) maxv;
  dval[3] = (double) nft;
  do {
    more = (fgetstr(card, 82, luin) >= 0);
    if (more && pcard(card, dval, 4, NULL) > 0) {
      nsym = (int) dval[0];
      maxm = (int) dval[1];
      maxv = (int) dval[2];
      nft = (int) dval[3];
      break;
    }
  } while (more);
  dcopy(10, &zero, 0, vpot, 1);
  while (more) {
    if (fgetstr(card, 82, luin) <= 0)
      break;
    if (pcard(card, dval, 2, NULL) < 2)
      break;
    i = (int) dval[0];
    if (i > 0 && i <= 9)
      vpot[i] += dval[1];
  }
  fclose(luin);
  dtaud = 360. / (double) nftc;
  dtau = dtaud * rad;
  if (cyclic) {
    dcopy(NFT, &zero, 0, fc, 1);
    pdbl = xcoord;
    for (iatm = 0; iatm < n3atm; ++iatm) {
      ii = iatm / 3;
      if ((iatm % 3) == 2 && istau[ii] != C0) {
        cmp = pdbl[ND3ATM]; tmp = pdbl[0] - cmp * tauval[0];
        cmp *= dtaud;
        pdbl[0] = tmp; k = NFT >> 1;
        for (itau = 1; itau < NFT; ++itau) {
          ii = itau;
          if (itau > k) ii -= NFT;
          /* setup linear phi */
          pdbl[itau * ND3ATM] = tmp + cmp * ii;
        }
      } else {
        dcopy(nstruct, pdbl, ND3ATM, fc, 1);
        fc[nstruct] = 0.;
        ftran(NFT, pdbl, ND3ATM, fc, -1); /* synthesize dist, theta, or phi */
      }
      ++pdbl;
    }
    pdbl = xcoord;
    for (itau = 0; itau < NFT; ++itau) {
      /*   get coordinates */
      /* fprintf(lulist,"coordinates for tau = %8.1f\n",itau * dtaud); */
      dcopy(nnatm, pdbl, 3, dist, 1);
      dcopy(nnatm, pdbl + 1, 3, theta, 1);
      dcopy(nnatm, pdbl + 2, 3, phi, 1);
      coords(NULL, nnatm, pdbl, wt, nx, dist, theta, phi);
      pdbl += ND3ATM;
    }
    pdbl = xcoord;
    for (iatm = 0; iatm < n3atm; ++iatm) {
      ftran(NFT, pdbl, ND3ATM, fc, 0);            /* analyze xcoord */
      ftran(NFT, &dcoord[iatm], ND3ATM, fc, 1);   /* synthesize derivative */
      ++pdbl;
    }
  } 
  /*   get primitive constants vs tau */
  ia = ii = iatm = 0;
  for (itau = 0; itau < nftc; ++itau) {
    gcalc(nnatm, wt, &xcoord[iatm], &dcoord[iatm], &gval[ii], &cval[ia]);
    dcopy(7, &zero, 0, &cval[ia + 4], 1);
    cval[ia + 3] = 1.;
    cval[ia + 7] = 1.;
    cval[ia + 11] = 1.;
    gval[ii + 8] = 0.;
    tmp = itau * dtau;
    gval[ii + 9] = cos(tmp) * vmutop;
    gval[ii + 10] = sin(tmp) * vmutop;
    ia += CINC;
    ii += NDCON;
    iatm += ND3ATM;
  }
  if (nsym == 0) {              /*   solve for PAS */
    ia = ii = 0;
    for (itau = 0; itau < nftc; ++itau) {
      /* diagonalize inertial matrix */
      iatm = itau * ND3ATM;
      pamsolv(nnatm, &xcoord[iatm], &gval[ii], &cval[ia + 3]);
      ia += CINC;
      ii += NDCON;
    }
    if (cyclic) {
      for (iatm = 0; iatm < n3atm; ++iatm) {      /* get new derivatives */
        ftran(NFT, &xcoord[iatm], ND3ATM, fc, 0);
        ftran(NFT, &dcoord[iatm], ND3ATM, fc, 1);
      }
    } else {
      for (iatm = 0; iatm < n3atm; ++iatm) {      /* get new derivatives */
        tmp = xcoord[iatm + ND3ATM] - xcoord[iatm];
        dcoord[iatm] = tmp; dcoord[iatm + ND3ATM] = tmp;
      }
    }
    ia = iatm = 0;
    for (itau = 0; itau < nftc; ++itau) {        /* get new values of C */
      gcalc(nnatm, wt, &xcoord[iatm], &dcoord[iatm], NULL, &cval[ia]);
      ia += CINC;
      iatm += ND3ATM;
    }
    angram = 0.;
  }
  /*   calculate pseudo potential terms */
  if (cyclic) {
    ftran(NFT, gval, NDCON, fc, 0);
    ftran(NFT, gval, NDCON, fc, 1);
    ftran(NFT, gder, 1, fc, 1);
    ftran(NFT, &gval[1], NDCON, fc, 0);
    ftran(NFT, &gder[NFT], 1, fc, 1);
  } else {
    tmp = gval[NDCON] - gval[0];
    gval[0] = tmp; gval[NDCON] = tmp;
    gder[0] = tmp; gder[1] = tmp;
    tmp = gval[NDCON + 1] - gval[1];
    gder[NFT] = tmp; gder[NFT + 1] = tmp;
  }
  fputs("    tau, pseudo potential,        G,", lulist);
  fputs("        Cz,       Cy,       Cx\n", lulist);
  ia = ii = 0;
  for (itau = 0; itau < nftc; ++itau) {
    tmp = gval[ii];
    tmp = 0.25 * (gder[itau + NFT] * tmp +
                  gval[ii + 1] * (gder[itau] + 0.25 * tmp * tmp));
    gval[ii] = tmp;
    tmp = itau * dtaud;
    if (tmp > 180.)
      tmp -= 360.;
    fprintf(lulist, "%7.1f %14.6E %13.2f %9.6f %9.6f %9.6f\n", tmp,
            gval[ii], gval[ii + 1], cval[ia], cval[ia + 1], cval[ia + 2]);
    gval[ii] *= -2.;
    ii += NDCON;
    ia += CINC;
  }
  angorg = 0.;
  if (nsym != 0 && cyclic) {              /*   solve for IAS */
    if (iamsolv(cval, fc) != 0) {
      strcpy(card, "iteration failed");
      fputs(card, lulist);
      fputs(card, luout);
      fputs(card, stdout);
      exit(EXIT_FAILURE);
    }
    rho = cval[0];
    if (rho < 0.) {
      angorg = 180.;
      rho = -rho;
    }
  }
  if (angram >= 360.) {
    dcopy(NDCON, gval, 1, gtmp, 1);
    /* transform moments of inertia to predict inertial direction */
    micvt(cval, gtmp, vmuz, vmuy, nsym);
    tmp = 2. * gtmp[7] / (gtmp[2] - 2. * gtmp[4]);
    angram = 0.5 * atan(tmp) / rad;
  }
  angram += angorg;
  tmp = angram * rad;
  ctau = cos(tmp);
  if (ctau < 1.) {
    stau = sin(tmp);
    ia = NFT * (CINC / 3);
    drot(ia, cval, 3, &cval[1], 3, ctau, stau);
    fprintf(lulist, "%s %9.3f, Ca = %14.10f, Cb = %14.10f\n",
            "Rotation angle for ab plane =", angram, cval[0], cval[1]);
  }
  fprintf(lulist, "rho = %12.9f\n", rho);
  fputs("    tau,  coordinates of Z,", lulist);
  fputs("                   coordinates of Y\n", lulist);
  ia = 3; iatm = 0;
  for (itau = 0; itau < nftc; ++itau) {
    tmp = itau * dtaud;
    if (tmp > 180.)
      tmp -= 360.;
    fprintf(luxyz,"%d\n tau = %7.1f\n", nxyz, tmp);
    for (i = 0; i < nnatm; ++i) {
      if (wt[i] < 0.5) continue;
      pdbl = &xcoord[iatm + 3 * i];
      for (ii = 0; ii < 3; ++ii) {
        dval[ii] = ddot(3, pdbl, 1, &cval[ia + ii], 3);
      }
      fprintf(luxyz,"%s %11.7f %11.7f %11.7f\n", xyztag[i],
               dval[2], dval[1], dval[0]);
    }
    fprintf(lulist, "%7.1f ", tmp);
    for (i = 0; i < 6; ++i) {
      fprintf(lulist, "%11.7f ", cval[ia + i]);
    }
    fputc('\n', lulist);
    ia += CINC; iatm += ND3ATM;
  }
  fclose(luxyz);
  fputs("    tau,  A-(B+C)/2,    (B+C)/2,    (B-C)/4,", lulist);
  fputs("        Dbc,        Dac,        Dab\n", lulist);
  ii = ia = 0;
  for (itau = 0; itau < nftc; ++itau) {
    /* transform moments of inertia */
    micvt(&cval[ia], &gval[ii], vmuz, vmuy, nsym);
    tmp = itau * dtaud;
    if (tmp > 180.)
      tmp -= 360.;
    fprintf(lulist, "%7.1f ", tmp);
    for (i = 2; i <= 7; ++i) {
      fprintf(lulist, "%11.2f ", gval[i + ii]);
    }
    fputc('\n', lulist);
    ii += NDCON;
    ia += CINC;
  }
  fclose(lulist);
  fprintf(luout, "%15.10f %4d %4d %4d %4d\n", cval[0], nsym, maxm, maxv,
          nft);
  /*  output results */
  tmp = ddot(10, vpot, 1, vpot, 1) + 1.e-6;
  cmp = (1.0e-7) * sqrt(tmp);
  for (ia = 0; ia < 10; ++ia) {
    if (fabs(vpot[ia]) > cmp) {
      fmtv[1] = '0' + (char) ia;
      fprintf(luout, "%8d %2d %2d %20.12E/ %s\n", 99, ia, 0, vpot[ia],
              fmtv);
    }
  }
  ia = NDCON;
  if (nsym != 0) {
    --ia;
    if (fabs(cval[1]) < 1e-10)
      --ia;
  }
  for (i = 0; i < ia; ++i) {
    ftran(nftc, &gval[i], NDCON, fc, 0);
    if (i == 0)
      fc[0] = 0.;
    if (i >= 17) {
      cmp = -2.;
      dscal(nftc, cmp, fc, 1);
    }
    cmp = (1.0e-7) * sqrt(ddot(nftc, fc, 1, fc, 1) + tmp);
    idpar = 0;
    for (ii = 0; ii < nftc; ++ii) {
      if (fabs(fc[ii]) > cmp) {
        mval = idpar;
        idpar = sparv[i].idv;
        if (mval != 0)
          idpar = -idpar;
        nval = (ii + 1) >> 1;
        if (ii != 0 && (ii & 1) == 0)
          nval = -nval;
        mval = 0;
        if (i == 1)
          mval = 2;
        if (i >= 17) 
          mval = 1;
        fprintf(luout, "%8d %2d %2d %20.12E/ %s\n", idpar, nval,
                mval, fc[ii], sparv[i].labl);
      }
    }
    tmp = 1.e-6;
  }
  return 0;
}                               /* MAIN */

int rdinp(luin, ndatm, natm, nx, wtx, dist, theta, phi, istau, xyztag)
FILE *luin;
int ndatm, *natm, *nx;
double *wtx, *dist, *theta, *phi;
char *istau;
char **xyztag;
{
#define BLNK (' ')
#define NCARD 9
  static struct swt {
    /*@null@*/ char *symb;
    double wt;
    char *xyztag;
  } swtv[] = { 
    {"C", 12.,"C "}, 
    {"H", 1.00782,"H "}, 
    {"O", 15.99491,"O "}, 
    {"N", 14.00307,"N "},
    {"F", 18.9984,"F "},
    {"H2", 2.01409,"H "},
    {"O18", 17.99915,"O "},
    {"C13", 13.00335,"C "},
    {"N15", 15.00011,"N "},
    {"Cl", 34.96885,"Cl"},
    {"Cl37", 36.96582,"Cl"},
    {"Br79", 78.9183,"Br"},
    {"Br81", 80.9163,"Br"},
    {"S", 31.9722,"S "},
    {"S34", 33.9688,"S "},
    {"P", 30.9744,"P "},
    {"I", 126.9004,"I "},
    {"D", 2.0140,"H "},
    {"X", -1.,"H "},
    {" ", 0.,"H "},
    {NULL, 0.," "}};
  struct swt *pswt;
  int *nxx;
  char *pstr;
  double dvec[NCARD];
  int ioff, k, noatm, nnatm, ich;

  noatm = nnatm = k = 0;
  nxx = nx;
  while (noatm < ndatm) {
    if (fgetstr(card, 82, luin) <= 0) {
      if (noatm != 0)
        break;
      continue;
    }
    for (pswt = swtv; (pstr = pswt->symb) != NULL; ++pswt) {
      for (k = 0; (ich = (int) pstr[k]) != 0 && ich == (int) card[k]; ++k);
      if (ich == 0 && card[k] == BLNK)
        break;                  /* match found */
    }
    if (pstr == NULL) {
      puts(card);
      puts("Weight symbol is not defined. (zero used)");
      pswt = swtv;
      k = 0;
    }
    pstr = &card[k];
    for (k = 0; k < NCARD; ++k)
      dvec[k] = 0.;
    if (pcard(pstr, dvec, NCARD, NULL) <= 0)
      break;
    ioff = 0;
    if (pswt->wt < -0.5) {
      wtx[noatm] = dvec[ioff++];
    } else {
      wtx[noatm] = pswt->wt; xyztag[noatm] = pswt->xyztag;
    }
    for (k = 0; k < 4; ++k) {
      nxx[k] = (int) dvec[ioff++];
      if (k > noatm)
        nxx[k] = 0;
      if (k != 0 && nxx[k] > noatm) {
        puts(card);
        puts("Atom designations are out of order.");
        return 0;
      }
    }
    dist[noatm] = dvec[ioff++];
    theta[noatm] = dvec[ioff++];
    phi[noatm] = dvec[ioff++];
    istau[noatm] = (char) dvec[ioff];
    ++noatm;
    if (nxx[0] == noatm)
      nnatm = noatm;
    if (nxx[0] > noatm)
      return 0;
    nxx += 4;
  }
  *natm = noatm;
  return nnatm;
}                               /* rdinp */

int coords(lu, noatm, coord, wt, nx, dist, theta, phi)
FILE *lu;
int noatm;
double *coord, *wt;
int *nx;
double *dist, *theta, *phi;
{
  double xcm[3], prcor[3], trans[9], xvectr[3], yvectr[3], cb, sb;
  double wtot, cth, rij, rik, ril, sth, tmp, sum, scale;
  int *nxx, i, j, iloop, ia, ib, ic, ii, iap, ibp, icp, ia3, ib3, ic3, ii3;

  nxx = nx;
  /*  THE COORDINATES OF ATOM 1 ARE SET AS THE ORIGIN */
  for (i = 0; i < 3; ++i) {
    coord[i] = 0.;
    coord[i + 3] = 0.;
    xcm[i] = 0.;
  }
  coord[3] = dist[1];
  xcm[0] = wt[1] * dist[1];
  wtot = wt[0] + wt[1];
  ia = ib = ic = ia3 = ii3 = -1;
  /*  THE X AND Y COORDINATE OF ATOM 3 ARE CALCULATED */
  if (noatm > 2) {
    nxx += 8;
    ia = nxx[1] - 1;
    ia3 = 3 * ia;
    cth = cos(theta[2] * rad);
    if (ia == 1)
      cth = -cth;
    coord[6] = coord[ia3] + dist[2] * cth;
    coord[7] = coord[ia3 + 1] + dist[2] * sqrt(1. - cth * cth);
    coord[8] = 0.;
    xcm[0] += wt[2] * coord[6];
    xcm[1] = wt[2] * coord[7];
    wtot += wt[2];
    ii3 = 9;
  }
  /*  THE COORDINATES OF ALL THE OTHER ATOMS ARE CALCULATED */
  ril = scale = 1.;
  for (iloop = 3; iloop < noatm; ++iloop) {
    iap = ia;
    ibp = ib;
    icp = ic;
    nxx += 4;
    ia = nxx[1] - 1;
    ib = nxx[2] - 1;
    ic = nxx[3] - 1;
    ia3 = 3 * ia;
    /*     THE TRANSFORMATION MATRIX IS OBTAINED */
    if (ia != iap || ib != ibp || ic != icp) {
      /*  THE COMPONENTS OF THE X PRIME DIRECTION ARE OBTAINED */
      rij = ril = 0.;
      ib3 = 3 * ib;
      ic3 = 3 * ic;
      for (i = 0; i < 3; ++i) {
        xvectr[i] = tmp = coord[i + ib3] - coord[i + ia3];
        rij += tmp * tmp;
        yvectr[i] = tmp = coord[i + ic3] - coord[i + ia3];
        ril += tmp * tmp;
      }
      rij = sqrt(rij);
      scale = 0.;
      for (i = 0; i < 3; ++i) {
        trans[i] = xvectr[i] / rij;
        scale += yvectr[i] * trans[i];
      }
      /*  THE COMPONENTS OF THE Y PRIME DIRECTION ARE OBTAINED */
      rik = ril - scale * scale;
      if (rik < 1e-8) {         /* Y PRIME PARALLEL TO X PRIME */
        for (i = 0; i < 3; ++i)
          yvectr[i] = 0.;
        i = 1;
        if (fabs(trans[1]) > 0.71)
          i = 2;
        yvectr[i] = 1.;
        scale = trans[i];
        ril = 1.;
        rik = ril - scale * scale;
      }
      rik = sqrt(rik);
      for (i = 0; i < 3; ++i) {
        trans[i + 3] = (yvectr[i] - scale * trans[i]) / rik;
      }
      /*  THE COMPONENTS OF THE Z DIRECTION ARE OBTAINED */
      trans[6] = trans[1] * trans[5] - trans[2] * trans[4];
      trans[7] = trans[2] * trans[3] - trans[0] * trans[5];
      trans[8] = trans[0] * trans[4] - trans[1] * trans[3];
    }
    /* THE SPHERICAL COORDINATES ARE TRANSFORMED TO CARTESIAN COORDINATES
       IN THE PRIME SYSTEM */
    if (nxx[0] == iloop + 1) {
      if (dist[iloop] < 0.) {
        tmp = dist[iloop] * cos(theta[iloop] * rad);
        cth = 4. * ril;
        cth = scale / sqrt(cth);
        sth = tmp * sqrt(0.5 - cth);
        cth = tmp * sqrt(cth + 0.5);
        sb = sin(phi[iloop] * rad);
        cb = cos(phi[iloop] * rad);
        prcor[0] = cth * cb - sth * sb;
        prcor[1] = sth * cb + cth * sb;
        prcor[2] = -dist[iloop] * sin(theta[iloop] * rad);
      } else {
        cth = cos(theta[iloop] * rad);
        tmp = dist[iloop] * sqrt(1. - cth * cth);
        prcor[0] = dist[iloop] * cth;
        prcor[1] = tmp * cos(phi[iloop] * rad);
        prcor[2] = tmp * sin(phi[iloop] * rad);
      }
      /*     THE PRIME COORDINATES ARE TRANSLATED AND ROTATED TO THE
         COORDINATE SYSTEM DESIRED */
      wtot += wt[iloop];
      for (i = 0; i < 3; ++i) {
        sum = 0.;
        ii = i;
        for (j = 0; j < 3; ++j) {
          sum += trans[ii] * prcor[j];
          ii += 3;
        }
        coord[ii3] = coord[i + ia3] + sum;
        xcm[i] += wt[iloop] * coord[ii3];
        ++ii3;
      }
    } else {                    /*  SPECIAL REVERSE CALCULATION OF DIST THETA AND PHI */
      rij = 0.;
      ic3 = 3 * (nxx[0] - 1);
      ii = 0;
      for (i = 0; i < 3; ++i) {
        coord[ii3] = coord[i + ic3];
        sum = 0.;
        for (j = 0; j < 3; ++j) {
          sum += trans[j + ii] * (coord[j + ic3] - coord[j + ia3]);
        }
        rij += sum * sum;
        prcor[i] = sum;
        ii += 3;
        ++ii3;
      }
      rij = sqrt(rij);
      phi[iloop] = 0.;
      if (dist[iloop] < 0.) {
        sth = prcor[2] / rij;
        dist[iloop] = -rij;
        theta[iloop] = asin(sth) / rad;
        if (fabs(sth) < 0.99999999999999) {
          cth = 4. * ril;
          sth = -sqrt(0.5 - cth);
          cth = -sqrt(cth + 0.5);
          sb = prcor[1] * cth - prcor[0] * sth;
          cb = prcor[0] * cth + prcor[1] * sth;
          phi[iloop] = atan2(sb, cb) / rad;
        }
      } else {
        dist[iloop] = rij;
        sth = prcor[0] / rij;
        theta[iloop] = acos(sth) / rad;
        if (fabs(sth) < 0.999999999999) {
          phi[iloop] = atan2(prcor[2], prcor[1]) / rad;
        }
      }
    }
  }
  if (wtot > 1.) {
    xcm[0] /= wtot;
    xcm[1] /= wtot;
    xcm[2] /= wtot;
  }
  ii3 = 0;
  for (iloop = 0; iloop < noatm; ++iloop) {
    for (i = 0; i < 3; ++i) {
      coord[ii3] -= xcm[i];
      ++ii3;
    }
  }
  if (lu == NULL) return 0;
  fprintf(lu, "WT,N,NA,NB,NC,DIST,THETA,PHI,X,  MASS= %10.6f\n", wtot);
  ii3 = 0;
  nxx = nx;
  for (iloop = 0; iloop < noatm; ++iloop) {
    fprintf(lu, "%10.6f", wt[iloop]);
    for (ii = 0; ii < 4; ++ii) {
      fprintf(lu, " %2d", *nxx);
      ++nxx;
    }
    fprintf(lu, " %8.4f %8.3f %8.3f", dist[iloop], theta[iloop],
            phi[iloop]);
    for (i = 0; i < 3; ++i) {
      fprintf(lu, " %9.5f", coord[ii3]);
      ++ii3;
    }
    fputc('\n', lu);
  }
  return 0;
}                               /* coords */

int gcalc(natm, wt, coord, dcoord, gval, cval)
int natm;
double *wt, *coord, *dcoord, *gval, *cval;
{
  static double tiny = 1.0e-30;
  static double conv = 505379.07;
  double mi[12], x[3], y, tmp;
  int iatm, i, ii, j, k, iatm3;

  y = tiny;
  dcopy(12, &zero, 0, mi, 1);
  iatm3 = 0;
  for (iatm = 0; iatm < natm; ++iatm) {
    /* z=0, y=1, x=2 */
    ii = 0;
    for (i = 0; i < 3; ++i) {
      j = i - 1;
      if (j < 0)
        j = 2;
      k = i + 1;
      if (k > 2)
        k = 0;
      tmp = wt[iatm] * coord[iatm3 + i];
      mi[k + 9] += tmp * dcoord[iatm3 + j];
      mi[j + 9] -= tmp * dcoord[iatm3 + k];
      for (j = 0; j < 3; ++j) {
        mi[ii] -= tmp * coord[iatm3 + j];
        ++ii;
      }
      y += wt[iatm] * dcoord[iatm3 + i] * dcoord[iatm3 + i];
    }
    iatm3 += 3;
  }
  dcopy(3, &mi[9], 1, x, 1);
  tmp = mi[0] + mi[4] + mi[8];
  mi[0] -= tmp;
  mi[4] -= tmp;
  mi[8] -= tmp;
  tmp = determ(mi, 3, 3, 4);
  dcopy(3, &mi[9], 1, cval, 1);
  if (gval) {
    y -= ddot(3, cval, 1, x, 1);
    tmp *= y;
    if (tmp < tiny)
      tmp = tiny;
    gval[0] = log(tmp);
    gval[1] = conv / y;
    gval[2] = conv * mi[0];
    gval[3] = conv * mi[4];
    gval[4] = conv * mi[8];
    gval[5] = conv * mi[5];
    gval[6] = conv * mi[2];
    gval[7] = conv * mi[1];
  }
  return 0;
}                               /* gcalc */

double determ(a, ndd, nr, nc)
double *a;
int ndd, nr, nc;
{ /* calculate inverse & determinant */
  /* use dumb Gauss-Jordon method with pivot */
  static double eps = 1.e-30;
  //   static double aa[36];
  static int nsav;
  /*@owned@*/ static int *iperm;
  double det, s;
  int i, ic, ir, iic, n, iswap;
  //   n = ndd * nr;
  //   dcopy(n, a, 1, aa, 1);
  n = (nc > nr) ? nr : nc;
  if (n > nsav) {
    if (nsav != 0) {
      free(iperm);
      iperm = NULL;
    }
    iperm = (int *) mallocq(n * sizeof(int));
    nsav = n;
    iperm[0] = 0;
  }
  det = 1.;
  iic = 0;
  for (ic = 0; ic < n; ++ic) {
    i = ic + iic;
    iswap = ic + (int) idamax(n - ic, &a[i], 1);
    iperm[ic] = iswap;
    if (iswap != ic) {
      dswap(nr, &a[iic], 1, &a[iswap * ndd], 1);
      det = -det;
    }
    s = a[i];
    a[i] = 1.;
    det *= s;
    if (fabs(s) < eps) {
      dcopy(nc, &zero, 0, &a[ic], ndd);
    } else {
      s = 1. / s;
      dscal(nc, s, &a[ic], ndd);
    }
    for (ir = 0; ir < nr; ++ir) {
      if (ir == ic)
        continue;
      i = ir + iic;
      s = -a[i];
      a[i] = 0.;
      daxpy(nc, s, &a[ic], ndd, &a[ir], ndd);
    }
    iic += ndd;
  }
  for (ic = n - 2; ic >= 0; --ic) {
    iswap = iperm[ic];
    if (iswap != ic)
      dswap(nc, &a[ic], ndd, &a[iswap], ndd);
  }
  /* test
     for (ic = 0; ic < nr; ++ic) {
     iic = 0;
     for (ir = 0; ir < nr; ++ir) {
     s = ddot(nr, &a[ic], ndd, &aa[iic], 1);
     if (ir == ic) s -= 1.;
     iic += ndd;
     }
     }
   */
  return det;
} /* determ */

int iamsolv(cval, fc)
double *cval, *fc;
{
  static double a[2 * NFT], ader[2 * NFT];
  static double eps = 1e-20;
  double rho, rhosq, ang, rho0, z, r, s, c;
  double *fcx, *b, *bder;
  int itr, k, kk, kz, nab, nc, nt, pivot;
  BOOL done;

  fcx = fc + NFT;
  b = a + NFT;
  bder = ader + NFT;
  nab = 2 * NFT;
  nc = CINC / 3;
  nt = nc * NFT;
  rho0 = 0.;
  for (itr = 0; itr < 100; ++itr) {
    kz = 0;
    z = cval[0];
    c = cval[1];
    s = cval[2];
    for (k = 1; k < NFT; ++k) {
      kz += CINC;
      z += cval[kz];
      c += cval[kz + 1];
      s += cval[kz + 2];
    }
    /*  apply tau-independent rotation */
    rhosq = s * s + c * c;
    if (rhosq > eps) {
      r = sqrt(rhosq);
      if (c < 0.)
        r = -r;
      c /= r;
      if (c < 1.) {
        s /= r;
        drot(nt, &cval[1], 3, &cval[2], 3, c, s);
      }
    } else {
      r = 0.;
    }
    rhosq += z * z;
    rho = sqrt(rhosq);
    if (rhosq > eps) {
      if (z < 0.)
        rho = -rho;
      c = z / rho;
      if (c < 1.) {
        s = r / rho;
        drot(nt, cval, 3, &cval[1], 3, c, s);
      }
    }
    if (itr == 0)
      rho0 = rho / NFT;
    r = ddot(NFT, &cval[1], CINC, &cval[1], CINC) +
        ddot(NFT, &cval[2], CINC, &cval[2], CINC);
    done = (r < eps);
    /* analyze Cz */
    ftran(NFT, cval, CINC, fc, 0);
    rho = fc[0];
    r = rho;
    if (done)
      r -= rho0;
    r = floor(r + 0.5);
    pivot = (int) r;
    if (pivot == 0) {
      fc[0] = 0;
      if (ddot(NFT, fc, 1, fc, 1) > eps)
        pivot = 1;
    }
    if (pivot != 0) {           /* rotate about z to make Cz constant */
      rho -= r;
      r *= pi2 / NFT;
      ftran(NFT, a, 1, fc, 2);  /* integrate Cz */
      kz = 0;
      for (k = 0; k < NFT; ++k) {
        ang = a[k] + k * r;
        c = cos(ang);
        s = sin(ang);
        drot(nc, &cval[kz + 2], 3, &cval[kz + 1], 3, c, s);
        cval[kz] = rho;
        kz += CINC;
      }
    }
    if (done)
      return 0;
    /*  analyze Cx, Cy */
    ftran(NFT, &cval[1], CINC, fc, 0);
    ftran(NFT, &cval[2], CINC, fcx, 0);
    r = rho * rho + eps;
    fc[0] *= rho / r;
    fcx[0] *= -rho / r;
    k = 0;
    for (kk = 1; kk < NFT; kk += 2) {
      r -= kk;
      ++k;
      c = fc[kk];
      s = fcx[kk + 1];
      fc[kk] = (rho * c - k * s) / r;
      fcx[kk + 1] = (k * c - rho * s) / r;
      c = fcx[kk];
      s = fc[kk + 1];
      fc[kk + 1] = (rho * s + k * c) / r;
      fcx[kk] = -(k * s + rho * c) / r;
    }
    /*  synthesize angles and derivatives */
    ftran(NFT, a, 1, fc, -1);
    ftran(NFT, ader, 1, fc, 1);
    ftran(NFT, b, 1, fcx, -1);
    ftran(NFT, bder, 1, fcx, 1);
    k = (int) idamax(nab, a, 1);
    r = 8. * fabs(a[k]);
    if (r > 1.) {               /* scale angles for linearity */
      r = 1. / r;
      dscal(nab, r, a, 1);
      dscal(nab, r, ader, 1);
    }
    /* do rotations */
    kz = 0;
    for (k = 0; k < NFT; ++k) {
      ang = a[k];
      c = cos(ang);
      s = sin(ang);
      drot(nc, &cval[kz], 3, &cval[kz + 1], 3, c, s);
      cval[kz + 2] += ader[k];
      ang = b[k];
      c = cos(ang);
      s = sin(ang);
      drot(nc, &cval[kz + 2], 3, &cval[kz], 3, c, s);
      cval[kz + 1] += bder[k];
      kz += CINC;
    }
  }
  return -1;
}                               /* iamsolv */

int pamsolv(natm, x, gval, tval)
int natm;
double *x, *gval, *tval;
{
  double mi[9], te[3], wk[3];
  int j, k, m, itmp;
  short isblk[4], ipasgn[3], iqnsep[3];

  for (j = 0; j < 4; ++j) {
    isblk[j] = (short) j;
  }
  for (k = 0; k < 3; ++k) {
    ipasgn[k] = 0;
    j = k - 1;
    if (j < 0)
      j = 2;
    m = k + 1;
    if (m > 2)
      m = 0;
    mi[j + k * 3] = gval[m + 5];
    mi[k + j * 3] = gval[m + 5];
    mi[k * 4] = gval[k + 2];
  }
  /*  diagonalize and find max overlap */
  hdiag(3, 3, mi, te, wk, iqnsep);
  if ((ordblk(3, 3, iqnsep, mi, te, isblk, wk, ipasgn) & 1) != 0) {
    for (j = 6; j < 9; ++j) {   /* make determinant positive */
      mi[j] = -mi[j];
    }
  }
  /*    apply transformations */
  itmp = 0;
  for (k = 0; k < natm; ++k) {
    dcopy(3, &x[itmp], 1, te, 1);
    for (j = 0; j < 3; ++j) {
      x[itmp] = ddot(3, &mi[j * 3], 1, te, 1);
      ++itmp;
    }
  }
  itmp = 0;
  for (j = 0; j < 3; ++j) {
    dcopy(3, &mi[itmp], 1, &tval[j], 3);
    itmp += 3;
  }
  return 0;
}                               /* pamsolv */

int micvt(cval, gval, vmuz, vmuy, flg)
double *cval, *gval, vmuz, vmuy;
int flg;
{
  double mi[9], tt[9], tv[3], tmp, ctmp, *mt;
  int j, k, m, itmp;

  ctmp = cval[0];
  if (flg > 0)
    cval[0] = 0.;

  for (k = 0; k < 3; ++k) {
    j = k - 1;
    if (j < 0)
      j = 2;
    m = k + 1;
    if (m > 2)
      m = 0;
    mi[j + k * 3] = gval[m + 5];
    mi[k + j * 3] = gval[m + 5];
    mi[k * 4] = gval[k + 2];
    tv[k] = gval[k + 8];
  }
  mt = &cval[3];
  for (j = 0; j < 3; ++j) {
    itmp = j;
    for (k = 0; k < 3; ++k) {
      tt[itmp] = ddot(3, &mi[j], 3, &mt[k], 3);
      itmp += 3;
    }
  }
  for (j = 0; j < 3; ++j) {
    itmp = 0;
    for (k = 0; k < 3; ++k) {
      mi[itmp + j] = ddot(3, &tt[itmp], 1, &mt[j], 3);
      itmp += 3;
    }
  }
  for (k = 0; k < 3; ++k) {
    j = k - 1;
    if (j < 0)
      j = 2;
    m = k + 1;
    if (m > 2)
      m = 0;
    tmp = gval[1] * cval[k];
    gval[k + 2] = mi[k * 4] + tmp * cval[k];
    gval[m + 5] = mi[j + k * 3] + tmp * cval[j];
    gval[k + 8] = ddot(3, &mt[k], 3, tv, 1);
    gval[k + 11] = mt[k] * vmuz;
    gval[k + 14] = mt[k + 3] * vmuy;
    gval[k + 17] = tmp;
  }
  tmp = 0.5 * (gval[3] + gval[4]);
  gval[3] = tmp;
  gval[2] -= tmp;
  gval[4] = 0.5 * (tmp - gval[4]);
  cval[0] = ctmp;
  return 0;
}                               /* micvt */

double xinterp(n, tauval, xtrans)
int n;
double *tauval, *xtrans;
{
  double ang, ang1, tt;
  int i, k, ii;
  tt = 0.;
  ii = 0;
  for (i = 0; i < n; ++i) {
    ang1 = tauval[i] - tauval[0];
    ang1 -= 360. * floor(ang1 / 360 + 0.5);
    tauval[i] = tauval[0] + ang1;
    ang = tauval[i];
    tt += ang1 * ang1;
    xtrans[ii] = 1.;
    ang *= rad; ang1 = ang;
    for (k = 1; k < n; ++k) {
      ++ii;
      if ((k & 1) != 0) {
        xtrans[ii] = cos(ang);
      } else {
        xtrans[ii] = sin(ang);
        ang += ang1;
      }
    }
    ++ii;
  }
  ang1 = determ(xtrans, n, n, n);
  if (n < 2) tt = 1.;
  tt = 1. / tt;
  return tt;
}                               /* xinterp */
