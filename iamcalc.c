/*  Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */
/*   Herbert M. Pickett, 20 March 1989 */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */ 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "calpgm.h"
#include "ftran.h"
#include "readopt.h"

#define MYDEBUG 0
#define MYDEBUG_EIG 0
#define NLABL 10
#define NFT   20
#define MPE  1 /* potntial egy multiplied by -0.5 */
#define MWN  2 /* wave number units */
#define MREV 4 /* primary operator has vx = 1 */
#define MV12 8 /* operator abs(nv) >= 25 */
#define MHAM 16 /* parameter is in hamiltonian */
#define ND2TOP 7
#define ODD(x) (((x) & 1) != 0)

typedef struct sxpec {
  /*@owned@*/ /*@null@*/ double *xpec;
  /* vector of expectation values */
  double xval; /* sum of parameter values used in Hamiltonian */
  int nval, mval, kdel, ixvib, isunit;
  /* for first element, nval = nsym, mval = max(m)/nsym, 
                        kdel = number of parameters in Hamiltonian 
                        ixvib = number of extra vibrations
                        isunit = number of parameter types */
} SXPEC;
typedef /*@out@*/ SXPEC *PSXPEC;

typedef struct spar {
  /*@dependent@*/ /*@null@*/ SXPEC *ptxpec, *ptxpecx;
  double parval, parerr;
  long idval;
  int  offset, ksym, flg;
  char labl[NLABL + 1];
} SIPAR;
typedef /*@out@ */ SIPAR *PSIPAR;

typedef struct seig {
  double rhok[9];
  /*@owned@*/ /*@notnull@*/ double *eigvec;
  /*@dependent@*/ int *morg;
  int irhok;
} SEIG;
typedef /*@out@*/ SEIG *PSEIG;

typedef struct s2top {
  int ix1, ix2, iqv, iqsig, inew;
} S2TOP;
 
static double zero = 0.;
static double eps = 1.e-35;
static double tiny = 1.e-12;
static double clight = 29979.2458;
static double rthalf =0.707;
static int morgsym[10];

int iamdiag( /*@null@*/ FILE * lu, int nvib, SXPEC ** pxpecv,
             /*@out@*/ SEIG * eig, /*@out@*/ double *egy,
             int ioff);
int expec(int ityp, int indx, PSXPEC *pxpecv, int nvib,
                 SEIG *eig, SEIG *eigp, int nxpec, int *isym);
int getparv(int nv, int mv, int iv, int nvib2, /*@out@*/ SIPAR *parp,
                   PSXPEC *pxpecv);
int ideqv(long *idtmp, /*@out@*/ double *fac);
int init2top(FILE *lu, int nvib2, int nsym, int nvib, int nvtot,
	     int nsig, short *kmap, S2TOP *vmap, 
             /*@out@*/int *sigma2, /*@out@*/int *ivtop);
int get2top(SIPAR *ppar,/*@out@*/ double *xpec2, int nvtot, int nft,
            int nsym, int *vmapsym, int *ivtop, int nvtop);
int wrpar(/*@null@*/FILE * lu, double *id, int nid, int isdip, int kflg,
          double *t, int nft, char *labl, int nfmt, double drhok);
int ftanal(int flg, int nft, double *xpec, /*@out@ */ double *fc);
int ksign(int kdel, int nft, double *xpec, /*@out@ */ double *fc);
void ftfix(int nftp2, int nftp, double *pxfc);
void dneg(int n, double *dx, int incx);
void sxpec_alloc( /*@out@*/ PSXPEC *pxpecv, const int nvec);
void spar_alloc(/*@out@*/ SIPAR *parwk, /*@out@*/ PSIPAR *par, const int nvec);
void seig_alloc(PSEIG eigv, const int nvec, const size_t nfill,
                       /*@dependent@*/ int *pmorgwk, const int nvib);

/*       program to calculate expectation values for iam basis */
/*  input file: name.iam */

/*  first line => rho, nsym, maxm, nvib, maxfc */
/*        abs(nsym)  = symmetry of basis (zero for no K series) */
/*                     (negative nsym for calculation at K = 0 only) */
/*        maxm       = maximum value of m for the basis */
/*        abs(nvib)  = number of vibrational states to include in the output*/
/*                     (neg nvib for no FC) */
/*        abs(maxfc) = maximum number of sample points in (rho * K) */
/*                     rounded up so mod(maxfc, abs(nsym)) = 0  */
/*                     (negative maxfc for sampling of K) */

/*  other lines => ityp, nf, nm, value */
/*        nf = coefficient in Fourier series in tau */
/*              (negative for sin, positive for cosine) */
/*        if |nf| >= 25 then nf = nf mod 25 and is special potential term */
/*        ityp, nm=          99,0 for V                      delta K=0 */
/*                           99,2 for F                      delta K=0 */
/*                         1099,0 for A - (B+C)/2            delta K=0 */
/*                          199,0 for (B+C)/2                delta K=0 */
/*                        40099,0 for (B-C)/4                delta K=2 */
/*                       210099,0 for D_bc                   delta K=2 */
/*                       410099,0 for D_ac                   delta K=1 */
/*                       610099,0 for D_ab                   delta K=1 */
/*                            1,0 for mu_a                   delta K=0 */
/*                            2,0 for mu_b                   delta K=1 */
/*                            3,0 for mu_c                   delta K=1 */
/*                       200099,1 for Pa                     delta K=0 */
/*                       400099,1 for Pb (RAM)               delta K=1 */
/*                       600099,1 for Pc (PAM)               delta K=1 */
/*               if nm >= 100, nmp = nm - 100 and define as dipole     */
/*  output file name.par */
/*       rho */
/*       V_n */
/*  
      Kinetic energy */
/*       expectation for ityp=2..10 (all Fourier series in tau combined) */

int main(int argc, char **argv);

int main(argc, argv)
     int argc;
     char *argv[];
{
#define LIN   0
#define LPAR  1
#define LOUT  
#define LINT  3
#define LBAK  4
#define NFILE 5
#define NRHOK 5
  static double dval[8], rhokv[NRHOK * 9], rhov[10], drhokv[9];
  static double *eigref[9];
  static int krhok[NRHOK], ilink[9];


  static const char *cext[NFILE] = { "iam", "par", "out", "int", "bak"};
  enum efile {eiam, epar, eout, eint, ebak};
  char *cfil[NFILE+1];
  static char card[82], title[82];
  FILE *luin, *lupar, *luout, *luint, *lubak, *lu;
  /*@owned@*/ SXPEC **pxpecv;
  /*@owned@*/ SXPEC *pxpecpwk;
  SXPEC *pxpecp;
  /*@owned@*/ SIPAR **par;
  /*@owned@*/ SIPAR *parwk;
  SIPAR *parp, *parq, *parr;
  SEIG *eigp, *eig, *eig0;
  /*@owned@*/ SEIG *eigv; 
  /*@owned@*/ S2TOP *vmap;
  /*@owned@*/ double *eigval, *eigvalk, *eigvalk2, *xfc, *rhokdif, *xpec2top;
  double *egy, *egyp, *egywk, *eigvec0, *eigvec1, *eigvec0p, *eigvec1p, *pdtmp;
  double *pxfc, *pxfc1, *eigvec;
  /*@owned@*/ double *txpec;
  /*@owned@*/ int *sigma, *isym, *pmorgwk, *ioffv, *ixeig, *ixpt0, *ivtop;
  /*@owned@*/ int *sigma2, *vmapsym;
  /*@owned@*/ short *kmap;
  int *pmorg;
  char *pstr, *pbgn;
  double rho, drhok, val, voff, vref, vmin, vmax, didtmp, valm[2], zval;
  double rhoang, ang, ang0, amp;
  size_t nl;
  long idtyp, idtmp, idtst, nnpar;
  /*@owned@*/ double *didtmpv;
  int nsize_p, iend, kflg, ipar, nvibsq, maxm, npar, nvib, ift, ncalc, nxvib;
  int ityp, nsym, ntyp, i, k, nbasis, negy, nvfac, nrhok, nft, nft2, nftp, m;
  int nvfac2, kd, kk, iv, ivp, ivv, ivvp, mv, nv, nkk, neig, nxpec, ldel;
  int isig, isigp, ivx, ivxp, nvs, nvp, kdel, kdel0, kmax, kmax0, ktest, ntyp0; 
  int neigp, nhalf, nsymb, isymtau, ksymtau, kdoff, nsin, mbasis, nsig, nvtot2;
  int kflg0, kflgm, nprn, isdip, nvtot, nvtotsq, nvtot0, nxrho, ncalcm, nvt2;
  int ixv, nxblks, nxdim, isigx, ixvv, nnblk, ndip, nvib2, nvt, nvtsq, nsymsq;
  int nvib0, mm, nsigsq, ii, nline2, vtmax, nvtop, ndtop, ndmtop, nneg;
  short iv1,iv2;
  char cval;
  BOOL hsym, negrho, nomorg, nofc, isoffd, headok, lsym, rezero;

  nsize_p = maxmem(&nl);
  rthalf = sqrt(0.5);
  filget(argc, argv, NFILE, cfil, cext);
  luin = fopenq(cfil[eiam], "r");
  luout = fopenq(cfil[eout], "w");
  kflg0 = 0; kflgm = -1;
  /* read input */
  if (fgetstr(card, 82, luin) <= 0)
    exit(EXIT_FAILURE);
  dval[0] = 0.;
  dval[1] = 1.;
  dval[2] = 30.;
  dval[3] = 1.;
  dval[4] = NFT;
  dval[5] = 0;
  dval[6] = -1.;
  pcard(card, dval, 6, NULL);
  rho = dval[0];
  nsym   = (int) dval[1];
  maxm   = (int) dval[2];
  nvib   = (int) dval[3];
  nft    = (int) dval[4];
  nline2 = (int) dval[5];
  ktest  = (int) dval[6];
  dcopy(10, &rho, 0, rhov,  1);
  ktest = ktest + ktest + 1;
  negrho = FALSE; nomorg = FALSE; nofc = FALSE;
  kmax = -1; kmax0 = nhalf = 0;
  if (nvib < 0) {
    nofc = TRUE; nvib = -nvib; ktest = 0;
    fputs(".par file does not use Fourier coefficients\n",luout);
  }
  if (nsym > 0) {
    nrhok = NRHOK;
    kflg = 100 + nsym;
    if (nft < 0) {
      kmax = -nft;
      if (kmax > 49) kmax = 49;
      kmax0 = kmax;
      nft = kmax + kmax + 1;
      drhok = rho; voff = 1. / (double) nft;
      nhalf = kmax;
      fputs("sampling at integer K\n",luout);
    } else {
      if (nft > 99) nft = 99;
      drhok = (double) nsym / nft;
      voff = rho / nsym; 
      negrho = (rho < 0.);
      fputs("sampling for periodic rho K\n",luout);
    }
  } else {
    nrhok = 1;
    nsym = -nsym;
    if (nsym == 0)
      nsym = 1;
    nft = 1;
    drhok = rho; voff = rho / nsym;
    kflg = 100;
    fputs("single sample in low rho K limit\n",luout);
  }
  lsym = (nsym > 1);
  if (ktest > nft) ktest = nft;
  ftfix(0, 0, &voff);
  for (k = 0; k < 9; ++k) {
    morgsym[k] = 0; ilink[k] = k;
  }
  nsymb = -10000;
  if ((nsym & 1) == 0) 
    nsymb = nsym >> 1;
  nsymsq = nsym * nsym;
  nl = (size_t) (nline2 * ND2TOP + 1);
  /* save 2-top info */
  kmap = (short *)mallocq (nl * sizeof(short));
  kmap[0] = (short) nline2; 
  ivv = 1; nsig = 1; vtmax = 0; nvib2 = 0; iv1 = 0; iv2 = 0;
  for (k = 0; k < nline2; ++k, ivv += ND2TOP){
    if (fgetstr(card, 82, luin) < 1) break;
    dcopy(ND2TOP, &zero, 0, dval,1);
    if (pcard(card, dval, ND2TOP, NULL) < 1) continue;
    for (kk = 0; kk < ND2TOP; ++kk) {
      kmap[ivv + kk] = (short) dval[kk];
    }
    /* find nsig */
    if (k == 0) {
      iv1 = kmap[ivv + 5]; iv2 = kmap[ivv + 6];
    } else if (kmap[ivv + 5] == iv1 && kmap[ivv + 6] == iv2) {
      ++nsig;
    }
    kk = (int) kmap[ivv] + 1;
    if (nvib2 < kk) nvib2 = kk;
    kk = (int) kmap[ivv + 5] + 1;
    if (vtmax < kk) vtmax = kk;
    kk = (int) kmap[ivv + 6] + 1;
    if (vtmax < kk) vtmax = kk;
  }
  nvt2 = 0; nvtot2 = 0; nvtop = 0;
  if (nvib2 > 0) {
    nvt2 = (nvib2 - 1) / nsig + 1; nvib2 = nvt2 * nsig;
    nvtot2 = nvib2; vtmax *= nsym;
    if (nvib < vtmax) nvib = vtmax;
  }
  /* scan for number of parameters */
  npar = 1;
  ndip = 0;
  ntyp = 1;
  nxvib = 0;
  isoffd = FALSE;
  nxrho = 0;
  ivx = 0;
  ivxp = 0;
  while (fgetstr(card, 82, luin) > 0) {
    if (pcard(card, dval, 5, NULL) < 4) continue;
    idtmp = (long) dval[0];
    nv = (int) dval[1];
    mv = (int) dval[2];
    if (idtmp < 0) idtmp = -idtmp;
    if (idtmp < 9L && mv < 100) {
      mv += 100; idtmp += 990;
    }
    iv = 0;
    if (mv >= 100) {
      ++ndip;
      ivv = (int)(idtmp / 10); idtmp -= ivv * 10;
    } else {
      ivv = (int)(idtmp % 100L); idtmp += 99 - ivv;
    }
    ivp = ivv / 10; iv = ivv - 10 * ivp;
    if (iv < ivp) {
      i = ivp; ivp = iv; iv = i;
    }
    if (ivv != 99) {
      if (iv >= 9) continue;
      if (nxvib < iv) nxvib = iv;
    }
    if (nv < 0) 
      nv = -nv;
    if (nv >= 25)
      nv -= 25;
    if (lsym && (nv % nsym) != 0) lsym = FALSE; 
    if (ivp != iv && idtmp == 99L && mv == 0) {
      isoffd = TRUE;
      if (nsymb > 1 && ilink[ivp] != ilink[iv]) {
        i = (nsym - nv) % nsym;
        if (i < 0) i += nsym;
        kk = ilink[iv];
        for (k = 1; k <= nxvib; ++k) {
          if (ilink[k] == kk) {
            ilink[k] = ilink[ivp]; 
            morgsym[k] = (morgsym[k] + i) % nsym;
          }
        }
      }
    }
    if (idtmp == 910099L) {
      if (ivp == iv && iv < 9) { 
        rhov[ivv] = dval[3]; ++nxrho;
      }
      continue;
    }
    ++npar; ++ntyp; 
    if (ivv == 99) {
      ++ivx; ++ivxp; 
    } else if (iv != ivp) {
      ++npar; ++ntyp; 
    }
    if (ideqv(&idtmp, &val) != 0){
      ++npar; ++ntyp; 
      if (ivv == 99) {
        ++ivx;
      } else if (iv != ivp) {
        ++npar; ++ntyp; 
      }
    }
  }
  npar += nxvib * ivx;
  ntyp += nxvib * (ivxp + nrhok);
  ++nxvib; 
  nxblks = nxvib; nxdim = 1;
  if (isoffd) {
    nxblks = 1; nxdim = nxvib;
  }
  if (nxrho > 0) {
    nxrho = 0;
    for (k = 1; k < nxvib; ++k) {
      if (fabs(rhov[k] - rhov[0]) > 1.e-10) {
        if (isoffd) {
          rhov[k] = rhov[0];
          fputs(" warning: rho values for ixv = ", luout);
          fprintf(luout, "%d set to rho for ixv = 0\n", k);
        } else {
          ++nxrho;
        }
      }
    }
  }
  if (kmax0 > 0) /* integer K sampling */
    dcopy(nxvib, rhov, 1, drhokv, 1);
  else
    dcopy(nxvib, &drhok, 0, drhokv, 1);
  if (npar < 2) {
    puts(" bad parameter file"); exit(EXIT_FAILURE);
  }
  rewind(luin);
  nftp = nft | 1; nsin = nft >> 1;
  for (k = 0; k < nrhok; ++k) { /* offsets for overlaps */
    kk = (k + 1) >> 1; /* i = 0,1,1,2,2 */
    i = kk;
    if (ODD(k)) i = -i;        
    krhok[k] = i;     /* i = 0,-1,1,-2,2 */
    if (nhalf > 0 && ODD(kk))
      ++i;/* i = 0,0,2,-2,2 */
    rhokv[k] = 0.5 * i * rhov[0];
    for (ixv = 1; ixv < nxvib; ++ixv) 
      rhokv[k + ixv * nrhok] = 0.5 * i * rhov[ixv];
  }
  if (maxm < 0) {
    maxm = -maxm; nomorg = TRUE;
  }
  maxm = (maxm - 1) / nsym + 1;
  mbasis = maxm + maxm + 1; nbasis = mbasis;
  nvt = (nvib - 1) / nsym + 1; 
  if (nvt > nbasis)
    nvt = nbasis;
  nvib = nvt * nsym;
  nvtsq = nvt * nvt;
  nbasis *= nxdim;
  nvib0 = nvib;
  if (nvib0 > nvib)
    nvib0 = nvib;
  if (nxvib != 2 || isoffd || nsym > 3)
    nvib2 = 0;
  nnblk = nxdim * nvib;
  nvtot = nvib * nxvib; nvtot0 = nvib0 * nxvib;
  nvibsq = nvib * nvib; nvtotsq = nvtot * nvtot;
  if (nvib2 <= 0) {
    nsig = nsym; nvtot2 = nvtot;
    nvtop = 1;
  } else {
    nsymb = -nsig; nvtot0 = nvib2;
    nvtop = ((nvt + 1) * (nvt + 2)) >> 1;
  }
  nvfac = 10;
  if (nvtot2 > 9)
    nvfac = 100;
  nvfac2 = nvfac * nvfac;
  nsigsq = nsig * nsig;
  nl = nvtot2 * sizeof(int);
  sigma2 = (int *)mallocq(nl);
  nl = (nvib2 + 1) * sizeof(S2TOP);
  vmap = (S2TOP *)mallocq(nl);
  ivtop = (int *)mallocq((size_t) nvtop * sizeof(int));
  nvtop = init2top(luout, nvib2, nsym, nvib, nvtot2, nsig, kmap, 
                   vmap, sigma2, ivtop);
  free(kmap);
  ndtop = 0; ndmtop = 0;
  /* allocate space for eigenvalues */
  negy = nft * nvtot;
  nl = (size_t) (negy + nvtot) * sizeof(double);
  eigval = (double *) mallocq(nl);
  eigval[0] = 0.;
  eigvalk = NULL;
  eigvalk2 = NULL;
  /* allocate space for eigenvector structures */
  neig = nft * nrhok; neigp = neig + 1;
  nl = (size_t) neigp * sizeof(SEIG);
  eigv = (SEIG *) mallocq(nl);
  /* allocate space for sigma */
  nl = (size_t) nvtot * sizeof(int);
  sigma = (int *) mallocq(nl);
  sigma[0] = 0;
  /* allocate space for morg + isym */
  isym = (int *) mallocq(nl);
  isym[0] = 0; 
  for (k = 1; k < nvtot; ++k)
    isym[k] = 0;
  nl *= (size_t) neigp;
  pmorgwk = (int *) mallocq(nl);
  pmorgwk[0] = 0;
  nl = (size_t) nbasis; 
  nl *= (size_t) nvtot * sizeof(double);  /* eigenvector size */
  /* allocate space for eigenvectors */
  seig_alloc(eigv, neig, nl, pmorgwk, nvtot);
  eig = &eigv[0]; kd = 0;
  for (ift = 0; ift < nft; ++ift) {
    k = ift; kk = nft - ift;
    if (k > kk) k -= nft;
    val = (double) k;
    for (kd = 0; kd < nrhok; ++kd) {
      for (ixv = 0; ixv < nxvib; ++ixv) {
        eig->rhok[ixv] = rhokv[kd + nrhok * ixv] + val * drhokv[ixv];
      }
      eig->irhok = kd; 
      ++eig;
    }
  }
  val = 0.5 * nsym;
  dcopy(nxvib, &val, 0, eig->rhok, 1);
  eig->irhok = 0;
  eig = NULL;
  vref = -0.5 + 0.5 * tiny / nsym; 
  for (k = 0; k < nsym; ++k) {  /* set up sigma values */
    kk = k + 1;
    isig = kk >> 1; 
    if (ODD(kk)) {
      isig = -isig;      /* isig = 0,+1,-1,.. */
    } 
    for (iv = k; iv < nvtot; iv += nsym)
      sigma[iv] = isig; 
    eig = &eigv[0];
    for (kk = 0; kk <= neig; ++kk) {
      /* set up morg values for each sigma value */
      ixvv = 0;
      for (ixv = 0; ixv < nxvib; ++ixv) {
        if (isoffd && ixv > 0) {
          kd = kd + morgsym[ixv] - morgsym[ixv - 1];
        } else {
          /* abs(rhok - morg) < 0.5 * nsym          */
          val = (eig->rhok[ixv] - isig) / nsym;
          kd = (int) floor(val + 0.5); 
          val -= kd; /* abs(val) <= 0.5 */
          if (val < vref) {
            --kd; val += 1.;
          }
          if (nomorg) {
            val += kd; kd = 0;
          }
          kd = isig + nsym * kd;
        }
        pmorg = eig->morg;
        for (iv = k; iv < nvib; iv += nsym) 
          pmorg[iv + ixvv] = kd;
        ixvv += nvib;
      }
      ++eig;
    }
  }
  eig = NULL;
  if (nvib2 == 0){
    for (k = 0; k < nvtot; ++k) sigma2[k] = sigma[k];
  }
  /* allocate space for ioffv + ixeig + ipt0 */
  ncalcm = neig * nsym  * nxblks;
  ncalc = ncalcm + nxblks;
  nl = (size_t)ncalc * sizeof(int);
  ioffv = (int *) mallocq(nl);
  ioffv[0] = 0;
  ixeig = (int *) mallocq(nl);
  ixeig[0] = 0;
  ixpt0 = (int *) mallocq(nl);
  ixpt0[0] = 0;
  nl = (size_t)ncalc * sizeof(double);
  rhokdif = (double *) mallocq(nl);
  rhokdif[0] = 0.;
  kk = 0;
  for (isigx = 0; isigx < nsym; ++isigx) {
    eig = &eigv[0];
    for (k = 0; k < neig; ++k) {
      ixvv = 0;
      for (ixv = 0; ixv < nxblks; ++ixv) {
        isig = isigx + ixvv;
        ioffv[kk] = isig;
        ixeig[kk] = k;
        ixpt0[kk] = kk;
        val = fabs(eig->rhok[ixv] - eig->morg[isig]);
        if (val < tiny) val = 0.;
        rhokdif[kk] = val / nsym;
        ++kk; ixvv += nvib;
      }
      ++eig;
    }
  }
  ixvv = 0;
  for (ixv = 0; ixv < nxblks; ++ixv) {
    ioffv[kk] = ixvv;
    ixeig[kk] = neig;
    ixpt0[kk] = kk;
    rhokdif[kk] = 0.5;
    ++kk; ixvv += nvib;
  }
  ixvv = 0;
  for (ixv = 0; ixv < nxblks; ++ixv) { /* sort rhokdif */
    for (i = ixv + nxblks; i < ncalc; i += nxblks) { 
      kk = i; vmin = rhokdif[i];
      for (k = i + nxblks; k < ncalc; k += nxblks) {
        val = rhokdif[k];
        if (val < vmin) {
          vmin = val; kk = k;
        }
      }  
      if (kk > i) {
        isig = ioffv[kk]; ioffv[kk] = ioffv[i]; ioffv[i] = isig; 
        kd   = ixeig[kk]; ixeig[kk] = ixeig[i]; ixeig[i] = kd;
        rhokdif[kk] = rhokdif[i]; rhokdif[i] = vmin;
      }
      kk = ixpt0[i - nxblks];
      val = rhokdif[kk]; 
      if (fabs(val - vmin) > tiny) {
        ixpt0[i] = i;
        continue;
      }
      ixpt0[i] = kk;
      ivp = eigv[ixeig[kk]].irhok;
      iv  = eigv[ixeig[ i]].irhok;
      if (ivp < iv) continue;
      if (ivp == iv && ioffv[kk] <= ioffv[i]) continue;
      isig = ioffv[kk]; ioffv[kk] = ioffv[i]; ioffv[i] = isig; 
      kd   = ixeig[kk]; ixeig[kk] = ixeig[i]; ixeig[i] = kd;
    }
    ixvv += nvib;
  }
  fgetstr(card, 82, luin);
  for (k = 0; k < nline2; ++k)
    fgetstr(card, 82, luin);
  /* allocate structures for parameters */
  nl = (size_t) npar * sizeof(SIPAR *);
  par = (SIPAR **) mallocq(nl);  /* allocate vector of parameter pointers */
  nl = (size_t) npar * sizeof(SIPAR);
  parwk = (SIPAR *) mallocq(nl); /* allocate parameter structures */
  par[0] = parwk;
  spar_alloc(parwk, par, npar);
  nl = (size_t) (ntyp + 1) * sizeof(SXPEC *);
  pxpecv = (SXPEC **) mallocq(nl);      /* allocate vector of expectation pointers */
  nl = (size_t) (ntyp + 1) * sizeof(SXPEC);
  pxpecpwk = (SXPEC *) mallocq(nl);     /* allocate expectation structures */
  pxpecv[0] = pxpecpwk;
  pxpecpwk->xpec = NULL;
  sxpec_alloc(pxpecv, ntyp);
  pxpecp = pxpecv[0];
  pxpecp->mval = maxm;
  pxpecp->nval = nsym;
  pxpecp->kdel = 0;
  ivx = nxvib;
  if (isoffd) ivx = -ivx; 
  pxpecp->ixvib = ivx;
  ntyp0 = ntyp;
  pxpecp->isunit = ntyp;
  pxpecp = NULL;
  par[0]->idval = 0; idtmp = 0; kdel0 = -1; 
  for (ipar = 0; ipar < npar; ++ipar) { /* read parameters */
    if (fgetstr(card, 82, luin) <= 0)
      break;
    dval[4] = 1.e+36;
    if (pcard(card, dval, 5, NULL) < 4) {
      --ipar; continue;
    }
    idtmp = (long) dval[0]; idtst = idtmp;
    nv = (int) dval[1];
    mv = (int) dval[2];
    if (idtst < 0) idtst = -idtst;
    if (idtst < 9L && mv < 100) {
      mv += 100; idtst += 990;
    } else if (idtst > 800000L && idtst <= 899099L) {
      k = (int) (idtst / 1000L) - 800;
      if (k > kmax) kmax = k;
      idtmp = 800099L;
    } else if (idtst == 910099L) {
      --ipar; continue;
    }
    if (mv >= 100) {
      ivv = (int)(idtst / 10); idtst -= ivv * 10;
    } else {
      ivv = (int)(idtst % 100); idtst += 99 - ivv;
    }
    ivp = ivv / 10; iv = ivv - 10 * ivp;
    if (ivp > iv) {
      ivv = ivp; ivp = iv; iv = ivv;
      ivv = iv + ivp * 10;
    }
    ivvp = ivp + iv * 10;
    if (iv >= 9 && ivv != 99) {
      --ipar; continue;
    }
    if (idtmp > 0) {
      idtmp = idtst;
    } else {
      idtmp = -idtst;
    }
    parp = par[ipar];
    parp->offset = 0;
    parp->parval = dval[3]; 
    parp->parerr = dval[4];
    parp->flg = 0;
    if (idtst == 99L && mv == 0 && nv > 0 && nv < 25 && ivp == iv) {
      parp->parval *= -0.5;
      parp->flg |= MPE;
    }
    if (parp->parerr < -eps) {
      parp->parerr *= -clight;
      parp->parval *= clight;
    }
    pstr = strchr(card, '/');
    if (pstr == NULL) {
      parp->labl[0] = '\0';
    } else {
      memcpy(parp->labl, pstr + 1, NLABL);
      parp->labl[NLABL] = '\0';
    }
    ivp = ivv;
    if (ivv == 99) 
      ivp = 0;
    parp->idval = idtmp;
    ntyp = getparv(nv, mv, ivp, nvib2, parp, pxpecv);
    kdel = ideqv(&idtmp, &val); 
    k = mv;
    if (nv < 0) ++k;
    if (nvib2 > 0 && nv <= -25) ++k;
    parp->ksym = k;
    if (ivv == 99) {
      for (ixv = 1; ixv < nxvib; ++ixv) {
        parr = par[++ipar];
        memcpy(parr, parp, sizeof(SIPAR));
        if (parr->idval > 0) parr->idval = -parr->idval;
        ivp = ixv * 11;
        ntyp = getparv(nv, mv, ivp, nvib2, parr, pxpecv);
      }
    } else if (ivv != ivvp) {
      parr = par[++ipar];
      memcpy(parr, parp, sizeof(SIPAR));
      if (parr->idval > 0) parr->idval = -parr->idval;
      ntyp = getparv(nv, mv, ivvp, nvib2, parr, pxpecv);
    }
    if (kdel != 0) {
      /* add parameters for complementary operators */
      parq = par[++ipar];
      memcpy(parq, parp, sizeof(SIPAR));
      parq->parval *= val;
      parq->parerr *= fabs(val);
      parq->idval = idtmp;
      parq->offset = 1;
      parq->ksym = k + 1;
      if (ivv == 99) {
        for (ixv = 1; ixv < nxvib; ++ixv) {
          parr = par[++ipar];
          memcpy(parr, parq, sizeof(SIPAR));
          if (parq->idval > 0) parq->idval = -parq->idval;
          ivp = ixv * 11;
          ntyp = getparv(nv, mv, ivp, nvib2, parr, pxpecv);
        }
      } else if (ivv != ivvp) {
        parr = par[++ipar];
        memcpy(parr, parq, sizeof(SIPAR));
        if (parr->idval > 0) parr->idval = -parr->idval;
        ntyp = getparv(nv, mv, ivvp, nvib2, parr, pxpecv);
      }
    }
    if (parp->idval < 0 && kdel != kdel0) 
      parp->idval = -parp->idval;
    kdel0 = kdel;
  }
  fclose(luin);
  if (kmax >= 0) {
    nl = (size_t)(kmax + 1) * sizeof(double);
    eigvalk = (double *) mallocq((size_t) nvtot * nl);
    eigvalk[0] = 0.;
    eigvalk2 = (double *) mallocq((size_t) nvtot2 * nl);
    eigvalk2[0] = 0.;
  }
  pxpecp = pxpecv[0];
  pxpecp->nval = nsym;
  parp = parq = NULL;
  /* allocate space for expectation matrices */
  nxpec = nft * nvtotsq;
  for (ityp = 1; ityp <= ntyp; ++ityp) {
    pxpecp = pxpecv[ityp];
    nv = nxpec;
    kdel = pxpecp->kdel;
    if (kdel > 0) {
      nv = nv << 1;
    } else if (kdel < 0) {
      pxpecp->kdel = 0;
    }
    if (pxpecp->xpec == NULL)
      pxpecp->xpec = (double *) mallocq((size_t)nv * sizeof(double));
    dcopy(nv, &zero, 0, pxpecp->xpec, 1);
  }
  pxpecp = NULL;
  nl = (size_t) (nsigsq * nftp + nftp + 1) * sizeof(double);
  xfc = (double *) mallocq(nl);
  nl = (size_t) nsigsq * sizeof(double);
  didtmpv = (double *) mallocq(nl);
  didtmpv[0] = 0.;
  fputs("ipar,   idpar, vv,  n,  m, kd, ip,  par\n", luout);
  if (par[0]->idval < 0)
    par[0]->idval = -par[0]->idval;
  for (ipar = 0; ipar < npar; ++ipar) { /*   print out parameters */
    parp = par[ipar];
    idtyp = parp->idval;
    pxpecp = parp->ptxpec;
    if ((parp->flg & MREV) != 0)
      pxpecp = parp->ptxpecx;
    if (pxpecp == NULL) continue;
    kdel = pxpecp->kdel;
    val = parp->parval;
    if ((parp->flg & MPE) != 0) {
      val *= -2.;
    }
    cval = ((parp->flg & MHAM) == 0)? ' ':'*';
    ixv = pxpecp->ixvib;
    fprintf(luout, "%4d %8ld %3d %3d %3d %3d %3d%c %15.5E %15.5E %s\n",
            ipar, idtyp, ixv, pxpecp->nval, pxpecp->mval, kdel,
            parp->offset, cval, val, parp->parerr, parp->labl);
    if (parp->offset != 0)
      parp->offset = nxpec;
    pxpecp = parp->ptxpec;
    if (pxpecp != NULL) {
      if (pxpecp->nval >= 25) {
        pxpecp->nval -= 25; parp->flg |= MV12;
      } else if (pxpecp->nval <= -25) {
        pxpecp->nval += 25; parp->flg |= MV12;
      }
    }
    pxpecp = parp->ptxpecx;
    if (pxpecp != NULL) {
      if (pxpecp->nval >= 25) {
        pxpecp->nval -= 25; parp->flg |= MV12;
      } else if (pxpecp->nval <= -25) {
        pxpecp->nval += 25; parp->flg |= MV12;
      }
    }
    if (idtyp > 0)
      idtyp = -idtyp;
    nkk = ipar + 1;
    for (iend = nkk; iend < npar; ++iend) {     /* sort related operators */
      parq = par[iend];
      idtmp = parq->idval;
      if (idtmp > 0)
        break;
      if (idtmp == idtyp) {
        for (k = iend; k > nkk; --k) {
          par[k] = par[k - 1];
        }
        par[nkk] = parq;
        break;
      }
    }
  }
  pxpecp = NULL;
  parp = parq = NULL;
  if (nxrho == 0) {
    fprintf(luout, " rho = %15.10f\n", rho);
  } else {
    for (i = 0; i < nxvib; ++i) 
      fprintf(luout, " rho (%1d) = %15.10f\n", i, rhov[i]);
  }
  egywk = &eigval[negy]; egy = egyp = egywk; eig = eigv;
  eigvec0 = eig->eigvec;
  for (ixv = 0; ixv < nxblks; ++ixv) {
    eigref[ixv] = NULL;
  }
  for (i = 0; i < ncalc; ++i) {
    /* calculate eigenvectors in order of rho K*/
    isig = ioffv[i]; 
    ixv = isig / nvib; isigx = isig - ixv * nvib;
    kk = ixeig[i]; 
    eig = &eigv[kk]; m = eig->morg[isig]; 
    val = eig->rhok[ixv] - m;
    if (val < 0.)
      rhokdif[i] = -rhokdif[i];
    lu = NULL; egy  = egywk;
    k = eig->irhok; 
    if (k == 0) {
      ift = kk / nrhok; 
      if (ift < nft)
        egy = &eigval[nvtot * ift];
      if (isigx == 0)
        lu = luout;
    }
    eigvec0 = eig->eigvec;
    iv = ixpt0[i]; 
    if (iv == i) {
      printf(" rhoK = %9.6f, %2d %2d \n", val, krhok[k], sigma[isigx]);
      egyp = egy;
      fflush(stdout);
      eig->rhok[ixv] = fabs(val); eig->morg[isig] = 0;
      iamdiag(lu, nvib, pxpecv, eig, egy, isig);
      eig->rhok[ixv] = val + m; eig->morg[isig] = m;
      eigvec1 = eigref[ixv];
      for (iv = 0; iv < nnblk; iv += nsym) {
        ivv = (iv + isig) * nbasis; ivvp = iv * nbasis;
        if (eigvec1 != NULL) {
          val = ddot(nbasis, &eigvec1[ivvp], 1, &eigvec0[ivv], 1);
          if (val < 0.) {
            dneg(nbasis, &eigvec0[ivv], 1);
          }
        }
      }
      eigref[ixv] = &eigvec0[isig * nbasis];
    } else {
      isigp = ioffv[iv];
      kk = ixeig[iv]; eigp = &eigv[kk];
      if (eigp->irhok == 0) {
        for (iv = 0; iv < nnblk; iv += nsym) {
          egy[iv + isig] = egyp[iv + isigp];
        }
      }
      eigvec1 = eigp->eigvec;
      for (iv = 0; iv < nnblk; iv += nsym) {
        ivv = (iv + isig) * nbasis; ivvp = (iv + isigp) * nbasis;
        dcopy(nbasis, &eigvec1[ivvp], 1, &eigvec0[ivv], 1);
      }
    }
  }
  ivv = iv = 0;
  hsym = FALSE; vmax = 0.; nneg = 0;
  eigp = &eigv[neig]; eigvec1p = eigp->eigvec;
  eig  = &eigv[0];    eigvec0p =  eig->eigvec;
  fputs("  v, sym, symmetry at 0 and mid rhoK\n", luout);
  for (ivp = ivvp = 0; ivp < nvtot; ++ivp, ivvp += nbasis) {
    if (sigma[ivp] == 0) {
      iv = ivp;
      /* determine symmetry for rhoK = 0.5 * nsym and 0 */
      eigvec0 = eigvec0p + ivvp; 
      eigvec1 = eigvec1p + ivvp;
      ivv = 0;
      valm[0] = 0.; valm[1] = 0.; vmax = 0;
      for (ixv = 0; ixv < nxdim; ++ixv) {
        if (ixv > 0) {
          eigvec0 += mbasis; eigvec1 += mbasis; ivv += nvib;
        }
        mv = -(eig->morg[ivv] << 1) / nsym; eigvec = eigvec0;
        for (kk = 0; kk <= 1; ++kk) {
          if (ODD(mv)) {
            mv = maxm + ((mv + 1) >> 1); mm = mv - 1;
            k = (mv > maxm)? maxm: mv;
            vref = 0.;
            for (m = 1; m <= k; ++m) { 
              /* mv - 1, mv */
              vref += eigvec[mv - m] * eigvec[mm + m];
            }
            nneg = k;
          } else {
            mv = maxm + (mv >> 1);
            k = (mv > maxm)? maxm: mv;
            vref = 0.5 * eigvec[mv] * eigvec[mv];
            for (m = 1; m <= k; ++m) {
              vref += eigvec[mv - m] * eigvec[mv + m];
            }
            nneg = k + 1;
          }
          valm[kk] += vref;
          if (kk != 0) break;
          nl = idamax(nneg, &eigvec[mv], 1) + mv;
          vref = eigvec[nl];
          if (fabs(vref) > fabs(vmax)) vmax = vref;
          mv = 1 - (eigp->morg[ivv] << 1) / nsym; eigvec = eigvec1;
        }
      }
      nneg = 0;
      if (vmax < 0.) nneg = 1;
      val = 2.* valm[0]; vref = 2. * valm[1]; 
      isym[iv] = 0;
      if (val < 0.)
        isym[iv] = 3;
      if (vref < 0.) 
        isym[iv] ^= 1;
      fprintf(luout, "%3d %3d  %9.6f %9.6f\n", iv, isym[iv] & 3, val, vref);
    } else {
      if (sigma[ivp] < 0) {
        isym[ivp - 1] = isym[iv] + 4;
        isym[ivp] = isym[iv] + 8;
      } else {
        isym[ivp] = isym[iv] + 12;
      }
    }
    if (ODD(isym[ivp])) {      /* correct for translation symmetry */
      hsym = TRUE;
      eig = &eigv[0];
      for (k = 0; k < neig; ++k) {
        kk = eig->morg[ivp] / nsym + nneg;
        if (ODD(kk)) {
          eigvec1 = eig->eigvec + ivvp;
          dneg(nbasis, eigvec1, 1);
        }
        ++eig;
      }
    } else if (nneg != 0) {
      eig = &eigv[0];
      for (k = 0; k < neig; ++k) {
        eigvec1 = eig->eigvec + ivvp;
        dneg(nbasis, eigvec1, 1);
        ++eig;
      }
    }
  }
  for (i = 1; i < ncalc; ++i) { /* fix negative rhoK */
    if (rhokdif[i] >= 0.) continue;
    isig = ioffv[i]; 
    ixv = isig / nvib; isigx = isig - ixv * nvib;
    kk   = ixeig[i]; 
    eig  = &eigv[kk]; 
    eigvec0 = eig->eigvec; 
    for (ii = 0; ii < nnblk; ii += nsym) {
      iv = ii + isig;
      eigvec = &eigvec0[iv * nbasis];
      if ((isym[iv] & 2) != 0) {
        dneg(nbasis, eigvec, 1);
      }
      ivv = isigx;
      for (ixv = 0; ixv < nxdim; ++ixv) {
        if (ixv > 0) {
          eigvec += mbasis;
        }
        if (morgsym[ixv] != 0) {
          mv = maxm; mm = mv - 1;
          for (m = 1; m <= mv; ++m) { 
            /* mv - 1, mv */
            vref = eigvec[mv - m];
            eigvec[mv - m] = eigvec[mm + m];
            eigvec[mm + m] = vref;
          }
        } else {
          mv = maxm;
          for (m = 1; m <= mv; ++m) {
            vref = eigvec[mv - m];
            eigvec[mv - m] = eigvec[mv + m];
            eigvec[mv + m] = vref;
          }
        }
      }
    }
  }
  egyp = NULL;
  eigvec1 = NULL; egy = NULL;
  if (eigvalk != NULL) { /* need energy vs K */
    egy = eigvalk; 
    kk = (kmax0 + 1) * nvtot;
    dcopy(kk, eigval, 1, egy, 1);
    if (kmax0 < kmax) {
      egy += kk;
      eig = &eigv[neig]; eig->irhok = 0; pmorg = eig->morg;
      for (kk = kmax0 + 1; kk <= kmax; ++kk) {
        for (ixv = 0; ixv < nxblks; ++ixv) {
          isig = ixv * nvib;    
          for (isigx = 0; isigx < nsym; ++isigx) {
            val = kk * rhov[ixv]; eig->rhok[ixv] = val; 
            val = (val - sigma[isigx]) / nsym; 
            pmorg[isig] = sigma[isigx] + nsym * (int)floor(val + 0.5);
            iamdiag(luout, nvib, pxpecv, eig, egy, isig); 
            ++isig;
          }
        }
        printf(" K = %3d, egy = %12.5e\n", kk, egy[0]);
        egy += nvib;
      }
    }
    if (eigvalk2 != NULL) {
      if (nvib2 > 0) {
        egy = eigvalk; egyp = eigvalk2; 
        for (kk = 0; kk <= kmax; ++kk) {
          for (k = 0; k < nvtot2; ++k) {
            iv = vmap[k].ix1; ivp = vmap[k].ix2;
            egyp[k] = egy[iv] + egy[ivp];
          }
          egy += nvtot; egyp += nvtot2;
        } 
        egy = NULL; egyp = NULL;
      } else {
        kk = nvtot * (kmax + 1);
        dcopy(kk, eigvalk, 1, eigvalk2, 1);
      }
    }
    free(eigvalk); eigvalk = NULL;
    if (eigvalk2 != NULL) {
      egy = eigvalk2;
      for (kk = 0; kk <= kmax; ++kk) {
        for (isig = 0; isig < nsig; ++isig) {
          if (sigma2[isig] < 0) {
            for (k = isig; k < nvtot2; k += nsym) {
              egy[k] = 0.5 * (egy[k - 1] - egy[k]);
              egy[k - 1] -= egy[k];
            }
          }
        }
        egy += nvtot2;
      }
      egy = NULL;
    }
  }
  egywk = NULL;
  iamdiag(NULL, 0, pxpecv, eigv, eigval, -1); /* free store */
#if MYDEBUG
  /*  test for good overlap */
  for (i = 1; i < ncalc; ++i) { /* sort rhokdif */
    kk = i; vmin = rhokdif[i];
    for (k = i + 1; k < ncalc; ++k) {
      val = rhokdif[k];
      if (val < vmin) {
        vmin = val; kk = k;
      }
    }  
    if (kk > i) {
      isig = ioffv[kk]; ioffv[kk] = ioffv[i]; ioffv[i] = isig; 
      kd = ixeig[kk]; ixeig[kk] = ixeig[i]; ixeig[i] = kd;
      rhokdif[kk] = rhokdif[i]; rhokdif[i] = vmin;
    }
  }
  eig = &eigv[0]; isig = 0;
  for (i = 1; i <= ncalc; ++i) {
    isigp = isig; eigp = eig;
    if (i < ncalc) {
      isig = ioffv[i];
      kk = ixeig[i]; val = rhokdif[i];
      if (kk == neig) continue;
    } else {
      isig = kk = 0; val = 0.;
    }
    eig = &eigv[kk];
    eigvec1 = eig->eigvec;
    eigvec0 = eigp->eigvec;
    for (iv = 0; iv < nvtot; iv += nsym) {
      ivv  = nbasis * (iv + isig);
      ivvp = nbasis * (iv + isigp);
      vref = ddot(nbasis, &eigvec0[ivvp], 1, &eigvec1[ivv], 1);
      if (ODD(isym[iv])) {
        kk = (eig->morg[iv] - eigp->morg[iv]) / nsym;
        if (ODD(kk))
          vref = -vref;
      }
      fprintf(luout, " rhoK = %9.6f,v = %2d, %10.6f\n", val, iv + isig, vref);
    }
  }
  eigvec1 = NULL;
#endif                          /* end test */
  for (ixv = 0; ixv < nxblks; ++ixv) eigref[ixv] = NULL;
  eigvec0 = NULL;
  if (nrhok == 1) {
    nft = 1; nftp = 1;
  }
  eig0 = &eigv[0];
  for (kk = 0; kk < nft; ++kk) {        /* calculate expectation values */
    eig = eigp = eig0;
    kd = 0;
    nkk = kk * nvtotsq;
    for (ityp = 1; ityp <= ntyp; ++ityp) {
      k = pxpecv[ityp]->kdel;
      if (kd < k) {
        kd = k; k = k << 1;
        if (k < nrhok) {
          eigp = eig0 + k; eig = eigp - 1;
        } else {
          eig = eigp = eig0;
        }
      }
      expec(ityp, nkk, pxpecv, nvib, eig, eigp, nxpec, isym);
    }
    eig0 += nrhok;
  }
  eig0 = NULL; eigp = NULL;
  free(ioffv); free(ixeig); free(rhokdif); free(ixpt0);
  for (i = 0; i <= neig; ++i) {
    eig = &eigv[i];
    if (eig->eigvec != NULL) {
      free(eig->eigvec);
      eig->eigvec = NULL;
    }
  }
  eig = eigp = NULL;
  free(pmorgwk);
  free(eigv);
  lupar = fopen(cfil[epar], "r");
  lubak = fopenq(cfil[ebak], "w");
  headok = FALSE;
  if (lupar != NULL) {
    if (fgetstr(card, 82, lupar) > 0) {
      fputs(card, lubak); fputc('\n', lubak);
      if (fgetstr(card, 82, lupar) > 0) {
        fputs(card, lubak); fputc('\n', lubak);
        mm = 0;
        for (;;) {/* read option cards */
          if (fgetstr(card, 82, lupar) < 0) break;
          fputs(card, lubak); fputc('\n', lubak); 
          if (readopt(card, dval, 8) <= 0) break;
          m = (int) dval[0];
          if (mm == 0) mm = m;
          if (dval[7] > -0.5) break;
        }
        headok = TRUE;
        if (nvtot2 != mm) 
          printf("%s (%3d) does not match .iam (%3d)\n",
                 "WARNING: # states in .par file", mm, nvtot2); 
      }
    }
    fclose(lupar);
  }
  luint = NULL;
  if (ndip >  0) 
    luint = fopen(cfil[eint], "r");
  if (luint != NULL){ 
    if (fgetstr(title, 82, luint) > 0 &&
        fgetstr(card, 82, luint) > 0) {
      fclose(luint);
      luint = fopen(cfil[eint], "w");
      if (luint != NULL) {
        fputs(title, luint); fputc('\n', luint); 
        fputs(card, luint); fputc('\n', luint);
      }
    } else if (luint != NULL) {
      fclose(luint); luint = NULL;
    }
  }  
  if (luint == NULL && ndip > 0) 
    luint = fopenq(cfil[eint], "w");
  nnpar = 0;
  nl = (size_t) (2 * nft) * sizeof(double);
  txpec = (double *) mallocq(nl);
  txpec[0] = 0.;
  /*   print out rho */
  if (kmax0 > 0 && !nofc) {
    val = 3. / nft; drhok = 1.;
  } else {
    val = fabs(3. * rho / nsym);
  }
  rhoang = val * acos(-0.5); /* 2*pi/3 */
  pstr = " 0. / rho";
  if (kmax0 > 0 && hsym) {
    hsym = FALSE;
    for (iv = 0; iv < nvtot; ++iv) isym[iv] &= 0x7ffe;
  }
  if (hsym && !nofc) 
    nxrho = 1;
  if (nxrho > 0) {
    ivv = 0;
    for (ixv = 0; ixv < nxvib; ++ixv) {
      val = rhov[ixv];
      for (iv = 0; iv < nvib; ++iv) {
        txpec[0] = val;
        if (ODD(isym[ivv]))
          txpec[0] = -txpec[0];
        didtmp = 9100. * (double) nvfac2 + ivv * (nvfac + 1);
        nnpar += wrpar(lubak, &didtmp, 1, 0, kflg0, txpec, 1, pstr, 0, drhok);
        ++ivv;
      }
    }
  } else {
    idtyp = 910099L;
    if (nvfac == 100)
      idtyp = idtyp * nvfac + 99;
    txpec[0] = val;
    didtmp = (double) idtyp;
    nnpar += wrpar(lubak, &didtmp, 1, 0, kflg0, txpec, 1, pstr, 0, drhok);
  }
  /* print out Fourier transform of energies */
  txpec[0] = par[0]->parval;
  idtyp = 0; didtmp = 0.;
  fputs("Fourier transform of energies\n", luout);
  wrpar(luout, &didtmp, 1, 0, kflg0, txpec, 1, "", 0, drhok);
  idtyp = 99;
  for (iv = 0; iv < nvtot2; ++iv) {
    if (nvib2 > 0) {
      ivx = vmap[iv].ix1; ivxp = vmap[iv].ix2 + nvib;
      for (k = 0; k < nft; ++k) {
        txpec[k] = eigval[ivx] + eigval[ivxp];
        ivx += nvtot; ivxp += nvtot;
      }
    } else {
      dcopy(nft, &eigval[iv], nvtot, txpec, 1);
    }
    ftanal(0, nft, txpec, xfc);
    didtmp = iv * (double)(nvfac + 1);
    wrpar(luout, &didtmp, 1, 0, kflgm, xfc, nftp, " Egy", nvfac2, drhok);
  }
  fputs("expectation values\n", luout);
  pstr = NULL; pbgn = NULL;
  rezero = !isoffd;
  if (nvib2 > 0) {
    ndmtop = nsymsq * (nvtop >> 1);
    vmapsym = (int *)mallocq((size_t) ndmtop * sizeof(int));
    vmapsym[0] = 0;
    for (k = 1; k < ndmtop; ++k) vmapsym[k] = 0;
    for (k = 0; k < nvtot2; ++k) {
      if (vmap[k].iqv < 0) continue;
      i = vmap[k].inew;
      vmapsym[i] = (vmap[k].iqsig & 3) + (vmap[k].iqv << 2);
    }     
    nvtotsq = ndmtop * ndmtop; nxpec = nvtotsq * nft;
  } else {
    vmapsym = (int *) mallocq(sizeof(int)); vmapsym[0] = 0;
  }
  xpec2top = (double *) mallocq(sizeof(double)); 
  xpec2top[0] = 0.;
  /* assemble expectation values */
  nprn = nftp; 
  strcpy(card, " 0.0  / egy");
  pbgn = strchr(card, '/') - 1; pstr = pbgn + 2;
  for (ipar = 0; ipar < npar; ipar = iend + 1) {
    zval = 0.; amp = 0.; kflgm = -1;
    parp = par[ipar];
    idtyp = parp->idval;
    isdip = 0; lu = lubak;
    if (idtyp < 0) {
      idtyp = -idtyp;
      if (idtyp < 99) {
        isdip = 1; lu = luint;
      }
    } else {
      if (idtyp < 99) {
        kflgm = -2; isdip = 1; lu = luint;
      }
      sprintf(card,"%12.4E  /%s", parp->parerr, parp->labl);
      pbgn = strchr(card, '/') - 1; pstr = pbgn + 2;
      txpec[0] = parp->parval;
      pxpecp = parp->ptxpec;
      if (pxpecp == NULL) break; 
      if (isdip > 0) {
        didtmp = 9900.; mv = 99;
        if (nvfac == 100)
          didtmp = didtmp * 100 + 99.;
      } else {
        if ((parp->flg & MPE) != 0)
          txpec[0] *= -2.;
        didtmp = (double)(idtyp / 100);
        nv = pxpecp->nval; mv = pxpecp->mval;
        if (nv != 0) {
          if (nv < 0) {
            nv = -nv;
            if ((parp->flg & MV12) != 0) nv += 25;
            i = (nv + 9) / 10;
          } else {
            if ((parp->flg & MV12) != 0) nv += 25;            
            i = (nv - 1) / 10;
          }
          nv += i * 10;
          didtmp += nv * 1.e+07;
        }
        if (nvfac == 100) {
          didtmp = didtmp * 100 + mv; mv = 99;
        } else {
          mv = mv * 10 + 9;
        }
      }
      nnpar += wrpar(lu, &didtmp, 1, isdip, kflg0, txpec, 1, card,
                     mv, drhok);
    }
    ksymtau = parp->ksym << 1;
    idtst = -idtyp;
    if (idtyp == 800099L) idtst = -99L;
    iend = ipar;
    for (i = ipar + 1; i < npar && idtst == par[i]->idval; ++i)
      iend = i;
    if (idtyp == 990099L)
      continue;
    if (idtyp == 800099L) {
      if (eigvalk2 == NULL) continue;
      for (iv = 0; iv < nvtot2; ++iv) {       
        val = 0.;
        egy = &eigvalk2[iv];
        if (sigma2[iv] >= 0 && iv > 0) {
          val = egy[0]; txpec[0] = val;
          didtmp = (double)(-iv * (nvfac + 1));
          nnpar += wrpar(lu, &didtmp, 1, 0, kflg0, txpec, 1, pbgn, 0, drhok);
        }
        for (kk = 1; kk <= kmax; ++kk) {
          egy += nvtot2;
          txpec[0] = (*egy) - val;
          didtmp = ((8000. + kk * 10.) * nvfac + iv) * nvfac + iv;
          didtmp = -didtmp;
          nnpar += wrpar(lu, &didtmp, 1, 0, kflg0, txpec, 1, pbgn, 0, drhok);
        }
      }
      continue;
    }
    pxpecp = parp->ptxpec; kdel = -2;
    if (pxpecp != NULL) kdel = pxpecp->kdel; 
    if (nvib2 > 0) {
      k = iend - ipar + 1;
      if (kdel > 0) ++k;
      if (k > ndtop) {
        free(xpec2top);
        ndtop = k;
        nl = (size_t) nxpec; 
        nl = nl * ndtop * sizeof(double);
        xpec2top = (double *)mallocq(nl);
      }
      pdtmp = xpec2top; xpec2top[0] = 0.;
      for (k = ipar; k <= iend; ++k) {
        if (k > ipar) pdtmp += nxpec;
        parq = par[k];
        get2top(parq, pdtmp, nvtot, nft, nsym, vmapsym, ivtop, nvtop);
      }
      pdtmp = NULL; 
    }
    voff = vref = 0.;
    for (ivx = 0; ivx < nvtot2; ivx += nsig) {
      for (ivxp = 0; ivxp <= ivx; ivxp += nsig) {
        if (nvib2 > 0) {
          isymtau = ksymtau ^ isym[vmap[ivx].ix1] ^ isym[vmap[ivxp].ix1] ^ 
            isym[vmap[ivx].ix2] ^ isym[vmap[ivxp].ix2];
        } else {
          isymtau = ksymtau ^ isym[ivx] ^ isym[ivxp];
        }
        pxfc = pxfc1 = xfc;
        xfc[0] = 0.;
        nvs = 0;
        nv = ivx + nsig - 1;
        for (iv = ivx; iv <= nv; ++iv) {
          nvp = ivxp + nsig - 1;
          if (nvp > iv) nvp = iv;
          for (ivp = ivxp; ivp <= nvp; ++ivp) {
            /* loops over sigma within torsion state */
            isig = sigma2[iv]; isigp = sigma2[ivp];            
            val = 0.;
            if (iv == ivp)
              val = voff;
            dcopy(nft, &val, 0, txpec, 1);
            for (k = ipar; k <= iend; ++k) {
              parp = par[k];
              pxpecp = parp->ptxpec;
              if (pxpecp == NULL) continue;
              val = parp->parval;
              idtmp = parp->idval;
              if (nvib2 > 0) {
                nl = (size_t) ndmtop * vmap[iv].inew + vmap[ivp].inew;
                nl += nxpec * (k - ipar);
                daxpy(nft, val, &xpec2top[nl], nvtotsq, txpec, 1);
              } else {
                pdtmp = pxpecp->xpec;
                if (pdtmp == NULL) continue;
                pdtmp += (iv + ivp * nvtot + parp->offset);
                daxpy(nft, val, pdtmp, nvtotsq, txpec, 1);
              }
            }
            pxpecp = NULL; pdtmp = NULL; 
            kdoff = 0;
            if (nofc) {
              nprn = ksign(kdel, nft, txpec, pxfc);
              nft2 = nft;
            } else {
              if (nhalf > 0 && ODD(kdel)) {
                kdoff = 1;
                txpec[nhalf] = 0.5 * (txpec[nhalf - 1] + txpec[nhalf + 1]);
              }
              nft2 = ftanal(isymtau, nft, txpec, pxfc);
              if (kdoff != 0) {
                ftfix(nft2, nftp, pxfc);
              }
            }
            if (negrho){
              dneg(nsin, &pxfc[2], 2);
            }
            parp = NULL;
            if (isig  == nsymb) isig = 0;
            if (isigp == nsymb) isigp = 0;
            ldel = isigp - isig;
            idtmp = ivp * nvfac + iv; 
            if (ldel < 0)
              idtmp = iv * nvfac + ivp;
            didtmp = (double) idtmp;
            if (isdip > 0) {
              didtmp = didtmp * 10 + (double)(idtyp % 10);
            } else if (idtyp > 99) {
              didtmp += nvfac2 * (double)(idtyp / 100);
            } else if (idtmp == 0) {    /* idtyp == 99 */
              voff = -txpec[0];
              xfc[0] += voff;
              for (k = 0; k < nft; ++k) {
                txpec[k] += voff;
              }
              vref = xfc[0];
              didtmp = -99.;
              if (nvfac == 100)
                didtmp = -9999.;
              nnpar += wrpar(lu, &didtmp, 1, 0, kflgm, xfc, 1, pbgn, 
                             nvfac2, drhok);
              didtmp = 0;
            }
            wrpar(luout, &didtmp, 1, isdip, kflg, txpec, nft, pbgn, 
                  nft, drhok);
            if (iv == ivp)
              pxfc[0] -= vref;
            if (ktest > 0 && nft2 == nft && dnrm2(nft, pxfc, 1) > tiny) {
              fputs("begin K prediction\n",luout);
              egy = txpec + nft;
              for (k = 0; k < ktest; ++k){
                kk = k;
                if ((k + k) >= ktest) kk -= ktest;
                ang0 = (double)kk;
                if (ODD(kdel)) ang0 += 0.5;
                ang0 *= rhoang; ang = ang0; 
                val = pxfc[0];
                for (kk = 1; kk < nft; kk += 2) {
                  val += pxfc[kk] * cos(ang) + pxfc[kk + 1] * sin(ang);
                  ang += ang0;
                }
                egy[k] = val;
              }
              pbgn[0] = '/'; pbgn[1] = ':';
              nnpar += wrpar(luout, &didtmp, 1, isdip, kflg, egy, ktest, 
                             pbgn, ktest, 1.);
              pbgn[0] = ' '; pbgn[1] = '/';
            }
            if (lsym && ldel != 0) continue;  
            didtmpv[nvs] = -didtmp;
            if (ldel == 0) {
              if (isig >= 0) { 
                pxfc1 = pxfc;
              } else {
                for (k = 0; k < nprn; ++k) {
                  pxfc[k] = 0.5 * (pxfc1[k] - pxfc[k]);
                  pxfc1[k] -= pxfc[k];
                }
                if (nvib2 == 0) {
                  k = isymtau & 2;
                  val = dnrm2(nprn, pxfc, 1) + tiny;
                  if (val > amp) amp = val;
                  val = 0;
                  for (; k < nprn; k += 2) {
                    val += fabs(pxfc[k]); 
                    if (rezero) 
                      pxfc[k] = 0.;
                    if (k == 0) --k;
                  }
                  val /= amp;
                  if (zval < val) zval = val;
                }
              }
              if (isig <= 0 && nvib2 == 0) {
                k = 2 - (isymtau & 2);
                val = dnrm2(nprn, pxfc, 1) + tiny;
                if (val > amp) amp = val;
                val = 0; 
                for (; k < nprn; k += 2) {
                  val += fabs(pxfc1[k]);
                  if (rezero) pxfc1[k] = 0.;
                  if (k == 0) --k;
                }
                val /= amp;
                if (zval < val) zval = val;
              }
            }
            if (iv < nvtot0) {
              ++nvs;
              pxfc += nprn;
            }
          }                     /*ivp loop */
        }                       /* iv loop */
        nnpar += wrpar(lu, didtmpv, nvs, isdip, kflgm, xfc, nprn, 
                       pbgn, -nvfac2, drhok);
        pxfc1 = NULL; pxfc = NULL; kflgm = -1;
      }                         /*ivxp loop */
    }                           /* ivx loop */
    if (nvib2 == 0)
      fprintf(luout,"%3d zero test = %e\n", ipar, zval);
  }                             /* par loop */
  if (luint != NULL) fclose(luint);
  fclose(lubak);
  fclose(luout);
  lubak = fopen(cfil[ebak], "r");
  if (lubak == NULL) exit(EXIT_FAILURE);
  lupar = fopen(cfil[epar], "w");
  if (lupar == NULL) {
    fclose(lubak);
    exit(EXIT_FAILURE);
  }
  if (fgetstr(card, 82, lubak) > 0) {
    fputs(card, lupar); fputc('\n', lupar); /* title */
    if (fgetstr(card, 82, lubak) > 0) {
      k = readint(card);
      if (headok) {
        fprintf(lupar,"%8ld%s\n", nnpar, &card[k]);
      } else {
        fputs(card, lupar); fputc('\n', lupar); 
      }
      for (;;) {/* read */
        if (fgetstr(card, 82, lubak) < 0) break;
        fputs(card, lupar); fputc('\n', lupar); 
      }
    }
  }
  fclose(lubak);
  fclose(lupar);
  parp = NULL;
  free(eigval);
  free(eigvalk2);
  free(didtmpv);
  free(xfc);
  ityp = 0;
  for (ityp = 1; ityp <= ntyp; ++ityp) {
    pxpecp = pxpecv[ityp];
    if (pxpecp->xpec != NULL)
      free(pxpecp->xpec);
    pxpecp->xpec = NULL;
  }
  pxpecp = NULL;
  free(pxpecv);
  pxpecv = NULL;
  free(pxpecpwk);
  pxpecpwk = NULL;
  free(xpec2top);
  free(vmapsym);
  free(ivtop);
  vmapsym = NULL;
  free(par);
  par = NULL;
  free(vmap); 
  free(sigma2);
  free(isym);
  free(sigma);
  return 0;
} /* main */

int iamdiag(lu, nvib, pxpecv, eig, egy, ioff)
     FILE *lu;
     int nvib;
     SXPEC **pxpecv;
     SEIG *eig;
     double *egy;
     int ioff;
{
  SXPEC *pxpec0;
  /*@only@*/ static double *t, *tv, *tegy, *twk;
  /*@only@*/ static short *iassgn, *idx;
  static short ibkptr[2];
  static double c, egy0, fval[9];
  static int init = 0;
  static int mbgn, nsym, nbasis, mbasis, nf, np;
  static int nsq, isdgn0, nxvib, nxdim, nxall;
  double *pt, *ptt;
  double vsym, val, rhokx, tmp, pmix, cv, sv;
  int i, iv, ivp, k, icol, nv, nlim, m, ivdgn, minv, isdgn, mm, ixvib;
  int iend, ixv, ixvp, isig, nvib0, kk, knt, kd, kmix, kkp;
  size_t nl;

  pxpec0 = *pxpecv;
  if (init == 0) { /* initialize */
    if (nbasis > 0) {
      free(iassgn);
      iassgn = NULL;
      free(idx);
      idx = NULL;
      free(tegy);
      tegy = NULL;
      free(twk);
      twk = NULL;
      free(t);
      t = NULL;
      free(tv);
      tv = NULL;
    }
    nbasis = 0; kkp = 0;
    if (nvib == 0) return 0;
    nxvib = pxpec0->ixvib; nxdim = 1; nxall = nxvib;
    if (nxvib < 0) {
      nxdim = -nxvib; nxall = nxdim;
    }
    nsym  = pxpec0->nval;
    m = pxpec0->mval;
    mbgn = -m * nsym;
    mbasis = m + m + 1; nbasis = mbasis * nxdim;
    ibkptr[0] = 0; ibkptr[1] = (short) nbasis;
    np = nbasis + 1;
    nl = (size_t) np * sizeof(short);
    iassgn = (short *) mallocq(nl);
    iassgn[0] = 0;
    idx = (short *) mallocq(nl);
    idx[0] = 0;
    nl = (size_t) nbasis * sizeof(double);
    tegy = (double *) mallocq(nl);
    twk = (double *) mallocq(nl);
    nsq = nbasis * nbasis;
    nl = (size_t) nsq * sizeof(double);
    t = (double *) mallocq(nl);
    if (nxvib > 1) nl *= nxvib;
    tv = (double *) mallocq(nl);
    dcopy(nxall, &zero, 0, fval, 1); 
    dcopy(nsq, &zero, 0, tv, 1); 
    for (k = 1; k < nxvib; ++k) 
      dcopy(nsq, &zero, 0, &tv[k * nsq], 1);
    nf = 0;
    minv = 1024; 
    iend = pxpec0->kdel;
    for (i = 1; i <= iend; ++i) {
      val = pxpecv[i]->xval;
      nv = pxpecv[i]->nval;
      iv = pxpecv[i]->ixvib;
      if (nv >= 25)
        continue;
      if (pxpecv[i]->mval == 0) {       /*  calculate potential terms */
        if (nv > 0 && nv < minv && (iv % 11) == 0)
          minv = nv;
        if (nv > 0) 
          val *= 0.5;
        k = nv / nsym;
        if (iv == 0) {
          pt = tv; ivp = 0;
        } else {
          ivp = iv / 10; iv -= ivp * 10; /* iv >= ivp */
          if (iv < ivp) continue;
          if (nxvib < 0) 
            ixv = (iv + ivp * nbasis) * mbasis;
          else
            ixv = nsq * iv;
          pt = &tv[ixv];
        }
        kd = 0;
        if (iv != ivp ) {
          if ((nsym & 1) == 0) {
            kk = morgsym[ivp] - morgsym[iv] - nv;
            k = kk / nsym;
            if (kk != k * nsym) continue;
            kd = (kk + nv + nv) / nsym;
          } else {
            kd = k; k = -k;
          }
        }
        for (knt = 0; knt < 2; ++knt) {
          ixv = k;
          if (k < 0) { 
            k = -k; ixv = k * nbasis;
          }
          for (kk = k; kk < mbasis; ++kk) {
            pt[ixv] += val; ixv += np;
          }
          if (nv ==  0 || iv == ivp) break;
          k = kd;
          nv = 0;
        }
        if (nv == 0) continue;
        val *= -2.; 
        ixv = 0;
        for (kk = 0; kk < mbasis; ++kk) {
          pt[ixv] += val; ixv += np;
        }
      } else {   /* isolate allowed F values */
        iv = iv / 11;
        fval[iv] += val;
      }
    }
    isdgn0 = minv;
    if (nsym >= 3 || minv < 3) isdgn0 = 0;
    c = 29979.2458;
  }
  /* begin normal computation */
  if (nvib <= 0) return -1;
  if (nxdim == 1) {
    ixvib = ioff / nvib; isig = ioff - nvib * ixvib;
    iend = ioff + nvib;
  } else {
    ixvib = 0; isig = ioff; iend = nvib * nxdim;
  }
  nvib0 = nvib / nsym;
  pt = tv;
  if (ixvib > 0) 
    pt = &tv[ixvib * nsq];
  dcopy(nsq, pt, 1, t, 1);
  /* calculate kinetic energy */
  rhokx = eig->morg[ioff] - eig->rhok[ixvib]; 
  isdgn = 0;
  if (init <= 0 || fabs(rhokx - 0.5 * nsym) < tiny)
    isdgn = isdgn0;
  kd = 0; kk = 0;
  iv = 0; nv = mbasis;  
  for (ixv = 0; ixv < nxdim; ++ixv) {
    val = fval[ixv + ixvib];
    m = mbgn;
    for (k = 0; k < nv; ++k) {
      tmp = m + rhokx;
      t[iv] += val * tmp * tmp;
#if MYDEBUG_EIG 
      if (lu != 0) {
        kkp = kk;
        for (i = kk; i < nbasis; ++i) { 
          pmix = t[iv + i - kkp];
          if (fabs(pmix) < tiny) continue;
          fprintf(lu, " %4d%4d%15.6e", i, kk, pmix);
          if (++kd > 4) {
            fputc('\n',lu); kd = 0;
          }
        }
      }
#endif
      m += nsym; iv += np; 
      ++kk;
    }
    rhokx += morgsym[ixv + 1] - morgsym[ixv];
  }
#if MYDEBUG_EIG 
  if (lu != 0) fputc('\n',lu);
#endif
  /* diagonalize */
  m = hdiag(nbasis, nbasis, t, tegy, twk, iassgn);
  ordblk(nbasis, nbasis, iassgn, t, tegy, ibkptr, twk, idx);
  m = pxpec0->mval; ivp = 0; nlim = nvib0; ixv = 0;
  for (iv = ioff; iv < iend; iv += nsym) { /* loop over states */
    val = tegy[ivp]; pt = &t[ivp * nbasis]; pmix = twk[ivp];
#if MYDEBUG_EIG 
    if (lu != 0) {
      kd = 0;
      fprintf(lu, "eignvector %3d\n", iv);
      for (i = 0; i < nbasis; ++i) { 
        tmp = pt[i];
        fprintf(lu, " %3d%15.6e", i, tmp);
        if (++kd > 4) {
          fputc('\n',lu); kd = 0;
        }
      }
      fputc('\n',lu);
    }
#endif
    pmix = 0.; kmix = 0; ptt = pt;
    for (ixvp = 0; ixvp < nxdim; ++ixvp) {
      if (ixvp > 0)
        ptt += mbasis;
      vsym = ddot(mbasis, ptt, 1, ptt, 1);
      if (vsym > pmix) {
        kmix = ixvp; pmix = vsym;
      }
    }
    ixvp = ixv; ++ivp;
    if (ivp >= nlim) {
      ++ixv; nlim = ivp + nvib0;
    }
    ivdgn = -1; 
    if (isdgn != 0 && nxvib > 0) {
      vsym = (1.e-10) * fabs(val); i = ivp;       
      /* find energy degeneracy */
      if (fabs(tegy[ivp] - val) < vsym) {
        ivdgn = ivp;
        ptt = &t[ivdgn * nbasis];
        mm = m;
        if (init > 0)
          --mm;
        vsym = 0.; sv = 0.;
        for (i = 1; i <= m; ++i) {
          vsym += pt[mm + i] * pt[m - i];
          sv  += pt[mm + i] * ptt[m - i] + ptt[mm + i] * pt[m - i];
        }
        i = iv;
        if (m == mm) {
          i += iv / isdgn;
          vsym = 2. * vsym + pt[m] * pt[m];
          sv += pt[m] * ptt[m]; 
        }
        cv = 1. - vsym;
        tmp = sqrt(cv * cv + sv * sv);
        cv /= tmp; sv /= tmp;
        if ((i & 1) == 0) {
          tmp = cv; cv = sv; sv = -tmp;
        }
        drot(nbasis, ptt, 1, pt, 1, cv, sv);
      }
    }
    if (init == 0) {
      init = -1;
      egy0 = val;
      if (lu != NULL) 
        fprintf(lu, 
                "  v,    rho*K, vmix, pmix, energy relative to %15.6f\n",
                egy0 / c);
    }
    val -= egy0;
    if (lu != NULL)
      fprintf(lu, "%3d %9.6f %2d %9.6f %15.6f\n", 
              iv, eig->rhok[ixvib], kmix, pmix, val / c);
    egy[iv] = val;
    icol = iv * nbasis;
    dcopy(nbasis, pt, 1, eig->eigvec + icol, 1);
  }
  init = 1;
  return 0;
}                               /* iamdiag */

int expec(ityp, indx, pxpecv, nvib, eigm, eigp, noff, isym)
     int ityp, indx;
     SXPEC **pxpecv;
     int nvib;
     SEIG *eigm, *eigp;
     int noff, *isym;
{
  static double facv[] = { 1.0, 0.5, 0.5, 0.25 };
  static double sumv[4];
  SXPEC *pxpecp, *pxpec0;
  SEIG *eigbra, *eigket;
  double *result;
  double *eigvbra, *eigvket;
  double pabra, paket, vtmp1, vtmp2, tmp, fac;
  int idel, nsym, i, nbasis, mbra, mket, maxm, isin, nxvib;
  int state, state0, nv, nv0, mv, m, nx, iv1, iv2, kd, isum, nvtot, nvcmp;
  int iv1s, iv2s, iv1d, iv2d, ivv, nn, mbasis, ixv1, ixv2, ixv1p, ixv2p;
  int nblk,ibgn1, ibgn2, ixv1m, ixv2m, iend1, iend2, nlim, iv1z, iv2z;
  int ivv0, isig1, isig2;
  /*    state=0 <eigm(iv1) |exp(-i*nval*tau)| eigp(iv2)> */
  /*    state=1 <eigm(iv1) |exp(+i*nval*tau)| eigp(iv2)> */
  /*    state=2 <eigp(iv1) |exp(-i*nval*tau)| eigm(iv2)> */
  /*    state=3 <eigp(iv1) |exp(+i*nval*tau)| eigm(iv2)> */
  pxpec0 = pxpecv[0];
  nsym = pxpec0->nval;
  maxm = pxpec0->mval;
  nxvib = pxpec0->ixvib;
  mbasis = maxm + maxm + 1; nbasis = mbasis;
  if (nxvib < 0)
    nbasis *= -nxvib;
  maxm *= nsym;
  pxpecp = pxpecv[ityp];
  result = pxpecp->xpec + indx;
  nv0 = pxpecp->nval;
  mv = pxpecp->mval;
  kd = pxpecp->kdel;
  ixv1 = pxpecp->ixvib;
  ixv2 = ixv1 / 10; ixv1 -= 10 * ixv2;
  if (nxvib > 0) {
    nblk = nvib; nvtot = nblk * nxvib;
    ibgn1 = ixv1 * nblk; ixv1p = 0; ixv1m = 0;
    ibgn2 = ixv2 * nblk; ixv2p = 0; ixv2m = 0;
  } else { 
    nblk = -nxvib * nvib; nvtot = nblk;
    ibgn1 = 0; ixv1p = ixv1 * mbasis; ixv1m = ixv1 * nvib; 
    ibgn2 = 0; ixv2p = ixv2 * mbasis; ixv2m = ixv2 * nvib;
  }
  iend1 = ibgn1 + nblk - 1;
  iend2 = ibgn2 + nblk - 1;
  isin = 0;
  if (nv0 < 0) {
    isin = 1;
    nv0 = -nv0;
  }
  state0 = 0;
  if (kd != 0)
    state0 = 2;
  if (nv0 != 0)
    ++state0;
  nvcmp = nv0 % nsym;
  fac = facv[state0];
  iv1z = ibgn1; iv2z = ibgn2; nv = nv0;
  /* loop through vibrational states */
  for (iv2s = ibgn2; iv2s <= iend2; ++iv2s) {
    iv2d = (int) isym[iv2s] >> 2;
    if (iv2d == 0) 
      iv2z = iv2s;
    isig2 = iv2s - iv2z;
    nlim = (iend1 < iv2s)? iend1: iv2s;
    ivv0  = iv2s * nvtot; 
    for (iv1s = ibgn1; iv1s <= nlim; ++iv1s) { /* iv1s <= iv2s */
      iv1d = (int) isym[iv1s] >> 2;
      if (iv1d == 0) 
        iv1z = iv1s;
      isig1 = iv1s - iv1z;
      if (nvcmp == 0 && isig1 != isig2) continue;
      eigbra = eigp; eigket = eigm;
      iv1 = iv1s; iv2 = iv2s;
      isum = 0;
      for (state = state0; state >= 0; --state) {
        if (state == 1) {
          if (nv0 == 0) 
            state = 0;
          iv1 = iv1z + isig2;
          iv2 = iv2z + isig1;
          eigbra = eigm;
          eigket = eigp;
        }
        nv = nv0;
        if ((state & 1) == 0)
          nv = -nv;
        sumv[state] = 0.;
        mbra = eigbra->morg[ixv1m + iv1 - iv1z];
        mket = eigket->morg[ixv2m + iv2 - iv2z];
        idel = nv + mket - mbra;
        if (nsym > 1) {
          if ((idel % nsym) == 0) { /*  check for good symmetry */
            idel /= nsym; 
          } else {
            continue; 
          }
        }
        if (idel > 0) {
          nx = mbasis - idel;
          if (nx <= 0) continue;
          eigvbra = eigbra->eigvec + (iv1 * nbasis + ixv1p + idel);
          eigvket = eigket->eigvec + (iv2 * nbasis + ixv2p);
          mbra += idel * nsym;
        } else {
          nx = mbasis + idel;
          if (nx <= 0) continue;
          eigvbra = eigbra->eigvec + (iv1 * nbasis + ixv1p);
          eigvket = eigket->eigvec + (iv2 * nbasis + ixv2p - idel);
          mket -= idel * nsym;
        }
        if (mv == 0) {
          tmp = ddot(nx, eigvbra, 1, eigvket, 1);
        } else {              /* mv > 0 */
          pabra = (mbra - maxm) - eigbra->rhok[ixv1];            
          paket = (mket - maxm) - eigket->rhok[ixv2];
          if (mv == 2 && nv == 0 && iv1s == 2 && iv2s == 2 && kd == 0)
            tmp = 0.;
          tmp = 0.;
          for (i = 0; i < nx; ++i) {
            vtmp1 = 0.5 * eigvbra[i] * eigvket[i]; vtmp2 = vtmp1;
            for (m = 0; m < mv; ++m) {
              vtmp1 *= pabra;
              vtmp2 *= paket;
            }
            tmp += vtmp1 + vtmp2;
            pabra += nsym; paket += nsym;
          }
        }
        sumv[state] = tmp;
        isum += nx;
      }  /* state loop */
      if (isum == 0) continue;
      if (isin != 0) {
        sumv[0] -= sumv[1];
        if (state0 == 3) sumv[2] -= sumv[3];
      } else if (nv != 0) {
        sumv[0] += sumv[1];
        if (state0 == 3) sumv[2] += sumv[3];
      }
      ivv = ivv0 + iv1s;
      if (state0 >= 2) {
        tmp = sumv[2] - sumv[0]; sumv[0] += sumv[2];
        if (isin != 0) tmp = -tmp;
        result[ivv + noff] = fac * tmp;
      }
      result[ivv] = fac * sumv[0];
    }   /* iv1 loop */
    nn = nlim - ibgn1; 
    ivv = iv2s  + nvtot * ibgn1; ivv0 += ibgn1;
    if (state0 >= 2)
      dcopy(nn, &result[ivv0 + noff], 1, &result[ivv + noff], nvtot);
    dcopy(nn, &result[ivv0], 1, &result[ivv], nvtot);
  }                             /* iv2 loop */
  return 0;
}                               /* expec */

int getparv(nv, mv, iv, nvib2, parp, pxpecv)
     int nv, mv, iv, nvib2;
     SIPAR *parp;
     SXPEC **pxpecv;
{
  static int n = 0;
  SXPEC *pxpecp, *pxpeci, *pxpec0;
  double val;
  long ityp;
  int i, iend, k, kd, nsym, md, knt, nknt, isunit, nlim;
  BOOL inham;

  ityp = parp->idval;
  if (ityp == 990099L) {
    mv = nv = 0;
  } else if (ityp == 800099L) {
    mv = nv = 0; parp->parval = 10000.;
  }
  pxpec0 = pxpecv[0];
  nsym = pxpec0->nval;
  md = pxpec0->mval;
  nlim = pxpec0->isunit;
  kd = ideqv(&ityp, &val); 
  if (mv >= 100) 
    mv -= 100;
  nknt = 1; isunit = 0;
  if (nv == 0 && mv == 0 && pxpec0->ixvib >= 0) {
    isunit = kd + 1;
  } else if (nv <= -25 || nv >= 25) { 
    isunit = -1;
  }
  if (nvib2 > 0) {
    nknt = 2;
    if (isunit != 0) {
      if (iv == 11) {
        parp->ptxpec = NULL; parp->ptxpecx = NULL;
        return n;
      } 
    }
  }
  for(knt = 0; knt < nknt; ++knt) {
    inham = FALSE;
    k = nv / nsym;
    /* check symmetry for inclusion in diagonalization */
    if (ityp == -99L && kd == 0 && nv >= 0 && nv < 25 && k <= (md + md)) {
      if ((iv % 11) == 0) {
        if (k * nsym == nv) {
          /* nv >= 0 and is a multiple of nsym */
          if (mv == 2 && nv == 0) {
            inham = TRUE;
          } else if (mv == 0 && nv >= 0) {
            inham = TRUE;
          } 
        }
      } else {
        if (mv == 0 && nv >= 0) {
          inham = TRUE;
        }
      } 
    }
    val = 0.; 
    if (inham && isunit == 0) {
      val = parp->parval; kd = -1;
      parp->flg |= MHAM;
    }
    ++n; iend = n;
    pxpecp = pxpecv[n];
    pxpecp->kdel = kd;
    pxpecp->mval = mv;
    pxpecp->nval = nv;
    pxpecp->xval = val;
    pxpecp->xpec = NULL;
    pxpecp->ixvib = iv;
    pxpecp->isunit = isunit;
    for (i = 1; i < n; ++i) {
      pxpeci = pxpecv[i];
      k = pxpeci->kdel - kd;
      if (k > 0) {                /* insert in list */
        for (k = n; k > i; --k) {
          pxpecv[k] = pxpecv[k - 1];
        }
        pxpecv[i] = pxpecp;
        iend = i;
        break;
      }
      if (k == 0 && pxpeci->nval == nv && pxpeci->mval == mv && 
          pxpeci->ixvib == iv) {
        /* full match */
        pxpecp = pxpeci;
        --n;
        if (kd < 0) pxpecp->xval += val;
        break;
      }
    } 
    if (n >= nlim) {
      printf("pexpecv too small :%d\n", nlim);
      exit(EXIT_FAILURE);
    }
    if (inham && pxpec0->kdel < iend)
      pxpec0->kdel = iend;
    if (knt > 0) {
      if (iv == 0) {
        parp->ptxpecx = parp->ptxpec;
        parp->ptxpec  = pxpecp;
        parp->flg |= MREV;
      } else {
        parp->ptxpecx = pxpecp;
      }
      break;
    }
    parp->ptxpec = pxpecp;
    iv = 11 - iv;
    if (isunit >= 0) {
      if (kd < 0) 
        kd = 0;
      if (isunit == 0) 
        isunit = 1 + kd;
      mv = 0; nv = 0;
    }
  }
  return n;
}                               /* getparv */

int ideqv(pidtmp, pfac)
     long *pidtmp;
     double *pfac;
{
#define IDIV 10000
  static struct {
    double fac;
    int ityp, kd;
  } cmpv[] = { 
    {-2.0,  4, 2},
    { 0.5, 21, 2}, 
    {-1.0, 40, 1},
    { 1.0, 60, 1},
    {-1.0, 61, 1},
    { 1.0, 41, 1}, 
    { 0.5, 22, 2},
    {-2.0, 12, 2},
    { 1.0,  0, 0},
  };
  long idtmp;
  int kdret, i, ii, key;

  kdret = 0;
  *pfac = 1.; 
  idtmp = (*pidtmp);
  if (idtmp < 0)
    idtmp = -idtmp;
  if (idtmp < 99) {
    key = (int) idtmp % 10;
    if (key == 2) {
      ++idtmp;
      kdret = 1;
      *pfac = -1.;
    } else if (key == 3) {
      --idtmp;
      kdret = 1;
    }
  } else {
    key = (int) ((idtmp / IDIV) % 100); 
    if (key != 0) {
      for (i = 0; (ii = cmpv[i].ityp) != 0; ++i) {
        if (ii == key) {
          kdret = cmpv[i].kd;
          *pfac = cmpv[i].fac;
          if (kdret != 0) {
            ii = i ^ 1;
            idtmp += IDIV * (long) (cmpv[ii].ityp - key);
          }
          break;
        }
      }
    }
  }
  *pidtmp = -idtmp;
  return kdret;
}                               /* ideqv */

int ftanal(iflg, nft, txpec, xfc)
     int iflg, nft;
     double *txpec, *xfc;
{
  int k, kk, nft2;
  dcopy(nft, txpec, 1, &txpec[nft], 1);
  iflg &= 1;
  nft2 = nft << 1;
  kk = 3;
  if (iflg != 0) {
    dneg(nft, &txpec[nft], 1);
    kk = 1;
  }
  ftran(nft2, txpec, 1, xfc, 0);
  for (k = 1; k < nft; k += 2) {
    xfc[k] = xfc[kk];
    xfc[k + 1] = xfc[kk + 1];
    kk += 4;
  }
  if (iflg != 0)
    nft = nft2;
  return nft;
}
int ksign(kdel, nft, txpec, xfc)
     int kdel, nft;
     double *txpec, *xfc;
{
  double dif,val0;
  int k, kp, km, iret;
  iret = nft - kdel;
  if ((kdel & 1) != 0) {
    val0 = 0.; ++iret; kp = 0;
  } else {
    val0 = txpec[0]; kp = 1;
  }
  km = nft - 1;
  xfc[0] = val0; 
  for (k = 1; k < iret; k += 2) {
    dif = 0.5 * (txpec[kp] - txpec[km]);
    xfc[k] = txpec[kp] - dif - val0;
    xfc[k + 1] = dif;
    --km; ++kp;
  }
  for (k = iret; k < nft; ++k) 
    xfc[k] = 0.;
  return iret;
}
void ftfix(int nft2, int nft, double *xfc)
{
  static double shalf, sbgn, chalf, cbgn;
  double c,s,vc,vs;
  int k;
  if (nft2 == 0) {
    vs = asin(1.) * (*xfc);
    sbgn = sin(vs); cbgn = cos(vs);
    vs *= 2.;
    shalf = sin(vs); chalf = cos(vs);
    return;
  }
  if (nft2 == nft) {
    s = shalf; c = chalf;
  } else {
    s = sbgn; c = cbgn;
  }
  for (k = 1; k < nft; k += 2){
    vc = xfc[k]; vs = xfc[k + 1];
    xfc[k + 1] = c * vs + s * vc;
    xfc[k]     = c * vc - s * vs;
    vs = s;
    s = c * shalf + vs * chalf;
    c = c * chalf - vs * shalf;
  }
}

void dneg(int n, double *dx, int incx)
{
  int k;
  if (n <= 0) return;
  *dx = -(*dx);
  for (k = 1; k < n; ++k) {
    dx += incx;
    *dx = -(*dx);
  }
}
void sxpec_alloc(PSXPEC *pxpecv, const int nvec)
{
  SXPEC *pxpec;
  int i;
  pxpec = pxpecv[0];
  pxpec->xpec = NULL;
  for (i = 1; i <= nvec; ++i) {
    pxpecv[i] = (++pxpec);
    pxpec->xpec = NULL;
    pxpec->nval = pxpec->mval = pxpec->kdel = 0;
    pxpec->xval = 0.;
  }
}
void spar_alloc(SIPAR *pswk, PSIPAR *par, const int nvec)
{
  SIPAR *pspar;
  int i;
  par[0] = pswk;
  pspar = par[0]; pspar->ptxpec = NULL; pspar->ptxpecx = NULL;
  for (i = 1; i < nvec; ++i) {
    par[i] = (++pspar); pspar->ptxpec = NULL; pspar->ptxpecx = NULL;
  }
}
void seig_alloc(SEIG *eigv, const int nvec, const size_t nfill,
                  int *pmorgwk, const int nvib)
{
  SEIG *eig;
  int *pmorg;
  int k;
  k = 0;
  eig = eigv;
  pmorg = pmorgwk;
  do {
    eig->eigvec = (double *) mallocq(nfill);
    eig->morg = pmorg;
    eig->irhok = 0;
    ++eig;
    pmorg += nvib;
  } while (++k <= nvec);
}
int init2top(FILE *lu, int nvib2, int nsym, int nvib, int nvtot,
              int nsig, short *kmap, S2TOP *vmap, int *sigma2, int *ivtop)
{  
  S2TOP *pmap, *pmapp;
  int i, iv, iv1, iv2, k, kk, nline, isig1, isig2, isig, isold, nvt, nvtop;
  sigma2[0] = 0;
  if (nvib2 <= 0) nvtot = 1;
  ivtop[0] = 0; k= 0;
  do {
    pmap = &vmap[k];
    sigma2[k] = 0; 
    pmap->iqv = -16; pmap->iqsig = 0;   
    pmap->ix1 = 0; pmap->ix2 = 0; pmap->inew = 0;
     ++k;
  } while (k < nvtot);
  if (nvib2 <= 0) return 0;
  pmap = NULL;
  nvt = nvib / nsym; nvtop = 0;
  nline = (int) kmap[0] * ND2TOP + 1;
  /* user input equivalent states */
  for (k = 1; k < nline; k += ND2TOP) {
    iv = (int) kmap[k]; 
    if (iv < 0 || iv >= nvtot) continue;
    pmap = &vmap[iv];
    isig1 = (int) kmap[k + 2]; 
    if (isig1 < 0) isig1 += nsym;
    if (isig1 >= nsym) continue;
    isig2 = (int) kmap[k + 3];
    if (isig2 < 0) isig2 += nsym;
    if (isig2 >= nsym) continue;
    i = (int) kmap[k + 1];
    if (isig1 == isig2) i = 0;
    if (i > 0) i = 2;
    else if (i < 0) i = 3;
    pmap->iqsig = i;
    if (isig1 > isig2) i ^= 1;
    if (i == 3) {
      kk = isig1; isig1 = isig2; isig2 = kk;
    }
    iv1 = (int) kmap[k + 5];
    if (iv1 < 0 || iv1 >= nvt) continue;
    iv2 = (int) kmap[k + 6];
    if (iv2 < 0 || iv2 >= nvt) continue;
    i = (int) kmap[k + 4];
    if (kmap[k + 5] == kmap[k + 6]) i = 0;
    if (i > 0) i = 2;
    else if (i < 0) i = 3;
    pmap->iqv = i;
    if (iv1 > iv2) i ^= 1;
    if (i == 3) {
      kk = iv1; iv1 = iv2; iv2 = kk;
    }
    pmap->ix1 = iv1 * nsym + isig1;
    pmap->ix2 = iv2 * nsym + isig2;
    for (kk = 0; kk < nvtop; kk += 2) {
      if (ivtop[kk] == iv1 && ivtop[kk + 1] == iv2) break;
    }
    if (kk == nvtop) {
      if (iv1 > iv2) {
        kk = iv1; iv1 = iv2; iv2 = kk;
      }
      ivtop[nvtop++] = iv1; ivtop[nvtop++] = iv2;
      if (iv1 != iv2) {
        ivtop[nvtop++] = iv2; ivtop[nvtop++] = iv1;
      }
    }
  }
  pmap = NULL;
  /* find default sigmas */
  for (k = 0; k <  nvtot; ++k) {
    pmap = &vmap[k];
    if (pmap->iqv >= 0) continue;
    iv = k / nsig; isig = k - iv * nsig;
    for (kk = isig; kk < nvtot; kk += nsig) {
      pmapp = &vmap[kk];
      if (pmapp->iqv < 0) continue;
      pmap->iqv += 8;
      pmap->iqsig  = pmapp->iqsig;
      pmap->ix1 += pmapp->ix1 % nsym;
      pmap->ix2 += pmapp->ix2 % nsym;
      break;
    }
  }
  pmap = NULL; pmapp = NULL;
  /* find default vibrations */
  for (k = 0; k < nvtot; ++k) {
    pmap = &vmap[k];
    if (pmap->iqv >= 0) continue;
    iv = k / nsig; i = iv * nsig;
    for (kk = 0; kk < nsig; ++kk) {
      pmapp = &vmap[i + kk];
      if (pmapp->iqv < 0) continue;
      pmap->iqv +=  pmapp->iqv + 8;
      pmap->ix1 += (pmapp->ix1 / nsym) * nsym;
      pmap->ix2 += (pmapp->ix2 / nsym) * nsym;
      break;
    }
  } 
  pmap = NULL; pmapp = NULL;
  /* find index for states */
  for (k = 0; k < nvtot; ++k) {
    pmap = &vmap[k]; pmap->inew = 0;
    if (pmap->iqv < 0) continue;
    isig1 = pmap->ix1; iv1 = isig1 / nsym; isig1 -= iv1 * nsym;
    isig2 = pmap->ix2; iv2 = isig2 / nsym; isig2 -= iv2 * nsym;
    for (kk = 0; kk < nvtop; kk += 2) {
      if (ivtop[kk] == iv1 && ivtop[kk + 1] == iv2) {
        pmap->inew = isig2 + nsym * (isig1 + nsym * (kk >> 1));
        break;
      }
    }
  }
  pmap = NULL;
  /* check symmetry */
  for (k = 0; k < nvtot; ++k) {
    pmap = &vmap[k];
    if (pmap->iqv < 0) continue;
    isig1 = pmap->ix1; iv1 = isig1 / nsym; isig1 -= iv1 * nsym;
    isig2 = pmap->ix2; iv2 = isig2 / nsym; isig2 -= iv2 * nsym;
    if (pmap->iqv == 3) {
      iv = pmap->inew - nsym * nsym; 
      for (kk = 0; kk < nvtot; ++kk) {
        pmapp = &vmap[kk];
        if (pmapp->inew == iv) {
          if (pmapp->iqv == 0) pmapp->iqv = 2;
          break;
        }
      }
    }
    if (pmap->iqsig == 3) {
      iv = pmap->inew + (nsym - 1) * (isig2 - isig1); 
      for (kk = 0; kk < nvtot; ++kk) {
        pmapp = &vmap[kk];
        if (pmapp->inew == iv) {
          if (pmapp->iqsig == 0) pmapp->iqsig = 2;
        }
      }
    }         
  } 
  pmap = NULL; pmapp = NULL;
  /* fill in sigma */
  isold = 0; kk = 0; iv2 = 0;
  for (k = 0; k < nsig; ++k){
    iv1 = iv2;
    pmap = &vmap[k];
    isig1 = pmap->ix1 % nsym; isig2 = pmap->ix2 % nsym;
    iv2 = isig1 + isig2;
    if (isig1 == 0 && isig2 == 0) {
      sigma2[k] = 0; isold = 0;
      continue;
    }
    if (isold != 0) {
      isold = 0;
      if (iv2 == nsym) {
        sigma2[k] = kk;
        continue;
      }
      if (((iv1 + iv2) % nsym) == 0) {
        sigma2[k] = -kk;
        continue;
      }
    }
    if (isold == 0) {
      if (isig1 !=  0 && isig2 != 0) ++kk;
      if (kk == 0) kk = 1; 
      sigma2[k] = kk; isold = 1;
    }
  }
  for (k = 0; k < nsig; ++k) {
    kk = sigma2[k];
    for (iv = nsig; iv < nvtot; iv += nsig)
      sigma2[k + iv] = kk;
  }
  if (lu == NULL) return nvtop;
  fputs("2 top mapping new v: iq old sig1 sig2, iqv vt, lz\n", lu);
  for (k = 0; k < nvtot; ++k){
    pmap = &vmap[k];
    i = pmap->iqv;
    if (i == 2) i = 1;
    if (i == 3) i = -1;
    isig = pmap->iqsig;
    if (isig == 2) isig = 1;
    if (isig == 3) isig = -1;
    isig1 = pmap->ix1; isig2 = pmap->ix2;
    iv2 = pmap->ix2 / nsym; iv1 = pmap->ix1 / nsym;
    fprintf(lu,"%3d: %2d %3d %3d, %2d %2d%2d, %3d %3d\n", k, 
            isig, isig1, isig2, i, iv1, iv2, sigma2[k], pmap->inew);
  }
  return nvtop;
} /* init2top */

int get2top(SIPAR *ppar,/*@out@*/ double *xpec2, int nvtot, int nft, int nsym,
            int *vmapsym, int *ivtop, int nvtop)
{
  SXPEC *pxpec1, *pxpec2;
  double *dxpec1, *dxpec2, *dxpec;
  double va1, va2, vb1, vb2, fac;
  int nsymsq, iv, ivp, iv1, iv1p, iv2, iv2p, isig1, isig1p, isig2, isig2p, ivv;
  int ix, ix1, ix2, ndm2, noff, noff2, ioff, isig, isunit1, isunit2;
  int kdel, ndtot, ndtot2, ift, nft2, iret, nvib, nv1, nv2, isigp;
  kdel = -2; 
  pxpec1 = ppar->ptxpec;\
  if (pxpec1 != NULL)  kdel = pxpec1->kdel;
  nsymsq = nsym * nsym; ndm2 = nsymsq * (nvtop >> 1); ndtot2 = ndm2 * ndm2;
  dxpec = xpec2; xpec2[0] = 0.;
  nft2 = nft; 
  if (kdel > 0) 
    nft2 = nft2 << 1;
  for (ift = 0; ift < nft2; ++ift) {
    if (ift > 0) dxpec += ndtot2;
    dcopy(ndtot2, &zero, 0, dxpec, 1);
  }
  if (pxpec1 == NULL) return 0;
  pxpec2 = ppar->ptxpecx; 
  if (pxpec2 == NULL) return 0;
  dxpec1 = pxpec1->xpec; dxpec2 = pxpec2->xpec;
  if (dxpec1 == NULL || dxpec2 == NULL) return 0;
  noff2 = ndtot2 * nft; ndtot = nvtot * nvtot; nvib = nvtot >> 1;
  noff  = ndtot  * nft; iret = 0;
  isunit1 = pxpec1->isunit; isunit2 = pxpec2->isunit;
  nv1 = pxpec1->nval; nv2 = pxpec2->nval; 
  /* load operators */
  ioff = 0;
  for (iv = 0; iv < ndm2; ++iv) {
    ivv = iv / nsymsq; isig = iv - ivv * nsymsq; 
    isig1 = isig  / nsym; isig2 = isig - isig1 * nsym;
    ivv = ivv << 1; 
    iv1 = isig1 + nsym * ivtop[ivv]; 
    iv2 = isig2 + nsym * ivtop[ivv + 1] + nvib;
    for (ivp = 0; ivp < ndm2; ++ivp) {
      ivv = ivp / nsymsq; isigp = ivp - ivv * nsymsq;
      isig1p = isigp  / nsym; isig2p = isigp - isig1p * nsym;
      if (nv1 == 0 && isig1p != isig1) continue;
      if (nv2 == 0 && isig2p != isig2) continue;
      ivv = ivv << 1; 
      iv1p = isig1p + nsym * ivtop[ivv]; 
      iv2p = isig2p + nsym * ivtop[ivv + 1] + nvib;
      if (isunit1 == 1 && iv1p != iv1) continue; 
      if (isunit2 == 1 && iv2p != iv2) continue; 
      ix1 = iv1 + nvtot * iv1p; ix2 = iv2 + nvtot * iv2p;
      va1 = dxpec1[ix1]; va2 = dxpec2[ix2];
      ix = iv + ndm2 * ivp;
      if (kdel > 0)
        ioff = (isig - isigp) * (ndm2 - 1) + noff2;
      for (ift = 0; ift < nft; ++ift) {
        if (ift > 0) {
          ix1 += ndtot; ix2 += ndtot; ix += ndtot2;
        }
        va1 = dxpec1[ix1]; va2 = dxpec2[ix2];
        if (kdel > 0) {
          vb1 = dxpec1[ix1 + noff]; vb2 = dxpec2[ix2 + noff];
          if (nv1 < 0) vb1 = -vb1;
          if (nv2 < 0) vb2 = -vb2;
          xpec2[ix] = (va1 + vb1) * (va2 + vb2);
          xpec2[ix + ioff] = (va1 - vb1) * (va2 - vb2);
        } else {
          xpec2[ix] = va1 * va2;
        }
      }
    }
  }
  /* sigma symmetry */ 
  for (iv = 0; iv < ndm2; ++iv) {
    if ((vmapsym[iv] & 3) == 0) continue; 
    ivv = iv / nsymsq; isig = iv - ivv * nsymsq;
    isig1 = isig  / nsym; isig2 = isig - isig1 * nsym;
    if (isig1 >= isig2) continue;
    ivv = iv + (nsym - 1) * (isig2 - isig1);
    ix1 = iv * ndm2; ix2 = ivv * ndm2; 
    for (ift = 0; ift < nft2; ++ift) {
      if (ift > 0) {
        ix1 += ndtot2; ix2 += ndtot2;
      }
      drot(ndm2, &xpec2[ix1], 1, &xpec2[ix2], 1, rthalf, rthalf);
    }
    ix1 = iv; ix2 = ivv; 
    for (ift = 0; ift < nft2; ++ift) {
      if (ift > 0) {
        ix1 += ndtot2; ix2 += ndtot2;
      }
      drot(ndm2, &xpec2[ix1], ndm2, &xpec2[ix2], ndm2, rthalf, rthalf);
    }
  }
  /* vibrational symmetry */ 
  for (iv = 0; iv < ndm2; ++iv) {
    if ((vmapsym[iv] & 0x0c) == 0) continue; 
    ivv = iv / nsymsq;
    ivv = ivv << 1; 
    iv1 = ivtop[ivv]; 
    iv2 = ivtop[ivv + 1];
    if (iv1 >= iv2) continue;
    iv1p = ivtop[ivv + 2]; iv2p = ivtop[ivv + 3];
    ivv = iv + nsymsq;
    ix1 = iv * ndm2; ix2 = ivv * ndm2;
    for (ift = 0; ift < nft2; ++ift) {
      if (ift > 0) {
        ix1 += ndtot2; ix2 += ndtot2;
      }
      drot(ndm2, &xpec2[ix1], 1, &xpec2[ix2], 1, rthalf, rthalf);
    }
    ix1 = iv; ix2 = ivv; 
    for (ift = 0; ift < nft2; ++ift) {
      if (ift > 0) {
        ix1 += ndtot2; ix2 += ndtot2;
      }
      drot(ndm2, &xpec2[ix1], ndm2, &xpec2[ix2], ndm2, rthalf, rthalf);
    }
  }
  if (kdel > 0) { /* final operator fixup */
    ivv = 0;
    if (ppar->offset != 0) ivv = 1;
    fac = 0.5;
    if (nv1 < 0 || nv2 < 0) fac = -fac;
    for (iv = 0; iv < ndm2; ++iv) {
      isig = iv % nsymsq;
      for (ivp = 0; ivp < ndm2; ++ivp) {
        isigp = ivp % nsymsq;
        ioff = (isig - isigp) * (ndm2 - 1) + noff2;
        ix = iv + ndm2 * ivp;
        for (ift = 0; ift < nft; ++ift) {
          if (ift > 0) {
            ix += ndtot2;
          }
          va1 = xpec2[ix]; vb1 = xpec2[ix + ioff];
          if (ivv == 0) {
            xpec2[ix] = fac * (va1 + vb1);
          } else  {
            xpec2[ix] = fac * (va1 - vb1);
          }
        }
      }
    }
  }
  ix = 5 * ndtot2; ioff = ix;
  if (kdel > 0) ioff += noff2;
  for (iv = 0; iv < ndm2; ++iv) {
    for (ivp = 0; ivp <= iv; ++ivp) {
      ix1 = iv + ndm2 * ivp; ix2 = ivp + ndm2 * iv;
      va1 = xpec2[ix1+ix];   va2 = xpec2[ix2+ix];
      vb1 = xpec2[ix1+ioff]; vb2 = xpec2[ix2+ioff];
      if (fabs(va1-va2) > 1.e-10) break;
      if (fabs(vb1-vb2) > 1.e-10) break;
    }
    if (ivp <= iv) break;
  }
  return iret;
} /* get2top */
  
int wrpar(lu, id, nid, isdip, kflg, t, nt, labl, nfmt, drhok)
     FILE *lu;
     double *id;
     int kflg;
     double *t;
     int nt;
     int nid;
     int isdip;
     char *labl;
     int nfmt;
     double drhok;
{
  /* id = parameter id */
  /* isdip > 0 for dipole, else parameter */
  /* kflg = flag indicating type of output */
  /* t = vector of Fourier coefficients */
  /* labl = parameter label */
  /* nvfac2 = 100 or 10000 (size of vib field) */
  /* nft = number of Fourier coefficients */
  static double cmp;
  static char fmt0[] = "  %14.0f %23.15E %s\n";
  static char fmt1[] = "%14.0f%02d %23.15E %s\n";
  static char fmt2[] = "%12.0f %10.6f %23.15E %s\n";
  double fidc, finc, fidd, fida, fids, tmp, fid, tval, *tv, *tvp;
  int k, kk, nsym, n, i, ires, nret;
  long idtmp, idd;
  if (lu == NULL) return 0;
  fid = id[0];
  tval = t[0];
  nret = 0;
  if (kflg > 0) {
    /* print to output file format */
    if (kflg <= 100) {
      if (fabs(tval) > cmp) {
        fprintf(lu, fmt0, fid, tval, labl);
      }
    } else {
      nsym = kflg - 100;
      if (dnrm2(nt, t, 1) > cmp) {
        for (k = 0; k < nt; ++k) {
          kk = k;
          if ((k + k) > nfmt) kk -= nfmt; 
          tmp = kk * drhok; ++nret;
          fprintf(lu, fmt2, fid, tmp, t[k], labl);
        }
      }
    }
    if (isdip > 0) nret = 0;
    return nret;
  } else if (kflg == 0) {
    /* if kflg == 0 then print head type file and update comparison */
    cmp = (fabs(tval) + tiny) * tiny;
    if (isdip == 0) {
      ++nret;
      if (nfmt > 0)
        fprintf(lu, fmt1, fid, nfmt, tval, labl);
      else
        fprintf(lu, fmt0, fid, tval, labl);
    }
    return nret;
  }
  /* if kflg < 0 then print par file format */
  finc = (double) (nfmt / 100);
  n = nt - 1;
  if (isdip > 0) {            /* dipole moments */
    finc = 10. * finc;  /* increment */
    fidc = finc * 70.;  /* base value */
    fidd = finc * 10;   /* sin - cos */
    fida = finc * 90;   /* carry */
    if (n > 198)
      n = 198;
  } else {                      /* parameters */
    finc *= 1.e7;
    fidc = 0.;
    fidd = finc * 10; 
    fida = fidd;
    if (n > 98)
      n = 98;
  }
  tv = t;

  tmp = 0.;
  for (i = 0; i < nid; ++i) {
    fid = id[i]; tval = (*tv);
    if (fabs(tval) > cmp) {
      if (kflg == -2) { 
        fid = fabs(fid);
        kflg = -1;
      }
      fprintf(lu, fmt0, fid, tval, labl);
      if (isdip == 0) ++nret;
    } else {
      *tv = 0.;
    }
    if (n > 0)
      tmp += dnrm2(n, &tv[1], 1);
    tv += nt;
  }
  if (tmp < tiny * nid)
    return nret;
  k = 0;
  if (isdip > 0) 
    ++k;
  for (kk = 1; kk <= n; kk += 2) {
    if ((++k) > 10) {
      fidc += fida;
      k = 1;
    }
    fidc += finc;
    tv = &t[kk];
    for (i = 0; i < nid; ++i) {
      idtmp = (long) id[i]; idd = idtmp / 100;
      ires = (int)(idtmp - 100 * idd);
      if (ires < 0) ires = - ires; 
      fid = fidc + (double)idd; tval = (*tv);
      fids = fid + fidd;
      if (fabs(tval) > cmp) {
        if (kflg == -2) { 
          fid = fabs(fid);
          kflg = -1;
        }
        fprintf(lu, fmt1, fid, ires, tval, labl); ++nret;
      } else {
        *tv = 0.;
      }
      tvp = tv + 1;
      tval = (*tvp);
      if (fabs(tval) > cmp) {
        if (kflg == -2) { 
          fids = fabs(fids);
          kflg = -1;
        }
        fprintf(lu, fmt1, fids, ires, tval, labl); ++nret;
      } else {
        *tvp = 0.;
      }
      tv += nt;
    }
  }
  if (isdip > 0) nret = 0;
  return nret;
}                               /* wrpar */
