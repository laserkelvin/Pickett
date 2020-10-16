/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   24 March  1999: Fixed sin scaling in dircos */
/*   25 March  1999: Fixed diag = 2 error in call to ordham*/
/*   27 March  1999: Changed diag = 3 to take care of l-doubled states */
/*   30 March  1999: Fixed sin scaling bug */
/*   21 April  1999: Fixed bug at end of hamx for KROLL case */
/*   22 April  1999: Fixed bug in Itot option */
/*   25 May    1999: fix bug when largest number of spins is not default */
/*   23 July   1999: fix bug in SETOPT that returned incorrect IQNFMT */
/*    7 Sept.  1999: fix bug in GETQQ for that filled ivs with short */
/*   29 Dec.   1999: improve ORDHAM */
/*   12 April  2000: negate Hamiltonian for kappa > 0 on first pass */
/*   12 April  2000: separate ORDHAM and FIXHAM */
/*   12 April  2000: subtract constant from energy for more precision*/
/*    1 August 2000: special code for unit matrix */
/*   12 Sept.  2000: fix up phase for cross-multiplet interaction */
/*    9 Oct.   2000: change dipole type 10 to double commutator with N^2 */
/*   29 Nov.   2000: redefine idpars and spar to put more calc. in initial */
/*    5 Feb.   2001: define augmented Euler and Fourier series */
/*    7 July   2001: fix bug for l doublet pairs in setopt */
/*   26 Sept.  2001: fix bug in idpari for oblate 6x and 2x interactions */
/*   12 March  2002: improve selections for off-diagonal l in dircos */
/*   27 June   2002: increase vibration limit to 359 */
/*   29 July   2002: fix bug in fixham (DIAG = 2,3)*/
/*   18 Aug.   2003: code cleanup, @comment@ is for splint */ 
/*   24 Nov.   2003: implement code for higher symmetry IAX options */
/*    3 Feb    2004: implement code for higher symmetry spins */
/*   19 May    2004: implement code for bcd storage of idpar and idip */
/*   20 July   2004: fix bug in code for bcd storage of idip */
/*   19 Jan.   2005: extended Euler operator */
/*   28 May    2005: fix bug in getmask for operator between l-doubled state */
/*   27 Sept.  2005: fix bug for 6-fold symmetry B-weights in getqn & setopt */
/*   17 Oct.   2005: fix bug in getmask for operator within l-doubled state */
/*   12 Dec.   2005: fix bug Itot basis and include mixed basis */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "calpgm.h"
#include "spinit.h"
#define ODD(i)  (((int)(i) & 1) != 0)
#define EVEN(i) (((int)(i) & 1) == 0)
#define ODD2(i) (((int)(i) & 2) != 0)
#define MOD(i,n) ((int)(i) % (int)(n))
#define C0 '\0'
#define MAXINT  12  /* maximum number of dipole types */
#define MAXVIB  999
#define MAXSPIN 9
#define MAXN_DIRCOS 359 /* MAX K in DIRCOS */
#define MAXVDEC 3
#define NDECPAR 5   /* number of digit pairs for idpar not including vib */ 
#define NSBCD   (2 * (NDECPAR + MAXVDEC + 1))
#define MAXII   20  /* I <= 10 */
#define MAXNS   20  /* 2 * NSPINS + 2 */
#define MNSQ     2  /* bit for operator which is different from previous by
                             only N*(N+1) */
#define MCOS_OK  4  /* bit for operator which uses cosine from previous */
#define MNOUNIT  8  /* bit for operator which is not a unit matrix */
#define MMASK   (0x0fe0) /* clear bits 0-4 */
#define MIDEN   16  /* bit for Identity operator under alpha */
#define MODD  0x040 /* bit for imaginary operator */
#define MSYM2 0x080 /* bit for I_alpha - I_-alpha */
#define MLZ   0x100 /* bit for Lz operator */
#define MINOQ    4  /* bit for NO Delta N = 0 (commutator) */
#define MELEC    2  /* bit for electric dipole (must be 2)*/
#define NDXCOM  8   /* XCOM size, XCOM enum: */
#define XDIM    0   /* 0: DIMENSION OF SUB-BLOCK */
#define XSYM    1   /* 1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
#define XNVAL   2   /* 2: 'N' QUANTUM NUMBER */
#define XKBGN   3   /* 3: BEGINNING K */
#define XLVAL   4   /* 4: L VALUE */
#define XVIB    5   /* 5: VIBRATIONAL QUANTA */
#define XISYM   6   /* 6: SPIN SYMMETRY */
#define XIQN    7   /* 6: SPIN SYMMETRY QUANTUM NUMBER */

/* Common Declarations */
typedef struct {  /* vibrational information */
  /*@dependent@*/ short *spt;
  int nspstat;
  short knmin[4], knmax, wt[5], ewt[2], lvqn, lvupper, gsym, nqn;
} SVIB;
static SVIB vinfo1;
/*@owned@*/ static SVIB *vinfo = &vinfo1;

typedef struct str_spar {  /* local parameter information */
  /*@null@*/ /*@owned@*/ struct str_spar *next;
  double zfac;
  int ip;
  unsigned int ipsym;
  short flags;
  signed char ksq, fc, msi1, msi2, mldel;
  unsigned char njq, mln, mld, mkdel, mins, msznz, euler, alpha;
} SPAR;
typedef /*@null@*/ /*@owned@*/ SPAR *PSPAR;

static PSPAR spar_head[MAXVIB];

static short sptzero[] = { 1, 0 };      /* default for no spins */

typedef struct struct_ssp {
  /*@null@*/ /*@owned@*/ /*@reldef@*/ struct struct_ssp *next;
  /*@notnull@*//*@owned@*/ short *sspt;
  int ssize, nitot;
} SSP;
static SSP ssp_head = { NULL, sptzero, 1 , 0};

static struct {   /* save data for getqn */
  int cblk, cnblk, csblk, cff, cwt[5];
} *cgetq, cgetv[2];

typedef struct { /* local dipole information */
  double fac;
  short  flg;
  signed char kd, ld, ldel, fc;
} SDIP;
/*@owned@*/ static SDIP *dipinfo;

static double zero = 0.;
static double zwk = 0.;
static double spfac[MAXII];
static double spfac2[MAXII];
static int zmoldv = 0;
static int zblkptr = 0;
static int zivs = 0;
static int zipder = 0;
static int revsym[] = { 0, 3, 2, 1};
static int isoddk[] = { 0, 1, 1, 0};
static int is_esym[MAXITOT];
static int lscom[MAXNS], iscom[MAXNS], jscom[MAXNS], ismap[MAXNS];
static int ixcom[NDXCOM], jxcom[NDXCOM];
static int itptr, itsym, nspin, nsqmax, ndmx, ndmax, nddip;
static short szero = 0;
static short zidx = 0;
static short zjdx = 0;
static short zibkptr = 0;
static short zikmin = 0;
static struct {
  int mxspin, idiag, nvib, nfit, nbkpj, ixz, nqnn, nqn, maxqn, vibfac;
  int parinit, maxblk, nitot, vibdec, esymdec, msmask, nqn0, iqfmt0;
  unsigned int msshft;
  BOOL lsym, esym, oblate, vibfmt;
} glob;
static char sbcd[NSBCD];
/* pointers to dynamically allocated arrays */
/*@owned@*/ static int *moldv;
/*@owned@*/ static int *blkptr;
/*@owned@*/ static int *ipder, *ivs;
/*@owned@*/ static double *wk;
/*@owned@*/ static short *idx, *jdx, *ibkptr, *ikmin;

int dclr(const int n1, const int n2, double *vec, const int ix);
int specop(const int neuler, BOOL * newblk, int *nsqj, int *ikq,
                  const int ksi, const int ksj, const int ni, const int nj,
                  const int ncos, double *wk, const short *ix,
                  const short *jx, const double par);
int specfc(const int ifc, const int iv, const int jv,
                  const int kdel, const int ksi, const int ksj,
                  const int ncos, double *wk, const short *ix,
                  const short *jx);
int sznzfix(const int sznz, const int ni, const int nj, int *ixcom,
                   int *jxcom, int *iscom, int *jscom);
int sznzop(const int ni, const int nj, const int ksi, const int ksj,
                  const int *iscom, const int *jscom, const int ncos,
                  double *wk, const short *ix, const short *jx);
unsigned int blksym(const int *ixcom, const int *jxcom);
int ordham(const int nn, double *egy, const short *isblk,
                  short *iswap);
int fixham(const int ndm, const int nn, double *t, double *egy,
                  double *p, const short *iswap);
BOOL kroll(const int nsizd, double *t, const int nsblk,
                  const short *sbkptr, const short *kmin);
int getqs(const int mvs, const int iff, const int nsiz, const int kbgn, 
               /*@out@*/ int *ixcom, /*@out@*/ int *iscom, /*@out@*/ int *iv);
int idpars(SPAR * pspar, /*@out@*/ int *ksq, /*@out@*/ int *itp,
               /*@out@*/ int *l, /*@out@*/ int *ld, /*@out@*/ int *kdel,
               /*@out@*/ int *ins, /*@out@*/ int *si1, /*@out@*/ int *si2,
               /*@out@*/ int *sznz, /*@out@*/ int *ifc, /*@out@*/ int *alpha, 
               /*@out@*/ int *ldel);
int getll(const int llf, const int ld, const int kd, const int lt, 
          const int si1, const int si2, int *lscom, const int *iscom, 
          const int *jscom);
int getmask(const int *xbra, const int *xket, const int kd, const int ldel, 
                   const int loff, const int alpha);
double rmatrx(const int ld, const int lv, const int *ixcom,
                     const int *jxcom);
int symnsq(const int inq, const int ins, const int *iscom,
                  const int *jscom, double *z);
int symksq(const int ikq, const int ksi, const int ksj, const int n,
                  double *wk, short *ix, short *jx);
int dpmake(const int nsize, double *dp, const double *t,
                  const int n, const double *wk, const short *ix,
                  const short *jx, const BOOL isunit);
int pasort(FILE * lu, const int npar, bcd_t *idpar,
                  const double *par);
int idpari(bcd_t *idval, int itp, /*@out@*/ SPAR * pspar);
int checksp(const BOOL first, int si1, int si2, const short *iiv1, 
            const short *iiv2, double *zfac);
int setwt(SVIB * pvinfo, const int ivib, const int iax,
                 const int iwtpl, const int iwtmn, double vsym);
int getwt(SVIB * pvinfo, const int isym, const int iispin, 
          /*@out@*/ int *ivwt);
int setgsym(const int gsym);
int getsp(const bcd_t *ispnx, SVIB *pvinfo);
void setsp(void);
int getqq(const int iblk, /*@out@*/ int *f, /*@out@*/ int *iwtb,
                 /*@out@*/ short *sbkptr, /*@out@*/ short *kmin,
                 /*@out@*/ int *vs);
int dircos(const int *xbra, const int *xket, const int ld,
                  const int kd, const int ncmax, /*@out@*/ double *direl,
                  /*@out@*/ short *ibra, /*@out@*/ short *iket, 
                  const int ifup, const int loff, int mask, 
                  /*@out@*/ BOOL *isunit);
int ffcal(const int nff, const int kff, /*@out@*/ double *ff);
/*****************************************************************************/
int hamx(iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump)
const int iblk, nsize, npar;
const BOOL ifdump;
const bcd_t *idpar;
const double *par;
double *egy, *t, *dedp, *pmix;
{
  /* .. PACKAGE FOR 98 INTERACTING VIBRATIONAL STATES WITH MULTI-SPIN */

  /*     calculate energies and derivatives for block IBLK */

  /*     on entry: */
  /*         IBLK= block number */
  /*         NSIZE= block size and dimension */
  /*         IDPAR= list of parameter identifiers ( element 0 is length ) */
  /*         PAR = list of parameter values */
  /*         NPARD= number of parameters with positive IDPAR */
  /*     on return: */
  /*         EGY = energy list */
  /*         T = eigenvector matrix */
  /*         DEDP = derivative of energy with respect to parameter */

  static double sqj[10];
  static int nold = -1;
  SPAR *spar_now;
  short *itau;
  double *pt, *dbar, pbar, sqnn, zpar, dn, ele, egy0;
  long ndm;
  int ipar, ierr, ncos, ispn, jspn, nsqj, iiwt[5], ndmd, alpha, mkd;
  int i, ii, lt, n, ibase, jbase, kbgni, kbgnj, nsblk, idflags;
  int kd, ld, nd, ni, nj, ivbase, jvbase, ivmin, ifc, iz, jz, npair, ldel;
  int si1, si2, iff, ijd, ikq, ins, neuler, sznz, ixx, jxx, kl, isgn;
  unsigned int ivsym, lastvv, ivcmp;
  BOOL oldpar, isegy0, roll, parskp, isunit, isneg, first;
  BOOL newblk, firstpar;

  if (ndmx <= 0) {
    puts("working vectors not allocated");
    exit(EXIT_FAILURE);
  }
  dbar = &wk[ndmx];
  ndm = nsize;
  ndmd = nsize + 1;
  isneg = (glob.oblate && !ifdump);
  cgetv[0].cblk = 0;
  cgetv[1].cblk = 0;
  /*     get F and sub-block structure */
  itau = &ikmin[glob.maxblk];
  nsblk = getqq(iblk, &iff, iiwt, ibkptr, ikmin, ivs);
  /* zero hamiltonian matrix */
  dclr(nsize, nsize, t, 1);
  /* zero derivative matrix */
  n = glob.nfit;
  dclr(nsize, n, dedp, 1);
  /* set up hamiltonian */
  oldpar = parskp = roll = isunit = firstpar = FALSE;
  ncos = sznz = 0; isgn = 1;
  pbar = zpar = 0.;
  first = TRUE;
  do {
    /* loop on wang blocks */
    lastvv = 0; newblk = TRUE; egy0 = 0.;
    for (ixx = 0; ixx < nsblk; ++ixx) {
      ibase = ibkptr[ixx];
      n = ibkptr[ixx + 1] - ibase;
      kbgni = ikmin[ixx];
      ispn = ivs[ixx];
      getqs(ispn, iff, n, kbgni, ixcom, iscom, &ivbase);
      ni = ixcom[XNVAL];
      /* set up itau */
      ld = ixcom[XLVAL];
      if (ld == 0)
        ld = (ni + ixcom[XSYM]) & 1;
      else
        ld = (ld > 0)? 1: 0;
      kd = kbgni + kbgni - ld;
      kd = kd + kd;
      if (ivbase != ixcom[XVIB])
        ++kd; /* upper l-doubled state */
      itau[ixx] = (short) kd;
      if (ni != nold) {         /* set up powers of N(N+1) */
        nold = ni;
        dn = (double) ni;
        sqnn = dn * (ni + 1);
        sqj[0] = 1;
        for (i = 1; i <= nsqmax; ++i) {
          sqj[i] = sqj[i - 1] * sqnn;
        }
      }
      for (jxx = 0; jxx <= ixx; ++jxx) { /*  set up sub-block quanta */
        ijd = ixx - jxx;
        jbase = ibkptr[jxx];
        n = ibkptr[jxx + 1] - jbase;
        kbgnj = ikmin[jxx];
        jspn = ivs[jxx];
        getqs(jspn, iff, n, kbgnj, jxcom, jscom, &jvbase);
        nj = jxcom[XNVAL];
        nd = ni - nj;
        if (nd > ndmax) continue;
        ivmin = (ivbase < jvbase) ? ivbase : jvbase;
        ivsym = (unsigned int) (ivbase + jvbase + ivmin * glob.vibfac) << 2;
        ivcmp = ivsym + 3; ivsym += blksym(ixcom, jxcom);
        if (lastvv != ivcmp) { 
          newblk = TRUE; lastvv = ivcmp;
        }
        isegy0 = FALSE; firstpar = TRUE;
        for (spar_now = spar_head[ivmin]; TRUE; spar_now = spar_now->next) {
          /* loop over parameters */
          idflags = 0;
          if (spar_now != NULL) {
            ivcmp = spar_now->ipsym;
            if (ODD(spar_now->euler)) {
              if (ivcmp > lastvv) continue;
              if (ivcmp == lastvv)
                idflags = (int) spar_now->flags;
            } else {
              if (ivcmp > ivsym) continue; 
              if (ivcmp == ivsym)
                idflags = (int) spar_now->flags;
            }
          }
          if ((idflags & MCOS_OK) == 0) {       /* new cosine or end */
            if (oldpar) {
              /*  add diagonal or sub-diagonal to Hamiltonian */
              /*  using composite parameter value PBAR */
              if (isegy0) {
                egy0 += pbar;
              } else {
                if (isneg)
                  pbar = -pbar;
                if (isunit) {
                  iz = idx[0];
                  jz = jdx[0];
                  pt = &t[jz * ndm];
                  for (i = 0; i < ncos; ++i) {
                    pt[iz] += pbar;
                    iz += ndmd;
                  }
                } else {
                  for (i = 0; i < ncos; ++i) {
                    iz = idx[i];
                    jz = jdx[i];
                    t[iz + jz * ndm] += pbar * wk[i];
                  }
                }
              }
              oldpar = FALSE;
            }
            if (idflags == 0) break; /* end loop over parameters */
            ncos = -1;
          } else { /* check if old cosine operator was zero */
            if (ncos == 0) continue;
          }
          if (spar_now == NULL) break;
          ipar = spar_now->ip;
          nsqj = (int) spar_now->njq;
          if ((idflags & MNSQ) != 0 && nd == 0) {
            /* sum parameters which are different only by N*(N+1) */
            if (parskp) continue; /* jump if mother operator was zero */
          } else {                /*  new type of operator */
            parskp = TRUE;
            if (sznz < 0)
              sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
            kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins,
                        &si1, &si2, &sznz, &ifc, &alpha, &ldel);
            if (sznz > 0) {
              sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
              if (sznz == 0) continue;    /* get next parameter */
            }
            mkd = getmask(ixcom, jxcom, kd, ldel, kl, alpha);
            if (mkd == 0) continue;
            /* test tensor compatability */
            npair = getll(0, ld, lt, kd, si1, si2, lscom, iscom, jscom);
            if (npair < 0) continue;
            /* get direction cosines and K dependent part of operator */
            if (ncos < 0) {
              ncos = dircos(ixcom, jxcom, ld, kd, ndmx, wk, idx, jdx,
                            ijd, kl, mkd, &isunit);
              if (ncos <= 0) {
                if (ncos == 0) continue;     /* get next parameter */
                printf("DIRCOS WORKING VECTOR TOO SHORT IN HAMX BY %d\n",
                       ncos);
                exit(EXIT_FAILURE);
              }
              if ((idflags & MNOUNIT) != 0) 
                isunit = FALSE;
              isgn = 1;
              if (isunit && wk[0] < 0.)
                isgn = -1;
              if (neuler != 0) {        /*Euler operator */
                ncos = specop(neuler, &newblk, &nsqj, &ikq, kbgni, kbgnj,
                              ni, nj, ncos, wk, idx, jdx, par[ipar]);
                if (ncos == 0) continue;
              }
              if (ifc != 0)
                specfc(ifc, ivbase, jvbase, kd, kbgni, kbgnj, ncos, 
                       wk, idx, jdx);
              if (sznz != 0)
                sznzop(ni, nj, kbgni, kbgnj, iscom, jscom, ncos, wk, idx,
                       jdx);
              if (ikq != 0) {   /* find K*K dependence (anti-commutated) */
                ncos = symksq(ikq, kbgni, kbgnj, ncos, wk, idx, jdx);
                if (ncos == 0) continue;
              }
              if (ibase != 0) {
                for (i = 0; i < ncos; ++i)
                  idx[i] = (short) (idx[i] + ibase);
              }
              if (jbase != 0) {
                for (i = 0; i < ncos; ++i)
                  jdx[i] = (short) (jdx[i] + jbase);
              }
              if (first) {
                isegy0 = FALSE;
                pbar = 0.;
                if (isunit && ixx == jxx) {     /* diagonal unit operator */
                  oldpar = TRUE;
                  if (ixx == 0) {
                    isegy0 = TRUE;
                  } else if (firstpar){
                    pbar = -egy0;
                    firstpar = FALSE;
                  }
                }
                if (ipar < 0) continue;     /* get next parameter */
              } else { /*  begin by setting up reduced expectation values */
                if (ipar < 0) continue;     /* get next parameter */
                dpmake(nsize, dbar, t, ncos, wk, idx, jdx, isunit);
              }
            }
            /* end for: get direction cosines and K dependent parts */
            if (si2 < 0 && nd == 0) continue;  /* commutator with N * N */         
            /* find K independent parts */
            zpar = rmatrx(ld, lt, ixcom, jxcom);
            /*  spherical tensor transformations  */
            tensor(&zpar, iscom, jscom, lscom, ismap, npair, alpha);
            if (sznz < 0)
              sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
            /*  find N*(N+1) and N.S dependence */
            symnsq(nsqj, ins, iscom, jscom, &zpar);
            parskp = (fabs(zpar) < 1e-30);
            if (parskp)
              continue;         /* get next parameter */
            if (si2 < 0) { /* commutator with N * N / 2 */
              zpar *= ni;
            }
          }
          if (ni == 0 && nsqj != 0) continue;  /* get next parameter */
          ele = zpar * sqj[nsqj] * spar_now->zfac;
          if (isgn < 0)
            ele = -ele;
          if (first) {          /*  sum into composite parameter */
            pbar += ele * par[ipar] ;
            oldpar = TRUE;
          } else {              /* sum into derivative matrix */
            i = ipder[ipar];
            if (i < 0) {
              i = -1 - i;
              ele *= par[ipar];
              i = ipder[i];
            }
            daxpy(nsize, ele, dbar, 1, &dedp[i * ndm], 1);
          }
        } /* end parameter loop */
      } /* end JXX loop */
    } /* end IXX loop */
    if (first) {
      if (ifdump) {
        for (i = 0; i < nsize; ++i) {
          pt = &t[i + i * ndm];
          *pt += egy0;
          egy[i] = (*pt);
          pmix[i] = 0.;
        }
        return 0;
      }
      if (glob.idiag < 0 || nsize <= 1) {
        nd = nsize + 1;
        dcopy(nsize, t, nd, egy, 1);
        for (i = 0; i < nsize; ++i) {
          pt = &t[i * ndm];
          dcopy(nsize, &zero, 0, pt, 1);
          pmix[i] = pt[i] = 1.;
        }
      } else {
        if (glob.idiag == 0) {
          roll = kroll(nsize, t, nsblk, ibkptr, ikmin);
        } else if (glob.idiag == 2) {
          nd = nsize + 1;
          dcopy(nsize, t, nd, egy, 1);
          ordham(nsize, egy, ibkptr, jdx);
        }
        /* diagonalize */
        ierr = hdiag(nsize, nsize, t, egy, pmix);
        if (ierr != 0) {
          printf(" diagonalization failure, err = %d in block %d, %d\n",
                 ierr, iblk, nsize);
          exit(EXIT_FAILURE);
        }
        if (glob.idiag == 0) {  /*  energy sort of subblock */
          ordblk(nsize, nsize, t, egy, ibkptr, pmix, idx);
        } else if (glob.idiag == 1) {   /* projection sort of full block */
          for (i = 0; i <= nsize; ++i) {
            jdx[i] = (short) i;
          }
          ordblk(nsize, nsize, t, egy, jdx, pmix, idx);
        } else if (glob.idiag == 2) {
          ordblk(nsize, nsize, t, egy, ibkptr, pmix, idx);
          fixham(nsize, nsize, t, egy, pmix, jdx);
        } else if (glob.idiag == 3) {
          jxx = 0;
          ivbase = -1;
          for (i = 0; i < nsblk; ++i) {
            ixx = ibkptr[i];
            jvbase = (int) ((unsigned) ivs[i] >> 2) - (itau[i] & 1);
            if (ivbase != jvbase) {
              jdx[jxx++] = (short) ixx;
              ivbase = jvbase;
            }
            wk[ixx] = (double) (itau[i] / 2);
            nd = ibkptr[i + 1] - 1;
            for (ii = ixx; ii < nd; ++ii) {
              wk[ii + 1] = wk[ii] + 4;
            }
          }
          jdx[jxx] = ibkptr[nsblk];
          ordblk(nsize, nsize, t, egy, jdx, pmix, idx);
          ordham(nsize, wk, jdx, idx);
          fixham(nsize, nsize, t, egy, pmix, idx);
        }
      }
      if (isneg) {              /* correct for oblate ordering */
        for (i = 0; i < nsize; ++i)
          egy[i] = -egy[i];
      }
      for (i = 0; i < nsize; ++i)       /* add energy offset */
        egy[i] += egy0;
      first = FALSE;
    } else {                    /* not first */
      if (roll) {               /*  go back and get expectation values */
        dcopy(nsize, &zero, 0, egy, 1);
        for (ipar = npar - 1; ipar >= 0; --ipar) {
          i = ipder[ipar];
          if (i >= 0) {
            ele = par[ipar];
            daxpy(nsize, ele, &dedp[i * ndm], 1, egy, 1);
          }
        }
      }
      break;
    }
  } while (!first);  /* first loop */
  i = 1;
  if (idpar[0] == (bcd_t)0)
    i = -1;     /* this makes idpar referenced */
  return i;
}  /* hamx */

int dclr(n1, n2, vec, ix)
const int n1, n2, ix;
double *vec;
{    /*  clear a N1*N2 block */
  static long nbig = 0x7ff0;
  long nsq;
  double *pvec;

  pvec = vec;
  nsq = n1 * (long) n2;
  while (nsq > nbig) {
    dcopy((int) nbig, &zero, 0, pvec, ix);
    pvec += nbig;
    nsq -= nbig;
  }
  dcopy((int) nsq, &zero, 0, pvec, ix);
  return 0;
} /* dclr */

int specop(neuler, newblk, nsqj, ikq, ksi, ksj, ni, nj, ncos, wk, ix, jx, par)
const int neuler, ksi, ksj, ni, nj, ncos;
int *nsqj, *ikq;
BOOL *newblk;
double *wk;
const double par;
const short *ix, *jx;
{
  /*     subroutine to make Euler transform operators */
  /*     on entry: */
  /*         neuler= type of operator */
  /*         oldblk= flag for new block (=0 initially) */
  /*         NSQJ=power of N*N+N +1 */
  /*         IKQ= power of K*K */
  /*         KSI,KSJ= starting values of K */
  /*         NCOS= length of vector WK */
  /*         WK= vector of a sub-diagonal */
  /*         IX,JX= relative indices of elements in WK */
  /*     on return: */
  /*         WK= modified vector */
#define NSPOP 5
  static double aden[NSPOP], bden[NSPOP];
  double akisq, akjsq, anisq, anjsq, da, dd, di, dj, xi, xj;
  double sqi, sqj, di0, dj0;
  int i, nkq, nkq0, nnq, nnq0, ki, kj, ist;
  unsigned int ipwr;

  if (*newblk) {
    for (i = 0; i < NSPOP; ++i) {
      aden[i] = 0.;
      bden[i] = 0.;
    }
    *newblk = FALSE;
  }
  ist = (neuler - 2) >> 1;
  i = ist;
  if (i >= NSPOP) i = 0;
  nkq0 = (*ikq);
  nnq0 = (*nsqj);
  if (ODD(neuler)) {
    if (nnq0 == 0 && nkq0 == 1) {
      aden[i] = par;
    } else if (nnq0 == 1 && nkq0 == 0) {
      bden[i] = par;
    }
    return 0;
  }
  dd = bden[i]; da = aden[i] - dd;
  sqi = (double)(ni * (ni + 1)); di0 = 1. + dd * sqi;
  if (ni == nj) {
    sqj = sqi; dj0 = di0;
  }else {
    sqj = (double)(nj * (nj + 1)); dj0 = 1. + dd * sqj;
  }
  for (i = 0; i < ncos; ++i) {
    nkq = nkq0;
    nnq = nnq0;
    xi = wk[i];
    ki = ksi + ((int) ix[i] << 1);
    akisq = (double) (ki * ki);
    di = di0 + da * akisq;
    akisq /= di;
    kj = ksj + ((int) jx[i] << 1);
    if (ki == kj && ni == nj) {
      if (nnq > 0) {            /* N*N+N-K*K part */
        ipwr = (unsigned int) nnq;
        anisq = sqi / di - akisq;
        for (;;) {
          if (ODD(ipwr)) xi *= anisq;
          if ((ipwr >>= 1) == 0) break;
          anisq *= anisq;
        }
      }
      if (nkq > 0) {            /* K*K part */
        ipwr = (unsigned int) nkq;
        for (;;) {
          if (ODD(ipwr)) xi *= akisq;
          if ((ipwr >>= 1) == 0) break;
          akisq *= akisq;
        }
      }
      if (ist != 0)
        xi /= di;               /* extra part */
      wk[i] = xi;
    } else {
      akjsq = (double) (kj * kj);
      dj = dj0 + da * akjsq;
      akjsq /= dj;
      xj = xi;
      if (nnq > 0) {            /* N*N+N-K*K part */
        ipwr = (unsigned int) nnq;
        anisq = sqi / di - akisq;
        anjsq = sqj / dj - akjsq;
        for (;;) {
          if (ODD(ipwr)) {
            xi *= anisq;
            xj *= anjsq;
          }
          if ((ipwr >>= 1) == 0) break;
          anisq *= anisq;
          anjsq *= anjsq;
        }
      }
      if (nkq > 0) {            /* K*K part */
        ipwr = (unsigned int) nkq;
        for (;;) {
          if (ODD(ipwr)) {
            xi *= akisq;
            xj *= akjsq;
          }
          if ((ipwr >>= 1) == 0) break;
          akisq *= akisq;
          akjsq *= akjsq;
        }
      }
      if (ist != 0) {
        xi /= di;
        xj /= dj;
      }                         /* extra part */
      wk[i] = 0.5 * (xi + xj);
    }
  }
  *nsqj = 0;
  *ikq = 0;
  return ncos;
}                               /* specop */

int specfc(ifc, iv, jv, kdel, ksi, ksj, ncos, wk, ix, jx)
const int ifc, iv, jv, kdel, ksi, ksj, ncos;
double *wk;
const short *ix, *jx;
{
  /*     subroutine to make Fourier coefficient operators */
  /*     on entry: */
  /*         ITY= type of operator */
  /*         IV= v', JV=v" */
  /*         KDEL= change in K */
  /*         KSI,KSJ= starting values of K */
  /*         NCOS= length of vector WK */
  /*         WK= vector of a sub-diagonal */
  /*         IX,JX= relative indices of elements in WK */
  /*     on return: */
  /*         WK= modified vector */
  static double pfac[MAXVIB];
  static short ipfac[MAXVIB];
  double di, dj, ang;
  int i, n, ki, kj, kd;

  if (ifc == 0) {
    if (ncos == 0) {
      for (i = 0; i < MAXVIB; ++i) {
        pfac[i] = 0.;
        ipfac[i] = 0;
      }
    } else {
      /*  multiply by pi/3 */
      pfac[iv] = wk[0] * acos(0.5);
      if (pfac[iv] < 0.) {
        pfac[iv] = -pfac[iv];
        ipfac[iv] = 1;
      }
    }
    return 0;
  }
  n = ifc;
  if (n < 0)
    n = -1 - n;
  ang = (double) n;
  if (ODD(ipfac[iv] + ipfac[jv])) ang -= 0.5;
  di = ang * pfac[iv];
  dj = ang * pfac[jv];
  for (i = 0; i < ncos; ++i) {
    ki = ksi + (ix[i] << 1);
    kj = ksj + (jx[i] << 1);
    kd = ki - kj;
    if (kd < 0)
      kd = -kd;
    if (kdel != kd)
      kj = -kj;
    if (n != 0) {
      ang = di * ki + dj * kj;
      if (ifc < 0) {
        wk[i] *= sin(ang);
      } else {
        wk[i] *= cos(ang);
      }
    } else {
      wk[i] *= 0.5 * (ki + kj);
    }
  }
  return 0;
}                               /* specfc */

int sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom)
const int sznz, ni, nj;
int *ixcom, *jxcom, *iscom, *jscom;
{
  if (sznz <= 0) {
    if (sznz <= -4) {
      jscom[nspin] = nj << 1;
      jxcom[XNVAL] = nj;
    } else {
      iscom[nspin] = ni << 1;
      ixcom[XNVAL] = ni;
    }
    return 0;
  } else {
    switch (sznz) {
    case 1:
      if (ni + nj == 0)
        return 0;
      break;
    case 2:
      if (ni <= 1)
        return 0;
      iscom[nspin] -= 2;
      --ixcom[XNVAL];
      break;
    case 3:
      if (ni <= 0)
        return 0;
      iscom[nspin] += 2;
      ++ixcom[XNVAL];
      break;
    case 4:
      if (nj <= 1)
        return 0;
      jscom[nspin] -= 2;
      --jxcom[XNVAL];
      break;
    case 5:
      if (nj <= 0)
        return 0;
      jscom[nspin] += 2;
      ++jxcom[XNVAL];
      break;
    default:
      return 0;
    }
    return -sznz;
  }
}                               /* sznzfix */

int sznzop(ni, nj, ksi, ksj, iscom, jscom, ncos, wk, ix, jx)
const int ni, nj, ksi, ksj, ncos;
const int *iscom, *jscom;
double *wk;
const short *ix, *jx;
{
  /*     subroutine to add SzNz */
  /*     on entry: */
  /*         ISCOM,JSCOM= spins */
  /*         NCOS= length of vector WK */
  /*         WK= vector of a sub-diagonal */
  /*         IX,JX= relative indices of elements in WK */
  /*     on return: */
  /*         WK= modified vector */
  int iflg, i, k, kk, nni, nnj, nni0, nnj0, jj, iis;
  double vali, valj, twk, sum;

  if (ncos <= 0)
    return 0;
  nni = ni << 1;
  nni0 = iscom[nspin];
  vali = (double) (nni0 + 1);
  nnj = nj << 1;
  nnj0 = jscom[nspin];
  valj = (double) (nnj0 + 1);
  if (nni != nni0) {
    iflg = 1;
    vali = sqrt(vali * (nni + 1));
  } else if (nnj != nnj0) {
    iflg = -1;
    valj = sqrt(valj * (nnj + 1));
  } else {
    iflg = 0;
  }
  if (iflg <= 0) {
    jj = jscom[0];
    iis = jscom[nspin + 1];
    k = ksj + ksj + nnj0 + nnj + jj + iis;
    if ((k & 3) != 0)
      valj = -valj;
    valj *= c6jj(nnj0, 2, nnj, iis, jj, iis);
  }
  if (iflg >= 0) {
    jj = iscom[0];
    iis = iscom[nspin + 1];
    k = ksi + ksi + nni0 + nni + jj + iis;
    if ((k & 3) != 0)
      vali = -vali;
    vali *= c6jj(nni, 2, nni0, iis, jj, iis);
  }
  for (i = 0; i < ncos; ++i) {
    twk = wk[i];
    sum = 0.;
    if (iflg <= 0) {
      k = ksj + ((int) jx[i] << 1);
      kk = k << 1;
      sum = valj * twk * k * c3jj(nnj0, 2, nnj, -kk, 0, kk);
    }
    if (iflg >= 0) {
      k = ksi + ((int) ix[i] << 1);
      kk = k << 1;
      sum += vali * twk * k * c3jj(nni, 2, nni0, -kk, 0, kk);
    }
    wk[i] = sum;
  }
  return 0;
}                               /* sznzop */

unsigned int blksym(ixcom, jxcom)
const int *ixcom, *jxcom;
{                               /* get block symmetry from state symmetry */
  return (unsigned int) (3 & (ixcom[XSYM] ^ jxcom[XSYM]));
}                               /* blksym */

int ordham(nn, egy, isblk, iswap)
const int nn;
double *egy;
const short *isblk;
short *iswap;
{
/* subroutine to order eigenvalues within sub-block like diagonal of */
/*      Hamiltonian */
  double tmp;
  int iblk, i, ii, iq, is, inext;

  iblk = 1;
  i = 0;
  inext = isblk[1];
  for (is = 1; is < nn; ++is) {
    if (is == inext) {
      inext = isblk[++iblk];
    } else {                    /* find min energy */
      tmp = egy[i];
      iq = i;
      for (ii = is; ii < inext; ++ii) {
        if (egy[ii] < tmp) {
          tmp = egy[ii];
          iq = ii;
        }
      }
      iswap[i] = (short) iq;
      if (iq > i) {
        egy[iq] = egy[i];
        egy[i] = tmp;
      }
    }
    i = is;
  }
  return 0;
}                               /* ordham */

int fixham(ndm, nn, t, egy, p, iswap)
const int ndm, nn;
double *t, *egy, *p;
const short *iswap;
{
/* subroutine to order eigenvalues within sub-block like diagonal of */
/*      Hamiltonian using permutation found with ordham */
  double tmp, *tp;
  long ndml;
  int i, iq, is;

  is = nn - 2;
  ndml = ndm;
  if (is >= 0) {
    tp = &t[is * ndml];
    for (i = is; i >= 0; --i) {
      iq = iswap[i];
      if (iq > i) {
        tmp = egy[iq];
        egy[iq] = egy[i];
        egy[i] = tmp;
        tmp = p[iq];
        p[iq] = p[i];
        p[i] = tmp;
        dswap(nn, tp, 1, &t[iq * ndml], 1);
      }
      tp -= ndm;
    }
  }
  return 0;
}                               /* fixham */

BOOL kroll(nsizd, t, nsblk, sbkptr, kmin)
const int nsizd, nsblk;
double *t;
const short *sbkptr, *kmin;
{
/* subroutine to make sure that diagonal elements of the Hamiltonian */
/*    are monotonically increasing (decreasing for oblate) */
/*    after K=KTHRSH */
  double vall, tlast, tmp, *ptmp;
  long ndmt;
  int ibgn, iend, i, k, n, ixx;
  BOOL roll;

  ndmt = nsizd + 1;
  roll = FALSE;
  for (ixx = 0; ixx < nsblk; ++ixx) {
    i = 5 - kmin[ixx];
    ibgn = sbkptr[ixx];
    if (i >= 2)
      ibgn += i >> 1;
    iend = sbkptr[ixx + 1] - 1;
    if (ibgn >= iend) continue;
    ptmp = &t[ibgn * ndmt];
    tmp = tlast = *ptmp;
    ++ibgn;
    for (i = ibgn; i <= iend; ++i) {
      ptmp += ndmt;
      tlast = tmp;
      tmp = (*ptmp);
      if (tmp < tlast) break;
    }
    if (i > iend) continue;
    k = kmin[ixx] + ((i - sbkptr[ixx]) << 1);
    printf(" roll-over at K = %d\n", k);
    roll = TRUE;
    vall = 1e15;
    for (; i <= iend; ++i) {
      n = i - 1;
      dcopy(n, &zero, 0, &t[i], nsizd);
      n = nsizd - i;
      dcopy(n, &zero, 0, &t[i * ndmt], 1);
      tlast += vall;
      *ptmp = tlast;
      ptmp += ndmt;
    }
  }
  return roll;
}                               /* kroll */

int getqs(im, iff, nsiz, kbgn, ixcom, iscom, iv)
const int im, iff, nsiz, kbgn;
int *ixcom, *iscom, *iv;
{
  /*     subroutine to get quanta corresponding to sub-block */
  /*     on entry: */
  /*         IM= position in sub-block mold */
  /*         IFF = F quantum number *2 */
  /*         NSIZ= size of sub-block */
  /*         KBGN= first K value */
  /*     on return: */
  /*         IXCOM = list of rotational containing */
  /*                 0: DIMENSION OF SUB-BLOCK */
  /*                 1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
  /*                 2: 'N' QUANTUM NUMBER */
  /*                 3: BEGINNING K */
  /*                 4: L VALUE */
  /*                 5: VIBRATIONAL QUANTA */
  /*                 5: ITOT SYMMETRY */
  /*                 5: EXTRA ITOT QUANTA */
  /*         ISCOM = list of angular momenta */
  /*         IV = vibrational quantum number   if L=0,1 */
  /*         IV = vibrational quantum number-1 if L=-1 */
  /*     find list of angular momenta */
  SVIB *pvinfo;
  int *isscom;
  short *iis, *jjs;
  int i, iiv, nset, nspinv, ispx, k, lupper;
  unsigned int mvs;

  mvs = (unsigned int) im;
  ispx = (int) (mvs >> glob.msshft);
  iiv = (int) (mvs >> 2) & glob.msmask;
  pvinfo = &vinfo[iiv];
  iis = pvinfo->spt;
  nset = iis[0];
  jjs = &iis[nset * ispx]; i = (int)(mvs & 3);
  iscom[0] = iff + (*jjs);
  ixcom[XNVAL] = iscom[0] >> 1; /* find N */
  ixcom[XDIM] = nsiz;
  ixcom[XSYM] = i;        /* find symmetry */
  if (nsiz < 0)
    setgsym((int) pvinfo->gsym);
  k = kbgn;
  if (k < 0) {
    k = pvinfo->knmin[i]; 
    i = pvinfo->knmin[3 - i];
    if (k > i) k = i;
    if (k < 0) k = 0;
  }
  ixcom[XKBGN] = k;
  ixcom[XVIB] = iiv;
  ixcom[XISYM] = 0;
  lupper = pvinfo->lvupper;
  if (ODD(lupper)) --iiv;
  *iv = iiv;
  i = pvinfo->lvqn;
  ixcom[XLVAL] = i;
  if (i != 0 || k != 0) i = 1;
  ixcom[XIQN] = i;
  if (nspin == 0)
    return 1;
  nspinv = nset - 1;
  if (nspinv > nspin)
    nspinv = nspin;
  isscom = &iscom[nspin];
  isscom[0] = iscom[0];
  for (i = 1; i < nspinv; ++i) {
    iscom[i - 1] = iff + jjs[i];
    isscom[i] = iis[i];
  }
  iscom[i - 1] = iff;
  isscom[i] = iis[i];
  while (i < nspin) {
    ++i; iscom[i - 1] = iff; isscom[i] = 0;
  }
  if (glob.nitot > 0) {
    if (glob.nitot >= 3) {
      if (itsym < nspinv) {
        i = jjs[itsym + 1]; iscom[itsym] = i;
        k = MOD(i >> 2, glob.nitot);
        i -= k << 2;
        if (is_esym[k] != 0) ++i;
        ixcom[XISYM] = k;
        ixcom[XIQN] |= i;
      } else {
        iscom[itsym] = 0;
      }
    }
    iscom[itptr] = jjs[itptr + 1];
  }
  return ispx;
}                               /* getqs */

int idpars(pspar, ksq, itp, ln, ld, kdel, ins, si1, si2, sznz, ifc, 
           alpha, ldel)
SPAR *pspar;
int *ksq, *itp, *ln, *ld, *kdel, *ins, *si1, *si2, *sznz, *ifc, *alpha, *ldel;
{
/*     subroutine to parse parameter identifiers to subfields */
/*     on entry: */
/*         IDPAR= 10*( parameter ID/1000) + ITP */
/*     on return: */
/*         KSQ= power of K*K */
/*         ITP= parameter subtype */
/*         LN= tensor order */
/*         LD= order of direction cosine ( L <= LD ) */
/*         KDEL= change in K */
/*         INS=power of N.S */
/*         SI1,SI2=spin identifiers */
  int iret;
  *ksq = (int) pspar->ksq;
  *itp = (int) pspar->euler;
  *ln = (int) pspar->mln;
  *ld = (int) pspar->mld;
  *kdel = (int) pspar->mkdel;
  *ins = (int) pspar->mins;
  *si1 = (int) pspar->msi1;
  *si2 = (int) pspar->msi2;
  *sznz = (int) pspar->msznz;
  *ifc = (int) pspar->fc;
  *alpha = (int)pspar->alpha;
  *ldel = (int)pspar->mldel;
  iret = (int)pspar->flags;
  return iret;
}                               /* idpars */

int getll(llf, ld, ln, kd, si1, si2, lscom, iscom, jscom)
const int llf, ld, ln, kd, si1, si2;
int lscom[];
const int iscom[], jscom[];
{
  /*   function to find sequence of L values for coupling of spins */
  /*     on entry: */
  /*         LLF = 2 * tensor order for F */
  /*         LD = order of direction cosine */
  /*         LN = order of N dependent tensor */
  /*         SI1, SI2 = position of spin values */
  /*         ISCOM,JSCOM = list of 2 * angular momenta */
  /*     on return: */
  /*         LSCOM = list of 2 * tensor orders for succesive couplings */
  /*                 followed tensor orders for spins  */
  /*     return true if delta QN not consistent with L */
  int *lsscom;
  int i, llj, lld, llmax, jdif, jsum, maxspin, ierr, nm2;

  /* check tensor sum rules for N */
  ierr = -1; llmax = 0;
  lld = ld << 1;
  llj = ln << 1;
  jdif = iscom[nspin] - jscom[nspin];
  if (jdif < 0)
    jdif = -jdif;
  if (lld >= llj) {
    if (llj < jdif)
      return ierr;
    llmax = kd << 1;
    if (llmax < lld)
      llmax = lld;  
  } else {
    if (lld < jdif)
      return ierr;
    llmax = kd << 1;
    if (llmax < llj)
      llmax = llj;
  }
  jsum = iscom[nspin] + jscom[nspin];
  if (llmax > jsum)
    return ierr;
  if (nspin == 0) {
    lscom[0] = llj; 
    return 0;
  }
  maxspin = si1; 
  if (si1 > itsym) 
    maxspin = nspin;
  lsscom = &lscom[nspin]; lsscom[0] = llj;
  for (i = 1; i <= nspin; ++i) {
    lsscom[i] = 0;
  }
  if (si1 > 0) {
    lsscom[si1] = 2;
    if (si2 > 0)
      lsscom[si2] += 2;
  }
  /* generate list of tensor orders, checking quantum numbers along the way */
  nm2 = nspin - 2;
  for (i = 0; i <= nm2; ++i) {
    if (i == itsym) { /* special for Itot */
      if (si1 <= itsym) {
        llj = 0;
      } else if (si2 <= itsym || si1 == si2) {
        llj = lsscom[si1];
      } 
      i = nm2;
    } else if (i >= maxspin) {  /* operator complete */
      llj = llf;
    } else {
      llj -= lsscom[i + 1];     /* reduce order by spin operator order */
      if (llj < 0)
        llj = -llj;
    }
    lscom[i] = llj;
    jdif = iscom[i] - jscom[i];
    if (jdif < 0)
      jdif = -jdif;
    if (llj < jdif)
      return ierr;
    jsum = iscom[i] + jscom[i];
    if (llj > jsum)
      return ierr;
  }
  lscom[nspin - 1] = llf;
  if (llf != 0)
    maxspin = nspin;
  return maxspin;
}  /* getll */

int getmask(const int *xbra, const int *xket, const int kd, const int ldel, 
            const int loff, const int alpha)
{
  /* CREATE MASK FOR DIRCOS */
  /* KD IS TOTAL K QUANTIUM CHANGE */
  int ldif, lbra, lket, mask, kbra, kket, nitot, kk, ksbra, ksket;
  mask = 7;
  nitot = glob.nitot;
  if (nitot >= 3) { /* check Itot symmetry */ 
    ksbra = xbra[XISYM]; ksket = xket[XISYM];
    kk = (xbra[XIQN] ^ xket[XIQN]) >> 1;
    if ((loff & MIDEN) != 0) { /* Identity operator */
      if (kk != 0) return 0; 
    } else if (EVEN(kk)) { /* A1|A1, A1|E, E|E, A2|A2, ..*/
      if ((loff & MSYM2) != 0) return 0;
    } else { /* kk is ODD */
      if ((loff & MSYM2) == 0) return 0; /* <A2| I_alpha |E> */
      if (is_esym[ksbra] == 0 && is_esym[ksket] == 0) 
        return 0; /*  <A1 |I_alpha| A2> == 0 */
    }
    mask = 0;
    kk = ksket - ksbra;
    if (alpha == 0) {
      if (kk == 0) mask = 5;
    } else if (alpha < nitot) {
      if (MOD(kk + alpha, nitot) == 0)
        mask = 1; 
      if (MOD(kk - alpha, nitot) == 0)  
        mask += 4;
    } else { /* quasi-diagonal */
      if (MOD(ksbra + ksket - alpha, nitot) == 0)
        mask = 2;
    }
    if (mask == 0) return mask;
  }
  /* mask bit 0 (val = 1): kbra = kket + kd, lbra = lket + ldel */
  /* mask bit 2 (val = 4): kbra = kket - kd, lbra = lket - ldel */
  /* mask bit 1 (val = 2): kbra = kd - kket, lbra = ldel - lket */
  lbra = xbra[XLVAL]; lket = xket[XLVAL]; 
  ldif = lbra - lket;
  if (ldel == 0) { /* operator diagonal in l */
    if (ldif != 0) 
      mask &= 2; /* clear bits 0 and 2 */
  } else { /* operator off-diagonal in l */
    if (ldif == 0) {
      mask &= 2; /* clear bits 0 and 2 */
    } else if (ldif != ldel) { 
      mask &= 6; /* clear bit 0 */
    } else if (-ldif != ldel) { 
      mask &= 3; /* clear bit 2 */
    }
  }
  if ((mask & 2) == 0) return mask; /* bit 1 not set */
  mask &= 5;  /* clear bit 1 temporarily */
  /* quasi-diagonal, KD = KBRA + KKET */
  if (lbra + lket != ldel) return mask; 
  kbra = xbra[XKBGN]; kket = xket[XKBGN]; 
  kk = kd - kbra - kket;
  if (kk < 0) return mask;
  if (EVEN(xket[XIQN])) kk -= 2;
  if (EVEN(xbra[XIQN])) kk -= 2;
  if (kk < 0) return mask;
  return (mask + 2);
} /* getmask */

double rmatrx(ld, lv, ixcom, jxcom)
const int ld, lv;
const int *ixcom, *jxcom;
{                               /* find  N corrections to reduced matrix elements */
  double dlsq, dnsq, dn, del, tmp;
  int nsum, ndif, ni, nj, lmin, lt;

  tmp = 1.;
  if (ld == lv)
    return tmp;
  ni = ixcom[XNVAL];
  nj = jxcom[XNVAL];
  nsum = ni + nj + 1;
  dn = del = (double) nsum;
  dnsq = dn * dn;
  if (ld > lv) {                /* contracting tensor order */
    lmin = lv;
    lt = ld;
    if (lmin == 0) {
      lmin = 1;
      tmp = ni * (double) (ni + 1) / dn;
    }
  } else {                      /* expanding tensor order */
    lmin = ld;
    lt = lv;
    if (lmin == 0) {
      lmin = 1;
      tmp = ni * (double) (ni + 1) * dn;
    }
  }
  if (lt >= nsum)
    return 0.;
  ndif = ni - nj;
  if (ndif != 0)
    del = (double) (ndif * ndif);
  while (lt > lmin) {
    dlsq = (double) (lt * lt);
    tmp *= 0.25 * (dnsq - dlsq);
    if (ndif != 0)
      tmp -= tmp * del / dlsq;
    --lt;
  }
  return sqrt(tmp);
}                               /* rmatrx */

int symnsq(inq, ins, iscom, jscom, zval)
const int inq, ins, *iscom, *jscom;
double *zval;
{
/*  subroutine to find powers of N*(N+1) and N.S */
/*     on entry: */
/*         INQ= power of N*(N+1) */
/*         INS= power of N.S */
/*         ISCOM,JSCOM= list of angular momenta *2 */
/*         ISPN = S *2 */
/*         Z= matrix multipling factor */
/*     on return: */
/*         Z=  matrix factor */
  double x, dotns, dotnsp, zr, zl;
  int modf, ispn, nni, nnj, j;
  unsigned int ipwr;

  nni = iscom[nspin];
  nnj = jscom[nspin];
  modf = nni - nnj;
  if (inq <= 0)
    modf = 0;
  if (modf == 0 && ins == 0)
    return 0;
  zl = zr = 0.5 * (*zval);
  if (modf != 0) {              /*  CORRECT ELE FOR COMMUTATION */
    if (nnj == 0) {
      zr = 0.;
    } else {
      x = nnj * (double) (nnj + 2) / (nni * (double) (nni + 2));
      ipwr = (unsigned int) inq;
      for (;;) {
        if (ODD(ipwr)) zr *= x;
        if ((ipwr >>= 1) == 0) break;
        x *= x;
      }
    }
  }
  if (ins > 0) {
    j = iscom[0];
    ispn = iscom[nspin + 1];
    x = (double) (ispn * (ispn + 2));
    dotns = 0.125 * ((j - nni) * (double) (j + nni + 2) - x);
    j = jscom[0];
    ispn = jscom[nspin + 1];
    x = (double) (ispn * (ispn + 2));
    dotnsp = 0.125 * ((j - nnj) * (double) (j + nnj + 2) - x);
    ipwr = (unsigned int) ins;
    for (;;) {
      if (ODD(ipwr)) {
        zl *= dotns;
        zr *= dotnsp;
      }
      if ((ipwr >>= 1) == 0) break;
      dotns *= dotns;
      dotnsp *= dotnsp;
    }
  }
  *zval = zl + zr;
  return 0;
}                               /* symnsq */

int symksq(ikq, ksi, ksj, n, wk, ix, jx)
const int ikq, ksi, ksj, n;
double *wk;
short *ix, *jx;
{
/*     subroutine to multiply right and left by a power of K*K */
/*     on entry: */
/*         IKQ= power of K*K */
/*         KSI,KSJ= starting values of K */
/*         N= length of vector WK */
/*         WK= vector of a sub-diagonal */
/*         IX,JX= relative indices of elements in WK */
/*         note: if IKQ < 0 then operator is only define for K=-1-IKQ */
/*     on return: */
/*         WK= modified vector */
  double akisq, akjsq, xi, xj;
  int kcmp, i, ki, kj, zflg;
  unsigned int nkq;

  kcmp = -1 - ikq;
  zflg = 0;
  for (i = 0; i < n; ++i) {
    xi = wk[i];
    ki = ksi + ((int) ix[i] << 1);
    kj = ksj + ((int) jx[i] << 1);
    if (ikq > 0) {
      nkq = (unsigned int) ikq;
      if (ki != kj) {
        zflg = n;
        xj = xi;
        akisq = (double) (ki * ki);
        akjsq = (double) (kj * kj);
        for (;;) {
          if (ODD(nkq)) {
            xi *= akisq;
            xj *= akjsq;
          }
          if ((nkq >>= 1) == 0) break;
          akisq *= akisq;
          akjsq *= akjsq;
        }
        wk[i] = 0.5 * (xi + xj);
      } else if (ki > 0) {
        zflg = n;
        akisq = (double) (ki * ki);
        for (;;) {
          if (ODD(nkq)) xi *= akisq;
          if ((nkq >>= 1) == 0) break;
          akisq *= akisq;
        }
        wk[i] = xi;
      } else {
        wk[i] = 0.;
      }
    } else if (ki == kcmp && ki == kj) {
      ix[0] = ix[i];
      jx[0] = jx[i];
      wk[0] = xi;
      return 1;
    }
  }
  return zflg;
}                               /* symksq */

int dpmake(nsize, dp, t, n, wk, idx, jdx, isunit)
const int nsize, n;
const BOOL isunit;
double *dp;
const double *t, *wk;
const short *idx, *jdx;
{ /*  find derivative contribution from sub-diagonal */
  const double *pt;
  double ele, tz;
  int i, k, iz, jz;

  if (isunit) {
    iz = idx[0];
    jz = jdx[0];
    if (iz != jz) {             /* off-diagonal unit matrix */
      for (k = 0; k < nsize; ++k) {
        dp[k] = 2. * ddot(n, &t[iz], 1, &t[jz], 1);
        iz += nsize;
        jz += nsize;
      }
    } else {                    /* diagonal unit matrix */
      for (k = 0; k < nsize; ++k) {
        pt = &t[iz];
        dp[k] = ddot(n, pt, 1, pt, 1);
        iz += nsize;
      }
    }
  } else {
    dcopy(nsize, &zero, 0, dp, 1);
    for (i = 0; i < n; ++i) {
      iz = idx[i];
      jz = jdx[i];
      ele = wk[i];
      if (iz != jz) {
        ele *= 2;
        for (k = 0; k < nsize; ++k) {
          dp[k] += t[iz] * t[jz] * ele;
          iz += nsize;
          jz += nsize;
        }
      } else {
        for (k = 0; k < nsize; ++k) {
          tz = t[iz];
          dp[k] += tz * tz * ele;
          iz += nsize;
        }
      }
    }
  }
  return 0;
}                               /* dpmake */

int tensor(zval, iscom, jscom, lscom, imap, npair, alpha)
double *zval;
const int iscom[], jscom[], lscom[], imap[];
int npair, alpha;
{
/*     subroutine to find spherical tensor spin coupling coefficients */
/*     on entry: */
/*         Z= initial multiplier */
/*         NS= number of spins */
/*         ISCOM,JSCOM= list of angular momenta *2 */
/*         LLCOM= list of tensor orders *2 */
/*         LSCOM = list of tensor orders for spins *2 */
/*     on return: */
/*         Z= modified multiplier */
  double zsq, z;
  int isgn, i, jjf, jji, nnf, llj, nni, lln, lls, ssf, ssi, ix, jj[9];

  if (npair <= 0)
    return 0;
  z = (*zval); isgn = 0; zsq = 1.;  
  for (i = 0; i < npair; ++i) {
    ix = imap[i + i];
    lln = lscom[ix]; nni = iscom[ix]; nnf = jscom[ix];
    ix = imap[i + i + 1];
    lls = lscom[ix]; ssi = iscom[ix]; ssf = jscom[ix];
    if (alpha >= 0 && ix == itptr) {
      ix  = itsym + nspin + 1;
      getzitot(&z, lls, iscom[ix], &lscom[ix], 
               &iscom[itsym], &jscom[itsym], alpha, glob.nitot);
      i = nspin - 1;
    }
    llj = lscom[i]; jji = iscom[i]; jjf = jscom[i];
    if (lln != 0) {
      if (llj == 0) {           /* T * U */
        isgn += jji + ssi + nnf;
        z *= c6jj(nni, lln, nnf, ssf, jji, ssi);
      } else if (lls == 0) {    /* U = 1 */
        isgn += jjf + nni + ssi + llj;
        zsq *= jji + 1;
        zsq *= jjf + 1;
        z *= c6jj(nni, llj, nnf, jjf, ssi, jji);
      } else {                  /* T X U */
        isgn += 2;
        zsq *= 2.5;
        ix = (jjf + 1) * (llj + 1);
        zsq *= ix;
        zsq *= jji + 1;
        jj[0] = nnf;
        jj[1] = lln;
        jj[2] = nni;
        jj[3] = ssf;
        jj[4] = lls;
        jj[5] = ssi;
        jj[6] = jjf;
        jj[7] = llj;
        jj[8] = jji;
        z *= c9jj(jj);
      }
    } else if (lls != 0) {      /* T = 1 */
      isgn += jji + nni + ssf + llj;
      zsq *= jji + 1;
      zsq *= jjf + 1;
      z *= c6jj(ssi, llj, ssf, jjf, nni, jji);
    }
  }
  if (ODD2(isgn)) z = -z;
  *zval = z * sqrt(zsq);
  return 0;
}  /* tensor */

int setint(lu, ifdiag, nsav, ndip, idip, isimag)
FILE *lu;
BOOL *ifdiag;
const int ndip;
int *nsav, *isimag;
bcd_t *idip;
{   /*   subroutine to initialize intensity data */
  SVIB *pvib1, *pvib2;
  short *iiv1, *iiv2;
  double zfac;
  unsigned int ijv;
  int i, icase, isym, jsym, ivdif, ii, iv1, iv2, lv1, lv2, ifac, lv, nsym;
  int si1, ibcd, ndecv, nbcd, ipty, iflg, kd, ld, ifc, imag, ldel, jj;
  bcd_t itmp;
  BOOL bad_dip;
  *ifdiag = ((glob.idiag + 1) != 0);
  *nsav = 1;
  free(dipinfo);
  dipinfo = (SDIP *) mallocq((size_t) ndip * sizeof(SDIP));
  nddip = ndip;
  dipinfo[0].fac  = 1.;
  dipinfo[0].flg  = (short) 0;
  dipinfo[0].kd   = C0;
  dipinfo[0].ld   = C0;
  dipinfo[0].ldel = C0;
  dipinfo[0].fc   = C0;
  ifac = glob.vibfac + 1; 
  ndecv = glob.vibdec; nbcd = (int) idip[0] & 0x7f;
  for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd) {
    ii = (int)idip[ibcd + ndecv + 1];
    si1 = (ii >> 4) & 0x0f;
    ii &= 0x0f;
    if (ndecv > 2) 
      ii = 100 * ii + bcd2i(idip[ibcd + 3]);
    if (ndecv > 1)
      ii = 100 * ii + bcd2i(idip[ibcd + 2]);
    itmp = idip[ibcd + 1];
    isym = (int) itmp & 0x0f;
    ii = 10 * ii + ((int)(itmp >> 4) & 0x0f);
    iv2 = ii / ifac;
    iv1 = ii - iv2 * ifac;
    icase = bcd2i(idip[ibcd + ndecv + 2]);
    ivdif = iv1 - iv2;
    if (ivdif < 0) {
      ii = iv1; iv1 = iv2; iv2 = ii;
    }
    pvib1 = &vinfo[iv1]; pvib2 = &vinfo[iv2];
    lv1 = pvib1->lvupper; lv2 = pvib2->lvupper;
    if (ODD(lv1)) --iv1;
    if (ODD(lv2)) --iv2;
    ipty = ((lv1 ^ lv2) & 2) >> 1;
    lv1 = pvib1->lvqn; lv2 = pvib2->lvqn;
    /* MAKE SURE VIBRATIONS ARE ORDERED */
    bad_dip = TRUE; ifc = 0; zfac = 1.;
    do { /* one pass */
      /*  CONVERT SYMMETRY */
      if (isym > 3 || icase > MAXINT) break;
      if (iv1 >= glob.nvib)  break;
      if (iv2 >= glob.nvib)  break;
      if ((lv1 + lv2) < 0) break;
      if (pvib1->spt    != pvib2->spt) break;
      if (pvib1->gsym   != pvib2->gsym) break;
      if (pvib1->ewt[0] != pvib2->ewt[0]) break;
      if (pvib1->ewt[1] != pvib2->ewt[1]) break;
      if (glob.oblate == FALSE)
        isym = revsym[isym];      
      iflg = 1;
      if (isym == 0) {
        ld = 0; zfac = 0.009274;
      } else {
        ld = 1; iflg |= MELEC;
      }
      imag = 0;
      kd = isoddk[isym]; /* ld = 0,1; kd = 0, 1, 1, 0 */
      if (icase > 0) {
        switch (icase) {
        case 1: /* ld = 2,2; kd = 0, 1, 1, 2 */
          ld = 2;
          if (isym == 3)
            kd = 2;
          break;
        case 2: /* ld = 2,2; kd = 2, 1, 1, 2 */
        case 11:
        case 12:
          ld = 2;
          if (kd == 0)
            kd = 2;
          if (icase == 11) {
            iflg |= MINOQ; ++imag; /* commutator */ 
          } else if (icase == 12) {
            ifc = -1;
            ++imag; isym = 3 - isym;
          }
          break;
        case 5: 
          if (isym != 0) {
            iflg |= MINOQ; ++imag; /* commutator */
          }
          break;
        case 6: /* ld = 0,3; kd = 2, 3, 3, 2 */
          if (isym != 0)
            ld = 3;
          kd += 2;
          break;
        case 7:
        case 8:
          ifc = icase * 10 + si1 - 70;
          if (ifc > 10) {  
            ifc = 9 - ifc;
            ++imag; isym = 3 - isym;
          }
          si1 = 0;
          break;
        case 10:
          iflg |= MINOQ; /* double commutator */ 
          break;
        }
        if ((iflg & MINOQ) != 0 && isym == 0) break;
      }
      if (icase == 11 && isym == 0) break;
      for (ii = 0; ii < 4; ++ii) {
        jsym = ii ^ isym; 
        if (pvib1->wt[ii] != pvib2->wt[jsym]) { 
          jj = pvib1->wt[4];
          if (jj == 0) break;
          if (pvib1->wt[ii] != pvib2->wt[jsym ^ jj]) break;
        }
      }
      if (ii < 4) break;
      if (pvib1->wt[4] != pvib2->wt[4]) break;
      imag += ld;
      ldel = lv1 - lv2;
      nsym = setgsym((int)pvib1->gsym);
      if (ldel != 0) {
        if (icase == 5 && isym == 0) break;
        if (ldel < 0)
          ldel = -ldel;
        if (ivdif < 0)  /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
      } else if (icase == 5 && isym == 0) {
        if (lv1 == 0 || si1 > 0) break;
        iflg |= MLZ; ld = 1;
      }
      if (ODD(imag)) {
        iflg |= MODD;
      } else if (isym > 0) {
        if (iv1 == iv2 && ldel == 0) break;
      }
      if (glob.nitot >= 3) {
        if (MOD(kd - ldel, glob.nitot) != 0) 
          break; /* require alpha == 0 */
        if (si1 <= itsym)  
          iflg |= MIDEN; /* identity operator */
      }
      lv = 1;
      if (si1 > 0) {
        /* check spin properties */
        if (ld == 0 && isym != 0) break;
        iiv1 = pvib1->spt; iiv2 = pvib2->spt;
        if (checksp(TRUE, si1, 0, iiv1, iiv2, &zfac) != 0) break;
        lv = iflg & MELEC; /*  simple electric dipole, lv = 2 */
      }
      if (ld == kd && ld < lv) 
        ld -= (ld - lv) & 0x7fe; /* subtract even number */
      isimag[i] = (isym > 0)? 1: 0;
      dipinfo[i].fac  = zfac;
      dipinfo[i].flg  = (short) iflg;
      dipinfo[i].kd   = (signed char) kd;
      dipinfo[i].ld   = (signed char) ld;
      dipinfo[i].ldel = (signed char) ldel;
      dipinfo[i].fc   = (signed char) ifc;
      bad_dip = FALSE;
    } while (FALSE);
    if (bad_dip) {
      putbcd(sbcd, NSBCD, &idip[ibcd]);
      fprintf(lu,
              " WARNING: dipole %3d %s has no matrix elements\n",
        (i + 1), sbcd);
      dipinfo[i].fac  = 0.;
      dipinfo[i].flg  = (short)0;
      dipinfo[i].kd   = C0;
      dipinfo[i].ld   = C0;
      dipinfo[i].ldel = C0;
      dipinfo[i].fc   = C0;
      iv2 = iv1 = glob.vibfac; isym = 0;
    }
    ii = iv2 * ifac + iv1;
    ijv = ((unsigned int) ii << 2) + (unsigned int) isym;
    if (ndecv < 3) {
      idip[ibcd + 5] = idip[ibcd + ndecv + 2];  
      idip[ibcd + 4] = idip[ibcd + ndecv + 1]; 
      idip[ibcd + 3] = (bcd_t)0;
    } else {
      idip[ibcd + 4] &= 0xf0;
      idip[ibcd + 3] = (bcd_t) ((ijv >> 16) & 0xff);
    }
    idip[ibcd + 2] = (bcd_t) (ijv >> 8);
    idip[ibcd + 1] = (bcd_t) ijv;
    if (si1 == 0) 
      idip[ibcd + 4] &= 0x0f;
  }
  return 0;
} /* setint */


int intens(iblk, isiz, jblk, jsiz, ndip, idip, dip, s)
const int iblk, isiz, jblk, jsiz, ndip;
const bcd_t *idip;
const double *dip;
double *s;
{
  /*  function to calculate intensities */
  /*     on entry: */
  /*         IBLK,JBLK= block numbers */
  /*         ISIZ,JSIZ= block sizes */
  /*         IDIP= dipole id codes , IDIP(0)= number of dipoles */
  /*         DIP= dipole values */
  /*     on return: */
  /*         S= transition dipole matrix in original basis */
  /*         INTENS returns true is S is not 0 */
  static int idipoff;
  SDIP *pdip;
  double dd;
  long ndms;
  unsigned int ijv;
  int ncos, ifup, isym, i, k, lv, n, icase, ibase, mkd, ibcd, nbcd;
  int jbase, kbgni, kbgnj, nblki, nblkj, ixtmp, iret, kd, ld, npair, ldel;
  int ix, jx, nx, si1, si2, iff, jff, nni, nnj, iwt[5], jwt[5], ndecv;
  int ivv, jvv, ixx, jxx, ifc, ksym, kl, ioff, joff, alpha, nitot, dipoff;
  BOOL isunit;
  bcd_t bijv1, bijv2, bijv3;

  if (ndmx <= 0) {
    puts("working vectors not allocated");
    exit(EXIT_FAILURE);
  }
  dipoff = idipoff;
  if (dipoff >= nddip) 
    dipoff = 0;
  idipoff = dipoff + ndip;
  iret = 0;
  cgetv[0].cblk = 0;
  cgetv[1].cblk = 0;
  /*     get quantum information */
  nblki = getqq(iblk, &iff, iwt, ibkptr, ikmin, ivs);
  ioff = glob.maxblk;
  nblkj = getqq(jblk, &jff, jwt, &ibkptr[ioff], &ikmin[ioff], &ivs[ioff]);
  ixx = iff - jff;
  if (ixx < -2 || ixx > 2) return 0;
  ixx = iff + jff;
  if (ixx < 2 || ODD(ixx)) return 0;
  if (iwt[0] != jwt[0]) return 0;
  if (iwt[3] != jwt[3]) return 0;
  nitot = glob.nitot;
  alpha = 0;
  bijv3 = (bcd_t) 0;
  /* clear dipole matrix */
  dclr(isiz, jsiz, s, 1);
  /* loop over sub-blocks */
  ndms = isiz; ndecv = glob.vibdec; nbcd = (int) idip[0] & 0x7f;
  for (ixx = 0; ixx < nblki; ++ixx) {
    ibase = ibkptr[ixx];
    n = ibkptr[ixx + 1] - ibase;
    kbgni = ikmin[ixx];
    ixtmp = ivs[ixx];
    getqs(ixtmp, iff, n, kbgni, ixcom, iscom, &ivv);
    nni = iscom[0];
    for (jxx = 0; jxx < nblkj; ++jxx) {
      joff = jxx + ioff;
      jbase = ibkptr[joff];
      n = ibkptr[joff + 1] - jbase;
      kbgnj = ikmin[joff];
      ixtmp = ivs[joff];
      getqs(ixtmp, jff, n, kbgnj, jxcom, jscom, &jvv);
      nnj = jscom[0];
      nx = ((nni > nnj) ? nni : nnj) >> 1;
      i = (ivv < jvv) ? ivv : jvv;
      ijv = (unsigned int) (ivv + jvv + i * glob.vibfac);
      ijv = (ijv << 2) + blksym(ixcom, jxcom);
      bijv1 = (bcd_t) ijv;
      bijv2 = (bcd_t)(ijv >> 8);
      if (ndecv == 3)
        bijv3 = (bcd_t)(ijv >> 16);
      for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd) {
        /*  set dipole for transition */
        if (bijv1 != idip[ibcd + 1]) continue;
        if (bijv2 != idip[ibcd + 2]) continue;
        if (bijv3 != idip[ibcd + 3]) continue;
        pdip = &dipinfo[dipoff + i];
        kl = pdip->flg;
        if ((kl & MINOQ) != 0 && nni == nnj) continue;
        isym = (int)(bijv1 & 3);
        si1 = (int) (idip[ibcd + 4] >> 4) & 0x0f; si2 = 0;
        icase = bcd2i(idip[ibcd + 5]);
        kd = (int) pdip->kd; ifc = (int) pdip->fc; ldel = (int) pdip->ldel;
        lv = 1; ld = (int) pdip->ld;
        if (si1 > 0) 
          lv = kl & MELEC;
        ksym = 0;
        if (nitot >= 3 && kd > 0) {
          alpha = nitot; ksym = 1;
        }
        do {
          if (ksym == 0) { 
            alpha = 0;
          }
          /*  find matrix elements */
          mkd = getmask(ixcom, jxcom, kd, ldel, kl, alpha);
          if (mkd == 0) continue;
          npair = getll(2, ld, lv, kd, si1, si2, lscom, iscom, jscom);
          if (npair < 0) continue;
          ifup = 1 - (kl & MELEC);
          ncos = dircos(ixcom, jxcom, ld, kd, ndmx, wk, idx, jdx,
                        ifup, kl, mkd, &isunit);
          if (ncos == 0) continue;
          if (ncos < 0) {
            puts("DIRCOS WORKING VECTOR TOO SHORT IN INTENS");
            exit(EXIT_FAILURE);
          }
          dd = dip[i] * pdip->fac;
          /* correct for special cases */
          switch (icase) {
          case 3:
            if (nni == nnj) {
              k = nx * (nx + 1);
            } else {
              k = nx * nx;
            }
            dd *= k;
            break;
          case 4:
            symksq(1, kbgni, kbgnj, ncos, wk, idx, jdx);
            break;
          case 5:
          case 11:
            if (ifup > 0) break;
            if (nni < nnj)
              dd = -dd;
            dd *= nx;
            break;
          case 9:
            if (nni == nnj) {
              k = nx * (nx + 1);
              dd *= k;
            } else {
              k = nx * nx;
              dd *= k + 1;
            }
            dd *= k;
            break;
          case 10:
            dd *= nx * nx;
            break;
          case 7:
          case 8:
          case 12:
            specfc(ifc, ivv, jvv, kd, kbgni, kbgnj, ncos, wk, idx, jdx);
            break;
          }
          /*  correct for reduced matrix of N */
          if (ld != lv)
            dd *= rmatrx(ld, lv, ixcom, jxcom);
    /*  couple dipoles through the spins */
          tensor(&dd, iscom, jscom, lscom, ismap, npair, alpha);
          for (k = 0; k < ncos; ++k) {
            ix = idx[k] + ibase;
            jx = jdx[k] + jbase;
            s[ix + jx * ndms] += dd * wk[k];
          }
          iret = iff + 1;
        } while (--ksym >= 0);
      }
    }
  }
  return iret;
}                               /* intens */

int setopt(luin, nfmt, itd, nbcd, namfil)
FILE *luin;
int *nfmt, *itd, *nbcd;
char *namfil;
{
  /* this subroutine reads in quantum quantities */
  /*     on entry: */
  /*         LUIN is the fortran unit number to read from */
  /*         NFMT is the maximum number of quanta */
  /*     on return: */
  /*         NFMT= number of quantum number formats for catalog */
  /*         ITD= 2 for linear, 3 for non-linear */
  /*         NAMFIL contains a file name on which to find parameter names */
  /*         SETOPT = number of cards read, -1 for end-of-file */
#define NVEC 10
#define NDECSPIN 11
#define MAXIAX 11
  static int isymv[] ={2,2,2,2,2,2,3,4,4,5,6,6};
  static int iaxv[]  ={2,1,2,3,4,5,2,4,5,2,4,5};
  SVIB *pvinfo, *pvinfom;
  double rvec[NVEC], rvec0[NVEC], vsym;
  size_t nl;
  int ivib, lopt, i, ii, j, k, nvibm, isym, iend, ivsym;
  int iwtpl, iwtmn, nt, ncards, knnmin, knnmax, si2, iax;
  int spinsgn, ewt0, ewt1, nspinqn;
  int ivwt[3];
  bcd_t bcdspin[NDECSPIN];
  char card[82], ctyp;

  cgetv[0].cblk = cgetv[1].cblk = 0;    /* reset store for getqn */
  cjjini();
  glob.lsym = TRUE; glob.esym = TRUE; glob.esymdec = 100;
  rvec0[0] = 1; 
  rvec0[1] = 0;
  rvec0[2] = (double)MAXN_DIRCOS;
  rvec0[3] = 0;
  rvec0[4] = 1;
  rvec0[5] = 1;
  rvec0[6] = 1;
  rvec0[7] = 0;
  rvec0[8] = 99;
  rvec0[9] = 0;
  ncards = nvibm = nt = si2 = nspinqn = 0;
  *itd = 2;
  bcdspin[0] = (bcd_t)NDECSPIN;
  /*     read option cards */
  do {
    if (fgetstr(card, 82, luin) <= 0)
      return -1;
    ++ncards;
    ctyp = C0;
    if (isalpha((int) card[0])) 
      ctyp = card[0];
    else if (isalpha((int) card[1]))
      ctyp = card[1];
    iend = getbcd(card, bcdspin, NDECSPIN);
    spinsgn = NEGBCD(bcdspin[0]);
    dcopy(NVEC, rvec0, 1, rvec, 1);
    rvec[7] = 0.;
    pcard(&card[iend], rvec, NVEC, NULL);
    lopt = (int) rvec[0];
    knnmin = (int) rvec[1];
    knnmax = (int) rvec[2];
    iax = (int) rvec[4];
    iwtpl = (int) rvec[5];
    iwtmn = (int) rvec[6];
    vsym = rvec[7];
    ewt0 = (int) rvec[8]; ewt1 = 0;
    if (vsym < -0.5) {
      ivsym = -1; vsym = 0.;
    } else if (vsym > 0.5 && vsym < 2.e+15) {
      ivsym = 1; 
    } else {
      ivsym = 0; vsym = 0.;
    }
    ivib = lopt;
    if (ivib < 0)
      ivib = -ivib;
    if (knnmax < 0 || knnmax > MAXN_DIRCOS)
      knnmax = MAXN_DIRCOS;
    if (knnmin > knnmax)
      knnmax = knnmin;
    if (knnmin != knnmax)
      *itd = 3;
    if (ncards == 1) {
      dcopy(NVEC, rvec, 1, rvec0, 1);
      strcpy(namfil, "sping.nam");
      if (ctyp != C0)
        namfil[4] = ctyp;
      glob.ixz = (int) rvec[3];
      glob.idiag = (int) rvec[9];
      if (glob.idiag > 3) {
        glob.idiag = 0;
      } else if (glob.idiag < 0) {
        glob.idiag = -1;
      }
      if (ivib > 80) {
        if (ivib > MAXVIB)
          ivib = MAXVIB;
        if (ivib > 80 && sizeof(int) == 2)
          ivib = 80;
      } else if (ivib <= 0) {
        ivib = 1;
      }
      if (glob.nvib > 1) {
        free(vinfo);
        vinfo = NULL;
      }
      vinfo1.spt = sptzero;
      if (ivib > 1) {
        nvibm = ivib - 1;
        nl = (size_t) ivib * sizeof(SVIB);
        vinfo = (SVIB *) mallocq(nl);
        memcpy(vinfo, &vinfo1, sizeof(SVIB));
      } else {
        vinfo = &vinfo1;
      }
      for (k = 0; k <= 4; ++k)
        vinfo->wt[k] = 0;
      vinfo->ewt[0] = vinfo->ewt[1] = 0;
      glob.nvib = ivib;
      if (ivib <= 9) {
        glob.vibfac =   9; glob.vibdec = 1; k = 4; /* 16 */
      } else if (ivib <= 99) {
        glob.vibfac =  99; glob.vibdec = 2; k = 7; /* 128 */
      } else {
        glob.vibfac = 999; glob.vibdec = 3; k = 10; /* 1024 */
      }
      glob.msshft = (unsigned int)(k + 2);
      /* msshft = shifts needed for vib and symmetry */
      glob.msmask = (int)(1 << (unsigned int) k) - 1;
      ivib = 0;
      glob.oblate = (lopt < 0);
      glob.nqnn = 3;
      if (spinsgn != 0)
        glob.nqnn = 1;
      if (ewt0 < 0) glob.esymdec = 1000;
    } else if (ivib >= glob.nvib) {
      /* if IVIB out of range read another card */
      continue;
    }
    pvinfo = &vinfo[ivib];
    k = iax;
    if (iax < 0)
      iax = -iax;
    if (iax > MAXIAX) 
      iax = 0;
    if (iax <= 3 && glob.oblate)
      iax = revsym[iax];
    isym = isymv[iax];
    iax = iaxv[iax];
    pvinfo->gsym = (short) (isym << 1);
    /*  set up minimum K values */
    pvinfo->knmax = (short) knnmax;
    pvinfo->knmin[0] = 0; 
    pvinfo->knmin[1] = 1;
    pvinfo->knmin[2] = 1;
    pvinfo->knmin[3] = 2;
    if (k < 0) { /* iax was negative */
      if (isym >= 3) {
        pvinfo->knmin[0] = -2;
        pvinfo->knmin[3] = -2;
      }
      pvinfo->gsym |= 1;
    }
    if (knnmax != 0 && glob.nqnn == 1)
      glob.nqnn = 2;
    /* set ewt */
    if (ewt0 < 0) 
      ewt0 = -ewt0;
    if (ewt0 >= glob.esymdec) {
      ewt1 = ewt0 / glob.esymdec; ewt0 -= glob.esymdec * ewt1;
    }
    if (isym < 3) { 
      ewt0 = glob.esymdec - 1; glob.esym = FALSE;
    } else if (isym == 4 && iax == 5) {/* 4-fold B entry has no E */
      ewt0 = -1;
    }
    if (isym > 3 && pvinfo->wt[4] < 4) {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short) 0;
      pvinfo->wt[4] |= 4;
    }
    if (isym == 6) { /* six-fold, iax = 4,5 */
      pvinfo->ewt[iax-4] = (short) ewt0;
    } else  if (ewt0 >= 0) {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short) ewt0;
    }
    pvinfo->lvupper = 0;
    if (ewt1 != 0) {
      if (ewt1 > 10) {
        pvinfo->lvupper = 2; 
        ewt1 -= 10;
      }
      if (isym >= 3 && ewt1 > isym)
        ewt1 = MOD(ewt1 - 1, isym) + 1;
      pvinfo->lvupper += (short) (ewt1 << 2);
      pvinfo->knmin[0] = 0;
      pvinfo->knmin[3] = 0; /* 0,1,1,0 */
    }
    pvinfo->lvqn = (short) ewt1; 
    for (i = 0; i <= 3; ++i) {
      k = pvinfo->knmin[i];
      if (k < 0 && knnmin == 0) 
        continue;
      while (k < knnmin)
        k += 2;
      pvinfo->knmin[i] = (short) k;
    }
    k = getsp(bcdspin, pvinfo);
    if (nspinqn < k)
      nspinqn = k;
    pvinfo->nqn = (short) k;
    if (ncards == 1) { /* fill in defaults with wt = 0 */
      pvinfom = vinfo;
      for (i = 1; i <= nvibm; ++i) {
        ++pvinfom;
        memcpy(pvinfom, vinfo, sizeof(SVIB));
      }
      ivib = -glob.nvib;
    }
    setwt(pvinfo, ivib, iax, iwtpl, iwtmn, vsym);       /* set up weights */
  } while (ivsym < 0);          /*   end of reading */
  setsp();
  /* set up format */
  if (glob.nqnn == 3)
    *itd = 3;
  glob.nqn = glob.nqnn;
  glob.vibfmt = (nvibm != 0); glob.iqfmt0 = 0;
  if (glob.vibfmt) {
    glob.nqn += 1;
    glob.iqfmt0 = 1000;
  }
  nt = glob.nqn; k = nt + nspinqn; glob.nqn = k;
  if (k > (*nfmt)) k = nt + 2;
  glob.maxqn = k; glob.nqn0 = nt;
  glob.iqfmt0 += 100 * nt + k;
  /*     make sure +/- l values come in pairs */
  pvinfo = pvinfom = vinfo;
  k = 0;
  for (i = 0; i <= nvibm; ++i) {
    pvinfo->nqn += (short) nt;
    pvinfo->wt[4] &= 3;
    j = pvinfo->wt[4];
    if (j != 0) {
      if (EVEN(pvinfo->gsym)) {
        pvinfo->wt[4] = 0; /* no Itot */
      } else if (pvinfo->wt[0] == pvinfo->wt[j] && 
                 pvinfo->wt[1] == pvinfo->wt[j ^ 1]) {
        pvinfo->wt[4] = 0;
      }
    }
    j = pvinfo->lvqn;
    if (k != 0) {  /* previous l was non-zero and first of pair */
      isym = (int)((unsigned int)pvinfo->gsym >> 1);
      if (isym >= 3 && MOD(j + k, isym) != 0) break; 
      if (k > j) {
        if (j == 0) break;
        pvinfom->lvqn = (short) (-j);
        if (pvinfo->knmin[0] == 0) pvinfo->knmin[0] = 2;
        if (pvinfo->knmin[3] == 0) pvinfo->knmin[3] = 2;
      } else {
        pvinfo->lvqn = (short) (-k); 
        if (pvinfom->knmin[0] == 0) pvinfom->knmin[0] = 2;
        if (pvinfom->knmin[3] == 0) pvinfom->knmin[3] = 2;
      }
      for (j = 0; j <= 4; ++j) {
        if (pvinfo->wt[j] != pvinfom->wt[j]) break;
      }
      if (j <= 4) break;
      if (pvinfo->gsym != pvinfom->gsym) break;
      if (pvinfo->ewt[0] != pvinfom->ewt[0])break;
      if (pvinfo->ewt[1] != pvinfom->ewt[1])break;
      pvinfo->lvupper += 1;
      j = 0;
    }
    k = j;
    pvinfom = pvinfo; ++pvinfo;
  }
  if (k != 0) {
    puts("L DOUBLETS SHOULD BE IN PAIRS");
    exit(EXIT_FAILURE);
  }
  k = 0;
  pvinfo = vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib) {
    nt = pvinfo->nspstat;
    for (isym = 0; isym < 4; ++isym) {
      if (getwt(pvinfo, isym, 0, ivwt) <= 0) continue;
      for (ii = 1; ii <= nt; ++ii) {
        if (getwt(pvinfo, isym, ii, ivwt) > 0)
          ++k;
      }
    }
    ++pvinfo;
  }
  if (glob.nbkpj > 0) {
    free(moldv);
    moldv = NULL;
    free(blkptr);
    blkptr = NULL;
  }
  glob.nbkpj = k;
  nl = (size_t) k * sizeof(int);
  moldv = (int *) mallocq(nl);
  moldv[0] = 0;
  nl += sizeof(int);
  blkptr = (int *) mallocq(nl);
  blkptr[0] = 0;
  k = 0;
  pvinfo = vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib) {
    nt = pvinfo->nspstat;
    for (isym = 0; isym < 4; ++isym) {
      if (getwt(pvinfo, isym, 0, ivwt) <= 0) continue;
      for (ii = 1; ii <= nt; ++ii) {
        if (getwt(pvinfo, isym, ii, ivwt) > 0) {
          blkptr[k] = k;
          moldv[k] = isym + (ivib << 2) + (ii << glob.msshft);
          ++k;
        }
      }
    }
    ++pvinfo;
  }
  blkptr[k] = k;
  pvinfo = NULL;
  *nbcd = (NDECPAR + 1) + glob.vibdec;
  if (sizeof(int) == 2 && glob.vibdec > 2) {
    puts("too many vibrations for 16-bit integer computer");
    ncards = 0;
  }
  *nfmt = 1;
  if (nspinqn != 0 && nvibm != 0) 
    *nfmt = setfmt(&k, -1);
  return ncards;
}  /* setopt */

int setfmt(iqnfmt, nfmt)
int *iqnfmt, nfmt;
{
  SVIB *pvinfo;
  short *jjs;
  int iqfmtq, iqfmt0, i, j, k, iv, si2, ff, nset, nspinv, nvib;
  iqfmt0 = 0;
  nvib = glob.nvib;
  if (nfmt > 0 && nfmt < nvib)
    nvib = nfmt;
  pvinfo = vinfo;
  for (iv = 0; iv < nvib; ++iv) {
    iqfmtq = glob.iqfmt0;
    i = setgsym((int) pvinfo->gsym); 
    jjs = pvinfo->spt;
    nset = jjs[0];
    nspinv = nset - 1;
    if (nspinv > nspin)
      nspinv = nspin;
    jjs += nset;
    ff = 200 - jjs[0];
    if ((int) pvinfo->nqn > glob.maxqn) {
      iqfmtq += 10 * (ff & 1) + 4000;
    } else {
      si2 = 0; k = glob.maxqn - glob.nqn0 - 1;
      for (i = 0; i < k; ++i) {
        if (i == itsym && itsym < itptr) { 
          i = itptr; j = 0; 
          si2 = si2 << 1;
        }
        if (i < nspinv) {
          j = jjs[i + 1];
          if (i < itsym) 
            j += ff;
        } else {
          j = ff;
        }
        si2 = (si2 + (j & 1)) << 1;
      }
      si2 += (ff & 1);
      if (si2 > 9)
        si2 &= 7; 
      if (itptr < nspin)
        iqfmtq += 2000;
      if (itsym < itptr)
        iqfmtq += 4000; 
      iqfmtq += 10 * si2;
    }
    if (iv == 0) iqfmt0 = iqfmtq;
    if (nfmt < 0) {
      if (iqfmtq != iqfmt0) return nvib;
    } else {
      iqnfmt[iv] = iqfmtq;
    }
    ++pvinfo;
  }
  if (nfmt <= 0) return 1;
  return glob.maxqn;
}

int setblk(lu, npar, idpar, par, nbkpf, negy)
FILE *lu;
const int npar;
int *nbkpf, *negy;
bcd_t *idpar;
const double *par;
{
  /* this subroutine initializes  quantum quantities */
  /*     on entry: */
  /*         LU= logical unit for listing */
  /*         IDPAR= parameter id */
  /*         NBKPF= maximum value of F */
  /*         NEGY = max dimensioned energy in main program */
  /*     on return: */
  /*         NBKPF= blocks per F */
  /*         NEGY= possibly downsized NEGY to account for local storage */
  static char csym[] = { 'x', 'c', 'b', 'a' };
  static int nvwt[] = {1,1,1,2,3,2,3};
  SVIB *pvinfo;
  SPAR *spar_now;
  short *jjs;
  double ai;
  size_t nl;
  int nodd, isym, jsym, ixzt, i, k, kl, lt, n, kdiag, nvibm, neven, nsize, alpha;
  int kd, ld, ii, jj, im, kk, ni, nj, iv, jv, ivmin, ix, lv, lx, npair, ldel;
  int knnmin, ntstat, si1, si2, iff, ikq, jjt, ins, neuler, sznz, iff0;
  int iwt, nsym, gsym, maxblk, maxsblk, itmp, nd, ifc, nf;
  unsigned int ivsym;
  int ivwt[3], jvwt[3];
  char ctmp;

  cgetv[0].cblk = cgetv[1].cblk = 0;    /* reset store for getqn */
  kdiag = glob.idiag;
  if (pasort(lu, npar, idpar, par) == 0)
    glob.idiag = -1;
  if (npar <= 0) {
    if (glob.nvib > 1)
      free(vinfo);
    vinfo = &vinfo1;
    glob.nvib = 0;
    if (glob.nbkpj > 0) {
      free(moldv); moldv = &zmoldv;
      free(blkptr); blkptr = &zblkptr;
      glob.nbkpj = 0;
    }
    if (glob.maxblk > 0) {
      free(ivs); ivs = &zivs;
      free(ikmin); ikmin = &zikmin;
      free(ibkptr); ibkptr = &zibkptr;
      glob.maxblk = 0;
    }
    if (ndmx > 0) {
      free(idx); idx = &zidx;
      free(jdx); jdx = &zjdx;
      free(wk); wk = &zwk;
      ndmx = 0;
    }
    return 0;
  }
  nvibm = glob.nvib - 1;
  n = (*nbkpf);
  if (n > 0) {
    iff0 = n << 1;
    pvinfo = vinfo;
    for (i = 0; i <= nvibm; ++i) {
      jjs = pvinfo->spt;
      k = jjs[0];
      /*  select largest N relative to F */
      ni = (iff0 - jjs[k]) >> 1;
      if ((int) pvinfo->knmax > ni)
        pvinfo->knmax = (short) ni;
      ++pvinfo;
    }
    pvinfo = NULL;
  }
  pvinfo = NULL;
  iff0 = 200;
  /*    find interactions */
  ntstat = glob.nbkpj;
  sznz = ndmax = 0;
  if (ntstat <= 0)
    return 0;
  nf = nspin - 1;
  for (ii = 0; ii < ntstat; ++ii) {
    im = moldv[ii];
    iff = iff0;
    kk = getqs(im, iff, -1, -1, ixcom, iscom, &iv);
    if (ODD(iscom[nspin])) {
      /* special for odd spin multiplicity (make N integer) */
      kk = getqs(im, --iff, 0, -1, ixcom, iscom, &iv);
    }
    isym = ixcom[XSYM];
    ni = ixcom[XNVAL];
    pvinfo = &vinfo[ixcom[XVIB]];
    gsym = pvinfo->gsym;
    nsym = setgsym(gsym);
    getwt(pvinfo, isym, kk, ivwt);
    for (jj = 0; jj < ii; ++jj) {
      jjt = moldv[jj];
      kk = getqs(jjt, iff, 0, -1, jxcom, jscom, &jv);
      jsym = jxcom[XSYM];
      nj = jxcom[XNVAL];
      nd = ni - nj;
      pvinfo = &vinfo[jxcom[XVIB]];
      if (pvinfo->gsym != (short) gsym) continue;
      getwt(pvinfo, jsym, kk, jvwt);
      iwt = nvwt[nsym];
      for (k = 0; k < iwt; ++k) {
        if (ivwt[k] != jvwt[k]) break;
      }
      if (k < iwt) continue;
      /* check for connection restrictions on N */
      ixzt = glob.ixz;
      if (ODD(ixzt) && iscom[nspin] != jscom[nspin]) continue;
      for (k = 0; k < nf; ++k) {
        if (k == itsym)
          k = itptr;
        ixzt = ixzt >> 1;
        /* check for matching spin multiplicity*/
        if (ODD(iscom[k] + jscom[k])) break; 
        /* check for connection restrictions on spin */
        if (ODD(ixzt) && iscom[k] != jscom[k]) break; 
      }
      if (k < nf)
        continue;
      ivmin = (iv < jv) ? iv : jv;
      itmp = iv + jv + ivmin * glob.vibfac;
      ivsym = blksym(ixcom, jxcom) + ((unsigned int)itmp << 2);
      for (spar_now = spar_head[ivmin]; spar_now != NULL;
           spar_now = spar_now->next) {
        if (spar_now->ipsym < ivsym) break;
        if (spar_now->ipsym > ivsym) continue;
        if (((int) spar_now->flags & MNSQ) != 0) continue;
        if (sznz < 0)
          sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
        kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
                    &sznz, &ifc, &alpha, &ldel);
        if (ODD(neuler)) continue;
        if (sznz > 0) {
          sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
        }
        if (getmask(ixcom, jxcom, kd, ldel, kl, alpha) == 0) continue;
        npair = getll(0, ld, lt, kd, si1, si2, lscom, iscom, jscom);
        if (npair < 0) continue;
        if (kd > 0)
          glob.idiag = kdiag;
        if (nd > ndmax)
          ndmax = nd;
        kd = blkptr[ii] - blkptr[jj];
        if (kd != 0) { /* II and JJ are coupled */
          if (kd > 0) {         /* rename ptr II to ptr JJ */
            kk = blkptr[ii];
            ix = jj;
          } else {              /* rename ptr JJ to ptr II */
            kk = blkptr[jj];
            ix = ii;
          }
          glob.idiag = kdiag;
          for (k = 0; k <= ii; ++k) {
            if (blkptr[k] == kk)
              blkptr[k] = blkptr[ix];
          }
        }
        break;
      } /* loop over parameters */
    } /* jj loop */
  } /* ii loop */
  /*  block connections now found */
  switch (glob.idiag) {
  case -1:
    fputs("NO DIAGONALIZATION NEEDED\n", lu);
    break;
  case 0:
    fputs("ENERGY SORT OF WANG SUB-BLOCKS\n", lu);
    break;
  case 1:
    fputs("EIGENVECTOR SORT OF STATES\n", lu);
    break;
  case 2:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN WANG SUB-BLOCK\n", lu);
    break;
  case 3:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN V AND SPIN SUB-BLOCKS\n", lu);
    break;
  }
  if (glob.oblate) {
    fputs("OBLATE ROTOR\n", lu);
    ctmp = csym[1];
    csym[1] = csym[3];
    csym[3] = ctmp;
  } else {
    fputs("PROLATE ROTOR\n", lu);
  }
  if (glob.nqnn == 1) {
    fputs("LINEAR MOLECULE QUANTA, K SUPPRESSED\n", lu);
  } else if (glob.nqnn == 2) {
    fputs("SYMMETRIC TOP QUANTA\n", lu);
  }
  fputs("    V KMIN KMAX WTPL WTMN ESYMWT NSYM SPINS\n", lu);
  pvinfo = vinfo; jj = 0;
  for (iv = 0; iv <= nvibm; ++iv) {
    knnmin = pvinfo->knmin[0];
    ii = pvinfo->knmin[1];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[2];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[3];
    if (ii < knnmin)
      knnmin = ii;
    if (knnmin < 0) 
      knnmin = 0;
    k  = pvinfo->gsym; isym = k >> 1;  
    ii = (int)(pvinfo->lvupper >> 2);
    if ((pvinfo->lvupper & 2) != 0) 
      ii += 10;
    if (ODD(k)) {
      isym = -isym; jj = 1;
    }
    ii = ii * glob.esymdec + pvinfo->ewt[0]; 
    fprintf(lu, " %4d %4d %4d %4d %4d %6d %4d", iv, knnmin, pvinfo->knmax, 
            pvinfo->wt[0], pvinfo->wt[3], ii, isym);
    jjs = pvinfo->spt;
    ni = jjs[0];
    for (i = 1; i < ni; ++i) {
      ai = jjs[i] * 0.5;
      fprintf(lu, "%5.1f", ai);
    }
    fputc('\n', lu);
    if (EVEN(isym)) { /* even */
      i = ii; ii -= pvinfo->ewt[0]; 
      if (isym != 4) ii += pvinfo->ewt[1];
      if (i != ii || pvinfo->wt[0] != pvinfo->wt[2] || 
                     pvinfo->wt[3] != pvinfo->wt[1]) {
        fprintf(lu, " %4d %4d %4d %4d %4d %6d    B\n", iv, knnmin, pvinfo->knmax, 
                pvinfo->wt[2], pvinfo->wt[1], ii);
      }
    }
    ++pvinfo;
  }
  if (jj != 0) {
    fprintf(lu, "I(TOT) IS LAST QUANTUM BEFORE F AND %s\n", 
      "IS INDICATED BY NEGATIVE VALUE OF NSYM ABOVE");
  }
  /* bubble sort block number */
  lx = ntstat;
  while (lx > 1) {
    lv = 0;
    im = 0;
    for (i = 1; i < lx; ++i) {
      if (blkptr[im] > blkptr[i] ||
          (blkptr[im] == blkptr[i] && moldv[im] > moldv[i])) {
        lv = i;
        jj = blkptr[i];
        blkptr[i] = blkptr[im];
        blkptr[im] = jj;
        jj = moldv[i];
        moldv[i] = moldv[im];
        moldv[im] = jj;
      }
      im = i;
    }
    lx = lv;
  }
  fputs("BLOCK - WT - SYM - V - TSP - N  -  other quanta  (rel. to F=0 )\n",
        lu);
  /* convert block label to pointer */
  i = blkptr[0];
  blkptr[0] = 0;
  n = nsize = neven = nodd = maxblk = maxsblk = 0;
  for (jj = 0; jj <= ntstat; ++jj) {
    if (blkptr[jj] != i) {
      ii = jj - blkptr[n++];
      if (maxblk < ii)
        maxblk = ii;
      blkptr[n] = jj;
      if (neven > nsize)
        nsize = neven;
      if (nodd > nsize)
        nsize = nodd;
      if (jj == ntstat) break;
      i = blkptr[jj];
      neven = 0;
      nodd = 0;
    }
    jjt = moldv[jj];
    jjt = getqs(jjt, 0, -1, 0, ixcom, iscom, &ii);
    isym = ixcom[XSYM];
    iv = ixcom[XVIB];
    pvinfo = &vinfo[iv];
    iwt = getwt(pvinfo, isym, jjt, ivwt);
    kk = pvinfo->knmin[isym];
    if (kk < 0) {
      ii = ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (is_esym[ii] < 0) kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of even block */
    neven += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    kk = pvinfo->knmin[3 - isym];
    if (kk < 0) {
      ii = ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (is_esym[ii] >= 0) kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of odd block */
    nodd += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    --jjt;
    fprintf(lu, "% 5d %4d    %c %4d %4d", n + 1, iwt, csym[isym], iv, jjt);
    ni = nspin;
    for (ii = 0; ii < nspin; ++ii) {
      ai = iscom[ni] * 0.5;
      fprintf(lu, " %5.1f", ai);
      if (ni == itsym && itsym < itptr) 
        ii = itptr;
      ni = ii;
    }
    fputc('\n', lu);
  }
  ++maxblk;
  glob.maxblk = maxblk;
  fprintf(lu, "Maximum Dimension for Hamiltonian = %d\n", nsize);
  fflush(lu);
  *nbkpf = n;
  glob.nbkpj = n;
  if (nsize < (*negy))
    *negy = nsize;
  if (ndmx > 0) {
    free(ivs);
    ivs = NULL;
    free(ikmin);
    ikmin = NULL;
    free(ibkptr);
    ibkptr = NULL;
    free(idx);
    idx = NULL;
    free(jdx);
    jdx = NULL;
    free(wk);
    wk = NULL;
    ndmx = 0;
  }
  ii = maxsblk + maxsblk - 1;
  ndmx = nsize;
  if (ndmx < ii)
    ndmx = ii;
  nl = (size_t) (nsize + ndmx) * sizeof(double);
  wk = (double *) mallocq(nl);
  wk[0] = zero;
  /* allocate space for sub-block information */
  nl = (size_t) (ndmx + 1) * sizeof(short);
  idx = (short *) mallocq(nl);
  idx[0] = 0;
  jdx = (short *) mallocq(nl);
  jdx[0] = 0;
  maxblk = maxblk + maxblk;
  nl = maxblk * sizeof(short);
  ibkptr = (short *) mallocq(nl);
  ibkptr[0] = 0;
  ikmin = (short *) mallocq(nl);
  ikmin[0] = 0;
  nl = (size_t) maxblk * sizeof(int);
  ivs = (int *) mallocq(nl);
  ivs[0] = 0;
  ii = glob.vibfac + 1;
  return (ii * ii);
}                               /* setblk */

int pasort(lu, npar, idpar, par)
FILE *lu;
const int npar;
bcd_t *idpar;
const double *par;
{
  /*     initialize spin factors and arrange operators */
  /*     on entry: */
  /*         NSIZE= block size and dimension */
  /*         IDPAR= list of parameter identifiers ( element 0 is length ) */
  static int idmy[NDXCOM], initl;
  static bcd_t idunit[NDECPAR];
  SVIB *pvib1, *pvib2;
  SPAR *spar_now, *spar_last, *spar_match, *spar_free;
  short *iiv1, *iiv2;
  double dtmp;
  size_t nl;
  int ifac, iv1d, ipar, iv12q, ikqp, insp, njqt, ityi, isym, neulerp, i, lt;
  int kk, ltp, si1, si2, iv1, iv2, lv1, lv2, kdp, ldp, iv12, ikq, alpha, gsym;
  int ins, neuler, sznz, sznzp, si1p, si2p, nt, kd, ld, ivdif, ibcd, ldel;
  int iret, idflags, ifc, ifcp, ilim, kl, klp, nitot, ndecv, nbcd, imag;
  int alphap, ldelp;
  unsigned int ivsym;
  bcd_t *idval;
  BOOL isused, first;

  if (glob.parinit > 0) {
    for (i = 0; i < glob.parinit; ++i) {
      spar_free = spar_head[i];
      while (spar_free != NULL) {
        spar_now = spar_free->next;
        free(spar_free);
        spar_free = spar_now;
      }
      spar_head[i] = spar_free;
    }
    spar_free = spar_now = NULL;
    glob.parinit = 0;
  }
  if (npar <= 0) {
    free(ipder);
    ipder = &zipder;
    initl = 0;
    return 0;
  }
  /*  find reduced spin matrix elements */
  spfac[0] = 0.; spfac2[0] = 0.; 
  spfac[1] = sqrt(1.5); spfac2[0] = 0.;
  nt = glob.mxspin;
  for (i = 2; i <= nt; ++i) {
    dtmp = 0.5 * i; /* dtmp = I */
    spfac[i] = sqrt(dtmp * (dtmp + 1.) * (i + 1));
    spfac2[i] = 0.25 * sqrt((dtmp + 1.5) / (dtmp - 0.5)) / dtmp;
  }
  /* initialize SPECFC and DIRCOS */
  specfc(0, 0, 0, 0, 0, 0, 0, &zero, &szero, &szero);
  dircos(idmy, idmy, 0, 0, 0, &zero, &szero, &szero, 0, MODD, 0, &first);
  ifac = glob.vibfac + 1; ndecv = glob.vibdec;
  ilim = ifac * ifac - 1;
  if (initl > 0) {
    free(ipder);
    ipder = NULL;
  }
  nl = (size_t) npar *sizeof(int);
  ipder = (int *) mallocq(nl);
  initl = npar;
  kk = 1;
  ityi = 1;
  ipder[0] = 0; ibcd = 0; nbcd = (int)idpar[0] & 0x7f;
  for (i = 1; i < npar; ++i) {
    ibcd += nbcd;
    if (NEGBCD(idpar[ibcd]) == 0) {
      ipder[i] = kk++;
      ityi = i + 1;
    } else {
      ipder[i] = -ityi;
    }
  }
  glob.nfit = kk;
  /* .. code parameter id  for power of N*(N+1) and equivalence */
  /* .. NJQ = POWER + 1  (first occurrance of power) */
  /* .. IPTYP= ITY + 1    (first occurrance of cosine) */
  /* ..      ITY = sequence number for hybrid operators */
  /* .. negate NJQ and IPTYP for successive occurrance */
  /* .. IP points to parameter */
  glob.parinit = glob.nvib;
  for (i = 0; i < glob.nvib; ++i) {
    spar_head[i] = NULL;
  }
  pvib1 = pvib2 = vinfo;
  nsqmax = ityi = iv1d = iret = iv1 = iv2 = nitot = gsym = 0; 
  ipar = npar - 1; ibcd = ipar * nbcd;
  spar_free = NULL;
  isused = FALSE; 
  while (ipar >= -1) {
    ityi = 0; first = TRUE;
    do {                        /* repeat until parameter subtypes exhausted */
      ivdif = 0;
      if (ipar < 0) {           /* make sure there is a unit operator */
        idval = idunit;
        iv12q = ilim;
      } else {
        idval = &idpar[ibcd + 1];
        iv12q = bcd2i(idval[0]);
        if (ndecv > 1)
          iv12q += bcd2i(idval[1]) * 100;
        if (ndecv > 2) 
          iv12q += bcd2i(idval[2]) * 10000;
        idval += ndecv;
        iv2 = iv12q / ifac; iv1 = iv12q - iv2 * ifac;
        ivdif = iv1 - iv2;
        if (ivdif < 0) {
          i = iv1; iv1 = iv2; iv2 = i;
          iv12q = iv2 * ifac + iv1;
        }
      }
      if (iv12q == ilim) {
        iv1 = iv2 = iv1d; ivdif = 0;
        iv12 = iv1d * (ifac + 1);
      } else {
        if (iv1 >= glob.nvib) break;
        iv12 = iv12q;
      }
      if (idval[1] == (bcd_t)0x91 && idval[0] == (bcd_t)0) {
        /* special for rho */
        if (iv1 == iv2) {
          isused = TRUE;
          dtmp = par[ipar];
          specfc(0, iv1, iv1, 0, 0, 0, 1, &dtmp, &szero, &szero);
        }
        break;
      }
      pvib1 = &vinfo[iv1]; pvib2 = &vinfo[iv2];
      gsym = pvib1->gsym;
      if ((short) gsym != pvib2->gsym) break;
      gsym = setgsym(gsym); nitot = glob.nitot;
      if (spar_free == NULL) {  /* allocate more structures */
        spar_free = (SPAR *) mallocq(sizeof(SPAR));
      }
      spar_now = spar_free;
      spar_now->flags = 0; spar_now->alpha = C0; spar_now->mldel = C0;
      ityi = idpari(idval, ityi, spar_now);
      if (ityi == 0) break;
      njqt = (int) spar_now->njq;
      if (njqt > nsqmax)
        nsqmax = njqt;
      isym = (int) spar_now->ipsym; 
      idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
             &sznz, &ifc, &alpha, &ldel);
      iiv1 = pvib1->spt; iiv2 = pvib2->spt;
      if (si1 > 0 && 
          checksp(first, si1, si2, iiv1, iiv2, &spar_now->zfac) != 0) break;
      if (sznz > 0) {           /*  setup for SzNz */
        if (checksp(first, 1, 0, iiv1, iiv2, &spar_now->zfac) != 0) break;
        spar_now->zfac *= 0.5;
      }
      first = FALSE;
      /* factor of 2 for B-C, etc */
      if (kd != 0 && isym == 0)
        spar_now->zfac *= 2.;
      /* setup flags */
      idflags = 1;
      imag = (kd > ld)? kd: ld;
      if (ifc < 0) 
        ++imag;
      if (si2 < 0)
        ++imag;  /* commutator is odd */
      if (pvib1->wt[4] != pvib2->wt[4]) break;
      lv1 = pvib1->lvqn; lv2 = pvib2->lvqn;
      if (lv1 + lv2 < 0) break;
      if (ODD(pvib1->lvupper)) {
        --iv12;
        --iv1;
      }
      if (ODD(pvib2->lvupper)) {
        iv12 -= ifac;
        --iv2;
      }
      /* setup for l-doubling */
      ldel = lv1 - lv2;
      if (ldel != 0) { /*  add flag for KL symmetry when off diagonal in l */
        if (ldel < 0) 
          ldel = -ldel;
        if (ivdif < 0)   /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
        glob.lsym = FALSE;
      } else if (lv1 != 0) {
        /*  change symmetry for Lz operator */
        i = ((pvib1->lvupper ^ pvib2->lvupper) & 2) >> 1;
        if (ODD(imag + i)){ 
          isym = 3 - isym; ++imag;
          idflags |= MLZ;
        }
      }
      spar_now->mldel = (signed char)ldel; 
      if (ODD(imag)) {
        if (iv1 == iv2 && ldel == 0) break;
        idflags |= MODD;
      }
      if (nitot < 3) { /* calculate alpha */
        if (gsym >= 3 && MOD(kd - ldel, gsym) != 0) 
          continue; /* matrix element breaks symmetry */
        alpha = 0;
      } else {
        alpha = MOD(ldel - kd, nitot);
        if (si1 <= itsym) {
          if (alpha != 0) continue;
          idflags |= MIDEN; /* identity spin operator */
        }
        if (alpha < 0)
          alpha += nitot;
        if (alpha == 0 || (alpha + alpha) == nitot) {
          if (ityi > 160) continue;
          if (ityi > 80 && kd == 0 && lv1 == 0 && lv2 == 0)
            continue;
        } else {
          glob.esym = FALSE;
          spar_now->zfac *= 0.5;
        }
        if (ityi >= 80) {
          i = ityi / 80; /* i = 3,2,1,0 */
          if (ODD2(i)){
            if (ODD(isym + imag)) 
              spar_now->zfac = -spar_now->zfac;
            idflags |= MSYM2; isym = 3 - isym;
          }
          if (ODD(i)){
            alpha += nitot; /* special for quasi-diagonal */
          }
        }
      }
      if (ikq != 0 || neuler != 0 || sznz != 0 || ifc != 0)
        idflags |= MNOUNIT; /* no unit matrix */
      spar_now->alpha = (unsigned char)alpha; 
      if (isym == 0 && iv1 == iv2 && kd != 0)
        iret = 1;
      ivsym = (unsigned int) isym + ((unsigned int)iv12 << 2);
      spar_now->ipsym = ivsym;
      spar_match = spar_last = NULL;
      for (spar_now = spar_head[iv2]; spar_now != NULL;
           spar_now = spar_now->next) {
        if (spar_now->ipsym <= ivsym) break;
        spar_last = spar_now;
      }
      if (EVEN(neuler)) {    /* not Euler denominator parameter */
        kl = idflags & MMASK;
        for ( /* no init */ ; spar_now != NULL; spar_now = spar_now->next) {
          if (spar_now->ipsym < ivsym) break;
          spar_last = spar_now;
          klp = idpars(spar_now, &ikqp, &neulerp, &ltp, &ldp, &kdp, &insp,
                       &si1p, &si2p, &sznzp, &ifcp, &alphap, &ldelp);
          klp &= MMASK;
          /* check for same K dependence */
          if (ikq == ikqp && ld == ldp && kd == kdp && kl == klp && 
              sznz == sznzp && ifc == ifcp && alpha == alphap && 
              ldel == ldelp && neuler == 0 && neulerp == 0) {
            /*  check for same operator except power of N*(N+1) */
            if (lt == ltp && ins == insp && si1 == si1p && si2 == si2p) {
              idflags |= MNSQ;
            } else {
              if ((idflags & MNSQ) != 0) break;
            }
            spar_match = spar_now;
          } else {              /* operator with different K dependence */
            if ((idflags & MNSQ) != 0) break;
          }
        }                       /* end loop over previous parameters */
      }
      isused = TRUE;
      if (spar_match != NULL) {
        if (ipar < 0) break;      /* already have a unit operator */
        idflags |= MCOS_OK;
        spar_last = spar_match;
      }
      spar_now = spar_free;
      spar_free = NULL;
      if (spar_last != NULL) {
        spar_now->next = spar_last->next;
        spar_last->next = spar_now;
      } else {
        spar_now->next = spar_head[iv2];
        spar_head[iv2] = spar_now;
      }
      spar_last = NULL;
      spar_now->ip = ipar;
      spar_now->flags = (short) idflags;
    } while (ityi > 1);
    i = ipar;
    if (iv12q == ilim) {
      ++iv1d;
      if (iv1d >= glob.nvib)
        --ipar;
    } else {
      --ipar; 
    }
    if (i > ipar) {
      if (!isused) {
        putbcd(sbcd, NSBCD, &idpar[ibcd]);
        fprintf(lu,
                " WARNING: parameter %6d %s has no matrix elements\n",
                (i + 1), sbcd);
      }
      iv1d = 0; ibcd -= nbcd;
      isused = FALSE;
    }
  }                             /* end ipar loop */
  if (spar_free != NULL)
    free(spar_free);
  if (glob.esym)
    glob.lsym = FALSE;
  return iret;
}                               /* pasort */

int idpari(idval, itp, pspar)
bcd_t *idval;
int itp;
SPAR *pspar;
{
  /* subroutine to parse parameter identifier for complex parameter type */
  /*     on entry: */
  /*         IDVAL= parameter ID / vibrational field */
  /*         ITP= value of sub-parameter type for hybrid operator */
  /*     on return: */
  /*         pspar modified */

  /*0  1  2   3   4   5   6   7  8   9 */
  static int itpnxt[10] = { 0, 0, 0,  2,  6,  6,  2,  0, 9,  0 };
  static int itpop[10] =  { 0, 0, 0, 10, 40, 40, 10,  1, 0,  1 };
  static int zfacv[10] =  { 1, 2, 6,  3, -8,  8, -6,  2, 4, -4 };
  int ins, ity, nsx, isy, itysav;
  int itmp, i, iphaz, ldt, lnt, kdt, ksq, si1, si2, exval;
  bcd_t ibtmp;

  iphaz = 0;
  if (itp > 0) {
    /* get current sub-parameter type from previous */
    if (itp < 10) {
      itp = itpnxt[itp];
      if (itp == 0)
        return 0;
    } else {
      iphaz = itp >> 4;
      itp = itpnxt[itp & 15];
    }
  }
  pspar->zfac = 1.;
  pspar->euler = C0;
  pspar->msznz = C0;
  pspar->fc = C0;
  exval = 0;
  ibtmp = idval[0]; nsx = (int)(ibtmp & 0x0f); ksq = (int)(ibtmp >> 4);
  ity = bcd2i(idval[1]); itysav = ity * 10 + ksq;
  ibtmp = idval[2]; ins = (int)(ibtmp & 0x0f); si1 = (int)(ibtmp >> 4);
  if (ins > 0 && nspin == 0)
    return 0;
  if (ins >= 5 && itsym == 0)
    return 0;
  ibtmp = idval[3]; si2 = (int)(ibtmp & 0x0f); itmp = (int)(ibtmp >> 4);
  if (si2 > si1) {
    /*  make sure SI1 > SI2 */
    i = si1; si1 = si2; si2 = i;
    if (si2 == 0) 
      si2 = -1;
  }
  if (si1 > nspin) return 0;
  if (ity > 90) {
    if (si1 != 0 || ins >= 5)
      return 0;
  }
  ibtmp = idval[4];
  if (ibtmp != (bcd_t)0) { 
    if (ibtmp >= (bcd_t) 0x60) 
      return 0;
    itmp += 10 * (int)(ibtmp & 0x0f);
    exval = (int) ibtmp & 0xf0;
  }
  if (itmp != 0) {
    i = (itmp - 1) / 10;
    itmp -= (i >> 1) * 10;
    if (ODD(i))
      itmp = 9 - itmp;
    pspar->fc = (char) itmp;
  }
  pspar->msi1 = (char) si1;
  pspar->msi2 = (char) si2;
  pspar->njq = (unsigned char) nsx;
  /*        exclude S * S and I * S from aa, bb, cc operators */
  if (itp == 2) {
    if (si2 > 0) {
      itp = itpnxt[itp];
      if (itp == 0 && iphaz == 0)
        return 0;
    }
  }
  if (itp != 0) {
    if (itp == 6 && si1 == 0) {
      itp = 8;
      pspar->njq = (unsigned char) (nsx + 1);
    }
  } else {
    if (iphaz != 0) {
      if (ins >= 5) 
        iphaz -= 1;
      else if (glob.nitot >= 3)
        iphaz -= 5;
    } else {
      if (ins >= 5) 
        iphaz = 4;      /* iphaz = 4,3,2,1,0 */
      if (glob.nitot >= 3)
        iphaz += 15;
    }
    if (ity == 0) {
      itp = 1;
    } else if (ity <= 3) {
      itp = ity;
      if (glob.oblate)
        itp = revsym[itp];
      itp += 2;
      /*  ITP=3,4,5 for A,B,C */
      if (itp == 3 && si1 == 0)
        itp = 7;
    } else {
      itp = 1;
      itmp = ity;
      if (ity >= 92)
        itmp = (ity - 92) >> 1;
      /*  negate odd powers of P+**2 */
      if (itmp <= 11 && EVEN(itmp))
        pspar->zfac = -(pspar->zfac);
    }
  }
  isy = ity / 20;
  if (isy > 3) {
    isy = 0;
  } else if (isy > 0) {
    if (glob.oblate) {
      itysav += (revsym[isy] - isy) * 200;
    } else {
      isy = revsym[isy];
    }
  }
  if ((int) pspar->fc < 0)
    isy = 3 - isy;              /* fix for sin coeff. */
  pspar->ipsym = (unsigned int) isy;
  if (itp > 1) {
    itysav += itpop[itp] - ity * 10;
    pspar->zfac *= 2. / (double) zfacv[itp];
  }
  itp += iphaz << 4;
  if (ins >= 5) {
    ins -= 5;
    if(iphaz >= 5) iphaz = MOD(iphaz, 5);
    pspar->msznz = (unsigned char) (iphaz + 1);
  }
  pspar->mins = (unsigned char) ins;
  ity = itysav / 10;
  ksq = itysav - ity * 10;
  ldt = kdt = 0;
  if (ity <= 3) {
    if (ity != 0)
      ldt = 2;
  } else if (ity < 12) {
    ldt = kdt = ity + ity - 6;
  } else if (ity < 20) {
    kdt = ity + ity - 22;
    ldt = kdt + 1;
  } else if (ity < 80) {
    isy = ity / 20;
    ldt = kdt = ity - isy * 20 + 1;
    if (ODD((isy >> 1) + ldt))
      --kdt;
  } else if (ity < 90) {
    /* K energies coded as 8n  with K=10*n+KSQ */
    if (exval != 0) 
      return 0;
    ksq = -1 - 10 * (ity - 80) - ksq;
  } else {
    if (exval != 0) 
      return 0;
    i = ity - 90;
    if (ODD(i)) {
      if (ksq > 1) { /* augmented Euler series */
        --i; ksq += 8;
        kdt = i;
      } else {
        pspar->ipsym = 3; /* symmetry = 3 to put first in list */
      }
    } else {
      kdt = i;
    }
    pspar->euler = (unsigned char) (i + 2);
  }
  pspar->ksq = (char) ksq;
  if (si1 == 0) {
    lnt = 0;
  } else if (si2 <= 0) {
    if (isy == 0 || ldt == 0) 
      si2 = 0;
    else if (si2 < 0)
      --ldt;
    lnt = 1;
  } else if (ldt < 2) {
    lnt = 0;
  } else {
    lnt = 2;
  }
  pspar->mln = (unsigned char) lnt;
  pspar->mkdel = (unsigned char) kdt;
  /*  try to use raising operator */
  if (ldt == kdt && ldt > lnt)
    ldt -= (ldt - lnt) & 0x7e; /* subtract even number */
  pspar->mld = (unsigned char) ldt;
  if (exval != 0) {
    /* check for Euler restrictions */
    ibtmp = (bcd_t) exval; exval = exval >> 3;
    if (idval[1] == (bcd_t) 0 && idval[2] == (bcd_t) 0 &&
        idval[3] == (bcd_t) 0 && idval[4] == ibtmp) {
      /* energy Euler term */
      if (idval[0] == (bcd_t) 0) 
        return 0; /* bad energy-like term */
    } else if (exval == 2) {
      exval = 12; /* operator Euler term */
    }
    pspar->euler = (unsigned char) exval;
  }
  return itp;
}                               /* idpari */

int checksp(const BOOL first, int si1, int si2, const short *iiv1,
            const short *iiv2, double *zfac)
{
  int ii, iip;
  ii = iip = 0;
  if (si2 > 0) {
    if ((short) si2 < iiv1[0])
      ii = iiv1[si2];
    if ((short) si2 < iiv2[0])
      iip = iiv2[si2];
    if (ii == iip) {
      if (ii == 0) return 2;
      if (si2 == si1) {
        if (ii < 2) return 3;
        *zfac *= spfac2[ii];
      }
      *zfac *= spfac[ii];
    } else { /* ii != iip */
      if (ODD(ii + iip)) return 4; /* check multiplicity */
      if (si2 > itsym) return 4;
    }
  }
  if (si2 != si1){
    ii = iip = 0;
    if ((short) si1 < iiv1[0])
      ii = iiv1[si1];
    if ((short) si1 < iiv2[0])
      iip = iiv2[si1];
    if (ii == iip) {
      if (ii == 0) return 5;
      *zfac *= spfac[ii];
    } else { /* ii != iip */
      if (ODD(ii + iip)) return 6; /* check multiplicity */
      if (si1 > itsym) return 6;
    }
  }
  if (first && si1 > itsym){
    iip = 0;
    if (glob.nitot >= 3) iip = si1 - itsym;
    if (si1 == si2) {
      if (iip > 1) return 1;
      setzitot(2, 0, 2, ii, glob.nitot); /* quadrupole */
    } else if (si2 > itsym) {
      if (iip > 2) return 1;
      setzitot(1, 1, 0, ii, glob.nitot); /* 2-spin product */
      setzitot(1, 1, 2, ii, glob.nitot); 
    } else {
      if (iip > 1) return 1;
      setzitot(1, 0, 1, ii, glob.nitot); /* 1-spin vector */
    }
  }  
  return 0;
} /* checksp */

int setwt(pvinfov, ivib, iax, iwtpl, iwtmn, vsym)
SVIB *pvinfov;
const int ivib, iax, iwtpl, iwtmn;
double vsym;
{
  /*  set weights */
  /*     on entry: */
  /*         IVIB= number of vibrational states * number of wang states */
  /*         IAX= axis of symmetry */
  /*         IWTPL= weight for even rotational states */
  /*         IWTMN= weight for odd rotational states */
  /*         VSYM= symmetry code for vibrational states */
  /*     on return: */
  /*         IWT= weights for vibration rotation states */
  SVIB *pvinfo;
  int n, nvsym, i, k;
  short ii, jj, itmp;

  n = 1; nvsym = 0;
  if (ivib < 0) {
    if (vsym > 0.5)
      nvsym = 1;
    n = -ivib;
  }
  pvinfo = pvinfov;
  for (i = 0; i < n; ++i) {
    ii = (short) iwtpl;
    jj = (short) iwtmn;
    if (nvsym != 0) { 
      k = (int) (fmod(vsym, 10.) + 0.5);
      vsym = (vsym - k) / 10.;
      if (vsym < 0.) 
        nvsym = 0;
      if (ODD(k)) {
        itmp = ii; ii = jj; jj = itmp;
      }
    }
    switch (iax) {
    case 1: /*  WEIGHTS FOR Z */
      pvinfo->wt[0] = ii; 
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = ii;
      pvinfo->wt[4] = 2;
      ii = jj;
      break;
    case 2: /* WEIGHTS FOR B */
      pvinfo->wt[0] = ii; 
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 3: /* WEIGHTS FOR X */
      pvinfo->wt[0] = ii; 
      pvinfo->wt[1] = ii;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 2;
      break;
    case 4: /* A1, A2 only */
      pvinfo->wt[0] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 5: /* B1, B2 only */
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[4] = 3;      
      break;
    } /* end switch */ 
    ++pvinfo;
  }
  return 0;
} /* setwt */

int getwt(pvinfo, isym, iispin, ivwt)
SVIB *pvinfo;
const int isym, iispin;
int *ivwt;
{
  int iwt, jsym, ndm, nsym, kodd;
  short *jjs; 
  if (iispin == 0){
    iwt = 0;
    if (pvinfo->knmax >= pvinfo->knmin[isym] ||
        pvinfo->knmax >= pvinfo->knmin[3 - isym])
      iwt = 1;
    return iwt;
  }
  jsym = isym; iwt = pvinfo->wt[4]; nsym = pvinfo->gsym >> 1;
  if (iwt != 0) {    /* resolve Itot symmetries */
    jjs = pvinfo->spt; ndm = jjs[0];
    if (ODD2(jjs[iispin * ndm + nspin - nsym + 1]))
      jsym ^= iwt;
  }
  ivwt[0] = pvinfo->wt[jsym]; iwt = ivwt[0];
  ivwt[1] = 0; ivwt[2] = 0;
  if (nsym >= 3) {
    kodd = isoddk[jsym]; 
    if (nsym == 4){
      if (kodd != 0) {
        ivwt[1] = pvinfo->ewt[0]; /* E state for odd K */
        iwt = ivwt[1];
      } else {
        ivwt[2] = pvinfo->wt[jsym ^ 2]; /* check B weights */
        if (iwt == 0) iwt = ivwt[2];
      }
    } else {
      ivwt[1] = pvinfo->ewt[kodd];
      if (nsym == 6 && kodd != 0) {
        ivwt[0] = 0;
        ivwt[2] = pvinfo->wt[jsym ^ 2]; /* check B weights */
        iwt = ivwt[2];
      }
      if (iwt == 0) iwt = ivwt[1];
    }
  }
  return iwt;
} /* getwt */

int setgsym(const int gsym)
{
  static int oldgsym = -1;
  static int nsym;
  int k, kk;
  if (gsym == oldgsym) return nsym;
  oldgsym = gsym;
  nsym = gsym >> 1;
  is_esym[0] = 0;
  for (k = 1, kk = nsym - 1; k < kk; ++k, --kk) {
    is_esym[k] = 1; is_esym[kk] = -1;
  }
  if (k == kk) is_esym[k] = 0;
  if (ODD(gsym)) {
    itptr = nspin - 2;
    itsym = nspin - nsym;
    glob.nitot = nsym;
  } else {
    itptr = itsym = nspin + nspin + 1;
    glob.nitot = 0;
  }
  /* set up nominal spin coupling map */
  for (k = 0; k < nspin; ++k) {
    ismap[k + k] = k - 1; /* last N,J,F1 .. */
    if (k < itsym) {
      ismap[k + k + 1] = k + nspin + 1; /* S,I1,I2 .. */
    } else {
      ismap[k + k + 1] = itptr; /* Itot */
      break;
    } 
  }
  ismap[0] = nspin; /* fix up position of N */
  return nsym;
} /* setgsym */

int getsp(ispnx, pvinfo)
const bcd_t *ispnx;
SVIB *pvinfo;
{
  /*  subroutine to set up spin structure                        */
  /*     on entry:                                               */
  /*         ISPNX = spin code                                   */
  /*         IVIB  = vibrational quantum number                  */
  /*     on return:                                              */
  /*         JJS array of spin combinations                      */
  /*         (this array is initialized mostly in setsp)         */
  /*                                                             */
  SSP *ssp_now;
  short *jjs;
  size_t nl;
  int i, ii, jj, nspinv, idec, nspnx, nsstat, icmp, nset, isbig;
  int nitot, iret;
  short iis[MAXSPIN+1];
  bcd_t itmp;

  jj = -1; nspnx = (int)(ispnx[0] & 0x7f);
  isbig = 0; itmp = ispnx[1];
  if ((int)(itmp & 0x0f) == 0) {
    jj = (int)(itmp >> 4);
    isbig = 1;
  }
  /*  decode spins */
  nspinv = 0;
  nl = 1; idec = isbig;
  for (i = 1; i <= MAXSPIN; ++i) {
    if (isbig == 0 && jj >= 0) {
      ii = jj; jj = -1;
    } else {
      ++idec; 
      if (idec >= nspnx) break;
      itmp = ispnx[idec];
      ii = (int)(itmp & 0x0f);
      if (isbig != 0) 
        ii = ii * 10 + jj;
      jj = (int)(itmp >> 4);
    }
    if (ii == 0) break;
    if (ii >= MAXII) {
      ii = 1;
    } else if (ii > 1) {
      nspinv = i;
    }
    nl *= (size_t) ii;
    --ii; 
    iis[i] = (short) ii;
    if (ii > glob.mxspin)
      glob.mxspin = ii;
  }
  nsstat = (int) nl;
  i = (int)(nsstat << glob.msshft);
  if (i < 0 || (size_t)(i >> glob.msshft) != nl) {
    puts("spin problem too big");
    exit(EXIT_FAILURE);
  }
  nitot = 0; iret = nspinv; i = pvinfo->gsym;
  if (ODD(i)) { 
    nitot = i >> 1;
    if (nitot > nspinv) {
      iret += nitot;
      if (iret > MAXSPIN)
        iret = MAXSPIN;
      for (i = nspinv + 1; i <= iret; ++i) 
        iis[i] = 0;
      nspinv = iret;     
    }
    if (nitot > 3) {
      iret -= (nitot - 3);
    }
  }
  pvinfo->nspstat = nsstat;
  nset = nspinv + 1;
  iis[0] = (short) nset;
  ssp_now = &ssp_head;
  do {                          /*  check previous combinations for match */
    icmp = 0;
    jjs = ssp_now->sspt;
    for (i = 0; i < nset; ++i) {
      icmp = iis[i] - jjs[i];
      if (icmp != 0) break;
    }
    if (icmp == 0) 
      icmp = nitot - ssp_now->nitot;
    /*  match obtained, so return */
    if (icmp == 0) {
      pvinfo->spt = jjs;
      return iret;
    }
    ssp_now = ssp_now->next;
  } while (ssp_now != NULL);
  nl = (size_t) (nsstat + 1) * (size_t) nset * sizeof(short);
  jjs = (short *) mallocq(nl);
  nl = sizeof(SSP);
  ssp_now = (SSP *) mallocq(nl);
  ssp_now->next = ssp_head.next;
  ssp_head.next = ssp_now;
  ssp_now->sspt = jjs;
  ssp_now->ssize = nsstat;
  ssp_now->nitot = nitot;
  jjs[0] = (short) nset;
  for (i = 1; i < nset; ++i) {
    jjs[i] = iis[i];
  }
  pvinfo->spt = jjs;
  return iret;
}                               /* getsp */

void setsp(void)
{
  /*  subroutine to set up spin structure                        */
  /*   let: jjs = array of spin combinations in struct ssp.sspt  */
  /*         jjs[0]  = NSET = NSPIN + 1 = size of the set        */
  /*                   {N,J,F1..F} or {N,J..Itot,F}              */
  /*         jjs[1]..jjs[NSPIN] = S, I1, I2..                    */
  /*         jjs[k*NSET]..jjs[NSPIN-1+k*NSET] = 2*(N-F),2*(J-F)..*/
  /*         jjs[NSPIN+k*NSET] = minimum value of 2*F +(1)       */
  /*                                                             */
  SSP *ssp_now;
  ITMIX *pitmix0;
  EITMIX *pitmix;
  int nset, nspinv, ispsum, nn, iv, itot, ibgn, iend, ns, ns0, ifhalf;
  int minf0, minft, jj, ival, ibase, i, itmp, knt, nlim, nitot, ii;
  int itsym0, itptr0;
  short *iis, *jjs, *jjs0;
  nspin = 0;
  ssp_now = &ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL) { /* loop over spin sets */
    iis = ssp_now->sspt;
    nspinv = iis[0] - 1;
    for (iv = nspinv; iv > nspin; --iv) {
      if (iis[0] > 0) {
        nspin = iv;
        break;
      }
    }
  }
  ssp_now = &ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL) { /* loop over spin sets */
    iis = ssp_now->sspt;
    nset = iis[0];
    nspinv = nset - 1;
    if (nspinv > nspin) 
      nspinv = nspin;
    ispsum = 0;
    for (iv = 1; iv <= nspinv; ++iv) {
      ispsum += iis[iv];
    }
    ifhalf = (ispsum & 1);
    iend = iis[nspinv];
    ival = iend + 1; ns = ival;
    pitmix0 = NULL; pitmix = NULL;
    nitot = ssp_now->nitot;
    if (nitot == 0) { 
      nlim = nspinv; ibgn = iend;
      itsym0 = itptr0 = nspin + nspin + 1;
    } else { /* Itot coupling */
      itptr0 = nspin - 2;
      itsym0 = nspin - nitot;
      nlim = itsym0 + 1; 
      ii = iend;
      for (i = nlim; i < nspinv; ++i) {
        if ((int) iis[i] != ii) break;
        ns *= ival;
      }
      if (i < nspin) {
        puts("spins under Itot should be the same");
        exit(EXIT_FAILURE);
      }
      pitmix0 = get_itmix(ii, nitot);
      if (nitot < 3)
        pitmix0 = NULL;
      iend = nitot * ii; ibgn = (iend & 1);
    }
    knt = ssp_now->ssize; ns0 = knt / ns; 
    jjs = iis + nset;
    for (nn = -ispsum; nn <= ispsum; nn += 2) { /* select N quanta */
      minf0 = -ifhalf;
      if (nn < minf0)
        minf0 = nn;             /* nn = 2 * (N - F) */
      if (pitmix0 != NULL) 
        pitmix = pitmix0->mix;
      for (itot = ibgn; itot <= iend; itot += 2) {
        ns = ns0 * (itot + 1);
        if (pitmix != NULL && itot > ibgn)
          pitmix += 1;
        for (iv = 0; iv < ns; ++iv) {   /* try all possible spins */
          jj = nn;
          jjs[0] = (short) nn;
          ival = iv;
          minft = minf0;
          for (i = 1; i < nlim; ++i) {
            ibase = iis[i] + 1;
            itmp = ival;
            ival /= ibase;
            itmp -= ival * ibase;
            jj += itmp - iis[i];        /* jj = J + N - S - 2F */
            if (jj < minft)
              minft = jj;
            jj += itmp;         /* jj = 2 * (J - F) */
            jjs[i] = (short) jj;
          }
          jj += ival - itot;
          if (jj < minft)
            minft = jj;
          jj += ival;
          if (jj == 0) {        /* this possibility gives F=0 */
            if (itptr0 < nspinv)
              jjs[itptr0 + 1] = (short) itot;
            jj = ifhalf - minft;
            minft = jj + (jj & 1);
            jjs[nspinv] = (short) minft;
            if (pitmix != NULL){
              jjs[itsym0 + 1] = pitmix->qsym[0];
              for (i = itsym0 + 2; i <= itptr0; ++i)
                jjs[i] = 0;
              nitot = pitmix->n; jjs0 = jjs;
              for (i = 1; i < nitot; ++i) {
                jjs += nset; --knt;
                for (jj = 0; jj <= nspinv; ++jj) 
                  jjs[jj] = jjs0[jj];
                jjs[itsym0 + 1] = pitmix->qsym[i];
              }
            }
            if (--knt <= 0) {
              ns0 = 0;
              break;
            }
            jjs += nset;
          }
        }                       /* end loop over all spins */
      }                         /* end loop over Itot */
    }                           /* end loop over N */
  }                             /* end loop over spin types */
}                               /* setsp */


int getqq(iblk, f, iwtb, sbkptr, kmin, vs)
const int iblk;
int *f, *iwtb, *vs;
short *sbkptr, *kmin;
{
  /* subroutine to get sub-block structure and major quantum numbers */
  /*     on entry: */
  /*         MXBLK=biggest block size */
  /*         IBLK= block number */
  /*     on return: */
  /*         GETQQ= number of sub-blocks */
  /*         F= F quantum number *2 */
  /*         IWTB= statistical weight */
  /*         SBKPTR= pointer to beginning index for sub-block */
  /*         KMIN= minimum K for each sub-block */
  /*         VS= vibration and symmetry for each sub-block */
  SVIB *pvinfo;
  short *jjs;
  int ibgn, iend, ksym, i, n, nsize, isblk, nsblk, k, kk, nn, iv;
  int mss, nset, ff;
  unsigned int mvs;

  isblk = iblk - 1;
  if (glob.nbkpj <= 0)
    return 0;
  ff = isblk / glob.nbkpj;
  isblk -= glob.nbkpj * ff;     /* form remainder */
  ff = ff + ff;
  /*  get information for sub-blocks */
  nsize = nsblk = 0;
  ibgn = blkptr[isblk];
  iend = blkptr[isblk + 1];
  for (i = ibgn; i < iend; ++i) {
    mvs = (unsigned int) moldv[i];
    ksym = (int) mvs & 3;
    iv = (int) (mvs >> 2) & glob.msmask;
    mss = (int) (mvs >> glob.msshft);
    pvinfo = &vinfo[iv];
    if (i == ibgn) setgsym((int)pvinfo->gsym);
    jjs = pvinfo->spt;
    nset = jjs[0];
    jjs += nset * mss;
    if (ff < jjs[nset - 1]) continue;
    nn = jjs[0] + ff;
    n = nn >> 1;
    if (ODD(n))
      ksym = 3 - ksym;
    k = pvinfo->knmin[ksym];
    if (k < 0) {
      kk = MOD(jjs[itsym + 1] >> 2, glob.nitot); /* get spin symmetry */
      k = ksym;
      if (is_esym[kk] < 0) k += 2; 
      k &= 2; 
    }
    kk = pvinfo->knmax;
    if (kk > n)
      kk = n;
    kk -= k;
    if (kk >= 0) {
      sbkptr[nsblk] = (short) nsize;
      kmin[nsblk] = (short) k;
      vs[nsblk] = (int) mvs;
      nsize += 1 + (kk >> 1);
      if (nsblk == 0) {
        *f = ff - (nn & 1);
        k = pvinfo->wt[4];
        if (k != 0 && ODD2(jjs[itsym + 1]))
          mvs ^= k;
        mvs &= 3;
        iwtb[0] = pvinfo->wt[mvs];
        iwtb[1] = pvinfo->ewt[isoddk[mvs]];
        iwtb[2] = pvinfo->wt[mvs ^ 2];
        iwtb[3] = pvinfo->gsym;
        iwtb[4] = pvinfo->nqn; 
      }
      ++nsblk;
    }
  }
  if (nsblk > 0)
    sbkptr[nsblk] = (short) nsize;
  return nsblk;
}                               /* getqq */

int getqn(iblk, indx, maxqn, iqn, idgn)
const int iblk, indx, maxqn;
short *iqn;
int *idgn;
{
  /*   subroutine to get quantum numbers */
  /*     on entry: */
  /*         IBLK= block number */
  /*         INDX= index within block */
  /*     on return: */
  /*         IQN= list of quantum numbers */
  /*         IDGN= degeneracy */
  /*         NCOD= seraching aid to find state within a Wang block */
  static int last;
  int ibgn, iend, iqsp, isym, i, k, nsblk, ioff, lower, iv, iqf;
  int n, ncod, ldgn, ix, kq, ivbase, lv;
  short mgsym, ngsym;

  /*  get sub-block info from store or GETQQ */
  /*     look at last used then oldest */
  if (iblk == cgetv[0].cblk) {
    cgetq = cgetv;
    last = 0;
    ioff = 0;
    nsblk = cgetq->cnblk;
  } else if (iblk == cgetv[1].cblk) {
    cgetq = &cgetv[1];
    last = 1;
    ioff = glob.maxblk;
    nsblk = cgetq->cnblk;
  } else if (last != 0) {
    cgetq = cgetv;
    nsblk = getqq(iblk, &cgetq->cff, cgetq->cwt, ibkptr, ikmin, ivs);
    if (nsblk == 0) {
      *idgn = 0;
      return glob.nqn;
    }
    last = 0;
    ioff = 0;
    cgetq->cnblk = nsblk;
    cgetq->cblk = iblk;
    cgetq->csblk = 0;
  } else {
    cgetq = &cgetv[1];
    ioff = glob.maxblk;
    nsblk = getqq(iblk, &cgetq->cff, cgetq->cwt, &ibkptr[ioff],
                  &ikmin[ioff], &ivs[ioff]);
    if (nsblk == 0) {
      *idgn = 0;
      return glob.nqn;
    }
    last = 1;
    cgetq->cnblk = nsblk;
    cgetq->cblk = iblk;
    cgetq->csblk = ioff;
  }
  /*  check for request of size */
  if (indx <= 0) {
    *idgn = ibkptr[nsblk + ioff];
    return glob.nqn;
  }
  /*  search for sub-block */
  ix = indx - 1;
  ibgn = cgetq->csblk;
  iend = nsblk + ioff;
  if (ix < ibkptr[ibgn])
    ibgn = ioff;
  for (i = ibgn + 1; i < iend; ++i) {
    if (ix < ibkptr[i]) break;
  }
  cgetq->csblk = ibgn = i - 1;
  /*  assemble quanta */
  ncod = ibkptr[i] - ibkptr[ibgn];
  ldgn = cgetq->cwt[0];
  ngsym = (short) setgsym(cgetq->cwt[3]);
  iqf = cgetq->cff;
  iqsp = getqs(ivs[ibgn], iqf, 0, 0, ixcom, iscom, &ivbase) - 1;
  iv = ixcom[XVIB];
  isym = ixcom[XSYM];
  n = ixcom[XNVAL];
  k = ix - ibkptr[ibgn];
  kq = ikmin[ibgn] + (k << 1);
  if (ngsym >= 3) {
    mgsym = (short)MOD(kq - ixcom[XLVAL] + ixcom[XISYM] + ngsym, ngsym);
    if (mgsym != 0) {
      if ((mgsym + mgsym) == ngsym) { /* B symmetry */
        ldgn = cgetq->cwt[2];
      } else { /* E symmetry */
        if (glob.esym)
          isym = 4;
        ldgn = cgetq->cwt[1];
      }
    }
  }
  lv = ixcom[XLVAL];
  if (lv != 0) { /* l-doubled state */
    if (glob.lsym)
      isym = 4;
    ncod = 1;
    if (kq == 0 && isym == 3) {
      if (iv == ivbase) 
        ++iv;
      else 
        --iv;
    }
  }
  if (glob.nqnn == 2) {
    if (ODD(isym))
      kq = -kq;
    ncod = 1;
  } else {
    if (isym <= 3)
      lower = (n + isym) & 1;
    else /* degenerate state */
      lower = (lv > 0)? 1: 0; 
    if (glob.oblate) {
      iqn[2] = (short) kq;
      kq = n - kq + lower;
      ncod = -ncod;
    } else {
      iqn[2] = (short) (n - kq + lower);
    }
  }
  iqn[0] = (short) n;
  iqn[1] = (short) kq;
  k = glob.nqnn;
  if (glob.vibfmt)
    iqn[k++] = (short) iv;
  if (cgetq->cwt[4] > maxqn) {
    iqn[k++] = (short) iqsp;
    iqn[k] = (short) ((iqf + 1) >> 1);
  } else {
    for (i = 0; k < maxqn; ++i) {
      iqn[k++] = (short) ((iscom[i] + 1) >> 1);
      if (i == itsym && ngsym > 3)
        i += ngsym - 3; 
    }
  }
  ldgn *= iqf + 1;
  if (ldgn < 0)
    ldgn = 0;
  *idgn = ldgn;
  return ncod;
}     /* getqn */

int dircos(xbra, xket, ld, kd, ncmax, direl, ibra, iket, ifup, loff, mask,
           isunit)
const int *xbra, *xket;
const int ld, kd, ncmax, ifup, loff;
int mask;
double *direl;
short *ibra, *iket;
BOOL *isunit;
{ /*  SUBROUTINE TO  CALCULATE DIRECTION COSINE ELEMENTS */
  /* XBRA,XKET INTEGER VECTORS CONTAINING: */
  /*           0: DIMENSION OF SUB-BLOCK */
  /*           1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
  /*           2: 'N' QUANTUM NUMBER */
  /*           3: BEGINNING K */
  /*           4: L VALUE */
  /*           5: VIBRATIONAL QUANTA */
  /*           6: SPIN SYMMETRY */
  /* L IS TENSOR ORDER OF COSINE */
  /* KD IS TOTAL K QUANTIUM CHANGE */
  /* DIREL IS A VECTOR OF DIRECTION COSINE ELEMENTS */
  /* NCMAX IS MAXIMUM LENGTH OF DIREL */
  /* IBRA,IKET INDICES OF MATRIX ELEMENTS IN DIREL (START WITH ZERO) */
  /* IFUP =0 IF NO UPPER TRIANGLE TO BE CALC. */
  /* (LOFF AND MSYM2) IF OPERATOR HAS ALTERNATE SYMMETRY UNDER ITOT */
  /* (LOFF AND MLZ) IF OPERATOR INCLUDES LZ */
  /* (LOFF AND MODD) IF OPERATOR IS ODD-ORDER */
  /*  L=0  KD=0    DIREL=1. */
  /* IFUP <0 IF ELECTRIC DIPOLE */
  /* ODD ORDER PARAMETERS ARE IMAGINARY AND */
  /* EVEN ARE REAL IF NOT ELECTRIC DIPOLE */
  /* FOR ELECTRIC DIPOLE ODD ORDER IS REAL AND EVEN ARE IMAGINARY  */
  /*  L=1  KD=0    DIREL=PHI(Z) */
  /*  L=1  KD=1    DIREL=PHI(X) OR PHI(Y) */
  /*  L=2  KD=0    DIREL=0.5*(3* PHI(Z)**2 -1.) */
  /*  L=2  KD=1    DIREL=PHI(Z)*PHI(X)+PHI(X)*PHI(Z)  OR */
  /*                      =PHI(Z)*PHI(Y)+PHI(Y)*PHI(Z) */
  /*  L=2  KD=2    DIREL=2*(PHI(X)**2 - PHI(Y)**2) OR */
  /*                      =PHI(X)*PHI(Y)+PHI(Y)*PHI(X) */
  /*     IF KD > L THEN (P-)**(KD-L) OPERATOR MULTIPLIED BY */
  /*               L,L OPERATOR */
  /*     WANG FUNCTIONS MULTIPLIED BY i**(SYMMETRY) TO MAKE ALL */
  /*               ODD ANGULAR MOMENTA IMAG. AND ALL EVEN ONES REAL */
  /* IFC = FOURIER SERIES, > 0 IS COS, < 0 IS SIN */
  /* MASK = 3-BIT MASK TO ENABLE MATRIX ELEMENT SUB-TYPE */
  /*     BEST RESULTS ON SPEED IF N FOR BRA CHANGES SLOWEST */
  /***********************************************************************/
#define NFAC_DIRCOS 20
  static double ff[MAXN_DIRCOS], fac[NFAC_DIRCOS];
  double ele, elex, eletmp;
  int kdel, kbra, kbgn, nbgn, kket, isgn;
  int ncase, isum, kbra0, kdel2, kket0, kket2, i, k;
  int n, ikdel, kkdel, isbra, nnbra, nsbra, isket, idif;
  int ii, il, kk, lld, ir, ks, nnket, nsket, isrev;
  int iphaz, isym, nbra, nket, kbra0x, kket0x;
  BOOL matsym;
 
  kbra0 = xbra[XKBGN]; kbra0x = xbra[XIQN] & 1; 
  kket0 = xket[XKBGN]; kket0x = xket[XIQN] & 1; 
  nsbra = xbra[XDIM]; isbra = xbra[XSYM]; 
  nsket = xket[XDIM]; isket = xket[XSYM];
  isym = (isbra ^ isket) & 3;
  if ((loff & MODD) != 0) {
    if (ncmax == 0) { /* initialize */
      fac[0] = 0.5;
      for (lld = 1; lld < NFAC_DIRCOS; ++lld) {
        ele = lld - 0.5;
        fac[lld] = fac[lld - 1] * sqrt(lld / ele);
      }
      fac[0] = 1. / fac[1];
      return 0;
    };
    isym = 3 - isym; 
  }
  /****** set up state phases, isum is power of i */
  isrev = isym >> 1;
    /* isym is odd if operator is imaginary */
  isum = 4 + isket - isbra + (isym & 1);
  if (ifup < 0) /* correct to make dipoles real */
    ++isum;
  if (ODD(isum) && xbra[XVIB] < xket[XVIB])
    isum += 2; /* imaginary vibration phase assumed */
  isgn = isum >> 1;
  if ((loff & MLZ) != 0 && xket[XLVAL] < 0)
    ++isgn;    /* Lz operator found */
  ncase = 0;
  /*********** TRY SCALAR ******************************/
  if (kd == 0 && ld == 0 && (mask & 5) != 0) {
    mask &= 2;
    ii = 0; ncase = nsket; n = nsbra;
    kdel = kbra0 - kket0; kk = kbra0;
    if (kdel != 0) {
      if (kdel > 0) {
        kdel = kdel >> 1; ncase -= kdel;
        if (ncase < 0) ncase = 0;
      } else {
        ii = (-kdel) >> 1; kdel = -ii; n -= ii;
        if (n < 0) n = 0;
        kk = kket0;
      }
    }
    if (ncase > n)
      ncase = n;
    if (ncase > 0) {
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (kk == 0) { /* kk = max beginning k */
        if (kbra0x == kket0x) {
          if (kbra0x == 0) mask = 0;
          *isunit = TRUE; kk = -1; 
        } else { 
          *isunit = FALSE; mask = 0; /* kk = 0 */
        }
      } else { /* kk > 0 */
        *isunit = TRUE; mask = 0;
      }
      ele = 1.;
      if (ODD(isgn+isrev))
        ele = -ele;
      for (k = 0; k < ncase; ++k) {
        direl[k] = ele;
        ibra[k] = (short) ii;
        iket[k] = (short) (ii + kdel);
        ++ii;
      }
      if (kk == 0) direl[0] *= fac[0];
    }
    if (mask == 0) 
      return ncase; /* no quasi-diagonal delta K = 0 */
  } /************** end Scalar *********/
  *isunit = FALSE;
  nbra = xbra[XNVAL]; nket = xket[XNVAL];
  kdel = kd;
  ikdel = kdel - ld;
  if (ikdel <= 0) {
    kkdel = kdel;
  } else {
    k = kbra0 + ikdel + ((nsbra - 1) << 1);
    ffcal(nbra, k, ff);
    kkdel = ld;
  }
  kdel2 = kkdel + kkdel;
  lld = ld + ld;
  nnbra = nbra + nbra;
  nnket = nket + nket;
  if (ld > 0) {
    ele = (double) (nnbra + 1);
    if (nnket != nnbra)
      ele = sqrt(ele * (nnket + 1));
    if (kkdel != 0)
      ele *= fac[ld];
    isgn += nbra + kket0 + kkdel;
  } else {
    ele = 1.;
    if (kdel != 0)
      ele = 0.5;
  }
  matsym = TRUE;
  k = mask & 5;
  if (kdel == 0){
    if (k != 0) 
      mask = (mask & 2) + 1;
    matsym = FALSE;
  } else {
    if (k == 1 || k == 4)
      matsym = FALSE;
  }
  if (matsym) {
    for (k = 1; k <= 5; ++k) {
      if (xbra[k] != xket[k]) {
        matsym = FALSE; break;
      }
    }
  }
  if (ODD(isgn + isrev))
    ele = -ele;
  for (iphaz = 1; iphaz <= 4; iphaz = iphaz << 1) {
    /* iphaz = 1,2,4 */
    if (mask < iphaz) break;
    if (iphaz != 2) {
      /* iphaz = 1: +KDEL, iphaz = 4: -KDEL */
      /*  START CALCULATION OF  <K ....K-KDEL> */
      /*     KBRA=KBRA0+ (IL+ i)*KINC   IL >= 0 */
      /*     KKET=KKET0+ (IR+ i)*KINC   IR >= 0, i=0,... */
      /*     KBRA=KKET+KDEL */
      if ((mask & iphaz) == 0) continue;
      if (iphaz == 4) {
        if (matsym) break;
        kdel = -kdel; kkdel = -kkdel;
      } 
      isum = kdel + kket0 - kbra0;        /* ISUM = (IL-IR)*KINC */
      if (isum >= 0) {
        il = isum >> 1; ir = 0;
        kbra = kdel + kket0; kket = kket0;
      } else {
        il = 0; ir = (-isum) >> 1;
        kbra = kbra0; kket = kbra0 - kdel;
      }
      n = nsket - ir; i = nsbra - il;
      if (n > i) 
        n = i;
      if (n <= 0) continue;
      nbgn = ncase; ncase += n; idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      kket2 = kket + kket;
      if (ld == 0) {
        for (n = nbgn; n < ncase; ++n) {
          direl[n] = ele;
          ibra[n] = (short) il; iket[n] = (short) ir;
          ++il; ++ir;
        }
      } else {
        for (n = nbgn; n < ncase; ++n) {
          kk = -kdel2 - kket2;
          direl[n] = ele * c3jj(nnbra, lld, nnket, kk, kdel2, kket2);
          kket2 += 4;
          ibra[n] = (short) il; iket[n] = (short) ir;
          ++il; ++ir;
        }
      }
      il = ir = 1;
      if (kbra == 0) il = kbra0x;
      if (kket == 0) ir = kket0x;
      if (ODD(il + ir)) /*  CORRECT FOR K=0 */
        direl[nbgn] *= fac[0];
      if (ikdel > 0) {         /* ADD ON P- OR P+ OPERATOR TO LEFT */
        kbgn = kbra;
        if (kdel > 0) 
          kbgn -= ikdel;
        for (n = nbgn; n < ncase; ++n) {
          eletmp = direl[n]; kk = kbgn + ikdel;
          for (k = kbgn; k < kk; ++k) {
            eletmp *= ff[k];
          }
          direl[n] = eletmp;
          kbgn += 2;
        }
      }
    } else { /* iphaz == 2 */
      /*  DO QUASI-DIAGONAL PART */
      /*   <KDEL-K ..... K>   K=1,...,KDEL-1 */
      /*       KBRA=KBRA0+ (IL - i)*KINC  IL >= n-1 */
      /*       KKET=KKET0+ (IR + i)*KINC  IR >= 0  , i=0 .. n-1 */
      /*       KDEL=KKET+KBRA */
      kdel2 = -kdel2;
      if (ODD(isrev)) ele = -ele;
      if ((mask & 2) == 0) continue;
      isum = kdel - kbra0 - kket0;
      if (isum < 0) continue;
      isum = isum >> 1; /* isum = (IL + IR) */
      ir = 0; /* IR = 0 */
      if (isum >= nsbra) 
        ir = isum - nsbra + 1;
      else if (kket0x == 0)
        ++ir;
      n = isum + 1;
      if (n > nsket)
        n = nsket;
      else if (kbra0x == 0)
        --n;
      n -= ir;
      if (n <= 0) continue;
      il = isum - ir;
      if (matsym) 
        n = (il - ir + 2) >> 1;
      kket = kket0 + (ir << 1);
      kbra = kdel - kket;
      nbgn = ncase;
      ncase += n;
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (ODD2(isbra + nnbra)) {
        elex = -ele;
      } else {
        elex = ele;
      }
      kket2 = kket + kket;
      for (n = nbgn; n < ncase; ++n) {
        ibra[n] = (short) il; iket[n] = (short) ir;
        eletmp = elex;
        if (ld != 0) {
          ks = -kdel2 - kket2;
          eletmp *= c3jj(nnbra, lld, nnket, ks, kdel2, kket2);
          kket2 += 4;
        }
        for (kk = 0; kk < ikdel; ++kk) {
          k = kbra - kk;
          if (k <= 0) {
            eletmp *= ff[-k];
          } else {
            eletmp *= ff[k - 1];
          }
        }
        direl[n] = eletmp;
        --il; ++ir;
        kbra -= 2;
      }  /* loop over elements */
    } 
  } /* loop over 3 cases */
  if (matsym && ncase > 0 && ifup != 0) {
    n = ncase - 1;
    if (ibra[n] != iket[n])
      ++n;
    k = ncase;
    ncase += n;
    idif = ncmax - ncase;
    if (idif < 0)
      return idif;
    for (i = 0; i < n; ++i, ++k) {
      direl[k] = direl[i];
      ibra[k] = iket[i];
      iket[k] = ibra[i];
    }
  }  
  return ncase;
}  /* dircos */

int ffcal(nff, kff, ff)
const int nff, kff;
double *ff;
{
  static int nlast = -1;
  static int klast = 0;
  static double sq;
  int k;

  if (nff != nlast) {
    nlast = nff;
    sq = (double) nff;
    sq *= nff + 1;
    k = 0;
  } else {
    if (klast >= kff)
      return 0;
    k = klast;
  }
  klast = kff;
  if (klast > nff)
    klast = nff;
  while (k < klast) {
    ff[k] = sqrt(sq);
    ++k;
    sq -= k + k;
  }
  return 0;
}                               /* ffcal */
