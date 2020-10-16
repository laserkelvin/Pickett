/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
#include <stdio.h>
#ifndef BOOL
typedef int bool;
#define BOOL bool
#endif
#ifndef TRUE
#define TRUE (BOOL) 1
#endif
#ifndef FALSE
#define FALSE (BOOL) 0
#endif
#ifndef NULL
#define NULL (void *)(0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
/**********common structures*********************************************/
#define NDHEAPF 1200000L /* suggested length of store for heap (calfit)*/
#define NDHEAPC 1200000L /* suggested length of store for heap (calcat)*/
#define MAXCAT 6  /* number of quanta in catalog file */
#define MAXQN 10  /* maximum number of quanta */
#define NDECDIP 6 /* number of digit pairs + 1 for idip */

typedef struct {
  double cfrq;   /* calc. frequency   */
  double xfrq;   /* expt. frequency   */
  /*@dependent@*/ double *dnudp; /* freq. derivatives */
  float  xwt;    /* weight of line    */
  float  xerr;   /* error  of line    */
  int    linku;  /* link to next upper state */
  int    linkl;  /* link to next lower state */
  int    ibu;    /* upper state block number */
  int    ibl;    /* lower state block number */
  short  inu;    /* upper state index number */
  short  inl;    /* lower state index number */
  short  bln;    /* blend flag               */
  short  qn[2*MAXQN]; /* quantum number input     */
} SXLINE; 
typedef unsigned char bcd_t;
#define NEGBCD(ivbcd) (int)(ivbcd & 0x80)
/************** BLAS interface ******************************************/
#include "cblas.h"
#define idamax cblas_idamax
#define dasum  cblas_dasum
#define daxpy  cblas_daxpy
#define dcopy  cblas_dcopy
#define ddot   cblas_ddot
#define dnrm2  cblas_dnrm2
#define drot   cblas_drot
#define drotg  cblas_drotg
#define dscal  cblas_dscal
#define dswap  cblas_dswap

/************** CNJJ interface ******************************************/
#include "cnjj.h"
/************** SLIB interface ***************************************** */
#include "slib.h"
/************** ULIB interface ******************************************/
int ordblk(const int ndm, const int nn, short *iqnsep, double *t, double *e,
            short *isblk, double *p, short *ip);
void etswap(const int ndm,const int nsize,const int ix1,const int ix2,
            double *t,/*@null@*/double *e, /*@null@*/double *q);
int hdiag(const int nm, const int nx, double *z, /*@out@*/ double *d, 
          /*@out@*/ double *e, /*@out@*/ short *iqnsep);
int triag(const int nm, const int nx, const int nz,
          double *z, double *d, double *e);
int getpar(FILE *lu, FILE *luout, /*@out@*/int *npar, /*@out@*/int* npartot,
           /*@out@*/ bcd_t *idpar,/*@out@*/double *par,
           /*@out@*/ double *erpar, /*@out@*/char *plbl, int plblen);
int getvar(FILE *lu, const int npar, /*@out@*/ double *var, bcd_t *idpar,
           double *erpar, const int flg);
int putvar(FILE *lu, const int npar, double *var, double * erpar);
double calerr(const int npar,const double *var,const double *derv);
int deflin(int iqnfmt, /*@out@*/ short *idqn);
int getlin(FILE *lu, const int nqn, short *idqn, /*@out@*/ short *iqn,
           /*@out@*/ double *xfreq, /*@out@*/ double *xerr,
           /*@out@*/ double *xwt, /*@out@*/ char *card, const int ncard);
int getbcd(const char *line, bcd_t *ivbcd, int nbcd);
int putbcd(/*@out@*/ char *line, int nlen, const bcd_t *ivbcd);
int bcd2i(bcd_t btmp);
bcd_t i2bcd(int i);
#include "catutil.h"
/************** SPINV or DPI interface *********************************/
int hamx(const int iblk, const int nsiz, const int npar, const bcd_t *idpar,
         const double *par, /*@out@*/ double *egy, /*@out@*/ double *t, 
         /*@out@*/ double *dedp, /*@out@*/ double *pmix, const BOOL ifdump);
int setint(FILE *lu, /*@out@*/BOOL *ifdiag, /*@out@*/int *nsav, const int ndip,
           bcd_t *idip, /*@out@*/ int *isimag);
int intens(const int iblk, const int isiz, const int jblk, const int jsiz,
           const int ndip, const bcd_t *idip, const double *dip,
           /*@out@*/double *s);
int getqn(const int iblk,const int indx,const int maxqn, /*@out@*/short *iqn, 
          /*@out@*/int *idgn);
int setopt(FILE *lu, /*@out@*/ int *nfmt, /*@out@*/ int *itd,
           /*@out@*/ int *nbcdpar, /*@out@*/char *namfil);
int setfmt(/*@out@*/ int *iqnfmt, int nfmt);
int setblk(FILE *lu, const int npar, bcd_t *idpar, const double *par,
           int *nblkpf, int *negy);
/************** SUBFIT interfaces ***********************************/
int getlbl(int npar, bcd_t *idpar, /*@out@*/ char *parlbl, char *fil, int idiv,
           int len);
int getblk(/*@out@*/ int *iblk, /*@out@*/ int *indx, short *iqnum,
           int nblkpf, int ipos, int nqn);
int filbak(char *flu, char *fbak);
int prcorr(FILE *lufit, int npar, double *cor, int ndcor, double *err);
/*@dependent@*/ SXLINE *lbufof(int iflg, int ipos);
void dnuadd(int npar, int nparx, int initl, int indx,
            int ifac, double *egy, double *egyder, int nsize,
            int line, double *par, double *fac);
double dnuget(int iflg, int npar, double f, int line, double *dvec);
int getdbk(int *link, /*@out@*/ int *iblk, /*@out@*/ int *indx, 
           /*@out@*/ int *initl, /*@out@*/ int *ifac);
int frqdat(int line, /*@out@*/ int *ibln, /*@out@*/ double *xfrq, 
        /*@out@*/ double *xwt, /*@out@*/ double *xerr, /*@out@*/ short *iqn);
int lnlink(int *prvblk, int nblk, int iblk, int line);
/************** SORTSUB interfaces ***********************************/
int sortn(char *inpname, char *outname, BOOL dokey);





