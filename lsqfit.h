int lsqfit(double *dk, int ndm, int nr, int nvec, 
           double *marqp, double *varv, /*@out@*/ double *dkold, 
           /*@out@*/ double *ediag,/*@out@*/ double *enorm, 
           /*@out@*/ int *iperm);
int jelim(double *dk, double *vec, int ndm, int nc, int nvec);
int dqrfac(double *dk, int ndm, int nc, int nr, int nvec, 
           /*@out@*/ double *wk, /*@out@*/ int *iperm);
int dqrsolv(double *dk, int ndm, int nrank, int nr, int nvec, 
            int *iperm);
int lmpar(int nr, int nc, int nvec, double *r, int ldr, /*@out@*/ double *x, 
          double *par, /*@out@*/ double *wk);

/* nc is the number of columns (number of experiments) */
/* nr is the number of rows (number of fitted parameters) */
/* nvec is the number of solution vectors */
/* nrank is the rank of the derivative matrix */
/* on input to lsqfit and dqrfac, */ 
/*            dk contains the derivatives in the first nr rows */
/*            and nvec rows of observations packed in succeeding rows */ 
/* on exit from lsqfit and dqrsolv, */ 
/*            dk contains the lower triangular pseudo inverse */
/*                     in the first nr rows  and nvec rows of solutions */ 
/* on entrance to jelim, dqrsolv, and lmpar, and exit from jelim and lmpar */
/*            dk contains a lower triangular matrix */
/*                     in the first nr rows */ 
/*                     and nvec rows of transformed observations */ 
