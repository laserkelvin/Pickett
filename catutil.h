int readqn(const char *qnstr, short *iqn, const int n);
void gupfmt(int igup, /*@out@*/ char *sgup);
int pcard(const char *card, /*@out@*/ double *val, const int nval,
          /*@null@*/ const int *fmtlen);
