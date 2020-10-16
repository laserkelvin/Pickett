#ifndef FILE
#include <stdio.h>
#endif
/************** SLIB interface ***************************************** */
int  chtime(char str[], const int n);
int  caldelay(int delay);
/*@out@*/ /*@only@*/ /*@notnull@*/ void *mallocq(size_t nl);
int  maxmem(/*@null@*/ /*@out@*/ size_t *nl);
int  filget(const int argc, char *argv[], const int nfile,
            /*@out@*/ char *cfil[], const char *cext[]);
/*@dependent@*/ FILE *fopenq(const char *fname, const char *opt)
                             /*@modifies fileSystem@*/ ;
int  fgetstr(/*@out@*/char buffer[], const int n, FILE *stream);
int  rqexit(const int ival);
void brkqr(int);
