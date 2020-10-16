#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include "slib.h"
#ifndef NULL
#define NULL (void *)(0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif


/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
/*   Revised version in c, 22 March 1999 */

/* potentially system dependent functions */

int chtime(str, N)
char str[];
const int N;
{ /* subroutine to return 24 characters with the time and date */
#define NCTIME 24
  char *buffer;
  time_t curtime;
  struct tm *loctime;
  int k, n;
  n = N;
  /* Get the current time. */
  curtime = time (NULL);
  /* Convert it to local time representation. */
  loctime = localtime (&curtime);
  /* copy the date and time in the standard format. */
  if (loctime != NULL) {
    buffer=asctime (loctime);
    k= (int) strlen(str);
    n -= NCTIME + 2;
    while(k < n) str[k++]=' ';
    str += n;
    memcpy(str,buffer, (size_t) 24);
    str[NCTIME] = '\n'; str[NCTIME + 1] = '\0';
  }
  return 0;
} /* chtime */

int caldelay(delay)
int delay;
{ 
  /* if delay is different, save time and return false */  
  /* else return true if time interval >= delay */
  static time_t tcmp = 0;
  static int old_delay = -1;

  if (delay == old_delay) {
    if (time(NULL) >= tcmp) {
      tcmp += delay; return 1;
    }
  } else {
    tcmp = time(NULL) + delay; old_delay = delay;
    if (delay <= 0) return 1;
  }
  return 0;
} /* caldelay */

void *mallocq(size_t nl)
{
  void *p;
  p = malloc((size_t) nl);
  if (p == NULL) {
    puts("memory allocation error");
    exit(EXIT_FAILURE);
  } 
  return p;
}
int maxmem(size_t *nll)
{
  /* On return 'nll' is the maximum double precision size that can be */ 
  /* allocated at once. This function returns the largest dimension of */
  /* a square double-precision matrix that can be allocated */
  size_t nl, nn, nroot;
  int k, n;
  n = (int) sizeof(size_t);
  nl = (size_t) 0x100 / sizeof(double) - 1;
  for (k = 1; k < n; ++k) nl = (nl << 8) + 0x0ff;
  if (nll) *nll = nl;
  nroot =  90; 
  for(k = 2; k < n; k += 2) nroot <<= 8; 
  do { /* iterative square root */
    nn = nroot; nroot = (nl /(nn + 1) + nn) >> 1;
  } while (nn != nroot);
  n = (int) sizeof(int);
  nl = 0x7fff;
  for(k = 2; k < n; ++k)  nl = (nl << 8) + 0x0ff;
  if (nroot > nl) nroot = nl;   
  return (int)nroot;
}
int filget(argc, argv,  nfile, cfil, cext)
const int   argc, nfile;
char *argv[];
char *cfil[];
const char *cext[];
{
#define NSTR 82
  char str[NSTR];
  char *master, *ext, *parg, *opt;
  /*@ dependent @*/ char *pstr;
  int  k, iext, len, iarg, knt, query;
  /*     Get NFILE names with extent given in CEXT */
  /*     Put names in CFIL */

  /*     gets file names from command line or interactively */
  /*     ARGC+1 is the number of command line arguments */
  /*     ARGV(I) is the I-th command line argument */

  /*     define character which delimits EXT field */

  for (iext=0; iext <= nfile; ++iext) cfil[iext]=NULL;
  knt = argc - 1; query = 0;
  if (knt <= 0) {
    query = 1; knt = nfile + 1;
  }
  master = NULL; opt = NULL; len = 0;
  for (iarg = 1; iarg <= knt; ++iarg) {
    if (query != 0) {
      parg = str;
      puts(" Enter file name ");
      if (fgetstr(parg, NSTR, stdin) <= 0) break;
    } else {
      parg = argv[iarg];
      if (parg == NULL) continue;
    }
    ext = strrchr(parg,'.');
    if (ext != NULL) {/* dot found */
      ++ext;
      for(iext = 0; iext < nfile; ++iext) {
        if (strcmp(cext[iext], ext) == 0) {
          k = (int) strlen(parg) + 1;
          pstr = (char *)malloc((size_t) k);
          if (pstr != NULL) memcpy(pstr, parg, (size_t)k);
          cfil[iext] = pstr;
          break;
        }
      }
    }
    if (parg[0] == '-') { /* set up options */
      if (opt != NULL) continue;
      k = (int) strlen(parg);
      opt = (char *)malloc((size_t) k);
      if (opt == NULL) break;
      --k;
      memcpy(opt, parg + 1, (size_t) k);
      opt[k] = '\0'; cfil[nfile] = opt;
      continue;
    }
    if (master == NULL) { /* set up default file name */
      k = (int) strlen(parg) + 1;
      if (ext != NULL) k -= (int)strlen(ext) + 1;
      if (k <= 0) continue;
      master = (char *)malloc((size_t) k);
      if (master == NULL) break; 
      len = k; --k;
      memcpy(master, parg, (size_t) k);
      master[k] = '.';
    }
  }
  if (master == NULL) return 0;
  for (iext = 0; iext < nfile; ++iext) {/* find default values */
    pstr = cfil[iext];
    if(pstr == NULL) {
      k = (int)strlen(cext[iext]) + 1;
      pstr = (char *) malloc((size_t)(len + k));
      if(pstr != NULL) { 
        memcpy(pstr, master, (size_t) len);
        memcpy(pstr + len, cext[iext], (size_t) k);
        cfil[iext] = pstr;
      }
    }
    if (pstr != NULL) puts(pstr);
  }
  free (master);
  return nfile;
} /* filget */

FILE *fopenq(const char *fname, const char *opt)
{
  FILE *hret;
  char *fnew, *path;
  if(opt[0] == 'l') {
    hret=fopen(fname, "r");
    if(hret == NULL) {
      path = getenv("SPECNAME");
      if (path != NULL){
        fnew = (char *)malloc(strlen(fname) + strlen(path) + 1);
        if (fnew != NULL) {
          strcpy(fnew,path); strcat(fnew,fname);
          hret = fopen(fnew,"r");
          free(fnew);
        }
      }else {
        puts("path variable SPECNAME can be set to path for names");
      }
    }
  } else {
    hret = fopen(fname,opt);
    if(hret == NULL) {
      printf(" Trouble opening %s\n", fname);
      exit(EXIT_FAILURE);
    }
  }
  return hret;
} /* fopenq */

int fgetstr(buffer, N, stream)
char *buffer;
const int N;
FILE *stream;
{ /* read first n bytes of line */
  char *pbuf;
  int k, n;
  n = N;
  if (stream != NULL && fgets(buffer,n,stream) == NULL) {
    buffer[0] = '\0';
    n = -1;
  } else {
    k = n = 0;
    for (pbuf = buffer; *pbuf != '\0' && *pbuf != '\n'; ++pbuf) {
      ++k;
      if (*pbuf != ' ' && !iscntrl(*pbuf)) n = k;
    }
    if(*pbuf == '\0') {
      do {
        k = fgetc(stream);
      } while (k != (int)'\n' && k != EOF); /* flush rest of line */
    }
    buffer[n] = '\0';
  }
  return n;
} /* fgetstr */

static int brkflg;

int rqexit(ival)
const int ival;
{
  if (ival != 0) {
    (void)signal(SIGINT, brkqr);
    brkflg = 0;
  }
  return brkflg;
} /* rqexit */

void brkqr(int i)
{
    if(i != 0) brkflg = 1;
} /* brkqr */













