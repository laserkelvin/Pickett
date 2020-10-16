/* slibg77.f -- translated by f2c (version 19970219).*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include "calpgm.h"

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */

static integer brkflg;

int chtime(str, n)
char *str;
int n;
{
/*     subroutine to return 24 characters with the time and date */
#include <time.h>
  char *buffer;
  time_t curtime;
  struct tm *loctime;
  int k;
  /* Get the current time. */
  time(&curtime);
  /* Convert it to local time representation. */
  loctime = localtime(&curtime);
  /* copy the date and time in the standard format. */
  buffer = asctime(loctime);
  k = strlen(str);
  while (k < n)
    str[k++] = ' ';
  str += n;
  strncpy(str, buffer, 24);
  str[24] = 0;
  return 0;
}                               /* chtime */

int filget(argc, argv, nfile, cfil, cext)
int argc, nfile;
char *argv[], *cfil[], *cext[];
{
  /* Local variables */
  char str[82], master[82], *pstr, *ext, *pext;
  int k, iext, len, iarg, query;
  int fgetstr();
/*     Get NFILE names with extent given in CEXT */
/*     Put names in CFIL */

/*     gets file names from command line or interactively */
/*     ARGC+1 is the number of command line arguments */
/*     ARGV(I) is the I-th command line argument */

/*     define character which delimits EXT field */

  for (iext = 0; iext < nfile; iext++)
    cfil[iext] = NULL;
  query = (argc <= 1);
  if (query)
    argc = nfile;
  len = 0;
  for (iarg = 1; iarg < argc; iarg++) {
    if (query) {
      puts(" Enter file name ");
      if (fgetstr(str, 82, stdin) <= 0)
        break;
      pstr = str;
    } else {
      pstr = argv[iarg];
    }
    ext = strrchr(pstr, '.');
    if (ext) {                  /* dot found */
      len = ext - pstr;
      strlwr(++ext);
      for (iext = 0; iext < nfile; iext++) {
        pext = cfil[iext];
        for (k = 0; pext[k]; k++) {
          if (pext[k] != ext[k])
            break;
        }
        if (pext[k]) {
          if (query) {
            pstr = malloc(strlen(str) + 1);
            if (pstr)
              strcpy(pstr, str);
          }
          cfil[iext] = pstr;
          break;
        }
      }
    }
    if (iarg == 1) {            /* set up default file name */
      if (len == 0)
        len = strlen(pstr);
      if (len > 80)
        len = 80;
      strncpy(master, pstr, len);
      master[len++] = '.';
    }
  }
/* find default values */
  for (iext = 0; iext < nfile; iext++) {
    if (cfil[iext] == NULL) {
      if (len == 1)
        return 0;
      k = strlen(cext[iext]) + 1;
      pstr = malloc(len + k);
      if (pstr) {
        strncpy(pstr, master, len);
        strncpy(pstr + len, cext[iext], k);
      }
      cfil[iext] = pstr;
    }
    puts(cfil[iext]);
  }
  return nfile;
}                               /* filget */

FILE *fopenq(fname, opt)
char *fname, *opt;
{
  FILE *hret;
  char *fnew, *path, *getenv();
  if (opt[0] = 'l') {
    hret = fopen(fname, "r");
    if (hret == NULL) {
      if (path = getenv("SPECNAME")) {
        fnew = malloc(strlen(fname) + strlen(path) + 1);
        strcpy(fnew, path);
        strcat(fnew, fname);
        hret = fopen(fnew, "r");
      } else {
        puts("path variable SPECNAME can be set to path for names");
      }
    }
  } else {
    hret = fopen(fname, opt);
    if (hret == NULL) {
      printf(" Trouble opening %s\n", fname);
      exit(1);
    }
  }
  return hret;
}                               /* fopenq */

int rqexit(ival)
int ival;
{
  if (ival) {
    signal(SIGINT, brkqr);
    brkflg = 0;
  }
  return brkflg;
}                               /* rqexit_ */

/* Subroutine */ void brkqr(int i)
{
  brkflg = -1;
}                               /* brkqr */
