#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calpgm.h"
#define NDILINE 130

int sortn(char *inpname, char *outname, int dokey);
int fillbuf(FILE * finp, char *tbuf, int nlim, char **pbuf);
int getkey(int *key, char *tbuf, char *wbuf);
int bufcmp(const char *b1, const char *b2, const int *key);
int memsort(int nline, const int *key, char **pbuf);
int memfind(int ibgn, int nbuf, const int *key, const char *buf0, 
            char **pbuf);

static int ndline;

int sortn(inpname, outname, dokey)
char *inpname, *outname;
BOOL dokey;
{
  FILE *finp, *fout, *ftmp, *fscr;
  int *key;
  char *buf0, *buf, **pbuf;
  int nline, k, kk, mxline, nn, ntot;
  BOOL init, scrflg;
  finp = fopen(inpname, "r");
  if (finp == NULL) {
    puts("no input file");
    return 1;
  }
  nn = 0x7fff;
  if (sizeof(size_t) > 2)
    nn = 16385 * 81;
  buf0 = (char *) malloc((size_t) nn);
  if (buf0 == NULL) {
    return 1;
  }
  key = (int *) malloc((size_t) NDILINE * sizeof(int));
  if (key == NULL) {
    free(buf0); return 1;
  }
  ndline = NDILINE;
  pbuf = &buf; pbuf[0] = buf0;
  nline = fillbuf(finp, buf0, 1, pbuf);
  if (nline <= 0) {
    free(buf0); free(key); return 1;
  }
  ntot = strlen(buf0); ndline = ntot + 1;
  if (ndline < 81) 
    ndline = 81;
  if (!dokey) {
    for (k = 0; k < 13; ++k)
      key[k] = k + k;
    for (k = 0; k < 40; ++k) {
      kk = k + 55;
      if (kk >= ntot) break;
      key[k + 13] = kk + kk;
    }
    key[k + 13] = -1;
  } else {
    key[0] = -1;
    getkey(key, buf0, &buf0[ndline]);
  }
  ntot = 1;
  mxline = nn / ndline;
  pbuf = (char **) malloc((size_t) mxline * sizeof(char *));
  if (pbuf == NULL) {
    free(buf0); free(key); return 1;
  }
  buf = buf0 + ndline; --mxline; 
  fout = fscr = NULL;
  init = TRUE;
  do {
    pbuf[0] = buf;
    nline = fillbuf(finp, buf, mxline, pbuf);
    ntot += nline;
    memsort(nline, key, pbuf);
    if (nline < mxline) {
      fclose(finp); finp = NULL;
      if (fscr != NULL) {
        fclose(fscr); fscr = fout;
      }
      fout = fopen(outname, "w");
    } else if (fout != NULL) {  /* MXLINE lines read without EOF after first pass */
      rewind(fscr);
      ftmp = fscr;
      fscr = fout;
      fout = ftmp;
    } else if (init) {          /* MXLINE lines read without EOF on first pass */
      fout = tmpfile();
      fscr = tmpfile();
    }
    scrflg = (init || fillbuf(fscr, buf0, 1, pbuf) != 0);
    k = 0;
    while (scrflg || k < nline) {       /* merge buffer with scratch file */
      if (scrflg && k < nline) {
        nn = memfind(k, nline, key, buf0, pbuf);
      } else {
        nn = nline;
      }
      while (k < nn) {
        fputs(pbuf[k], fout); fputc('\n', fout);
        ++k;
      }
      if (scrflg) {
        fputs(buf0, fout); fputc('\n', fout);
        scrflg = (!init && fillbuf(fscr, buf0, 1, pbuf) != 0);
      }
    }
    printf("sorted %8d lines\n", ntot);
    fflush(stdout);
    rewind(fout);
    init = FALSE;
  } while (finp != NULL);
  if (fout != NULL)
    fclose(fout);
  if (fscr != NULL)
    fclose(fscr);
  buf = NULL;
  free(buf0); free(key); free(pbuf);
  return 0;
}

int fillbuf(finp, tbuf, nlim, pbuf)
FILE *finp;
char *tbuf, **pbuf;
int nlim;
{
  int i, iret, ich;
  iret = 0;
  for (;;) {
    if (fgets(tbuf, ndline, finp) == NULL)
      break;
    i = strlen(tbuf) - 1;
    if (tbuf[i] == '\n') {
      tbuf[i] = '\0';
    } else {
      do {
        ich = getc(finp);
      } while (ich != '\n' && ich != EOF);      /* flush rest of line */
    }
    i += 2;
    if (i < ndline) 
      memset(&tbuf[i], 0, ndline - i);
    if (++iret == nlim)
      break;
    tbuf += ndline;
    pbuf[iret] = tbuf;
  }
  return iret;
}

int bufcmp(b1, b2, key)
const char *b1, *b2;
const int *key;
{
  int iret, i, k, ib1, ib2, irev;
  irev = 0;
  for (i = 0; (k = key[i]) >= 0; i++) {
    if ((k & 1) != 0)
      irev = 0;
    k = k >> 1;
    ib1 = b1[k] & 0x7f;
    if ((char) ib1 == '-') { 
      ib1 = -1; irev = 0;
    }
    ib2 = b2[k] & 0x7f;
    if ((char) ib2 == '-') {
      ib2 = -1; irev = 0;
    }
    iret = ib1 - ib2;
    if (iret != 0) {
      if (irev != 0) iret = -iret;
      return iret;
    }
    if (ib1 < 0) irev = 1;
  }
  iret = (int)(b1 - b2);
  return iret;
}

int memsort(nline, key, pbuf)
int nline;
const int *key;
char **pbuf;
{                               /* heapsort */
  int left, right, i, j;
  char *pcmp;
  left = nline >> 1;
  right = nline - 1;
  while (right > 0) {
    if (left > 0) {
      pcmp = pbuf[--left];
      i = left;
      j = i + i + 1;
    } else {
      pcmp = pbuf[right];
      pbuf[right--] = pbuf[0];
      i = 0;
      j = 1;
    }
    while (j <= right) {
      if (j < right && bufcmp(pbuf[j], pbuf[j + 1], key) < 0)
        j++;
      if (bufcmp(pcmp, pbuf[j], key) > 0)
        break;
      pbuf[i] = pbuf[j];
      i = j;
      j = i + i + 1;
    }
    pbuf[i] = pcmp;
  }
  return 0;
}

int memfind(ibgn, nbuf, key, buf0, pbuf) 
int ibgn, nbuf;
const int *key;
const char *buf0;
char **pbuf;
{
  /* return position n such that buf0 > pbuf[k] for k = ibgn, .., (n - 1) */
  int left, right, k;
  if (bufcmp(buf0, pbuf[ibgn], key) <= 0) 
    return ibgn;
  left = ibgn + 1; right = nbuf;
  while (left < right) {
    k = (left + right) >> 1;
    if (bufcmp(buf0, pbuf[k], key) > 0) {
      left = k + 1;
    } else {
      right = k;
    }
  }
  return left;
}

int getkey(key, tbuf, wbuf)
int key[];
char tbuf[], wbuf[];
{
  int k, i, init, klast;
  char cmp;
  printf("type numbers for key fields (1= most significant)\n%s\n", tbuf);
  fflush(stdout);
  fgets(wbuf, ndline, stdin);
  i = 0;
  for (cmp = '1'; cmp <= '9'; ++cmp) {
    init = 1; klast = -2;
    for (k = 0; wbuf[k] != '\0'; ++k) {
      if (wbuf[k] == cmp) {
        if ((k - klast) > 1) init = 1;
        key[i++] = k + k + init;
        klast = k;
      }
      init = 0;
    }
  }
  key[i] = -1;
  return i;
}
