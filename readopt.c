#include "catutil.h"
#include <ctype.h>
#ifndef NULL
#define NULL (void *)(0)
#endif

int readopt(char *line, double *dval, int nval)
{
  int kbgn, iret;
  char c;
  kbgn = -1;
  for(;;) {
    /* ignore leading whitespace */
    c = line[++kbgn];
    if (c >= '0' && c <= '9') break;
    if (c ==',' || kbgn > 1024 || c == '\0')
      return -1;
  } 
  while (kbgn < 1024) {
    /* find first trailing whitespace */
    c = line[++kbgn];
    if (c < '0' || c > '9') break;
  } 
  if (c != '\0' || kbgn < 1024) {
    iret = pcard(&line[kbgn], dval, nval, NULL);
    return kbgn;
  }
  return 0;
}
 
int readint(char *buf)
{
  int state, ii, k;
  state = 0;
  for (k = 0; buf[k] != '\0'; ++k) {
    ii = isdigit(buf[k]);
    if (state == 0 && ii != 0) {
      state = 1;
    } else if (state != 0 && ii == 0) {
      return k;
    }
  }
  return 0;
}
