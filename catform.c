#include <stdio.h>
#include <stdlib.h>

main()
{
#define NLIM 1000
#define BIG  99999999.9999
#define CLIGHT 29979.2458
  int cl, i, moltag, nline, iver, cnvt, ll, nmol, first, all;
  char *molname, *entval, *entnam, *catdir(), *fmakeword(), *makeword();
  char fqmin[14], fqmax[14], line[81], tbuf[22];
  double fmin, fmax, strlim, f, err, str, stof();
  float qrot[7];
  entval = getenv("REQUEST_METHOD");
  if (entval == 0 || strcmp(entval, "POST")) {
    printf("Content-type: text/html%c%c", 10, 10);
    beginform(stdout);
    moltag = 0;
    while (nxtdir(&moltag)) {
      molname = catdir(moltag, &nline, qrot, &iver);
      printf("<OPTION>%d %s %c", moltag, molname, 10);
    }
    endform(stdout);
    exit(0);
  } else {
    cnvt = nmol = fqmin[0] = fqmax[0] = 0;
    fmin = fmax = 0.;
    strlim = -500.;
    nline = 0;
    printf("Content-type: text/plain%c%c", 10, 10);
    cl = atoi(getenv("CONTENT_LENGTH"));
    while (cl && !feof(stdin)) {
      entval = fmakeword(stdin, '&', &cl);
      plustospace(entval);
      unescape_url(entval);
      entnam = makeword(entval, '=');
      /* printf("%s = %s\n",entnam,entval); */
      if (strcmp(entnam, "MinNu") == 0) {
        fmin = atof(entval);
      } else if (strcmp(entnam, "MaxNu") == 0) {
        fmax = atof(entval);
      } else if (strcmp(entnam, "UnitNu") == 0) {
        if (cnvt = strcmp(entval, "GHz")) {
          fmin = fmin * CLIGHT;
          fmax = fmax * CLIGHT;
        } else {
          fmin = fmin * 1000.;
          fmax = fmax * 1000.;
        }
        if (fmin > BIG)
          fmin = BIG;
        else if (fmin < 0.)
          fmin = 0.;
        if (fmax > BIG)
          fmax = BIG;
        else if (fmax < 0.)
          fmax = 0.;
        sprintf(fqmin, "%13.4f", fmin);
        sprintf(fqmax, "%13.4f", fmax);
      } else if (strcmp(entnam, "StrLim") == 0) {
        strlim = atof(entval);
      } else if (strcmp(entnam, "Mol") == 0) {
        all = (strcmp(entval, "All") == 0);
        if (all) {
          moltag = 0;
          nxtdir(&moltag);
        } else {
          moltag = atoi(entval);
        }
        while (moltag) {
          first = 1;
          if ((ll = catfrq(moltag, fqmin, line)) > 0) {
            while (strcmp(fqmax, line) >= 0) {
              if (stof(line + 21, 8) > strlim) {
                if (first) {
                  molname = catdir(moltag, &i, qrot, &iver);
                  printf("             %6d %s\n", moltag, molname);
                  nmol++;
                }
                nline++;
                first = 0;
                if (cnvt) {
                  f = stof(line, 13);
                  err = stof(line + 13, 8);
                  f = f / CLIGHT;
                  err = err / CLIGHT;
                  sprintf(tbuf, "%13.5f%8.5f", f, err);
                  for (i = 0; i < 21; i++)
                    line[i] = tbuf[i];
                }
                puts(line);
                if (nline >= NLIM) {
                  printf("%s %d lines.%s\n",
                         "THIS form is currently limited to",
                         nline, " Please limit your search.");
                  exit(1);
                }
              }
              if (catrd(moltag, ++ll, line))
                break;
            }
          }
          if (all) {
            if (nxtdir(&moltag) == 0)
              exit(1);
          } else {
            moltag = 0;
          }
        }
      }
    }
    if (nmol == 0)
      puts("No lines were found for your search criteria.");
  }
}

double stof(str, len)
char *str;
int len;
{
  int neg, dp, last, i;
  char c;
  double f;
  static double dpwr[6] = { 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001 };
  neg = 0;
  dp = len;
  last = -1;
  for (i = 0; i < len; i++) {
    c = str[i];
    if (c == '-') {
      neg = 1;
    } else if (c == '.') {
      dp = i;
    } else if (c >= '0' && c <= '9') {
      if (last >= 0) {
        f = f * 10. + (c - '0');
      } else {
        f = c - '0';
      }
      last = i;
    } else if (last >= 0) {
      break;
    }
  }
  if (last < 0)
    return 0.;
  if (neg)
    f = -f;
  dp = last - dp;
  while (dp) {
    i = dp;
    if (i > 6)
      i = 6;
    f = f * dpwr[i - 1];
    dp = dp - i;
  }
  return f;
}
