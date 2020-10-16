#include <string.h>
#include <stdlib.h>
#include "calpgm.h"

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 Mar 1989 */
/*   Revised version in c, 22 March 1999 */

/**************************************************************************/
char buf[82];

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *finp, *fout;
  static double clight = 29979.2458;
  static int fmt[]={13,8};
  double frq, err, dval[2];
  if (argc < 3) {
    puts("usage: cnvir input_file output_file");
    exit(1);
  }
  finp = fopen(argv[1],"r"); fout = fopen(argv[2],"w");
  if (finp == NULL || fout == NULL) {
    puts("bad file");
    exit(1);
  }
  while (fgetstr(buf,82,finp) > 0) {
    pcard(buf, dval, 2, fmt);
    frq = dval[0] / clight; err = dval[1] / clight;
    fprintf(fout,"%13.7f%8.7f%s\n", frq, err, buf + 21);
  }
}