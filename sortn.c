#include <stdio.h>
#include <stdlib.h>
int sortn(char *impname, char *outname, int dokey);

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *finp;
  char *outnam;
  if (argc < 2) {
    puts("sortn infile [outfile]");
    exit(1);
  }
  finp = fopen(argv[1], "r");
  if (finp == NULL) {
    puts("no input file");
    exit(1);
  }
  if (argc >= 3) {
    outnam = argv[2];
  } else {
    outnam = argv[1];
  }
  return sortn(argv[1], outnam, 1);
}
