#include <time.h>
int strdate_(char *pbuffer, int plen)
{
  char *buffer;
  time_t curtime;
  struct tm *loctime;
  int i;

  /* Get the current time. */
  time(&curtime);

  /* Convert it to local time representation. */
  loctime = localtime(&curtime);

  /* Print out the date and time in the standard format. */
  buffer = asctime(loctime);
  for (i = 0; i < plen; i++) {
    if (*buffer) {
      pbuffer[i] = *buffer;
      buffer++;
    } else {
      pbuffer[i] = ' ';
    }
  }
  return 0;
}
