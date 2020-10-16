#include <stdio.h>

#define LF 10
#define CR 13

void getword(word, line, stop)
char *word, *line, stop;
{
  int x = 0, y;

  for (x = 0; ((line[x]) && (line[x] != stop)); x++)
    word[x] = line[x];

  word[x] = '\0';
  if (line[x])
    ++x;
  y = 0;

  while (line[y++] = line[x++]);
}

char *makeword(line, stop)
char *line, stop;
{
  int x = 0, y;
  char *word = (char *) malloc(sizeof(char) * (strlen(line) + 1));

  for (x = 0; ((line[x]) && (line[x] != stop)); x++)
    word[x] = line[x];

  word[x] = '\0';
  if (line[x])
    ++x;
  y = 0;

  while (line[y++] = line[x++]);
  return word;
}

char *fmakeword(f, stop, cl)
FILE *f;
char stop;
int *cl;
{
  int wsize;
  char *word;
  int ll;

  wsize = 102400;
  ll = 0;
  word = (char *) malloc(sizeof(char) * (wsize + 1));

  while (1) {
    word[ll] = (char) fgetc(f);
    if (ll == wsize) {
      word[ll + 1] = '\0';
      wsize += 102400;
      word = (char *) realloc(word, sizeof(char) * (wsize + 1));
    }
    --(*cl);
    if ((word[ll] == stop) || (feof(f)) || (!(*cl))) {
      if (word[ll] != stop)
        ll++;
      word[ll] = '\0';
      return word;
    }
    ++ll;
  }
}

char x2c(what)
char *what;
{
  register char digit;

  digit =
      (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A') + 10 : (what[0] - '0'));
  digit *= 16;
  digit +=
      (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A') + 10 : (what[1] - '0'));
  return (digit);
}

unescape_url(url)
char *url;
{
  register int x, y;

  for (x = 0, y = 0; url[y]; ++x, ++y) {
    if ((url[x] = url[y]) == '%') {
      url[x] = x2c(&url[y + 1]);
      y += 2;
    }
  }
  url[x] = '\0';
}

plustospace(str)
char *str;
{
  register int x;

  for (x = 0; str[x]; x++)
    if (str[x] == '+')
      str[x] = ' ';
}

int rind(s, c)
char *s, c;
{
  register int x;
  for (x = strlen(s) - 1; x != -1; x--)
    if (s[x] == c)
      return x;
  return -1;
}

int getline(s, n, f)
char *s;
int n;
FILE *f;
{
  register int i = 0;

  while (1) {
    s[i] = (char) fgetc(f);

    if (s[i] == CR)
      s[i] = fgetc(f);

    if ((s[i] == 0x4) || (s[i] == LF) || (i == (n - 1))) {
      s[i] = '\0';
      return (feof(f) ? 1 : 0);
    }
    ++i;
  }
}

send_fd(f, fd)
FILE *f, *fd;
{
  int num_chars = 0;
  char c;

  while (1) {
    c = fgetc(f);
    if (feof(f))
      return;
    fputc(c, fd);
  }
}

int ind(s, c)
char *s, c;
{
  register int x;

  for (x = 0; s[x]; x++)
    if (s[x] == c)
      return x;

  return -1;
}

escape_shell_cmd(cmd)
char *cmd;
{
  register int x, y, l;

  l = strlen(cmd);
  for (x = 0; cmd[x]; x++) {
    if (ind("&;`'\"|*?~<>^()[]{}$\\\n", cmd[x]) != -1) {
      for (y = l + 1; y > x; y--)
        cmd[y] = cmd[y - 1];
      l++;                      /* length has been increased */
      cmd[x] = '\\';
      x++;                      /* skip the character */
    }
  }
}
