/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
#ifndef CNJJ_H
int cjjini (void);
double c3jj (const int  j1,const int  j2,const int  j3,
             const int  m1,const int  m2,const int  m3);
double c6jj (const int   j1,const int   j2,const int   j3,
             const int  jj1,const int  jj2,const int  jj3);
int    deljj(const int j1,const int j2,const int j3,double *val);
double c9jj (int *jj);
#define CNJJ_H
#endif
