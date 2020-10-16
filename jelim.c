int jelim(t, vec, d, ndm, n)
double *t, *vec, *d;
int ndm, n;
{				/*   do orthogonal transformations to rotate vector VEC into matrix T */
  /*   NDM  = dimensioned column length of T */
  /*   IBGN = first element of VEC */
  /*   N    = length of VEC */

  static double tiny = 1.e-37;
  double c, s, r, rr, as, ac, q, *td;
  int nn;

  q = 1.;
  for (nn = n; nn > 0; --nn, t += ndm) {
    s = (*vec);
    as = fabs(s);
    c = (*t);
    ac = fabs(c);
    if (as > ac) {
      if (q < 1.) {
	dscal(nn + 1, q, vec, 1);
	q = 1.;
      }
      dswap(nn + 1, t, 1, vec, 1);
      as = ac;
      s = (*vec);
      c = (*t);
    }
    td = t;
    ++t;
    ++vec;
    ++d;
    if (as < tiny)
      continue;
    s = -s / c;
    daxpy(nn, s, t, 1, vec, 1);
    s *= q;
    rr = 1. + s * s;
    s *= -q / rr;
    daxpy(nn, s, vec, 1, t, 1);
    if (rr > 1.) {
      r = sqrt(rr);
      q /= r;
      dscal(nn + 1, r, td, 1);
    }
  }
  *d = 1.;
  return 0;
}				/* jelim */
