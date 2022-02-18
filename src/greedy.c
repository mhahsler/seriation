/*
 * seriation - Infrastructure for seriation
 * Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */


#include <R.h>
#include <Rdefines.h>

/* greedy endpoint ordering based on arbitrary similarities.
 * this is trivial.
 *
 * input is a lower triangular distance matrix. returns the
 * merge tree), the corresponding order, and the height (see
 * hclust).
 *
 * note that the height need not be monotonically increasing!
 *
 * (C) ceeboo 2005
 */

typedef struct { double v; int i; } MDS;

static MDS minDist(double *x, int j, int *c, int *p, int n) {

  R_xlen_t i, k, l;
  double v;

  MDS m;

  m.v = R_PosInf;
  l = 0;
  for (k = 0; k < n; k++) {
    i = c[k];
    if (i > j)
      v = x[i+p[j]];
    else
      v = x[j+p[i]];
    if (v < m.v) {
      m.v = v;
      m.i = i;
      l = 1;
    }
    else if (v == m.v) {
      if (unif_rand() > (double) l/(l+1))
        m.i = i;
      l++;
    }
  }

  return m;
}

/* swap */

static void swap(int *x1, int *x2) {

  int x = *x1;

  *x1 = *x2;
  *x2 = x;
}

SEXP order_greedy(SEXP R_dist) {

  R_xlen_t n, i, j, h, k;
  int *left, *right, *order, *c, *p;

  double *x, *height;

  MDS l, ll, r, rr = {0.0, 0};


  SEXP R_obj;

  n =  1 + (int) sqrt(2 * LENGTH(R_dist));

  if (LENGTH(R_dist) != n*(n-1)/2)
    error("order_greedy: \"dist\" invalid length");

  PROTECT(R_obj = NEW_LIST(3));

  SET_ELEMENT(R_obj, 0, allocMatrix(INTSXP, n-1, 2));	/* merge */
  SET_ELEMENT(R_obj, 1, NEW_INTEGER(n));		/* order */
  SET_ELEMENT(R_obj, 2, NEW_NUMERIC(n-1));		/* height */

  left   = INTEGER(VECTOR_ELT(R_obj, 0));
  right  = INTEGER(VECTOR_ELT(R_obj, 0))+n-1;
  order  = INTEGER(VECTOR_ELT(R_obj, 1));
  height =    REAL(VECTOR_ELT(R_obj, 2));

  x = REAL(R_dist);			    /* distance matrix */

  GetRNGstate();

  p = Calloc(n-1, int);		    /* column pointers */
  c = Calloc(n, int);

  for (k = 0; k < n-1; k++) {
    c[k] = k;			    /* candidate leaves */
    p[k] = k*(n-1)-k*(k+1)/2-1;
    order[k] = k;			    /* here backreference */
  }
  c[k] = k;
  order[k] = k;

  i = (int) (unif_rand() * n);	    /* initial leaf */
  h = l.i = ll.i = r.i = rr.i = i;

  for (k = 0; k < n-1; k++) {
    swap(c+order[h], c+n-k-1);
    swap(order+h, order+c[order[h]]);

    if (ll.i == h)
      ll = minDist(x, l.i, c, p, n-k-1);
    if (k == 0)
      rr = ll;
    else if (rr.i == h)
      rr = minDist(x, r.i, c, p, n-k-1);

    if (!R_FINITE(ll.v) || !R_FINITE(rr.v)) {
      Free(c); Free(p);
      error("order_greedy: non-finite values");
    }

    if (ll.v < rr.v) {
      l = ll;
      h = l.i;
      left[k] = -h-1;
      right[k] = k;
      height[k] = l.v;
    }
    else {
      r = rr;
      h = r.i;
      left[k] = k;
      right[k] = -h-1;
      height[k] = r.v;
    }
  }
  left[0] = -i-1;

  /* in each step a leaf was merged. so, we can simply
   * descend the tree and place it on the next left
   * or right position.
   */

  i = 0;
  j = n-1;
  for (k = n-2; k >= 0; k--)
  if (left[k] > 0)
    order[j--] = -right[k];
  else
    order[i++] = -left[k];
  order[j] = -right[0];

  Free(c);
  Free(p);

  PutRNGstate();

  UNPROTECT(1);

  return R_obj;
}

/**/
