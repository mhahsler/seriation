/*
 * seriation - Infrastructure for seriation
 * Copyright (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdbool.h>

#include "lt.h"


/*
 * path length can be found in optimal.c
 */

/*
 * least-squares criterion
 */

SEXP least_squares_criterion(SEXP R_dist, SEXP R_order) {

  double sum = 0.0;
  int p = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  double *dist = REAL(R_dist);
  double x = 0.0;
  SEXP R_out;

  /* since d is symmetric we only need to sum up half the matrix */
  for (int i = 1; i <= p; i++) {
    for (int j = 1; j < i; j++) {
      x = (dist[LT_POS(p, o[i-1], o[j-1])] - abs(i-j));
      sum += x*x;
    }
  }
  sum *= 2.0;

  PROTECT(R_out = allocVector(REALSXP, 1));
  REAL(R_out)[0] = sum;
  UNPROTECT(1);

  return(R_out);
}

/*
 * inertia criterion
 */

SEXP inertia_criterion(SEXP R_dist, SEXP R_order) {

  double sum = 0.0;
  int p = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  double *dist = REAL(R_dist);
  int x = 0;
  SEXP R_out;

  /* since d ist symmetric we only need to sum up half the matrix */
  for (int i = 1; i <= p; i++) {
    for (int j = 1; j < i; j++) {
      x = abs(i-j);
      sum += dist[LT_POS(p, o[i-1], o[j-1])] * x*x;
    }
  }
  sum *= 2.0;

  PROTECT(R_out = allocVector(REALSXP, 1));
  REAL(R_out)[0] = sum;
  UNPROTECT(1);

  return(R_out);
}

/*
 * Anti-Robinson Events
 */

SEXP ar(SEXP R_dist, SEXP R_order, SEXP R_which) {

  /*
   * which indicates the weighing scheme
   * 1 ... no weighting (i)
   * 2 ... abs. deviations (s)
   * 3 ... weighted abs. deviations (w)
   */

  int p = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  double *dist = REAL(R_dist);
  int which = INTEGER(R_which)[0];

  double sum = 0.0;
  double d_ij = 0.0;
  double d_ik = 0.0;

  SEXP R_out;


  /* sum_i=1^p sum_j<k<i I(d_ij < d_ik) * weight */
  for (int i = 3; i <= p; i++) {
    for(int k = 2; k < i; k++) {
      d_ik = dist[LT_POS(p, o[i-1], o[k-1])];

      for(int j = 1; j < k; j++) {
        d_ij = dist[LT_POS(p, o[i-1], o[j-1])];

        if(d_ij < d_ik) {
          if(which == 1) {
            sum++;
          }else if(which == 2) {
            sum += fabs(d_ij - d_ik);
          }else if(which == 3)
            sum += abs(o[j-1]-o[k-1]) * fabs(d_ij - d_ik);
        }
      }
    }
  }

  /* sum_i=1^p sum_i<j<k I(d_ij > d_ik) * weight */
  for (int i = 1; i < (p-1); i++) {
    for(int j = i+1; j < p; j++) {
      d_ij = dist[LT_POS(p, o[i-1], o[j-1])];
      for(int k = j+1; k <= p; k++) {
        d_ik = dist[LT_POS(p, o[i-1], o[k-1])];

        if(d_ij > d_ik) {
          if(which == 1) {
            sum++;
          }else if(which == 2) {
            sum += fabs(d_ij - d_ik);
          }else if(which == 3)
            sum += abs(o[j-1]-o[k-1]) * fabs(d_ij - d_ik);
        }
      }
    }
  }

  PROTECT(R_out = allocVector(REALSXP, 1));
  REAL(R_out)[0] = sum;
  UNPROTECT(1);

  return(R_out);
}

/*
 * Relative Generalized Anti-Robinson Events
 */

SEXP rgar(SEXP R_dist, SEXP R_order, SEXP R_w, SEXP R_relative) {

  int n = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  int relative = INTEGER(R_relative)[0];
  double *dist = REAL(R_dist);
  /* w is in [2, n-1] (window size) */
  int w = INTEGER(R_w)[0];

  double d_ij = 0.0;
  double d_ik = 0.0;
  int ar = 0;     /* AR events */
  int total = 0;  /* total number of possible AR events */
  int i, j, k;

  SEXP R_out;

  /* sum_i=1^n sum_{(i-w)<=j<k<i} I(d_ij < d_ik) */
  for (i = 3; i <= n; i++) {
    /* Rprintf("i1=%d\n", i); */
    for(k = MAX(i-w+1, 2); k < i; k++) {
      /* Rprintf("k1=%d\n", k); */
      d_ik = dist[LT_POS(n, o[i-1], o[k-1])];
      for(j = MAX(i-w, 1); j < k; j++) {
        /* Rprintf("j1=%d\n\n", j); */
        d_ij = dist[LT_POS(n, o[i-1], o[j-1])];

        total++;
        if(d_ij < d_ik) ar++;
      }
    }
  }

  /* sum_i=1^n sum_i<j<k<=(i+w) I(d_ij > d_ik) * weight */
  for (i = 1; i <= (n-2); i++) {
    /* Rprintf("i2=%d\n", i); */
    for(j = i+1; j <= MIN(i+w-1, n-1); j++) {
      /* Rprintf("j2=%d\n", j); */
      d_ij = dist[LT_POS(n, o[i-1], o[j-1])];
      for(k = j+1; k <= MIN(i+w, n); k++) {
        /* Rprintf("k2=%d\n\n", k); */
        d_ik = dist[LT_POS(n, o[i-1], o[k-1])];

        total++;
        if(d_ij > d_ik) ar++;
      }
    }
  }

  /* Note: total = (2/3-n)*w + n*w^2 - 2/3*w^3 */

  PROTECT(R_out = allocVector(REALSXP, 1));
  if(relative) REAL(R_out)[0] = (double) ar / (double) total;
  else REAL(R_out)[0] = (double) ar;
  UNPROTECT(1);

  return(R_out);
}

/*
 * Gradient Measure
 */

SEXP gradient(SEXP R_dist, SEXP R_order, SEXP R_which) {

  /*
   * which indicates the weighing scheme
   * 1 ... no weighting
   * 2 ... weighted
   */

  int p = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  double *dist = REAL(R_dist);
  int which = INTEGER(R_which)[0];

  double sum = 0.0;
  double d_ij = 0.0;
  double d_ik = 0.0;
  double d_kj = 0.0;
  double diff;

  SEXP R_out;
  int i, k, j;


  /* sum_i<k<j(f(d_ik,d_ij)) */
  for (i = 1; i <= p-2; i++){
    for(k = i+1; k <= p-1; k++) {
      d_ik = dist[LT_POS(p, o[i-1], o[k-1])];

      for(j = k+1; j <= p; j++) {
        d_ij = dist[LT_POS(p, o[i-1], o[j-1])];

        /* diff = d_ik - d_ij; seems to be wrong in the book*/
        diff = d_ij - d_ik;

        if(which > 1) {
          /* weighted */
          sum += diff;
        }else{
          /* unweighted */
          if(diff > 0) sum += 1.0;
          else if(diff < 0) sum -= 1.0;
        }

        /* second sum */
        d_kj = dist[LT_POS(p, o[k-1], o[j-1])];

        /* diff = d_kj - d_ij; seems to be wrong in the book*/
        diff = d_ij - d_kj;


        if(which > 1) {
          /* weighted */
          sum += diff;
        }else{
          /* unweighted */
          if(diff > 0) sum += 1.0;
          else if(diff < 0) sum -= 1.0;
        }

      }
    }
  }


  PROTECT(R_out = allocVector(REALSXP, 1));
  REAL(R_out)[0] = sum;
  UNPROTECT(1);

  return(R_out);
}

/*
 * Lazy Path length (see Earle and Hurley 2015)
 */
SEXP lazy_path_length(SEXP R_dist, SEXP R_order) {

  double tour_length = 0.0;
  SEXP R_tour_length;
  double segment;
  bool posinf = false;
  bool neginf = false;

  int *order = INTEGER(R_order);
  int n = INTEGER(getAttrib(R_dist, install("Size")))[0];

  double *dist = REAL(R_dist);

  if (n != LENGTH(R_order))
    error("length of distance matrix and tour do not match");

  for (int i = 1; i <= n-1; i++) {
    segment = dist[LT_POS(n, order[i-1], order[i])];

    // check Inf
    if (segment == R_PosInf) posinf = true;
    else if (segment == R_NegInf) neginf = true;
    else tour_length +=  (n-i) * segment;
  }

  // do not close tour!

  // inf
  if (posinf && neginf) tour_length = NA_REAL;
  else if (posinf) tour_length = R_PosInf;
  else if (neginf) tour_length = R_NegInf;

  // create R object
  PROTECT(R_tour_length = NEW_NUMERIC(1));
  REAL(R_tour_length)[0] = tour_length;
  UNPROTECT(1);

  return R_tour_length;
}


/*
 * Banded Anti-Robinson Form (see Earle and Hurley, 2015)
 */
SEXP bar(SEXP R_dist, SEXP R_order, SEXP R_b) {

  int n = INTEGER(getAttrib(R_dist, install("Size")))[0];
  int *o = INTEGER(R_order);
  double *dist = REAL(R_dist);
  /* 1 <= b < n */
  int b = INTEGER(R_b)[0];

  double ar = 0;
  int i, j;

  SEXP R_out;

  /* sum_{|i-j|<=b} (b+1-|i-j|) d_{ij} */
  for (i = 1; i <= n-1; i++) {
    for (j = i+1; j <= MIN(i+b, n); j++) {
      ar += (b+1-abs(i-j)) * dist[LT_POS(n, o[i-1], o[j-1])];
    }
  }

  // create R object
  PROTECT(R_out = NEW_NUMERIC(1));
  REAL(R_out)[0] = ar;
  UNPROTECT(1);

  return R_out;
}
