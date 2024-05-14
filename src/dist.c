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

#include "lt.h"

/*
 * Reorder a dist object with a given order
 * Beware: all checking and attribute stuff has to be done in the R wrapper
 */

SEXP reorder_dist(SEXP R_dist, SEXP R_order) {

  SEXP R_dist_out;

  int n = INTEGER(getAttrib(R_dist, install("Size")))[0];
  R_xlen_t n_out = LENGTH(R_order);
  int *o = INTEGER(R_order);

  PROTECT(R_dist_out = allocVector(REALSXP, n_out*(n_out-1)/2));

  double *dist = REAL(R_dist);
  double *dist_out = REAL(R_dist_out);

  for (int i = 1; i <= n_out; i++) {
    for (int j = (i+1); j <=n_out; j++) {

      if(o[i-1] == o[j-1]) dist_out[LT_POS(n_out, i, j)] = 0.0;
      else dist_out[LT_POS(n_out, i, j)] =
        dist[LT_POS(n, o[i-1], o[j-1])];
    }
  }

  UNPROTECT(1);
  return R_dist_out;
}

