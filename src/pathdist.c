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
#include <Rdefines.h>

#include "lt.h"

/* Calculate the path distance for iVAT */
/* Note this changes A! */

/* FIXME: INF and NA */
SEXP pathdist_floyd(SEXP R_x) {
  int *dimX = INTEGER( GET_DIM(R_x) );
  R_xlen_t i, j, k, n = dimX[0];
  SEXP R_y;
  double *x = REAL(R_x);
  double *y;

  PROTECT(R_y = allocMatrix(REALSXP, dimX[0], dimX[1]));
  y = REAL(R_y);

  /* initialize y with paths of length 1 */
  for(i=0; i<n*n; i++) y[i] = x[i];

  /* dynamic programming */
  for (k=0; k<n; k++)
	  for (i=0; i<n; i++)
	    for (j=0; j<n; j++)
		    if (MAX(y[i+n*k], y[k+n*j]) < y[i+n*j])
          y[i+n*j] = MAX(y[i+n*k], y[k+n*j]);

  UNPROTECT(1);
  return R_y;
}
