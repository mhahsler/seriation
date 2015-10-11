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
#include <Rinternals.h>
#include <math.h>

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

    /* since d ist symmetric we only need to sum up half the matrix */
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

SEXP rgar(SEXP R_dist, SEXP R_order, SEXP R_w) {

    int p = INTEGER(getAttrib(R_dist, install("Size")))[0];
    int *o = INTEGER(R_order);
    double *dist = REAL(R_dist);
    /* w ... window size [p,1] */
    int w = INTEGER(R_w)[0];
    
    double d_ij = 0.0;
    double d_ik = 0.0;
    int sum = 0;
    int count = 0;
    
    SEXP R_out;

    /* sum_i=1^p sum_(i-w)<=j<k<i I(d_ij < d_ik) */
    for (int i = 3; i <= p; i++) {
        for(int k = MAX(i-w+1, 2); k < i; k++) {
            d_ik = dist[LT_POS(p, o[i-1], o[k-1])];
            
            for(int j = MAX(i-w, 1); j < k; j++) {
                d_ij = dist[LT_POS(p, o[i-1], o[j-1])];
                
		count++;
                if(d_ij < d_ik) {
		    sum++;
                }    
            }
        }
    }
                    
    /* sum_i=1^p sum_i<j<k<=(i+w) I(d_ij > d_ik) * weight */
    for (int i = 1; i < (p-1); i++) {
        for(int j = i+1; j < MIN(i+w-1, p); j++) {
            d_ij = dist[LT_POS(p, o[i-1], o[j-1])];
            for(int k = j+1; k <= MIN(i+w, p); k++) {
                d_ik = dist[LT_POS(p, o[i-1], o[k-1])];

		count++;
                if(d_ij > d_ik) {
                    sum++;
                }    
            }
        }
    }

    PROTECT(R_out = allocVector(REALSXP, 1));
    REAL(R_out)[0] = (double) sum / (double) count; 
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
