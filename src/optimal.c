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

/* compute the length of an order, i.e. the sum of
 * the edge weights along the path defined by the
 * order.
 *
 * note that the order is a tour with the leg between
 * the first and the last city omitted.
 *
 * ceeboo 2005
 */

static double orderLength(double *x, int *o, int n) {

    double v, z;
    R_xlen_t i, j, k;

    z = 0;	/* path length */
    i = o[0];
    for (k = 0; k < n-1; k++) {
	j = o[k+1];
	if (i > j)
	   v = x[i+j*(n-1)-j*(j+1)/2-1];
	else
	   if (i == j)
	      return NA_REAL;
           else
	      v = x[j+i*(n-1)-i*(i+1)/2-1];
	if (!R_FINITE(v))
	   return NA_REAL;
	z += v;
	i = j;
    }

    return z;
}

/* R wrapper
 */

SEXP order_length(SEXP R_dist, SEXP R_order) {

    R_xlen_t n, k;
    int *o;

    SEXP R_obj;

    n = LENGTH(R_order);

    if (LENGTH(R_dist) != n * (n - 1) / 2)
       error("order_length: length of \"dist\" and \"order\" do not match");

    o = Calloc(n, int);

    for (k = 0; k < n; k++)		/* offset to C indexing */
	o[k] = INTEGER(R_order)[k]-1;

    PROTECT(R_obj = NEW_NUMERIC(1));

    REAL(R_obj)[0] = orderLength(REAL(R_dist), o, n);
    Free(o);

    UNPROTECT(1);

    return R_obj;
}

/* check validity of a merge tree representation */

int checkRmerge(int *x, int n) {

    R_xlen_t k;
    int v;

    if (x[0] > 0 || x[n-1] > 0)     /* initial merge */
       return 0;

    for (k = 0; k < 2*(n-1); k++) {
        v = x[k];
        if (v < -n || v > n-1)
           return 0;
	if (v > 0 && v > k+1)
	   return 0;
    }

    return 1;
}

/* Z. Bar-Joseph, E. D. Demaine, D. K. Gifford, and T. Jaakkola.
 * (2001) Fast Optimal Leaf Ordering for Hierarchical Clustering.
 * Bioinformatics, Vol. 17 Suppl. 1, pp. 22-29.
 *
 * this implementation builds on the improvements of a more recent paper
 * available at the website of Bar-Joseph!
 *
 * as input we exepct a matrix with the distances in the lower triangle,
 * a merge tree, i.e. two arrays holding n-1 indexes of the left and right
 * subtrees (or leaves) merged at the kth step (for details see dist and
 * hclust).
 *
 * returns a list with a matrix (merge) and two vectors (order and length).
 *
 * The algorithm has the following stages:
 *
 * 1) find a leaf ordering consistent with the supplied merge tree.
 *    the order of the leaves of a tree consists of the order of the
 *    leaves in the left subtree followed by the order of the leaves
 *    in the right subtree.
 *
 * note that the tree (leaf) indexes must have an offset of one because
 * the leaves are coded as negative numbers. subtrees are referenced by
 * their position in the merge sequence (see hclust). this sucks!
 *
 * we compute for each left and right subtree the offset of the leftmost
 * leaf in the total order of leaves, and the number of leaves in both
 * trees, i.e. in the parent tree.
 *
 * 2) recursively compute for each pair of outer endpoints, i.e. a left
 *    endpoint from the left subtree and a right endpoint from the right
 *    subtree the length of the optimal ordering of the leaves.
 *
 * the temporary tables are stored in the lower triangle as well as the
 * similarities. the lengths of the best linear orderings are stored in
 * the upper triangle.
 *
 * for the improved computations at the root the diagonal is used as
 * storage for temporary results.
 *
 * the time complexity of finding all the partial optimal leaf orderings
 * is O(n^3).
 *
 * the suggested improvement based on early termination of the search is
 * currently not implemented. however, ties are broken randomly.
 *
 * 3) recursively find the total optimal leaf ordering.
 *
 * 4) find the merge tree corresponding to the optimal ordering.
 *
 * fixme: using similarities would allow a remapping of non-finite
 *	  values to zero and thus sanitizing of overflows. also for
 *	  missing values this would be a more user friendly approach.
 *
 * (C) ceeboo 2005
 */

static int calcAllOrder(double *x, int *e, int *oi, int *ok, int *oj,
				           int  ci, int  ck, int  cj, int n) {

    R_xlen_t i, ii, j, jj, k, kk, h = 0, l;
    double s, z;

    for (i = 0; i < ci; i++) {
	ii = oi[i];
	for (j = 0; j < cj; j++) {
	    jj = oj[j];
	    l = 0;
	    z = R_PosInf;
	    for (k = 0; k < ck; k++) {
		kk = ok[k];
		if (ii > kk)
		   s  = x[kk+ii*n];
		else
		   s  = x[ii+kk*n];
		if (kk > jj)
		   s += x[kk+jj*n];
		else
		   s += x[jj+kk*n];
		if (s < z) {
		   z = s;
		   h = kk;
		   l = 1;
		}
		else if (s == z) {
			if (unif_rand() > (double) l/(l+1))
			   h = kk;
			l++;
		     }
	    }
	    if (!R_FINITE(z))
	       return 0;		    /* error */

	    if (ii > jj)
	       x[jj+ii*n] = z;
	    else
	       x[ii+jj*n] = z;
	    e[ii+jj*n] = h;
	}
    }
    return 1;
}

static int calcEndOrder(double *x, int *e, int *oi, int *ok,
				           int  ci, int  ck, int n) {

    R_xlen_t i, ii, k, kk, h = 0, l;
    double s, z;

    for (i = 0; i < ci; i++) {
	ii = oi[i];
	l = 0;
	z = R_PosInf;
	for (k = 0; k < ck; k++) {
	    kk = ok[k];
	    if (ii > kk)
		s = x[kk+ii*n];
	    else
		s = x[ii+kk*n];
	    if (s < z) {
		z = s;
		h = kk;
		l = 1;
	    }
	    else if (s == z) {
		    if (unif_rand() > (double) l/(l+1))
		       h = kk;
		    l++;
		 }
	}
	if (!R_FINITE(z))
	   return 0;

	x[ii+ii*n] = z;
	e[ii+ii*n] = h;
    }
    return 1;
}

static int debug = FALSE;

SEXP order_optimal(SEXP R_dist, SEXP R_merge) {

    R_xlen_t n, i, ii, j, jj, k, kk, h, a = 0, b = 0;
    int cl = 0, cll = 0, clr = 0, cr = 0, crl = 0, crr = 0;
    int *l, *r, *c, *e;
    int *left, *right, *o, *ol = 0, *oll = 0, *olr = 0, *or = 0, *orl = 0, *orr = 0;

    double s, z, zz;
    double *x;

    SEXP R_obj;

    n = 1 + (int) sqrt(2 * LENGTH(R_dist));

    if (LENGTH(R_dist) < 3 || LENGTH(R_dist) != n*(n-1)/2)
       error("order_optimal: invalid length");

    if (LENGTH(GET_DIM(R_merge)) != 2)
       error("order_optimal: \"merge\" invalid");

    if (INTEGER(GET_DIM(R_merge))[0] != n-1)
       error("order_optimal: \"dist\" and \"merge\" do not conform");

    if (!checkRmerge(INTEGER(R_merge), n))
       error("order_optimal: \"merge\" invalid");

    /* copy similarities into lower triangle */

    x = Calloc(n*n, double);	    /* data + part order lengths + temporary */

    k = 0;
    for (i = 0; i < n-1; i++)
	for (j = i+1; j < n; j++) {
	    z = REAL(R_dist)[k++];
	    if (!R_FINITE(z)) {
	       Free(x);
	       error("order_optimal: \"dist\" invalid");
	    }
	    else
	       x[j+i*n] = z;
	}

    PROTECT(R_obj = NEW_LIST(3));   /* result list */

    SET_ELEMENT(R_obj, 0, duplicate(R_merge));	/* merge */
    SET_ELEMENT(R_obj, 1, NEW_INTEGER(n));	/* order */
    SET_ELEMENT(R_obj, 2, NEW_NUMERIC(1));	/* length */

    left  = INTEGER(VECTOR_ELT(R_obj, 0));
    right = INTEGER(VECTOR_ELT(R_obj, 0))+n-1;
    o	  = INTEGER(VECTOR_ELT(R_obj, 1));

    GetRNGstate();

    l = Calloc(n,   int);	/* offset of leftmost leaf of left tree */
    r = Calloc(n,   int);	/* offset of leftmost leaf of right tree;
				 * reverse mapping of order */
    c = Calloc(n-1, int);	/* number of leaves in a tree */

    e = Calloc(n*n, int);	/* inner endpoints */

    /* for each tree count the number of leaves.
     */

    for (k = 0; k < n-1; k++) {
	if (left[k] > 0)
	   c[k] += c[left[k]-1];
	else
	   c[k]  = 1;
	if (right[k] > 0)
	   c[k] += c[right[k]-1];
	else
	   c[k] += 1;
    }

    /* backpropagate the counts to obtain the current
     * leaf order and the offset of the leftmost leaf
     * of the left and right subtree.
     */

    for (k = n-2; k >= 0; k--) {
	if (left[k] > 0) {
	   h = l[k] + c[left[k]-1];
	   if (right[k] > 0)
	      l[right[k]-1] = h;
	   else
	      o[h] = -right[k]-1;
	   l[left[k]-1] = l[k];
	}
	else {
	   h = l[k] + 1;
	   if (right[k] > 0)
	      l[right[k]-1] = h;
	   else
	      o[h] = -right[k]-1;
	   o[l[k]] = -left[k]-1;
	}
	r[k] = h;
    }

    /* determine for each subtree the optimal order
     * for each pair of left and right endpoints
     * (leaves). this is done in the order provided
     * by the merge tree.
     */

    for (k = 0; k < n-1; k++) {

	ol = o + l[k];		/* order of left subtree */
	or = o + r[k];		/* order of right subtree */

	cl = r[k] - l[k];	/* number of leaves in left subtree */
	cr = c[k] - cl;		/* number of leaves in right subtree */

	if (cl > 1) {		/* a left tree */
	   h = left[k]-1;

	   oll = o + l[h];
	   olr = o + r[h];

	   cll = r[h] - l[h];
	   clr = c[h] - cll;
	}
	else {			/* a left leaf */
	   oll = olr = ol;
	   cll = clr = cl;
	}
	if (cr > 1) {		/* a right tree */
	   h = right[k]-1;

	   orl = o + l[h];
	   orr = o + r[h];

	   crl = r[h] - l[h];
	   crr = c[h] - crl;
	}
	else {			/* a right leaf */
	   orl = orr = or;
	   crl = crr = cr;
	}

	if (k == n-2)		/* optimized search at the root */
	   break;

	/* compute temporary sums for all endpoints */

	if (!calcAllOrder(x, e, oll, olr, or, cll, clr, cr, n)) {
	   Free(x); Free(r); Free(l); Free(c); Free(e);
	   error("order_optimal: non-finite values");
	}

	if (olr != oll)
	   if (!calcAllOrder(x, e, olr, oll, or, clr, cll, cr, n)) {
	      Free(x); Free(r); Free(l); Free(c); Free(e);
	      error("order_optimal: non-finite values");
	   }

	/* copy temporary sums to lower triangle */

	for (i = 0; i < cl; i++) {
	    ii = ol[i];
	    for (j = 0; j < cr; j++) {
		jj = or[j];
		if (ii > jj)
		   x[ii+jj*n] = x[jj+ii*n];
		else
		   x[jj+ii*n] = x[ii+jj*n];
	    }
	}

	/* compute best orders for all endpoints */

	if (!calcAllOrder(x, e, orl, orr, ol, crl, crr, cl, n)) {
	   Free(x); Free(r); Free(l); Free(c); Free(e);
	   error("order_optimal: non-finite values");
	}

	if (orr != orl)
	   if (!calcAllOrder(x, e, orr, orl, ol, crr, crl, cl, n)) {
	      Free(x); Free(r); Free(l); Free(c); Free(e);
	      error("order_optimal: non-finite values");
	   }

	/* now that we know both endpoints we can store
	 * the inner endpoint from the left tree at the
	 * correct addresse.
	 */

	for (i = 0; i < cr; i++) {
	    ii = or[i];
	    for (j = 0; j < cl; j++) {
	        jj = ol[j];
		kk = e[ii+jj*n];
		if (ii > jj)
		   x[ii+jj*n] = (double) e[jj+kk*n];
		else
	           x[jj+ii*n] = (double) e[jj+kk*n];
	    }
	}

	/* copy back */

	for (i = 0; i < cl; i++) {
	    ii = ol[i];
	    for (j = 0; j < cr; j++) {
		jj = or[j];
		if (ii > jj)
		   e[ii+jj*n] = (int) x[ii+jj*n];
		else
		   e[ii+jj*n] = (int) x[jj+ii*n];
	    }
	}
    }

    /* find the best linear order for each endpoint
     * of the left and right subtree of the root
     */

    if (!calcEndOrder(x, e, oll, olr, cll, clr, n)) {
       Free(x); Free(r); Free(l); Free(c); Free(e);
       error("order_optimal: non-finite values");
    }

    if (olr != oll)
       if (!calcEndOrder(x, e, olr, oll, clr, cll, n)) {
	  Free(x); Free(r); Free(l); Free(c); Free(e);
	  error("order_optimal: non-finite values");
       }

    if (!calcEndOrder(x, e, orl, orr, crl, crr, n)) {
       Free(x); Free(r); Free(l); Free(c); Free(e);
       error("order_optimal: non-finite values");
    }

    if (orr != orl)
       if (!calcEndOrder(x, e, orr, orl, crr, crl, n)) {
	  Free(x); Free(r); Free(l); Free(c); Free(e);
	  error("order_optimal: non-finite values");
       }

    /* find the best linear order at the root */

    k = 0;
    z = R_PosInf;
    for (i = 0; i < cl; i++) {
	ii = ol[i];
	zz = x[ii+ii*n];
	for (j = 0; j < cr; j++) {
	    jj = or[j];
	    s  = zz + x[jj+jj*n];
	    if (ii > jj)
	       s += x[ii+jj*n];
	    else
	       s += x[jj+ii*n];
	    if (s < z) {
	       z = s;
	       a = ii;
	       b = jj;
	       k = 1;
	    }
	    else if (s == z) {
		    if (unif_rand() > (double) k/(k+1)) {
		       a = ii;
		       b = jj;
		    }
		    k++;
		 }
	}
	if (!R_FINITE(z)) {
	   Free(x); Free(r); Free(l); Free(c); Free(e);
	   error("order_optimal: non-finite values");
	}
    }
    REAL(VECTOR_ELT(R_obj, 2))[0] = z;	/* set length */

    /* the order can be found by double recursion.
     * for this we use a stack, one for the left
     * and one for the right endpoints.
     */

    l[0] = b;		    /* push endpoints of right tree on the stack*/
    r[0] = e[b+b*n];

    i = e[a+a*n];	    /* start with endpoints of left tree */
    j = a;

    h = 0;
    k = 1;
    while (h < n) {
	if (i == j) {	    /* backtrack */
	   o[h++] = i;
	   k--;
	   if (k < 0)
	      break;
	   i = l[k];	    /* pop endpoints */
	   j = r[k];
	}
	else {
	   l[k] = e[j+i*n]; /* push endpoints of right tree on the stack */
	   r[k] = j;
	   k++;
	   j = e[i+j*n];    /* recurse left tree */
	}
    }

    /* adjust the merge tree to the optimal order
     *
     * 1) for each pair of leaves from a left and right
     *    subtree the order relation is the same. thus,
     *    use the leftmost leaves as representatives.
     *
     * 2) if the order is reversed we must swap the
     *    subtrees at the parent.
     */

    for (k = 0; k < n; k++)	/* reverse mapping of optimal order */
        r[o[k]] = k;

    for (k = 0; k < n-1; k++) {
	if (left[k] > 0)	/* left leaf in left subtree */
	   i = l[left[k]-1];
	else
	   i = -left[k]-1;
	if (right[k] > 0)	/* left leaf in right subtree */
	   j = l[right[k]-1];
	else
	   j = -right[k]-1;
	if (r[i] > r[j]) {	/* swap the subtrees */
	          h = right[k];
	   right[k] = left[k];
	    left[k] = h;
	}
	l[k] = i;		/* left leaf in parent tree */
    }

    for (k = 0; k < n; k++)	/* offset to R indexing */
	o[k]++;

    if (debug) {
       i = e[a+a*n];
       j = e[b+b*n];

       if (i > j)
          x[j+i*n] = z;
       else
          x[i+j*n] = z;

       for (k = 0; k < n-1; k++) {
	   if (left[k] > 0)
	      l[k] = l[left[k]-1];
	   else
	      l[k] = -left[k]-1;
	   if (right[k] > 0)
	      r[k] = r[right[k]-1];
	   else
	      r[k] = -right[k]-1;

	   i = l[k];
	   j = r[k];
	   if (i > j)
	      z = x[j+i*n];
	   else
	      z = x[i+j*n];


	   // left and right are int
	   // k, i and j are R_xlen_t which is typedefed to ptrdiff_t so we cast to int
	   Rprintf(" %3i | %4i %4i | %3i %3i | %f\n", (int) k+1, left[k], right[k],
						      (int) i+1, (int) j+1, z);
       }
    }

    Free(x);
    Free(l);
    Free(r);
    Free(c);
    Free(e);

    PutRNGstate();

    UNPROTECT(1);

    return R_obj;
}

/**/
