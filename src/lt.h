/* LT_POS to access a lower triangle matrix by C. Buchta
 * modified by M. Hahsler
 * n ... number of rows/columns
 * i,j ... column and row index (starts with 1)
 */

/*  R_xlen_t is for long vector support */

#ifndef LT_POS
#define LT_POS(n, i, j)					\
  (i)==(j) ? 0 : (i)<(j) ? (R_xlen_t)(n) * ((i) - 1) - (R_xlen_t)(i)*((i)-1)/2 + (j)-(i) -1	\
        : (R_xlen_t)(n)*((j)-1) - (R_xlen_t)(j)*((j)-1)/2 + (i)-(j) -1
#endif


/* M_POS to access matrix column-major order by i and j index (starts with 1)
 * n is the number of rows
 */

#ifndef M_POS
#define M_POS(n, i, j) ((i)+(R_xlen_t)(n)*(j))
#endif


/*
 * MIN/MAX
 */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


