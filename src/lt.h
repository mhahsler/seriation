
/* LT_POS to access a lower triangle matrix by C. Buchta
 * modified by M. Hahsler
 * n ... number of rows/columns
 * i,j ... column and row index (starts with 1)
 */

#ifndef LT_POS
#define LT_POS(n, i, j)					\
  (i)==(j) ? 0 : (i)<(j) ? n*((i)-1) - (i)*((i)-1)/2 + (j)-(i) -1	\
        : n*((j)-1) - (j)*((j)-1)/2 + (i)-(j) -1
#endif


/* M_POS to access matrix column-major order by i and j index (starts with 1) */

#ifndef M_POS
#define M_POS(n, i, j) ((i)+(n)*(j))
#endif


/*
 * MIN/MAX
 */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


