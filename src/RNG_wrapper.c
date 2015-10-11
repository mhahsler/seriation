/* FORTRAN Wrapper for R RNG */

#include <R.h>

void F77_SUB(getrngstate)(void) { GetRNGstate(); }
void F77_SUB(putrngstate)(void) { PutRNGstate(); }

/* Note: R's unif_rand returns 0<=x<=1 while FORTRAN's RAND returns 0<=x<1 */
void F77_SUB(unifrand)(float* x) { 
  do{
    *x = (float) unif_rand();
  }while(*x >= 1.0 || *x <0.0);
}
