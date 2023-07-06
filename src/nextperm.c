#include <stdio.h>

void    swap(double *, int, int);
void 		permNext(double *, int *);
void 		isMon(double *, double *, int *, int *);

void
swap(double *x, int i, int j)
{
  float           temp;
  temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}


void
permNext(double *x, int *nn)
{
  int             i, j, r, s, n = *nn;
  i = n - 1;
  while (x[i - 1] >= x[i])
    i--;
  if (i == 0)
    return;
  j = n;
  while (x[j - 1] <= x[i - 1])
    j--;
  swap(x, i - 1, j - 1);
  j = n;
  i++;
  while (i < j) {
    swap(x, i - 1, j - 1);
    j--;
    i++;
  }
}

void
isMon(double *x, double *y, int *nn, int *what)
{
  int             n = *nn, i, j;
  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++)
      if (((x[i] - x[j]) * (y[i] - y[j])) <= 0)
        *what = 0;
}
