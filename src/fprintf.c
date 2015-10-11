#include <R.h>

void F77_NAME(fprintf) (const char *format, int *nchar,
	double *d1, double *d2)
{

    int nc = *nchar;
    int k;
    char formatc[nc+2]; 

    for (k = 0; k < nc; k++)
	formatc[k] = format[k];
    formatc[nc] = '\0';

    Rprintf(formatc, *d1, *d2);
    Rprintf("\n");
}

