#include <stdio.h>
#include <math.h>
#include "bisnwt.h"

double fexp (double x, void *prm) {
   return exp(x)-2;
}

double dfexp (double x, void *prm) {
   return exp(x);
}

int main (void) {
   double a=-5, b=1, dlt=.1, arr, tol=1e-12;
   int maxit=10;
   bisnwt(a,b,&arr,&dlt,tol,maxit,&fexp,&dfexp,NULL);
   return 0;
}
