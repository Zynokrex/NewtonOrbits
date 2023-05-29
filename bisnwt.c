#include <stdio.h>
#include "bisnwt.h"
#include <math.h>


int bisnwt (double a, double b, double *arr,
      double *dlt, double tol, int maxit,
      double (*f)(double,void*), double (*df)(double,void*), void *prm) {
    double c,xn,xn1,fa,fc,fx,dfx;
    while(1){
        while (fabs(a - b) > *dlt) {
            c = (a + b) / 2;
            fa = (*f)(a, prm);
            fc = (*f)(c, prm);
            if ((fa * fc) <= 0) b = c;
            else a = c;
        }
        if (*dlt <= tol){
            *arr = c;
            return -1;
        }
        else{
            xn = c;
            for (int i = 1; i <= maxit; i++) {
                fx = (*f)(xn, prm);
                dfx = (*df)(xn, prm);
                xn1 = xn - (fx / dfx);
                if (fabs(xn - xn1) < tol) {
                    *arr = xn1;
                    return i;
                }
                xn=xn1;
            }
            *dlt = *dlt / 2;
        }
    }
}