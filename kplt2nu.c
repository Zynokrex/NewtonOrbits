#include <stdio.h>
#include <math.h>
#include "bisnwt.h"


double fexp(double x, void * prm){
    return x-((double *)prm)[1]*sin(x)-((double *)prm)[0];
}
double dfexp(double x, void * prm){
    return 1-((double *)prm)[1]*cos(x);
}
/*
 * Compilar:
	gcc -o kplt2nu -g -Wall kplt2nu.c bisnwt.c -lm
 */

int main (int argc, char *argv[]) {
   double e, T, M0,tf,ti,t,tp,a,b,E,v,dlt,tol,nuacos;
   double prm[2];
   int nt,maxit,k;



   if (argc<6
         || sscanf(argv[1], "%lf", &e)!=1
         || sscanf(argv[2], "%lf", &T)!=1
         || sscanf(argv[3], "%lf", &M0)!=1
         || sscanf(argv[4], "%lf", &tf)!=1
         || sscanf(argv[5], "%d", &nt)!=1
      ) {
      fprintf(stderr,"%s e T M0 tf nt\n", argv[0]);
      return -1;
   }
   prm[0]=M0;
   prm[1]=e;
   t=0;
   k=0;
   ti=tf/nt;
   tol=1e-12;
   maxit=10;
   while(t<tf){
       dlt=2.5;
       a=prm[0]-M_PI;
       b=prm[0]+M_PI;
       bisnwt(a,b,&E,&dlt,tol,maxit,fexp,dfexp,prm);
       nuacos=acos((prm[1]-cos(E))/(prm[1]*cos(E)-1));

       if(E>(k+1)*M_PI) k++;
       if(k%2==0){
           v=nuacos + (k*M_PI);
       }
       else{
           v=(k*M_PI)+(M_PI-nuacos);
       }

       printf("%.16g %.16g %.16g\n",t,prm[0],v); //t M v

       tp=t-((prm[0]/(2*M_PI))*T);
       t=t+ti;//temps
       prm[0]=2*M_PI*((t-tp)/T);//M
   }
   return 0;

}
