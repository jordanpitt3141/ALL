 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void getufromG(double *h, double *G, double *bhbc, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int bnBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc);
 extern void evolvewrapCONST(double *G, double *h,double *bhbc,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend,double *hMbeg1 , double *hMend1,double *wMbeg1, double *wMend1,double *GMbeg1, double *GMend1,double *uMbeg1, double *uMend1,double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhbc , double theta, double *hhbc, double *whbc,double *Ghbc,double *ubc, double *x, double t);
 extern void RKstep(double *a, double *b, int n);

 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);

 extern double hALLW(double *hbc,int n,double dx) ;
 extern double uhALLW(double *hbc,double *ubc,int n,double dx) ;
 extern double HamilW(double *hbc,double *ubc,int n,double dx) ;




 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void getufromG(double *h, double *G, double *bhbc, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int bnBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc);
 extern void evolvewrapCONST(double *G, double *h,double *bhbc,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend,double *hMbeg1 , double *hMend1,double *wMbeg1, double *wMend1,double *GMbeg1, double *GMend1,double *uMbeg1, double *uMend1,double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhbc , double theta, double *hhbc, double *whbc,double *Ghbc,double *ubc, double *x, double t);
 extern void RKstep(double *a, double *b, int n);

 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);

 extern double hALLW(double *hbc,int n,double dx) ;
 extern double uhALLW(double *hbc,double *ubc,int n,double dx) ;
 extern double HamilW(double *hbc,double *ubc,int n,double dx) ;

