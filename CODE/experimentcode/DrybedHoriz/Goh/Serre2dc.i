 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void getufromG(double *h, double *G, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double theta, double dx , int n, int m, int nGhBC,int unBC, int nGhbc, int nubc, double *u, double *hhbc,double *Ghbc);
 extern void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh);
 extern void RKstep(double *a, double *b, int n);
 extern void evolvewrapperconsistenttime(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);
 extern double evolvewrapperADAP(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);

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
 extern void getufromG(double *h, double *G, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double theta, double dx , int n, int m, int nGhBC,int unBC, int nGhbc, int nubc, double *u, double *hhbc,double *Ghbc);
 extern void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh);
 extern void RKstep(double *a, double *b, int n);
 extern void evolvewrapperconsistenttime(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);
 extern double evolvewrapperADAP(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);

 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);

 extern double hALLW(double *hbc,int n,double dx) ;
 extern double uhALLW(double *hbc,double *ubc,int n,double dx) ;
 extern double HamilW(double *hbc,double *ubc,int n,double dx) ;
