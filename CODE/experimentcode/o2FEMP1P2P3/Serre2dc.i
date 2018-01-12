 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern void getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern void evolvewrapForcing(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc, double *x, double h1, double h2, double h3, double u1, double u2, double u3, double b1, double b2, double b3);
 extern void evolvewrap(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc);
 extern void evolveBCChange(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double *hMbeg1 , double *hMend1,double *wMbeg1 , double *wMend1,double *GMbeg1, double *GMend1 ,double *uMbeg1 , double *uMend1 , double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc);
 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern void getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern void evolvewrapForcing(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc, double *x, double h1, double h2, double h3, double u1, double u2, double u3, double b1, double b2, double b3);
 extern void evolvewrap(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc);
 extern void evolveBCChange(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double *hMbeg1 , double *hMend1,double *wMbeg1 , double *wMend1,double *GMbeg1, double *GMend1 ,double *uMbeg1 , double *uMend1 , double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc);
