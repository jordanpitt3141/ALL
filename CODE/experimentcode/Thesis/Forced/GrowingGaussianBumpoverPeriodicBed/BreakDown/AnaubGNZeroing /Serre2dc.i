 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x); 
 extern void ReconLin( double *q, double *qMbeg, double *qMend, int n, int nBC, int nbc, double theta, double *qbc);
 extern void ReconQuart( double *q, double *qMbeg, double *qMend, int n,int nMBC, int nbcBC, int nbc, double *qbc, double dx);
 extern void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9);
 extern void getufromG(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, double *ubc);

 extern int *RegSplit(double *h,int n);
 extern void getufromGsplit(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, int *RegIndx, int m1 ,int RegIndxLen, double *ubc);
 extern int readfrom2DmemINT(int *x,int i,int j, int m);
 extern void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc,double *bbc);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);



 %}
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);

 extern void ReconLin( double *q, double *qMbeg, double *qMend, int n, int nBC, int nbc, double theta, double *qbc);
 extern void ReconQuart( double *q, double *qMbeg, double *qMend, int n,int nMBC, int nbcBC, int nbc, double *qbc, double dx);
 extern void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9);
 extern void getufromG(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, double *ubc);

 extern int *RegSplit(double *h,int n);
 extern void getufromGsplit(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, int *RegIndx, int m1 ,int RegIndxLen, double *ubc);
 extern int readfrom2DmemINT(int *x,int i,int j, int m);
 extern void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc,double *bbc);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);
