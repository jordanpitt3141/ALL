 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x); 

 extern void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc, int CnBC , int Cnbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc,double *bbc, double *bMbc, double *dbMbc, double *ddbCbc);

 extern void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc, int CnBC, int Cnbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);



 %}
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);

 extern void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc, int CnBC , int Cnbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc,double *bbc, double *bMbc, double *dbMbc, double *ddbCbc);

 extern void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc, int CnBC, int Cnbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);
