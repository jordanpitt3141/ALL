 /* Serre2.i */
 %module Serre2
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *u);
 extern void edgevaluesSplit(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duMbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duMend, double *uEend, double *ddbCend, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double *hMbc, double *GMbc, double *wMbc, double *bMbc, double *duMbc, double *uEbc, double *ddbCbc, double dx, double theta);
 extern void evolvewrapBC(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duEbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duEend, double *uEend, double *ddbCend, double *hMbeg1, double *GMbeg1, double  *wMbeg1, double *duEbeg1, double *uEbeg1, double *hMend1, double *GMend1, double  *wMend1, double *duEend1, double *uEend1, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double dx, double dt, double g, double theta);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);
 %}

 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double minmod(double a, double b, double c);
 extern void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *u);
 extern void edgevaluesSplit(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duMbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duMend, double *uEend, double *ddbCend, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double *hMbc, double *GMbc, double *wMbc, double *bMbc, double *duMbc, double *uEbc, double *ddbCbc, double dx, double theta);
 extern void evolvewrapBC(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duEbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duEend, double *uEend, double *ddbCend, double *hMbeg1, double *GMbeg1, double  *wMbeg1, double *duEbeg1, double *uEbeg1, double *hMend1, double *GMend1, double  *wMend1, double *duEend1, double *uEend1, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double dx, double dt, double g, double theta);

 extern double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double Gall(double *x,double *G,int n, int nBC,double dx);
