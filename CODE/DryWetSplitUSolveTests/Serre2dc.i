 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);

 extern void deallocPy(void *x);
 extern int *RegSplit(double *h,int n);
 extern int readfrom2DmemINT(int *x,int i,int j, int m);
 extern void getufromGsplit(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc, double *bedhbc);






 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);

 extern void deallocPy(void *x);
 extern int *RegSplit(double *h,int n);
 extern int readfrom2DmemINT(int *x,int i,int j, int m);

 extern void getufromGsplit(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc, double *bedhbc);
