 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);

 extern void evolvewrapperconsistenttimeForced(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend,double *hMbeg1 , double *hMend1,double *GMbeg1 , double *GMend1,double *uMbeg1, double *uMend1, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp, double *x, double t, double *ubc1, double *ubcp1);


 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);

 extern void evolvewrapperconsistenttimeForced(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend,double *hMbeg1 , double *hMend1,double *GMbeg1 , double *GMend1,double *uMbeg1, double *uMend1, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp, double *x, double t, double *ubc1, double *ubcp1);
