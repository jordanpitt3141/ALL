 /* Serre3.i */
 %module Serre3
 %{
 /* Put header files here or function declarations like below */
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void TDMA(double *a, double *b, double *c, double *d, int n, double *x);
 extern void PENT(double *e, double *a, double *d, double *c,double *f, double *B, int n, double *x);
 extern void midpt2ca(double *qm, double qabeg, double qaend , double dx , int n, double *qa); 
 extern void ca2midpt(double *qa, double qabeg, double qaend, double dx, int n,double *qm);
 extern void ufromGh(double *G, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *u);
 extern void Gfromuh(double *u, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *G);
 extern double phikm(double r);
 extern double phikp(double r);
 extern void weightsum(double a,double *x, double b, double *y, int n, double *z);
 extern void evolvewrap(double *Ga, double *ha, double *Gabeg, double *Gaend, double *habeg, double *haend, double *hmbeg, double *hmend, double *uabeg, double *uaend, double *umbeg, double *umend, double *Gabeg1, double *Gaend1, double *habeg1, double *haend1, double *hmbeg1, double *hmend1, double *uabeg1, double *uaend1, double *umbeg1, double *umend1, double *Gabeg2, double *Gaend2, double *habeg2, double *haend2, double *hmbeg2, double *hmend2, double *uabeg2, double *uaend2, double *umbeg2, double *umend2, int nfcBC, int nGsBC, double g, double dx, double dt, int n, int nBCa, int nBCm, double t, double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);
 %} 
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern double *mallocPy(int n);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void TDMA(double *a, double *b, double *c, double *d, int n, double *x);
 extern void PENT(double *e, double *a, double *d, double *c,double *f, double *B, int n, double *x);
 extern void midpt2ca(double *qm, double qabeg, double qaend , double dx , int n, double *qa); 
 extern void ca2midpt(double *qa, double qabeg, double qaend, double dx, int n,double *qm);
 extern void ufromGh(double *G, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *u);
 extern void Gfromuh(double *u, double *h,double *hbeg,double *hend,double *ubeg,double *uend, double dx , int n, int nBC, double *G);
 extern double phikm(double r);
 extern double phikp(double r);
 extern void weightsum(double a,double *x, double b, double *y, int n, double *z);
 extern void evolvewrap(double *Ga, double *ha, double *Gabeg, double *Gaend, double *habeg, double *haend, double *hmbeg, double *hmend, double *uabeg, double *uaend, double *umbeg, double *umend, double *Gabeg1, double *Gaend1, double *habeg1, double *haend1, double *hmbeg1, double *hmend1, double *uabeg1, double *uaend1, double *umbeg1, double *umend1, double *Gabeg2, double *Gaend2, double *habeg2, double *haend2, double *hmbeg2, double *hmend2, double *uabeg2, double *uaend2, double *umbeg2, double *umend2, int nfcBC, int nGsBC, double g, double dx, double dt, int n, int nBCa, int nBCm, double t, double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx);
