 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh);
 extern void RKstep(double *a, double *b, int n);
 extern void evolvewrapperconsistenttime(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);
 extern double evolvewrapperADAP(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp,double *MollF, double Molls, int m2, int Molllen);
 extern void convolveC(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result);
 extern void convolveS(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result);
 extern void Mollifunc(double *C, double *x, double D, double e, int n);

 extern void convlv0ff(double *data1, unsigned long n, double *respns1, unsigned long m,int isign, double *ans1);
 extern void realft0ff(double *data1, unsigned long n, int isign);
 extern void getufromG(double *h, double *G, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double theta, double dx , int n, int m, int nGhBC,int unBC, int nGhbc, int nubc, double *u, double *hhbc,double *Ghbc, double *MollF, double Molls, int m2, int Molllen);


 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh);
 extern void RKstep(double *a, double *b, int n);
 extern void evolvewrapperconsistenttime(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp);
 extern double evolvewrapperADAP(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp,double *MollF, double Molls, int m2, int Molllen);
 extern void convolveC(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result);
 extern void convolveS(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result);
 extern void Mollifunc(double *C, double *x, double D, double e, int n);

 extern void convlv0ff(double *data1, unsigned long n, double *respns1, unsigned long m,int isign, double *ans1);
 extern void realft0ff(double *data1, unsigned long n, int isign);
 extern void getufromG(double *h, double *G, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double theta, double dx , int n, int m, int nGhBC,int unBC, int nGhbc, int nubc, double *u, double *hhbc,double *Ghbc, double *MollF, double Molls, int m2, int Molllen);
