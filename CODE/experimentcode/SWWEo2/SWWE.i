 /* SWWE.i */
 %module SWWE
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
 extern void evolvewrap(double *uh, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uhMbeg , double *uhMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
 extern void evolvewrapChangeBC(double *u, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend,double *hMbeg1, double *hMend1,double *wMbeg1, double *wMend1,double *uMbeg1, double *uMend1, double *bMbeg1, double *bMend1, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
 extern void evolvewrapIncomwavetoDir(double *u, double *h,double *bed,double WG1ht , double WG1ht1 , double *hMend, double *wMend , double *uMend, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
 %}
 extern double *mallocPy(int n);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx); 
 extern double hall(double *x,double *h,int n, int nBC,double dx);
 extern double uhall(double *x,double *h,double *u,int n, int nBC,double dx);
 extern void evolvewrap(double *u, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
 extern void evolvewrapChangeBC(double *u, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend,double *hMbeg1, double *hMend1,double *wMbeg1, double *wMend1,double *uMbeg1, double *uMend1, double *bMbeg1, double *bMend1, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
 extern void evolvewrapIncomwavetoDir(double *u, double *h,double *bed,double WG1ht , double WG1ht1 , double *hMend, double *wMend , double *uMend, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc);
