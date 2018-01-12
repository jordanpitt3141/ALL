 /* matrixNR.i */
 %module matrixNR
 %{
 /* Put header files here or function declarations like below */
 
 extern double *mallocPy(int n);
 extern double **malloc22Py(int r, int c);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void writeto2mem(double **x, int i, int j , double f);
 extern double readfrom2Dmem(double **x,int i,int j);
 extern unsigned long *mallocLongPy(int n);
 extern void banmul(double **a, unsigned long n, int m1, int m2, double x0[], double b0[]);
 extern double bandec(double **a, unsigned long n, int m1, int m2, double **al,unsigned long indx0[]);
 extern void banbks(double **a, unsigned long n, int m1, int m2, double **al, unsigned long indx0[], double b0[]);
 extern void getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern void getufromGBAND(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern double **dmatrix(long nrl, long nrh, long ncl, long nch);

 %}
 extern double *mallocPy(int n);
 extern double **malloc22Py(int r, int c);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void writeto2mem(double **x, int i, int j , double f);
 extern double readfrom2Dmem(double **x,int i,int j);
 extern unsigned long *mallocLongPy(int n);
 extern void banmul(double **a, unsigned long n, int m1, int m2, double x0[], double b0[]);
 extern double bandec(double **a, unsigned long n, int m1, int m2, double **al,unsigned long indx0[]);
 extern void banbks(double **a, unsigned long n, int m1, int m2, double **al, unsigned long indx0[], double b0[]);
 extern void getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern void getufromGBAND(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend, double *bMbeg, double *bMend, double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbhc, double *u, double *hhbc, double *whbc,double *Ghbc, double *bedhbc);
 extern double **dmatrix(long nrl, long nrh, long ncl, long nch);
