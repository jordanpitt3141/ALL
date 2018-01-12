 /* matrix1.i */
 %module matrix1
 %{
 /* Put header files here or function declarations like below */
 
 extern double *mallocPy(int n);
 extern double *malloc2Py(int n, int m);
 extern double **malloc22Py(int r, int c);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void banmul(double *a, int n, int m1, int m2, double *x, double *b); 
 extern void writeto2mem(double **x, int i, int j , double f);
 extern double readfrom2Dmem(double **x,int i,int j);
 extern void banmul2D(double **a, unsigned long n, int m1, int m2, double *x, double *b);
 extern void banLUP(double **a, unsigned long n, int m1, int m2);
 extern unsigned long *mallocLongPy(int n);

 %}
 extern double *mallocPy(int n);
 extern double *malloc2Py(int n, int m);
 extern void conc(double *a , double *b, double *c,int n,int m ,int k, double *d);
 extern void writetomem(double *x, int i , double f);
 extern double readfrommem(double *x,int i);
 extern void deallocPy(double *x);
 extern void banmul(double *a, int n, int m1, int m2, double *x, double *b);
 extern double **malloc22Py(int r, int c);
 extern void writeto2mem(double **x, int i, int j , double f);
 extern double readfrom2Dmem(double **x,int i,int j);
 extern void banmul2D(double **a, unsigned long n, int m1, int m2, double *x, double *b);
 extern void banLUP(double **a, unsigned long n, int m1, int m2);
 extern unsigned long *mallocLongPy(int n);
