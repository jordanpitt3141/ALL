 /* Serre2dc.i */
 %module Serre2dc
 %{
 /* Put header files here or function declarations like below */
 extern int readfrom2DmemINT(int *x,int i,int j, int m);
 extern void deallocMem(void* array);
 extern int *append(int *A, int n, int m, int *B);
 extern int *RegSplit(int n);





 %}


 extern int readfrom2DmemINT(int *x,int i,int j, int m);
 extern void deallocMem(void* array);
 extern int *append(int *A, int n, int m, int *B);
 extern int *RegSplit(int n);

