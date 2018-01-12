#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20


void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

double *mallocPy(int n)
{
    double *x = malloc(n*sizeof(double));
    return x;
}

unsigned long *mallocLongPy(int n)
{
    unsigned long *x = malloc(n*sizeof(unsigned long));
    return x;
}

double *malloc2Py(int n, int m)
{

    double *arr = (double *)malloc(n * m * sizeof(double));
    return arr;
}

double **malloc22Py(int r, int c)
{

    int i;
    double **arr = (double **)malloc(r * sizeof(double *));
    for (i=0; i<r; i++)
    {
         arr[i] = (double *)malloc(c * sizeof(double));
    }
    return arr;
}


void writetomem(double *x, int i , double f)
{
    x[i] = f;

}

void writeto2mem(double **x, int i, int j , double f)
{
    x[i][j] = f;

}


void banmul(double *a, int n, int m1, int m2, double *x, double *b)
{
    int m = m1 + m2 + 1;
    int i,j,si,ei;
    for (i=0;i<n;i++) {
    
        si = 0;
        ei = m1 + m2 + 1;
        if(i < m1 + 1)
        {
            si = m1 - i; 
        }
        
        if(i > n - m2)
        {
            ei = m1 + n - i;
        }
        b[i]=0.0;
        for (j=si ;j<ei;j++) 
        {
            //printf("%d | %d | %d | %d | %f | %f \n",i,j,si,ei,a[i*(m) + j],x[i + j-m1]);
            b[i] += a[i*(m) + j]*x[i -m1 + j];
        }
        
    }


}

long LMIN(long a, long b)
{
    return (b >= a)*a + (b < a)*b;

}

long LMAX(long a, long b)
{
    return (b >= a)*b + (b < a)*a;

}

void banmul2D(double **a, unsigned long n, int m1, int m2, double *x, double *b)
{
    unsigned long i,j,k,tmploop;
    for (i=0;i<n;i++) {
        k=i-m1;
        tmploop=LMIN(m1+m2+1,n-k);
        b[i]=0.0;
        //printf("%lu | %d | %lu | %lu \n", i, m1+m2+1, n - k, tmploop);
        //printf("%lu \n", LMAX(0,-k));
        for (j=LMAX(0,-k);j< tmploop;j++) 
        {
            //printf("%lu | %lu | %f | %f \n", i,j,a[i][j],x[j+k]);
            b[i] += a[i][j]*x[j+k];
        }
    }
   
        
 


}

void banLUP(double **a, unsigned long n, int m1, int m2)
{

    unsigned long i,j,k,l;
    int mm;
    double dum;
    mm=m1+m2+1;

    for (i = 0; i < n ; i ++)
        // Find largest element in a columnb (now a diagonal in reduced format)
        k = i - m1;
        tmploop=LMIN(m1+m2+1,n-k);
 

    /*
    *d=1.0;
    l=m1;
    for (k=0;k<n;k++) 
    {
        //    For each row...
        dum=a[k][1];
        i=k;
        if (l < n) l++;
        for (j=k;j<l;j++) {
        //Find the pivot element.
            if (fabs(a[j][1]) > fabs(dum)) {
                dum=a[j][1];
                i=j;
            }
        }

        indx[k]=i;
        if (dum == 0.0) a[k][1]=TINY;
        //Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in
        //some applications).
        if (i != k) 
        {
            //Interchange rows.
            *d = -(*d);
            for (j=0;j<mm;j++) SWAP(a[k][j],a[i][j])
        }

        for (i=k;i<l;i++) 
        {
            //Do the elimination.
            dum=a[i][1]/a[k][1];
            al[k][i-k]=dum;
            for (j=1;j<mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
            a[i][mm]=0.0;
        }
    }
    */


}


double readfrommem(double*x,int i)
{
    return x[i];
}

double readfrom2Dmem(double **x,int i,int j)
{
    return x[i][j];
}

void deallocPy(double *x)
{
    free(x);
}

int main()
{
    printf("h");
    return 1;
}
