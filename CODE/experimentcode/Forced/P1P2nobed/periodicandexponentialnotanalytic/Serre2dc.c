#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

const double i24 = 1.0/24.0;
const double i12 = 1.0/12.0;
const double i16 = 1.0/16.0;
const double i3 = 1.0/3.0;
const double i6 = 1.0/6.0;
const double i8 = 1.0/8.0;
const double i48 = 1.0/48.0;
const double i5 = 1.0/5.0;
const double i42 = 1.0/42.0;
const double i280 = 1.0/280.0;
const double i1680 = 1.0/1680.0;
const double i560 = 1.0/560.0;
const double i140 = 1.0/140.0;
const double i210 = 1.0/210.0;
const double i105 = 1.0 / 105.0;
const double i84 = 1.0/84.0;
const double i70 = 1.0/ 70.0;
const double i80 = 1.0/80.0;
const double i28 = 1.0/28.0;
const double i15 = 1.0/15.0;
const double i10 = 0.1;
const double i9 = 1.0 / 9.0;
const double i120 = 1.0/ 120.0;
const double i20 = 1.0/ 20.0;
const double i35 = 1.0/35.0;
const double i14 = 1.0/ 14.0;
const double i420 = 1.0/420.0;
const double i1344 = 1.0/ 1344.0;
const double i6720 = 1.0 / 6720.0;
const double i2240 = 1.0 / 2240.0;
const double i320 = 1.0 / 320.0;
const double i336 = 1.0 / 336.0;
const double i112 = 1.0 / 122.0;
const double i240 = 1.0 / 240.0;
const double i60 = 1.0 / 60.0;
const double i30 = 1.0 / 30.0;


#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20
#define div_0 1.0e-7
#define ep_h 1.0e-12

#define diffSTAB 1.0e-5

#define NR_END 1


#define FREE_ARG char*

// ############################################################################ START OF NUMERICAL RECIPES CODE #########################################################################

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


long LMIN(long a, long b)
{
    return (b >= a)*a + (b < a)*b;

}

long LMAX(long a, long b)
{
    return (b >= a)*b + (b < a)*a;

}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


void banmul(double **a, unsigned long n, int m1, int m2, double x0[], double b0[])
{
// x, b are not zero offset, as is a
    double *x, *b;
    x = x0 -1;
    b = b0 - 1;
    unsigned long i,j,k,tmploop;
    for (i=1;i<=n;i++) 
    {
        k=i-m1-1;
        tmploop=LMIN(m1+m2+1,n-k);
        b[i]=0.0;

        for (j=LMAX(1,1-k);j<=tmploop;j++) b[i] += a[i][j]*x[j+k];
    }

}

double bandec(double **a, unsigned long n, int m1, int m2, double **al,unsigned long indx0[])
{

    double d = 1.0;
    unsigned long *indx;
    indx = indx0 -1;

    unsigned long i,j,k,l;
    int mm;
    double dum;

    mm=m1+m2+1;
    l=m1;


    for (i=1;i<=m1;i++) 
    {
        //Rearrange the storage a bit.
        for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];

        l--;

        for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
    }

    l=m1;
    for (k=1;k<=n;k++) 
    {
        //For each row...
        dum=a[k][1];
        i=k;
        if (l < n) l++;

        for (j=k+1;j<=l;j++) 
        {
        //Find the pivot element.
            if (fabs(a[j][1]) > fabs(dum)) 
            {
                dum=a[j][1];
                i=j;
            }
        }

        indx[k]=i;

        if ( fabs(dum) <=  div_0) a[k][1]=div_0;
        //Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in   some applications).

        if (i != k) 
        {
        //Interchange rows.
            d = -(d);
            for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
        }

        for (i=k+1;i<=l;i++) 
        {
            //Do the elimination.
            dum=a[i][1]/a[k][1];
            al[k][i-k]=dum;
            for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
            a[i][mm]=0.0;
        }
    }

    return d;

}


void banbks(double **a, unsigned long n, int m1, int m2, double **al, unsigned long indx0[], double b0[])
{

/*Given the arrays a , al , and indx as returned from bandec , and given a right-hand side vector
b[1..n] , solves the band diagonal linear equations A · x = b. The solution vector x overwrites
b[1..n] . The other input arrays are not modified, and can be left in place for successive calls
with different right-hand sides.*/

    unsigned long i,k,l;
    int mm;
    double dum;

    unsigned long *indx;
    indx = indx0 -1;

    double *b;
    b = b0 - 1;

    mm=m1+m2+1;
    l=m1;

    for (k=1;k<=n;k++) 
    {
        //Forward substitution, unscrambling the permuted rows as we go.
        i=indx[k];

        if (i != k) SWAP(b[k],b[i])

        if (l < n) l++;

        for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
    }

    l=1;

    for (i=n;i>=1;i--) 
    {
        //Backsubstitution.
        dum=b[i];
        for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];

        b[i]=dum/a[i][1];

        if (l < mm) l++;
    }
}


// ############################################################################ END OF NUMERICAL RECIPES CODE #########################################################################


// ####################################################################################################### START OF  CODE REQUIRED TO INTERFACE WITH PYTHON ##############################################################
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

double **malloc22Py(int r, int c)
{

   return dmatrix(1,r,1, c);

}


void writetomem(double *x, int i , double f)
{
    x[i] = f;

}

void writeto2mem(double **x, int i, int j , double f)
{
    x[i][j] = f;

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

// ####################################################################################################### END OF  CODE REQUIRED TO INTERFACE WITH PYTHON ##############################################################


/*
 $$$$$$$$$$$$$$$$$$$$$$$$$$$ Hamiltonian Stuff $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/

double interpquarticval(double *coeff,double xj,double x)
{    
    return coeff[0]*(x -xj)*(x -xj)*(x -xj)*(x -xj) + coeff[1]*(x -xj)*(x -xj)*(x -xj)
    + coeff[2]*(x -xj)*(x -xj) + coeff[3]*(x -xj)+ coeff[4];
}  
  
double interpquarticgrad(double *coeff,double xj,double x)
{
    
    return 4*coeff[0]*(x -xj)*(x -xj)*(x -xj) + 3*coeff[1]*(x -xj)*(x -xj)
    + 2*coeff[2]*(x -xj) + coeff[3];
}    
void interpquartcoeff(double *q,double *coeff,int j,double dx)
{
    double idx = 1.0/dx;

    coeff[0] = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2]);
    coeff[1] = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2]);
    coeff[2] = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2]);
    coeff[3] = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2]);
    coeff[4] = q[j];
}

    
double HankEnergyacrosscell(double *x,double *h,double *u,double g,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    double fgpux = interpquarticgrad(ucoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu*fgpu + g*fgph*fgph + i3*(fgph*fgph*fgph)*fgpux*fgpux;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    double sgpux = interpquarticgrad(ucoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu*sgpu + g*sgph*sgph + i3*(sgph*sgph*sgph)*sgpux*sgpux;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
}

double hacrosscell(double *x,double *h,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    
    double fgpe = fgph;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    
    double sgpe = sgph;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    
    double tgpe = tgph;

	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
}

double uhacrosscell(double *x,double *h,double *u,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
}
    
double HankEnergyall(double *x,double *h,double *u,double g,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + HankEnergyacrosscell(x,h,u,g,i,dx);
	}
    return 0.5*sum1; 

}

double hall(double *x,double *h,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + hacrosscell(x,h,i,dx);
	}
    return sum1; 

}

double uhall(double *x,double *h,double *u,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + uhacrosscell(x,h,u,i,dx);
	}
    return sum1; 

}


// New Interpolation

double interplinearval(double *coeff,double xj,double x)
{    
    return coeff[0]*(x -xj) + coeff[1];
}  
      
void interplinearcoeff(double *hbc,double *coeff,int j,double dx)
{
    double idx = 1.0/dx;

    coeff[0] = idx*(hbc[j +1] - hbc[j - 1]);
    coeff[1] = 0.5*(hbc[j +1] + hbc[j - 1]);
}

double interpquadval(double *coeff,double xj,double x)
{    
    return coeff[0]*(x -xj)*(x - xj) + coeff[1]*(x -xj) + coeff[2];
}  
      
void interpquadcoeff(double *ubc,double *coeff,int j,double dx)
{
    double idx = 1.0/dx;

    coeff[0] = 2*idx*idx*(ubc[j +1] - 2*ubc[j] + ubc[j - 1]);
    coeff[1] = idx*(ubc[j +1] - ubc[j - 1]);
    coeff[2] = ubc[j];
}


double minmod(double a, double b, double c)
{
    if((a > 0) && (b>0) && (c>0))
    {
        return fmin(a,fmin(b,c));
    }
    else if((a < 0) && (b<0) && (c<0))
    {
        return fmax(a,fmax(b,c));
    }
        return 0.0;
}

double hacrosscellW(double *hbc,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));

    double mass;
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,j,dx);
    
    //intgxjph = 0.5*hcoeff[0]*0.25*dx*dx + hcoeff[1]*0.5*dx;
    //intgxjmh = 0.5*hcoeff[0]*0.25*dx*dx - hcoeff[1]*0.5*dx;

    mass = 2*hcoeff[1]*0.5*dx;

	free(hcoeff);
    
    return mass;
}

double uhacrosscellW(double *hbc,double *ubc,int uj, int hj,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));
    double *ucoeff = malloc(3*sizeof(double));

    double mass;
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,hj,dx);
    interpquadcoeff(ubc,ucoeff,uj,dx);
    
    /*intgxjph = 0.25*hcoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) + i3*(hcoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx) + 0.5*(hcoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)
               + hcoeff[1]*ucoeff[2]*(0.5*dx);


    intgxjmh = 0.25*hcoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) - i3*(hcoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx) + 0.5*(hcoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)
           - hcoeff[1]*ucoeff[2]*(0.5*dx);*/

    mass = 2*(i3*(hcoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[1]*ucoeff[2]*(0.5*dx));

	free(hcoeff);
    free(ucoeff);
    
    return mass;
}

double HamilacrosscellW(double *hbc,double *ubc,int uj, int hj,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));
    double *ucoeff = malloc(3*sizeof(double));

    double u2h, h2, h3ux2;
    
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,hj,dx);
    interpquadcoeff(ubc,ucoeff,uj,dx);


    /*u2hjph = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             + i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;

    u2hjmh = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             - i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             - i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             - hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;*/

    u2h  =2*( i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
           + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx)); 
    
    //h2jph = i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx);

    //h2jmh = -i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) - hcoeff[1]*hcoeff[1]*(0.5*dx);

    h2  = 2*( i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx)) ;

    /*h3ux2jph = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);


    h3ux2jmh = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - i5*(4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            - hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);*/

    h3ux2 =  2*(i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
              + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
              + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx) );

	free(hcoeff);
    free(ucoeff);
    
    return 0.5*(u2h + 9.81*h2 + i3*h3ux2);
}

double HamilKEHacrosscellW(double *hbc,double *ubc,int uj, int hj,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));
    double *ucoeff = malloc(3*sizeof(double));

    double u2h, h2, h3ux2;
    
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,hj,dx);
    interpquadcoeff(ubc,ucoeff,uj,dx);


    /*u2hjph = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             + i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;

    u2hjmh = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             - i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             - i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             - hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;*/

    u2h  =2*( i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
           + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx)); 
    
    //h2jph = i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx);

    //h2jmh = -i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) - hcoeff[1]*hcoeff[1]*(0.5*dx);

    h2  = 2*( i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx)) ;

    /*h3ux2jph = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);


    h3ux2jmh = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - i5*(4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            - hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);*/

    h3ux2 =  2*(i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
              + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
              + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx) );

	free(hcoeff);
    free(ucoeff);
    
    return 0.5*(u2h);
}

double HamilKEVacrosscellW(double *hbc,double *ubc,int uj, int hj,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));
    double *ucoeff = malloc(3*sizeof(double));

    double u2h, h2, h3ux2;
    
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,hj,dx);
    interpquadcoeff(ubc,ucoeff,uj,dx);


    /*u2hjph = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             + i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;

    u2hjmh = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             - i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             - i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             - hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;*/

    u2h  =2*( i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
           + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx)); 
    
    //h2jph = i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx);

    //h2jmh = -i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) - hcoeff[1]*hcoeff[1]*(0.5*dx);

    h2  = 2*( i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx)) ;

    /*h3ux2jph = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);


    h3ux2jmh = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - i5*(4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            - hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);*/

    h3ux2 =  2*(i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
              + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
              + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx) );

	free(hcoeff);
    free(ucoeff);
    
    return 0.5*(i3*h3ux2);
}

double HamilGPEacrosscellW(double *hbc,double *ubc,int uj, int hj,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(2*sizeof(double));
    double *ucoeff = malloc(3*sizeof(double));

    double u2h, h2, h3ux2;
    
	

    //jth cell
    interplinearcoeff(hbc,hcoeff,hj,dx);
    interpquadcoeff(ubc,ucoeff,uj,dx);


    /*u2hjph = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             + i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;

    u2hjmh = i6*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)
             - i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             + (0.5*hcoeff[0]*ucoeff[0]*ucoeff[2] + 0.25*hcoeff[0]*ucoeff[1]*ucoeff[1] + 0.5*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
             - i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
             + (0.5*hcoeff[0]*ucoeff[2]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[2])*(0.5*dx)*(0.5*dx)   
             - hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx) ;*/

    u2h  =2*( i5*(2*hcoeff[0]*ucoeff[0]*ucoeff[1] + hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
           + i3*(2*hcoeff[0]*ucoeff[1]*ucoeff[2] + 2*hcoeff[1]*ucoeff[0]*ucoeff[2] + hcoeff[1]*ucoeff[1]*ucoeff[1])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + hcoeff[1]*ucoeff[2]*ucoeff[2]*(0.5*dx)); 
    
    //h2jph = i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx);

    //h2jmh = -i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[0]*hcoeff[1]*(0.5*dx)*(0.5*dx) - hcoeff[1]*hcoeff[1]*(0.5*dx);

    h2  = 2*( i3*hcoeff[0]*hcoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx) + hcoeff[1]*hcoeff[1]*(0.5*dx)) ;

    /*h3ux2jph = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);


    h3ux2jmh = 2*i3*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[0]*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - i5*(4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (0.25*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[1]*ucoeff[1] + 3*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 3*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx) 
            - (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx)  
            + (3*0.5*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 2*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1])*(0.5*dx)*(0.5*dx)
            - hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx);*/

    h3ux2 =  2*(i5*( 4*hcoeff[0]*hcoeff[0]*hcoeff[0]*ucoeff[0]*ucoeff[1] + 12*hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[0]*ucoeff[0] )*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)*(0.5*dx)  
              + (hcoeff[0]*hcoeff[0]*hcoeff[1]*ucoeff[1]*ucoeff[1] + 4*hcoeff[0]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[1] + 4*i3*hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[0]*ucoeff[0])*(0.5*dx)*(0.5*dx)*(0.5*dx) 
              + hcoeff[1]*hcoeff[1]*hcoeff[1]*ucoeff[1]*ucoeff[1]*(0.5*dx) );

	free(hcoeff);
    free(ucoeff);
    
    return 0.5*(9.81*h2);
}


double hALLW(double *hbc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + hacrosscellW(hbc,3*i + 1,dx);

    }


    return sum; 
}

double uhALLW(double *hbc,double *ubc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + uhacrosscellW(hbc,ubc,2*i + 1, 3*i + 1,dx);

    }


    return sum; 
}


double HamilW(double *hbc,double *ubc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + HamilacrosscellW(hbc,ubc,2*i + 1, 3*i + 1,dx);

    }


    return sum; 
}

double HamilKEHW(double *hbc,double *ubc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + HamilKEHacrosscellW(hbc,ubc,2*i + 1, 3*i + 1,dx);

    }


    return sum; 
}

double HamilKEVW(double *hbc,double *ubc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + HamilKEVacrosscellW(hbc,ubc,2*i + 1, 3*i + 1,dx);

    }


    return sum; 
}

double HamilGPEW(double *hbc,double *ubc,int n,double dx) 
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere
    int i;
    double sum = 0;

    for(i = 1; i < n + 1 ; i ++)
    {

        sum = sum + HamilGPEacrosscellW(hbc,ubc,2*i + 1, 3*i + 1,dx);

    }


    return sum; 
}

// Solver Code
void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void RKstep(double *a, double *b, int n)
{
    int i;
    for(i =0;i <n;i++)
    {
        a[i] = 0.5*(a[i] + b[i]);

    }

}

void getufromG(double *h, double *G, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double theta, double dx , int n, int m, int nGhBC,int unBC, int nGhbc, int nubc, double *u, double *hhbc,double *Ghbc)
{
    //trying to join different B.C systems...., also try incorporating into  C code that already has beds, like o2bedfix.

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now

    double idx = 1.0/ dx;
    //double *uais = malloc((m-2)*sizeof(double));
    //double *ubis = malloc((m-1)*sizeof(double));
    //double *ucis = malloc((m)*sizeof(double));
    //double *udis = malloc((m-1)*sizeof(double));
    //double *ueis = malloc((m-2)*sizeof(double));

    double **AL, **Ared;
    unsigned long *indx;

    Ared = dmatrix(1, m, 1, 5);
    AL = dmatrix(1, m, 1, 5);
    indx = mallocLongPy(m);

    double *nGis = malloc((m)*sizeof(double));
    //double *utemp = malloc((m)*sizeof(double));

    int i,j;

    double dGib,dGim,dGif,dhib,dhim,dhif,dGi,dhi,Gjphm,Gjmhp,hjphm,hjmhp;

    double Gintia11,Gintia21,Gintia31,uhintia11,uhintia12,uhintia13,uhintia21,uhintia22,uhintia23,uhintia31,uhintia32,uhintia33;
    double h3uxintia11,h3uxintia12,h3uxintia13,h3uxintia21,h3uxintia22,h3uxintia23,h3uxintia31,h3uxintia32,h3uxintia33;
    double LHSa11,LHSa12,LHSa13,LHSa21,LHSa22,LHSa23,LHSa31,LHSa32,LHSa33;

    // Zero out the diagonal arrays
    for(i = 0; i < m ; i++)
    {
        Ared[i + 1][1] = 0;
        Ared[i + 1][2] = 0;
        Ared[i + 1][3] = 0;
        Ared[i + 1][4] = 0;
        Ared[i + 1][5] = 0;
        
        AL[i + 1][1] = 0;
        AL[i + 1][2] = 0;
        AL[i + 1][3] = 0;
        AL[i + 1][4] = 0;
        AL[i + 1][5] = 0;

        nGis[i] = 0;
    
    }



    j = 3;
    for (i =1;i < n-1 ; i++)
    {
        // Reconstruct G and h
        dGib = (G[i] - G[i-1]);
        dGim = 0.5*(G[i+1] - G[i-1]);
        dGif = (G[i+1] - G[i]);

        dhib = (h[i] - h[i-1]);
        dhim = 0.5*(h[i+1] - h[i-1]);
        dhif = (h[i+1] - h[i]);

        dGi = minmod(theta*dGib, dGim, theta*dGif);
        dhi = minmod(theta*dhib, dhim, theta*dhif);

        hjphm=  h[i] + 0.5*dhi; 
        hjmhp=  h[i] - 0.5*dhi;

        Gjphm= G[i] + 0.5*dGi;
        Gjmhp= G[i] - 0.5*dGi; 



        // G integral (RHS)
        Gintia11 = i6*dx*(Gjmhp);
        Gintia21 = i6*dx*(2*Gjmhp + 2*Gjphm);
        Gintia31 = i6*dx*(Gjphm);
               
        //uh integral
        uhintia11 = (1.0 *i60)*dx*(7*hjmhp + hjphm);
        uhintia12 = (1.0 *i60)*dx*(4*hjmhp );
        uhintia13 = (1.0 *i60)*dx*(-hjmhp - hjphm);
        
        uhintia21 = (1.0 *i60)*dx*(4*hjmhp);
        uhintia22 = (1.0 *i60)*dx*(16*hjmhp + 16*hjphm);
        uhintia23 = (1.0 *i60)*dx*(4*hjphm);
        
        uhintia31 = (1.0 *i60)*dx*(-hjmhp - hjphm);
        uhintia32 = (1.0 *i60)*dx*(4*hjphm);
        uhintia33 = (1.0 *i60)*dx*(hjmhp + 7*hjphm);
        
        //h3ux
        
        h3uxintia11 = (2.0*i3)*idx*((79.0*i120)*hjmhp*hjmhp*hjmhp +  (39.0*i120)*hjmhp*hjmhp*hjphm + (3.0*i24)*hjmhp*hjphm*hjphm + (7.0*i120)*hjphm*hjphm*hjphm);        
        h3uxintia12 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
        h3uxintia13 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);
        
        
        h3uxintia21 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
        h3uxintia22 = (2.0*i3)*idx*((14.0*i15)*hjmhp*hjmhp*hjmhp +  (6.0*i15)*hjmhp*hjmhp*hjphm + (6.0*i15)*hjmhp*hjphm*hjphm + (14.0*i15)*hjphm*hjphm*hjphm);        
        h3uxintia23 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);
        
        h3uxintia31 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);       
        h3uxintia32 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);        
        h3uxintia33 = (2.0*i3)*idx*((7.0*i120)*hjmhp*hjmhp*hjmhp +  (3.0*i24)*hjmhp*hjmhp*hjphm + (39.0*i120)*hjmhp*hjphm*hjphm + (79.0*i120)*hjphm*hjphm*hjphm);
        
        
        // LHS 
        
        LHSa11 = uhintia11 + h3uxintia11;
        LHSa12 = uhintia12 + h3uxintia12; 
        LHSa13 = uhintia13 + h3uxintia13;
        LHSa21 = uhintia21 + h3uxintia21;
        LHSa22 = uhintia22 + h3uxintia22; 
        LHSa23 = uhintia23 + h3uxintia23;
        LHSa31 = uhintia31 + h3uxintia31;
        LHSa32 = uhintia32 + h3uxintia32;
        LHSa33 = uhintia33 + h3uxintia33;

        Ared[j-1 + 3][1] = Ared[j-1 + 3][1] + LHSa31;

        Ared[j-1 +2][2] = Ared[j-1 + 2][2] + LHSa21;
        Ared[j + 2][2] = Ared[j + 2][2] + LHSa32;

        Ared[j-1 + 1][3] = Ared[j-1 + 1][3] + LHSa11;
        Ared[j + 1][3] = Ared[j + 1][3] + LHSa22;
        Ared[j+1 + 1][3] = Ared[j+1 + 1][3] + LHSa33;

        Ared[j-1 + 1][4] = Ared[j-1 + 1][4] + LHSa12;
        Ared[j + 1][4] = Ared[j + 1][4] + LHSa23;
        
        Ared[j-1 + 1 ][5] = Ared[j-1 + 1][5] + LHSa13;


        nGis[j-1] = nGis[j-1] + Gintia11;
        nGis[j] = nGis[j] + Gintia21; 
        nGis[j+1] = nGis[j+1] + Gintia31;


        hhbc[3*i + (nGhBC)] = hjmhp;
        hhbc[3*i+1 + (nGhBC)] =  h[i];
        hhbc[3*i+2 + (nGhBC)] = hjphm; 

        Ghbc[3*i + (nGhBC)] = Gjmhp;
        Ghbc[3*i+1 + (nGhBC)] = G[i];
        Ghbc[3*i+2 + (nGhBC)] = Gjphm; 
      
        
        j = j + 2 ;       


    }


    // first 
    i = 0;
    j = 1;

    // Get it from interior
    // Reconstruct w

    // Reconstruct G and h
    hjphm= hhbc[i + 2*nGhBC];
    hjmhp= hMbeg[nGhBC-1];

    Gjphm= Ghbc[i + 2*nGhBC];
    Gjmhp= GMbeg[nGhBC-1];


    // G integral (RHS)
    Gintia11 = dx*i6*(Gjmhp);
    Gintia21 = dx*i6*(2*Gjmhp + 2*Gjphm);
    Gintia31 = dx*i6*(Gjphm);
    
    
    //uh integral
    uhintia11 = (1.0 *i60)*dx*(7*hjmhp + hjphm);
    uhintia12 = (1.0 *i60)*dx*(4*hjmhp );
    uhintia13 = (1.0 *i60)*dx*(-hjmhp - hjphm);
    
    uhintia21 = (1.0 *i60)*dx*(4*hjmhp);
    uhintia22 = (1.0 *i60)*dx*(16*hjmhp + 16*hjphm);
    uhintia23 = (1.0 *i60)*dx*(4*hjphm);
    
    uhintia31 = (1.0 *i60)*dx*(-hjmhp - hjphm);
    uhintia32 = (1.0 *i60)*dx*(4*hjphm);
    uhintia33 = (1.0 *i60)*dx*(hjmhp + 7*hjphm);
    
    //h3ux
    
    h3uxintia11 = (2.0*i3)*idx*((79.0*i120)*hjmhp*hjmhp*hjmhp +  (39.0*i120)*hjmhp*hjmhp*hjphm + (3.0*i24)*hjmhp*hjphm*hjphm + (7.0*i120)*hjphm*hjphm*hjphm);        
    h3uxintia12 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
    h3uxintia13 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);
    
    
    h3uxintia21 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
    h3uxintia22 = (2.0*i3)*idx*((14.0*i15)*hjmhp*hjmhp*hjmhp +  (6.0*i15)*hjmhp*hjmhp*hjphm + (6.0*i15)*hjmhp*hjphm*hjphm + (14.0*i15)*hjphm*hjphm*hjphm);        
    h3uxintia23 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);
    
    h3uxintia31 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);       
    h3uxintia32 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);        
    h3uxintia33 = (2.0*i3)*idx*((7.0*i120)*hjmhp*hjmhp*hjmhp +  (3.0*i24)*hjmhp*hjmhp*hjphm + (39.0*i120)*hjmhp*hjphm*hjphm + (79.0*i120)*hjphm*hjphm*hjphm);
    
    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11;
    LHSa12 = uhintia12 + h3uxintia12; 
    LHSa13 = uhintia13 + h3uxintia13;
    LHSa21 = uhintia21 + h3uxintia21;
    LHSa22 = uhintia22 + h3uxintia22; 
    LHSa23 = uhintia23 + h3uxintia23;
    LHSa31 = uhintia31 + h3uxintia31;
    LHSa32 = uhintia32 + h3uxintia32;
    LHSa33 = uhintia33 + h3uxintia33;


    Ared[j-1 + 3][1] = Ared[j-1 + 3][1] + LHSa31;

    Ared[j-1 +2][2] = Ared[j-1 + 2][2] + LHSa21;
    Ared[j + 2][2] = Ared[j + 2][2] + LHSa32;

    Ared[j-1 + 1][3] = 1;
    Ared[j + 1][3] = Ared[j + 1][3] + LHSa22;
    Ared[j+1 + 1][3] = Ared[j+1 + 1][3] + LHSa33;

    Ared[j-1 + 1][4] = 0;
    Ared[j + 1][4] = Ared[j + 1][4] + LHSa23;
    
    Ared[j-1 + 1 ][5] = 0;

    nGis[j-1] = uMbeg[unBC-1];
    nGis[j] = nGis[j] + Gintia21; 
    nGis[j+1] = nGis[j+1] + Gintia31;     

    hhbc[3*i + (nGhBC)] = hjmhp;
    hhbc[3*i+1 + (nGhBC)] = 0.5*(hjmhp +  hjphm);
    hhbc[3*i+2 + (nGhBC)] = hjphm; 


    Ghbc[3*i + (nGhBC)] = Gjmhp;
    Ghbc[3*i+1 + (nGhBC)] = 0.5*(Gjmhp +  Gjphm);
    Ghbc[3*i+2 + (nGhBC)] = Gjphm;
  

// last
    j = m-2;
    i = n-1;

    // Get it from interior 
    hjphm= hMend[0];
    hjmhp= hhbc[nGhbc - nGhBC -4];

    Gjphm= GMend[0];
    Gjmhp= Ghbc[nGhbc - nGhBC -4];

    // G integral (RHS)
    Gintia11 = dx*i6*(Gjmhp);
    Gintia21 = dx*i6*(2*Gjmhp + 2*Gjphm);
    Gintia31 = dx*i6*(Gjphm);
    
    
    //uh integral
    uhintia11 = (1.0 *i60)*dx*(7*hjmhp + hjphm);
    uhintia12 = (1.0 *i60)*dx*(4*hjmhp );
    uhintia13 = (1.0 *i60)*dx*(-hjmhp - hjphm);
    
    uhintia21 = (1.0 *i60)*dx*(4*hjmhp);
    uhintia22 = (1.0 *i60)*dx*(16*hjmhp + 16*hjphm);
    uhintia23 = (1.0 *i60)*dx*(4*hjphm);
    
    uhintia31 = (1.0 *i60)*dx*(-hjmhp - hjphm);
    uhintia32 = (1.0 *i60)*dx*(4*hjphm);
    uhintia33 = (1.0 *i60)*dx*(hjmhp + 7*hjphm);
    
    //h3ux
    
    h3uxintia11 = (2.0*i3)*idx*((79.0*i120)*hjmhp*hjmhp*hjmhp +  (39.0*i120)*hjmhp*hjmhp*hjphm + (3.0*i24)*hjmhp*hjphm*hjphm + (7.0*i120)*hjphm*hjphm*hjphm);        
    h3uxintia12 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
    h3uxintia13 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);
    
    
    h3uxintia21 = (2.0*i3)*idx*((-23.0*i30)*hjmhp*hjmhp*hjmhp +  (-3.0*i10)*hjmhp*hjmhp*hjphm + (-3.0*i30)*hjmhp*hjphm*hjphm + (-1.0*i6)*hjphm*hjphm*hjphm);        
    h3uxintia22 = (2.0*i3)*idx*((14.0*i15)*hjmhp*hjmhp*hjmhp +  (6.0*i15)*hjmhp*hjmhp*hjphm + (6.0*i15)*hjmhp*hjphm*hjphm + (14.0*i15)*hjphm*hjphm*hjphm);        
    h3uxintia23 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);
    
    h3uxintia31 = (2.0*i3)*idx*((13.0*i120)*hjmhp*hjmhp*hjmhp +  (-3.0*i120)*hjmhp*hjmhp*hjphm + (-3.0*i120)*hjmhp*hjphm*hjphm + (13.0*i120)*hjphm*hjphm*hjphm);       
    h3uxintia32 = (2.0*i3)*idx*((-1.0*i6)*hjmhp*hjmhp*hjmhp +  (-3.0*i30)*hjmhp*hjmhp*hjphm + (-3.0*i10)*hjmhp*hjphm*hjphm + (-23.0*i30)*hjphm*hjphm*hjphm);        
    h3uxintia33 = (2.0*i3)*idx*((7.0*i120)*hjmhp*hjmhp*hjmhp +  (3.0*i24)*hjmhp*hjmhp*hjphm + (39.0*i120)*hjmhp*hjphm*hjphm + (79.0*i120)*hjphm*hjphm*hjphm);
    
    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11;
    LHSa12 = uhintia12 + h3uxintia12; 
    LHSa13 = uhintia13 + h3uxintia13;
    LHSa21 = uhintia21 + h3uxintia21;
    LHSa22 = uhintia22 + h3uxintia22; 
    LHSa23 = uhintia23 + h3uxintia23;
    LHSa31 = uhintia31 + h3uxintia31;
    LHSa32 = uhintia32 + h3uxintia32;
    LHSa33 = uhintia33 + h3uxintia33;
       

    Ared[j-1 + 3][1] = 0;

    Ared[j-1 +2][2] = Ared[j-1 + 2][2] + LHSa21;
    Ared[j + 2][2] = 0;

    Ared[j-1 + 1][3] = Ared[j-1 + 1][3] + LHSa11;
    Ared[j + 1][3] = Ared[j + 1][3] + LHSa22;
    Ared[j+1 + 1][3] = 1;

    Ared[j-1 + 1][4] = Ared[j-1 + 1][4] + LHSa12;
    Ared[j + 1][4] = Ared[j + 1][4] + LHSa23;
    
    Ared[j-1 + 1 ][5] = Ared[j-1 + 1][5] + LHSa13;

    
    nGis[j-1] = nGis[j-1] + Gintia11; 
    nGis[j] = nGis[j] + Gintia21; 
    nGis[j+1] = uMend[0];

    hhbc[3*i + (nGhBC)] = hjmhp;
    hhbc[3*i+1 + (nGhBC)] = 0.5*(hjmhp + hjphm);
    hhbc[3*i+2 + (nGhBC)] = hjphm; 

    Ghbc[3*i + (nGhBC)] = Gjmhp;
    Ghbc[3*i+1 + (nGhBC)] = 0.5*(Gjmhp + Gjphm);
    Ghbc[3*i+2 + (nGhBC)] = Gjphm;


    double d;
    d = bandec(Ared, m, 2,2, AL ,indx);

    banbks(Ared, m,2,2, AL, indx, nGis);



    memcpy(u,uMbeg, (unBC-1)*sizeof(double));
    memcpy(u+(unBC-1),nGis,m*sizeof(double));
    memcpy(u+nubc - (unBC-1),uMend + 1,(unBC-1)*sizeof(double));

    // B.C stuff
    for(i=0;i < nGhBC;i++)
    {
        // Assuming bed constant in ghost cells
        hhbc[i] = hMbeg[i];
        hhbc[nGhbc-nGhBC + i] = hMend[i];
        Ghbc[i] = GMbeg[i];
        Ghbc[nGhbc-nGhBC +i] =GMend[i];
    }

    free_dmatrix(Ared, 1, m, 1, 5);
    free_dmatrix(AL, 1, m, 1, 5);
    free(indx);
    free(nGis);

}


double barhxt(double x, double t)
{

    return 10.0*x - exp(0.5*t)*cos(2*x)/2;

}

double barGxt(double x, double t)
{
    return (0.25*exp(0.5*t)*sin(2*x)*sin(2*x) 
        + 5.0*sin(2*x))*exp(0.5*t) 
        + 4*(3.33333333333333*exp(0.5*t)*sin(2*x)*sin(2*x)*sin(2*x) 
        - 0.25*exp(1.0*t)*sin(2*x)*sin(2*x)*cos(2*x)*cos(2*x)
        - 0.125*exp(1.0*t)*cos(2*x)*cos(2*x)*cos(2*x)*cos(2*x) 
        + 25.0*sin(2*x)*sin(2*x))*exp(1.0*t) 
        + 4*(75.0*exp(0.5*t)*sin(2*x)*sin(2*x) 
        + 5.0*exp(1.0*t)*sin(2*x)*sin(2*x)*sin(2*x) 
        - 0.25*exp(1.5*t)*sin(2*x)*sin(2*x)*cos(2*x)*cos(2*x) 
        - 0.125*exp(1.5*t)*cos(2*x)*cos(2*x)*cos(2*x)*cos(2*x) 
        + 500.0*sin(2*x))*exp(0.5*t)/3;
}


double Fluxh(double x, double t)
{
    return 20.0*exp(0.5*t)*cos(2*x) + 1.0*exp(1.0*t)*sin(2*x)*cos(2*x);
}

double FluxG(double x, double t, double g)
{
    return 490.5*t 
        + (1343.33333333333*exp(1.0*t)*cos(2*x) 
        + 534.0*exp(1.5*t)*sin(2*x)*cos(2*x) 
        + 60.0*exp(2.0*t)*sin(2*x)*sin(2*x)*cos(2*x) 
        + 2.13333333333333*exp(2.5*t)*sin(2*x)*sin(2*x)*sin(2*x)*cos(2*x))*cos(2*x) 
        - 2.66666666666667*(1000.0*exp(1.0*t) + 200.0*exp(1.5*t)*sin(2*x)
        + 15.0*exp(2.0*t)*sin(2*x)*sin(2*x) + 0.4*exp(2.5*t)*sin(2*x)*sin(2*x)*sin(2*x))*sin(2*x)*sin(2*x) 
        + 196.2*exp(0.5*t)*sin(2*x) + 4.905*exp(1.0*t)*sin(2*x)*sin(2*x);
}

double massRHS(double x, double t, double a, double b, double c, double d, double e)
{
    return b*exp(c*t)*exp(e*t)*cos(b*x)*cos(d*x) + c*exp(c*t)*sin(b*x) - d*(a + exp(c*t)*sin(b*x))*exp(e*t)*sin(d*x);
}

double momeRHS(double x, double t, double a, double b, double c, double d, double e, double g)
{
    return b*c*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*sin(d*x)*cos(b*x) 
            + 2*b*c*d*(a + exp(c*t)*sin(b*x))*exp(2*c*t)*exp(e*t)*sin(b*x)*sin(d*x)*cos(b*x) 
            - 2*b*d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(2*e*t)*sin(d*x)*sin(d*x)*cos(b*x)
            + b*d*e*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*sin(d*x)*cos(b*x) 
            - b*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*(d*d*exp(2*e*t)*sin(d*x)*sin(d*x) 
            + d*d*exp(2*e*t)*cos(d*x)*cos(d*x) + d*e*exp(e*t)*sin(d*x))*exp(c*t)*cos(b*x) 
            - b*exp(c*t)*exp(2*e*t)*cos(b*x)*cos(d*x)*cos(d*x) 
            + c*d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*sin(b*x)*cos(d*x) 
            - 4*d*d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(2*e*t)*sin(d*x)*cos(d*x)/3 
            + 2*d*(a + exp(c*t)*sin(b*x))*exp(2*e*t)*sin(d*x)*cos(d*x) 
            - d*(b*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*sin(d*x)*cos(b*x) 
            + d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(e*t)*cos(d*x)/3 
            + (a + exp(c*t)*sin(b*x))*exp(e*t)*cos(d*x))*exp(e*t)*sin(d*x)
            + (-b*b*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*sin(b*x)*sin(d*x) 
            + 2*b*b*d*(a + exp(c*t)*sin(b*x))*exp(2*c*t)*exp(e*t)*sin(d*x)*cos(b*x)*cos(b*x) 
            + 2*b*d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(c*t)*exp(e*t)*cos(b*x)*cos(d*x) 
            + b*exp(c*t)*exp(e*t)*cos(b*x)*cos(d*x) - d*d*d*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*(a + exp(c*t)*sin(b*x))*exp(e*t)*sin(d*x)/3 
            - d*(a + exp(c*t)*sin(b*x))*exp(e*t)*sin(d*x))*exp(e*t)*cos(d*x);

}

double hxt(double x, double t)
{
    return  exp(0.5*t)*sin(2*x) + 10.0;

}

double Gxt(double x, double t)
{
    return  4*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*exp(0.5*t)*cos(2*x)/3 + 4*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*exp(1.0*t)*sin(2*x)*cos(2*x) + (exp(0.5*t)*sin(2*x) + 10.0)*exp(0.5*t)*cos(2*x);

}

double Fh(double x, double t)
{
    return (exp(0.5*t)*sin(2*x) + 10.0)*exp(0.5*t)*cos(2*x);
}

double FG(double x, double t)
{
    return -8*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*exp(1.0*t)*sin(2*x)*exp(1.0*t)*sin(2*x)/3 
            + 4.905*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0) 
            + (4*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*exp(0.5*t)*cos(2*x)/3 
            + 4*(exp(0.5*t)*sin(2*x) + 10.0)*(exp(0.5*t)*sin(2*x) + 10.0)*exp(1.0*t)*sin(2*x)*cos(2*x) 
            + (exp(0.5*t)*sin(2*x) + 10.0)*exp(0.5*t)*cos(2*x))*exp(0.5*t)*cos(2*x);
}

//include BCs
void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh, double *x, double t)
{
    double idx = 1.0 / dx;  
	int i;
    double her,Ger,uer,duer,hel,Gel,uel,duel,fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG;
	double uai,ubi,uaip1,ubip1;
    double hbarn1,hbarn, hfiph,hfimh,Gbarn1,Gbarn, Gfiph,Gfimh;

    // i = -1
    i = -1;

    //Figure out how the interior numberings work.
    

    uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
    ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

    uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
    ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

    her = hhbc[3*(i+1) +nGhBC];
    Ger = Ghbc[3*(i+1) +nGhBC];
    uer  = ubc[2*i + unBC + 1];
    duer = -uaip1*(dx) + ubip1; //(2*0.5*dx)

    hel = hhbc[3*(i) +nGhBC +2];
    Gel = Ghbc[3*(i) +nGhBC +2];
    uel  = ubc[2*i + unBC + 1];
    duel = uai*(dx) + ubi; //special formula

    fhel = uel*hel;
    fher = uer*her;

    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

    sqrtghel = sqrt(g* hel);
    sqrtgher = sqrt(g* her);

    sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

    isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);	

    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));

    fih = foh;
    fiG = foG;
    for(i = 0;i < n;i++)
    {

        uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
        ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

        uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
        ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

        her = hhbc[3*(i+1) +nGhBC];
	    Ger = Ghbc[3*(i+1) +nGhBC];
	    uer  = ubc[2*i + unBC + 1];
	    duer = -uaip1*(dx) + ubip1; //(2*0.5*dx)

	    hel = hhbc[3*(i) +nGhBC +2];
	    Gel = Ghbc[3*(i) +nGhBC +2];
	    uel  = ubc[2*i + unBC + 1];
	    duel = uai*(dx) + ubi; //special formula


	    fhel = uel*hel;
	    fher = uer*her;

	
	    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
	    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

        //printf("%d | %e | %e \n",i,fGel,fGer);

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);	

        //hbarn1 = hxt(x[i], t + dt);
        //hbarn = hxt(x[i], t );
       // hfiph = Fh(x[i] + 0.5*dx, t + 0.5*dt);
        //hfimh = Fh(x[i] - 0.5*dx, t + 0.5*dt);


        //Gbarn1 = Gxt(x[i], t + dt);
        //Gbarn = Gxt(x[i], t );
        //Gfiph = FG(x[i] + 0.5*dx, t + 0.5*dt);
        //Gfimh = FG(x[i] - 0.5*dx, t + 0.5*dt);


	    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));

        //newh[i] = hbarn1 + (hhbc[3*(i) +nGhBC + 1] - hbarn) - dt*idx*((foh - hfiph) - (fih - hfimh));
        //newG[i] = Gbarn1 + (Ghbc[3*(i) +nGhBC + 1] - Gbarn) - dt*idx*((foG - Gfiph) - (fiG - Gfimh));

        newh[i] = (hhbc[3*(i) +nGhBC + 1]) - dt*idx*((foh ) - (fih));
        newG[i] =(Ghbc[3*(i) +nGhBC + 1]) - dt*idx*((foG ) - (fiG ));

        //printf("AVG %d | h bar err : %1.15f | G bar err : %1.15f \n ",i,hhbc[3*(i) +nGhBC + 1] - hbarn , Ghbc[3*(i) +nGhBC + 1] - Gbarn);
        //printf("FLUX left %d | h f err : %1.15f | G f err : %1.15f \n ",i, (foh - hfiph) , (foG - Gfiph) );
        //printf("FLUX right %d | h f err : %1.15f | G f err : %1.15f \n \n ",i, (fih - hfimh), (fiG - Gfimh) );


	    fih = foh;
	    fiG = foG;

    }

}

double dmaxarray( double *q, int n)
{
    int i;
    double max = fabs(q[0]);
    for (i = 1 ; i < n ; i++)
    {
        if(fabs(q[i]) > max) max = fabs(q[i]);
    }

    return max;

}



void forcing(double *G, double *h, double *x, double t, double dt, double dx, int n, double g)
{

    int i = 0;
    double hbarn1, hbarn,hfiph,hfimh,Gbarn1,Gbarn,Gfiph,Gfimh;

    for (i = 0 ; i < n ; i++)
    {

        hbarn1 = (barhxt(x[i] + 0.5*dx,t + dt) - barhxt(x[i] - 0.5*dx,t + dt) ) /dx;
        hbarn = (barhxt(x[i] + 0.5*dx,t) - barhxt(x[i] - 0.5*dx,t) ) /dx;
        hfiph = (Fluxh(x[i] + 0.5*dx,t+dt) - Fluxh(x[i] + 0.5*dx,t)) /dt;
        hfimh = (Fluxh(x[i] - 0.5*dx,t+dt) - Fluxh(x[i] - 0.5*dx,t)) /dt;


        Gbarn1 = (barGxt(x[i] + 0.5*dx,t + dt) - barGxt(x[i] - 0.5*dx,t + dt) ) /dx;
        Gbarn = (barGxt(x[i] + 0.5*dx,t) - barGxt(x[i] - 0.5*dx,t) ) /dx;
        Gfiph = (FluxG(x[i] + 0.5*dx,t+dt,g) - FluxG(x[i] + 0.5*dx,t,g)) /dt;
        Gfimh = (FluxG(x[i] - 0.5*dx,t+dt,g) - FluxG(x[i] - 0.5*dx,t,g)) /dt;


        h[i] = h[i] + (hbarn1-  hbarn + dt/dx*(hfiph - hfimh) );

        G[i] = G[i] + (Gbarn1-  Gbarn + dt/dx*(Gfiph - Gfimh) );
    }


}

void evolvewrapperconsistenttimeForced(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend,double *hMbeg1 , double *hMend1,double *GMbeg1 , double *GMend1,double *uMbeg1, double *uMend1, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp, double *x, double t, double a, double b, double c, double d, double e)
{

    int u_length = 2*n + 1;

    //nBCs is number of cells to define ghost cells
    getufromG(h,G,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);
    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,dt,n, theta ,Gp, hp,x,t);

    forcing(Gp,hp, x, t, dt, dx,n,g);

    getufromG(hp,Gp,hMbeg1,hMend1,GMbeg1,GMend1,uMbeg1,uMend1,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);
    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,dt,n, theta ,Gpp, hpp,x,t + dt);

    forcing(Gpp,hpp, x, t + dt, dt, dx,n,g);

    int i;
    for(i=0;i<n;i++)
    {
        h[i] = 0.5*(h[i] + hpp[i]);
        G[i] = 0.5*(G[i] + Gpp[i]);
        
    }

    //forcing(G,h, x, t, dt, dx,n,g);

}

void evolvewrapperconsistenttime(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp)
{

    int u_length = 2*n + 1;

    //nBCs is number of cells to define ghost cells
    getufromG(h,G,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    double Mass, Mome, Ener;

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);


    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,dt,n, theta ,Gp, hp,h,0);
//void evolve(double *Ghbc, double *hhbc, double *ubc, int nGhhbc, int nubc, int nGhBC, int unBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh)
    getufromG(hp,Gp,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);


    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,dt,n, theta ,Gpp, hpp,h,0);

    int i;
    for(i=0;i<n;i++)
    {
        h[i] = 0.5*(h[i] + hpp[i]);
        G[i] = 0.5*(G[i] + Gpp[i]);
        
    }


    getufromG(h,G,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);

}
/*
double evolvewrapperADAP(double *G, double *h,double *hMbeg , double *hMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double g, double dx, double dt, int n, int nGhBC, int unBC, int nGhhbc, int nubc, double theta, double *hhbc,double *Ghbc,double *ubc, double *Gp, double *hp, double *Gpp, double *hpp)
{

    int u_length = 2*n + 1;
    int i;



    //nBCs is number of cells to define ghost cells
    getufromG(h,G,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    double Mass, Mome, Ener, Energ1, Energ2, Energ3;

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);

    Energ1 = HamilKEVW(hhbc,ubc,n,dx);
    Energ2 = HamilKEHW(hhbc,ubc,n,dx);
    Energ3 = HamilGPEW(hhbc,ubc,n,dx);

    printf("1 : %f | %f | %f (%f, %f, %f) \n", Mass, Mome, Ener,Energ1,Energ2,Energ3);


    double hmax,umax;
    hmax = dmaxarray(hhbc,nGhhbc);
    umax = dmaxarray(ubc, nubc);

    double ndt;

    ndt = 0.5/ (umax + sqrt(g*hmax)) * dx;

    //printf("%f | %f | %f | %f \n ",hmax,umax,dt,ndt);

    if(ndt > dt) ndt = dt;



    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,ndt,n, theta ,Gp, hp);
    getufromG(hp,Gp,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);



    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);

    Energ1 = HamilKEVW(hhbc,ubc,n,dx);
    Energ2 = HamilKEHW(hhbc,ubc,n,dx);
    Energ3 = HamilGPEW(hhbc,ubc,n,dx);

    printf("2 : %f | %f | %f (%f, %f, %f) \n", Mass, Mome, Ener,Energ1,Energ2,Energ3);

    evolve(Ghbc, hhbc, ubc,nGhhbc,nubc,nGhBC,unBC,g, dx,ndt,n, theta ,Gpp, hpp);

    getufromG(hpp,Gpp,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);

    Energ1 = HamilKEVW(hhbc,ubc,n,dx);
    Energ2 = HamilKEHW(hhbc,ubc,n,dx);
    Energ3 = HamilGPEW(hhbc,ubc,n,dx);

    printf("3 : %f | %f | %f (%f, %f, %f) \n", Mass, Mome, Ener,Energ1,Energ2,Energ3);

    for(i=0;i<n;i++)
    {
        h[i] = 0.5*(h[i] + hpp[i]);
        G[i] = 0.5*(G[i] + Gpp[i]);
        
    }

    getufromG(h,G,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,theta,dx ,n,u_length ,nGhBC,unBC,nGhhbc,nubc,ubc,hhbc,Ghbc);

    Ener = HamilW(hhbc,ubc,n,dx);
    Mome = uhALLW(hhbc,ubc,n,dx);
    Mass = hALLW(hhbc,n,dx);

    Energ1 = HamilKEVW(hhbc,ubc,n,dx);
    Energ2 = HamilKEHW(hhbc,ubc,n,dx);
    Energ3 = HamilGPEW(hhbc,ubc,n,dx);

    printf("F : %f | %f | %f (%f, %f, %f) \n", Mass, Mome, Ener,Energ1,Energ2,Energ3);

    return ndt;

}
*/




int main()
{
    printf("h");
    return 1;
}
