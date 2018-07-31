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

const double E = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069;


#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-60
#define div_0 1.0e-15

#define htol 1.0e-15


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

//Conservation

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

double Gacrosscell(double *x,double *G,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *Gcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(G,Gcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(Gcoeff,x[j],fgp);
    
    double fgpe = fgph;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(Gcoeff,x[j],sgp);    
    double sgpe = sgph;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(Gcoeff,x[j],tgp);
    
    double tgpe = tgph;

	free(Gcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}

double hacrosscell(double *x,double *h,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    
    double fgpe = fgph;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);    
    double sgpe = sgph;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    
    double tgpe = tgph;

	free(hcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
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
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}


    
double HankEnergyacrosscell(double *x,double *h,double *u, double *b,double g,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *ucoeff = malloc(5*sizeof(double));
	double *hcoeff = malloc(5*sizeof(double));
	double *bcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(u,ucoeff,j,dx);
    interpquartcoeff(h,hcoeff,j,dx);
    interpquartcoeff(b,bcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    double fgpux = interpquarticgrad(ucoeff,x[j],fgp);
	double fgpb = interpquarticval(bcoeff,x[j],fgp);
	double fgpbx = interpquarticval(bcoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu*fgpu + 2*g*fgph*(0.5*fgph + fgpb) + i12*(fgph*fgph*fgph)*fgpux*fgpux 
				+ fgph*(fgpu*fgpbx - 0.5*fgpux*fgph)*(fgpu*fgpbx - 0.5*fgpux*fgph);
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    double sgpux = interpquarticgrad(ucoeff,x[j],sgp);
	double sgpb = interpquarticval(bcoeff,x[j],sgp);
	double sgpbx = interpquarticval(bcoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu*sgpu + 2*g*sgph*(0.5*sgph + sgpb) + i12*(sgph*sgph*sgph)*sgpux*sgpux 
				+ sgph*(sgpu*sgpbx - 0.5*sgpux*sgph)*(sgpu*sgpbx - 0.5*sgpux*sgph);

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
	double tgpb = interpquarticval(bcoeff,x[j],tgp);
	double tgpbx = interpquarticval(bcoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + 2*g*tgph*(0.5*tgph + tgpb) + i12*(tgph*tgph*tgph)*tgpux*tgpux 
				+ tgph*(tgpu*tgpbx - 0.5*tgpux*tgph)*(tgpu*tgpbx - 0.5*tgpux*tgph);

	free(ucoeff);
	free(hcoeff);
	free(bcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
}
    
double HankEnergyall(double *x,double *h,double *u, double *b,double g,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + HankEnergyacrosscell(x,h,u,b,g,i,dx);
		//printf("i : %d || x : %f || h : %f || u : %f \n",i,x[i],h[i],u[i]);
	}
    return 0.5*sum1; 

}

double uhall(double *x,double *h,double *u,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + uhacrosscell(x,h,u,i,dx);
		//printf("i : %d || x : %f || h : %f || u : %f \n",i,x[i],h[i],u[i]);
	}
    return sum1; 

}

double hall(double *x,double *h,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + hacrosscell(x,h,i,dx);
		//printf("i : %d || x : %f || h : %f || u : %f \n",i,x[i],h[i],u[i]);
	}
    return sum1; 

}


double Gall(double *x,double *G,int n, int nBC,double dx)
{
    double sum1 = 0.0;
	int i;
	for(i = nBC; i < n - nBC;i++)
	{
       sum1 = sum1 + hacrosscell(x,G,i,dx);
		//printf("i : %d || x : %f || h : %f || u : %f \n",i,x[i],h[i],u[i]);
	}
    return sum1; 

}


//GENERAL STUFF
void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    int i;
    for(i = 0; i < n;i++)
	{
       d[i] = a[i];
	}
    for(i = n; i < n + m;i++)
	{
       d[i] = b[i-n];
	}
    for(i = n + m; i < n + m + k;i++)
	{
       d[i] = c[i-n - m];
	}
    //memcpy(d,a,n*sizeof(double));
    //memcpy(d+n,b,m*sizeof(double));
    //memcpy(d+n+m,c,k*sizeof(double));
}


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

void deallocPy(void *x)
{
    free(x);
}

int readfrom2DmemINT(int *x,int i,int j, int m)
{
    return x[i*m + j];
}
//HMMM memory




double ht(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
    return a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a5 + (a2*(-a3 - a2*t + x))/a4);
}

double Gt(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
    return ((6*a1*a2*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/pow(a4,2) + 
      (2*a2*a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/pow(a4,2) + 
      (3*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(a5 + (a2*(-a3 - a2*t + x))/a4))/a4 - 
      (6*pow(a1,2)*a6*pow(E,2*a5*t + a7*t - (3*pow(-a3 - a2*t + x,2))/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,2)*(a5 + (a2*(-a3 - a2*t + x))/a4))/pow(a4,2) - 
      (3*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2)*(a5 + (a2*(-a3 - a2*t + x))/a4))/pow(a4,2) + 
      (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(a7 + (a2*(-a3 - a2*t + x))/a4))/a4 - 
      (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,2)*(a7 + (a2*(-a3 - a2*t + x))/a4))/pow(a4,2) - 
      (3*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2)*(a5 + a7 + (2*a2*(-a3 - a2*t + x))/a4))/pow(a4,2))/3. + 
   a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*(a5 + (a2*(-a3 - a2*t + x))/a4)*
    (1 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*cos(a9*x))/a4 + pow(a8,2)*pow(a9,2)*pow(cos(a9*x),2) - 
      (a8*pow(a9,2)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a9*x))/2.) + 
   a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*(a7 + (a2*(-a3 - a2*t + x))/a4)*
    (1 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*cos(a9*x))/a4 + pow(a8,2)*pow(a9,2)*pow(cos(a9*x),2) - 
      (a8*pow(a9,2)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a9*x))/2.) + 
   a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*
    ((a1*a2*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*cos(a9*x))/a4 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*(a5 + (a2*(-a3 - a2*t + x))/a4)*cos(a9*x))/a4 - 
      (a1*a8*pow(a9,2)*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a5 + (a2*(-a3 - a2*t + x))/a4)*sin(a9*x))/2.);
}

double Fluxhdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6,double a7,double a8, double a9, double g)
{
    return -((a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*(-a3 - a2*t + x))/a4) - (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x))/a4;
}

double FluxGdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
    return (-4*pow(a6,2)*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(3.*pow(a4,2)) - 
   (a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*g*(-a3 - a2*t + x))/a4 + 
   (2*a1*pow(a6,2)*pow(E,a5*t + 2*a7*t - (3*pow(-a3 - a2*t + x,2))/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/pow(a4,3) + 
   (4*pow(a6,2)*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(3.*pow(a4,3)) - 
   (pow(a6,2)*a8*a9*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*cos(a9*x))/a4 + 
   (2*a1*pow(a6,2)*a8*a9*pow(E,a5*t + 2*a7*t - (3*pow(-a3 - a2*t + x,2))/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,2)*cos(a9*x))/pow(a4,2) + 
   (2*pow(a6,2)*a8*a9*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2)*cos(a9*x))/pow(a4,2) + 
   (pow(a6,2)*a8*pow(a9,2)*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x)*sin(a9*x))/a4 - 
   (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*(((a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3))/a4 - 
           (3*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2))/pow(a4,2) - 
           (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,2))/pow(a4,2))/3. + 
        a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*
         (1 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*cos(a9*x))/a4 + pow(a8,2)*pow(a9,2)*pow(cos(a9*x),2) - 
           (a8*pow(a9,2)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a9*x))/2.)))/a4 + 
   a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(((-9*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/pow(a4,2) - 
         (3*a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/pow(a4,2) + 
         (6*pow(a1,2)*a6*pow(E,2*a5*t + a7*t - (3*pow(-a3 - a2*t + x,2))/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,3))/pow(a4,3) + 
         (9*a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/pow(a4,3) + 
         (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/pow(a4,3))/3. - 
      (a1*a6*pow(E,a5*t + a7*t - pow(-a3 - a2*t + x,2)/a4)*(-a3 - a2*t + x)*(1 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*cos(a9*x))/a4 + pow(a8,2)*pow(a9,2)*pow(cos(a9*x),2) - 
           (a8*pow(a9,2)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a9*x))/2.))/a4 - 
      (a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x)*
         (1 - (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*cos(a9*x))/a4 + pow(a8,2)*pow(a9,2)*pow(cos(a9*x),2) - 
           (a8*pow(a9,2)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a9*x))/2.))/a4 + 
      a6*pow(E,a7*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*
       (-((a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*cos(a9*x))/a4) - (a8*pow(a9,3)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*cos(a9*x))/2. + 
         (a1*a8*a9*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*pow(-a3 - a2*t + x,2)*cos(a9*x))/pow(a4,2) + (3*a1*a8*pow(a9,2)*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4))*(-a3 - a2*t + x)*sin(a9*x))/(2.*a4) - 
         2*pow(a8,2)*pow(a9,3)*cos(a9*x)*sin(a9*x)));
}

double SourceG(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
    return a8*a9*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*g*cos(a9*x) + (pow(a6,2)*a8*pow(a9,2)*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*pow(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)),2)*
      (-a3 - a2*t + x)*sin(a9*x))/(2.*a4) + pow(a6,2)*pow(a8,2)*pow(a9,3)*pow(E,2*a7*t - pow(-a3 - a2*t + x,2)/a4)*(a0 + a1*pow(E,a5*t - pow(-a3 - a2*t + x,2)/(2.*a4)))*cos(a9*x)*sin(a9*x);
}


double hFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double phi  = x - a2*t  ;
      return a0 + a1*pow(E, (-pow(phi - a3,2)/(2*a4)))*pow(E,a5*t);
}

double hxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double phi  = x - a2*t  ;
      return -a1/a4*(phi - a3)*pow(E, (-pow(phi - a3,2)/(2*a4)))*pow(E,a5*t);
}

double uFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double phi  = x - a2*t  ;
      return a6*pow(E, (-pow(phi - a3,2)/(2*a4)))*pow(E,a7*t);
}

double uxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double phi  = x - a2*t  ;
      return -a6/a4*(phi - a3)*pow(E, (-pow(phi - a3,2)/(2*a4)))*pow(E,a7*t);
}

double uxxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double phi  = x - a2*t  ;
      return -a6/(a4*a4)*(a4 - pow(phi - a3,2))*pow(E, (-pow(phi - a3,2)/(2*a4)))*pow(E,a7*t);
}

double bFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      return a8*sin(a9*x);
}

double bxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      return a8*a9*cos(a9*x);
}

double bxxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      return -a8*a9*a9*sin(a9*x);
}

double wFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
      double h = hFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
      double b = bFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
      return h + b;
}

double GFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9, double g)
{
     double u = uFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
     double h = hFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
     double hx = hxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
     double bx = bxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
     double bxx = bxxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    double ux = uxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    double uxx = uxxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

     return u*h*(1 + hx*bx + 0.5*h*bxx + bx*bx) - h*h*hx*ux - h*h*h/3.0*uxx;
}


//include BCs
void evolveForce(double g, double dx, double dt, int n, double *newG, double *newh, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9)
{
    double idx = 1.0 / dx;  
	int i;
    double her,Ger,dber,uer,duer,hel,Gel,dbel,uel,duel,fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG,th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
	double wil,wir,wip1l,bip1l,bil,bir,nbi,hihm,hihp,uai,ubi,uaip1,ubip1;
    double himhp,bedai,bedbi,bedci,bedaip1,bedbip1,bedcip1;

    double hS, GS;



    // i = -1
    i = -1;

    her = hFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    Ger = GFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    dber = bxFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    uer  = uFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    duer = uxFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

    hel = hFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    Gel = GFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    dbel = bxFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    uel  = uFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
    duel = uxFunc(x[0] - 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

    fhel = uel*hel;
    fher = uer*her;

    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

    sqrtghel = sqrt(g* hel);
    sqrtgher = sqrt(g* her);

    sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

    isrmsl = 0.0;

    if(sr - sl > TINY) isrmsl = 1.0 / (sr - sl);	

    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));

    fih = foh;
    fiG = foG;

    for(i = 0;i < n;i++)
    {

        her = hFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        Ger = GFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        dber = bxFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        uer  = uFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        duer = uxFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

        hel = hFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        Gel = GFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        dbel = bxFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        uel  = uFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
        duel = uxFunc(x[i] + 0.5*dx,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

	    fhel = uel*hel;
	    fher = uer*her;

	    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
	    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr - sl > TINY) isrmsl = 1.0 / (sr - sl);	

	    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));


        //centerted values
        th = hFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
		tu = uFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
		tux = uxFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
		tbx = bxFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
		tbxx = bxxFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);
		
		sourcec = -g*th*tbx -  0.5*th*th*tu*tux*tbxx + th*tu*tu*tbx*tbxx ;

        hS = ht(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

        GS = Gt(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

        //printf("hchange : %d | %e \n ",i,- dt*idx*(foh - fih));
        //printf("Gchange : %d | %e \n ",i,- dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + dx*sourcec));


        /*if(isnan(foh))
        {
            printf("Right Values %d : %e | %e | %e | %e | %e \n",i,her,Ger,dber,uer,duer);
            printf("Right Fluxes %d : %e | %e | %e | %e | %e \n",i,fher,fGer,sqrtgher,sr,isrmsl);
            printf("Left  Values %d : %e | %e | %e | %e | %e \n",i,hel,Gel,dbel,uel,duel);
             printf("Left Fluxes %d : %e | %e | %e | %e | %e \n",i,fhel,fGel,sqrtghel,sl,isrmsl);
        } */
		newh[i] = hFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) - dt*idx*(foh - fih) + dt*hS;
		newG[i] = GFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g)- dt*idx*(foG - fiG) + dt*idx*(dx*sourcec)+ dt*GS;

        fih = foh;
        fiG = foG;

    }

}



void evolvewrapForcingANA(double *h, double *G,int n,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9)
{
    double idx = 1.0 /dx;
    //n is number of cells
    int u_length = 2*n + 1;
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

    //Solve for u given wet/dry regions
    evolveForce(g, dx,dt,n,Gp, hp,x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);

    //Solve for u given wet/dry regions
    evolveForce(g, dx,dt,n,Gpp, hpp,x,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9); 


    int i;
    for(i=0;i<n;i++)
    {
        G[i] =0.5*(GFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + Gpp[i]);
        h[i] =0.5*(hFunc(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + hpp[i]);
        
    }

    //Reg2 = RegSplit(h,n);
    //Solve for u given wet/dry regions
    //getufromGsplit(h,G,bed,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,wMbeg,wMend,bMbeg,bMend,theta,dx ,n,u_length ,nGhBC,unBC,bnBC,nGhhbc,nubc,nbhc,ubc,hhbc,Ghbc,whbc,bedhbc,Reg2);


    free(hp);
    free(Gp);
    free(Gpp);
    free(hpp);
}

int main()
{
    printf("h");
    return 1;
}
