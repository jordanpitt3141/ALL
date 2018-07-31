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
#define TINY 1.0e-20
#define div_0 1.0e-15

#define hbase 1.0e-3
#define htol 1.0e-6


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

void deallocPy(void *x)
{
    free(x);
}

int readfrom2DmemINT(int *x,int i,int j, int m)
{
    return x[i*m + j];
}
//HMMM memory

int *append(int *A, int n, int m, int *B)
{
    int j;
    int *NA;

    if(n == 0)
    {
       //NA = mallocFlat2DInt(1, 5 );

       NA = (int*) malloc((1*5)*sizeof(int)); 
    }
    else
    {
       NA = (int*) realloc(A,(n + 1)*m*sizeof(int));  
    }

    for (j=0;j < m;j++)
    {
        //printf("%d  | %d | %d \n",j,n,m);
        NA[n*m +j] = B[j];
    }
    return NA;

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


int *RegSplit(double *h,int n)
{
    int arr[5];
    int ini;
    int Csi,Cregval,regvali,Cei,i,j;
    int *Reg,*Reg1;

    ini = 0;    
    Csi = 0;
    if (h[0] < htol)
    {
        Cregval = 0;
    }
    else
    {
        Cregval = 1;
    }

    for (i = 0; i<n;i++)
    {
        if (h[i] < htol)
        {
            regvali = 0;
        }
        else
        {
            regvali = 1;
        }
        if(regvali != Cregval)
        {
            Cei = i-1;
            arr[0] = Csi;
            arr[1] = Cei;
            arr[2] = (Cei - Csi);
            arr[3] = Cregval;

            Reg1 = append(Reg, ini, 5, arr);

            Csi = i ;
            Cregval =regvali ;
            ini = ini + 1;

            Reg = (int*) malloc((ini*5)*sizeof(int)); 
            memcpy(Reg, Reg1, (ini*5)*sizeof(int));

            for(j =0;j < ini;j++)
            {
               Reg[j*5 +4] = ini;
            }

            free(Reg1);

        }
        else if(regvali == Cregval && i==n-1 )
        {
            Cei = i;
            arr[0] = Csi;
            arr[1] = Cei;
            arr[2] = (Cei - Csi);
            arr[3] = Cregval;

            Reg1 = append(Reg, ini, 5, arr);

            ini = ini + 1;

            Reg = (int*) malloc((ini*5)*sizeof(int)); 
            memcpy(Reg, Reg1, (ini*5)*sizeof(int));
            

            for(j =0;j < ini;j++)
            {
               Reg[j*5 +4] = ini;
            }

            free(Reg1);
            
        }
    }
    
    return Reg;
}




// ####################################################################################################### END OF  CODE REQUIRED TO INTERFACE WITH PYTHON ##############################################################


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

void getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc, double *bedhbc)
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

    double wi,wip1,wim1,dwib,dwim,dwif,dwi,wjphm,wjmhp, bedbegmiddle,bedai,bedbi,bedci,beddi,bedendmiddle,bjmh,bjms,bjps,bjph;
    double h2bxuvxa11,h2bxuvxa12,h2bxuvxa13,h2bxuvxa21,h2bxuvxa22,h2bxuvxa23,h2bxuvxa31,h2bxuvxa32,h2bxuvxa33;
    double h2bxuxva11,h2bxuxva12,h2bxuxva13,h2bxuxva21,h2bxuxva22,h2bxuxva23,h2bxuxva31,h2bxuxva32,h2bxuxva33;
    double hbx2uva11,hbx2uva12,hbx2uva13,hbx2uva21,hbx2uva22,hbx2uva23,hbx2uva31,hbx2uva32,hbx2uva33;

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



    // first 
    i = 0;
    j = 1;

    // Get it from interior
    // Reconstruct w

    // reconstruct bed
    bedai = i6*idx*idx*idx*(bed[i+3] -3*bed[i+2] +3*bed[i+1] - bed[i]);
    bedbi = 0.5*idx*idx*(-bed[i+3] + 4*bed[i+2] - 5*bed[i+1] + 2*bed[i]);
    bedci = i6*idx*(2*bed[i+3] - 9*bed[i+2] + 18*bed[i+1] - 11*bed[i]);
    beddi = bed[i];

    bjmh = -bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) - bedci*(0.5*dx) + beddi;
    bjms = -bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) - bedci*(dx*i6) + beddi;
    bjps = bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) + bedci*(dx*i6) + beddi;
    bjph = bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) + bedci*(0.5*dx) + beddi;

    dGib = (G[i] - GMbeg[nGhBC-2]);
    dGim = 0.5*(G[i+1] - GMbeg[nGhBC-2]);
    dGif = (G[i+1] - G[i]);

    dhib = (h[i] - hMbeg[nGhBC-2]);
    dhim = 0.5*(h[i+1] - hMbeg[nGhBC-2]);
    dhif = (h[i+1] - h[i]);

    wi = h[i] + bed[i];
    wip1 = h[i+1] + bed[i+1];
    wim1 = hMbeg[nGhBC-2] + (bedai*(-dx)*(-dx)*(-dx) + bedbi*(-dx)*(-dx) + bedci*(-dx) + beddi);
    dwib = (wi - wim1);
    dwim = 0.5*(wip1 - wim1);
    dwif = (wip1 - wi);

    dGi = minmod(theta*dGib, dGim, theta*dGif);
    dhi = minmod(theta*dhib, dhim, theta*dhif);
    dwi = minmod(theta*dwib, dwim, theta*dwif);

    wjphm= wi + 0.5*dwi;
    wjmhp= wi - 0.5*dwi; 

    hjphm=  h[i] + 0.5*dhi; 
    hjmhp=  h[i] - 0.5*dhi;  

    Gjphm= G[i] + 0.5*dGi;
    Gjmhp= G[i] - 0.5*dGi;

    hhbc[3*i] = hjmhp;
    hhbc[3*i+1] = h[i];
    hhbc[3*i+2] = hjphm; 


    Ghbc[3*i] = Gjmhp;
    Ghbc[3*i+1] = G[i];
    Ghbc[3*i+2] = Gjphm;

    whbc[3*i] = wjmhp;
    whbc[3*i+1] =  wi;
    whbc[3*i+2] = wjphm; 

    bedhbc[4*i - 1] = 0.5*(bjmh + bedhbc[4*i - 1]);
    bedhbc[4*i] = bedhbc[4*i - 1];
    bedhbc[4*i+1] = bjms;
    bedhbc[4*i+2] = bjps;
    bedhbc[4*i+3] = bjph; 

    hjphm= hjphm + (hbase)/ (hjphm + htol); 
    hjmhp= hjmhp + (hbase)/ (hjmhp + htol);



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

    //h2bxuvx
    h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

    h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

    h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

    //h2bxuxv
     h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

    h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

    h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

    h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

    h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

    h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

    h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

    h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

    h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

    //hbx2uv

    hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph);

    hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

    hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

    hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

    hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph);

    hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

    hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

    hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

    hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph);
    
    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuvxa11 + h2bxuxva11 + hbx2uva11;
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuvxa12 + h2bxuxva12 + hbx2uva12; 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuvxa13 + h2bxuxva13 + hbx2uva13;
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuvxa21 + h2bxuxva21 + hbx2uva21;
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuvxa22 + h2bxuxva22 + hbx2uva22; 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuvxa23 + h2bxuxva23 + hbx2uva23;
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuvxa31 + h2bxuxva31 + hbx2uva31;
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuvxa32 + h2bxuxva32 + hbx2uva32;
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuvxa33 + h2bxuxva33 + hbx2uva33;


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

        wi = h[i] + bed[i];
        wip1 = h[i+1] + bed[i+1];
        wim1 = h[i-1] + bed[i-1];
        dwib = (wi - wim1);
        dwim = 0.5*(wip1 - wim1);
        dwif = (wip1 - wi);

        dGi = minmod(theta*dGib, dGim, theta*dGif);
        dhi = minmod(theta*dhib, dhim, theta*dhif);
        dwi = minmod(theta*dwib, dwim, theta*dwif);

        wjphm= wi + 0.5*dwi;
        wjmhp= wi - 0.5*dwi; 

        hjphm=  h[i] + 0.5*dhi; 
        hjmhp=  h[i] - 0.5*dhi;  

        Gjphm= G[i] + 0.5*dGi;
        Gjmhp= G[i] - 0.5*dGi;  


        // reconstruct bed
        if(i == 1) 
        {
            bedbegmiddle = -bMbeg[0]*i16 + 9*bMbeg[1]*i16 + 9*bMbeg[2]*i16 - bMbeg[3]*i16;

            bedai = i12*idx*idx*idx*(-bedbegmiddle + 2*bed[i-1] -2*bed[i+1] + bed[i+2]);
            bedbi = i6*idx*idx*(bedbegmiddle - bed[i-1] - bed[i+1] + bed[i+2]);
            bedci = i12*idx*(bedbegmiddle - 8*bed[i-1] + 8*bed[i+1] - bed[i+2]);
            beddi = -bedbegmiddle*i6 + 2*i3*bed[i-1] + 2*i3*bed[i+1] - bed[i+2]*i6;
        }
        else if(i == n-2)
        {
            bedendmiddle = -bMend[0]*i16 + 9*bMend[1]*i16 + 9*bMend[2]*i16 - bMend[3]*i16;

            bedai = i12*idx*idx*idx*(-bed[i-2] + 2*bed[i-1] -2*bed[i+1] + bedendmiddle);
            bedbi = i6*idx*idx*(bed[i-2] - bed[i-1] - bed[i+1] + bedendmiddle);
            bedci = i12*idx*(bed[i-2] - 8*bed[i-1] + 8*bed[i+1] - bedendmiddle);
            beddi = -bed[i-2]*i6 + 2*i3*bed[i-1] + 2*i3*bed[i+1] - bedendmiddle*i6;

        }
        else
        {
            bedai = i12*idx*idx*idx*(-bed[i-2] + 2*bed[i-1] -2*bed[i+1] + bed[i+2]);
            bedbi = i6*idx*idx*(bed[i-2] - bed[i-1] - bed[i+1] + bed[i+2]);
            bedci = i12*idx*(bed[i-2] - 8*bed[i-1] + 8*bed[i+1] - bed[i+2]);
            beddi = -bed[i-2]*i6 + 2*i3*bed[i-1] + 2*i3*bed[i+1] - bed[i+2]*i6;

        }

        bjmh = -bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) - bedci*(0.5*dx) + beddi;
        bjms = -bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) - bedci*(dx*i6) + beddi;
        bjps = bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) + bedci*(dx*i6) + beddi;
        bjph = bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) + bedci*(0.5*dx) + beddi;  


        hhbc[3*i] = hjmhp;
        hhbc[3*i+1] = h[i];
        hhbc[3*i+2] = hjphm; 


        Ghbc[3*i] = Gjmhp;
        Ghbc[3*i+1] = G[i];
        Ghbc[3*i+2] = Gjphm;

        whbc[3*i] = wjmhp;
        whbc[3*i+1] =  wi;
        whbc[3*i+2] = wjphm; 

        bedhbc[4*i-1] = 0.5*(bjmh + bedhbc[4*i - 1]);
        bedhbc[4*i] = bedhbc[4*i-1];
        bedhbc[4*i+1] = bjms;
        bedhbc[4*i+2] = bjps;
        bedhbc[4*i+3] = bjph; 

    hjphm= hjphm + (hbase)/ (hjphm); 
    hjmhp= hjmhp + (hbase)/ (hjmhp);





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

        //h2bxuvx
        h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

        h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

        h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

        h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

        h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

        h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

        h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

        h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

        h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

        //h2bxuxv
         h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

        h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

        h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

        h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

        h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

        h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

        h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

        h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

        h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

        //hbx2uv

        hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph);

        hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

        hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

        hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

        hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph);

        hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

        hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

        hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

        hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph);
        
        // LHS 
        
        LHSa11 = uhintia11 + h3uxintia11 + h2bxuvxa11 + h2bxuxva11 + hbx2uva11;
        LHSa12 = uhintia12 + h3uxintia12 + h2bxuvxa12 + h2bxuxva12 + hbx2uva12; 
        LHSa13 = uhintia13 + h3uxintia13 + h2bxuvxa13 + h2bxuxva13 + hbx2uva13;
        LHSa21 = uhintia21 + h3uxintia21 + h2bxuvxa21 + h2bxuxva21 + hbx2uva21;
        LHSa22 = uhintia22 + h3uxintia22 + h2bxuvxa22 + h2bxuxva22 + hbx2uva22; 
        LHSa23 = uhintia23 + h3uxintia23 + h2bxuvxa23 + h2bxuxva23 + hbx2uva23;
        LHSa31 = uhintia31 + h3uxintia31 + h2bxuvxa31 + h2bxuxva31 + hbx2uva31;
        LHSa32 = uhintia32 + h3uxintia32 + h2bxuvxa32 + h2bxuxva32 + hbx2uva32;
        LHSa33 = uhintia33 + h3uxintia33 + h2bxuvxa33 + h2bxuxva33 + hbx2uva33;

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

      
        
        j = j + 2 ;       


    }


  

// last
    j = m-2;
    i = n-1;

   // Get it from interior

    // reconstruct bed

    bedai = i6*idx*idx*idx*(-bed[i-3] +3*bed[i-2] -3*bed[i-1] + bed[i]); 
    bedbi = 0.5*idx*idx*(-bed[i-3] + 4*bed[i-2] - 5*bed[i-1] + 2*bed[i]);
    bedci = i6*idx*(-2*bed[i-3] + 9*bed[i-2] - 18*bed[i-1] + 11*bed[i]);
    beddi = bed[i];

    bjmh = -bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) - bedci*(0.5*dx) + beddi;
    bjms = -bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) - bedci*(dx*i6) + beddi;
    bjps = bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) + bedci*(dx*i6) + beddi;
    bjph = bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) + bedci*(0.5*dx) + beddi;


    // Reconstruct w
    wi = h[i] + bed[i]; 

    dGib = (G[i] - G[i-1]);
    dGim = 0.5*(GMend[1] - G[i-1]);
    dGif = (GMend[1] - G[i]);

    dhib = (h[i] - h[i-1]);
    dhim = 0.5*(hMend[1] - h[i-1]);
    dhif = (hMend[1] - h[i]);

    wi = h[i] + bed[i];
    wip1 = hMend[1] + (bedai*(dx)*(dx)*(dx) + bedbi*(dx)*(dx) + bedci*(dx) + beddi);
    wim1 = h[i-1] + bed[i-1];
    dwib = (wi - wim1);
    dwim = 0.5*(wip1 - wim1);
    dwif = (wip1 - wi);

    dGi = minmod(theta*dGib, dGim, theta*dGif);
    dhi = minmod(theta*dhib, dhim, theta*dhif);
    dwi = minmod(theta*dwib, dwim, theta*dwif);

    wjphm= wi + 0.5*dwi;
    wjmhp= wi - 0.5*dwi; 

    hjphm=  h[i] + 0.5*dhi; 
    hjmhp=  h[i] - 0.5*dhi;  

    Gjphm= G[i] + 0.5*dGi;
    Gjmhp= G[i] - 0.5*dGi;  

    hhbc[3*i] = hjmhp;
    hhbc[3*i+1] = h[i];
    hhbc[3*i+2] = hjphm; 


    Ghbc[3*i] = Gjmhp;
    Ghbc[3*i+1] = G[i];
    Ghbc[3*i+2] = Gjphm;

    whbc[3*i] = wjmhp;
    whbc[3*i+1] =  wi;
    whbc[3*i+2] = wjphm; 

    bedhbc[4*i-1] = 0.5*(bjmh + bedhbc[4*i - 1]);
    bedhbc[4*i] = bedhbc[4*i-1];
    bedhbc[4*i+1] = bjms;
    bedhbc[4*i+2] = bjps;
    bedhbc[4*i+3] = 0.5*(bjph + bedhbc[4*i+4]);
    bedhbc[4*i+4] = bedhbc[4*i+3];


    hjphm= hjphm + (hbase)/ (hjphm+ htol); 
    hjmhp= hjmhp + (hbase)/ (hjmhp+ htol);



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

    //h2bxuvx
    h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

    h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

    h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

    h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

    h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

    //h2bxuxv
     h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph);

    h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph);

    h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph);

    h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph);

    h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph);

    h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph);

    h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph);

    h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph);

    h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph);

    //hbx2uv

    hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph);

    hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

    hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

    hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph);

    hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph);

    hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

    hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph);

    hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph);

    hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph);
    
    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuvxa11 + h2bxuxva11 + hbx2uva11;
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuvxa12 + h2bxuxva12 + hbx2uva12; 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuvxa13 + h2bxuxva13 + hbx2uva13;
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuvxa21 + h2bxuxva21 + hbx2uva21;
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuvxa22 + h2bxuxva22 + hbx2uva22; 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuvxa23 + h2bxuxva23 + hbx2uva23;
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuvxa31 + h2bxuxva31 + hbx2uva31;
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuvxa32 + h2bxuxva32 + hbx2uva32;
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuvxa33 + h2bxuxva33 + hbx2uva33;
       

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




    double d;
    d = bandec(Ared, m, 2,2, AL ,indx);

    banbks(Ared, m,2,2, AL, indx, nGis);



    memcpy(u,nGis,m*sizeof(double));

    free_dmatrix(Ared, 1, m, 1, 5);
    free_dmatrix(AL, 1, m, 1, 5);
    free(indx);
    free(nGis);

}

void getufromGsplit(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc, double *bedhbc)
{
   double idx = 1.0 / dx;
   int i,j,hbegI,bbegI;
   double bedbegmiddle,bedendmiddle, bedai, bedbi, bedci, beddi,bjmh,bjms,bjph,bjps;
   int hIbeg,hIend,uIbeg,uIend,bIbeg,bIend;


   int *RegInd = RegSplit(h,n); 
   int m1 = 5;
   int RegIndLen = RegInd[m1-1];


    //First Do beg initial conditions
    // B.C stuff

    bedhbc[0] = bMbeg[nbBC-4];
    bedhbc[1] = bMbeg[nbBC-3];
    bedhbc[2] = bMbeg[nbBC-2];
    bedhbc[3] = bMbeg[nbBC-1];

    i = n-1;
    bedhbc[4*(i+1) + (nbBC)] = bMend[0];
    bedhbc[4*(i+1) + 1 + (nbBC)] = bMend[1];
    bedhbc[4*(i+1) + 2 + (nbBC)] = bMend[2];
    bedhbc[4*(i+1) + 3 + (nbBC)] = bMend[3];

    for(i=0;i < nGhBC;i++)
    {
        // Assuming bed constant in ghost cells
        whbc[i] = wMbeg[i];
        whbc[nGhbc-nGhBC + i] = wMend[i];
        hhbc[i] = hMbeg[i];
        hhbc[nGhbc-nGhBC + i] = hMend[i];
        Ghbc[i] = GMbeg[i];
        Ghbc[nGhbc-nGhBC +i] = GMend[i];
    }

    // Get middle and far edge from this
    memcpy(u,uMbeg, (unBC-1)*sizeof(double));
    memcpy(u+nubc - (unBC-1),uMend + 1,(unBC-1)*sizeof(double));

    //Ensure thatall the zeroing has been done, so we can read the boundary conditions from the neighbours
    for (i =0;i<RegIndLen;i++)
    {
        if(RegInd[i*m1 +  3] != 1)
        {
           for(j=RegInd[i*m1]; j <= RegInd[i*m1 +  1] ; j++)
           {
                u[(unBC-1) + 2*j] = 0;
                u[(unBC-1) + 2*j +1] = 0;

                hhbc[nGhBC + 3*j] = 0;
                hhbc[nGhBC + 3*j + 1] = 0;
                hhbc[nGhBC + 3*j + 2] = 0;

                Ghbc[nGhBC + 3*j] = 0;
                Ghbc[nGhBC + 3*j + 1] = 0;
                Ghbc[nGhBC + 3*j + 2] = 0;


                // reconstruct bed
                if(j == 1) 
                {
                    bedbegmiddle = -bMbeg[0]*i16 + 9*bMbeg[1]*i16 + 9*bMbeg[2]*i16 - bMbeg[3]*i16;

                    bedai = i12*idx*idx*idx*(-bedbegmiddle + 2*bed[j-1] -2*bed[j+1] + bed[j+2]);
                    bedbi = i6*idx*idx*(bedbegmiddle - bed[j-1] - bed[j+1] + bed[j+2]);
                    bedci = i12*idx*(bedbegmiddle - 8*bed[j-1] + 8*bed[j+1] - bed[j+2]);
                    beddi = -bedbegmiddle*i6 + 2*i3*bed[j-1] + 2*i3*bed[j+1] - bed[j+2]*i6;
                }
                else if(j == 0)
                {
                    bedbegmiddle = -bMbeg[0]*i16 + 9*bMbeg[1]*i16 + 9*bMbeg[2]*i16 - bMbeg[3]*i16;

                    bedai = i6*idx*idx*idx*(-bedbegmiddle + 3*bed[j] - 3*bed[j + 1] + bed[j + 2]);
                    bedbi = 0.5*idx*idx*(bedbegmiddle  - 2*bed[j] + bed[j + 1]);
                    bedci = idx*(-bedbegmiddle/3.0 - bed[j]/2.0 + bed[j + 1] - bed[j + 2]/6.0);
                    beddi = bed[j];

                }
                else if(j == n-2)
                {
                    bedendmiddle = -bMend[0]*i16 + 9*bMend[1]*i16 + 9*bMend[2]*i16 - bMend[3]*i16;

                    bedai = i12*idx*idx*idx*(-bed[j-2] + 2*bed[j-1] -2*bed[j+1] + bedendmiddle);
                    bedbi = i6*idx*idx*(bed[j-2] - bed[j-1] - bed[j+1] + bedendmiddle);
                    bedci = i12*idx*(bed[j-2] - 8*bed[j-1] + 8*bed[j+1] - bedendmiddle);
                    beddi = -bed[j-2]*i6 + 2*i3*bed[j-1] + 2*i3*bed[j+1] - bedendmiddle*i6;

                }
                else if(j == n-1)
                {
                    bedendmiddle = -bMend[0]*i16 + 9*bMend[1]*i16 + 9*bMend[2]*i16 - bMend[3]*i16;

                    bedai = i6*idx*idx*idx*(bedendmiddle - 3*bed[j] + 3*bed[j -1]- bed[j-2]);
                    bedbi = 0.5*idx*idx*(bedendmiddle - 2*bed[j] + bed[j -1]);
                    bedci = idx*(bedendmiddle/3.0 + bed[j]/2.0 - bed[j -1] + bed[j -2]/6.0);
                    beddi = bed[j];

                }
                else
                {
                    bedai = i12*idx*idx*idx*(-bed[j-2] + 2*bed[j-1] -2*bed[j+1] + bed[j+2]);
                    bedbi = i6*idx*idx*(bed[j-2] - bed[j-1] - bed[j+1] + bed[j+2]);
                    bedci = i12*idx*(bed[j-2] - 8*bed[j-1] + 8*bed[j+1] - bed[j+2]);
                    beddi = -bed[j-2]*i6 + 2*i3*bed[j-1] + 2*i3*bed[j+1] - bed[j+2]*i6;

                }


                bjmh = -bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) - bedci*(0.5*dx) + beddi;
                bjms = -bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) - bedci*(dx*i6) + beddi;
                bjps = bedai*(dx*i6)*(dx*i6)*(dx*i6) + bedbi*(dx*i6)*(dx*i6) + bedci*(dx*i6) + beddi;
                bjph = bedai*(0.5*dx)*(0.5*dx)*(0.5*dx) + bedbi*(0.5*dx)*(0.5*dx) + bedci*(0.5*dx) + beddi;

                if(j == n-1)
                {

                    bedhbc[nbBC + 4*j -1 ] =0.5*( bjmh + bedhbc[nbBC + 4*j -1 ]);
                    bedhbc[nbBC + 4*j] = bedhbc[nbBC + 4*j -1 ];
                    bedhbc[nbBC + 4*j +1] = bjms;
                    bedhbc[nbBC + 4*j +2] = bjps;
                    bedhbc[nbBC + 4*j + 3] = 0.5*(bjph + bedhbc[nbBC + 4*j + 4]);
                    bedhbc[nbBC + 4*j + 4] = bedhbc[nbBC + 4*j + 3] ;

                }
                else
                {

                    bedhbc[nbBC + 4*j -1 ] =0.5*( bjmh + bedhbc[nbBC + 4*j -1 ]);
                    bedhbc[nbBC + 4*j] = bedhbc[nbBC + 4*j -1 ];
                    bedhbc[nbBC + 4*j +1] = bjms;
                    bedhbc[nbBC + 4*j +2] = bjps;
                    bedhbc[nbBC + 4*j + 3] = bjph;

                }

                whbc[nGhBC + 3*j] =  bjmh;
                whbc[nGhBC + 3*j + 1] = beddi;
                whbc[nGhBC + 3*j + 2] =  bjph;
                
           }
        }
    }

    for (i =0;i<RegIndLen;i++)
    {
        if(RegInd[i*m1 +  3] == 1)
        {
            hIbeg = nGhBC + 3*RegInd[i*m1];
            hIend = nGhBC + 3*(RegInd[i*m1 + 1]) + 3;

            uIbeg = (unBC-1) + 2*RegInd[i*m1];
            uIend = (unBC-1) + 2*(RegInd[i*m1 + 1]) + 2;
            
            bIbeg = nbBC + 4*RegInd[i*m1];
            bIend = nbBC + 4*(RegInd[i*m1 + 1]) + 4;
            // give h + offset for h
            // use hhbc to get hbeg
            // w

            getufromG(h + RegInd[i*m1] ,G + RegInd[i*m1] ,bed + RegInd[i*m1],hhbc+hIbeg -nGhBC,hhbc+hIend,Ghbc+hIbeg -nGhBC,Ghbc+hIend,u+uIbeg -unBC,u+uIend,whbc+hIbeg -nGhBC,whbc+hIend,bedhbc + bIbeg - nbBC,bedhbc + bIend,theta,dx , (RegInd[i*m1 + 2] +1) , 2*(RegInd[i*m1+2]+1) + 1,nGhBC,unBC,nbBC,nGhbc,nubc,nbbc,u + uIbeg,hhbc+ hIbeg ,Ghbc+hIbeg ,whbc+hIbeg ,bedhbc+bIbeg );
        }
    }

    free(RegInd);


}




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


//include BCs
void evolveForce(double *Ghbc, double *hhbc, double *whbc, double *ubc,double *bedhbc, int nGhhbc, int nubc, int nbhc, int nGhBC, int unBC, int bnBC, double g, double dx, double dt, int n, double theta, double *newG, double *newh, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9)
{
    double idx = 1.0 / dx;  
	int i;
    double her,Ger,dber,uer,duer,hel,Gel,dbel,uel,duel,fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG,th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
	double wil,wir,wip1l,bip1l,bil,bir,nbi,hihm,hihp,uai,ubi,uaip1,ubip1;
    double himhp,bedai,bedbi,bedci,bedaip1,bedbip1,bedcip1;
    double dui,duif,duim,duib,duip1,duip1f,duip1m,duip1b,ujphm,ujphp;

    double hS, GS;



    // i = -1
    i = -1;

    //Figure out how the interior numberings work.
    
    wil = whbc[3*i + nGhBC];
    wir = whbc[3*i + nGhBC + 2];

    // bil and bir
    bil = wil - hhbc[3*i + nGhBC];
    bir = wir - hhbc[3*i + nGhBC + 2];

    wip1l = whbc[3*(i+1) + nGhBC];
    //wip1r = whbc[3*(i+1) + nGhBC + 2];

    // bil and bir
    bip1l = wip1l - hhbc[3*(i+1) +nGhBC];
    //bip1r = wip1r - hhbc[3*(i+1) +nGhBC +2];

    // new bed and height
    nbi = fmax(bip1l, bir);
	hihm = fmax(0, wir - nbi);
	hihp = fmax(0, wip1l - nbi);

    uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
    ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

    bedai = 0.5*9*idx*idx*idx*(-bedhbc[4*i + bnBC]+ 3*bedhbc[4*i + bnBC+1] - 3*bedhbc[4*i + bnBC+2] + bedhbc[4*i + bnBC +3]);
    bedbi = 0.25*9*idx*idx*(bedhbc[4*i + bnBC] - bedhbc[4*i + bnBC+1] - bedhbc[4*i + bnBC+2] + bedhbc[4*i + bnBC+3]);
    bedci = i8*idx*(bedhbc[4*i + bnBC] - 27*bedhbc[4*i + bnBC+1] + 27*bedhbc[4*i + bnBC+2]- bedhbc[4*i + bnBC+3]);
    //beddi = -bedhbc[4*i + bnBC]*i16 + 9*bedhbc[4*i + bnBC+1]*i16 + 9*bedhbc[4*i + bnBC+2]*i16 - bedhbc[4*i + bnBC]*i16;

    uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
    ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

    bedaip1 = 0.5*9*idx*idx*idx*(-bedhbc[4*(i+1) + bnBC]+ 3*bedhbc[4*(i+1) + bnBC+1] - 3*bedhbc[4*(i+1) + bnBC+2] + bedhbc[4*(i+1) + bnBC +3]);
    bedbip1 = 0.25*9*idx*idx*(bedhbc[4*(i+1) + bnBC] - bedhbc[4*(i+1) + bnBC+1] - bedhbc[4*(i+1) + bnBC+2] + bedhbc[4*(i+1) + bnBC+3]);
    bedcip1 = i8*idx*(bedhbc[4*(i+1) + bnBC] - 27*bedhbc[4*(i+1) + bnBC+1] + 27*bedhbc[4*(i+1) + bnBC+2]- bedhbc[4*(i+1) + bnBC+3]);
    //beddip1 = -bedhbc[4*(i+1) + bnBC]*i16 + 9*bedhbc[4*(i+1) + bnBC+1]*i16 + 9*bedhbc[4*(i+1) + bnBC+2]*i16 - bedhbc[4*(i+1) + bnBC]*i16;



    her = hihp;
    Ger = Ghbc[3*(i+1) +nGhBC];
    //ber = bip1l;
    dber = 3*bedaip1*(0.5*dx)*(0.5*dx) - 2*bedbip1*(0.5*dx) + bedcip1;
    uer  = ubc[2*i + unBC + 1];
    duer = -uaip1*(dx) + ubip1; //(2*0.5*dx)

    hel = hihm;
    Gel = Ghbc[3*(i) +nGhBC +2];
    //bel = bir;
    dbel = 3*bedai*(0.5*dx)*(0.5*dx) + 2*bedbi*(0.5*dx) + bedci;
    uel  = ubc[2*i + unBC + 1];
    duel = uai*(dx) + ubi; //special formula

    fhel = uel*hel;
    fher = uer*her;

    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

    //printf("%d | %e | %e \n",i,fGel,fGer);

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
    himhp = hihp;
    for(i = 0;i < n;i++)
    {
        //Figure out how the interior numberings work.
        
	    wil = whbc[3*i + nGhBC];
	    wir = whbc[3*i + nGhBC + 2];

	    // bil and bir
	    bil = wil - hhbc[3*i + nGhBC];
	    bir = wir - hhbc[3*i + nGhBC + 2];

	    wip1l = whbc[3*(i+1) + nGhBC];
	    //wip1r = whbc[3*(i+1) + nGhBC + 2];

	    // bil and bir
	    bip1l = wip1l - hhbc[3*(i+1) +nGhBC];
	    //bip1r = wip1r - hhbc[3*(i+1) +nGhBC +2];

        // new bed and height
        nbi = fmax(bip1l, bir);
		hihm = fmax(0, wir - nbi);
		hihp = fmax(0, wip1l - nbi);

        bedai = 0.5*9*idx*idx*idx*(-bedhbc[4*i + bnBC]+ 3*bedhbc[4*i + bnBC+1] - 3*bedhbc[4*i + bnBC+2] + bedhbc[4*i + bnBC +3]);
        bedbi = 0.25*9*idx*idx*(bedhbc[4*i + bnBC] - bedhbc[4*i + bnBC+1] - bedhbc[4*i + bnBC+2] + bedhbc[4*i + bnBC+3]);
        bedci = i8*idx*(bedhbc[4*i + bnBC] - 27*bedhbc[4*i + bnBC+1] + 27*bedhbc[4*i + bnBC+2]- bedhbc[4*i + bnBC+3]);
        //beddi = -bedhbc[4*i + bnBC]*i16 + 9*bedhbc[4*i + bnBC+1]*i16 + 9*bedhbc[4*i + bnBC+2]*i16 - bedhbc[4*i + bnBC]*i16;

        uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
        ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

        bedaip1 = 0.5*9*idx*idx*idx*(-bedhbc[4*(i+1) + bnBC]+ 3*bedhbc[4*(i+1) + bnBC+1] - 3*bedhbc[4*(i+1) + bnBC+2] + bedhbc[4*(i+1) + bnBC +3]);
        bedbip1 = 0.25*9*idx*idx*(bedhbc[4*(i+1) + bnBC] - bedhbc[4*(i+1) + bnBC+1] - bedhbc[4*(i+1) + bnBC+2] + bedhbc[4*(i+1) + bnBC+3]);
        bedcip1 = i8*idx*(bedhbc[4*(i+1) + bnBC] - 27*bedhbc[4*(i+1) + bnBC+1] + 27*bedhbc[4*(i+1) + bnBC+2]- bedhbc[4*(i+1) + bnBC+3]);
        //beddip1 = -bedhbc[4*(i+1) + bnBC]*i16 + 9*bedhbc[4*(i+1) + bnBC+1]*i16 + 9*bedhbc[4*(i+1) + bnBC+2]*i16 - bedhbc[4*(i+1) + bnBC]*i16;

        her = hihp;
	    Ger = Ghbc[3*(i+1) +nGhBC];
	    //ber = bip1l;
	    dber = 3*bedaip1*(0.5*dx)*(0.5*dx) - 2*bedbip1*(0.5*dx) + bedcip1;
	    uer  = ubc[2*i + unBC + 1];
	    duer = -uaip1*(dx) + ubip1; //(2*0.5*dx)

	    hel = hihm;
	    Gel = Ghbc[3*(i) +nGhBC +2];
	    //bel = bir;
	    dbel = 3*bedai*(0.5*dx)*(0.5*dx) + 2*bedbi*(0.5*dx) + bedci;
	    uel  = ubc[2*i + unBC + 1];
	    duel = uai*(dx) + ubi; //special formula

	    fhel = uel*hel;
	    fher = uer*her;
        //printf("h : %d | %e | %e \n",i,fhel,fher);

	
	    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
	    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;
        //printf("G : %d | %e | %e \n",i,fGel,fGer);

        //printf("%d | %e | %e \n",i,fGel,fGer);

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);	

	    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));


        //centerted values
        th = hhbc[3*i + nGhBC+ 1];
		tu = ubc[2*i + unBC];
		tux = ubi;
		tbx = idx*(bir - bil);
		tbxx =2*bedbi;
		
		sourcer = g*0.5*(hihm*hihm - hhbc[3*(i) +nGhBC +2]*hhbc[3*(i) +nGhBC +2]);
		sourcec = -g*th*tbx -  0.5*th*th*tu*tux*tbxx + th*tu*tu*tbx*tbxx ;
		sourcel = g*0.5*(hhbc[3*(i) +nGhBC]*hhbc[3*(i) +nGhBC] - himhp*himhp);

        //printf("hchange : %d | %e \n ",i,- dt*idx*(foh - fih));
        //printf("Gchange : %d | %e \n ",i,- dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + dx*sourcec));
    
		newh[i] = hhbc[3*(i) +nGhBC + 1] - dt*idx*(foh - fih);
		newG[i] = Ghbc[3*(i) +nGhBC + 1] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + dx*sourcec);

		newh[i] = newh[i]*(newh[i] > htol);
		newG[i] = newG[i]*(newh[i] > htol);


	    fih = foh;
	    fiG = foG;
        himhp = hihp;

    }
        

}


void evolvewrapForcing(double *G, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *GMbeg , double *GMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend,double *hMbeg1, double *hMend1,double *wMbeg1, double *wMend1,double *GMbeg1, double *GMend1,double *uMbeg1, double *uMend1, double g, double dx, double dt, int n, int nGhBC, int unBC,int bnBC, int nGhhbc, int nubc, int nbhc , double theta, double *hhbc, double *whbc,double *Ghbc,double *bedhbc,double *ubc, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9)
{
    double idx = 1.0 /dx;
    //n is number of cells
    int u_length = 2*n + 1;
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

    //Solve for u given wet/dry regions
    getufromGsplit(h,G,bed,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,wMbeg,wMend,bMbeg,bMend,theta,dx ,n,u_length ,nGhBC,unBC,bnBC,nGhhbc,nubc,nbhc,ubc,hhbc,Ghbc,whbc,bedhbc);
    //Get fluxes and sources
    evolveForce(Ghbc, hhbc,whbc, ubc, bedhbc,nGhhbc,nubc,nbhc,nGhBC,unBC,bnBC,g, dx,dt,n, theta ,Gp, hp,x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);

    //Solve for u given wet/dry regions
    getufromGsplit(hp,Gp,bed,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,wMbeg,wMend,bMbeg,bMend,theta,dx ,n,u_length ,nGhBC,unBC,bnBC,nGhhbc,nubc,nbhc,ubc,hhbc,Ghbc,whbc,bedhbc);
    //Get fluxes and sources
    evolveForce(Ghbc, hhbc,whbc, ubc, bedhbc,nGhhbc,nubc,nbhc,nGhBC,unBC,bnBC,g, dx,dt,n, theta ,Gpp, hpp,x,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9); 


    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gpp[i]);
        h[i] = 0.5*(h[i] + hpp[i]);
        //h[i] = h[i]*(h[i] > htol);
        //G[i] = G[i]*(h[i] > htol);
        
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
