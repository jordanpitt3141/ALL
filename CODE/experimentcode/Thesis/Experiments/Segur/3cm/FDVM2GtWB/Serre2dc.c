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
#define div_0 1.0e-20

#define htol 1.0e-6
#define hbase 1.0e-3

#define RegTol 30


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


//Region Split


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

int *RegCondense(int *Reg, int m)
{
    int arr[5];
    int i,j,Cregval,regval,Csi,Cei;
    int ini = 0;
    int *Reg2,*Reg1;

    for(i = 0; i < m;i++)
    {
        Csi = Reg[5*i];
        Cregval = Reg[5*i + 3];
        for(j = i; j < m;j++)
        {
            regval = Reg[5*j + 3];
            if(Cregval != regval)
            {
                i = j-1;
                break;

            }
        }

        Cei = Reg[5*(i) + 1];
        arr[0] = Csi;
        arr[1] = Cei;
        arr[2] = Cei - Csi;
        arr[3] = Cregval;
        arr[4] = ini;

        Reg1 = append(Reg2, ini, 5, arr);
        ini = ini + 1;

        Reg2 = (int*) malloc((ini*5)*sizeof(int)); 
        memcpy(Reg2, Reg1, (ini*5)*sizeof(int));
        free(Reg1);

        
    }

    for(j =0;j < ini;j++)
    {
       Reg2[j*5 +4] = ini;
    }

    return Reg2;

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
            Cei = i;
            if((Cei - Csi) > RegTol)
            {
                arr[0] = Csi;
                arr[1] = Cei;
                arr[2] = (Cei - Csi);
                arr[3] = Cregval;
            }
            else
            {
                arr[0] = Csi;
                arr[1] = Cei;
                arr[2] = (Cei - Csi);
                arr[3] = 0;
            }

            Reg1 = append(Reg, ini, 5, arr);

            Csi = i ;
            Cregval =regvali ;
            ini = ini + 1;

            Reg = (int*) malloc((ini*5)*sizeof(int)); 
            memcpy(Reg, Reg1, (ini*5)*sizeof(int));

            free(Reg1);

        }
        else if(regvali == Cregval && i==n-1 )
        {
            Cei = i+1;
            arr[0] = Csi;
            arr[1] = Cei;
            arr[2] = (Cei - Csi);
            arr[3] = Cregval;

            Reg1 = append(Reg, ini, 5, arr);

            ini = ini + 1;

            Reg = (int*) malloc((ini*5)*sizeof(int)); 
            memcpy(Reg, Reg1, (ini*5)*sizeof(int));

            free(Reg1);
            
        }
    }

    for(j =0;j < ini;j++)
    {
       Reg[j*5 +4] = ini;
    }


    Reg1 = RegCondense(Reg,ini);
    Reg = (int*) malloc((ini*5)*sizeof(int)); 
    memcpy(Reg, Reg1, (ini*5)*sizeof(int));

    free(Reg1);
    //now condense
    
    return Reg;
}

void getufromGFD(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *u)
{
	double idx = 1.0 / dx;


    double **AL, **Ared;
    unsigned long *indx;

    Ared = dmatrix(1, n, 1, 3);
    AL = dmatrix(1, n, 1, 3);
    indx = mallocLongPy(n);

    double *nGis = malloc((n)*sizeof(double));

	double th,thx,tbx,tbxx,D,ai,bi,ci;
	int i;


	for(i = 1; i < n-1; i++)
	{
		th = h[i]*( (h[i] + hbase) / (h[i] + htol));
		thx = 0.5*idx*(h[i+1] - h[i-1]);
		tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

		D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx;
		
		ai = -i3*idx*idx*th*th + 0.5*idx*th*thx;
		bi = D + 2.0*i3*idx*idx*th*th;
		ci = -i3*idx*idx*th*th - 0.5*idx*th*thx;


        Ared[i + 1][1] = ai;
        Ared[i + 1][2] = bi;
        Ared[i + 1][3] = ci;
        nGis[i] = G[i] / ( th + htol);
	}

	//Boundary 
	//i = 0
	i = 0;
	th = h[i]*( (h[i] + hbase) / (h[i] + htol));
	thx = 0.5*idx*(h[i+1] - h[i-1]);
	tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

	D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx;
	
	ai = -i3*idx*idx*th*th + 0.5*idx*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th;
	ci = -i3*idx*idx*th*th - 0.5*idx*th*thx;

    Ared[i + 1][2] = bi;
    Ared[i + 1][3] = ci;
    nGis[i] = G[i]/ ( th + htol) - u0*ai;

	i = n-1;
	
	th = h[i]*( (h[i] + hbase) / (h[i] + htol));
	thx = 0.5*idx*(h[i+1] - h[i-1]);
	tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

	D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx;
	
	ai = -i3*idx*idx*th*th + 0.5*idx*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th;
	ci = -i3*idx*idx*th*th - 0.5*idx*th*thx;


    Ared[i + 1][1] = ai;
    Ared[i + 1][2] = bi;
    nGis[i] = G[i]/ ( th + htol)  - u1*ci;

	
    double d;
    d = bandec(Ared, n, 1,1, AL ,indx);

    banbks(Ared, n,1,1, AL, indx, nGis);


    memcpy(u,nGis,n*sizeof(double));

    free_dmatrix(Ared, 1, n, 1, 3);
    free_dmatrix(AL, 1,n, 1, 3);
    free(indx);
    free(nGis);

}

void uSolveandRecon(double *h, double *G, double *b,double *hMbeg, double *hMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend, double *duMbeg, double *duMend,int n, int hnBC , int hnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *ubc, double *dubc,int *RegIndex ,int m,int RegLen)
{
    int i,j;
    double ujm1,ujp1;
    double idx = 1.0 / dx;

    for(i = 0;i < unBC;i++)
    {
        ubc[i] = uMbeg[i];
        ubc[unbc - unBC + i] = uMend[i];

        dubc[i] = duMbeg[i];
        dubc[unbc - unBC + i] = duMend[i];
        
    }

    double *u = malloc((n)*sizeof(double));
    getufromGFD(h,G,b,uMbeg[unBC-2],uMend[1],hMbeg[hnBC-2],hMend[1],bMbeg[hnBC-2],bMend[1],dx,n,u);

    for(i =0;i<RegLen;i++)
    {
        if(RegIndex[i*m +  3] == 0)
        {
           for(j=RegIndex[i*m]; j < RegIndex[i*m +  1] ; j++)
           {


                if(j == RegIndex[i*m])
                {
		            if(j == 0)
		            {
		                ujm1 = uMbeg[unBC-2];
		            }
		            else
		            {
		                ujm1 = u[j-1];
		            }
		            
		            ubc[unBC + 2*(j)] = 0;
		            ubc[unBC + 2*(j) + 1] = 0;

                    dubc[unBC + 2*(j) - 1] =idx*(0 - ujm1);
                    dubc[unBC + 2*(j)] = 0;
                    dubc[unBC + 2*(j) + 1] = 0;
                }
                else if(j == RegIndex[i*m +  1] - 1)
                {
		            if(j == n-1)
		            {
		                ujp1 = uMend[1];
		            }
		            else
		            {
		                ujp1 = u[j+1];
		            }


					ubc[unBC + 2*(j) - 1] = 0;
		            ubc[unBC + 2*(j)] = 0;
		            
                    dubc[unBC + 2*(j) - 1] =0;
                    dubc[unBC + 2*(j)] = 0;
                    dubc[unBC + 2*(j) + 1] = idx*(ujp1 - 0);
                }
                else
                {
                    ubc[unBC + 2*(j)] = 0;
                    ubc[unBC + 2*(j) + 1] = 0;

                    dubc[unBC + 2*(j)] = 0;
                    dubc[unBC + 2*(j) + 1] = 0;
                }
           }
        }

        if(RegIndex[i*m +  3] == 1)
        {
           for(j=RegIndex[i*m]; j < RegIndex[i*m +  1] ; j++)
           {

                if(j == 0)
                {
                    ujm1 = uMbeg[unBC-2];
                }
                else
                {
                    ujm1 = u[j-1];
                }

                if(j == n-1)
                {
                    ujp1 = uMend[1];
                }
                else
                {
                    ujp1 = u[j+1];
                }


                if(j == RegIndex[i*m])
                {
                    ubc[unBC + 2*(j) - 1] = 0.5*(u[j] + ujm1);
                    ubc[unBC + 2*(j)] = u[j];
                    ubc[unBC + 2*(j) + 1] = 0.5*(u[j] + ujp1);

                    dubc[unBC + 2*(j) - 1] =idx*(u[j] - ujm1);
                    dubc[unBC + 2*(j)] = 0.5*idx*(ujp1 - ujm1);
                    dubc[unBC + 2*(j) + 1] = idx*(ujp1 - u[j]);


                }
                else
                {
                    ubc[unBC + 2*(j)] = u[j];
                    ubc[unBC + 2*(j) + 1] = 0.5*(u[j] + ujp1);

                    dubc[unBC + 2*(j)] = 0.5*idx*(ujp1 - ujm1);
                    dubc[unBC + 2*(j) + 1] = idx*(ujp1 - u[j]);
                }
           }
        }

    }


    free(u);

}

void ReconLin( double *q, double *qMbeg, double *qMend, int n, int nBC, int nbc, double theta, double *qbc)
{
    double dqim,dqib,dqif,dqi,qimh,qiph;
    int i;


    for(i = 0;i < nBC;i++)
    {
        qbc[i] = qMbeg[i];
        qbc[nbc - nBC + i] = qMend[i];
        
    }

    i = 0;
    dqib = (q[i] - qMbeg[nBC - 2]);
    dqim = 0.5*(q[i+1] - qMbeg[nBC - 2]);
    dqif = (q[i+1] - q[i]);

    dqi = minmod(theta*dqib,dqim,theta*dqif);

    qimh = q[i] - 0.5*dqi;
    qiph = q[i] + 0.5*dqi;

    qbc[nBC + 3*(i)] = qimh;
    qbc[nBC + 3*(i) + 1] = q[i];
    qbc[nBC + 3*(i) + 2] = qiph;

    for(i = 1;i < n -1;i++)
    {
        dqib = (q[i] - q[i -1]);
        dqim = 0.5*(q[i+1] - q[i -1]);
        dqif = (q[i+1] - q[i]);

        dqi = minmod(theta*dqib,dqim,theta*dqif);

        qimh = q[i] - 0.5*dqi;
        qiph = q[i] + 0.5*dqi;

        qbc[nBC + 3*(i)] = qimh;
        qbc[nBC + 3*(i) + 1] = q[i];
        qbc[nBC + 3*(i) + 2] = qiph;

    }

    i = n-1;

    dqib = (q[i] - q[i -1]);
    dqim = 0.5*(qMend[1] - q[i -1]);
    dqif = (qMend[1] - q[i]);

    dqi = minmod(theta*dqib,dqim,theta*dqif);

    qimh = q[i] - 0.5*dqi;
    qiph = q[i] + 0.5*dqi;

    qbc[nBC + 3*(i)] = qimh;
    qbc[nBC + 3*(i) + 1] = q[i];
    qbc[nBC + 3*(i) + 2] = qiph;

}

void ReconZero(double *qbc, int n, int nBC, int nbc, int *RegIndex, int m, int RegLen)
{
   int i,j;
   for(i =0;i<RegLen;i++)
    {
        if(RegIndex[i*m +  3] == 0)
        {
           for(j=RegIndex[i*m]; j < RegIndex[i*m +  1] ; j++)
           {

                qbc[nBC + 3*(j)] = 0;
                qbc[nBC + 3*(j) + 1] = 0;
                qbc[nBC + 3*(j) + 2] = 0;
           }
        }

    }
}

void Reconbfromwh(double *hbc,double *wbc, int nbc,double *bbc)
{
    int i;
    
   for(i =0;i<nbc;i++)
   {
        bbc[i] = wbc[i] - hbc[i];

   }
    

}


//include BCs
void evolveForce(double *hbc, double *Gbc, double *wbc, double *bMbc, double *dbMbc, double *ddbCbc, double *uFDbc,double *duFDbc,double g, double dx, double dt, int n, int hnBC, int hnbc, int unBC, int unbc, int CnBC, int Cnbc, double *FjhH, double *FjhG, double *SjG)
{
    double idx = 1.0 / dx;  
	int i;
    double her,Ger,dber,uer,duer,hel,Gel,dbel,uel,duel,fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG,th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
	double wil,wir,wip1l,bip1l,bil,bir,nbi,hihm,hihp,hir,hip1l,hil;
    double himhp;

    double hS, GS;


    // i = -1
    i = -1;
    
    wip1l = wbc[hnBC + 3*(i+1)];
    wir = wbc[hnBC + 3*(i) + 2];
    wil = wbc[hnBC + 3*(i)];

    hip1l = hbc[hnBC + 3*(i+1)];
    hir = hbc[hnBC + 3*(i) + 2];
    hil = hbc[hnBC + 3*(i)];

    bip1l = wip1l - hip1l;
    bir = wir - hir;
    bil = wil - hil;

    nbi = fmax(bip1l, bir);
    hihm = fmax(0, wir - nbi);
    hihp = fmax(0, wip1l - nbi);  

    her = hihp ;
    Ger = Gbc[hnBC + 3*(i+1)];
    dber = dbMbc[hnBC + 3*(i+1)];
    uer  = uFDbc[unBC + 2*(i) + 1];
    duer = duFDbc[unBC + 2*(i) + 1];


    hel = hihm ;
    Gel = Gbc[hnBC + 3*(i) + 2];
    dbel = dbMbc[hnBC + 3*(i) + 2];;
    uel  = uFDbc[unBC + 2*(i) + 1];
    duel = duFDbc[unBC + 2*(i) + 1];

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

    FjhH[0] = foh;
    FjhG[0] = foG;
    
    himhp = hihp;

    for(i = 0;i < n;i++)
    {
    
        wip1l = wbc[hnBC + 3*(i+1)];
		wir = wbc[hnBC + 3*(i) + 2];
		wil = wbc[hnBC + 3*(i)];

		hip1l = hbc[hnBC + 3*(i+1)];
		hir = hbc[hnBC + 3*(i) + 2];
		hil = hbc[hnBC + 3*(i)];

		bip1l = wip1l - hip1l;
		bir = wir - hir;
		bil = wil - hil;

		nbi = fmax(bip1l, bir);
		hihm = fmax(0, wir - nbi);
		hihp = fmax(0, wip1l - nbi);  

        her = hihp ;
        Ger = Gbc[hnBC + 3*(i+1)];
        dber = dbMbc[hnBC + 3*(i+1)];
        uer  = uFDbc[unBC + 2*(i) + 1];
        duer = duFDbc[unBC + 2*(i) + 1];


        hel = hihm;
        Gel = Gbc[hnBC + 3*(i) + 2];
        dbel = dbMbc[hnBC + 3*(i) + 2];;
        uel  = uFDbc[unBC + 2*(i) + 1];
        duel = duFDbc[unBC + 2*(i) + 1];

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
        th = hbc[hnBC + 3*(i) + 1];
		tu =  uFDbc[2*i + unBC];
		tux = duFDbc[2*i + unBC];
		tbx = dbMbc[hnBC + 3*(i) + 1];
		tbxx = ddbCbc[CnBC + (i)];
		
		sourcer = g*0.5*(hihm*hihm - hir*hir);
		sourcec = -g*th*tbx -  0.5*th*th*tu*tux*tbxx + th*tu*tu*tbx*tbxx ;
		sourcel = g*0.5*(hil*hil - himhp*himhp);

        FjhH[i + 1] = foh;
        FjhG[i + 1] = foG;
        SjG[i] = dx*sourcec + sourcer + sourcel;
        
        himhp = hihp;

    }

}

void Fluxsplit(double *hbc, double *Gbc, double *wbc, double *bMbc, double *dbMbc, double *ddbCbc, double *uFDbc,double *duFDbc,double g, double dx, double dt, int n, int hnBC, int hnbc, int unBC, int unbc, int CnBC, int Cnbc, int *RegIndx, int m1 ,int RegIndxLen, double *FjhH, double *FjhG, double *SjG)
{

    int hIbeg,uIbeg,uIend,bIbeg,FIbeg,SIbeg,CIbeg;
    double hS, GS;
    //Dry First

    int i,j,i1;
    for (i =0;i<RegIndxLen;i++)
    {

        //printf("%d : %d | %d | %d | %d  \n",i,RegIndx[i*m1],RegIndx[i*m1 +  1],RegIndx[i*m1 +  2],RegIndx[i*m1 +  3]    );
        if(RegIndx[i*m1 +  3] == 0)
        {

           if(RegIndx[i*m1] == 0)
           {
               FjhH[0] = 0;
               FjhG[0] = 0;
                
           }

           if(RegIndx[i*m1 +  1] == n)
           {
               FjhH[n] = 0;
               FjhG[n] = 0;
                
           }
                
           //First and last dry cell have no source
           
           SjG[RegIndx[i*m1]] = 0;


           for(j=RegIndx[i*m1] + 1; j < RegIndx[i*m1 +  1]; j++)
           {
                FjhH[j] = 0;
                FjhG[j] = 0;
                SjG[j] = 0;

           }

        }

        if(RegIndx[i*m1 +  3] == 1)
        {
            //printf("%d : %d | %d | %d | %d  \n",i,RegIndx[i*m1],RegIndx[i*m1 +  1],RegIndx[i*m1 +  2],RegIndx[i*m1 +  3]    );
            hIbeg = hnBC + 3*(RegIndx[i*m1] - 1);
            uIbeg = (unBC-1) + 2*(RegIndx[i*m1] - 1);
            CIbeg = CnBC + (RegIndx[i*m1] - 1);
            FIbeg = RegIndx[i*m1];
            SIbeg = RegIndx[i*m1];

           // printf("%d | %d | %d | %d | %d  \n \n",hIbeg,uIbeg,bIbeg,FIbeg,SIbeg   );



            evolveForce(hbc+ hIbeg,Gbc+ hIbeg,wbc+ hIbeg,bMbc+ hIbeg,dbMbc+ hIbeg,ddbCbc+ CIbeg,uFDbc+ uIbeg,duFDbc+ uIbeg,g,dx,dt,RegIndx[i*m1 +  2],hnBC,hnbc,unBC,unbc,CnBC,Cnbc,FjhH + FIbeg,FjhG + FIbeg ,SjG + SIbeg);
        }
    }

}

void FluxEvo(double *hbc, double *Gbc,double g, double dx, double dt, int n, int hnBC, double *FjhH, double *FjhG, double *SjG, double *hp, double *Gp)
{
    double hS,GS;
    int i;
    double idx = 1.0 / dx;
    for (i = 0; i < n; i++)
    {

        hp[i] = hbc[hnBC + 3*(i) + 1] - dt*idx*(FjhH[i+1] - FjhH[i]);
        Gp[i] = Gbc[hnBC + 3*(i) + 1] - dt*idx*(FjhG[i+1] - FjhG[i]) + dt*idx*(SjG[i]);
        //Gp[i] = Gbc[hnBC + 3*(i) + 1] + dt*idx*(SjG[i]) + dt*SjGa[i];

        //printf("h Flux  : %d | %e | %e | %e  \n ",i,- dt*idx*(FjhH[i+1] - FjhH[i]),dt*Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), - dt*idx*(FjhH[i+1] - FjhH[i]) + dt*Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g));
        //printf("G Flux  : %d | %e | %e | %e  \n ",i,- dt*idx*(FjhG[i+1] - FjhG[i]),dt*FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), - dt*idx*(FjhG[i+1] - FjhG[i]) + dt*FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g));
        //printf("G Source: %d | %e | %e | %e  \n ",i,dt*idx*(SjG[i]),dt*SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), dt*idx*(SjG[i])+ dt*SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g) );

    }

}

void FluxWrap(double *hbc, double *Gbc, double *wbc, double *bMbc, double *dbMbc, double *ddbCbc, double *uFDbc, double *duFDbc,double g, double dx, double dt, int n, int hnBC, int hnbc, int unBC, int unbc,int CnBC,int Cnbc, int *RegIndx, int m1 ,int RegIndxLen, double *hp, double *Gp)
{
    double *FjhH = malloc((n+1)*sizeof(double));
    double *FjhG = malloc((n+1)*sizeof(double));
    double *SjG = malloc((n)*sizeof(double));
    
    Fluxsplit(hbc,Gbc,wbc,bMbc,dbMbc,ddbCbc,uFDbc,duFDbc,g,dx,dt,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,RegIndx,m1,RegIndxLen,FjhH,FjhG,SjG);

    FluxEvo(hbc,Gbc,g,dx,dt,n,hnBC,FjhH,FjhG,SjG,hp,Gp);

    free(FjhH);
    free(FjhG);
    free(SjG);
}

void Stagehbed(int n, double *h, double *b, double *w, int *RegionIndex, int m, int RegIndLen)
{
    int i,j;
    for (i =0;i<RegIndLen;i++)
    {
        //printf("Reg %d : %d | %d | %d | %d : %d \n",i,RegInd[i*m1],RegInd[i*m1+1],RegInd[i*m1+2],RegInd[i*m1+3],RegInd[i*m1+4]);
        if(RegionIndex[i*m +  3] == 0)
        {
           for(j=RegionIndex[i*m]; j < RegionIndex[i*m +  1] ; j++)
           {

                w[j] = b[j];
           }
        }
        else if(RegionIndex[i*m +  3] == 1)
        {

           for(j=RegionIndex[i*m]; j < RegionIndex[i*m +  1] ; j++)
           {
                w[j] = h[j] + b[j];
           }

        }
    }

}


void Recondb(double *bbc, double *dbMbeg, double *dbMend, double dx, int n, int hnBC, int hnbc, double *dbbc  )
{
    int i;
    double idx = 1.0 / dx;

    for(i = 0;i < hnBC;i++)
    {
        dbbc[i] = dbMbeg[i];
        dbbc[hnbc - hnBC + i] = dbMend[i];
        
    }

    for(i = 0;i < n ;i++)
    {
        dbbc[hnBC + 3*(i)] = idx*(bbc[hnBC + 3*i + 1] - bbc[hnBC + 3*(i - 1) + 1]);
        dbbc[hnBC + 3*(i) + 1] = idx*(bbc[hnBC + 3*(i) + 2] - bbc[hnBC + 3*(i)]);
        dbbc[hnBC + 3*(i) + 2] = idx*(bbc[hnBC + 3*(i+1) + 1] - bbc[hnBC + 3*i + 1]);

    }

}

void Reconddb(double *b, double *ddbMbeg, double *ddbMend, double *bMbeg  , double *bMend, double dx, int n,int hnBC,int CnBC, int Cnbc, double *ddbbc  )
{
    int i;
    double idx = 1.0 / dx;

    for(i = 0;i < CnBC;i++)
    {
        ddbbc[i] = ddbMbeg[i];
        ddbbc[Cnbc - CnBC + i] = ddbMend[i];
        
    }
    
    i = 0;
    ddbbc[i +CnBC ] = idx*idx*(b[i+ 1] - 2*b[i] + bMbeg[hnBC - 2] );

    for(i = 1;i < n-1 ;i++)
    {
        ddbbc[i +CnBC ] = idx*idx*(b[i+ 1] - 2*b[i] + b[i-1]);

    }

    i = n-1;
    ddbbc[i +CnBC ] = idx*idx*(bMend[1] - 2*b[i] + b[i-1] );

}


void ReconWrap(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int unBC , int unbc, int CnBC , int Cnbc,double theta,double dx, double dt, double g, double *w, double *hbc, double *Gbc, double *wbc , double *bMbc, double *dbMbc, double *ddbCbc,int *RegInd ,int m1,int RegIndLen)
{

    // calculate stage
    Stagehbed(n,h,b,w,RegInd,m1,RegIndLen);

    //Reconstruct
    ReconLin(h, hMbeg, hMend,n,hnBC,hnbc,theta,hbc);
    ReconLin(G, GMbeg, GMend,n,hnBC,hnbc,theta,Gbc);
    ReconLin(w, wMbeg, wMend,n,hnBC,hnbc,theta,wbc);

    ReconZero(hbc,n,hnBC,hnbc,RegInd,m1,RegIndLen);
    ReconZero(Gbc,n,hnBC,hnbc,RegInd,m1,RegIndLen);

    Reconbfromwh(hbc,wbc,hnbc,bMbc);
    Recondb(bMbc,dbMbeg,dbMend,dx,n,hnBC,hnbc,dbMbc);

    Reconddb(b,ddbCbeg,ddbCend,bMbc,bMbc + 3*(n) + hnBC,dx,n,hnBC,CnBC,Cnbc,ddbCbc); 

}

void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend, double *duMbeg, double *duMend,int n, int hnBC , int hnbc, int unBC , int unbc, int CnBC , int Cnbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc, double *bMbc, double *dbMbc, double *ddbCbc,double *uFDbc, double *duFDbc)
{

    double *w = malloc((n)*sizeof(double));
    int u_length = 2*n + 1;
    int *RegInd = RegSplit(h,n); 
    int m1 = 5;
    int RegIndLen = RegInd[m1-1];

    ReconWrap(h,G,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,dbMbeg,dbMend,ddbCbeg,ddbCend,uMbeg,uMend,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,w,hbc,Gbc,wbc,bMbc,dbMbc,ddbCbc,RegInd,m1,RegIndLen);
    uSolveandRecon(h,G,b,hMbeg,hMend,bMbc,bMbc + 3*(n) + hnBC,uMbeg,uMend,duMbeg,duMend,n,hnBC,hnbc, unBC, unbc,theta,dx,dt,g,uFDbc,duFDbc,RegInd,m1,RegIndLen);

    free(RegInd);
    free(w);
}

void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *dbMbeg, double *dbMend,double *ddbCbeg, double *ddbCend, double *uMbeg, double *uMend, double *duMbeg, double *duMend,int n, int hnBC , int hnbc, int unBC , int unbc, int CnBC, int Cnbc,double theta,double dx, double dt, double g)
{
    double idx = 1.0 /dx;
    //n is number of cells
    int u_length = 2*n + 1;
    double *w = malloc((n)*sizeof(double));
    double *hp = malloc((n)*sizeof(double));
    double *wp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

    double *Gbc = malloc((hnbc)*sizeof(double));
    double *hbc = malloc((hnbc)*sizeof(double));
    double *wbc = malloc((hnbc)*sizeof(double));

    double *bMbc = malloc((hnbc)*sizeof(double));
    double *dbMbc = malloc((hnbc)*sizeof(double));

    double *ddbCbc = malloc((Cnbc)*sizeof(double));

    double *ubc = malloc((unbc)*sizeof(double));
    double *uFDbc = malloc((unbc)*sizeof(double));
    double *duFDbc = malloc((unbc)*sizeof(double));

    double *Gpbc = malloc((hnbc)*sizeof(double));
    double *hpbc = malloc((hnbc)*sizeof(double));
    double *wpbc = malloc((hnbc)*sizeof(double));


    double *upbc = malloc((unbc)*sizeof(double));
    double *upFDbc = malloc((unbc)*sizeof(double));
    double *dupFDbc = malloc((unbc)*sizeof(double));



    //Dry Regions
    int *RegInd = RegSplit(h,n); 
    int m1 = 5;
    int RegIndLen = RegInd[m1-1];

    ReconWrap(h,G,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,dbMbeg,dbMend,ddbCbeg,ddbCend,uMbeg,uMend,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,w,hbc,Gbc,wbc,bMbc,dbMbc,ddbCbc,RegInd,m1,RegIndLen);

    uSolveandRecon(h,G,b,hMbeg,hMend,bMbc,bMbc + 3*(n) + hnBC,uMbeg,uMend,duMbeg,duMend,n,hnBC,hnbc, unBC, unbc,theta,dx,dt,g,uFDbc,duFDbc,RegInd,m1,RegIndLen);

    FluxWrap(hbc,Gbc,wbc,bMbc,dbMbc,ddbCbc,uFDbc,duFDbc,g,dx,dt,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,RegInd,m1 ,RegIndLen,hp,Gp);

 
    int *RegIndp  = RegSplit(hp,n); 
    m1 = 5;
    RegIndLen = RegIndp[m1-1];

    ReconWrap(hp,Gp,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,dbMbeg,dbMend,ddbCbeg,ddbCend,uMbeg,uMend,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,wp,hpbc,Gpbc,wpbc,bMbc,dbMbc,ddbCbc,RegIndp,m1,RegIndLen);

    uSolveandRecon(hp,Gp,b,hMbeg,hMend,bMbc,bMbc + 3*(n) + hnBC,uMbeg,uMend,duMbeg,duMend,n,hnBC,hnbc, unBC, unbc,theta,dx,dt,g,upFDbc,dupFDbc,RegIndp,m1,RegIndLen);

    FluxWrap(hpbc,Gpbc,wpbc,bMbc,dbMbc,ddbCbc,upFDbc,dupFDbc,g,dx,dt,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,RegIndp,m1 ,RegIndLen,hpp,Gpp);

    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gpp[i]);  
        h[i] = 0.5*(h[i] + hpp[i]);
        
    }


    free(w);
    free(wp);
    free(hp);
    free(Gp);
    free(Gpp);
    free(hpp);

    free(Gbc);
    free(hbc);
    free(wbc);
    free(bMbc);
    free(ubc);

    free(Gpbc);
    free(hpbc);
    free(wpbc);
    free(upbc);

    free(dbMbc);
    free(ddbCbc);
    free(uFDbc);
    free(duFDbc);
    free(upFDbc);
    free(dupFDbc);


    free(RegInd);
    free(RegIndp);
}


int main()
{
    printf("h");
    return 1;
}
