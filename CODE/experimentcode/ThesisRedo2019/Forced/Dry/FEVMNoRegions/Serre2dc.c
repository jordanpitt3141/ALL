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
const double sq2 = 1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641572735013846230912297024924836055850737212644121497099935831413222665927505592755799950501152782060571470109559971605970274534596862014728517418640889198609552329230484308714321450839762603627995251407989687253396546331808829640620615258352395054745750287759961729835575220337531857011354374603408498;

const double E = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069;

 
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20
#define div_0 1.0e-20

#define htol 1.0e-12
#define hbase 1.0e-8

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




void getufromG(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, double *ubc)
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

    bjmh = bbc[3*i];
    bjms = bbc[3*i +1];
    bjps = bbc[3*i + 2];
    bjph = bbc[3*i + 3];


	hjphm= (hbc[3*(i) + 2]*hbc[3*(i) + 2] + hbase) / (hbc[3*(i) + 2]);
	hjmhp=  (hbc[3*(i)]*hbc[3*(i) ] + hbase) / (hbc[3*(i)]);

    Gjphm= Gbc[3*(i) + 2];
    Gjmhp= Gbc[3*(i)];



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

    nGis[j-1] = uMbeg[2];
    nGis[j] = nGis[j] + Gintia21; 
    nGis[j+1] = nGis[j+1] + Gintia31;     



    j = 3;
    for (i =1;i < n-1 ; i++)
    {
        bjmh = bbc[3*i];
        bjms = bbc[3*i +1];
        bjps = bbc[3*i + 2];
        bjph = bbc[3*i + 3];


		hjphm= (hbc[3*(i) + 2]*hbc[3*(i) + 2] + hbase) / (hbc[3*(i) + 2]);
		hjmhp=  (hbc[3*(i)]*hbc[3*(i) ] + hbase) / (hbc[3*(i)]);

        Gjphm= Gbc[3*(i) + 2];
        Gjmhp= Gbc[3*(i)];


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

    bjmh = bbc[3*i];
    bjms = bbc[3*i +1];
    bjps = bbc[3*i + 2];
    bjph = bbc[3*i + 3];


	hjphm= (hbc[3*(i) + 2]*hbc[3*(i) + 2] + hbase) / (hbc[3*(i) + 2]);
	hjmhp=  (hbc[3*(i)]*hbc[3*(i) ] + hbase) / (hbc[3*(i)]);

    Gjphm= Gbc[3*(i) + 2];
    Gjmhp= Gbc[3*(i)];

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


   
    memcpy(ubc,nGis,m*sizeof(double));

    free_dmatrix(Ared, 1, m, 1, 5);
    free_dmatrix(AL, 1, m, 1, 5);
    free(indx);
    free(nGis);

}

void getufromGsplit(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, int *RegIndx, int m1 ,int RegIndxLen, double *ubc)
{

    int hIbeg,uIbeg,uIend,bIbeg;
    // B.C's
    memcpy(ubc,uMbeg,unBC*sizeof(double));
    memcpy(ubc + unbc- unBC,uMend,unBC*sizeof(double));

    //Dry First

    hIbeg = hnBC;

    uIbeg = (unBC-1);
    uIend = (unBC-1) + 2*n;
    
    bIbeg = (bnBC-1);

    //getufromG(hbc + hIbeg, Gbc + hIbeg , bbc + bIbeg,ubc+(uIbeg+ 1) - unBC, ubc+uIend,dx , n, 2*n + 1, ubc + uIbeg);


   int i,j;
    for (i =0;i<RegIndxLen;i++)
    {
        //printf("Reg %d : %d | %d | %d | %d : %d \n",i,RegInd[i*m1],RegInd[i*m1+1],RegInd[i*m1+2],RegInd[i*m1+3],RegInd[i*m1+4]);
        if(RegIndx[i*m1 +  3] == 0)
        {
           for(j=RegIndx[i*m1]; j < RegIndx[i*m1 +  1] ; j++)
           {
                if(j == RegIndx[i*m1])
                {
                    ubc[unBC + 2*(j) - 1] = 0;
                    ubc[unBC + 2*(j)] = 0;
                    ubc[unBC + 2*(j) + 1] = 0;
                }
                else
                {
                    ubc[unBC + 2*(j)] = 0;
                    ubc[unBC + 2*(j) + 1] = 0;
                }
           }
        }
    }

    
    for (i =0;i<RegIndxLen;i++)
    {

        if(RegIndx[i*m1 +  3] == 1)
        {
            hIbeg = hnBC + 3*RegIndx[i*m1];

            uIbeg = (unBC-1) + 2*RegIndx[i*m1];
            uIend = (unBC-1) + 2*(RegIndx[i*m1 + 1]);
            
            bIbeg = (bnBC-1) + 3*RegIndx[i*m1];

            getufromG(hbc + hIbeg, Gbc + hIbeg , bbc + bIbeg,ubc+(uIbeg+ 1) - unBC, ubc+uIend,dx , RegIndx[i*m1 + 2], 2* RegIndx[i*m1 + 2] + 1, ubc + uIbeg);

        }
    }
    
}




void ReconQuart( double *q, double *qMbeg, double *qMend, int n,int nMBC, int nbcBC, int nbc, double *qbc, double dx)
{
    double qMBCmid0,qMBCmid1,qMBCmidn,qMBCmidnp1,qai,qbi,qci,qdi,qjmh,qjms,qjps,qjph;
    int i;
    double idx = 1.0 /dx;


    for(i = 0;i < nbcBC;i++)
    {
        qbc[i] = qMbeg[(nMBC- nbcBC) + i];
        qbc[nbc - nbcBC + i] = qMend[i];
    }

    // i = 0
    i = 0;
    qMBCmid0 = -qMbeg[0]*i16 + 9*qMbeg[1]*i16 + 9*qMbeg[2]*i16 - qMbeg[3]*i16;
    qMBCmid1 = -qMbeg[3]*i16 + 9*qMbeg[4]*i16 + 9*qMbeg[5]*i16 - qMbeg[6]*i16;

    qai = i12*idx*idx*idx*(-qMBCmid0 + 2*qMBCmid1 -2*q[i+1] + q[i+2]);
    qbi = i6*idx*idx*(qMBCmid0 - qMBCmid1 - q[i+1] + q[i+2]);
    qci = i12*idx*(qMBCmid0 - 8*qMBCmid1 + 8*q[i+1] - q[i+2]);
    qdi = -qMBCmid0*i6 + 2*i3*qMBCmid1 + 2*i3*q[i+1] - q[i+2]*i6;

    qjmh = -qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) - qci*(0.5*dx) + qdi;
    qjms = -qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) - qci*(dx*i6) + qdi;
    qjps = qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) + qci*(dx*i6) + qdi;
    qjph = qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) + qci*(0.5*dx) + qdi;  

    qbc[3*i + nbcBC- 1] = 0.5*(qjmh + qbc[3*i + nbcBC - 1]);
    qbc[3*i + nbcBC] = qjms;
    qbc[3*i + nbcBC +1] = qjps;
    qbc[3*i + nbcBC + 2] = qjph; 

    i = 1;
    qMBCmid0 = -qMbeg[0]*i16 + 9*qMbeg[1]*i16 + 9*qMbeg[2]*i16 - qMbeg[3]*i16;

    qai = i12*idx*idx*idx*(-qMBCmid1 + 2*q[i-1] -2*q[i+1] + q[i+2]);
    qbi = i6*idx*idx*(qMBCmid1 - q[i-1] - q[i+1] + q[i+2]);
    qci = i12*idx*(qMBCmid1 - 8*q[i-1] + 8*q[i+1] - q[i+2]);
    qdi = -qMBCmid1*i6 + 2*i3*q[i-1] + 2*i3*q[i+1] - q[i+2]*i6;

    qjmh = -qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) - qci*(0.5*dx) + qdi;
    qjms = -qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) - qci*(dx*i6) + qdi;
    qjps = qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) + qci*(dx*i6) + qdi;
    qjph = qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) + qci*(0.5*dx) + qdi;  

    qbc[3*i + nbcBC- 1] = 0.5*(qjmh + qbc[3*i + nbcBC- 1]);
    qbc[3*i + nbcBC] = qjms;
    qbc[3*i + nbcBC +1] = qjps;
    qbc[3*i + nbcBC + 2] = qjph; 


    for(i = 2;i < n -2;i++)
    {
        qai = i12*idx*idx*idx*(-q[i-2] + 2*q[i-1] -2*q[i+1] + q[i+2]);
        qbi = i6*idx*idx*(q[i-2] - q[i-1] - q[i+1] + q[i+2]);
        qci = i12*idx*(q[i-2] - 8*q[i-1] + 8*q[i+1] - q[i+2]);
        qdi = -q[i-2]*i6 + 2*i3*q[i-1] + 2*i3*q[i+1] - q[i+2]*i6;

        qjmh = -qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) - qci*(0.5*dx) + qdi;
        qjms = -qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) - qci*(dx*i6) + qdi;
        qjps = qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) + qci*(dx*i6) + qdi;
        qjph = qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) + qci*(0.5*dx) + qdi;  

        qbc[3*i + nbcBC - 1] = 0.5*(qjmh + qbc[3*i + nbcBC - 1]);
        qbc[3*i + nbcBC] = qjms;
        qbc[3*i + nbcBC +1] = qjps;
        qbc[3*i + nbcBC + 2] = qjph; 

    }

    i = n-2;
    qMBCmidn = -qMend[0]*i16 + 9*qMend[1]*i16 + 9*qMend[2]*i16 - qMend[3]*i16;
 

    qai = i12*idx*idx*idx*(-q[i-2] + 2*q[i-1] -2*q[i+1] + qMBCmidn);
    qbi = i6*idx*idx*(q[i-2] - q[i-1] - q[i+1] + qMBCmidn);
    qci = i12*idx*(q[i-2] - 8*q[i-1] + 8*q[i+1] - qMBCmidn);
    qdi = -q[i-2]*i6 + 2*i3*q[i-1] + 2*i3*q[i+1] - qMBCmidn*i6;

    qjmh = -qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) - qci*(0.5*dx) + qdi;
    qjms = -qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) - qci*(dx*i6) + qdi;
    qjps = qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) + qci*(dx*i6) + qdi;
    qjph = qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) + qci*(0.5*dx) + qdi;  

    qbc[3*i + nbcBC- 1] = 0.5*(qjmh + qbc[3*i + nbcBC - 1]);
    qbc[3*i + nbcBC] = qjms;
    qbc[3*i + nbcBC +1] = qjps;
    qbc[3*i + nbcBC + 2] = qjph; 

    i = n-1;
    qMBCmidn = -qMend[0]*i16 + 9*qMend[1]*i16 + 9*qMend[2]*i16 - qMend[3]*i16;
    qMBCmidnp1 = -qMend[3]*i16 + 9*qMend[4]*i16 + 9*qMend[5]*i16 - qMend[6]*i16;

    qai = i12*idx*idx*idx*(-q[i-2] + 2*q[i-1] -2*qMBCmidn + qMBCmidnp1);
    qbi = i6*idx*idx*(q[i-2] - q[i-1] - qMBCmidn + qMBCmidnp1);
    qci = i12*idx*(q[i-2] - 8*q[i-1] + 8*qMBCmidn - qMBCmidnp1);
    qdi = -q[i-2]*i6 + 2*i3*q[i-1] + 2*i3*qMBCmidn - qMBCmidnp1*i6;

    qjmh = -qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) - qci*(0.5*dx) + qdi;
    qjms = -qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) - qci*(dx*i6) + qdi;
    qjps = qai*(dx*i6)*(dx*i6)*(dx*i6) + qbi*(dx*i6)*(dx*i6) + qci*(dx*i6) + qdi;
    qjph = qai*(0.5*dx)*(0.5*dx)*(0.5*dx) + qbi*(0.5*dx)*(0.5*dx) + qci*(0.5*dx) + qdi;  

    qbc[3*i + nbcBC- 1] = 0.5*(qjmh + qbc[3*i + nbcBC - 1]);
    qbc[3*i + nbcBC] = qjms;
    qbc[3*i + nbcBC +1] = qjps;
    qbc[3*i + nbcBC + 2] = 0.5*(qjph + qbc[3*i + nbcBC + 2]); 



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



double ht(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return (a1*a2*(-a3 - a2*t + x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)));
}

double Gt(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return ((9*a1*a2*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/a4)) + 
      (3*a2*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) - 
      (6*pow(a1,2)*a2*a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) - 
      (9*a1*a2*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
      (a2*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))))/3. + 
   (a1*a2*a5*(-a3 - a2*t + x)*(1 - (a1*a6*a7*(-a3 - a2*t + x)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) - 
        (a6*pow(a7,2)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a7*x))/2.))/(a4*pow(E,pow(-a3 - a2*t + x,2)/a4)) + 
   (a2*a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x)*(1 - (a1*a6*a7*(-a3 - a2*t + x)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) - 
        (a6*pow(a7,2)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a7*x))/2.))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + 
   (a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*((a1*a2*a6*a7*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) - 
        (a1*a2*a6*a7*pow(-a3 - a2*t + x,2)*cos(a7*x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) - (a1*a2*a6*pow(a7,2)*(-a3 - a2*t + x)*sin(a7*x))/(2.*a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))))/
    pow(E,pow(-a3 - a2*t + x,2)/(2.*a4));
}

double Fluxhdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6,double a7, double g)
{
    return -((a1*a5*(-a3 - a2*t + x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/a4))) - (a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)));
}

double FluxGdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return (-4*pow(a5,2)*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(3.*pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
   (a1*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*g*(-a3 - a2*t + x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + 
   (2*a1*pow(a5,2)*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) + 
   (4*pow(a5,2)*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(3.*pow(a4,3)*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
   (pow(a5,2)*a6*a7*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/a4)) + 
   (2*a1*pow(a5,2)*a6*a7*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,2)*cos(a7*x))/(pow(a4,2)*pow(E,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) + 
   (2*pow(a5,2)*a6*a7*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2)*cos(a7*x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/a4)) + 
   (pow(a5,2)*a6*pow(a7,2)*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x)*sin(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
   (a5*(-a3 - a2*t + x)*(((a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) - 
           (3*a1*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
           (a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,2))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))))/3. + 
        (a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*(1 - (a1*a6*a7*(-a3 - a2*t + x)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) - 
             (a6*pow(a7,2)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a7*x))/2.))/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + 
   (a5*(((-9*a1*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
           (3*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + 
           (6*pow(a1,2)*a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) + 
           (9*a1*a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,pow(-a3 - a2*t + x,2)/a4)) + 
           (a5*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))))/3. - 
        (a1*a5*(-a3 - a2*t + x)*(1 - (a1*a6*a7*(-a3 - a2*t + x)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) - 
             (a6*pow(a7,2)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a7*x))/2.))/(a4*pow(E,pow(-a3 - a2*t + x,2)/a4)) - 
        (a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x)*(1 - (a1*a6*a7*(-a3 - a2*t + x)*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) - 
             (a6*pow(a7,2)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*sin(a7*x))/2.))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + 
        (a5*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-((a1*a6*a7*cos(a7*x))/(a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))) - (a6*pow(a7,3)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*cos(a7*x))/2. + 
             (a1*a6*a7*pow(-a3 - a2*t + x,2)*cos(a7*x))/(pow(a4,2)*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) + (3*a1*a6*pow(a7,2)*(-a3 - a2*t + x)*sin(a7*x))/(2.*a4*pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))) - 
             2*pow(a6,2)*pow(a7,3)*cos(a7*x)*sin(a7*x)))/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4))))/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4));
}

double SourceG(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return a6*a7*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*g*cos(a7*x) + (pow(a5,2)*a6*pow(a7,2)*pow(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x)*sin(a7*x))/
    (2.*a4*pow(E,pow(-a3 - a2*t + x,2)/a4)) + (pow(a5,2)*pow(a6,2)*pow(a7,3)*(a0 + a1/pow(E,pow(-a3 - a2*t + x,2)/(2.*a4)))*cos(a7*x)*sin(a7*x))/pow(E,pow(-a3 - a2*t + x,2)/a4);
}


double hFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double phi  = x - a2*t  ;
      return a0 + a1*pow(E, (-pow(phi - a3,2)/(2*a4)));
}

double hxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double phi  = x - a2*t  ;
      return -a1/a4*(phi - a3)*pow(E, (-pow(phi - a3,2)/(2*a4)));
}

double uFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double phi  = x - a2*t  ;
      return a5*pow(E, (-pow(phi - a3,2)/(2*a4)));
}

double uxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double phi  = x - a2*t  ;
      return -a5/a4*(phi - a3)*pow(E, (-pow(phi - a3,2)/(2*a4)));
}

double uxxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double phi  = x - a2*t  ;
      return -a5/(a4*a4)*(a4 - pow(phi - a3,2))*pow(E, (-pow(phi - a3,2)/(2*a4)));
}

double bFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      return a6*sin(a7*x);
}

double bxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      return a6*a7*cos(a7*x);
}

double bxxFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      return -a6*a7*a7*sin(a7*x);
}

double wFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
      double h = hFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
      double b = bFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
      return h + b;
}

double GFunc(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
     double u = uFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
     double h = hFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
     double hx = hxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
     double bx = bxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
     double bxx = bxxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
    double ux = uxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);
    double uxx = uxxFunc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g);

     return u*h*(1 + hx*bx + 0.5*h*bxx + bx*bx) - h*h*hx*ux - h*h*h/3.0*uxx;
}





//include BCs
void evolveForce(double *hbc, double *Gbc, double *wbc, double *bbc, double *ubc,double g, double dx, double dt, int n, int hnBC, int hnbc, int bnBC, int bnbc, int unBC, int unbc, double *FjhH, double *FjhG, double *SjG,double *SjHa, double *SjGa, double *x, double t, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)
{
    double idx = 1.0 / dx;  
	int i;
    double her,Ger,dber,uer,duer,hel,Gel,dbel,uel,duel,fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG,th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
	double wil,wir,wip1l,bip1l,bil,bir,nbi,hihm,hihp,uai,ubi,uaip1,ubip1,hir,hip1l,hil;
    double himhp,bedai,bedbi,bedci,bedaip1,bedbip1,bedcip1,beddi,beddip1;

    double hS, GS;



    // i = -1
    i = -1;

    uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
    ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

    uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
    ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

    bedai = 0.5*9*idx*idx*idx*(-bbc[3*i + bnBC- 1]+ 3*bbc[3*i + bnBC] - 3*bbc[3*i + bnBC +1] + bbc[3*i + bnBC + 2]);
    bedbi = 0.25*9*idx*idx*(bbc[3*i + bnBC- 1] - bbc[3*i + bnBC] - bbc[3*i + bnBC +1] + bbc[3*i + bnBC + 2]);
    bedci = i8*idx*(bbc[3*i + bnBC- 1] - 27*bbc[3*i + bnBC] + 27*bbc[3*i + bnBC +1]- bbc[3*i + bnBC + 2]);
    beddi = -bbc[3*i + bnBC- 1]*i16 + 9*bbc[3*i + bnBC]*i16 + 9*bbc[3*i + bnBC +1]*i16 - bbc[3*i + bnBC + 2]*i16;

    bedaip1 = 0.5*9*idx*idx*idx*(-bbc[3*(i+1) + bnBC -1]+ 3*bbc[3*(i+1) + bnBC] - 3*bbc[3*(i+1) + bnBC +1] + bbc[3*(i+1) + bnBC + 2]);
    bedbip1 = 0.25*9*idx*idx*(bbc[3*(i+1) + bnBC- 1] - bbc[3*(i+1) + bnBC] - bbc[3*(i+1) + bnBC+1] + bbc[3*(i+1) + bnBC + 2]);
    bedcip1 = i8*idx*(bbc[3*(i+1) + bnBC- 1] - 27*bbc[3*(i+1) + bnBC] + 27*bbc[3*(i+1) + bnBC +1]- bbc[3*(i+1) + bnBC + 2]);
    beddip1 = -bbc[3*(i+1) + bnBC- 1]*i16 + 9*bbc[3*(i+1) + bnBC]*i16 + 9*bbc[3*(i+1) + bnBC +1]*i16 - bbc[3*(i+1) + bnBC + 2]*i16;


    her = hbc[hnBC + 3*(i+1)] ;
    Ger = Gbc[hnBC + 3*(i+1)];
    dber = 3*bedaip1*(0.5*dx)*(0.5*dx) - 2*bedbip1*(0.5*dx) + bedcip1;
    uer  = ubc[unBC + 2*(i) + 1];
    duer = -uaip1*(dx) + ubip1;


    hel = hbc[hnBC + 3*(i) + 2] ;
    Gel = Gbc[hnBC + 3*(i) + 2];
    dbel = 3*bedai*(0.5*dx)*(0.5*dx) + 2*bedbi*(0.5*dx) + bedci;
    uel  = ubc[unBC + 2*(i) + 1];
    duel = uai*(dx) + ubi;

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

    for(i = 0;i < n;i++)
    {
        wir = wbc[hnBC + 3*(i) + 2];
        wil = wbc[hnBC + 3*(i)];

        hir = hbc[hnBC + 3*(i) + 2];
        hil = hbc[hnBC + 3*(i)];
   

        bir = wir - hir;
        bil = wil - hil;     
        
        uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
        ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

        uaip1 =2*idx*idx*(ubc[2*(i+1) + unBC - 1] - 2*ubc[2*(i+1) + unBC] + ubc[2*(i+1) + unBC + 1]);
        ubip1 =idx*(-ubc[2*(i+1) + unBC - 1]+ ubc[2*(i+1) + unBC + 1]);

        bedai = 0.5*9*idx*idx*idx*(-bbc[3*i + bnBC- 1]+ 3*bbc[3*i + bnBC] - 3*bbc[3*i + bnBC +1] + bbc[3*i + bnBC + 2]);
        bedbi = 0.25*9*idx*idx*(bbc[3*i + bnBC- 1] - bbc[3*i + bnBC] - bbc[3*i + bnBC +1] + bbc[3*i + bnBC + 2]);
        bedci = i8*idx*(bbc[3*i + bnBC- 1] - 27*bbc[3*i + bnBC] + 27*bbc[3*i + bnBC +1]- bbc[3*i + bnBC + 2]);
        beddi = -bbc[3*i + bnBC- 1]*i16 + 9*bbc[3*i + bnBC]*i16 + 9*bbc[3*i + bnBC +1]*i16 - bbc[3*i + bnBC + 2]*i16;

        bedaip1 = 0.5*9*idx*idx*idx*(-bbc[3*(i+1) + bnBC -1]+ 3*bbc[3*(i+1) + bnBC] - 3*bbc[3*(i+1) + bnBC +1] + bbc[3*(i+1) + bnBC + 2]);
        bedbip1 = 0.25*9*idx*idx*(bbc[3*(i+1) + bnBC- 1] - bbc[3*(i+1) + bnBC] - bbc[3*(i+1) + bnBC+1] + bbc[3*(i+1) + bnBC + 2]);
        bedcip1 = i8*idx*(bbc[3*(i+1) + bnBC- 1] - 27*bbc[3*(i+1) + bnBC] + 27*bbc[3*(i+1) + bnBC +1]- bbc[3*(i+1) + bnBC + 2]);
        beddip1 = -bbc[3*(i+1) + bnBC- 1]*i16 + 9*bbc[3*(i+1) + bnBC]*i16 + 9*bbc[3*(i+1) + bnBC +1]*i16 - bbc[3*(i+1) + bnBC + 2]*i16;


        her = hbc[hnBC + 3*(i+1)] ;
        Ger = Gbc[hnBC + 3*(i+1)];
        dber = 3*bedaip1*(0.5*dx)*(0.5*dx) - 2*bedbip1*(0.5*dx) + bedcip1;
        uer  = ubc[unBC + 2*(i) + 1];
        duer = -uaip1*(dx) + ubip1;;


        hel = hbc[hnBC + 3*(i) + 2] ;
        Gel = Gbc[hnBC + 3*(i) + 2];
        dbel = 3*bedai*(0.5*dx)*(0.5*dx) + 2*bedbi*(0.5*dx) + bedci;
        uel  = ubc[unBC + 2*(i) + 1];
        duel = uai*(dx) + ubi;

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
		tu =  ubc[2*i + unBC];
		tux = ubi;
		tbx = bedci;
		tbxx =2*bedbi;
		
		//sourcer = g*0.5*(hihm*hihm - hir*hir);
		sourcec = -g*th*tbx -  0.5*th*th*tu*tux*tbxx + th*tu*tu*tbx*tbxx ;
		//sourcel = g*0.5*(hil*hil - himhp*himhp);

        FjhH[i + 1] = foh;
        FjhG[i + 1] = foG;
        SjG[i] = dx*sourcec;

        hS = ht(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g);
        GS = Gt(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g);

        SjHa[i] = hS;
        SjGa[i] = GS;


    }

}

void Fluxsplit(double *hbc, double *Gbc, double *wbc, double *bbc, double *ubc,double g, double dx, double dt, int n, int hnBC, int hnbc, int bnBC, int bnbc, int unBC, int unbc, int *RegIndx, int m1 ,int RegIndxLen, double *FjhH, double *FjhG, double *SjG,double *SjHa, double *SjGa, double *x, double t, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)
{

    int hIbeg,uIbeg,uIend,bIbeg,FIbeg,SIbeg;
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

            hS = ht(x[RegIndx[i*m1]],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + Fluxhdiff(x[RegIndx[i*m1]],t,a0,a1,a2,a3,a4,a5,a6,a7,g);
            GS = Gt(x[RegIndx[i*m1]],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + FluxGdiff(x[RegIndx[i*m1]],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + SourceG(x[RegIndx[i*m1]],t,a0,a1,a2,a3,a4,a5,a6,a7,g);

            SjHa[RegIndx[i*m1]] = hS;
            SjGa[RegIndx[i*m1]] = GS;

           for(j=RegIndx[i*m1] + 1; j < RegIndx[i*m1 +  1]; j++)
           {
                FjhH[j] = 0;
                FjhG[j] = 0;
                SjG[j] = 0;
                SjHa[j] = 0;
                SjGa[j] = 0;

           }

            hS = ht(x[RegIndx[i*m1 + 1] - 1],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + Fluxhdiff(x[RegIndx[i*m1 + 1] - 1],t,a0,a1,a2,a3,a4,a5,a6,a7,g);
            GS = Gt(x[RegIndx[i*m1 + 1] - 1],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + FluxGdiff(x[RegIndx[i*m1 + 1] - 1],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + SourceG(x[RegIndx[i*m1 + 1] - 1],t,a0,a1,a2,a3,a4,a5,a6,a7,g);

            SjHa[RegIndx[i*m1 + 1] - 1] = hS;
            SjGa[RegIndx[i*m1 + 1] - 1] = GS;
        }

        if(RegIndx[i*m1 +  3] == 1)
        {
            //printf("%d : %d | %d | %d | %d  \n",i,RegIndx[i*m1],RegIndx[i*m1 +  1],RegIndx[i*m1 +  2],RegIndx[i*m1 +  3]    );
            hIbeg = hnBC + 3*(RegIndx[i*m1] - 1);
            uIbeg = (unBC-1) + 2*(RegIndx[i*m1] - 1);
            bIbeg = (bnBC-1) + 3*(RegIndx[i*m1] - 1);
            FIbeg = RegIndx[i*m1];
            SIbeg = RegIndx[i*m1];

           // printf("%d | %d | %d | %d | %d  \n \n",hIbeg,uIbeg,bIbeg,FIbeg,SIbeg   );



            evolveForce(hbc+ hIbeg,Gbc+ hIbeg,wbc+ hIbeg,bbc + bIbeg,ubc+ uIbeg,g,dx,dt,RegIndx[i*m1 +  2],hnBC,hnbc,bnBC,bnbc,unBC,unbc,FjhH + FIbeg,FjhG + FIbeg ,SjG + SIbeg,SjHa+ SIbeg,SjGa+ SIbeg,x+ SIbeg,t,a0,a1,a2,a3,a4,a5,a6,a7);
        }
    }

}

void FluxEvo(double *hbc, double *Gbc,double g, double dx, double dt, int n, int hnBC, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double *FjhH, double *FjhG, double *SjG, double *hp, double *Gp,double *SjHa, double *SjGa)
{
    double hS,GS;
    int i;
    double idx = 1.0 / dx;
    for (i = 0; i < n; i++)
    {

        hp[i] = hbc[hnBC + 3*(i) + 1] - dt*idx*(FjhH[i+1] - FjhH[i])  + dt*SjHa[i];
        Gp[i] = Gbc[hnBC + 3*(i) + 1] - dt*idx*(FjhG[i+1] - FjhG[i]) + dt*idx*(SjG[i]) + dt*SjGa[i];
        //Gp[i] = Gbc[hnBC + 3*(i) + 1] + dt*idx*(SjG[i]) + dt*SjGa[i];

        //printf("h Flux  : %d | %e | %e | %e  \n ",i,- dt*idx*(FjhH[i+1] - FjhH[i]),dt*Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), - dt*idx*(FjhH[i+1] - FjhH[i]) + dt*Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g));
        //printf("G Flux  : %d | %e | %e | %e  \n ",i,- dt*idx*(FjhG[i+1] - FjhG[i]),dt*FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), - dt*idx*(FjhG[i+1] - FjhG[i]) + dt*FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g));
        //printf("G Source: %d | %e | %e | %e  \n ",i,dt*idx*(SjG[i]),dt*SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g), dt*idx*(SjG[i])+ dt*SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,g) );

    }

}

void FluxWrap(double *hbc, double *Gbc, double *wbc, double *bbc, double *ubc,double g, double dx, double dt, int n, int hnBC, int hnbc, int bnBC, int bnbc, int unBC, int unbc, int *RegIndx, int m1 ,int RegIndxLen, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double *hp, double *Gp)
{
    double *FjhH = malloc((n+1)*sizeof(double));
    double *FjhG = malloc((n+1)*sizeof(double));
    double *SjG = malloc((n)*sizeof(double));

    double *SjHa = malloc((n)*sizeof(double));
    double *SjGa = malloc((n)*sizeof(double));
    
    Fluxsplit(hbc,Gbc,wbc,bbc,ubc,g,dx,dt,n,hnBC,hnbc,bnBC,bnbc,unBC,unbc,RegIndx,m1,RegIndxLen,FjhH,FjhG,SjG,SjHa,SjGa,x,t,a0,a1,a2,a3,a4,a5,a6,a7);

    FluxEvo(hbc,Gbc,g,dx,dt,n,hnBC,x,t,a0,a1,a2,a3,a4,a5,a6,a7,FjhH,FjhG,SjG,hp,Gp,SjHa,SjGa);

    free(FjhH);
    free(FjhG);
    free(SjG);

    free(SjHa);
    free(SjGa);
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

void ReconWrap(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *w, double *hbc, double *Gbc, double *wbc, double *bbc,int *RegInd ,int m1,int RegIndLen)
{

    // calculate stage
    Stagehbed(n,h,b,w,RegInd,m1,RegIndLen);

    //Reconstruct
    ReconLin(h, hMbeg, hMend,n,hnBC,hnbc,theta,hbc);
    ReconLin(G, GMbeg, GMend,n,hnBC,hnbc,theta,Gbc);
    ReconLin(w, wMbeg, wMend,n,hnBC,hnbc,theta,wbc);

    ReconZero(hbc,n,hnBC,hnbc,RegInd,m1,RegIndLen);
    ReconZero(Gbc,n,hnBC,hnbc,RegInd,m1,RegIndLen);

    ReconQuart(b,bMbeg,bMend,n,bnMBC, bnBC, bnbc,bbc, dx);

}

void ReconandSolve(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *Gbc,double *hbc,double *wbc,double *ubc,double *bbc)
{

    double *w = malloc((n)*sizeof(double));
    int u_length = 2*n + 1;
    int *RegInd = RegSplit(h,n); 
    int m1 = 5;
    int RegIndLen = RegInd[m1-1];

    ReconWrap(h,G,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,bMbeg,bMend,uMbeg,uMend,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,w,hbc,Gbc,wbc,bbc,RegInd,m1,RegIndLen);
    getufromGsplit(hbc, Gbc, bbc,uMbeg,uMend,dx ,n,u_length ,hnBC,bnBC,hnbc,bnbc,unBC,unbc,RegInd,m1,RegIndLen,ubc);

    free(RegInd);
    free(w);
}

void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7)
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

    double *ubc = malloc((unbc)*sizeof(double));

    double *Gpbc = malloc((hnbc)*sizeof(double));
    double *hpbc = malloc((hnbc)*sizeof(double));
    double *wpbc = malloc((hnbc)*sizeof(double));

    double *upbc = malloc((unbc)*sizeof(double));

    double *bbc = malloc((bnbc)*sizeof(double));


    //Dry Regions
    int *RegInd = RegSplit(h,n); 
    int m1 = 5;
    int RegIndLen = RegInd[m1-1];

    ReconWrap(h,G,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,bMbeg,bMend,uMbeg,uMend,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,w,hbc,Gbc,wbc,bbc,RegInd,m1,RegIndLen);

    getufromGsplit(hbc, Gbc, bbc,uMbeg,uMend,dx ,n,u_length ,hnBC,bnBC,hnbc,bnbc,unBC,unbc,RegInd,m1,RegIndLen,ubc);

    FluxWrap(hbc,Gbc,wbc,bbc,ubc,g,dx,dt,n,hnBC,hnbc,bnBC,bnbc,unBC,unbc,RegInd,m1 ,RegIndLen,x,t,a0,a1,a2,a3,a4,a5,a6,a7,hp,Gp);

 
    int *RegIndp  = RegSplit(hp,n); 
    m1 = 5;
    RegIndLen = RegIndp[m1-1];

    ReconWrap(hp,Gp,b,hMbeg,hMend,GMbeg,GMend,wMbeg,wMend,bMbeg,bMend,uMbeg,uMend,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,wp,hpbc,Gpbc,wpbc,bbc,RegIndp,m1,RegIndLen);

    getufromGsplit(hpbc, Gpbc, bbc,uMbeg,uMend,dx ,n,u_length ,hnBC,bnBC,hnbc,bnbc,unBC,unbc,RegIndp,m1,RegIndLen,upbc);

    FluxWrap(hpbc,Gpbc,wpbc,bbc,upbc,g,dx,dt,n,hnBC,hnbc,bnBC,bnbc,unBC,unbc,RegIndp,m1 ,RegIndLen,x,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,hpp,Gpp);

    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gpp[i]);  
        h[i] = 0.5*(h[i] + hpp[i]);
        
    }

    //Reg2 = RegSplit(h,n);
    //Solve for u given wet/dry regions
    //getufromGsplit(h,G,bed,hMbeg,hMend,GMbeg,GMend,uMbeg,uMend,wMbeg,wMend,bMbeg,bMend,theta,dx ,n,u_length ,nGhBC,unBC,bnBC,nGhhbc,nubc,nbhc,ubc,hhbc,Ghbc,whbc,bedhbc,Reg2);

    free(w);
    free(wp);
    free(hp);
    free(Gp);
    free(Gpp);
    free(hpp);

    free(Gbc);
    free(hbc);
    free(wbc);
    free(bbc);
    free(ubc);

    free(Gpbc);
    free(hpbc);
    free(wpbc);
    free(upbc);

    free(RegInd);
    //free(RegIndp);
}


int main()
{
    printf("h");
    return 1;
}
