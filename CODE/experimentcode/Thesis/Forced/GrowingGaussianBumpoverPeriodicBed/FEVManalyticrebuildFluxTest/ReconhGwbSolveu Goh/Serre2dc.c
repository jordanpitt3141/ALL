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

void getufromG(double *hbc, double *Gbc, double *bbc, double *uMbeg, double *uMend, double dx , int n, int m, int hnBC, int bnBC, int hnbc, int bnbc, int unBC, int unbc, double *ubc)
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

    double Gva11,Gva21,Gva31,uva11,uva12,uva13,uva21,uva22,uva23,uva31,uva32,uva33;
    double h2uxvxa11,h2uxvxa12,h2uxvxa13,h2uxvxa21,h2uxvxa22,h2uxvxa23,h2uxvxa31,h2uxvxa32,h2uxvxa33;
    double LHSa11,LHSa12,LHSa13,LHSa21,LHSa22,LHSa23,LHSa31,LHSa32,LHSa33;
    double Gohjphm,Gohjmhp;

    double wi,wip1,wim1,dwib,dwim,dwif,dwi,wjphm,wjmhp, bedbegmiddle,bedai,bedbi,bedci,beddi,bedendmiddle,bjmh,bjms,bjps,bjph;
    double hbxuvxa11,hbxuvxa12,hbxuvxa13,hbxuvxa21,hbxuvxa22,hbxuvxa23,hbxuvxa31,hbxuvxa32,hbxuvxa33;
    double hbxuxva11,hbxuxva12,hbxuxva13,hbxuxva21,hbxuxva22,hbxuxva23,hbxuxva31,hbxuxva32,hbxuxva33;
    double bx2uva11,bx2uva12,bx2uva13,bx2uva21,bx2uva22,bx2uva23,bx2uva31,bx2uva32,bx2uva33;

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

    bjmh = bbc[3*i + bnBC- 1];
    bjms = bbc[3*i + bnBC];
    bjps = bbc[3*i + bnBC + 1];
    bjph = bbc[3*i + bnBC + 2];


    hjphm=  hbc[hnBC + 3*(i) + 2]; 
    hjmhp=  hbc[hnBC + 3*(i)];  

    Gjphm= Gbc[hnBC + 3*(i) + 2];
    Gjmhp= Gbc[hnBC + 3*(i)];

    Gohjphm= Gjphm / hjphm;
    Gohjmhp= Gjmhp /  hjmhp;


   // G integral (RHS)
    Gva11 = dx*i6*(Gohjmhp);
    Gva21 = dx*i6*(2*Gohjmhp + 2*Gohjphm);
    Gva31 = dx*i6*(Gohjphm);

           
    //u integral
    uva11 = 0.5*dx*((4.0/15));
    uva12 = 0.5*dx*((2.0/15));
    uva13 = 0.5*dx*((-1.0/15));

    uva21 = 0.5*dx*((2.0/15));
    uva22 = 0.5*dx*((16.0/15));
    uva23 = 0.5*dx*((2.0/15));

    uva31 = 0.5*dx*((-1.0/15));
    uva32 = 0.5*dx*((2.0/15));
    uva33 = 0.5*dx*((4.0/15));
    
    //h2ux
    
    h2uxvxa11 = 2*i3*idx*((23.0/30)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
    h2uxvxa12 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
    h2uxvxa13 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);

    h2uxvxa21 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
    h2uxvxa22 = 2*i3*idx*((16.0/15)*hjmhp*hjmhp+(4.0/15)*hjphm*hjmhp+(4.0/15)*hjmhp*hjphm+(16.0/15)*hjphm*hjphm);
    h2uxvxa23 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);

    h2uxvxa31 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
    h2uxvxa32 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);
    h2uxvxa33 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(23.0/30)*hjphm*hjphm);

    //hbxuvx
    hbxuvxa11 = -2*0.5*idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
    hbxuvxa12 = -2*0.5*idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
    hbxuvxa13 = -2*0.5*idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);

    hbxuvxa21 = -2*0.5*idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
    hbxuvxa22 = -2*0.5*idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
    hbxuvxa23 = -2*0.5*idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);

    hbxuvxa31 = -2*0.5*idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);
    hbxuvxa32 = -2*0.5*idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);
    hbxuvxa33 = -2*0.5*idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

    //hbxuxv
    hbxuxva11 = -2*0.5*idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
    hbxuxva12 = -2*0.5*idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
    hbxuxva13 = -2*0.5*idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);

    hbxuxva21 = -2*0.5*idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
    hbxuxva22 = -2*0.5*idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
    hbxuxva23 = -2*0.5*idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);

    hbxuxva31 = -2*0.5*idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);
    hbxuxva32 = -2*0.5*idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);
    hbxuxva33 = -2*0.5*idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

    //bx2uv

    bx2uva11 = 2*idx*((1777.0/1680)*bjmh*bjmh+(-207.0/140)*bjms*bjmh+(297.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-207.0/140)*bjmh*bjms+(243.0/112)*bjms*bjms+(-243.0/280)*bjps*bjms+(99.0/560)*bjph*bjms+(297.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/560)*bjps*bjps+(-27.0/280)*bjph*bjps+(-23.0/210)*bjmh*bjph+(99.0/560)*bjms*bjph+(-27.0/280)*bjps*bjph+(7.0/240)*bjph*bjph);

    bx2uva12 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);

    bx2uva13 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);

    bx2uva21 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);

    bx2uva22 = 2*idx*((13.0/42)*bjmh*bjmh+(-9.0/35)*bjms*bjmh+(-9.0/70)*bjps*bjmh+(8.0/105)*bjph*bjmh+(-9.0/35)*bjmh*bjms+(27.0/14)*bjms*bjms+(-54.0/35)*bjps*bjms+(-9.0/70)*bjph*bjms+(-9.0/70)*bjmh*bjps+(-54.0/35)*bjms*bjps+(27.0/14)*bjps*bjps+(-9.0/35)*bjph*bjps+(8.0/105)*bjmh*bjph+(-9.0/70)*bjms*bjph+(-9.0/35)*bjps*bjph+(13.0/42)*bjph*bjph);

    bx2uva23 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);

    bx2uva31 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);

    bx2uva32 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);

    bx2uva33 = 2*idx*((7.0/240)*bjmh*bjmh+(-27.0/280)*bjms*bjmh+(99.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-27.0/280)*bjmh*bjms+(243.0/560)*bjms*bjms+(-243.0/280)*bjps*bjms+(297.0/560)*bjph*bjms+(99.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/112)*bjps*bjps+(-207.0/140)*bjph*bjps+(-23.0/210)*bjmh*bjph+(297.0/560)*bjms*bjph+(-207.0/140)*bjps*bjph+(1777.0/1680)*bjph*bjph);
    
    // LHS 
    
    LHSa11 = uva11 + h2uxvxa11 + hbxuvxa11 + hbxuxva11 + bx2uva11;
    LHSa12 = uva12 + h2uxvxa12 + hbxuvxa12 + hbxuxva12 + bx2uva12;
    LHSa13 = uva13 + h2uxvxa13 + hbxuvxa13 + hbxuxva13 + bx2uva13;

    LHSa21 = uva21 + h2uxvxa21 + hbxuvxa21 + hbxuxva21 + bx2uva21;
    LHSa22 = uva22 + h2uxvxa22 + hbxuvxa22 + hbxuxva22 + bx2uva22;
    LHSa23 = uva23 + h2uxvxa23 + hbxuvxa23 + hbxuxva23 + bx2uva23;

    LHSa31 = uva31 + h2uxvxa31 + hbxuvxa31 + hbxuxva31 + bx2uva31;
    LHSa32 = uva32 + h2uxvxa32 + hbxuvxa32 + hbxuxva32 + bx2uva32;
    LHSa33 = uva33 + h2uxvxa33 + hbxuvxa33 + hbxuxva33 + bx2uva33;


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
    nGis[j] = nGis[j] + Gva21; 
    nGis[j+1] = nGis[j+1] + Gva31;     



    j = 3;
    for (i =1;i < n-1 ; i++)
    {
        bjmh = bbc[3*i + bnBC- 1];
        bjms = bbc[3*i + bnBC];
        bjps = bbc[3*i + bnBC + 1];
        bjph = bbc[3*i + bnBC + 2];


        hjphm=  hbc[hnBC + 3*(i) + 2]; 
        hjmhp=  hbc[hnBC + 3*(i)];  

        Gjphm= Gbc[hnBC + 3*(i) + 2];
        Gjmhp= Gbc[hnBC + 3*(i)];

        Gohjphm= Gjphm / hjphm;
        Gohjmhp= Gjmhp /  hjmhp;


       // G integral (RHS)
        Gva11 = dx*i6*(Gohjmhp);
        Gva21 = dx*i6*(2*Gohjmhp + 2*Gohjphm);
        Gva31 = dx*i6*(Gohjphm);

               
        //u integral

        uva11 = 0.5*dx*((4.0/15));
        uva12 = 0.5*dx*((2.0/15));
        uva13 = 0.5*dx*((-1.0/15));

        uva21 = 0.5*dx*((2.0/15));
        uva22 = 0.5*dx*((16.0/15));
        uva23 = 0.5*dx*((2.0/15));

        uva31 = 0.5*dx*((-1.0/15));
        uva32 = 0.5*dx*((2.0/15));
        uva33 = 0.5*dx*((4.0/15));
            
        //h2ux
        
        h2uxvxa11 = 2*i3*idx*((23.0/30)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
        h2uxvxa12 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
        h2uxvxa13 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);

        h2uxvxa21 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
        h2uxvxa22 = 2*i3*idx*((16.0/15)*hjmhp*hjmhp+(4.0/15)*hjphm*hjmhp+(4.0/15)*hjmhp*hjphm+(16.0/15)*hjphm*hjphm);
        h2uxvxa23 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);

        h2uxvxa31 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
        h2uxvxa32 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);
        h2uxvxa33 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(23.0/30)*hjphm*hjphm);

        //hbxuvx
        hbxuvxa11 = -idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
        hbxuvxa12 = -idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
        hbxuvxa13 = -idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);

        hbxuvxa21 = -idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
        hbxuvxa22 = -idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
        hbxuvxa23 = -idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);

        hbxuvxa31 = -idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);
        hbxuvxa32 = -idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);
        hbxuvxa33 = -idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

        //hbxuxv
        hbxuxva11 = -idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
        hbxuxva12 = -idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
        hbxuxva13 = -idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);

        hbxuxva21 = -idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
        hbxuxva22 = -idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
        hbxuxva23 = -idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);

        hbxuxva31 = -idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);
        hbxuxva32 = -idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);
        hbxuxva33 = -idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

        //bx2uv

        bx2uva11 = 2*idx*((1777.0/1680)*bjmh*bjmh+(-207.0/140)*bjms*bjmh+(297.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-207.0/140)*bjmh*bjms+(243.0/112)*bjms*bjms+(-243.0/280)*bjps*bjms+(99.0/560)*bjph*bjms+(297.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/560)*bjps*bjps+(-27.0/280)*bjph*bjps+(-23.0/210)*bjmh*bjph+(99.0/560)*bjms*bjph+(-27.0/280)*bjps*bjph+(7.0/240)*bjph*bjph);
        bx2uva12 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);
        bx2uva13 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);

        bx2uva21 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);
        bx2uva22 = 2*idx*((13.0/42)*bjmh*bjmh+(-9.0/35)*bjms*bjmh+(-9.0/70)*bjps*bjmh+(8.0/105)*bjph*bjmh+(-9.0/35)*bjmh*bjms+(27.0/14)*bjms*bjms+(-54.0/35)*bjps*bjms+(-9.0/70)*bjph*bjms+(-9.0/70)*bjmh*bjps+(-54.0/35)*bjms*bjps+(27.0/14)*bjps*bjps+(-9.0/35)*bjph*bjps+(8.0/105)*bjmh*bjph+(-9.0/70)*bjms*bjph+(-9.0/35)*bjps*bjph+(13.0/42)*bjph*bjph);
        bx2uva23 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);

        bx2uva31 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);
        bx2uva32 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);
        bx2uva33 = 2*idx*((7.0/240)*bjmh*bjmh+(-27.0/280)*bjms*bjmh+(99.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-27.0/280)*bjmh*bjms+(243.0/560)*bjms*bjms+(-243.0/280)*bjps*bjms+(297.0/560)*bjph*bjms+(99.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/112)*bjps*bjps+(-207.0/140)*bjph*bjps+(-23.0/210)*bjmh*bjph+(297.0/560)*bjms*bjph+(-207.0/140)*bjps*bjph+(1777.0/1680)*bjph*bjph);
        
        // LHS 
        
        LHSa11 = uva11 + h2uxvxa11 + hbxuvxa11 + hbxuxva11 + bx2uva11;
        LHSa12 = uva12 + h2uxvxa12 + hbxuvxa12 + hbxuxva12 + bx2uva12;
        LHSa13 = uva13 + h2uxvxa13 + hbxuvxa13 + hbxuxva13 + bx2uva13;

        LHSa21 = uva21 + h2uxvxa21 + hbxuvxa21 + hbxuxva21 + bx2uva21;
        LHSa22 = uva22 + h2uxvxa22 + hbxuvxa22 + hbxuxva22 + bx2uva22;
        LHSa23 = uva23 + h2uxvxa23 + hbxuvxa23 + hbxuxva23 + bx2uva23;

        LHSa31 = uva31 + h2uxvxa31 + hbxuvxa31 + hbxuxva31 + bx2uva31;
        LHSa32 = uva32 + h2uxvxa32 + hbxuvxa32 + hbxuxva32 + bx2uva32;
        LHSa33 = uva33 + h2uxvxa33 + hbxuvxa33 + hbxuxva33 + bx2uva33;

        Ared[j-1 + 3][1] = Ared[j-1 + 3][1] + LHSa31;

        Ared[j-1 +2][2] = Ared[j-1 + 2][2] + LHSa21;
        Ared[j + 2][2] = Ared[j + 2][2] + LHSa32;

        Ared[j-1 + 1][3] = Ared[j-1 + 1][3] + LHSa11;
        Ared[j + 1][3] = Ared[j + 1][3] + LHSa22;
        Ared[j+1 + 1][3] = Ared[j+1 + 1][3] + LHSa33;

        Ared[j-1 + 1][4] = Ared[j-1 + 1][4] + LHSa12;
        Ared[j + 1][4] = Ared[j + 1][4] + LHSa23;
        
        Ared[j-1 + 1 ][5] = Ared[j-1 + 1][5] + LHSa13;


        nGis[j-1] = nGis[j-1] + Gva11;
        nGis[j] = nGis[j] + Gva21; 
        nGis[j+1] = nGis[j+1] + Gva31;

      
        
        j = j + 2 ;       


    }


  

// last
    j = m-2;
    i = n-1;

   // Get it from interior

    // reconstruct bed

    bjmh = bbc[3*i + bnBC- 1];
    bjms = bbc[3*i + bnBC];
    bjps = bbc[3*i + bnBC + 1];
    bjph = bbc[3*i + bnBC + 2];


    hjphm=  hbc[hnBC + 3*(i) + 2]; 
    hjmhp=  hbc[hnBC + 3*(i)];  

    Gjphm= Gbc[hnBC + 3*(i) + 2];
    Gjmhp= Gbc[hnBC + 3*(i)];

    Gohjphm= Gjphm / hjphm;
    Gohjmhp= Gjmhp /  hjmhp;


    // G integral (RHS)
    Gva11 = dx*i6*(Gohjmhp);
    Gva21 = dx*i6*(2*Gohjmhp + 2*Gohjphm);
    Gva31 = dx*i6*(Gohjphm);

           
    //u integral

    uva11 = 0.5*dx*((4.0/15));
    uva12 = 0.5*dx*((2.0/15));
    uva13 = 0.5*dx*((-1.0/15));

    uva21 = 0.5*dx*((2.0/15));
    uva22 = 0.5*dx*((16.0/15));
    uva23 = 0.5*dx*((2.0/15));

    uva31 = 0.5*dx*((-1.0/15));
    uva32 = 0.5*dx*((2.0/15));
    uva33 = 0.5*dx*((4.0/15));
        
    //h2ux
    
    h2uxvxa11 = 2*i3*idx*((23.0/30)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
    h2uxvxa12 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
    h2uxvxa13 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);

    h2uxvxa21 = 2*i3*idx*((-13.0/15)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-1.0/5)*hjphm*hjphm);
    h2uxvxa22 = 2*i3*idx*((16.0/15)*hjmhp*hjmhp+(4.0/15)*hjphm*hjmhp+(4.0/15)*hjmhp*hjphm+(16.0/15)*hjphm*hjphm);
    h2uxvxa23 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);

    h2uxvxa31 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(-1.0/60)*hjphm*hjmhp+(-1.0/60)*hjmhp*hjphm+(1.0/10)*hjphm*hjphm);
    h2uxvxa32 = 2*i3*idx*((-1.0/5)*hjmhp*hjmhp+(-2.0/15)*hjphm*hjmhp+(-2.0/15)*hjmhp*hjphm+(-13.0/15)*hjphm*hjphm);
    h2uxvxa33 = 2*i3*idx*((1.0/10)*hjmhp*hjmhp+(3.0/20)*hjphm*hjmhp+(3.0/20)*hjmhp*hjphm+(23.0/30)*hjphm*hjphm);

    //hbxuvx
    hbxuvxa11 = -idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
    hbxuvxa12 = -idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
    hbxuvxa13 = -idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);

    hbxuvxa21 = -idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
    hbxuvxa22 = -idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
    hbxuvxa23 = -idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);

    hbxuvxa31 = -idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);
    hbxuvxa32 = -idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);
    hbxuvxa33 = -idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

    //hbxuxv
    hbxuxva11 = -idx*((677.0/840)*hjmhp*bjmh+(137.0/1680)*hjphm*bjmh+(-57.0/56)*hjmhp*bjms+(-39.0/560)*hjphm*bjms+(15.0/56)*hjmhp*bjps+(-3.0/560)*hjphm*bjps+(-47.0/840)*hjmhp*bjph+(-11.0/1680)*hjphm*bjph);
    hbxuxva12 = -idx*((-209.0/210)*hjmhp*bjmh+(-37.0/420)*hjphm*bjmh+(9.0/7)*hjmhp*bjms+(9.0/140)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(2.0/21)*hjmhp*bjph+(37.0/420)*hjphm*bjph);
    hbxuxva13 = -idx*((53.0/280)*hjmhp*bjmh+(11.0/1680)*hjphm*bjmh+(-15.0/56)*hjmhp*bjms+(3.0/560)*hjphm*bjms+(33.0/280)*hjmhp*bjps+(39.0/560)*hjphm*bjps+(-11.0/280)*hjmhp*bjph+(-137.0/1680)*hjphm*bjph);

    hbxuxva21 = -idx*((17.0/42)*hjmhp*bjmh+(11.0/140)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(3.0/14)*hjphm*bjms+(-27.0/70)*hjmhp*bjps+(-51.0/140)*hjphm*bjps+(19.0/420)*hjmhp*bjph+(1.0/14)*hjphm*bjph);
    hbxuxva22 = -idx*((-10.0/21)*hjmhp*bjmh+(-13.0/105)*hjphm*bjmh+(3.0/7)*hjmhp*bjms+(6.0/35)*hjphm*bjms+(6.0/35)*hjmhp*bjps+(3.0/7)*hjphm*bjps+(-13.0/105)*hjmhp*bjph+(-10.0/21)*hjphm*bjph);
    hbxuxva23 = -idx*((1.0/14)*hjmhp*bjmh+(19.0/420)*hjphm*bjmh+(-51.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(3.0/14)*hjmhp*bjps+(-9.0/140)*hjphm*bjps+(11.0/140)*hjmhp*bjph+(17.0/42)*hjphm*bjph);

    hbxuxva31 = -idx*((-137.0/1680)*hjmhp*bjmh+(-11.0/280)*hjphm*bjmh+(39.0/560)*hjmhp*bjms+(33.0/280)*hjphm*bjms+(3.0/560)*hjmhp*bjps+(-15.0/56)*hjphm*bjps+(11.0/1680)*hjmhp*bjph+(53.0/280)*hjphm*bjph);
    hbxuxva32 = -idx*((37.0/420)*hjmhp*bjmh+(2.0/21)*hjphm*bjmh+(-9.0/140)*hjmhp*bjms+(-27.0/70)*hjphm*bjms+(9.0/140)*hjmhp*bjps+(9.0/7)*hjphm*bjps+(-37.0/420)*hjmhp*bjph+(-209.0/210)*hjphm*bjph);
    hbxuxva33 = -idx*((-11.0/1680)*hjmhp*bjmh+(-47.0/840)*hjphm*bjmh+(-3.0/560)*hjmhp*bjms+(15.0/56)*hjphm*bjms+(-39.0/560)*hjmhp*bjps+(-57.0/56)*hjphm*bjps+(137.0/1680)*hjmhp*bjph+(677.0/840)*hjphm*bjph);

    //bx2uv

    bx2uva11 = 2*idx*((1777.0/1680)*bjmh*bjmh+(-207.0/140)*bjms*bjmh+(297.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-207.0/140)*bjmh*bjms+(243.0/112)*bjms*bjms+(-243.0/280)*bjps*bjms+(99.0/560)*bjph*bjms+(297.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/560)*bjps*bjps+(-27.0/280)*bjph*bjps+(-23.0/210)*bjmh*bjph+(99.0/560)*bjms*bjph+(-27.0/280)*bjps*bjph+(7.0/240)*bjph*bjph);
    bx2uva12 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);
    bx2uva13 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);

    bx2uva21 = 2*idx*((587.0/1680)*bjmh*bjmh+(-39.0/112)*bjms*bjmh+(3.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-39.0/112)*bjmh*bjms+(243.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(33.0/560)*bjph*bjms+(3.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(81.0/560)*bjps*bjps+(-3.0/560)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(33.0/560)*bjms*bjph+(-3.0/560)*bjps*bjph+(-79.0/1680)*bjph*bjph);
    bx2uva22 = 2*idx*((13.0/42)*bjmh*bjmh+(-9.0/35)*bjms*bjmh+(-9.0/70)*bjps*bjmh+(8.0/105)*bjph*bjmh+(-9.0/35)*bjmh*bjms+(27.0/14)*bjms*bjms+(-54.0/35)*bjps*bjms+(-9.0/70)*bjph*bjms+(-9.0/70)*bjmh*bjps+(-54.0/35)*bjms*bjps+(27.0/14)*bjps*bjps+(-9.0/35)*bjph*bjps+(8.0/105)*bjmh*bjph+(-9.0/70)*bjms*bjph+(-9.0/35)*bjps*bjph+(13.0/42)*bjph*bjph);
    bx2uva23 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);

    bx2uva31 = 2*idx*((-127.0/1680)*bjmh*bjmh+(99.0/1120)*bjms*bjmh+(-9.0/560)*bjps*bjmh+(11.0/3360)*bjph*bjmh+(99.0/1120)*bjmh*bjms+(-81.0/560)*bjms*bjms+(81.0/1120)*bjps*bjms+(-9.0/560)*bjph*bjms+(-9.0/560)*bjmh*bjps+(81.0/1120)*bjms*bjps+(-81.0/560)*bjps*bjps+(99.0/1120)*bjph*bjps+(11.0/3360)*bjmh*bjph+(-9.0/560)*bjms*bjph+(99.0/1120)*bjps*bjph+(-127.0/1680)*bjph*bjph);
    bx2uva32 = 2*idx*((-79.0/1680)*bjmh*bjmh+(-3.0/560)*bjms*bjmh+(33.0/560)*bjps*bjmh+(-11.0/1680)*bjph*bjmh+(-3.0/560)*bjmh*bjms+(81.0/560)*bjms*bjms+(-81.0/560)*bjps*bjms+(3.0/560)*bjph*bjms+(33.0/560)*bjmh*bjps+(-81.0/560)*bjms*bjps+(243.0/560)*bjps*bjps+(-39.0/112)*bjph*bjps+(-11.0/1680)*bjmh*bjph+(3.0/560)*bjms*bjph+(-39.0/112)*bjps*bjph+(587.0/1680)*bjph*bjph);
    bx2uva33 = 2*idx*((7.0/240)*bjmh*bjmh+(-27.0/280)*bjms*bjmh+(99.0/560)*bjps*bjmh+(-23.0/210)*bjph*bjmh+(-27.0/280)*bjmh*bjms+(243.0/560)*bjms*bjms+(-243.0/280)*bjps*bjms+(297.0/560)*bjph*bjms+(99.0/560)*bjmh*bjps+(-243.0/280)*bjms*bjps+(243.0/112)*bjps*bjps+(-207.0/140)*bjph*bjps+(-23.0/210)*bjmh*bjph+(297.0/560)*bjms*bjph+(-207.0/140)*bjps*bjph+(1777.0/1680)*bjph*bjph);
    // LHS 
    
    LHSa11 = uva11 + h2uxvxa11 + hbxuvxa11 + hbxuxva11 + bx2uva11;
    LHSa12 = uva12 + h2uxvxa12 + hbxuvxa12 + hbxuxva12 + bx2uva12;
    LHSa13 = uva13 + h2uxvxa13 + hbxuvxa13 + hbxuxva13 + bx2uva13;

    LHSa21 = uva21 + h2uxvxa21 + hbxuvxa21 + hbxuxva21 + bx2uva21;
    LHSa22 = uva22 + h2uxvxa22 + hbxuvxa22 + hbxuxva22 + bx2uva22;
    LHSa23 = uva23 + h2uxvxa23 + hbxuvxa23 + hbxuxva23 + bx2uva23;

    LHSa31 = uva31 + h2uxvxa31 + hbxuvxa31 + hbxuxva31 + bx2uva31;
    LHSa32 = uva32 + h2uxvxa32 + hbxuvxa32 + hbxuxva32 + bx2uva32;
    LHSa33 = uva33 + h2uxvxa33 + hbxuvxa33 + hbxuxva33 + bx2uva33;
       

    Ared[j-1 + 3][1] = 0;

    Ared[j-1 +2][2] = Ared[j-1 + 2][2] + LHSa21;
    Ared[j + 2][2] = 0;

    Ared[j-1 + 1][3] = Ared[j-1 + 1][3] + LHSa11;
    Ared[j + 1][3] = Ared[j + 1][3] + LHSa22;
    Ared[j+1 + 1][3] = 1;

    Ared[j-1 + 1][4] = Ared[j-1 + 1][4] + LHSa12;
    Ared[j + 1][4] = Ared[j + 1][4] + LHSa23;
    
    Ared[j-1 + 1 ][5] = Ared[j-1 + 1][5] + LHSa13;

    
    nGis[j-1] = nGis[j-1] + Gva11; 
    nGis[j] = nGis[j] + Gva21; 
    nGis[j+1] = uMend[0];




    double d;
    d = bandec(Ared, m, 2,2, AL ,indx);

    banbks(Ared, m,2,2, AL, indx, nGis);


   
    memcpy(ubc + (unBC - 1),nGis,m*sizeof(double));
    memcpy(ubc,uMbeg,unBC*sizeof(double));
    memcpy(ubc + unbc- unBC,uMend,unBC*sizeof(double));

    free_dmatrix(Ared, 1, m, 1, 5);
    free_dmatrix(AL, 1, m, 1, 5);
    free(indx);
    free(nGis);

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
void evolveForce(double *hbc, double *Gbc, double *wbc, double *bbc, double *ubc,double g, double dx, double dt, int n, int hnBC, int hnbc, int bnBC, int bnbc, int unBC, int unbc, double *newG, double *newh, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9)
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

    fih = foh;
    fiG = foG;
    for(i = 0;i < n;i++)
    {

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

        hS = ht(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + Fluxhdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

        GS = Gt(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + FluxGdiff(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + SourceG(x[i],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

		newh[i] = hbc[hnBC + 3*(i) + 1] - dt*idx*(foh - fih) + dt*hS;
		newG[i] = Gbc[hnBC + 3*(i) + 1] - dt*idx*(foG - fiG) + dt*idx*(dx*sourcec)+ dt*GS;

        fih = foh;
        fiG = foG;

    }

}

void Stagehbed(int n, double *h, double *b, double *w)
{
    int i;
    for(i = 0 ; i < n; i++)
    {
        w[i] = h[i] + b[i];
    }

}


void evolvewrapForcingANA(double *h, double *G, double *b,double *hMbeg, double *hMend,double *GMbeg, double *GMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend, double *uMbeg, double *uMend,int n, int hnBC , int hnbc, int bnBC, int bnMBC , int bnbc, int unBC , int unbc,double theta,double dx, double dt, double g, double *x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9)
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

    double *bbc = malloc((bnbc)*sizeof(double));

    // calculate stage
    Stagehbed(n,h,b,w);

    //Reconstruct
    ReconLin(h, hMbeg, hMend,n,hnBC,hnbc,theta,hbc);
    ReconLin(G, GMbeg, GMend,n,hnBC,hnbc,theta,Gbc);
    ReconLin(w, wMbeg, wMend,n,hnBC,hnbc,theta,wbc);

    ReconQuart(b,bMbeg,bMend,n,bnMBC, bnBC, bnbc,bbc, dx);

    getufromG(hbc, Gbc, bbc,uMbeg, uMend,dx , n, u_length ,hnBC,bnBC,hnbc, bnbc, unBC, unbc, ubc);

    evolveForce(hbc,Gbc,wbc,bbc,ubc,g, dx,dt,n,hnBC,hnbc,bnBC,bnbc,unBC,unbc,Gp, hp,x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9); 

    Stagehbed(n,hp,b,wp);

    //Solve for u given wet/dry regions
    ReconLin(hp, hMbeg, hMend,n,hnBC,hnbc,theta,hbc);
    ReconLin(Gp, GMbeg, GMend,n,hnBC,hnbc,theta,Gbc);
    ReconLin(wp, wMbeg, wMend,n,hnBC,hnbc,theta,wbc);

    getufromG(hbc, Gbc, bbc,uMbeg, uMend,dx , n, u_length ,hnBC,bnBC,hnbc, bnbc, unBC, unbc, ubc);

    evolveForce(hbc,Gbc,wbc,bbc,ubc,g, dx,dt,n,hnBC,hnbc,bnBC,bnbc,unBC,unbc,Gpp, hpp,x,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9); 

    int i;
    for(i=0;i<n;i++)
    {
        G[i] =0.5*(G[i] + Gpp[i]);
        h[i] =0.5*(h[i]+ hpp[i]);
        
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
}


int main()
{
    printf("h");
    return 1;
}
