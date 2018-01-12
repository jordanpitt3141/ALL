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
#define div_0 1.0e-15
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

void RKstep(double *a, double *b, int n)
{
    int i;
    for(i =0;i <n;i++)
    {
        a[i] = 0.5*(a[i] + b[i]);

    }

}


void getufromG(double *h, double *G, double u0, double u1, double h0, double h1, double dx , int n, double *u)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    //double *a = malloc((n-1)*sizeof(double));
    //double *b = malloc(n*sizeof(double));
    //double *c = malloc((n-1)*sizeof(double));

    double **AL, **Ared;
    unsigned long *indx;

    Ared = dmatrix(1, n, 1, 3);
    AL = dmatrix(1, n, 1, 3);
    indx = mallocLongPy(n);
    
    double *nGis = malloc((n)*sizeof(double));

    int i;
    double thx;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);
        


        Ared[i + 1][1] =  -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
        Ared[i + 1][2] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
        Ared[i + 1][3] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

        nGis[i] = G[i];



    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);
        
    Ared[i + 1][1] =  -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
    Ared[i + 1][2] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    Ared[i + 1][3] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;


    nGis[i] = G[i] - u0*(-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx);

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);
        
    Ared[i + 1][1] =  -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
    Ared[i + 1][2] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    Ared[i + 1][3] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    nGis[i] = G[i] - u1*(-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx);

    double d;
    d = bandec(Ared, n, 1,1, AL ,indx);

    banbks(Ared, n,1,1, AL, indx, nGis);

    memcpy(u,nGis,n*sizeof(double));




    free_dmatrix(Ared, 1, n, 1, 3);
    free_dmatrix(AL, 1, n, 1, 3);
    free(indx);
    free(nGis);

}

void getGfromu(double *h, double *u, double u0, double u1, double h0, double h1, double dx , int n, double *G)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;

    int i;
    double thx;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);

        G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
                (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
                (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];
    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u0 +
            (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u[i+1];

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);

    G[i] = (-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx)*u[i-1] +
            (h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i])*u[i] +
            (-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx)*u1;
}

void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}

void evolveBC(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs,double *nG, double *nh, double *nu)
{ 
    //maybe an error in calculating G
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int j = nBCs -1;


    double *Gb = malloc(nBC*sizeof(double));
    double *Ge = malloc(nBC*sizeof(double));
    double *ub = malloc(nBC*sizeof(double));
    double *ue = malloc(nBC*sizeof(double));
    double *hb = malloc(nBC*sizeof(double));
    double *he = malloc(nBC*sizeof(double));
    double *u = malloc(n*sizeof(double));
    //ADD BC we get h,G need to solve for u then add in boundaries for G using u and h
    getufromG(h,G, u0[j], u1[0], h0[j], h1[0],dx,n,u);

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    int i = -1;
    int k = nBC -1;
    j = nBCs -1;


    double hx = 0.5*idx*(h[i+1] - h0[j-1]);

    double ai =0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    double bi = h0[j] + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j];
    double ci = -0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j];
    Gb[k] = ai*u0[j-1]
                 + bi*u0[j]
                 + ci*u[i+1];
    ub[k] = u0[j];
    hb[k] = h0[j];


    for(k = k-1; k > -1 ; k--)
    {
        j--;
        hx = 0.5*idx*(h0[j+1] - h0[j-1]);

        Gb[k] = (0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j-1]
                    + (h0[j] + 2*ithree*idx*idx*h0[j]*h0[j]*h0[j])*u0[j]
                    + (-0.5*idx*h0[j]*h0[j]*hx - idx*idx*ithree*h0[j]*h0[j]*h0[j])*u0[j+1];

        ub[k] = u0[j];
        hb[k] = h0[j];
    }

    //back end
    i = n;
    k = 0;
    j = 0;

    hx = 0.5*idx*(h1[j+1] - h[i-1]);
    Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u[i-1]
            + (h1[j] + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
            + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];

    ue[k] = u1[j];
    he[k] = h1[j];

    for(k = k+1; k < nBC ; k++)
    {
        j++;
        hx = 0.5*idx*(h1[j+1] - h1[j-1]);
        Ge[k] = (0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j-1]
                    + (h1[j] + 2*ithree*idx*idx*h1[j]*h1[j]*h1[j])*u1[j]
                    + (-0.5*idx*h1[j]*h1[j]*hx - idx*idx*ithree*h1[j]*h1[j]*h1[j])*u1[j+1];
        ue[k] = u1[j];
        he[k] = h1[j];

    }

    //bring them all together

    conc(Gb,G,Ge,nBC,n,nBC,nG);
    conc(hb,h,he,nBC,n,nBC,nh);
    conc(ub,u,ue,nBC,n,nBC,nu);

    free(Gb);
    free(Ge);
    free(hb);
    free(he);
    free(ub);
    free(ue);
    free(u);    

}

void evolve(double *G, double *h, double *u, double g, double dx, double dt, int nBC, int n,double *nG, double *nh)
{
    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

    double ue = 0.5*(u[i+1]+u[i]);

    //calculate values at right of i cell
    double hir = h[i];
    double Gir = G[i];
    double uir = u[i];

    //calculate values at left of i+1 cell
    double hip1l = h[i+1];
    double Gip1l = G[i+1];
    double uip1l = u[i + 1];

    //right force

    double duer = idx*(u[i+1] - u[i]);
    double duel = idx*(u[i+1] - u[i]);

    double sqrtghel = sqrt(g* hir);
    double sqrtgher = sqrt(g* hip1l);

    double sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
    double sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

    double felh = uir*hir;
    double felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duel*duel;
    double ferh = uip1l*hip1l;
    double ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*ithree*hip1l*hip1l*hip1l*duer*duer;

    double isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
    double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));

    double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

    double fih = foh;
    double fiG = foG;
    for (i = nBC ; i < n +nBC;i++)
    {
        //i right

        ue = 0.5*(u[i+1]+u[i]);
        hir = h[i];
        Gir = G[i];
        uir = u[i];

        //i+1 left

        hip1l = h[i+1];
        Gip1l = G[i+1];
        uip1l = u[i + 1];


        duer = idx*(u[i+1] - u[i]);
        duel = idx*(u[i] - u[i- 1]);

        sqrtghel = sqrt(g*hir);
        sqrtgher = sqrt(g*hip1l);

        sl = fmin(0,fmin(uir - sqrtghel, uip1l - sqrtgher));
        sr = fmax(0,fmax(uir + sqrtghel, uip1l + sqrtgher));

        felh = uir*hir;
        felG = Gir*uir + 0.5*g*hir*hir - 2*ithree*hir*hir*hir*duel*duel;
        ferh = uip1l*hip1l;
        ferG = Gip1l*uip1l + 0.5*g*hip1l*hip1l -2*ithree*hip1l*hip1l*hip1l*duer*duer;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);
   
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        //source term
        nh[i -nBC] = h[i] -dt*idx*(foh - fih);
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG);

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

double evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs)
{
    //first allocate memory for BC variables
    double *Gbc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    double *ubc = malloc((n + 2*nBC)*sizeof(double));

    evolveBC(G,h,h0,h1,u0,u1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc);

    double hmax,umax;
    hmax = dmaxarray(hbc,n + 2*nBC);
    umax = dmaxarray(ubc,n + 2*nBC);

    double ndt;

    ndt = 0.5 / (umax + sqrt(g*hmax)) * dx;

    //printf("%f | %f | %f | %f \n ",hmax,umax,dt,ndt);

    if(ndt > dt) ndt = dt;

    evolve(Gbc,hbc,ubc,g,dx,ndt,nBC,n,G,h);


    free(Gbc);
    free(hbc);
    free(ubc);

    return ndt;

}








int main()
{
    printf("h");
    return 1;
}
