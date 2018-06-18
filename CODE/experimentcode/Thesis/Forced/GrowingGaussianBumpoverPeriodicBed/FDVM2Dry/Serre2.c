#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

const double i24 = 1.0/24.0;
const double i48 = 1.0/48.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20
#define div_0 1.0e-8

#define htol 1.0e-4


#define NR_END 1


#define FREE_ARG char*

const double E = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069;


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


double *mallocPy(int n)
{
    double *x = malloc(n*sizeof(double));
    return x;
}


int main()
{
    printf("h");
    return 1;
}


// APPLICATION SPECIFIC
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
	else
	{
        return 0.0;
	}
}

void TDMA(double *a, double *b, double *c,double *d, int n, double *x)
{
	//tridiag matrix solver for [0abc0][x] = [d]
	double *alpha = malloc(n*sizeof(double));
	double *beta = malloc(n*sizeof(double));
	
	//int n will be size of b,d,x

	alpha[0] = c[0] / b[0];
	beta[0] = d[0] / b[0];
	int i;
	double m;
	for(i = 1; i < n -1;i++)
	{
		m = 1.0 / (b[i] - a[i-1]*alpha[i-1]);
		alpha[i] = c[i]*m;
		beta[i] = (d[i] - a[i-1]*beta[i-1])*m;
       
	}

	m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2]);
	beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m;

	x[n-1] = beta[n-1];

	for(i = n-2; i > -1; i--)
	{
		x[i] = beta[i] - alpha[i]*x[i+1];
       
	}

    free(alpha);
    free(beta);		
}


void getufromG(double *h, double *G,double *bed, double u0, double u1, double h0, double h1, double b0, double b1,double dx, int n,double *u)
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
		th = h[i];
		thx = 0.5*idx*(h[i+1] - h[i-1]);
		tbx = 0.5*idx*(bed[i+1] - bed[i-1]);
		tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1]);

		D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
		
		ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
		bi = D + 2.0*i3*idx*idx*th*th*th;
		ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;


        Ared[i + 1][1] = ai;
        Ared[i + 1][2] = bi;
        Ared[i + 1][3] = ci;
        nGis[i] = G[i];

	}

	//Boundary 
	//i = 0
	i = 0;
	th = h[i];
	thx = 0.5*idx*(h[i+1] - h0);
	tbx = 0.5*idx*(bed[i+1] - b0);
	tbxx = idx*idx*(bed[i+1] -2*bed[i]+ b0);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;

    Ared[i + 1][2] = bi;
    Ared[i + 1][3] = ci;
    nGis[i] = G[i] - u0*ai;

	i = n-1;
	
	th = h[i];
	thx = 0.5*idx*(h1- h[i-1]);
	tbx = 0.5*idx*(b1 - bed[i-1]);
	tbxx = idx*idx*(b1 -2*bed[i]+ bed[i-1]);

	D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx;
	
	ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx;
	bi = D + 2.0*i3*idx*idx*th*th*th;
	ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx;


    Ared[i + 1][1] = ai;
    Ared[i + 1][2] = bi;
    nGis[i] = G[i] - u1*ci;

	
    double d;
    d = bandec(Ared, n, 1,1, AL ,indx);

    banbks(Ared, n,1,1, AL, indx, nGis);


    memcpy(u,nGis,n*sizeof(double));

    free_dmatrix(Ared, 1, n, 1, 3);
    free_dmatrix(AL, 1,n, 1, 3);
    free(indx);
    free(nGis);

}

void edgevaluesSplit(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duEbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duEend, double *uEend, double *ddbCend, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double *hMbc, double *GMbc, double *wMbc, double *bMbc, double *duEbc, double *uEbc, double *ddbCbc, double dx, double theta)
{
// Take all h,G, bed and calculate all the edge values required for evolve, with split
    int nMbegI,nMendI,nEbegI,nEendI,nCbegI,nCendI;


    double bedjm1, hjm1, Gjm1,bedjp1, hjp1, Gjp1, uNp1,uNm1;                 


   double idx = 1.0 /dx;
   double wj, wjm1,wjp1,dwjf,dwjm,dwjb,dwjlim,dhjf,dhjm,dhjb,dhjlim,dGjf,dGjm,dGjb,dGjlim;
   int i,j;

   int *RegInd = RegSplit(h,n); 
   int m1 = 5;
   int RegIndLen = RegInd[m1-1];

    


    //put in BC
    for(i=0;i < nMBC;i++)
    {
        wMbc[i] = wMbeg[i];
        wMbc[nMbc-nMBC + i] = wMend[i];

        bMbc[i] = bMbeg[i];
        bMbc[nMbc-nMBC + i] = bMend[i];

        hMbc[i] = hMbeg[i];
        hMbc[nMbc-nMBC + i] = hMend[i];

        GMbc[i] = GMbeg[i];
        GMbc[nMbc-nMBC +i] = GMend[i];
    }

    for(i=0;i < nEBC;i++)
    {
        uEbc[i] = uEbeg[i];
        uEbc[nEbc-nEBC + i] = uEend[i];

        duEbc[i] = duEbeg[i];
        duEbc[nEbc-nEBC + i] = duEend[i];

    }

    for(i=0;i < nCBC;i++)
    {
        ddbCbc[i] = ddbCbeg[i];
        ddbCbc[nCbc-nCBC + i] = ddbCend[i];

    }

    //Ensure thatall the zeroing has been done, so we can read the boundary conditions from the neighbours
    for (i =0;i<RegIndLen;i++)
    {
        if(RegInd[i*m1 +  3] != 1)
        {

           // dry region
           for(j=RegInd[i*m1]; j <= RegInd[i*m1 +  1] ; j++)
           {


                // everything except bed and w is zero
                //u[(unBC-1) + 2*j] = 0;
                //u[(unBC-1) + 2*j +1] = 0;

                hMbc[nMBC + 3*j] = 0;
                hMbc[nMBC + 3*j + 1] = 0;
                hMbc[nMBC + 3*j + 2] = 0;

                GMbc[nMBC + 3*j] = 0;
                GMbc[nMBC + 3*j + 1] = 0;
                GMbc[nMBC + 3*j + 2] = 0;

                duEbc[nMBC + 2*j] = 0;
                duEbc[nMBC + 2*j + 1] = 0;

                uEbc[nEBC + 2*j] = 0 ;
                uEbc[nEBC + 2*j + 1] = 0;


                // possible undefined: bed[j+1] , bed[j-1]

                if(j == 0)
                {
                    bedjm1 = bMbeg[1];
                } 
                else
                {
                    bedjm1 = bed[j-1];
                    
                }

                if(j == n-1)
                {
                    bedjp1 = bMend[1];
                } 
                else
                {
                    bedjp1 = bed[j+1];
                    
                }



                wj = bed[j];
                wjp1 = bedjp1;
                wjm1 = bedjm1;
                dwjf = idx*(wjp1 - wj) ;
                dwjb = idx*(wj - wjm1) ;
                dwjm = 0.5*idx*(wjp1 - wjm1) ;
                dwjlim = minmod(theta*dwjf,dwjm,theta*dwjb);

                wMbc[nMBC + 3*j] = wj - 0.5*dx*dwjlim ;
                wMbc[nMBC + 3*j + 1] = wj;
                wMbc[nMBC + 3*j + 2] = wj + 0.5*dx*dwjlim ;

                bMbc[nMBC + 3*j] = wMbc[nMBC + 3*j];
                bMbc[nMBC + 3*j + 1] = wMbc[nMBC + 3*j + 1];
                bMbc[nMBC + 3*j + 2] = wMbc[nMBC+ 3*j + 2];


                ddbCbc[nCBC + j] = idx*idx*(bedjp1 - 2*bed[j] + bedjm1) ;

                

                
           }
        }

    }

    //have to zero first so we can use the BCs

    for (i =0;i<RegIndLen;i++)
    {
        if(RegInd[i*m1 +  3] == 1)
        {
            //wet region

            // where will BCS be             
            nMbegI = nMBC  + 3*RegInd[i*m1];
            nMendI = nMBC  + 3*(RegInd[i*m1 + 1] + 1);

            nEbegI = (nEBC -1) + 2*RegInd[i*m1];
            nEendI = (nEBC-1) + 2*(RegInd[i*m1 + 1] + 1);

            nCbegI = nCBC  + RegInd[i*m1];
            nCendI = nCBC  + (RegInd[i*m1 + 1] + 1);


            // have u solved
            // perform reconstruction
            double *uN = malloc((RegInd[i*m1 + 2] + 1)*sizeof(double));
            getufromG(h + RegInd[i*m1],G + RegInd[i*m1],bed + RegInd[i*m1],uEbc[nEbegI - 1],uEbc[nEendI+1], hMbc[nMbegI - 1],hMbc[nMendI+1], bMbc[nMbegI - 1],bMbc[nMendI+1],dx,RegInd[i*m1 + 2] + 1,uN);

           for(j=RegInd[i*m1]; j <= RegInd[i*m1 +  1] ; j++)
           {
                // all these j+1 and j possibly undefined, seperate cases if 
                if(j == 0)
                {
                    bedjm1 = bMbeg[1];
                    hjm1 = hMbeg[1];
                    Gjm1 = GMbeg[1];
                } 
                else
                {
                    bedjm1 = bed[j-1];
                    hjm1 = h[j-1];
                    Gjm1 = G[j - 1];
                    
                    
                }

                if(j == n-1)
                {
                    bedjp1 = bMend[1];
                    hjp1 = hMend[1];
                    Gjp1 = GMend[1];
                } 
                else
                {
                    bedjp1 = bed[j+1];
                    hjp1 = h[j+1];
                    Gjp1 = G[j + 1];
     
                }


                if(j == RegInd[i*m1])
                {
                    uNm1 = uEbc[nEbegI - 1];
                }
                else
                {
                    uNm1 = uN[j - RegInd[i*m1] - 1];
                }

                if(j == RegInd[i*m1+1])
                {
                    uNp1 = uEbc[nEendI + 1];
                }
                else
                {
                    uNp1 = uN[j - RegInd[i*m1] + 1];
                }

 
                //printf("%d : %f | %f | %f \n",j,uNm1,uN[j - RegInd[i*m1]],uNp1);
    

                //Reconstruct w and h -> b
                wj = h[j] + bed[j];
                wjp1 = hjp1 + bedjp1;
                wjm1 = hjm1 + bedjm1;
                dwjf = idx*(wjp1 - wj) ;
                dwjb = idx*(wj - wjm1) ;
                dwjm = 0.5*idx*(wjp1 - wjm1) ;
                dwjlim = minmod(theta*dwjf,dwjm,theta*dwjb);

                wMbc[nMBC + 3*j] = wj - 0.5*dx*dwjlim ;
                wMbc[nMBC + 3*j + 1] = wj;
                wMbc[nMBC + 3*j + 2] = wj + 0.5*dx*dwjlim ;

                dhjf = idx*(hjp1 - h[j]) ;
                dhjb = idx*(h[j] - hjm1) ;
                dhjm = 0.5*idx*(hjp1 - hjm1) ;
                dhjlim = minmod(theta*dhjf,dhjm,theta*dhjb);

                hMbc[nMBC + 3*j] = h[j] - 0.5*dx*dhjlim ;
                hMbc[nMBC + 3*j + 1] = h[j];
                hMbc[nMBC + 3*j + 2] = h[j] + 0.5*dx*dhjlim ;

                bMbc[nMBC + 3*j] = wMbc[nMBC + 3*j] - hMbc[nMBC + 3*j];
                bMbc[nMBC + 3*j + 1] = bed[j];
                bMbc[nMBC + 3*j + 2] = wMbc[nMBC + 3*j + 2] - hMbc[nMBC + 3*j + 2];

                dGjf = idx*(Gjp1 - G[j]) ;
                dGjb = idx*(G[j] - Gjm1) ;
                dGjm = 0.5*idx*(Gjp1 - Gjm1) ;
                dGjlim = minmod(theta*dGjf,dGjm,theta*dGjb);

                GMbc[nMBC + 3*j] = G[j] - 0.5*dx*dGjlim ;
                GMbc[nMBC + 3*j + 1] = G[j];
                GMbc[nMBC + 3*j + 2] = G[j] + 0.5*dx*dGjlim ;


                duEbc[nMBC + 2*j] = idx*0.5*(uNp1 - uNm1);
                duEbc[nMBC + 2*j + 1] = idx*(uNp1 -  uN[j - RegInd[i*m1]]);


                uEbc[nEBC + 2*j] = uN[j - RegInd[i*m1]] ;
                uEbc[nEBC + 2*j + 1] =  0.5*(uN[j - RegInd[i*m1]] +  uNp1) ;

                ddbCbc[nCBC + j] = idx*idx*(bedjp1 - 2*bed[j] + bedjm1) ;

           }
          free(uN);
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


void evolveDRY(double *hMbc, double *GMbc, double *wMbc, double *bMbc, double *duEbc, double *uEbc, double *ddbCbc, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double dx, double dt, double g, double theta, double *hblank, double *Gblank, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9)
{

    double idx = 1.0 / dx;  
	int i;

	double th,tu,tux,tbx,tbxx,sourcer,sourcel,sourcec;
    double hip1l,hir,Gip1l,Gir,wip1l,wir,bip1l,bir,ue,due,hil,bil,hS,GS;
    double her,Ger,ber,dber,uer,duer,hel,Gel,bel,dbel,uel,duel;
    double fhel,fher,fGel,fGer,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,foG,fih,fiG,himhp,hihm,hihp,nbi;

    //i=0, ghost cell caculate Fi+1/2
    i = 0;
    hip1l = hMbc[nMBC + 3*i];
    hir = hMbc[nMBC + 3*i -1];

    Gip1l = GMbc[nMBC + 3*i];
    Gir = GMbc[nMBC + 3*i -1];

    wip1l = wMbc[nMBC + 3*i];
    wir = wMbc[nMBC + 3*i -1];

    bip1l = bMbc[nMBC + 3*i];
    bir = bMbc[nMBC + 3*i -1];

    ue = uEbc[(nEBC-1) + 2*i];
    due = duEbc[(nEBC-1) + 2*i];

    nbi = fmax(bip1l, bir);
	hihm = fmax(0, wir - nbi);
	hihp = fmax(0, wip1l - nbi);

	her = hihp;
	Ger = Gip1l;
	ber = bip1l;
	dber = idx*(bMbc[nMBC + 3*(i+1)] - bip1l);
	uer  = ue;
	duer = due;

	hel = hihm;
	Gel = Gir;
	bel = bir;
	dbel = idx*(bMbc[nMBC + 3*i -2] -  bMbc[nMBC + 3*(i)  +1]);
	uel  = ue;
	duel = due;

	fhel = uel*hel;
	fher = uer*her;
	
	fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
	fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

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



    for(i =nCBC; i < n + nCBC;i++)
    {

        hip1l = hMbc[nMBC + 3*i];
        hir = hMbc[nMBC + 3*i -1];

        hil = hMbc[nMBC + 3*(i-1)];

        Gip1l = GMbc[nMBC + 3*i];
        Gir = GMbc[nMBC + 3*i -1];

        wip1l = wMbc[nMBC + 3*i];
        wir = wMbc[nMBC + 3*i -1];

        bip1l = bMbc[nMBC + 3*i];
        bir = bMbc[nMBC + 3*i -1];
        bil = bMbc[nMBC + 3*(i-1)];

        ue = uEbc[(nEBC-1) + 2*i];
        due = duEbc[(nEBC-1) + 2*i];

        nbi = fmax(bip1l, bir);
	    hihm = fmax(0, wir - nbi);
	    hihp = fmax(0, wip1l - nbi);

	    her = hihp;
	    Ger = Gip1l;
	    ber = bip1l;
	    dber = idx*(bMbc[nMBC + 3*(i+1)] - bip1l);
	    uer  = ue;
	    duer = due;

	    hel = hihm;
	    Gel = Gir;
	    bel = bir;
	    dbel = idx*(bir -  bMbc[nMBC + 3*(i-1)  - 1]);
	    uel  = ue;
	    duel = due;

	    fhel = uel*hel;
	    fher = uer*her;
	
	    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel;
	    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber;

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);	

	    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));


		th = hMbc[nMBC + 3*i -2];
		tu = uEbc[(nEBC-1) + 2*i -1];
		tux = duEbc[(nEBC-1) + 2*i -1];
		tbx = idx*(bir- bil);
		tbxx = ddbCbc[i];
		
		sourcer = g*0.5*(hihm*hihm - hir*hir);
		sourcec = -g*th*tbx -  0.5*th*th*tu*tux*tbxx + th*tu*tu*tbx*tbxx ;
		sourcel = g*0.5*(hil*hil - himhp*himhp);

        hS = ht(x[i- nCBC],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + Fluxhdiff(x[i- nCBC],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

        GS = Gt(x[i- nCBC],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + FluxGdiff(x[i- nCBC],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g) + SourceG(x[i- nCBC],t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g);

       
		hblank[i - nCBC] = hMbc[nMBC + 3*i -2] - dt*idx*(foh - fih) + dt*hS;
		Gblank[i - nCBC] = GMbc[nMBC + 3*i -2] - dt*idx*(foG -fiG) + dt*idx*(sourcer+sourcel + dx*sourcec) + dt*GS;

	    fih = foh;
	    fiG = foG;
	    himhp = hihp;

    }

}



void evolvewrapBC(double *h, double *G,double *bed, double *hMbeg, double *GMbeg, double  *wMbeg, double *bMbeg, double *duEbeg, double *uEbeg, double *ddbCbeg, double *hMend, double *GMend, double  *wMend, double *bMend, double *duEend, double *uEend, double *ddbCend, double *hMbeg1, double *GMbeg1, double  *wMbeg1, double *duEbeg1, double *uEbeg1, double *hMend1, double *GMend1, double  *wMend1, double *duEend1, double *uEend1, int nMBC, int nEBC, int nCBC, int n,int nMbc, int nEbc, int nCbc, double dx, double dt, double g, double theta, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7,double a8, double a9)
{ 

//UNFINISHED~~~
    double *hp = malloc((n)*sizeof(double));
    double *Gp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *Gpp = malloc((n)*sizeof(double));

	double *hMbc = malloc((nMbc)*sizeof(double));
	double *GMbc = malloc((nMbc)*sizeof(double));
	double *wMbc = malloc((nMbc)*sizeof(double));
	double *bMbc = malloc((nMbc)*sizeof(double));

	double *uEbc = malloc((nEbc)*sizeof(double));
	double *duEbc = malloc((nEbc)*sizeof(double));
    
    double *ddbCbc = malloc((nCbc)*sizeof(double));


    //update h' and G' 
    edgevaluesSplit(h,G,bed, hMbeg,GMbeg,wMbeg,bMbeg,duEbeg,uEbeg,ddbCbeg,hMend,GMend,wMend,bMend,duEend,uEend,ddbCend,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,hMbc,GMbc,wMbc,bMbc,duEbc,uEbc,ddbCbc,dx,theta);
	evolveDRY(hMbc,GMbc,wMbc,bMbc,duEbc,uEbc,ddbCbc,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta,hp,Gp,x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);

    edgevaluesSplit(hp,Gp,bed, hMbeg1,GMbeg1,wMbeg1,bMbeg,duEbeg1,uEbeg1,ddbCbeg,hMend1,GMend1,wMend1,bMend,duEend1,uEend1,ddbCend,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,hMbc,GMbc,wMbc,bMbc,duEbc,uEbc,ddbCbc,dx,theta);      
	evolveDRY(hMbc,GMbc,wMbc,bMbc,duEbc,uEbc,ddbCbc,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta,hpp,Gpp,x,t+ dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);


	int i;
	for(i = 0; i < n ; i++)
	{
		h[i] =0.5*(h[i] +hpp[i]);
		G[i] =0.5*(G[i] +Gpp[i]);
	}
	
    free(hp);
    free(Gp);
    free(hpp);
    free(Gpp);


    free(hMbc);
    free(GMbc);
    free(wMbc);
    free(bMbc);

    free(uEbc);
    free(duEbc);
    free(ddbCbc);

}




