#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//just do a simple first order Scheme, with no bed variations
const double i24 = 1.0/24.0;
const double i48 = 1.0/48.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;

const double E = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069;


void TDMA(double *a, double *b, double *c, double *d, int n, double *x)
{
	double *alpha = malloc(n*sizeof(double));
	double *beta = malloc(n*sizeof(double));

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];


    int i;
    double m;
    for(i=1; i<n-1; i++)
    {
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1]);
        alpha[i] = c[i]*m;
        beta[i] = (d[i] - a[i-1]*beta[i-1])*m;   
        
    }

    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2]);
    beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m;

    x[n-1] = beta[n-1];

    for (i=n-2; i > -1; i--)
    {
        x[i] = beta[i] - alpha[i]*x[i+1];  
   
    }

    free(alpha);
    free(beta);

}

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
    double fgp = 0.5*dx*sqrt(3.0/5.0) + x[j];
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
    double tgp = -0.5*dx*sqrt(3.0/5.0) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe);
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

void getufromG(double *h, double *G, double u0, double u1, double h0, double h1, double dx , int n, double *u)
{
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));

    int i;
    double thx;


    for (i =1;i < n-1 ; i++)
    {
        thx = 0.5*idx*(h[i+1] - h[i-1]);
        
        a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
        b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
        c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    }

    //Boundaries
    i = 0;

    thx = 0.5*idx*(h[i+1] - h0);
        
    //a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*th*th*thx
    b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    double tmpG1 = G[i];

    G[i] = tmpG1 - u0*(-ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx);

    i = n-1;

    thx = 0.5*idx*(h1 - h[i-1]);
        
    a[i-1] = -ithree*idx*idx*h[i]*h[i]*h[i] + 0.5*idx*h[i]*h[i]*thx;
    b[i] = h[i] + 2.0*ithree*idx*idx*h[i]*h[i]*h[i];
    //c[i] = -ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx;

    double tmpG2 = G[i];
    G[i] = tmpG2 - u1*(-ithree*idx*idx*h[i]*h[i]*h[i] - 0.5*idx*h[i]*h[i]*thx);

    TDMA(a,b,c,G,n,u);

    G[0] = tmpG1;
    G[n-1] = tmpG2;


    free(a);
    free(b);
    free(c);

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

double sech(double x)
{
    return 2.0/ (exp(x) + exp(-x));
}


double ht(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return -(a1*a2*a3*pow(sech(a2*(-(a3*t) + x)),2));
}


double Gt(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return (a0 + a1*tanh(a2*(-(a3*t) + x)))*(a4 + a5*tanh(a2*(-(a3*t) + x)))*((a1*a2*a3*a6*pow(a7,2)*pow(sech(a2*(-(a3*t) + x)),2)*sin(a7*x))/2. + 
      2*a1*pow(a2,2)*a3*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2)*tanh(a2*(-(a3*t) + x))) - 
   a2*a3*a5*pow(sech(a2*(-(a3*t) + x)),2)*(a0 + a1*tanh(a2*(-(a3*t) + x)))*(1 + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) + a1*a2*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2) - 
      (a6*pow(a7,2)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) - a1*a2*a3*pow(sech(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x)))*
    (1 + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) + a1*a2*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2) - (a6*pow(a7,2)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) + 
   (6*pow(a1,2)*pow(a2,3)*a3*a5*pow(sech(a2*(-(a3*t) + x)),6)*(a0 + a1*tanh(a2*(-(a3*t) + x))) - 
      18*a1*pow(a2,3)*a3*a5*pow(sech(a2*(-(a3*t) + x)),4)*tanh(a2*(-(a3*t) + x))*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2) - 2*pow(a2,3)*a3*a5*pow(sech(a2*(-(a3*t) + x)),4)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3) + 
      4*pow(a2,3)*a3*a5*pow(sech(a2*(-(a3*t) + x)),2)*pow(tanh(a2*(-(a3*t) + x)),2)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3))/3.;
}

double Fluxhdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a, double g)
{
    return a2*a5*pow(sech(a2*(-(a3*t) + x)),2)*(a0 + a1*tanh(a2*(-(a3*t) + x))) + a1*a2*pow(sech(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x)));
}

double FluxGdiff(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return a1*a2*g*pow(sech(a2*(-(a3*t) + x)),2)*(a0 + a1*tanh(a2*(-(a3*t) + x))) + pow(a2,2)*pow(a5,2)*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),4)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2) - 
   2*a1*pow(a2,3)*pow(a5,2)*pow(sech(a2*(-(a3*t) + x)),6)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2) + 
   (8*pow(a2,3)*pow(a5,2)*pow(sech(a2*(-(a3*t) + x)),4)*tanh(a2*(-(a3*t) + x))*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3))/3. + 
   2*a1*pow(a2,2)*a5*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),4)*(a0 + a1*tanh(a2*(-(a3*t) + x)))*(a4 + a5*tanh(a2*(-(a3*t) + x))) - 
   a2*a5*a6*pow(a7,2)*pow(sech(a2*(-(a3*t) + x)),2)*sin(a7*x)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x))) - 
   2*pow(a2,2)*a5*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2)*tanh(a2*(-(a3*t) + x))*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x))) + 
   a2*a5*pow(sech(a2*(-(a3*t) + x)),2)*((a0 + a1*tanh(a2*(-(a3*t) + x)))*(a4 + a5*tanh(a2*(-(a3*t) + x)))*
       (1 + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) + a1*a2*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2) - (a6*pow(a7,2)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) + 
      (-3*a1*pow(a2,2)*a5*pow(sech(a2*(-(a3*t) + x)),4)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2) + 2*pow(a2,2)*a5*pow(sech(a2*(-(a3*t) + x)),2)*tanh(a2*(-(a3*t) + x))*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3))/3.) + 
   (a4 + a5*tanh(a2*(-(a3*t) + x)))*((a0 + a1*tanh(a2*(-(a3*t) + x)))*(a4 + a5*tanh(a2*(-(a3*t) + x)))*
       (-2*pow(a6,2)*pow(a7,3)*cos(a7*x)*sin(a7*x) - (3*a1*a2*a6*pow(a7,2)*pow(sech(a2*(-(a3*t) + x)),2)*sin(a7*x))/2. - 2*a1*pow(a2,2)*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2)*tanh(a2*(-(a3*t) + x)) - 
         (a6*pow(a7,3)*cos(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) + a2*a5*pow(sech(a2*(-(a3*t) + x)),2)*(a0 + a1*tanh(a2*(-(a3*t) + x)))*
       (1 + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) + a1*a2*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2) - (a6*pow(a7,2)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) + 
      a1*a2*pow(sech(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x)))*(1 + pow(a6,2)*pow(a7,2)*pow(cos(a7*x),2) + a1*a2*a6*a7*cos(a7*x)*pow(sech(a2*(-(a3*t) + x)),2) - 
         (a6*pow(a7,2)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))))/2.) + (-6*pow(a1,2)*pow(a2,3)*a5*pow(sech(a2*(-(a3*t) + x)),6)*(a0 + a1*tanh(a2*(-(a3*t) + x))) + 
         18*a1*pow(a2,3)*a5*pow(sech(a2*(-(a3*t) + x)),4)*tanh(a2*(-(a3*t) + x))*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2) + 2*pow(a2,3)*a5*pow(sech(a2*(-(a3*t) + x)),4)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3) - 
         4*pow(a2,3)*a5*pow(sech(a2*(-(a3*t) + x)),2)*pow(tanh(a2*(-(a3*t) + x)),2)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),3))/3.);
}

double SourceG(double x, double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7, double g)
{
    return a6*a7*g*cos(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x))) - (a2*a5*a6*pow(a7,2)*pow(sech(a2*(-(a3*t) + x)),2)*sin(a7*x)*pow(a0 + a1*tanh(a2*(-(a3*t) + x)),2)*(a4 + a5*tanh(a2*(-(a3*t) + x))))/2. + 
   pow(a6,2)*pow(a7,3)*cos(a7*x)*sin(a7*x)*(a0 + a1*tanh(a2*(-(a3*t) + x)))*pow(a4 + a5*tanh(a2*(-(a3*t) + x)),2);
}


void evolve(double *G, double *h, double *u, double g, double dx, double dt, int nBC, int n,double *nG, double *nh, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7)
{
    double hS, GS;

    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = nBC - 1;

    double ue = 0.5*(u[i+1]+u[i]);

    //calculate values at right of i cell
    double hir = h[i];
    double Gir = G[i];
    double uir = ue;

    //calculate values at left of i+1 cell
    double hip1l = h[i+1];
    double Gip1l = G[i+1];
    double uip1l = ue;

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
        uir = ue;

        //i+1 left

        hip1l = h[i+1];
        Gip1l = G[i+1];
        uip1l = ue;


        duer = idx*(u[i+1] - u[i]);
        duel = idx*(u[i+1] - u[i]);

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


        hS = ht(x[i-nBC],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + Fluxhdiff(x[i-nBC],t,a0,a1,a2,a3,a4,a5,a6,a7,g);

        GS = Gt(x[i-nBC],t,a0,a1,a2,a3,a4,a5,a6,a7,g) + FluxGdiff(x[i-nBC],t,a0,a1,a2,a3,a4,a5,a6,a7,g);

   
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir));

        //source term
        nh[i -nBC] = h[i] -dt*idx*(foh - fih) + dt*hS;
        nG[i -nBC] = G[i] -dt*idx*(foG -fiG) + dt*GS;

        fih = foh;
        fiG = foG;  
 
    }
    
   
}

void evolvewrap(double *G, double *h, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs, double *x,double t,double a0,double a1,double a2, double a3, double a4, double a5, double a6, double a7)
{
    //first allocate memory for BC variables
    double *Gbc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));
    double *ubc = malloc((n + 2*nBC)*sizeof(double));

    evolveBC(G,h,h0,h1,u0,u1,g,dx,dt,nBC,n,nBCs,Gbc,hbc,ubc);
    evolve(Gbc,hbc,ubc,g,dx,dt,nBC,n,G,h,x,t,a0,a1,a2,a3,a4,a5,a6,a7);


    free(Gbc);
    free(hbc);
    free(ubc);

}


double *mallocPy(int n)
{
    double *x = malloc(n*sizeof(double));
    return x;
}

void writetomem(double *x, int i , double f)
{
    x[i] = f;

}

double readfrommem(double*x,int i)
{
    return x[i];
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