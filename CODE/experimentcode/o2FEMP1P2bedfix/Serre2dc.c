#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//put in slope limiting.

const double i24 = 1.0/24.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;
const double i48 = 1.0/48.0;
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

void PENT(double *e, double *a, double *d, double *c,double *f, double *B, int n, double *x)
{
    //works but is destructive on inputs
    int i;
    double m;
    for(i=1; i<n-1; i++)
    {
        m = a[i-1] /d[i-1];
        
        d[i] = d[i] - m*c[i-1];
        c[i] = c[i] - m*f[i-1];
        B[i] = B[i] - m*B[i-1];

        m = e[i-1] /d[i-1];
        a[i] = a[i] - m*c[i-1];
        d[i+1] = d[i+1] - m*f[i-1];
        B[i+1] = B[i+1] - m*B[i-1];
    }

    m = a[n-2] / d[n-2];
    d[n-1] = d[n-1] - m*c[n-2];
    x[n-1] = (B[n-1] - m*B[n-2]) / d[n-1];
    x[n-2] = (B[n-2] - c[n-2]*x[n-1]) / d[n-2];

    for (i=n-3; i > -1; i--)
    {
        x[i] = (B[i] - f[i]*x[i+2] - c[i]*x[i+1]) / d[i];  
   
    }
}

void getufromG(double *h, double *G, double *bed, double *hebeg, double *heend, double *Gebeg, double *Geend,double *uebeg, double *ueend, double *bedbeg, double *bedend, double theta, double dx , int n, int m, int nBC,int nhbc , double *u, double *hhbc,double *Ghbc, double *bedhbc)
{
    //trying to join different B.C systems...., also try incorporating into  C code that already has beds, like o2bedfix.

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now

    double idx = 1.0 / dx;
    double *uais = malloc((m-2)*sizeof(double));
    double *ubis = malloc((m-1)*sizeof(double));
    double *ucis = malloc((m)*sizeof(double));
    double *udis = malloc((m-1)*sizeof(double));
    double *ueis = malloc((m-2)*sizeof(double));

    double *nGis = malloc((m)*sizeof(double));

    int i,j;

    double dGib,dGim,dGif,dhib,dhim,dhif,dGi,dhi,Gjphm,Gjmhp,hjphm,hjmhp,bedai,bedbi,bedci,bjmh,bjph,bj;

    double Gintia11,Gintia21,Gintia31,uhintia11,uhintia12,uhintia13,uhintia21,uhintia22,uhintia23,uhintia31,uhintia32,uhintia33;
    double h3uxintia11,h3uxintia12,h3uxintia13,h3uxintia21,h3uxintia22,h3uxintia23,h3uxintia31,h3uxintia32,h3uxintia33;
    double h2bxuxva11,h2bxuxva12,h2bxuxva13,h2bxuxva21,h2bxuxva22,h2bxuxva23,h2bxuxva31,h2bxuxva32,h2bxuxva33;
    double h2bxuvxa11,h2bxuvxa12,h2bxuvxa13,h2bxuvxa21,h2bxuvxa22,h2bxuvxa23,h2bxuvxa31,h2bxuvxa32,h2bxuvxa33;
    double hbxbxuva11,hbxbxuva12,hbxbxuva13,hbxbxuva21,hbxbxuva22,hbxbxuva23,hbxbxuva31,hbxbxuva32,hbxbxuva33;
    double LHSa11,LHSa12,LHSa13,LHSa21,LHSa22,LHSa23,LHSa31,LHSa32,LHSa33;

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

        Gjphm= G[i] + 0.5*dGi;
        Gjmhp= G[i] - 0.5*dGi;

        hjphm= h[i] + 0.5*dhi;
        hjmhp= h[i] - 0.5*dhi;

        // reconstruct bed

        bedai = 0.5*idx*idx*(bed[i +1] - 2*bed[i] + bed[i  -1]);
        bedbi = 0.5*idx*(bed[i +1] - bed[i -1]);
        bedci = bed[i];

        bjmh = bedai*(0.5*dx)*(0.5*dx) - bedbi*0.5*dx + bedci;
        bjph = bedai*(0.5*dx)*(0.5*dx) + bedbi*0.5*dx + bedci;
        bj = bed[i];

        // G integral (RHS)
        Gintia11 = dx/6.0*(Gjmhp);
        Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm);
        Gintia31 = dx/6.0*(Gjphm);
        
        
        //uh integral
        uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm);
        uhintia12 = (1.0 /60.0)*dx*(4*hjmhp );
        uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
        
        uhintia21 = (1.0 /60.0)*dx*(4*hjmhp);
        uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm);
        uhintia23 = (1.0 /60.0)*dx*(4*hjphm);
        
        uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
        uhintia32 = (1.0 /60.0)*dx*(4*hjphm);
        uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm);
        
        //h3ux
        
        h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm);     
        h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm);        
        h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm);
        
        
        h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm) ;       
        h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm);       
        h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);
        
        h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm) ;       
        h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);        
        h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm);
        
        //h2 bx ux v

        h2bxuxva11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh 
                        +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj 
                        +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph ) ;   

        h2bxuxva12  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                        +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph );  

        h2bxuxva13  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
                        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
                        +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph );

        h2bxuxva21  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
                        +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph ) ;    

        h2bxuxva22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
                        +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
                        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;                         

        h2bxuxva23  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
                        +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ) ;
        
        h2bxuxva31  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
                +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
                +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph ) ;
        
        h2bxuxva32  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
                +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
                +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );
        
        h2bxuxva33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
                +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj
                +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ; 
        
        //h2 bx u vx
        
        h2bxuvxa11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh
                +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj
                +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph );  
        
        h2bxuvxa12  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
            +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph );
        
        h2bxuvxa13  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
            +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
            +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph);  

        h2bxuvxa21  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
            +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
            +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph ); 

        h2bxuvxa22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
            +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
            +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;

        h2bxuvxa23  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
            +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );         

        h2bxuvxa31  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
            +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
            +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph ) ;

        h2bxuvxa32  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
            +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
            +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ); 

        h2bxuvxa33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj 
            +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ;


        //h bx bx u v
        
        hbxbxuva11  = 2*idx*(     ((337.0/840)*hjmhp + (31.0/840)*hjphm )*bjmh*bjmh
                              + 2*((-1.0/2)*hjmhp + (-3.0/70)*hjphm )*bjmh*bj
                              + 2*((83.0/840)*hjmhp + (1.0/168)*hjphm )*bjmh*bjph
                              +   ((22.0/35)*hjmhp + (2.0/35)*hjphm )*bj*bj
                              + 2*((-9.0/70)*hjmhp + (-1.0/70)*hjphm )*bj*bjph
                              +   ((5.0/168)*hjmhp + (1.0/120)*hjphm )*bjph*bjph ) ;

        hbxbxuva12  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh
                              + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj
                              + 2*((1.0/42)*hjmhp + (0.0)*hjphm )*bjmh*bjph
                              +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj
                              + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph
                              +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

        hbxbxuva13  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                              + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj 
                              + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                              +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                              + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                              +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;
        
        hbxbxuva21  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh 
                              + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj 
                              + 2*((1.0/42)*hjmhp + (0)*hjphm )*bjmh*bjph 
                              +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj 
                              + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph 
                              +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

        hbxbxuva22  = 2*idx*(     ((2.0/7)*hjmhp + (2.0/15)*hjphm )*bjmh*bjmh 
                              + 2*((-8.0/35)*hjmhp + (-8.0/105)*hjphm )*bjmh*bj 
                              + 2*((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bjmh*bjph 
                              +   ((32.0/105)*hjmhp + (32.0/105)*hjphm )*bj*bj 
                              + 2*((-8.0/105)*hjmhp + (-8.0/35)*hjphm )*bj*bjph 
                              +   ((2.0/15)*hjmhp + (2.0/7)*hjphm )*bjph*bjph ) ;

        hbxbxuva23  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh 
                              + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj 
                              + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph 
                              +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj 
                              + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph 
                              +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ) ;

        hbxbxuva31  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                              + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj
                              + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                              +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                              + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                              +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;

        hbxbxuva32  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh
                              + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj
                              + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph
                              +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj
                              + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph
                              +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ); 
        
        hbxbxuva33  = 2*idx*(     ((1.0/120)*hjmhp + (5.0/168)*hjphm )*bjmh*bjmh
                              + 2*((-1.0/70)*hjmhp + (-9.0/70)*hjphm )*bjmh*bj
                              + 2*((1.0/168)*hjmhp + (83.0/840)*hjphm )*bjmh*bjph
                              +   ((2.0/35)*hjmhp + (22.0/35)*hjphm )*bj*bj
                              + 2*((-3.0/70)*hjmhp + (-1.0/2)*hjphm )*bj*bjph
                              +   ((31.0/840)*hjmhp + (337.0/840)*hjphm )*bjph*bjph );        

        // LHS 
        
        LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbxbxuva11;  
        LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbxbxuva12; 
        LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbxbxuva13; 
        LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbxbxuva21; 
        LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbxbxuva22; 
        LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbxbxuva23; 
        LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbxbxuva31; 
        LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbxbxuva32;  
        LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbxbxuva33; 
        

        uais[j-1] = uais[j-1] + LHSa31;
        
        ubis[j-1] = ubis[j-1] + LHSa21;
        ubis[j] = ubis[j] + LHSa32;
        
        ucis[j-1] = ucis[j-1] + LHSa11;
        ucis[j] = ucis[j] + LHSa22;
        ucis[j+1] = ucis[j+1] + LHSa33;
        
        udis[j-1] = udis[j-1] + LHSa12;
        udis[j] = udis[j] + LHSa23;
        
        ueis[j-1] = ueis[j-1] + LHSa13;
        
        
        nGis[j-1] = nGis[j-1] + Gintia11;  
        nGis[j] = nGis[j] + Gintia21;
        nGis[j+1] = nGis[j+1] + Gintia31; 

        hhbc[3*i+1] = hjmhp;
        hhbc[3*i+2] = h[i];
        hhbc[3*i+3] = hjphm; 

        Ghbc[3*i+1] = Gjmhp;
        Ghbc[3*i+2] = G[i];
        Ghbc[3*i+3] = Gjphm;

        bedhbc[3*i+1] = bjmh;
        bedhbc[3*i+2] = bj;
        bedhbc[3*i+3] = bjph;
        
        
        j = j + 2 ;       


    }


    // first 
    i = 0;
    j = 1;


    // Reconstruct G and h
    Gjphm= Gebeg[nBC-1] + 2*(G[0] - Gebeg[nBC-1]);
    Gjmhp= Gebeg[nBC-1];

    hjphm= hebeg[nBC-1] + 2*(h[0] - hebeg[nBC-1]);
    hjmhp=hebeg[nBC-1];

    //Reconstruct bed
    bedai = 0.5*idx*idx*(bed[i +1] - 2*bed[i] + bedbeg[nBC-1]);
    bedbi = 0.5*idx*(bed[i +1] - bedbeg[nBC-1]);
    bedci = bed[i];

    bjmh = bedai*(0.5*dx)*(0.5*dx) - bedbi*0.5*dx + bedci;
    bjph = bedai*(0.5*dx)*(0.5*dx) + bedbi*0.5*dx + bedci;
    bj = bed[i];

    // G integral (RHS)
    Gintia11 = dx/6.0*(Gjmhp);
    Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm);
    Gintia31 = dx/6.0*(Gjphm);
    
    
    //uh integral
    uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm);
    uhintia12 = (1.0 /60.0)*dx*(4*hjmhp );
    uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
    
    uhintia21 = (1.0 /60.0)*dx*(4*hjmhp);
    uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm);
    uhintia23 = (1.0 /60.0)*dx*(4*hjphm);
    
    uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
    uhintia32 = (1.0 /60.0)*dx*(4*hjphm);
    uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm);
    
    //h3ux
    
    h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm);     
    h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm);        
    h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm);
    
    
    h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm) ;       
    h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm);       
    h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);
    
    h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm) ;       
    h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);        
    h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm);
    
    //h2 bx ux v

    h2bxuxva11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh 
                    +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj 
                    +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph ) ;   

    h2bxuxva12  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                    +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                    +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph );  

    h2bxuxva13  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
                    +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
                    +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph );

    h2bxuxva21  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
                    +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                    +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph ) ;    

    h2bxuxva22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
                    +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
                    +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;                         

    h2bxuxva23  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                    +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
                    +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ) ;
    
    h2bxuxva31  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
            +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
            +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph ) ;
    
    h2bxuxva32  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
            +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );
    
    h2bxuxva33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj
            +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ; 
    
    //h2 bx u vx
    
    h2bxuvxa11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh
            +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj
            +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph );  
    
    h2bxuvxa12  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
        +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph );
    
    h2bxuvxa13  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
        +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
        +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph);  

    h2bxuvxa21  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
        +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph ); 

    h2bxuvxa22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
        +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;

    h2bxuvxa23  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );         

    h2bxuvxa31  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
        +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph ) ;

    h2bxuvxa32  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
        +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ); 

    h2bxuvxa33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj 
        +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ;


    //h bx bx u v
    
    hbxbxuva11  = 2*idx*(     ((337.0/840)*hjmhp + (31.0/840)*hjphm )*bjmh*bjmh
                          + 2*((-1.0/2)*hjmhp + (-3.0/70)*hjphm )*bjmh*bj
                          + 2*((83.0/840)*hjmhp + (1.0/168)*hjphm )*bjmh*bjph
                          +   ((22.0/35)*hjmhp + (2.0/35)*hjphm )*bj*bj
                          + 2*((-9.0/70)*hjmhp + (-1.0/70)*hjphm )*bj*bjph
                          +   ((5.0/168)*hjmhp + (1.0/120)*hjphm )*bjph*bjph ) ;

    hbxbxuva12  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh
                          + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj
                          + 2*((1.0/42)*hjmhp + (0.0)*hjphm )*bjmh*bjph
                          +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj
                          + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph
                          +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

    hbxbxuva13  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                          + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj 
                          + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                          +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                          + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                          +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;
    
    hbxbxuva21  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh 
                          + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj 
                          + 2*((1.0/42)*hjmhp + (0)*hjphm )*bjmh*bjph 
                          +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj 
                          + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph 
                          +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

    hbxbxuva22  = 2*idx*(     ((2.0/7)*hjmhp + (2.0/15)*hjphm )*bjmh*bjmh 
                          + 2*((-8.0/35)*hjmhp + (-8.0/105)*hjphm )*bjmh*bj 
                          + 2*((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bjmh*bjph 
                          +   ((32.0/105)*hjmhp + (32.0/105)*hjphm )*bj*bj 
                          + 2*((-8.0/105)*hjmhp + (-8.0/35)*hjphm )*bj*bjph 
                          +   ((2.0/15)*hjmhp + (2.0/7)*hjphm )*bjph*bjph ) ;

    hbxbxuva23  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh 
                          + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj 
                          + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph 
                          +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj 
                          + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph 
                          +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ) ;

    hbxbxuva31  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                          + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj
                          + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                          +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                          + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                          +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;

    hbxbxuva32  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh
                          + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj
                          + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph
                          +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj
                          + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph
                          +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ); 
    
    hbxbxuva33  = 2*idx*(     ((1.0/120)*hjmhp + (5.0/168)*hjphm )*bjmh*bjmh
                          + 2*((-1.0/70)*hjmhp + (-9.0/70)*hjphm )*bjmh*bj
                          + 2*((1.0/168)*hjmhp + (83.0/840)*hjphm )*bjmh*bjph
                          +   ((2.0/35)*hjmhp + (22.0/35)*hjphm )*bj*bj
                          + 2*((-3.0/70)*hjmhp + (-1.0/2)*hjphm )*bj*bjph
                          +   ((31.0/840)*hjmhp + (337.0/840)*hjphm )*bjph*bjph );        

    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbxbxuva11;  
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbxbxuva12; 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbxbxuva13; 
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbxbxuva21; 
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbxbxuva22; 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbxbxuva23; 
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbxbxuva31; 
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbxbxuva32;  
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbxbxuva33; 
    

    uais[j-1] = uais[j-1] + LHSa31;
    
    ubis[j-1] = ubis[j-1] + LHSa21;
    ubis[j] = ubis[j] + LHSa32;
    
    ucis[j-1] = 1;
    ucis[j] = ucis[j] + LHSa22;
    ucis[j+1] = ucis[j+1] + LHSa33;
    
    udis[j-1] = 0;
    udis[j] = udis[j] + LHSa23;
    
    ueis[j-1] = 0;
    
    
    nGis[j-1] = uebeg[nBC-1];  
    nGis[j] = nGis[j] + Gintia21;
    nGis[j+1] = nGis[j+1] + Gintia31; 

    hhbc[3*i+1] = hjmhp;
    hhbc[3*i+2] = h[i];
    hhbc[3*i+3] = hjphm; 

    Ghbc[3*i+1] = Gjmhp;
    Ghbc[3*i+2] = G[i];
    Ghbc[3*i+3] = Gjphm;

    bedhbc[3*i+1] = bjmh;
    bedhbc[3*i+2] = bj;
    bedhbc[3*i+3] = bjph;
  

// last
    j = m-2;
    i = n-1;

    Gjmhp= 2*G[i] - Geend[0];
    Gjphm= Geend[0];

    hjmhp= 2*h[i] - heend[0];
    hjphm= heend[0];

    // reconstruct bed

    bedai = 0.5*idx*idx*(bedend[0] - 2*bed[i] + bed[i  -1]);
    bedbi = 0.5*idx*(bedend[0]- bed[i -1]);
    bedci = bed[i];

    bjmh = bedai*(0.5*dx)*(0.5*dx) - bedbi*0.5*dx + bedci;
    bjph = bedai*(0.5*dx)*(0.5*dx) + bedbi*0.5*dx + bedci;
    bj = bed[i];

    // G integral (RHS)
    Gintia11 = dx/6.0*(Gjmhp);
    Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm);
    Gintia31 = dx/6.0*(Gjphm);
    
    
    //uh integral
    uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm);
    uhintia12 = (1.0 /60.0)*dx*(4*hjmhp );
    uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
    
    uhintia21 = (1.0 /60.0)*dx*(4*hjmhp);
    uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm);
    uhintia23 = (1.0 /60.0)*dx*(4*hjphm);
    
    uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm);
    uhintia32 = (1.0 /60.0)*dx*(4*hjphm);
    uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm);
    
    //h3ux
    
    h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm);     
    h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm);        
    h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm);
    
    
    h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm) ;       
    h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm);       
    h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);
    
    h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm) ;       
    h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm);        
    h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm);
    
    //h2 bx ux v

    h2bxuxva11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh 
                    +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj 
                    +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph ) ;   

    h2bxuxva12  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                    +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                    +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph );  

    h2bxuxva13  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
                    +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
                    +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph );

    h2bxuxva21  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
                    +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
                    +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph ) ;    

    h2bxuxva22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
                    +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
                    +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;                         

    h2bxuxva23  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
                    +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
                    +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ) ;
    
    h2bxuxva31  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
            +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
            +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph ) ;
    
    h2bxuxva32  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
            +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );
    
    h2bxuxva33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
            +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj
            +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ; 
    
    //h2 bx u vx
    
    h2bxuvxa11  = -idx*(((83.0/168)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (3.0/280)*hjphm*hjphm)*bjmh
            +   ((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bj
            +   ((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjph );  
    
    h2bxuvxa12  = -idx*(((23.0/70)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjmh
        +   ((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjph );
    
    h2bxuvxa13  = -idx*(((-47.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (23.0/840)*hjphm*hjphm)*bjmh
        +   ((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bj
        +   ((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjph);  

    h2bxuvxa21  = -idx*(((-127.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
        +   ((26.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bj
        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bjph ); 

    h2bxuvxa22  = -idx*(((-34.0/105)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-2.0/35)*hjphm*hjphm)*bjmh
        +   ((8.0/21)*hjmhp*hjmhp + 2*(16.0/105)*hjmhp*hjphm + (8.0/21)*hjphm*hjphm)*bj
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bjph ) ;

    h2bxuvxa23  = -idx*(((13.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (-29.0/210)*hjphm*hjphm)*bjmh
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(2.0/35)*hjmhp*hjphm + (26.0/35)*hjphm*hjphm)*bj
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bjph );         

    h2bxuvxa31  = -idx*(((31.0/280)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (-1.0/168)*hjphm*hjphm)*bjmh
        +   ((-29.0/210)*hjmhp*hjmhp + 2*(1.0/210)*hjmhp*hjphm + (13.0/210)*hjphm*hjphm)*bj
        +   ((23.0/840)*hjmhp*hjmhp + 2*(-3.0/280)*hjmhp*hjphm + (-47.0/840)*hjphm*hjphm)*bjph ) ;

    h2bxuvxa32  = -idx*(((-1.0/210)*hjmhp*hjmhp + 2*(-1.0/35)*hjmhp*hjphm + (-1.0/210)*hjphm*hjphm)*bjmh
        +   ((-2.0/35)*hjmhp*hjmhp + 2*(-8.0/105)*hjmhp*hjphm + (-34.0/105)*hjphm*hjphm)*bj
        +   ((13.0/210)*hjmhp*hjmhp + 2*(11.0/105)*hjmhp*hjphm + (23.0/70)*hjphm*hjphm)*bjph ); 

    h2bxuvxa33  = -idx*(((-1.0/168)*hjmhp*hjmhp + 2*(1.0/168)*hjmhp*hjphm + (31.0/280)*hjphm*hjphm)*bjmh
        +   ((-1.0/210)*hjmhp*hjmhp + 2*(-13.0/210)*hjmhp*hjphm + (-127.0/210)*hjphm*hjphm)*bj 
        +   ((3.0/280)*hjmhp*hjmhp + 2*(47.0/840)*hjmhp*hjphm + (83.0/168)*hjphm*hjphm)*bjph ) ;


    //h bx bx u v
    
    hbxbxuva11  = 2*idx*(     ((337.0/840)*hjmhp + (31.0/840)*hjphm )*bjmh*bjmh
                          + 2*((-1.0/2)*hjmhp + (-3.0/70)*hjphm )*bjmh*bj
                          + 2*((83.0/840)*hjmhp + (1.0/168)*hjphm )*bjmh*bjph
                          +   ((22.0/35)*hjmhp + (2.0/35)*hjphm )*bj*bj
                          + 2*((-9.0/70)*hjmhp + (-1.0/70)*hjphm )*bj*bjph
                          +   ((5.0/168)*hjmhp + (1.0/120)*hjphm )*bjph*bjph ) ;

    hbxbxuva12  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh
                          + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj
                          + 2*((1.0/42)*hjmhp + (0.0)*hjphm )*bjmh*bjph
                          +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj
                          + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph
                          +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

    hbxbxuva13  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                          + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj 
                          + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                          +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                          + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                          +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;
    
    hbxbxuva21  = 2*idx*(     ((13.0/70)*hjmhp + (4.0/105)*hjphm )*bjmh*bjmh 
                          + 2*((-22.0/105)*hjmhp + (-4.0/105)*hjphm )*bjmh*bj 
                          + 2*((1.0/42)*hjmhp + (0)*hjphm )*bjmh*bjph 
                          +   ((8.0/35)*hjmhp + (0)*hjphm )*bj*bj 
                          + 2*((-2.0/105)*hjmhp + (4.0/105)*hjphm )*bj*bjph 
                          +   ((-1.0/210)*hjmhp + (-4.0/105)*hjphm )*bjph*bjph ); 

    hbxbxuva22  = 2*idx*(     ((2.0/7)*hjmhp + (2.0/15)*hjphm )*bjmh*bjmh 
                          + 2*((-8.0/35)*hjmhp + (-8.0/105)*hjphm )*bjmh*bj 
                          + 2*((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bjmh*bjph 
                          +   ((32.0/105)*hjmhp + (32.0/105)*hjphm )*bj*bj 
                          + 2*((-8.0/105)*hjmhp + (-8.0/35)*hjphm )*bj*bjph 
                          +   ((2.0/15)*hjmhp + (2.0/7)*hjphm )*bjph*bjph ) ;

    hbxbxuva23  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh 
                          + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj 
                          + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph 
                          +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj 
                          + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph 
                          +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ) ;

    hbxbxuva31  = 2*idx*(     ((-31.0/840)*hjmhp + (-1.0/120)*hjphm )*bjmh*bjmh
                          + 2*((3.0/70)*hjmhp + (1.0/70)*hjphm )*bjmh*bj
                          + 2*((-1.0/168)*hjmhp + (-1.0/168)*hjphm )*bjmh*bjph
                          +   ((-2.0/35)*hjmhp + (-2.0/35)*hjphm )*bj*bj
                          + 2*((1.0/70)*hjmhp + (3.0/70)*hjphm )*bj*bjph
                          +   ((-1.0/120)*hjmhp + (-31.0/840)*hjphm )*bjph*bjph ) ;

    hbxbxuva32  = 2*idx*(     ((-4.0/105)*hjmhp + (-1.0/210)*hjphm )*bjmh*bjmh
                          + 2*((4.0/105)*hjmhp + (-2.0/105)*hjphm )*bjmh*bj
                          + 2*((0)*hjmhp + (1.0/42)*hjphm )*bjmh*bjph
                          +   ((0)*hjmhp + (8.0/35)*hjphm )*bj*bj
                          + 2*((-4.0/105)*hjmhp + (-22.0/105)*hjphm )*bj*bjph
                          +   ((4.0/105)*hjmhp + (13.0/70)*hjphm )*bjph*bjph ); 
    
    hbxbxuva33  = 2*idx*(     ((1.0/120)*hjmhp + (5.0/168)*hjphm )*bjmh*bjmh
                          + 2*((-1.0/70)*hjmhp + (-9.0/70)*hjphm )*bjmh*bj
                          + 2*((1.0/168)*hjmhp + (83.0/840)*hjphm )*bjmh*bjph
                          +   ((2.0/35)*hjmhp + (22.0/35)*hjphm )*bj*bj
                          + 2*((-3.0/70)*hjmhp + (-1.0/2)*hjphm )*bj*bjph
                          +   ((31.0/840)*hjmhp + (337.0/840)*hjphm )*bjph*bjph );        

    // LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbxbxuva11;  
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbxbxuva12; 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbxbxuva13; 
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbxbxuva21; 
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbxbxuva22; 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbxbxuva23; 
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbxbxuva31; 
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbxbxuva32;  
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbxbxuva33; 
    

    uais[j-1] = 0;
    
    ubis[j-1] = ubis[j-1] + LHSa21;
    ubis[j] = 0;
    
    ucis[j-1] = ucis[j-1] + LHSa11;
    ucis[j] = ucis[j] + LHSa22;
    ucis[j+1] = 1;
    
    udis[j-1] = udis[j-1] + LHSa12;
    udis[j] = udis[j] + LHSa23;
    
    ueis[j-1] = ueis[j-1] + LHSa13;
    
    
    nGis[j-1] = nGis[j-1] + Gintia11;  
    nGis[j] = nGis[j] + Gintia21;
    nGis[j+1] = ueend[0]; 

    hhbc[3*i+1] = hjmhp;
    hhbc[3*i+2] = h[i];
    hhbc[3*i+3] = hjphm; 

    Ghbc[3*i+1] = Gjmhp;
    Ghbc[3*i+2] = G[i];
    Ghbc[3*i+3] = Gjphm;

    bedhbc[3*i+1] = bjmh;
    bedhbc[3*i+2] = bj;
    bedhbc[3*i+3] = bjph;

    PENT(uais,ubis,ucis,udis,ueis,nGis,m,u);  

    free(uais);
    free(ubis);
    free(ucis);
    free(udis);
    free(ueis);
    free(nGis);

}

void getGfromuFD(double *h, double *u, double u0, double u1, double h0, double h1, double dx , int n, double *G)
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

void evolve(double *Gbc, double *hbc, double *ubc, double g, double dt, double dx, double theta,int nCells,double *nG, double *nh)
{

double idx = 1.0 / dx;


// counts what edge we are one on xubc, so we only count edges once, and we want to start on a real cell edge
int i = 1;

// counts what edge we are on by its left index
int j = 0;

// counts what cell we are in on the real cells

int k = -1;


// so we start on our cell left cell edge, want to have h,u,G anbd ux at this point.

//u is given from its list:

double ue = ubc[i];

//h and G similar, but we use different index
double hel = hbc[j];
double Gel = Gbc[j];

// we calculate derivative of u to the left at the edge (backwards)

double duel = idx*(ubc[i] - ubc[i-1]);

// now we calculate the right side of the esdge ffor h, G, du
double her = hbc[j+1];
double Ger = Gbc[j + 1];

double duer = idx*(ubc[i+1] - ubc[i]);

double sqrtghel = sqrt(g*hel);
double sqrtgher = sqrt(g*her);

double sl = fmin(0,fmin(ue - sqrtghel, ue - sqrtgher));
double sr = fmax(0,fmax(ue + sqrtghel, ue + sqrtgher));

double felh = ue*hel;
double felG = Gel*ue + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
double ferh = ue*her;
double ferG = Ger*ue + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

double isrmsl = 0.0;

if(sr != sl) isrmsl = 1.0 / (sr - sl);

double foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel));
double foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel));

double fih = foh;
double fiG = foG;

// increment our numbers, ok , so our cell and edge count only go up by 1, lets loop over cells midpoints, since we are calculating them on the way out
    i = 2;
    j = 3;
    for (k =0;k < nCells ; k++)
    {

        // so we want to calculate the flux across the i + 1/2 boundary in the kth cell, since we from our previous calculation have the flux across the i - 1/2 


        ue = ubc[i];

        //h and G similar, but we use different index
        hel = hbc[j];
        Gel = Gbc[j];

        // we calculate derivative of u to the left at the edge (backwards)

        duel = idx*(ubc[i] - ubc[i-1]);

        // now we calculate the right side of the esdge ffor h, G, du
        her = hbc[j+1];
        Ger = Gbc[j + 1];

        duer = idx*(ubc[i+1] - ubc[i]);

        sqrtghel = sqrt(g*hel);
        sqrtgher = sqrt(g*her);

        printf("h[%d] : %1.20f | G[%d] : %1.20f | u[%d] : %1.20f\n", j-1,hbc[j-1],j-1,Gbc[j-1],i,0.5*(ubc[i] + ubc[i - 1]));

        sl = fmin(0,fmin(ue - sqrtghel, ue - sqrtgher));
        sr = fmax(0,fmax(ue + sqrtghel, ue + sqrtgher));

        felh = ue*hel;
        felG = Gel*ue + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
        ferh = ue*her;
        ferG = Ger*ue + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);

        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel));
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel));

        nh[k] = hbc[j -1] -dt*idx*(foh - fih);
        nG[k] = Gbc[j -1] -dt*idx*(foG -fiG);

        fih = foh;
        fiG = foG;
        i = i +1;
        j = j + 3;

    }
    printf("%d | %d | %d \n", k,i,j);


}

void evolvewrap(double *G, double *h,double *bed,double *hebeg , double *heend ,double *Gebeg , double *Geend,double *uebeg , double *ueend, double *bedbeg, double *bedend, double g, double dx, double dt, int n, int nBCs, double theta, double *ubc, double *hhbc, double *Ghbc,double *bedhbc)
{

    //n is number of cells

    //nBCs is number of cells to define ghost cells
    getufromG(h,G,bed,hebeg,heend,Gebeg,Geend,uebeg,ueend,bedbeg,bedend,theta,dx ,n, 2*n + 1,2,3*n +2,ubc,hhbc,Ghbc,bedhbc);
    //evolve(Ghbc, hhbc, ubc,g, dt,dx, theta, n ,Gp, hp);


    //getufromG(hp,Gp, hebeg,heend,Gebeg,Geend,uebeg,ueend,theta,dx ,n, n+ 3 ,3*n +2,2, upbc,hphbc,Gphbc);
    //evolve(Gphbc, hphbc, upbc,g, dt,dx, theta, n,Gpp, hpp);

    //int i;
    //for(i=0;i<n;i++)
    //{
        //G[i] = 0.5*(G[i] + Gpp[i]);
        //h[i] = 0.5*(h[i] + hpp[i]);
        //G[i] = Gp[i];
        //h[i] = hp[i];
    //}

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
