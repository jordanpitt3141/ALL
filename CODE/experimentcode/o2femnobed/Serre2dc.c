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

void getufromG(double *h, double *G, double *hebeg, double *heend, double *Gebeg, double *Geend,double *uebeg, double *ueend, double theta, double dx , int n, int nu, int nhbc, int nubc , double *u, double *hhbc, double *Ghbc)
{

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now
    //need more B.Cs, lets make it handle any number properly, implying we have more than the necessary ones.

    //n is number of cells
    //nu is number of points in our final 'u' so its (n - 1 + 2*nubc)
    //nhbc is the number of BCs for h and G
    //nubc is number of BCs for u

    double idx = 1.0 / dx;
    double *a = malloc((n-2)*sizeof(double));
    double *b = malloc((n-1)*sizeof(double));
    double *c = malloc((n-2)*sizeof(double));

    double *d = malloc((n-1)*sizeof(double));

    double *uf = malloc((n-1)*sizeof(double));

    int i =0;

    double Gpimh, Gmiph, Gpiph,Gmip1h,hpimh,hmiph,hpiph,hmip1h;
    double dGip1b,dGip1m,dGip1f,dGip1,dhip1b, dhip1m, dhip1f,dhip1;

    Gmiph= Gebeg[nhbc-1] + 2*(G[0] - Gebeg[nhbc-1]);
    Gpimh= Gebeg[nhbc-1];

    hmiph= hebeg[nhbc-1] + 2*(h[0] - hebeg[nhbc-1]);
    hpimh= hebeg[nhbc-1];

    hhbc[0] = hebeg[nhbc-2];
    hhbc[1] = hpimh;
    hhbc[2] = h[0];
    hhbc[3] = hmiph; 

    Ghbc[0] = Gebeg[nhbc-2];
    Ghbc[1] = Gpimh;
    Ghbc[2] = G[0];
    Ghbc[3] = Gmiph; 

    //gradient of right cell

    dGip1b = (G[i+1] - G[i]);
    dGip1m = 0.5*(G[i+2] - G[i]);
    dGip1f = (G[i+2] - G[i+1]);

    dhip1b = (h[i+1] - h[i]);
    dhip1m = 0.5*(h[i+2] - h[i]);
    dhip1f = (h[i+2] - h[i+1]);

    dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);
    dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);

    Gmip1h= G[i+1] + 0.5*dGip1;
    Gpiph= G[i+1] - 0.5*dGip1;

    hmip1h= h[i+1] + 0.5*dhip1;
    hpiph= h[i+1] - 0.5*dhip1;

    b[i] = (hpimh + 3*hmiph + 3*hpiph  + hmip1h) +  idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph +  hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h);
    c[i] = (hpiph + hmip1h) - idx*idx*(hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h);

    d[i] =  2*Gpimh + 4*Gmiph + 4*Gpiph + 2*Gmip1h - uebeg[nubc - 1]*((hpimh + hmiph) -idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph) );

    
    Gmiph = Gmip1h;
    Gpimh = Gpiph;

    hmiph = hmip1h;
    hpimh = hpiph;
    

    for (i =1;i < n-2 ; i++)
    {

        //only need to reconstruct right cell
        dGip1b = (G[i+1] - G[i]);
        dGip1m = 0.5*(G[i+2] - G[i]);
        dGip1f = (G[i+2] - G[i+1]);

        dhip1b = (h[i+1] - h[i]);
        dhip1m = 0.5*(h[i+2] - h[i]);
        dhip1f = (h[i+2] - h[i+1]);

        dGip1 = minmod(theta*dGip1b, dGip1m, theta*dGip1f);
        dhip1 = minmod(theta*dhip1b, dhip1m, theta*dhip1f);

        Gmip1h= G[i+1] + 0.5*dGip1;
        Gpiph= G[i+1] - 0.5*dGip1;

        hmip1h= h[i+1] + 0.5*dhip1;
        hpiph= h[i+1] - 0.5*dhip1;

        a[i - 1] = (hpimh + hmiph) -idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph);
        b[i] = (hpimh + 3*hmiph + 3*hpiph  + hmip1h) +  idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph +  hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h);
        c[i] = (hpiph + hmip1h) - idx*idx*(hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h);

        d[i] =  2*Gpimh + 4*Gmiph + 4*Gpiph + 2*Gmip1h;

        hhbc[3*i+1] = hpimh;
        hhbc[3*i+2] = h[i];
        hhbc[3*i+3] = hmiph; 

        Ghbc[3*i+1] = Gpimh;
        Ghbc[3*i+2] = G[i];
        Ghbc[3*i+3] = Gmiph;

        
        Gmiph = Gmip1h;
        Gpimh = Gpiph;

        hmiph = hmip1h;
        hpimh = hpiph;

    }

    i = n-2;

    Gmip1h= Geend[0];
    Gpiph= Geend[0] - 2*(Geend[0] - G[i+1]);

    hmip1h= heend[0];
    hpiph= heend[0] - 2*(heend[0] - h[i+1]);

    a[i - 1] = (hpimh + hmiph) -idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph);
    b[i] = (hpimh + 3*hmiph + 3*hpiph  + hmip1h) +  idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph +  hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h);

    d[i] =  2*Gpimh + 4*Gmiph + 4*Gpiph + 2*Gmip1h - ueend[0]*((hpiph + hmip1h) - idx*idx*(hpiph*hpiph*hpiph + hpiph*hpiph*hmip1h +  hpiph*hmip1h*hmip1h  + hmip1h*hmip1h*hmip1h));

    hhbc[3*i+1] = hpimh;
    hhbc[3*i+2] = h[i];
    hhbc[3*i+3] = hmiph; 

    Ghbc[3*i+1] = Gpimh;
    Ghbc[3*i+2] = G[i];
    Ghbc[3*i+3] = Gmiph;

    hhbc[3*(i+1) +1] = hpiph;
    hhbc[3*(i+1)+2] = h[i+1];
    hhbc[3*(i+1)+3] = hmip1h; 

    hhbc[3*(i+1)+4] = heend[1]; 

    Ghbc[3*(i+1)+1] = Gpiph;
    Ghbc[3*(i+1)+2] = G[i+1];
    Ghbc[3*(i+1)+3] = Gmip1h; 

    Ghbc[3*(i+1)+4] = Geend[1]; 

    TDMA(a, b,c, d, n-1, uf);

    memcpy(u+nubc,uf,(n-1)*sizeof(double));

    for (i =0;i < nubc ; i++)
    {
        u[i] = uebeg[nubc -1 -i];
        u[nu - 1 - i] = ueend[i];

    }  

    free(a);
    free(b);
    free(c);
    free(d);
    free(uf);

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

void evolve(double *G, double *h, double *u, double g, double dx, double dt,int nBcell,int nubc,double *nG, double *nh)
{
    // G and h now contain, cell centre and edges
    //u is now just cell edges

    //they contain the virtual cells as well

    //n cells

    //modifies nh and nG to give the new values of h and G after a single time step
    double idx = 1.0 / dx;
    double ithree = 1.0 / 3.0;
    int i = 0,j = 0;

    double sqrtghel,sqrtgher,sl,sr,felh,ferh,felG,ferG,isrmsl,foh,fih,foG,fiG;
    double hir,Gir,uir,hip1l,Gip1l,uip1l,duel, duer;

    fih = 0;
    fiG = 0;
    
    //this is a boundary cell (cell = -1), so is a counter of ghost cells as well
    i = 0;
    j = nubc - 1;

    hir = h[0];
    Gir = G[0];
    uir = u[j];

    hip1l = h[1];
    Gip1l = G[1];
    uip1l = u[j];

    duel = idx*(1.5*u[j] - 2*u[j-1] + 0.5*u[j-2]);
    duer = idx*(-1.5*u[j] + 2*u[j+1] - 0.5*u[j+2]);

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

    fih = foh;
    fiG = foG; 

    j = nubc;
//so what is the starting i

    for (i = 1 ; i < nBcell-1;i++)
    {
        hir = h[3*i];
        Gir = G[3*i];
        uir = u[j];

        hip1l = h[3*i+1];
        Gip1l = G[3*i+1];
        uip1l = u[j];

        duel = idx*(1.5*u[j] - 2*u[j-1] + 0.5*u[j-2]);
        duer = idx*(-1.5*u[j] + 2*u[j+1] - 0.5*u[j+2]);
        //printf("%d | hir : %f | Gir : %f | uir : %f \n", i,h[3*i], G[3*i],u[j]);
        //printf("%d | hip1l : %f | Gip1l : %f | uip1l : %f \n", i,h[3*i +1], G[3*i+1],u[j]);
        //printf("%d | duel : %f | duer : %f \n",i, duel, duer);

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

        nh[i - 1] = h[3*i -1] -dt*idx*(foh - fih);
        nG[i - 1] = G[3*i - 1] -dt*idx*(foG -fiG);

        fih = foh;
        fiG = foG; 
        j = j + 1;
 
    }

}


void evolvewrap(double *G, double *h,double *hebeg , double *heend ,double *Gebeg , double *Geend,double *uebeg , double *ueend, double g, double dx, double dt, int n, int nhbc, int nubc, double theta)
{

    //n is number of cells

    //nhbc is size of hebeg and heend (and G as well)

    //nubc is size of uebeg


    //first allocate memory for BC variables

    //counting cells here, we have (n-2) interior, cells (alle xcept most left and most right)
    //each of the interior cells has 3 points, the edge ones have2 points each, and the remaining nhbc ones come from hebeg
    int nhBC =  (3*(n -2) + 2*2 + 2*nhbc);

    // we count the interior cells, which we solve for plus, the boundary cells we set   
    int nuBC = (n-1) + 2*nubc;


    
    double *Ghbc = malloc(nhBC*sizeof(double));
    double *hhbc = malloc(nhBC*sizeof(double));
    double *ubc = malloc(nuBC*sizeof(double));


    double *Gp = malloc(n*sizeof(double));
    double *hp = malloc(n*sizeof(double));

    getufromG(h,G, hebeg,heend,Gebeg,Geend,uebeg,ueend,theta,dx ,n, nuBC ,nhbc,nubc, ubc,hhbc,Ghbc);

    evolve(Ghbc,hhbc,ubc,g,dx,dt,n + 2,nubc,Gp, hp);

    
    getufromG(hp,Gp, hebeg,heend,Gebeg,Geend,uebeg,ueend,theta,dx ,n, nuBC ,nhbc,nubc, ubc,hhbc,Ghbc);
    evolve(Ghbc,hhbc,ubc,g,dx,dt,n + 2,nubc ,Gp, hp);

    int i;
    for(i=0;i<n;i++)
    {
        G[i] = 0.5*(G[i] + Gp[i]);
        h[i] = 0.5*(h[i] + hp[i]);
        //G[i] = Gp[i];
        //h[i] = hp[i];
    }
    /*
    int i;
    for(i=0;i<nhBC;i++)
    {
        printf("hhbc[%d] : %f | Ghbc[%d] : %f \n ", i, hhbc[i],i, Ghbc[i]);
    }

    for(i=0;i<nuBC;i++)
    {
        printf("ubc[%d] : %f  \n", i, ubc[i]);
    }
    */


    free(Ghbc);
    free(hhbc);
    free(ubc);
    free(Gp);
    free(hp);

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
