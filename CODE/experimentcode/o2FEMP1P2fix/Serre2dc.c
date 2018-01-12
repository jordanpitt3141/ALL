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

void getufromG(double *h, double *G, double *hebeg, double *heend, double *Gebeg, double *Geend,double *uebeg, double *ueend, double theta, double dx , int n, int nu, int nhbc, int nBC , double *u, double *hhbc, double *Ghbc)
{

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now

    double idx = 1.0 / dx;
    double *a = malloc((n-2)*sizeof(double));
    double *b = malloc((n-1)*sizeof(double));
    double *c = malloc((n-2)*sizeof(double));
    double *d = malloc((n-1)*sizeof(double));

    double *uf = malloc((n-1)*sizeof(double));

    int i =0;

    double Gpimh, Gmiph, Gpiph,Gmip1h,hpimh,hmiph,hpiph,hmip1h;
    double dGip1b,dGip1m,dGip1f,dGip1,dhip1b, dhip1m, dhip1f,dhip1;
    double dGib,dGim,dGif,dGi,dhib, dhim, dhif,dhi;

    Gmiph= Gebeg[nBC-1] + 2*(G[0] - Gebeg[nBC-1]);
    Gpimh= Gebeg[nBC-1];

    hmiph= hebeg[nBC-1] + 2*(h[0] - hebeg[nBC-1]);
    hpimh= hebeg[nBC-1];

    hhbc[0] = hebeg[nBC-2];
    hhbc[1] = hpimh;
    hhbc[2] = h[0];
    hhbc[3] = hmiph; 

    Ghbc[0] = Gebeg[nBC-2];
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

    d[i] =  2*Gpimh + 4*Gmiph + 4*Gpiph + 2*Gmip1h - uebeg[nBC - 1]*((hpimh + hmiph) -idx*idx*(hpimh*hpimh*hpimh + hpimh*hpimh*hmiph +  hpimh*hmiph*hmiph  + hmiph*hmiph*hmiph) );
    

    for (i =1;i < n-2 ; i++)
    {

        dGib = (G[i] - G[i-1]);
        dGim = 0.5*(G[i+1] - G[i-1]);
        dGif = (G[i+1] - G[i]);

        dhib = (h[i] - h[i-1]);
        dhim = 0.5*(h[i+1] - h[i-1]);
        dhif = (h[i+1] - h[i]);

        dGi = minmod(theta*dGib, dGim, theta*dGif);
        dhi = minmod(theta*dhib, dhim, theta*dhif);

        Gmiph= G[i] + 0.5*dGi;
        Gpimh= G[i] - 0.5*dGi;

        hmiph= h[i] + 0.5*dhi;
        hpimh= h[i] - 0.5*dhi;

        //reconstruct both cells
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

        


    }

    i = n-2;

    dGib = (G[i] - G[i-1]);
    dGim = 0.5*(G[i+1] - G[i-1]);
    dGif = (G[i+1] - G[i]);

    dhib = (h[i] - h[i-1]);
    dhim = 0.5*(h[i+1] - h[i-1]);
    dhif = (h[i+1] - h[i]);

    dGi = minmod(theta*dGib, dGim, theta*dGif);
    dhi = minmod(theta*dhib, dhim, theta*dhif);

    Gmiph= G[i] + 0.5*dGi;
    Gpimh= G[i] - 0.5*dGi;

    hmiph= h[i] + 0.5*dhi;
    hpimh= h[i] - 0.5*dhi;

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

    memcpy(u+2,uf,(n-1)*sizeof(double));
    u[0] = uebeg[nBC-2];
    u[1] = uebeg[nBC-1];
    u[nu -2] = ueend[0];
    u[nu -1] = ueend[1];   

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

void evolvewrap(double *G, double *h,double *hebeg , double *heend ,double *Gebeg , double *Geend,double *uebeg , double *ueend, double g, double dx, double dt, int n, int nBCs, double theta, double *ubc, double *hhbc, double *Ghbc)
{

    //n is number of cells

    //nBCs is number of cells to define ghost cells
    double *Gphbc = malloc((3*n +2)*sizeof(double));
    double *hphbc = malloc((3*n +2)*sizeof(double));
    double *upbc = malloc(((n+ 3))*sizeof(double));

    double *Gp = malloc(n*sizeof(double));
    double *hp = malloc(n*sizeof(double));

    double *Gpp = malloc(n*sizeof(double));
    double *hpp = malloc(n*sizeof(double));

    
    getufromG(h,G, hebeg,heend,Gebeg,Geend,uebeg,ueend,theta,dx ,n, n+ 3 ,3*n +2,2, ubc,hhbc,Ghbc);
    evolve(Ghbc, hhbc, ubc,g, dt,dx, theta, n ,Gp, hp);


    //getufromG(hp,Gp, hebeg,heend,Gebeg,Geend,uebeg,ueend,theta,dx ,n, n+ 3 ,3*n +2,2, upbc,hphbc,Gphbc);
    //evolve(Gphbc, hphbc, upbc,g, dt,dx, theta, n,Gpp, hpp);

    int i;
    for(i=0;i<n;i++)
    {
        //G[i] = 0.5*(G[i] + Gpp[i]);
        //h[i] = 0.5*(h[i] + hpp[i]);
        G[i] = Gp[i];
        h[i] = hp[i];
    }


    free(Gp);
    free(hp);
    free(Gpp);
    free(hpp);
    free(Gphbc);
    free(hphbc);
    free(upbc);
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
