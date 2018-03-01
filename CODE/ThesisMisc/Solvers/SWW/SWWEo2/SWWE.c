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
    double idx = 1.0*idx;

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
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
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
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    double tgpux = interpquarticgrad(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
}

double hacrosscell(double *x,double *h,int j,double dx)
{
    //so we have h,u at midpoints
    //epsilon and sigma are everywhere

	double *hcoeff = malloc(5*sizeof(double));
	

    //jth cell
    interpquartcoeff(h,hcoeff,j,dx);
    
    //first gauss point
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    
    double fgpe = fgph;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    
    double sgpe = sgph;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    
    double tgpe = tgph;

	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
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
    double fgp = 0.5*dx*sqrt(3.0*i5) + x[j];
    double fgph = interpquarticval(hcoeff,x[j],fgp);
    double fgpu = interpquarticval(ucoeff,x[j],fgp);
    
    double fgpe = fgph*fgpu;
        
    //second gauss point
    double sgp = x[j];
    double sgph = interpquarticval(hcoeff,x[j],sgp);
    double sgpu = interpquarticval(ucoeff,x[j],sgp);
    
    double sgpe = sgph*sgpu;

    //third gauss point
    double tgp = -0.5*dx*sqrt(3.0*i5) + x[j];
    double tgph = interpquarticval(hcoeff,x[j],tgp);
    double tgpu = interpquarticval(ucoeff,x[j],tgp);
    
    double tgpe = tgph*tgpu;

	free(ucoeff);
	free(hcoeff);
    
    return 0.5*dx*( (5.0*i9)*fgpe + (8.0*i9)*sgpe + (5.0*i9)*tgpe);
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

void calculateedgesIncom(double *h, double *uh, double *bed, double WG1ht, double *hMend, double *wMend, double *uhMend, double *bMend, double theta, double dx , int n, int nBC, int nbc, double *uhbc, double *hhbc, double *whbc, double *bhbc)
{
    //trying to join different B.C systems...., also try incorporating into  C code that already has beds, like o2bedfix.

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now
    int i;
    double idx = 1.0/ dx;

    double duhib,duhim,duhif,dhib,dhim,dhif,duhi,dhi,uhjphm,uhjmhp,hjphm,hjmhp,bjphm,bjmhp;

    double wim1, wi,wip1,dwim,dwif,dwi,dwib,wjmhp,wjphm;

    for (i =1;i < n-1 ; i++)
    {
        // Reconstruct G and h
        dhib = (h[i] - h[i-1]);
        dhim = 0.5*(h[i+1] - h[i-1]);
        dhif = (h[i+1] - h[i]);

        duhib = (uh[i] - uh[i-1]);
        duhim = 0.5*(uh[i+1] - uh[i-1]);
        duhif = (uh[i+1] - uh[i]);

        wi = h[i] + bed[i];
        wip1 = h[i+1] + bed[i+1];
        wim1 = h[i-1] + bed[i-1];
        dwib = (wi - wim1);
        dwim = 0.5*(wip1 - wim1);
        dwif = (wip1 - wi);

        duhi = minmod(theta*duhib, duhim, theta*duhif);
        dhi = minmod(theta*dhib, dhim, theta*dhif);
        dwi = minmod(theta*dwib, dwim, theta*dwif);

        uhjphm= uh[i] + 0.5*duhi;
        uhjmhp= uh[i] - 0.5*duhi;

        hjphm= h[i] + 0.5*dhi;
        hjmhp= h[i] - 0.5*dhi;   

        wjphm= wi + 0.5*dwi;
        wjmhp= wi - 0.5*dwi; 

        bjphm = wjphm - hjphm;
        bjmhp = wjmhp - hjmhp;  

 


        hhbc[3*i + (nBC)] = hjmhp;
        hhbc[3*i+1 + (nBC)] = h[i];
        hhbc[3*i+2 + (nBC)] = hjphm; 

        whbc[3*i + (nBC)] = wjmhp;
        whbc[3*i+1 + (nBC)] = wi;
        whbc[3*i+2 + (nBC)] = wjphm;

        uhbc[3*i + (nBC)] = uhjmhp;
        uhbc[3*i+1 + (nBC)] = uh[i];
        uhbc[3*i+2 + (nBC)] = uhjphm;

        bhbc[3*i + (nBC)] = bjmhp;
        bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp + bjphm);
        bhbc[3*i+2 + (nBC)] = bjphm;     


    }


    // first 
    i = 0;
    // Get it from interior
    // Reconstruct w
    wi = h[i] + bed[i];

    // Reconstruct G and h

    hjphm= hhbc[nBC +3];
    hjmhp= WG1ht;

    wjphm=  hjphm + bed[i];
    wjmhp=  hjmhp + bed[i];   

    bjphm = wjphm - hjphm;
    bjmhp = wjmhp - hjmhp;
    

    hhbc[3*i + (nBC)] = hjmhp;
    hhbc[3*i+1 + (nBC)] = 0.5*(hjmhp +hjphm);
    hhbc[3*i+2 + (nBC)] = hjphm; 

    whbc[3*i + (nBC)] = wjmhp;
    whbc[3*i+1 + (nBC)] = 0.5*(wjmhp +wjphm);
    whbc[3*i+2 + (nBC)] = wjphm;
 
    bhbc[3*i + (nBC)] = bjmhp;
    bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp +bjphm);
    bhbc[3*i+2 + (nBC)] = bjphm; 


    double uhgrad = idx*(uhbc[7] - uhbc[4]);
    uhbc[1] = 2*uh[0] - uhbc[7];
    uhbc[0] = -0.5*dx*(uhgrad) + uhbc[1];
    uhbc[2] = 0.5*dx*(uhgrad) + uhbc[1];
    uhbc[3] = uhbc[2];
    uhbc[4] = uh[0];
    uhbc[5] = 0.5*(uh[0] + uh[1]); 
     


// last
    i = n-1;

    // Get it from interior
    // Reconstruct w
    wi = h[i] + bed[i];
 
    wjphm= wMend[0];
    wjmhp= whbc[nbc - nBC -4]; 

    uhjphm= uhMend[0];
    uhjmhp= uhbc[nbc - nBC -4];

    hjphm= hMend[0];
    hjmhp= hhbc[nbc - nBC -4];

    bjphm = wjphm - hjphm;
    bjmhp = wjmhp - hjmhp;

    hhbc[3*i + (nBC)] = hjmhp;
    hhbc[3*i+1 + (nBC)] = 0.5*(hjmhp +hjphm);
    hhbc[3*i+2 + (nBC)] = hjphm; 

    whbc[3*i + (nBC)] = wjmhp;
    whbc[3*i+1 + (nBC)] = 0.5*(wjmhp +wjphm);
    whbc[3*i+2 + (nBC)] = wjphm;

    uhbc[3*i + (nBC)] = uhjmhp;
    uhbc[3*i+1 + (nBC)] = 0.5*(uhjmhp +uhjphm);
    uhbc[3*i+2 + (nBC)] = uhjphm;

    bhbc[3*i + (nBC)] = bjmhp;
    bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp +bjphm);
    bhbc[3*i+2 + (nBC)] = bjphm;

    // B.C stuff
    for(i=0;i < nBC;i++)
    {

        //ending
        whbc[nbc-nBC + i] = wMend[i];
        hhbc[nbc-nBC + i] = hMend[i];
        uhbc[nbc-nBC +i] = uhMend[i];
        bhbc[nbc-nBC +i] = bMend[i];

        //beginning
        whbc[(nBC -1) - i] = 2*whbc[(nBC -1) - i + 1] - whbc[(nBC -1) - i + 2];
        hhbc[(nBC -1) - i] = 2*hhbc[(nBC -1) - i + 1] - hhbc[(nBC -1) - i + 2];
        bhbc[(nBC -1) - i] = 2*bhbc[(nBC -1) - i + 1] - bhbc[(nBC -1) - i + 2];
    }


}


void calculateedges(double *h, double *uh, double *bed, double *hMbeg, double *hMend, double *wMbeg, double *wMend,double *uhMbeg, double *uhMend, double *bMbeg, double *bMend, double theta, double dx , int n, int nBC, int nbc, double *uhbc, double *hhbc, double *whbc, double *bhbc)
{
    //trying to join different B.C systems...., also try incorporating into  C code that already has beds, like o2bedfix.

    // n is length of h, G
    //nu is length of u (should be n + 1, to count edges of n cells)

    //enforcing B.Cs at cell edges now
    int i;

    double duhib,duhim,duhif,dhib,dhim,dhif,duhi,dhi,uhjphm,uhjmhp,hjphm,hjmhp,bjphm,bjmhp;

    double wim1, wi,wip1,dwim,dwif,dwi,dwib,wjmhp,wjphm;

    for (i =1;i < n-1 ; i++)
    {
        // Reconstruct G and h
        dhib = (h[i] - h[i-1]);
        dhim = 0.5*(h[i+1] - h[i-1]);
        dhif = (h[i+1] - h[i]);

        duhib = (uh[i] - uh[i-1]);
        duhim = 0.5*(uh[i+1] - uh[i-1]);
        duhif = (uh[i+1] - uh[i]);

        wi = h[i] + bed[i];
        wip1 = h[i+1] + bed[i+1];
        wim1 = h[i-1] + bed[i-1];
        dwib = (wi - wim1);
        dwim = 0.5*(wip1 - wim1);
        dwif = (wip1 - wi);

        duhi = minmod(theta*duhib, duhim, theta*duhif);
        dhi = minmod(theta*dhib, dhim, theta*dhif);
        dwi = minmod(theta*dwib, dwim, theta*dwif);

        uhjphm= uh[i] + 0.5*duhi;
        uhjmhp= uh[i] - 0.5*duhi;

        hjphm= h[i] + 0.5*dhi;
        hjmhp= h[i] - 0.5*dhi;   

        wjphm= wi + 0.5*dwi;
        wjmhp= wi - 0.5*dwi; 

        bjphm = wjphm - hjphm;
        bjmhp = wjmhp - hjmhp;  

 


        hhbc[3*i + (nBC)] = hjmhp;
        hhbc[3*i+1 + (nBC)] = h[i];
        hhbc[3*i+2 + (nBC)] = hjphm; 

        whbc[3*i + (nBC)] = wjmhp;
        whbc[3*i+1 + (nBC)] = wi;
        whbc[3*i+2 + (nBC)] = wjphm;

        uhbc[3*i + (nBC)] = uhjmhp;
        uhbc[3*i+1 + (nBC)] = uh[i];
        uhbc[3*i+2 + (nBC)] = uhjphm;

        bhbc[3*i + (nBC)] = bjmhp;
        bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp + bjphm);
        bhbc[3*i+2 + (nBC)] = bjphm;     


    }


    // first 
    i = 0;
    // Get it from interior
    // Reconstruct w
    wi = h[i] + bed[i];
 
    wjphm= whbc[nBC +3];
    wjmhp= wMbeg[nBC-1];  

    // Reconstruct G and h

    hjphm= hhbc[nBC +3];
    hjmhp= hMbeg[nBC-1];

    uhjphm= uhbc[nBC +3];
    uhjmhp= uhMbeg[nBC-1];

    bjphm = wjphm - hjphm;
    bjmhp = wjmhp - hjmhp;
    

    hhbc[3*i + (nBC)] = hjmhp;
    hhbc[3*i+1 + (nBC)] = 0.5*(hjmhp +hjphm);
    hhbc[3*i+2 + (nBC)] = hjphm; 

    whbc[3*i + (nBC)] = wjmhp;
    whbc[3*i+1 + (nBC)] = 0.5*(wjmhp +wjphm);
    whbc[3*i+2 + (nBC)] = wjphm;

    uhbc[3*i + (nBC)] = uhjmhp;
    uhbc[3*i+1 + (nBC)] = 0.5*(uhjmhp +uhjphm);
    uhbc[3*i+2 + (nBC)] = uhjphm;
 
    bhbc[3*i + (nBC)] = bjmhp;
    bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp +bjphm);
    bhbc[3*i+2 + (nBC)] = bjphm; 

// last
    i = n-1;

    // Get it from interior
    // Reconstruct w
    wi = h[i] + bed[i];
 
    wjphm= wMend[0];
    wjmhp= whbc[nbc - nBC -4]; 

    uhjphm= uhMend[0];
    uhjmhp= uhbc[nbc - nBC -4];

    hjphm= hMend[0];
    hjmhp= hhbc[nbc - nBC -4];

    bjphm = wjphm - hjphm;
    bjmhp = wjmhp - hjmhp;

    hhbc[3*i + (nBC)] = hjmhp;
    hhbc[3*i+1 + (nBC)] = 0.5*(hjmhp +hjphm);
    hhbc[3*i+2 + (nBC)] = hjphm; 

    whbc[3*i + (nBC)] = wjmhp;
    whbc[3*i+1 + (nBC)] = 0.5*(wjmhp +wjphm);
    whbc[3*i+2 + (nBC)] = wjphm;

    uhbc[3*i + (nBC)] = uhjmhp;
    uhbc[3*i+1 + (nBC)] = 0.5*(uhjmhp +uhjphm);
    uhbc[3*i+2 + (nBC)] = uhjphm;

    bhbc[3*i + (nBC)] = bjmhp;
    bhbc[3*i+1 + (nBC)] = 0.5*(bjmhp +bjphm);
    bhbc[3*i+2 + (nBC)] = bjphm;

    // B.C stuff
    for(i=0;i < nBC;i++)
    {
        whbc[i] = wMbeg[i];
        whbc[nbc-nBC + i] = wMend[i];
        hhbc[i] = hMbeg[i];
        hhbc[nbc-nBC + i] = hMend[i];
        uhbc[i] = uhMbeg[i];
        uhbc[nbc-nBC +i] = uhMend[i];
        bhbc[i] = bMbeg[i];
        bhbc[nbc-nBC +i] = bMend[i];

    }


}

void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
    //replace with memcpy for performance?
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}



//include BCs
void evolve(double *uhbc, double *hhbc, double *whbc,double *bhbc, int nbc, int nBC, double g, double dx, double dt, int n, double theta, double *newh, double *newuh)
{
    double idx = 1.0 / dx;  
	int i;
    double her,uer,hel,uel,fhel,fher,fuhel,fuher,sqrtghel,sqrtgher,sl,sr,isrmsl,foh,fouh,fih,fiuh,th,tbx,sourcer,sourcel,sourcec;
	double wir,wip1l,bip1l,bil,bir,nbi,hihm,hihp;
    double himhp;

    // i = -1
    i = -1;

    //Figure out how the interior numberings work.
    
    //wil = whbc[3*i + nBC];
    wir = whbc[3*i + nBC + 2];

    // bil and bir
    bil = bhbc[3*i + nBC];
    bir = bhbc[3*i + nBC + 2];

    wip1l = whbc[3*(i+1) + nBC];
    bip1l = bhbc[3*(i+1) + nBC];

    // new bed and height
    nbi = fmax(bip1l, bir);
	hihm = fmax(0, wir - nbi);
	hihp = fmax(0, wip1l - nbi);


    her = hihp;
    uer  = uhbc[3*(i+1) + nBC] / hihp;

    hel = hihm;
    uel  = uhbc[3*(i) +nBC +2] / hihm;

    fhel = uhbc[3*(i) +nBC +2];
    fher = uhbc[3*(i+1) + nBC];

    fuhel = uel*fhel + 0.5*g*hel*hel;
    fuher = uer*fher + 0.5*g*her*her;

    //printf("%d | %e | %e \n",i,fGel,fGer);

    sqrtghel = sqrt(g* hel);
    sqrtgher = sqrt(g* her);

    sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

    isrmsl = 0.0;

    if(sr != sl) isrmsl = 1.0 / (sr - sl);	

    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
    fouh =isrmsl*( sr*fuhel - sl*fuher + sl*sr*(fher -fhel));

    fih = foh;
    fiuh = fouh;
    himhp = hihp;

    for(i = 0;i < n;i++)
    {
        //Figure out how the interior numberings work.
        
	    //wil = whbc[3*i + nBC];
	    wir = whbc[3*i + nBC + 2];

	    // bil and bir
	    bil = bhbc[3*i + nBC];
	    bir = bhbc[3*i + nBC + 2];

	    wip1l = whbc[3*(i+1) + nBC];
	    bip1l = bhbc[3*(i+1) + nBC];

        // new bed and height
        nbi = fmax(bip1l, bir);
		hihm = fmax(0, wir - nbi);
		hihp = fmax(0, wip1l - nbi);

        her = hihp;
        uer  = uhbc[3*(i+1) + nBC]/ hihp;

        hel = hihm;
        uel  = uhbc[3*(i) +nBC +2] / hihm;


	    fhel = uhbc[3*(i) +nBC +2];
	    fher = uhbc[3*(i+1) + nBC];

	
	    fuhel = uel*fhel + 0.5*g*hel*hel;
	    fuher = uer*fher + 0.5*g*her*her;

        //printf("%d | %e | %e \n",i,fGel,fGer);

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl) isrmsl = 1.0 / (sr - sl);	

	    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
	    fouh =isrmsl*( sr*fuhel - sl*fuher + sl*sr*(fher - fhel));



        //centerted values
        th = hhbc[3*i + nBC+ 1];
		tbx = idx*(bir - bil);
		
		sourcer = g*0.5*(hihm*hihm - hhbc[3*(i) +nBC +2]*hhbc[3*(i) +nBC +2]);
		sourcec = -g*th*tbx;
		sourcel = g*0.5*(hhbc[3*(i) +nBC]*hhbc[3*(i) +nBC] - himhp*himhp);

		newh[i] = hhbc[3*(i) +nBC + 1] - dt*idx*(foh - fih);
		newuh[i]= uhbc[3*(i) +nBC + 1] - dt*idx*(fouh - fiuh) + dt*idx*(sourcer+sourcel + dx*sourcec);


	    fih = foh;
	    fiuh = fouh;
        himhp = hihp;

    }

}



void evolvewrap(double *uh, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uhMbeg , double *uhMend, double *bMbeg, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc)
{

    //n is number of cells
    double *hp = malloc((n)*sizeof(double));
    double *uhp = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *uhpp = malloc((n)*sizeof(double));

    //nBCs is number of cells to define ghost cells
    calculateedges(h,uh,bed,hMbeg,hMend,wMbeg,wMend,uhMbeg,uhMend,bMbeg,bMend,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hp,uhp);
    calculateedges(hp,uhp,bed,hMbeg,hMend,wMbeg,wMend,uhMbeg,uhMend,bMbeg,bMend,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hpp,uhpp);
    int i;
    for(i=0;i<n;i++)
    {
        uh[i] = 0.5*(uh[i] + uhpp[i]);
        h[i] = 0.5*(h[i] + hpp[i]);
        //printf("hpp | %d | %f \n",i,hpp[i]);
        //printf("hp  | %d | %f \n",i,hp[i]);
        //G[i] = Gp[i];
        //h[i] = hp[i];
        
    }

    free(hp);
    free(uhp);
    free(uhpp);
    free(hpp);
}

void evolvewrapChangeBC(double *u, double *h,double *bed,double *hMbeg , double *hMend,double *wMbeg , double *wMend,double *uMbeg , double *uMend, double *bMbeg, double *bMend,double *hMbeg1, double *hMend1,double *wMbeg1, double *wMend1,double *uMbeg1, double *uMend1, double *bMbeg1, double *bMend1, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc)
{

    //n is number of cells
    double *hp = malloc((n)*sizeof(double));
    double *up = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *upp = malloc((n)*sizeof(double));

    //nBCs is number of cells to define ghost cells
    calculateedges(h,u,bed,hMbeg,hMend,wMbeg,wMend,uMbeg,uMend,bMbeg,bMend,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hp,up);
    calculateedges(hp,up,bed,hMbeg1,hMend1,wMbeg1,wMend1,uMbeg1,uMend1,bMbeg1,bMend1,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hpp,upp);
    int i;
    for(i=0;i<n;i++)
    {
        u[i] = 0.5*(u[i] + upp[i]);
        h[i] = 0.5*(h[i] + hpp[i]);
        //printf("hpp | %d | %f \n",i,hpp[i]);
        //printf("hp  | %d | %f \n",i,hp[i]);
        //G[i] = Gp[i];
        //h[i] = hp[i];
        
    }

    free(hp);
    free(up);
    free(upp);
    free(hpp);
}

void evolvewrapIncomwavetoDir(double *u, double *h,double *bed,double WG1ht , double WG1ht1 , double *hMend, double *wMend , double *uMend, double *bMend, double g, double dx, double dt, int n, int nBC, int nbc, double theta, double *hhbc, double *whbc,double *bhbc,double *uhbc)
{

    //n is number of cells
    double *hp = malloc((n)*sizeof(double));
    double *up = malloc((n)*sizeof(double));
    double *hpp = malloc((n)*sizeof(double));
    double *upp = malloc((n)*sizeof(double));

    //nBCs is number of cells to define ghost cells
    calculateedgesIncom(h,u,bed,WG1ht,hMend,wMend,uMend,bMend,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hp,up);
    calculateedgesIncom(hp,up,bed,WG1ht1,hMend,wMend,uMend,bMend,theta,dx,n,nBC,nbc,uhbc,hhbc,whbc,bhbc);
    evolve(uhbc,hhbc,whbc,bhbc,nbc,nBC,g,dx,dt,n,theta,hpp,upp);
    int i;
    for(i=0;i<n;i++)
    {
        u[i] = 0.5*(u[i] + upp[i]);
        h[i] = 0.5*(h[i] + hpp[i]);
        //printf("hpp | %d | %f \n",i,hpp[i]);
        //printf("hp  | %d | %f \n",i,hp[i]);
        //G[i] = Gp[i];
        //h[i] = hp[i];
        
    }

    free(hp);
    free(up);
    free(upp);
    free(hpp);
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
