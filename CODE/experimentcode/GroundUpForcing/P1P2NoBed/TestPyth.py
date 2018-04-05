from scipy import *
from pylab import plot

def minmod(a,b,c):
    if a>0 and b > 0 and c>0:
        return min(a,b,c)
    elif a <0 and b< 0 and c<0:
        return max(a,b,c)
    else:
        return 0

def hGDef(x,t,a,b,c,d,e):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = a + sin(b*x[i])*exp(c*t)
        u = cos(d*x[i])*exp(e*t)
        uxi = -d*exp(e*t)*sin(d*x[i]) 
        hxi = b*exp(c*t)*cos(b*x[i])
        uxxi = -d**2*exp(e*t)*cos(d*x[i])
        G[i] = u*h[i] - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
    
    return h,G

def uDef(x,t,a,b):
    n = len(x)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        u[i] = cos(a*x[i])*exp(b*t)
    
    return u       


def Reconstruction(q,qbeg,qend,theta,FluxGhostCells):
    qext = concatenate((qbeg,q,qend))
    n1 = len(qbeg)
    n = len(q)
    n2 = len(qext)
    qR = []
    for i in range(n1 - FluxGhostCells ,n + n1 + FluxGhostCells):
        qib = (qext[i] - qext[i-1])
        qim = 0.5*(qext[i+1] - qext[i-1])
        qif = (qext[i+1] - qext[i])
        
        dqi = minmod(theta*qib,qim,theta*qif)
        
        qimhp = qext[i] - 0.5*dqi
        qiphm = qext[i] + 0.5*dqi
        
        qR.append(qimhp)
        qR.append(qext[i])
        qR.append(qiphm)
    return qR

def FluxCalc(hR,GR,uR,FluxGhostCells,n,dt,dx):
    idx = 1.0/dx
    #i includes Ghost cells
    for i in range(FluxGhostCells, n + FluxGhostCells):
        
        #hi = hR[3*i + 1]
        #Gi = GR[3*i + 1]
        #ui = uR[2*i + 1]
        
        #FluxIn
        uai = 2*idx*idx*(uR[2*i] - 2*uR[2*i + 1] + uR[2*(i+1)])
        ubi = idx*(-uR[2*i] + uR[2*(i+1)])
        
        uaim1 = 2*idx*idx*(uR[2*(i-1)] - 2*uR[2*(i-1) + 1] + uR[2*(i)])
        ubim1 = idx*(-uR[2*(i-1)] + uR[2*(i)])
        
        her = hR[3*i]
        Ger = GR[3*i]
        uer = uR[2*i]
        duer = -uai*(dx) + ubi;
        
        hel = hR[3*i-1]
        Gel = GR[3*i-1]
        uel = uR[2*i]
        duel = uaim1*(dx) + ubim1;
        
        fhel = uel*hel;
        fher = uer*her;

        fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
        fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl):
            isrmsl = 1.0 / (sr - sl);	

        fih =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
        fiG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));
        
        #Flux out
        
        uaip1 = 2*idx*idx*(uR[2*(i+1)] - 2*uR[2*(i+1) + 1] + uR[2*(i+2)])
        ubip1 = idx*(-uR[2*(i+1)] + uR[2*(i+2)])
        
        her = hR[3*(i+1)]
        Ger = GR[3*(i+1)]
        uer = uR[2*(i+1)]
        duer = -uaip1*(dx) + ubip1;
        
        hel = hR[3*(i+1)-1]
        Gel = GR[3*(i+1)-1]
        uel = uR[2*(i+1)]
        duel = uai*(dx) + ubi;
        
        fhel = uel*hel;
        fher = uer*her;

        fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel;
        fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer;

        sqrtghel = sqrt(g* hel);
        sqrtgher = sqrt(g* her);

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher));
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher));

        isrmsl = 0.0;

        if(sr != sl):
            isrmsl = 1.0 / (sr - sl);	

        foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel));
        foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel));
        
        nh.append(hR[3*i + 1] - dt*idx*(foh - fih))
        nG.append(GR[3*i + 1] - dt*idx*(foG - fiG))
    return nh,nG


g = 9.81

FluxGhostCells = 1
HalfWidth = 2
dx = 0.01 
dt = 0.1 
t = 0      
x = arange(-HalfWidth,HalfWidth+0.1*dx,dx) 
xbeg = array([x[0] - 2*dx,x[0] - dx])      
xend = array([x[-1] + dx,x[-1] + 2*dx]) 

xu = arange(-HalfWidth - (FluxGhostCells+0.5)*dx ,HalfWidth + (FluxGhostCells+0.5)*dx  +0.1*dx,0.5*dx) 

xbc =  concatenate((xbeg,x,xend))  

n = len(x)

a0 = 10
a1 = 2
a2 = 0.5
a3 = 3
a4 = 0.7

theta = 2

h,G = hGDef(x,t,a0,a1,a2,a3,a4)
u = uDef(x,t,a3,a4)
hbeg,Gbeg = hGDef(xbeg,t,a0,a1,a2,a3,a4)
hend,Gend = hGDef(xend,t,a0,a1,a2,a3,a4)


hR = Reconstruction(h,hbeg,hend,theta,FluxGhostCells)
GR = Reconstruction(G,Gbeg,Gend,theta,FluxGhostCells)

uR = uDef(xu,t,a3,a4)

FluxCalc(hR,GR,uR,FluxGhostCells,n,dt,dx)