from scipy import *
from numpy.linalg import norm
from pylab import plot
import os

i3 = 1.0/ 3.0


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


def hGSol(x,t,a0,a1,g):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        phi = x[i] - c*t;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        tanhkphi = sechkphi*((exp(k*phi) - exp(-k*phi))/2.0)
        hdx = -2*a1*k*tanhkphi*sechkphi*sechkphi
        hdxx = a1*(4*k*k*tanhkphi*tanhkphi*sechkphi*sechkphi - 2*k*k*sechkphi*sechkphi*sechkphi*sechkphi)
        h[i] = a0 + a1*sechkphi*sechkphi
        u =  c* ((h[i] - a0) / h[i])
        G[i] = u*h[i] - i3*h[i]*h[i]*h[i]*(a0*c*(h[i]*hdxx - 2*hdx*hdx)/(h[i]*h[i]*h[i])) - h[i]*h[i]*hdx*(a0*c*hdx/(h[i]*h[i]))
         
    
    return h,G

def uSol(x,t,a0,a1,g):
    n = len(x)
    u = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    
    for i in range(n):
        phi = x[i] - c*t;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        h = a0 + a1*sechkphi*sechkphi
        u[i] =  c* ((h - a0) / h)
         
    
    
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

def FluxCalc(hR,GR,uR,g,FluxGhostCells,n,dt,dx, nh, nG):
    idx = 1.0/dx
    #i includes Ghost cells
    for i in range(FluxGhostCells, n + FluxGhostCells):
        
        #hi = hR[3*i + 1]
        #Gi = GR[3*i + 1]
        #ui = uR[2*i + 1]
        
        #print(hi,Gi,ui)
        
        #FluxIn
        uai = 2*idx*idx*(uR[2*i] - 2*uR[2*i + 1] + uR[2*i + 2])
        ubi = idx*(-uR[2*i] + uR[2*i + 2])
        
        uaim1 = 2*idx*idx*(uR[2*(i-1)] - 2*uR[2*(i-1) + 1] + uR[2*(i-1) + 2])
        ubim1 = idx*(-uR[2*(i-1)] + uR[2*(i-1) + 2])
        
        herL = hR[3*i]
        GerL = GR[3*i]
        uerL = uR[2*i]
        duerL = -uai*(dx) + ubi;
        
        helL = hR[3*i-1]
        GelL = GR[3*i-1]
        uelL = uR[2*i]
        duelL = uaim1*(dx) + ubim1;
        
        fhelL = uelL*helL;
        fherL = uerL*herL;

        fGelL = GelL*uelL + 0.5*g*helL*helL - 2*i3*helL*helL*helL*duelL*duelL;
        fGerL = GerL*uerL + 0.5*g*herL*herL - 2*i3*herL*herL*herL*duerL*duerL;

        sqrtghelL = sqrt(g*helL);
        sqrtgherL = sqrt(g*herL);

        slL = min(0,uelL - sqrtghelL, uerL - sqrtgherL);
        srL = max(0,uelL + sqrtghelL, uerL + sqrtgherL);

        isrmslL = 0.0;

        if(srL != slL):
            isrmslL = 1.0 / (srL - slL);	

        fih =isrmslL*( srL*fhelL - slL*fherL + slL*srL*(herL - helL));
        fiG =isrmslL*( srL*fGelL - slL*fGerL + slL*srL*(GerL - GelL));
        
        #Flux out
        
        uaip1 = 2*idx*idx*(uR[2*(i+1)] - 2*uR[2*(i+1) + 1] + uR[2*(i+1) + 2])
        ubip1 = idx*(-uR[2*(i+1)] + uR[2*(i+1) + 2])
        
        herR = hR[3*(i+1)]
        GerR = GR[3*(i+1)]
        uerR = uR[2*(i+1)]
        duerR = -uaip1*(dx) + ubip1;
        
        helR = hR[3*(i+1)-1]
        GelR = GR[3*(i+1)-1]
        uelR = uR[2*(i+1)]
        duelR = uai*(dx) + ubi;
        
        fhelR = uelR*helR;
        fherR = uerR*herR;

        fGelR = GelR*uelR + 0.5*g*helR*helR - 2*i3*helR*helR*helR*duelR*duelR;
        fGerR = GerR*uerR + 0.5*g*herR*herR - 2*i3*herR*herR*herR*duerR*duerR;

        sqrtghelR = sqrt(g* helR);
        sqrtgherR = sqrt(g* herR);

        slR = min(0,uelR - sqrtghelR, uerR - sqrtgherR);
        srR = max(0,uelR + sqrtghelR, uerR + sqrtgherR);

        isrmslR = 0.0;

        if(srR != slR):
            isrmslR = 1.0 / (srR - slR);	

        foh =isrmslR*( srR*fhelR - slR*fherR + slR*srR*(herR - helR));
        foG =isrmslR*( srR*fGelR - slR*fGerR + slR*srR*(GerR - GelR));
        
        nh[i - FluxGhostCells] = hR[3*i + 1] - dt*idx*(foh - fih)
        nG[i - FluxGhostCells] = GR[3*i + 1] - dt*idx*(foG - fiG)

def RKStep(h,G,un,unp1,hbeg,hend,Gbeg,Gend,hbeg1,hend1,Gbeg1,Gend1,g,theta,FluxGhostCells,dx,dt,nh,nG):
    n = len(h)
    hp = zeros(n)
    Gp = zeros(n)    
    hpp = zeros(n)
    Gpp = zeros(n) 
    
    hR = Reconstruction(h,hbeg,hend,theta,FluxGhostCells)
    GR = Reconstruction(G,Gbeg,Gend,theta,FluxGhostCells)
    
    FluxCalc(hR,GR,un,g,FluxGhostCells,n,dt,dx,hp,Gp)
 
    hpR = Reconstruction(hp,hbeg1,hend1,theta,FluxGhostCells)
    GpR = Reconstruction(Gp,Gbeg1,Gend1,theta,FluxGhostCells)   
    
    FluxCalc(hpR,GpR,unp1,g,FluxGhostCells,n,dt,dx,hpp,Gpp)
    
    for i in range(n):
        nh[i] = 0.5*(h[i] + hpp[i])
        nG[i] = 0.5*(G[i] + Gpp[i])
    
    

g = 9.81
wdir = "./Diagnostics/Flux1/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 

for i in range(1,15):
    FluxGhostCells = 1
    startx = -100
    endx = 100
    dx = 10.0 / (2**i)
    l = 0.5 / (2*sqrt(g*2))
    dt = l*dx
    t = 0 
    endt = 0.01  
    x = arange(startx,endx+0.1*dx,dx) 
    xbeg = array([x[0] - 2*dx,x[0] - dx])      
    xend = array([x[-1] + dx,x[-1] + 2*dx]) 
    
    xu = arange(startx - (FluxGhostCells+0.5)*dx ,endx + (FluxGhostCells+0.5)*dx  +0.1*dx,0.5*dx) 
    
    xbc =  concatenate((xbeg,x,xend))  
    
    n = len(x)
    
    a0 = 1
    a1 = 1
    
    theta = 2
    
    h,G = hGSol(x,t,a0,a1,g)
    u = uSol(x,t,a0,a1,g)
    nh = zeros(n)
    nG = zeros(n)
    
    while t < endt:
        print(t)
        hbeg,Gbeg = hGSol(xbeg,t,a0,a1,g)
        hend,Gend = hGSol(xend,t,a0,a1,g)
        
        hbeg1,Gbeg1 = hGSol(xbeg,t+dt,a0,a1,g)
        hend1,Gend1 = hGSol(xend,t+dt,a0,a1,g)
        
        un = uSol(xu,t,a0,a1,g)
        unp1 = uSol(xu,t + dt,a0,a1,g)
        
        RKStep(h,G,un,unp1,hbeg,hend,Gbeg,Gend,hbeg1,hend1,Gbeg1,Gend1,g,theta,FluxGhostCells,dx,dt,nh,nG)
        
        t = t + dt
    
    hF,GF = hGSol(x,t,a0,a1,g)
        
    hnorm = norm(nh - hF)/ norm(hF)
    Gnorm = norm(nG - GF)/ norm(GF)
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)
        
    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
        file1.write(s)