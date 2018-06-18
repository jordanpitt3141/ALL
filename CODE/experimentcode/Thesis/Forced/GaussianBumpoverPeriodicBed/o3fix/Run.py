from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from Serre3 import *
from numpy.linalg import norm
from time import time
from scipy.special import erf

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def TDMApy(a,b,c,d):
    n = len(d)
    alpha = []
    beta = []
    x = [0]*n
    
    alpha.append((1.0*c[0])/b[0])
    beta.append((1.0*d[0])/b[0] )  
 
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha.append(c[i]* m)
        beta.append((d[i] - a[i-1]*beta[i-1]) * m)
        
    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta.append((d[n-1] - a[n-2]*beta[n-2]) * m)  

    x[n-1] = beta[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
 
    return array(x)
    
#Tested from CFD forum post
def pentadiagsolve(e,a,d,c,f,B):
    n = len(d)
    X = zeros(n)
    
    for i in range(1,n-1):
        xmult = float(a[i-1]) / d[i-1]
        
        d[i] = d[i] - xmult*c[i-1]
        c[i] = c[i] - xmult*f[i-1]
        B[i] = B[i] - xmult*B[i-1]
        
        xmult = float(e[i-1]) /d[i-1]
        a[i] = a[i] - xmult*c[i-1]
        d[i+1] = d[i+1] - xmult*f[i-1]
        B[i+1] = B[i+1] - xmult*B[i-1]
        
    xmult = float(a[n-2]) / d[n-2]
    d[n-1] = d[n-1] - xmult*c[n-2]
    X[n-1] = (B[n-1] - xmult*B[n-2]) / float(d[n-1])
    X[n-2] = (B[n-2] - c[n-2]*X[n-1]) / float(d[n-2])
    
    for i in range(n-3,-1,-1):
        X[i] = (B[i] - f[i]*X[i+2] - c[i]*X[i+1])/float(d[i])
        
    return X    
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
    
def midpointtocellaverages(mq,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    idx = 1.0/dx
    i24 = 1.0 / 24.0
    n = len(mq)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    for i in range(1,n-1):
        ai = -i24
        bi = 26*i24
        ci = -i24

        a[i-1] = ai
        b[i] = bi
        c[i] = ci
    
    #i = 0
    i = 0
    ai =0.0 #-i24
    bi =1.0 #26*i24
    ci =0.0 #-i24

    b[i] = bi
    c[i] = ci
    
    #mq[i] = mq[i] - ai*qbeg[0]
    
    #i = 0
    i = n-1
    ai =0.0# -i24
    bi =1.0# 26*i24
    ci =0.0# -i24

    a[i-1] = ai
    b[i] = bi
    
    #mq[i] = mq[i] - ci*qend[0]
    
    q = TDMApy(a,b,c,mq)
    
    return q
    
def cellaveragestomidpoints(q,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    i24 = 1.0 / 24.0
    n = len(q)
    mq = zeros(n)
    for i in range(1,n-1):
        #iterate  over the cell midpoints, there are 2 edge values for each (except the first and last cell)
        
        #variables
        #ai = (q[i+1] - 2*q[i] + q[i-1])*0.5*idx*idx
        #bi = (q[i+1] - q[i-1])*0.5*idx
        ci = i24*(-q[i+1] + 26*q[i]  -q[i-1])
        mq[i] = ci
    
    #i = 0
    i = 0
    ci = q[i] #i24*(-q[i+1] + 26*q[i] - qbeg[0])
    mq[i] = ci
    
    #i = n-1
    i = n-1
    ci = q[i]#i24*(-qend[0] + 26*q[i] - q[i-1])
    mq[i] = ci 
    
    return mq
    
def hI(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    return a0*x - a1*sqrt(a4)*sqrt(pi/2.)*erf((a3 + a2*t - x)/(sqrt(2)*sqrt(a4)))
    
def hA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7):
    hxi1 = hI(xi1,t,a0,a1,a2,a3,a4,a5,a6,a7)
    hxi2 = hI(xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    hAv = (hxi2 - hxi1)/(xi2 - xi1)
    return hAv

def uI(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    return - a5*sqrt(a4)*sqrt(pi/2.)*erf((a3 + a2*t - x)/(sqrt(2)*sqrt(a4)))

def uA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7):
    uxi1 = uI(xi1,t,a0,a1,a2,a3,a4,a5,a6,a7)
    uxi2 = uI(xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    uAv = (uxi2 - uxi1)/(xi2 - xi1)
    return uAv

    
def GI(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    return (a5*(2*pow(e,((a3 + a2*t)*x)/a4)*pow(a1*pow(e,((a3 + a2*t)*x)/a4) + a0*pow(e,(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2))/(2.*a4)),3)*(-a3 - a2*t + x) + \
       3*a1*pow(a4,1.5)*pow(e,(2*(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2)))/a4)*sqrt(pi)*erf((-a3 - a2*t + x)/sqrt(a4)) + \
       3*a0*pow(a4,1.5)*pow(e,(2*(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2)))/a4)*sqrt(2*pi)*erf((-a3 - a2*t + x)/(sqrt(2)*sqrt(a4)))))/(6.*a4*pow(e,(2*(pow(a3 + a2*t,2) + pow(x,2)))/a4))

def GA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7):
    Gxi1 = GI(xi1,t,a0,a1,a2,a3,a4,a5,a6,a7)
    Gxi2 = GI(xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    GAv = (Gxi2 - Gxi1)/ (xi2 - xi1)
    return GAv
    
def ForcedbedA(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    ha = zeros(n)
    ua = zeros(n)
    Ga = zeros(n)
    
    for i in range(n):
        xi1 = x[i] - 0.5*dx
        xi2 = x[i] + 0.5*dx
        ha[i] = hA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
        ua[i] = uA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
        Ga[i] = GA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)

    return ha,ua,Ga
    
    
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
           
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        
        G[i] = u[i]*h[i] - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
        
    return h,u,G
 


def solveGfromuh(u,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of u,h and gives midpoint values of G    
    
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(u)
    
    G = zeros(n)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
        
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
    
    G[i] = ai*ubeg[-2] + bi*ubeg[-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]

    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*ubeg[-1] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*uend[0]
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*uend[0] + ei*uend[1]

    return G
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u   
    
  
    #Forcing Problem    
wdir = "../../../../../../../data/raw/Forced/FDVM3NoBed/GaussBedAll/EvotBest3/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(4,15):
    g =9.81

    a0 = 1
    a1 = 0.2
    a2 = 1.3
    a3 = 0.4
    a4 = 1.5
    a5 = 0.1
    a6 = 0
    a7 = 0.1
    
    width = 20.0
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (a5 + sqrt(g*(a0 + a1)))
    dt = l*dx
    startx = -width/2
    endx = width/2
    startt = 0.0
    endt = 1
            
    
    t = startt
    
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
    
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    
    xmbeg = arange(startx - niBC*dx ,x[0], dx)
    xmend = arange(x[-1] + dx ,x[-1] + (niBC+ 0.1)*dx, dx)
    

    ts = []
    
    n = len(x)  
    theta = 2
    
    gap = int(1.0/dt)
                    
    hm,um,Gm = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hmbeg,umbeg,Gmbeg = ForcedbedM(xmbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hmend,umend,Gmend = ForcedbedM(xmend,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    cnBC = niBC - nGsBC
    
    umbc = concatenate([umbeg[-cnBC:],um,umend[:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[:cnBC]])  
    Gmbc = concatenate([Gmbeg[-cnBC:],Gm,Gmend[:cnBC]]) 
    
    #calculate averages
    Gabc = midpointtocellaverages(Gmbc,dx)
    habc = midpointtocellaverages(hmbc,dx)
    uabc = midpointtocellaverages(umbc,dx)
    
    
    Gabeg = Gabc[:cnBC]
    Ga = Gabc[cnBC:-cnBC]
    Gaend = Gabc[-cnBC:]
    habeg = habc[:cnBC]
    ha = habc[cnBC:-cnBC]
    haend = habc[-cnBC:]
    uabeg = uabc[:cnBC]
    ua = uabc[cnBC:-cnBC]
    uaend = uabc[-cnBC:]
    
    

    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
    
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
    
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
    
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
    
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
    
    u_c = mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
    x_c = copyarraytoC(x)

    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7)
        t = t + dt
        ts.append(t)
        print(t)
        
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
    
    GaC = copyarrayfromC(Ga_c,n)
    haC = copyarrayfromC(ha_c,n)
    
    GC = copyarrayfromC(G_c,n)
    hC = copyarrayfromC(h_c,n)


    ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)


    uC = copyarrayfromC(u_c,n)
    
    haA,uaA,GaA = ForcedbedA(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hmA,umA,GmA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    hnorm = norm(haC - haA, ord=2)/ norm(haA, ord=2)
    unorm = norm(uC - umA, ord=2)/ norm(umA, ord=2)
    Gnorm = norm(GaC - GaA, ord=2)/ norm(GaA, ord=2)
    
    
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hnorm)
        file1.write(s)
        

    
    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Gnorm)
        file1.write(s)   
    
    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",unorm)
        file1.write(s) 
 
    deallocPy(ha_c)
    deallocPy(Ga_c)

    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c) 
    
    deallocPy(hmbeg_c)
    deallocPy(umbeg_c)
    deallocPy(hmend_c)
    deallocPy(umend_c)
 
    deallocPy(habeg_c)
    deallocPy(Gabeg_c)
    deallocPy(uabeg_c)
    deallocPy(haend_c)
    deallocPy(Gaend_c)
    deallocPy(uaend_c)


    
"""
##Accuracy Test
### Soliton Accuracy ################
wdir = "../../../../../../../data/raw/Forced/FDVM3NoBed/GaussBedAll/Soliton1/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity','Hamiltonian Difference'])
    
for k in range(6,21):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    #g = 1.0
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -250.0
    endx = 250.0 + dx
    startt = 0.0
    endt = 1 + dt
    
    print(dx,dt)
        
    szoomx = startx
    ezoomx = endx
        
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
        
        
    wdatadir = wdir+ str(k) + "/" 
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    gap =int(10.0/dt)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    t0 = 0.0   
    hm,um = solitoninit(n,a0,a1,g,x,t0,dx)
                
    umbeg = um[0]*ones(niBC)
    umend = um[-1]*ones(niBC)
    hmbeg = hm[0]*ones(niBC)
    hmend = hm[-1]*ones(niBC)    
                
    #calculate G midpoint
    cnBC = niBC - nGsBC
        
    umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
    Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
    #calculate averages
    Gabc = midpointtocellaverages(Gmbc,dx)
    habc = midpointtocellaverages(hmbc,dx)
    uabc = midpointtocellaverages(umbc,dx)
        
    #so we can just go from here with Ga ang ha?
    Gabeg = Gabc[0:cnBC]
    Ga = Gabc[cnBC:-cnBC]
    Gaend = Gabc[-cnBC:]
    habeg = habc[0:cnBC]
    ha = habc[cnBC:-cnBC]
    haend = habc[-cnBC:]
    uabeg = uabc[0:cnBC]
    ua = uabc[cnBC:-cnBC]
    uaend = uabc[-cnBC:]
    
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
    
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
    
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
    
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
    
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
    
    u_c = mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])  
    
    x_c = copyarraytoC(x)
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*niBC)
    ubc_c = mallocPy(n + 2*niBC)
    Evals = []
    
    for i in range(1,len(t)):
        
        if(i % gap == 0 or i ==1):
            
            ca2midpt(ha_c,dx,n,h_c)
            ca2midpt(Ga_c,dx,n,G_c)
            ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
            
            conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
            conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)
            Eval =1
            #Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
            
            Evals.append(Eval)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            
            c = sqrt(g*(a0 + a1))
            htrue = zeros(n)
            utrue = zeros(n)
            for j in range(n):             
                he = soliton(x[j],t[i],g,a0,a1)
                htrue[j] = he
                utrue[j] = c* ((he - a0) / he) 
                
            s = wdatadir + "saveoutputts" + str(i) + ".txt"
            print t[i]
            print(h[3],G[3]) 
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
                writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
                 
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC,x_c,t[i-1],0.0,0.0,1,1,1,0.0,1,1)
    
        #evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
        print t[i]
        print(h[3],G[3]) 
    
        
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
    ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
    
    conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
    conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
    Eval = 1
    
    Evals.append(Eval)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[i],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
        
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
     
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)
    normHamdiff = (Evals[-1] - Evals[0])/ Evals[0]

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi),str(normHamdiff)])
        
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normhdiffi)
        file1.write(s)
    
    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normudiffi)
        file1.write(s) 
"""

