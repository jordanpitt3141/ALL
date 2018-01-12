# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time
from scipy.interpolate import interp1d
from scipy import signal
from scipy import sqrt
from numpy.fft import fft

def DrybedANA(h1,x,t,g):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    
    
    for i in range(n):
         if(x[i] >= -t*sqrt(g*h1) and x[i] <= 2*t*sqrt(g*h1) ):
             u[i] = 2.0 / 3.0 *(sqrt(g*h1) + x[i] / t)
             h[i] = 4.0 / (9.0 * g) *(sqrt(g*h1) - 0.5*x[i] / t)**2
             ux = 2.0 / 3.0 *(1.0 / t)
             uxx = 0
             hx = 2.0 / (9.0 * g * t*t) *(x[i] - 2*t*sqrt(g*h1))
             G[i] = u[i]*h[i] - h[i]*h[i]*hx*ux
         elif(x[i] < -t*sqrt(g*h1)):
             h[i] = h1
             
    return h,u, G
    

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

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    

#FD solution 

#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 

def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def Dambreak(h0,h1,x0,x):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        if(x[i] > x0):
            h[i] = h0
        else:
            h[i] = h1
    
    return h,u,G,b
    
def DambreakS(h0,h1,x0,x,diffuse):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        h[i] = h0 + 0.5*(h1 - h0)*(1 + tanh(diffuse*(x0 - x[i])))
    
    return h,u,G,b
    
def DamNreakDRYANA(h1,x,t,g):
    n = len(x)
    bed = zeros(n)
    h, u,G = DrybedANA(h1,x,t,g)
    G1 = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    
    return h,u,G,G1


def solitoninitGana(a0,a1,g,x,t0,bot,dx):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    ux = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    i3 = 1.0/ 3.0
    for i in range(n):
        phi = x[i] - c*t0;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        tanhkphi = sechkphi*((exp(k*phi) - exp(-k*phi))/2.0)
        hdx = -2*a1*k*tanhkphi*sechkphi*sechkphi
        hdxx = a1*(4*k*k*tanhkphi*tanhkphi*sechkphi*sechkphi - 2*k*k*sechkphi*sechkphi*sechkphi*sechkphi)
        bx[i] = bot
        h[i] = a0 + a1*sechkphi*sechkphi
        u[i] =  c* ((h[i] - a0) / h[i])
        ux[i] = (a0*c*hdx/(h[i]*h[i]))
        G[i] = u[i]*h[i] - i3*h[i]*h[i]*h[i]*(a0*c*(h[i]*hdxx - 2*hdx*hdx)/(h[i]*h[i]*h[i])) - h[i]*h[i]*hdx*(a0*c*hdx/(h[i]*h[i]))
         
    
    return h,u,G,bx,ux 
    
def SolitonMass(a,b,a0,a1,g):
    c = sqrt(g*(a0 + a1)) 
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    
    return a0*(b - a) + a1/k*(tanh(k*b) - tanh(k*a))
    
def SolitonMomentum(a,b,a0,a1,g):
    c = sqrt(g*(a0 + a1)) 
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    
    return c*a1/k*(tanh(k*b) - tanh(k*a))
  

def sech(x) :
    return (2./(exp(x) + exp(-x)))
    
    
def SolE(xbeg,xend):
    AB = 21.0068*tanh(0.555719*xbeg) - 19.2569*arctanh(0.641689*tanh(0.555719*xbeg))
    AE = 21.0068*tanh(0.555719*xend) - 19.2569*arctanh(0.641689*tanh(0.555719*xend))
    
    BB = 9.81*(xbeg) + tanh(0.555719*xbeg)*(2.88329*sech(0.555719*xbeg)**2 + 30.4805)
    BE =9.81*(xend) + tanh(0.555719*xend)*(2.88329*sech(0.555719*xend)**2 + 30.4805)
    
    CB = 307.641*(tanh(0.555719*xbeg)*(0.049539 - 0.00937224*sech(0.555719*xbeg)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xbeg))))
    CE = 307.641*(tanh(0.555719*xend)*(0.049539 - 0.00937224*sech(0.555719*xend)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xend))))

    
    A = AE - AB  
    B = BE - BB    
    C = CE - CB
    
    print(A,B,C)


    #1527.68293
    return 0.5*(A + B + C)


#Solver
#So our solver solves the analytic soliton problem with second order accuracy in h, u and G. 
"""
wdir = "../../../../../data/raw/DryTest/Solver/SolitonENERGY/theta1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "norms.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx','theta','hnorm', 'Gnorm', 'unorm', 'Rhnorm', 'RGnorm', 'dunorm'])    


for j in range(9):        

    a0 = 1
    a1 = 0.7
    g = 9.81
    t0 = 0
    bot = 0
    
    dx = 1.0 / 2**j
    l =  0.5 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -50
    endx = 250 + 0.9*dx
    startt = 0.0
    endt = 1 + dt
            
    szoomx = startx
    ezoomx = endx
    
    t0 = 0
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    ts = []
    
    
    n = len(x)  
    
    theta = 2
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    
    idx = 1.0 / dx
    
    
    h,u,G,bx,ux = solitoninitGana(a0,a1,g,x,t0,bot,dx)
    
    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)
    uMend = u[-1]*ones(unBC) 
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    
    hp_c = mallocPy(n)
    Gp_c = mallocPy(n)
    hpp_c = mallocPy(n)
    Gpp_c = mallocPy(n)  
    
    ct = 0
    while ct < endt:
        
        evolvewrapperconsistenttime(G_c, h_c,hMbeg_c , hMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,g,dx, dt,n,GhnBC,unBC,nGhhbc,nubc,theta, hhbc_c,Ghbc_c,ubc_c,Gp_c,hp_c, Gpp_c,hpp_c)
        ct = ct + dt
        print(ct)

    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n)     
    getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
    
    Mass = hALLW(hhbc_c,n,dx) 
    Momentum = uhALLW(hhbc_c,ubc_c,n,dx)
    Hamil = HamilW(hhbc_c,ubc_c,n,dx)
    
    MassA = SolitonMass(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    MomentumA = SolitonMomentum(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    EA = SolE(x[0] - 0.5*dx,x[-1] + 0.5*dx)
    
    MassRel = abs(Mass - MassA) / abs(MassA)
    MomeRel = abs(Momentum - MomentumA) / abs(MomentumA)
    EnerRel = abs(Hamil - EA)/ abs(EA)
    
    ubcC = copyarrayfromC(ubc_c,nubc)
    uCti = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)

    hA,uA,GA,bx_ta,ux_ta = solitoninitGana(a0,a1,g,x,t0 + ct,bot,dx)
    
    unorm = norm(uA - uCti,ord =1) / norm(uA,ord=1)
    hnorm = norm(hA - hC,ord =1) / norm(hA,ord=1)
    Gnorm = norm(GA - GC,ord =1) / norm(GA,ord=1)
    
    s = wdir + "Mass.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MassRel)
        file1.write(s)

    s = wdir + "Mome.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MomeRel)
        file1.write(s)

    s = wdir + "Ener.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",EnerRel)
        file1.write(s)

    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
        file1.write(s)         

"""




## This FEM reconstructs the soliton problem (with analytic G) with second order accuracy for h (G) at centres and edges and u and du at the edges
"""
wdir = "../../../../../data/raw/DryTest/FEMENERGY/Soliton/theta1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "norms.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx','theta','hnorm', 'Gnorm', 'unorm', 'Rhnorm', 'RGnorm', 'dunorm'])    


for j in range(12):        

    a0 = 1
    a1 = 0.7
    g = 9.81
    t0 = 0
    bot = 0
    
    dx = 1.0 / 2**j
    l =  0.5 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -50
    endx = 250 + 0.9*dx
    startt = 0.0
    endt = 50
            
    szoomx = startx
    ezoomx = endx
    
    t0 = 0
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    ts = []
    
    
    n = len(x)  
    
    theta = 1
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    
    idx = 1.0 / dx
    
    
    h,u,G,bx,ux = solitoninitGana(a0,a1,g,x,t0,bot,dx)
    
    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)
    uMend = u[-1]*ones(unBC) 
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    
           
    getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)

    Mass = hALLW(hhbc_c,n,dx) 
    Momentum = uhALLW(hhbc_c,ubc_c,n,dx)
    Hamil = HamilW(hhbc_c,ubc_c,n,dx)
    
    MassA = SolitonMass(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    MomentumA = SolitonMomentum(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    EA = SolE(x[0] - 0.5*dx,x[-1] + 0.5*dx)


    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    ubcC = copyarrayfromC(ubc_c,nubc)
    uCti = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    

    
    #Calculate u gradients
    du = []
    xdu = []
    for i in range(n):
        
        uai =2*idx*idx*(ubcC[2*i + unBC - 1] - 2*ubcC[2*i + unBC] + ubcC[2*i + unBC + 1])
        ubi =idx*(-ubcC[2*i + unBC - 1]+ ubcC[2*i + unBC + 1])
        duiph = uai*(dx) + ubi;
        duimh = -uai*(dx) + ubi;
        du.append(duimh)
        du.append(duiph)
        xdu.append(x[i] - 0.5*dx)
        xdu.append(x[i] + 0.5*dx) 
        
    hh,hu,hG,hbx,hux = solitoninitGana(a0,a1,g,xdu,t0,bot,dx)
    
    xhbc = []
    xubc = []
    for i in range(len(xG)):
        if(i ==0):           
            xubc.append(xG[i] - 0.5*dx)
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)           
        else:
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)
            
        xhbc.append(xG[i] - 0.5*dx)
        xhbc.append(xG[i])
        xhbc.append(xG[i] + 0.5*dx)
            
    xubc = array(xubc)    
    xhbc = array(xhbc)
    
    Rh,Ru,RG,Rbx,Rux = solitoninitGana(a0,a1,g,xhbc,t0,bot,dx)
    
    unorm = norm(u - uCti,ord =1) / norm(u,ord=1)
    hnorm = norm(h - hC,ord =1) / norm(h,ord=1)
    Gnorm = norm(G - GC,ord =1) / norm(G,ord=1)
    
    MassRel = abs(Mass - MassA) / abs(MassA)
    MomeRel = abs(Momentum - MomentumA) / abs(MomentumA)
    EnerRel = abs(Hamil - EA)/ abs(EA)
    # derivatives and reconstructions
    
    rhnorm = norm(Rh - hhbcC,ord =1) / norm(Rh,ord=1)
    rGnorm = norm(RG - GhbcC,ord =1) / norm(RG,ord=1)
    
    dunorm = norm(hux - du,ord =1) / norm(hux,ord=1)
    
    s = wdir + "norms.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx),str(theta),str(hnorm), str(Gnorm), str(unorm),str(rhnorm),str(rGnorm), str(dunorm), str(MassRel), str(MomeRel), str(EnerRel)])  


    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
        file1.write(s)         

    s = wdir + "rh.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rhnorm)
        file1.write(s)

    s = wdir + "rG.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rGnorm)
        file1.write(s) 
        
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "du.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",dunorm)
        file1.write(s) 
        
    s = wdir + "Mass.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MassRel)
        file1.write(s) 
        
    s = wdir + "Momentum.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MomeRel)
        file1.write(s) 
        
    s = wdir + "Energy.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",EnerRel)
        file1.write(s) 
 
#add some deallocs
 
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(x_c)
    deallocPy(u_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
"""    
  

#Dry bed test
st = time()
h0 = 0.0
h1 = 0.3
x0 = 0
g = 9.81

t0 = 0

dx = 0.01
l =  0.5 / sqrt(g*(h1 + h0))
dt = l*dx
startx = -50
endx = 50 + 0.9*dx
startt = t0
endt = 5 + t0 
szoomx = startx
ezoomx = endx

diffuse = 0.5
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
ts = []


n = len(x)  

theta = 1

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC

idx = 1.0 / dx

#FEM handles dry dam-break with 0 height and 0 velocity well          
#h,u,G,b = Dambreak(h0,h1,x0,x)


h,u,G,b = Dambreak(h0,h1,x0,x)

#h,u,G,b = DambreakS(h0,h1,0,x,diffuse)
#h,u,G,G1 =DamNreakDRYANA(h1,x,t0,g)


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = 0*ones(unBC) 
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)


hp_c = mallocPy(n)
Gp_c = mallocPy(n)
hpp_c = mallocPy(n)
Gpp_c = mallocPy(n)


ct = startt
while ct < endt:
    
    #evolvewrapperconsistenttime(G_c, h_c,hMbeg_c , hMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,g,dx, dt,n,GhnBC,unBC,nGhhbc,nubc,theta, hhbc_c,Ghbc_c,ubc_c,Gp_c,hp_c, Gpp_c,hpp_c)
    dt = evolvewrapperADAP(G_c, h_c,hMbeg_c , hMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,g,dx, dt,n,GhnBC,unBC,nGhhbc,nubc,theta, hhbc_c,Ghbc_c,ubc_c,Gp_c,hp_c, Gpp_c,hpp_c)
       
    
    ct = ct + dt
    
    if(dt < 10**-8):
        break
    print(ct)


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 

hF,uF,GF,G1F =DamNreakDRYANA(h1,x,ct,g)  
#getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
#ubcC = copyarrayfromC(ubc_c,nubc)
#ufC = ubcC[unBC:-unBC:2]


deallocPy(h_c)
deallocPy(G_c)
deallocPy(hp_c)
deallocPy(Gp_c)
deallocPy(hpp_c)
deallocPy(Gpp_c)
deallocPy(u_c)

deallocPy(ubc_c)
deallocPy(hhbc_c)
deallocPy(Ghbc_c)

deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
et = time()
tt = et - st
