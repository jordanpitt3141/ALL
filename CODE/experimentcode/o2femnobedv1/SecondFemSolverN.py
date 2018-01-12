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
from numpy.linalg import norm
from time import time

"""
#Recon good, so B.C?    
def ReconTest(h,G,h0,h1,h2,G0,G1,G2,uev0,uev1,theta,dx):
    n = len(h)
    idx = 1.0 / dx 
    
    Gpjmhs = zeros(n)
    Gmjphs = zeros(n)
    Gpjphs = zeros(n)
    Gmjp1hs = zeros(n)
    hpjmhs = zeros(n)
    hmjphs = zeros(n)
    hpjphs = zeros(n)
    hmjp1hs = zeros(n)
        
    
    for j in range(1,n-2):
        #First we should reconstruct all the ones we need for this step
        Gpjmh = Recon2P(G,0,0,0,j-1,theta,dx)
        Gmjph = Recon2M(G,0,0,0,j,theta,dx)
        Gpjph = Recon2P(G,0,0,0,j,theta,dx)
        Gmjp1h = Recon2M(G,0,0,0,j+1,theta,dx)
        
        hpjmh = Recon2P(h,0,0,0,j-1,theta,dx)
        hmjph = Recon2M(h,0,0,0,j,theta,dx)
        hpjph = Recon2P(h,0,0,0,j,theta,dx)
        hmjp1h = Recon2M(h,0,0,0,j+1,theta,dx)
        
        Gpjmhs[j] = Gpjmh
        Gmjphs[j] = Gmjph
        Gpjphs[j] = Gpjph
        Gmjp1hs[j] = Gmjp1h
        hpjmhs[j] = hpjmh
        hmjphs[j] = hmjph
        hpjphs[j] = hpjph
        hmjp1hs[j] = hmjp1h

    #This gives coefficients from one equation, but we have many equations
    
    #j = 0
    j = 0
    
    Gpjmh = Recon2P(G,G0,G[0],G[1],j-1,theta,dx)
    Gmjph = Recon2M(G,G0,0,0,j,theta,dx)    
    Gpjph = Recon2P(G,0,0,0,j,theta,dx)
    Gmjp1h = Recon2M(G,0,0,0,j+1,theta,dx)
        
    hpjmh = Recon2P(h,h0,h[0],h[1],j-1,theta,dx)
    hmjph = Recon2M(h,h0,0,0,j,theta,dx)
    hpjph = Recon2P(h,0,0,0,j,theta,dx)
    hmjp1h = Recon2M(h,0,0,0,j+1,theta,dx)
    
    Gpjmhs[j] = Gpjmh
    Gmjphs[j] = Gmjph
    Gpjphs[j] = Gpjph
    Gmjp1hs[j] = Gmjp1h
    hpjmhs[j] = hpjmh
    hmjphs[j] = hmjph
    hpjphs[j] = hpjph
    hmjp1hs[j] = hmjp1h
    

    
    #j = n-2
    j = n-2
    
    Gpjmh = Recon2P(G,0,0,0,j-1,theta,dx)
    Gmjph = Recon2M(G,0,0,0,j,theta,dx)
    Gpjph = Recon2P(G,0,0,G1,j,theta,dx)
    Gmjp1h = Recon2M(G,0,0,G1,j+1,theta,dx)
    
    hpjmh = Recon2P(h,0,0,0,j-1,theta,dx)
    hmjph = Recon2M(h,0,0,0,j,theta,dx)
    hpjph = Recon2P(h,0,0,h1,j,theta,dx)
    hmjp1h = Recon2M(h,0,0,h1,j+1,theta,dx)
    
    Gpjmhs[j] = Gpjmh
    Gmjphs[j] = Gmjph
    Gpjphs[j] = Gpjph
    Gmjp1hs[j] = Gmjp1h
    hpjmhs[j] = hpjmh
    hmjphs[j] = hmjph
    hpjphs[j] = hpjph
    hmjp1hs[j] = hmjp1h
    
    
    #j = n-1
    j = n-1
    
    Gpjmh = Recon2P(G,0,0,G1,j-1,theta,dx)
    Gmjph = Recon2M(G,0,0,G1,j,theta,dx)
    Gpjph = Recon2P(G,0,G1,G2,j,theta,dx)
    Gmjp1h = Recon2M(G,G[-1],G1,G2,j+1,theta,dx)
    
    hpjmh = Recon2P(h,0,0,h1,j-1,theta,dx)
    hmjph = Recon2M(h,0,0,h1,j,theta,dx)
    hpjph = Recon2P(h,0,h1,h2,j,theta,dx)
    hmjp1h = Recon2M(h,h[-1],h1,h2,j+1,theta,dx)
    
    Gpjmhs[j] = Gpjmh
    Gmjphs[j] = Gmjph
    Gpjphs[j] = Gpjph
    Gmjp1hs[j] = Gmjp1h
    hpjmhs[j] = hpjmh
    hmjphs[j] = hmjph
    hpjphs[j] = hpjph
    hmjp1hs[j] = hmjp1h
    
    #NU = TDMApy(uais,ubis,ucis,nGis)
    
    return Gpjmhs,Gmjphs,Gpjphs,Gmjp1hs,hpjmhs,hmjphs,hpjphs,hmjp1hs
"""

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

def minmodpy(a, b, c):
    if((a > 0) and (b>0) and (c>0)):
        return min(a,b,c)
    elif((a < 0) and (b<0) and (c<0)):
        return max(a,b,c)
    else:
        return 0.0

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
    

def interpquarticvalPy(aj,bj,cj,dj,ej,xj,x):
    
    return aj*(x -xj)*(x -xj)*(x -xj)*(x -xj) + bj*(x -xj)*(x -xj)*(x -xj) \
    + cj*(x -xj)*(x -xj) + dj*(x -xj)+ ej
    
def interpquarticgradPy(aj,bj,cj,dj,ej,xj,x):
    
    return 4*aj*(x -xj)*(x -xj)*(x -xj) + 3*bj*(x -xj)*(x -xj) \
    + 2*cj*(x -xj) + dj
    
def interpquartcoeffPy(q,j,dx):
    i24 = 1.0 / 24.0
    i12 = 1.0 / 12.0
    idx = 1.0/dx
    aj = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2])
    bj = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2])
    cj = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2])
    dj = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2])
    ej = q[j]
    
    return aj,bj,cj,dj,ej


#FD solution 

def getGfromupy(h,u,u0,u1,h0,h1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        
        D = th
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
            
    D = th
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
        
    D = th
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 
   
def testsolSin(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = sin(3*xp)
        h[i] = sin(10*xp) + 3
        G[i] = sin(3*xp)*(sin(10*xp) +3) - 30*(sin(10*xp) + 3)**2 *cos(10*xp)*cos(3*xp) + 3*(sin(10*xp) +3)**3*sin(3*xp)
        
    return h,u,G

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a    

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a
    
def testsolSol(a0,a1,g,x):
    n = len(x)
    u = zeros(n)
    ux = zeros(n)
    h = zeros(n)
    G = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        xp = x[i]
        
        
        h[i] = a0 + a1*sech2(k*xp)
        u[i] = c*(1 - float(a0)/h[i])
        
        hx = -2*a1*k*tanh(k*xp)*sech2(k*xp)
        hxx = 2*a1*k*k*(cosh(2*k*xp) -2)*sech2(k*xp)*sech2(k*xp)
        
        ux[i] = a0*c*hx / (h[i]*h[i])
        
        G[i] = u[i]*h[i] - a0*c*hx*hx - 1.0/3*(a0*c*(h[i]*hxx - 2*hx*hx))
        
    return h,u,G,ux
    
def testsolSec(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = 10*xp
        h[i] = 3*xp + 5
        G[i] = u[i]*h[i] - h[i]**2*(3)*(10)
        
    return h,u,G
    
def Recon2M(q,qim1M,qiM,qip1M,i,theta,dx):
    n = len(q)
    
    if(i == 0):
        qi = q[i]
        qip1 = q[i+1]
        qim1 = qim1M
    elif(i == n-1):
         qi = q[i]
         qip1 = qip1M
         qim1 = q[i-1]
    elif(i < 0 or i > n-1):
        qi = qiM
        qip1 = qip1M
        qim1 = qim1M
    else:
        qi = q[i]
        qip1 = q[i+1]
        qim1 = q[i-1]
    

    
    #_j+1/2 ^-
    dqib = qi - qim1
    dqim = 0.5*(qip1 - qim1)
    dqif = qip1 - qi
    
    dqi = minmodpy(theta*dqib,theta*dqif,dqim)
    
    qir = qi + 0.5*dqi
    
    return qir
    
def Recon2P(q,qiM,qip1M,qip2M,i,theta,dx):
    n = len(q)
    if(i == n-2):
        qi = q[i]
        qip1 = q[i+1]
        qip2 = qip2M
    elif(i == n-1):
         qi = q[i]
         qip1 = qip1M
         qip2 = qip2M
    elif(i < 0 or i > n-1):
        qi = qiM
        qip1 = qip1M
        qip2 = qip2M
    else:
        qi = q[i]
        qip1 = q[i+1]
        qip2 = q[i+2]
    
    
    #_j+1/2 ^-
    dqip1b = qip1 - qi
    dqip1m = 0.5*(qip2 - qi)
    dqip1f = qip2 - qip1
    
    dqip1 = minmodpy(theta*dqip1b,theta*dqip1f,dqip1m)
    
    qip1l = qip1- 0.5*dqip1
    
    return qip1l
    
    
def FEM2Assem(h,G,hbeg,hend,Gbeg,Gend,umh,uph,theta,dx):
    n = len(h)
    idx = 1.0 / dx 
    
    
    uais = zeros(n-2)
    ubis = zeros(n-1)
    ucis = zeros(n-2)
    nGis = zeros(n-1)
        
    
    for j in range(1,n-2):
        #First we should reconstruct all the ones we need for this step
        Gpjmh = Recon2P(G,0,0,0,j-1,theta,dx)
        Gmjph = Recon2M(G,0,0,0,j,theta,dx)
        Gpjph = Recon2P(G,0,0,0,j,theta,dx)
        Gmjp1h = Recon2M(G,0,0,0,j+1,theta,dx)
        
        hpjmh = Recon2P(h,0,0,0,j-1,theta,dx)
        hmjph = Recon2M(h,0,0,0,j,theta,dx)
        hpjph = Recon2P(h,0,0,0,j,theta,dx)
        hmjp1h = Recon2M(h,0,0,0,j+1,theta,dx)
        
        #u coefficients
        #lets do uai for j-1/2 ()
        # lets for ubi for j+1/2
        #lets use uci for j+3/2
        
        #split into parts
        uaipar1 = hpjmh + hmjph
        ubipar1 = hpjmh + 3*hmjph + 3*hpjph  + hmjp1h
        ucipar1 = hpjph + hmjp1h
        
        
        uaipar2 = -idx*idx*(hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3)
        ubipar2 = idx*idx*( hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3 +  hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)    
        ucipar2 = -idx*idx*(hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)
        
        uai = uaipar1 + uaipar2
        ubi = ubipar1 + ubipar2
        uci = ucipar1 + ucipar2
        #G coefficients
        #lets do gai for _j-1/2 ^+
        # lets for gbi for _j+1/2 ^-
        # lets for gci for _j+1/2 ^+
        #lets use gdi for _j+3/2 ^-
        
        gai = 2
        gbi = 4
        gci = 4
        gdi = 2
        
        nGi = gai*Gpjmh + gbi*Gmjph + gci*Gpjph + gdi*Gmjp1h
        
        uais[j-1] = uai
        ubis[j] = ubi
        ucis[j] = uci
    
        nGis[j] = nGi
    
    #This gives coefficients from one equation, but we have many equations
    
    #j = 0
    j = 0
    
    Gpjmh = Recon2P(G,Gbeg[-1],G[0],G[1],j-1,theta,dx)
    Gmjph = Recon2M(G,Gbeg[-1],0,0,j,theta,dx)    
    Gpjph = Recon2P(G,0,0,0,j,theta,dx)
    Gmjp1h = Recon2M(G,0,0,0,j+1,theta,dx)
        
    hpjmh = Recon2P(h,hbeg[-1],h[0],h[1],j-1,theta,dx)
    hmjph = Recon2M(h,hbeg[-1],0,0,j,theta,dx)
    hpjph = Recon2P(h,0,0,0,j,theta,dx)
    hmjp1h = Recon2M(h,0,0,0,j+1,theta,dx)
    
    uaipar1 = hpjmh + hmjph
    ubipar1 = hpjmh + 3*hmjph + 3*hpjph  + hmjp1h
    ucipar1 = hpjph + hmjp1h
        
        
    uaipar2 = -idx*idx*(hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3)
    ubipar2 = idx*idx*( hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3 +  hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)    
    ucipar2 = -idx*idx*(hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)
        
    uai = uaipar1 + uaipar2
    ubi = ubipar1 + ubipar2
    uci = ucipar1 + ucipar2

        
    gai = 2
    gbi = 4
    gci = 4
    gdi = 2
        
    nGi = gai*Gpjmh + gbi*Gmjph + gci*Gpjph + gdi*Gmjp1h
    
    ubis[j] = ubi
    ucis[j] = uci
    
    nGis[j] = nGi - uai*umh  
    
    #j = n-2
    j = n-2
    
    Gpjmh = Recon2P(G,0,0,0,j-1,theta,dx)
    Gmjph = Recon2M(G,0,0,0,j,theta,dx)
    Gpjph = Recon2P(G,0,0,Gend[0],j,theta,dx)
    Gmjp1h = Recon2M(G,0,0,Gend[0],j+1,theta,dx)
    
    hpjmh = Recon2P(h,0,0,0,j-1,theta,dx)
    hmjph = Recon2M(h,0,0,0,j,theta,dx)
    hpjph = Recon2P(h,0,0,hend[0],j,theta,dx)
    hmjp1h = Recon2M(h,0,0,hend[0],j+1,theta,dx)
    
    #u coefficients
    #lets do uai for j-1/2 ()
    # lets for ubi for j+1/2
    #lets use uci for j+3/2
    
    #split into parts
    uaipar1 = hpjmh + hmjph
    ubipar1 = hpjmh + 3*hmjph + 3*hpjph  + hmjp1h
    ucipar1 = hpjph + hmjp1h
    
    
    uaipar2 = -idx*idx*(hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3)
    ubipar2 = idx*idx*( hpjmh**3 + hpjmh**2*hmjph +  hpjmh*hmjph**2  + hmjph**3 +  hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)    
    ucipar2 = -idx*idx*(hpjph**3 + hpjph**2*hmjp1h +  hpjph*hmjp1h**2  + hmjp1h**3)
    
    uai = uaipar1 + uaipar2
    ubi = ubipar1 + ubipar2
    uci = ucipar1 + ucipar2
    #G coefficients
    #lets do gai for _j-1/2 ^+
    # lets for gbi for _j+1/2 ^-
    # lets for gci for _j+1/2 ^+
    #lets use gdi for _j+3/2 ^-
    
    gai = 2
    gbi = 4
    gci = 4
    gdi = 2
    
    nGi = gai*Gpjmh + gbi*Gmjph + gci*Gpjph + gdi*Gmjp1h
    
    uais[j-1] = uai
    ubis[j] = ubi

    nGis[j] = nGi - uci*uph
    
    NU = TDMApy(uais,ubis,ucis,nGis)
    NUs = concatenate((array([umh]),NU,array([uph])))
    
    return NUs


   
#Soliton test
dxs = []
unormsC = []
GnormsC = []
hnormsC = []
uxrnormsC = []
uxlnormsC = []
wdir = "../../../../data/raw/FEMtesto2C1SOLThighnonAMP1M/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)
#Sin test, very little error good
for i in range(11,12):   
    a0 = 1.0
    a1 = 1.0
    dx = 10.0/2**i
    l = 0.1
    dt = l*dx
    startx = -70
    endx = 70.0 + dx
    startt = 0.0
    endt = 30 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            
    g = 9.81
    theta = 1.2
    
    gap = int(1.0/dt)
    niBC = 3
            
    xbc,t = makevar(startx - niBC*dx,endx +niBC*dx,dx,startt,endt,dt)
    
    x = xbc[niBC:-niBC]
    
    
    m = len(x)
            
    hbc,ubc,Gbc,uxbc = testsolSol(a0,a1,g,xbc)
    
    GbcPy = getGfromupy(hbc,ubc,0,0,a0,a0,dx)
    
    h = hbc[niBC:-niBC]
    hbeg = hbc[:niBC]
    hend = hbc[-niBC:]
    
    G = Gbc[niBC:-niBC]
    Gbeg = Gbc[:niBC]
    Gend = Gbc[-niBC:]
    
    u = ubc[niBC:-niBC]
    ubeg = ubc[:niBC]
    uend = ubc[-niBC:]
    
    xmh = x[0] - 0.5*dx
    xm1h = x[0] - 1.5*dx
    xph = x[-1] + 0.5*dx
    xp1h = x[-1] + 1.5*dx
    
    umh = ubc[0]
    uph = ubc[-1]
    um1h = ubc[0]
    up1h = ubc[-1]
    
    hph = hbc[-1]
    Gph = Gbc[-1]
    
    hmh = hbc[0]
    Gmh = Gbc[0]
    
    
    Gebeg = [Gmh,Gmh]
    hebeg = [hmh,hmh]
    uebeg = [um1h,umh]
    
    Geend = [Gph,Gph]
    heend = [hph,hph]
    ueend = [uph,up1h]
    
    uFEMpy = FEM2Assem(h,G,hbeg,hend,Gbeg,Gend,umh,uph,theta,dx)
    
    uFEMpycomp = concatenate((array([um1h]),uFEMpy,array([up1h])))
    
    
    
    niBC = 2
    nhbc = 2
    nubc = 3
    nBCs = 2
    nhBC =  (3*(m -2) + 2*2 + 2*nhbc)
    nuBC = (m-1) + 2*nubc
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    u_c = mallocPy(nuBC)
    
    hebeg_c = copyarraytoC(hebeg)
    Gebeg_c = copyarraytoC(Gebeg)
    uebeg_c = copyarraytoC(uebeg)
    
    heend_c = copyarraytoC(heend)
    Geend_c = copyarraytoC(Geend)
    ueend_c = copyarraytoC(ueend)
    
    hhbc_c = mallocPy(nhBC)
    Ghbc_c = mallocPy(nhBC)
    

    
    getufromG(h_c,G_c,hebeg_c,heend_c,Gebeg_c,Geend_c,uebeg_c,ueend_c,theta,dx ,m, nuBC, nhbc,nubc ,u_c, hhbc_c, Ghbc_c);
    
    uFEMC = copyarrayfromC(u_c,nuBC)
    hRecC = copyarrayfromC(hhbc_c,nhBC)
    GRecC = copyarrayfromC(Ghbc_c,nhBC)
    
    xhbc = []
    xhbc.append(xmh)
    xhbc.append(xmh)
    xhbc.append(x[0])
    xhbc.append(x[0] + 0.5*dx)
    for i in range(1,m-1):
        xhbc.append(x[i] - 0.5*dx)
        xhbc.append(x[i])
        xhbc.append(x[i] + 0.5*dx)
    
    xhbc.append(x[-1] - 0.5*dx)
    xhbc.append(x[-1])
    xhbc.append(xph)    
    xhbc.append(xph) 
    
    hhbc,uhbc,Ghbc,uxhbc = testsolSol(a0,a1,g,xhbc)
    
    xubc,t = makevar(xmh - 2*dx,xph + 2.9*dx,dx,startt,endt,dt)
    
    uhhbc,uuhbc,uGhbc,uuxhbc = testsolSol(a0,a1,g,xubc)
    
    dxs.append(dx)
    
    uFEMduer = (-1.5*array(uFEMC[1:-2]) +  2*array(uFEMC[2:-1]) - 0.5*array(uFEMC[3:] )) / dx
    uFEMduel = (1.5*array(uFEMC[2:-1]) -  2*array(uFEMC[1:-2]) + 0.5*array(uFEMC[0:-3] )) / dx
    
    unormC =  norm(uFEMC -uuhbc,ord=1) / norm(uuhbc,ord=1) 
    
    hnormC =  norm(hRecC -hhbc,ord=1) / norm(hhbc,ord=1) 
    
    GnormC =  norm(GRecC -Ghbc,ord=1) / norm(Ghbc,ord=1) 
    
    uxrnormC =  norm(uFEMduer -uuxhbc[1:-2],ord=1) / norm(uuxhbc[1:-2],ord=1) 
    uxlnormC =  norm(uFEMduel -uuxhbc[2:-1],ord=1) / norm(uuxhbc[2:-1],ord=1) 

    uxrnormsC.append(uxrnormC)
    uxlnormsC.append(uxlnormC)
    unormsC.append(unormC)
    hnormsC.append(hnormC)
    GnormsC.append(GnormC)
               


"""    
#Sin test, very little error good
dxs = []
unormsC = []
GnormsC = []
hnormsC = []
wdir = "../../../data/raw/FEMtesto2C1/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)
#Sin test, very little error good
for i in range(5,20):   
    dx = 10.0/2**i
    l = 0.1
    dt = l*dx
    startx = 0
    endx = 10.0 + dx
    startt = 0.0
    endt = 30 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
            
    g = 9.81
    theta = 1.4
    
    gap = int(1.0/dt)
    niBC = 3
            
    xbc,t = makevar(startx - niBC*dx,endx +niBC*dx,dx,startt,endt,dt)
    
    x = xbc[niBC:-niBC]
    
    
    m = len(x)
            
    hbc,ubc,Gbc = testsolSin(xbc)
    
    h = hbc[niBC:-niBC]
    hbeg = hbc[:niBC]
    hend = hbc[-niBC:]
    
    G = Gbc[niBC:-niBC]
    Gbeg = Gbc[:niBC]
    Gend = Gbc[-niBC:]
    
    u = ubc[niBC:-niBC]
    ubeg = ubc[:niBC]
    uend = ubc[-niBC:]
    
    xmh = x[0] - 0.5*dx
    xm1h = x[0] - 1.5*dx
    xph = x[-1] + 0.5*dx
    xp1h = x[-1] + 1.5*dx
    umh = sin(3*xmh)
    uph = sin(3*xph)
    um1h = sin(3*xm1h)
    up1h = sin(3*xp1h)
    
    hph = sin(10*xph) + 3
    Gph = sin(3*xph)*(sin(10*xph) +3) - 30*(sin(10*xph) + 3)**2 *cos(10*xph)*cos(3*xph) + 3*(sin(10*xph) +3)**3*sin(3*xph)
    
    hmh = sin(10*xmh) + 3
    Gmh = sin(3*xmh)*(sin(10*xmh) +3) - 30*(sin(10*xmh) + 3)**2 *cos(10*xmh)*cos(3*xmh) + 3*(sin(10*xmh) +3)**3*sin(3*xmh)
    
    
    Gebeg = [Gmh,Gmh]
    hebeg = [hmh,hmh]
    uebeg = [um1h,umh]
    
    Geend = [Gph,Gph]
    heend = [hph,hph]
    ueend = [uph,up1h]
    
    uFEMpy = FEM2Assem(h,G,hbeg,hend,Gbeg,Gend,umh,uph,theta,dx)
    
    uFEMpycomp = concatenate((array([um1h]),uFEMpy,array([up1h])))
    
    
    
    niBC = 2
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    u_c = mallocPy(m + 3)
    
    hebeg_c = copyarraytoC(hebeg)
    Gebeg_c = copyarraytoC(Gebeg)
    uebeg_c = copyarraytoC(uebeg)
    
    heend_c = copyarraytoC(heend)
    Geend_c = copyarraytoC(Geend)
    ueend_c = copyarraytoC(ueend)
    
    hhbc_c = mallocPy(3*m +2)
    Ghbc_c = mallocPy(3*m+2)
    
    getufromG(h_c,G_c,hebeg_c,heend_c,Gebeg_c,Geend_c,uebeg_c,ueend_c,theta,dx , m, m + 3,3*m + 2,niBC, u_c,hhbc_c,Ghbc_c)
    
    uFEMC = copyarrayfromC(u_c,m + 3)
    hRecC = copyarrayfromC(hhbc_c,3*m +2)
    GRecC = copyarrayfromC(Ghbc_c,3*m + 2)
    
    xhbc = []
    xhbc.append(xmh)
    xhbc.append(xmh)
    xhbc.append(x[0])
    xhbc.append(x[0] + 0.5*dx)
    for i in range(1,m-1):
        xhbc.append(x[i] - 0.5*dx)
        xhbc.append(x[i])
        xhbc.append(x[i] + 0.5*dx)
    
    xhbc.append(x[-1] - 0.5*dx)
    xhbc.append(x[-1])
    xhbc.append(xph)    
    xhbc.append(xph) 
    
    hhbc,uhbc,Ghbc = testsolSin(xhbc)
    
    xubc,t = makevar(xmh - dx,xph + 1.9*dx,dx,startt,endt,dt)
    
    uhhbc,uuhbc,uGhbc = testsolSin(xubc)
    
    dxs.append(dx)
    
    unormC =  norm(uFEMC -uuhbc,ord=1) / norm(uuhbc,ord=1) 
    
    hnormC =  norm(hRecC -hhbc,ord=1) / norm(hhbc,ord=1) 
    
    GnormC =  norm(GRecC -Ghbc,ord=1) / norm(Ghbc,ord=1) 

    unormsC.append(unormC)
    hnormsC.append(hnormC)
    GnormsC.append(GnormC)
        
s = wdir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unormsC[i])
        file1.write(s)
        
s = wdir + "h.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hnormsC[i])
        file1.write(s)


s = wdir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",GnormsC[i])
        file1.write(s)
        
"""

#second order test
"""
dx = 0.0001
l = 0.1
dt = l*dx
startx = 0
endx = 10.0 + dx
startt = 0.0
endt = 30 + (dt*0.9)  
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
theta = 1.4

gap = int(1.0/dt)
niBC = 3
        
xbc,t = makevar(startx - niBC*dx,endx +niBC*dx,dx,startt,endt,dt)

x = xbc[niBC:-niBC]


m = len(x)
        
hbc,ubc,Gbc = testsolSec(xbc)

h = hbc[niBC:-niBC]
hbeg = hbc[:niBC]
hend = hbc[-niBC:]

G = Gbc[niBC:-niBC]
Gbeg = Gbc[:niBC]
Gend = Gbc[-niBC:]

u = ubc[niBC:-niBC]
ubeg = ubc[:niBC]
uend = ubc[-niBC:]

xmh = x[0] - 0.5*dx
xph = x[-1] + 0.5*dx
umh = 10*xmh 
uph = 10*xph

uFEM = FEM2Assem(h,G,hbeg,hend,Gbeg,Gend,umh,uph,theta,dx)

xhbc,t = makevar(xmh,xph+0.9*dx,dx,startt,endt,dt)

hhbc,uhbc,Ghbc = testsolSec(xhbc)
"""
