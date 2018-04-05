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
   
def testsolSin(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    bed = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = sin(3*xp)
        h[i] = sin(10*xp) + 3
        bed[i] = sin(7*xp)
        G[i] = u[i]*h[i] - 30*(h[i])**2*cos(10*xp)*cos(3*xp) + 3*(h[i])**3*sin(3*xp) \
                +   u[i]*h[i]*10*cos(10*xp)*7*cos(7*xp) + 0.5*u[i]*h[i]*h[i]*(-49*sin(7*xp))  + u[i]*h[i]*(7*cos(7*xp))**2      
    return h,bed,u,G

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
    
def Quad(q,qmh,qph,i,dx):
    n = len(q)
    idx = 1.0 / dx
    if(i == n-1):
        qim1 = q[i-1]
        qi = q[i]
        qip1 = qph
    elif(i == 0):
        qim1 = qmh
        qi = q[i]
        qip1 = q[i +1]
    else:
        qim1 = q[i-1]
        qi = q[i]
        qip1 = q[i+1]
        
    ai = 0.5*idx*idx*(qip1 - 2*qi + qim1)
    bi = 0.5*idx*(qip1 - qim1)
    ci = qi
    
    qim1o2 = ai*(0.5*dx)**2 - bi*0.5*dx + ci
    qip1o2 = ai*(0.5*dx)**2 + bi*0.5*dx + ci
    
    return qim1o2,qi,qip1o2

def Cubic(q,i,dx):
    n = len(q)
    idx = 1.0 / dx
    i3 = 1.0/3.0
    i6 = i3*0.5
    i12 = i6*0.5
    
    
    #at boundaries
    if(i >= n-2):
        #Right edge do left one
        ai = i6*idx*idx*idx*(-q[i-3] +3*q[i-2] -3*q[i-1] + q[i]) 
        bi = 0.5*idx*idx*(-q[i-3] + 4*q[i-2] - 5*q[i-1] + 2*q[i])
        ci = i6*idx*(-2*q[i-3] + 9*q[i-2] - 18*q[i-1] + 11*q[i])
        di = q[i]

    elif(i <= 1):
        #Left edge do right one
        ai = i6*idx*idx*idx*(q[i+3] -3*q[i+2] +3*q[i+1] - q[i])
        bi = 0.5*idx*idx*(-q[i+3] + 4*q[i+2] - 5*q[i+1] + 2*q[i])
        ci = i6*idx*(2*q[i+3] - 9*q[i+2] + 18*q[i+1] - 11*q[i])
        di = q[i]

    else:
        #Interior do centered one       
        ai = i12*idx*idx*idx*(-q[i-2] + 2*q[i-1] -2*q[i+1] + q[i+2])
        bi = i6*idx*idx*(q[i-2] - q[i-1] - q[i+1] + q[i+2])
        ci = i12*idx*(q[i-2] - 8*q[i-1] + 8*q[i+1] - q[i+2])
        di = -q[i-2]*i6 + 2*i3*q[i-1] + 2*i3*q[i+1] - q[i+2]*i6
    
    qim1o2 = -ai*(0.5*dx)**3 + bi*(0.5*dx)**2 - ci*(0.5*dx) + di
    qim1o6 = -ai*(dx/6.0)**3 + bi*(dx/6.0)**2 - ci*(dx/6.0) + di
    qip1o2 = ai*(0.5*dx)**3 + bi*(0.5*dx)**2 + ci*(0.5*dx) + di
    qip1o6 = ai*(dx/6.0)**3 + bi*(dx/6.0)**2 + ci*(dx/6.0) + di
    
    return qim1o2,qim1o6,qip1o6,qip1o2        
    
    
    
#### CODES wrong not giving correct answers, tomorrow work on paper more    
def FEM2Assem(h,b,G,hbeg,hend,bbeg,bend,Gbeg,Gend,umh,uph,theta,dx):
    n = len(h)
    m = 2*n + 1
    idx = 1.0 / dx 
    
    
    uais = zeros(m-2)
    ubis = zeros(m-1)
    ucis = zeros(m)
    udis = zeros(m-1)
    ueis = zeros(m-2)
    
    nGis = zeros(m)
    
        
    j = 3
    for i in range(1,n-1):
        #First we should reconstruct all the ones we need for this step

        Gjmhp = Recon2P(G,0,0,0,i-1,theta,dx)      
        Gjphm = Recon2M(G,0,0,0,i,theta,dx)
        
        hjmhp = Recon2P(h,0,0,0,i-1,theta,dx)
        hjphm = Recon2M(h,0,0,0,i,theta,dx)
        
        bjmh,bjms,bjps,bjph = Cubic(b,i,dx)

        # G integral (RHS)
        Gintia11 = dx/6.0*(Gjmhp)
        Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm)
        Gintia31 = dx/6.0*(Gjphm)
        
        
        #uh integral
        uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm)
        uhintia12 = (1.0 /60.0)*dx*(4*hjmhp )
        uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
        
        uhintia21 = (1.0 /60.0)*dx*(4*hjmhp)
        uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm)
        uhintia23 = (1.0 /60.0)*dx*(4*hjphm)
        
        uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
        uhintia32 = (1.0 /60.0)*dx*(4*hjphm)
        uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm)
        
        #h3ux
        
        h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm)        
        h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
        h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)
        
        
        h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
        h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm)        
        h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)
        
        h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)        
        h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)        
        h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm)
        
        #h2 bx ux v
        h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
        h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
        h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)

        h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
        h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
        h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)

        h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)
        h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)
        h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)

        #h2 bx u vx
        h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
        h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
        h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)

        h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
        h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
        h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)

        h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)
        h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)
        h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)
                                        

        #h bx bx u v
        hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph)
        hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
        hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)


        hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
        hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph)
        hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)

        hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)                
        hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)
        hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph)       

        ## LHS 
        
        LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbx2uva11
        LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbx2uva12 
        LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbx2uva13
        LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbx2uva21
        LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbx2uva22 
        LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbx2uva23
        LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbx2uva31
        LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbx2uva32 
        LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbx2uva33
        
        uais[j-1] = uais[j-1] + LHSa31
        
        ubis[j-1] = ubis[j-1] + LHSa21
        ubis[j] = ubis[j] + LHSa32
        
        ucis[j-1] = ucis[j-1] + LHSa11
        ucis[j] = ucis[j] + LHSa22
        ucis[j+1] = ucis[j+1] + LHSa33
        
        udis[j-1] = udis[j-1] + LHSa12
        udis[j] = udis[j] + LHSa23
        
        ueis[j-1] = ueis[j-1] + LHSa13
        
        
        nGis[j-1] = nGis[j-1] + Gintia11  
        nGis[j] = nGis[j] + Gintia21 
        nGis[j+1] = nGis[j+1] + Gintia31 
        
        
        j = j + 2
        
    #First point
    j = 1
    i = 0
    
    Gjmhp = Recon2P(G,Gbeg[-1],G[0],G[1],i-1,theta,dx)
    Gjphm = Recon2M(G,Gbeg[-1],0,0,i,theta,dx)    
        
    hjmhp = Recon2P(h,hbeg[-1],h[0],h[1],i-1,theta,dx)
    hjphm = Recon2M(h,hbeg[-1],0,0,i,theta,dx)
    
    bjmh,bjms,bjps,bjph = Cubic(b,i,dx)
    
    
    
    # G integral (RHS)
    Gintia11 = dx/6.0*(Gjmhp)
    Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm)
    Gintia31 = dx/6.0*(Gjphm)
    
    
    #uh integral
    uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm)
    uhintia12 = (1.0 /60.0)*dx*(4*hjmhp )
    uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
    
    uhintia21 = (1.0 /60.0)*dx*(4*hjmhp)
    uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm)
    uhintia23 = (1.0 /60.0)*dx*(4*hjphm)
    
    uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
    uhintia32 = (1.0 /60.0)*dx*(4*hjphm)
    uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm)
    
    #h3ux
    
    h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm)        
    h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
    h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)
    
    
    h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
    h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm)        
    h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)
    
    h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)        
    h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)        
    h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm)
    
    #h2 bx ux v
    h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
    h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
    h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)

    h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
    h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
    h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)

    h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)
    h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)
    h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)

    #h2 bx u vx
    h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
    h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
    h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)

    h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
    h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
    h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)

    h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)
    h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)
    h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)
                                    

    #h bx bx u v
    hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph)
    hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
    hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)


    hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
    hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph)
    hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)

    hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)                
    hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)
    hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph)       

    ## LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbx2uva11
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbx2uva12 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbx2uva13
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbx2uva21
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbx2uva22 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbx2uva23
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbx2uva31
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbx2uva32 
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbx2uva33
    
    
    uais[j-1] = uais[j-1] + LHSa31
    
    ubis[j-1] = ubis[j-1] + LHSa21
    ubis[j] = ubis[j] + LHSa32
    
    ucis[j-1] = 1
    ucis[j] = ucis[j] + LHSa22
    ucis[j+1] = ucis[j+1] + LHSa33
    
    udis[j-1] = 0
    udis[j] = udis[j] + LHSa23
    
    ueis[j-1] = 0
    
    
    nGis[j-1] = umh
    nGis[j] = nGis[j] + Gintia21 
    nGis[j+1] = nGis[j+1] + Gintia31 
    
    
    #last point
    j = m-2
    i = n-1
    
    Gjmhp = Recon2P(G,0,0,0,i-1,theta,dx)
    Gjphm = Recon2M(G,0,0,0,i,theta,dx)
    
    hjmhp = Recon2P(h,0,0,0,i-1,theta,dx)
    hjphm = Recon2M(h,0,0,0,i,theta,dx)
    
    bjmh,bjms,bjps,bjph = Cubic(b,i,dx)
    
    
    # G integral (RHS)
    Gintia11 = dx/6.0*(Gjmhp)
    Gintia21 = dx/6.0*(2*Gjmhp + 2*Gjphm)
    Gintia31 = dx/6.0*(Gjphm)
    
    
    #uh integral
    uhintia11 = (1.0 /60.0)*dx*(7*hjmhp + hjphm)
    uhintia12 = (1.0 /60.0)*dx*(4*hjmhp )
    uhintia13 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
    
    uhintia21 = (1.0 /60.0)*dx*(4*hjmhp)
    uhintia22 = (1.0 /60.0)*dx*(16*hjmhp + 16*hjphm)
    uhintia23 = (1.0 /60.0)*dx*(4*hjphm)
    
    uhintia31 = (1.0 /60.0)*dx*(-hjmhp - hjphm)
    uhintia32 = (1.0 /60.0)*dx*(4*hjphm)
    uhintia33 = (1.0 /60.0)*dx*(hjmhp + 7*hjphm)
    
    #h3ux
    
    h3uxintia11 = (2.0/3.0)*idx*((79.0/120)*hjmhp*hjmhp*hjmhp +  (39.0/120)*hjmhp*hjmhp*hjphm + (3.0/24)*hjmhp*hjphm*hjphm + (7.0/120)*hjphm*hjphm*hjphm)        
    h3uxintia12 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
    h3uxintia13 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)
    
    
    h3uxintia21 = (2.0/3.0)*idx*((-23.0/30)*hjmhp*hjmhp*hjmhp +  (-3.0/10)*hjmhp*hjmhp*hjphm + (-3.0/30)*hjmhp*hjphm*hjphm + (-1.0/6)*hjphm*hjphm*hjphm)        
    h3uxintia22 = (2.0/3.0)*idx*((14.0/15)*hjmhp*hjmhp*hjmhp +  (6.0/15)*hjmhp*hjmhp*hjphm + (6.0/15)*hjmhp*hjphm*hjphm + (14.0/15)*hjphm*hjphm*hjphm)        
    h3uxintia23 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)
    
    h3uxintia31 = (2.0/3.0)*idx*((13.0/120)*hjmhp*hjmhp*hjmhp +  (-3.0/120)*hjmhp*hjmhp*hjphm + (-3.0/120)*hjmhp*hjphm*hjphm + (13.0/120)*hjphm*hjphm*hjphm)        
    h3uxintia32 = (2.0/3.0)*idx*((-1.0/6)*hjmhp*hjmhp*hjmhp +  (-3.0/30)*hjmhp*hjmhp*hjphm + (-3.0/10)*hjmhp*hjphm*hjphm + (-23.0/30)*hjphm*hjphm*hjphm)        
    h3uxintia33 = (2.0/3.0)*idx*((7.0/120)*hjmhp*hjmhp*hjmhp +  (3.0/24)*hjmhp*hjmhp*hjphm + (39.0/120)*hjmhp*hjphm*hjphm + (79.0/120)*hjphm*hjphm*hjphm)
    
    #h2 bx ux v
    h2bxuxva11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
    h2bxuxva12 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
    h2bxuxva13 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)

    h2bxuxva21 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
    h2bxuxva22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
    h2bxuxva23 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)

    h2bxuxva31 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)
    h2bxuxva32 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)
    h2bxuxva33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)

    #h2 bx u vx
    h2bxuvxa11 = -idx*((31.0/42)*hjmhp*hjmhp*bjmh+(19.0/280)*hjphm*hjmhp*bjmh+(19.0/280)*hjmhp*hjphm*bjmh+(23.0/1680)*hjphm*hjphm*bjmh+(-537.0/560)*hjmhp*hjmhp*bjms+(-33.0/560)*hjphm*hjmhp*bjms+(-33.0/560)*hjmhp*hjphm*bjms+(-3.0/280)*hjphm*hjphm*bjms+(39.0/140)*hjmhp*hjmhp*bjps+(-3.0/280)*hjphm*hjmhp*bjps+(-3.0/280)*hjmhp*hjphm*bjps+(3.0/560)*hjphm*hjphm*bjps+(-97.0/1680)*hjmhp*hjmhp*bjph+(1.0/560)*hjphm*hjmhp*bjph+(1.0/560)*hjmhp*hjphm*bjph+(-1.0/120)*hjphm*hjphm*bjph)
    h2bxuvxa12 = -idx*((71.0/210)*hjmhp*hjmhp*bjmh+(1.0/15)*hjphm*hjmhp*bjmh+(1.0/15)*hjmhp*hjphm*bjmh+(1.0/84)*hjphm*hjphm*bjmh+(-3.0/20)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(9.0/70)*hjphm*hjphm*bjms+(-3.0/14)*hjmhp*hjmhp*bjps+(-6.0/35)*hjphm*hjmhp*bjps+(-6.0/35)*hjmhp*hjphm*bjps+(-27.0/140)*hjphm*hjphm*bjps+(11.0/420)*hjmhp*hjmhp*bjph+(2.0/105)*hjphm*hjmhp*bjph+(2.0/105)*hjmhp*hjphm*bjph+(11.0/210)*hjphm*hjphm*bjph)
    h2bxuvxa13 = -idx*((-19.0/280)*hjmhp*hjmhp*bjmh+(-23.0/1680)*hjphm*hjmhp*bjmh+(-23.0/1680)*hjmhp*hjphm*bjmh+(-43.0/1680)*hjphm*hjphm*bjmh+(33.0/560)*hjmhp*hjmhp*bjms+(3.0/280)*hjphm*hjmhp*bjms+(3.0/280)*hjmhp*hjphm*bjms+(3.0/28)*hjphm*hjphm*bjms+(3.0/280)*hjmhp*hjmhp*bjps+(-3.0/560)*hjphm*hjmhp*bjps+(-3.0/560)*hjmhp*hjphm*bjps+(-21.0/80)*hjphm*hjphm*bjps+(-1.0/560)*hjmhp*hjmhp*bjph+(1.0/120)*hjphm*hjmhp*bjph+(1.0/120)*hjmhp*hjphm*bjph+(19.0/105)*hjphm*hjphm*bjph)

    h2bxuvxa21 = -idx*((-193.0/210)*hjmhp*hjmhp*bjmh+(-8.0/105)*hjphm*hjmhp*bjmh+(-8.0/105)*hjmhp*hjphm*bjmh+(-1.0/84)*hjphm*hjphm*bjmh+(171.0/140)*hjmhp*hjmhp*bjms+(9.0/140)*hjphm*hjmhp*bjms+(9.0/140)*hjmhp*hjphm*bjms+(0)*hjphm*hjphm*bjms+(-27.0/70)*hjmhp*hjmhp*bjps+(0)*hjphm*hjmhp*bjps+(0)*hjmhp*hjphm*bjps+(-9.0/140)*hjphm*hjphm*bjps+(1.0/12)*hjmhp*hjmhp*bjph+(1.0/84)*hjphm*hjmhp*bjph+(1.0/84)*hjmhp*hjphm*bjph+(8.0/105)*hjphm*hjphm*bjph)
    h2bxuvxa22 = -idx*((-41.0/105)*hjmhp*hjmhp*bjmh+(-3.0/35)*hjphm*hjmhp*bjmh+(-3.0/35)*hjmhp*hjphm*bjmh+(-4.0/105)*hjphm*hjphm*bjmh+(12.0/35)*hjmhp*hjmhp*bjms+(3.0/35)*hjphm*hjmhp*bjms+(3.0/35)*hjmhp*hjphm*bjms+(3.0/35)*hjphm*hjphm*bjms+(3.0/35)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(12.0/35)*hjphm*hjphm*bjps+(-4.0/105)*hjmhp*hjmhp*bjph+(-3.0/35)*hjphm*hjmhp*bjph+(-3.0/35)*hjmhp*hjphm*bjph+(-41.0/105)*hjphm*hjphm*bjph)
    h2bxuvxa23 = -idx*((8.0/105)*hjmhp*hjmhp*bjmh+(1.0/84)*hjphm*hjmhp*bjmh+(1.0/84)*hjmhp*hjphm*bjmh+(1.0/12)*hjphm*hjphm*bjmh+(-9.0/140)*hjmhp*hjmhp*bjms+(0)*hjphm*hjmhp*bjms+(0)*hjmhp*hjphm*bjms+(-27.0/70)*hjphm*hjphm*bjms+(0)*hjmhp*hjmhp*bjps+(9.0/140)*hjphm*hjmhp*bjps+(9.0/140)*hjmhp*hjphm*bjps+(171.0/140)*hjphm*hjphm*bjps+(-1.0/84)*hjmhp*hjmhp*bjph+(-8.0/105)*hjphm*hjmhp*bjph+(-8.0/105)*hjmhp*hjphm*bjph+(-193.0/210)*hjphm*hjphm*bjph)

    h2bxuvxa31 = -idx*((19.0/105)*hjmhp*hjmhp*bjmh+(1.0/120)*hjphm*hjmhp*bjmh+(1.0/120)*hjmhp*hjphm*bjmh+(-1.0/560)*hjphm*hjphm*bjmh+(-21.0/80)*hjmhp*hjmhp*bjms+(-3.0/560)*hjphm*hjmhp*bjms+(-3.0/560)*hjmhp*hjphm*bjms+(3.0/280)*hjphm*hjphm*bjms+(3.0/28)*hjmhp*hjmhp*bjps+(3.0/280)*hjphm*hjmhp*bjps+(3.0/280)*hjmhp*hjphm*bjps+(33.0/560)*hjphm*hjphm*bjps+(-43.0/1680)*hjmhp*hjmhp*bjph+(-23.0/1680)*hjphm*hjmhp*bjph+(-23.0/1680)*hjmhp*hjphm*bjph+(-19.0/280)*hjphm*hjphm*bjph)
    h2bxuvxa32 = -idx*((11.0/210)*hjmhp*hjmhp*bjmh+(2.0/105)*hjphm*hjmhp*bjmh+(2.0/105)*hjmhp*hjphm*bjmh+(11.0/420)*hjphm*hjphm*bjmh+(-27.0/140)*hjmhp*hjmhp*bjms+(-6.0/35)*hjphm*hjmhp*bjms+(-6.0/35)*hjmhp*hjphm*bjms+(-3.0/14)*hjphm*hjphm*bjms+(9.0/70)*hjmhp*hjmhp*bjps+(3.0/35)*hjphm*hjmhp*bjps+(3.0/35)*hjmhp*hjphm*bjps+(-3.0/20)*hjphm*hjphm*bjps+(1.0/84)*hjmhp*hjmhp*bjph+(1.0/15)*hjphm*hjmhp*bjph+(1.0/15)*hjmhp*hjphm*bjph+(71.0/210)*hjphm*hjphm*bjph)
    h2bxuvxa33 = -idx*((-1.0/120)*hjmhp*hjmhp*bjmh+(1.0/560)*hjphm*hjmhp*bjmh+(1.0/560)*hjmhp*hjphm*bjmh+(-97.0/1680)*hjphm*hjphm*bjmh+(3.0/560)*hjmhp*hjmhp*bjms+(-3.0/280)*hjphm*hjmhp*bjms+(-3.0/280)*hjmhp*hjphm*bjms+(39.0/140)*hjphm*hjphm*bjms+(-3.0/280)*hjmhp*hjmhp*bjps+(-33.0/560)*hjphm*hjmhp*bjps+(-33.0/560)*hjmhp*hjphm*bjps+(-537.0/560)*hjphm*hjphm*bjps+(23.0/1680)*hjmhp*hjmhp*bjph+(19.0/280)*hjphm*hjmhp*bjph+(19.0/280)*hjmhp*hjphm*bjph+(31.0/42)*hjphm*hjphm*bjph)
                                    

    #h bx bx u v
    hbx2uva11 = 2*idx*((1333.0/1344)*hjmhp*bjmh*bjmh+(443.0/6720)*hjphm*bjmh*bjmh+(-3141.0/2240)*hjmhp*bjms*bjmh+(-171.0/2240)*hjphm*bjms*bjmh+(1161.0/2240)*hjmhp*bjps*bjmh+(27.0/2240)*hjphm*bjps*bjmh+(-145.0/1344)*hjmhp*bjph*bjmh+(-11.0/6720)*hjphm*bjph*bjmh+(-3141.0/2240)*hjmhp*bjmh*bjms+(-171.0/2240)*hjphm*bjmh*bjms+(4617.0/2240)*hjmhp*bjms*bjms+(243.0/2240)*hjphm*bjms*bjms+(-1863.0/2240)*hjmhp*bjps*bjms+(-81.0/2240)*hjphm*bjps*bjms+(387.0/2240)*hjmhp*bjph*bjms+(9.0/2240)*hjphm*bjph*bjms+(1161.0/2240)*hjmhp*bjmh*bjps+(27.0/2240)*hjphm*bjmh*bjps+(-1863.0/2240)*hjmhp*bjms*bjps+(-81.0/2240)*hjphm*bjms*bjps+(891.0/2240)*hjmhp*bjps*bjps+(81.0/2240)*hjphm*bjps*bjps+(-27.0/320)*hjmhp*bjph*bjps+(-27.0/2240)*hjphm*bjph*bjps+(-145.0/1344)*hjmhp*bjmh*bjph+(-11.0/6720)*hjphm*bjmh*bjph+(387.0/2240)*hjmhp*bjms*bjph+(9.0/2240)*hjphm*bjms*bjph+(-27.0/320)*hjmhp*bjps*bjph+(-27.0/2240)*hjphm*bjps*bjph+(131.0/6720)*hjmhp*bjph*bjph+(13.0/1344)*hjphm*bjph*bjph)
    hbx2uva12 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
    hbx2uva13 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)


    hbx2uva21 = 2*idx*((103.0/336)*hjmhp*bjmh*bjmh+(3.0/70)*hjphm*bjmh*bjmh+(-183.0/560)*hjmhp*bjms*bjmh+(-3.0/140)*hjphm*bjms*bjmh+(3.0/112)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(-11.0/1680)*hjmhp*bjph*bjmh+(0)*hjphm*bjph*bjmh+(-183.0/560)*hjmhp*bjmh*bjms+(-3.0/140)*hjphm*bjmh*bjms+(243.0/560)*hjmhp*bjms*bjms+(0)*hjphm*bjms*bjms+(-81.0/560)*hjmhp*bjps*bjms+(0)*hjphm*bjps*bjms+(3.0/80)*hjmhp*bjph*bjms+(3.0/140)*hjphm*bjph*bjms+(3.0/112)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-81.0/560)*hjmhp*bjms*bjps+(0)*hjphm*bjms*bjps+(81.0/560)*hjmhp*bjps*bjps+(0)*hjphm*bjps*bjps+(-3.0/112)*hjmhp*bjph*bjps+(3.0/140)*hjphm*bjph*bjps+(-11.0/1680)*hjmhp*bjmh*bjph+(0)*hjphm*bjmh*bjph+(3.0/80)*hjmhp*bjms*bjph+(3.0/140)*hjphm*bjms*bjph+(-3.0/112)*hjmhp*bjps*bjph+(3.0/140)*hjphm*bjps*bjph+(-1.0/240)*hjmhp*bjph*bjph+(-3.0/70)*hjphm*bjph*bjph)
    hbx2uva22 = 2*idx*((101.0/420)*hjmhp*bjmh*bjmh+(29.0/420)*hjphm*bjmh*bjmh+(-6.0/35)*hjmhp*bjms*bjmh+(-3.0/35)*hjphm*bjms*bjmh+(-3.0/28)*hjmhp*bjps*bjmh+(-3.0/140)*hjphm*bjps*bjmh+(4.0/105)*hjmhp*bjph*bjmh+(4.0/105)*hjphm*bjph*bjmh+(-6.0/35)*hjmhp*bjmh*bjms+(-3.0/35)*hjphm*bjmh*bjms+(27.0/28)*hjmhp*bjms*bjms+(27.0/28)*hjphm*bjms*bjms+(-27.0/35)*hjmhp*bjps*bjms+(-27.0/35)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(-3.0/28)*hjphm*bjph*bjms+(-3.0/28)*hjmhp*bjmh*bjps+(-3.0/140)*hjphm*bjmh*bjps+(-27.0/35)*hjmhp*bjms*bjps+(-27.0/35)*hjphm*bjms*bjps+(27.0/28)*hjmhp*bjps*bjps+(27.0/28)*hjphm*bjps*bjps+(-3.0/35)*hjmhp*bjph*bjps+(-6.0/35)*hjphm*bjph*bjps+(4.0/105)*hjmhp*bjmh*bjph+(4.0/105)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(-3.0/28)*hjphm*bjms*bjph+(-3.0/35)*hjmhp*bjps*bjph+(-6.0/35)*hjphm*bjps*bjph+(29.0/420)*hjmhp*bjph*bjph+(101.0/420)*hjphm*bjph*bjph)
    hbx2uva23 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)

    hbx2uva31 = 2*idx*((-443.0/6720)*hjmhp*bjmh*bjmh+(-13.0/1344)*hjphm*bjmh*bjmh+(171.0/2240)*hjmhp*bjms*bjmh+(27.0/2240)*hjphm*bjms*bjmh+(-27.0/2240)*hjmhp*bjps*bjmh+(-9.0/2240)*hjphm*bjps*bjmh+(11.0/6720)*hjmhp*bjph*bjmh+(11.0/6720)*hjphm*bjph*bjmh+(171.0/2240)*hjmhp*bjmh*bjms+(27.0/2240)*hjphm*bjmh*bjms+(-243.0/2240)*hjmhp*bjms*bjms+(-81.0/2240)*hjphm*bjms*bjms+(81.0/2240)*hjmhp*bjps*bjms+(81.0/2240)*hjphm*bjps*bjms+(-9.0/2240)*hjmhp*bjph*bjms+(-27.0/2240)*hjphm*bjph*bjms+(-27.0/2240)*hjmhp*bjmh*bjps+(-9.0/2240)*hjphm*bjmh*bjps+(81.0/2240)*hjmhp*bjms*bjps+(81.0/2240)*hjphm*bjms*bjps+(-81.0/2240)*hjmhp*bjps*bjps+(-243.0/2240)*hjphm*bjps*bjps+(27.0/2240)*hjmhp*bjph*bjps+(171.0/2240)*hjphm*bjph*bjps+(11.0/6720)*hjmhp*bjmh*bjph+(11.0/6720)*hjphm*bjmh*bjph+(-9.0/2240)*hjmhp*bjms*bjph+(-27.0/2240)*hjphm*bjms*bjph+(27.0/2240)*hjmhp*bjps*bjph+(171.0/2240)*hjphm*bjps*bjph+(-13.0/1344)*hjmhp*bjph*bjph+(-443.0/6720)*hjphm*bjph*bjph)                
    hbx2uva32 = 2*idx*((-3.0/70)*hjmhp*bjmh*bjmh+(-1.0/240)*hjphm*bjmh*bjmh+(3.0/140)*hjmhp*bjms*bjmh+(-3.0/112)*hjphm*bjms*bjmh+(3.0/140)*hjmhp*bjps*bjmh+(3.0/80)*hjphm*bjps*bjmh+(0)*hjmhp*bjph*bjmh+(-11.0/1680)*hjphm*bjph*bjmh+(3.0/140)*hjmhp*bjmh*bjms+(-3.0/112)*hjphm*bjmh*bjms+(0)*hjmhp*bjms*bjms+(81.0/560)*hjphm*bjms*bjms+(0)*hjmhp*bjps*bjms+(-81.0/560)*hjphm*bjps*bjms+(-3.0/140)*hjmhp*bjph*bjms+(3.0/112)*hjphm*bjph*bjms+(3.0/140)*hjmhp*bjmh*bjps+(3.0/80)*hjphm*bjmh*bjps+(0)*hjmhp*bjms*bjps+(-81.0/560)*hjphm*bjms*bjps+(0)*hjmhp*bjps*bjps+(243.0/560)*hjphm*bjps*bjps+(-3.0/140)*hjmhp*bjph*bjps+(-183.0/560)*hjphm*bjph*bjps+(0)*hjmhp*bjmh*bjph+(-11.0/1680)*hjphm*bjmh*bjph+(-3.0/140)*hjmhp*bjms*bjph+(3.0/112)*hjphm*bjms*bjph+(-3.0/140)*hjmhp*bjps*bjph+(-183.0/560)*hjphm*bjps*bjph+(3.0/70)*hjmhp*bjph*bjph+(103.0/336)*hjphm*bjph*bjph)
    hbx2uva33 = 2*idx*((13.0/1344)*hjmhp*bjmh*bjmh+(131.0/6720)*hjphm*bjmh*bjmh+(-27.0/2240)*hjmhp*bjms*bjmh+(-27.0/320)*hjphm*bjms*bjmh+(9.0/2240)*hjmhp*bjps*bjmh+(387.0/2240)*hjphm*bjps*bjmh+(-11.0/6720)*hjmhp*bjph*bjmh+(-145.0/1344)*hjphm*bjph*bjmh+(-27.0/2240)*hjmhp*bjmh*bjms+(-27.0/320)*hjphm*bjmh*bjms+(81.0/2240)*hjmhp*bjms*bjms+(891.0/2240)*hjphm*bjms*bjms+(-81.0/2240)*hjmhp*bjps*bjms+(-1863.0/2240)*hjphm*bjps*bjms+(27.0/2240)*hjmhp*bjph*bjms+(1161.0/2240)*hjphm*bjph*bjms+(9.0/2240)*hjmhp*bjmh*bjps+(387.0/2240)*hjphm*bjmh*bjps+(-81.0/2240)*hjmhp*bjms*bjps+(-1863.0/2240)*hjphm*bjms*bjps+(243.0/2240)*hjmhp*bjps*bjps+(4617.0/2240)*hjphm*bjps*bjps+(-171.0/2240)*hjmhp*bjph*bjps+(-3141.0/2240)*hjphm*bjph*bjps+(-11.0/6720)*hjmhp*bjmh*bjph+(-145.0/1344)*hjphm*bjmh*bjph+(27.0/2240)*hjmhp*bjms*bjph+(1161.0/2240)*hjphm*bjms*bjph+(-171.0/2240)*hjmhp*bjps*bjph+(-3141.0/2240)*hjphm*bjps*bjph+(443.0/6720)*hjmhp*bjph*bjph+(1333.0/1344)*hjphm*bjph*bjph)       

    ## LHS 
    
    LHSa11 = uhintia11 + h3uxintia11 + h2bxuxva11 + h2bxuvxa11 + hbx2uva11
    LHSa12 = uhintia12 + h3uxintia12 + h2bxuxva12 + h2bxuvxa12 + hbx2uva12 
    LHSa13 = uhintia13 + h3uxintia13 + h2bxuxva13 + h2bxuvxa13 + hbx2uva13
    LHSa21 = uhintia21 + h3uxintia21 + h2bxuxva21 + h2bxuvxa21 + hbx2uva21
    LHSa22 = uhintia22 + h3uxintia22 + h2bxuxva22 + h2bxuvxa22 + hbx2uva22 
    LHSa23 = uhintia23 + h3uxintia23 + h2bxuxva23 + h2bxuvxa23 + hbx2uva23
    LHSa31 = uhintia31 + h3uxintia31 + h2bxuxva31 + h2bxuvxa31 + hbx2uva31
    LHSa32 = uhintia32 + h3uxintia32 + h2bxuxva32 + h2bxuvxa32 + hbx2uva32 
    LHSa33 = uhintia33 + h3uxintia33 + h2bxuxva33 + h2bxuvxa33 + hbx2uva33   
   
    uais[j-1] = 0
    
    ubis[j-1] = ubis[j-1] + LHSa21
    ubis[j] = 0
    
    ucis[j-1] = ucis[j-1] + LHSa11
    ucis[j] = ucis[j] + LHSa22
    ucis[j+1] = 1
    
    udis[j-1] = udis[j-1] + LHSa12
    udis[j] = udis[j] + LHSa23
    
    ueis[j-1] = ueis[j-1] + LHSa13
    
    
    nGis[j-1] = nGis[j-1] + Gintia11  
    nGis[j] = nGis[j] + Gintia21 
    nGis[j+1] = uph
    
    
    
    #This gives coefficients from one equation, but we have many equations
    #NUs = diag(uais,k=-2) + diag(ubis,k=-1) + diag(ucis,k=0) + diag(udis,k=1) + diag(ueis,k=2)
    
    ue = pentadiagsolve(uais,ubis,ucis,udis,ueis,nGis)
    
    return ue,uais,ubis,ucis,udis,ueis,nGis
               


#Python Code
"""
#Sin test, very little error good
dxs = []
unorms = []
dunorms = []
wdir = "../../../data/raw/FEMTESTSinbed/"  
if not os.path.exists(wdir):
    os.makedirs(wdir)
for i in range(5,15):
    dx = 1.0/2**i
    l = 0.1
    dt = l*dx
    startx = -10
    endx = 10 + 0.9*dx
    startt = 0.0
    endt = 30 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            
            
    g = 9.81
    theta = 1.4
    
    gap = int(1.0/dt)
    nBC = 3
            
    xbc,t = makevar(startx - nBC*dx,endx +nBC*dx,dx,startt,endt,dt)
    
    x = xbc[nBC:-nBC]
    
    
    n = len(x)
    idx = 1.0 / dx
            
    hbc,bbc,ubc,Gbc = testsolSin(xbc)
    
    h = hbc[nBC:-nBC]
    hbeg = hbc[:nBC]
    hend = hbc[-nBC:]
    
    G = Gbc[nBC:-nBC]
    Gbeg = Gbc[:nBC]
    Gend = Gbc[-nBC:]
    
    u = ubc[nBC:-nBC]
    ubeg = ubc[:nBC]
    uend = ubc[-nBC:]
    
    b = bbc[nBC:-nBC]
    bbeg = bbc[:nBC]
    bend = bbc[-nBC:]
    
    xmh = x[0] - 0.5*dx
    xph = x[-1] + 0.5*dx
    umh = sin(3*xmh)
    uph = sin(3*xph)
    
    xebeg,t0beg = makevar(xmh -nBC*dx,xmh + 0.1*dx ,dx,startt,endt,dt)
    xeend,t0beg = makevar(xph,xph + nBC*dx + 0.1*dx ,dx,startt,endt,dt)
    
    xh,t0 = makevar(xmh,xph+0.1*dx ,0.5*dx,startt,endt,dt)
    
    hebeg,bebeg,uebeg,Gebeg = testsolSin(xebeg)
    heend,beend,ueend,Geend = testsolSin(xeend)
    
    he,be,ue,Ge = testsolSin(xh)
    
    
    Gpy  = getGfromupy(h,u,b,ubeg[-1],uend[0],hbeg[-1],hend[0],bbeg[-1],bend[0],dx)
    
    ueFEMpy,uais,ubis,ucis,udis,ueis,nGis  = FEM2Assem(h,b,G,hbeg,hend,bbeg,bend,Gbeg,Gend,umh,uph,theta,dx)
    uFEMpy = ueFEMpy[1:-1:2]
    
    
    normu = norm(array(ueFEMpy) - array(ue),ord=1) / norm(ue,ord=1)

    
    dxs.append(dx)
    unorms.append(normu)
  
s = wdir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unorms[i])
        file1.write(s)
"""


#C Code
 
#Sin test, very little error good
dxs = []
unorms = []
hnorms = []
Gnorms = []
bnorms = []
dunorms = []
wdir = "../../../data/raw/FEMTESTSinbedC/"  
if not os.path.exists(wdir):
    os.makedirs(wdir)
for i in range(7,8):
    dx = 1.0
    l = 0.1
    dt = l*dx
    startx = 0
    endx = 2*dx + 0.9*dx
    startt = 0.0
    endt = 30 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            
            
    g = 9.81
    theta = 1.4
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    
    
    n = len(x)
    idx = 1.0 / dx
            
    h,b,u,G = testsolSin(x)
    
    #Make x's for hhbc, ubc and bedhbc (ubc is xh)
    xhbc = []
    xubc = []
    xbhbc = []
    for i in range(len(xG)):
        if(i ==0 or i == len(xG) -1):
            xhbc.append(xG[i] - 0.5*dx)
            xhbc.append(xG[i] + 0.5*dx)
        else:
            xhbc.append(xG[i] - 0.5*dx)
            xhbc.append(xG[i])
            xhbc.append(xG[i] + 0.5*dx) 
        xubc.append(xG[i] - 0.5*dx)
        xubc.append(xG[i])
        xubc.append(xG[i] + 0.5*dx)
        xbhbc.append(xG[i] - 0.5*dx)
        xbhbc.append(xG[i] - dx/6.0)
        xbhbc.append(xG[i] + dx/6.0)
        xbhbc.append(xG[i] + 0.5*dx)
            
    xubc = array(xubc)    
    xhbc = array(xhbc)
    xbhbc = array(xbhbc)
    
    xuMbeg = xubc[:unBC]
    xuMend = xubc[-unBC:]
    xGhMbeg = xubc[:GhnBC]
    xGhMend = xubc[-GhnBC:]
    xbMbeg = xubc[:bnBC]
    xbMend = xubc[-bnBC:]
    
    hMbeg,b_ta,u_ta,GMbeg = testsolSin(xGhMbeg)
    
    bmbeg = sin
    
    
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(b)
    u_c = mallocPy(n)
    ubc_c = mallocPy(2*n + 1)
    hhbc_c = mallocPy(3*n)
    Ghbc_c = mallocPy(3*n)
    bedhbc_c = mallocPy(4*n)
    
    #hMbeg_c = copyarraytoC(hMbeg)
    #hMend_c = copyarraytoC(hMend)  
    #GMbeg_c = copyarraytoC(GMbeg)
    #GMend_c = copyarraytoC(GMend) 
    #uMbeg_c = copyarraytoC(uMbeg)
    #uMend_c = copyarraytoC(uMend)
    
"""    
    evolvewrap(G_c,h_c,bed_c,hebeg_c,heend_c ,Gebeg_c,Geend_c,uebeg_c,ueend_c,bedbeg_c,bedend_c,g,dx,dt,n,nBC,theta,ubc_c, hhbc_c,Ghbc_c,bedhbc_c)
    ubcC = copyarrayfromC(ubc_c,2*n + 1)
    hhbcC = copyarrayfromC(hhbc_c,3*n)
    GhbcC = copyarrayfromC(Ghbc_c,3*n)
    bedhbcC = copyarrayfromC(bedhbc_c,4*n)
    
    normG = norm(GhbcC - GhbcA,ord=1) / norm(GhbcA,ord=1)
    normh = norm(hhbcC - hhbcA,ord=1) / norm(hhbcA,ord=1)
    normb = norm(bedhbcC - bbhbc,ord=1) / norm(bbhbc,ord=1)
    normu = norm(ubcC - ue,ord=1) / norm(ue,ord=1)

    
    dxs.append(dx)
    unorms.append(normu)
    hnorms.append(normh)
    Gnorms.append(normG)
    bnorms.append(normb)
    
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bedhbc_c)
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(u_c)
    deallocPy(hebeg_c)
    deallocPy(Gebeg_c)
    deallocPy(uebeg_c)
    deallocPy(bedbeg_c)
    deallocPy(heend_c)
    deallocPy(Geend_c)
    deallocPy(ueend_c)
    deallocPy(bedend_c)    
  
s = wdir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unorms[i])
        file1.write(s)
        
s = wdir + "b.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",bnorms[i])
        file1.write(s)
        
s = wdir + "h.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hnorms[i])
        file1.write(s)
        
s = wdir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",Gnorms[i])
        file1.write(s)
"""




     
