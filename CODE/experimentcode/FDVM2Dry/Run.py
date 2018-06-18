# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe

from scipy.optimize import bisect

#gives exact up to linears, so is second order accurate huzzah    
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

    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def makeX(sx,ex,dx): 
    x = arange(sx, ex, dx)
    return x 

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(a0,a1,g,x,t0,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    w = zeros(n)
    b = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        w[i] = h[i]
        u[i] =  c* (1 - a0 / h[i])
        #G[i] = 2.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**4*h[i] + h[i]*u[i] - 4.0/3*a0*a1**2*c*k**2*sech(k*(x[i] - c*t0))**4*tanh(k*(x[i] - c*t0))**2 - 4.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**2*h[i]*tanh(k*(x[i] - c*t0))**2
    G = getGfromupy(h,u,0.0,0.0,a0,a0,dx)
    return h,u,G,w,b



sdir = "./tests/Solitont3/"
if not os.path.exists(sdir):
    os.makedirs(sdir)

for ki in range(6,20):
    theta = 1.2

    a0 = 1.0
    a1 = 0.7
        
    g = 9.81
    sx = -50
    ex = 50
    dx = 100.0/ 2**ki
    
    st = 0
    et = 1
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    
    nMBC = 3
    nEBC = 3
    nCBC = 1
    
    x = makeX(sx,ex + 0.1*dx,dx)
    n= len(x)
    xMbeg  = array([x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx])
    xMend  = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xCbc = concatenate(([x[0] - dx], x, [x[-1] + dx]))
    
    nMbc = 3*n + 2*nMBC
    nEbc = 2*n - 1 + 2*nEBC
    nCbc = n + 2*nCBC
    
    
    h,u,G,w,b = solitoninit(a0,a1,g,x,st,dx)
    
    hMbeg = h[0]*ones(nMBC)
    wMbeg = hMbeg
    GMbeg = zeros(nMBC)
    bMbeg = zeros(nMBC)
    
    uEbeg = zeros(nEBC)
    duEbeg = zeros(nEBC)
    
    ddbCbeg = zeros(nCBC)
    
    hMend = h[-1]*ones(nMBC)
    wMend = hMend
    GMend = zeros(nMBC)
    bMend = zeros(nMBC)
    
    uEend = zeros(nEBC)
    duEend = zeros(nEBC)
    
    ddbCend = zeros(nCBC)    
    
      

    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    b_c = copyarraytoC(b)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    uEbeg_c = copyarraytoC(uEbeg)
    uEend_c = copyarraytoC(uEend)
    duEbeg_c = copyarraytoC(duEbeg)
    duEend_c = copyarraytoC(duEend)
    ddbCbeg_c = copyarraytoC(ddbCbeg)
    ddbCend_c = copyarraytoC(ddbCend)
    
    hMbc_c = mallocPy(nMbc)
    GMbc_c = mallocPy(nMbc)
    wMbc_c = mallocPy(nMbc)
    bMbc_c = mallocPy(nMbc)
    
    duEbc_c = mallocPy(nEbc)
    uEbc_c = mallocPy(nEbc)
    
    ddbCbc_c = mallocPy(nCbc)
    
    
    t = 0.0
    #Just an FEM solve here
    while t < et: 
        evolvewrapBC(h_c,G_c,b_c,hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c,ddbCend_c,hMbeg_c,GMbeg_c,wMbeg_c,duEbeg_c,uEbeg_c,hMend_c,GMend_c,wMend_c,duEend_c,uEend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta)
        t = t + dt
        print(t)
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
    uEbcC = copyarrayfromC(uEbc_c,nEbc)
    uC = uEbcC[nEBC:-nEBC:2]
    
    hA,uA,GA,wA,bA = solitoninit(a0,a1,g,x,t,dx)
    
    hL1 = norm(hC - hA,ord=1)/ norm(hA,ord=1)
    GL1 = norm(GC - GA,ord=1)/ norm(GA,ord=1)    
    uL1 = norm(uC - uA,ord=1)/ norm(uA,ord=1)
    
    deallocPy(hMbc_c)
    deallocPy(GMbc_c)
    deallocPy(wMbc_c)
    deallocPy(bMbc_c)
    deallocPy(duEbc_c)
    deallocPy(uEbc_c)
    deallocPy(ddbCbc_c)
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(b_c)
    
    deallocPy(hMbeg_c)
    deallocPy(hMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
    deallocPy(GMbeg_c)
    deallocPy(GMend_c)
    deallocPy(bMbeg_c)
    deallocPy(bMend_c)
    deallocPy(uEbeg_c)
    deallocPy(uEend_c)
    deallocPy(duEbeg_c)
    deallocPy(duEend_c)
    deallocPy(ddbCbeg_c)
    deallocPy(ddbCend_c)


    
    s = sdir + "hL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",hL1)
            file1.write(s)

    s = sdir + "GL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",GL1)
            file1.write(s)

    s = sdir + "uL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uL1)
            file1.write(s)