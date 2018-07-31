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

def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    ux = zeros(n)
    bxx = zeros(n)
    
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        ux[i] = uxi
        bxx[i] = bxxi
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w,ux,bxx

sdir = "./tests/FDsolverDryT/"
if not os.path.exists(sdir):
    os.makedirs(sdir)

for ki in range(8,9):

    a0 = 0.0
    a1 = 0.2
    a2 = 1.3
    a3 = 0.4
    a4 = 1.5
    a5 = 0.6
    a6 = 0.01
    a7 = 0.2
        
    g = 9.81
    sx = -50
    ex = 50
    dx = float(abs(ex - sx))/ 2**ki
    
    st = 0
    et = 0
    dt = 0
    
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
    
    
    h,u,G,b,w,uxta,bxxta = ForcedbedM(x,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    hMbeg,uEbeg,GMbeg,bMbeg,wMbeg,duEbeg,bxxta = ForcedbedM(xMbeg,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    ddbCbeg = [bxxta[1]]
    hMend,uEend,GMend,bMend,wMend,duEend,bxxta = ForcedbedM(xMend,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    ddbCend = [bxxta[1]]
      

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
    
    
    theta = 1.2
    
    xMbc = []
    xEbc = []
    for i in range(len(xCbc)):
        if(i ==0):           
            xEbc.append(xCbc[i] - 0.5*dx)
            xEbc.append(xCbc[i])
            xEbc.append(xCbc[i] + 0.5*dx)            
        else:
            xEbc.append(xCbc[i])
            xEbc.append(xCbc[i] + 0.5*dx)
            
        xMbc.append(xCbc[i] - 0.5*dx)
        xMbc.append(xCbc[i])
        xMbc.append(xCbc[i] + 0.5*dx)
            
    xMbc = array(xMbc)    
    xEbc = array(xEbc)
    xCbc = array(xCbc)
    
    
    
    edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
    
    hMbcC = copyarrayfromC(hMbc_c,nMbc)
    GMbcC = copyarrayfromC(GMbc_c,nMbc)
    wMbcC = copyarrayfromC(wMbc_c,nMbc)
    bMbcC = copyarrayfromC(bMbc_c,nMbc)
    duEbcC = copyarrayfromC(duEbc_c,nEbc)
    
    uEbcC = copyarrayfromC(uEbc_c,nEbc)
    
    ddbCbcC = copyarrayfromC(ddbCbc_c,nCbc)
    
    hMbcA,uMbcA,GMbcA,bMbcA,wMbcA,uxMbcA,bxxMbcA = ForcedbedM(xMbc,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hEbcA,uEbcA,GEbcA,bEbcA,wEbcA,uxEbcA,bxxEbcA = ForcedbedM(xEbc,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hCbcA,uCbcA,GCbcA,bCbcA,wCbcA,uxCbcA,bxxCbcA = ForcedbedM(xCbc,st,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    hL1 = norm(hMbcC - hMbcA,ord=1)/ norm(hMbcA,ord=1)
    GL1 = norm(GMbcC - GMbcA,ord=1)/ norm(GMbcA,ord=1)
    wL1 = norm(wMbcC - wMbcA,ord=1)/ norm(wMbcA,ord=1)
    bL1 = norm(bMbcC - bMbcA,ord=1)/ norm(bMbcA,ord=1)
    
    uL1 = norm(uEbcC - uEbcA,ord=1)/ norm(uEbcA,ord=1)
    duL1 = norm(duEbcC - uxEbcA,ord=1)/ norm(uxEbcA,ord=1)
    
    ddbL1 = norm(ddbCbcC - bxxCbcA,ord=1)/ norm(bxxCbcA,ord=1)
    
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
            
    s = sdir + "wL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",wL1)
            file1.write(s)
            
    s = sdir + "bL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",bL1)
            file1.write(s)

    s = sdir + "uL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uL1)
            file1.write(s)
            
    s = sdir + "duL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",duL1)
            file1.write(s)

    s = sdir + "ddbL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",ddbL1)
            file1.write(s)