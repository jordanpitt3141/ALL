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

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


    
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):

        h[i] = a0 + a1*sin(a2*x[i])*exp(a3*t)
        u[i] = a4*cos(a5*x[i])*exp(a6*t)
        b[i] = a7*sin(a8*x[i])
        w[i] = h[i] + b[i]
        
        hxi = a1*a2*cos(a2*x[i])*exp(a3*t)
        uxi = -a4*a5*sin(a5*x[i])*exp(a6*t)

        uxxi = -a4*a5**2*cos(a5*x[i])*exp(a6*t)
        
        bxi = a7*a8*cos(a8*x[i])
        bxxi = -a7*a8**2*sin(a8*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w
    
#G has error that flattens it out, mainly from the boundary


#Forcing Problem    
wdir = "../../../../../../../data/raw/Forced/P1P2P3BedFEM/PerandExp/SmollGrowth0p1/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(15):

    """
    a0 = 1.1
    a1 = 0.2
    a2 = 0.3
    a3 = 0.4
    a4 = 0.5
    a5 = 0.6
    a6 = 0.1
    a7 = 2.0
    a8 = 0.9
"""
    

    a0 = 10
    a1 = 0.1
    a2 = 0.2
    a3 = 0.1
    a4 = 0.3
    a5 = 0.4
    a6 = 0.1
    a7 = -0.1
    a8 = 0.1

    
    width = 50
    
    g = 9.81
    
    startx = -width/2
    endx = width/2
    startt = 0.0
    endt =  0.01
    
    dx = width / (2.0)**(j)
    l =  0.5 / (a4*exp(a6*endt) + sqrt(g*(a0 + a1*exp(a3*endt))))
    dt = l*dx
            
    szoomx = startx
    ezoomx = endx
    
    t = startt
            
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
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    
    xbMend = array([x[-1] + 0.5*dx, x[-1] + 5*dx/6.0, x[-1] + 7*dx/6.0, x[-1] + 1.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - 7*dx/6.0,x[0] - 5*dx/6.0 , x[0] -0.5*dx])
                
        
    
    h,u,G,b,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    bMbeg = a7*sin(a8*xbMbeg)
    bMend = a7*sin(a8*xbMend)
    
    hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t+dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    hMend1,uMend1,GMend1,bta,wMend1 = ForcedbedM(xhuMend,t+dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    b_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
   
 
    hMbeg1_c = copyarraytoC(hMbeg1)
    hMend1_c = copyarraytoC(hMend1)
    wMbeg1_c = copyarraytoC(wMbeg1)
    wMend1_c = copyarraytoC(wMend1)
    GMbeg1_c = copyarraytoC(GMbeg1)
    GMend1_c = copyarraytoC(GMend1) 
    uMbeg1_c = copyarraytoC(uMbeg1)
    uMend1_c = copyarraytoC(uMend1)
   
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bhbc_c = mallocPy(nbhbc)
    
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg1_c,hMend1_c,wMbeg1_c,wMend1_c,GMbeg1_c,GMend1_c,uMbeg1_c,uMend1_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8)
        t = t + dt
        ts.append(t)
        print(t)
        
        copywritearraytoC(hMbeg1,hMbeg_c)
        copywritearraytoC(hMend1,hMend_c)
        copywritearraytoC(uMbeg1,uMbeg_c)
        copywritearraytoC(uMend1,uMend_c)
        copywritearraytoC(GMbeg1,GMbeg_c)
        copywritearraytoC(GMend1,GMend_c)
        
        hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t+dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
        hMend1,uMend1,GMend1,bta,wMend1 = ForcedbedM(xhuMend,t+dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
        
        copywritearraytoC(hMbeg1,hMbeg1_c)
        copywritearraytoC(hMend1,hMend1_c)
        copywritearraytoC(uMbeg1,uMbeg1_c)
        copywritearraytoC(uMend1,uMend1_c)
        copywritearraytoC(GMbeg1,GMbeg1_c)
        copywritearraytoC(GMend1,GMend1_c)
        
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    hA,uA,GA,bA,wA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,g,dx)
    
    hnorm = norm(hC - hA, ord=2)/ norm(hC, ord=2)
    unorm = norm(uC - uA, ord=2)/ norm(uC, ord=2)
    Gnorm = norm(GC - GA, ord=2)/ norm(GC, ord=2)
    
    
    
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
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)


"""
#Soliton Problem

wdir = "../../../../../../data/raw/Forced/P1P2P3BedFEM/GaussBedO/Soltest/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(20):
    g =9.81
    a0 = 1.0
    a1 = 1.0
    
    width = 50
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (sqrt(g*(a0 + a1)))
    dt = l*dx
    startx = -width/2
    endx = width/2 + 0.9*dx
    startt = 0.0
    endt = 0.1
            
    szoomx = startx
    ezoomx = endx
    
    t = 0
            
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
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
                
        
    h,u,G,b = solitoninit(n,a0,a1,g,x,startt,0,dx)
    w = h + b
    
    
    
    print(t)
    
    hMbeg = a0*ones(GhnBC)
    hMend = a0*ones(GhnBC)
    
    wMbeg = a0*ones(GhnBC)
    wMend = a0*ones(GhnBC)
    
    uMbeg = zeros(GhnBC)
    uMend = zeros(GhnBC)
    
    GMbeg = zeros(GhnBC)
    GMend = zeros(GhnBC)
    
    bMbeg = zeros(bnBC)
    bMend = zeros(bnBC)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    b_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
   
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bhbc_c = mallocPy(nbhbc)
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t)
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    hA,uA,GA,bA = solitoninit(n,a0,a1,g,x,t,0,dx)
    wA = hA + bA
    
    hnorm = norm(hC - hA, ord=2)/ norm(hC, ord=2)
    unorm = norm(uC - uA, ord=2)/ norm(uC, ord=2)
    Gnorm = norm(GC - GA, ord=2)/ norm(GC, ord=2)
    
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
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
"""