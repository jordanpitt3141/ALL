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

   
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        u[i] = a6*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)
        b[i] = a8*sin(a9*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        uxi = -a6/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)

        uxxi = -a6/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)*exp(a7*t)
        
        bxi = a8*a9*cos(a9*x[i]) 
        bxxi = -a8*a9**2*sin(a9*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

#Forcing Problem    
wdir = "../../../../data/2018/raw/Thesis/TestReconQuartLong/Wet/FEVM2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(3,16):
    
    wdirji = wdir + str(ki) + "/"
    if not os.path.exists(wdirji):
        os.makedirs(wdirji)
    
    a8 = 1.0
    a9 = 2*pi/50.0
    
    width = 2*(2*pi/a9)
        
    a0 = 1.0
    a1 = 0.5
    a2 =  ((2*pi) / a9)/10.0
    a3 = -pi/2.0/a9 -width/4.0
    a4 = width/2**6
    a5 = 0.0
    a6 = a1
    a7 = 0.0

    
    g = 9.81
    
    startx = -pi/2.0/a9 -width
    endx = -pi/2.0/a9 +width
    startt = 0.0
    endt = 0.1#(2*pi/a9) / a2
    
    dx = width / (2.0)**(ki)
    l =  0.5 / (a6*exp(a7*endt) + sqrt(g*(a0 + a1*exp(a5*endt))))
    dt = l*dx
            
    szoomx = startx
    ezoomx = endx
    
    t = startt
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
            
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate(([x[0] - dx],x,[x[-1]+ dx]))
    
    xbc = []
    for i in range(len(xG)):
        if i == 0:
            xbc.append(xG[i] - 0.5*dx)
            xbc.append(xG[i] - dx/6.0)
            xbc.append(xG[i] + dx/6.0)
            xbc.append(xG[i] + 0.5*dx)
        else:
            xbc.append(xG[i] - dx/6.0)
            xbc.append(xG[i] + dx/6.0)
            xbc.append(xG[i] + 0.5*dx)
        
    
    xMbeg = [x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx]
    xMend = [x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx]
    n = len(x)
    nMBC = 7
    nbcBC = 4
    nbc = 3*n + 1 + 2*(nbcBC -1)
 
    hbc_c =  mallocPy(nbc)
    ubc_c =  mallocPy(nbc)
    Gbc_c =  mallocPy(nbc)
    bbc_c =  mallocPy(nbc)
    wbc_c =  mallocPy(nbc)
    
    h,u,G,b,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hMbeg,uMbeg,GMbeg,bMbeg,wMbeg = ForcedbedM(xMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hMend,uMend,GMend,bMend,wMend = ForcedbedM(xMend,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    h_c = copyarraytoC(h)
    u_c = copyarraytoC(u)
    G_c = copyarraytoC(G)
    b_c = copyarraytoC(b)
    w_c = copyarraytoC(w)
    
    hMbeg_c = copyarraytoC(hMbeg)
    uMbeg_c = copyarraytoC(uMbeg)
    GMbeg_c = copyarraytoC(GMbeg)
    bMbeg_c = copyarraytoC(bMbeg)
    wMbeg_c = copyarraytoC(wMbeg)
    
    hMend_c = copyarraytoC(hMend)
    uMend_c = copyarraytoC(uMend)
    GMend_c = copyarraytoC(GMend)
    bMend_c = copyarraytoC(bMend)
    wMend_c = copyarraytoC(wMend)
    
    ReconQuart(h_c, hMbeg_c, hMend_c, n,nMBC,nbcBC,nbc, hbc_c,dx)
    ReconQuart(w_c, wMbeg_c, wMend_c, n,nMBC,nbcBC,nbc, wbc_c,dx)
    ReconQuart(G_c, GMbeg_c, GMend_c, n,nMBC,nbcBC,nbc, Gbc_c,dx)
    ReconQuart(u_c, uMbeg_c, uMend_c, n,nMBC,nbcBC,nbc, ubc_c,dx)
    ReconQuart(b_c, bMbeg_c, bMend_c, n,nMBC,nbcBC,nbc, bbc_c,dx)


    hbcC = copyarrayfromC(hbc_c,nbc)
    ubcC = copyarrayfromC(ubc_c,nbc)
    GbcC = copyarrayfromC(Gbc_c,nbc)
    bbcC = copyarrayfromC(bbc_c,nbc)
    wbcC = copyarrayfromC(wbc_c,nbc)
    
    
    hAbc,uAbc,GAbc,bAbc,wAbc= ForcedbedM(xbc,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)

    hnorm = norm(hbcC -hAbc, ord=1)/ norm(hAbc, ord=1)
    Gnorm = norm(GbcC -GAbc, ord=1)/ norm(GAbc, ord=1)
    unorm = norm(ubcC -uAbc, ord=1)/ norm(uAbc, ord=1)
    bnorm = norm(bbcC -bAbc, ord=1)/ norm(bAbc, ord=1)
    wnorm = norm(wbcC -wAbc, ord=1)/ norm(wAbc, ord=1)

    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hnorm)
        file1.write(s)
    
    s = wdir + "GL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",unorm)
        file1.write(s)
    
    s = wdir + "bL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",bnorm)
        file1.write(s)       
 
    s = wdir + "wL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",wnorm)
        file1.write(s)  


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