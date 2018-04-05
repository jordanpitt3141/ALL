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

def minmodpy(a, b, c):
    if((a > 0) and (b>0) and (c>0)):
        return min(a,b,c)
    elif((a < 0) and (b<0) and (c<0)):
        return max(a,b,c)
    else:
        return 0.0
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    #return y1  + (xi)*(y2 - y1)/(x2 - x1)  
    return y1  + (xi)*(y2 - y0)/(x2 - x0)  

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
   


def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    b = []
    
    cxv = x[0] - dx
    b.append(cos(a5*(cxv - 0.5*dx)))
    b.append(cos(a5*(cxv - dx/6.0)))
    bavg= (sin(a5*(cxv + 0.5*dx)) - sin(a5*(cxv - 0.5*dx)))/(a5*dx)
    b.append(bavg)
    b.append(cos(a5*(cxv + dx/6.0)))
    b.append(cos(a5*(cxv + 0.5*dx)))

    
    for i in range(n):
        b.append(cos(a5*(x[i] - dx/6.0)))
        bavg= (sin(a5*(x[i] + 0.5*dx)) - sin(a5*(x[i] - 0.5*dx)))/(a5*dx)
        b.append(bavg)
        b.append(cos(a5*(x[i] + dx/6.0)))
        b.append(cos(a5*(x[i] + 0.5*dx)))         
        
        
        h[i] = a0 + sin(a1*x[i])*exp(a2*t)
        u[i] = cos(a3*x[i])*exp(a4*t)
        bi = cos(a5*x[i])
        w[i] = h[i] + bi
        
        hxi = a1*exp(a2*t)*cos(a1*x[i])
        uxi = -a3*exp(a4*t)*sin(a3*x[i]) 
        bxi = -a5*sin(a5*x[i]) 
        
        uxxi = -a3**2*exp(a4*t)*cos(a3*x[i])
        bxxi = -a5**2*cos(a5*x[i])
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
    
    cxv = x[-1] + dx
    b.append(cos(a5*(cxv - dx/6.0)))
    bavg= (sin(a5*(cxv + 0.5*dx)) - sin(a5*(cxv - 0.5*dx)))/(a5*dx)
    b.append(bavg)
    b.append(cos(a5*(cxv + dx/6.0)))
    b.append(cos(a5*(cxv + 0.5*dx)))  
    
    return h,u,G,b,w
 

"""
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    b = []
    
    cxv = x[0] - dx
    b.append(0)
    b.append(0)
    bavg= 0
    b.append(bavg)
    b.append(0)
    b.append(0)

    
    for i in range(n):
        b.append(0)
        bavg= 0
        b.append(0)
        b.append(0)
        b.append(0)         
        
        
        h[i] = a0 + sin(a1*x[i])*exp(a2*t)
        u[i] = cos(a3*x[i])*exp(a4*t)
        bi = 0
        w[i] = h[i] + bi
        
        hxi = a1*exp(a2*t)*cos(a1*x[i])
        uxi = -a3*exp(a4*t)*sin(a3*x[i]) 
        bxi = 0
        
        uxxi = -a3**2*exp(a4*t)*cos(a3*x[i])
        bxxi = 0
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
    
    cxv = x[-1] + dx
    b.append(0)
    bavg= 0
    b.append(0)
    b.append(0)
    b.append(0)  
    
    return h,u,G,b,w
"""  
 
#Spatial Steps

wdir = "../../../../../data/raw/Forced/P1P2P3BedFEM/ALL2N1/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(20):
    a0 = 10.0
    a1 = 2.0
    a2 = 0.5
    a3 = 3.0
    a4 = 0.7
    a5 = 5
    
    width = 1
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (1 + sqrt(10*11))
    dt = l*dx
    startx = -width/2
    endx = width/2 + 0.9*dx
    startt = 0.0
    endt = 0.01
            
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
    bnBC = 5
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n - 1 + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - dx/6.0,x[0] + dx/6.0 , x[0] -0.5*dx])
    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xbMend = array([x[-1] + 0.5*dx, x[-1] -  dx/6.0, x[-1] + dx/6.0, x[-1] + 1.5*dx])
    
    xbhbc = []
    for i in range(len(xG)):  
        if i == 0:
            xbhbc.append(xG[i] - 0.5*dx)
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
        else:
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
            
        
    
    h,u,G,bhbc,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMend1 ,uMend1 ,GMend1 ,bta,wMend1 = ForcedbedM(xhuMend ,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    bhbc_c = copyarraytoC(bhbc)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
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
    
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:
        evolvewrapCONST(G_c,h_c,bhbc_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,hMbeg1_c,hMend1_c,wMbeg1_c,wMend1_c,GMbeg1_c,GMend1_c,uMbeg1_c,uMend1_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,ubc_c,x_c,t)
            
            
        copywritearraytoC(hMbeg1,hMbeg_c)
        copywritearraytoC(hMend1,hMend_c)
        copywritearraytoC(uMbeg1,uMbeg_c)
        copywritearraytoC(uMend1,uMend_c)
        copywritearraytoC(GMbeg1,GMbeg_c)
        copywritearraytoC(GMend1,GMend_c)
        copywritearraytoC(wMbeg1,wMbeg_c)
        copywritearraytoC(wMend1,wMend_c)
        
        hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
        hMend1 ,uMend1 ,GMend1 ,bta,wMend1 = ForcedbedM(xhuMend ,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
        
        copywritearraytoC(hMbeg1,hMbeg1_c)
        copywritearraytoC(hMend1,hMend1_c)
        copywritearraytoC(uMbeg1,uMbeg1_c)
        copywritearraytoC(uMend1,uMend1_c)
        copywritearraytoC(GMbeg1,GMbeg1_c)
        copywritearraytoC(GMend1,GMend1_c)
        copywritearraytoC(wMbeg1,wMbeg1_c)
        copywritearraytoC(wMend1,wMend1_c)
        
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, bhbc_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c)
    
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    hA,uA,GA,bhbc,wA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx)
    
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
    
    deallocPy(hMbeg1_c)
    deallocPy(GMbeg1_c)
    deallocPy(uMbeg1_c)
    deallocPy(hMend1_c)
    deallocPy(GMend1_c)
    deallocPy(uMend1_c)
    deallocPy(wMbeg1_c)
    deallocPy(wMend1_c)




"""
wdir = "../../../../../data/raw/Forced/P1P2P3BedFEM/Temp3/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

largej = 18
for j in range(largej):
    a0 = 10.0
    a1 = 2.0
    a2 = 0.5
    a3 = 3.0
    a4 = 0.7
    a5 = 5
    
    width = 1
    
    g = 9.81
    
    dx = 0.01
    l =  0.01
    dt = 0.01 / (2**j)
    startx = -width/2
    endx = width/2 + 0.9*dx
    startt = 0.0
    endt = 0.01
            
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
    bnBC = 5
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n - 1 + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - dx/6.0,x[0] + dx/6.0 , x[0] -0.5*dx])
    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xbMend = array([x[-1] + 0.5*dx, x[-1] -  dx/6.0, x[-1] + dx/6.0, x[-1] + 1.5*dx])
    
    xbhbc = []
    for i in range(len(xG)):  
        if i == 0:
            xbhbc.append(xG[i] - 0.5*dx)
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
        else:
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
            
        
    
    h,u,G,bhbc,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMend1 ,uMend1 ,GMend1 ,bta,wMend1 = ForcedbedM(xhuMend ,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    bhbc_c = copyarraytoC(bhbc)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
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
    
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:
        evolvewrapCONST(G_c,h_c,bhbc_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,hMbeg1_c,hMend1_c,wMbeg1_c,wMend1_c,GMbeg1_c,GMend1_c,uMbeg1_c,uMend1_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,ubc_c,x_c,t)
            
            
        copywritearraytoC(hMbeg1,hMbeg_c)
        copywritearraytoC(hMend1,hMend_c)
        copywritearraytoC(uMbeg1,uMbeg_c)
        copywritearraytoC(uMend1,uMend_c)
        copywritearraytoC(GMbeg1,GMbeg_c)
        copywritearraytoC(GMend1,GMend_c)
        copywritearraytoC(wMbeg1,wMbeg_c)
        copywritearraytoC(wMend1,wMend_c)
        
        hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
        hMend1 ,uMend1 ,GMend1 ,bta,wMend1 = ForcedbedM(xhuMend ,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
        
        copywritearraytoC(hMbeg1,hMbeg1_c)
        copywritearraytoC(hMend1,hMend1_c)
        copywritearraytoC(uMbeg1,uMbeg1_c)
        copywritearraytoC(uMend1,uMend1_c)
        copywritearraytoC(GMbeg1,GMbeg1_c)
        copywritearraytoC(GMend1,GMend1_c)
        copywritearraytoC(wMbeg1,wMbeg1_c)
        copywritearraytoC(wMend1,wMend1_c)
        
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, bhbc_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c)
    
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    href,uref,Gref,bhbc,wA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx)
    hC = array(hC)    
    GC = array(GC) 
    uC = array(uC) 
    

    if(j == largej):
        href = hC
        uref = uC
        Gref = GC

    
    hnorm = norm(hC - href, ord=2)/ norm(href, ord=2)
    unorm = norm(uC - uref, ord=2)/ norm(uref, ord=2)
    Gnorm = norm(GC - Gref, ord=2)/ norm(Gref, ord=2)
    
    
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dt," ",hnorm)
        file1.write(s)
    
    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dt," ",Gnorm)
        file1.write(s)   
    
    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dt," ",unorm)
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
    
    deallocPy(hMbeg1_c)
    deallocPy(GMbeg1_c)
    deallocPy(uMbeg1_c)
    deallocPy(hMend1_c)
    deallocPy(GMend1_c)
    deallocPy(uMend1_c)
    deallocPy(wMbeg1_c)
    deallocPy(wMend1_c)
"""


#FEM check
"""
wdir = "../../../../../data/raw/Forced/P1P2P3BedFEM/FEMt1N12345/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

j = 16

for j in range(30):
    
    a0 = 10.0
    a1 = 2.0
    a2 = 0.5
    a3 = 3.0
    a4 = 0.7
    a5 = 5
    
    width = 10
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.01
    dt = l*dx
    startx = -width/2
    endx = width/2 + 0.9*dx
    startt = 0.0
    endt = 0
            
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
    bnBC = 5
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n - 1 + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - dx/6.0,x[0] + dx/6.0 , x[0] -0.5*dx])
    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xbMend = array([x[-1] + 0.5*dx, x[-1] -  dx/6.0, x[-1] + dx/6.0, x[-1] + 1.5*dx])
    
    xbhbc = []
    for i in range(len(xG)):  
        if i == 0:
            xbhbc.append(xG[i] - 0.5*dx)
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
        else:
            xbhbc.append(xG[i] - dx/6.0)
            xbhbc.append(xG[i])
            xbhbc.append(xG[i] + dx/6.0)
            xbhbc.append(xG[i] + 0.5*dx)
            
        
    
    h,u,G,bhbc,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,g,dx)
    bM = cos(a5*x)
    
    
    print(t)
    hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMbeg1,uMbeg1,GMbeg1,bta,wMbeg1 = ForcedbedM(xhuMbeg,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,a0,a1,a2,a3,a4,a5,g,dx)
    hMend1 ,uMend1 ,GMend1 ,bta,wMend1 = ForcedbedM(xhuMend ,t + dt,a0,a1,a2,a3,a4,a5,g,dx)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    bhbc_c = copyarraytoC(bhbc)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
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
    
    
    
    t = 0.0
    ndt = dt
    ts.append(t)
    #Just an FEM solve here
    while t < endt:
        
            
        t = t + ndt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, bhbc_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c)
     
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
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
    
    hhbc,uta,Ghbc,bta,whbc = ForcedbedM(xhbc,t,a0,a1,a2,a3,a4,a5,g,dx)
    hta,ubc,Gt,bta,wt = ForcedbedM(xubc,t,a0,a1,a2,a3,a4,a5,g,dx)
    
    hnorm = norm(hhbc - hhbcC, ord=2)/ norm(hhbc, ord=2)
    wnorm = norm(whbc - whbcC, ord=2)/ norm(whbc, ord=2)
    Gnorm = norm(Ghbc - GhbcC, ord=2)/ norm(Ghbc, ord=2)
    unorm = norm(ubc - ubcC, ord=2)/ norm(ubc, ord=2)


    s = wdir + "w.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",wnorm)
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
    
    deallocPy(hMbeg1_c)
    deallocPy(GMbeg1_c)
    deallocPy(uMbeg1_c)
    deallocPy(hMend1_c)
    deallocPy(GMend1_c)
    deallocPy(uMend1_c)
    deallocPy(wMbeg1_c)
    deallocPy(wMend1_c)
"""