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
             
    return h,u,G
    

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

def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        phi = x[i] - c*t0;
        bx[i] = bot
        h[i] = a0 + a1*(2./(exp(k*phi) + exp(-k*phi)))*(2./(exp(k*phi) + exp(-k*phi)))
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,u,G,bx 


    
def ForcedbedM(x,t,c,a0,a1,a2,a3,b1,b2,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - c*t  
        
        
        h[i] = a0 + a1*exp(-(phi - a2)**2/(2*a3))
        u[i] =  b1*exp(-(phi - a2)**2/(2*a3))
        b[i] = sin(b2*x[i])
        w[i] = h[i] + b[i]
        
        uxi = -b1*(-2*a2 - 2*c*t + 2*x[i])*exp(-(-a2 - c*t + x[i])**2/(2*a3))/(2*a3)
        hxi = -a1*(-2*a2 - 2*c*t + 2*x[i])*exp(-(-a2 - c*t + x[i])**2/(2*a3))/(2*a3)
        uxxi = b1*(-1 + (a2 + c*t - x[i])**2/a3)*exp(-(a2 + c*t - x[i])**2/(2*a3))/a3
        
        bxi = b2*cos(b2*x[i]) 
        bxxi = -b2**2*sin(b2*x[i])
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w
    

#Forcing Problem    
wdir = "../../../../../../data/raw/Forced/P1P2P3BedFEM/GaussBedO/RT12/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(10,11):
    g =9.81
    c = 2.0
    a0 = 0.0
    a1 = 1.0
    a2 = 2.0
    a3 = 3.0
    b1 = 0.5
    b2 = 0.1
    
    width = 50
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (1 + sqrt(10*11))
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
    bnBC = 5
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - dx/6.0,x[0] + dx/6.0 , x[0] -0.5*dx])
    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xbMend = array([x[-1] + 0.5*dx, x[-1] -  dx/6.0, x[-1] + dx/6.0, x[-1] + 1.5*dx])
                
        
    
    h,u,G,b,w = ForcedbedM(x,t,c,a0,a1,a2,a3,b1,b2,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,c,a0,a1,a2,a3,b1,b2,g,dx)
    hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,c,a0,a1,a2,a3,b1,b2,g,dx)
    bMbeg = sin(b2*xbMbeg)
    bMend = sin(b2*xbMend)
    
    
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
    
    hA,uA,GA,bhbc,wA = ForcedbedM(x,t,c,a0,a1,a2,a3,b1,b2,g,dx)
    
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