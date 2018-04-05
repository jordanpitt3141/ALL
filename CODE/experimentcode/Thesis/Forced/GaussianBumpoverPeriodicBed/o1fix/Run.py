# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre1 import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os 
from numpy.linalg import norm 

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
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


    
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
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
    
#Forcing Problem    
wdir = "../../../../../../../data/raw/Forced/FDVM1nobed/GaussBedAll/RT2/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(15):
    g =9.81

    a0 = 1
    a1 = 0.2
    a2 = 1.3
    a3 = 0.4
    a4 = 1.5
    a5 = 0.6
    a6 = 0
    a7 = 0
    
    width = 50
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (a5 + sqrt(g*(a0 + a1)))
    dt = l*dx
    startx = -width/2
    endx = width/2 
    startt = 0.0
    endt = 0.1
            
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
    
    nBC = 3
    nBCs = 4
        
    idx = 1.0 / dx              
        
    h,u,G = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    #bM = cos(a5*x)
    hbeg = h[0]*ones(nBCs)
    hend = h[-1]*ones(nBCs)
    ubeg = zeros(nBCs)
    uend = zeros(nBCs)
    
    
    print(t)
    
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    hbeg_c = copyarraytoC(hbeg)
    hend_c = copyarraytoC(hend)
    ubeg_c = copyarraytoC(ubeg)
    uend_c = copyarraytoC(uend)
       
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrap(G_c,h_c,hbeg_c,hend_c, ubeg_c,uend_c,g,dx,dt,nBC,n,nBCs,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7)
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c,G_c,ubeg[-1], uend[0], hbeg[-1],hend[0],dx ,n,u_c)
    
    uC = copyarrayfromC(u_c,n)
    
    hA,uA,GA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
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

    
    deallocPy(hbeg_c)
    deallocPy(ubeg_c)
    deallocPy(hend_c)
    deallocPy(uend_c)