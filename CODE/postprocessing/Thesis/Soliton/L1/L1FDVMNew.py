# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from numpy import ones

from numpy import tanh

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


def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a
#Functions for the total Energy, h , uh, G and H.

def SolitonMass(a0,a1,k,xb,xe):    
    return a0*(xe - xb) + a1*(tanh(k*xe) - tanh(k*xb)) / k
    

def SolitonMome(a0,a1,c,k,xb,xe):    
    return a1*c*(tanh(k*xe) - tanh(k*xb)) / k
    
def SolitonG(a0,a1,c,k,xb,xe):
    return a1*c / (3*k) *( (3 + 2*a0**2*k**2*sech(k*xe)**2 + 2*a0*a1*k**2*sech(k*xe)**4)*tanh(k*xe) \
   -(3 + 2*a0**2*k**2*sech(k*xb)**2 + 2*a0*a1*k**2*sech(k*xb)**4)*tanh(k*xb) )

def Solitonghsquare(a0,a1,c,k,x):
     return g/ (12*k)*sech(k*x)**3 *(9*a0**2*k*x*cosh(k*x) + 3*a0**2*k*x*cosh(3*k*x) + 4*a1*(3*a0 + 2*a1 + (3*a0 + a1)*cosh(2*k*x))*sinh(k*x))

def Solitonhusquare(a0,a1,c,k,x):
     return sqrt(a1)*c**2*( -a0*arctanh( sqrt(a1)*tanh(k*x) / sqrt(a0 + a1) )/ sqrt(a0 + a1)  + sqrt(a1)*tanh(k*x))/k

def Solitonhcubeddusquare(a0,a1,c,k,x):
     return (2*a0**2*c**2*k*(a0 + 2*a1 + a0*cosh(2*k*x))*sech(k*x)**2) *(-3*a0*sqrt(a0 + a1)*arctanh( sqrt(a1)*tanh(k*x) / sqrt(a0 + a1)) + sqrt(a1)*(3*a0 + a1 - a1*sech(k*x)**2)*tanh(k*x)) \
      / (9*sqrt(a1)*(a0 + a1*sech(k*x)**2))
    
    
def SolitonHam(a0,a1,c,k,xb,xe):
    
    ghsqInt = Solitonghsquare(a0,a1,c,k,xe) -Solitonghsquare(a0,a1,c,k,xb)
    husqInt = Solitonhusquare(a0,a1,c,k,xe) -Solitonhusquare(a0,a1,c,k,xb)
    hcubedusq = Solitonhcubeddusquare(a0,a1,c,k,xe) - Solitonhcubeddusquare(a0,a1,c,k,xb)
    
    return 0.5*(ghsqInt + husqInt + hcubedusq)
    


def midpointtoca(h,dx):
    n = len(h)
    b = zeros(n)
    c = zeros(n)
    a = zeros(n)
    i24 = 1.0/24

    for i in range(n): 
        a[i-1] = -i24
        b[i] = 26*i24
        c[i] = -i24
        
    #i =0
    i = 0;
    b[i] = 1.0;
    c[i] = 0.0;

    #i=n-1
    i = n-1;
    a[i-1] = 0.0;
    b[i] = 1.0; 
    
    return TDMApy(a,b,c,h)
    
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
    


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* (1 - a0 / h[i])
        G[i] = 2.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**4*h[i] + h[i]*u[i] - 4.0/3*a0*a1**2*c*k**2*sech(k*(x[i] - c*t0))**4*tanh(k*(x[i] - c*t0))**2 - 4.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**2*h[i]*tanh(k*(x[i] - c*t0))**2
    
    return h,u,G

a0 = 1.0
a1 = 0.7
g = 9.81
k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
c = sqrt(g*(a0 + a1))


meth = "FDVM1"
#meth = "FDVM2"
#meth = "FDVM3"

#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

wdirb = "/home/jp/Documents/PhD/project/data/2018/raw/Thesis/SolitonAgn/" + meth+ "/" 
L1hs = []
L1us = []
L1Gs = []
dxs=  []

sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/SolitonFDVMRedo/" +meth+"/L1/"

if not os.path.exists(sdir):
        os.makedirs(sdir)
        
        
for ki in range(6,18):
        
        
    wdir = wdirb + str(ki) + "/"
     

    s = wdir + "outlast.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         G = []

         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                x.append(float(row[3]))
                h.append(float(row[4]))
                G.append(float(row[5]))
                u.append(float(row[6]))
             j = j + 1 


    t0 = 0
      
    n = len(x)
    niBC = 4
    startx = x[0]
    endx = x[-1]
        
    hA,uA,GA = solitoninit(n,a0,a1,g,x,t,dx)
    
    L1h  = norm(hA - h,ord=1)/ norm(hA,ord=1)
    L1u  = norm(uA - u,ord=1)/ norm(uA,ord=1)
    L1G  = norm(GA - G,ord=1)/ norm(GA,ord=1)
    
    L1hs.append(L1h)
    L1us.append(L1u)
    L1Gs.append(L1G)
    dxs.append(dx)
    



n= len(dxs)

s = sdir + "L1h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L1hs[i])
        file1.write(s) 

s = sdir + "L1u.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L1us[i])
        file1.write(s)          
    
s = sdir + "L1G.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L1Gs[i])
        file1.write(s)  

     
