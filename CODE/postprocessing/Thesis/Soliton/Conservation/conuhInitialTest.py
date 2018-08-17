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
from Hamil import *

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


wdirord = "o3"
sdirord = "testinitialcomp"

#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

wdirb = "../../../../../data/raw/Solnon0p7NEW/"+wdirord+"/"
Mns = []
Pns = []
Mis = []
Pis = []
dxs = []
Ens = []
Eis = []
Gns = []
Gis = []

sdir = "../../../../../data/ThesisPost/Soliton/"+sdirord+"/C1/"

if not os.path.exists(sdir):
        os.makedirs(sdir)
        
        
for ki in range(6,18):
        
        
    wdir = wdirb + "/"+ str(ki) + "/"
     
    """
    s = wdir + "outlast.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         he = []
         ue = []
         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                x.append(float(row[3]))
                h.append(float(row[4]))
                u.append(float(row[6]))
                ai = float(float(row[7]))
             j = j + 1  
    """
      
    #plot(x,h)
    #plot(x,he)
    dx = 100.0 / 2**ki
    #dx = 0.1
    t0 = 0
      
    
    niBC = 4
    startx = -250
    endx = 250
    
    xbc = arange(startx - niBC*dx ,endx + niBC*dx + 0.1*dx,dx)
    
    #u0 = u[0]*ones(niBC)
    #u1 = u[-1]*ones(niBC)   
    #h0 = h[0]*ones(niBC)
    #h1 = h[-1]*ones(niBC)
    #G0 = G[0]*ones(niBC)
    #G1 = G[-1]*ones(niBC)

    hbc,ubc,Gbc = solitoninit(len(xbc),a0,a1,g,xbc,t0,dx)
    n = len(xbc ) - 2*niBC
    
    #hbc =  concatenate([h0,h,h1])
    #ubc =  concatenate([u0,u,u1])
    #Gbc =  concatenate([G0,G,G1])
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
    Gbc_c = copyarraytoC(Gbc)
    
    #hi,ui = solitoninit(n,1,1,9.81,x,0,dx)

    En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    Pn = uhall(xbc_c,hbc_c,ubc_c,n + 2*niBC,niBC,dx)
    Mn = hall(xbc_c,hbc_c,n + 2*niBC,niBC,dx)
    Gcn = Gall(xbc_c,Gbc_c,n + 2*niBC,niBC,dx)
    
    xbeg = -250 - 0.5*dx
    xend = 250 + 0.5*dx

    Mi = SolitonMass(a0,a1,k,xbeg,xend)
    Pi = SolitonMome(a0,a1,c,k,xbeg,xend)
    Gci =  SolitonG(a0,a1,c,k,xbeg,xend)
    Ei = SolitonHam(a0,a1,c,k,xbeg,xend)
    
    
    Pns.append(Pn)
    Pis.append(Pi)
    Mns.append(Mn)
    Mis.append(Mi)
    Gns.append(Gcn)
    Gis.append(Gci)
    dxs.append(dx)
    Ens.append(En)
    Eis.append(Ei)

Gns = array(Gns)
Gis = array(Gis) 
Ens = array(Ens)
Eis = array(Eis)    
Mns = array(Mns)
Pns = array(Pns)
Mis = array(Mis)
Pis = array(Pis)

relerrP = abs(Pis - Pns)
relerrM = abs(Mis - Mns)/ abs(Mis)
relerrG = abs(Gis - Gns)/ abs(Gis)
relerrE = abs(Eis - Ens)/ abs(Eis)



n= len(dxs)

s = sdir + "conh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",relerrM[i])
        file1.write(s) 

s = sdir + "conG.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",relerrG[i])
        file1.write(s)          
    
s = sdir + "conuh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",relerrP[i])
        file1.write(s)  
        
s = sdir + "conH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",relerrE[i])
        file1.write(s) 
       