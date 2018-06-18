# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan

This code calculates the dispersion error for our FDVM
"""
from scipy import *
from numpy import matrix
from scipy.linalg import eigvals
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

I = sqrt(-1)
   

def Rpo3(k,dx):
    Rp1 = (5 + 2*exp(-I*dx*k) - exp(I*dx*k) )*exp(I*dx*k)/6.0
    return Rp1

def Rmo3(k,dx):
    Rm1 = (5 - exp(-I*dx*k) +2*exp(I*dx*k) )/6.0
    return Rm1
    
def Rpo2(k,dx):
    Rp1 = exp(I*dx*k)*(1 - I*sin(dx*k)/2)
    return Rp1

def Rmo2(k,dx):
    Rm1 = 1 + I*sin(dx*k)/2
    return Rm1
    
def Rpo1(k,dx):
    Rp1 = exp(I*dx*k)
    return Rp1

def Rmo1(k,dx):
    Rm1 = 1
    return Rm1
    
def M3(k,dx):
    return (26 - 2*cos(k*dx))/24.0
    
def MA(k,dx):
    return k*dx/ (2*sin(k*dx/2.0))


def GA(H,k,M,dx):
    Gden = H + H**3/3*k**2
    Gnum = exp(I*k*dx/2.0)
    G1= Gnum*M/Gden
    return G1
    
def G2(H,k,dx):
    Gden = -H**3*(2*cos(dx*k) - 2)/(3*dx**2) + H
    Gnum = (1 + exp(I*k*dx))/2.0
    G1= Gnum/Gden
    return G1
 
def G4(H,k,M,dx):
    Gden = -H**3*(32*cos(dx*k) - 2*cos(2*dx*k) - 30)/(36*dx**2) + H
    Gnum = (-exp(-I*k*dx) + 9*exp(I*k*dx) - exp(2*I*k*dx) + 9)/16.0
    G1 = Gnum*M / Gden
    return G1  
    
def GNFEM(H,Rm,Rp,k,dx):
    GRec = dx/6.0*(Rm + Rp)
    GRHSp1 = 4*cos(k*dx/2.0) - 2*cos(k*dx) + 8
    GRHSp2 = 14 - 16*cos(k*dx/2.0) + 2*cos(k*dx)
    uInt = H*dx/30.0*GRHSp1 + H**3/(9*dx)*GRHSp2
    GFEM1 = GRec / uInt

    return GFEM1  
      
    
#The difference-shift operator (comes out of difference of fluxes)    
def D(k,dx):
    return 1 - exp(-I*k*dx)
 

#calculation of the flux contributions 
def Fnn(g,H,Rms,Rps):
    return -sqrt(H*g)*(Rps - Rms)/2.0
 
def FnG(H,Gs):
    return H*Gs
    
def FGG(g,H,Rms,Rps):
    return -sqrt(H*g)*(Rps - Rms)/2.0
 
def FGn(g,H,Rms,Rps):
    return H*g*(Rms + Rps)/2.0
    

def FnnA():
    return 0
 
def FnGA(k,H):
    return 3*I*k / (3 + H**2*k**2)
    
def FGGA():
    return 0
 
def FGnA(k,H,g):
    return I*g*k*H
    
    
    
def o1(k,g,h,dx,dt, l1A, l2A):
    from scipy import log

    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    
    G1 = G2(h,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    if (lams[0].imag > 0):
        l1 = abs(lams[0] - l1A) / abs(l1A)
        l2 = abs(lams[1] - l2A) / abs(l2A)
    else:
        l1 = abs(lams[1] - l1A) / abs(l1A)
        l2 = abs(lams[0] - l2A) / abs(l2A)
    
    return l1,l2
    

def o2(k,g,h,dx,dt, l1A, l2A):
    from scipy import log

    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    G1 = G2(h,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    if (lams[0].imag > 0):
        l1 = abs(lams[0] - l1A) / abs(l1A)
        l2 = abs(lams[1] - l2A) / abs(l2A)
    else:
        l1 = abs(lams[1] - l1A) / abs(l1A)
        l2 = abs(lams[0] - l2A) / abs(l2A)
    
    return l1,l2

def oFEM2(k,g,h,dx,dt, l1A, l2A):
    from scipy import log

    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    G1 = GNFEM(h,Rm1,Rp1,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    if (lams[0].imag > 0):
        l1 = abs(lams[0] - l1A) / abs(l1A)
        l2 = abs(lams[1] - l2A) / abs(l2A)
    else:
        l1 = abs(lams[1] - l1A) / abs(l1A)
        l2 = abs(lams[0] - l2A) / abs(l2A)
    
    return l1,l2


def o3(k,g,h,dx,dt, l1A, l2A):
    from scipy import log

    #calculate the elementary contribution factors
    M =   M3(k,dx)
    
    
    Rp1 = Rpo3(k,dx)
    Rm1 = Rmo3(k,dx)
    
    G1 = G4(h,k,M,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    
    if (lams[0].imag > 0):
        l1 = abs(lams[0] - l1A) / abs(l1A)
        l2 = abs(lams[1] - l2A) / abs(l2A)
    else:
        l1 = abs(lams[1] - l1A) / abs(l1A)
        l2 = abs(lams[0] - l2A) / abs(l2A)
        
    
    return l1,l2
    

from scipy import pi



h = 1.0
k = 0.5
g = 9.81

FnnAv = FnnA()
FnGAv =FnGA(k,h)
FGGAv =  FGGA()
FGnAv = FGnA(k,h,g)

FmatA = matrix([[FnnAv , FnGAv], [FGnAv, FGGAv]])
lams = eigvals(FmatA)
l1A = lams[0]
l2A = lams[1]

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),1,num=1000)   
n = len(dxs)

erro1l1 = zeros(n)
erro1l2 = zeros(n)

erro2l1 = zeros(n)
erro2l2 = zeros(n)

errof2l1 = zeros(n)
errof2l2 = zeros(n)

erro3l1 = zeros(n)
erro3l2 = zeros(n)

for i in range(n): 

    dx = dxs[i]
    Cr = 0.5
    l = Cr/(sqrt(g*h))
    dt = l*dx
    
    erro1l1[i],erro1l2[i] = o1(k,g,h,dx,dt,l1A,l2A)
    erro2l1[i],erro2l2[i] = o2(k,g,h,dx,dt,l1A,l2A)
    errof2l1[i],errof2l2[i] = oFEM2(k,g,h,dx,dt,l1A,l2A)
    erro3l1[i],erro3l2[i] = o3(k,g,h,dx,dt,l1A,l2A)


plot(dxs,erro1l2,'b',label="o1")
plot(dxs,erro2l2,'r',label="o2")
plot(dxs,errof2l2,'g',label="o2FEM")
plot(dxs,erro3l2,'black',label="o3")
legend()




    