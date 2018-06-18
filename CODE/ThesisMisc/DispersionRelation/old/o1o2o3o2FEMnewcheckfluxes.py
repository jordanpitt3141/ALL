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
    
    
    
def o1(k,g,h,dx,dt, FGNA, FGGA, FNNA,FNGA):
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
    
    relerrGn = abs(FGn1 - FGNA)/ abs(FGNA)
    relerrGG = abs(FGG1 - FGGA)
    relerrnn = abs(Fnn1 - FNNA)
    relerrnG = abs(FnG1 - FNGA)/ abs(FNGA)
    
    return relerrnn,relerrnG,relerrGn,relerrGG 
    

def o2(k,g,h,dx,dt, FGNA, FGGA, FNNA,FNGA):
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
    
    relerrGn = abs(FGn1 - FGNA)/ abs(FGNA)
    relerrGG = abs(FGG1 - FGGA)
    relerrnn = abs(Fnn1 - FNNA)
    relerrnG = abs(FnG1 - FNGA)/ abs(FNGA)
    
    return relerrnn,relerrnG,relerrGn,relerrGG

def oFEM2(k,g,h,dx,dt, FGNA, FGGA, FNNA,FNGA):
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
    
    relerrGn = abs(FGn1 - FGNA)/ abs(FGNA)
    relerrGG = abs(FGG1 - FGGA)
    relerrnn = abs(Fnn1 - FNNA)
    relerrnG = abs(FnG1 - FNGA)/ abs(FNGA)
    
    return relerrnn,relerrnG,relerrGn,relerrGG


def o3(k,g,h,dx,dt, FGNA, FGGA, FNNA,FNGA):
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
    
    relerrGn = abs(FGn1 - FGNA)/ abs(FGNA)
    relerrGG = abs(FGG1 - FGGA)
    relerrnn = abs(Fnn1 - FNNA)
    relerrnG = abs(FnG1 - FNGA)/ abs(FNGA)
    
    return relerrnn,relerrnG,relerrGn,relerrGG
    

from scipy import pi



h = 1.0
k = 0.5
g = 9.81

FnnAv = FnnA()
FnGAv =FnGA(k,h)
FGGAv =  FGGA()
FGnAv = FGnA(k,h,g)

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),1,num=1000)   
n = len(dxs)

erro1nn = zeros(n)
erro1nG = zeros(n)
erro1Gn = zeros(n)
erro1GG = zeros(n)

erro2nn = zeros(n)
erro2nG = zeros(n)
erro2Gn = zeros(n)
erro2GG = zeros(n)

errof2nn = zeros(n)
errof2nG = zeros(n)
errof2Gn = zeros(n)
errof2GG = zeros(n)

erro3nn = zeros(n)
erro3nG = zeros(n)
erro3Gn = zeros(n)
erro3GG = zeros(n)


for i in range(n): 

    dx = dxs[i]
    Cr = 0.5
    l = Cr/(sqrt(g*h))
    dt = l*dx
    
    erro1nn[i],erro1nG[i],erro1Gn[i],erro1GG[i] = o1(k,g,h,dx,dt, FGnAv, FGGAv, FnnAv,FnGAv)
    erro2nn[i],erro2nG[i],erro2Gn[i],erro2GG[i] = o2(k,g,h,dx,dt, FGnAv, FGGAv, FnnAv,FnGAv)
    errof2nn[i],errof2nG[i],errof2Gn[i],errof2GG[i] = oFEM2(k,g,h,dx,dt, FGnAv, FGGAv, FnnAv,FnGAv)
    erro3nn[i],erro3nG[i],erro3Gn[i],erro3GG[i] = o3(k,g,h,dx,dt, FGnAv, FGGAv, FnnAv,FnGAv)

"""
plot(dxs,erro1nn,'b',label="o1")
plot(dxs,erro2nn,'r',label="o2")
plot(dxs,errof2nn,'g',label="o2FEM")
plot(dxs,erro3nn,'black',label="o3")
legend()

plot(dxs,erro1nG,'b',label="o1")
plot(dxs,erro2nG,'r',label="o2")
plot(dxs,errof2nG,'g',label="o2FEM")
plot(dxs,erro3nG,'black',label="o3")
legend()

plot(dxs,erro1GG,'b',label="o1")
plot(dxs,erro2GG,'r',label="o2")
plot(dxs,errof2GG,'g',label="o2FEM")
plot(dxs,erro3GG,'black',label="o3")
legend()
"""
plot(dxs,erro1Gn,'b',label="o1")
plot(dxs,erro2Gn,'r',label="o2")
plot(dxs,errof2Gn,'g',label="o2FEM")
plot(dxs,erro3Gn,'black',label="o3")
legend()




    