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
    Rp1 = exp(I*dx*k)*(1 - I*sin(dx*k)/2.0)
    return Rp1

def Rmo2(k,dx):
    Rm1 = 1 + I*sin(dx*k)/2.0
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
    
def G2(U,H,k,dx):
    Gden = -H**3*(2*cos(dx*k) - 2)/(3*dx**2) + H
    Gnum = (1 + exp(I*k*dx))/2.0
    GG = Gnum/Gden
    Gn = -U*GG 
    return GG,Gn
 
def G4(U,H,k,M,dx):
    Gden = -H**3*(32*cos(dx*k) - 2*cos(2*dx*k) - 30)/(36*dx**2) + H
    Gnum = (-exp(-I*k*dx) + 9*exp(I*k*dx) - exp(2*I*k*dx) + 9)/16.0
    GG = Gnum*M / Gden
    Gn = -U*GG
    return GG,Gn
    
def GNFEM(U,H,Rm,Rp,k,dx):
    GRec = dx*(Rm + Rp)/6.0
    GRHSp1 = 4*cos(k*dx/2.0) - 2*cos(k*dx) + 8 
    GRHSp2 = 14 - 16*cos(k*dx/2.0) + 2*cos(k*dx)
    uInt = H*dx/30.0*GRHSp1 + H**3/(9*dx)*GRHSp2
    GG = GRec / uInt
    Gn = -U*GG

    return GG,Gn 
      
    
#The difference-shift operator (comes out of difference of fluxes)    
def D(k,dx):
    return 1 - exp(-I*k*dx)
 

#calculation of the flux contributions 
def Fnn(U,g,H,Gn,Rms,Rps):
    return H*Gn + U/2.0*(Rms + Rps) -sqrt(H*g)*(Rps - Rms)/2.0
 
def FnG(H,Gs):
    return H*Gs
    
def FGG(U,g,H,GG,Rms,Rps):
    return U*H*GG + U/2.0*(Rms + Rps)  -sqrt(H*g)*(Rps - Rms)/2.0
 
def FGn(U,g,H,Gn,Rms,Rps):
    return U*H*Gn + g*H/2.0*(Rms + Rps)  +U*sqrt(H*g)*(Rms - Rps)/2.0
    
    
    
def o1(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    
    GG,Gn = G2(U,h,k,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2)
    
    return wn1,wn2    
    

def o2(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    GG,Gn = G2(U,h,k,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(h,GG)
    
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2 + 0.5*dt*dt*l2*l2)
    
    return wn1,wn2

def oFEM2(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    GG,Gn = GNFEM(U,h,Rm1,Rp1,k,dx)
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2 + 0.5*dt*dt*l2*l2)
    
    return wn1,wn2


def o3(U,k,g,h,dx,dt):

    #calculate the elementary contribution factors
    M =   M3(k,dx)
    
    
    Rp1 = Rpo3(k,dx)
    Rm1 = Rmo3(k,dx)
    
    GG,Gn = G4(U,h,k,M,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1 - dt*dt*dt*l1*l1*l1/6.0)
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2  - dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2
    
    

    
def wactual(U,k,g,h0):
    w1 =U*k + k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    w2 =U*k - k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,w2

sdir = "../../../../../data/ThesisRAW/FullDispRelCr1o2PicsHyp1t1/"    

if not os.path.exists(sdir):
    os.makedirs(sdir)

us = [1000,100,10,1,0.1,0.01,0.001,0]
Hs = [100,10,1,0.1,0.01,0.001,0.0001]
ks = [100,10,1,0.1,0.01]

for u in us:
    for h in Hs:
            #    for k in ks:


    
            #Dispersion Error 
            g = 1
            k = 1
            
            wdir = "../../../../../data/ThesisRAW/FullDispRelCr1o2/U" + str(u) + "/h" + str(h) + "/k"+ str(k) + "/"
            
            if not os.path.exists(wdir):
                os.makedirs(wdir)
                
            
            w1,w2 = wactual(u,k,g,h)
            
            
            #this then measures the dispersion relation from the highest resolution soliton problem until 0.5
            dxs = linspace(10**(-5),1.0/k,num=1000)   
            n = len(dxs)
            w1s = zeros(n)
            w2s = zeros(n) 
            
            o21s = zeros(n,dtype=complex)
            o22s = zeros(n,dtype=complex)
            of21s = zeros(n,dtype=complex)
            of22s = zeros(n,dtype=complex)
            o11s = zeros(n,dtype=complex)
            o12s = zeros(n,dtype=complex)
            o31s = zeros(n,dtype=complex)
            o32s = zeros(n,dtype=complex)
            of31s = zeros(n,dtype=complex)
            of32s = zeros(n,dtype=complex)
            
            erro21 = zeros(n)
            erro22 = zeros(n)
            errof21 = zeros(n)
            errof22 = zeros(n)
            erro11 = zeros(n)
            erro12 = zeros(n)
            erro31 = zeros(n)
            erro32 = zeros(n)
            errof31 = zeros(n)
            errof32 = zeros(n)
            
            erro1A = zeros(n)
            erro2A = zeros(n)
            erro2fA = zeros(n)
            erro3A = zeros(n)
            
            for i in range(n):  
            
                v1 = w1/k
                v2 = w2/k
                dx = dxs[i]
                Cr = 0.5
                l = Cr/ (abs(u) + sqrt(g*h))
                dt = l*dx
                        
                o11,o12 = o1(u,k,g,h,dx,dt)
                o21,o22 = o2(u,k,g,h,dx,dt)
                of21,of22 = oFEM2(u,k,g,h,dx,dt)
                o31,o32 = o3(u,k,g,h,dx,dt)
            
            
                #Somehow this process does not guarantee that o11 and o12 are either always w1 or w2
                #so we check that o11 is positive or not, and then assign o11s and o12s accordingly so that
                # o11s[i] corresponds to the positive analytic dispersion w1
            
                if u == 0 :
                    if(o11.real > 0):
                        o11s[i] = o11
                        o12s[i] = o12
                    else:
                        o12s[i] = o11
                        o11s[i] = o12
                        
                    if(o21.real > 0):
                        o21s[i] = o21
                        o22s[i] = o22
                    else:
                        o22s[i] = o21
                        o21s[i] = o22
                        
                    if(of21.real > 0):
                        of21s[i] = of21
                        of22s[i] = of22
                    else:
                        of22s[i] = of21
                        of21s[i] = of22
                 
                    if(o31.real > 0):
                        o31s[i] = o31
                        o32s[i] = o32
                    else:
                        o32s[i] = o31
                        o31s[i] = o32
                
                else:
                    o11s[i] = o11
                    o12s[i] = o12
                    
                    o21s[i] = o21
                    o22s[i] = o22
                
                    
                    of21s[i] = of21
                    of22s[i] = of22
                
                    o31s[i] = o31
                    o32s[i] = o32
            
            
                
                erro11[i] = abs(w1+ o12s[i])  /abs(w1)
                #erro12[i] = abs(w2 - o12s[i]) / abs(w2)
                
                erro21[i] = abs(w1+ o22s[i])  /abs(w1)
                #erro22[i] = abs(w2 - o22s[i]) / abs(w2)
                
                errof21[i] = abs(w1+ of22s[i])  /abs(w1)
                #errof22[i] = abs(w2 - of22s[i]) / abs(w2)
                
                erro31[i] = abs(w1 + o32s[i])  /abs(w1)
                #erro32[i] = abs(w2 - o32s[i]) / abs(w2)

            
            plot(k*dxs,erro11,label="o1 N")
            #plot(k*dxs,erro1A,label="o1 A")
            plot(k*dxs,erro21,label="o2 N")
            #plot(k*dxs,erro2A,label="o2 A")
            plot(k*dxs,errof21,label="o2FEM N")
            #plot(k*dxs,erro2fA,label="o2FEM A")
            plot(k*dxs,erro31,label="o3 N")
            #plot(k*dxs,erro3A,label="o3 A")
            legend(loc=2)
            s = sdir + "U" + str(u) +  "H" + str(h) +  "K" + str(k) + ".png"
            savefig(s, bbox_inches='tight')
            clf()


