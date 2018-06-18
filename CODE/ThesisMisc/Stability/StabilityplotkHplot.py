# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from scipy.linalg import norm, eig
from numpy import roots
import os

I = sqrt(-1)


#given a matrix [[a,b],[c,d]] calculates the eigenvalues of it
def eigenvaluecalc2by2(a,b,c,d):
    T = a + d
    D = a*d - b*c
    
    l1 = T/2.0 + sqrt((T*T)/4.0 - D)
    l2 = T/2.0 - sqrt((T*T)/4.0 - D)
    
    return l1,l2
 
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
    return 1 - e**(-I*k*dx)
 #calculation of the flux contributions  
def Fnn(U,g,H,Gn,Rms,Rps):
    if(U > sqrt(g*H)):
        return H*Gn + U*Rms
    else:
        return H*Gn + U/2.0*(Rms + Rps) -sqrt(H*g)*(Rps - Rms)/2.0
 
def FnG(U,g,H,Gs):
    if(U > sqrt(g*H)):
        return H*Gs
    else:
        return H*Gs
    
def FGG(U,g,H,GG,Rms,Rps):
    if(U > sqrt(g*H)):
        return U*Rms + U*H*GG
    else:
        return U*H*GG + U/2.0*(Rms + Rps)  -sqrt(H*g)*(Rps - Rms)/2.0
 
def FGn(U,g,H,Gn,Rms,Rps):
    if(U > sqrt(g*H)):
        return U*H*Gn + g*H*Rms
    else:
        return U*H*Gn + g*H/2.0*(Rms + Rps)  +U*sqrt(H*g)*(Rms - Rps)/2.0


def SWWEFnn(U,g,H,Rm1,Rp1):
    if(U > sqrt(g*H)):
        return U*Rm1
    else:
        return U*(Rm1 + Rp1)/2.0 - sqrt(g*H)/2.0*(Rp1 - Rm1)

def SWWEFnv(U,g,H,Rm1,Rp1):
    if(U > sqrt(g*H)):
        return H*Rm1
    else:
        return U*H*(Rm1 - Rp1)/(2.0*sqrt(g*H))  + H/2.0*(Rp1 + Rm1)    
 
def SWWEFvn(U,g,H,Rm1,Rp1):
    if(U > sqrt(g*H)):
        return g*Rm1
    else:
        return U*g*(Rm1 - Rp1) /(2.0*sqrt(g*H))  + g/2.0*(Rp1 + Rm1)  

def SWWEFvv(U,g,H,Rm1,Rp1):
    if(U > sqrt(g*H)):
        return U*Rm1
    else:
        return U*(Rm1 + Rp1)/2.0 - sqrt(g*H)/2.0*(Rp1 - Rm1)
   
def SWWEo2(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    
    
    Fnn1 = SWWEFnn(U,g,h,Rm1,Rp1)
    Fnv1 = SWWEFnv(U,g,h,Rm1,Rp1)
    Fvn1 = SWWEFvn(U,g,h,Rm1,Rp1)
    Fvv1 = SWWEFvv(U,g,h,Rm1,Rp1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fnn1 = D1/(dx)*Fnn1
    Fnv1 = D1/(dx)*Fnv1
    Fvn1 = D1/(dx)*Fvn1
    Fvv1 = D1/(dx)*Fvv1
    
    F = matrix([[Fnn1,Fnv1], [Fvn1,Fvv1]])
    
    M = eye(2) -dt*F + (dt*F)**2/2.0 
    
    lam , v = eig(M)
    return max(abs(lam))   
    
def o1(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    
    GG,Gn = G2(U,h,k,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(U,g,h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    F = matrix([[Fnn1,FnG1], [FGn1,FGG1]])
    
    M = eye(2) -dt*F
    
    lam , v = eig(M)
    return max(abs(lam))


def o2(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    GG,Gn = G2(U,h,k,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(U,g,h,GG)
    
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1

    
    F = matrix([[Fnn1,FnG1], [FGn1,FGG1]])
    
    M = 0.5*(2*eye(2) - 2*dt*F + dt*dt*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))

def o2FEM(U,k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    GG,Gn = GNFEM(U,h,Rm1,Rp1,k,dx)
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(U,g,h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    F = matrix([[Fnn1,FnG1], [FGn1,FGG1]])
    
    M = 0.5*(2*eye(2) - 2*dt*F + dt*dt*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))

def o3(U,k,g,h,dx,dt):

    #calculate the elementary contribution factors
    M =   M3(k,dx)
    
    
    Rp1 = Rpo3(k,dx)
    Rm1 = Rmo3(k,dx)
    
    GG,Gn = G4(U,h,k,M,dx)
    
    

    Fgn1 = FGn(U,g,h,Gn,Rm1,Rp1)
    Fgg1 = FGG(U,g,h,GG,Rm1,Rp1)
    Fnn1 = Fnn(U,g,h,Gn,Rm1,Rp1)
    Fng1 = FnG(U,g,h,GG)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    F = matrix([[Fnn1,FnG1], [FGn1,FGG1]])
    
    M = (eye(2) - dt*F + 0.5*dt*dt*F*F - (1.0/6.0)*dt*dt*dt*F*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))
    
def naiveFD(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    a = -2*im*dt*idx*u*sin(k*dx)
    b = -2*im*dt*idx*h*sin(k*dx)
    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    F = matrix([[a,b,1,0],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    lam , v = eig(F)
    return max(abs(lam))
    
    
def LWN(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    
    a = 1 -dt*idx*( 0.5*h*c*im*sin(k*dx) - u*im*sin(k*dx) -u*u*dt*idx*(cos(k*dx) - 1))
    b = -dt*idx*( 0.5*h*im*sin(k*dx)*(d+1) -u*h*dt*idx*(cos(k*dx) - 1))
    
    e = -dt*idx*h*0.5*im*sin(k*dx)
    
    F = matrix([[a,b,0,e],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    lam , v = eig(F)
    return max(abs(lam))


tol = 10**(-14)
tol2 = 10**-6

Cr = 0.5

dxs = [0.001,0.01,0.1,1]
us = [100,1000,10000]
Hs = linspace(tol2,5,num=50)  
ks = linspace(tol2,5,num=50)  
g = 9.81
for dx in dxs :
    
    sdir = "../../../../../data/ThesisRAW/StabilityFixed/Cr0p5PicskH1/"   + "dx"+str(dx) +"/"

    if not os.path.exists(sdir):
        os.makedirs(sdir)
        
    for u in us:
        
        
        mks = []
        mhs = []
        pks = []
        phs = []
        qks = []
        qhs = []
        rks = []
        rhs = []
        nks = []
        nhs = []
        sks = []
        shs = []
        
        mksn = []
        mhsn = []
        pksn = []
        phsn = []
        qksn = []
        qhsn = []
        rksn = []
        rhsn = []
        nksn = []
        nhsn = []
        sksn = []
        shsn = []
        
        for h in Hs:
            for k in ks:
                l= Cr/ (u + sqrt(g*h))
                dt = l*dx
                
                m = naiveFD(k,g,u,h,dx,dt)
                q = LWN(k,g,u,h,dx,dt)
                p = o1(u,k,g,h,dx,dt)
                r = o2(u,k,g,h,dx,dt)
                n = o3(u,k,g,h,dx,dt)
                s = o2FEM(u,k,g,h,dx,dt)
                
                if m <= 1 + tol:
                    mks.append(k)
                    mhs.append(h)
                else:
                    mksn.append(k)
                    mhsn.append(h)                
    
                if q <= 1 + tol:
                    qks.append(k)
                    qhs.append(h)
                else:
                    qksn.append(k)
                    qhsn.append(h)  
                    
                if r <= 1 + tol:
                    rks.append(k)
                    rhs.append(h) 
                else:
                    rksn.append(k)
                    rhsn.append(h)  
                    
                if p <= 1 + tol:
                    pks.append(k)
                    phs.append(h)
                else:
                    pksn.append(k)
                    phsn.append(h)  
                    
                if n <= 1 + tol:
                    nks.append(k)
                    nhs.append(h)
                else:
                    nksn.append(k)
                    nhsn.append(h)  
                    
                if s <= 1 + tol:
                    sks.append(k)
                    shs.append(h)
                else:
                    sksn.append(k)
                    shsn.append(h)  
    
        xlabel("k")
        ylabel("H")
        plot(pks,phs,'.',label="stable")
        plot(pksn,phsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "o1U" + str(u) + ".png"
        savefig(s, bbox_inches='tight')
        clf()
        
        xlabel("k")
        ylabel("H")
        plot(rks,rhs,'.',label="stable")
        plot(rksn,rhsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "o2U" + str(u) + ".png"
        savefig(s, bbox_inches='tight')
        clf()
        
        xlabel("k")
        ylabel("H")
        plot(sks,shs,'.',label="stable")
        plot(sksn,shsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "o2FEMU" + str(u) + ".png"
        savefig(s, bbox_inches='tight')
        clf()    
        
        xlabel("k")
        ylabel("H")
        plot(nks,nhs,'.',label="stable")
        plot(nksn,nhsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "o3U" + str(u) + ".png"
        savefig(s, bbox_inches='tight')
        clf() 
    
        xlabel("k")
        ylabel("H")
        plot(mks,mhs,'.',label="stable")
        plot(mksn,mhsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "nFDU" + str(u)  + ".png"
        savefig(s, bbox_inches='tight')
        clf() 
        
        xlabel("k")
        ylabel("H")
        plot(qks,qhs,'.',label="stable")
        plot(qksn,qhsn,'.r',label="unstable")
        legend(loc=2)
        s = sdir + "LWU" + str(u) + ".png"
        savefig(s, bbox_inches='tight')
        clf() 



    
    