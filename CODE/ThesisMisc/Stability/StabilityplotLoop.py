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
    
    F = matrix([[Fnn1,FnG1], [FGn1,FGG1]])
    
    M = eye(2) - dt*F
    
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
    Fng1 = FnG(h,GG)
    
    
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
    Fng1 = FnG(h,GG)
    
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
    Fng1 = FnG(h,GG)
    
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
    

def wactual(k,g,h0,u0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return u0*k+ w1,u0*k-w1


sdir = "../../../../../data/ThesisRAW/StabilityCr0p5PicsHyp1t2/"  

if not os.path.exists(sdir):
    os.makedirs(sdir)

us = [100,10,1,0.1,0.01,0.001,0]
Hs = [100,10,1,0.1,0.01,0.001,0.0001]
ks = [100,10,1,0.1,0.01]
g = 9.81

for u in us:
    #for h in Hs:
        h = sqrt(u**2/g)*2
        for k in ks:
        
            
            w1,w2 = wactual(k,g,h,u)
            
            #calculate the phase speeds by the relation omega/k
            v1 = w1/k
            v2 = w2/k
            
            wdir = "../../../../data/ThesisRAW/Stability/h" + str(h) + "/k"+ str(k) + "/"
            
            if not os.path.exists(wdir):
                os.makedirs(wdir)
            
            
            
            #this then measures the dispersion relation from the highest resolution soliton problem until 0.5
            dxs = linspace(10**(-15),1.0/k,num=1000)  
            n1 = len(dxs)
            
            Ss= []
            Ms= []
            Ns = []
            Ps = []
            Qs = []
            Rs = []
            dts = []
            for i in range(n1):
                dx = dxs[i]
                Cr = 0.5
                l = Cr/ (u + sqrt(g*h))
                dt = l*dx
                
                m = naiveFD(k,g,u,h,dx,dt)
                q = LWN(k,g,u,h,dx,dt)
                p = o1(u,k,g,h,dx,dt)
                r = o2(u,k,g,h,dx,dt)
                n = o3(u,k,g,h,dx,dt)
                s = o2FEM(u,k,g,h,dx,dt)
            
                Ms.append(m)
                Qs.append(q)
                Ns.append(n) 
                Ps.append(p)
                Rs.append(r)
                Ss.append(s)
                dts.append(dt)
            
            """
                s = wdir + "Naive.dat"
                with open(s,'a') as file1:
                    s ="%3.8f%5s%1.15f\n" %(k*dx," ",m)
                    file1.write(s) 
             
                s = wdir + "LW.dat"
                with open(s,'a') as file1:
                    s ="%3.8f%5s%1.15f\n" %(k*dx," ",q)
                    file1.write(s) 
            
            NormMS = norm(Ms,ord=1) / n1
            NormQS = norm(Qs,ord=1) / n1
            
            s = wdir + "Norms.dat"
            with open(s,'w') as file1:
                s ="u | h | k | Norm Naive \n"
                file1.write(s) 
                s ="%3.8f%5s%3.8f%5s%3.8f%5s%1.20f\n \n" %(u," ",h," ",k," ",NormMS)
                file1.write(s) 
                s ="u | h | k | Norm LW \n"
                file1.write(s) 
                s ="%3.8f%5s%3.8f%5s%3.8f%5s%1.20f\n" %(u," ",h," ",k," ",NormQS)
                file1.write(s) 
            """
            print('max Qs:')
            print(max(Qs))
            print('max Ms:')
            print(max(Ms))
            print()
            print()
            
            
            xlabel("dx")
            ylabel("Growth Factor")
            plot(k*dxs,Ps,label="o1")
            plot(k*dxs,Rs,label="o2")
            plot(k*dxs,Ss,label="o2FEM")
            plot(k*dxs,Ns,label="o3")
            plot(k*dxs,Ms,label="Naive FD")
            plot(k*dxs,Qs,label="Lax Wendroff")
            ylim([0.8,1.1])  
            legend(loc='lower left')
            s = sdir + "U" + str(u) +  "H" + str(h) +  "K" + str(k) + ".png"
            savefig(s, bbox_inches='tight')
            clf()

    
    