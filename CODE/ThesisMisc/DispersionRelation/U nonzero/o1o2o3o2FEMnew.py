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
    
    #calculate eigenvalues
    lams = eigvals(F)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2 + 0.5*dt*dt*l2*l2)
    
    return wn1,wn2   
    
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
    Fng1 = FnG(U,g,h,GG)
    
    
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
    Fng1 = FnG(U,g,h,GG)
    
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
    Fng1 = FnG(U,g,h,GG)
    
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
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2 - dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2
    
def naiveFD(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    a = -2*im*dt*idx*u*sin(k*dx)
    b = -2*im*dt*idx*h*sin(k*dx)
    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    F = matrix([[a,b,1,0],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    #calculate eigenvalues
    lams = eigvals(F)
    l1 = lams[0]
    l2 = lams[1]
    l3 = lams[2]
    l4 = lams[3]

    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(l1)
    wn2 = 1.0/(I*dt)*log(l2)
    wn3 = 1.0/(I*dt)*log(l3)
    wn4 = 1.0/(I*dt)*log(l4)
    
    return wn1,wn2,wn3,wn4
    
    
def LWN(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    
    a = 0.5*(2 + dt**2*( - 4*u**2*sin(dx*k/2.0)**2 / (dx**2) + 3*(exp(2*im*dx*k) -1)*g*h/ (6*dx*dx + 8*h**2*sin(dx*k/2.0)**2)) \
     - 2*im*dt*u*sin(dx*k)/dx - 3*im*dt*dt*exp(-im*dx*k)*g*h*sin(k*dx) / (3*dx*dx + 4*h*h*sin(dx*k/2.0)**2))
     
    b = dt/(2*dx*dx)*h*(dt*u*(cos(2*dx*k) + 2*cos(k*dx) - 3) - im*dx*sin(dx*k))
    
    e = -dt*idx*h*0.5*im*sin(k*dx)
    
    F = matrix([[a,b,0,e],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    #calculate eigenvalues
    lams = eigvals(F)
    l1 = lams[0]
    l2 = lams[1]
    l3 = lams[2]
    l4 = lams[3]

    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(l1)
    wn2 = 1.0/(I*dt)*log(l2)
    wn3 = 1.0/(I*dt)*log(l3)
    wn4 = 1.0/(I*dt)*log(l4)
    
    return wn1,wn2,wn3,wn4
    

def close(alist,anumber):
    
    nlen = len(alist)
    si = 0

    for i in range(nlen):
        if not isnan(alist[i]):
            si = i
            result = alist[i]
            break


    for i in range(si,nlen) :
        if(abs(anumber - alist[i]) < abs(anumber - result)):
            result = alist[i]
    return result

    
def wactual(U,k,g,h0):
    w1 =U*k + k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    w2 =U*k - k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,w2

def SWWEwactual(U,k,g,h0):
    w1 =U*k + k*sqrt(g*h0)
    w2 =U*k - k*sqrt(g*h0)
    return w1,w2

#Dispersion Error 
u = 0.0
h = 1
k = 2.5
k = pi/10.0
wavel = (2.0*pi)/k
g = 9.81

wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/FullDispRel/U" + str(u) + "/h" + str(h) + "/k"+ str(k) + "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    

w1,w2 = wactual(u,k,g,h)
Sw1,Sw2 = SWWEwactual(u,k,g,h)


#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-15),0.2*wavel,num=500)   
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 

oD21s = zeros(n,dtype=complex)
oD22s = zeros(n,dtype=complex)
oW21s = zeros(n,dtype=complex)
oW22s = zeros(n,dtype=complex)
o21s = zeros(n,dtype=complex)
o22s = zeros(n,dtype=complex)
So21s = zeros(n,dtype=complex)
So22s = zeros(n,dtype=complex)
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
erroD21 = zeros(n)
erroD22 = zeros(n)
erroW21 = zeros(n)
erroW22 = zeros(n)
errSo21 = zeros(n)
errSo22 = zeros(n)
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
            
    o1v = o1(u,k,g,h,dx,dt)
    o2v= o2(u,k,g,h,dx,dt)
    oD2v = naiveFD(k,g,u,h,dx,dt)
    oW2v = LWN(k,g,u,h,dx,dt)
    of2v = oFEM2(u,k,g,h,dx,dt)
    o3v = o3(u,k,g,h,dx,dt)
    So2v = SWWEo2(u,k,g,h,dx,dt)


    #Somehow this process does not guarantee that o11 and o12 are either always w1 or w2
    #so we check that o11 is positive or not, and then assign o11s and o12s accordingly so that
    # o11s[i] corresponds to the positive analytic dispersion w1

    o11s[i] = close(o1v,-w1)
    o12s[i] = close(o1v,-w2)
    
    o21s[i] = close(o2v,-w1)
    o22s[i] = close(o2v,-w2)
    
    of21s[i] = close(of2v,-w1)
    of22s[i] = close(of2v,-w2)
    
    oD21s[i] = close(oD2v,-w1)
    oD22s[i] = close(oD2v,-w2)
    
    oW21s[i] = close(concatenate((oW2v,-array(oW2v))),-w1)
    oW22s[i] = close(concatenate((oW2v,-array(oW2v))),-w2)
    
    o31s[i] = close(o3v,-w1)
    o32s[i] = close(o3v,-w2)   

    
    erro11[i] = abs(-w1 - o11s[i])  /abs(w1)   
    erro21[i] = abs(-w1  - o21s[i])  /abs(w1)
    errof21[i] = abs(-w1 - of21s[i])  /abs(w1)    
    erro31[i] = abs(-w1 - o31s[i])  /abs(w1)
    
    erroD21[i] = abs(-w1  - oD21s[i])  /abs(w1)
    
    erroW21[i] = abs(-w1  - oW21s[i])  /abs(w1)
    
#
#    s = wdir + "o1.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",erro11[i])
#        file1.write(s)   
#        
#    s = wdir + "o2.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",erro21[i])
#        file1.write(s) 
#        
#    s = wdir + "oD2.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",erroD21[i])
#        file1.write(s) 
#        
#    s = wdir + "oW2.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",erroW21[i])
#        file1.write(s) 
#        
#    s = wdir + "o2FEM.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",errof21[i])
#        file1.write(s) 
#
#    s = wdir + "o3.dat"
#    with open(s,'a') as file1:
#        s ="%3.8f%5s%1.15f\n" %(dxs[i]/wavel," ",erro31[i])
#        file1.write(s)  

  



plot(dxs/wavel,erro11,label="o1 N")
plot(dxs/wavel,erro21,label="o2 N")
plot(dxs/wavel,errof21,label="o2FEM N")
plot(dxs/wavel,erro31,label="o3 N")
plot(dxs/wavel,erroD21,label="oD2 N")
plot(dxs/wavel,erroW21,label="oW2 N")
legend(loc=2)

"""
plot(k*dxs,erro11,label="o1 N")
#plot(k*dxs,erro1A,label="o1 A")
plot(k*dxs,erro21,label="o2 N")
plot(k*dxs,erroD21,label="oD2 N")
plot(k*dxs,erroW21,label="oW2 N")
#plot(k*dxs,errSo21,label="So2 N")
#plot(k*dxs,erro2A,label="o2 A")
plot(k*dxs,errof21,label="o2FEM N")
#plot(k*dxs,erro2fA,label="o2FEM A")
plot(k*dxs,erro31,label="o3 N")
#plot(k*dxs,erro3A,label="o3 A")
legend(loc=2)
"""



## G Error   
"""
h = 1.0
k = 0.5
g = 9.81

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),10,num=1000)   
n = len(dxs)

erro1 = zeros(n)
erro2 = zeros(n)
errof2 = zeros(n)
erro3 = zeros(n)

for i in range(n): 

    dx = dxs[i]
    v1 = w1/k
    v2 = w2/k
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    Rp = Rpo2(k,dx)
    Rm = Rmo2(k,dx)
    
    M = MA(k,dx)
    M3v = M3(k,dx)
    
    Gav = GA(h,k,M,dx)
    G2v = G2(h,k,dx)
    G4v = G4(h,k,M3v,dx)
    GFv = GNFEM(h,Rm,Rp,k,dx)
    
    erro1[i] = abs(Gav - G2v) / abs(Gav)
    erro2[i] = abs(Gav - G2v) / abs(Gav)
    erro3[i] = abs(Gav - G4v) / abs(Gav)
    errof2[i] = abs(Gav - GFv) / abs(Gav)


plot(k*dxs,erro1,'b',label="o1")
plot(k*dxs,erro2,'r',label="o2")
plot(k*dxs,errof2,'g',label="o2FEM")
plot(k*dxs,erro3,'black',label="o3")
legend()
"""

"""
# Fng Error

h = 1.0
k = 0.5
g = 9.81

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),1.0/k,num=1000)   
n = len(dxs)

erro1 = zeros(n)
erro2 = zeros(n)
errof2 = zeros(n)
erro3 = zeros(n)

for i in range(n): 

    dx = dxs[i]
    v1 = w1/k
    v2 = w2/k
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    Rp = Rpo2(k,dx)
    Rm = Rmo2(k,dx)
    
    M2 = 1
    M3v =   M3(k,dx)

    G2v = G2(h,k,dx)
    G4v = G4(h,k,dx)
    GFv = GNFEM(h,Rm,Rp,k,dx)
    
    Fv = 3*I*k / (3 + h**2*k**2)
    F2v = h* G2v*M2*(1- exp(-I*k*dx))/dx
    F2fv = h* GFv*M2*(1- exp(-I*k*dx))/dx
    F4v = h* G4v*M2*(1- exp(-I*k*dx))/dx
    
    
    
    erro2[i] = abs(Fv - F2v) / abs(Fv)
    errof2[i] = abs(Fv - F2fv) / abs(Fv)
    erro3[i] = abs(Fv - F4v) / abs(Fv)
"""
    

    