# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *



def qoff(k,w,xoff,toff,x,t,dx,dt):
    return exp(I*k*xoff*dx)*exp(I*w*toff*dt)

def qdiffA(k,w,x,t,derivvartup):
    dinit = exp(I*k*x + I*w*t)
    n = len(derivvartup)
    d = dinit
    for i in range(n):
        d = diff(d,derivvartup[i])
    
    dcoeff = simplify(d/dinit)
    return dcoeff


def dddFD(k,w,x,t,dx,dt):
    qjp2 = qoff(k,w,2,0,x,t,dx,dt) 
    qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
    qj = qoff(k,w,0,0,x,t,dx,dt) 
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
    qjm2 = qoff(k,w,-2,0,x,t,dx,dt)
    
    Term = simplify( expand(simplify(qjp2 - qjm2),trig=True) + simplify( - 2*qjp1 + 2*qjm1))/dx/dx/dx/2
    
    
    return Term

def ddFD(k,w,x,t,dx,dt):
    qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
    qj = qoff(k,w,0,0,x,t,dx,dt) 
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
    
    Term = simplify(qjm1 - 2*qj + qjp1)/dx/dx
    
    
    return Term

def dFD(k,w,x,t,dx,dt):
    qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
    qj = qoff(k,w,0,0,x,t,dx,dt) 
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
    
    Term = simplify(-qjm1 + qjp1)/dx/2
    
    
    return Term


def NaiveNumMass(k,w,x,t,U,H,dx,dt):
    
    d2 = dFD(k,w,x,t,dx,dt)
    dd2 = ddFD(k,w,x,t,dx,dt)
    ddd2 = dddFD(k,w,x,t,dx,dt)
    
    etan = simplify(-2*dt*(U*d2))
    upsn = simplify(-2*dt*(H*d2))
    
    etanm1 = 1
    upsnm1 = 0
    
    return (simplify(etanm1),simplify(etan),simplify(upsnm1),simplify(upsn))


def NaiveNumMome(k,w,x,t,U,H,dx,dt):
    
    d2 = dFD(k,w,x,t,dx,dt)
    dd2 = ddFD(k,w,x,t,dx,dt)
    ddd2 = dddFD(k,w,x,t,dx,dt)
    
    RHScoeff = 1 - H**2/3*dd2
    
    etan = simplify(-g*2*dt*d2 / (RHScoeff))
    upsn = simplify(2*dt*(-U*d2 + H**2/3*(U*ddd2))/RHScoeff)
    
    etanm1 = 0
    upsnm1 = 1
    
    return (simplify(etanm1),simplify(etan),simplify(upsnm1),simplify(upsn))    

def LWNumMass(k,w,x,t,U,H,dx,dt,A00,A01):
    qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
    qj = qoff(k,w,0,0,x,t,dx,dt) 
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
    
    MomeCoeff = NaiveNumMome(k,w,x,t,U,H,dx,dt)
    
    hnphjphe = simplify(qjp1 + qj)/2 - dt/(2*dx)*(U*simplify(qjp1 - qj))
    hnphjphu = - dt/(2*dx)*(H*simplify(qjp1 - qj))
    
    hnphjmhe = simplify(qj + qjm1)/2 - dt/(2*dx)*(U*simplify(qj - qjm1))
    hnphjmhu = - dt/(2*dx)*(H*simplify(qj - qjm1))
    
    print( - dt/dx*(U*simplify(hnphjphu- hnphjmhu)))    
    
    unphjphu = (MomeCoeff[3] + 1)
    #unphjphunm1 = 1/4
    unphjphh = MomeCoeff[1]/4
    
    #h contribution not including u update
    hnp1je = qj - dt/dx*(U*simplify(hnphjphe - hnphjmhe))
    hnp1je1  = collect(expand(hnp1je),U*dt/dx)
    hnp1je1 = hnp1je1.subs(exp(I*dx*k)/2 + exp(-I*dx*k)/2, cos(dx*k))
    hnp1je1 = hnp1je1.subs(-exp(I*dx*k)/2 + exp(-I*dx*k)/2, -I*sin(dx*k))
    
    #addh from u caclulation
    hnp1je1 = hnp1je1 -dt/dx*H*simplify((1 - exp(-I*k*dx))*(exp(I*k*dx) +1))*unphjphh   
    
    Masscoeffhn = hnp1je1
    Masscoeffhnm1 = 0 
    Masscoeffunm1 =-dt/dx*H*simplify((1 - exp(-I*k*dx))*(exp(I*k*dx) +1))/4
    
    Masscoeffunp1 = - dt/dx*(U*simplify(hnphjphu- hnphjmhu))
    Masscoeffunp1 = Masscoeffunp1.subs((exp(I*dx*k) - 1)**2*exp(-I*dx*k),2*(cos(k*dx) -1))
    Masscoeffun =-dt/dx*H*simplify((1 - exp(-I*k*dx))*(exp(I*k*dx) +1))*unphjphu/4 +Masscoeffunp1  
    
    return (simplify(Masscoeffhnm1),simplify(Masscoeffhn),simplify(Masscoeffunm1),simplify(Masscoeffun))
    
    
    

x = Symbol('x')
t = Symbol('t')
dx = Symbol('dx')
dt = Symbol('dt')

U = Symbol('U')
H = Symbol('H')
h = Symbol('h')
u = Symbol('u')

k = Symbol('k')
w = Symbol('w')

g = Symbol('g')

TaylorSeriesTermNum= 10
#Elliptic Equation

# qxx  -> q | qxx = Cq
d2 = dFD(k,w,x,t,dx,dt)
dd2 = ddFD(k,w,x,t,dx,dt)
ddd2 = dddFD(k,w,x,t,dx,dt)

MassCoeff = NaiveNumMass(k,w,x,t,U,H,dx,dt)
MomeCoeff = NaiveNumMome(k,w,x,t,U,H,dx,dt)

LWMassCoeff = LWNumMass(k,w,x,t,U,H,dx,dt,MomeCoeff[1],MomeCoeff[3])
