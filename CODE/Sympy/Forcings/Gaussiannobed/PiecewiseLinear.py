# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from sympy.plotting import plot


def hxt(x,t,c,a0,a1,a2,a3):
    return a0 + a1*exp(-((x - c*t)  - a2)**2/(2*a3))

def uxt(x,t,c,a0,a1,a2,a3,b1):
    return b1*exp(-((x - c*t)  - a2)**2/(2*a3))

def bxt(x,b2):
    return sin(b2*x)
    
def Gxt(h,u,b):
    ux = diff(u,x)
    hx = diff(h,x)
    bx = diff(b,x)
    bxx = diff(bx,x)
    return u*h*(1 + hx*bx + h*bxx/2 + bx*bx) - diff(h**3/3*ux,x)

def hFluxx(h,u):
    #First Part
    hF = diff(h*u,x)
    
    return hF

def GFluxx(h,u,b,G,g):
    #First Part
    ux = diff(u,x)
    bx = diff(b,x)
    bxx = diff(bx,x)

    GF = diff(G*u + g/2.0*h**2 - 2*h**3*ux*ux/3.0  + h*h*u*ux*bx,x)

    
    return GF

def Source(h,u,b,G,g):
    ux = diff(u,x)
    bx = diff(b,x)
    bxx = diff(bx,x)
    
    SourceG = -(h*h*u*ux*bxx/2 - h*u*u*bx*bxx + g*h*bx)
    
    return SourceG
    
    
    
g =9.81
c = 2.0
a0 = 0.0
a1 = 1.0
a2 = 2.0
a3 = 3.0
b1 = 0.5
b2 = 0.1


x = Symbol('x')
t = Symbol('t')

a0 = Symbol('a0')
a1 = Symbol('a1')
a2 = Symbol('a2')
a3 = Symbol('a3')
a4 = Symbol('a4')
a5 = Symbol('a5')
a6 = Symbol('a6')


h = hxt(x,t,c,a0,a1,a2,a3)
u = uxt(x,t,c,a0,a1,a2,a3,b1)
b = bxt(x,b2)


G = Gxt(h,u,b)


#Derivative in t
ht = diff(h,t)

Gt = diff(G,t)

hF = hFluxx(h,u)

GF = GFluxx(h,u,b,G,g)

GS = Source(h,u,b,G,g)






