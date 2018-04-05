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

    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def hFluxx(h,u):
    #First Part
    hF = diff(h*u,x)
    
    return hF

def GFluxx(h,u,G,g):
    #First Part
    ux = diff(u,x)
    GF = diff(G*u + g/2.0*h*h - 2*h*h*h*ux*ux/3,x)

    
    return GF
    
g =9.81
c = 2.0
a0 = 0.0
a1 = 1.0
a2 = 2.0
a3 = 3.0
b1 = 0.5


x = Symbol('x')
t = Symbol('t')
"""
a0 = Symbol('a0')
a1 = Symbol('a1')
a2 = Symbol('a2')
a3 = Symbol('a3')
b1 = Symbol('b1')
c = Symbol('c')
"""

h = hxt(x,t,c,a0,a1,a2,a3)
u = uxt(x,t,c,a0,a1,a2,a3,b1)

G = Gxt(h,u)


#Derivative in t
ht = diff(h,t)

Gt = diff(G,t)

hF = hFluxx(h,u)

GF = GFluxx(h,u,G,g)






