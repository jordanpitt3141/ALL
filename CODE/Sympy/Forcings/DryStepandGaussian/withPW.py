# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from sympy.plotting import plot


def hxt(x,t,x0,h0,c,l):
    nh = Piecewise(
        (0,x <= (x0 + c*t) - l),
        (0,x >= (x0 + c*t) + l),
        (h0,True)) 
    return nh

def uxt(x,t,x0,c,a0,a1,a2):
    nu1 = a0*exp(-((x - c*t)  - a1)**2/(2*a2))
    return nu1

    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def hFluxx(h,u):
    #First Part
    F1 = diff(h*u,x)
    #Second Part

    return F1

def GFluxx(h,u,G,g):
    #First Part
    
    ux = diff(u,x)
    F1 = diff(G*u + g/2.0*h*h - 2*h*h*h*ux*ux/3,x)
    
    return F1
    
g =9.81

l = 10.0
x0 =1.0
h0 = 1.0       
g =9.81
c = 2.0
a0 = h0
a1 = x0
a2 = l/4

x = Symbol('x')
t = Symbol('t')


h = hxt(x,t,x0,h0,c,l)
u = uxt(x,t,x0,c,a0,a1,a2)

G = Gxt(h,u)


ht = diff(h,t)


#G
Gt = diff(G,t)


hF =  hFluxx(h,u)

GF= GFluxx(h,u,G,g)




