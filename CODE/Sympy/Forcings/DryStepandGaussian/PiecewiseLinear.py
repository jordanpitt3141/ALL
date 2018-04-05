# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from sympy.plotting import plot


def hxt(x,t,x0,h0,c,l):
    nh1 = 0
    nh2 = h0
    nh3 = 0
    b1 = (x0 + c*t) - l/2
    b2 = (x0 + c*t) + l/2
    return nh1,nh2,nh3,b1,b2

def uxt(x,t,x0,c,a0,a1,a2):
    nu1 = a0*exp(-((x - c*t)  - a1)**2/(2*a2))
    return nu1

    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def hFluxx(h1,h2,h3,u):
    #First Part
    F1 = diff(h1*u,x)
    #Second Part
    F2 = diff(h2*u,x)
    
    F3 = diff(h3*u,x)

    return F1, F2,F3

def GFluxx(h1,h2,h3,u,G1,G2,G3,g):
    #First Part
    ux = diff(u,x)
    F1 = diff(G1*u + g/2.0*h1*h1 - 2*h1*h1*h1*ux*ux/3,x)
    #Second Part

    F2 = diff(G2*u + g/2.0*h2*h2 - 2*h2*h2*h2*ux*ux/3,x)
    
    F3 = diff(G3*u + g/2.0*h3*h3 - 2*h3*h3*h3*ux*ux/3,x)
    
    return F1, F2,F3
    
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


h1,h2,h3,b1,b2= hxt(x,t,x0,h0,c,l)
u = uxt(x,t,x0,c,a0,a1,a2)

G1 = Gxt(h1,u)
G2 = Gxt(h2,u)
G3 = Gxt(h3,u)

#Derivative in t

#h
h1t = diff(h1,t)
h2t = diff(h2,t)
h3t = diff(h3,t)

#G
G1t = diff(G1,t)
G2t = diff(G2,t)
G3t = diff(G3,t)

hF1,hF2, hF3=  hFluxx(h1,h2,h3,u)

GF1, GF2, GF3 = GFluxx(h1,h2,h3,u,G1,G2,G3,g)




