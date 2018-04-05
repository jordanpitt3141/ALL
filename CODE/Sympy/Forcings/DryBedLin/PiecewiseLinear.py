# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from sympy.plotting import plot


def hxt(x,t,x0,h0,c):
    nh1 = h0
    nh2 = (-h0 / (2*c*t))*(x - (x0 - c*t)) + h0
    nh3 = 0
    b1 = (x0 - c*t)
    b2 = (x0 - c*t,x0 + c*t)
    b3 = (x0 + c*t)
    return nh1,nh2,nh3,b1,b2,b3

def uxt(x,t,x0,h0,c):
    nu1 = 0
    nu2 = (h0 / (2*c*t))*(x - (x0 - c*t))
    nu3 = 0
    return nu1,nu2,nu3

    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def hFluxx(h1,h2,h3,u1,u2,u3):
    #First Part
    F1 = diff(h1*u1,x)
    #Second Part
    F2 = diff(h2*u2,x)
    #Third Part
    F3 = diff(h3*u3,x)
    
    return F1, F2, F3

def GFluxx(h1,h2,h3,u1,u2,u3,G1,G2,G3,g):
    #First Part
    u1x = diff(u1,x)
    F1 = diff(G1*u1 + g/2.0*h1*h1 - 2*h1*h1*h1*u1x*u1x/3,x)
    #Second Part
    u2x = diff(u2,x)
    F2 = diff(G2*u2 + g/2.0*h2*h2 - 2*h2*h2*h2*u2x*u2x/3,x)
    #Third Part
    u3x = diff(u3,x)
    F3 = diff(G3*u3 + g/2.0*h3*h3 - 2*h3*h3*h3*u3x*u3x/3,x)
    
    return F1, F2, F3
    
g =9.81

x0 = 0
h0 = 1
c = 3

x = Symbol('x')
t = Symbol('t')


h1,h2,h3,b1,b2,b3 = hxt(x,t,x0,h0,c)
u1,u2,u3 = uxt(x,t,x0,h0,c)

G1 = Gxt(h1,u1)
G2 = Gxt(h2,u2)
G3 = Gxt(h3,u3)

#Derivative in t

#h
h1t = diff(h1,t)
h2t = diff(h2,t)
h3t = diff(h3,t)

#G
G1t = diff(G1,t)
G2t = diff(G2,t)
G3t = diff(G3,t)

hF1,hF2,hF3 = hFluxx(h1,h2,h3,u1,u2,u3)

GF1, GF2, GF3 = GFluxx(h1,h2,h3,u1,u2,u3,G1,G2,G3,g)




