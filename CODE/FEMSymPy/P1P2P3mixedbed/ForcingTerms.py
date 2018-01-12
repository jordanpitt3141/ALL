# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv  
from sympy import Matrix, solve_linear_system

def QuadraticInitial(a,b,c,x):   

    return a*x**2 + b*x + c
  



h1 = Symbol('h1')
h2 = Symbol('h2')
h3 = Symbol('h3')

u1 = Symbol('u1')
u2 = Symbol('u2')
u3 = Symbol('u3')

b1 = Symbol('b1')
b2 = Symbol('b2')
b3 = Symbol('b3')

g = Symbol('g')

x = Symbol('x')


h = QuadraticInitial(h1,h2,h3,x)
u = QuadraticInitial(u1,u2,u3,x)
b = QuadraticInitial(b1,b2,b3,x)

hx = diff(h,x)
ux = diff(u,x)
bx = diff(b,x)

hxx = diff(hx,x)
uxx = diff(ux,x)
bxx = diff(bx,x)

G = u*h*(1 + hx*bx + h*bxx/2 + bx*bx) - diff(h*h*h*ux/3,x)

F1 = -diff(u*h,x)

F2 = -diff(G*u + g*h*h/2 - 2*h/3*h*h*ux**2 + h**2*u*ux*bx ,x) - h*h*u*ux*bxx/2 + h*u*u*bx*bxx - g*h*bx

F2S = - h*h*u*ux*bxx/2 + h*u*u*bx*bxx - g*h*bx
F2F = -diff(G*u + g*h*h/2 - 2*h/3*h*h*(ux)**2 + h**2*u*ux*bx ,x) 
F2Fp1 = -diff(G*u,x) 
F2Fp2 = -diff(g*h*h/2,x) 
F2Fp3 = -diff(- 2*h/3*h*h*(ux)**2 ,x)
F2Fp4 = -diff( h**2*u*ux*bx ,x)
 

F2nodx = -(G*u + g*h*h/2 - 2*h/3*h*h*(ux)**2 + h**2*u*ux*bx)






