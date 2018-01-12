#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 10:31:20 2017

@author: jp
"""

from sympy import *
from scipy import arange
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

Hamiltonian = []
Int1s = []
Int2s = []
Int2p1s = []
Int2p2s = []
Int2p3s = []

ts = arange(0,0.1,0.00001)
for t in ts:
    #t = i
    a = -1000
    b = 1000
    
    h1 = 10.0
    g = 9.81
    
    x = Symbol('x')
    
    
    Int1 = 0.5*integrate(g*h1*h1,(x,a,-t*sqrt(g*h1)))
    
    u2 = 2.0/3.0 * (sqrt(g*h1) + x/t)
    h2 =  4/(9*g) * (sqrt(g*h1) - x/(2*t))**2
    
    Int2 = 0.5*integrate(h2*u2*u2 + h2*h2*h2/3.0*diff(u2,x)*diff(u2,x) + g*h2*h2,(x,-t*sqrt(g*h1),2*t*sqrt(g*h1)))
 
    
    Int2p1 = 0.5*integrate(h2*u2*u2,(x,-t*sqrt(g*h1),2*t*sqrt(g*h1)))
    Int2p2 = 0.5*integrate(h2*h2*h2/3.0*diff(u2,x)*diff(u2,x),(x,-t*sqrt(g*h1),2*t*sqrt(g*h1)))
    Int2p3 = 0.5*integrate(g*h2*h2,(x,-t*sqrt(g*h1),2*t*sqrt(g*h1)))

    Hamiltonian.append(Int1 + Int2)
    Int1s.append(Int1)
    Int2s.append(Int2)
    Int2p1s.append(Int2p1)
    Int2p2s.append(Int2p2)
    Int2p3s.append(Int2p3)
    print(-t*sqrt(g*h1) , 2*t*sqrt(g*h1) , diff(u2,x))