# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hDef(x,t,x0,c,h0):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        if x[i] < x0 - c*t:
            h[i] = h0
            u[i] = 0
        elif(x[i] >=  x0 - c*t and x[i] <x0 + c*t):
            h[i] = (-h0 / (2*c*t))*(x[i] - (x0 - c*t)) + h0
            u[i] = (h0 / (2*c*t))*(x[i] - (x0 - c*t))
        elif(x[i] >=  x0 + c*t):
            h[i] = 0
            u[i] = 0
    return h,u

        
g =9.81
x0 = 10.0
c = 5.0
h0 = 1.0
t = 10.0

x = arange(-100,100,0.1)

h,u = hDef(x,t,x0,c,h0)







