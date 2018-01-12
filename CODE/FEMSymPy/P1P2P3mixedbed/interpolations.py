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

def cubicCent(a,b,c,d,f,g,h,i,x):
    system = Matrix( (((-2*x)**3,(-2*x)**2,(-2*x),1,f), ((-1*x)**3,(-1*x)**2,(-1*x),1,g), \
                  ((x)**3,(x)**2,x,1,h), ((2*x)**3,(2*x)**2,(2*x),1,i)))

    return solve_linear_system(system, a, b,c,d)
 
def cubicCentOverCell(a,b,c,d,f,g,h,i,x):
    system = Matrix( (((-x/2)**3,(-x/2)**2,(-x/2),1,f), ((-x/6)**3,(-x/6)**2,(-x/6),1,g), \
                  ((x/6)**3,(x/6)**2,x/6,1,h), ((x/2)**3,(x/2)**2,(x/2),1,i)))

    return solve_linear_system(system, a, b,c,d)

def cubicLeft(a,b,c,d,f,g,h,i,x):
    system = Matrix( (((-3*x)**3,(-3*x)**2,(-3*x),1,f), ((-2*x)**3,(-2*x)**2,(-2*x),1,g), \
                  ((-x)**3,(-x)**2,-x,1,h), (0,0,0,1,i)))

    return solve_linear_system(system, a, b,c,d)   

def cubicRight(a,b,c,d,f,g,h,i,x):
    system = Matrix( (((3*x)**3,(3*x)**2,(3*x),1,f), ((2*x)**3,(2*x)**2,(2*x),1,g), \
                  ((x)**3,(x)**2,x,1,h), (0,0,0,1,i)))

    return solve_linear_system(system, a, b,c,d)  

def quadCent(a,b,c,f,g,h,x):
    system = Matrix( (((-x)**2,(-x),1,f), (0,0,1,g), \
                  ((x)**2,x,1,h)))

    return solve_linear_system(system, a, b,c)
 
def quadCentCell(a,b,c,f,g,h,x):
    system = Matrix( (((-x/2)**2,(-x/2),1,f), (0,0,1,g), \
                  ((x/2)**2,x/2,1,h)))

    return solve_linear_system(system, a, b,c)    



a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')

f = Symbol('f')
g = Symbol('g')
h = Symbol('h')
i = Symbol('i')

x = Symbol('x')


CubicCent = cubicCent(a,b,c,d,f,g,h,i,x)
QuadCent = quadCent(a,b,c,f,g,h,x)

CubicRight = cubicRight(a,b,c,d,f,g,h,i,x)
CubicLeft = cubicLeft(a,b,c,d,f,g,h,i,x)

CubicCell = cubicCentOverCell(a,b,c,d,f,g,h,i,x)
QuadCell = quadCentCell(a,b,c,f,g,h,x)







