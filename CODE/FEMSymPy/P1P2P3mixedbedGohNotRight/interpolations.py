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

def cubicCentInterior(a,b,c,d,f,g,h,i,x):
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




a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')

qjm3 = Symbol('qjm3')
qjm2 = Symbol('qjm2')
qjm1 = Symbol('qjm1')
qj = Symbol('qj')
qjp1 = Symbol('qjp1')
qjp2 = Symbol('qjp2')
qjp3 = Symbol('qjp3')

qjmh = Symbol('qjmh')
qjms = Symbol('qjms')
qjps = Symbol('qjps')
qjph = Symbol('qjph')

x = Symbol('x')

#Cubic over jm2,jm1,jp1,jp2
CubicCent = cubicCentInterior(a,b,c,d,qjm2,qjm1,qjp1,qjp2,x)

#calculate bj from b values in cell
CubicoverCell = cubicCentOverCell(a,b,c,d,qjmh ,qjms,qjps ,qjph ,x)

#Cubic over j,jp1,jp2,jp3
CubicRightInt = cubicRight(a,b,c,d,qjp3,qjp2,qjp1,qj,x)


#Cubic over j,jm1,jm2,jm3
CubicLeftInt = cubicLeft(a,b,c,d,qjm3,qjm2,qjm1,qj,x)






