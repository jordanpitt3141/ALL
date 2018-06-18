# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def lineardiscP(x,x0,x1):
    return ((x - S(x1)) / (S(x0) - S(x1)) ,S(x0), S(x1))

def lineardiscM(x,x0,x1):
    return ((x - S(x0)) / (S(x1) - S(x0)) ,S(x0), S(x1))

        
def quadcontEM(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x0),S(x2))

def quadcontEP(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x2),S(x0))


def quadcontM(x,xm1,x0,x1):
    fexp = (x - S(xm1))*(x - S(x1)) / ((S(x0) - S(xm1))*(S(x0) - S(x1)))

    return (fexp,S(xm1),S(x1))

def phideriv(x,phi):
        return (diff(phi[0],x),phi[1],phi[2])        

x = Symbol('x')

xjm3o2 = "-2"
xjm1 = "-3/2"
xjm1o2 = "-1"
xj = "0"
xjp1o2 = "1"
xjp1 = "3/2"
xjp3o2= "2"

phijm1o2R = quadcontEP(x,xjm1o2,xj,xjp1o2)
phij = quadcontM(x,xjm1o2,xj,xjp1o2)
phijp1o2L = quadcontEP(x,xjp1o2,xj,xjm1o2)
phijp1o2R = quadcontEM(x,xjp1o2,xjp1,xjp3o2)
phijp1 = quadcontM(x,xjp1o2,xjp1,xjp3o2)
phijp3o2L = quadcontEM(x,xjp3o2,xjp1,xjp1o2)


#derivatives
dphijm1o2R = phideriv(x,phijm1o2R)
dphij = phideriv(x,phij )
dphijp1o2L = phideriv(x,phijp1o2L)
dphijp1o2R = phideriv(x,phijp1o2R)
dphijp1 = phideriv(x,phijp1 )
dphijp3o2L = phideriv(x,phijp3o2L)


wjm1o2p = lineardiscP(x,xjm1o2,xjp1o2) 
wjp1o2m = lineardiscM(x,xjm1o2,xjp1o2) 
wjp1o2p = lineardiscP(x,xjp1o2,xjp3o2) 
wjp3o2m = lineardiscM(x,xjp1o2,xjp3o2) 

#Gmatrix

Gvec00 = integrate(wjm1o2p[0]*phijm1o2R[0] ,(x,xjm1o2,xjp1o2))
Gvec01 = integrate(wjp1o2m[0]*phijm1o2R[0] ,(x,xjm1o2,xjp1o2))

Gvec10 = integrate(wjm1o2p[0]*phij[0] ,(x,xjm1o2,xjp1o2))
Gvec11 = integrate(wjp1o2m[0]*phij[0] ,(x,xjm1o2,xjp1o2))

Gvec20 = integrate(wjm1o2p[0]*phijp1o2L[0] ,(x,xjm1o2,xjp1o2))
Gvec21 = integrate(wjp1o2m[0]*phijp1o2L[0] ,(x,xjm1o2,xjp1o2))
Gvec22 = integrate(wjp1o2p[0]*phijp1o2R[0] ,(x,xjp1o2,xjp3o2))
Gvec23 = integrate(wjp3o2m[0]*phijp1o2R[0] ,(x,xjp1o2,xjp3o2))

Gvec30 = integrate(wjp1o2p[0]*phijp1[0] ,(x,xjp1o2,xjp3o2))
Gvec31 = integrate(wjp3o2m[0]*phijp1[0] ,(x,xjp1o2,xjp3o2))

Gvec40 = integrate(wjp1o2p[0]*phijp3o2L[0] ,(x,xjp1o2,xjp3o2))
Gvec41 = integrate(wjp3o2m[0]*phijp3o2L[0] ,(x,xjp1o2,xjp3o2))

#PhiMatrix

Phivec0 =integrate(phijp1o2L[0]*phijm1o2R[0] ,(x,xjm1o2,xjp1o2))
Phivec1 =integrate(phijp1o2L[0]*phij[0] ,(x,xjm1o2,xjp1o2))
Phivec2 =integrate(phijp1o2L[0]*phijp1o2L[0] ,(x,xjm1o2,xjp1o2))
Phivec3 =integrate(phijp1o2R[0]*phijp1o2R[0] ,(x,xjp1o2,xjp3o2))
Phivec4 =integrate(phijp1o2R[0]*phijp1[0] ,(x,xjp1o2,xjp3o2))
Phivec5 =integrate(phijp1o2R[0]*phijp3o2L[0] ,(x,xjp1o2,xjp3o2))

#dPhiMatrix
dPhivec0 =integrate(dphijp1o2L[0]*dphijm1o2R[0] ,(x,xjm1o2,xjp1o2))
dPhivec1 =integrate(dphijp1o2L[0]*dphij[0] ,(x,xjm1o2,xjp1o2))
dPhivec2 =integrate(dphijp1o2L[0]*dphijp1o2L[0] ,(x,xjm1o2,xjp1o2))
dPhivec3 =integrate(dphijp1o2R[0]*dphijp1o2R[0] ,(x,xjp1o2,xjp3o2))
dPhivec4 =integrate(dphijp1o2R[0]*dphijp1[0] ,(x,xjp1o2,xjp3o2))
dPhivec5 =integrate(dphijp1o2R[0]*dphijp3o2L[0] ,(x,xjp1o2,xjp3o2))