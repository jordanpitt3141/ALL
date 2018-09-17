# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time



def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


    
def ForcedbedM(x,t,a0,a1,a2,a3,a4):
    n = len(x)
    h = zeros(n)

    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))

       
    return h
    
def SWWETransform(h,h0,htol):
    n = len(x)
    ht = zeros(n)
    for i in range(n):
        """
        if(h[i] > htol):
            ht[i] = h[i]*((h[i] + h0) / (h[i] + htol))
        else:
            ht[i] = 0
        """
        ht[i] = h[i]*((h[i] + h0) / (h[i] + htol))
    return ht

def LiftTran(h,h0):
    n = len(x)
    ht = zeros(n)
    for i in range(n):
        ht[i] = (h[i] + h0)
    return ht
    
def LiftSWWETran(h,h0,htol):
    n = len(x)
    ht = zeros(n)
    for i in range(n):
        ht[i] =  ((h[i] + h0)**2 / (h[i] + htol))
    return ht
    
def LiftSWWETranSep(h,h0,htol):
    n = len(x)
    ht = zeros(n)
    for i in range(n):
        ht[i] =  ((h[i]**2 + h0**2) / (h[i] + htol))
    return ht

def LiftSWWETranNew(h,h0,htol):
    n = len(x)
    ht = zeros(n)
    st = zeros(n)
    k =  1
    for i in range(n):
        st[i] = (1.0 / (1 + exp(2*k*(x[i] - htol))))
        ht[i] = h[i] + h0*(1 / (1 + exp(2*k*(h[i] - htol))))
    return ht,st

def LiftMax(h,h0):
    n = len(x)
    ht = zeros(n)
    for i in range(n):
        ht[i] = max(h[i],h0)
    return ht

a0 =0
a1 = 10**-6
a2 = 0
a3 = 0
a4 = 5


sx = -40
ex = 40
dx = 0.01
x,t = makevar(sx,ex,dx,0,0,1)

t0 = 0

h = ForcedbedM(x,t0,a0,a1,a2,a3,a4)
hcut = h*(h > 10**-12)
u1 = ForcedbedM(x,t0,a0,a1,a2,a3,a4)
htSWWE = SWWETransform(hcut,10**-8,10**-12) 
htSWWEnorm = norm(htSWWE - hcut, ord = 1)/ norm(hcut, ord = 1)

G1 = u1*hcut

u1A = G1 / hcut
u1tA = G1 / htSWWE

"""
htLift = LiftTran(h,10**-3)
htLiftSWWE = LiftSWWETran(h,10**-16,10**-30)
htLiftSWWESep = LiftSWWETranSep(h,10**-16,10**-30)
htLiftSWWENew,st = LiftSWWETranNew(h,10**-4,10**-2)
hliftmax = LiftMax(h,10**-4)

htSWWEnorm = norm(htSWWE - h, ord = 1)/ norm(h, ord = 1)
htLiftnorm = norm(htLift - h, ord = 1)/ norm(h, ord = 1)
htLiftSWWEnorm = norm(htLiftSWWE  - h, ord = 1)/ norm(h, ord = 1)
"""