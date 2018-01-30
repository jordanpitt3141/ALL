# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:25:13 2018

@author: jp
"""

h0 = 1
h1 = 1.8
a = 0
b = 1000
g = 9.81

from scipy import *

def Hamiltonian(h0,h1,g,a,b,alpha):
    return 0.25*g*(h0**2 - h1**2 + alpha*((h1 - h0)**2)*tanh((a-b)/(2*alpha)))
    
def HamilDB(alpha,dx):
    return 10.3986*(1000 + dx) - 0.7848*((2.0/alpha)*tanh(alpha * (500.0+ 0.5*dx)))    