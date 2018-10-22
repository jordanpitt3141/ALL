# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 10:28:35 2018

@author: jp
"""
from scipy import *
from pylab import plot, legend, xlabel, ylabel
x = linspace(0,2,50)

y = tanh(x)
z = x

plot(x,y, label="Tanh(x)")
plot(x,x, label="x")
#plot(x,y - x , label="Tanh(x) - x")
legend(loc='upper left')
xlabel('x')
ylabel('y')
