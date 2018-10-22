# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 10:28:35 2018

@author: jp
"""
from scipy import *
from pylab import plot, legend, xlabel, ylabel, xlim, ylim
x = linspace(0,1000,1000)

y = 1*(x > 500) + 1.8*(x <= 500)

plot(x,y)
#plot(x,y - x , label="Tanh(x) - x")
xlabel('x(m)')
ylabel('h(m)')
ylim([0.8,2])
