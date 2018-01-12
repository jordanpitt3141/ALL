# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 12:02:43 2017

@author: jp
"""

from scipy import *
from pylab import plot

x = linspace(0.025,0.075,num=30)
x1 = linspace(0.2025,0.2275,num=30)

z = linspace(-1,1,num=1000)

y = 0.8 + 0.05*sin((x - 0.05) * (2*pi/0.05))

y1 = 0.8 + 0.05*sin((x1 - 0.2025) * (2*pi/0.05))

for i in range(len(x)):
   print (x1[i],y1[i])