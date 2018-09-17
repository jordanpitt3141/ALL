# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 09:50:15 2018

@author: jp
"""
from scipy import *
import csv
import os

def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b
  
wdir = "/home/jp/Documents/PhD/project/data/2018/raw/Presentation/ForcedDryBed/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

dx = 0.1
  
a6= 1.0
a7 = 2*pi/50.0

width = 2*(2*pi/a7)
    
a0 = 0.0
a1 = 0.5
a2 =  ((2*pi) / a7)/10.0
a3 = -pi/2.0/a7 -width/4.0
a4 = width/2**6
a5 = a1


g = 9.81

startx = -pi/2.0/a7 -width
sx= startx
endx = -pi/2.0/a7 +width
ex = endx
startt = 0.0
endt = (2*pi/a7) / a2
et = endt

ts = [0,2.5,5.0,7.5,10.0]
x = arange(startx,endx +0.1*dx, dx)
n = len(x)

for t in ts:
    
    h,u,G,b = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    s = wdir + "h" + str(t) + ".dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i])
            file1.write(s)

    s = wdir + "u" + str(t) + ".dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",u[i])
            file1.write(s)
            
    s = wdir + "G" + str(t) + ".dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",G[i])
            file1.write(s)
            
    s = wdir + "b" + str(t) + ".dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
            file1.write(s)