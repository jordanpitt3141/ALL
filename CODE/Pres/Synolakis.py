# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 09:50:15 2018

@author: jp
"""
from scipy import *
import csv
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 

def cot(x):
    return 1.0/ tan(x)
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton(x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
    
def ForcedbedM(x,t,beta,a0,a1,x0):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        
        if (x[i] <= cot(beta)):
            b[i] = - x[i]*tan(beta)
            
        else:
            b[i] = -1
        
        if b[i] >= 0:
            h[i] = 0
            u[i] = 0
            w[i] = b[i]
        else:
            w[i] = soliton(x[i] - x0,t,g,a0,a1) - a0
            h[i] = w[i] - b[i]
            u[i]= -c* (1 - a0/ (w[i] + a0))
        
    b0 =- (x[0] - dx)*tan(beta)         
    G = getGfromupy(h,u,b,0,0,0,a0,b0,-1,dx)     

    return h,u,G,b,w
  
wdir = "/home/jp/Documents/PhD/project/data/2018/raw/Presentation/Syn/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

dx = 0.1
  
g = 1

H = 0
d = 1
x1 = 38.5
x0 = 19.85
beta = arctan(1.0/x0)

startx = -30
sx = startx
endx = 100
ex = endx

startt = 0.0
st = startt
endt = 71
et = endt
dx = 0.05
l =  0.1
dt = l*dx

t = startt

ts = [0]
x = arange(startx,endx +0.1*dx, dx)
n = len(x)

for t in ts:
    
    h,u,G,b,w = ForcedbedM(x,t,beta,d,H,x1)
    

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

    s = wdir + "w" + str(t) + ".dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",w[i])
            file1.write(s)