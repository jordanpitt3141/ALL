# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2dc import *
from scipy import *
import csv
import os
import time

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
        

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(a0,a1,g,x,t0,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* (1 - a0 / h[i])
        G[i] = 2.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**4*h[i] + h[i]*u[i] - 4.0/3*a0*a1**2*c*k**2*sech(k*(x[i] - c*t0))**4*tanh(k*(x[i] - c*t0))**2 - 4.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**2*h[i]*tanh(k*(x[i] - c*t0))**2
        
    return h,u,G
    
#Soliton Problem    
wdirb = "/home/jp/Documents/PhD/project/data/2019/TimeSoliton/FDVM2/"

if not os.path.exists(wdirb):
    os.makedirs(wdirb)

TotTimes = []

trials = 10
for nt in range(trials):
    
    ki = 11
    start_time = time.time()
    wdir = wdirb + str(ki) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    c = sqrt(g*(a0 + a1))
    
    
    startx = -250
    endx = 250
    
    startt = 0.0
    endt = 50
    
    
    dx = 100.0/ 2**ki
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    
    t = startt
    theta = 1.2
    
    x = arange(startx,endx +0.1*dx, dx)
    
    n = len(x)
    nBC = 3
    nBCs = 4
    
    h,u,G = solitoninit(a0,a1,g,x,0,dx)
    
    u0a = zeros(nBCs)
    u1a = zeros(nBCs)    
    h0a = a0*ones(nBCs)
    h1a = a0*ones(nBCs)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    h0_c  = copyarraytoC(h0a)
    h1_c  = copyarraytoC(h1a)
    u0_c  = copyarraytoC(u0a)
    u1_c  = copyarraytoC(u1a)
    u_c = mallocPy(n)
    
    stepcount = 0
    while t < endt :        
        evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs,theta)
        print (t)
        
        t = t + dt
        stepcount = stepcount + 1
    
    getufromG(h_c,G_c,u0a[-1],u1a[0],h0a[-1],h1a[0], dx ,n,u_c)
    uF = copyarrayfromC(u_c,n)
    GF = copyarrayfromC(G_c,n)
    hF = copyarrayfromC(h_c,n)  
    

    s = wdir +  "outlast.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'h', 'G' , 'u(m/s)' ])        
                   
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t),str(x[j]), str(hF[j]) , str(GF[j]) , str(uF[j])])
            
    end_time = time.time()   
    TotTime = end_time - start_time   
    
    s = wdir +  "TimeInfo.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['st' ,'et','et-st',"stepcount" ,'(et-  st) /stepcount ' ])        
        writefile2.writerow([str(start_time),str(end_time),str(TotTime),str(stepcount),str(TotTime / stepcount)])
        
    TotTimes.append(TotTime)

    deallocPy(G_c)
    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(u_c)
    
avgTotTime = sum(TotTimes)/len(TotTimes)
avgTotTimeperStep = avgTotTime/stepcount

s = wdirb +  "TimeInfo.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['# trials''avg time','step' ,'avg time per step'])        
    writefile2.writerow([str(len(TotTimes)),str(avgTotTime),str(stepcount),str(avgTotTimeperStep)])
           
