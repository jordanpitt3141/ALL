# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
import csv
import os
import time
from pylab import plot

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

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])


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
  


#Soliton  
wdirb = "/home/jp/Documents/PhD/project/data/2019/TimeSoliton/FEVM/" 
if not os.path.exists(wdirb):
    os.makedirs(wdirb)

TotTimes = []

trials = 11
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
    
    h,u,G = solitoninit(a0,a1,g,x,startt,dx)

    
    n = len(x)  
    
    GhnBC = 3
    unBC = 3
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    
    uMbeg = zeros(GhnBC)
    uMend = zeros(GhnBC)
    hMbeg = a0*ones(GhnBC)
    hMend = a0*ones(GhnBC) 
    GMbeg = zeros(GhnBC)
    GMend = zeros(GhnBC)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend)
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
       
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    
    
    t = 0.0
    stepcount = 0  
    #Just an FEM solve here
    while t < endt: 
        evolvewrapForcing(G_c,h_c,hMbeg_c,hMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,nGhhbc,nubc,theta,hhbc_c,Ghbc_c,ubc_c)
        t = t + dt
        print(t)
        stepcount = stepcount + 1
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
    
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    

    s = wdir +  "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'h', 'G' , 'u(m/s)' ])        
                   
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t),str(x[j]), str(hC[j]) , str(GC[j]) , str(uC[j])])
    

    end_time = time.time()   
    TotTime = end_time - start_time   
    
    s = wdir +  "TimeInfo.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['st' ,'et','et-st',"stepcount" ,'(et-  st) /stepcount ' ])        
        writefile2.writerow([str(start_time),str(end_time),str(TotTime),str(stepcount),str(TotTime / stepcount)])
        
    TotTimes.append(TotTime)
    
    
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)


avgTotTime = sum(TotTimes)/len(TotTimes)
avgTotTimeperStep = avgTotTime/stepcount

s = wdirb +  "TimeInfo.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['# trials''avg time','step' ,'avg time per step'])        
    writefile2.writerow([str(len(TotTimes)),str(avgTotTime),str(stepcount),str(avgTotTimeperStep)])