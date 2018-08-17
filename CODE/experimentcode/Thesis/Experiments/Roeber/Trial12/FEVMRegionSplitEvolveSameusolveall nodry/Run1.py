# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time
from scipy import signal

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

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


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

def MollifyFunc(C,x,e):
    ie = 1.0 / e
    if(abs(ie*x) <1):
        return ie*C*exp(1.0/(abs(ie*x)**2 - 1))
    else:
        return 0

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton(x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)    

def Roeberflume(x,t,a0,a1,g,x0,xexp,bedexp,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    MollPy = zeros(n)
    
    D = 1.0/0.444994
    eta =3*dx
    rublen = 600
    
    for i in range(n):
        
        MollPy[i] = MollifyFunc(D,x[i] - x[n/2],eta)
        
        if(x[i] <= xexp[1]):
            bed[i] = bedexp[1]
        elif(xexp[1] < x[i] < xexp[-1]):
            j = [ nin for nin, nv in enumerate(xexp) if nv>=x[i] ][0]
            bed[i] = bedexp[j-1] + ( (bedexp[j] - bedexp[j-1]) / 0.05)*(x[i] - xexp[j-1])
            
        elif(x[i] >= xexp[-1]):
            bed[i] = bedexp[-1]

    b0 = bed[0]*ones(rublen)
    b1 = bed[-1]*ones(rublen)
    bedbc = signal.convolve(concatenate((b0,bed,b1)), concatenate((zeros(rublen),MollPy,zeros(rublen))), mode='same') / sum(MollPy)
    bed = bedbc[rublen:-rublen]
        
    for i in range(n):
        if (x[i] > x0 -50 and x[i]<x0 + 50):
            h[i] = soliton(x[i] - x0,t,g,a0,a1)
            u[i] = sqrt(g*(a0 + a1))*(1 - a0 / h[i])
        else:
            h[i] = max(0.0 - bed[i],0.0)
       
            
    G = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    return h,u,G,bed

def BejiEdge(x,hc0,vc0,ft):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    hb = 2.5
    

    
    i = n-1
    et = ft
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) 
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    #linear extrapolation
    i = n - 2
    et = 2*ft - (hc0 - hb)
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) 
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    for i in range(n-3,-1,-1):
        et = 2*(eta[i+1]) - (eta[i+2])
        
        h1 = hb + et
        c1 = sqrt(g*(hb+ et)) 
        ut = (c1*et) / (h1)
        eta[i] = et
        v[i] = ut
        
    i = -1
    et = 2*(eta[i+1]) - (eta[i+2])
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) 
    ut = (c1*et) / (h1)
    e0 = et
    v0 = ut 
    h0 = hb + e0
    
    hv = hb+ eta

    
    G = getGfromupy(hv,v,bed,v0,vc0,h0,hc0,0,0,dx)  
    return hv,v,G

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)
    
def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    #return y1  + (xi)*(y2 - y1)/(x2 - x1)  
    return y1  + (xi)*(y2 - y0)/(x2 - x0) 

  
def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var

def find_nearestidx(array1, value1):
    array2 = asarray(array1)
    idx = (abs(array2 - value1)).argmin()
    return idx



#Roeber Problem    
expdir = "../../../../../../../../data/Experimental/HIreef/Trial12/"
wdir = "../../../../../../../../data/raw/Thesis/Experiment/Roeber/Trial12/FEVM/t1/"  


if not os.path.exists(wdir):
    os.makedirs(wdir)

s = expdir + "bed.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    xexp = []
    bedexp = []
    for row in readfile:       
            xexp.append(float(row[0]))
            bedexp.append(float(row[4]))
            
s = expdir + "WG1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG1exp.append(float(row[4]))
            
s = expdir + "WG2.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG2exp = []
    for row in readfile:       
            WG2exp.append(float(row[4]))
            
s = expdir + "WG3.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG3exp = []
    for row in readfile:       
            WG3exp.append(float(row[4]))
            
s = expdir + "WG4.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG4exp = []
    for row in readfile:       
            WG4exp.append(float(row[4]))
            
s = expdir + "WG5.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG5exp = []
    for row in readfile:       
            WG5exp.append(float(row[4]))
            
s = expdir + "WG6.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG6exp = []
    for row in readfile:       
            WG6exp.append(float(row[4]))

s = expdir + "WG7.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG7exp = []
    for row in readfile:       
            WG7exp.append(float(row[4]))
 
s = expdir + "WG8.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG8exp = []
    for row in readfile:       
            WG8exp.append(float(row[4]))
            
s = expdir + "WG9.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG9exp = []
    for row in readfile:       
            WG9exp.append(float(row[4]))
            
s = expdir + "WG10.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG10exp = []
    for row in readfile:       
            WG10exp.append(float(row[4]))
            
s = expdir + "WG11.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG11exp = []
    for row in readfile:       
            WG11exp.append(float(row[4]))
            
s = expdir + "WG12.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG12exp = []
    for row in readfile:       
            WG12exp.append(float(row[4]))

s = expdir + "WG13.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG13exp = []
    for row in readfile:       
            WG13exp.append(float(row[4]))

s = expdir + "WG14.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WG14exp = []
    for row in readfile:       
            WG14exp.append(float(row[4]))
            
s = expdir + "wavegauagelocations.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WGloc = []
    for row in readfile:       
            WGloc.append(round(float(row[0]),1))  


g = 9.81
sr = 0.02
dt = sr/ (2**3)
dx = 0.1 / (2**1)
l = 1.0 / (5 + sqrt(g*4))

theta = 1.0
sx = -200
ex = 100
st = texp[0]
et = 35.2
hb = 2.5



GhnBC = 3
unBC = 3
bnBC = 4

x = arange(sx,ex +0.1*dx, dx)

wg1i = find_nearestidx(x, WGloc[0])
wg2i = find_nearestidx(x, WGloc[1])
wg3i = find_nearestidx(x, WGloc[2])
wg4i = find_nearestidx(x, WGloc[3])
wg5i = find_nearestidx(x, WGloc[4])
wg6i = find_nearestidx(x, WGloc[5])
wg7i = find_nearestidx(x, WGloc[6])
wg8i = find_nearestidx(x, WGloc[7])
wg9i = find_nearestidx(x, WGloc[8])
wg10i = find_nearestidx(x, WGloc[9])
wg11i = find_nearestidx(x, WGloc[10])
wg12i = find_nearestidx(x, WGloc[11])
wg13i = find_nearestidx(x, WGloc[12])
wg14i = find_nearestidx(x, WGloc[13])

n = len(x) 

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

xbMend = array([x[-1] + 0.5*dx, x[-1] + 5*dx/6.0, x[-1] + 7*dx/6.0, x[-1] + 1.5*dx])
xbMbeg = array([x[0] - 1.5*dx, x[0] - 7*dx/6.0,x[0] - 5*dx/6.0 , x[0] -0.5*dx])


tts = [] 
ftcs = []
nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nwg8s = []
nwg9s = []
nwg10s = []
nwg11s = []
nwg12s = []
nwg13s = []
nwg14s = []

a0 = hb
a1 = 0.75
x0 = -139.26

tij = 0.0
ij = 0
h,u,G,b = Roeberflume(x,st,a0,a1,g,x0,xexp,bedexp,dx)

hMbeg = h[0]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)
GMend = G[-1]*ones(GhnBC)
wbMbeg = b[0]*ones(GhnBC)
wbMend = b[-1]*ones(GhnBC)
wMend = h[-1]*ones(GhnBC) + wbMend
wMbeg = h[0]*ones(GhnBC) + wbMbeg
bMbeg = b[0]*ones(bnBC)
bMend = b[-1]*ones(bnBC)



h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)
b_c = copyarraytoC(b)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)
   

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)

hbwg1 = hMbeg[1]    
     
nwg1s.append(hbwg1)
nwg2s.append(hbwg1)  
nwg3s.append(hbwg1)  
nwg4s.append(hbwg1)  
nwg5s.append(hbwg1)  
nwg6s.append(hbwg1)  
nwg7s.append(hbwg1)  
nwg8s.append(hbwg1)
nwg9s.append(hbwg1)  
nwg10s.append(hbwg1)  
nwg11s.append(hbwg1)  
nwg12s.append(hbwg1)  
nwg13s.append(hbwg1)  
nwg14s.append(hbwg1) 

hbwg2 = CELLRECON(readfrommem(h_c,wg2i -1),readfrommem(h_c,wg2i) ,readfrommem(h_c,wg2i + 1),x[wg2i -1],x[wg2i],x[wg2i + 1],WGloc[1] - x[wg2i])
hbwg3 = CELLRECON(readfrommem(h_c,wg3i -1),readfrommem(h_c,wg3i) ,readfrommem(h_c,wg3i + 1),x[wg3i -1],x[wg3i],x[wg3i + 1],WGloc[2] - x[wg3i])
hbwg4 = CELLRECON(readfrommem(h_c,wg4i -1),readfrommem(h_c,wg4i) ,readfrommem(h_c,wg4i + 1),x[wg4i -1],x[wg4i],x[wg4i + 1],WGloc[3] - x[wg4i])
hbwg5 = CELLRECON(readfrommem(h_c,wg5i -1),readfrommem(h_c,wg5i) ,readfrommem(h_c,wg5i + 1),x[wg5i -1],x[wg5i],x[wg5i + 1],WGloc[4] - x[wg5i])
hbwg6 = CELLRECON(readfrommem(h_c,wg6i -1),readfrommem(h_c,wg6i) ,readfrommem(h_c,wg6i + 1),x[wg6i -1],x[wg6i],x[wg6i + 1],WGloc[5] - x[wg6i])
hbwg7 = CELLRECON(readfrommem(h_c,wg7i -1),readfrommem(h_c,wg7i) ,readfrommem(h_c,wg7i + 1),x[wg7i -1],x[wg7i],x[wg7i + 1],WGloc[6] - x[wg7i])
hbwg8 = CELLRECON(readfrommem(h_c,wg8i -1),readfrommem(h_c,wg8i) ,readfrommem(h_c,wg8i + 1),x[wg8i -1],x[wg8i],x[wg8i + 1],WGloc[7] - x[wg8i])
hbwg9 = CELLRECON(readfrommem(h_c,wg9i -1),readfrommem(h_c,wg9i) ,readfrommem(h_c,wg9i + 1),x[wg9i -1],x[wg9i],x[wg9i + 1],WGloc[8] - x[wg9i])
hbwg10 = CELLRECON(readfrommem(h_c,wg10i -1),readfrommem(h_c,wg10i) ,readfrommem(h_c,wg10i + 1),x[wg10i -1],x[wg10i],x[wg10i + 1],WGloc[9] - x[wg10i])
hbwg11 = CELLRECON(readfrommem(h_c,wg11i -1),readfrommem(h_c,wg11i) ,readfrommem(h_c,wg11i + 1),x[wg11i -1],x[wg11i],x[wg11i + 1],WGloc[10] - x[wg11i])
hbwg12 = CELLRECON(readfrommem(h_c,wg12i -1),readfrommem(h_c,wg12i) ,readfrommem(h_c,wg12i + 1),x[wg12i -1],x[wg12i],x[wg12i + 1],WGloc[11] - x[wg12i])
hbwg13 = CELLRECON(readfrommem(h_c,wg13i -1),readfrommem(h_c,wg13i) ,readfrommem(h_c,wg13i + 1),x[wg13i -1],x[wg13i],x[wg13i + 1],WGloc[12] - x[wg13i])
hbwg14 = CELLRECON(readfrommem(h_c,wg14i -1),readfrommem(h_c,wg14i) ,readfrommem(h_c,wg14i + 1),x[wg14i -1],x[wg14i],x[wg14i + 1],WGloc[13] - x[wg14i])
 


t = 0.0
tts.append(t)
#Just an FEM solve here
while t < et: 
    
    evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t,0,0,0,0,0,0,0,0,0,0)
    
    getufromGsplit(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
        
    uc0 = readfrommem(ubc_c,unBC)
    hc0 = readfrommem(h_c,0) 
    
    if(isnan(hc0)) :
        break
    
    wg1h = CELLRECON(readfrommem(h_c,wg1i -1),readfrommem(h_c,wg1i) ,readfrommem(h_c,wg1i + 1),x[wg1i -1],x[wg1i],x[wg1i + 1],WGloc[0] - x[wg1i])
    wg2h = CELLRECON(readfrommem(h_c,wg2i -1),readfrommem(h_c,wg2i) ,readfrommem(h_c,wg2i + 1),x[wg2i -1],x[wg2i],x[wg2i + 1],WGloc[1] - x[wg2i])
    wg3h = CELLRECON(readfrommem(h_c,wg3i -1),readfrommem(h_c,wg3i) ,readfrommem(h_c,wg3i + 1),x[wg3i -1],x[wg3i],x[wg3i + 1],WGloc[2] - x[wg3i])
    wg4h = CELLRECON(readfrommem(h_c,wg4i -1),readfrommem(h_c,wg4i) ,readfrommem(h_c,wg4i + 1),x[wg4i -1],x[wg4i],x[wg4i + 1],WGloc[3] - x[wg4i])
    wg5h = CELLRECON(readfrommem(h_c,wg5i -1),readfrommem(h_c,wg5i) ,readfrommem(h_c,wg5i + 1),x[wg5i -1],x[wg5i],x[wg5i + 1],WGloc[4] - x[wg5i])
    wg6h = CELLRECON(readfrommem(h_c,wg6i -1),readfrommem(h_c,wg6i) ,readfrommem(h_c,wg6i + 1),x[wg6i -1],x[wg6i],x[wg6i + 1],WGloc[5] - x[wg6i])
    wg7h = CELLRECON(readfrommem(h_c,wg7i -1),readfrommem(h_c,wg7i) ,readfrommem(h_c,wg7i + 1),x[wg7i -1],x[wg7i],x[wg7i + 1],WGloc[6] - x[wg7i])
    wg8h = CELLRECON(readfrommem(h_c,wg8i -1),readfrommem(h_c,wg8i) ,readfrommem(h_c,wg8i + 1),x[wg8i -1],x[wg8i],x[wg8i + 1],WGloc[7] - x[wg8i])
    wg9h = CELLRECON(readfrommem(h_c,wg9i -1),readfrommem(h_c,wg9i) ,readfrommem(h_c,wg9i + 1),x[wg9i -1],x[wg9i],x[wg9i + 1],WGloc[8] - x[wg9i])
    wg10h = CELLRECON(readfrommem(h_c,wg10i -1),readfrommem(h_c,wg10i) ,readfrommem(h_c,wg10i + 1),x[wg10i -1],x[wg10i],x[wg10i + 1],WGloc[9] - x[wg10i])
    wg11h = CELLRECON(readfrommem(h_c,wg11i -1),readfrommem(h_c,wg11i) ,readfrommem(h_c,wg11i + 1),x[wg11i -1],x[wg11i],x[wg11i + 1],WGloc[10] - x[wg11i])
    wg12h = CELLRECON(readfrommem(h_c,wg12i -1),readfrommem(h_c,wg12i) ,readfrommem(h_c,wg12i + 1),x[wg12i -1],x[wg12i],x[wg12i + 1],WGloc[11] - x[wg12i])
    wg13h = CELLRECON(readfrommem(h_c,wg13i -1),readfrommem(h_c,wg13i) ,readfrommem(h_c,wg13i + 1),x[wg13i -1],x[wg13i],x[wg13i + 1],WGloc[12] - x[wg13i])
    wg14h = CELLRECON(readfrommem(h_c,wg14i -1),readfrommem(h_c,wg14i) ,readfrommem(h_c,wg14i + 1),x[wg14i -1],x[wg14i],x[wg14i + 1],WGloc[13] - x[wg14i])
          

    nwg1s.append(wg1h - hbwg1 + hbwg1) 
    
    nwg2s.append(wg2h - hbwg2 + hbwg1)  
    nwg3s.append(wg3h - hbwg3 + hbwg1) 
    nwg4s.append(wg4h - hbwg4 + hbwg1) 
    nwg5s.append(wg5h - hbwg5 + hbwg1) 
    nwg6s.append(wg6h - hbwg6 + hbwg1) 
    nwg7s.append(wg7h - hbwg7 + hbwg1) 
    nwg8s.append(wg8h - hbwg8 + hbwg1)  
    nwg9s.append(wg9h - hbwg9 + hbwg1) 
    nwg10s.append(wg10h - hbwg10 + hbwg1) 
    nwg11s.append(wg11h - hbwg11 + hbwg1) 
    nwg12s.append(wg12h - hbwg12 + hbwg1) 
    nwg13s.append(wg13h - hbwg13 + hbwg1) 
    nwg14s.append(wg14h - hbwg14 + hbwg1)  

    
    t = t + dt
    tts.append(t)
    
  
    print(t)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 

getufromGsplit(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bhbcC = copyarrayfromC(bhbc_c,nbhbc)

nn = len(tts)

s = wdir + "nWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writefile2.writerow(["t","nWG1","nWG2","nWG3","nWG4","nWG5","nWG6","nWG7","nWG8","nWG9","nWG10","nWG11","nWG12","nWG13","nWG14"])  

    for j in range(nn):
        writefile2.writerow([str(tts[j]),str(nwg1s[j]),str(nwg2s[j]), str(nwg3s[j]), str(nwg4s[j]),str(nwg5s[j]),str(nwg6s[j]), str(nwg7s[j]),str(nwg8s[j]),str(nwg9s[j]), str(nwg10s[j]), str(nwg11s[j]),str(nwg12s[j]),str(nwg13s[j]), str(nwg14s[j])])  