# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe

from scipy.optimize import bisect

#gives exact up to linears, so is second order accurate huzzah    
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

    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    

def DingFlume(x,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):

        if(0 <= x[i] < 6):
            bed[i] = 0.0
            h[i] = 0.4            
        elif(6 <= x[i] <= 12):
            bed[i] = 0.05*(x[i] - 6)
            h[i] = 0.4 - bed[i]
        elif(12 < x[i] <= 14):
            bed[i] = 0.3
            h[i] = 0.1
        elif(14 < x[i] <= 17):
            bed[i] = 0.3 - 0.1*(x[i] - 14)
            h[i] = 0.4 - bed[i]
        elif(17 < x[i] <= 18.95):
            bed[i] = 0.0
            h[i] = 0.4 - bed[i]
        else:
            bed[i] = 0.0
            h[i] = 0.4  - bed[i]

    G = getGfromupy(h,u,bed,0,0,0.4,h[-1],0,bed[-1],dx)
    
    return h,u,G,bed

def BejiEdge(x,hc0,vc0,ft):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    dv = zeros(n)
    hb = 0.4
    

    
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
    
    dv[-1] = (v[-2] - vc0)/dx
    dv[-2] = (v[-1] - v[0])/dx
    dv[-3] = (v[-2] - v0)/dx

    
    G = getGfromupy(hv,v,bed,v0,vc0,h0,hc0,0,0,dx)  
    return hv,v,G,dv

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


#Beji Problem    
exp = "sh"
expdir = "../../../../../../../data/Experimental/Data 1994 Paper/CSV/"
wdir = "../../../../../../../data/raw/Thesis/Experiment/Beji/"+str(exp)+"/FDVM/r2/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

wg1loc = 5.7
wg2loc = 10.5
wg3loc = 12.5
wg4loc = 13.5
wg5loc = 14.5
wg6loc = 15.7
wg7loc = 17.3


a0 = 0.01
k = 0.1

### WAVE LENGTH

g = 9.81
sr = 0.039312
#dt = sr/ 2**5
#dx = (0.1/2.0**4)

dt = sr/ 2**5 
dx = (0.1/2.0**4) 

theta = 1.2
sx = wg1loc + 0.5*dx
ex = 150
st = 0
et = 60

hb = 0.4


nMBC = 3
nEBC = 3
nCBC = 1

x = arange(sx,ex +0.1*dx, dx)

n = len(x) 

nMbc = 3*n + 2*nMBC
nEbc = 2*n - 1 + 2*nEBC
nCbc = n + 2*nCBC




xMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])


tts = [] 
ftcs = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nwg8s = []
nwg1s = []


tij = 0.0
ij = 0
h,u,G,b = DingFlume(x,dx)

uEend = u[-1]*ones(nEBC)
duEend = zeros(nEBC)

hMend = h[-1]*ones(nMBC)
GMend = G[-1]*ones(nMBC)
wMend = h[-1]*ones(nMBC) + b[-1]*ones(nMBC)
bMbeg = zeros(nMBC)
bMend = b[-1]*ones(nMBC)
ddbCbeg = zeros(nCBC)
ddbCend = zeros(nCBC)


wg2i =  find_nearestidx(x, wg2loc)
wg3i = find_nearestidx(x, wg3loc)
wg4i = find_nearestidx(x, wg4loc)
wg5i = find_nearestidx(x, wg5loc)
wg6i = find_nearestidx(x, wg6loc)
wg7i = find_nearestidx(x, wg7loc)

s = expdir + exp + ".csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    ts = [0.0]
    rs = [0.0]
    wg1s = [0.0]
    wg2s = [0.0]
    wg3s = [0.0]
    wg4s = [0.0]
    wg5s = [0.0]
    wg6s = [0.0]
    wg7s = [0.0]

    j = -1
    for row in readfile:   
        if (j >= 0):
            ts.append((j + 1)*sr)
            rs.append(float(row[0]))
            wg1s.append(float(row[1]))
            wg2s.append(float(row[2]))
            wg3s.append(float(row[3]))
            wg4s.append(float(row[4]))
            wg5s.append(float(row[5]))
            wg6s.append(float(row[6]))
            wg7s.append(float(row[7]))
        j = j + 1


hMbeg,uEbeg,GMbeg,duEbeg = BejiEdge(xMbeg,h[0],u[0],0)
wMbeg = hMbeg

ftcs.append(0)

ct = dt
mp = int(ct/sr)
ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp])

hMbegdt,uEbegdt,GMbegdt,duEbegdt =  BejiEdge(xMbeg,h[0],u[0],ftc) 
wMbegdt = hMbegdt

ftcs.append(ftc)

h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)
b_c = copyarraytoC(b)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend)
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)
uEbeg_c = copyarraytoC(uEbeg)
uEend_c = copyarraytoC(uEend)
duEbeg_c = copyarraytoC(duEbeg)
duEend_c = copyarraytoC(duEend)
ddbCbeg_c = copyarraytoC(ddbCbeg)
ddbCend_c = copyarraytoC(ddbCend)
   
hMbegdt_c = copyarraytoC(hMbegdt)
wMbegdt_c = copyarraytoC(wMbegdt) 
uEbegdt_c = copyarraytoC(uEbegdt)  
duEbegdt_c = copyarraytoC(duEbegdt) 
GMbegdt_c = copyarraytoC(GMbegdt)
   

duEbc_c = mallocPy(nEbc)
uEbc_c = mallocPy(nEbc)
hMbc_c = mallocPy(nMbc)
wMbc_c = mallocPy(nMbc)
GMbc_c = mallocPy(nMbc)
bMbc_c = mallocPy(nMbc)
ddbCbc_c = mallocPy(nCbc)

nwg1s.append(0)
nwg2s.append(0)  
nwg3s.append(0)  
nwg4s.append(0)  
nwg5s.append(0)  
nwg6s.append(0)  
nwg7s.append(0)  

hbwg1 = hMbeg[1]

hbwg2 = CELLRECON(readfrommem(h_c,wg2i -1),readfrommem(h_c,wg2i) ,readfrommem(h_c,wg2i + 1),x[wg2i -1],x[wg2i],x[wg2i + 1],wg2loc - x[wg2i])
hbwg3 = CELLRECON(readfrommem(h_c,wg3i -1),readfrommem(h_c,wg3i) ,readfrommem(h_c,wg3i + 1),x[wg3i -1],x[wg3i],x[wg3i + 1],wg3loc - x[wg3i])
hbwg4 = CELLRECON(readfrommem(h_c,wg4i -1),readfrommem(h_c,wg4i) ,readfrommem(h_c,wg4i + 1),x[wg4i -1],x[wg4i],x[wg4i + 1],wg4loc - x[wg4i])
hbwg5 = CELLRECON(readfrommem(h_c,wg5i -1),readfrommem(h_c,wg5i) ,readfrommem(h_c,wg5i + 1),x[wg5i -1],x[wg5i],x[wg5i + 1],wg5loc - x[wg5i])
hbwg6 = CELLRECON(readfrommem(h_c,wg6i -1),readfrommem(h_c,wg6i) ,readfrommem(h_c,wg6i + 1),x[wg6i -1],x[wg6i],x[wg6i + 1],wg6loc - x[wg6i])
hbwg7 = CELLRECON(readfrommem(h_c,wg7i -1),readfrommem(h_c,wg7i) ,readfrommem(h_c,wg7i + 1),x[wg7i -1],x[wg7i],x[wg7i + 1],wg7loc - x[wg7i])


t = 0.0
tts.append(t)
#Just an FEM solve here
while t < et: 
    
    evolvewrapBC(h_c,G_c,b_c,hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c,ddbCend_c,hMbegdt_c,GMbegdt_c,wMbegdt_c,duEbegdt_c,uEbegdt_c,hMend_c,GMend_c,wMend_c,duEend_c,uEend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta)
   
    edgevaluesSplit(h_c,G_c,b_c, hMbegdt_c,GMbegdt_c,wMbegdt_c,bMbeg_c,duEbegdt_c,uEbegdt_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
        
    uc0 = readfrommem(uEbc_c,nEBC)
    hc0 = readfrommem(h_c,0) 
    
    wg2h = CELLRECON(readfrommem(h_c,wg2i -1),readfrommem(h_c,wg2i) ,readfrommem(h_c,wg2i + 1),x[wg2i -1],x[wg2i],x[wg2i + 1],wg2loc - x[wg2i])
    wg3h = CELLRECON(readfrommem(h_c,wg3i -1),readfrommem(h_c,wg3i) ,readfrommem(h_c,wg3i + 1),x[wg3i -1],x[wg3i],x[wg3i + 1],wg3loc - x[wg3i])
    wg4h = CELLRECON(readfrommem(h_c,wg4i -1),readfrommem(h_c,wg4i) ,readfrommem(h_c,wg4i + 1),x[wg4i -1],x[wg4i],x[wg4i + 1],wg4loc - x[wg4i])
    wg5h = CELLRECON(readfrommem(h_c,wg5i -1),readfrommem(h_c,wg5i) ,readfrommem(h_c,wg5i + 1),x[wg5i -1],x[wg5i],x[wg5i + 1],wg5loc - x[wg5i])
    wg6h = CELLRECON(readfrommem(h_c,wg6i -1),readfrommem(h_c,wg6i) ,readfrommem(h_c,wg6i + 1),x[wg6i -1],x[wg6i],x[wg6i + 1],wg6loc - x[wg6i])
    wg7h = CELLRECON(readfrommem(h_c,wg7i -1),readfrommem(h_c,wg7i) ,readfrommem(h_c,wg7i + 1),x[wg7i -1],x[wg7i],x[wg7i + 1],wg7loc - x[wg7i])
    
  
    
    nwg1s.append((readfrommem(hMbegdt_c,2)  - hbwg1))
    nwg2s.append((wg2h - hbwg2))  
    nwg3s.append((wg3h - hbwg3)) 
    nwg4s.append((wg4h - hbwg4)) 
    nwg5s.append((wg5h - hbwg5)) 
    nwg6s.append((wg6h - hbwg6)) 
    nwg7s.append(wg7h - hbwg7) 
    
    copywritearraytoC(hMbegdt,hMbeg_c)
    copywritearraytoC(GMbegdt,GMbeg_c)
    copywritearraytoC(uEbegdt,uEbeg_c)
    copywritearraytoC(duEbegdt,duEbeg_c)
    copywritearraytoC(hMbegdt,wMbeg_c)
    
    t = t + dt
    tts.append(t)
    
    ct = t + dt
    mp = int(ct/sr)
    ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct -ts[mp])
    if(t <= et ):
        ftcs.append(ftc)
    
    
    hMbegdt,uEbegdt,GMbegdt,duEbegdt =  BejiEdge(xMbeg,hc0,uc0,ftc) 
    copywritearraytoC(hMbegdt,hMbegdt_c)
    copywritearraytoC(hMbegdt,wMbegdt_c)
    copywritearraytoC(GMbegdt,GMbegdt_c)
    copywritearraytoC(uEbegdt,uEbegdt_c)
    copywritearraytoC(duEbegdt,duEbegdt_c)
  
    print(t)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 
edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
uEbcC = copyarrayfromC(uEbc_c,nEbc)
uC = uEbcC[nEBC:-nEBC:2]


nn = len(tts)
ne = len(ts)
s = wdir + "eWG1.dat"
with open(s,'w') as file1:
    for j in range(ne):
        ss ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg1s[j]/100.0)
        file1.write(ss)

s = wdir + "eWG2.dat"
with open(s,'w') as file1:
    for j in range(ne):
        ss ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg2s[j]/100.0)
        file1.write(ss)

s = wdir + "eWG3.dat"
with open(s,'w') as file1:
    for j in range(ne):
        s ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg3s[j]/100.0)
        file1.write(s)

s = wdir + "eWG4.dat"
with open(s,'w') as file1:
    for j in range(ne):
        s ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg4s[j]/100.0)
        file1.write(s)

s = wdir + "eWG5.dat"
with open(s,'w') as file1:
    for j in range(ne):
        s ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg5s[j]/100.0)
        file1.write(s)

s = wdir + "eWG6.dat"
with open(s,'w') as file1:
    for j in range(ne):
        s ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg6s[j]/100.0)
        file1.write(s)

s = wdir + "eWG7.dat"
with open(s,'w') as file1:
    for j in range(ne):
        s ="%3.8f%5s%1.15f\n" %(ts[j]," ",wg7s[j]/100.0)
        file1.write(s)


s = wdir + "nWG1.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg1s[j])
        file1.write(ss)


s = wdir + "nWG2.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg2s[j])
        file1.write(ss)

s = wdir + "nWG3.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg3s[j])
        file1.write(ss)
        
s = wdir + "nWG4.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg4s[j])
        file1.write(ss)

s = wdir + "nWG5.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg5s[j])
        file1.write(ss)

s = wdir + "nWG6.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg6s[j])
        file1.write(ss)

s = wdir + "nWG7.dat"
with open(s,'w') as file1:
    for j in range(nn):
        ss ="%3.8f%5s%1.15f\n" %(tts[j]," ",nwg7s[j])
        file1.write(ss)