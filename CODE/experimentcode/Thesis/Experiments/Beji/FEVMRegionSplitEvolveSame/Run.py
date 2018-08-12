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



#Beji Problem    
exp = "sl"
expdir = "/home/jp/Documents/PhD/project/data/Experimental/Beji/Out/"
wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Beji/"+str(exp)+"/FEVM/r1/"  

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
dt = sr/ 2**5 
dx = (0.1/2.0**4) 

theta = 1.2
sx = wg1loc + 0.5*dx
ex = 150
st = 0
et = 60

hb = 0.4


GhnBC = 3
unBC = 3
bnBC = 4

x = arange(sx,ex +0.1*dx, dx)

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

hMend = h[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)
GMend = G[-1]*ones(GhnBC)
wMend = h[-1]*ones(GhnBC) + b[-1]*ones(GhnBC)
bMbeg = zeros(bnBC)
bMend = b[-1]*ones(bnBC)




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


hMbeg,uMbeg,GMbeg = BejiEdge(xhuMbeg,h[0],u[0],0)
wMbeg = hMbeg

ftcs.append(0)

ct = dt
mp = int(ct/sr)
ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp])

hMbegdt,uMbegdt,GMbegdt =  BejiEdge(xhuMbeg,h[0],u[0],ftc) 
wMbegdt = hMbegdt

ftcs.append(ftc)

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
   
hMbegdt_c = copyarraytoC(hMbegdt)
wMbegdt_c = copyarraytoC(wMbegdt) 
uMbegdt_c = copyarraytoC(uMbegdt)  
GMbegdt_c = copyarraytoC(GMbegdt)
   

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)

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
    
    evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbegdt_c,hMend_c,wMbegdt_c,wMend_c,GMbegdt_c,GMend_c,uMbegdt_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t,0,0,0,0,0,0,0,0,0,0)
    
    getufromGsplit(h_c, G_c, b_c,hMbegdt_c,hMend_c,GMbegdt_c,GMend_c,uMbegdt_c,uMend_c,wMbegdt_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
        
    uc0 = readfrommem(ubc_c,unBC)
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
    copywritearraytoC(uMbegdt,uMbeg_c)
    copywritearraytoC(hMbegdt,wMbeg_c)
    
    t = t + dt
    tts.append(t)
    
    ct = t + dt
    mp = int(ct/sr)
    ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct -ts[mp])
    if(t <= et ):
        ftcs.append(ftc)
    
    
    hMbegdt,uMbegdt,GMbegdt =  BejiEdge(xhuMbeg,hc0,uc0,ftc) 
    copywritearraytoC(hMbegdt,hMbegdt_c)
    copywritearraytoC(hMbegdt,wMbegdt_c)
    copywritearraytoC(GMbegdt,GMbegdt_c)
    copywritearraytoC(uMbegdt,uMbegdt_c)
  
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


deallocPy(h_c)
deallocPy(G_c)
deallocPy(u_c)

deallocPy(ubc_c)
deallocPy(hhbc_c)
deallocPy(whbc_c)
deallocPy(Ghbc_c)
deallocPy(bhbc_c)

deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
