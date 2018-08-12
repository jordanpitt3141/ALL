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

def BejiEdgeO(x,hc0,vc0,ft,dx):
    idx = 1.0  /dx
    n = len(x)
    eta = zeros(n)
    v = zeros(n)
    dv = zeros(n)
    hb = 2.46
    bed = -hb*ones(n)
    

    
    i = n-1
    et = 0.5*(ft +(hc0 - hb) )
    h1 = hb + et
    c1 = sqrt(g*(hb + et)) 
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    #linear extrapolation
    i = n - 2
    et = ft
    h1 = hb + et
    c1 = sqrt(g*(hb + et)) 
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    for i in range(n-3,-1,-1):
        et = 2*(eta[i+1]) - (eta[i+2])
        
        h1 = hb + et
        c1 = sqrt(g*(hb + et)) 
        ut = (c1*et) / (h1)
        eta[i] = et
        v[i] = ut
        
    i = -1
    et = 2*(eta[i+1]) - (eta[i+2])
    h1 = hb + et
    c1 = sqrt(g*(hb + et)) 
    ut = (c1*et) / (h1)
    e0 = et
    v0 = ut 
    h0 = hb + e0
    
    hv = hb+ eta
    
    dv[2] = idx*(vc0 -v[1] )
    dv[1] = idx*(v[2] -v[0] )
    dv[0] = idx*(v[1] -v0 )
    
    G = getGfromupy(hv,v,bed,v0,vc0,h0,hc0,-hb,-hb,dx)  
    return hv,v,dv,G

def BejiEdge(x,hc0,vc0,ft):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    hb = 0.4
    
    #SH
    k = 3.06218

    
    i = n-1
    et = ft
    #c1 = sqrt(g*(hb+ et))
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496 #sqrt(9.81*hb) * sqrt(3.0 / (k*k*hb*hb+ 3))
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    #linear extrapolation
    i = n - 2
    et = 2*ft - (hc0 - hb)
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496
    ut = (c1*et) / (h1)
    #c1 = sqrt(g*(hb+ et))
    #ut = (c1*et) / (hb + et)
    eta[i] = et
    v[i] = ut
    
    for i in range(n-3,-1,-1):
        et = 2*(eta[i+1]) - (eta[i+2])
        #c1 = sqrt(g*(hb+ et))
        #ut = (c1*et) / (hb + et)
        
        h1 = hb + et
        c1 = sqrt(g*(hb+ et)) #1.641496
        ut = (c1*et) / (h1)
        eta[i] = et
        v[i] = ut
        
    i = -1
    et = 2*(eta[i+1]) - (eta[i+2])
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496
    ut = (c1*et) / (h1)
    #c1 = sqrt(g*(hb+ et))
    #ut = (c1*et) / (hb + et)
    e0 = et
    v0 = ut 
    h0 = hb + e0
    
    hv = hb+ eta

    dv[2] = idx*(vc0 -v[1] )
    dv[1] = idx*(v[2] -v[0] )
    dv[0] = idx*(v[1] -v0 )
    
    G = getGfromupy(hv,v,bed,v0,vc0,h0,hc0,0,0,dx)  
    return hv,v,dv,D

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    return y1  + (xi)*(y2 - y0)/(x2 - x0) 
    
def MollifyFunc(C,x,e):
    ie = 1.0 / e
    if(abs(ie*x) <1):
        return ie*C*exp(1.0/(abs(ie*x)**2 - 1))
    else:
        return 0

def Roeberflume(x,xexp,bedexp,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    MollPy = zeros(n)
    
    D = 1.0/0.444994
    eta = 10*dx
    rublen = 600
    
    for i in range(n):
        
        MollPy[i] = MollifyFunc(D,x[i] - x[n/2],eta)
        
        if(x[i] <= xexp[0]):
            bed[i] = bedexp[1]
        elif(xexp[0] < x[i] < xexp[-1]):
            j = [ nin for nin, nv in enumerate(xexp) if nv>=x[i] ][0]
            bed[i] = bedexp[j-1] + ( (bedexp[j] - bedexp[j-1]) / 0.05)*(x[i] - xexp[j-1])
            
        elif(x[i] >= xexp[-1]):
            bed[i] = bedexp[-1]

    b0 = bed[0]*ones(rublen)
    b1 = bed[-1]*ones(rublen)
    bedbc = signal.convolve(concatenate((b0,bed,b1)), concatenate((zeros(rublen),MollPy,zeros(rublen))), mode='same') / sum(MollPy)
    bed = bedbc[rublen:-rublen]
        
    for i in range(n):
        h[i] = max(0.0 - bed[i],0)
            
    G = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    
    return h,u,G,bed

  
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

expdir = "/home/jp/Documents/PhD/project/data/Experimental/Roeber/Out/Trial8/"
#wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Roeber/Trial8/FDVM/r1/"  
wdir = "/home/jp/Documents/PhD/project/data/test/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = expdir + "bed.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    xexp = []
    bedexp = []
    for row in readfile:       
            xexp.append(float(row[0]))
            bedexp.append(float(row[1]))
            
s = expdir + "WGs.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    WG2exp = []
    WG3exp = []
    WG4exp = []
    WG5exp = []
    WG6exp = []
    WG7exp = []
    WG8exp = []
    WG9exp = []
    WG10exp = []
    WG11exp = []
    WG12exp = []
    WG13exp = []
    WG14exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG1exp.append(float(row[1]))    
            WG2exp.append(float(row[2]))
            WG3exp.append(float(row[3]))    
            WG4exp.append(float(row[4]))            
            WG5exp.append(float(row[5]))    
            WG6exp.append(float(row[6]))
            WG7exp.append(float(row[7]))    
            WG8exp.append(float(row[8]))
            WG9exp.append(float(row[9]))    
            WG10exp.append(float(row[10]))            
            WG11exp.append(float(row[11]))    
            WG12exp.append(float(row[12]))
            WG13exp.append(float(row[13]))    
            WG14exp.append(float(row[14]))
            
s = expdir + "WGLoc.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    WGloc = []
    for row in readfile:       
            WGloc.append(round(float(row[0]),1))  


g = 9.81
sr = 0.02
dt = sr/ (2**3)
dx = 0.1 / (2**2)
l = 1.0 / (5 + sqrt(g*4))

theta = 1.2
sx = WGloc[0] + dx
ex = 200
st = 0
et = 28

hb = 2.46



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

t = st
        
n = len(x)
hnBC = 3
hnbc = 3*n + 2*hnBC
unBC = 3
unbc = 2*n + 1 + 2*(unBC -1)
CnBC = 1
Cnbc = n + 2*CnBC 

   
niBC = 4

xhMbeg =array( [x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx])
xhMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
  
xCbeg =array( [x[0] - dx])
xCend = array([x[-1] + dx])

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
  

theta = 1.2

tij = 0
ij = 0
h,u,G,b = Roeberflume(x,xexp,bedexp,dx)

hMend = h[-1]*ones(hnBC)
uMend = u[-1]*ones(unBC)
GMend = G[-1]*ones(hnBC)
bMbeg = b[0]*ones(hnBC)
bMend = b[-1]*ones(hnBC)
wMend = h[-1]*ones(hnBC) + bMend


hMbeg,uMbeg,duMbeg,GMbeg = BejiEdge(xhMbeg,h[0],u[0],0,dx)
wMbeg = hMbeg +bMbeg


dbMbeg = zeros(hnBC)
dbMend = zeros(hnBC)
 
duMend = zeros(unBC)
   
ddbCbeg = zeros(CnBC)
ddbCend = zeros(CnBC)

ftcs.append(0)

ct = dt
mp = int(ct/sr)
ftc = lineinterp(WG1exp[mp] - hb ,WG1exp[mp + 1] - hb,texp[mp],texp[mp + 1],ct - texp[mp])

hMbegdt,uMbegdt,duMbegdt,GMbegdt =  BejiEdge(xhMbeg,h[0],u[0],ftc,dx) 
wMbegdt = hMbegdt +bMbeg

print(ct,ftc,mp, WG1exp[mp], WG1exp[mp + 1])
print(hMbegdt)
print(uMbegdt)
print(duMbegdt)
print(GMbegdt)
 
ftcs.append(ftc)


uMbegdt_c = copyarraytoC(uMbegdt)
duMbegdt_c = copyarraytoC(duMbegdt)
hMbegdt_c = copyarraytoC(hMbegdt)
wMbegdt_c = copyarraytoC(wMbegdt)
GMbegdt_c = copyarraytoC(GMbegdt)


uMbeg_c = copyarraytoC(uMbeg)
hMbeg_c = copyarraytoC(hMbeg)
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg)

dbMbeg_c = copyarraytoC(dbMbeg)
ddbCbeg_c = copyarraytoC(ddbCbeg)
duMbeg_c = copyarraytoC(duMbeg)

uMend_c = copyarraytoC(uMend)
hMend_c = copyarraytoC(hMend)
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend)

dbMend_c = copyarraytoC(dbMend)
ddbCend_c = copyarraytoC(ddbCend)
duMend_c = copyarraytoC(duMend)

bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)

h_c = copyarraytoC(h)
b_c = copyarraytoC(b)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)

hbc_c =  mallocPy(hnbc)
bMbc_c =  mallocPy(hnbc)
dbMbc_c =  mallocPy(hnbc)
ddbCbc_c =  mallocPy(Cnbc)
wbc_c =  mallocPy(hnbc)
ubc_c =  mallocPy(unbc)
uFDbc_c =  mallocPy(unbc)
duFDbc_c =  mallocPy(unbc)
Gbc_c =  mallocPy(hnbc)

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

   




tts.append(t)
#Just an FEM solve here
while t < et:  
    
    #print(hMbegdt,uMbegdt,duMbegdt,GMbegdt )
    
    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,hMbegdt_c,hMend_c,GMbegdt_c,GMend_c,wMbegdt_c,wMend_c,uMbegdt_c,uMend_c,duMbegdt_c,duMend_c,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g)
    
    #ReconandSolve(h_c,G_c,b_c,hMbegdt_c,hMend_c,GMbegdt_c,GMend_c,wMbegdt_c,wMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbegdt_c,uMend_c,duMbegdt_c,duMend_c,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bMbc_c,dbMbc_c,ddbCbc_c, uFDbc_c,duFDbc_c)    
    
    uc0 = readfrommem(uFDbc_c,unBC)
    hc0 = readfrommem(h_c,0) 
    
    if(isnan(hc0)) :
        break
    
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
          

    nwg1s.append((readfrommem(hMbegdt_c,2)  - hbwg1) + hbwg1)
    
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

    copywritearraytoC(hMbegdt,hMbeg_c)
    copywritearraytoC(GMbegdt,GMbeg_c)
    copywritearraytoC(uMbegdt,uMbeg_c)
    copywritearraytoC(duMbegdt,duMbeg_c)
    copywritearraytoC(wMbegdt,wMbeg_c) 
    
    t = t + dt
    tts.append(t)
    
    ct = t + dt
    mp = int(ct/sr)
    ftc = lineinterp(WG1exp[mp] - hb ,WG1exp[mp + 1] - hb,texp[mp],texp[mp + 1],ct - texp[mp])
         
    hMbegdt,uMbegdt,duMbegdt,GMbegdt =  BejiEdge(xhMbeg,hc0,uc0,ftc,dx) 
    wMbegdt = hMbegdt + bMbeg   
    
    #print(ct,ftc,mp, WG1exp[mp], WG1exp[mp + 1])
    #print(hMbegdt)
    #print(uMbegdt)
    #print(duMbegdt)
    #print(GMbegdt)
    
    copywritearraytoC(hMbegdt,hMbegdt_c)
    copywritearraytoC(wMbegdt,wMbegdt_c)
    copywritearraytoC(GMbegdt,GMbegdt_c)
    copywritearraytoC(uMbegdt,uMbegdt_c)
    copywritearraytoC(duMbegdt,duMbegdt_c)

    print(t)



ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bMbc_c,dbMbc_c,ddbCbc_c, uFDbc_c,duFDbc_c)    

wbcC = copyarrayfromC(wbc_c,hnbc)  
bMbcC = copyarrayfromC(bMbc_c,hnbc)
dbMbcC = copyarrayfromC(dbMbc_c,hnbc)
ddbCbcC = copyarrayfromC(ddbCbc_c,Cnbc)
hbcC = copyarrayfromC(hbc_c,hnbc)  
ubcC = copyarrayfromC(ubc_c,unbc)  
uFDbcC = copyarrayfromC(uFDbc_c,unbc) 
duFDbcC = copyarrayfromC(duFDbc_c,unbc) 
GbcC = copyarrayfromC(Gbc_c,hnbc)  

hiC = copyarrayfromC(h_c,n)
GiC = copyarrayfromC(G_c,n) 
uiC = uFDbcC[unBC:-unBC:2]
wiC = hiC + b




deallocPy(h_c)
deallocPy(G_c)
deallocPy(b_c)
deallocPy(x_c)

deallocPy(wbc_c)
deallocPy(bMbc_c)
deallocPy(dbMbc_c)
deallocPy(ddbCbc_c)
deallocPy(hbc_c)
deallocPy(ubc_c)
deallocPy(uFDbc_c)
deallocPy(duFDbc_c)
deallocPy(Gbc_c)



deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(wMbeg_c)
deallocPy(uMbeg_c)
deallocPy(dbMbeg_c)
deallocPy(ddbCbeg_c)
deallocPy(duMbeg_c)

deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(wMend_c)
deallocPy(uMend_c)
deallocPy(dbMend_c)
deallocPy(ddbCend_c)
deallocPy(duMend_c)
