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
from scipy import signal

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
    eta = 3*dx
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
    

def BejiEdge(x,hc0,vc0,ft):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    dv = zeros(n)
    hb = 2.46
    

    
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
    


def minmod(a,b,c):
    
    if(a>0 and b> 0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0

def ReconEdge(qj,qjm1,qjp1,theta,idx):
    
    dqjf = idx*(qjp1 - qj) ;
    dqjb = idx*(qj - qjm1) ;
    dqjm = 0.5*idx*(qjp1 - qjm1) ;
    dqjlim = minmod(theta*dqjf,dqjm,theta*dqjb);
    
    qjmhp = qj - 0.5*dx*dqjlim ;
    qjphm = qj + 0.5*dx*dqjlim ;
    
    return qjmhp,qjphm
    
    
    

   
def BedSlopes(b,h,x,theta,dx):
    idx = 1.0 / dx
    n = len(x)
    bxe = []
    xbxe = []
    bxm = []
    xbxm = []
    bxxm = []
    xbxxm = []
    
    for j in range(3,n-3):
        wjmhp , wjphm = ReconEdge(h[j] + b[j],h[j-1] + b[j-1],h[j+1] + b[j+1],theta,idx)
        hjmhp , hjphm = ReconEdge(h[j],h[j-1],h[j+1],theta,idx)
        bjmhp = wjmhp - hjmhp;
        bjphm = wjphm - hjphm;
    
        wjm3hp , wjmhm = ReconEdge(h[j-1] + b[j-1],h[j-2] + b[j-2],h[j] + b[j],theta,idx)
        hjm3hp , hjmhm = ReconEdge(h[j-1],h[j-2],h[j],theta,idx)
        bjm3hp = wjm3hp - hjm3hp;
        bjmhm = wjmhm - hjmhm;   
        
        wjphp , wjp3hm = ReconEdge(h[j+1] + b[j+1],h[j] + b[j],h[j+2] + b[j+2],theta,idx)
        hjphp , hjp3hm = ReconEdge(h[j+1],h[j],h[j+2],theta,idx)
        bjphp = wjphp - hjphp;
        bjp3hm = wjp3hm - hjp3hm;
        
        wjp3hp , wjp5hm = ReconEdge(h[j+2] + b[j+2],h[j+1] + b[j+1],h[j+3] + b[j+3],theta,idx)
        hjp3hp , hjp5hm = ReconEdge(h[j+2],h[j+ 1],h[j+3],theta,idx)
        bjp3hp = wjp3hp - hjp3hp;
        bjp5hm = wjp5hm - hjp5hm;
        
        dbjphp = idx*(bjp3hp - bjphp)
        dbjphm = idx*(bjphm - bjmhm)
        
        bxe.append(dbjphm)
        bxe.append(dbjphp)
        xbxe.append(x[j] + 0.5*dx)
        xbxe.append(x[j] + 0.5*dx)
        
        dbj = idx*(bjphm - bjmhp)
        bxm.append(dbj)
        xbxm.append(x[j])
        
        bxxm.append((b[j+1] - 2*b[j] + b[j-1])*idx**2)
        xbxxm.append(x[j])
        
    
    
    return bxe,xbxe,bxm,xbxm,bxxm,xbxxm


expdir = "/home/jp/Documents/PhD/project/data/Experimental/Roeber/Out/Trial8/"
#wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Roeber/Trial8/FDVM/r3/"  
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
sx = WGloc[0] + 0.5*dx
ex = 100
st = 20
et = 33.55

hb = 2.46

nMBC = 3
nEBC = 3
nCBC = 1

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


nMbc = 3*n + 2*nMBC
nEbc = 2*n - 1 + 2*nEBC
nCbc = n + 2*nCBC

xMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

tts = [] 
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

h,u,G,b = Roeberflume(x,xexp,bedexp,dx)

bxe,xbxe,bxm,xbxm,bxxm,xbxxm = BedSlopes(b,h,x,theta,dx)

uEend = u[-1]*ones(nEBC)
duEend = zeros(nEBC)

hMend = h[-1]*ones(nMBC)
GMend = G[-1]*ones(nMBC)
wMend = h[-1]*ones(nMBC) + b[-1]*ones(nMBC)
bMbeg = b[0]*ones(nMBC)
bMend = b[-1]*ones(nMBC)
ddbCbeg = zeros(nCBC)
ddbCend = zeros(nCBC)

ct = st
hMbeg,uEbeg,GMbeg,duEbeg = BejiEdge(xMbeg,h[0],u[0],0)
wMbeg = hMbeg + bMbeg

ct =st +  dt
mp = int(ct/sr)
ftc = lineinterp(WG1exp[mp] - hb ,WG1exp[mp + 1] - hb,texp[mp],texp[mp + 1],ct - texp[mp])

hMbegdt,uEbegdt,GMbegdt,duEbegdt =  BejiEdge(xMbeg,h[0],u[0],ftc) 
wMbegdt = hMbegdt + bMbeg
         
         
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


t = st
tts.append(t)
#Just an FEM solve here
while t < et: 
    
    evolvewrapBC(h_c,G_c,b_c,hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c,ddbCend_c,hMbegdt_c,GMbegdt_c,wMbegdt_c,duEbegdt_c,uEbegdt_c,hMend_c,GMend_c,wMend_c,duEend_c,uEend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta)
   
    edgevaluesSplit(h_c,G_c,b_c, hMbegdt_c,GMbegdt_c,wMbegdt_c,bMbeg_c,duEbegdt_c,uEbegdt_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
        
    uc0 = readfrommem(uEbc_c,nEBC)
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
    copywritearraytoC(uEbegdt,uEbeg_c)
    copywritearraytoC(duEbegdt,duEbeg_c)
    copywritearraytoC(wMbegdt,wMbeg_c)    
    
    t = t + dt
    tts.append(t)
    
    ct = t + dt
    mp = int(ct/sr)
    ftc = lineinterp(WG1exp[mp] - hb ,WG1exp[mp + 1] - hb,texp[mp],texp[mp + 1],ct - texp[mp])
    
    hMbegdt,uEbegdt,GMbegdt,duEbegdt =  BejiEdge(xMbeg,hc0,uc0,ftc) 
    
    wMbegdt = hMbegdt + bMbeg
    
    copywritearraytoC(hMbegdt,hMbegdt_c)
    copywritearraytoC(wMbegdt,wMbegdt_c)
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

s = wdir + "nWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writefile2.writerow(["t","nWG1","nWG2","nWG3","nWG4","nWG5","nWG6","nWG7","nWG8","nWG9","nWG10","nWG11","nWG12","nWG13","nWG14"])  

    for j in range(nn):
        writefile2.writerow([str(tts[j]),str(nwg1s[j]),str(nwg2s[j]), str(nwg3s[j]), str(nwg4s[j]),str(nwg5s[j]),str(nwg6s[j]), str(nwg7s[j]),str(nwg8s[j]),str(nwg9s[j]), str(nwg10s[j]), str(nwg11s[j]),str(nwg12s[j]),str(nwg13s[j]), str(nwg14s[j])])  

