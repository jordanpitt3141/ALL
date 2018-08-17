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


def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    G = zeros(n)
    h = ones(n)*h1
    bed = zeros(n)
    for i in range(n):
        if (x[i] <0 - dx and x[i] > -2*b - dx):
            h[i] = h0
    
    w = h + bed

    return h,u,G,bed,w
  
def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    return y1  + (xi)*(y2 - y0)/(x2 - x0) 
    
def find_nearestidx(array1, value1):
    array2 = asarray(array1)
    idx = (abs(array2 - value1)).argmin()
    return idx

#Forcing Problem    
wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Segur/3cm/FDVM/"
#wdir = "/home/jp/Documents/PhD/project/data/Tests/Syn/FDVM/"
if not os.path.exists(wdir):
    os.makedirs(wdir)

WGLocdivh = [0,50,100,150,200]
WGloc = [0,5,10,15,20]


g = 9.81

tl = 60.0
blength = 0.61
h0 = 0.07
h1 = 0.1
g = 9.81
theta = 1.2

startx = -tl
sx = startx
endx = tl
ex = endx

startt = 0.0
st = startt
endt = 50
et = endt
dx = 0.01
Cr = 0.5
l = Cr / sqrt(g*h1)
dt = l*dx

t = st
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
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
  

theta = 1.2
 

h,u,G,b,w = experiment1(x,blength,h0,h1,dx)

wg1i = find_nearestidx(x, WGloc[0])
wg2i = find_nearestidx(x, WGloc[1])
wg3i = find_nearestidx(x, WGloc[2])
wg4i = find_nearestidx(x, WGloc[3])
wg5i = find_nearestidx(x, WGloc[4])

tts = []
nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []


hMbeg = h1*ones(hnBC)
hMend = h1*ones(hnBC)

wMbeg = h1*ones(hnBC)
wMend = h1*ones(hnBC)

bMbeg = zeros(hnBC)
bMend = zeros(hnBC)

GMbeg = zeros(hnBC)
GMend = zeros(hnBC)

uMbeg = zeros(unBC)
uMend = zeros(unBC)

dbMbeg = zeros(hnBC)
dbMend = zeros(hnBC)
 
duMbeg = zeros(unBC)
duMend = zeros(unBC)
   
ddbCbeg = zeros(CnBC)
ddbCend = zeros(CnBC)


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
   
xbegC = arange(sx - niBC*dx,sx,dx)
xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 

u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)
b0C = b[0]*ones(niBC)
b1C = b[-1]*ones(niBC)

xbcC =  concatenate([xbegC,x,xendC])
bbcC =  concatenate([b0C,b,b1C])
xbcC_c = copyarraytoC(xbcC)
bbcC_c = copyarraytoC(bbcC)


hbcC =  concatenate([h0C,h,h1C])
ubcC =  concatenate([u0C,u,u1C])
GbcC =  concatenate([G0C,G,G1C])

hbcC_c = copyarraytoC(hbcC)
ubcC_c = copyarraytoC(ubcC)
GbcC_c = copyarraytoC(GbcC)

Eni = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
Pni = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
Mni = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
Gni = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)

deallocPy(hbcC_c)
deallocPy(ubcC_c)
deallocPy(GbcC_c)



tts.append(t)
nwg1s.append(h1)    
nwg2s.append(h1)  
nwg3s.append(h1) 
nwg4s.append(h1) 
nwg5s.append(h1) 
#Just an FEM solve here
while t < endt:  
            
    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g);      
    t = t + dt
    
    wg1h = CELLRECON(readfrommem(h_c,wg1i -1),readfrommem(h_c,wg1i) ,readfrommem(h_c,wg1i + 1),x[wg1i -1],x[wg1i],x[wg1i + 1],WGloc[0] - x[wg1i])
    wg2h = CELLRECON(readfrommem(h_c,wg2i -1),readfrommem(h_c,wg2i) ,readfrommem(h_c,wg2i + 1),x[wg2i -1],x[wg2i],x[wg2i + 1],WGloc[1] - x[wg2i])
    wg3h = CELLRECON(readfrommem(h_c,wg3i -1),readfrommem(h_c,wg3i) ,readfrommem(h_c,wg3i + 1),x[wg3i -1],x[wg3i],x[wg3i + 1],WGloc[2] - x[wg3i])
    wg4h = CELLRECON(readfrommem(h_c,wg4i -1),readfrommem(h_c,wg4i) ,readfrommem(h_c,wg4i + 1),x[wg4i -1],x[wg4i],x[wg4i + 1],WGloc[3] - x[wg4i])
    wg5h = CELLRECON(readfrommem(h_c,wg5i -1),readfrommem(h_c,wg5i) ,readfrommem(h_c,wg5i + 1),x[wg5i -1],x[wg5i],x[wg5i + 1],WGloc[4] - x[wg5i])
    
    
    nwg1s.append(wg1h)    
    nwg2s.append(wg2h)  
    nwg3s.append(wg3h) 
    nwg4s.append(wg4h) 
    nwg5s.append(wg5h) 
    tts.append(t)
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





u0C = uiC[0]*ones(niBC)
u1C = uiC[-1]*ones(niBC)   
h0C = hiC[0]*ones(niBC)
h1C = hiC[-1]*ones(niBC)
G0C = GiC[0]*ones(niBC)
G1C = GiC[-1]*ones(niBC)

hbcC1 =  concatenate([h0C,hiC,h1C])
ubcC1 =  concatenate([u0C,uiC,u1C])
GbcC1 =  concatenate([G0C,GiC,G1C])

hbcC_c = copyarraytoC(hbcC1)
ubcC_c = copyarraytoC(ubcC1)
GbcC_c = copyarraytoC(GbcC1)

En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)

deallocPy(hbcC_c)
deallocPy(ubcC_c)
deallocPy(GbcC_c)


s = wdir +  "outListLast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
    n = len(x)         
    for j in range(n):
        writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
        
s = wdir +  "outSingLast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G" ,"Energyi", "Massi", "Momentumi", "Gi"  ])   
    writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni)  ]) 
    
  
nt = len(tts)

s = wdir +  "outWGList.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["t(s)","WG1s","WG2s" ,"WG3s" ,"WG4s" ,"WG5s"  ])        
    n = len(x)         
    for j in range(nt):
        writefile2.writerow([str(tts[j]),str(nwg1s[j]),str(nwg2s[j]),str(nwg3s[j]),str(nwg4s[j]),str(nwg5s[j])])
        
s = wdir + "wG1Num.dat"
with open(s,'a') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[0]/h1 ," ",1.5*((nwg1s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG2Num.dat"
with open(s,'a') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[1]/h1 ," ",1.5*((nwg2s[j] - h1)/h1))
        file1.write(s)
        
s = wdir + "wG3Num.dat"
with open(s,'a') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[2]/h1 ," ",1.5*((nwg3s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG4Num.dat"
with open(s,'a') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[3]/h1 ," ",1.5*((nwg4s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG5Num.dat"
with open(s,'a') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[4]/h1 ," ",1.5*((nwg5s[j] - h1)/h1))
        file1.write(s)

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
