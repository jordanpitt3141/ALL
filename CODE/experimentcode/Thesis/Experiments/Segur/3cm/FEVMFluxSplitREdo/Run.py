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
    
def HamilA(g,a0,a1,x0,x1,x2,x3):
    Hamil = 0.5*g*(a0**2*(x1 - x0) + a1**2*(x2 - x1) + a0**2*(x3 - x2))
    return Hamil

def MassA(g,a0,a1,x0,x1,x2,x3):
    Hamil = (a0*(x1 - x0) + a1*(x2 - x1) + a0*(x3 - x2))
    return Hamil

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
        if (x[i] <0  and x[i] > -2*b):
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
wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Segur/3cm/FEVM/"

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
endt = 0.5
et = endt
dx = 0.01
Cr = 0.5
l = Cr / sqrt(g*h1)
dt = l*dx


t = startt
x = arange(startx,endx +0.1*dx, dx)
n = len(x)

hnBC = 3
hnbc = 3*n + 2*hnBC
bnMBC = 7
bnBC = 4
bnbc = 3*n + 1 + 2*(bnBC -1)
unBC = 3
unbc = 2*n + 1 + 2*(unBC -1)



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




niBC = 4
xbegC = arange(sx - niBC*dx,sx,dx)
xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 


b0C = b[0]*ones(niBC)
b1C = b[-1]*ones(niBC)


xbcC1 =  concatenate([xbegC,x,xendC])
bbcC1 =  concatenate([b0C,b,b1C])
xbcC1_c = copyarraytoC(xbcC1)
bbcC1_c = copyarraytoC(bbcC1)

u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)


hbcC1 =  concatenate([h0C,h,h1C])
ubcC1 =  concatenate([u0C,u,u1C])
GbcC1 =  concatenate([G0C,G,G1C])

hbcC1_c = copyarraytoC(hbcC1)
ubcC1_c = copyarraytoC(ubcC1)
GbcC1_c = copyarraytoC(GbcC1)

Eni = HankEnergyall(xbcC1_c,hbcC1_c,ubcC1_c,bbcC1_c,g,n + 2*niBC,niBC,dx)
Pni = uhall(xbcC1_c,hbcC1_c,ubcC1_c,n + 2*niBC,niBC,dx)
Mni = hall(xbcC1_c,hbcC1_c,n + 2*niBC,niBC,dx)
Gni = Gall(xbcC1_c,GbcC1_c,n + 2*niBC,niBC,dx)

deallocPy(hbcC1_c)
deallocPy(ubcC1_c)
deallocPy(GbcC1_c)

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
   

ubc_c = mallocPy(unbc)
hbc_c = mallocPy(hnbc)
wbc_c = mallocPy(hnbc)
Gbc_c = mallocPy(hnbc)
bbc_c = mallocPy(bnbc)

tts.append(t)
nwg1s.append(h1)    
nwg2s.append(h1)  
nwg3s.append(h1) 
nwg4s.append(h1) 
nwg5s.append(h1) 

t = 0.0
#Just an FEM solve here
while t < endt: 
    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g);
    
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
    
    t = t + dt
    print(t)


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 

ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bbc_c)    

ubcC = copyarrayfromC(ubc_c,unbc)  
uC = ubcC[unBC:-unBC:2]
hbcC = copyarrayfromC(hbc_c,hnbc)
wbcC = copyarrayfromC(wbc_c,hnbc)
GbcC = copyarrayfromC(Gbc_c,hnbc)
bbcC = copyarrayfromC(bbc_c,bnbc)

hbcC1 =  concatenate([h0C,hC,h1C])
ubcC1 =  concatenate([u0C,uC,u1C])
GbcC1 =  concatenate([G0C,GC,G1C])

hbcC1_c = copyarraytoC(hbcC1)
ubcC1_c = copyarraytoC(ubcC1)
GbcC1_c = copyarraytoC(GbcC1)

En = HankEnergyall(xbcC1_c,hbcC1_c,ubcC1_c,bbcC1_c,g,n + 2*niBC,niBC,dx)
Pn = uhall(xbcC1_c,hbcC1_c,ubcC1_c,n + 2*niBC,niBC,dx)
Mn = hall(xbcC1_c,hbcC1_c,n + 2*niBC,niBC,dx)
Gn = Gall(xbcC1_c,GbcC1_c,n + 2*niBC,niBC,dx)


x0 = x[0] - 0.5*dx
x1 = -2*blength
x2 = 0
x3 = x[-1] + 0.5*dx
MA = MassA(g,h1,h0,x0,x1,x2,x3)
HA = HamilA(g,h1,h0,x0,x1,x2,x3)
"""
s = wdir +  "outListLast.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
    n = len(x)         
    for j in range(n):
        writefile2.writerow([str(x[j]), str(hC[j]) , str(GC[j]) , str(uC[j]),str(b[j]),str(hC[j] + b[j])])
        
s = wdir +  "outSingLast.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G" ,"Energyi", "Massi", "Momentumi", "Gi"  ])   
    writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni)  ]) 
    
  
nt = len(tts)

s = wdir +  "outWGList.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["t(s)","WG1s","WG2s" ,"WG3s" ,"WG4s" ,"WG5s"  ])        
    n = len(x)         
    for j in range(nt):
        writefile2.writerow([str(tts[j]),str(nwg1s[j]),str(nwg2s[j]),str(nwg3s[j]),str(nwg4s[j]),str(nwg5s[j])])
        
s = wdir + "wG1Num.dat"
with open(s,'w') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[0]/h1 ," ",1.5*((nwg1s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG2Num.dat"
with open(s,'w') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[1]/h1 ," ",1.5*((nwg2s[j] - h1)/h1))
        file1.write(s)
        
s = wdir + "wG3Num.dat"
with open(s,'w') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[2]/h1 ," ",1.5*((nwg3s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG4Num.dat"
with open(s,'w') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[3]/h1 ," ",1.5*((nwg4s[j] - h1)/h1))
        file1.write(s)

s = wdir + "wG5Num.dat"
with open(s,'w') as file1:
    for j in range(nt):
        s ="%3.8f%5s%1.15f\n" %(tts[j]*sqrt(g/h1) - WGloc[4]/h1 ," ",1.5*((nwg5s[j] - h1)/h1))
        file1.write(s)
"""



deallocPy(h_c)
deallocPy(G_c)
deallocPy(u_c)

deallocPy(ubc_c)
deallocPy(hbc_c)
deallocPy(wbc_c)
deallocPy(Gbc_c)
deallocPy(bbc_c)

deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
