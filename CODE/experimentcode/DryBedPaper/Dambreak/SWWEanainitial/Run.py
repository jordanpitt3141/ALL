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

  
def Dambreak(h0,h1,x0,x):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    
    for i in range(n):
        
        if (x[i] < x0):
            h[i] = h1
        else:
            h[i] = h0
    
    return h,u,G,b,h
    
def DrybedSWWANA(h1,x,t,g):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    b = zeros(n)
    
    
    for i in range(n):
         if(x[i] >= -t*sqrt(g*h1) and x[i] <= 2*t*sqrt(g*h1) ):
             u[i] = 2.0 / 3.0 *(sqrt(g*h1) + x[i] / t)
             h[i] = 4.0 / (9.0 * g) *(sqrt(g*h1) - 0.5*x[i] / t)**2
             ux = 2.0 / 3.0 *(1.0 / t)
             uxx = 0
             hx = 2.0 / (9.0 * g * t*t) *(x[i] - 2*t*sqrt(g*h1))
             G[i] = u[i]*h[i] - h[i]*h[i]*hx*ux
         elif(x[i] < -t*sqrt(g*h1)):
             h[i] = h1
             
    return h,u, G,b,h

#Forcing Problem    
wdir = "/home/jp/Documents/PhD/project/data/DryBedPaper/Dambreak/SWWEanateq1/50s/dx0p01/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)


g = 9.81
h1 = 1.0
h0 = 0.0
x0 = 0


startx = -400
sx = startx
endx = 400
ex = endx

startt = 1.0
st = startt
endt = 50.0
et = endt
dx = 0.01
l =  0.01
dt = l*dx

t = startt


x = arange(startx,endx +0.1*dx, dx)

xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

xbMbeg = [x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx]
xbMend = [x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx]
 

theta = 1.2

h,u,G,b,w = DrybedSWWANA(h1,x,startt,g)

hMbeg,uMbeg,GMbeg,bta,wMbeg = DrybedSWWANA(h1,xhuMbeg,startt,g)
hMend ,uMend ,GMend ,bta,wMend = DrybedSWWANA(h1,xhuMend,startt,g)

hta,uta,Gta,bMbeg,wta = DrybedSWWANA(h1,xbMbeg,startt,g)
hta,uta,Gta,bMend,wta = DrybedSWWANA(h1,xbMend,startt,g)



n = len(x)
hnBC = 3
hnbc = 3*n + 2*hnBC
bnMBC = 7
bnBC = 4
bnbc = 3*n + 1 + 2*(bnBC -1)
unBC = 3
unbc = 2*n + 1 + 2*(unBC -1)


niBC = 4
xbegC = arange(sx - niBC*dx,sx,dx)
xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 


b0C = b[0]*ones(niBC)
b1C = b[-1]*ones(niBC)
u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)

xbcC =  concatenate([xbegC,x,xendC])
bbcC =  concatenate([b0C,b,b1C])
hbcC =  concatenate([h0C,h,h1C])
ubcC =  concatenate([u0C,u,u1C])
GbcC =  concatenate([G0C,G,G1C])

xbcC_c = copyarraytoC(xbcC)
bbcC_c = copyarraytoC(bbcC)
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

s = wdir +  "outList" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
               
    for j in range(n):
        writefile2.writerow([str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(b[j]),str(w[j])])
        
s = wdir +  "outSing" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum"])   
    writefile2.writerow([str(dx),str(dt),str(t),str(Eni),str(Mni),str(Pni),str(Gni) ]) 

t = startt
#Just an FEM solve here
while t < endt: 
    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g);
    t = t + dt
    print(t)


hSWWE,uSWWE,GSWWE,bSWWE,wSWWE = DrybedSWWANA(h1,x,t,g)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 

ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bbc_c)    

ubcC = copyarrayfromC(ubc_c,unbc)  
uC = ubcC[unBC:-unBC:2]
hbcC = copyarrayfromC(hbc_c,hnbc)
wbcC = copyarrayfromC(wbc_c,hnbc)
GbcC = copyarrayfromC(Gbc_c,hnbc)
bbcC = copyarrayfromC(bbc_c,bnbc)


u0Cn = uC[0]*ones(niBC)
u1Cn = uC[-1]*ones(niBC)   
h0Cn = hC[0]*ones(niBC)
h1Cn = hC[-1]*ones(niBC)
G0Cn = GC[0]*ones(niBC)
G1Cn = GC[-1]*ones(niBC)

hbcC =  concatenate([h0Cn,hC,h1Cn])
ubcC =  concatenate([u0Cn,uC,u1Cn])
GbcC =  concatenate([G0Cn,GC,G1Cn])

hbcC_c = copyarraytoC(hbcC)
ubcC_c = copyarraytoC(ubcC)
GbcC_c = copyarraytoC(GbcC)

En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)

Eerr = abs(En- Eni)/ abs(Eni)
Perr = abs(Pn- Pni)
Gerr = abs(Gn- Gni)
Merr = abs(Mn- Mni)/ abs(Mni)



s = wdir +  "outList" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
               
    for j in range(n):
        writefile2.writerow([str(x[j]), str(hC[j]) , str(GC[j]) , str(uC[j]),str(b[j]),str(wC[j])])
        
s = wdir +  "outSing" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G","Energy0", "Mass0", "Momentum0", "G0"   ])   
    writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni) ]) 



deallocPy(hbcC_c)
deallocPy(ubcC_c)
deallocPy(GbcC_c)






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

