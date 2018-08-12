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

def cot(x):
    return 1.0/ tan(x)
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton(x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
   
def MollifyFunc(C,x,e):
    ie = 1.0 / e
    if(abs(ie*x) <1):
        return ie*C*exp(1.0/(abs(ie*x)**2 - 1))
    else:
        return 0   
   
def ForcedbedM(x,t,beta,a0,a1,x0):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    
    MollPy = zeros(n)
    
    D = 1.0/0.444994
    eta = 3*dx
    rublen = 600
    
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        
        MollPy[i] = MollifyFunc(D,x[i] - x[n/2],eta)
        
        if (x[i] <= cot(beta)):
            b[i] = - x[i]*tan(beta)
            
        else:
            b[i] = -1
        
        if b[i] >= 0:
            h[i] = 0
            u[i] = 0
            w[i] = b[i]
        else:
            w[i] = soliton(x[i] - x0,t,g,a0,a1) - a0
            h[i] = w[i] - b[i]
            u[i]= -c* (1 - a0/ (w[i] + a0))

    
       
    b0 =- (x[0] - dx)*tan(beta) 
    
    G = getGfromupy(h,u,b,0,0,0,a0,b0,-1,dx)     

    return h,u,G,b,w
  
def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var

#Forcing Problem    
wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Synolakis/H0p0185/FDVMNoMull/"
#wdir = "/home/jp/Documents/PhD/project/data/Tests/Syn/FDVM/"
if not os.path.exists(wdir):
    os.makedirs(wdir)



g = 1

H = 0.0185
d = 1
x1 = 38.5
x0 = 19.85
beta = arctan(1.0/x0)

startx = -30
sx = startx
endx = 150
ex = endx

startt = 0.0
st = startt
endt = 71
et = endt
dx = 0.05
l =  0.1
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
  

theta = 1
 

h,u,G,b,w = ForcedbedM(x,t,beta,d,H,x1)

hMbeg,uMbeg,GMbeg,bMbeg,wMbeg = ForcedbedM(xhMbeg,t,beta,d,H,x1)
hMend,uMend,GMend,bMend,wMend = ForcedbedM(xhMend,t,beta,d,H,x1)

dbMbeg = (1.0/x0)*ones(hnBC)
dbMend = ones(hnBC)
 
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
b0C = - tan(beta)*xbegC
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




ts = [30,40,50,60,65,70]
#Just an FEM solve here
while t < endt:  

    if close(t,ts,dt):
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
        
        uh = h *u
        uhiC = array(uiC)*array(hiC)
        
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
        
        s = wdir +  "outList" + str(t)+"s.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
                       
            for j in range(n):
                writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
                
        s = wdir +  "outSing" + str(t)+"s.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G"  ])   
            writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn) ]) 
            
            
        s = wdir + "w" +str(t)+ ".dat"
        with open(s,'a') as file1:
            for j in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[j]," ",wiC[j])
                file1.write(s)

        s = wdir + "b" +str(t)+ ".dat"
        with open(s,'a') as file1:
            for j in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[j]," ",b[j])
                file1.write(s)

        
    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g);
    t = t + dt
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

uh = h *u
uhiC = array(uiC)*array(hiC)


hnorm = norm(hiC -h, ord=1)/ norm(h, ord=1)
Gnorm = norm(GiC -G, ord=1)
unorm = norm(uiC -u, ord=1)
uhnorm = norm(uhiC -uh, ord=1)



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


s = wdir +  "outList" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
               
    for j in range(n):
        writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
        
s = wdir +  "outSing" + str(t)+"s.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G"  ])   
    writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn) ]) 
    
    
s = wdir + "w" +str(t)+ ".dat"
with open(s,'a') as file1:
    for j in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[j]," ",wiC[j])
        file1.write(s)

s = wdir + "b" +str(t)+ ".dat"
with open(s,'a') as file1:
    for j in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[j]," ",b[j])
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