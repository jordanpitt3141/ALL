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
    
def ForcedbedM(x,t,beta,a0,a1,x0):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        
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
wdir = "/home/jp/Documents/PhD/project/master/FigureData/Presentation/CTAC/Synolakis/NumericalLonger2/" 
if not os.path.exists(wdir):
    os.makedirs(wdir)


wdirD = wdir + "Data/"
if not os.path.exists(wdirD):
        os.makedirs(wdirD)
        
wdirT = wdir + "Texs/"
if not os.path.exists(wdirT):
        os.makedirs(wdirT)

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
endt = 120.0
et = endt
dx = 0.05
l =  0.1
dt = l*dx

t = startt


x = arange(startx,endx +0.1*dx, dx)

xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

xbMbeg = [x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx]
xbMend = [x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx]
 

theta = 1.2

h,u,G,b,w = ForcedbedM(x,t,beta,d,H,x1)

hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,beta,d,H,x1)
hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend,t,beta,d,H,x1)

hta,uta,Gta,bMbeg,wta = ForcedbedM(xhuMbeg,t,beta,d,H,x1)
hta,uta,Gta,bMend,wta = ForcedbedM(xhuMend,t,beta,d,H,x1)

b = b + 1
w = w + 1
bMbeg = bMbeg + 1
bMend = bMend + 1
wMbeg = wMbeg + 1
wMend = wMend + 1

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


b0C = - tan(beta)*xbegC
b1C = b[-1]*ones(niBC)


xbcC =  concatenate([xbegC,x,xendC])
bbcC =  concatenate([b0C,b,b1C])
xbcC_c = copyarraytoC(xbcC)
bbcC_c = copyarraytoC(bbcC)

u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)

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

t = 0.0
j = 0
#Just an FEM solve here

fileDn = wdirD + "b.dat"
with open(fileDn,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
        file1.write(s) 
        
while t < endt: 
    
    if j%20 == 0:
        
        timesn = "{:3.1f}".format(t)
        
        hiC = copyarrayfromC(h_c,n)
        wiC = hiC + b
                

        fileDn = wdirD + "ht="+ timesn + ".dat"
        with open(fileDn,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.20f\n" %(x[i]," ",wiC[i])
                file1.write(s) 
                
        fileTn = wdirT + "h" + str(j/20) +  ".tex"
        #Make A test file
        with open(fileTn,'w') as file1:
        
            sheader= " \\documentclass[]{standalone} \n \\usepackage{pgfplots} \n \\usepgfplotslibrary{fillbetween} \n \\usepackage{tikz} \n" \
            + " \\usepackage{amsmath} \n \\usepackage{pgfplots}  \n \\usepackage{sansmath} \n \\sansmath \n \\usetikzlibrary{calc} " \
            + "\\pgfplotsset{compat = newest, every axis plot post/.style={line join=round}, label style={font=\\Huge},every tick label/.append style={font=\\Huge} }\n"
            file1.write(sheader) 
            
            sdocs  = "\\begin{document}  \n	\\begin{tikzpicture}\n"   
            file1.write(sdocs) 
            
            saxisoptions = "\\fontfamily{cmss} \n	\\begin{axis}[ \n width=40cm, \n height = 20cm, \n	xtick={-10,0,10,20,30,40,50,60}, \n	ytick = {0.9,0.95,1,1.05,1.1}, \n" \
        	+ "xmin=-10, \n xmax=60, \n ymin =0.9, \n ymax = 1.1,\n  xlabel=$x'$, \n	ylabel=$z'$ ]\n"
            file1.write(saxisoptions) 
            
            snode = "\\node[label={\Huge$t'="+timesn+"$}] at (axis cs:55,1.07) {}; \n"
            file1.write(snode) 
            
            splot = "	\\addplot [name path=s,blue,no markers,ultra thick] table {../Data/ht="+timesn+".dat }; \n \\addplot [name path=b,brown!60!black,no markers,ultra thick] table {../Data/b.dat };  " \
        	+ "\\path[name path=a] (axis cs:-30,-0.1) -- (axis cs:150,-0.1);\n"
            file1.write(splot) 
            
            sfillplots = "\\addplot [thick,color=brown!60!black,fill=brown!60!black, fill opacity=0.3] fill between[of=b and a]; \n " \
        	+ "\\addplot [thick,color=blue,	fill=blue,fill opacity=0.3] fill between[of=s and b]; \n"
         
            file1.write(sfillplots) 
         
            send = "\\end{axis} \n \\end{tikzpicture} \n \\end{document}"
            file1.write(send) 

    evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g);
    t = t + dt
    j = j+1
    print(t)



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
