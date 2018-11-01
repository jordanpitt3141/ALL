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

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(a0,a1,g,x,t0,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    w = zeros(n)
    b = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        w[i] = h[i]
        u[i] =  c* (1 - a0 / h[i])
        G[i] = 2.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**4*h[i] + h[i]*u[i] - 4.0/3*a0*a1**2*c*k**2*sech(k*(x[i] - c*t0))**4*tanh(k*(x[i] - c*t0))**2 - 4.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**2*h[i]*tanh(k*(x[i] - c*t0))**2
        
    return h,u,G,w,b

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    


#Soliton  
wdir = "/home/jp/Documents/PhD/project/master/FigureData/Presentation/CTAC/Soliton/Numerical/" 
if not os.path.exists(wdir):
    os.makedirs(wdir)


wdirD = wdir + "Data/"
if not os.path.exists(wdirD):
        os.makedirs(wdirD)
        
wdirT = wdir + "Texs/"
if not os.path.exists(wdirT):
        os.makedirs(wdirT)

ki= 10


a0 = 1.0
a1 = 0.7
g = 9.81
k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
c = sqrt(g*(a0 + a1))


startx = -50
endx = 250

startt = 0.0
endt = 50.1


dx = 100.0/ 2**ki
Cr = 0.5
l = 1.0 / (sqrt(g*(a0 + a1)))
dt = Cr*l*dx

t = startt
theta = 1.2


x = arange(startx,endx +0.1*dx, dx)
xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

xbMend = array([x[-1] + 0.5*dx, x[-1] + 5*dx/6.0, x[-1] + 7*dx/6.0, x[-1] + 1.5*dx])
xbMbeg = array([x[0] - 1.5*dx, x[0] - 7*dx/6.0,x[0] - 5*dx/6.0 , x[0] -0.5*dx])


h,u,G,w,b = solitoninit(a0,a1,g,x,startt,dx)


hMbeg,uMbeg,GMbeg,wMbeg,bta = solitoninit(a0,a1,g,xhuMbeg,startt,dx)
hMend ,uMend ,GMend,wMend,bta = solitoninit(a0,a1,g,xhuMend,startt,dx)

hta,uta,Gta,wta,bMbeg = solitoninit(a0,a1,g,xbMbeg,startt,dx)
hta,uta,Gta,wta,bMend = solitoninit(a0,a1,g,xbMend,startt,dx)

n = len(x)  

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

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


t = 0.0
#Just an FEM solve here
j = 0
while t < endt: 
    evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t,0,0,0,0,0,0,0,0,0,0)
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    if j%10 == 0:
        

        timesn = "{:1.1f}".format(t)

        fileDn = wdirD + "ht="+ timesn + ".dat"
        with open(fileDn,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.20f\n" %(x[i]," ",hC[i])
                file1.write(s) 
                
        fileTn = wdirT + "h" + str(j/10) +  ".tex"
        #Make A test file
        with open(fileTn,'w') as file1:
        
            sheader= " \\documentclass[]{standalone} \n \\usepackage{pgfplots} \n \\usepgfplotslibrary{fillbetween} \n \\usepackage{tikz} \n" \
            + " \\usepackage{amsmath} \n \\usepackage{pgfplots}  \n \\usepackage{sansmath} \n \\sansmath \n \\usetikzlibrary{calc} " \
            + "\\pgfplotsset{compat = newest, every axis plot post/.style={line join=round}, label style={font=\\Huge},every tick label/.append style={font=\\Huge} }\n"
            file1.write(sheader) 
            
            sdocs  = "\\begin{document}  \n	\\begin{tikzpicture}\n"   
            file1.write(sdocs) 
            
            saxisoptions = "\\fontfamily{cmss} \n	\\begin{axis}[ \n width=40cm, \n height = 20cm, \n	xtick={-50,0,50,100,150,200,250}, \n	ytick = {0,0.5,1,1.5,2}, \n" \
        	+ "xmin=-50, \n xmax=250, \n ymin =-0.1, \n ymax = 2.0,\n  xlabel=$x(m)$, \n	ylabel=$z(m)$ ]\n"
            file1.write(saxisoptions) 
            
            snode = "\\node[label={\Huge$t="+timesn+"s$}] at (axis cs:230,1.8) {}; \n"
            file1.write(snode) 
            
            splot = "	\\addplot [name path=s,blue,no markers,ultra thick] table {../Data/ht="+timesn+".dat }; \n \\draw[brown!80!black,ultra thick,name path=b] (-50, 0  ) -- (250,0); \n " \
        	+ "\\path[name path=a] (axis cs:-50,-0.1) -- (axis cs:250,-0.1);\n"
            file1.write(splot) 
            
            sfillplots = "\\addplot [thick,color=brown!60!black,fill=brown!60!black, fill opacity=0.3] fill between[of=b and a]; \n " \
        	+ "\\addplot [thick,color=blue,	fill=blue,fill opacity=0.3] fill between[of=s and b]; \n"
         
            file1.write(sfillplots) 
         
            send = "\\end{axis} \n \\end{tikzpicture} \n \\end{document}"
            file1.write(send) 
    
    t = t + dt
    j = j + 1
    print(t)



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


