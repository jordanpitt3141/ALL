import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

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
             
    return h,u,G,b,h


#wdir = "/home/jp/Documents/PhD/project/data/DryBedPaper/DambreakConvergence/Smoothalpha0p1/15/"
wdir = "/home/jp/Documents/PhD/project/data/DryBedPaper/DambreakConvergence/NewSmallTime/Normal/16/"

h1 = 1
g = 9.81
    
s = wdir + "outList10.0s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    x = []
    h = []
    u = []
    G = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            x.append(float(row[0]))
            h.append(float(row[1]))
            G.append(float(row[2]))
            u.append(float(row[3]))
    
        j = j + 1

n = len(x)

s = wdir + "outSing10.0s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    row1 = next(readfile)
    row1 = next(readfile)
    
    t = float(row1[2])

hSWWE,uSWWE,GSWWE,bSWWE,wSWWE = DrybedSWWANA(h1,x,t,g)


"""
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",Merr[i])
        file1.write(s) 

s = sdir + "H.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",Eerr[i])
        file1.write(s) 

s = sdir + "uh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",Perr[i])
        file1.write(s) 

s = sdir + "G.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",Gerr[i])
        file1.write(s) 
"""