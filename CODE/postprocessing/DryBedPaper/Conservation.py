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


wdir = "/home/jp/Documents/PhD/project/data/DryBedPaper/DambreakConvergence/HintonLumpuhAlpha0p1/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/DryBedPaper/HintonLumpuhAlpha0p1/C1/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

h1 = 1
g = 9.81
    
s = wdir + "Conservation.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    dxs = []
    Eerr = []
    Merr = []
    Gerr = []
    Perr = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            dxs.append(float(row[0]))
            Eerr.append(float(row[2]))
            Merr.append(float(row[3]))
            Perr.append(float(row[4]))
            Gerr.append(float(row[5]))
    
        j = j + 1

    


n = len(dxs)

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