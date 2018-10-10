import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os



alpha = "1"

wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/bigsmooth/o2/11/" +alpha+"/"

sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/DamBreak/" +alpha+"/"

if not os.path.exists(sdir):
        os.makedirs(sdir)

s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    h = []
    G = []
    u = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            G.append(float(row[5]))
            u.append(float(row[6]))
            diff = float(row[8])
    
        j = j + 1

    n = len(x)        
    x = array(x)
    h = array(h)
    G = array(G)
    u = array(u)

gap = 8
sxl = 350
exl = 650
sxi = int(sxl / dx   )
exi = int(exl / dx)

xn = x[sxi:exi:gap]
hn = h[sxi:exi:gap]

nn = len(xn)

s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(nn):
        s ="%3.8f%5s%1.20f\n" %(xn[i]," ",hn[i])
        file1.write(s) 

   