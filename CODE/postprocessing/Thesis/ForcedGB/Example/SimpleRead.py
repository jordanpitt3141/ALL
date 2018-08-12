import csv
from numpy.linalg import norm
from scipy import *
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones

def IterateIn(list1,word):
    for l in list1:
        if word in l:
            break
    return l

meth1 = "FEVM2"
num = "10"

wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/Forced/Wet/" +meth1+ "/" +num+ "/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/Forced/Wet/" +meth1+"/Ex"+num+"/"


if not os.path.exists(sdir):
        os.makedirs(sdir)

ts = "10.0"
outfiles=os.listdir(wdir)
sname = IterateIn(outfiles,"Sing"+ts)
lname = IterateIn(outfiles,"List"+ts)  

s = wdir + lname
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    b = []
    w = []
    h = []
    u = []
    G = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            x.append(float(row[0]))
            h.append(float(row[1]))
            G.append(float(row[2]))
            u.append(float(row[3]))
            b.append(float(row[4]))
            w.append(float(row[5]))
                  
        j = j + 1

    n = len(x)        
    x = array(x)
    b = array(b)
    w = array(w)
    h = array(h)
    u = array(u)    
    G = array(G)  

dx = x[1] - x[0]
incr = 1
x0 = -50
x1 = 25

x0i = int((x0 - x[0]) / dx)
x1i = int((x1 - x[0]) / dx)

rx = x[x0i:x1i:incr]
rh = h[x0i:x1i:incr]
rG = G[x0i:x1i:incr]
ru = u[x0i:x1i:incr]
rb = b[x0i:x1i:incr]
rw = w[x0i:x1i:incr]

n= len(rx)

s = sdir + "Stage"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(rx[i]," ",rw[i])
        file1.write(s)
    
s = sdir + "Bed"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(rx[i]," ",rb[i])
        file1.write(s)

s = sdir + "h"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(rx[i]," ",rh[i])
        file1.write(s)

s = sdir + "u"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(rx[i]," ",ru[i])
        file1.write(s)

s = sdir + "G"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(rx[i]," ",rG[i])
        file1.write(s)
