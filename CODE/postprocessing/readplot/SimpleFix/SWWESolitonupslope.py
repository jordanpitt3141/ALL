import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os


wdir = "../../../../../data/raw/swdata/lsolslope2NTA/"
sdir = "../../../../../data/postprocessing/solslopeN/o2/100s/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

xbeg = 174
gap = 1
xend = 179
g = 9.81
         
#time = 0.0995978291745
filen = 20*500
#filen = 1

s = wdir + "saveoutputts" +str(filen)+ ".txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    bed = []
    u = []
    he = []
    ue = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            h.append(float(row[3]))
            u.append(float(row[5]))
            bed.append(float(row[6]))

            
            
                
        j = j + 1
x = arange(-200,200 + dx, dx)
xbegi = int((xbeg - x[0])/dx)    
xendi = int((xend - x[0])/dx)   
h = array(h[xbegi:xendi:gap])
b = array(bed[xbegi:xendi:gap])
x = array(x[xbegi:xendi:gap])
u = array(u[xbegi:xendi:gap])


n = len(x)
s = sdir + "SWstage.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i] + b[i])
        file1.write(s)
        
s = sdir + "SWh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "SWu.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)

s = sdir + "SWbed.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",b[i])
        file1.write(s)
  