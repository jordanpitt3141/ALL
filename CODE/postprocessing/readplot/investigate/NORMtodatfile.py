import csv
import os
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


wdir = "../../../../../data/raw/SolitondxfixTIME/"
sdir = "../../../../../data/postprocessing/SolitondxFEM/"

timeinsecs = 30

gap = 8
g = 9.81
         
#filen = 200*2560
#s = wdir + "saveoutputts" + str(int(filen)) + ".txt"
 
#s = wdir + "saveoutputtslast.txt"
 
if not os.path.exists(sdir):
    os.makedirs(sdir)


s = wdir +"savenorms.csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    dx = []
    L1h = []
    L1u = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx.append(float(row[0]))
            L1h.append(float(row[1]))
            L1u.append(float(row[2]))

            

            
            
                
        j = j + 1

n = len(dx)        
   
s = sdir + "L1h.dat"

with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dx[i]," ",L1h[i])
        file1.write(s)
        
s = sdir + "L1u.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dx[i]," ",L1u[i])
        file1.write(s)