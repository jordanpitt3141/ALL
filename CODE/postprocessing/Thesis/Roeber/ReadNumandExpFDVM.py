import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

meth = 'FDVM'

expdir = "/home/jp/Documents/PhD/project/data/Experimental/Roeber/Out/Trial8/"  

#wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Roeber/Trial8/"+meth+"/r1/" 
wdir = "/home/jp/Documents/PhD/project/data/ThesisRaws/Experiment/Roeber/Trial8/o2bedBC/r5/" 
sdirb = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/Roeber/"+meth+"/"  
sdirexp =sdirb + "Exp/"
sdirnum =sdirb + "Num/"

if not os.path.exists(sdirexp):
   os.makedirs(sdirexp)
   
if not os.path.exists(sdirnum):
   os.makedirs(sdirnum)

nWG = 14  
texp = []
WGsall = []
for i in range(nWG):
    WGsall.append([]) 

s = expdir + "WGs.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for row in readfile:       
            texp.append(float(row[0]))
            for i in range(nWG):
                WGsall[i].append(float(row[i + 1]))    

tn = []
WGsalln = []
for i in range(nWG):
    WGsalln.append([])    
    
s = wdir + "nWGs.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    j = -1
    for row in readfile:    
        if(j >= 0):
            tn.append(float(row[0]))
            for i in range(nWG):
                WGsalln[i].append(float(row[i + 1]))  
        j =j + 1

    
dt = tn[1] -tn[0]
incr = 10
st = 20
et = 80

sti =int(st/dt) 
eti =int(et/dt) 

rtn = tn[sti : eti :incr ]
rWGsalln = []
for i in range(nWG):
    rWGsalln.append(WGsalln[i][sti : eti :incr ]) 
    




for i in range(nWG):
    
    n = len(texp)
    s = sdirexp + "WG" +str(i +1)+".dat"
    with open(s,'w') as file1:
        for j in range(n):
            ss ="%3.8f%5s%1.15f\n" %(texp[j]," ",WGsall[i][j])
            file1.write(ss)
     
    n = len(rtn)       
    s = sdirnum + "WG" +str(i +1)+".dat"
    with open(s,'w') as file1:
        for j in range(n):
            ss ="%3.8f%5s%1.15f\n" %(rtn[j]," ",rWGsalln[i][j])
            file1.write(ss)