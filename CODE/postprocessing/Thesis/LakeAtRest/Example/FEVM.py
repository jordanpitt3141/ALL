import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

def IterateIn(list1,word):
    for l in list1:
        if word in l:
            break
    return l

meth = "FDVM2WB"
numexp = "10"

wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/LakeAtRest/Dry/" +meth+"/" +numexp+"/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/LakeAtRest/" +meth+"/Ex" +numexp+"/"

if not os.path.exists(sdir):
    os.makedirs(sdir)
    

outfiles=os.listdir(wdir)
sname = IterateIn(outfiles,"Sing")
lname = IterateIn(outfiles,"List")

s = wdir + lname 
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    w = []
    b = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            x.append(float(row[0]))
            b.append(float(row[4]))
            w.append(float(row[5]))

    
        j = j + 1


        
x = array(x)
w = array(w)
b = array(b)
   

n = len(x)
s = sdir + "w.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",w[i])
        file1.write(s) 
        
s = sdir + "b.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
        file1.write(s) 
