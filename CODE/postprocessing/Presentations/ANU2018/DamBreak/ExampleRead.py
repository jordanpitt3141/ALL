import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os





wdirbase = "/home/jp/Documents/PhD/project/data/ThesisRaw/bigsmooth/"

sdir = "/home/jp/Documents/PhD/project/master/presentations/ANU2018/Diagrams/SteepGradients/ExpPlots/Data/"

#dxis = [9]
dxis = [5,7,9,11]

for dxi in dxis:
    
    method = "o3"
    dxstr = str(dxi)
    diffstr = "20"
    
    wdir = wdirbase + method + "/"  + dxstr + "/" +diffstr+ "/"
    
    if not os.path.exists(sdir):
            os.makedirs(sdir)
    
    s = wdir + "outlast.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        G = []    
        h = []
        u = []
        b = []
        Gt = []    
        ht = []
        ut = []
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
                b.append(float(row[7]))
             
                
                
    
                
                
                    
            j = j + 1

    print(len(x))
    startx = 300
    endx = 700
    
    plotdx = 0.1
    
    
    sxi = int(startx/dx)
    exi = int(endx/dx)
    zdx = max(int(plotdx/dx),1)
    
    xn = x[sxi:exi:zdx]
    hn = h[sxi:exi:zdx]


    n = len(xn)
    s = sdir +dxstr + "h.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(xn[i]," ",hn[i])
            file1.write(s) 


