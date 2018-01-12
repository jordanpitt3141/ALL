import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


wdir = "../../../../../data/raw/solslopelarger10p1R/o2/"

timeinsecs = 30

gap = 1
g = 9.81

mult = 175
         
#time = 0.0995978291745
filen = 400*mult

#s = wdir + "saveoutputts" + str(int(filen)) + ".txt"
 
#s = wdir + "saveoutputtslast.txt"
s = wdir + "out" +str(filen)+ ".txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    bed = []
    u = []
    he = []
    ue = []
    x = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            #he.append(float(row[8]))
            #ue.append(float(row[9]))
            bed.append(float(row[7]))
            #diffuse = float(row[8])

            
                
        j = j + 1

    n = len(x)        
    x = array(x)
    u = array(u)
    h = array(h)
    bed = array(bed)
    
wdir = "../../../../../data/raw/swdata/lsolslope2NTAsteep/"

timeinsecs = 30

gap = 1
g = 9.81
         
#time = 0.0995978291745
filen = 20*mult

#s = wdir + "saveoutputts" + str(int(filen)) + ".txt"
 
#s = wdir + "saveoutputtslast.txt"
s = wdir + "saveoutputts" +str(filen)+ ".txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    hSW = []
    bedSW = []
    uSW = []
    xSW = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            hSW.append(float(row[3]))
            uSW.append(float(row[5]))
            bedSW.append(float(row[6]))

            
                
        j = j + 1
    
    xSW = arange(-200,200+dx,dx)
    uSW = array(uSW)
    hSW = array(hSW)
    bedSW = array(bedSW)
    
   