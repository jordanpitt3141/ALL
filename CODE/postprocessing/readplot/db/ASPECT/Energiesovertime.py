import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones



wdir = "../../../../../../data/raw/P1P2P3/Dambreak10to1/"

g = 9.81
h1 = 10
h0 = 1
#from 0 to 168
Masss = [3200.0]
Momes = [0]
Momecorr = [0]
Energs = [148131.0]



ts = [0]
for i in range(1):
    s = wdir + "out"+str(i)+".txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        h = []
        u = []
        x = []
        j = -1
        for row in readfile:       
            if (j >= 0):
                x.append(float(row[0]))
                t =float(row[1])
                theta = float(row[2])
                h.append(float(row[3]))
                u.append(float(row[4]))
                Massi = float(row[6])
                Momei = float(row[7])
                Energi = float(row[8])
                Massc = float(row[9])
                Momec = float(row[10])
                Energc = float(row[11])
                        
            j = j + 1
    
        n = len(x)        
        x = array(x)
        u = array(u)
        h = array(h)
        
    MomentumCorrection = g*t*0.5*(h1**2 -h0**2)
    
    Momecorr.append(MomentumCorrection)
    Masss.append(Massc)
    Momes.append(Momec)
    Energs.append(Energc)
    ts.append(t)
   