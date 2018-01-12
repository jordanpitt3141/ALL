import csv
import os
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


wdir = "../../../../../data/raw/P1P2P3/Roeberf/"
sdir = "../../../../../data/postprocessing/P1P2P3/Roeberf/"
 
if not os.path.exists(sdir):
    os.makedirs(sdir)


s = wdir +"eWGs.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    et = []
    eWG1 = []
    eWG2 = []
    eWG3 = []
    eWG4 = []
    eWG5 = []
    eWG6 = []
    eWG7 = []
    eWG8 = []
    eWG9 = []
    eWG10 = []
    eWG11 = []
    eWG12 = []
    eWG13 = []
    eWG14 = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            et.append(float(row[0]))
            eWG1.append(float(row[1]))
            eWG2.append(float(row[2]))
            eWG3.append(float(row[3]))
            eWG4.append(float(row[4]))
            eWG5.append(float(row[5]))
            eWG6.append(float(row[6]))            
            eWG7.append(float(row[7]))
            eWG8.append(float(row[8]))
            eWG9.append(float(row[9]))
            eWG10.append(float(row[10]))            
            eWG11.append(float(row[11]))
            eWG12.append(float(row[12]))            
            eWG13.append(float(row[13]))
            eWG14.append(float(row[14]))                
        j = j + 1


s = wdir +"nWGs.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nt = []
    nWG1 = []
    nWG2 = []
    nWG3 = []
    nWG4 = []
    nWG5 = []
    nWG6 = []
    nWG7 = []
    nWG8 = []
    nWG9 = []
    nWG10 = []
    nWG11 = []
    nWG12 = []
    nWG13 = []
    nWG14 = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            nt.append(float(row[0]))
            nWG1.append(float(row[1]))
            nWG2.append(float(row[2]))
            nWG3.append(float(row[3]))
            nWG4.append(float(row[4]))
            nWG5.append(float(row[5]))
            nWG6.append(float(row[6]))            
            nWG7.append(float(row[7]))
            nWG8.append(float(row[8]))
            nWG9.append(float(row[9]))
            nWG10.append(float(row[10]))            
            nWG11.append(float(row[11]))
            nWG12.append(float(row[12]))            
            nWG13.append(float(row[13]))
            nWG14.append(float(row[14]))                
        j = j + 1