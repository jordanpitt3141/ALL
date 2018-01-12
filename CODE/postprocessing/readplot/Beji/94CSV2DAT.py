import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdatadir = "../../../../../data/raw/P1P2P3/Beji/"
sdatadir = "../../../../../data/postprocessing/Beji94FEM/o2/"
exp = "sl"
wdir = wdatadir + exp+ "/"
sexpdir = sdatadir + exp + "/"

endt = 60
begt = 40  
nts = []
nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []

s = wdir + "NumWaveGauge.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            nts.append(float(row[0]))
            nwg1s.append(float(row[1]))
            nwg2s.append(float(row[2]))
            nwg3s.append(float(row[3]))
            nwg4s.append(float(row[4]))
            nwg5s.append(float(row[5]))
            nwg6s.append(float(row[6]))
            nwg7s.append(float(row[7]))
            
            
         j = j + 1
         

ets = []
ewg1s = []
ewg2s = []
ewg3s = []
ewg4s = []
ewg5s = []
ewg6s = []
ewg7s = []        
s = wdir + "WaveGauge.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            ets.append(float(row[0]))
            ewg1s.append(float(row[1]))
            ewg2s.append(float(row[2]))
            ewg3s.append(float(row[3]))
            ewg4s.append(float(row[4]))
            ewg5s.append(float(row[5]))
            ewg6s.append(float(row[6]))
            ewg7s.append(float(row[7]))
            
            
         j = j + 1

ets0 = ets[1]

endtei = int(endt/ets[1]) + 1 
begtei = int(begt/ets[1]) - 1
ets = array(ets[begtei:endtei])
ewg1s = array(ewg1s[begtei:endtei])*100
ewg2s = array(ewg2s[begtei:endtei])*100
ewg3s = array(ewg3s[begtei:endtei])*100
ewg4s = array(ewg4s[begtei:endtei])*100
ewg5s = array(ewg5s[begtei:endtei])*100
ewg6s = array(ewg6s[begtei:endtei])*100
ewg7s = array(ewg7s[begtei:endtei])*100

        
ExpCom = []
ExpCom.append(ets)
ExpCom.append(ewg1s)
ExpCom.append(ewg2s)
ExpCom.append(ewg3s)
ExpCom.append(ewg4s)
ExpCom.append(ewg5s)
ExpCom.append(ewg6s)
ExpCom.append(ewg7s)


mult = int(ets0/nts[1])
endtni = int(endt/nts[1]) + 1
begtni = int(begt/nts[1]) - 1 

nts = array(nts[begtni:endtni:mult])
nwg1s = (array(nwg1s[begtni:endtni:mult])-0.4)*100
nwg2s = array(nwg2s[begtni:endtni:mult])*100
nwg3s = array(nwg3s[begtni:endtni:mult])*100
nwg4s = array(nwg4s[begtni:endtni:mult])*100
nwg5s = array(nwg5s[begtni:endtni:mult])*100
nwg6s = array(nwg6s[begtni:endtni:mult])*100
nwg7s = array(nwg7s[begtni:endtni:mult])*100
 
 
NumCom = []
NumCom.append(nts)
NumCom.append(nwg1s)
NumCom.append(nwg2s)
NumCom.append(nwg3s)
NumCom.append(nwg4s)
NumCom.append(nwg5s)
NumCom.append(nwg6s)
NumCom.append(nwg7s)


nc = len(NumCom)

for j in range(1,nc):
    sdir = sexpdir +"WaveGauge" + str(j) + "/"
    if not os.path.exists(sdir):
        os.makedirs(sdir)
    nn = len(nts)    
    s = sdir + "Numerical.dat"
    with open(s,'w') as file1:
        for i in range(nn):
            s ="%3.8f%5s%1.15f\n" %(NumCom[0][i]," ",NumCom[j][i])
            file1.write(s)
    ne = len(ets)        
    s = sdir + "Experimental.dat"
    with open(s,'w') as file1:
        for i in range(ne):
            s ="%3.8f%5s%1.15f\n" %(ExpCom[0][i]," ",ExpCom[j][i])
            file1.write(s)
 