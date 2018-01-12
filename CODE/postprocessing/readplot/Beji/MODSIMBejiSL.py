import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdatadir = "../../../../../data/raw/Beji94/o2/"
Stevedatadir = "../../../../../data/raw/swdata/BCTestSL/"
sdatadir = "../../../../../data/postprocessing/Beji94ChrisMODSIMSL/o2/"
exp = "sl"
wdir = wdatadir + exp+ "/"
sexpdir = sdatadir + exp + "/"

SrI =0.039312
mult = 1

Begtime = 49
Endtime = 59

nts1 = []
nwg1s1 = []
nwg2s1 = []
nwg3s1 = []
nwg4s1 = []
nwg5s1 = []
nwg6s1 = []
nwg7s1 = []

s = Stevedatadir + "NumWaveGauge.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            nts1.append(float(row[0]))
            nwg1s1.append(float(row[1])*100)
            nwg2s1.append(float(row[2])*100)
            nwg3s1.append(float(row[3])*100)
            nwg4s1.append(float(row[4])*100)
            nwg5s1.append(float(row[5])*100)
            nwg6s1.append(float(row[6])*100)
            nwg7s1.append(float(row[7])*100)
            
            
         j = j + 1


Numgap = int(mult*SrI /nts1[1] )
BegtimeNumi = int(Begtime/nts1[1]) -1
EndtimeNumi = int(Endtime/nts1[1]) +1

         
nts1 = array(nts1[BegtimeNumi:EndtimeNumi:Numgap])
nwg1s1 = array(nwg1s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg2s1 = array(nwg2s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg3s1 = array(nwg3s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg4s1 = array(nwg4s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg5s1 = array(nwg5s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg6s1 = array(nwg6s1[BegtimeNumi:EndtimeNumi:Numgap])
nwg7s1 = array(nwg7s1[BegtimeNumi:EndtimeNumi:Numgap])
 
NumCom1 = []
NumCom1.append(nts1)
NumCom1.append(nwg1s1)
NumCom1.append(nwg2s1)
NumCom1.append(nwg3s1)
NumCom1.append(nwg4s1)
NumCom1.append(nwg5s1)
NumCom1.append(nwg6s1)
NumCom1.append(nwg7s1)
         



        
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
            nwg1s.append(float(row[2]))
            nwg2s.append(float(row[3]))
            nwg3s.append(float(row[4]))
            nwg4s.append(float(row[5]))
            nwg5s.append(float(row[6]))
            nwg6s.append(float(row[7]))
            nwg7s.append(float(row[8]))
            
            
         j = j + 1


Numgap = int(mult*SrI /nts[1] )
BegtimeNumi = int(Begtime/nts[1]) -1
EndtimeNumi = int(Endtime/nts[1]) +1

         
nts = array(nts[BegtimeNumi:EndtimeNumi:Numgap])
nwg1s = array(nwg1s[BegtimeNumi:EndtimeNumi:Numgap])
nwg2s = array(nwg2s[BegtimeNumi:EndtimeNumi:Numgap])
nwg3s = array(nwg3s[BegtimeNumi:EndtimeNumi:Numgap])
nwg4s = array(nwg4s[BegtimeNumi:EndtimeNumi:Numgap])
nwg5s = array(nwg5s[BegtimeNumi:EndtimeNumi:Numgap])
nwg6s = array(nwg6s[BegtimeNumi:EndtimeNumi:Numgap])
nwg7s = array(nwg7s[BegtimeNumi:EndtimeNumi:Numgap])
 
 
NumCom = []
NumCom.append(nts)
NumCom.append(nwg1s)
NumCom.append(nwg2s)
NumCom.append(nwg3s)
NumCom.append(nwg4s)
NumCom.append(nwg5s)
NumCom.append(nwg6s)
NumCom.append(nwg7s)

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

Egap = mult
BegtimeExpi = int(Begtime/ets[1]) -1
EndtimeExpi = int(Endtime/ets[1]) +1
 
ets = array(ets[BegtimeExpi:EndtimeExpi:Egap])
ewg1s = array(ewg1s[BegtimeExpi:EndtimeExpi:Egap])
ewg2s = array(ewg2s[BegtimeExpi:EndtimeExpi:Egap])
ewg3s = array(ewg3s[BegtimeExpi:EndtimeExpi:Egap])
ewg4s = array(ewg4s[BegtimeExpi:EndtimeExpi:Egap])
ewg5s = array(ewg5s[BegtimeExpi:EndtimeExpi:Egap])
ewg6s = array(ewg6s[BegtimeExpi:EndtimeExpi:Egap])
ewg7s = array(ewg7s[BegtimeExpi:EndtimeExpi:Egap])

        
ExpCom = []
ExpCom.append(ets)
ExpCom.append(ewg1s)
ExpCom.append(ewg2s)
ExpCom.append(ewg3s)
ExpCom.append(ewg4s)
ExpCom.append(ewg5s)
ExpCom.append(ewg6s)
ExpCom.append(ewg7s)
nc = len(NumCom)


for j in range(1,nc):
    sdir = sexpdir +"WaveGauge" + str(j) + "/"
    if not os.path.exists(sdir):
        os.makedirs(sdir)
    nSW= len(NumCom1[0])    
    s = sdir + "NumericalSWW.dat"
    with open(s,'w') as file1:
        for i in range(nSW):
            s ="%3.8f%5s%1.15f\n" %(NumCom1[0][i]," ",NumCom1[j][i]*100)
            file1.write(s)
    nn = len(nts)    
    s = sdir + "NumericalSerre.dat"
    with open(s,'w') as file1:
        for i in range(nn):
            s ="%3.8f%5s%1.15f\n" %(NumCom[0][i]," ",NumCom[j][i]*100)
            file1.write(s)
    ne = len(ets)        
    s = sdir + "Experimental.dat"
    with open(s,'w') as file1:
        for i in range(ne):
            s ="%3.8f%5s%1.15f\n" %(ExpCom[0][i]," ",ExpCom[j][i]*100)
            file1.write(s) 

