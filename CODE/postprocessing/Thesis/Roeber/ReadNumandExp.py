import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os



expdir = "../../../../../data/Experimental/HIreef/Trial8/"

wdir = "../../../../../data/raw/Thesis/Experiment/Roeber/Trial8/FEVM/run4/"
sdirb = "../../../../../data/ThesisPost/Experiment/Roeber/Trial8/FEVM/"
sdirexp =sdirb + "Exp/"
sdirnum =sdirb + "Num/"

if not os.path.exists(sdirexp):
   os.makedirs(sdirexp)
   
if not os.path.exists(sdirnum):
   os.makedirs(sdirnum)
   
texpsall = []
WGsall = []

for i in range(1,15):
    WGs = []
    texps =[]
    s = expdir + "WG" + str(i)+".txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in readfile: 
                texps.append(float(row[0]))
                WGs.append(float(row[1]))
    texpsall.append(texps)
    WGsall.append(WGs)
    

s = wdir + "nWGs.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    t = []    
    nwg1 = []
    nwg2 = []
    nwg3 = []
    nwg4 = []
    nwg5 = []
    nwg6 = []
    nwg7 = []
    nwg8 = []
    nwg9 = []
    nwg10 = []
    nwg11 = []
    nwg12 = []
    nwg13 = []
    nwg14 = []
    j = -1
    for row in readfile:  
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            t.append(float(row[0]))
            nwg1.append(float(row[1]))
            nwg2.append(float(row[2]))
            nwg3.append(float(row[3]))
            nwg4.append(float(row[4]))
            nwg5.append(float(row[5]))
            nwg6.append(float(row[6]))
            nwg7.append(float(row[7]))     
            nwg8.append(float(row[8]))
            nwg9.append(float(row[9]))
            nwg10.append(float(row[10]))
            nwg11.append(float(row[11]))
            nwg12.append(float(row[12]))
            nwg13.append(float(row[13]))
            nwg14.append(float(row[14])) 
        j = j + 1

dt = t[1] -t[0]
incr = 10
st = 20
et = 80

sti =int(st/dt) 
eti =int(et/dt) 

rt = t[sti : eti :incr ]
rnwg1 = nwg1[sti : eti :incr ]
rnwg2 = nwg2[sti : eti :incr ]
rnwg3 = nwg3[sti : eti :incr ]
rnwg4 = nwg4[sti : eti :incr ]
rnwg5 = nwg5[sti : eti :incr ]
rnwg6 = nwg6[sti : eti :incr ]
rnwg7 = nwg7[sti : eti :incr ]
rnwg8 = nwg8[sti : eti :incr ]
rnwg9 = nwg9[sti : eti :incr ]
rnwg10 = nwg10[sti : eti :incr ]
rnwg11 = nwg11[sti : eti :incr ]
rnwg12 = nwg12[sti : eti :incr ]
rnwg13 = nwg13[sti : eti :incr ]
rnwg14 = nwg14[sti : eti :incr ]

n = len(rt)
s = sdirnum + "nWG1.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg1[j])
        file1.write(ss)

s = sdirnum + "nWG2.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg2[j])
        file1.write(ss)

s = sdirnum + "nWG3.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg3[j])
        file1.write(ss)

s = sdirnum + "nWG4.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg4[j])
        file1.write(ss)

s = sdirnum + "nWG5.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg5[j])
        file1.write(ss)

s = sdirnum + "nWG6.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg6[j])
        file1.write(ss)

s = sdirnum + "nWG7.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg7[j])
        file1.write(ss)

s = sdirnum + "nWG8.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg8[j])
        file1.write(ss)
        
s = sdirnum + "nWG9.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg9[j])
        file1.write(ss)

s = sdirnum + "nWG10.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg10[j])
        file1.write(ss)

s = sdirnum + "nWG11.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg11[j])
        file1.write(ss)

s = sdirnum + "nWG12.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg12[j])
        file1.write(ss)

s = sdirnum + "nWG13.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg13[j])
        file1.write(ss)

s = sdirnum + "nWG14.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnwg14[j])
        file1.write(ss)

n = len(texpsall[0])
for i in range(0,14):
    
    s = sdirexp + "WG" +str(i +1)+".dat"
    with open(s,'w') as file1:
        for j in range(n):
            ss ="%3.8f%5s%1.15f\n" %(texpsall[i][j]," ",WGsall[i][j])
            file1.write(ss)
