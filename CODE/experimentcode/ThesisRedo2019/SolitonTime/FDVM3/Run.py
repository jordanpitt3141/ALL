from scipy import *
import csv
import os
from Serre3 import *
import time

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b 
    

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
 
def SolitonMass(a0,a1,k,xb,xe):    
    return a0*(xe - xb) + a1*(tanh(k*xe) - tanh(k*xb)) / k
    
def SolitonG(a0,a1,c,k,xb,xe):
    return a1*c / (3*k) *( (3 + 2*a0**2*k**2*sech(k*xe)**2 + 2*a0*a1*k**2*sech(k*xe)**4)*tanh(k*xe) \
   -(3 + 2*a0**2*k**2*sech(k*xb)**2 + 2*a0*a1*k**2*sech(k*xb)**4)*tanh(k*xb) )
 
    
def solitonAveragest0(a0,a1,g,x,dx):
    n = len(x)
    ha = zeros(n)
    Ga = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    idx = 1.0 / dx
    for i in range(n):
        ha[i] = idx*SolitonMass(a0,a1,k,x[i] - 0.5*dx,x[i] + 0.5*dx)
        Ga[i] = idx*SolitonG(a0,a1,c,k,x[i] - 0.5*dx,x[i] + 0.5*dx)
        
    return ha,Ga


#Soliton Problem    
wdirb = "/home/jp/Documents/PhD/project/data/2019/TimeSoliton/FDVM3/"

if not os.path.exists(wdirb):
    os.makedirs(wdirb)

TotTimes = []

trials = 10
for nt in range(trials):
    
    ki = 10
    start_time = time.time()
    wdir = wdirb + str(ki) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    c = sqrt(g*(a0 + a1))
    
    
    startx = -250
    endx = 250
    
    startt = 0.0
    endt = 50
    
    
    dx = 100.0/ 2**ki
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    
    t = startt
    
    x = arange(startx,endx +0.1*dx, dx)
    
    n = len(x)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
    
    #hm,um,Gm = solitoninit(a0,a1,g,x,0,dx)
    
    umbeg = zeros(niBC)
    umend = zeros(niBC)
    hmbeg = a0*ones(niBC)
    hmend = a0*ones(niBC) 
    
    #need averages
    ha,Ga = solitonAveragest0(a0,a1,g,x,dx)
    
    cnBC = niBC - nGsBC
    
    Gabeg = zeros(cnBC)
    Gaend = zeros(cnBC)
    uabeg = zeros(cnBC)
    uaend = zeros(cnBC)
    habeg = a0*ones(cnBC)
    haend = a0*ones(cnBC)
    
    
    #C arrays
    
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
    
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
    
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
    
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
    
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
    
    
    stepcount = 0
    while t < endt :        
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)

        print (t)
        
        t = t + dt
        stepcount = stepcount + 1
 
    um_c = mallocPy(n)
    Gm_c = mallocPy(n)
    hm_c = mallocPy(n)
    
    ca2midpt(ha_c,dx,n,hm_c)
    ca2midpt(Ga_c,dx,n,Gm_c)
    ufromGh(Gm_c,hm_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, um_c)
    
    umF = copyarrayfromC(um_c,n)
    GmF = copyarrayfromC(Gm_c,n)
    hmF = copyarrayfromC(hm_c,n)


    GaF = copyarrayfromC(Ga_c,n)
    haF = copyarrayfromC(ha_c,n)    


    s = wdir +  "outlast.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'h', 'G' , 'u(m/s)','ha', 'Ga' ])        
                   
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t),str(x[j]), str(hmF[j]) , str(GmF[j]) , str(umF[j]), str(haF[j]), str(GaF[j])])

    end_time = time.time() 
    
    deallocPy(Ga_c)
    deallocPy(Gabeg_c)
    deallocPy(Gaend_c)
    deallocPy(ha_c)
    deallocPy(habeg_c)
    deallocPy(haend_c)
    deallocPy(uabeg_c)
    deallocPy(uaend_c)

    deallocPy(hmbeg_c)
    deallocPy(hmend_c)
    deallocPy(umbeg_c)
    deallocPy(umend_c)    
    
    deallocPy(hm_c)
    deallocPy(um_c)
    deallocPy(Gm_c)  

    
    TotTime = end_time - start_time           

    s = wdir +  "TimeInfo.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['st' ,'et','et-st',"stepcount" ,'(et-  st) /stepcount ' ])        
        writefile2.writerow([str(start_time),str(end_time),str(TotTime),str(stepcount),str(TotTime / stepcount)])
        
    TotTimes.append(TotTime)
    
avgTotTime = sum(TotTimes)/len(TotTimes)
avgTotTimeperStep = avgTotTime/stepcount

s = wdirb +  "TimeInfo.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['# trials''avg time','step' ,'avg time per step'])        
    writefile2.writerow([str(len(TotTimes)),str(avgTotTime),str(stepcount),str(avgTotTimeperStep)])
