# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def lineardiscP(x,x0,x1):
    return ((x - S(x1)) / (S(x0) - S(x1)) ,S(x0), S(x1))

def lineardiscM(x,x0,x1):
    return ((x - S(x0)) / (S(x1) - S(x0)) ,S(x0), S(x1))

        
def quadcontE(x,xm2,xm1,x0,x1,x2):
    fexp = (x - S(xm2))*(x - S(xm1)) / ((S(x0) - S(xm2))*(S(x0) - S(xm1)))
    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (fexp,S(xm2),S(x0),sexp,S(x0),S(x2))


def quadcontM(x,xm1,x0,x1):
    fexp = (x - S(xm1))*(x - S(x1)) / ((S(x0) - S(xm1))*(S(x0) - S(x1)))

    return (fexp,S(xm1),S(x1))

def phiMderiv(x,phi):
    return (diff(phi[0],x),phi[1],phi[2])

def phiEderiv(x,phi):
    return (diff(phi[0],x),phi[1],phi[2],diff(phi[3],x),phi[4],phi[5]  )

x = Symbol('x')


xjm3o2 = "-1"
xjm1 = "-2/3"
xjm1o2 = "-1/3"
xj = "0"
xjp1o2 = "1/3"
xjp1 = "2/3"
xjp3o2 = "1"


phijm1o2 = quadcontE(x,xjm3o2,xjm1,xjm1o2,xj,xjp1o2)
phij = quadcontM(x,xjm1o2,xj,xjp1o2)
phijp1o2 = quadcontE(x,xjm1o2,xj,xjp1o2,xjp1,xjp3o2)


#derivatives
dphijm1o2 = phiEderiv(x,phijm1o2)
dphij = phiMderiv(x,phij )
dphijp1o2 = phiEderiv(x,phijp1o2)


wjm1o2p = lineardiscP(x,xjm1o2,xjp1o2) 
wjp1o2m = lineardiscM(x,xjm1o2,xjp1o2) 


phiEs = [phijm1o2,phijp1o2]
phiMs = [phij]

lphiEs = [phijm1o2,phijp1o2]
lphiMs = [phij]

ldphiEs = [phijm1o2,phijp1o2]
ldphiMs = [phij]

dphiEs = [dphijm1o2,dphijp1o2]
dphiMs = [dphij]

ws = [wjm1o2p,wjp1o2m]
wsNAME = ["wjm1o2p","wjp1o2m"]
phiEsNAME = ["phijm1o2","phijp1o2"]
phiMsNAME = ["phij"]
lphiEsNAME = ["phijm1o2","phijp1o2"]
lphiMsNAME = ["phij"]
dphiEsNAME = ["dphijm1o2","dphijp1o2"]
dphiMsNAME = ["dphij"]
ldphiEsNAME = ["dphijm1o2","dphijp1o2"]
ldphiMsNAME = ["dphij"]

## Integral of G

#First cell cjm1
nwis = len(ws)
nphiEs = len(phiEs)
nlphiEs = len(lphiEs)
nphiMs = len(phiMs)
nlphiMs = len(lphiMs)
ndphiEs = len(dphiEs)
ndphiMs = len(dphiMs)
nldphiEs = len(ldphiEs)
nldphiMs = len(ldphiMs)

Gintegrallist = []
for i in range(nwis):
    iw = ws[i][0]
    iwN = wsNAME[i]
    lbi = ws[i][1]
    ubi = ws[i][2]
    for j in range(nlphiEs):
        if(lphiEs[j][1] == lbi and lphiEs[j][2] == ubi):
            iphi = lphiEs[j][0]
            iphiN = lphiEsNAME[j]
            Intiphiiw = integrate(iphi*iw,(x,lbi,ubi))
            Gintegrallist.append( (iw,iphi ,lbi,ubi,Intiphiiw,iwN,iphiN ) )
            #print(lbi ,ubi , iphi,iw,Intiphiiw   )
            #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )
            
        elif(lphiEs[j][4] == lbi and lphiEs[j][5] == ubi):
            iphi = lphiEs[j][3]
            iphiN = lphiEsNAME[j]
            Intiphiiw = integrate(iphi*iw,(x,lbi,ubi))
            Gintegrallist.append(  (iw,iphi ,lbi,ubi,Intiphiiw,iwN,iphiN ))
            #print(lbi ,ubi , iphi,iw,Intiphiiw   )
            #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )
            
    for j in range(nlphiMs):
        if(lphiMs[j][1] == lbi and lphiMs[j][2] == ubi):
            iphi = lphiMs[j][0]
            iphiN = lphiMsNAME[j]
            Intiphiiw = integrate(iphi*iw,(x,lbi,ubi))
            Gintegrallist.append(  (iw,iphi ,lbi,ubi,Intiphiiw,iwN,iphiN ) )
            #print(lbi ,ubi , iphi,iw,Intiphiiw   )
            #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )

s = "GintegralsE.csv"
n = len(Gintegrallist)
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['W Term' ,'Phi Term','lower bound of integral','upper bound of integral', 'Value of Integral', 'Name of W Term' , 'name of phi Term'])        
               
     for i in range(n):
         writefile2.writerow([str(Gintegrallist[i][0]),str(Gintegrallist[i][1]),str(Gintegrallist[i][2]),str(Gintegrallist[i][3]),str(Gintegrallist[i][4]),str(Gintegrallist[i][5]),str(Gintegrallist[i][6])])  


s = "GintegralsE.txt"
with open(s,'w') as file1:
    s = "%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s\n" %('W Term',' ' ,'Phi Tem',' ','lbi',' ','ubi',' ', 'vi',' ', 'W name',' ' , 'Phi name')
    file1.write(s) 
    for i in range(n):
        s = "%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s\n" %(str(Gintegrallist[i][0])," ",str(Gintegrallist[i][1])," ",str(Gintegrallist[i][2])," ",str(Gintegrallist[i][3])," ",str(Gintegrallist[i][4])," ",str(Gintegrallist[i][5])," ",str(Gintegrallist[i][6]))
        file1.write(s)    
        
        
#### uh integrals
uhintegrallist = []       
for i in range(nwis):
    iw = ws[i][0]
    iwN = wsNAME[i]
    lbi = ws[i][1]
    ubi = ws[i][2]
    for j in range(nphiEs):
        if(phiEs[j][1] == lbi and phiEs[j][2] == ubi):
            jphiN = phiEsNAME[j]
            jphi = phiEs[j][0]
            for k in range(nlphiEs):
                if(lphiEs[k][1] == lbi and lphiEs[k][2] == ubi):
                    kphi = lphiEs[k][0]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
                elif(lphiEs[k][4] == lbi and lphiEs[k][5] == ubi):
                    kphi = lphiEs[k][3]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
            for k in range(nlphiMs):
                if(lphiMs[k][1] == lbi and lphiMs[k][2] == ubi):
                    kphi = lphiMs[k][0]
                    kphiN = lphiMsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iphi,iw,Intiphiiw   )
                    #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )
        if(phiEs[j][4] == lbi and phiEs[j][5] == ubi):
            jphiN = phiEsNAME[j]
            jphi = phiEs[j][3]
            for k in range(nlphiEs):
                if(lphiEs[k][1] == lbi and lphiEs[k][2] == ubi):
                    kphi = lphiEs[k][0]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
                elif(lphiEs[k][4] == lbi and lphiEs[k][5] == ubi):
                    kphi = lphiEs[k][3]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
            for k in range(nlphiMs):
                if(lphiMs[k][1] == lbi and lphiMs[k][2] == ubi):
                    kphi = lphiMs[k][0]
                    kphiN = lphiMsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iphi,iw,Intiphiiw   )
                    #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )  
                    
    for j in range(nphiMs):
        if(phiMs[j][1] == lbi and phiMs[j][2] == ubi):
            jphiN = phiMsNAME[j]
            jphi = phiMs[j][0]
            
            for k in range(nlphiEs):
                if(lphiEs[k][1] == lbi and lphiEs[k][2] == ubi):
                    kphi = lphiEs[k][0]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
                elif(lphiEs[k][4] == lbi and lphiEs[k][5] == ubi):
                    kphi = lphiEs[k][3]
                    kphiN = lphiEsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    
            for k in range(nlphiMs):
                if(lphiMs[k][1] == lbi and lphiMs[k][2] == ubi):
                    kphi = lphiMs[k][0]
                    kphiN = lphiMsNAME[k]
                    Intiphiiw = integrate(iw*jphi*kphi,(x,lbi,ubi))
                    uhintegrallist.append((iw , jphi ,kphi ,lbi,ubi,Intiphiiw,iwN,jphiN,kphiN ))
                    #print(lbi ,ubi , iw ,jphi,kphi,Intiphiiw   )
                    #print(lbi ,ubi ,iwN, jphiN,kphiN,Intiphiiw   )
                    #Gintegrallist.append( (iphi , iw ,lbi,ubi,Intiphiiw,iphiN,iwN ) )
                    #print(lbi ,ubi , iphi,iw,Intiphiiw   )
                    #print(lbi ,ubi , iphiN,iwN,Intiphiiw   )  


s = "uhintegralsE.csv"
n = len(uhintegrallist)
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['W Term' ,'first Phi Term','second Phi Term','lower bound of integral','upper bound of integral', 'Value of Integral', 'Name of W Term' , 'name of first phi Term', 'name of second phi term'])        
               
     for i in range(n):
         writefile2.writerow([str(uhintegrallist[i][0]),str(uhintegrallist[i][1]),str(uhintegrallist[i][2]),str(uhintegrallist[i][3]),str(uhintegrallist[i][4]),str(uhintegrallist[i][5]),str(uhintegrallist[i][6]),str(uhintegrallist[i][7]),str(uhintegrallist[i][8])])  


s = "uhintegralsE.txt"
with open(s,'w') as file1:
    s = "%30s%5s%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s%5s%10s\n" %('W Term',' ' ,'1 Phi Tem',' ' ,'2 Phi Tem',' ','lbi',' ','ubi',' ', 'vi',' ', 'W name',' ' , '1Phi name',' ' , '2Phi name')
    file1.write(s) 
    for i in range(n):
        s = "%30s%5s%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s%5s%10s\n" %(str(uhintegrallist[i][0])," ",str(uhintegrallist[i][1])," ",str(uhintegrallist[i][2])," ",str(uhintegrallist[i][3])," ",str(uhintegrallist[i][4])," ",str(uhintegrallist[i][5])," ",str(uhintegrallist[i][6])," ",str(uhintegrallist[i][7])," ",str(uhintegrallist[i][8]))
        file1.write(s) 
                    
#### h^3 u_x phi_x integrals
h2uxintegrallist = [] 

for i in range(nwis):
    iw = ws[i][0]
    iwN = wsNAME[i]
    lbi = ws[i][1]
    ubi = ws[i][2]
    for l in range(nwis):
        if(ws[l][1] == lbi and ws[l][2] == ubi):
            lw = ws[l][0]
            lwN = wsNAME[l]
            for m in range(nwis):
                if(ws[m][1] == lbi and ws[m][2] == ubi):
                    mw = ws[m][0]
                    mwN = wsNAME[m]
                    
                    for j in range(ndphiEs):
                        if(dphiEs[j][1] == lbi and dphiEs[j][2] == ubi):
                            jdphiN = dphiEsNAME[j]
                            jdphi = dphiEs[j][0]
                            for k in range(nldphiEs):
                                if(ldphiEs[k][1] == lbi and ldphiEs[k][2] == ubi):
                                    kdphi = ldphiEs[k][0]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
           
                                    
                                elif(ldphiEs[k][4] == lbi and ldphiEs[k][5] == ubi):
                                    kdphi = ldphiEs[k][3]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                            for k in range(nldphiMs):
                                if(ldphiMs[k][1] == lbi and ldphiMs[k][2] == ubi):
                                    kdphi = ldphiMs[k][0]
                                    kdphiN = ldphiMsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                        elif(dphiEs[j][4] == lbi and dphiEs[j][5] == ubi):
                            jdphiN = dphiEsNAME[j]
                            jdphi = dphiEs[j][3]
                            for k in range(nldphiEs):
                                if(ldphiEs[k][1] == lbi and ldphiEs[k][2] == ubi):
                                    kdphi = ldphiEs[k][0]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                                elif(ldphiEs[k][4] == lbi and ldphiEs[k][5] == ubi):
                                    kdphi = ldphiEs[k][3]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                            for k in range(nldphiMs):
                                if(ldphiMs[k][1] == lbi and ldphiMs[k][2] == ubi):
                                    kdphi = ldphiMs[k][0]
                                    kdphiN = ldphiMsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                    for j in range(ndphiMs):
                        if(dphiMs[j][1] == lbi and dphiMs[j][2] == ubi):
                            jdphiN = dphiMsNAME[j]
                            jdphi = dphiMs[j][0]
                            
                            for k in range(nldphiEs):
                                if(ldphiEs[k][1] == lbi and ldphiEs[k][2] == ubi):
                                    kdphi = ldphiEs[k][0]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                                elif(ldphiEs[k][4] == lbi and ldphiEs[k][5] == ubi):
                                    kdphi = ldphiEs[k][3]
                                    kdphiN = ldphiEsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)
                                    
                            for k in range(nldphiMs):
                                if(ldphiMs[k][1] == lbi and ldphiMs[k][2] == ubi):
                                    kdphi = ldphiMs[k][0]
                                    kdphiN = ldphiMsNAME[k]
                                    Intiphiiw = integrate(iw*lw*mw*jdphi*kdphi,(x,lbi,ubi))
                                    h2uxintegrallist.append((iw,lw,mw, jdphi ,kdphi ,lbi,ubi,Intiphiiw,iwN,lwN,mwN,jdphiN,kdphiN ))
                                    print(iwN,lwN,mwN,jdphiN,kdphiN,lbi,ubi)

s = "h3uxintegralsE.csv"
n = len(h2uxintegrallist)
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['first W Term', 'second W term', 'third W term' ,'first Phi Term','second Phi Term','lower bound of integral','upper bound of integral', 'Value of Integral', 'Name of first W Term', 'Name of second W Term', 'Name of third W Term' , 'name of first phi Term', 'name of second phi term'])        
               
     for i in range(n):
         writefile2.writerow([str(h2uxintegrallist[i][0]),str(h2uxintegrallist[i][1]),str(h2uxintegrallist[i][2]),str(h2uxintegrallist[i][3]),str(h2uxintegrallist[i][4]),str(h2uxintegrallist[i][5]),str(h2uxintegrallist[i][6]),str(h2uxintegrallist[i][7]),str(h2uxintegrallist[i][8]),str(h2uxintegrallist[i][9]),str(h2uxintegrallist[i][10]),str(h2uxintegrallist[i][11]),str(h2uxintegrallist[i][12])])  


s = "h3uxintegralsE.txt"
with open(s,'w') as file1:
    s = "%30s%5s%30s%5s%30s%5s%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s%5s%10s%5s%10s%5s%10s\n" %('1 W Term',' ','2 W Term',' ','3 W Term',' ','1 Phi Tem',' ' ,'2 Phi Tem',' ','lbi',' ','ubi',' ', 'vi',' ', '1W name',' ', '2W name',' ', '3W name',' ' , '1Phi name',' ' , '2Phi name')
    file1.write(s) 
    for i in range(n):
        s = "%30s%5s%30s%5s%30s%5s%30s%5s%30s%5s%5s%5s%5s%5s%5s%5s%10s%5s%10s%5s%10s%5s%10s%5s%10s\n" %(str(h2uxintegrallist[i][0])," ",str(h2uxintegrallist[i][1])," ",str(h2uxintegrallist[i][2])," ",str(h2uxintegrallist[i][3])," ",str(h2uxintegrallist[i][4])," ",str(h2uxintegrallist[i][5])," ",str(h2uxintegrallist[i][6])," ",str(h2uxintegrallist[i][7])," ",str(h2uxintegrallist[i][8])," ",str(h2uxintegrallist[i][9])," ",str(h2uxintegrallist[i][10])," ",str(h2uxintegrallist[i][11])," ",str(h2uxintegrallist[i][12]))
        file1.write(s) 
