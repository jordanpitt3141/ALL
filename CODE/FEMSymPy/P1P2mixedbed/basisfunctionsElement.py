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

        
def quadcontEM(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x0),S(x2))

def quadcontEP(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x2),S(x0))


def quadcontM(x,xm1,x0,x1):
    fexp = (x - S(xm1))*(x - S(x1)) / ((S(x0) - S(xm1))*(S(x0) - S(x1)))

    return (fexp,S(xm1),S(x1))

def phideriv(x,phi):
        return (diff(phi[0],x),phi[1],phi[2])        


def Finbounds(phi,lb,ub):
    return phi[1] == lb and phi[2] == ub


def ExpressionConvert(s,hyperlist,x,Cterm,Ctermname,Clb,Cub,filen):
    # hpyerlist has is a list of tuples of (char,functionlist,functionname)
    n = len(hyperlist)
    if (len(s) == 0):
        #can do the interesting stuff here, write out to file
        Intv= integrate(Cterm,(x,Clb,Cub))
        print(Cterm, Ctermname, Clb, Cub)
        writefile2.writerow([str(Cterm),str(Clb),str(Cub), str(Intv)] + Ctermname)        

    elif(Cterm is None):        
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                for j in range(m):
                    ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][0],[hyperlist[i][2][j]],Clb,Cub,filen)
    else:
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                for j in range(m):
                    ExpressionConvert(s[:-1],hyperlist,x,Cterm*hyperlist[i][1][j][0],[hyperlist[i][2][j]] + Ctermname,Clb,Cub,filen)
                    
       

    
            
           
    
    

x = Symbol('x')


xjm1o2 = "-1"
xj = "0"
xjp1o2 = "1"

phijm1o2 = quadcontEM(x,xjm1o2,xj,xjp1o2)
phij = quadcontM(x,xjm1o2,xj,xjp1o2)
phijp1o2 = quadcontEP(x,xjp1o2,xj,xjm1o2)


#derivatives
dphijm1o2 = phideriv(x,phijm1o2)
dphij = phideriv(x,phij )
dphijp1o2 = phideriv(x,phijp1o2)


wjm1o2p = lineardiscP(x,xjm1o2,xjp1o2) 
wjp1o2m = lineardiscM(x,xjm1o2,xjp1o2) 

wsNAME = ["wjm1o2p","wjp1o2m"]
phisNAME = ["phijm1o2","phij","phijp1o2"]
dphisNAME = ["dphijm1o2","dphij","dphijp1o2"]
ws = [wjm1o2p,wjp1o2m]
phis = [phijm1o2,phij,phijp1o2] 
dphis = [dphijm1o2,dphij,dphijp1o2]

hyperlist = []
hyperlist.append( ('w',ws,wsNAME))
hyperlist.append( ('p',phis,phisNAME))
hyperlist.append( ('d',dphis,dphisNAME))




s = "Gv.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'P Term 1'])        
     ExpressionConvert('wp',hyperlist,x,None,None,-1,1,writefile2)      
   
s = "huv.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1' , 'P Term 1', 'P Term 2'])        
     ExpressionConvert('wpp',hyperlist,x,None,None,-1,1,writefile2)  

s = "h3uxvx.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'W Term 3' , 'D Term 1', 'D Term 2'])        
     ExpressionConvert('wwwdd',hyperlist,x,None,None,-1,1,writefile2)  

s = "h2bxuxv.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'D Term 1', 'D Term 2', 'P Term 1'])        
     ExpressionConvert('wwddp',hyperlist,x,None,None,-1,1,writefile2)  

s = "h2bxuvx.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'D Term 1', 'P Term 1', 'D Term 2'])        
     ExpressionConvert('wwdpd',hyperlist,x,None,None,-1,1,writefile2)  

    
s = "hbx2uv.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1' , 'D Term 1', 'D Term 2', 'P Term 1', 'P Term 2'])        
     ExpressionConvert('wddpp',hyperlist,x,None,None,-1,1,writefile2) 

