# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hxt(x,t,a,b,c,d,e,f):
    return (a*x**2 + b*x + c)*(d*t**2 + e*t + f)

def uxt(x,t,a,b,c,d,e,f):
    return (a*x**2 + b*x + c)*(d*t**2 + e*t + f)
    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def Masseq(h,u,x,t):
    
    return diff(h,t) + diff(h*u,x)

def Momeeq(h,u,g,x,t):
    ux = diff(u,x)
    uxx = diff(ux,x)
    uxt = diff(ux,t)
    phi = ux**2 - u*uxx - uxt
    
    return diff(h*u,t) + diff(u**2*h + g/2*h**2 + h**3/3*phi ,x)

def MomeGeq(h,u,G,g,x,t):
    ux = diff(u,x)
    
    return diff(G,t) + diff(G*u + g/2*h**2 - 2*h**3/3*ux*ux ,x)

def RHSFF(H,x,t):
    return integrate(integrate(H,x),t)
        

def PythonToCexp(str1):
    n = len(str1)
    #strn = "return "
    strn= ""
    i = 0
    while i < n:    
        #assuming nor bracketin
        if(str1[i] == "*" and str1[i+1]=="*"):
            
            strsym = ""
            for j1 in range(i-1,-1,-1):
                
                if (str1[j1] == " " or str1[j1] == "+" or str1[j1] == "*"):
                    break
                else:
                    strsym= strsym+ str1[j1]
    
            strnum = ""
            for j2 in range(i+2,n):
                
                if (str1[j2] == " " or str1[j2] == "+" or str1[j2] == "*"):
                    break
                else:
                    strnum= strnum+ str1[j2]
                                    
            strnum = int(strnum)        
            
            strn = strn + (strnum- 1)*("*"+strsym)
            i = j2 
            
        else:
            strn = strn + str1[i]
            i = i+ 1
    
    #return strn + ";"
    return strn

x = Symbol('x')
t = Symbol('t')
g = Symbol('g')

a0 = Symbol('a0')
a1 = Symbol('a1')
a2= Symbol('a2')
a3 = Symbol('a3')
a4 = Symbol('a4')
a5 = Symbol('a5')

b0 = Symbol('b0')
b1 = Symbol('b1')
b2 = Symbol('b2')
b3 = Symbol('b3')
b4 = Symbol('b4')
b5 = Symbol('b5')


hxt = hxt(x,t,1.1,1.2,1.3,1.4,1.5,1.6)
uxt = uxt(x,t,0.11,0.12,0.13,0.14,0.15,0.16)

Gxt = Gxt(hxt,uxt)

barhxt = integrate(hxt,x)

barGxt = integrate(Gxt,x)

Fhint = integrate(hxt*uxt,t)

uxtx = diff(uxt,x)

FGint = integrate(Gxt*uxt + g/2.0*hxt**2 - 2*hxt**3*uxtx*uxtx/3.0,t)

str1 = str(barhxt)
str2 = str(barGxt)
str3 = str(Fhint)
str4 = str(FGint)

str1n = PythonToCexp(str1)
str2n = PythonToCexp(str2)
str3n = PythonToCexp(str3)
str4n = PythonToCexp(str4)

#no original coefficients are small, so if number comes out at floating point error ish, then the expressions are the same
print("hbar correct?")
print(simplify(S(str1n) -barhxt).subs(x,1).subs(t,1))
print("Gbar correct?")
print(simplify(S(str2n) -barGxt).subs(x,1).subs(t,1))
print("h flux correct?")
print(simplify(S(str3n) -Fhint).subs(x,1).subs(t,1))
print("G flux correct?")
print(simplify(S(str4n) -FGint).subs(x,1).subs(t,1))



print("Code")
print
print
print('h bar')
print(str1n)
print
print
print('G bar')
print(str2n)
print
print
print('h Flux')
print(str3n)
print
print
print('G Flux')
print(str4n)

        
        #read symbol infront
        #read number after
    #print(str1[i])


