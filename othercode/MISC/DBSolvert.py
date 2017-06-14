# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:41:33 2016

@author: jp
"""

h1 = 1.8
h0 = 1.0
h2 = 1.36898

g = 9.81
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks

sdir = "./paper/"
	
def analytical_sol(x,h0,h1,h2,t):
    n = len(x)    # number of cells

    u = zeros(n)
    h = zeros(n)
    S2 = ((2.0*h2)/(h2 - h0)) *(sqrt(g*h1) - sqrt(g*h2))        
    u2 = 2*(sqrt(g*h1) - sqrt(g*h2))    
    
    for i in range(n):
        # Calculate Analytical Solution at time t > 0
        u3 = 2.0/3.0*(sqrt(g*h1)+x[i]/t)
        h3 = 4.0/(9.0*g)*(sqrt(g*h1)-x[i]/(2.0*t))*(sqrt(g*h1)-x[i]/(2.0*t))
        
        
        if ( x[i] <= -sqrt(g*h1)*t):
            u[i] = 0.0
            h[i] = h1
        elif ( x[i] > -(sqrt(g*h1)*t) and x[i] <= t*(u2  - sqrt(g*h2)) ):
            u[i] = u3
            h[i] = h3
        elif ( x[i] > t*(u2  - sqrt(g*h2)) and x[i] < t*S2 ):
            u[i] = u2
            h[i] = h2
        elif ( x[i] >= t*S2 ):
            u[i] = 0.0
            h[i] = h0
            
    print(-sqrt(g*h1)*t + 500 )
    print(-sqrt(g*h1)*t + 500 , t*(u2  - sqrt(g*h2)) + 500)
    print(t*(u2  - sqrt(g*h2)) + 500 ,  t*S2 + 500 )
        
    return h , u*h, u
     
x = arange(-500, 500, 0.1)

h,uh,u = analytical_sol(x,h0,h1,h2,30)

plot(x + 500,h,'k')
ylabel("$h$ ($m$)")
xlabel("$x$ ($m$)")
ylim ([1,2])
xlim([0,1000])
n = len(x)        
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i] + 500," ",h[i])
        file1.write(s)
