# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from SWWE import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time
from scipy.interpolate import interp1d
from scipy import signal
from scipy import sqrt
from numpy.fft import fft
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

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

def minmodpy(a, b, c):
    if((a > 0) and (b>0) and (c>0)):
        return min(a,b,c)
    elif((a < 0) and (b<0) and (c<0)):
        return max(a,b,c)
    else:
        return 0.0
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)



def Dambreak(x,x0,h0,h1):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    b = zeros(n)
    
    for i in range(n):
        if (x[i] < x0):
            h[i] = h0
        else:
            h[i] = h1
            
    return h,u*h,b

def SineWave(A,f,l,hb,x,t):
    il = 1.0 / l
    return hb + A*sin(2*pi*(x*il + f*t))

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    #return y1  + (xi)*(y2 - y1)/(x2 - x1)  
    return y1  + (xi)*(y2 - y0)/(x2 - x0) 


def soloverslope(x,a0,a1,solbeg,solend,slopbeg,slopend,g):
    
	"""
	soloverslope : set up initial conditions for the bed
	    
	Input:
		    
		x   :   array of cell centres
		a0  :   deep height, also still water level for soliton
		a1  :   soliton height above the still water
		solbeg  :  the beginning of the range over which the soliton scheme is defined
		solend  :  the end of the range over which the soliton scheme is defined
		slopbeg :  the beginning of the slope
		g       :  the acceleration due to gravity
	
	Output:
	    
		h    :   array of heights at cell centres
		G    :   array of G values at cell centres
		bed  :   array of bed heights at cell centres
	
	"""
	
	n = len(x)
	h = zeros(n)
	u = zeros(n)
	bed = zeros(n)
    
    
	#speed of the soliton
	c = sqrt(g*(a0 + a1))

	for i in range(n):
		
		#This is the range over which we define a soliton with a bed of 0 beneath it, which is smaller than the constant depth area
		if (x[i] > solbeg and x[i] < solend):
			bed[i] = 0
			h[i] = soliton(x[i],0,g,a0,a1)
			u[i] =  c* ((h[i] - a0) / h[i])

		#We still have a constant bed, but now the height is constant as well
		elif(x[i] <= slopbeg):
			bed[i] = 0
			h[i] = a0 - bed[i]
			u[i] = 0

		#This is the region in which the bed has a linear slope with a constant stage
		elif(x[i] > slopbeg and x[i] <= slopend):
			bed[i] = 0.02*(x[i] - slopbeg)
			h[i] = a0 - bed[i]
			u[i] = 0

		#After the slope the bed is constant as is the stage
		elif(x[i] >= slopend):
			bed[i] = 0.99
			h[i] = a0 - bed[i]
			u[i] = 0  
    
	#return h,G,bed
	return h,u,bed
 
 
def BejiEdge(hm1o2,h0,u0,u1,hb,g,dx,b,hpeak):
    hMbeg = zeros(3)
    uMbeg = zeros(3)
 
    #i=-1/2
    hMbeg[2] = hm1o2  
    uMbeg[2] = sqrt(g*hMbeg[2])*(1  -  hb / hMbeg[2])   
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = sqrt(g*hMbeg[1])*(1  -  hb / hMbeg[1])  
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = sqrt(g*hMbeg[0])*(1  -  hb / hMbeg[0])   
    

    
    return hMbeg,uMbeg

def DingFlume(x,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):

        if(0 <= x[i] <= 6):
            bed[i] = 0.0
            h[i] = 0.4            
        elif(6 < x[i] <= 12):
            bed[i] = 0.05*(x[i] - 6)
            h[i] = 0.4 - bed[i]
        elif(12 < x[i] <= 14):
            bed[i] = 0.3
            h[i] = 0.1
        elif(14 < x[i] <= 17):
            bed[i] = 0.3 - 0.1*(x[i] - 14)
            h[i] = 0.4 - bed[i]
        elif(17 < x[i] <= 18.95):
            bed[i] = 0.0
            h[i] = 0.4 - bed[i]
        elif(18.95 < x[i] <= 23.95):
            bed[i] = 0.04*(x[i] - 18.95)
            h[i] = 0.4  - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            h[i] = 0.4  - bed[i]
        else:
            bed[i] = 0.0
            h[i] = 0.4  - bed[i]
    return h,u,bed

def analytical_sol(x,g,h0,h1,h2,t):
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
            
    print(-sqrt(g*h1)*t)
    print(-sqrt(g*h1)*t , t*(u2  - sqrt(g*h2)))
    print(t*(u2  - sqrt(g*h2)) ,  t*S2)
        
    return h , u*h, u

def Roeberflume(x,xexp,bedexp,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        if(x[i] <= xexp[0]):
            bed[i] = bedexp[1]
            h[i] = 0.0 - bed[i] 
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
                
        elif(xexp[0] < x[i] < xexp[-1]):
            j = [ nin for nin, nv in enumerate(xexp) if nv>=x[i] ][0]
            bed[i] = bedexp[j-1] + ( (bedexp[j] - bedexp[j-1]) / 0.05)*(x[i] - xexp[j-1])
            h[i] = 0.0 - bed[i]
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
            
        elif(x[i] >= xexp[-1]):
            bed[i] = bedexp[-1]
            h[i] = 0.0 - bed[i]
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
            
    return h,u,u,bed
 
"""
## Dambreak accuracy
h1 = 1.0
h0 = 1.8
h2 = 1.36898
x0 = 500

g = 9.81  

normhs = []
normuhs = []
normus = []
dxs = []

wdatadir = "../../data/raw/SWWE/DBANAconshortNEW/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)

s = wdatadir + "norms.csv"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'norm h','norm u' ,"norm uh"])   
  
for it in range(1,3):    
    dx = 1.0 / (2**it)
    l =  0.01
    dt = l*dx
    startx = 300
    endx = 700 + 0.9*dx
    startt = 0.0
    endt = 30 + (dt*0.1)  
            
    szoomx = startx
    ezoomx = endx
    
    wdir = wdatadir + str(it) + "/"
    if not os.path.exists(wdir):
        os.makedirs(wdir)
            
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    
    n = len(x)  
    
    theta = 1.2
    
    gap = int(1.0/dt)
    
    nBC = 3
    nbc = 3*n + 2*(nBC)
    
    idx = 1.0 / dx
            
    h,uh,b = Dambreak(x,x0,h0,h1)   
    
    SWWh,SWWuh,SWWu = analytical_sol(x- 500,g,h1,h0,h2,30)    
    
    hMbeg = h[0]*ones(nBC)
    uhMbeg = uh[0]*ones(nBC)
    hMend = h[-1]*ones(nBC)
    uhMend = uh[-1]*ones(nBC) 
    bMbeg = b[0]*ones(nBC)
    bMend = b[-1]*ones(nBC)
    
    wMbeg = hMbeg + bMbeg
    wMend = hMend + bMend
        
    h_c = copyarraytoC(h)
    uh_c = copyarraytoC(uh)
    b_c = copyarraytoC(b)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend) 
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend) 
    uhMbeg_c = copyarraytoC(uhMbeg)
    uhMend_c = copyarraytoC(uhMend)
    
    
    hhbc_c = mallocPy(nbc)
    whbc_c = mallocPy(nbc)
    uhbc_c = mallocPy(nbc)
    bhbc_c = mallocPy(nbc)
    
    
    #Just an FEM solve here
    for i in range(1,len(t)):
        evolvewrap(uh_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,uhMbeg_c,uhMend_c,bMbeg_c,bMend_c,g,dx,dt,n,nBC,nbc,theta,hhbc_c,whbc_c,bhbc_c,uhbc_c)
        print(t[i])
    
    
    hC = copyarrayfromC(h_c,n)
    uhC = copyarrayfromC(uh_c,n)     
    bC = copyarrayfromC(b_c,n) 

    s = wdir + "outlast.csv"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'h', 'uh' , 'u','h Ana', 'uh Ana', 'u Ana' ])        
               
        for k in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[k]) ,str(hC[k]) , str(uhC[k]) , str(uhC[k]/hC[k]) , str(SWWh[k]), str(SWWuh[k]), str(SWWu[k])])      


    
    normh = norm(SWWh - hC,ord=1)/ norm(SWWh,ord=1)
    normuh = norm(SWWuh - array(uhC),ord=1)/ norm(SWWuh,ord=1)
    normu = norm(SWWu - array(uhC)/array(hC),ord=1)/ norm(SWWu,ord=1)


    s = wdatadir + "norms.csv"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx),str(normh),str(normu),str(normuh)  ])     
    
    
    normhs.append(normh)
    normuhs.append(normuh)
    normus.append(normu)
    dxs.append(dx)
    
    
    deallocPy(h_c)
    deallocPy(b_c)
    deallocPy(uh_c)
    
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(uhbc_c)
    deallocPy(bhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(uhMbeg_c)
    deallocPy(hMend_c)
    deallocPy(uhMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
    deallocPy(bMbeg_c)
    deallocPy(bMend_c)    
 

s = wdatadir + "norms.csv"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    writefile2.writerow(['dx' ,'norm h','norm uh' ,"normu" ])        
           
    for k in range(len(dxs)):
        writefile2.writerow([str(dxs[k]),str(normhs[k]),str(normuhs[k]),str(normus[k])])   
    
##Timing

BEGT = time.time()
h1 = 1.0
h0 = 1.8
x0 = 500

g = 9.81
dx = 0.1
l =  0.01
dt = l*dx
startx = 300
endx = 700 + 0.9*dx
startt = 0.0
endt = 30 + (dt*0.9)  
        
szoomx = startx
ezoomx = endx
        

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  

theta = 1.2

gap = int(1.0/dt)

nBC = 3
nbc = 3*n + 2*(nBC)

idx = 1.0 / dx
        
h,u,b = Dambreak(x,x0,h0,h1)     

hMbeg = h[0]*ones(nBC)
uMbeg = u[0]*ones(nBC)
hMend = h[-1]*ones(nBC)
uMend = u[-1]*ones(nBC) 
bMbeg = b[0]*ones(nBC)
bMend = b[-1]*ones(nBC)

wMbeg = hMbeg + bMbeg
wMend = hMend + bMend
    
h_c = copyarraytoC(h)
u_c = copyarraytoC(u)
b_c = copyarraytoC(b)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)


hhbc_c = mallocPy(nbc)
whbc_c = mallocPy(nbc)
uhbc_c = mallocPy(nbc)
bhbc_c = mallocPy(nbc)


#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrap(u_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,dt,n,nBC,nbc,theta,hhbc_c,whbc_c,bhbc_c,uhbc_c)
    print(t[i])


hC = copyarrayfromC(h_c,n)
uC = copyarrayfromC(u_c,n)     
bC = copyarrayfromC(b_c,n)      

hhbcC = copyarrayfromC(hhbc_c,nbc)
whbcC = copyarrayfromC(whbc_c,nbc)
uhbcC = copyarrayfromC(uhbc_c,nbc)
bhbcC = copyarrayfromC(bhbc_c,nbc)


deallocPy(h_c)
deallocPy(b_c)
deallocPy(u_c)

deallocPy(hhbc_c)
deallocPy(whbc_c)
deallocPy(uhbc_c)
deallocPy(bhbc_c)

deallocPy(hMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(uMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
deallocPy(bMbeg_c)
deallocPy(bMend_c)
ENDT = time.time()
print(ENDT - BEGT)
"""

#Dingemans Flume
"""
wdatadir = "../../../../data/raw/SWWE/Beji/" 
expdir = "../../../../data/Experimental/Data 1994 Paper/CSV/"

exp = "sl"

g = 9.81
Cr = 0.5
l = Cr / (sqrt(g*(0.43) ))
sr = 0.039312
dt = sr/ (2**5)
dx = (0.1/2.0**4)

theta = 2
startx = 5.7
endx = 200
startt = 0
endt = 100 + dt  

hb = 0.4

wdir = wdatadir + exp + "/"


if not os.path.exists(wdir):
    os.makedirs(wdir)

        
szoomx = startx
ezoomx = endx

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.4

gap = int(1.0/dt)


nBC = 3
nbc = 3*n + 2*(nBC)


idx = 1.0 / dx
        
h,uh,b = DingFlume(x,dx) 


hMend = h[-1]*ones(nBC)
bMend = b[-1]*ones(nBC)
uMend = uh[-1]*ones(nBC)
wMend = hMend + bMend

wg2i = int((10.5 - startx) / dx ) + 1 #good one
wg3i = int((12.5 - startx) / dx ) #G
wg4i = int((13.5 - startx) / dx ) #G
wg5i = int((14.5 - startx) / dx ) #
wg6i = int((15.7 - startx) / dx )
wg7i = int((17.3 - startx) / dx )

s = expdir + exp + ".csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    ts = [0.0]
    rs = [0.0]
    wg1s = [0.0]
    wg2s = [0.0]
    wg3s = [0.0]
    wg4s = [0.0]
    wg5s = [0.0]
    wg6s = [0.0]
    wg7s = [0.0]
    j = -1
    for row in readfile:   
        if (j >= 0):
            ts.append((j + 1)*sr)
            rs.append(float(row[0]))
            wg1s.append(float(row[1]))
            wg2s.append(float(row[2]))
            wg3s.append(float(row[3]))
            wg4s.append(float(row[4]))
            wg5s.append(float(row[5]))
            wg6s.append(float(row[6]))
            wg7s.append(float(row[7]))
        j = j + 1

nwg1s = [hb]
nwg2s = [h[wg2i]]
nwg3s = [h[wg3i]]
nwg4s = [h[wg4i]]
nwg5s = [h[wg5i]]
nwg6s = [h[wg6i]]
nwg7s = [h[wg7i]]


ct = dt
mp = int(ct/sr)
ftc0 = hb
ftc1 = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp]) + hb


    
h_c = copyarraytoC(h)
uh_c = copyarraytoC(uh)
b_c = copyarraytoC(b)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend) 
uMend_c = copyarraytoC(uMend)

uhbc_c = mallocPy(nbc)
hhbc_c = mallocPy(nbc)
whbc_c = mallocPy(nbc)
bhbc_c = mallocPy(nbc)

wg2im1h = readfrommem(h_c,wg2i - 1) 
wg2ih = readfrommem(h_c,wg2i) 
wg2ip1h = readfrommem(h_c,wg2i + 1) 
hbwg2 = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],10.5 - x[wg2i])   

wg3im1h = readfrommem(h_c,wg3i - 1) 
wg3ih = readfrommem(h_c,wg3i) 
wg3ip1h = readfrommem(h_c,wg3i + 1) 
hbwg3 = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],12.5 - x[wg3i])   

wg4im1h = readfrommem(h_c,wg4i - 1) 
wg4ih = readfrommem(h_c,wg4i) 
wg4ip1h = readfrommem(h_c,wg4i + 1) 
hbwg4 = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],13.5 - x[wg4i])  

wg5im1h = readfrommem(h_c,wg5i - 1) 
wg5ih = readfrommem(h_c,wg5i) 
wg5ip1h = readfrommem(h_c,wg5i + 1) 
hbwg5 = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],14.5 - x[wg5i])  

wg6im1h = readfrommem(h_c,wg6i - 1) 
wg6ih = readfrommem(h_c,wg6i) 
wg6ip1h = readfrommem(h_c,wg6i + 1) 
hbwg6 = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],15.7 - x[wg6i])  

wg7im1h = readfrommem(h_c,wg7i - 1) 
wg7ih = readfrommem(h_c,wg7i) 
wg7ip1h = readfrommem(h_c,wg7i + 1) 
hbwg7 = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],17.3 - x[wg7i]) 


h0s = [hb]
h1s = [hb]
u0s = [0]
u1s = [0]

#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrapIncomwavetoDir(uh_c, h_c,b_c,ftc0, ftc1 , hMend_c, wMend_c , uMend_c,bMend_c,g,dx,dt,n,nBC,nbc, theta, hhbc_c,whbc_c,bhbc_c,uhbc_c)
    nwg1s.append(ftc0)
    
    wg2im1h = readfrommem(h_c,wg2i - 1) 
    wg2ih = readfrommem(h_c,wg2i) 
    wg2ip1h = readfrommem(h_c,wg2i + 1) 
    wg2h = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],10.5 - x[wg2i])   
    
    wg3im1h = readfrommem(h_c,wg3i - 1) 
    wg3ih = readfrommem(h_c,wg3i) 
    wg3ip1h = readfrommem(h_c,wg3i + 1) 
    wg3h = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],12.5 - x[wg3i])   
    
    wg4im1h = readfrommem(h_c,wg4i - 1) 
    wg4ih = readfrommem(h_c,wg4i) 
    wg4ip1h = readfrommem(h_c,wg4i + 1) 
    wg4h = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],13.5 - x[wg4i])  
    
    wg5im1h = readfrommem(h_c,wg5i - 1) 
    wg5ih = readfrommem(h_c,wg5i) 
    wg5ip1h = readfrommem(h_c,wg5i + 1) 
    wg5h = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],14.5 - x[wg5i])  
    
    wg6im1h = readfrommem(h_c,wg6i - 1) 
    wg6ih = readfrommem(h_c,wg6i) 
    wg6ip1h = readfrommem(h_c,wg6i + 1) 
    wg6h = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],15.7 - x[wg6i])  
    
    wg7im1h = readfrommem(h_c,wg7i - 1) 
    wg7ih = readfrommem(h_c,wg7i) 
    wg7ip1h = readfrommem(h_c,wg7i + 1) 
    wg7h = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],17.3 - x[wg7i]) 
    
    nwg2s.append((wg2h - hbwg2))  
    nwg3s.append((wg3h - hbwg3)) 
    nwg4s.append((wg4h - hbwg4)) 
    nwg5s.append((wg5h - hbwg5)) 
    nwg6s.append((wg6h - hbwg6)) 
    nwg7s.append(wg7h - hbwg7) 
    
    h0t = readfrommem(h_c,0) 
    h1t = readfrommem(h_c,1) 
    uh0t = readfrommem(uh_c,0) 
    uh1t = readfrommem(uh_c,1)
    
    h0s.append(h0t)
    h1s.append(h1t)
    u0s.append(uh0t/h0t)
    u1s.append(uh1t/h1t)
    
    ct = t[i] +dt
    mp = int(ct/sr)
    ftc0 = ftc1 
    ftc1 = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp]) + hb
   
    print(t[i])

hC = copyarrayfromC(h_c,n)
uhC = copyarrayfromC(uh_c,n)   
hhbcC = copyarrayfromC(hhbc_c,nbc)  
whbcC = copyarrayfromC(whbc_c,nbc)
bhbcC = copyarrayfromC(bhbc_c,nbc)    
uhbcC = copyarrayfromC(uhbc_c,nbc)      

s = wdir + "NumWaveGauge.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["Time(s)","hts(m)","WG1(m)","WG2(m)", "WG3(m)","WG4(m)","WG5(m)","WG6(m)","WG7(m)","WG8(m)"]) 
    
    for j in range(len(t)):
        writefile.writerow([str(t[j]), str(nwg1s[j]), str(nwg2s[j]), str(nwg3s[j]), str(nwg4s[j]), str(nwg5s[j]), str(nwg6s[j]),str(nwg7s[j])]) 

s = wdir +  "h0u0.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['x','t' ,'h0','u0'])        
                   
     for j in range(len(h0s)):
         writefile2.writerow([str(x[0]), str(t[j]), str(h0s[j]),u0s[j]])


s = wdir +  "h1u1.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['x','t' ,'h0','u0'])        
                   
     for j in range(len(h0s)):
         writefile2.writerow([str(x[1]), str(t[j]), str(h1s[j]),u1s[j]])
         

   
deallocPy(h_c)
deallocPy(b_c)
deallocPy(uh_c) 
"""


#Roeber Data

nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nwg8s = []
nwg9s = []
nwg10s = []
nwg11s = []
nwg12s = []
nwg13s = []
nwg14s = []
nts = []


wdir = "../../../../data/raw/SWWE/RoeberWGtrial8/" 

expdir = "../../../../data/Experimental/HIreef/Trial8/"

s = expdir + "bed.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    xexp = []
    bedexp = []
    for row in readfile:       
            xexp.append(float(row[0]))
            bedexp.append(float(row[1]))
            
s = expdir + "WG1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG1exp.append(float(row[1]))
            
s = expdir + "WG2.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG2exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG2exp.append(float(row[1]))

s = expdir + "WG3.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG3exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG3exp.append(float(row[1]))

s = expdir + "WG4.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG4exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG4exp.append(float(row[1]))
            
s = expdir + "WG5.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG5exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG5exp.append(float(row[1]))

s = expdir + "WG6.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG6exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG6exp.append(float(row[1]))

s = expdir + "WG7.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG7exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG7exp.append(float(row[1]))
            
s = expdir + "WG8.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG8exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG8exp.append(float(row[1]))

s = expdir + "WG9.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG9exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG9exp.append(float(row[1]))

s = expdir + "WG10.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG10exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG10exp.append(float(row[1]))
            
s = expdir + "WG11.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG11exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG11exp.append(float(row[1]))

s = expdir + "WG12.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG12exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG12exp.append(float(row[1]))

s = expdir + "WG13.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG13exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG13exp.append(float(row[1]))

s = expdir + "WG14.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG14exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG14exp.append(float(row[1]))

a1 = 1.2
a0 = 2.46
hb = 2.46
g = 9.81
sr = 0.02
dt = sr/ (2**4)
dx = 0.05

theta = 1.2
startx = 17.6 + 0.5*dx #17.6
endx = 1000
startt = 0
endt = 119.96 
xWG2 = 28.6040 
xWG3 = 35.9060
xWG4 = 40.5780  
xWG5 = 44.2530  
xWG6 = 46.0930  
xWG7 = 48.2330  
xWG8 = 50.3730 
xWG9 = 54.4060 
xWG10 = 58.0500 
xWG11 = 61.7000
xWG12 = 65.3800  
xWG13 = 72.7200 
xWG14 = 80.0300  



if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.2


nBC = 3
nbc = 3*n + 2*(nBC)

        
h1,u1,uh1,b1 = Roeberflume(x,xexp,bedexp,dx) 
h = h1[0]*ones(n)
b = b1[0]*ones(n)
uh = zeros(n)
u = zeros(n)

hMend = h[-1]*ones(nBC)
uhMend = uh[-1]*ones(nBC)  
wMend = (h[-1] + b[-1])*ones(nBC)
bMend = b[-1]*ones(nBC)

ct = dt
mp = int(ct/sr)
ftc0 = hb
ftc1 = lineinterp(WG1exp[mp],WG1exp[mp + 1],texp[mp],texp[mp + 1],ct - texp[mp])


    
h_c = copyarraytoC(h)
uh_c = copyarraytoC(uh)
b_c = copyarraytoC(b)
bed_c = copyarraytoC(b)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend)
uhMend_c = copyarraytoC(uhMend)


hhbc_c = mallocPy(nbc)
whbc_c = mallocPy(nbc)
uhbc_c = mallocPy(nbc)
bhbc_c = mallocPy(nbc)

wg2i = int((xWG2  - startx) / dx ) + 1 #good one
wg3i = int((xWG3  - startx) / dx ) #G
wg4i = int((xWG4  - startx) / dx ) #G
wg5i = int((xWG5 - startx) / dx ) #
wg6i = int((xWG6  - startx) / dx )
wg7i = int((xWG7 - startx) / dx )
wg8i = int((xWG8  - startx) / dx )
wg9i = int((xWG9 - startx) / dx )
wg10i = int((xWG10 - startx) / dx )
wg11i = int((xWG11 - startx) / dx ) - 1
wg12i = int((xWG12  - startx) / dx )
wg13i = int((xWG13 - startx) / dx )
wg14i = int((xWG14 - startx) / dx )

hbwg1 = h[0]
bbwg1 = b[0]

wg2im1h = readfrommem(h_c,wg2i - 1) 
wg2ih = readfrommem(h_c,wg2i) 
wg2ip1h = readfrommem(h_c,wg2i + 1) 
hbwg2 = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])   

wg3im1h = readfrommem(h_c,wg3i - 1) 
wg3ih = readfrommem(h_c,wg3i) 
wg3ip1h = readfrommem(h_c,wg3i + 1) 
hbwg3 = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])   

wg4im1h = readfrommem(h_c,wg4i - 1) 
wg4ih = readfrommem(h_c,wg4i) 
wg4ip1h = readfrommem(h_c,wg4i + 1) 
hbwg4 = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  

wg5im1h = readfrommem(h_c,wg5i - 1) 
wg5ih = readfrommem(h_c,wg5i) 
wg5ip1h = readfrommem(h_c,wg5i + 1) 
hbwg5 = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  

wg6im1h = readfrommem(h_c,wg6i - 1) 
wg6ih = readfrommem(h_c,wg6i) 
wg6ip1h = readfrommem(h_c,wg6i + 1) 
hbwg6 = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  

wg7im1h = readfrommem(h_c,wg7i - 1) 
wg7ih = readfrommem(h_c,wg7i) 
wg7ip1h = readfrommem(h_c,wg7i + 1) 
hbwg7 = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) 

wg8im1h = readfrommem(h_c,wg8i - 1) 
wg8ih = readfrommem(h_c,wg8i) 
wg8ip1h = readfrommem(h_c,wg8i + 1) 
hbwg8 = CELLRECON(wg8im1h,wg8ih,wg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])   

wg9im1h = readfrommem(h_c,wg9i - 1) 
wg9ih = readfrommem(h_c,wg9i) 
wg9ip1h = readfrommem(h_c,wg9i + 1) 
hbwg9 = CELLRECON(wg9im1h,wg9ih,wg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])   

wg10im1h = readfrommem(h_c,wg10i - 1) 
wg10ih = readfrommem(h_c,wg10i) 
wg10ip1h = readfrommem(h_c,wg10i + 1) 
hbwg10 = CELLRECON(wg10im1h,wg10ih,wg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  

wg11im1h = readfrommem(h_c,wg11i - 1) 
wg11ih = readfrommem(h_c,wg11i) 
wg11ip1h = readfrommem(h_c,wg11i + 1) 
hbwg11 = CELLRECON(wg11im1h,wg11ih,wg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  

wg12im1h = readfrommem(h_c,wg12i - 1) 
wg12ih = readfrommem(h_c,wg12i) 
wg12ip1h = readfrommem(h_c,wg12i + 1) 
hbwg12 = CELLRECON(wg12im1h,wg12ih,wg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  

wg13im1h = readfrommem(h_c,wg13i - 1) 
wg13ih = readfrommem(h_c,wg13i) 
wg13ip1h = readfrommem(h_c,wg13i + 1) 
hbwg13 = CELLRECON(wg13im1h,wg13ih,wg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i]) 

wg14im1h = readfrommem(h_c,wg14i - 1) 
wg14ih = readfrommem(h_c,wg14i) 
wg14ip1h = readfrommem(h_c,wg14i + 1) 
hbwg14 = CELLRECON(wg14im1h,wg14ih,wg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i]) 


bwg2im1h = readfrommem(bed_c,wg2i - 1) 
bwg2ih = readfrommem(bed_c,wg2i) 
bwg2ip1h = readfrommem(bed_c,wg2i + 1) 
bwg2 = CELLRECON(bwg2im1h,bwg2ih,bwg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])   

bwg3im1h = readfrommem(bed_c,wg3i - 1) 
bwg3ih = readfrommem(bed_c,wg3i) 
bwg3ip1h = readfrommem(bed_c,wg3i + 1) 
bwg3 = CELLRECON(bwg3im1h,bwg3ih,bwg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])   

bwg4im1h = readfrommem(bed_c,wg4i - 1) 
bwg4ih = readfrommem(bed_c,wg4i) 
bwg4ip1h = readfrommem(bed_c,wg4i + 1) 
bwg4 = CELLRECON(bwg4im1h,bwg4ih,bwg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  

bwg5im1h = readfrommem(bed_c,wg5i - 1) 
bwg5ih = readfrommem(bed_c,wg5i) 
bwg5ip1h = readfrommem(bed_c,wg5i + 1) 
bwg5 = CELLRECON(bwg5im1h,bwg5ih,bwg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  

bwg6im1h = readfrommem(bed_c,wg6i - 1) 
bwg6ih = readfrommem(bed_c,wg6i) 
bwg6ip1h = readfrommem(bed_c,wg6i + 1) 
bwg6 = CELLRECON(bwg6im1h,bwg6ih,bwg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  

bwg7im1h = readfrommem(bed_c,wg7i - 1) 
bwg7ih = readfrommem(bed_c,wg7i) 
bwg7ip1h = readfrommem(bed_c,wg7i + 1) 
bwg7 = CELLRECON(bwg7im1h,bwg7ih,bwg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) 

bwg8im1h = readfrommem(bed_c,wg8i - 1) 
bwg8ih = readfrommem(bed_c,wg8i) 
bwg8ip1h = readfrommem(bed_c,wg8i + 1) 
bwg8 = CELLRECON(bwg8im1h,bwg8ih,bwg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])   

bwg9im1h = readfrommem(bed_c,wg9i - 1) 
bwg9ih = readfrommem(bed_c,wg9i) 
bwg9ip1h = readfrommem(bed_c,wg9i + 1) 
bwg9 = CELLRECON(bwg9im1h,bwg9ih,bwg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])   

bwg10im1h = readfrommem(bed_c,wg10i - 1) 
bwg10ih = readfrommem(bed_c,wg10i) 
bwg10ip1h = readfrommem(bed_c,wg10i + 1) 
bwg10 = CELLRECON(bwg10im1h,bwg10ih,bwg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  

bwg11im1h = readfrommem(bed_c,wg11i - 1) 
bwg11ih = readfrommem(bed_c,wg11i) 
bwg11ip1h = readfrommem(bed_c,wg11i + 1) 
bwg11 = CELLRECON(bwg11im1h,bwg11ih,bwg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  

bwg12im1h = readfrommem(bed_c,wg12i - 1) 
bwg12ih = readfrommem(bed_c,wg12i) 
bwg12ip1h = readfrommem(bed_c,wg12i + 1) 
bwg12 = CELLRECON(bwg12im1h,bwg12ih,bwg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  

bwg13im1h = readfrommem(bed_c,wg13i - 1) 
bwg13ih = readfrommem(bed_c,wg13i) 
bwg13ip1h = readfrommem(bed_c,wg13i + 1) 
bwg13 = CELLRECON(bwg13im1h,bwg13ih,bwg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i]) 

bwg14im1h = readfrommem(bed_c,wg14i - 1) 
bwg14ih = readfrommem(bed_c,wg14i) 
bwg14ip1h = readfrommem(bed_c,wg14i + 1) 
bwg14 = CELLRECON(bwg14im1h,bwg14ih,bwg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i]) 

nwg1s.append(hbwg1)
nwg2s.append(hbwg2 + bwg2 + hb)
nwg3s.append(hbwg3 + bwg3 + hb)
nwg4s.append(hbwg4 + bwg4 + hb)
nwg5s.append(hbwg5 + bwg5 + hb)
nwg6s.append(hbwg6 + bwg6 + hb)
nwg7s.append(hbwg7 + bwg7 + hb)
nwg8s.append(hbwg8 + bwg8 + hb)
nwg9s.append(hbwg9 + bwg9 + hb)
nwg10s.append(hbwg10 + bwg10 + hb)
nwg11s.append(hbwg11 + bwg11 + hb)
nwg12s.append(hbwg12 + bwg12 + hb)
nwg13s.append(hbwg13 + bwg13 + hb)
nwg14s.append(hbwg14 + bwg14 + hb)


h0s = [hb]
h1s = [hb]
u0s = [0]
u1s = [0]
#Just an FEM solve here
for i in range(1,len(t)):
    #evolveBCChange(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    evolvewrapIncomwavetoDir(uh_c, h_c,b_c,ftc0, ftc1 , hMend_c, wMend_c , uhMend_c,bMend_c,g,dx,dt,n,nBC,nbc, theta, hhbc_c,whbc_c,bhbc_c,uhbc_c)
    wg2im1h = readfrommem(h_c,wg2i - 1) 
    wg2ih = readfrommem(h_c,wg2i) 
    wg2ip1h = readfrommem(h_c,wg2i + 1) 
    nwg2s.append(CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])  + bwg2 + hb )
    
    wg3im1h = readfrommem(h_c,wg3i - 1) 
    wg3ih = readfrommem(h_c,wg3i) 
    wg3ip1h = readfrommem(h_c,wg3i + 1) 
    nwg3s.append(CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])  + bwg3 + hb)   
    
    wg4im1h = readfrommem(h_c,wg4i - 1) 
    wg4ih = readfrommem(h_c,wg4i) 
    wg4ip1h = readfrommem(h_c,wg4i + 1) 
    nwg4s.append(CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  + bwg4 + hb)  
    
    wg5im1h = readfrommem(h_c,wg5i - 1) 
    wg5ih = readfrommem(h_c,wg5i) 
    wg5ip1h = readfrommem(h_c,wg5i + 1) 
    nwg5s.append(CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  + bwg5 + hb)  
    
    wg6im1h = readfrommem(h_c,wg6i - 1) 
    wg6ih = readfrommem(h_c,wg6i) 
    wg6ip1h = readfrommem(h_c,wg6i + 1) 
    nwg6s.append(CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  + bwg6 + hb)  
    
    wg7im1h = readfrommem(h_c,wg7i - 1) 
    wg7ih = readfrommem(h_c,wg7i) 
    wg7ip1h = readfrommem(h_c,wg7i + 1) 
    nwg7s.append(CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i])  + bwg7 + hb )
    
    wg8im1h = readfrommem(h_c,wg8i - 1) 
    wg8ih = readfrommem(h_c,wg8i) 
    wg8ip1h = readfrommem(h_c,wg8i + 1) 
    nwg8s.append(CELLRECON(wg8im1h,wg8ih,wg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])  + bwg8 + hb)   
    
    wg9im1h = readfrommem(h_c,wg9i - 1) 
    wg9ih = readfrommem(h_c,wg9i) 
    wg9ip1h = readfrommem(h_c,wg9i + 1) 
    nwg9s.append(CELLRECON(wg9im1h,wg9ih,wg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])  + bwg9 + hb)   
    
    wg10im1h = readfrommem(h_c,wg10i - 1) 
    wg10ih = readfrommem(h_c,wg10i) 
    wg10ip1h = readfrommem(h_c,wg10i + 1) 
    nwg10s.append(CELLRECON(wg10im1h,wg10ih,wg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  + bwg10 + hb ) 
    
    wg11im1h = readfrommem(h_c,wg11i - 1) 
    wg11ih = readfrommem(h_c,wg11i) 
    wg11ip1h = readfrommem(h_c,wg11i + 1) 
    nwg11s.append(CELLRECON(wg11im1h,wg11ih,wg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  + bwg11 + hb)  
    
    wg12im1h = readfrommem(h_c,wg12i - 1) 
    wg12ih = readfrommem(h_c,wg12i) 
    wg12ip1h = readfrommem(h_c,wg12i + 1) 
    nwg12s.append(CELLRECON(wg12im1h,wg12ih,wg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  + bwg12 + hb) 
    
    wg13im1h = readfrommem(h_c,wg13i - 1) 
    wg13ih = readfrommem(h_c,wg13i) 
    wg13ip1h = readfrommem(h_c,wg13i + 1) 
    nwg13s.append(CELLRECON(wg13im1h,wg13ih,wg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i])  + bwg13 + hb) 
    
    wg14im1h = readfrommem(h_c,wg14i - 1) 
    wg14ih = readfrommem(h_c,wg14i) 
    wg14ip1h = readfrommem(h_c,wg14i + 1) 
    nwg14s.append(CELLRECON(wg14im1h,wg14ih,wg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i])  + bwg14 + hb) 

    
    nwg1s.append(ftc1)    
    
    ct = t[i] +dt
    mp = int(ct/sr)
    ftc0 = ftc1
    ftc1 = lineinterp(WG1exp[mp],WG1exp[mp + 1],texp[mp],texp[mp + 1],ct - texp[mp])
    
    h0t = readfrommem(h_c,0) 
    h1t = readfrommem(h_c,1) 
    uh0t = readfrommem(uh_c,0) 
    uh1t = readfrommem(uh_c,1)
    
    h0s.append(h0t)
    h1s.append(h1t)
    u0s.append(uh0t/h0t)
    u1s.append(uh1t/h1t)
    
 
    print(t[i])
    


hC = copyarrayfromC(h_c,n)
uhC = copyarrayfromC(uh_c,n)  
uC = array(uhC) / array(hC)


s = wdir +  "h0u0.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['x','t' ,'h0','u0'])        
                   
     for j in range(len(h0s)):
         writefile2.writerow([str(x[0]), str(t[j]), str(h0s[j]),u0s[j]])


s = wdir +  "h1u1.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['x','t' ,'h0','u0'])        
                   
     for j in range(len(h0s)):
         writefile2.writerow([str(x[1]), str(t[j]), str(h1s[j]),u1s[j]])
         
s = wdir +  "WGE.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['t','1','2','3','4','5','6','7','8','9','10','11','12','13','14' ])        
                   
     for j in range(len(texp)):
         writefile2.writerow([str(texp[j]), str(WG1exp[j]), str(WG2exp[j]), str(WG3exp[j]), str(WG4exp[j]), str(WG5exp[j]), str(WG6exp[j]), str(WG7exp[j]), str(WG8exp[j]), str(WG9exp[j]), str(WG10exp[j]), str(WG11exp[j]), str(WG12exp[j]), str(WG13exp[j]), str(WG14exp[j])])


s = wdir +  "WGN.txt"
with open(s,'w') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['t','1','2','3','4','5','6','7','8','9','10','11','12','13','14' ])        
                   
     for j in range(len(t)):
         writefile2.writerow([str(t[j]), str(nwg1s[j]), str(nwg2s[j]), str(nwg3s[j]), str(nwg4s[j]), str(nwg5s[j]), str(nwg6s[j]), str(nwg7s[j]), str(nwg8s[j]), str(nwg9s[j]), str(nwg10s[j]), str(nwg11s[j]), str(nwg12s[j]), str(nwg13s[j]), str(nwg14s[j])])



"""
wdir = "../../../../data/raw/solslopelargerSWWE60s/o2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
#Soliton Slope
a0 = 1.0
a1 = 0.01

g = 9.81
dx = 0.01
Cr = 0.5
l = 0.1
dt = l*dx
theta = 2.0
startx = -100
endx = 250.0 + dx
startt = 0.0
endt = 60 + dt  


        
szoomx = startx
ezoomx = endx
        

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  

theta = 1.2

gap = int(1.0/dt)

nBC = 3
nbc = 3*n + 2*(nBC)

idx = 1.0 / dx
        
h,u,b = soloverslope(x,a0,a1,-300,100,100,149.5,g)    

hMbeg = h[0]*ones(nBC)
uMbeg = u[0]*ones(nBC)
hMend = h[-1]*ones(nBC)
uMend = u[-1]*ones(nBC) 
bMbeg = b[0]*ones(nBC)
bMend = b[-1]*ones(nBC)

wMbeg = hMbeg + bMbeg
wMend = hMend + bMend
    
h_c = copyarraytoC(h)
u_c = copyarraytoC(u)
b_c = copyarraytoC(b)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)


hhbc_c = mallocPy(nbc)
whbc_c = mallocPy(nbc)
uhbc_c = mallocPy(nbc)
bhbc_c = mallocPy(nbc)


#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrap(u_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,dt,n,nBC,nbc,theta,hhbc_c,whbc_c,bhbc_c,uhbc_c)
    print(t[i])




hC = copyarrayfromC(h_c,n)
uC = copyarrayfromC(u_c,n)     
bC = copyarrayfromC(b_c,n)      

hhbcC = copyarrayfromC(hhbc_c,nbc)
whbcC = copyarrayfromC(whbc_c,nbc)
uhbcC = copyarrayfromC(uhbc_c,nbc)
bhbcC = copyarrayfromC(bhbc_c,nbc)


s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','Eval' , 'initial Eval',"cell midpoint" ,'h', 'u','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]), str(hC[j]) , str(uC[j]) ,str(b[j])])

deallocPy(h_c)
deallocPy(b_c)
deallocPy(u_c)

deallocPy(hhbc_c)
deallocPy(whbc_c)
deallocPy(uhbc_c)
deallocPy(bhbc_c)

deallocPy(hMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(uMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
deallocPy(bMbeg_c)
deallocPy(bMend_c)
"""