# -*- coding: utf-8 -*-
"""
Created on Mon Apr 07 17:31:45 2014

@author: Jordan
"""


from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os

def minmod(a,b,c):
    
    if (a>0 and b>0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0
        
def solveucon(con):
    
    n = len(con)
    u = []
    for i in range(n):
        if (con[i][0] == 0.0):
            u.append(0.0)
        else:
            u.append(con[i][1] / con[i][0])
    
    return u
            
            
    
    
def evolve(con,bed,dt,dx,g,beta,u0,u1,h0,h1,b0,b1):
    #get averages
    idx = 1.0 / dx  
    n = len(con)
    
    ncon = zeros((n,2))
    #add in ghost cells
    
    #beginning   
    conlb = zeros((3,2))
    conrb = zeros((3,2))
    bbeg = b0[1:4]
    bend = b1[0:3]
    
    
    conlb[0][0] = h0[-1]
    conlb[1][0] = h0[-2]
    conlb[2][0] = h0[-3]
    conlb[0][1] = u0[-1]*h0[-1]
    conlb[1][1] = u0[-2]*h0[-2]
    conlb[2][1] = u0[-3]*h0[-3]
    
    conrb[0][0] = h1[0]
    conrb[1][0] = h1[1]
    conrb[2][0] = h1[2]
    conrb[0][1] = u1[0]*h1[0]
    conrb[1][1] = u1[1]*h1[1]
    conrb[2][1] = u1[2]*h1[2] 
    
    con = concatenate([conlb, con, conrb])
    bed = concatenate([bbeg,bed,bend])
    
    #do normal stuff
    i = 2
    #define the stage
    wi = con[i][0] + bed[i]
    wip1 = con[i+1][0] + bed[i+1]
    wip2 = con[i+2][0] + bed[i+2]
    wip3 = con[i+3][0] + bed[i+3]
    wim1 = con[i-1][0] + bed[i-1]
    wim2 = con[i-2][0] + bed[i-2]
        
    #reconstruct i left and right values
    #gradients
    dwib = wi - wim1
    dwif = wip1 - wi
    dwim = 0.5*(wip1 - wim1)
    dhib = con[i][0] - con[i-1][0]
    dhif = con[i+1][0] - con[i][0]
    dhim = 0.5*(con[i+1][0] - con[i-1][0])
    duhib = con[i][1] - con[i-1][1]
    duhif = con[i+1][1] - con[i][1]
    duhim = 0.5*(con[i+1][1] - con[i-1][1])

    dwi = minmod(beta*dwib, beta*dwif, dwim)
    dhi = minmod(beta*dhib, beta*dhif, dhim)
    duhi = minmod(beta*duhib, beta*duhif, duhim)
        
    #i left
    uhil = con[i][1] - 0.5*duhi
    hil = con[i][0] - 0.5*dhi
    wil = wi - 0.5*dwi
    bil = wil - hil
        
    #i right
    uhir = con[i][1] + 0.5*duhi
    hir = con[i][0] + 0.5*dhi
    wir = wi + 0.5*dwi
    bir = wir - hir
        
    #reconstruct i+1 left and right values
    dwip1b = wip1 - wi
    dwip1f = wip2 - wip1
    dwip1m = 0.5*(wip2 - wi)
    dhip1b = con[i+1][0] - con[i][0]
    dhip1f = con[i+2][0] - con[i+1][0]
    dhip1m = 0.5*(con[i+2][0] - con[i][0])
    duhip1b = con[i+1][1] - con[i][1]
    duhip1f = con[i+2][1] - con[i+1][1]
    duhip1m = 0.5*(con[i+2][1] - con[i][1])

    dwip1 = minmod(beta*dwip1b, beta*dwip1f, dwip1m)
    dhip1 = minmod(beta*dhip1b, beta*dhip1f, dhip1m)
    duhip1 = minmod(beta*duhip1b, beta*duhip1f, duhip1m)

    #i+1 left
    uhip1l = con[i+1][1] - 0.5*duhip1
    hip1l = con[i+1][0] - 0.5*dhip1
    wip1l = wip1 - 0.5*dwip1
    bip1l = wip1l - hip1l

    #i+1 right
    uhip1r = con[i+1][1] + 0.5*duhip1
    hip1r = con[i+1][0] + 0.5*dhip1
    wip1r = wip1 + 0.5*dwip1 
    bip1r = wip1r - hip1r
        
    #reconstruct i+2 left values
    dwip2b = wip2 - wip1
    dwip2f = wip3 - wip2
    dwip2m = 0.5*(wip3 - wip1)
    dhip2b = con[i+2][0] - con[i+1][0]
    dhip2f = con[i+3][0] - con[i+2][0]
    dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
    duhip2b = con[i+2][1] - con[i+1][1]
    duhip2f = con[i+3][1] - con[i+2][1]
    duhip2m = 0.5*(con[i+3][1] - con[i+1][1])

    dwip2 = minmod(beta*dwip2b, beta*dwip2f, dwip2m)
    dhip2 = minmod(beta*dhip2b, beta*dhip2f, dhip2m)
    duhip2 = minmod(beta*duhip2b, beta*duhip2f, duhip2m)

    #i+2 left
    uhip2l = con[i+2][1] - 0.5*duhip2
    hip2l = con[i+2][0] - 0.5*dhip2
    wip2l = wip2 - 0.5*dwip2
    bip2l = wip2l - hip2l
        
    #reconstruct i-1 right values
    dwim1b = wim1 - wim2
    dwim1f = wi - wim1
    dwim1m = 0.5*(wi - wim2)
    dhim1b = con[i-1][0] - con[i-2][0]
    dhim1f = con[i][0] - con[i-1][0]
    dhim1m = 0.5*(con[i][0] - con[i-2][0])
    duhim1b = con[i-1][1] - con[i-2][1]
    duhim1f = con[i][1] - con[i-1][1]
    duhim1m = 0.5*(con[i][1] - con[i-2][1])
        
    dwim1 = minmod(beta*dwim1b, beta*dwim1f, dwim1m)
    dhim1 = minmod(beta*dhim1b, beta*dhim1f, dhim1m)
    duhim1 = minmod(beta*duhim1b, beta*duhim1f, duhim1m)
                
    #i-1 right
    uhim1r = con[i-1][1] + 0.5*duhim1
    him1r = con[i-1][0] + 0.5*dhim1
    wim1r = wim1 + 0.5*dwim1
    bim1r = wim1r - him1r
    
    nb = max(bip1l,bir)
    hihm = max(0,wir-nb)
    hihp = max(0,wip1l-nb)
    
    her = hihp
    uher = uhip1l
        
    hel = hihm
    uhel = uhir
    
    if(hip1l == 0):
        uer = 0.0
    else:
        uer = uher / hip1l
    if(hir == 0):
        uel = 0.0
    else:
        uel = uhel / hir
    
        
    #I think these need to be altered, its so close to a reasonable solution that it must be a small reason like this:        
    sqrtghel = sqrt(g*hel)
    sqrtgher = sqrt(g*her)
    sl = min(0,uel - sqrtghel, uer - sqrtgher)
    sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
    felh = uel*hel
    felG = uel*uel*hel + 0.5*g*hel*hel
    ferh = uer*her
    ferG = uer*uer*her + 0.5*g*her*her
               
        
    if (sl >= 0):
        foh = felh
        foG = felG
    elif(sr <= 0):
        foh = ferh
        foG = ferG
    else:
        srmsl = 1.0 / (sr - sl)
        foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
        foG = srmsl*(sr*felG - sl*ferG + sl*sr*(uher - uhel ))
    
    fih = foh
    fiG = foG
    himhp = hihp
    
    #update new values
    him1r = hir
    wim1r = wir
    bim1r = bir
    uhim1r = uhir
        
    hil = hip1l
    wil = wip1l
    bil = bip1l
    uhil = uhip1l
        
    hir = hip1r
    wir = wip1r
    bir = bip1r
    uhir = uhip1r
    
    hip1l = hip2l
    wip1l = wip2l
    bip1l = bip2l
    uhip1l = uhip2l      
        
    
    dhip1 = dhip2
    dwip1 = dwip2
    duhip1 = duhip2
    
    for i in range(3,len(con)-3):
        
        #define the stage
        wip1 = con[i+1][0] + bed[i+1]
        wip2 = con[i+2][0] + bed[i+2]
        wip3 = con[i+3][0] + bed[i+3]
        

        #i+1 right
        uhip1r = con[i+1][1] + 0.5*duhip1
        hip1r = con[i+1][0] + 0.5*dhip1
        wip1r = wip1 + 0.5*dwip1 
        bip1r = wip1r - hip1r
        
             
        #reconstruct i+2 left values
        dwip2b = wip2 - wip1
        dwip2f = wip3 - wip2
        dwip2m = 0.5*(wip3 - wip1)
        dhip2b = con[i+2][0] - con[i+1][0]
        dhip2f = con[i+3][0] - con[i+2][0]
        dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
        duhip2b = con[i+2][1] - con[i+1][1]
        duhip2f = con[i+3][1] - con[i+2][1]
        duhip2m = 0.5*(con[i+3][1] - con[i+1][1])

        dwip2 = minmod(beta*dwip2b, beta*dwip2f, dwip2m)
        dhip2 = minmod(beta*dhip2b, beta*dhip2f, dhip2m)
        duhip2 = minmod(beta*duhip2b, beta*duhip2f, duhip2m)

        #i+2 left
        uhip2l = con[i+2][1] - 0.5*duhip2
        hip2l = con[i+2][0] - 0.5*dhip2
        wip2l = wip2 - 0.5*dwip2
        bip2l = wip2l - hip2l
    
        nb = max(bip1l,bir)
        hihm = max(0,wir-nb)
        hihp = max(0,wip1l-nb)

        her = hihp
        uher = uhip1l
        #ber = bip1l
        
        hel = hihm
        uhel = uhir
        #bel = bir
        
        if(hip1l == 0):
            uer = 0.0
        else:
            uer = uher / hip1l
        if(hir == 0):
            uel = 0.0
        else:
            uel = uhel / hir
        
        #calculate the source term
        th = con[i][0]
        tbx = (bil - bir)
        
        sourcer = g*0.5*(hihm*hihm - hir*hir)
        sourcec = g*th*tbx     
        sourcel = g*0.5*(hil*hil - himhp*himhp)
                             
        
        sqrtghel = sqrt(g*hel)
        sqrtgher = sqrt(g*her)
        sl = min(0,uel - sqrtghel, uer - sqrtgher)
        sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
        felh = uel*hel
        felG = uel*uel*hel + 0.5*g*hel*hel
        ferh = uer*her
        ferG = uer*uer*her + 0.5*g*her*her
  
        if (sl >= 0):
           foh = felh
           foG = felG
        elif(sr <= 0):
           foh = ferh
           foG = ferG
        else:
           srmsl = 1.0 / (sr - sl)
           foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
           foG = srmsl*(sr*felG - sl*ferG + sl*sr*(uher - uhel ))
        
        
        
        ncon[i-3][0] = con[i][0] - dt*idx*(foh - fih)
        ncon[i-3][1] = con[i][1] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + sourcec)
        
                
        fih = foh
        fiG = foG
        himhp = hihp
        
        #update new values
        him1r = hir
        wim1r = wir
        bim1r = bir
        uhim1r = uhir
        
        hil = hip1l
        wil = wip1l
        bil = bip1l
        uhil = uhip1l
        
        hir = hip1r
        wir = wip1r
        bir = bip1r
        uhir = uhip1r
        
        hip1l = hip2l
        wip1l = wip2l
        bip1l = bip2l
        Gip1l = Gip2l
        uhip1l = uhip2l       
        
        
        dhip1 = dhip2
        dwip1 = dwip2
        duhip1 = duhip2
 
    return ncon
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t  
    
def dambreak(n,divn,normheight,damheight,bot,u0,u1,h0,h1,b0,b1,dx):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    dampos = n / divn
    
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(0,dampos-1):
        con[i][0]= damheight
        bx[i] = bot
    
    for i in range(dampos-1,n):
        con[i][0]= normheight
        bx[i] = bot
        
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
            
    return con,bx  
    
def flatjump(n,stage,addheight,addarea,bot,dx):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    jumppos = n / 2
    
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(n):
        bx[i] = bot
        con[i][0] = stage - bx[i]
        
        if( i  > jumppos - addarea and i  < jumppos + addarea):
            con[i][0]= con[i][0] + addheight
            
    
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
    return con,bx  

def linearbump(n,stage,vel,normbot,slope,addareal,addarear,dx):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    bx = zeros(n)
    h = zeros(n)
    u = zeros(n)
    
   
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(n):
        u[i] = vel
        bx[i] = normbot
        
        if( i  >= addareal and i  <=  addarear):
            bx[i] = normbot + slope*(i - addareal)
        if(i > addarear):
            bx[i] = normbot + slope*(addarear- addareal)
        
        con[i][0] = stage - bx[i]
            
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
            
    return con,bx
    
def parabolicbump(n,stage,vel,normbot,slope,top,addareal,addarear,dx):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    bx = zeros(n)
    h = zeros(n)
    u = zeros(n)
    
   
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(n):
        u[i] = vel
        bx[i] = normbot
        
        midpoint = (addareal + addarear)/2 
        
        if( i  >= addareal and i  <=  addarear):
            bx[i] = top + slope*(i - midpoint)*(i - midpoint)
        
        con[i][0]= stage - bx[i]
            
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
            
    return con,bx
    
def linearbumpjump(n,stage,vel,normbot,slope,baddareal,baddarear,haddareal,haddarear,hadd,dx):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
   
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(n):
        u[i] = vel
        bx[i] = normbot
        
        if( i  >= baddareal and i  <=  baddarear):
            bx[i] = normbot + slope*(i - baddareal)
        if(i > baddarear):
            bx[i] = normbot + slope*(baddarear- baddareal)
        
        con[i][0]= stage - bx[i]
        
        if( i  >= haddareal and i  <=  haddarear):
            con[i][0] = con[i][0] + hadd
            
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
            
    return con,bx
    
def flatlake(n,stage,normbot,vel,func,c1,c2,haddl,haddr,hadd):
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        u[i] = vel
        bx[i] =  normbot + c1*func(c2*(i - n/2))
        con[i][0] = stage - bx[i]
        
        if(i >= haddl and i <= haddr):
            con[i][0] = con[i][0]+hadd
        
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
    return con,bx
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a


def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3*a1) / (2*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        con[i][0] = he
        u[i] = c* ((he - a0) / he)
    
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
    return con,bx
    
def dingemans(x,dx):
    n = len(x)
    nbot = -0.4
    hbot = -0.1
    stage = 0.0
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    #very start
    for i in range(n):
        bx[i] = nbot        
        if(x[i] >= 6.0 and x[i] <= 12.0):
            bx[i] = nbot + ((hbot - nbot)/ 6.0)*(x[i] - 6.0)
        if(x[i] > 12.0 and x[i] < 14.0):
            bx[i] = hbot
        if(x[i] >= 14.0 and x[i] <= 17.0):
            bx[i] = hbot + ((nbot - hbot)/ 3.0)*(x[i] - 14.0)
        
        con[i][0] = stage - bx[i]
        
    
    for i in range(n):
        con[i][1] = u[i]*con[i][0]
       

    return con,bx

def eta(x,t,T,a,l):
    il = 1.0 /l
    t1 = (x*il - t/T)
    ft = a*cos(2*pi*t1)
    st = pi*il*a*a*cos(4*pi*t1)
    tt = - pi*pi*il*il*a*a*a*0.5*(cos(2*pi*t1) - cos(6*pi*t1))
    
    return ft + st + tt
    
    

#set it up so its exact floating point
dx = 0.05
l = 0.02
dt = l*dx
startx = 0.0
endx = 30.0
startt = 0.0
endt = 20*dt

szoomx = startx
ezoomx = endx

wdir = "../../../data/raw/SWWtest/test/"
if not os.path.exists(wdir):
    os.makedirs(wdir)


g = 9.8
ithree = 1.0 / 3
idx = 1/ dx

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

stage = 10.0
u0 = 0.0
gap = 20
gapbig = gap * 25
"""
g1u = []
g2u = []
g3u = []
g4u = []
g5u = []
g6u = []
g7u = []
g8u = []
g9u = []
g10u = []
g1h = []
g2h = []
g3h = []
g4h = []
g5h = []
g6h = []
g7h = []
g8h = []
g9h = []
g10h = []
"""

#con,bed = dambreak(n,2,h1,h0,0,u0,u1,h0,h1,0.0,0.0,dx)
#con,bed = flatjump(n,stage,2,100,0,dx)
#con,bed = linearbump(n,stage,u0,1.0,-0.004,n/2 - 50,n/2 + 50,dx)
#con,bed = linearbumpjump(n,stage,u0,1.0,0.002,n/2-200,5*n/8+200,n/3 - 20,n/3 + 20,0.3,dx)
#con,bed = parabolicbump(n,stage,u0,1.0,-0.000062,2.0,15*n/32,17*n/32,dx)
con,bed = flatlake(n,stage,3,u0,sin,1.5,0.01575,n/2 - 20,n/2 + 20,0.0)

#con,bed = solitoninit(n,stage,a,g,x,t[0],0,dx)

#con, bed = dingemans(x,dx)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([0.0,0.0,0.0,0.0])
u1 = array([0.0,0.0,0.0,0.0])

h0 = array([stage,stage,stage,stage]) - b0
h1 = array([stage,stage,stage,stage]) - b1

for i in range(len(t)):
   
        
        
    if(i % gapbig == 0):
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file1:
            writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile.writerow(['dt(s)','#timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
            for j in range(n):
             
                 writefile.writerow([str(dt),str(i), str(con[j][0]) , str(con[j][1]) , str(con[j][1] / con[j][0])])  
             
        
    con1 = evolve(con ,bed,dt,dx,g,1.2,u0,u1,h0,h1,b0,b1)
    con2 = evolve(con1,bed,dt,dx,g,1.2,u0,u1,h0,h1,b0,b1)
    con = 0.5*(con + con2)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dt(s)','#timestep', 'height(m)', 'G' , 'u(m/s)' ])        
       
    for j in range(n):
     
         writefile.writerow([str(dt),str(i), str(con[j][0]) , str(con[j][1]) , str(con[j][1] / con[j][0])]) 

"""    
s = wdir + "gauges.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['1h','1u','2h','2u','3h','3u','4h','4u','5h','5u','6h','6u','7h','7u','8h','8u','9h','9u','10h','10u',])        
               
    for j in range(len(t)):
             
        writefile.writerow([g1h[j],g1u[j],g2h[j],g2u[j],g3h[j],g3u[j],g4h[j],g4u[j],g5h[j],g5u[j],g6h[j],g6u[j],g7h[j],g7u[j],g8h[j],g8u[j],g9h[j],g9u[j],g10h[j],g10u[j]])    
"""

    
    
##READ FILE
"""
s = wdir + "saveoutputts" + str(0) + ".txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    acon = zeros((n,2))
    j = -1
    for row in readfile:       
        if (j >= 0):
            dt = float(row[0])
            i = int(row[1])
            acon[j][0] = float(row[2])
            acon[j][1] = float(row[3])
        
            
        j = j + 1
"""


