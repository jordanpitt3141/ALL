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
    
def splitcon(con):
    
    n = len(con)
    u = []
    h = []
    for i in range(n):
        h.append(con[i][0])
        if (con[i][0] == 0.0):
            u.append(0.0)
        else:
            u.append(con[i][1] / con[i][0])
    
    return u,h
            
def evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt): 
    n = len(con)
    ncon = zeros((n,2))
    #update h' and G'    
    conp = evolve(con,bed,dt,dx,g,beta,u0,u1,h0,h1,b0,b1)
        
    #update h'' and G''
    conpp = evolve(conp,bed,dt,dx,g,beta,u0,u1,h0,h1,b0,b1)
    
    ncon = 0.5 * (con + conpp)
    
    return ncon            
    
    
def evolve(con,bed,dt,dx,g,beta,u0,u1,h0,h1,b0,b1):
    #get averages
    idx = 1.0 / dx  
    n = len(con)
    
    ncon = zeros((n,2))
    #add in ghost cells
    
    #beginning   
    conlb = zeros((2,2))
    conrb = zeros((2,2))
    bbeg = b0[2:4]
    bend = b1[0:2]
    
    
    conlb[0][0] = h0[-2]
    conlb[1][0] = h0[-3]
    conlb[0][1] = u0[-2]*h0[-2]
    conlb[1][1] = u0[-3]*h0[-3]
    
    conrb[0][0] = h1[0]
    conrb[1][0] = h1[1]
    conrb[0][1] = u1[0]*h1[0]
    conrb[1][1] = u1[1]*h1[1]
    
    con = concatenate([conlb, con, conrb])
    bed = concatenate([bbeg,bed,bend])
    
    #do normal stuff
    i = 1
    #define the stage
    wi = con[i][0] + bed[i]
    wip1 = con[i+1][0] + bed[i+1]
    wip2 = con[i+2][0] + bed[i+2]
    wim1 = con[i-1][0] + bed[i-1]
        
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
    feluh = uel*uel*hel + 0.5*g*hel*hel
    ferh = uer*her
    feruh = uer*uer*her + 0.5*g*her*her
               
        
    if (sl >= 0):
        foh = felh
        fouh = feluh
    elif(sr <= 0):
        foh = ferh
        fouh = feruh
    else:
        srmsl = 1.0 / (sr - sl)
        foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
        fouh = srmsl*(sr*feluh - sl*feruh + sl*sr*(uer*her - uel*hel ))
    
    fih = foh
    fiuh = fouh
    himhp = hihp
        
    hil = hip1l
    bil = bip1l
        
    dwi = dwip1
    dhi = dhip1
    duhi = duhip1 
    
    for i in range(2,len(con)-2):
        
        #define the stage
        wi = con[i][0] + bed[i]
        wip1 = con[i+1][0] + bed[i+1]
        wip2 = con[i+2][0] + bed[i+2]
            
        #reconstruct i left and right values            
        #i right
        uhir = con[i][1] + 0.5*duhi
        hir = con[i][0] + 0.5*dhi
        wir = wi + 0.5*dwi
        bir = wir - hir
            
        #reconstruct i+1 left values
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
        feluh = uel*uel*hel + 0.5*g*hel*hel
        ferh = uer*her
        feruh = uer*uer*her + 0.5*g*her*her
  
        if (sl >= 0):
            foh = felh
            fouh = feluh
        elif(sr <= 0):
            foh = ferh
            fouh = feruh
        else:
            srmsl = 1.0 / (sr - sl)
            foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
            fouh = srmsl*(sr*feluh - sl*feruh + sl*sr*(uer*her - uel*hel ))      
        
        
        ncon[i-2][0] = con[i][0] - dt*idx*(foh - fih)
        ncon[i-2][1] = con[i][1] - dt*idx*(fouh - fiuh) + dt*idx*(sourcer+sourcel + sourcec)
        
                
        fih = foh
        fiuh = fouh
        himhp = hihp
        
        hil = hip1l
        bil = bip1l
        
        dwi = dwip1
        dhi = dhip1
        duhi = duhip1 
 
    return ncon
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t  
    

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    con = zeros((n,2))
    bx = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        ue =  c* ((he - a0) / he)
        
        con[i][0] = he
        con[i][1] = ue*he
         
    return con,bx
    
def solitoninit1(n,stage,a0,a1,g,x,t0,bot,dx):
    con = zeros((n,2))
    bx = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        ue =  c* ((he - a0) / he)
        
        con[i][0] = stage + (he-a0)
        con[i][1] = ue*he
         
    return con,bx

def flatlake(x,dx,stage,normbot,vel,func,c1,c2,haddl,haddr,hadd):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    for i in range(n):
        bx[i] =  normbot + c1*func(c2*(x[i] - x[n/2]))
        con[i][0] = stage - bx[i]
        
        if(i >= haddl and i <= haddr):
            con[i][0] = con[i][0]+hadd
            
        con[i][1] = vel*con[i][0]
        
    return con,bx
    
def compflatlake(x,dx,stage,normbot,vel,func,c1,c2,haddl,haddr,hadd):
    n = len(x)
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        u[i] = vel
        bx[i] =  normbot + c1*func(c2*(x[i] - x[n/2]))
        h[i] = stage - bx[i]
    return h,u
    
def flatlakedisc(x,dx,stage,bot1,bot2,xlb,xrb):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    for i in range(n):
        con[i][1] = 0.0
        
        if(x[i] >= xlb and x[i] <= xrb):
            bx[i] = bot2
        else:
            bx[i] = bot1
            
        con[i][0] = stage - bx[i]
        
    return con,bx
    
def compflatlakedisc(x,dx,stage,bot1,bot2,xlb,xrb):
    n = len(x)
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        u[i] = 0.0
        
        if(x[i] >= xlb and x[i] <= xrb):
            bx[i] = bot2
        else:
            bx[i] = bot1
            
        h[i] = stage - bx[i]

    return h,u
    

def dambreak(x,xc,hf,hl,bot):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    for i in range(n):
       
        if(x[i] >= xc):
            con[i][0] = hl
        else:
            con[i][0] = hf
            
        bx[i] = bot
        
    return con ,bx
    
def flowoverbump(x,stage,bot,bc,bm,ba,u0,u1,h0,h1):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    bg = (bm - bot)/ (ba*ba)
    
    for i in range(n):
        if(x[i] >= bc - ba and x[i] <= bc + ba):
            bx[i] = bm - bg*(x[i] - bc)*(x[i] - bc)
            
        con[i][0] = stage - bx[i]
        con[i][1] = u0*con[i][0]
            
        
    return con,bx
    
def solitonupslope(x,stage,bot1,bot2,bxl,bxr,a1,g,t0,sxl,sxr):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    a0 = stage - bot1
    
    c = (bot1 - bot2) / (bxl - bxr)
    a = bot1 - c*bxl
    d = sqrt(g*(a0 + a1))
    
    for i in range(n):
       
        if(x[i] >= bxl and x[i] <= bxr):
            bx[i] = a + c*x[i]
        elif(x[i] < bxl):
            bx[i] = bot1
        else:
            bx[i] = bot2

        con[i][0] = stage - bx[i]
        
        if(x[i] >= sxl and x[i] <= sxr):
            sxm = 0.5*(sxl + sxr)
            he = soliton((x[i] - sxm),t0,g,a0,a1)
            ue =  d* ((he - a0) / he)
        
            con[i][0] = con[i][0] + (he - a0)
            con[i][1] = ue*con[i][0]
            
        
    return con,bx
    
def soloverslope(x,a0,a1,slopbeg,slopend,topbed,g):
    """
    soloverslope : initial conditions for soliton traveling up a slope experiment
    
        Input:
            x       : array of cell centres
            a0      : stage of water
            a1      : amplitude of soliton
            slopbeg : location of start of slope
            slopend : location of end of slope
            topbed  : height of bed at the end of the slope (starting height is 0)
            g       : acceleration due to gravity
            
            
        Output:
            arrays of height (h) and velocity (u) of water and bed profile (bed) at cell centres
    """
	
    n = len(x)
    con = zeros((n,2))
    bed = zeros(n)

    #speed of the soliton
    c = sqrt(g*(a0 + a1))

    for i in range(n):
		
      #This is the range over which we define a soliton with a bed of 0 beneath it, which is smaller than the constant depth area  
      if (x[i] < slopbeg):
          bed[i] = 0
          con[i][0] = soliton(x[i],0,g,a0,a1)
          con[i][1] =  c* ((con[i][0] - a0) / con[i][0])

      #This is the region in which the bed has a linear slope with a constant stage
      elif(x[i] >= slopbeg and x[i] <= slopend):
          bed[i] = (float(topbed) / (slopend - slopbeg))*(x[i] - slopbeg)
          con[i][0] = a0 - bed[i]
          con[i][1] = 0

      #After the slope the bed is constant as is the stage
      elif(x[i] > slopend):
          bed[i] = topbed
          con[i][0] = a0 - bed[i]
          con[i][1] = 0
   
    return con,bed
    
    
def jumpupslope(x,stage,bot1,bot2,bxl,bxr,ha,g,t0,sxl,sxr):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    c = (bot1 - bot2) / (bxl - bxr)
    a = bot1 - c*bxl
    
    for i in range(n):
       
        if(x[i] >= bxl and x[i] <= bxr):
            bx[i] = a + c*x[i]
        elif(x[i] < bxl):
            bx[i] = bot1
        else:
            bx[i] = bot2

        con[i][0] = stage - bx[i]
        
        if(x[i] >= sxl and x[i] <= sxr):       
            con[i][0] = con[i][0] + ha
            
        
    return con,bx
    
def flatbump(x,stage,bot,hadd,uadd,xr,xl):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    
    for i in range(n):
        bx[i] = bot            
        con[i][0] = stage - bx[i]
        
        if(x[i] >= xr and x[i] <= xl):
            con[i][0] = con[i][0] +  hadd
            con[i][1] = uadd*con[i][0]
        
    return con, bx

###### DAM BREAK #############
"""
#set it up so its exact floating point
dx = 0.5
l = 0.01
dt = l*dx
startx = 0.0
endx = 1000.0
startt = 0.0
endt = 30.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../swdata/dambreak/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
ui = 0.0
bot = 0.0
hf = 1.8
hl = 1.0
gap = 20
gapbig = gap * 5



con,bed = dambreak(x,500,hf,hl,bot)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([ui,ui,ui,ui])
u1 = array([ui,ui,ui,ui])
    
h0 = array([hf,hf,hf,hf])
h1 = array([hl,hl,hl,hl])

for i in range(1,len(t)):      
        
    if(i % gapbig == 0 or i ==1):
        u = solveucon(con)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])         
               
            for j in range(n):
             
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
u = solveucon(con)
with open(s,'a') as file2:
    
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dt(s)','timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
     for j in range(n):
             
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
 """
  




######## FLOW OVER BUMP ##############
"""
#set it up so its exact floating point
dx = 0.1
l = 0.02
dt = l*dx
startx = 0.0
endx = 600.0
startt = 0.0
endt = 60.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../swdata/bump/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
uf = 1.5
ul = 1.5
hf = 3.0
hl = 3.0
bot = 0.0
stage = 3.0
gap = 20
gapbig = gap * 5



con, bed = flowoverbump(x,stage,bot,300,1.5,50,uf,ul,hf,hl)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([uf,uf,uf,uf])
u1 = array([ul,ul,ul,ul])
    
h0 = array([hf,hf,hf,hf])
h1 = array([hl,hl,hl,hl])

for i in range(1,len(t)):      
        
    if(i % gapbig == 0 or i ==1):
        u = solveucon(con)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])         
               
            for j in range(n):
             
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
u = solveucon(con)
with open(s,'a') as file2:
    
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dt(s)','timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
     for j in range(n):
             
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""


### SMOOTH FLAT LAKE ####################################
"""
from numpy.linalg import norm

wdir = "../swdata/smfl/"
beta = 2.0

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','L1-norm Difference Height', 'L1-norm Difference Velocity'])

nxs = [10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
#do accuracy test
for k in range(len(nxs)):
    dx = nxs[k]
    l = 0.02
    dt = l*dx
    startx = 0.0
    endx = 1000.0
    startt = 0
    endt = 30.0 + dt
    
    
    szoomx = startx
    ezoomx = endx
        
    g = 10.0
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    ui = 0.0
    c1 = 2.0
    c2 = 0.02
    normbot = 2.0
    stage = 10.0
    gap = 20
    gapbig = gap * 25
        
    x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
    x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

    
    
    con,bed = flatlake(x,dx,stage,normbot,0.0,sin,c1,c2,0,0,0.0)
    
    #change bed and h0
    u0 = array([ui,ui,ui,ui])
    u1 = array([ui,ui,ui,ui])
    
    b0 = normbot + c1*sin(c2*(x0 - x[n/2]))
    b1 = normbot + c1*sin(c2*(x1 - x[n/2]))
    
    h0 = array([stage,stage,stage,stage]) - b0
    h1 = array([stage,stage,stage,stage]) - b1
    
    

    for i in range(1,len(t)): 
                       
        con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
        print t[i]
        print con[5]

                
    h = []
    for j in range(n):
        h = append(h,con[j][0])
    
    u = solveucon(con)
            
    htrue,utrue = compflatlake(x,dx,stage,normbot,0.0,sin,c1,c2,0,0,0.0)

    normhdiffi = norm(h - htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1)
    
    
    s = wdir + "save"+ str(k)+".txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow(['time','dx','Height Approximate', 'Velocity Approximate', 'Bed', 'Height Exact', 'Velocity Exact'])        
               
        for j in range(n):
             
            writefile.writerow([str(t[-1]),str(dx),str(h[j]), str(u[j]), str(bed[j]), str(htrue[j]), str(utrue[j])])
            
    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)    
                           
        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)])
"""

##############Discontinuous Flate Lake#################################
"""
from numpy.linalg import norm

wdir = "../swdata/dfl/"
beta = 2.0

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','L1-norm Difference Height', 'L1-norm Difference Velocity'])

nxs = [10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
#do accuracy test
for k in range(16,len(nxs)):
    dx = nxs[k]
    l = 0.02
    dt = l*dx
    startx = 0.0
    endx = 1000.0
    startt = 0
    endt = 30.0 + dt
    
    
    szoomx = startx
    ezoomx = endx
        
    g = 10.0
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    ui = 0.0
    xlb = 400
    xrb = 600
    bot1 = 3.0
    bot2 = 5.0
    stage = 10.0
    gap = 20
    gapbig = gap * 25
        
    x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
    x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

    
    
    con,bed = flatlakedisc(x,dx,stage,bot1,bot2,xlb,xrb)
    
    #change bed and h0
    u0 = array([ui,ui,ui,ui])
    u1 = array([ui,ui,ui,ui])
    
    b0 = array([bot1,bot1,bot1,bot1])
    b1 = array([bot1,bot1,bot1,bot1])
    
    h0 = array([stage,stage,stage,stage]) - b0
    h1 = array([stage,stage,stage,stage]) - b1
    
    

    for i in range(1,len(t)):
                       
        con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
        print t[i]
        print con[5]

                
    h = []
    for j in range(n):
        h = append(h,con[j][0])
    
    u = solveucon(con)
            
    htrue,utrue = compflatlakedisc(x,dx,stage,bot1,bot2,xlb,xrb)

    normhdiffi = norm(h - htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1)
    
    
    s = wdir + "save"+ str(k)+".txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow(['time','dx','Height Approximate', 'Velocity Approximate', 'Bed', 'Height Exact', 'Velocity Exact'])        
               
        for j in range(n):
             
            writefile.writerow([str(t[-1]),str(dx),str(h[j]), str(u[j]), str(bed[j]), str(htrue[j]), str(utrue[j])])
            
    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)    
                           
        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)])           
"""


#########PLOT TO CHECK################
"""
        if(i % gap == 0 or i == 1):
            ni = i / gap
            h = []
            G = []
            for j in range(n):
                h = append(h,con[j][0])
                G = append(G,con[j][1])
        
            u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        
            plot(x,h+bed,'b', label="1")
            plot(x,bed,'g', label="2")
            xlim([startx,endx])
            ylim([-0.1,11.0])
            title("Flow over Bump")
            xlabel("Distance (m)")
            ylabel("Water Height (m)")
    
            s = wdir + "height" + str(ni) + ".png"
            
            savefig(s, bbox_inches='tight')        
            clf()

            
            plot(x,u,'r', label="1")
            xlim([startx,endx])
            ylim([-1.0,1.0])
            title("Flow over Bump")
            xlabel("Distance (m)")
            ylabel("Velocity (m/s)")
    
            s = wdir + "velocity" + str(ni) + ".png"
    
            savefig(s, bbox_inches='tight')
            clf()
"""

##### Deep Water ########
"""
#set it up so its exact floating point
dx = 0.1
l = 0.1
dt = l*dx
startx = 0.0
endx = 300.0
startt = 0.0
endt = 45.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../swdata/deep/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
stage = 1.0
bot = 0.0
hadd = 0.01
uadd= 0.0
xr = 147.5
xl = 152.5
gap = 20
gapbig = gap * 5

con, bed = flatbump(x,stage,bot,hadd,uadd,xr,xl)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([0.0,0.0,0.0,0.0])
u1 = array([0.0,0.0,0.0,0.0])
    
h0 = stage - b0
h1 = stage - b1

for i in range(1,len(t)):      
        
    if(i % gapbig == 0 or i ==1):
        u = solveucon(con)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])         
               
            for j in range(n):
             
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
u = solveucon(con)
with open(s,'a') as file2:
    
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dt(s)','timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
     for j in range(n):
             
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""

##### Deep Water Soliton ########
"""
#set it up so its exact floating point
dx = 0.2
l = 0.02
dt = l*dx
startx = -500.0 
endx = 1500.0
startt = 0.0
endt = 100.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../swdata/soliton/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
stage = 10.0
a0 = 10.0
a1 = 1.0
bot = 0.0
hadd = 1.0
uadd= 0.0
xr = 48
xl = 52
gap = 20
gapbig = gap * 5

#con, bed = solitoninit1(n,stage,a0,a1,g,x,t[0],bot,dx)
con,bed = solitoninit(n,a0,a1,g,x,t[0],0,dx)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([0.0,0.0,0.0,0.0])
u1 = array([0.0,0.0,0.0,0.0])
    
h0 = stage - b0
h1 = stage - b1

for i in range(1,len(t)):      
        
    if(i % gapbig == 0 or i ==1):
        u = solveucon(con)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])         
               
            for j in range(n):
             
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
u = solveucon(con)
with open(s,'a') as file2:
    
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dt(s)','timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
     for j in range(n):
             
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""

##### Soliton over linear slope ########
a0 = 1.0
a1 = 0.01

g = 9.81
dx = 0.01
l = 0.1
dt = l*dx
beta = 2.0
startx = -200
endx = 200.0 + dx
startt = 0.0
endt = 100 + dt  

wdir = "../../../data/raw/swdata/lsolslopedxlower/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
m = len(t)

gap = 20
    
#slope properties
#beginning of slope
slopbeg = 100

#end of slope
slopend =149.5

#end height of bed
topbed = 0.99 

con,bed = soloverslope(x,a0,a1,slopbeg,slopend,topbed,g)



h = []
G = []
for j in range(n):
    h = append(h,con[j][0])
    G = append(G,con[j][1])
u = solveucon(con)


x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([0.0,0.0,0.0,0.0])
u1 = array([0.0,0.0,0.0,0.0])
    
h0 = a0 - b0
h1 = a0 - b1

for i in range(1,len(t)):      
        
    if(i % gap == 0 or i ==1):
        u = solveucon(con)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])         
               
            for j in range(n):
             
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]
    
s = wdir + "saveoutputtslast.txt"
u,h = splitcon(con)
with open(s,'a') as file2:
    
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dt(s)','timestep', 'height(m)', 'G' , 'u(m/s)' ])        
               
     for j in range(n):
             
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
