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
            
def evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,u0h,h0h,beta,dx,dt): 
    n = len(con)
    ncon = zeros((n,2))
    #update h' and G'    
    conp = evolve(con,bed,dt,dx,g,beta,u0,u1,h0,h1,b0,b1)
        
    #update h'' and G''
    conpp = evolve(conp,bed,dt,dx,g,beta,u0h,u1,h0h,h1,b0,b1)
    
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
  
def flat(x):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    for i in range(n):        
        con[i][0] = 0.4
        con[i][1] = 0
         
    return con,bx

def BejiFlume(x):
    """
    BejiFlume: gives the bed profile given the spatial array x
    
        Input:
            x : spatial grid
            
        Output:
            h : array of h(x_i) which is the water depth at the x[i] grid points
            u : array of u(x_i) which is the water velocity at the x[i] grid points
            bed: array of b(x_i) which is the bed profile at the x[i] grid points
    """
    
    # intialise our variables as lists same size as x
    n = len(x)
    bed = zeros(n)
    con = zeros((n,2))
    
    
    for i in range(n):

        #Define the bed profile with if statements as is piecewise.
        # h at x[i] is just defined so that h + bed is constant at 0.4m
        if(0 <= x[i] <= 6):
            bed[i] = 0.0
            con[i][0] = 0.4
        elif(6 < x[i] <= 12):
            bed[i] = 0.05*(x[i] - 6)
            con[i][0] = 0.4 - bed[i]
        elif(12 < x[i] <= 14):
            bed[i] = 0.3
            con[i][0] = 0.1
        elif(14 < x[i] <= 17):
            bed[i] = 0.3 - 0.1*(x[i] - 14)
            con[i][0] = 0.4 - bed[i]
        elif(17 < x[i] <= 18.95):
            bed[i] = 0.0
            con[i][0] = 0.4 - bed[i]
        elif(18.95 < x[i] <= 23.95):
            bed[i] = (0.2/5.0)*(x[i] - 18.95)
            con[i][0] = 0.4 - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            con[i][0] = 0.4 - bed[i]
        else:
            bed[i] = 0.0
            con[i][0] = 0.4  - bed[i]
            
    return con,bed
    
def BejiEdge(x,hc0,vc0,ft):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    hb = 0.4

    
    i = n-1
    et = ft
    #c1 = sqrt(g*(hb+ et))
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496 #sqrt(9.81*hb) * sqrt(3.0 / (k*k*hb*hb+ 3))
    ut = (c1*et) / (h1)
    eta[i] = et
    v[i] = ut
    
    #linear extrapolation
    i = n - 2
    et = 2*ft - (hc0 - hb)
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496
    ut = (c1*et) / (h1)
    #c1 = sqrt(g*(hb+ et))
    #ut = (c1*et) / (hb + et)
    eta[i] = et
    v[i] = ut
    
    for i in range(n-3,-1,-1):
        et = 2*(eta[i+1]) - (eta[i+2])
        #c1 = sqrt(g*(hb+ et))
        #ut = (c1*et) / (hb + et)
        
        h1 = hb + et
        c1 = sqrt(g*(hb+ et)) #1.641496
        ut = (c1*et) / (h1)
        eta[i] = et
        v[i] = ut
        
    i = -1
    et = 2*(eta[i+1]) - (eta[i+2])
    h1 = hb + et
    c1 = sqrt(g*(hb+ et)) #1.641496
    ut = (c1*et) / (h1)
    #c1 = sqrt(g*(hb+ et))
    #ut = (c1*et) / (hb + et)
    e0 = et
    v0 = ut 
    h0 = hb + e0
    
    hv = hb+ eta
    return hv,v

    
def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)


wdir = "../../../data/raw/swdata/BCTestSL/"
expdir = "../../../data/Experimental/Data 1994 Paper/CSV/"
exp = "sl"

if not os.path.exists(wdir):
    os.makedirs(wdir)
"""
sr = 0.039312

g = 9.81
dx = 0.05
l = 0.1
dt = sr / 4.0
beta = 2.0
startx = 5.7 + dx
endx = 100.0 + dx
startt = 0.0
endt = 60 + dt  

tts = []
hts = []    
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nwg1s = []

#Experimental data
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

    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
m = len(t)

gap = 20
    

con,bed = BejiFlume(x)



h = []
G = []
for j in range(n):
    h = append(h,con[j][0])
    G = append(G,con[j][1])
u = solveucon(con)

nBC = 4
x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = bed[0] * ones(nBC)
b1 = bed[-1] * ones(nBC)
u0 = zeros(nBC)
u1 = zeros(nBC)
h0 = h[0]*ones(nBC)
h1 = h[-1]*ones(nBC)

ct = dt
mp = int(ct/sr)
ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp])

h0h ,u0h = BejiEdge(x0,h[0],u[0],ftc) 

wg2i = int((10.5 - startx) / dx ) #good one
wg3i = int((12.5 - startx) / dx ) #G
wg4i = int((13.5 - startx) / dx ) #G
wg5i = int((14.5 - startx) / dx ) + 1 #
wg6i = int((15.7 - startx) / dx ) + 1
wg7i = int((17.3 - startx) / dx )
nwg1s.append(0.4)
nwg2s.append(0.4)  
nwg3s.append(0.4)  
nwg4s.append(0.4)  
nwg5s.append(0.4)  
nwg6s.append(0.4)  
nwg7s.append(0.4) 


for i in range(1,len(t)):        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,u0h,h0h,beta,dx,dt)
    
    
    ct = t[i]
    mp = int(ct/sr)
    ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct -ts[mp])
    h0,u0 = BejiEdge(x0,con[0][0],con[0][1]/con[0][0],ftc) 
    
    ct = t[i] + dt
    mp = int(ct/sr)
    ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct -ts[mp])    
    h0h,u0h = BejiEdge(x0,con[0][0],con[0][1]/con[0][0],ftc) 
    
    nwg1s.append((h0[-1] + b0[-1]))
    nwg2s.append((con[wg2i][0] + bed[wg2i]))
    nwg3s.append((con[wg3i][0] + bed[wg3i]))
    nwg4s.append((con[wg4i][0] + bed[wg4i])) 
    nwg5s.append((con[wg5i][0] + bed[wg5i]))
    nwg6s.append((con[wg6i][0] + bed[wg6i]))
    nwg7s.append((con[wg7i][0] + bed[wg7i]))
    
    print t[i]
    print con[200]


h = []
G = []
for j in range(n):
    h = append(h,con[j][0])
    G = append(G,con[j][1])
u = solveucon(con)
"""

s = wdir + "NumWaveGauge.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["Time(s)","WG1(m)","WG2(m)", "WG3(m)","WG4(m)","WG5(m)","WG6(m)","WG7(m)","WG8(m)"]) 
    
    for j in range(len(t)):
        writefile.writerow([str(t[j]), str((nwg1s[j] - 0.4)/100.0), str((nwg2s[j] - 0.4)/100.0), str((nwg3s[j] - 0.4)/100.0), str((nwg4s[j] - 0.4)/100.0), str((nwg5s[j] - 0.4)/100.0), str((nwg6s[j] - 0.4)/100.0),str((nwg7s[j] - 0.4)/100.0)]) 


s = wdir + "WaveGauge.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["Time(s)","WG1(m)","WG2(m)", "WG3(m)","WG4(m)","WG5(m)","WG6(m)","WG7(m)","WG8(m)"]) 
    
    for j in range(len(ts)):
        writefile.writerow([str(ts[j]),str(wg1s[j]), str(wg2s[j]), str(wg3s[j]), str(wg4s[j]), str(wg5s[j]), str(wg6s[j]),str(wg7s[j])]) 
