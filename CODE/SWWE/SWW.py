# -*- coding: utf-8 -*-
"""
Created on Mon Apr 07 17:31:45 2014

@author: Jordan
"""


from scipy import *
from pylab import plot, show, legend,ylim,xlim
    

def minmod(a,b,c):
    
    if (a>0 and b>0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0

def evolveGTMVLHPL(con,dt,dx,g,beta):
    #get averages
    idx = 1.0 / dx  
    ig = 1.0 / g
    i16 = 1.0 / 16
    avcon = con
    ncon = zeros((n,2))
    
    avcon[0] = avcon[1]
    avcon[-1] = avcon[-2]
    
    ncon[0] = avcon[0]
    ncon[-1] = avcon[-2]
    
    #flux through right
    #calculate the temporary gradients
    utr1 =  avcon[1][1] - avcon[0][1]   
    utr2 =  avcon[2][1] - avcon[1][1]
    utr3 =  0.5*(avcon[2][1] - avcon[0][1])
    htr1 =  avcon[1][0] - avcon[0][0]   
    htr2 =  avcon[2][0] - avcon[1][0]
    htr3 =  0.5*(avcon[2][0] - avcon[0][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    #flux through right
    ul = ugrad*0.5 + avcon[1][1]
    hl = hgrad*0.5 + avcon[1][0]
        
        
    #calculate the temporary gradients
    utr1 =  avcon[2][1] - avcon[1][1]   
    utr2 =  avcon[3][1] - avcon[2][1]
    utr3 =  0.5*(avcon[3][1] - avcon[1][1])
    htr1 =  avcon[2][0] - avcon[1][0]   
    htr2 =  avcon[3][0] - avcon[2][0]
    htr3 =  0.5*(avcon[3][0] - avcon[1][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    ur = -ugrad*0.5 + avcon[2][1]
    hr = -hgrad*0.5 + avcon[2][0]
        
    sqrtghl = sqrt(g*hl)
    sqrtghr = sqrt(g*hr)
        
    um =0.5*(ur + ul) + sqrtghl - sqrtghr
    hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
    sqrtghm = sqrt(g*hm)
        
    sl = min(ul - sqrtghl, um - sqrtghm)
    sr = max(ur + sqrtghr, um + sqrtghm)
        
    fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
    fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
    if (sl >= 0):
            
        fo = fl 
            
    elif(sr <= 0):
         fo = fr
            
    else:
         srmsl = 1.0 / (sr - sl)
         fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[2] -avcon[1] ))
            

    #calculate the temporary gradients
    utr1 =  avcon[0][1] - avcon[0][1]   
    utr2 =  avcon[1][1] - avcon[0][1]
    utr3 =  0.5*(avcon[1][1] - avcon[0][1])
    htr1 =  avcon[0][0] - avcon[0][0]   
    htr2 =  avcon[1][0] - avcon[0][0]
    htr3 =  0.5*(avcon[1][0] - avcon[0][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    #flux through left
    ul = ugrad*0.5 + avcon[0][1]
    hl = hgrad*0.5 + avcon[0][0]
        
        
    #calculate the temporary gradients
    utr1 =  avcon[1][1] - avcon[0][1]   
    utr2 =  avcon[2][1] - avcon[1][1]
    utr3 =  0.5*(avcon[2][1] - avcon[0][1])
    htr1 =  avcon[1][0] - avcon[0][0]   
    htr2 =  avcon[2][0] - avcon[1][0]
    htr3 =  0.5*(avcon[2][0] - avcon[0][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    ur = -ugrad*0.5 + avcon[1][1]
    hr = -hgrad*0.5 + avcon[1][0]
        
    sqrtghl = sqrt(g*hl)
    sqrtghr = sqrt(g*hr)
        
    um =0.5*(ur + ul) + sqrtghl - sqrtghr
    hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
    sqrtghm = sqrt(g*hm)
        
    sl = min(ul - sqrtghl, um - sqrtghm)
    sr = max(ur + sqrtghr, um + sqrtghm)
        
    fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
    fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
    if (sl >= 0):
        fi = fl 
    elif(sr <= 0):
         fi = fr            
    else:
        srmsl = 1.0 / (sr - sl)
        fi = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[1] -avcon[0] ))
            
    ncon[1] = avcon[1] - dt*idx*(fo - fi)
    
    fi = fo
    
    for i in range(2, len(avcon) - 2):
        #DO INTERIOR FIRST AS TEST OF METHOD
        #ONLY NEED TO CHANGE INTERFACE VALUES

        #calculate the temporary gradients
        utr1 =  avcon[i][1] - avcon[i-1][1] 
        utr2 =  avcon[i+1][1] - avcon[i][1]
        utr3 =  0.5*(avcon[i+1][1] - avcon[i-1][1])
        htr1 =  avcon[i][0] - avcon[i-1][0]   
        htr2 =  avcon[i+1][0] - avcon[i][0]
        htr3 =  0.5*(avcon[i+1][0] - avcon[i-1][0])
        
        ugrad = minmod(beta*utr1,beta*utr2,utr3)
        hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
        #print(ugrad,hgrad)
        
        #flux through right
        ul = ugrad*0.5 + avcon[i][1]
        hl = hgrad*0.5 + avcon[i][0]
        
        
        #calculate the temporary gradients
        utr1 =  avcon[i+1][1] - avcon[i][1]   
        utr2 =  avcon[i+2][1] - avcon[i+1][1]
        utr3 =  0.5*(avcon[i+2][1] - avcon[i][1])
        htr1 =  avcon[i+1][0] - avcon[i][0]   
        htr2 =  avcon[i+2][0] - avcon[i+1][0]
        htr3 =  0.5*(avcon[i+2][0] - avcon[i][0])
        
        ugrad = minmod(beta*utr1,beta*utr2,utr3)
        hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
        ur = -ugrad*0.5 + avcon[i+1][1]
        hr = -hgrad*0.5 + avcon[i+1][0]
        
        sqrtghl = sqrt(g*hl)
        sqrtghr = sqrt(g*hr)
        
        um =0.5*(ur + ul) + sqrtghl - sqrtghr
        hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
        sqrtghm = sqrt(g*hm)
        
        sl = min(ul - sqrtghl, um - sqrtghm)
        sr = max(ur + sqrtghr, um + sqrtghm)
        
        fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
        fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
        if (sl >= 0):
            fo = fl 
        elif(sr <= 0):
             fo = fr
        else:
            srmsl = 1.0 / (sr - sl)
            fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[i+1] -avcon[i] ))
            
            
        ncon[i] = avcon[i] - dt*idx*(fo - fi)
        fi = fo
    
    #flux through right
    #calculate the temporary gradients
    utr1 =  avcon[-2][1] - avcon[-3][1]   
    utr2 =  avcon[-1][1] - avcon[-2][1]
    utr3 =  0.5*(avcon[-1][1] - avcon[-3][1])
    htr1 =  avcon[-2][0] - avcon[-3][0]   
    htr2 =  avcon[-1][0] - avcon[-2][0]
    htr3 =  0.5*(avcon[-1][0] - avcon[-3][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    #flux through right
    ul = ugrad*0.5 + avcon[-2][1]
    hl = hgrad*0.5 + avcon[-2][0]
        
        
    #calculate the temporary gradients
    utr1 =  avcon[-1][1] - avcon[-2][1]   
    utr2 =  avcon[-1][1] - avcon[-1][1]
    utr3 =  0.5*(avcon[-1][1] - avcon[-2][1])
    htr1 =  avcon[-1][0] - avcon[-2][0]   
    htr2 =  avcon[-1][0] - avcon[-1][0]
    htr3 =  0.5*(avcon[-1][0] - avcon[-2][0])
        
    ugrad = minmod(beta*utr1,beta*utr2,utr3)
    hgrad = minmod(beta*htr1,beta*htr2,htr3)
        
    ur = -ugrad*0.5 + avcon[-1][1]
    hr = -hgrad*0.5 + avcon[-1][0]
        
    sqrtghl = sqrt(g*hl)
    sqrtghr = sqrt(g*hr)
        
    um =0.5*(ur + ul) + sqrtghl - sqrtghr
    hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
    sqrtghm = sqrt(g*hm)
        
    sl = min(ul - sqrtghl, um - sqrtghm)
    sr = max(ur + sqrtghr, um + sqrtghm)
        
    fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
    fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
    if (sl >= 0):
            
        fo = fl 
            
    elif(sr <= 0):
         fo = fr
            
    else:
         srmsl = 1.0 / (sr - sl)
         fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[-1] -avcon[-2] ))
         
    ncon[-2] = avcon[-2] - dt*idx*(fo - fi)

    return ncon


def evolveGTMVLH(con,dt,dx,g):
    #get averages
    idx = 1.0 / dx  
    ig = 1.0 / g
    i16 = 1.0 / 16
    avcon = con
    ncon = zeros((n,2))
    
    avcon[0] = avcon[1]
    avcon[-1] = avcon[-2]
    
    ncon[0] = avcon[0]
    ncon[-1] = avcon[-2]
    
    #flux through right
    ul = avcon[1][1]
    hl = avcon[1][0]
    ur = avcon[2][1]
    hr = avcon[2][0]
        
    sqrtghl = sqrt(g*hl)
    sqrtghr = sqrt(g*hr)
        
    um =0.5*(ur + ul) + sqrtghl - sqrtghr
    hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
    sqrtghm = sqrt(g*hm)
        
    sl = min(ul - sqrtghl, um - sqrtghm)
    sr = max(ur + sqrtghr, um + sqrtghm)
        
    fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
    fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
    if (sl >= 0):
            
        fo = fl 
            
    elif(sr <= 0):
         fo = fr
            
    else:
         srmsl = 1.0 / (sr - sl)
         fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[2] -avcon[1] ))
            
    ul = avcon[0][1]
    hl = avcon[0][0]
    ur = avcon[1][1]
    hr = avcon[1][0]
        
    sqrtghl = sqrt(g*hl)
    sqrtghr = sqrt(g*hr)
        
    um =0.5*(ur + ul) + sqrtghl - sqrtghr
    hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
    sqrtghm = sqrt(g*hm)
        
    sl = min(ul - sqrtghl, um - sqrtghm)
    sr = max(ur + sqrtghr, um + sqrtghm)
        
    fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
    fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
    if (sl >= 0):
        fi = fl 
    elif(sr <= 0):
         fi = fr            
    else:
        srmsl = 1.0 / (sr - sl)
        fi = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[1] -avcon[0] ))
            
    ncon[1] = avcon[1] - dt*idx*(fo - fi)
    
    fi = fo
    
    for i in range(2, len(avcon) - 1):
        #flux through right
        ul = avcon[i][1]
        hl = avcon[i][0]
        ur = avcon[i+1][1]
        hr = avcon[i+1][0]
        
        sqrtghl = sqrt(g*hl)
        sqrtghr = sqrt(g*hr)
        
        um =0.5*(ur + ul) + sqrtghl - sqrtghr
        hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
        sqrtghm = sqrt(g*hm)
        
        sl = min(ul - sqrtghl, um - sqrtghm)
        sr = max(ur + sqrtghr, um + sqrtghm)
        
        fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
        fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
        if (sl >= 0):
            fo = fl 
        elif(sr <= 0):
             fo = fr
        else:
            srmsl = 1.0 / (sr - sl)
            fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[i+1] -avcon[i] ))
            
            
        ncon[i] = avcon[i] - dt*idx*(fo - fi)
        fi = fo

    return ncon


def evolveGTMVL(con,dt,dx,g):
    #get averages
    idx = 1.0 / dx  
    ig = 1.0 / g
    i16 = 1.0 / 16
    avcon = con
    ncon = zeros((n,2))
    
    avcon[0] = avcon[1]
    avcon[-1] = avcon[-2]
    
    ncon[0] = avcon[0]
    ncon[-1] = avcon[-2]
    
    for i in range(1, len(avcon) - 1):
        #flux through right
        ul = avcon[i][1]
        hl = avcon[i][0]
        ur = avcon[i+1][1]
        hr = avcon[i+1][0]
        
        sqrtghl = sqrt(g*hl)
        sqrtghr = sqrt(g*hr)
        
        um =0.5*(ur + ul) + sqrtghl - sqrtghr
        hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
        sqrtghm = sqrt(g*hm)
        
        sl = min(ul - sqrtghl, um - sqrtghm)
        sr = max(ur + sqrtghr, um + sqrtghm)
        
        fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
        fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
        if (sl >= 0):
            
            fo = fl 
            
        elif(sl < 0 and sr > 0):
                        
            srmsl = 1.0 / (sr - sl)
            fo = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[i+1] -avcon[i] ))
            
        elif(sr <= 0):
            fo = fr
            
        ul = avcon[i-1][1]
        hl = avcon[i-1][0]
        ur = avcon[i][1]
        hr = avcon[i][0]
        
        sqrtghl = sqrt(g*hl)
        sqrtghr = sqrt(g*hr)
        
        um =0.5*(ur + ul) + sqrtghl - sqrtghr
        hm = i16*ig*(ul + 2*sqrtghl - ur + sqrtghr)*(ul + 2*sqrtghl - ur + sqrtghr)
        sqrtghm = sqrt(g*hm)
        
        sl = min(ul - sqrtghl, um - sqrtghm)
        sr = max(ur + sqrtghr, um + sqrtghm)
        
        fl = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl]) 
        fr = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])
        
        if (sl >= 0):
            
            fi = fl 
            
        elif(sl < 0 and sr > 0):
                        
            srmsl = 1.0 / (sr - sl)
            fi = srmsl*(sr*fl - sl*fr + sl*sr*(avcon[i] -avcon[i-1] ))
            
        elif(sr <= 0):
            fi = fr
            
        ncon[i] = avcon[i] - dt*idx*(fo - fi)
        
       
    
    
    
    return ncon

def evolveGTM(con,dt,dx,g):
    #get averages
    idx = 1.0 / dx  
    avcon = con
    ncon = zeros((n,2))
    
    """
    for i in range(1, len(con) - 1 ):
        avcon[i] = 0.5* (con[i] + con[i+1])
        
    """
        
    avcon[0] = avcon[1]
    avcon[-1] = avcon[-2]
    
    ncon[0] = avcon[0]
    ncon[-1] = avcon[-2]
    
    
        
    #calculate flux using these averages as the left and right
    for i in range(1, len(avcon) - 1):
        
        #flux through right
        ul = avcon[i][1]
        hl = avcon[i][0]
        ur = avcon[i+1][1]
        hr = avcon[i+1][0]
        
        sqrthl = sqrt(hl)
        sqrthr = sqrt(hr)
        
        uavg = (sqrthl * ul + sqrthr*ur) / (sqrthl + sqrthr)
        cavg = sqrt(g*0.5*(hl + hr))
        icavg = 1 / cavg
        havg = sqrthl*sqrthr
        
        l1 = uavg - cavg
        l2 = uavg + cavg
        
        a1 = 0.5* (hr  - hl - havg*icavg*(ur - ul))
        a2 = 0.5* (hr  - hl + havg*icavg*(ur - ul))
        
        e1 = array([1, l1])
        e2 = array([1, l2])
        
        
        ful = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl])  
        fur = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])  
        fr = 0.5* (ful + fur) - 0.5*(abs(l1)*a1*e1 + abs(l2)*a2*e2)
        
        #flux through left
        
        ul = avcon[i-1][1]
        hl = avcon[i-1][0]
        ur = avcon[i][1]
        hr = avcon[i][0]
        
        sqrthl = sqrt(hl)
        sqrthr = sqrt(hr)
        
        uavg = (sqrthl * ul + sqrthr*ur) / (sqrthl + sqrthr)
        cavg = sqrt(g*0.5*(hl + hr))
        icavg = 1 / cavg
        havg = sqrthl*sqrthr
        
        l1 = uavg - cavg
        l2 = uavg + cavg
        
        a1 = 0.5* (hr  - hl - havg*icavg*(ur - ul))
        a2 = 0.5* (hr  - hl + havg*icavg*(ur - ul))
        
        e1 = array([1, l1])
        e2 = array([1, l2])
        
        
        ful = array([ul*hl, ul*ul*hl + 0.5*g*hl*hl])  
        fur = array([ur*hr, ur*ur*hr + 0.5*g*hr*hr])  
        fl = 0.5* (ful + fur) - 0.5*(abs(l1)*a1*e1 + abs(l2)*a2*e2)
        
        ncon[i] = avcon[i] - dt*idx*(fr - fl)
        
    return ncon
    
    
        
def evolveUP(con,dt,dx,g):
       
    idx = 1.0 / dx    
    
    #main loop
	
    ncon = zeros((n,2))
	
    #ghost cell at front
 
    ncon[0] = array([con[0][0], con[1][1]])
	
    sqrt2 = sqrt(2)
    isqrt2 = 1 / sqrt2
    for i in range(1,len(con) -1):
        #all flow left to right
        h = con[i][0]
        uh =con[i][1]
        u = (uh*1.0) / h
        c = sqrt(abs(g*h))
        G = array([[0 ,1],[-(u*u) + 0.5*g*h,2*u]])
        Da = array([[abs(u + c*isqrt2) ,0],[0,abs(u - c*isqrt2)]])
        P = array([[1 ,1],[u + c*isqrt2,u - c*isqrt2]])
        Ga = dot(P,dot(Da,inv(P)))
        F = dot(G,con[i])
        Fa = dot(Ga,con[i])
        
        h1 = con[i+1][0]
        uh1 =con[i+1][1]
        u1 = (uh1*1.0) / h1
        c1 = sqrt(abs(g*h1))
        G1 = array([[0 ,1],[-(u1*u1) + 0.5*g*h1,2*u1]])
        Da1 = array([[abs(u1 + c1*isqrt2) ,0],[0,abs(u1 - c1*isqrt2)]])
        P1 = array([[1 ,1],[u1 + c1*isqrt2,u1 - c1*isqrt2]])
        Ga1 = dot(P1,dot(Da1,inv(P1)))
        F1 = dot(G1,con[i+1])
        Fa1 = dot(Ga1,con[i+1])
                
        fluxi = 0.5*(F + F1) + 0.5*(Fa - Fa1)
        
        h1 = con[i-1][0]
        uh1 =con[i-1][1]
        u1 = (uh1*1.0) / h1
        c1 = sqrt(abs(g*h1))
        G1 = array([[0 ,1],[-(u1*u1) + 0.5*g*h1,2*u1]])
        Da1 = array([[abs(u1 + c1*isqrt2) ,0],[0,abs(u1 - c1*isqrt2)]])
        P1 = array([[1 ,1],[u1 + c1*isqrt2,u1 - c1*isqrt2]])
        Ga1 = dot(P1,dot(Da1,inv(P1)))
        F1 = dot(G1,con[i-1])
        Fa1 = dot(Ga1,con[i-1])
        
        fluxo = 0.5*(F + F1) + 0.5*(Fa1 - Fa)
        
        ncont = con[i] - dt*idx*(fluxi-fluxo)
        ncon[i] = ncont
        
    #upwindscheme
        
    #ghost cells   
    ncon[-1] = array([con[-1][0], con[-2][1]])
    
    return ncon
    
    
def makevar(sx,ex,dx,st,et,dt):
    
    x = arange(startx,endx,dx)
    t = arange(startt,endt,dt)
    
    return x,t    
    
def dambreak(n,divn,normheight,damheight):
    #returns the initial values for the dam breka problem
    
    con = zeros((n,2))
    
    bx = zeros(n)
    
    dampos = n / divn
    
    #initial conditions on h
    # dam break, set all heights behind dam wall to value
    for i in range(1,dampos-1):
        con[i][0] = damheight
    
    for i in range(dampos-1,len(con)-1):
        con[i][0]  = normheight
        
    con[0][0] = damheight
    con[-1][0] = normheight
    con[0][1] = 0.0
    con[-1][1] = 0.0
    
    return con,bx
    
    
    


#set it up so its exact floating point
dx = 0.1
l = 0.1
dt = l*dx
startx = 0.0
endx = 10.0
startt = 0.0
endt = 30

g = 9.8

from math import sqrt

sqrt2 = sqrt(2)

x,t = makevar(startx,endx,dx,startt,endt,dt)


n = len(x)

con,bx = dambreak(n,2,1,2)
con1,bx = dambreak(n,2,1,2)

#u,h,bx = Flat(n,2)
#u,h,bx = sine(x,n,2.0,0.5)

#con = evolveGTM(con,dt,dx,g)


#should be t
for ts in t:
    #con = evolveUP(con,dt,dx,g) 
    con = evolveGTMVLHPL(con,dt,dx,g,1.0)
    con1 = evolveGTMVLHPL(con1,dt*0.5,dx,g,1.0)
    cont = evolveGTMVLHPL(con1,dt,dx,g,1.0)
    con1 = 0.5*(cont + con1)
"""
    h = []
    if (int(ts/dt) % 100 == 0):
        for i in range(n):
            h = append(h,con[i][0])
        
        plot(x,h)
        show()
    
"""

h = []
h1 = []
for i in range(n):
    h = append(h,con[i][0])
    h1 = append(h1,con1[i][0])
    
    
plot(x,h1,label="2")
plot(x,h, label="1")
legend()
show()




