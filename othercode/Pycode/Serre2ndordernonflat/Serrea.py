from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv

def minmod(a,b,c):
    if (a>0 and b>0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0
          

#algorithm for solving tridiagonal systems from wiki page    
def TDMA(a,b,c,d):
    n = len(d)
    alpha = []
    beta = []
    x = [0]*n
    
    alpha.append((1.0*c[0])/b[0])
    beta.append((1.0*d[0])/b[0] )  
 
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha.append(c[i]* m)
        beta.append((d[i] - a[i-1]*beta[i-1]) * m)
        
    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta.append((d[n-1] - a[n-2]*beta[n-2]) * m)  

    x[n-1] = beta[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
 
    return array(x)

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

#gives exact up to linears, so is second order accurate huzzah
def getufromG(con,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
    n = len(con)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    
    G = zeros(n)
    
    for i in range(n):
        G[i] = con[i][1]
    
    
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        
        a[i-1] = ai
        b[i] =  bi
        c[i] = ci
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    c[i] = ci
    b[i] = bi
    
    G[i] = G[i] - u0*ai
    
    #i = n-1
    i = n-1
    
    
    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    a[i-1] = ai
    b[i] = bi
    G[i] = G[i] - u1*ci
    
    u = TDMA(a,b,c,G)
        
    return u 

#gives exact up to linears, so is second order accurate huzzah    
def getGfromu(con,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(con)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 
    


def evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt): 
    n = len(con)
    ncon = zeros((n,2))
    #update h' and G'    
    conp = evolve(con,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt)
        
    #update h'' and G''
    conpp = evolve(conp,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt)
    
    ncon = 0.5 * (con + con)
    
    return ncon
    
    
def evolve(con,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt):
    #get averages
    idx = 1.0 / dx  
    ithree = 1.0 / 3.0
    n = len(con)
    
    ncon = zeros((n,2))
    
    #calculate u at time step
    u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    
    #boundaries
    #beginning
    #i=-1
    i = -1
    h = h0[-1]
    hx = 0.5*idx*(con[i+1][0] - h0[-2])
    bx = 0.5*idx*(bed[i+1] - b0[-2])
    bxx = idx*idx*(bed[i+1] - 2*b0[-1] + b0[-2])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb3 = ai*u0[-2] + bi*u0[-1] +ci*u[i+1]
    
    #i=-2
    h = h0[-2]
    hx = 0.5*idx*(h0[-1] - h0[-3])
    bx = 0.5*idx*(b0[-1] - b0[-3])
    bxx = idx*idx*(b0[-1] - 2*b0[-2] + b0[-3])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb2 = ai*u0[-3] + bi*u0[-2] +ci*u0[-1]
    
    #i=-3
    h = h0[-3]
    hx = 0.5*idx*(h0[-2] - h0[-4])
    bx = 0.5*idx*(b0[-2] - b0[-4])
    bxx = idx*idx*(b0[-2] - 2*b0[-3] + b0[-4])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb1 = ai*u0[-4] + bi*u0[-3] +ci*u0[-2]
    
    #i = n
    i = n
    h = h1[0]
    hx = 0.5*idx*(h1[1] - con[i-1][0])
    bx = 0.5*idx*(b1[1] - bed[i-1])
    bxx = idx*idx*(b1[1] - 2*b1[0] + bed[i-1])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge1 = ai*u[i-1] + bi*u1[0] + ci*u1[1]
    
    #i = n+1
    h = h1[1]
    hx = 0.5*idx*(h1[2] - h1[0])
    bx = 0.5*idx*(b1[2] - b1[0])
    bxx = idx*idx*(b1[2] - 2*b1[1] + b1[0])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge2 = ai*u1[0] + bi*u1[1] +ci*u1[2]
    
    #i = n+2    
    h = h1[2]
    hx = 0.5*idx*(h1[3] - h1[1])
    bx = 0.5*idx*(b1[3] - b1[1])
    bxx = idx*idx*(b1[3] - 2*b1[2] + b1[1])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge3 = ai*u1[1] + bi*u1[2] +ci*u1[3]
    
    
    conlb = zeros((3,2))
    conrb = zeros((3,2))
    ubeg = u0[1:4]
    uend = u1[0:3]
    bbeg = b0[1:4]
    bend = b1[0:3]
    
    
    conlb[0][0] = h0[-3]
    conlb[1][0] = h0[-2]
    conlb[2][0] = h0[-1]
    conlb[0][1] = gb1
    conlb[1][1] = gb2
    conlb[2][1] = gb3
    
    conrb[0][0] = h1[0]
    conrb[1][0] = h1[1]
    conrb[2][0] = h1[2]
    conrb[0][1] = ge1
    conrb[1][1] = ge2
    conrb[2][1] = ge3
    
    con = concatenate([conlb, con, conrb])
    bed = concatenate([bbeg,bed,bend])
    u = concatenate([ubeg,u,uend])
        
    #do normal stuff 
        
    #i = 2
    i = 2
    #define the stage
    wi = con[i][0] + bed[i]
    wip1 = con[i+1][0] + bed[i+1]
    wip2 = con[i+2][0] + bed[i+2]
    wip3 = con[i+3][0] + bed[i+3]
    wim1 = con[i-1][0] + bed[i-1]
    wim2 = con[i-2][0] + bed[i-2] 
        
    #reconstruct common values first
        
    #i left and right
        
    #gradients
    dwib = (wi - wim1)
    dwif = (wip1 - wi)
    dwim = 0.5*(wip1 - wim1)
    dhib = (con[i][0] - con[i-1][0])
    dhif = (con[i+1][0] - con[i][0])
    dhim = 0.5*(con[i+1][0] - con[i-1][0])
    dGib = (con[i][1] - con[i-1][1])
    dGif = (con[i+1][1] - con[i][1])
    dGim = 0.5*(con[i+1][1] - con[i-1][1])
    duib = (u[i] - u[i-1])
    duif = (u[i+1] - u[i])
    duim = 0.5*(u[i+1] - u[i-1])
        
    #limiting
    dwi = minmod(beta*dwib,beta*dwif,dwim)
    dhi = minmod(beta*dhib,beta*dhif,dhim)
    dGi = minmod(beta*dGib,beta*dGif,dGim)
    dui = minmod(beta*duib,beta*duif,duim)
        
    #reconstruct right
    hir = con[i][0] + 0.5*dhi
    wir = wi + 0.5*dwi
    Gir = con[i][1] + 0.5*dGi
    uir = u[i] + 0.5*dui
    bir = wir - hir
        
    #reconstruct left
    hil = con[i][0] - 0.5*dhi
    wil = wi - 0.5*dwi
    Gil = con[i][1] - 0.5*dGi
    uil = u[i] - 0.5*dui
    bil = wil - hil
        
    #only left of i+1 common but do both
        
    #gradients
    dwip1b = (wip1 - wi)
    dwip1f = (wip2 - wip1)
    dwip1m = 0.5*(wip2 - wi)
    dhip1b = (con[i+1][0] - con[i][0])
    dhip1f = (con[i+2][0] - con[i+1][0])
    dhip1m = 0.5*(con[i+2][0] - con[i][0])
    dGip1b = (con[i+1][1] - con[i][1])
    dGip1f = (con[i+2][1] - con[i+1][1])
    dGip1m = 0.5*(con[i+2][1] - con[i][1])
    duip1b = (u[i+1] - u[i])
    duip1f = (u[i+2] - u[i+1])
    duip1m = 0.5*(u[i+2] - u[i])
        
    #limiting
    dwip1 = minmod(beta*dwip1b,beta*dwip1f,dwip1m)
    dhip1 = minmod(beta*dhip1b,beta*dhip1f,dhip1m)
    dGip1 = minmod(beta*dGip1b,beta*dGip1f,dGip1m)
    duip1 = minmod(beta*duip1b,beta*duip1f,duip1m)
        
    #reconstruct right
    hip1r = con[i+1][0] + 0.5*dhip1
    wip1r = wip1 + 0.5*dwip1
    Gip1r = con[i+1][1] + 0.5*dGip1
    uip1r = u[i+1] + 0.5*duip1
    bip1r = wip1r - hip1r
        
    #reconstruct left
    hip1l = con[i+1][0] - 0.5*dhip1
    wip1l = wip1 - 0.5*dwip1
    Gip1l = con[i+1][1] - 0.5*dGip1
    uip1l = u[i+1] - 0.5*duip1
    bip1l = wip1l - hip1l
        
        
    #only right of i-1
    #i-1  right
        
    #gradients
    dwim1b = (wim1 - wim2)
    dwim1f = (wi - wim1)
    dwim1m = 0.5*(wi - wim2)
    dhim1b = (con[i-1][0] - con[i-2][0])
    dhim1f = (con[i][0] - con[i-1][0])
    dhim1m = 0.5*(con[i][0] - con[i-2][0])
    dGim1b = (con[i-1][1] - con[i-2][1])
    dGim1f = (con[i][1] - con[i-1][1])
    dGim1m = 0.5*(con[i][1] - con[i-2][1])
    duim1b = (u[i-1] - u[i-2])
    duim1f = (u[i] - u[i-1])
    duim1m = 0.5*(u[i] - u[i-2])
        
    #limiting
    dwim1 = minmod(beta*dwim1b,beta*dwim1f,dwim1m)
    dhim1 = minmod(beta*dhim1b,beta*dhim1f,dhim1m)
    dGim1 = minmod(beta*dGim1b,beta*dGim1f,dGim1m)
    duim1 = minmod(beta*duim1b,beta*duim1f,duim1m)
        
    #reconstruct right
    him1r = con[i-1][0] + 0.5*dhim1
    wim1r = wim1 + 0.5*dwim1
    Gim1r = con[i-1][1] + 0.5*dGim1
    uim1r = u[i-1] + 0.5*duim1
    bim1r = wim1r - him1r
        
    #reconstruct i+2 left
    
    #gradients
    dwip2b = (wip2 - wip1)
    dwip2f = (wip3 - wip2)
    dwip2m = 0.5*(wip3 - wip1)
    dhip2b = (con[i+2][0] - con[i+1][0])
    dhip2f = (con[i+3][0] - con[i+2][0])
    dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
    dGip2b = (con[i+2][1] - con[i+1][1])
    dGip2f = (con[i+3][1] - con[i+2][1])
    dGip2m = 0.5*(con[i+3][1] - con[i+1][1])
    duip2b = (u[i+2] - u[i+1])
    duip2f = (u[i+3] - u[i+2])
    duip2m = 0.5*(u[i+3] - u[i+1])
        
    #limiting
    dwip2 = minmod(beta*dwip2b,beta*dwip2f,dwip2m)
    dhip2 = minmod(beta*dhip2b,beta*dhip2f,dhip2m)
    dGip2 = minmod(beta*dGip2b,beta*dGip2f,dGip2m)
    duip2 = minmod(beta*duip2b,beta*duip2f,duip2m)
                
    #reconstruct left
    hip2l = con[i+2][0] - 0.5*dhip2
    wip2l = wip2 - 0.5*dwip2
    Gip2l = con[i+2][1] - 0.5*dGip2
    uip2l = u[i+2] - 0.5*duip2
    bip2l = wip2l - hip2l
    
    #calculate forces              
        
    #right force i
    nbi = max(bip1l,bir)
    hihm = max(0,wir-nbi)
    hihp = max(0,wip1l-nbi)

    her = hihp
    Ger = Gip1l
    uer = uip1l
        
    hel = hihm
    Gel = Gir
    uel = uir
        
    duer = idx*(uip2l - uip1l)
    dber = idx*(bip2l - bip1l)
            
    duel = idx*(uir - uim1r)
    dbel = idx*(bir - bim1r)
        
    sqrtghel = sqrt(g*hel)
    sqrtgher = sqrt(g*her)
    sl = min(0,uel - sqrtghel, uer - sqrtgher)
    sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
    felh = uel*hel
    felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
    ferh = uer*her
    ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
  
          
    if(sr == 0 and sl ==0):
        foh = 0.0
        foG = 0.0
    else:
        isrmsl = 1.0 / (sr - sl)
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))
    
    fih = foh
    fiG = foG
    himhp = hihp
        
    him1r = hir
    wim1r = wir
    bim1r = bir
    Gim1r = Gir
    uim1r = uir
        
    hil = hip1l
    wil = wip1l
    bil = bip1l
    Gil = Gip1l
    uil = uip1l
        
    hir = hip1r
    wir = wip1r
    bir = bip1r
    Gir = Gip1r
    uir = uip1r
    
    hip1l = hip2l
    wip1l = wip2l
    bip1l = bip2l
    Gip1l = Gip2l
    uip1l = uip2l       
        
    
    dhip1 = dhip2
    dwip1 = dwip2
    duip1 = duip2
    dGip1 = dGip2
    
    for i in range(3,len(con)-3):
        #update both forces at same time

        #define the stage
        wi = con[i][0] + bed[i]
        wip1 = con[i+1][0] + bed[i+1]
        wip2 = con[i+2][0] + bed[i+2]
        wip3 = con[i+3][0] + bed[i+3]
        wim1 = con[i-1][0] + bed[i-1]
        wim2 = con[i-2][0] + bed[i-2] 
        
        #reconstruct common values first
                
        #only left of i+1 common but do both        
        
        #reconstruct right
        hip1r = con[i+1][0] + 0.5*dhip1
        wip1r = wip1 + 0.5*dwip1
        Gip1r = con[i+1][1] + 0.5*dGip1
        uip1r = u[i+1] + 0.5*duip1
        bip1r = wip1r - hip1r
        
        
        #reconstruct i+2 left
        
        #gradients
        dwip2b = (wip2 - wip1)
        dwip2f = (wip3 - wip2)
        dwip2m = 0.5*(wip3 - wip1)
        dhip2b = (con[i+2][0] - con[i+1][0])
        dhip2f = (con[i+3][0] - con[i+2][0])
        dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
        dGip2b = (con[i+2][1] - con[i+1][1])
        dGip2f = (con[i+3][1] - con[i+2][1])
        dGip2m = 0.5*(con[i+3][1] - con[i+1][1])
        duip2b = (u[i+2] - u[i+1])
        duip2f = (u[i+3] - u[i+2])
        duip2m = 0.5*(u[i+3] - u[i+1])
        
        #limiting
        dwip2 = minmod(beta*dwip2b,beta*dwip2f,dwip2m)
        dhip2 = minmod(beta*dhip2b,beta*dhip2f,dhip2m)
        dGip2 = minmod(beta*dGip2b,beta*dGip2f,dGip2m)
        duip2 = minmod(beta*duip2b,beta*duip2f,duip2m)
                
        #reconstruct left
        hip2l = con[i+2][0] - 0.5*dhip2
        wip2l = wip2 - 0.5*dwip2
        Gip2l = con[i+2][1] - 0.5*dGip2
        uip2l = u[i+2] - 0.5*duip2
        bip2l = wip2l - hip2l


                
        #calculate forces              
        
        #right force i
        nbi = max(bip1l,bir)
        hihm = max(0,wir-nbi)
        hihp = max(0,wip1l-nbi)

        her = hihp
        Ger = Gip1l
        uer = uip1l
        
        hel = hihm
        Gel = Gir
        uel = uir
        
        duer = idx*(uip2l - uip1l)
        dber = idx*(bip2l - bip1l)
            
        duel = idx*(uir - uim1r)
        dbel = idx*(bir - bim1r)
        
        sqrtghel = sqrt(g*hel)
        sqrtgher = sqrt(g*her)
        sl = min(0,uel - sqrtghel, uer - sqrtgher)
        sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
        felh = uel*hel
        felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
        ferh = uer*her
        ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
  
        if(sr == 0 and sl ==0):
            foh = 0.0
            foG = 0.0
        else:
            isrmsl = 1.0 / (sr - sl)
            foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
            foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))
        
        
        #calculate the source term
        th = con[i][0]
        tu = u[i]
        tux = (uil - uir)
        tbx = (bil - bir)
        tbxx = idx*idx*(bed[i+1] - 2*bed[i] + bed[i-1])
        
        sourcer = g*0.5*(hihm*hihm - hir*hir)
        sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx       
        sourcel = g*0.5*(hil*hil - himhp*himhp)
        
        
        ncon[i-3][0] = con[i][0] - dt*idx*(foh - fih)
        ncon[i-3][1] = con[i][1] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + sourcec)
        
        fih = foh
        fiG = foG
        himhp = hihp
        
        him1r = hir
        wim1r = wir
        bim1r = bir
        Gim1r = Gir
        uim1r = uir
        
        hil = hip1l
        wil = wip1l
        bil = bip1l
        Gil = Gip1l
        uil = uip1l
        
        hir = hip1r
        wir = wip1r
        bir = bip1r
        Gir = Gip1r
        uir = uip1r
    
        hip1l = hip2l
        wip1l = wip2l
        bip1l = bip2l
        Gip1l = Gip2l
        uip1l = uip2l       
        
    
        dhip1 = dhip2
        dwip1 = dwip2
        duip1 = duip2
        dGip1 = dGip2
        

 
    return ncon
        
    
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
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        ue =  c* ((he - a0) / he)
        
        con[i][0] = he
        u[i] = ue
         
    G = getGfromu(con,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    for i in range(len(con)):
        con[i][1] = G[i]
    return con,bx
    
def solitoninit1(n,stage,a0,a1,g,x,t0,bot,dx):
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        ue =  c* ((he - a0) / he)
        
        con[i][0] = stage + (he-a0)
        u[i] = ue
         
    G = getGfromu(con,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    for i in range(len(con)):
        con[i][1] = G[i]
    return con,bx

def flatlake(x,dx,stage,normbot,vel,func,c1,c2,haddl,haddr,hadd):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        u[i] = vel
        bx[i] =  normbot + c1*func(c2*(x[i] - x[n/2]))
        con[i][0] = stage - bx[i]
        
        if(i >= haddl and i <= haddr):
            con[i][0] = con[i][0]+hadd
            
        b0 = normbot + c1*func(c2*(x[0] - dx - x[n/2]))
        b1 = normbot + c1*func(c2*(x[-1] + dx - x[n/2]))
        
    G = getGfromu(con,u,bx,u[0],u[-1],stage - b0,stage-b1,b0,b1,dx)
    for i in range(n):
        con[i][1] = G[i]
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
    u = zeros(n)
    
    for i in range(n):
        u[i] = 0.0
        
        if(x[i] >= xlb and x[i] <= xrb):
            bx[i] = bot2
        else:
            bx[i] = bot1
            
        con[i][0] = stage - bx[i]
            
        b0 = bot1
        b1 = bot1
        
    G = getGfromu(con,u,bx,u[0],u[-1],stage - b0,stage-b1,b0,b1,dx)
    for i in range(n):
        con[i][1] = G[i]
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
    u = zeros(n)
    
    for i in range(n):
       
        if(x[i] >= xc):
            con[i][0] = hl
        else:
            con[i][0] = hf
            
        bx[i] = bot
         
    G = getGfromu(con,u,bx,u[0],u[-1],hf,hl,bot,bot,dx)
    for i in range(n):
        con[i][1] = G[i]
        
    return con ,bx
    
def flowoverbump(x,stage,bot,bc,bm,ba,u0,u1,h0,h1):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    bg = (bm - bot)/ (ba*ba)
    
    for i in range(n):
        u[i] = u0
       
        if(x[i] >= bc - ba and x[i] <= bc + ba):
            bx[i] = bm - bg*(x[i] - bc)*(x[i] - bc)
            
        con[i][0] = stage - bx[i]
            
            

         
    G = getGfromu(con,u,bx,u0,u1,h0,h1,bot,bot,dx)
    for i in range(n):
        con[i][1] = G[i]
        
    return con,bx
    
def solitonupslope(x,stage,bot1,bot2,bxl,bxr,a1,g,t0,sxl,sxr):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
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
            u[i] = ue
            
    
    G = getGfromu(con,u,bx,0.0,0.0,stage - bot1,stage - bot2,bot1,bot2,dx)
    for i in range(n):
        con[i][1] = G[i]
        
    return con,bx
    
def jumpupslope(x,stage,bot1,bot2,bxl,bxr,ha,g,t0,sxl,sxr):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
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
    
    G = getGfromu(con,u,bx,0.0,0.0,stage - bot1,stage - bot2,bot1,bot2,dx)
    for i in range(n):
        con[i][1] = G[i]
        
    return con,bx
    
def flatbump(x,stage,bot,hadd,uadd,xr,xl):
    n = len(x)
    con = zeros((n,2))
    bx = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        bx[i] = bot            
        con[i][0] = stage - bx[i]
        
        if(x[i] >= xr and x[i] <= xl):
            con[i][0] = con[i][0] +  hadd
            u[i] = u[i] + uadd

    G = getGfromu(con,u,bx,0.0,0.0,stage-bot,stage-bot,bot,bot,dx)
    for i in range(n):
        con[i][1] = G[i]
        
    return con, bx
    
def powerfunction(r,n):
    if ( r >= 0 and r<= 1):
        return (1-r)**n
    else:
        return 0
    
    
    
def flowoverbumpChris(x,stage,center,width,height,vel,l):
    n = len(x)
    con = zeros((n,2))
    u = zeros(n)
    bed = zeros(n)
    
    for i in range(n): 
        r = abs(x[i] - center) / width
        bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
        con[i][0] = stage - bed[i]
        u[i] = vel
       
    G = getGfromu(con,u,bed,vel,vel,con[0][0],con[-1][0],bed[0],bed[-1],dx)
    for i in range(n):
        con[i][1] = G[i]   
    
    
    return con,bed


############FLOW OVER BUMP
#set it up so its exact floating point

wdir = "../../../data/bumpChris/o2/"

stage = 1.0
center = 1000.0
width = 150
height = 0.5
el = 4
vel = 2

g = 9.81
dx = 0.1
Cr = 0.5
l = Cr / (2 + sqrt(g*stage) )
dt = l*dx
beta = 1.0
startx = 0.0
endx = 2000.0 + dx
startt = 0.0
endt = 2*dt  

if not os.path.exists(wdir):
    os.makedirs(wdir)

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 1.2

gapbig = int(0.5/dt)

con,bed = flowoverbumpChris(x,stage,center,width,height,vel,l)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
ui = vel
u0 = array([ui,ui,ui,ui])
u1 = array([ui,ui,ui,ui])
    
h0 = array([con[0][0],con[0][0],con[0][0],con[0][0]])
h1 = array([con[-1][0],con[-1][0],con[-1][0],con[-1][0]])

for i in range(1,len(t)):
      
    if(i % gapbig == 0 or i ==1):
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()

h = zeros(n)
G = zeros(n)
for i in range(n):
    h[i] = con[i][0]
    G[i] = con[i][1]







    
   
    

###### DAM BREAK #############
"""
#set it up so its exact floating point
dx = 0.1
l = 0.01
dt = l*dx
startx = 0.0
endx = 1000.0
startt = 0.0
endt = 30.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../data/dambreakf/"

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
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
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

wdir = "../data/bump/"

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
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""



##### Soliton over linear slope ########
"""
#set it up so its exact floating point
dx = 250.0
l = 0.004
dt = l*dx
startx = -250000
endx = 750000.0
startt = 0.0
endt = 8000.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../data/lsolslope2/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
stage = 1000.0
bot1 = 0.0
bot2 = 900.0
bxl = 300000.0
bxr = 450000.0
a1 = 2.0
sxl = 100000.0
sxr = 300000.0
gap = 1
gapbig = gap * 20



con, bed = solitonupslope(x,stage,bot1,bot2,bxl,bxr,a1,g,t[0],sxl,sxr)
#con, bed = jumpupslope(x,stage,bot1,bot2,bxl,bxr,a1,g,t[0],sxl,sxr)


x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
u0 = array([0.0,0.0,0.0,0.0])
u1 = array([0.0,0.0,0.0,0.0])
    
h0 = stage - b0
h1 = stage - b1

h = []
G = []
for j in range(n):
    h = append(h,con[j][0])
    G = append(G,con[j][1])
u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        
plot(x,h+bed,'b', label="1")
plot(x,bed,'g', label="2")
xlim([startx,endx])
ylim([-0.1,1010.0])
title("Flow over Bump")
xlabel("Distance (m)")
ylabel("Water Height (m)")
    
s = wdir + "height"  + ".png"
            
savefig(s, bbox_inches='tight')        
clf()

plot(x,h+bed,'b', label="1")
plot(x,bed,'g', label="2")
xlim([startx,endx])
ylim([990,1010])
title("Flow over Bump")
xlabel("Distance (m)")
ylabel("Water Height (m)")
    
s = wdir + "heightz"  + ".png"
            
savefig(s, bbox_inches='tight')        
clf()

            
plot(x,u,'r', label="1")
xlim([startx,endx])
ylim([-1.0,1.0])
title("Flow over Bump")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
    
s = wdir + "velocity" + ".png"
    
savefig(s, bbox_inches='tight')
clf()

for i in range(1,len(t)):
      
    if(i % gapbig == 0 or i ==1):
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
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

wdir = "../data/deep/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
stage = 1.0
bot = 0.0
hadd = 0.01
uadd = 0.0
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
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""

##### Deep Water Soliton ########
"""
#set it up so its exact floating point
dx = 0.1
l = 0.006
dt = l*dx
startx = -200.0
endx = 500.0
startt = 0.0
endt = 30.0 + dt

szoomx = startx
ezoomx = endx

wdir = "../data/deep/"

g = 10.0

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 2.0
stage = 100.0
a0 = 10.0
a1 = 1.0
bot = 0.0
gap = 20
gapbig = gap * 5

con,bed = solitoninit1(n,stage,a0,a1,g,x,t[0],bot,dx)

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
        u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()
"""


### SMOOTH FLAT LAKE ####################################
"""
from numpy.linalg import norm

wdir = "../data/smfl/"
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
    G = []
    for j in range(n):
        h = append(h,con[j][0])
        G = append(G,con[j][1])
    
    u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
            
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

wdir = "../data/dfl/"
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
    G = []
    for j in range(n):
        h = append(h,con[j][0])
        G = append(G,con[j][1])
    
    u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
            
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










### SOLITON ACCURACY TEST! ##############################
"""
from numpy.linalg import norm
hs = []
us = []
hts = []
uts = []
xs = []
dxs = []
beds = [] 
normhdiff = []
normudiff = [] 
 
wdir = "../data/soliton/"
beta = 2.0

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity'])

nxs = [10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
#do accuracy test
for k in range(2,len(nxs)):
    dx = nxs[k]
    l = 0.02
    dt = l*dx
    startx = -500.0
    endx = 1500.0 + dx
    startt = 0
    endt = 100.0 + dt
    
    szoomx = startx
    ezoomx = endx
        
    g = 10.0
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    ui = 0.0
    stage = 10.0
    a0 = 10.0
    a1 = 1.0
    gap = 20
    gapbig = gap * 25
        
    x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
    x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

    b0 = 0.0*x0
    b1 = 0.0*x1
    u0 = array([ui,ui,ui,ui])
    u1 = array([ui,ui,ui,ui])
    
    h0 = array([stage,stage,stage,stage]) - b0
    h1 = array([stage,stage,stage,stage]) - b1
    
    con,bed = solitoninit(n,a0,a1,g,x,t[0],0,dx)

    for i in range(1,len(t)): 
                
        con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
        print t[i]
        print con[5]

                
    h = []
    G = []
    for j in range(n):
        h = append(h,con[j][0])
        G = append(G,con[j][1])
    
    u = getufromG(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
            
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],(t[-1]),g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)
    
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

