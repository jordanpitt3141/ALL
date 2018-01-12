# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
import os
from subprocess import call

diffw = "11"
wdirord = "o3"

#wdirords = ["o3","FDcent","grim","o2","o1"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]



wdirords = ["o3","FDcent","grim","o2","o1"]
ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]


#wdirords = ["o3","FDcent","o2","o1"]
#ordsubtup = [[6,7],[5,6], [6,8], [6,7]]

diffws = [12]
dxs = [11]
#diffws = [6]

for k in dxs:
    for jp in diffws:
        
        diffw = str(jp)

        ylims = [[0.0,2],[1.0,1.8],[1.3,1.45],[1.3,1.45],[1.3,1.45]]
        xlims = [[0,1000],[299,701],[499,561],[524,566],[527,537]]
        #gaps = [[2,4,6,8,10,12,14],[1,1,2,6,10,14,14],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]]
        #zoom = [True,True,False,False,False]
        #zoomints = [[500,560],[500,560],[500,560],[500,560],[530,540]]
        #zoomgap = [2,1,1,1,1] 
        
        gaps = [[30,30,30,30,30], \
                [20,20,20,20,20], \
                [10,10,10,10,10], \
                [2,2,2,2,2], \
                [1,1,1,1,1]]
        zoom = [False,True,False,False,False]
        zoomints = [[350,650],[525,575],[350,650],[350,650],[350,650]]
        zoomgap = [[8,8,8,8,8], [5,5,5,5,5] , [2,2,2,2,2], [2,2,2,2,2],[2,2,2,2,2],[2,2,2,2,2]]  
        
        
        for l in range(len(ylims)):
        #for l in range(1):
            
            cylim = ylims[l]
            cxlim = xlims[l]
            
            sdir = "../../../../../../data/postprocessing/smoothdball/AMRd/"+str(diffw)+ "/"
            if not os.path.exists(sdir):
                os.makedirs(sdir)
            
            sdirf = sdir +str(l) + "/"
            if not os.path.exists(sdirf):
                os.makedirs(sdirf)
            
            for ip in range(len(wdirords)):
                wdirord = wdirords[ip]
                gap = gaps[l][ip]
                
                wdir = "../../../../../../data/raw/Joebigsmooth/"  +wdirord +"/" + str(k) + "/" + diffw + "/"
                
            
                     
                s = wdir + "outlast.txt"
                print(s)
                with open(s,'r') as file1:
                     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     h = []
                     u = []
                     x = []
                     j = -1
                     for row in readfile:       
                         if (j >= 0):
                            dx = float(row[0])
                            dt = float(row[1])
                            t = float(row[2])
                            x.append(float(row[3]))
                            h.append(float(row[4]))
                            u.append(float(row[ordsubtup[ip][0]])) #5 for FDcent, 6 otherwise 
                            #beta = float(row[ordsubtup[ip][1]]) #could be 8 as well for o2
                         j = j + 1
                     igap = int(gap)
                     if zoom[l]:
                         xbeg = int(cxlim[0]/dx)
                         xend = int(cxlim[1]/dx)
                         zbeg = int(zoomints[l][0]/dx)
                         zend = int(zoomints[l][1]/dx)
                         print(l,k,zend,xend,igap)
                         #zoomgapi = int(zoomgap[l][k-numbeg])
                         zoomgapi = int(zoomgap[l][ip])
                         xt = x[xbeg:zbeg:igap] + x[zbeg:zend:zoomgapi] + x[zend:xend:igap]
                         ht = h[xbeg:zbeg:igap] + h[zbeg:zend:zoomgapi] + h[zend:xend:igap]
                     else:
                         xbeg = int(cxlim[0]/dx)
                         xend = int(cxlim[1]/dx)
                         xt = x[xbeg:xend:igap]
                         ht = h[xbeg:xend:igap]
                     x = array(xt)
                     h = array(ht)     
               
                beta = jp               
                m = len(x)
                ap = 1.736397786
                #h2 = 1.36898
                const = ap*ones(m)
                s = str(dx) 
                plot(x,h ,label=s)
                plot(x,const,"--k" ,label="a ref")
                
                ylim(cylim)
                xlim(cxlim)
                #eyticks = [h2]
                #yticks(list(yticks()[0]) + eyticks)               
                xlabel("$x$ ($m$)")
                ylabel("$h$ ($m$)")
                
                n = len(xt)
                s = sdirf  + "mod" + str(wdirord) + "h.dat"
                with open(s,'w') as file1:
                    for i in range(n):
                        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
                        file1.write(s)
            s = "Dam Break: " + wdirord + " diff = " + str(beta)
            title(s)
            s = sdir +str(l)+".png"       
            savefig(s, bbox_inches='tight')
            legend()
            s = sdir +str(l)+ "leg.png"       
            savefig(s, bbox_inches='tight')        
            clf()
                #legend()
