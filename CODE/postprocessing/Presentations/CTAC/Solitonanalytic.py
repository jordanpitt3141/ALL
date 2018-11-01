from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os   


def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
 
 
def solitoninit(a0,a1,g,x,t0):
    n = len(x)
    h = zeros(n)
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)        
    return h

a0 = 1
a1 = 0.7
g = 9.81    

sx = -50
ex = 250
dx = 0.1

st = 0
dt = 0.1
et = 50

wdirb = "/home/jp/Documents/PhD/project/master/FigureData/Presentation/CTAC/Soliton/Analytic/"

wdirD = wdirb + "Data/"
if not os.path.exists(wdirD):
        os.makedirs(wdirD)
        
wdirT = wdirb + "Texs/"
if not os.path.exists(wdirT):
        os.makedirs(wdirT)


#make Data files

x = arange(sx,ex+dx,dx)
t = arange(st,et+dt,dt)
n = len(x)
nt = len(t)

for j in range(nt):
  
    ts = t[j]
    h = solitoninit(a0,a1,g,x,ts)  
    
    timesn = "{:02.1f}".format(ts)

    fileDn = wdirD + "ht="+ timesn + ".dat"
    with open(fileDn,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i])
            file1.write(s) 
            
    fileTn = wdirT + "h" + str(j) +  ".tex"
    #Make A test file
    with open(fileTn,'w') as file1:
    
        sheader= " \\documentclass[]{standalone} \n \\usepackage{pgfplots} \n \\usepgfplotslibrary{fillbetween} \n \\usepackage{tikz} \n" \
        + " \\usepackage{amsmath} \n \\usepackage{pgfplots}  \n \\usepackage{sansmath} \n \\sansmath \n \\usetikzlibrary{calc} " \
        + "\\pgfplotsset{compat = newest, every axis plot post/.style={line join=round}, label style={font=\\Huge},every tick label/.append style={font=\\Huge} }\n"
        file1.write(sheader) 
        
        sdocs  = "\\begin{document}  \n	\\begin{tikzpicture}\n"   
        file1.write(sdocs) 
        
        saxisoptions = "\\fontfamily{cmss} \n	\\begin{axis}[ \n width=40cm, \n height = 20cm, \n	xtick={-50,0,50,100,150,200,250}, \n	ytick = {0,0.5,1,1.5,2}, \n" \
    	+ "xmin=-50, \n xmax=250, \n ymin =-0.1, \n ymax = 2.0,\n  xlabel=$x(m)$, \n	ylabel=$z(m)$ ]\n"
        file1.write(saxisoptions) 
        
        snode = "\\node[label={\Huge$t="+timesn+"s$}] at (axis cs:230,1.8) {}; \n"
        file1.write(snode) 
        
        splot = "	\\addplot [name path=s,blue,no markers,ultra thick] table {../Data/ht="+timesn+".dat }; \n \\draw[brown!80!black,thick,name path=b,ultra thick] (-50, 0  ) -- (250,0); \n " \
    	+ "\\path[name path=a] (axis cs:-50,-0.1) -- (axis cs:250,-0.1);\n"
        file1.write(splot) 
        
        
        sfillplots = "\\addplot [thick,color=brown!60!black,fill=brown!60!black, fill opacity=0.3] fill between[of=b and a]; \n " \
    	+ "\\addplot [thick,color=blue,	fill=blue,fill opacity=0.3] fill between[of=s and b]; \n"
     
        file1.write(sfillplots) 
     
        send = "\\end{axis} \n \\end{tikzpicture} \n \\end{document}"
        file1.write(send) 

