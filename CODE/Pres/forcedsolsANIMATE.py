# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 09:50:15 2018

@author: jp
"""

from scipy import *
import csv
import os

def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b
  
wdir = "/home/jp/Documents/PhD/project/data/2018/raw/Presentation/ForcedDryBed/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

dx = 0.1
  
a6= 1.0
a7 = 2*pi/50.0

width = 2*(2*pi/a7)
    
a0 = 0.0
a1 = 0.5
a2 =  ((2*pi) / a7)/10.0
a3 = -pi/2.0/a7 -width/4.0
a4 = width/2**6
a5 = a1


g = 9.81

startx = -pi/2.0/a7 -width
sx= startx
endx = -pi/2.0/a7 +width
ex = endx
startt = 0.0
endt = (2*pi/a7) / a2
et = endt


#ts = [0,2.5,5.0,7.5,10.0]
x = arange(startx,endx +0.1*dx, dx)
n = len(x)
    
h,u,G,b = ForcedbedM(x,0,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)


import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-50, 25), ylim=(-1, 1.5))
ax.plot(x, b,'-',color='tan')
ax.fill_between(x, -1, b,color='tan')
plt.xlabel('$x$ ($m$)')
plt.ylabel('$z$ ($m$)')
plt.xticks([-50,-25,0,25])
plt.yticks([-1,-0.5,0,0.5,1.0,1.5])
line, = ax.plot([], [],'-b', lw=2)
time_text = ax.text(0.8, 0.05, '', transform=ax.transAxes)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line,time_text

# animation function.  This is called sequentially
def animate(i):
    h,u,G,b = ForcedbedM(x,i*0.02,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    time_text.set_text('Time = %.1f' % (i*0.02))
    line.set_data(x, h + b)
    return line,time_text

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=500, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('Final.mp4', fps=60, extra_args=['-vcodec', 'libx264'])

plt.show()