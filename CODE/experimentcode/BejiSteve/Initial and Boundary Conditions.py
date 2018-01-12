# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from scipy import *
import csv
    
def makevar(sx,ex,dx,st,et,dt): 
    """ 
    makevar: creates my arrays of my spatial variable x and my temporal variable t
    
        Input:
            sx: starting location of spatial interval
            ex: end location of spatial interval
            dx : spacing of spatial grid
        
            st: starting time
            et: ending time
            dt : spacing of temporal grid
        Output:
        
           x : array of my cell centres, created from arange , will be an interval [sx, ex] with a spacing of dx
           t : array of my times at which I have my numerical solution, will be an interval [st,et] with a spacing of dt
    """
    
    # ex +0.9*dx because arange doesnt include end points, so we extend it to ensure, ex is in the interval at least
    x = arange(sx, ex +0.9*dx, dx)
    t = arange(st, et +0.9*dt, dt)
    
    return x,t 

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
    h = zeros(n)
    u = zeros(n)
    
    
    for i in range(n):
            
        #Define the bed profile with if statements as is piecewise.
        # h at x[i] is just defined so that h + bed is constant at 0.4m
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
            bed[i] = (0.2/5.0)*(x[i] - 18.95)
            h[i] = 0.4 - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            h[i] = 0.4 - bed[i]
        else:
            bed[i] = 0.0
            h[i] = 0.4  - bed[i]
            
    return h,u,bed

   
def lineinterp(y0,y1,x0,x1,x):
    """
    lineinterp: gives linear interpolation of f(x) between (x0,y0) and (x1,y1) at x
    
        Input:
            y0 : evaluation of function at x0, y0 = f(x0)
            y1 : evaluation of function at x1, y1 = f(x1)
            x0 : location of function evaluation
            x1 : location of function evaluation
            
        Output:
            Linear interpolated value
    """
    
    return y0  + (y1 - y0)/(x1 - x0)*(x - x0)
  


# \Delta x, as in the paper, want to divide 0.1m sections nicely
dx = (0.1/2.0**4)

#sr is the sample rate of the wave gauges, importantly giving the rate at which our B.Cs are defined
sr = 0.039312

#want our time steps to be nice dividers of the sample rate
dt = sr/ (2**5)

#starting location, will depend on B.Cs, importantly we want to define the incoming waves at 5.7m
#so in this set up we use another cell centered at 5.7m as a B.C
#importantly this means with our dx that all wvae gauges and the important features of the bed are all at cell centres and thus properly defined
startx = 5.7 + dx

#end of tank, just set large so over the time period we have no boundary effects
endx = 300

#start time and end time
startt = 0
endt = 60

#g acceleration due to gravity
g = 9.81

#theta, limiter
theta = 1.2

#get our spatial array x of cell centres, and our temporal array of evaluation times
x,t = makevar(startx,endx,dx,startt,endt,dt)

#get our water depth 'h', fluid velocity 'u' and bed term 'bed'
h,u,bed = BejiFlume(x)


# This is the location of the experimental files, which contains all the wave gauage data in csv format
expdir = "./94 Paper CSV/"

# This is the experiment we are doing it, this represents the high frequency sine wave
exp = "sh"

#Now we read the wave gauge data
s = expdir + exp + ".csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
  
    #  calculated_times is an array of the calculated times, from using the sample rate of wave gauages
    #   we will use calculated_times since it isnt rounded unlike read_times which was  
  
    # read_times is an array of the listed times from the file
    # wg4s represents the array containing the wave gauage data for wave gauage 4, converted to metres
    calculated_times = [0.0]
    read_times = [0.0]
    wg1s = [0.0]
    wg2s = [0.0]
    wg3s = [0.0]
    wg4s = [0.0]
    wg5s = [0.0]
    wg6s = [0.0]
    wg7s = [0.0]
    
    #j = -1, ensures we dont read the header, probably some smarter way to do it with csv.reader
    #in this case, I also use j to calculate the time
    j = -1
    for row in readfile:  
        
        if (j >= 0):
            calculated_times.append((j + 1)*sr)
            read_times.append(float(row[0]))
            
            #wave gauge data was in centimetres, so I convert to metres
            wg1s.append(float(row[1])/100.0)
            wg2s.append(float(row[2])/100.0)
            wg3s.append(float(row[3])/100.0)
            wg4s.append(float(row[4])/100.0)
            wg5s.append(float(row[5])/100.0)
            wg6s.append(float(row[6])/100.0)
            wg7s.append(float(row[7])/100.0)
        j = j + 1
        
# So we now have our initial conditions h,u and bed        
# We also have read in the wave gauage data
        
# I calculated the boundary conditions in the following way
# Since we need boundary conditions at all t and not all of them have values defined
# From wave gauage 1 , I use linear interpolation to get the wave heights at any time, here is an example

#In this example, I calculate the wave height given by linear interpolation of the wave gauage at current_time (which is 10s in this case)     
current_time = 10

#Now want to locate where current_time is in calculated_times, since calculated_times is just some n \times sr, then int(current_time / sr) gives the
# location of current_time in calculated_times, importantly if   time_lowerindex = int(current_time / sr)  then we know 
# calculated_times[time_lowerindex] less than or equal to time_lowerindex  and  time_lowerindex less than or equal to calculated_times[time_lowerindex + 1]
time_lowerindex = int(current_time/sr)

#now we calculate interpolated wave gauage at current_time in the interval  [calculated_times[time_lowerindex] , calculated_times[time_lowerindex + 1]]
interpolated_wave_height = lineinterp(wg1s[time_lowerindex],wg1s[time_lowerindex + 1], \
                                      calculated_times[time_lowerindex],calculated_times[time_lowerindex + 1],current_time)
                                      
#To calculate the water depth at 5.7, you need to add in the water depth 0.4m
interpolated_water_depth = 0.4 +  interpolated_wave_height