import random
from scipy import *
random.seed(1)

#Constants
i3 = 1.0/3.0
t2i3 = 2*i3



#Functions
def RandList(a,b,n):
    Out = zeros(n)
    
    for i in range(n):
        Out[i] = random.uniform(a, b) 
    
    avg = sum(Out)/ n
    
    Out =Out / avg
    return Out
    

a = 0.0
b = 1.0
n = 100

A = RandList(a,b,n)
B = RandList(a,b,n)

C= i3*A + t2i3*B  
D = (A + B)/2.0 

c = i3*sum(A)+ t2i3*sum(B)
c1 = (sum(A)+ 2*sum(B))/3.0
d = 0.5*sum(A)+ 0.5*sum(B)