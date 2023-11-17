# -*- coding: utf-8 -*-
"""
Loading diagrams due to wing weight 
"""

import scipy as sp
import numpy as np 
import matplotlib.pyplot as plt

### three components: dry Ww, fuel, engine 
b = 66.9 #m
Ww = 38229.5/2 #kg
Wf = (125407 + 522.9)/2 #kg
Weng = 6033 #kg 
grav = 9.81 #m/s^2

x_vals = np.arange(0, b/2, 0.5)

def closest(lst, val): #finding item in list closest to val
    lst = np.asarray(lst)
    idx = (np.abs(lst-val)).argmin()
    
    return lst[idx]

def drydistr(x):
    b = 66.9 #m
    Ww = 38229.5/2 #kg
    grav = 9.81 #m/s^2
    
    c = 2*Ww*grav/0.5*b
    a = -b/0.5*b 
    f = a*x + c 
    return f 

def fueldistr(x):
    b = 66.9 #m
    Wf = (125407 + 522.9)/2 #kg
    grav = 9.81 #m/s^2
    g = 0 
    
    if x < b/4: 
        g = Wf*grav/(b/4)
    if x > b/4:
        g = 0 
    return g

##def engdistr(x): #figure out later
   # b = 66.9 #m
   # Weng = 6033 #kg 
   # grav = 9.81 #m/s^2
    
def loaddistr(x):
    f = drydistr(x) + fueldistr(x)
    return f

sheardistr = [] 
momentdistr = []

def sheardistribution(x):
    b = 66.9 #m
    shear, error = sp.integrate.quad(loaddistr, x, b/2)
    return shear 

def momentdistribution(x):
    b = 66.9 #m 
    moment, error = sp.integrate.quad(sheardistribution, x, b/2)
    return moment 

for element in x_vals:
    y = sheardistribution(element)
    sheardistr.append(y)

for element in x_vals:
    momentdistr.append(momentdistribution(element))
    
""" checking loading distributions """   
fuelload = []
for element in x_vals: 
    fuelload.append(fueldistr(element))
    
dryload = []
for element in x_vals:
    dryload.append(drydistr(element))
    
load = []
for element in x_vals:
    load.append(loaddistr(element))

#plt.plot(x_vals, dryload)
#plt.plot(x_vals, fuelload)
plt.plot(x_vals, load, label = "load")

""" shear & moment plots """
#plt.plot(x_vals,sheardistr, label = "shear")
#plt.plot(x_vals, momentdistr, label = "moment")
plt.legend()
plt.show()
    
## test test test 
    