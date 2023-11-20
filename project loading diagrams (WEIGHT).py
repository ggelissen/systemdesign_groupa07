# -*- coding: utf-8 -*-
"""
Loading diagrams due to wing weight 
"""
import scipy as sp
import numpy as np 
import matplotlib.pyplot as plt

### assumed coordinate system: x to the right, y downwards 

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

""" coordinate system """


""" loading distribution """

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
    
def cts_loaddistr(x):
    f = drydistr(x) + fueldistr(x)
    return f
    
 
fuelload = []
for element in x_vals: 
    fuelload.append(fueldistr(element))
    
dryload = []
for element in x_vals:
    dryload.append(drydistr(element))
    
engload = []
for element in x_vals:
    if element != closest(x_vals, (b/2)*0.35):
        engload.append(0)
    if element == closest(x_vals, (b/2)*0.35):
        engload.append(Weng*grav)
    
cts_load = []
for element in x_vals:
    cts_load.append(cts_loaddistr(element))
    
load = np.array(cts_load) + np.array(engload) ## ALL LOADING DUE TO WEIGHT

""" shear and moment calculations"""
cts_sheardistr = [] 
def sheardistribution(x):
    b = 66.9 #m
    shear, error = sp.integrate.quad(cts_loaddistr, x, b/2)
    return shear 

for element in x_vals:
    cts_sheardistr.append(sheardistribution(element))

engshear = []
for element in x_vals: 
    if element <= closest(x_vals, (b/2)*0.35):
        engshear.append(Weng*grav)
   # if element == closest(x_vals, (b/2)*0.35):
    #    engshear.append(-Weng*grav)
    if element > closest(x_vals, (b/2)*0.35):
        engshear.append(0)
    
#shear = np.array(cts_sheardistr).astype(float) + np.array(engshear).astype(float)
shear = np.array(cts_sheardistr) + np.array(engshear)
        

cts_momentdistr = []
def momentdistribution(x):
    b = 66.9 #m 
    moment, error = sp.integrate.quad(sheardistribution, x, b/2)
    return -moment ## since the moment is the negative integral of shear

for element in x_vals:
    cts_momentdistr.append(momentdistribution(element))
    
engmoment = []
### assuming Msh = Ma - P*g (where p is the load, in this case Weng*grav, and Ma is p*L with L being the moment arm)
for element in x_vals:
    if element <= closest(x_vals, (b/2)*0.35):
        engmoment.append(Weng*grav*(element - (0.35*0.5*b)))
    if element > closest(x_vals, (b/2)*0.35):
        engmoment.append(0)


moment = np.array(cts_momentdistr) + np.array(engmoment)

""" plots """
#plt.plot(x_vals, load, label = "load")
#plt.plot(x_vals, engshear)

plt.plot(x_vals,shear, label = "shear*EI")
plt.plot(x_vals, moment, label = "moment*EI")
plt.legend()
plt.show()
    