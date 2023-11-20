# -*- coding: utf-8 -*-
"""
Loading diagrams due to wing weight 
"""
import scipy as sp
import numpy as np 
import matplotlib.pyplot as plt

## assumed coordinate system: x forward (pointing to nose), y to the right (pointing to wingtip), z downwards

## three components: dry Ww, fuel, engine 
b = 67 #m
Ww = 38229.5/2 #kg
Wf = (125407 + 522.9)/2 #kg
Weng = 6033 + 3127.482/2 #kg 
grav = 9.81 #m/s^2

y_vals = np.arange(0, b/2, 0.5)

def closest(lst, val): #finding item in list closest to val
    lst = np.asarray(lst)
    idx = (np.abs(lst-val)).argmin()
    
    return lst[idx]

""" coordinate system """
## finish? 

""" loading distribution """
## functions to define all loading distribution (ie decreasing triangular shape for dry, const for fuel)
def drydistr(y):
    b = 66.9 #m
    Ww = 38229.5/2 #kg
    grav = 9.81 #m/s^2
    
    c = 2*Ww*grav/0.5*b
    a = -b/0.5*b 
    f = a*y + c 
    return f 

def fueldistr(y):
    b = 66.9 #m
    Wf = (125407 + 522.9)/2 #kg
    grav = 9.81 #m/s^2
    g = 0 
    
    if y < b/4: 
        g = Wf*grav/(b/4)
    if y > b/4:
        g = 0 
    return g
    
def cts_loaddistr(y): #not really necessary, could be written differently/more efficient
    f = drydistr(y) + fueldistr(y) 
    return f
    
fuelload = []
for element in y_vals: 
    fuelload.append(fueldistr(element))
    
dryload = []
for element in y_vals:
    dryload.append(drydistr(element))
    
engload = []
for element in y_vals:
    if element != closest(y_vals, (b/2)*0.35):
        engload.append(0)
    if element == closest(y_vals, (b/2)*0.35):
        engload.append(Weng*grav)
    
cts_load = []
for element in y_vals:
    cts_load.append(cts_loaddistr(element))
    
load = np.array(cts_load) + np.array(engload) ## ALL LOADING DUE TO WEIGHT

""" shear and moment calculations"""
cts_sheardistr = [] ## integrating the actual functions, later adding the shear due to point load (eng)
def sheardistribution(y):
    b = 66.9 #m
    shear, error = sp.integrate.quad(cts_loaddistr, y, b/2)
    return shear 

for element in y_vals:
    cts_sheardistr.append(sheardistribution(element))

engshear = []
## assuming the shear is a constants function due to a point force
for element in y_vals: 
    if element <= closest(y_vals, (b/2)*0.35): 
        engshear.append(Weng*grav)
    if element > closest(y_vals, (b/2)*0.35):
        engshear.append(0)
    
shear = np.array(cts_sheardistr) + np.array(engshear) ##ALL SHEAR CAUSED BY WEIGHT
        

cts_momentdistr = [] ## same reasoning as for the shear 
def momentdistribution(y):
    b = 66.9 #m 
    moment, error = sp.integrate.quad(sheardistribution, y, b/2)
    return -moment # since the moment is the negative integral of shear

for element in y_vals:
    cts_momentdistr.append(momentdistribution(element))
    
engmoment = []
## assuming Msh = Ma - P*g (where p is the load, in this case Weng*grav, and Ma is p*L with L being the moment arm)
for element in y_vals:
    if element <= closest(y_vals, (b/2)*0.35):
        engmoment.append(Weng*grav*(element - (0.35*0.5*b)))
    if element > closest(y_vals, (b/2)*0.35):
        engmoment.append(0)


moment = np.array(cts_momentdistr) + np.array(engmoment) ##ALL MOMENTS CAUSED BY WEIGHT

""" plots """
#plt.plot(x_vals, load, label = "load")
#plt.plot(x_vals, engshear)

plt.plot(y_vals,shear, label = "shear*EI")
plt.plot(y_vals, moment, label = "moment*EI")
plt.legend()
plt.show()
    