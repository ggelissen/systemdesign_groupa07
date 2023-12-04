# -*- coding: utf-8 -*-
"""
Loading diagrams due to wing weight 
"""
import scipy as sp
import numpy as np 
import matplotlib.pyplot as plt

## assumed coordinate system: x forward (pointing to nose), y to the right (pointing to wingtip), z downwards
## this coord sys. implies weight is pos, as it points to earth, so shear is pos & moment is neg

## three components: dry Ww, fuel, engine 
b = 66.9 #m
Ww = 38229.5/2 #kg
Wf = (125407 + 522.9)/2 #kg
Weng = 6033 #kg 
grav = 9.81 #m/s^2

y_vals = np.arange(0, b/2, 0.5)

def closest(lst, val): #finding entry in list closest to val
    lst = np.asarray(lst)
    idx = (np.abs(lst-val)).argmin()
    
    return lst[idx]

""" loading distribution """
## functions to define all loading distribution (ie decreasing triangular shape for dry, const for fuel)

def drydistr(y, b, WW):
    grav = 9.81 #m/s^2
    
    if y == 0:
        f = 0
    if y != 0 :
        c = 2*Ww*grav/0.5*b
        a = -b/0.5*b 
        f = a*y + c 
    return f 

def fueldistr(y, b, Wf):
    grav = 9.81 #m/s^2
    if y < b/4 and y > 0: 
        g = Wf*grav/(b/4)
    if y > b/4 or y == 0:
        g = 0 
    return g
    
fuelload = []
for element in y_vals: 
    if element == 0:
        fuelload.append(0)
    if element != 0 :
        fuelload.append(fueldistr(element, b, Wf))
    
dryload = []
for element in y_vals:
    dryload.append(drydistr(element,b,Ww))
    
engload = []
for element in y_vals:
    if element != closest(y_vals, (b/2)*0.35):
        engload.append(0)
    if element == closest(y_vals, (b/2)*0.35):
        engload.append(Weng*grav)
    
load = np.array(dryload) + np.array(fuelload) + np.array(engload) ## ALL LOADING DUE TO WEIGHT
"""function loading to be integrated """
def cts_loaddistr(y):    
    if y == 0:
        f = g = 0  
    if y != 0 :
        c = 2*Ww*grav/0.5*b
        a = -b/0.5*b 
        f = a*y + c     
        
    if y <= b/4 and y > 0: 
        g = Wf*grav/(b/4)
    if y > b/4:
        g = 0 
    return f + g ##f is structural weight, g is fuel weight

""" shear and moment calculations"""
"""cts_sheardistr = [] ## integrating the actual functions, later adding the shear due to point load (eng)
def sheardistribution(y):
    shear, error = sp.integrate.quad(cts_loaddistr, y, b/2)
    return shear 

for element in y_vals:
    if element == 0:
        cts_sheardistr.append(0)
    if element > 0:
        cts_sheardistr.append(sheardistribution(element))

engshear = []
## assuming the shear is a constants function due to a point force
for element in y_vals: 
    if element <= closest(y_vals, (b/2)*0.35) and element > 0: 
        engshear.append(Weng*grav)
    if element > closest(y_vals, (b/2)*0.35) or element == 0:
        engshear.append(0)
    
shear = np.array(cts_sheardistr) + np.array(engshear) ##ALL SHEAR CAUSED BY WEIGHT
        

cts_momentdistr = [] ## same reasoning as for the shear 
def momentdistribution(y):
    b = 66.9 #m 
    moment, error = sp.integrate.quad(sheardistribution, y, b/2)
    return -moment # since the moment is the negative integral of shear

for element in y_vals:
    if element == 0:
        cts_momentdistr.append(0)
    if element > 0:
        cts_momentdistr.append(momentdistribution(element))
    
engmoment = []
## assuming Msh = Ma - P*g (where p is the load, in this case Weng*grav, and Ma is p*L with L being the moment arm)
for element in y_vals:
    if element <= closest(y_vals, (b/2)*0.35) and element > 0:
        engmoment.append(Weng*grav*(element - (0.35*0.5*b)))
    if element > closest(y_vals, (b/2)*0.35) or element == 0:
        engmoment.append(0)


moment = np.array(cts_momentdistr) + np.array(engmoment) """ ##ALL MOMENTS CAUSED BY WEIGHT

""" plots """

## loading plots (force, weight)
plt.subplot(221)
plt.plot(y_vals[1:], load[1:])
plt.xlabel('Spanwise location [m]')
plt.ylabel('Weight [N]')
plt.title("total load")

plt.subplot(222)
plt.plot(y_vals[1:], dryload[1:])
plt.xlabel('Spanwise location [m]')
plt.ylabel('Weight [N]')
plt.title("Structural weight")


plt.subplot(223)
plt.plot(y_vals, fuelload)
plt.axis([0, 33, 0, 40000]) 
plt.xlabel('Spanwise location [m]')
plt.ylabel('Weight [N]')
plt.title("Fuel weight")

plt.subplot(224)
plt.plot(y_vals, engload)
plt.xlabel('Spanwise location [m]')
plt.ylabel('Weight [N]')
plt.title("Engine weight")

plt.subplots_adjust(wspace=0.6, hspace=0.7)
"""

## shear and moment plots 
plt.subplot(121)
plt.plot(y_vals,shear)
plt.xlabel('Spanwise location [m]')
plt.ylabel('Shear [N]')
plt.title("shear*EI")

plt.subplot(122)
plt.plot(y_vals, moment)
plt.xlabel('Spanwise location [m]')
plt.ylabel('Moment [Nm]')
plt.title("moment*EI") 

plt.subplots_adjust(wspace=0.45)"""




plt.show()

    