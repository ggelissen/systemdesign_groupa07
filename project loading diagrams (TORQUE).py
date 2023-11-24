# -*- coding: utf-8 -*-
"""
Torque diagrams due to engine weight & thrust
"""

import numpy as np 
import matplotlib.pyplot as plt

b = 66.9 #m
Weng = 6630 #kg
grav  = 9.81 #m/s^2 
T = 346944 #N
h = 2.127 #heigh of eng to center wing
sweep05 = 31.49 #deg, sweep at half chord 
c035 = 10.0232 #m
engcenter = 1.43 #wrt leading edge
taper = 0.28
d2r = np.pi/180

sweep_LE = 37.12 #deg
cr = 13.4 #m 

y_vals = np.arange(0, b/2, 0.5)

""" Torque values """
#sweepengcen = np.arctan(np.tan(sweep_LE*d2r)-(engcenter/c035)*(2*cr*(1-taper)/b))/d2r #sweep at engine center
Tw = Weng*grav*((c035/2)-engcenter)

Tt = np.cos(sweep05*d2r)*T*h

def closest(lst, val): #finding item in list closest to val
    lst = np.asarray(lst)
    idx = (np.abs(lst-val)).argmin()
    
    return lst[idx]

torque = []
for element in y_vals :
    if element <= closest(y_vals, (b/2)*0.35) and element >= 0:
        torque.append(Tw+Tt)
    if element > closest(y_vals, (b/2)*0.35): #or element == 0:
        torque.append(0)

plt.plot(y_vals, torque)
plt.show()

#------------------------------------------------------------------------------------------------------------------
"""
twist diagram
"""
import scipy as sp
t_1=float(input('Enter t1: '))
t_11=float(input('Enter t11: '))
t_2=float(input('Enter t2: '))
t_3=float(input('Enter t3: '))
t_4=float(input('Enter t4: '))
L_4=float(input('Enter L4: '))
G=26000000000
integrands=[]
def chord(y):
   c_root = 13.4
   c_tip = 3.8
   half_span = 33.45
   m = (c_tip - c_root) / half_span
   return c_root + (m * y)

import Torsional_stiffness
for i in range(len(y_vals)):
   h_length = float(0.5*chord(y_vals[i]))
   l_up=float(0.5004620365222*chord(y_vals[i]))
   l_low=float(0.5003958533001*chord(y_vals[i]))
   L_1 = float(0.1082 * chord(y_vals[i]))
   L_2= float(0.0668 * chord(y_vals[i]))
   L_3 = float(0.35 * chord(y_vals[i]))
   c=float(y_vals[i])
   T=float(torque[i])
   if y_vals[i] <= L_4:
      integrands.append(Torsional_stiffness.torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,2,T))
   else:
      integrands.append(Torsional_stiffness.torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,1,T))

integral_values = sp.integrate.cumtrapz(integrands, x=y_vals, initial=0)
integral_values = integral_values*180/math.pi
plt.plot(y_vals, integral_values)
plt.show()
