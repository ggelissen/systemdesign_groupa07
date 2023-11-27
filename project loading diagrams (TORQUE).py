# -*- coding: utf-8 -*-
"""
Torque diagrams due to engine weight & thrust
"""

import numpy as np 
import matplotlib.pyplot as plt
from LiftWeight_ShearBendingDiagrams import Cmc4lst
from LiftWeight_ShearBendingDiagrams import ylst 

## assumed coordinate system: x forward (pointing to nose), y to the right (pointing to wingtip), z downwards

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

rho = 1.225 #kg/m^3
V = 242 #m/s
S = 574.3 #m^2

sweep_LE = 37.12 #deg
cr = 13.4 #m 

y_vals = np.arange(0, b/2, 0.5)

""" Torque values """
#sweepengcen = np.arctan(np.tan(sweep_LE*d2r)-(engcenter/c035)*(2*cr*(1-taper)/b))/d2r #sweep at engine center
Tw = Weng*grav*((c035/2)-engcenter)

Tt = np.cos(sweep05*d2r)*T*h

def closest(lst, val): #finding entry in list closest to val
    lst = np.asarray(lst)
    idx = (np.abs(lst-val)).argmin()
    
    return lst[idx]

torque = []
for element in ylst:
    if element <= closest(ylst, (b/2)*0.35) and element > 0:
        torque.append(-Tw+Tt)
    if element > closest(ylst, (b/2)*0.35) or element == 0:
        torque.append(0)
        
moment = []
def chord(y):
    return cr - cr*(1-taper)*(y/(b/2))
    
for i in range(len(ylst)):
    moment.append(Cmc4lst[i]*0.5*rho*(chord(ylst[i])**2)*V**2) ## m = (1/2)Cm*rho*c^2*V^2

print(len(torque))
total = np.array(moment) + np.array(torque) ## questioning if the 

plt.plot(ylst, total)
plt.show()