from csv import reader
import scipy as sp
from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np

rho = 0.324438
v = 258.9704
q = 0.5*rho*(v**2)
halfspan = 33.5
n = 2.5
b = 67  # m
Ww = 38229.5 / 2  # kg
Wf = (125407 + 522.9) / 2  # kg
Weng = 6033  # kg
grav = 9.81  # m/s^2
T = 0
h = 2.127 
engcenter = 1.43 
taper = 0.28
d2r = np.pi/180
S = 574.3 #m^2
sweep_LE = 37.12 #deg
cr = 13.4 #m 

CL0 = 0.04647
CL10 = 0.97586
CLD = 0.5785

# Create arrays for the values in the CSV file
ylst0 = np.array([])
chordlst0 = np.array([])
Ailst0 = np.array([])
Cllst0 = np.array([])
ICdlst0 = np.array([])
Cmc4lst0 = np.array([])

ylst10 = np.array([])
chordlst10 = np.array([])
Ailst10 = np.array([])
Cllst10 = np.array([])
ICdlst10 = np.array([])
Cmc4lst10 = np.array([])

# Read the CSV file
data0 = []
data10 = []
with open('A07csv0.csv', 'r') as file:
    csv_reader = reader(file)
    for row in csv_reader:
        data0.append(row)

with open('A07csv10.csv', 'r') as file2:
    csv_reader = reader(file2)
    for row in csv_reader:
        data10.append(row)

# Append correct values from csv_reader to arrays
for row in data0[51:81]:      # Range can be adjusted here!
    ylst0 = np.append(ylst0, float(row[0]))
    chordlst0 = np.append(chordlst0, float(row[1]))
    Ailst0 = np.append(Ailst0, float(row[2]))
    Cllst0 = np.append(Cllst0, float(row[3]))
    ICdlst0 = np.append(ICdlst0, float(row[5]))
    Cmc4lst0 = np.append(Cmc4lst0, float(row[7]))

for row in data10[51:81]:      # Range can be adjusted here!
    ylst10 = np.append(ylst10, float(row[0]))
    chordlst10 = np.append(chordlst10, float(row[1]))
    Ailst10 = np.append(Ailst10, float(row[2]))
    Cllst10 = np.append(Cllst10, float(row[3]))
    ICdlst10 = np.append(ICdlst10, float(row[5]))
    Cmc4lst10 = np.append(Cmc4lst10, float(row[7]))

# Functions to interpolate the values
def yCl(y, Cl):
    return sp.interpolate.interp1d(y,Cl,kind='cubic',fill_value="extrapolate")

def ychord(y, chord):
    return sp.interpolate.interp1d(y,chord,kind='cubic',fill_value="extrapolate")

def yICd(y, ICd):
    return sp.interpolate.interp1d(y,ICd,kind='cubic',fill_value="extrapolate")

def yCmc4(y, cmc4):
    return sp.interpolate.interp1d(y,cmc4,kind='cubic',fill_value="extrapolate")

# Define set of values for y
yvalues = np.arange(0, halfspan, 0.01)
yCl_result0 = yCl(ylst0, Cllst0)
ychord_result0 = ychord(ylst0, chordlst0)
yICd_result0 = yICd(ylst0, ICdlst0)
yCmc4_result0 = yCmc4(ylst0, Cmc4lst0)

yCl_result10 = yCl(ylst10, Cllst10)
ychord_result10 = ychord(ylst10, chordlst10)
yICd_result10 = yICd(ylst10, ICdlst10)
yCmc4_result10 = yCmc4(ylst10, Cmc4lst10)

# finding entry in list closest to val
def closest(lst, val):  
    lst = np.asarray(lst)
    idx = (np.abs(lst - val)).argmin()
    return lst[idx]

# functions to define all loading distribution (ie decreasing triangular shape for dry, const for fuel)
engload = []
for element in yvalues:
    if element != closest(yvalues, (b / 2) * 0.35):
        engload.append(0)
    if element == closest(yvalues, (b / 2) * 0.35):
        engload.append(Weng * grav)
def cts_loaddistr(y):
    f=0
    g=0
    #if y == 0:
    #   f = g = 0
    if y >= 0:
        c = 4 * Ww * grav / b
        a = (-1 * (Ww * grav * 8)) / (b ** 2)
        f = a * y + c
    if y <= b / 4 and y >= 0:
        h = 8 * Wf * (1 - 1/2.56) * grav / b
        d = 8 * Wf * grav / (2.56 * b)
        m = (d - h) / (b / 4)
        g = h + m * y
    if y > b / 4:
        g = 0
    return f + g  #f is structural weight, g is fuel weight

# Angle of Attack
alpha_sin = (CLD-CL0)/(CL10-CL0) * np.sin(10*np.pi/180)
alpha = np.arcsin(alpha_sin)

# Functions to calculate distributed load an point load
def Ldistribution0(x):
    return yCl_result0(x) * q * ychord_result0(x)
liftdistributionlst0 = np.array([])

def Ldistribution10(x):
    return yCl_result10(x) * q * ychord_result10(x)
liftdistributionlst10 = np.array([])

def LdistributionD(x):
    return (Ldistribution0(x) + ((CLD - CL0)/(CL10 - CL0)) * (Ldistribution10(x) - Ldistribution0(x))) * np.cos(alpha)
liftdistributionlst = np.array([])
for element in yvalues:
    liftdistributionlst = np.append(liftdistributionlst, (LdistributionD(element)*n) - cts_loaddistr(element))   
#plt.plot(yvalues, liftdistributionlst, "r")
#plt.xlabel('Spanwise location [m]')
#plt.ylabel('Deflection [m]')
#plt.title('Deflection Graph')
#plt.show()
liftdistributionlst = np.flip(liftdistributionlst)        
def sheardistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    sheardistributionlst = np.flip(shear)
    for i in range(len(yvalues)):
        if yvalues[i] <= (b / 2) * 0.35:
            sheardistributionlst[i] = sheardistributionlst[i] - Weng * grav * (b / 2) * 0.35
    for i in range(len(yvalues)):
        if yvalues[i] <= 5.8:
            sheardistributionlst[i] = sheardistributionlst[i] - 11383.7 / 2 * grav * 5.8
    return sheardistributionlst

def momentdistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    sheardistributionlst = np.flip(shear)
    for i in range(len(yvalues)):
        if yvalues[i] <= (b / 2) * 0.35:
            sheardistributionlst[i] = sheardistributionlst[i] - Weng * grav * (b / 2) * 0.35
    for i in range(len(yvalues)):
        if yvalues[i] <= 5.8:
            sheardistributionlst[i] = sheardistributionlst[i] - 11383.7 / 2 * grav * 5.8
    sheardistributionlst = np.flip(sheardistributionlst)
    moment = integrate.cumtrapz(sheardistributionlst, y, initial=0)
    momentdistributionlst = -1 * np.flip(moment)
    return momentdistributionlst


def chord(y):
    return cr - cr*(1-taper)*(y/(b/2))

Tw = Weng*grav*((chord(0.35*b/2)/2)+engcenter)
Tt = T*h*np.cos(sweep_LE*d2r)

torque = []
for element in yvalues:
    if element <= closest(yvalues, (b/2)*0.35) and element >= 0:
        torque.append(Tt-Tw)
    if element > closest(yvalues, (b/2)*0.35):
        torque.append(0)
        
       
moment = []
  
for i in range(len(yvalues)):
    moment.append(yCmc4_result0(yvalues[i])*0.5*rho*chord(yvalues[i])*S*v**2) ## M = (1/2)Cm*rho*c*S*V^2

total = np.array(moment) + np.array(torque)


sheardist = sheardistribution(yvalues)
sheardist[0] = 0
momentdist = momentdistribution(yvalues)

def momentone(y):
    i=np.where(yvalues == y)[0]
    return momentdist[i]
        

def y_moment(y, M):
    return sp.interpolate.interp1d(y, M, kind='cubic', fill_value="extrapolate")

momentfunction = y_moment(yvalues, momentdistribution(yvalues))

plt.subplot(1,3,1)
plt.plot(yvalues, sheardist, "b")
plt.xlabel('Spanwise location [m]')
plt.ylabel('Shear [N]')
plt.title('Shear Distribution')

plt.subplot(1,3,2)
plt.plot(yvalues,momentdist, "g")
plt.xlabel('Spanwise location [m]')
plt.ylabel('Moment [Nm]')
plt.title('Moment Distribution')


plt.subplot(1,3,3)
plt.plot(yvalues, total)
plt.xlabel('Spanwise location [m]')
plt.ylabel('Torque [Nm]')
plt.title('Torque Distribution')
plt.show()

##plt.plot(yvalues, deflection(yvalues), "r")
##plt.xlabel('Spanwise location [m]')
##plt.ylabel('Deflection [m]')
##plt.title('Deflection Graph')

plt.subplots_adjust(wspace=0.45)
#plt.show()
