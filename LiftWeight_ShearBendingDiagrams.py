from csv import reader
import scipy as sp
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np

# Read the CSV file
data = []
with open('A07csv3.csv', 'r') as file:
    csv_reader = reader(file)
    for row in csv_reader:
        data.append(row)

# Defined constants
rho = 1.225
v = 10
q = 0.5*rho*(v**2)
halfspan = 33.5
centroid = 14.4486

# Create arrays for the values in the CSV file
ylst = np.array([])
chordlst = np.array([]) 
Ailst = np.array([]) 
Cllst = np.array([]) 
ICdlst = np.array([]) 
Cmc4lst = np.array([]) 

# Append correct values from csv_reader to arrays
for row in data[51:81]:      # Range can be adjusted here!
    ylst = np.append(ylst, float(row[0]))
    chordlst = np.append(chordlst, float(row[1]))
    Ailst = np.append(Ailst, float(row[2]))
    Cllst = np.append(Cllst, float(row[3]))
    ICdlst = np.append(ICdlst, float(row[5]))
    Cmc4lst = np.append(Cmc4lst, float(row[7]))

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
yvalues = np.arange(0, halfspan, 0.5)

yCl_result = yCl(ylst, Cllst)
ychord_result = ychord(ylst, chordlst)
yICd_result = yICd(ylst, ICdlst)
yCmc4_result = yCmc4(ylst, Cmc4lst)

# Functions to calculate distributed load an point load
def Ldistribution(x):
    return yCl_result(x) * q * ychord_result(x)

def pointload():
    totallift, _ = sp.integrate.quad(Ldistribution, 0, halfspan, limit=1000)
    return totallift

def moment():
    return pointload() * centroid

# Functions to define shear and moment distributions
def sheardistribution(y):
    estimateshear , _ = sp.integrate.quad(Ldistribution, y, halfspan, limit=1000)
    return estimateshear - pointload() if y == 0 else estimateshear

sheardistributionlst = np.array([])
for element in yvalues:
    sheardistributionlst = np.append(sheardistributionlst, sheardistribution(element))

def momentdistribution(z):
    estimatemoment , _ = sp.integrate.quad(sheardistribution, z, halfspan, limit=1000)
    return -1*estimatemoment + moment() if z == 0 else -1*(estimatemoment)

momentdistributionlst = np.array([])
for element in yvalues:
    momentdistributionlst = np.append(momentdistributionlst, momentdistribution(element))

# Plot shear and moment distributions
plt.subplot(1,2,1)
plt.plot(yvalues, sheardistributionlst, 'b')
plt.xlabel('Spanwise location [m]')
plt.ylabel('Shear [N]')
plt.title('Shear distribution')

plt.subplot(1,2,2)
plt.plot(yvalues, momentdistributionlst, 'g')
plt.xlabel('Spanwise location [m]')
plt.ylabel('Moment [Nm]')
plt.title('Moment distribution')

plt.subplots_adjust(wspace=0.5)
plt.show()