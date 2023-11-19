from csv import reader
import scipy as sp
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np

# Read the CSV file
data = []
with open('A07csv2.csv', 'r') as file:
    csv_reader = reader(file)
    for row in csv_reader:
        data.append(row)

# Defined constants
rho = 1.225
v = 10
q = 0.5*rho*(v**2)
halfspan = 33.5

# Create arrays for the values in the CSV file
ylst = np.empty((1, 1), float)
chordlst = np.empty((1, 1), float)
Ailst = np.empty((1, 1), float)
Cllst = np.empty((1, 1), float)
ICdlst = np.empty((1, 1), float)
Cmc4lst = np.empty((1, 1), float)

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
yvalues = np.arange(0, halfspan, 1)
yCl_result = yCl(ylst, Cllst)
ychord_result = ychord(ylst, chordlst)
yICd_result = yICd(ylst, ICdlst)
yCmc4_result = yCmc4(ylst, Cmc4lst)

# Functions to calculate distributed load an point load
def Ldistribution(x):
    global q
    return yCl_result(x) * q * ychord_result(x)

def pointload(x):
    global halfspan
    totallift, error = sp.integrate.quad(Ldistribution, 0, halfspan, limit=1000)
    if x == 0:
        return totallift
    else:
        return 0

#def moment(x):
#    return pointload(x) * 

# Functions to define shear and moment distributions
def sheardistribution(y):
    estimateshear,errorshear = sp.integrate.quad(Ldistribution, y , 33.5, limit=1000)
    return estimateshear - pointload(y)

sheardistributionlst = np.array((1, 1), float)
for element in yvalues:
    sheardistributionlst = np.append(sheardistributionlst, sheardistribution(element))

def momentdistribution(z):
    estimatemoment,errormoment = sp.integrate.quad(sheardistribution,z,33.5, limit=1000)
    return estimatemoment

momentdistributionlst = np.array((1, 1), float)
for element in yvalues:
    momentdistributionlst = np.append(momentdistributionlst, momentdistribution(element))

# origin = [0, 0]
# shearpointload = [0, sheardistributionlst[0]]
# momentpointload = [0, momentdistributionlst[0]]

# Plot shear and moment distributions
plt.subplot(1)
plt.plot(yvalues, sheardistributionlst, 'b')
plt.subplot(2)
plt.plot(yvalues, momentdistributionlst, 'b')

plt.show()
