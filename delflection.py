import numpy as np
from csv import reader
import scipy as sp
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import integrate





# Read the CSV file
data = []
with open('A07csv3.csv', 'r') as file:
    csv_reader = reader(file)
    for row in csv_reader:
        data.append(row)

# Defined constants
rho = 1.225
v = 242
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

# Define set of values for y
yvalues = np.arange(0, halfspan, 0.5)

yCl_result = yCl(ylst, Cllst)
ychord_result = ychord(ylst, chordlst)

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

def momentdistribution(z):
    estimatemoment , _ = sp.integrate.quad(sheardistribution, z, halfspan, limit=1000)
    return -1*estimatemoment + moment() if z == 0 else -1*(estimatemoment)







def calculate_moment_of_inertia(t_1, w_u1, w_d1, A1, n_str1, y):
    # Constants
    c_root = 13.4
    c_tip = 3.8
    half_span = 33.45
    m = (c_tip - c_root) / half_span

    # Wing geometry
    c = c_root + (m * y)
    f_spar = 0.1082 * c
    r_spar = 0.0668 * c

    # Dimension multipliers
    t = t_1 * c
    w_u = w_u1 * c
    w_d = w_d1 * c

    # Centroids and areas
    x_centroids = np.array([0, 0.5 * c, 0.5 * c, 0.5 * c * 0.5])
    z_centroids = np.array([0.5 * f_spar, (0.0417 * c) + (0.5 * 0.0450 * c), f_spar - ((0.0665 - 0.0450) * c * 0.5), (0.0417 - 0.0218) * c * 0.5])
    l_parts = np.array([f_spar, r_spar, np.sqrt((0.5 * c)**2 + ((0.0665 - 0.0450) * c)**2), np.sqrt((0.5 * c)**2 + ((0.0417 - 0.0218) * c)**2)])
    t_parts = np.array([t, t, w_u, w_d])
    l_x = np.array([t, t, np.sqrt((0.5 * c)**2 + ((0.0665 - 0.0450) * c)**2),np.sqrt((0.5 * c)**2 + ((0.0417 - 0.0218) * c)**2) ])
    h_x = np.array([f_spar, r_spar, w_u, w_d])

    # Moment of inertia calculation
    I_x = np.zeros(4)
    centroid_x = np.sum(x_centroids * (l_parts * t_parts)) / np.sum(l_parts * t_parts)
    centroid_z = np.sum(z_centroids * (l_parts * t_parts)) / np.sum(l_parts * t_parts)
    for i in range(4):
        I_x[i] = (l_x[i] * h_x[i]**3 / 12) + (l_parts[i] * t_parts[i]) * ((z_centroids[i] - centroid_z)**2)

    # Stringer contributions
    z_str1 = f_spar - ((0.0665 - 0.0450) * c * 0.5)
    z_str2 = (0.0417 - 0.0218) * c * 0.5
    for i in range(n_str1):
        I_x = np.append(I_x, A1 * ((z_str1 - centroid_z)**2))
        I_x = np.append(I_x, A1 * ((z_str2 - centroid_z)**2))

    return np.sum(I_x)


# Input values
t_1 = float(input('Enter the spar thickness: '))
w_u1 = float(input('Enter the thickness of upper skin: '))
w_d1 = float(input('Enter the thickness of lower skin: '))
A1 = float(input('Enter the area of the stringers: '))
n_str1 = int(input('Enter the number of stringers: '))
y = float(input('Enter the spanwise position: '))

# Calculate moment of inertia
#moment_of_inertia = calculate_moment_of_inertia(t_1, w_u1, w_d1, A1, n_str1, y)
#print("Moment of Inertia:", moment_of_inertia)

def load_integrand(y):
    return (momentdistribution(y) * -1) / (calculate_moment_of_inertia(t_1, w_u1, w_d1, A1, n_str1,y) * 69*10**9)

def deflection(y):
    def load(y):
        return integrate.quad(load_integrand, 0, y)[0]

    deflection_result = integrate.quad(load, 0, y)[0]
    return deflection_result

deflectionlst = np.vectorize(deflection)

print(deflectionlst(yvalues))