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
v = 258.9704
q = 0.5*rho*(v**2)
halfspan = 33.5
centroid = 14.4486
n = 2.5

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
    return yCl_result(x) * q * ychord_result(x) * n
liftdistributionlst = np.array([])
for element in yvalues:
    liftdistributionlst = np.append(liftdistributionlst, Ldistribution(element))

def sheardistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    sheardistributionlst = np.flip(shear)
    return sheardistributionlst

def momentdistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    moment = integrate.cumtrapz(shear, yvalues, initial=0)
    momentdistributionlst = -1 * np.flip(moment)
    return momentdistributionlst

def calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y):
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

    # Multiple spars
    l_spar1 = np.zeros(n_spar - 1)
    x_spar1 = np.zeros(n_spar - 1)
    z_spar1 = np.zeros(n_spar - 1)
    t_spar = np.zeros(n_spar - 1)
    m_up =  (0.045 - 0.0665)/ 0.5
    m_down = (0.0417 - 0.0218)/ 0.5
    l_moi = np.zeros(n_spar - 1)
    h_moi = np.zeros(n_spar - 1)
    if n_spar > 2:
        for i in range(n_spar - 2):
            l_spar1[i] = c * ((0.0665 - (i * m_up * 0.5 /(n_spar - 1))) + (0.0417 - (i * m_down * 0.5 /(n_spar - 1))))
            x_spar1[i] = c * i * 0.5 / (n_spar - 1)     
            z_spar1[i] = (c * (i * m_down * 0.5 /(n_spar - 1))) + (l_spar1[i]) * 0.5
            t_spar[i] = t
            l_moi[i] = t
            h_moi[i] = c * ((0.0665 - (i * m_up * 0.5 /(n_spar - 1))) + (0.0417 - (i * m_down * 0.5 /(n_spar - 1))))

    # Centroids and areas
    x_centroid = np.array([0, 0.5 * c, 0.5 * c, 0.5 * c * 0.5])
    x_centroids = np.concatenate((x_centroid, x_spar1))  # Taking into account spars
    z_centroid = np.array([0.5 * f_spar, (0.0417 * c) + (0.5 * 0.0668 * c), f_spar - ((0.0665 - 0.0450) * c * 0.5), (0.0417 - 0.0218) * c * 0.5])
    z_centroids = np.concatenate((z_centroid, z_spar1))
    l_part = np.array([f_spar, r_spar, np.sqrt((0.5 * c)**2 + ((0.0665 - 0.0450) * c)**2), np.sqrt((0.5 * c)**2 + ((0.0417 - 0.0218) * c)**2)])
    l_parts = np.concatenate((l_part, l_spar1))
    t_part = np.array([t, t, w_u, w_d])
    t_parts = np.concatenate((t_part, t_spar))
    l_x1 = np.array([t, t, np.sqrt((0.5 * c)**2 + ((0.0665 - 0.0450) * c)**2),np.sqrt((0.5 * c)**2 + ((0.0417 - 0.0218) * c)**2) ])
    l_x = np.concatenate((l_x1, l_moi))
    h_x1 = np.array([f_spar, r_spar, w_u, w_d])
    h_x = np.concatenate((h_x1, h_moi))

    # Moment of inertia calculation
    I_x = np.zeros(n_spar + 2)
    centroid_x = np.sum(x_centroids * (l_parts * t_parts)) / np.sum(l_parts * t_parts)
    centroid_z = np.sum(z_centroids * (l_parts * t_parts)) / np.sum(l_parts * t_parts)
    for i in range(n_spar + 2):
        I_x[i] = (l_x[i] * h_x[i]**3 / 12) + (l_parts[i] * t_parts[i]) * ((z_centroids[i] - centroid_z)**2)

    # Stringer contributions
    z_str1 = f_spar - ((0.0665 - 0.0450) * c * 0.5)
    z_str2 = (0.0417 - 0.0218) * c * 0.5
    for i in range(n_str1):
        I_x = np.append(I_x, A1 * ((z_str1 - centroid_z)**2))
        I_x = np.append(I_x, A1 * ((z_str2 - centroid_z)**2))

    return np.sum(I_x)

# Input values
n_spar = int(input('Enter the number of spars: '))
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
    momentdistributionlst = momentdistribution(yvalues)
    integrandlst = np.array([])
    n = 0
    for element in y:
        momentvalues = momentdistributionlst[n]
        n = n+1
        integrandlst = np.append(integrandlst, (momentvalues *- 1)/(calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1,element) * 69*10**9))
    return integrandlst
0.0
def deflection(y):
    def load(y):
        return integrate.cumtrapz(load_integrand(y), y, initial = 0)
    deflection_result = integrate.cumtrapz(load(y), y, initial = 0)
    return deflection_result

plt.plot(yvalues, deflection(yvalues), "b")
plt.xlabel('Spanwise location [m]')
plt.ylabel('Deflection [m]')
plt.title('Deflection Graph')
plt.show()
