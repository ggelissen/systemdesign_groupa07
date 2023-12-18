import math
import numpy as np
import matplotlib.pyplot as plt
from Moment import momentone

# Input Variables
K = 4                                           # Clamped on both ends
E = 68.9E9                                      # Pa     
# L = float(input("Enter rib spacing: "))       # m  
y = 1.2 # float(input("Enter spanwise position: "))   # m
M = momentone(y)                           # Nm  

A = 0.001875             #float(input("Cross-sectional area of the stringer (m^2): "))
length =  0.06           #float(input("Length of the stringer (mm): ")) / 1000
width = 0.04             #float(input("Width of the stringer (mm): ")) / 1000
thickness = 0.025        #float(input("Thickness of the stringer (mm): ")) / 1000



def centroid_x_stringer(a, b, t):
    xA1 = (b / 2) * (b * t)
    xA2 = (t / 2) * (a - t) * t
    A12 = (b * t) + (a - t) * t
    return (xA1 + xA2) / A12

def centroid_z_stringer(a, b, t):
    zA1 = (t / 2) * (b * t)
    zA2 = ((a + t) / 2) * (a - t) * t
    A12 = (b * t) + (a - t) * t
    return (zA1 + zA2) / A12



def momentofinertia_xx_stringer(a, b, t):
    I1 = (1 / 12) * b *(t**3) + (b * t) * (centroid_z_stringer(a, b, t) - t/2)**2
    I2 = (1 /12) * t * (a - t)**3 + (a - t) * t * (centroid_z_stringer(a, b, t) - (a + t) / 2)**2
    return I1 + I2

def momentofinertia_zz_stringer(a, b, t):
    I1 = (1 / 12) * t * (b**3) + (b * t) * (centroid_x_stringer(a, b, t) - b / 2)**2
    I2 = (1 / 12) * (a - t) * (t**3) + (a - t) * t * (centroid_x_stringer(a, b, t) - t / 2)**2
    return I1 + I2

def momentofinertia_xz_stringer(a, b, t):
    I1 = (b * t) * (centroid_x_stringer(a, b, t) - b / 2) * -1 * (centroid_z_stringer(a, b, t) - t / 2)
    I2 = (a - t) * t * (centroid_x_stringer(a, b, t) - t / 2) * -1 * (centroid_z_stringer(a, b, t) - (a + t) / 2)
    return I1 + I2


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
    t = t_1 
    w_u = w_u1 
    w_d = w_d1 

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
            l_spar1[i] = c * ((0.0665 - ((i + 1) * m_up * 0.5 /(n_spar - 1))) + (0.0417 - ((i + 1) * m_down * 0.5 /(n_spar - 1))))
            x_spar1[i] = c * (i+1) * 0.5 / (n_spar - 1)     
            z_spar1[i] = (c * ((i + 1) * m_down * 0.5 /(n_spar - 1))) + (l_spar1[i]) * 0.5
            t_spar[i] = t
            l_moi[i] = t
            h_moi[i] = c * ((0.0665 - ((i + 1) * m_up * 0.5 /(n_spar - 1))) + (0.0417 - ((i + 1) * m_down * 0.5 /(n_spar - 1))))

    # Centroids and areas
    x_centroid = np.array([0, 0.5 * c, 0.5 * c * 0.5, 0.5 * c * 0.5])
    x_centroids = np.concatenate((x_centroid, x_spar1))  # Taking into account spars
    z_centroid = np.array([0.5 * f_spar, (0.0417 * c) + (0.5 * 0.0668 * c), f_spar, 0])#f_spar - ((0.0665 - 0.0450) * c * 0.5), (0.0417 - 0.0218) * c * 0.5])
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
    z_str1 = f_spar #- ((0.0665 - 0.0450) * c * 0.5)
    z_str2 = 0 #(0.0417 - 0.0218) * c * 0.5
    for i in range(n_str1):
        I_x = np.append(I_x, A1 * ((z_str1 - centroid_z)**2))
        I_x = np.append(I_x, A1 * ((z_str2 - centroid_z)**2))

    return np.sum(I_x), centroid_z, (f_spar - centroid_z)

Ixx, z_down, z_up = calculate_moment_of_inertia(3, 0.02, 0.025, 0.025, 0.001875, 18, y)


def columnbuckling_stringer(a, b, t, L):
    global K, E, A
    return (K * (np.pi**2) * E * momentofinertia_xx_stringer(a, b, t)) / (L**2 * A)

print('The critical buckling stress is: ',columnbuckling_stringer(length, width, thickness, 1.65) )

def bendingstress_stringer(y, z):
    return (-M * z) / calculate_moment_of_inertia(3, 0.02, 0.025, 0.025, 0.0001875, 18, y)[0]

print('The buckling stress is : ', bendingstress_stringer(y,z_up))

def margin_of_safety_column(L):
    global y, length, width, thickness, z_down
    return columnbuckling_stringer(length, width, thickness, L) / bendingstress_stringer(y, z_down)

# Iterate over rib spacing to find optimal value
for rib_space in np.arange(0.25, 33.5, 0.1):
    if margin_of_safety_column(rib_space) < 1:
        print(f"Rib Spacing: {rib_space-0.1} m")
        if y < 33.5 - (rib_space-0.1):
            y += rib_space-0.1
        else:
            break
    else:
        continue

'''
#--------------------------------------------------------------
#only works if M and Ixx vary. Otherwise it will be straight line 
yvalues = np.arange(0, 33.5, 0.01)

stress_up = []
stress_down = []
yield_stress_tension = []
yield_stress_compress = []

load_factor = float(input("Load Factor: "))
moment=[]
for i in range(len(yvalues)):
    #moment.append(momentfunction(yvalues[i]))
    #moment=momentdistribution(yvalues)
    stress_up.append(bendingstress_stringer(z_up))
    stress_down.append(bendingstress_stringer(z_down))
    yield_stress_tension.append(276000000)
    yield_stress_compress.append(-241000000)
stress_up = np.array(stress_up[:-1])
stress_down = np.array(stress_down[:-1]) * -1
yield_stress_compress = np.array(yield_stress_compress[:-1])
yield_stress_tension = np.array(yield_stress_tension[:-1])
stress_up = stress_up * load_factor
stress_down = stress_down * load_factor

if load_factor>=0:
   compressive_stress = stress_up
   tension_stress = stress_down
   margin_of_safety_compressive = yield_stress_compress / stress_up
   #i_for_compress = np.argmax(margin_of_safety_compressive > 2)
   margin_of_safety_tensional = yield_stress_tension / stress_down
   #i_for_tension = np.argmax(margin_of_safety_tensional > 2)
else:
   compressive_stress = stress_down
   tension_stress = stress_up
   margin_of_safety_compressive = yield_stress_compress / stress_down
   #i_for_compress = np.argmax(margin_of_safety_compressive > 2)
   margin_of_safety_tensional = yield_stress_tension / stress_up
   #i_for_tension = np.argmax(margin_of_safety_tensional > 2)

#check the value
found_compress = False
for i in range(len(compressive_stress)):
    if compressive_stress[i] < -241000000:
        print(f"Value smaller than -241000000 found in compress")
        found_compress = True
        break
found_tension = False
for i in range(len(tension_stress)):
    if tension_stress[i] > 267000000:
        print(f"Value larger than 267000000 found in tension")
        found_tension = True
        break
        
#plt.plot(yvalues, moment)
#plt.xlabel("span")
#plt.ylabel("compress stress")
#plt.title("q")
#plt.show()

plt.plot(yvalues[:-1], compressive_stress)
plt.plot(yvalues[:-1], yield_stress_compress)
plt.xlabel("span")
plt.ylabel("compress stress")
#plt.title("q")
plt.show()

plt.plot(yvalues[:-1], tension_stress)
plt.plot(yvalues[:-1], yield_stress_tension)
plt.xlabel("span")
plt.ylabel("tension stress")
#plt.title("q")
plt.show()

plt.plot(yvalues[:-1], margin_of_safety_compressive)
plt.plot([0,33.5], [1,1],linestyle='--', label='Dashed Line')
plt.xlabel("span")
plt.ylabel("margin of safety compress")
#plt.title("q")
plt.show()

plt.plot(yvalues[:-1], margin_of_safety_compressive)
plt.plot([0,33.5], [1,1],linestyle='--', label='Dashed Line')
plt.xlabel("span")
plt.ylabel("margin of safety tension")
#plt.title("q")
plt.show()
'''