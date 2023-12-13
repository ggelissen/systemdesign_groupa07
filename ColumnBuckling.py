import math
import numpy as np
import matplotlib.pyplot as plt

# Input Variables
K = 4               # Clamped on both ends
E = 68.9E9        # Pa     
L = 33.5            # m  
M = -7.9E8 / 46        # Nm  

A = float(input("Cross-sectional area of the stringer (m^2): "))
length = float(input("Length of the stringer (mm): ")) / 1000
width = float(input("Width of the stringer (mm): ")) / 1000
thickness = float(input("Thickness of the stringer (mm): ")) / 1000



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

    return np.sum(I_x), centroid_z, (f_spar - centroid_z)

Ixx, z_down, z_up = calculate_moment_of_inertia(3, 0.03, 0.03, 0.03, 0.004, 23, 33.5)


def columnbuckling_stringer(a, b, t):
    global K, E, L, A
    return (K * (np.pi**2) * E * momentofinertia_xx_stringer(a, b, t)) / (L**2 * A)

print(columnbuckling_stringer(length, width, thickness))

def bendingstress_stringer(z):
    return (M * z) / Ixx

print(bendingstress_stringer(z_down))
print(bendingstress_stringer(z_up))

#diff = columnbuckling_stringer(length, width, thickness) - bendingstress_stringer()
#print(diff)

