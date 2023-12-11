import math
import numpy as np
import matplotlib.pyplot as plt
from delflection import calculate_moment_of_inertia

# Input Variables
K = 4               # Clamped on both ends
E = 68.9E9        # Pa     
L = 33.5            # m  
M = -3.1E8 / 46        # Nm  

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


Ixx = calculate_moment_of_inertia(3, 0.03, 0.03, 0.03, 0.004, 23, 33.5)


def columnbuckling_stringer(a, b, t):
    global K, E, L, A
    return (K * (np.pi**2) * E * momentofinertia_xx_stringer(a, b, t)) / (L**2 * A)

print(columnbuckling_stringer(length, width, thickness))

def bendingstress_stringer(a, b, t):
    return (M * ) / Ixx

print(bendingstress_stringer(length, width, thickness))

diff = columnbuckling_stringer(length, width, thickness) - bendingstress_stringer(length, width, thickness)
print(diff)

