import numpy as np
from scipy import integrate
from Main_CSV import momentdistribution

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
    return (momentdistribution(y) * -1) / (calculate_moment_of_inertia(t_1, w_u1, w_d1, A1, n_str1,y) * E)

def deflection(y):
    def load(y):
        return integrate.quad(load_integrand, 0, y)[0]

    deflection_result = integrate.quad(load, 0, y)[0]
    return deflection_result