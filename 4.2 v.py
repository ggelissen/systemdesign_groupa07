import math
import numpy as np
import matplotlib.pyplot as plt


t_1 = float(input('Enter the spar thickness: '))
w_u1 = float(input('Enter the thickness of upper skin: '))
w_d1 = float(input('Enter the thickness of lower skin: '))
A1 = float(input('Enter the area of the stringers: '))
n_str1 = int(input('Enter the number of stringers: '))

def I_xfinal(y):
    c_root = 13.4
    c_tip = 3.8
    half_span = 31.95
    m = (c_tip - c_root) / half_span

    c = c_root + (m * y)
    #c = float(input('Enter the chord length: '))

    #Contant throughout geometry
    f_spar = 0.1082 * c
    r_spar = 0.0668 * c
    h_length = 0.5 * c


    #Values to be inputted
    t = t_1
    w_u = w_u1
    w_d = w_d1
    A = A1
    n_str = n_str1

    X_c = 0
    Z_c = 0

    # Semi-constant values
    x_4 = 0.5 * c * 0.5
    x_3 = 0.5 * c * 0.5
    x_2 = 0.5 * c
    x_1 = 0

    x_centroids = np.array([x_1, x_2, x_3, x_4])
    Q_xcentroid = np.zeros(4)

    z_up = f_spar - ((0.0665 - 0.0450) * c * 0.5)
    z_low = (0.0417 - 0.0218) * c * 0.5
    z_1 = 0.5 * f_spar
    z_2 = (0.0417 * c) + (0.5 * 0.0450 * c)
    z_str1 = z_up  # assumption that it's almost horizontal
    z_str2 = z_low  # same horizontal assumption

    l_up = math.sqrt((0.5 * c)**2 + ((0.0665 - 0.0450) * c)**2)
    l_low = math.sqrt((0.5 * c)**2 + ((0.0417 - 0.0218) * c)**2)

    z_centroids = np.array([z_1, z_2, z_up, z_low])
    Q_zcentroid = np.zeros(4)

    l_parts = np.array([f_spar, r_spar, l_up, l_low])
    t_parts = np.array([t, t, w_u, w_d])
    Area_parts = np.zeros(4)

    # Function to obtain product of centroid and area
    def Q_x(l, t, d):
        A = l * t
        Q = l * t * d
        return Q, A

    # Obtain Qx and A for all parts w/o stringer
    for i in range(4):
        Q_xcentroid[i] = Q_x(l_parts[i], t_parts[i], x_centroids[i])[0]
        Area_parts[i] = Q_x(l_parts[i], t_parts[i], x_centroids[i])[1]

    # Obtain Qz for all parts w/o stringer
    for i in range(4):
        Q_zcentroid[i] = Q_x(l_parts[i], t_parts[i], z_centroids[i])[0]

    # Qz with stringer
    for i in range(n_str):
        Q1 = A * z_str1
        Q2 = A * z_str2
        Q_zcentroid = np.append(Q_zcentroid, [Q1])
        Q_zcentroid = np.append(Q_zcentroid, [Q2])
        Area_parts = np.append(Area_parts, [A])
        Area_parts = np.append(Area_parts, [A])

    # Qx with stringer
    for i in range(n_str):
        x = i * ((0.5 * c) / (n_str - 1))
        Q1 = A * x  # accounting for both lower and upper
        Q_xcentroid = np.append(Q_xcentroid, [Q1])
        Q_xcentroid = np.append(Q_xcentroid, [Q1])

    # Function to obtain centroid of the whole system
    def centroid(Q, A):
        d = sum(Q) / sum(A)
        return d

    X_c = centroid(Q_xcentroid, Area_parts)
    Z_c = centroid(Q_zcentroid, Area_parts)
    #print(Z_c)

    I = 0

    def A_moi(l,h,d):
        I_xx = (l * h**3)/12  # Assuming top almost rectangular
        par_thm = (l * h) * (d**2)
        I = I_xx + par_thm
        return I

    l_x = np.array([t, t, l_up, l_low])
    h_x = np.array([f_spar, r_spar, w_u, w_d])

    def I_str(n,A):
        I_1 = n * (A * ((z_str1 - Z_c) **2))
        I_2 = n * (A * ((z_str2 - Z_c)**2))
        I = I_1 + I_2
        return I

    I_x = np.zeros(4)

    for i in range(4):
        z_centroids[i] = z_centroids[i] - Z_c

    for i in range(4):
        I_x[i] = A_moi(l_x[i], h_x[i], z_centroids[i])

    I_x = np.append(I_x, I_str(n_str,A))

    I_xx = sum(I_x)
    return I_xx
