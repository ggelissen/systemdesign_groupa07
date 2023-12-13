import numpy as np 
import scipy as sp
import math
import matplotlib.pyplot as plt


b = 67
Weng = 6630 
grav  = 9.81 
T = 0 
h = 2.127 
engcenter = 1.43 
taper = 0.28
d2r = np.pi/180
rho = 0.324438
V = 258.9704 
S = 574.3 
sweep_LE = 37.12 
cr = 13.4 

yvalues = np.zeros(70)
ylst = np.zeros(70)

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

    return np.sum(I_x)

#upper
alpha = math.acos(0.5/0.5004620365222)
#lower
beta = math.acos(0.5/0.5003958533001)
# Find Am
def enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2):
   #thickness of front spar, thickness of the rear spar, thickness of upper plate, thickness of lower plate, width, length of front spar, length of rear spar
   return (h_length - t_1/2 - t_11/2) * (L_1 + L_2 - 2 * t_2 - 2 * t_3) / 2
#Find integral
def line_integral(t_1,t_11,t_2,t_3,L_1,L_2,option):
   #front spar line integral
   I_1 = (L_1 - t_1 * math.tan(alpha) / 2 - t_2 * math.cos(alpha) / 2 - t_1 * math.tan(beta) / 2 - t_3 * math.cos(beta) / 2)/t_1
   #rear spar line integral
   I_2 = (L_2 - t_11 * math.tan(alpha) / 2 - t_2 * math.cos(alpha) / 2 - t_11 * math.tan(beta) / 2 - t_3 * math.cos(beta) / 2)/t_11
   #upper part line integral
   I_3 = (l_up - t_1 / (2 * math.cos(alpha)) - t_11 / (2 * math.cos(alpha)))/t_2
   #lower part line integral
   I_4 = (l_low - t_1 / (2 *math.cos(beta)) - t_11 / (2 * math.cos(beta)))/t_3
   #Add all together
   I = I_1 + I_2 + I_3 + I_4
   integrals = {
      1: I,
      2: I_1,
      3: I_2
   }
   return integrals.get(option, None)
#Find torsional constant
def torsional_constant_1(t_1,t_11,t_2,t_3,L_1,L_2):
   return 4 * enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2) ** 2 / line_integral(t_1,t_11,t_2,t_3,L_1,L_2,1)
#Find torsional stiffness
def torsional_stiffness_single_cell(t_1,t_11,t_2,t_3,L_1,L_2,y,G):
   #y is the distance from wing root, G is shear modulus
   return G * torsional_constant_1(t_1,t_11,t_2,t_3,L_1,L_2)/y
#-----------------------------------------------------------------------------------------------------------------------------
#Find torsional stiffness distribution for double-cell wing box k(y) over the half wing span by St Venant’s torsion constant
#Find length of the middle spar
def length_of_middle_spar(L_1,L_3,t_2,t_3):
   #L_3 is the distance from the front spar to midline of middle spar
   return L_1 - L_3 * math.tan(alpha) - L_3 * math.tan(beta) - t_2 / (math.cos(alpha) * 2) - t_3 /( math.cos(alpha) * 2 )
def enclosed_area_2(t_1,t_2,t_3,L_1,L_3):
   #L_3 is the distance from the front spar to midline of middle spar
   return (L_3 - t_1 / 2) * (L_1 - t_1 * math.tan(alpha) / 2 - t_2 * math.cos(alpha) / 2 - t_1 * math.tan(beta) / 2 - t_3 * math.cos(beta) / 2 + length_of_middle_spar(L_1,L_3,t_2,t_3))/2
def enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3):
   return enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2) - enclosed_area_2(t_1,t_2,t_3,L_1,L_3)
def rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #t_4 is the thickness of middle spar
   coeff =  np.array([line_integral(t_1,t_11,t_2,t_3,L_1,L_2,2) + (L_3 / math.cos(beta) - t_1 / (2 * math.cos(beta)))/t_3 + length_of_middle_spar(L_1,L_3,t_2,t_3) /t_4 + (L_3 / math.cos(alpha) - t_1 / (2 * math.cos(alpha)))/t_2 , -length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4])
   coeff = coeff / (2 * enclosed_area_2(t_1,t_2,t_3,L_1,L_3) * G)
   coeff = np.append(coeff, -1.)
   return coeff
def rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   coeff =  np.array([-length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 , line_integral(t_1,t_11,t_2,t_3,L_1,L_2,3) + (l_up - L_3 / math.cos(alpha) - t_11 / (2 * math.cos(alpha)))/t_2 + length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 + (l_low - L_3 / math.cos(beta) - t_11 / (2 * math.cos(beta)))/t_3])
   coeff = coeff / (2 * enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3) * G)
   coeff = np.append(coeff, -1.)
   return coeff
def rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   matrix = np.array([[2*enclosed_area_2(t_1,t_2,t_3,L_1,L_3),2*enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3),0.],rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G),rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)])
   righthandside = np.array([1.,0.,0.])
   solution = np.linalg.solve(matrix,righthandside)
   return solution[2]
def torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   G=G
   return 1/(rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)*G)
def torsional_stiffness_double_cell(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,y,G):
   return G * torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)/y
#-----------------------------------------------------------------------------------------------------------------------------
def torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,option,T):
   G=G
   T_1 = T / (G * torsional_constant_1(t_1,t_11,t_2,t_3,L_1,L_2))
   T_2 = T / (G * torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G))
   integrand = {
      1: T_1,
      2: T_2,
   }
   return integrand.get(option, None)

G=26000000000 #Pa

def chord_calc(y):
   c_root = 13.4
   c_tip = 3.8
   half_span = 33.45
   m = (c_tip - c_root) / half_span
   return c_root + (m * y)

torsional_stiffness=[]
integrands=[]
for i in range(len(ylst)):
    h_length = float(0.5*chord_calc(yvalues[i]))
    l_up=float(0.5004620365222*chord_calc(yvalues[i]))
    l_low=float(0.5003958533001*chord_calc(yvalues[i]))
    L_1 = float(0.1082 * chord_calc(yvalues[i]))
    L_2= float(0.0668 * chord_calc(yvalues[i]))
    L_3 = float(spac * chord_calc(yvalues[i]))
    t_1 = t_01 
    t_11 = t_011 
    t_2 = t_02 
    t_3 = t_03
    t_4 = t_04
    T=float(total[i])
    if ylst[i] <= L_4:
        integrands.append(torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,2,T))
    elif ylst[i] > L_4:
        integrands.append(torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,1,T))
    if ylst[i] <= L_4 and ylst[i] != 0:
        torsional_stiffness.append(torsional_stiffness_double_cell(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,ylst[i],G))
    elif ylst[i] > L_4 and ylst[i] != 0:
        torsional_stiffness.append(torsional_stiffness_single_cell(t_1,t_11,t_2,t_3,L_1,L_2,ylst[i],G))

integral_values = sp.integrate.cumtrapz(integrands, x=ylst, initial=0)
integral_values = integral_values*180/math.pi

plt.plot(ylst,torsional_stiffness)
plt.title("Torsional Stiffness")
plt.ylabel("Torsional Stiffness [Nm/rad]")
plt.xlabel("Spanwise location [m]")
plt.show()

plt.plot(ylst, integral_values)
plt.title("Twist Angle Graph")
plt.ylabel("Twist angle [deg]")
plt.xlabel("Spanwise location [m]")
plt.show()

