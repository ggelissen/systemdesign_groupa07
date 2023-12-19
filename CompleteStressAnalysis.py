from csv import reader
import scipy as sp
from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np

# Critical Loading Conditions
v = 258.9704
n = 1.5

rho = 1.225
q = 0.5*rho*(v**2)
halfspan = 33.5
safetyfactor = 1.5
b = 67 
Ww = 38229.5 / 2 
Wf = (125407 + 522.9) / 2 
Weng = 6033 
Wlg = 11383.7 / 2
grav = 9.81 
T = 0
h = 2.127
engcenter = 1.43
taper = 0.28
d2r = np.pi/180
S = 574.3
sweep_LE = 37.12
cr = 13.4
ct = 3.8

# Aerodynamic Properties
CL0 = 0.04647
CL10 = 0.97586
CLD = 0.174647669

# Material / Structural Properties
K = 4
E = 68.9 * 10 ** 9
G = 26 * 10 ** 9
poisson = 0.33
t_c = 0.113
k_v = 2
sigmayield_tens = 276 * 10 ** 6
sigmayield_comp = -276 * 10 ** 6

# Stringer Properties
A = 0.001875             #float(input("Cross-sectional area of the stringer (m^2): "))
length =  0.06           #float(input("Length of the stringer (mm): ")) / 1000
width = 0.04             #float(input("Width of the stringer (mm): ")) / 1000
thickness = 0.025        #float(input("Thickness of the stringer (mm): ")) / 1000


## ------------------- Data Import ------------------- ##


# Create arrays for the values in the CSV file
ylst0 = np.array([])
chordlst0 = np.array([])
Ailst0 = np.array([])
Cllst0 = np.array([])
ICdlst0 = np.array([])
Cmc4lst0 = np.array([])

ylst10 = np.array([])
chordlst10 = np.array([])
Ailst10 = np.array([])
Cllst10 = np.array([])
ICdlst10 = np.array([])
Cmc4lst10 = np.array([])

# Read the CSV file
data0 = []
data10 = []
with open('A07csv0.csv', 'r') as file:
    csv_reader = reader(file)
    for row in csv_reader:
        data0.append(row)

with open('A07csv10.csv', 'r') as file2:
    csv_reader = reader(file2)
    for row in csv_reader:
        data10.append(row)

# Append correct values from csv_reader to arrays
for row in data0[51:81]:      # Range can be adjusted here!
    ylst0 = np.append(ylst0, float(row[0]))
    chordlst0 = np.append(chordlst0, float(row[1]))
    Ailst0 = np.append(Ailst0, float(row[2]))
    Cllst0 = np.append(Cllst0, float(row[3]))
    ICdlst0 = np.append(ICdlst0, float(row[5]))
    Cmc4lst0 = np.append(Cmc4lst0, float(row[7]))

for row in data10[51:81]:      # Range can be adjusted here!
    ylst10 = np.append(ylst10, float(row[0]))
    chordlst10 = np.append(chordlst10, float(row[1]))
    Ailst10 = np.append(Ailst10, float(row[2]))
    Cllst10 = np.append(Cllst10, float(row[3]))
    ICdlst10 = np.append(ICdlst10, float(row[5]))
    Cmc4lst10 = np.append(Cmc4lst10, float(row[7]))

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
yvalues = np.arange(0, halfspan, 0.01)
yCl_result0 = yCl(ylst0, Cllst0)
ychord_result0 = ychord(ylst0, chordlst0)
yICd_result0 = yICd(ylst0, ICdlst0)
yCmc4_result0 = yCmc4(ylst0, Cmc4lst0)

yCl_result10 = yCl(ylst10, Cllst10)
ychord_result10 = ychord(ylst10, chordlst10)
yICd_result10 = yICd(ylst10, ICdlst10)
yCmc4_result10 = yCmc4(ylst10, Cmc4lst10)



## ------------------- Loading Analysis ------------------- ##



# finding entry in list closest to val
def closest(lst, val):
    lst = np.asarray(lst)
    idx = (np.abs(lst - val)).argmin()
    return lst[idx]

# functions to define all loading distribution (ie decreasing triangular shape for dry, const for fuel)

def cts_loaddistr(y):
    f=0
    g=0
    #if y == 0:
    #   f = g = 0
    if y >= 0:
        c = 4 * Ww * grav / b
        a = (-1 * (Ww * grav * 8)) / (b ** 2)
        f = a * y + c
    if y <= b / 4 and y >= 0:
        h = 8 * Wf * (1 - 1/2.56) * grav / b
        d = 8 * Wf * grav / (2.56 * b)
        m = (d - h) / (b / 4)
        g = h + m * y
    if y > b / 4:
        g = 0
    return f + g  #f is structural weight, g is fuel weight


# Angle of Attack
alpha_sin = (CLD-CL0)/(CL10-CL0) * np.sin(10*np.pi/180)
alpha = np.arcsin(alpha_sin)

# Functions to calculate distributed load an point load
def Ldistribution0(x):
    return yCl_result0(x) * q * ychord_result0(x)
liftdistributionlst0 = np.array([])

def Ldistribution10(x):
    return yCl_result10(x) * q * ychord_result10(x)
liftdistributionlst10 = np.array([])

def LdistributionD(x):
    return (Ldistribution0(x) + ((CLD - CL0)/(CL10 - CL0)) * (Ldistribution10(x) - Ldistribution0(x))) * np.cos(alpha)
liftdistributionlst = np.array([])
for element in yvalues:
    liftdistributionlst = np.append(liftdistributionlst, ((LdistributionD(element)) - cts_loaddistr(element)) * n * safetyfactor)
liftdistributionlst = np.flip(liftdistributionlst)    

def sheardistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    sheardistributionlst = np.flip(shear)
    for i in range(len(yvalues)):
        if yvalues[i] >= (b / 2) * 0.35 and yvalues[i] <= halfspan:
            sheardistributionlst[i] = sheardistributionlst[i] - Weng * grav 
    for i in range(len(yvalues)):
        if yvalues[i] >= 5.8 and yvalues[i] <= halfspan:
            sheardistributionlst[i] = sheardistributionlst[i] - 11383.7 / 2 * grav 
    return sheardistributionlst

def y_shear(y, S):
    return sp.interpolate.interp1d(y, S, kind='cubic', fill_value="extrapolate")
shearfunction = y_shear(yvalues, sheardistribution(yvalues))

def momentdistribution(y):
    shear = integrate.cumtrapz(liftdistributionlst, y, initial=0)
    sheardistributionlst = np.flip(shear)
    for i in range(len(yvalues)):
        if yvalues[i] >= (b / 2) * 0.35 and yvalues[i] <= halfspan:
            sheardistributionlst[i] = sheardistributionlst[i] - Weng * grav 
    for i in range(len(yvalues)):
        if yvalues[i] >= 5.8 and yvalues[i] <= halfspan:
            sheardistributionlst[i] = sheardistributionlst[i] - 11383.7 / 2 * grav 
    sheardistributionlst = np.flip(sheardistributionlst)
    moment = integrate.cumtrapz(sheardistributionlst, y, initial=0)
    momentdistributionlst = -1 * np.flip(moment)
    return momentdistributionlst


def y_moment(y, M):
    return sp.interpolate.interp1d(y, M, kind='cubic', fill_value="extrapolate")
momentfunction = y_moment(yvalues, momentdistribution(yvalues))

def chord(y):
    return cr - cr*(1-taper)*(y/(b/2))

sheardist = sheardistribution(yvalues)
sheardist[0] = 0

momentdist = momentdistribution(yvalues)
'''
def momentone(y):
    i=np.where(yvalues == y)[0]
    return momentdist[i]
'''
Tw = Weng*grav*((chord(0.35*b/2)/2)+engcenter)
Tt = T*h*np.cos(sweep_LE*d2r)

torque = []
for element in yvalues:
    if element <= closest(yvalues, (b/2)*0.35) and element >= 0:
        torque.append(Tt-Tw)
    if element > closest(yvalues, (b/2)*0.35):
        torque.append(0)
        
moment = []
for i in range(len(yvalues)):
    moment.append(yCmc4_result10(yvalues[i])*0.5*rho*chord(yvalues[i])*S*v**2)

total_torque = np.array(moment) + np.array(torque)

def y_torque(y, T):
    return sp.interpolate.interp1d(y, T, kind='cubic', fill_value="extrapolate")
torquefunction = y_torque(yvalues, total_torque)



## ------------------- Moment of Inertia ------------------- ##



def calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y):
    # Wing geometry
    m = (ct - cr) / halfspan
    c = cr + (m * y)
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
            x_spar1[i] = c * (i + 1) * 0.5 / (n_spar - 1)     
            z_spar1[i] = (c * ((i + 1) * m_down * 0.5 /(n_spar - 1))) + (l_spar1[i]) * 0.5
            t_spar[i] = t
            l_moi[i] = t
            h_moi[i] = c * ((0.0665 - ((i + 1) * m_up * 0.5 /(n_spar - 1))) + (0.0417 - ((i + 1) * m_down * 0.5 /(n_spar - 1))))

    # Centroids and areas
    x_centroid = np.array([0, 0.5 * c, 0.5 * c, 0.5 * c * 0.5])
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


## ------------------- Buckling Analysis | Web Buckling ------------------- ##


tf = 0.03 #float(input("Enter the thickness of the spar [mm]: "))*(10**-3)
tr = tf
tm = tf
tsk = 0.03 #float(input("Enter the thickness of the skin [mm]: "))*(10**-3)
#nr = int(input("Enter the amount of ribs: "))
ns = 3 #int(input("Enter the amount of spars: "))
stringerarea = 0.001875 # float(input("Enter the cross-sectional area of the stringer [m^2]: "))
stringernumber = 18 #int(input("Enter the number of stringers [#]: "))
y_value = float(input("Enter spanwise position: "))
ribspacing = float(input("Enter rib spacing: "))
z_down = calculate_moment_of_inertia(ns, tf, tsk, tsk, stringerarea, stringernumber, y_value)[1]
z_up = calculate_moment_of_inertia(ns, tf, tsk, tsk, stringerarea, stringernumber, y_value)[2]



a = ribspacing ## assuming equal spacing

kslst_web = [15, 13, 11.8, 11, 10.5, 9.8, 9.7, 9.6]
a_hlst_web = [1, 1.2, 1.5, 1.7, 2, 2.5, 3, 3.5]

interpolation_function_web = sp.interpolate.interp1d(a_hlst_web, kslst_web, kind = 'cubic', fill_value = 'extrapolate')
a_h_interp_web = np.linspace(1, 3.5, 251)
ks_interp_web = interpolation_function_web(a_h_interp_web)

distr_ks_web = []
def webbuckling_crit(t, y):
    tolerance = 1e-10
    h = 0.1082 * chord(y)
    
    current_a_h = round(a/h,2)
    index_match = np.where(np.abs(a_h_interp_web - current_a_h) < tolerance)[0]
    if index_match.size > 0 and current_a_h < 3.5:
        ks = float(ks_interp_web[index_match[0]])
        distr_ks_web.append(ks)
    else:
        ks = 9.6
        distr_ks_web.append(ks)
    return ((E*ks*np.pi**2)/(12*(1 - poisson**2)))*(t/h)**2

def webbuckling_avg(shear, y):
    h = 0.1082 * chord(y) 
    return shear/((tf + tm + tr)*h)

def webbuckling_torque(torque, y):
    wf = (chord(y)/2) - tf - (tm/2)
    wr = (chord(y)/2) - tr - (tm/2)
    Af = chord(y)*t_c*wf
    Ar = chord(y)*t_c*wr
    
    matr_A = np.array([[(chord(y)/2)/(Af*tsk*G) + (chord(y)*t_c)/(2*Af*tf*G) + (chord(y)*t_c)/(2*Af*tm*G), -(chord(y)*t_c)/(2*Af*tm*G), -1], 
    [-(chord(y)*t_c)/(2*Ar*tm*G), (chord(y)/2)/(Ar*tsk*G) + (chord(y)*t_c)/(2*Ar*tr*G) + (chord(y)*t_c)/(2*Ar*tm*G), -1], [2*Af, 2*Ar, 0]])
    matr_B= np.array([0, 0, torque])
    matr_C = np.linalg.solve(matr_A, matr_B)
    
    tau_F = matr_C[0]/tf
    tau_R = matr_C[1]/tr
    tau_M = (matr_C[0] + matr_C[1])/tm
    return tau_F, tau_M, tau_R

tau_cr_f = []
tau_cr_m = []
tau_cr_r = []
tau_avg = []
tau_torque_f = []
tau_torque_m = []
tau_torque_r = []

i = 0
for value in yvalues:
    tau_cr_f.append(webbuckling_crit(tf, value))
    tau_cr_m.append(webbuckling_crit(tm, value))
    tau_cr_r.append(webbuckling_crit(tr, value))
    
    tau_avg.append(webbuckling_avg(shearfunction(value), value))
    tau_torque_f.append(webbuckling_torque(total_torque[i], value)[0])
    tau_torque_m.append(webbuckling_torque(total_torque[i], value)[1])
    tau_torque_r.append(webbuckling_torque(total_torque[i], value)[2])
    
    i += 1 




## ------------------- Buckling Analysis | Skin Buckling ------------------- ##



kclst_skin = [15, 12, 10.7, 8.6, 8, 7.8, 7.4, 7.3, 7.2]
a_hlst_skin = [0.7, 0.8, 0.9, 1.5, 2, 2.5, 3, 3.5, 4]

interpolation_function_skin = sp.interpolate.interp1d(a_hlst_skin, kclst_skin, kind='cubic', fill_value='extrapolate')
a_h_interp_skin = np.linspace(0.7, 4, 331)
ks_interp_skin = interpolation_function_skin(a_h_interp_skin)

distr_kc_skin = []
def skinbuckling_crit(t, y):
    tolerance = 1e-10
    h = chord(y) * 0.5
    if a <= h:
        current_a_h = round(h / a, 2)  # h is the chord, a is the rib pitch
    elif h < a:
        current_a_h = round(a / h, 2)
    index_match = np.where(np.abs(a_h_interp_skin - current_a_h) < tolerance)[0]
    if index_match.size > 0 and current_a_h < 4:
        kc = float(ks_interp_skin[index_match[0]])
        distr_kc_skin.append(kc)
    else:
        kc = 7.2
        distr_kc_skin.append(kc)
    if a <= h:  
        return ((E * 7 * np.pi ** 2) / (12 * (1 - poisson ** 2))) * (t / a) ** 2
    elif h < a:  
        return ((E * 7 * np.pi ** 2) / (12 * (1 - poisson ** 2))) * (t / h) ** 2

def skinbuckling_bendingstress(y, z):
    return (-momentfunction(y) * z) / calculate_moment_of_inertia(ns, tf, tsk, tsk, stringerarea, stringernumber, y)[0]


momentlst = momentdistribution(yvalues)
sigma_cr_lst = []
sigma_lst = []
margin_of_safety_skin = []

for i in range(3350):
    moment_span_location = -1*momentlst[i]
    sigma_cr_lst.append(skinbuckling_crit(tsk, i/100))
    sigma_lst.append(skinbuckling_bendingstress(moment_span_location, i/100))
    margin_of_safety_skin.append(skinbuckling_bendingstress(moment_span_location, i/100)/skinbuckling_crit(tsk, i/100))



## ------------------- Buckling Analysis | Column Buckling ------------------- ##



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

def columnbuckling_crit(a, b, t, L):
    global K, E, A
    return (K * (np.pi**2) * E * momentofinertia_xx_stringer(a, b, t)) / (L**2 * A)

def columnbuckling_bendingstress(y, z):
    return (-momentfunction(y)* z) / calculate_moment_of_inertia(ns, tf, tsk, tsk, stringerarea, stringernumber, y)[0]

def margin_of_safety_column(y, z, L):
    global length, width, thickness
    return columnbuckling_crit(length, width, thickness, L) / columnbuckling_bendingstress(y, z)




## ------------------- Buckling Analysis | Compression/Tension ------------------- ##



stress_up = []
stress_down = []
yield_stress_tension = []
yield_stress_compress = []

#load_factor = float(input("Load Factor: "))

def compressiontension_crit(y, z_comp, z_tens):
    for i in range(len(y)):
        stress_up.append(columnbuckling_bendingstress(z_comp))
        stress_down.append(columnbuckling_bendingstress(z_tens))
        yield_stress_tension.append(sigmayield_tens)
        yield_stress_compress.append(sigmayield_comp)
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

    return compressive_stress, tension_stress

'''
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
'''


## ------------------- Printing Statements ------------------- ##


# Web Buckling
print("\n Web Buckling")
print(f"Critical stress: {np.round(webbuckling_crit(tf, y_value))}")
print(f"True stress: {np.round(webbuckling_avg(shearfunction(y_value), y_value) * k_v + webbuckling_torque(torquefunction(y_value), y_value)[0])}")

# Skin Buckling
print("\n Skin Buckling")
print(f"Critical stress: {np.round(skinbuckling_crit(tsk, y_value))}")
print(f"True stress: {np.round(skinbuckling_bendingstress(y_value, z_up))}")

# Column Buckling
print("\n Column Buckling")
print(f"Critical stress: {np.round(columnbuckling_crit(length, width, thickness, ribspacing))}")
print(f"True stress: {np.round(columnbuckling_bendingstress(y_value, z_up))}")

# Compression
print("\n Compression")
print(f"Critical stress: {sigmayield_comp}")
print(f"True stress: {np.round(columnbuckling_bendingstress(y_value, z_up))}")

# Tension
print("\n Tension")
print(f"Critical stress: {sigmayield_tens}")
print(f"True stress: {np.round(columnbuckling_bendingstress(y_value, z_down))} \n")