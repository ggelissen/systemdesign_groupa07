#Merge this file under the Rimaz's code
#Find torsional stiffness distribution for single-cell wing box k(y) over the half wing span by St Venant’s torsion constant
#Find constant angle
#import math
#upper
alpha = math.acos(0.5/0.25092279689)
#lower
beta = math.acos(0.5/0.2507907693676)
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
   return L_1 - L_3 * math.sin(alpha) - L_3 * math.sin(beta) - t_2 * math.cos(alpha) / 2 - t_3 * math.cos(alpha) / 2
def enclosed_area_2(t_1,t_2,t_3,L_1,L_3):
   #L_3 is the distance from the front spar to midline of middle spar
   return (L_3 - t_1 / 2) *Length_of_middle_spar(L_1,L_3,t_2,t_3)/2
def enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3):
   return Enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2) - Enclosed_area_2(t_1,t_2,t_3,L_1,L_3)
def rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #t_4 is the thickness of middle spar
   coeff =  [line_integral(t_1,t_11,t_2,t_3,L_1,L_2,2)/t_1 + (L_3 / math.cos(beta) - t_1 / (2 * math.cos(beta))/t_3 + length_of_middle_spar(L_1,L_3,t_2,t_3) /t_4 + (L_3 / math.cos(alpha) - t_1 / (2 * math.cos(alpha)))/t_2 , -length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4]
   coeff.append(-2 * enclosed_area_2(t_1,t_2,t_3,L_1,L_3) * G)
   return coeff
def rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   coeff =  [-length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 , line_integral(t_1,t_11,t_2,t_3,L_1,L_2,3)/t_11 + (l_up - L_3 / math.cos(alpha) - t_11 / (2 * math.cos(alpha))/t_2 + length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 + (l_low - L_3 / math.cos(beta) - t_11 / (2 * math.cos(beta))/t_3]
   coeff.append(-2 * enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3) * G)
   return coeff
def rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #import numpy as np
   matrix = np.array([[2*enclosed_area_2(t_1,t_2,t_3,L_1,L_3),2*enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3),0.],rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G),rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)])
   righthandside = np.array([1.,0.,0.])
   solution = np.linalg.solve(matrix,righthandside)
   return solution[2]
def torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   return 1/(rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)*G)
def torsional_stiffness_double_cell(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,y,G):
   return G * torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)/y
#-----------------------------------------------------------------------------------------------------------------------------
#import scipy as sp
def torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,option):
   T_1 = T / (G * torsional_constant_1(t_1,t_11,t_2,t_3,L_1,L_2))
   T_2 = T / (G * torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G))
   integrand = {
      1: T_1,
      2: I_2,
   }
   return integrand.get(option, None)
def twist_angle(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,y,L_4):
   #L_4 length of the double cell wing box
   twist_angle_at_L_4 = sp.integrate.quad(torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,2),0,L_4)
   if y > L_4:
      angle_diff = sp.integrate.quad(torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,1),L_4,y)
      theta = angle_diff[0] + twist_angle_at_L_4
      return theta
   else:
      theta = sp.integrate.quad(torque_over_GJ(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G,2),0,y)
      return theta[0]


