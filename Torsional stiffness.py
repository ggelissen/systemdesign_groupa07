#Merge this file under the Rimaz code
#Find torsional stiffness distribution for single-cell wing box k(y) over the half wing span by St Venant’s torsion constant
# Find Am
def enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2):
   #thickness of front spar, thickness of the rear spar, thickness of upper plate, thickness of lower plate, width, length of front spar, length of rear spar
   Am = (h_length - t_1/2 - t_11/2) * (L_1 + L_2 - 2 * t_2 - 2 * t_3) / 2
   return Am
#Find integral
def line_integral(t_1,t_11,t_2,t_3,L_1,L_2,option):
   #front spar integral
   I_1 = (L_1 - (t_2 + t_3)/2)/t_1
   #rear spar integral
   I_2 = (L_2 - (t_2 + t_3)/2)/t_11
   #upper part integral
   I_3 = (l_up - t_1 * l_up / (2 * h_length) - t_11 * l_up / (2 * h_length))/t_2
   #lower part integral
   I_4 = (l_low - t_1 * l_low / (2 * h_length) - t_11 * l_low / (2 * h_length))/t_3
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
   J = 4 * enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2) ** 2 / line_integral(t_1,t_11,t_2,t_3,L_1,L_2,1)
   return J
#Find torsional stiffness
def torsional_stiffness_single_cell(t_1,t_11,t_2,t_3,L_1,L_2,y,G):
   #y is the distance from wing root, G is shear modulus
   k = G * torsional_constant_1(t_1,t_11,t_2,t_3,L_1,L_2)/y
   return k

#-----------------------------------------------------------------------------------------------------------------------------
#Find torsional stiffness distribution for double-cell wing box k(y) over the half wing span by St Venant’s torsion constant
#Find length of the middle spar
def length_of_middle_spar(L_1,L_3,t_2,t_3):
   #L_3 is the distance from the front spar to midline of middle spar
   L = L_1 * L_3 / h_length - t_2 * l_up / (2 * h_length) - t_3 * l_low / (2 * h_length)
   return L
# Find Am
def enclosed_area_2(t_1,t_2,t_3,L_1,L_3):
   #L_3 is the distance from the front spar to midline of middle spar
   Am = (L_3 - t_1 / 2) *Length_of_middle_spar(L_1,L_3,t_2,t_3)/2
   return Am
def enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3):
   Am = Enclosed_area_1(t_1,t_11,t_2,t_3,L_1,L_2) - Enclosed_area_2(t_1,t_2,t_3,L_1,L_3)
   return Am
def rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #t_4 is the thickness of middle spar
   coeff =  [line_integral(t_1,t_11,t_2,t_3,L_1,L_2,2)/t_1 + ((L_3 * l_low / h_length)-t_1 * l_low/(2 * h_length))/t_3 + length_of_middle_spar(L_1,L_3,t_2,t_3) /t_4 + ((L_3*l_up/h_length)-t_1 * l_up/(2 * h_length))/t_2 , -length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4]
   coeff.append(-2 * enclosed_area_2(t_1,t_2,t_3,L_1,L_3) * G)
   return coeff
def Rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #t_4 is the thickness of middle spar
   coeff =  [-length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 , line_integral(t_1,t_11,t_2,t_3,L_1,L_2,3)/t_11 + (l_up-(L_3 * l_up / h_length)-t_11 * l_up/(2 * h_length))/t_2 + length_of_middle_spar(L_1,L_3,t_2,t_3)/t_4 +(l_low-(L_3 * l_low / h_length)-t_11 * l_low/(2 * h_length))/t_3]
   coeff.append(-2 * enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3) * G)
   return coeff
def rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   #import numpy as np
   matrix = np.array([[2*enclosed_area_2(t_1,t_2,t_3,L_1,L_3),2*enclosed_area_3(t_1,t_11,t_2,t_3,L_1,L_2,L_3),0.],rate_of_twist_1(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G),rate_of_twist_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)])
   righthandside = np.array([1.,0.,0.])
   solution = np.linalg.solve(matrix,righthandside)
   return solution[3]
def torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G):
   J = 1/(rate_of_twist_value(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)*G)
   return J
def torsional_stiffness_double_cell(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,y,G):
   #y is the distance from wing root, G is shear modulus
   k = G * torsional_constant_2(t_1,t_11,t_2,t_3,t_4,L_1,L_2,L_3,G)/y
   return k
