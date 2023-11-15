#Find torsional stiffness distribution for single-cell wing box k(y) over the half wing span by St Venant’s torsion constant
# Find Am
def Enclosed_Area(t_1,t_2,t_3,w,L_1,L_2):
   #thickness of spar, thickness of upper plate, thickness of lower plate, width, front spar, rear spar
   Am = (w - t_1) * (L_1 + L_2 - 2 * t_2 - 2 * t_3) / 2
   return Am
#Find integral
def line_integral(t_1,t_2,t_3,w,L_1,L_2,L_3,L_4):
   #thickness of spar, thickness of upper plate, thickness of lower plate, width, front spar, rear spar
   #front spar integral
   I_1 = (L_1 - (t_1 + t_2)/2)/t_1
   #rear spar integral
   I_2 = (L_2 - (t_1 + t_2)/2)/t_1
   #upper part integral
   I_3 = (L_3)/t_2
   #lower part integral
   I_4 = (L_4)/t_3
   #Add all together
   I = I_1 + I_2 + I_3 + I_4
   return I
#Find torsional constant distribution
def Torsional_constant(t_1,t_2,t_3,w,L_1,L_2,L_3,L_4):
   J = 4 * Enclosed_Area(t_1,t_2,t_3,w,L_1,L_2) ** 2 / line_integral(t_1,t_2,t_3,w,L_1,L_2,L_3,L_4)
   return J
#Find torsional stiffness distribution
def Torsional_stiffness(t_1,t_2,t_3,w,L_1,L_2,L_3,L_4,y,G):
   #y is the distance from wing root, G is shear modulus
   k = G * Torsional_constant(t_1,t_2,t_3,w,L_1,L_2,L_3,L_4)/y
   return k
