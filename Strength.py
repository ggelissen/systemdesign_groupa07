import numpy as np
from LiftWeight_ShearBendingDiagrams import momentdistribution
from delflection import calculate_moment_of_inertia
from ColumnBuckling import bendingstress_stringer

#def stress_distribution_skin(moment, n_spar, t_1, w_u1, w_d1, A1, n_str1, y, c):
#   return moment * (0.1082 * c - calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y)[1]) / calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y)[0]
#def chord(y):
#   c_root = 13.4
#   c_tip = 3.8
#   half_span = 33.45
#   m = (c_tip - c_root) / half_span
#   return c_root + (m * y)
#
yvalues = np.arange(0, 33.5, 0.01)

stress=[]
n_spar=
n_str1=
for i in range(len(yvalues)):
   t_1=
   w_u1=
   w_d1=
   A1=
   stress.append(bendingstress_stringer(calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, yvalues[i])[2]))
