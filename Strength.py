import numpy as np
import math
from LiftWeight_ShearBendingDiagrams import momentdistribution
from delflection import calculate_moment_of_inertia

def stress_distribution_skin(n_spar, t_1, w_u1, w_d1, A1, n_str1, y,c):
   return momentdistribution(y) * (0.1082 * c - calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y)[1]) / calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, y)[0]
