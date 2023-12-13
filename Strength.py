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

yvalues = np.arange(0, 33.5, 0.01)

stress_up = []
stress_down = []
yield_stres = []
n_spar =
n_str1 =
safety_factor =
load_factor =

for i in range(len(yvalues)):
   t_1 =
   w_u1 =
   w_d1 =
   A1 =
   stress_up.append(bendingstress_stringer(calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, yvalues[i])[2]))
   stress_down.append(bendingstress_stringer(calculate_moment_of_inertia(n_spar, t_1, w_u1, w_d1, A1, n_str1, yvalues[i])[1]))
   yield_stress.append(276000000)
stress_up = np.array(stress_up[:-1])
stress_down = np.array(stress_down[:-1])
yield_stress = np.array(yield_stress[:-1])
stress_up = stress_up * load_factor * safety_factor
stress_down = stress_down * load_factor * safety_factor

if load_factor>=0:
   compressive_stress = stress_up
   tension_stress = stress_down
   margin_of_safety_compressive = - yield_stress / stress_up
   i_for_compress = np.argmax(margin_of_safety_compressive > 2)
   margin_of_safety_tensional = yield_stress / stress_down
   i_for_tension = np.argmax(margin_of_safety_tensional > 2)
else:
   compressive_stress = stress_down
   tension_stress = stress_up
   margin_of_safety_compressive = - yield_stress / stress_down
   i_for_compress = np.argmax(margin_of_safety_compressive > 2)
   margin_of_safety_tensional = yield_stress / stress_up
   i_for_tension = np.argmax(margin_of_safety_tensional > 2)
   
plt.plot(yvalues, compressive_stress)
plt.plot(yvalues, yield_stress)
plt.xlabel("span")
plt.ylabel("compress stress")
#plt.title("q")
plt.show()

plt.plot(yvalues, tension_stress)
plt.plot(yvalues, yield_stress)
plt.xlabel("span")
plt.ylabel("tension stress")
#plt.title("q")
plt.show()

plt.plot(yvalues[:i_for_compress], margin_of_safety_compressive[:i_for_compress])
plt.plot([0,33.5], [1,1],linestyle='--', label='Dashed Line')
plt.xlabel("span")
plt.ylabel("margin of safety compress")
#plt.title("q")
plt.show()

plt.plot(yvalues[:i_for_tension], margin_of_safety_compressive[:i_for_tension])
plt.plot([0,33.5], [1,1],linestyle='--', label='Dashed Line')
plt.xlabel("span")
plt.ylabel("margin of safety tension")
#plt.title("q")
plt.show()
