import csv
import scipy as sp
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import integrate

import numpy as np
file = open('A07csv2.csv')
type(file)

y = []
chord = []
Ai = []
Cl = []
ICd = []
Cmc4 = []
csvreader = list(csv.reader(file))

            ###########EDIT THE VALUES IN THE BRACKET BELOW############
for row in csvreader[51:81]:
    for i in range(0, len(row)):
        row[i] = float(row[i])
    print(row)
    y.append(row[0])
    chord.append(row[1])
    Ai.append(row[2])
    Cl.append(row[3])
    ICd.append(row[5])
    Cmc4.append(row[7])

print(y)
print(Cl)


yCl = sp.interpolate.interp1d(y,Cl,kind='cubic',fill_value="extrapolate")
yICd = sp.interpolate.interp1d(y,ICd,kind='cubic',fill_value="extrapolate")
yCmc4 = sp.interpolate.interp1d(y,Cmc4,kind='cubic',fill_value="extrapolate")

yvalues = np.arange(0, 33.5, 0.01)

plt.figure()
plt.plot(yvalues, yCmc4(yvalues))
plt.show()