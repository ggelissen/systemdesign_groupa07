import csv
import scipy as sp
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import integrate

import numpy as np
file = open('A07csv2.csv')
type(file)

rho = 1.225
v = 10
q = 0.5*rho*(v**2)
ylst = []
chordlst = []
Ailst = []
Cllst = []
ICdlst = []
Cmc4lst = []
csvreader = list(csv.reader(file))

            ###########EDIT THE VALUES IN THE BRACKET BELOW############
for row in csvreader[51:81]:
    for i in range(0, len(row)):
        row[i] = float(row[i])
    ylst.append(row[0])
    chordlst.append(row[1])
    Ailst.append(row[2])
    Cllst.append(row[3])
    ICdlst.append(row[5])
    Cmc4lst.append(row[7])

def yCl(y, Cl, yvalues):
    yCl = sp.interpolate.interp1d(y,Cl,kind='cubic',fill_value="extrapolate")
    return yCl(yvalues)
def yICd(y,ICd, yvalues):
    yICd = sp.interpolate.interp1d(y,ICd,kind='cubic',fill_value="extrapolate")
    return yICd(yvalues)
def yCmc4(y, cmc4, yvalues):
    yCmc4 = sp.interpolate.interp1d(y,cmc4,kind='cubic',fill_value="extrapolate")
    return yCmc4(yvalues)
def ychord(y, chord, yvalues):
    ychord = sp.interpolate.interp1d(y,chord,kind='cubic',fill_value="extrapolate")
    return ychord(yvalues)

yvalues = np.arange(0, 33.5, 0.5)
def Ldistribution(y):
    global q
    return yCl(ylst, Cllst, y)*q*ychord(ylst, chordlst, y)


sheardistributionlst = []
def sheardistribution(xlst):
    for element in xlst:
        estimateshear,errorshear = sp.integrate.quad(Ldistribution,element,33.5)
        sheardistributionlst.append(estimateshear)
    return sheardistributionlst

momentdistributionlst = []
def momentdistribution(xlst):
    for element in xlst:
        estimatemoment,errormoment = sp.integrate.quad(sheardistribution,element,33.5)
        momentdistributionlst.append(estimatemoment)
    return momentdistributionlst

plt.figure()
plt.plot(yvalues, momentdistribution(yvalues))
plt.show()
