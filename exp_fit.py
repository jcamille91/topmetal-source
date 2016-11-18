import pyplot
import numpy as np
import h5py
from scipy.optimize import curve_fit

def func(x, a, b, c, d):
    return a*np.exp(-c*(x-b))+d

popt, pcov = curve_fit(func, x, y, [100,400,0.001,0])
print popt

plot(x,y)
x=linspace(400,6000,10000)
plot(x,func(x,*popt))
show()

