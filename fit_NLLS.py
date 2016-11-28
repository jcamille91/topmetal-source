import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import util

# REAL NOISEY DATA. guess parameters: [0.008, 0, 1/35., baseline]

# duration 14420-14470. 14430 exactly is the peak
# 8mv amplitude, offset from start to peak is about 9 samples, time constant is about 1/35 samples
start = 14420
peak = 14430
end = 14480

data = util.pull('../data_TM1x1/out22_dmux.h5', 233)
baseline = np.mean(data[:1000])
stddev = np.std(data[:1000])
data = data[peak:end]

x = np.linspace(0, len(data)-1, len(data))
y = data

# FAKE IDEAL DATA. guess parameters [0.01, 0, 1/40., 0.8]

# nPt = 2**30/2 		   # for TM1x1, 1 GB, 2 byte datapoints.
# frameSize = 4*72**2
# mStart = 913
# nFrame = (nPt-mStart)/frameSize


# a = np.ones(500, dtype=np.float)*0.8
# exp = 0.01*np.exp(-(1/40.)*np.linspace(0, 500,500))
# a += exp
# a += np.random.normal(0,0.001,500)

# x = np.linspace(0, len(a)-1, len(a))
# y = a

def func(x, amp, tau, offset):
    return amp*np.exp(-tau*(x))+offset

popt, pcov = curve_fit(func, x, y, [0.008, 1/35., baseline], sigma=stddev)
print popt

#def func(x, amp, delay, tau, offset):
#   return amp*np.exp(-tau*(x-delay))+offset

# popt, pcov = curve_fit(func, x, y, [0.008, 0, 1/35., baseline], sigma=stddev)
# print popt

fig, axis = plt.subplots(1,1)
axis.step(x,y)
axis.scatter(x,func(x,*popt), marker = 'o')
plt.show()
