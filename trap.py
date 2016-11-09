import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import convolve

exp = np.exp((np.arange(300)/300.)[::-1])
exp = 0.01*exp/exp[0]
epts = np.arange(len(exp))

index = np.arange(2000)
step = np.ones(2000) 
step[:500] = 0.828 
step[500:550] = np.linspace(0.828, 0.978, 50)
step[550:] = 0.978

step[1300:1600] += exp 
step[1400:1700] += exp
step[1600:] = 0.978
step += np.random.normal(0, 0.002, 2000)


def plotter(data):
   plot.step(np.arange(len(data)), data)
   plot.grid()
   plot.show()
   
def trapfilter(pulse, n_avg, n_gap):
   out = np.zeros(len(pulse))
   for i in np.arange(len(pulse) - (2*n_avg+n_gap)):
      avg1 = 0
      avg2 = 0 
      for j in np.arange(n_avg): 
         avg1 += pulse[i+j]
         avg2 += pulse[i+j+n_avg+n_gap]          
      out[i] = (avg2 - avg1)/(n_avg)
   return out

def trigger(Y, threshold):

   # image derivative
   kernel = [1,0,-1]
   dY = convolve(Y, kernel, 'valid')
   # this limits derivative flips to only a specific value
   S = np.sign(dY)
   dS = convolve(S, kernel, 'valid')
   candidates = np.where(dY > 0)[0] + (len(kernel)-1)
   
   peaks = sorted(set(candidates).intersection(np.where(dS == 2)[0] + 1))
   peaks = np.array(peaks)[Y[peaks] > threshold]
   return peaks

# simple baseline correction
# step -=  np.mean(step[:20])
filt = trapfilter(step, 10, 10)

