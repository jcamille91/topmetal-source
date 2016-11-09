import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm

dump = 0 
frameSize = 4096
nPix = 4096

###			 [1]    [2]    [3]    [4]  ###
###command line args: filename event channel pixel ###

#choose file and event by command line
file_name = sys.argv[1]
event = sys.argv[2]

#open the file
with h5py.File(file_name,'r') as hf:
   d = hf.get(event)
   data = np.array(d)

#choose channel and pixel by command line
channel = int(sys.argv[3])
pixel = int(sys.argv[4])

chpix = (channel)*nPix + (pixel) # take ch 0-7, pixels 0-5183 (72**2 - 1)

samples = np.linspace(0, len(data[0])-1, len(data[0]))

if dump :
	for i in range(0,200000):
	   print i, np_data[channel][i]


plt.step(samples, data[chpix])
axes = plt.gca()
plt.grid()
plt.show()  
