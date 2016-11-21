import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm

nPix = 72**2

# for TM1X1 sensor, event and channel are always zero, so we'll
# fix them for now. use them as arguments again if necessary. 

###			  [1]	 [2]  [3]
### command line args: filename pixel dump


#choose file and event by command line
file_name = sys.argv[1]
# event = sys.argv[2]
event = 'C0'

#open the file
with h5py.File(file_name,'r') as hf:
   d = hf.get(event)
   data = np.array(d)

#choose channel and pixel by command line
# channel = int(sys.argv[3])
channel = 0
pixel = int(sys.argv[2])

chpix = (channel)*nPix + (pixel) # take ch 0-7, pixels 0-5183 (72**2 - 1)

samples = np.linspace(0, len(data[0])-1, len(data[0]))

if sys.argv[3] :
	for i in np.arange(0,len(data[channel])):
	   print i, data[channel][i]


plt.step(samples, data[chpix])
axes = plt.gca()
plt.grid()
plt.show()  
