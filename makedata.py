import numpy as np
import h5py
import sys

# define constants related to the data
dt = 6*10**-7
nPt = 256*1024**2/16
frameSize = 4*72**2
nch = 8
mChSpl = 4
mStart = 10418
nFrame = (nPt-mStart)/frameSize
adcX = (-1)*2**-15 # 16 bit ADC resolution
pixel = 2345
channel = 2

# ramp function to be input to the demux
ramp = np.linspace(int(-0.8/adcX), int((1-0.8)/adcX), nFrame, dtype=np.uint16)

# file to be written to
file = '/Users/josephcamilleri/notebook/topmetal/data3/dummy.h5'

# open hdf5 file for read/write
f = h5py.File(file, "r+")
C_raw = f["C0"]
C = C_raw[channel]


# input the fake data so it's as if it had been time multiplexed
# like real data is by the system ADC every clock cycle

for i in np.arange(nFrame):
   index_pixel = mStart + mChSpl*pixel + frameSize*i
   for j in np.arange(mChSpl):
      index_datapoint = index_pixel + j
      C[index_datapoint] = ramp[i]

# assigning to the hdf5 file a numpy array of the same shape is the only
# way to write. so first we copied the hdf5 object to a np array, 
# did stuff to it, then read it back into the hdf5 object.

C_raw[channel] = C

f.flush()
f.close()
