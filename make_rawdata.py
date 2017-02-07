import numpy as np
import h5py
import sys

# define constants related to the data
#dt = 6*10**(-7)       # for TM1x8
dt = 3.2*10**(-8)      # for TM1x1
# nPt = 256*1024**2/16 # for TM1x8, 256MB, 2 byte datapoints. 
nPt = 2**30/2 		   # for TM1x1, 1 GB, 2 byte datapoints.

frameSize = 4*(72**2)
mChSpl = 4
mStart = 913
nFrame = (nPt-mStart)/frameSize
adcX = (-1)*2**-15 # 16 bit ADC resolution
pixel = 2093
channel = 0

# ramp function to be input to the demux
# ramp = np.linspace(int(-0.8/adcX), int((1-0.8)/adcX), nFrame, dtype=np.uint16)


# step function to test the trapezoidal filter algorithm
a = np.ones(nFrame/2, dtype=np.int16)*(-11111)
b = np.ones(nFrame/2, dtype=np.int16)*(-12121)
step = np.concatenate([a,b])

# file to be written to
file = '/Users/josephcamilleri/notebook/topmetal/data_TM1x1/dtest.h5'

# open hdf5 file for read/write
f = h5py.File(file, "r+")
C_raw = f["C0"]
C = C_raw[channel]


# input the fake data so it's as if it had been time multiplexed
# like real data is by the system ADC every clock cycle

# every pixel in a frame gets 4 (mChSpl many) datapoints.
# offset by the pixel location and number of frames deep
for i in np.arange(nFrame):
   index_pixel = mStart + mChSpl*pixel + frameSize*i
   for j in np.arange(mChSpl):
      index_datapoint = index_pixel + j
      C[index_datapoint] = step[i]

# assigning to the hdf5 file a numpy array of the same shape is the only
# way to write. so first we copied the hdf5 object to a np array, 
# did stuff to it, then read it back into the hdf5 object.

C_raw[channel] = C

f.flush()
f.close()
