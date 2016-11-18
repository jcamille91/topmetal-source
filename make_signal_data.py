import numpy as np
import h5py

# define constants related to the data
# nPt = 256*1024**2/16 # for TM1x8, 256MB, 2 byte datapoints. 
nPt = 2**30/2 		   # for TM1x1, 1 GB, 2 byte datapoints.
frameSize = 4*72**2
mStart = 913
nFrame = (nPt-mStart)/frameSize
pixel = 4558

a = np.ones(nFrame, dtype=np.float)*0.8
exp = 0.01*np.exp(-(1/40.)*np.linspace(0, 500,500))
a[nFrame/2:(nFrame/2)+len(exp)] += exp


# file to be written to
file = '/Users/josephcamilleri/notebook/topmetal/data_TM1x1/exp.h5'

# open hdf5 file for read/write
f = h5py.File(file, "r+")
C_raw = f["C0"]
C = C_raw[pixel]

C = a

# assigning to the hdf5 file a numpy array of the same shape is the only
# way to write. so first we copied the hdf5 object to a np array, 
# did stuff to it, then read it back into the hdf5 object.

C_raw[pixel] = C

f.flush()
f.close()
