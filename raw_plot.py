import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt

###          [1]    [2]    [3]
### input filename event channel
###        "abc.h5" "C0"  "0-7"

text = 0
dump = 0
plot = 1
frameSize = 4*(72**2)   # 64*(72**2)
file_name = sys.argv[1]
event = sys.argv[2]
channel = int(sys.argv[3])

with h5py.File(file_name,'r') as hf:
   data = hf.get(event)
   np_data = np.array(data)
   if text:
      print('List of arrays in this file: \n', hf.keys())
      print('file attributes: \n', hf.attrs.keys())
      print('waveform attributes: \n', hf.attrs.get('Waveform Attributes'))
      print('Shape of the dataset: \n', np_data.shape)

samples = np.linspace(0, len(np_data[0])-1, len(np_data[0]))

if dump :
	for i in range(0,200000):
	   print i, np_data[channel][i]

if plot :
	plt.step(samples[:200000], np_data[channel][:200000])
	axes = plt.gca()
	plt.grid()
	plt.show()  
