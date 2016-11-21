# some useful functions for testing things that are used frequently.

import numpy as np
import h5py 
from ctypes import *
import matplotlib.pyplot as plt

 # demux and pre_shaper are pre-processing functions, applied to all data.
def demux(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def smooth(infile, outfile, l, k, M):
   # apply trapezoidal filter with 'rough parameters' to smooth the dataset so we can search for peaks.
   lib = CDLL("shaper.so")
   lib.shaper.argtypes = [c_char_p, c_char_p, c_ulong, c_ulong, c_double]
   lib.shaper(c_char_p(infile), c_char_p(outfile), c_ulong(l), c_ulong(k), c_double(M))

#def find_peaks(infile, outfile):
   # do a first check for peaks in the dataset. After finding peaks, should create a list of 
   # 'event' objects that will be modified as the data is further processed.
    #return array

#def fit_pulse():
	# for each peak, fit the pulse so we can extract the time constant and apply the correct shaping filter.
	

#def event_shaper(infile, pixel, evt_pt, tau):
   # apply trapezoidal filter to data on event-by-event basis
   # extract the pulse height and width.

def push(infile, pixel):

def pull(infile, pixel):
	# retrieve pixel signal data into numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf:
		d = hf.get(event)
		data = np.array(d)

	return data[pixel]

def plot(data):
	plt.step(np.arange(len(data)), data)
	plt.grid()
	plt.show()

#def img_plot(infile, point):

def signal_plot(infile, pixel, dump):
	event = 'C0'
	channel = 0
	#open the file
	with h5py.File(file_name,'r') as hf:
	   d = hf.get(event)
	   data = np.array(d)

	#choose channel and pixel 
	chpix = (channel)*nPix + (pixel) # take ch 0-7, pixels 0-5183 (72**2 - 1)

	samples = np.linspace(0, len(data[0])-1, len(data[0]))

	if dump :
		for i in np.arange(0,len(data[channel])):
		   print i, data[channel][i]


	plt.step(samples, data[chpix])
	axes = plt.gca()
	plt.grid()
	plt.show()  

#def raw_plot():
