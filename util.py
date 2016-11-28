# some useful functions for testing things that are used frequently.

import numpy as np
import h5py 
from ctypes import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt, convolve

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

def peakdet_cwt(infile, threshold):
   # do a first check for peaks in the dataset. After finding peaks, should create a list of 
   # 'event' objects that will be modified as the data is further processed.
   width = np.array([1,10,20,30,40,50])
   vector = pull(infile, 233)
   candidates = find_peaks_cwt(vector,width)
   fig, axis = plt.subplots(1,1)	
   plot(vector, axis)
   axis.scatter(candidates, vector[candidates], marker='o', color='g', s=40)
   plt.show()
   return candidates

def img_derivative(infile, pixel):
	Y = (-1)*pull(infile,pixel)
	#Obtaining derivative
	kernel = [1, 0, -1]
	dY = convolve(Y, kernel, 'valid') 

	#Checking for sign-flipping
	S = np.sign(dY)
	ddS = convolve(S, kernel, 'valid')

	#These candidates are basically all negative slope positions
	#Add one since using 'valid' shrinks the arrays
	candidates = np.where(dY < 0)[0] + (len(kernel) - 1)

	#Here they are filtered on actually being the final such position in a run of
	#negative slopes
	peaks = sorted(set(candidates).intersection(np.where(ddS == 2)[0] + 1))

	plt.step(np.arange(len(Y)), Y)

	#If you need a simple filter on peak size you could use:
	alpha = -0.003
	peaks = np.array(peaks)[Y[peaks] < alpha]

	plt.scatter(peaks, Y[peaks], marker='x', color='g', s=40)
	plt.show()
#def fit_pulse():
	# for each peak, fit the pulse so we can extract the time constant and apply the correct shaping filter.
	

#def event_shaper(infile, pixel, evt_pt, tau):
   # apply trapezoidal filter to data on event-by-event basis
   # extract the pulse height and width.

#def push(infile, pixel):

def pull(infile, pixel):
	# retrieve pixel signal data into numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf:
		d = hf.get(event)
		data = np.array(d)

	return data[pixel]


def plot(data, axis):
	axis.step(np.arange(len(data)), data)

#def img_plot(infile, point):

def signal_plot(infile, pixel):
	event = 'C0'
	channel = 0
	nPix = 72**2
	dump = 0

	#open the file
	with h5py.File(infile,'r') as hf:
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
