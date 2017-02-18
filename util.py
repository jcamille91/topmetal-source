# some useful functions that are used frequently.

# for importing data from HDF5 files into Numpy arrays
import h5py 

# Ctypes and Numpy support for calling C functions defined in shared libraries.
# These functions can be slow when implemented in Python.
from ctypes import * 
import numpy.ctypeslib as npct

# plotting tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

# math/matrices, statistics, and fitting
import numpy as np
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt, convolve, savgol_filter

# annoying/harmless warning when scipy tries to use some external library 
# for least squares fitting. https://github.com/scipy/scipy/issues/5998
import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")  

# named tuples can be convenient for returning items from functions without needing
# to remember indices/ordering.
from collections import namedtuple

### a note on the two types of python pointer functionalities used here ###
# 1. use Ctypes 'POINTER' functionality for defining/building data structures
# to use between C and Python.
# 2. use Numpy.Ctypeslib's 'ndpointer' to easily pass and return
# pointer-arrays between C and python for quick analysis of data in python environment.

# define pointer types for passing arrays between C and Python
float_ptr = npct.ndpointer(c_float, ndim=1, flags='CONTIGUOUS')
sizet_ptr = npct.ndpointer(c_ulong, ndim=1, flags='CONTIGUOUS')
double_ptr = npct.ndpointer(c_double, ndim=1, flags='CONTIGUOUS')

c_ulong_p = POINTER(c_ulong)
c_double_p = POINTER(c_double)

# make a peak Ctypes structure for trapezoidal filtering. "Structure" is from Ctypes.
class peaks_t(Structure):
	_fields_ = [("nPk", c_ulong), ("LEFT", c_ulong_p), ("RIGHT", c_ulong_p),
				("l", c_ulong_p), ("k", c_ulong_p), ("M", c_double_p)]

def	pk_hdl(nPk, left, right, l, k, M):

	# this function converts numpy arrays into the appropriate 
	# pointers for the peaks_t structure
	nPk = c_ulong(nPk)
	left_p = left.ctypes.data_as(c_ulong_p)
	right_p = right.ctypes.data_as(c_ulong_p)
	l_p = l.ctypes.data_as(c_ulong_p)
	k_p = k.ctypes.data_as(c_ulong_p)
	M_p = M.ctypes.data_as(c_double_p)
	return peaks_t(nPk, left_p, right_p, l_p, k_p, M_p)

def demux(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def demuxD(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_dbl.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def demuxF(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_flt.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def smooth(infile, outfile, l, k, M):
   # apply trapezoidal filter with 'rough parameters' to smooth the dataset so we can search for peaks.
   lib = CDLL("shaper.so")
   lib.shaper.argtypes = [c_char_p, c_char_p, c_ulong, c_ulong, c_double]
   lib.shaper(c_char_p(infile), c_char_p(outfile), c_ulong(l), c_ulong(k), c_double(M))

def savgol_gsl(data, order, der, window):
	# BROKEN 2/13/2017
	# get array with really small first value, then values aren't total bogus but definitely incorrect...

	# apply savitsky-golay 'least squares' filter to a numpy array, return the
	# numpy array for quick analysis.

	lib = CDLL("filters.so")
	lib.savgol_np.restype = c_int
						  # 	  in 		out        length  order   der   window
	lib.savgol_np.argtypes = [double_ptr, double_ptr, c_ulong, c_int, c_int, c_int]
	filt = np.empty_like(data)
	ret = lib.savgol_np(data, filt, c_ulong(len(data)), c_int(order), c_int(der), c_int(window))
	return np.array(filt)

def peak_det(data,threshold, minsep):

	# 'threshold' in volts above the average signal value
	# 'minsep' minimum number of samples between peaks. 
	#  peaks closer than minsep are discarded.

	# sign == 1 does positive data with peaks
	# sign == -1 does negative data with valleys

	# GSL based savgol filter is broken... values aren't ridiculously wrong... but
	# don't follow trend or baseline of the data...
	
	# GSL_Sign, Yuan's convolution by FFT
	trig = namedtuple('trig', 'mean dY S ddS cds peaks toss pkm')
	sign = 1
	threshold = 0.005
	minsep = 5
	# infile = '../data_TM1x1/demuxdouble.h5'
	# fig, axis = plt.subplots(1,1)
	# fig2, axis2 = plt.subplots(1,1)
	# raw = pull_one(infile, 1321)
	# #data = savgol_gsl(rawdata, 4, 0, 15)
	# filt = savgol_scipy(raw,15,4)
	# rcut = raw[9000:9300]
	# fcut = filt[9000:9300]


	# first implement for positive values

	# Y = filt
	# plot(raw, axis)

	Y = data
	filt = savgol_scipy(data,15,4)
	mean = np.mean(filt[:2000])
	kernel = [1,0,-1]
	# derivative
	dY = convolve(Y, kernel, 'valid') # note: each convolution cuts length of array by len(kernel)-1
	# normalize to 1. this array then tells us which direction the data is going, without indication of how quick.
	S = np.sign(dY)
	# the second derivative of the normalized derivative.
	# should only have values for peaks and valleys, where the value of the derivative changes.
	ddS = convolve(S, kernel, 'valid')

	# first, find all of the positive derivative values. going up the peak.
	# this returns indices of possible candidates. we want to offset by two because
	# the convolution cuts the array length by len(kernel)-1
	if (sign == 1) :
		candidates = np.where(dY > 0)[0] + (len(kernel)-1)
	elif (sign == -1) :
		candidates = np.where(dY < 0)[0] + (len(kernel)-1)

	peaks = sorted(set(candidates).intersection(np.where(ddS == -sign*2)[0] + 1))
	alpha = mean + sign*threshold
	# currently peaks is a set, cast as numpy array with a condition that we're 
	# above the minimum threshold for a peak to be valid.
	if (sign == 1) :
		peaks = np.array(peaks)[Y[peaks] > alpha]
	elif (sign == -1) :
		peaks = np.array(peaks)[Y[peaks] < alpha]

	# remove peaks within the minimum separation
	
	minsep = 5
	toss = np.array([0])
	for i in np.arange(len(peaks)-1) :
		# if the peaks are closer than the minimum separation and the second peak is
		# larger than the first, throw out the first peak. 
		if ((peaks[i+1]-peaks[i]) < minsep) :
		#if ((peaks[i+1]-peaks[i]) < minsep) and (Y[peaks[i+1]] < Y[peaks[i]]) :
			toss = np.append(toss, i+1)

	# remove junk element we left to initialize array.
	toss = np.delete(toss, 0)
	peaks_minsep = np.delete(peaks, toss)
		
	# axis.scatter(peaks, Y[peaks], marker='x', color='r', s=40)
	# axis.scatter(peaks_minsep, Y[peaks_minsep], marker='o', color='g', s=40)
	# fig.show()

	return peaks_minsep
	#return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=peaks, toss=toss, pkm=peaks_minsep)


def savgol_scipy(array, npt, order):
	
	out = savgol_filter(array, npt, order)
	return out

def shaper_np(data, l, k, M):
	# apply trapezoidal filter to a numpy array, return the numpy array 
	# for quick analysis / plotting.

	# import the library
	lib = CDLL("shaper.so")
	lib.trapezoid.restype = None
						  # 	  in 		out   	 length  	l  		k 		  M
	lib.trapezoid.argtypes = [float_ptr, float_ptr, c_ulong, c_ulong, c_ulong, c_double]
	# allocate an array to hold output.
	filt = np.empty_like(data)
	lib.trapezoid(data, filt, c_ulong(len(data)), c_ulong(l), c_ulong(k), c_double(M))
	return np.array(filt)

def peakdet_cwt(data, axis):
   # do a first check for peaks in the dataset. After finding peaks, should create a list of 
   # 'event' objects that will be modified as the data is further processed.
   width = np.array([1,10,20,30,40,50])
   candidates = find_peaks_cwt(data,width)	
   axis.scatter(candidates, data[candidates], marker='o', color='r', s=40)

   return candidates

def img_der(data, axis):
	
	Y = (-1)*data
	mean = np.mean(Y[0:100])
	#Obtaining derivative
	kernel = [1, 0, -1]
	dY = convolve(Y, kernel, 'valid') 

	# Checking for sign-flipping
	# returns 1,-1, or 0 depending on the sign of the derivative values.
	# normalizes derivative values to single magnitude.
	S = np.sign(dY)
	ddS = convolve(S, kernel, 'valid')

    # all points where the slope is negative, determined by first derivative.
	candidates = np.where(dY < 0)[0] + (len(kernel) - 1)

	# Here they are filtered on actually being the final such position in a run of
	# negative slopes

	# set is unordered collection of unique elements.

	# get the candidates with the same indices as the locations where the second derivative is
	# equal to two. then put them into order.
	peaks = sorted(set(candidates).intersection(np.where(ddS == 2)[0] + 1))
	
	#for i in np.arange(len(peaks)-1)
	#	if abs(peaks[i]-peaks[i+1])

	# plt.step(np.arange(len(Y)), Y)

	# simple filter on peak size 
	alpha = mean - 0.004
	#alpha = 0
	# make an array out of peaks with the condition that they pass our threshold alpha
	peaks = np.array(peaks)[Y[peaks] < alpha]

	axis.scatter(peaks, -1*Y[peaks], marker='x', color='g', s=40)
	
	return peaks
	#return [dY,S,ddS,candidates, peaks] 



#def push(infile, pixel):

def pull_one(infile, pixel):
	#infile = '../data_TM1x1/out22_dmux'
	# retrieve pixel signal data into 1D numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d, dtype=np.float64)
	return data[pixel]

def pull_all(infile):
	
	# retrieve signal data for all 5184 pixels into 2D numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d, dtype=np.float32)

	return data

def plot(data, axis):
	axis.step(np.arange(len(data)), data)

def plotter(data):
	plt.step(np.arange(len(data)), data)
	plt.show()

#def img_plot(infile, point):

def signal_plot(infile, pixel):
	event = 'C0'
	channel = 0
	nPix = 72**2
	dump = 0

	# open the file for read
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

