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
from matplotlib import axes
ax_class = axes.Axes

# math/matrices, statistics, and fitting
import numpy as np
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt, convolve, savgol_filter

# "harmless" warning when scipy tries to use some external library 
# for least squares fitting. https://github.com/scipy/scipy/issues/5998
import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")  

# named tuples are convenient for returning items from functions without needing
# to specify indices/ordering.
from collections import namedtuple
wfm = namedtuple('wfm', 'avg rms data') #  waveform info.
trig = namedtuple('trig', 'mean dY S ddS cds peaks toss pkm') # each step of peak detection algorithm.

### a note on the two types of python pointer functionalities used here ###
# 1. use Ctypes 'POINTER' functionality for defining/building data structures
# to use between C and Python.
# 2. use Numpy.Ctypeslib's 'ndpointer' to easily pass and return
# data arrays between C and python for quick analysis in python environment.

float_ptr = npct.ndpointer(c_float, ndim=1, flags='CONTIGUOUS')
sizet_ptr = npct.ndpointer(c_ulong, ndim=1, flags='CONTIGUOUS')
double_ptr = npct.ndpointer(c_double, ndim=1, flags='CONTIGUOUS')

c_ulong_p = POINTER(c_ulong)
c_double_p = POINTER(c_double)

# make a peak Ctypes structure for trapezoidal filtering. "Structure" is from Ctypes.
# for more general specification of where to begin/end filter changes.
class peaks_t(Structure):
	_fields_ = [("nPk", c_ulong), ("LEFT", c_ulong_p), ("RIGHT", c_ulong_p),
				("l", c_ulong_p), ("k", c_ulong_p), ("M", c_double_p)]
# #class peak()
# 	### stuff to put inside this class. ###
# 	location
# 	M
# 	chisq
# 	height

def	peaks_handle(nPk, left, right, l, k, M):

	# this function converts numpy arrays into the appropriate 
	# pointers for the peaks_t structure
	nPk = c_ulong(nPk)
	left_p = left.ctypes.data_as(c_ulong_p)
	right_p = right.ctypes.data_as(c_ulong_p)
	l_p = l.ctypes.data_as(c_ulong_p)
	k_p = k.ctypes.data_as(c_ulong_p)
	M_p = M.ctypes.data_as(c_double_p)
	return peaks_t(nPk, left_p, right_p, l_p, k_p, M_p)

	# simple three parameter decaying exponential model for fitting.
def model_func(x, A, l, off) :
	return ( A * np.exp( -l * (x) ) ) + off



def get_wfm(file, ch, npt, plt) :

	""" get data for a channel and calculate its average and root mean square. 
	npt defaults to max value (25890 for current dataset) for zero input or too large of input."""
	data = pull_one(file, ch)

	avg = np.mean(data[:npt])
	rms = np.std(data[:npt])

	if plt == True :
		print 'average = ', avg, 'Volts'
		print 'sigma = ', rms, 'Volts RMS'
		plotter(data)


	return wfm(avg=avg, rms=rms, data=data)

def get_wfm_all(file, npt) :

	""" get data for 72*72 sensor and calculate each channel's average and root mean square voltage. 
	npt defaults to max value (25890 for current dataset) for zero input or too large of input."""

	nch = 72**2
	avg = np.zeros(nch)
	rms = np.zeros(nch)

	data = pull_all(file)
	length = len(data[0])

	if ((npt == False) or (npt > length)) :
		print 'set calculation length to raw data array length =', length 
		npt = length
	for i in np.arange(nch) :
		avg = np.mean(data[i])
		rms = np.std(data[i])

	return wfm(avg=avg, rms=rms, data=data)


def demux(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def smooth(infile, outfile, l, k, M):
   # apply trapezoidal filter with 'rough parameters' to smooth the dataset so we can search for peaks.
   lib = CDLL("shaper.so")
   lib.shaper.argtypes = [c_char_p, c_char_p, c_ulong, c_ulong, c_double]
   lib.shaper(c_char_p(infile), c_char_p(outfile), c_ulong(l), c_ulong(k), c_double(M))

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

def get_peaks(data, mean, threshold, minsep, sgwin, sgorder):

	"""'threshold' in volts above the average signal value
	'minsep' minimum number of samples between peaks. 
	peaks closer than minsep are discarded."""

	# sign = 1 does positive data with peaks, sign = -1 does negative data with valleys
	sign = 1

	Y = data
	filt = savgol_scipy(data, sgwin, sgorder)
	kernel = [1, 0, -1]
	
	# get derivative
	dY = convolve(Y, kernel, 'valid') # note: each convolution cuts length of array by len(kernel)-1
	
	# normalize derivative to 1. three values: increasing 1, decreasing -1, constant 0.
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
	alpha = mean + (sign * threshold)

	if (sign == 1) :
		peaks = np.array(peaks)[Y[peaks] > alpha]
	elif (sign == -1) :
		peaks = np.array(peaks)[Y[peaks] < alpha]

	# remove peaks within the minimum separation... can do this in a smarter way.
	#minsep = 5
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

	# use the 'trig' namedtuple for debugging / accessing each step of the peak detection.
	# return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=peaks, toss=toss, pkm=peaks_minsep)
	return peaks_minsep

def fit_tau(data, peaks, fudge, fit_length, ax) :
	
	# input data and peak locations. return time constants and chi-square values
	# for each corresponding peak so we can filter the peaks.
	tau = np.zeros(len(peaks), dtype=c_double)
	chisq = np.zeros(len(peaks), dtype=c_double)
	avg = np.mean(data[:2000])
	rms = np.std(data[:2000]) # measure rms of dataset for the lsq fit.
	#rms = 0.003 # can fix this to a few mV if we are getting weird results.
	
	# assumption: data is modeled by a 3 parameter decaying exponential.
	M = 3

	# 			   PARAMETER FORMAT: ([A, l, off]
	#				 MIN 					MAX
	bounds = ([0.0, 1.0/100, avg], [0.03, 1.0/10, avg+0.3])
	guess = [0.008, 1.0/35, avg]
	
	# use fit_length if the peaks are far away. if they are close, 
	# fit only as much data as is available between peaks.
	npk = len(peaks)
	end = len(data)-1
	peaks = np.append(peaks, end) 
	for j in np.arange(npk) :
		if peaks[j+1]-peaks[j] < fit_length : 
			yi = data[peaks[j]:peaks[j+1]-fudge]
		else :								  
			yi = data[peaks[j]:peaks[j]+fit_length]
		N=len(yi)
		xi = np.arange(N)	
		par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
				             sigma=np.ones(N)*rms, absolute_sigma=True, check_finite=False, \
				             bounds=bounds, method='trf')

		f_xi = model_func(xi, *par)
		Xsq, P = chisquare(f_obs=yi, f_exp=f_xi, ddof=N-M)
		
		tau[j] = 1.0/par[1]
		chisq[j] = Xsq
	
		if isinstance(ax, ax_class):
			ax.scatter(xi+peaks[j],model_func(xi,*par), marker = 'o')



	return tau
	#return (tau, chisq)

def shaper_peaks(data, peaks, l, k, M, offset):


	### this is broken for an array of only one peak. and zero. ####


	npk = len(peaks)

	# for now just fix l,k.
	l = np.ones(npk, dtype=c_ulong)*100
	k = np.ones(npk, dtype=c_ulong)*20

	# put peak locations into left/right for the filter.
	# offset left/rights of peak with 'off'. use negative sign
	# to offset before the peak, positive for after the peak.
	LEFT = np.zeros(npk, dtype=c_ulong)
	RIGHT = np.zeros(npk, dtype=c_ulong)

	for i in np.arange(npk-1):
		LEFT[i] = peaks[i]
		RIGHT[i] = peaks[i+1]
		
	# LEFT -=	offset
	# RIGHT -= offset
	LEFT[0] = 0
	LEFT[npk-1] = peaks[npk-1]
	RIGHT[npk-1] = len(data)

	print 'l', l
	print 'k', k
	print 'M', M
	print 'number of peaks =', npk
	print 'LEFT = ', LEFT
	print 'RIGHT = ', RIGHT

	PEAK = peaks_handle(npk, LEFT, RIGHT, l, k, M)
	out = np.empty_like(data)
	lib = CDLL("shaper.so")
	lib.shaper_peaks.restype = None
	lib.shaper_peaks.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t)]
	lib.shaper_peaks(data, out, c_ulong(len(data)), byref(PEAK))

	return np.array(out)

def savgol_scipy(array, npt, order):
	
	out = savgol_filter(array, npt, order)
	return out

def peakdet_cwt(data, axis):
   # do a first check for peaks in the dataset. After finding peaks, should create a list of 
   # 'event' objects that will be modified as the data is further processed.
   width = np.array([1,10,20,30,40,50])
   candidates = find_peaks_cwt(data, width)	
   axis.scatter(candidates, data[candidates], marker='o', color='r', s=40)

   return candidates

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
		data = np.array(d, dtype=np.float64)

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

def demuxD(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_dbl.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def demuxF(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_flt.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))


