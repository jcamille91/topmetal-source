# some useful functions that are used frequently.

### 2.7 vs 3

# 2.7 xrange becomes range in 3
# 2.7 raw_input() becomes input in 3()

# test for slow pieces of code
import time

# for importing data from HDF5 files into Numpy arrays
import h5py 

# Ctypes and Numpy support for calling C functions defined in shared libraries.
# These functions can be slow when implemented in Python.
from ctypes import * 
import numpy.ctypeslib as npct

# plotting tools. so far, sticking to matplotlib.
import matplotlib
matplotlib.use('TkAgg') # matplotlib backend for plotting interactively. 
						# configure this before importing pylab/pyplot.
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import axes
ax_obj = axes.Axes
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LinearSegmentedColormap

# Intel's OpenCV functions for image processing.
import cv2

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
# to specify and remember indices/ordering.
from collections import namedtuple
wfm = namedtuple('wfm', 'avg rms data') #  waveform data + info.
trig = namedtuple('trig', 'mean dY S ddS cds peaks toss pkm') # each step of peak detection algorithm.
PEAKS = namedtuple('PEAKS', 'none few mid many pkch') # for testing peak detection performance of all 72x72 sensor channels

class Peak(object):

	def __init__(self)
		self.dude = 0 

### a note on the two types of python pointer functionalities used here ###
# 1. use Ctypes 'POINTER' functionality for defining/building data structures
# to use between C and Python.
# 2. use Numpy.Ctypeslib's 'ndpointer' to easily pass and return
# pointer-arrays between C and python for quick analysis in python environment.

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

	''' this function converts numpy arrays into the appropriate 
		pointers for the peaks_t structure '''

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


def get_wfm_one(infile, ch, npt, plt) :

	''' get data for a channel and calculate its average and root mean square. 
	npt defaults to max value (25890 for current dataset) for zero input or too large of input.'''
	data = pull_one(infile, ch)

	avg = np.mean(data[:npt])
	rms = np.std(data[:npt])

	if plt == True :
		print 'average = ', avg, 'Volts'
		print 'sigma = ', rms, 'Volts RMS'
		plotter(data)


	return wfm(avg=avg, rms=rms, data=data)

def get_wfm_all(infile, npt) :

	''' get data for 72*72 sensor and calculate each channel's average and root mean square voltage. 
	npt defaults to max value (25890 for current dataset) for zero input or too large of input.'''
	dead = 3
	nch = 72**2
	avg = np.zeros(nch)
	rms = np.zeros(nch)

	data = pull_all(infile)
	length = len(data[0])

	if ((npt == False) or (npt > length)) :
		print 'set calculation length to raw data array length =', length 
		npt = length
	for i in xrange(dead,nch) : # leave channels 0-2 with avg = 0, rms = 0.
		avg[i] = np.mean(data[i]) # they have no signal info, only for identifying demux
		rms[i] = np.std(data[i]) # frames in demux algorithm.

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

	''' 
	input:
	- data : numpy array with peaks to look for.
	- mean: average value of 'data'.
	- threshold : minimum acceptable value (volts) above the average signal level for a peak.
	- minsep : minimum number of samples between peaks.  peaks closer than minsep are discarded.

	*** sgwin and sgorder are filter parameters for smoothing the data with a savitsky-golay filter
	before peak detection. derivative based peak detection on noisy signals requires smoothing. ***

	- sgwin : window of savitsky-golay filter. smaller window picks up smaller features
	and vice-versa. window is number of adjacent points for polynomial fitting.
	- sgorder : order of savitsky-golay polynomial for fitting to data.
	
	return:
	- pkm : return array of peak locations given as indices of the original input array 'data'.
	- trig(optional) : namedtuple for testing each step of peak detection if desired.
	'''

	# sign = 1 does positive data with peaks, sign = -1 does negative data with valleys
	sign = 1

	# smooth data for peak detection. if using non-noisy/ideal/fake data, can skip this step.
	filt = savgol_scipy(data, sgwin, sgorder)
	Y = filt

	# calculate derivative
	kernel = [1, 0, -1]	
	dY = convolve(Y, kernel, 'valid') # note: each convolution cuts length of array by len(kernel)-1
	
	# normalize derivative to one. three values/meanings: 1 is increasing, -1 is decreasing, 0 is constant.
	S = np.sign(dY)
	
	# the second derivative of the normalized derivative.
	# should only have non-zero values for peaks and valleys, where the value of the derivative changes.
	ddS = convolve(S, kernel, 'valid')

	# first, find all of the positive derivative values. going up the peak.
	# this returns indices of possible candidates. we want to offset by two because
	# the convolution cuts the array length by len(kernel)-1
	if (sign == 1) :
		candidates = np.where(dY > 0)[0] + (len(kernel)-1)
	elif (sign == -1) :
		candidates = np.where(dY < 0)[0] + (len(kernel)-1)

	pk = sorted(set(candidates).intersection(np.where(ddS == -sign*2)[0] + 1))
	alpha = mean + (sign * threshold)

	if (sign == 1) :
		pk = np.array(pk)[Y[pk] > alpha]
	elif (sign == -1) :
		pk = np.array(pk)[Y[pk] < alpha]

	# remove peaks within the minimum separation... can do this smarter.
	toss = np.array([])
	for i in xrange(len(pk)-1) :
		# if the peaks are closer than the minimum separation and the second peak is
		# larger than the first, throw out the first peak. 
		if ((pk[i+1]-pk[i]) < minsep) :
		#if ((peaks[i+1]-peaks[i]) < minsep) and (Y[peaks[i+1]] < Y[peaks[i]]) :
			toss = np.append(toss, i+1)

	pkm = np.delete(pk, toss)

	# cons = 5 # consecutively increasing values preceeding a peak
	# for j in xrange(len(pkm))
	# 	for k in xrange(cons)

	# use the 'trig' namedtuple for debugging / accessing each step of the peak detection.
	#return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=pk, toss=toss, pkm=pkm)
	return pkm

def fit_tau(data, avg, rms, peaks, fudge, fit_length, ax) :
	'''
	apply non-linear least squares fitting to one or multiple peaks on a given channel. use 'ax' to optionally
	superimpose fit as a scatterplot onto input data plot.

	inputs:
	- data : array containing pulses to fit
	- avg : baseline of the signal
	- rms : std dev of the signal
	- peaks : array of peak locations
	- fudge :  used in the case that there are consecutive peaks closer than the specified
		fit length. 'fudge' moves current fit 'fudge' many data points back from next peak,
		to try and ensure the following peak doesn't disturb the current peak's fit. 
	- fit_length : number of points data points used for least squares fit.
		no more than about 3x expected tau is typical.
	- ax : supply a an axis from fig, ax = matplotlib.pyplot.subplots(y,y) to 
		superimpose a scatterplot of the fitted data onto the original data.

	return: 
	- tau : array of fitted tau values. same length as number of input peaks.
	'''
	tau = np.zeros(len(peaks), dtype=np.float64)
	chisq = np.zeros(len(peaks), dtype=np.float64)
	Q = np.zeros(len(peaks), dtype=np.float64)
	P = np.zeros(len(peaks), dtype=np.float64)

	# assumption: data is modeled by a 3 parameter decaying exponential.
	M = 3

	# 			   PARAMETER FORMAT: ([A, l, off]
	#				 MIN 					MAX
	bounds = ([0.0, 1.0/100, avg-0.3], [0.03, 1.0/10, avg+0.3])
	guess = [0.008, 1.0/35, avg]
	
	# use fit_length if the peaks are far away. if they are close, 
	# fit only as much data as is available between peaks.

	# exact peak location is fuzzy in presence of noise, so
	# include a 'fudge factor' to improve fit for consecutive peaks.

	npk = len(peaks)
	end = len(data)-1 
	peaks = np.append(peaks, end) 

	for j in xrange(npk) :
		if peaks[j+1]-peaks[j] < fit_length : 
			yi = data[peaks[j]:peaks[j+1]-fudge]
		else :								  
			yi = data[peaks[j]:peaks[j]+fit_length]
		N=len(yi)
		xi = np.arange(N)

		if rms == 0 :
			par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
				              check_finite=False, bounds=bounds, method='trf')
		else :
			par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
				             sigma=np.ones(N)*rms, absolute_sigma=True, check_finite=False, \
				             bounds=bounds, method='trf')

		# par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, \
		# 		             sigma=np.ones(N)*rms, absolute_sigma=True, check_finite=False, \
		# 		             bounds=bounds, method='trf')


		f_xi = model_func(xi, *par)
		Xsq, pval = chisquare(f_obs=yi, f_exp=f_xi, ddof=N-M)
		
		tau[j] = 1.0/par[1]
		# chisq[j] = Xsq
		# P[j] = pval
		# Q[j] = 1-pval

		# if axis object is provided, add the fit as a scatter plot to the axis.
		if isinstance(ax, ax_obj) :
			ax.scatter(xi+peaks[j], model_func(xi,*par), marker = 'o')



	return tau
	#return (tau, chisq, Q, P)

def shaper_single(data, l, k, M, baseline):
	# apply trapezoidal filter to a numpy array, return the numpy array 
	# for quick analysis / plotting.

	# import the library
	lib = CDLL("shaper.so")
	lib.shaper_single.restype = None
						  # 	     in 		out   	  length  	 l  		k 		  M 	baseline
	lib.shaper_single.argtypes = [double_ptr, double_ptr, c_ulong, c_ulong, c_ulong, c_double, c_double]
	# allocate an array to hold output.
	filt = np.empty_like(data)
	lib.shaper_single(data, filt, c_ulong(len(data)), c_ulong(l), c_ulong(k), c_double(M), c_double(baseline))
	return np.array(filt)

def shaper_multi(data, peaks, l, k, M, offset, baseline):

	''' apply trapezoidal filter to data. at each peak location provided, switch
	to the desired filter parameters, l, k, M. also provide the average baseline of
	the input data. l,k,M, should have as many elements as there are peaks. '''

	npk = len(peaks)

	# for now just fix l,k.
	l_arr = np.ones(npk, dtype = c_ulong)*l
	k_arr = np.ones(npk, dtype = c_ulong)*k

	LR = pk2LR(peaks, offset, len(data)-1)

	# print 'l', l
	# print 'k', k
	# print 'M', M
	# print 'number of peaks =', npk
	# print 'LEFT = ', LR[0]
	# print 'RIGHT = ', LR[1]

	PEAK = peaks_handle(npk, LR[0], LR[1], l_arr, k_arr, M)
	out = np.empty_like(data)
	lib = CDLL("shaper.so")
	lib.shaper_multi.restype = None
	lib.shaper_multi.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t), c_double]
	lib.shaper_multi(data, out, c_ulong(len(data)), byref(PEAK), c_double(baseline))

	return np.array(out)

def pk2LR(peaks, offset, end) :
	
	''' put peak locations into arrays of left and rights for trapezoidal shaper.
	 apply desired offset. Both arrays are the same length as input 'peaks' array.
	
	inputs:
	- peaks : numpy array of peak locations.
	- offset : '-' shifts L/R back, '+' shifts L/R forward 
	0- end : length of entire dataset. 25890 for out22.h5 after
	demux.

	return:
	- LEFT and RIGHT : beginning and end points for each set of 
	trapezoidal filter parameters l,k, and M.'''

	# probably a better way to do this function, it 
	# feels sort of clumsy.

	npk = len(peaks)
	LEFT = np.zeros(npk)
	RIGHT = np.zeros(npk)

	for i in xrange(npk-1):
		LEFT[i]  = peaks[i]   + offset
		RIGHT[i] = peaks[i+1] + offset
		
	LEFT[0] = 0
	LEFT[npk-1] = peaks[npk-1] + offset
	RIGHT[npk-1] = end

	# trapezoidal filter uses size_t, or c_ulong, as its datatype
	# for left and right. they index locations possibly larger than int allows.
	LEFT = np.array(LEFT, dtype = c_ulong)
	RIGHT = np.array(RIGHT, dtype = c_ulong)

	return (LEFT, RIGHT)

def savgol_scipy(data, npt, order):
	''' applying savitsky-golay filter to smooth raw data. Currently,
	this is used to aid derivative based peak-detection.

	input:
	- data : raw data as a 1D numpy array.
	- sgwin : window of savitsky-golay filter. smaller window picks up smaller features
	and vice-versa. window is number of adjacent points for polynomial fitting.
	- sgorder : order of savitsky-golay polynomial for fitting to data.

	return:
	- out : 1D numpy array of smoothed data, same length as input.
	'''

	out = savgol_filter(data, npt, order)
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
	''' retrieve signal data for a single pixel into a 1D numpy array.
	input:
	- infile : string specifying desired .h5/.hdf5 file to read.
	- pixel : desired pixel to retrieve. choose 0-5183. (72**2 total)
	'''

	#infile = '../data_TM1x1/out22_dmux' # file for current dataset.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d, dtype=np.float64)
	return data[pixel]

def pull_all(infile):
	
	''' retrieve signal data for all 5184 pixels into 2D numpy array.
	input:
	- infile : string specifying desired .h5/.hdf5 file to read.
	'''

	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d, dtype=np.float64)

	return data

def close_figs(): 
	'''
	closes all active matplotlib.pyplot figures
	'''

	plt.close("all")

def plot(data, axis):
	''' 
	1 dimensional 'step' plot.
	inputs:
	- data: 1D numpy array of data 
	- axis: supply a pyplot axis object to plot to. For use as a quick and dirty
	 plot in ipython, supply 0 for axis to generate a standalone fig, axis plot.

	'''
	if isinstance(axis, ax_obj) :	# if an axis is supplied, plot to it.
		axis.step(np.arange(len(data)), data)

	else :	# no axis, make a quick standalone plot.
		plt.step(np.arange(len(data)), data)
		plt.show()

def plot_multich(data, ch):
	
	'''
	This function plots a series of channels to 
	observe a trend or difference. In ipython, press enter to 
	progess to the next channel.

	input: 
	- data : (5184 x wfm length) 2D numpy array with waveform data for sensor.
	- ch : 1D numpy array of channels to be plotted. 
	'''

	nch = len(ch)
	wfmlen = len(data[0])
	x = np.arange(wfmlen)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_title("testing plot")

	ax.step(x, np.zeros(wfmlen))
	fig.show()
	ax.figure.canvas.draw()

	for i in xrange(nch) :
		raw_input("press enter for next channel")
		plt_ch = ch[i]
		ax.cla()
		ax.set_title("channel no. %i" % plt_ch)
		ax.step(x, data[plt_ch])
		ax.figure.canvas.draw()


def hist_plot(data, nbins, end, axis) :
	''' 
	1 dimensional histogram plot.
	inputs:
	- data: 1D numpy array of data 
	- nbins: number of desired bins.
	- end: cutoff point for histogram x-axis
	- axis: supply a pyplot axis object to plot to. For use as a quick and dirty
	 plot in ipython, supply 0 for axis to generate a standalone fig,axis plot.

	'''
	if isinstance(axis, ax_obj) : # axis supplied
		axis.set_title("Sensor Noise Histogram")
		axis.hist(data, nbins)
		axis.set_xlabel('Volts RMS')
		axis.set_ylabel('# of channels (5181 total)')
		axis.set_title('Sensor Noise')
		axis.set_xlim() # x limits, y limits
		axis.set_ylim()
		axis.grid(True)

	else : 			# no axis supplied, make standalone plot.
		fig = plt.figure(1)
		axis = fig.add_subplot(111)
		axis.set_title("Sensor Noise Histogram")
		axis.hist(data, nbins)
		axis.set_xlabel('Volts RMS')
		axis.set_ylabel('# of channels (5181 total)')
		axis.set_title('Sensor Noise')
		axis.set_xlim() # x limits, y limits
		axis.set_ylim()
		axis.grid(True)
		fig.show()


def locate_pixel(x, y) :
	''' Take x and y location of square pixel array and convert to location in linear 1D array
	input:
	- x, y : cartesian coordinates of pixel in square array.
	'''
	dim = 72
	out=(y*72)+x
	print out
	return out

def arr2square(data) :

	''' reshape a 1D array into a set of rows and columns (2D array).
	input:
	- data : 1D array for reshaping. 'row' determines shape of square 2D array.
	'''

	row = 72 # for a 72x72 pixel array
	data_2d = np.reshape(data, (row, -1)) # convert to square matrix
	return data_2d

def pixel_tri_value(data) :


	''' provide 72X72 data array on channel status. 0 is a channel with no found peaks.
	1 is a channel exceeding the max number of peaks. 0.5 is 'normal' channels.
	
	input:
	1. data: array of length 5184 to be reshaped into 72x72 pixel picture.
	output:
	plots 72x72 pixel picture, with a three bin color map. can customize color map
	if desired. '''

	row = 72
	data_2d = np.reshape(data, (row, -1)) # convert to square matrix
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	ax.set_title("pixel_status")

	# make a simple custom 3 value color map (Red, Green, Blue)

	###			LEGEND			###
	# 0 = BLACK, 0.5 = GREEN, 1.0 = RED
	#    NONE  		 OK			 BUSY        

	colors = [(0,0,0), (0,1,0), (1,0,0)]
	nbins = 3
	cmap_name = 'pix_test'
	mycolormap = LinearSegmentedColormap.from_list(
			   cmap_name, colors, N=nbins)

	im = ax.imshow(data_2d, cmap=mycolormap, vmin=0.0, vmax=1.0)
	fig.show()
	im.axes.figure.canvas.draw()

def pixel_status(data):


	''' Plot number of peaks found on each of the 5184 channels onto a 2D
	pixel array.

	input:
	- data: 1D numpy array of length 5184 to be reshaped into 72x72 pixel picture.
	'''

	row = 72
	data_2d = np.reshape(data, (row, -1)) # convert to square matrix
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	ax.set_title("pixel_status")     

	mn = 0
	mx = np.amax(data)

	im = ax.imshow(data_2d, cmap=cm.hot, vmin=mn, vmax=mx)
	fig.show()
	im.axes.figure.canvas.draw()


def pixelate_single(data, sample):
	''' Take 5184 channel x 25890 data point array and plot desired points in
	time as 5184 pixel array.
	input:
	- data : 2D numpy array (5184 x 25890)
	- sample : desired point in sample space to plot 0-25889
	'''


	# dark = min(data)
	# bright = max(data)
	timestep = (4*72**2)*(3.2*10**-8)
	row = 72
	data_2d = np.reshape(data[:,sample], (row, -1)) # convert to square matrix

	fig, ax = plt.subplots()

	# make bounds between -1 mV and 5 mV.
	im = ax.imshow(data_2d, cmap=cm.RdYlBu_r, vmin=-0.001, vmax=0.005)
	fig.colorbar(im)
	ax.grid(True)
	fig.show()


	plt.close()

def pixelate_multi(data, start, stop, stepsize):

	''' plot successive pixelated images of the 72*72 sensor.
	input:
	- data : a (5184 x number of time samples) numpy array.
	- start and stop : specify desired points in time to plot.
	- stepsize : decides how many of the time points to plot. if '1',
	all of the data is plotted. if '10' for example, each successive
	10th point in time is plotted. stepsize must divide into integers. 
	'''

	sample_time = (4*72**2)*(3.2*10**-8)
	row = 72
	
	a = np.arange(start, stop, stepsize)
	npt = len(a)

	# get the data into (72 x 72 x npt) array
	data_2d = np.zeros((72, 72, npt))

	for i in xrange(npt) :
		data_2d[:,:,i] = np.reshape(data[:,a[i]], (row,-1)) # convert to square matrix


	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	ax.set_title("topmetal data")
	
	#im = ax.imshow(np.zeros((72,72)), cmap=cm.viridis, vmin=-0.001, vmax=0.015)
	#im = ax.imshow(np.zeros((72,72)), cmap=cm.RdYlBu_r, vmin=-0.001, vmax=0.015)
	im = ax.imshow(np.zeros((72,72)), cmap=cm.jet, vmin=0.0, vmax=0.007)
	fig.show()
	im.axes.figure.canvas.draw()

	tstart = time.time()
	for j in xrange(npt) :
		t = j*stepsize*sample_time
		ax.set_title("Time elapsed: %f seconds" % t)
		im.set_data(data_2d[:,:,j])
		im.axes.figure.canvas.draw()

def text_dump(infile, pixel):
	

	''' outputs waveform values textually to the terminal command line
	input:
	- infile : string name for file to retrieve data from
	- pixel : choose a pixel 0-5183 to dump waveform data
	'''

	event = 'C0'
	channel = 0
	nPix = 72**2
	dump = 0

	# open the file for read
	with h5py.File(infile,'r') as hf:
	   d = hf.get(event)
	   data = np.array(d)

	# choose channel and pixel 
	chpix = (channel)*nPix + (pixel) # take ch 0-7, pixels 0-5183 (72**2 - 1)

	samples = np.linspace(0, len(data[0])-1, len(data[0]))

	for i in xrange(len(data[channel])):
	   print i, data[channel][i]


def demuxD(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_dbl.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def demuxF(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_flt.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))


