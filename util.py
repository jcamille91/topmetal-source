# some useful functions

# test for slow pieces of code
import time

# for importing data from HDF5 files into Numpy arrays
import h5py 

import pickle
# saving Sensor objects after doing all the analysis, so they can 
# quickly be re-tested without redoing entire filtering procedure.
#### sensor objects with data are large, and there is some issue
####  using pickle for large objects. it is capable of serializingl large objects
#### there's still a bug associated with larger sizes however...


# Ctypes and Numpy support for calling C functions defined in shared libraries.
# These functions can be slow when implemented in Python.

# Ctypes datatypes
from ctypes import c_ulong, c_int, c_double, c_float
# C equivalents =  size_t,  int,   double,   float

# CDLL used to load desired shared library. 
# POINTER creates C pointer object. byref() for passing pointers.
# ctypes Structure object allows creation of C 'struct'.
from ctypes import CDLL, POINTER, byref, Structure

import numpy.ctypeslib as npct

# plotting tools. so far, sticking to matplotlib.
import matplotlib
#matplotlib.use('TkAgg') # matplotlib backend for plotting interactively. 
						# configure this before importing pylab/pyplot.
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import axes
ax_obj = axes.Axes
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Circle

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

def save_object(obj, outfile):
	''' save an object (serialize it) with pickle library.
	input:
	-obj : instance of a python class.
	-outfile: string name, including path, for .pkl file.
	string must end in ".pkl"
	'''
	with open(outfile, 'wb') as output:
		pickle.dump(obj, output)

def read_object(infile):
	with open(infile, 'rb') as read:
		return pickle.load(read)

class Sensor(object):
	
	''' This class defines some basic features of sensor. Creates a
	list of 'Pixel' objects, which each contain a list of 'Peak' objects.
	'''

	# channels with charge possibly trapped. 
	# Slowly rising baselines with magnitude much larger than signal we want to measure.
	CHARGED  = [268,269,340,718,719,805,806,1047,1048,1617,1618,2037,2038,
    2188,2189,2260,2396,2397, 2539, 2540, 2768, 2769, 4868]
	
	# exceptionally noisy channels, larger than signals we are looking for.
	NOISY = [1898, 1899, 1970]

	# these channels are tied to a potential, don't contain any signal.
	# used to identify each frame when demuxing data.
	DEMUXCH = [0,1,2]
	
	# these non-standard channels yield bad results when submitted to analysis routine.
	# just filter them separately, or maybe even set some to 0 because of very negative values
	# that mess up 2d pixel images.
	BAD = CHARGED + NOISY + DEMUXCH
	BAD.sort()

	def __init__(self, infile = '/Users/josephcamilleri/notebook/topmetal/data_TM1x1/demuxdouble.h5'):
		''' 
		input:
		-select : choose to analyze all the channels, or only some. 
		the string 'all' or a list of desired channels are acceptable input.
		-infile : string input, hdf5 file (.h5/.hdf5) where raw data is stored.
		*** currently this is the already demuxed file. we can add demux functionality to this python code. ***
		
		can find data in this folder:
		/Users/josephcamilleri/notebook/topmetal/data_TM1x1

		Sensor attributes:
		-row : dimension of the square array
		-nch : total number of pixels
		-dead : channels used for frame identification in demux. contain no signal data.
		-pix : list of 'Pixel' objects containg pixel raw data
		 and attributes of each pixel.
		-daq_event : name of hdf5 event to be analyzed from 'infile' (currently 25890 datapoints per pixel)
		-daq_length : number of datapoints per pixel in a dataset (daq_event)
		-sensor_channel : in case of a sensor with multiple 72x72 arrays, use this to select one.
		*** currently setup for analysis of just one 72x72 sensor. ***
		-tsample : sampling period for a single pixel. (npt per frame * master sampling period)
		'''
		self.select = None 
		self.infile = infile
		self.row = 72
		self.nch = 72**2
		self.dead = 3
		self.pix = []
		self.daq_event = 'C0' 
		self.daq_length = None
		self.sensor_channel = 0
		self.tsample = (4*72**2)*(3.2*10**-8)

		# demux DAQ parameters
		self.demuxoutfile = '/Users/josephcamilleri/notebook/topmetal/data_TM1x1/demuxdouble2.h5'
		print ("Analyze data from", infile, "for Topmetal sensor with" , self.nch, "channels")

	def demux(rawinfile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
		''' 
		C implemented function for de-multiplexing raw data read from the TCP protocol.
		the pixel array is time multiplexed, so each frame 
		(4 time samples for each of the 5184 pixels) is sequential in memory.
		To get this data in time-order (an array of 5184 channels x total number of samples),
		we need to find the first frame of data acquisition 'mStart' and input parameters
		for the data acquisition.

		input:

		(file reading and writing)
		-infile : input string for .h5/.hdf5 file containing data to be de-multiplexed.
			using "h5dump -A infile.h5" gives information on datatypes and data array shape.
		-outfile : input string for .h5/.hdf5 file to be created containing 
			array of raw data (5184 x daq_length).

		(DAQ parameters)
		-mStart : When the sensor acquires data, there is some dead time before
		everything is synchronized and taking meaningful data. Pixels 0-2 are tied to a potential
		far apart from signal data, so they can easily identify the beginning of a frame. that is, mStart
		is the first of four time samples of pixel 0.
		-mChLen : this value refers to the number of time samples per pixel per frame. in this case,
		there are four time samples per pixel per frame.
		-mNCh : this refers to the number of channels to be output, also referred to as 'pixels'. This topmetal
		sensor is a 72x72 square, so there are 5184 pixels.
		-mChOff, mChSpl, mChLen : mChLen is the numeber of time samples per frame.
		with these four samples, mChSpl are averaged as a data point. mChOff is the number of
		samples offset from the first sample to begin averaging.
		example: mChOff = 1, mChLen = 4, mChSpl = 2. the first sample is skipped, the 2nd and 3rd are averaged,
		and the fourth is ignored.
		'''

		# DAQ parameters for this particular dataset.
		rawinfile = '../data_TM1x1/out22.h5'
		outdmuxfile = '../data_TM1x1/demuxdouble.h5'
		pixel = 2093
		mStart = 913
		mChLen = 4
		mNCh = 72**2
		mChOff = 1
		mChSpl = 2
		frameSize = (72**2)*4

		# call the C function with ctypes
		lib = CDLL("demux.so")
		lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
		lib.demux(c_char_p(self.infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))


	def load_pixels(self, npt = 25890, cut=[]):
		
		'''retrieve the waveform data from .hdf5/.h5 file for each pixel 
		and calculate each pixel's average and root mean square voltage. 
		npt defaults to length of daq event (25890 for current dataset) 
		for zero input or a value exceeding dataset length.
		input:
		-npt :  number of points starting from beginning of dataset to calculate
		the average voltage and RMS value. defaults to length of dataset.
		-cut : a list with two elements: beginning and endpoint of desired data.
		defaults to empty list, which uses the entire dataset.
		*** nothing has been done to make 'cut' work with the analysis chain yet.
		if may not be practical because the savitsky-golay filter and trapezoidal filter will
		give different results if used on the whole dataset versus only a chunk of it.
		'''

		with h5py.File(self.infile,'r') as hf: # open file for read, then close
			d = hf.get(self.daq_event)
			data = np.array(d, dtype=np.float64)
		self.cut = cut
		
		if not cut: # do all of the data

			self.daq_length = len(data[0])


			# iterate over a list of 'Pixel' objects, supply data, avg, and rms.
			if ((npt == False) or (npt >= self.daq_length)) :
				print('set avg/rms calculation length to raw data array length =', self.daq_length)
				self.pix = [Pixel(i, data[i], np.mean(data[i]), np.std(data[i])) for i in range(self.nch)]
			
			else :
				print('set avg/rms calculation length to =', npt, '.')
				print('data array length is ', self.daq_length)
				self.pix = [Pixel(i, data[i], np.mean(data[i][:npt]), np.std(data[i][:npt])) for i in range(self.nch)]
		
		else : # do a portion of the dataset.

			self.daq_length = cut[1]-cut[0]


			# iterate over a list of 'Pixel' objects, supply data, avg, and rms.
			if ((npt == False) or (npt >= self.daq_length)) :
				print('set avg/rms calculation length to raw data array length =', self.daq_length)
				self.pix = [Pixel(i, data[i][cut[0]:cut[1]], np.mean(data[i]), np.std(data[i])) for i in range(self.nch)]
			
			else :
				print('set avg/rms calculation length to =', npt, '.')
				print('data array length is ', self.daq_length)
				self.pix = [Pixel(i, data[i][cut[0]:cut[1]], np.mean(data[i][:npt]), np.std(data[i][:npt])) for i in range(self.nch)]


		### do we want to set average and rms values to 0 for the dead channels or not?
		### depends on how we will do sensor noise calculation later.
		# for j in range(self.dead) :
		# 	self.pix[j].av

		print(self.infile, "contains", self.daq_length, "datapoints per channel with timestep=", self.tsample, "seconds")
		print(npt, "points used in calculation for each channel's baseline and rms in Volts.")


	def input_fakedata(self):
		'''
		use this function to make data to input to the analysis chain.
		if using data without noise, be sure to change the RMS input to
		the fitting procedure.
		'''
		# create list to store data. normally 5184 channels are placed in this list,
		# but one or just several are fine for experimenting with fake data.
		self.pix = []

		# define number of points for dataset.
		# make sure there are enough points to accomodate 
		# whatever features are defined below.
		self.daq_length = 10000

		# specify M=tau, peak location, amplitude, and the number of points to define
		# the peak.
		M_i = np.array([30, 20, 40], dtype=np.float64)
		M_loc = np.array([1000, 5000, 8000])
		M_a = np.array([0.02, 0.015, 0.011])
		npk = len(M_i)
		M_npt = res
		baseline = 0.8


		if mV_noise == 0 :
			noise = 0
		else :
			noise = np.random.normal(loc = 0, scale = mV_noise*.001, size=data_len)

		# build an array with exponentially decaying pulses just defined above.
		exp = np.ones(data_len, dtype=np.float64)*baseline
		for i in range(npk) :
			exp[M_loc[i]:M_loc[i]+M_npt] += M_a[i]*np.exp(-(1.0/M_i[i])*np.linspace(0, M_npt-1, M_npt)) 
		exp += noise



	def noise_hist(self, nbins=2000, end=0.003, axis=0):
		''' make histogram of noise across 5184 channels. 
		1 dimensional histogram plot.
		inputs:
		- data: 1D numpy array of data 
		- nbins: number of desired bins.
		- end: cutoff point for histogram x-axis
		- axis: supply a pyplot axis object to plot to. For use as a quick and dirty
		  plot, supply 0 for axis to generate a standalone fig,axis plot.
		'''

		# build a list of all the RMS noise values, excluding the dead channels.
		rmsvalues = [self.pix[i].rms for i in range(self.dead, self.nch)]

		if isinstance(axis, ax_obj) : # axis supplied
			axis.set_title("Sensor Noise Histogram")
			axis.hist(rmsvalues, nbins)
			axis.set_xlabel('Volts RMS')
			axis.set_ylabel('# of channels (5181 total)')
			axis.set_title('Sensor Noise')
			axis.set_xlim(0,end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig = plt.figure(1)
			axis = fig.add_subplot(111)
			axis.set_title("Sensor Noise Histogram")
			axis.hist(rmsvalues, nbins)
			axis.set_xlabel('Volts RMS')
			axis.set_ylabel('# of channels (5181 total)')
			axis.set_title('Sensor Noise')
			axis.set_xlim(0, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()



	def analyze(self, simple = False, select = [], fake = [], noisethresh = 0.002,
				minsep = 50, threshold = 0.006,  	   			   # peak det
			    fit_length = 300, fudge = 20,   	   			   # lsq fit
			    l = 150, k = 20, M_def = 40, shaper_offset = 20): # shaper
		''' analysis chain: 
		1. find peaks 
		2. fit peaks 
		3. apply trapezoidal shaper with fitted tau to dataset
		
		input:
		-simple : set 'True' to just filter all the channels with the default filter parameters.
		 no peak detection / exponential decay fitting.
		-select : provide python list of desired channels. provide empty list to analyze all channels.
		-noisethresh : threshold in volts to classify channels. channels above this value are labeled with
		a string  as 'noisy'.

		(peak detection)
		-minsep : minimum separation from peaks. if two are more peaks are within minsep,
		only the first peak candidate is saved.
		-threshold : minimum value for accepted peaks. candidate peaks
		with max value smaller than the threshold are rejected.

		(tau fitting)
		-fit_length : number of datapoints to fit. generally (2-3) expected tau in length.
		-fudge : if consecutive peaks are closer than fit_length, the procedure fits the 
		 smaller space between the two consecutive peaks. In this case, 'fudge' many points 
		 are removed from the end of the fit, to try and ensure the tail of the fit 
		 does not catch the rise of the next peak.

		(trapezoidal filter)
		-l, k, M : trapezoidal filter parameters. 
		 l-k specifies the flat-top duration.
		 l specifies the delay of the filter
		 k specifies the rise time to the flat top
		 M specifies the time constant of the pulse to be filtered.
		-shaper_offset : shaper transitions parameters (l, k, M) at each
		 peak location. offset moves the transition fwd/bkwd (+/-) 
		 the specified number of indices. this is useful to avoid 
		 transitioning parameters after the peak has already begun
		 (it won't yield a trapezoidal response) 

		'''
		# these 'Pixel' class attributes are used by methods for analysis.
		# main three 'Pixel' methods are : peak_det(), fit_pulses(), filter_peaks().
		Pixel.daq_length = self.daq_length 
		Pixel.noisethresh = 0.002

		Pixel.threshold = threshold
		Pixel.minsep = minsep

		Pixel.fudge = fudge
		Pixel.fit_length = fit_length

		Pixel.l = l # fixing l and k shaper parameters
		Pixel.k = k
		Pixel.M_def = M_def # default M is no peaks are found.
		Pixel.shaper_offset = shaper_offset
		

		# do analysis on some fake data generated in software. the list 'fake' should
		#  only be defined if Sensor.input_fakedata() was executed.
		if fake : 
			for i in self.pix:
				i.peak_det()
				i.fit_pulses()
				i.filter_peaks()

		if simple : # do simple analysis, filter everything with default l,k,M.
			for i in range(self.nch) :
				self.pix[i].filter_peaks

		elif not select : # if list is empty, analyze all pixels
			
			# remove 'BAD' channels from regular analysis
			ch = [i for i in range(self.nch)]
			ch = list(set(ch) - set(Sensor.BAD))

			# regular analysis
			for i in ch:
				print('channel no. %i' % i)
				self.pix[i].peak_det()
				self.pix[i].fit_pulses()
				self.pix[i].filter_peaks()

			# just default filter bad channels
			for i in Sensor.BAD :
				print('bad channel no. %i' % i)
				self.pix[i].filter_peaks()

		else :
			for i in select : # analyze only pixels in the 'select' list.
				self.pix[i].peak_det()
				self.pix[i].fit_pulses()
				self.pix[i].filter_peaks()

	def label_events(self):

		# just a quick and dirty function to manually input events of
		# interest for future reference 
		# all taken from out22.h5

		# the length of an event should be ~ l+k, determined by the trapezoidal
		# filter parameters.
		self.trap_length = Pixel.l + Pixel.k

		# circle shaped events; use select_circle(x, y, radius, index)
		self.alpha_events = [
		Event(10, 29, 10, 3520), #
		# 4800, next 320, 470, 500, 900
		Event(15, 9, 9, 6510), #
		Event(45, 31, 10, 6580), #
		Event(11, 30, 10, 6845), #
		Event(18, 45, 15, 7020), #
		Event(25, 46, 13, 7505), #
		Event(18, 45, 15, 7020), #
		Event(23, 10, 10, 7865), #
		Event(19, 10, 10, 8510), #
		Event(43, 29, 12, 8910), #
		Event(28, 16, 11, 9095), #
		# CHECKME Event(18, 28, 13, 9545), # this one moves a little bit but is focused well
		Event(36, 29, 12, 10610), #
		Event(32, 10, 10, 10910), #
		Event(39, 14, 12, 11030), # lopsided but focused, rectangle is better here
		Event(42, 17, 12, 12460), #
		Event(22, 48, 12, 12850), #
		Event(9, 28, 9, 13060), #
		Event(42, 17, 12, 12460), #
		Event(18, 28, 10, 13700), #
		Event(20, 9, 9, 14425), # elliptical and focused
		Event(28, 9, 9, 15140), #
		Event(18, 23, 9, 15475), #
		Event(40, 42, 14, 16825), # lots of dark pixels here
		# CHECKME Event(12, 7, 7, 17600), # elliptical and a little broken up
		Event(40, 14, 10, 18345), # nice circle
		Event(42, 30, 10, 18755), # nice circle but not isolated
		Event(41, 32, 9, 20370), # moves a bit and is a little bit separated
		Event(40, 30, 11, 20550), # fat elliptical circle
		# CHECK ME Event(24, 42, 11, 21140), # circle, lots of darker pixels
		# CHECK ME Event(24, 42, 11, 21375), # pretty similar in location/consistency to previous event
		Event(40, 52, 10, 21830), # small ellipse
		Event(20, 11, 10, 22270), # starts circular becomes elliptical
		Event(40, 35, 10, 22675), # big messy circle
		Event(12, 8, 8, 23060), # elliptical blip
		Event(7, 25, 7, 23220), # two disjointed blobs, probably can't capture both with a circle
		Event(16, 7, 7, 23400), #
		Event(40, 12, 10, 23545), # circle, but not totally isolated from other signal
		Event(31, 13, 12, 23705), #
		#CHECK ME Event(26, 6, 6, 23775), # very small circle, not isolated
		Event(25, 54, 13, 23890), # kind of scattered, not a closed circle
		Event(12, 30, 12, 24175), # scattered circle, in happens very close in time to two other events
		Event(13, 29, 12, 24940), # big elliptical blob, focused


		]
		








		# # 'blobs' some  
		# ev = Event(45, 15, 15, 2135) #
		# ev = Event(38, 46, 16, 3020) #
		# self.select_rectangle([17,37] , [19,60]) #
		# ev = Event(50, 31, 14, 16440) # this one definitely moves...
		# ev = Event(24, 43, 14, 16590) # big elliptical, focused
		# ev = Event(22, 38, 9, 19180) # big elliptical
		# ev = Event(1,1,1, 21520) # pretty messy, hard to discren a shape, but definitely an event
		# # moving

		# ev = Event(27, 36, 8, 1800) #
		# self.alpha_events = [ev1, ev2, ev3, ev4]

	
	def vsum_hist(self, show = True, lr = (0, 300), nbins=20, axis = 0):

		'''a word on this measurement: for the experimental setup, each alpha event should deposit all of its
		energy into ionizing air. nearly all charge due to this event should be picked up by the sensor. Therefore,
		if we add all the consituent voltages of an event, we should get a measure of the alpha particle's energy
		spectrum (a "sharp" peak).

		This function makes takes selections of pixels within circles
		to define events. then, over all frames constituting an event, voltage is summed.
		
		input:
		-show : if True, step through a pixel image of each event with a circle enclosing the region of interest.
		'''

		ring = [] # contains info for each as a tuple so we can easily print later.
		self.alphaE = [] # contains voltage summation for a single event.

		# get location for every event. get selection of pixels based on this location for 
		# each event. store the voltage summation for histogram.
		for ev in self.alpha_events :
			

			circle = self.select_circle(ev.x, ev.y, ev.radius)
			# we'll use generator-expressions over list comprehensions here since we
			# don't need to store all the constituent values to be summed.
			vsum = sum(self.pix[i].filt[j] for i in circle for j in range(ev.index, ev.index+self.trap_length)) 
			#vsum = np.sum(np.array([self.pix[i].filt[j] for i in circle for j in range(frames[0], frames[1])]))
			#valcheck = np.array([self.pix[i].filt[j] for i in circle for j in range(frames[0], frames[1])])
			self.alphaE.append(vsum)
			ring.append((ev.index+self.trap_length/2, ev.x, ev.y, ev.radius, vsum))


		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(x=self.alphaE, bins=nbins, range=lr)
			axis.set_xlabel('Volts, summed over event pixels and frames')
			axis.set_ylabel('counts')
			axis.set_title('alpha energy spectrum')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig = plt.figure(1)
			axis = fig.add_subplot(111)
			axis.hist(x=self.alphaE, bins=nbins, range=lr)
			axis.set_xlabel('Volts, summed over event pixels and frames')
			axis.set_ylabel('counts')
			axis.set_title('alpha energy spectrum')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()	 

		# show the region of interest (possible event)
		# use this to verify we are taking data from desired channels.
		# a = np.zeros(5184)
		# a[self.selection] = 1
		# self.pixelate_single(sample = 0, arr=a)

		# show the the middle most frame of events of interest.
		if show :

			fig2, ax2 = plt.subplots(1,1)

			# p is a tuple (midpoint, x, y, r)
			for p in ring:
				input("press enter to show next event:")	
				ax2.cla()
				ax2.set_title('frame no. %i coordinate: (%i, %i) radius: %i | voltage sum = %f volts' % p)
				self.pixelate_single(sample = int(p[0]), arr=[], axis = ax2)
				# add a circle 'artist' to the event we are analyzing
				circ = plt.Circle((p[1], p[2]), p[3], color = 'r', fill=False, linewidth = 1.5, alpha=1)
				ax2.add_artist(circ)
				fig2.show()

	def select_ellipse(self, x, y, a, b, angle) :

		'''make an elliptical selection on the grid.
			
		input:
		-x : coordinate x
		-y : coordinate y
		-a : x axis length defining ellipsoid
		-b : y axis length defining ellipsoid
		-angle : angle of rotation with respect to the x axis. 
		y axis is 90 degrees, -x is 180 degrees, etc. etc.
		'''

		# caluclate the constants
		angle =  angle-np.pi
		sin = np.sin(angle)
		sin2 = (np.sin(angle))**2
		cos = np.cos(angle)
		cos2 = (np.cos(angle))**2
		a2 = a**2
		b2 = b**2

		# make a rectangle based on the major axis
		# needs to be able to account for rotation
		if a > b :
			major = a
		else :
			major = b

		rect = [((i, j), i+j*self.row) for j in range(y-major, y+major) for i in range(x-major,x+major)]

		# cut out the ellipse. this is just an ellipse distance equation with an angle included..
		selection = [rect[i][1] for i in range(len(rect)) if \
		(((((rect[i][0][0]-x)*cos+(rect[i][0][1]-y)*sin2)**2)/a2) + \
		((((rect[i][0][0]-x)*sin-(rect[i][0][1]-y)*cos2)**2)/b2)) < 1]

		return selection

	def select_circle(self, x, y, r) :

		# circular selection.
		rect = [((i, j), i+j*self.row) for j in range(y-r, y+r) for i in range(x-r,x+r)]
		# linear location in 5184 element array of a 'radius' sized rectangle
		#rectangle = [((self.row)*i)+j for j in range(x_0-radius, x_0+radius) for i in range(y_0-radius, y_0+radius)]
		# list of ((x,y), i) tuples -> (xy coordinates, linear index) 
		#xy = [(self.pix[i].loc, i) for i in rectangle]
		# if the element is inside of the circle's radius, include it. We're cutting a circle out of a rectangle
		selection = [rect[i][1] for i in range(len(rect)) if np.sqrt((rect[i][0][0]-x)**2 + (rect[i][0][1]-y)**2) < r]
		
		return selection

		# retrieve all voltage values for the selected pixels and frames.
		# values = np.array([self.pix[i].filt[j] for i in self.selection for j in frames])
	def select_rectangle(self, lr, tb) :

		'''make a rectangular selection of pixels.
		input:
		-lr : the left and right most points as a two element list.
		-tb : the top and bottom most points as a two element list.

		'''

		return [i+j*self.row for j in range(tb[0], tb[1]) for i in range(lr[0], lr[1])]

	def select_test(self, rect=[],circle=[], ellipse=[], arr=[]) :

		if not (rect and circle and ellipse and arr):
			print('supply list of arguments for one of the four selection types:')
			print('rect=[], circle=[], ellipse=[], or a list of custom values "arr".')

	def selection(self, x_0, y_0, radius, frames, select = [], nbins = 2000, alpha = 1, axis = 0) :

		### make histograms for selection of pixels and frames.
		### fit the histogram. Yuan used sum of two gaussian functions, let's
		### see how they unfitted histogram looks first.
		### plot total fit, then plot the contituent fits.

		''' generate a voltage histogram of signal selection inside of a defined circle.
			alternatively, provide a list of pixels to use as a selection. 
			
			*** cartesian coordinates -> origin is top left corner. ***
			*** left is increasing x, down is increasing y. 	    ***
			
			input:
			-x_0, y_0 : center point of a circle to define an event.
			-radius : distance from center to define circle enclosing event.
			-frames : python list oftime points / frames to use in calculation.
			**note: even a single frame input must be input as a list e.g. [12] 
			can choose single point, or multiple. frames do not necessarily need to be contiguous.
			-select : use select to supply a list of desired pixels to plot instead of default circle.
			-nbins : number of bins for voltage histogram.
			-axis : option to supply Axes  for plot with multiple steps.
		'''

		if not select :
			# do the circular selection.

			# linear location in 5184 element array
			self.rectangle = [((self.row)*i)+j for j in range(x_0-radius, x_0+radius) for i in range(y_0-radius, y_0+radius)]
			# list of ((x,y), i) tuples -> (xy coordinates, linear index) 
			self.xy = [(self.pix[i].loc, i) for i in self.rectangle]
			# if the element is inside of the circle's radius, include it.
			self.selection = [self.xy[i][1] for i in range(len(self.xy)) if np.sqrt((self.xy[i][0][0]-x_0)**2 + (self.xy[i][0][1]-y_0)**2) < radius]
			
			# retrieve all voltage values for the selected pixels and frames.
			values = np.array([self.pix[i].filt[j] for i in self.selection for j in frames])

		else :
			# use the selected pixels

			# numpy array of list comprehension of all filtered data to add to histogram.
			self.selection = select
			values = np.array([self.pix[i].filt[j] for i in self.selection for j in frames])
		


		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(values, nbins)
			axis.set_xlabel('Volts RMS')
			axis.set_ylabel('# of channels (5181 total)')
			axis.set_title('Sensor Noise')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig = plt.figure(1)
			axis = fig.add_subplot(111)
			axis.hist(values, nbins)
			axis.set_xlabel('Volts')
			axis.set_ylabel('# of datapoints')
			axis.set_title('filtered signal values')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()	 

		# show the region of interest (possible event)
		# use this to verify we are taking data from desired channels.
		# a = np.zeros(5184)
		# a[self.selection] = 1
		# self.pixelate_single(sample = 0, arr=a)

		# show the frame that with an event of interest.
		fig2, ax2 = plt.subplots(1,1)
		if (len(frames) == 1) :
			self.pixelate_single(sample = frames, arr=[], axis = ax2)
		elif (len(frames) > 1) :
			self.pixelate_single(sample = frames[0], arr=[], axis = ax2)
		# add a circle to the event we are analyzing
		circle = plt.Circle((x_0, y_0), radius, color = 'r', fill=False, linewidth = 1.5, alpha=alpha)
		ax2.add_artist(circle)
		fig2.show()

	def signal_hist(self, nbins = 10000, begin = -0.03, end = 0.03, axis=0) :
		''' generate a histogram with the range of signal values we are dealing with.
		'''
		# use numpy.ravel() for 1d view of 2d array.
		# use numpy.flatten() for 1d copy of 2d array.
		filt1d = self.filt.ravel() 

		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(filt1d, nbins)
			axis.set_xlabel('Volts RMS')
			axis.set_ylabel('# of channels (5181 total)')
			axis.set_title('Sensor Noise')
			axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig = plt.figure(1)
			axis = fig.add_subplot(111)
			axis.hist(filt1d, nbins)
			axis.set_xlabel('Volts')
			axis.set_ylabel('# of datapoints')
			axis.set_title('filtered signal values')
			axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()

	def tau_hist(self, nbins = 10000, begin = -0.03, end = 0.03, axis=0) :
		''' generate a histogram with the range of tau values we are dealing with.
		'''
		# use numpy.ravel() for 1d view of 2d array.
		# use numpy.flatten() for 1d copy of 2d array.
		filt1d = self.filt.ravel() 

		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(filt1d, nbins)
			axis.set_xlabel('Volts RMS')
			axis.set_ylabel('# of channels (5181 total)')
			axis.set_title('Sensor Noise')
			axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig = plt.figure(1)
			axis = fig.add_subplot(111)
			axis.hist(filt1d, nbins)
			axis.set_xlabel('Volts')
			axis.set_ylabel('# of datapoints')
			axis.set_title('filtered signal values')
			axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()

	def reform_data(self) :
		'''
		after analyzing the each Pixel's data and getting 
		filtered result, we can slam it all into one big 
		two-dimensional array for making 2d pixellated images.
		'''

		self.filt = np.array([self.pix[i].filt for i in range(self.nch)])

	def pixelate_multi(self, start, stop, stepsize, stepthru = False, vmin = 0, vmax = 0.007) :

		''' plot successive pixelated images of the 72*72 sensor. Set 'stepthru' to 'True'
		to step through each image by keyboard input.

		input:
		-data : a (5184 x number of time samples) numpy array.
		-start and stop : specify desired points in time to plot.
		-stepsize : decides how many of the time points to plot. if '1',
		all of the data is plotted. if '10' for example, each successive
		10th point in time is plotted. stepsize must divide into integers. 
		-stepthru : 
		'''

		a = np.arange(start, stop, stepsize)
		nplot = len(a)

		# get the data into (72 x 72 x npt) array
		data_2d = np.zeros((self.row, self.row, nplot))

		for i in range(nplot) :
			data_2d[:,:,i] = np.reshape(self.filt[:,a[i]], (self.row,-1)) # convert to square matrix

		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		ax.set_title("topmetal data")

		im = ax.imshow(np.zeros((self.row,self.row)), cmap=cm.jet, vmin = vmin, vmax = vmax)
		
		fig.show()
		im.axes.figure.canvas.draw()

		# quickly stream a series of images
		if not stepthru :

			#tstart = time.time()
			for j in range(nplot) :
				#t = j*stepsize*self.tsample
				ax.set_title("Frame %i" % a[j])
				#ax.set_title("Time elapsed: %f seconds" % t)
				im.set_data(data_2d[:,:,j])
				im.axes.figure.canvas.draw()

		# step through series of images by keyboard input 'enter'
		else : 

			for j in range(nplot) :
				input("press enter for next frame")
				#t = j*stepsize*self.tsample
				ax.set_title("Frame %i" % a[j])
				#ax.set_title("Time elapsed: %f seconds" % t)
				im.set_data(data_2d[:,:,j])
				im.axes.figure.canvas.draw()
		# close the figure
		plt.close(fig)

	def pixelate_single(self, sample, arr = [], vmin=-0.001, vmax=0.007, axis = 0):
		''' Take 5184 channel x 25890 data point array and plot desired points in
		time as 5184 pixel array.
		
		input:
		-sample : desired point in sample space to plot 0-25889
		-arr : input array of length (self.row)**2 to plot.
		can use this for quick tests.
		-vmin , vmax : minimum and maximum values for the colormap.
		-axis: supply an axis to impose this pixel plot to that axis.
		if no axis is supplied (default), then the function generates a standalone image.

		'''

		if isinstance(axis, ax_obj) :

			# default case, plot data by specifying sample in dataset.
			if not len(arr) :
				data_2d = np.reshape(self.filt[:,sample], (self.row, -1)) # convert to square matrix
				# make value bounds for the plot and specify the color map.
				im = axis.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			
			# if array is input, plot this image instead.
			else :
				data_2d = np.reshape(arr, (self.row, -1))
				im = axis.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			axis.grid(True)

		else :	
			fig, ax = plt.subplots()

			# default case, plot data by specifying sample in dataset.
			if not len(arr) :
				data_2d = np.reshape(self.filt[:,sample], (self.row, -1)) # convert to square matrix
				# make value bounds for the plot and specify the color map.
				im = ax.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			
			# if array is input, plot this image instead.
			else :
				data_2d = np.reshape(arr, (self.row, -1))
				im = ax.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
		
			fig.colorbar(im)
			ax.grid(True)
			fig.show()


	def plot_waveform(self, pixels, choose, avg=True, fit = False, lr = None) :
		'''
		This function can plot a series of channels. When the user
		presses 'Enter', the next channel's plot is replaces the prior.
		User can plot raw data, filtered data, or both with the 
		peak least squares fits.

		input: 
		-pixels : 1D np.array/list containing Pixel indices to plot. 
		-choose : 'd' plots raw data, 'f' plots filtered data.
			'b' plots the raw data and filtered data.
		-avg : impose a line of the average voltage onto the raw data plot.
		-fit : fit == 'True' superimposes scatter plot of fits 
			for a channel's peaks onto the same axis plotting the raw data.
		-lr : input a tuple to slice the plot. for example, if the dataset has
		25890 points, lr = (2000,4000) would only plot points 2000 through 4000.
		lr defaults to plotting the entire dataset. 
		'''
		if not lr :
			lr = (0, self.daq_length)

		x = np.arange(lr[0], lr[1])
		datalen = len(x)

		# setup
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_title("multich plot")

		ax.step(x, np.zeros(datalen))

		fig.show()
		ax.figure.canvas.draw()

		# plot raw data, optionally superimpose fits.
		if (choose == 'b') :
			
			fig2 = plt.figure()
			ax2 = fig2.add_subplot(111)
			ax2.step(x, np.zeros(datalen))
			fig2.show()
			ax2.figure.canvas.draw()

			print('plotting raw data and filtered data')
			for i in pixels :
				input("press enter for next channel")
				ax.cla()
				ax.set_title('raw data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))

				ax.step(x, self.pix[i].data[lr[0]:lr[1]])
				if avg :
					ax.step(x, np.ones(datalen)*self.pix[i].avg)


				if fit :
					for j in range(self.pix[i].npk) :
						ax.scatter(self.pix[i].peaks[j].fit_pts + self.pix[i].peaks[j].index, 
							model_func(self.pix[i].peaks[j].fit_pts, *self.pix[i].peaks[j].fit_par), 
							marker='o')
				ax2.cla()
				ax2.set_title('filtered data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))
				ax2.step(x, self.pix[i].filt[lr[0]:lr[1]])
				

				ax.figure.canvas.draw()
				ax2.figure.canvas.draw()

		if (choose == 'd') :
			print('plotting raw data')
			for i in pixels :
				input("press enter for next channel")
				ax.cla()
				ax.set_title('raw data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))
				ax.step(x, self.pix[i].data[lr[0]:lr[1]])
				if avg :
					ax.step(x, np.ones(datalen)*self.pix[i].avg)

				if fit :
					for j in range(self.pix[i].npk) :
						ax.scatter(self.pix[i].peaks[j].fit_pts + self.pix[i].peaks[j].index, 
							model_func(self.pix[i].peaks[j].fit_pts, *self.pix[i].peaks[j].fit_par), 
							marker='o', color = 'r')

				ax.figure.canvas.draw()

		# plot filtered data.
		elif (choose == 'f') :
			print('plotting filtered data')
			for i in pixels :
				input("press enter for next channel")
				ax.cla()
				ax.set_title('raw data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))
				ax.step(x, self.pix[i].filt[lr[0]:lr[1]])
				ax.figure.canvas.draw()

		else :
			print('unknown input: set "choose" = "f" for filtered data, "d" for raw data, or "b" for both.')

	def pixelate_tri_val(self) :
		''' provide 72X72 data array on channel status. 0 is a channel with no found peaks.
		1 is a channel exceeding the max number of peaks. 0.5 is 'normal' channels.
		'''

		
		data_2d = np.reshape(self.filt, (self.row, -1)) # convert to square matrix
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

	def pixel_status(self):


		''' Plot number of peaks found on each of the 5184 channels onto a 2D
		pixel array.
		'''

		# list comprehension of the number of peaks for each of the 5184 channels.
		pkperpix = np.array([self.pix[i].npk for i in range(self.nch)])

		data_2d = np.reshape(pkperpix, (self.row, -1)) # convert to square matrix
		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		ax.set_title("pixel_status")     

		mn = 0
		mx = np.amax(pkperpix)

		im = ax.imshow(data_2d, cmap=cm.hot, vmin=mn, vmax=mx)
		fig.show()
		im.axes.figure.canvas.draw()

	def prepare_openCV(self) :
		a = 'blah'
	
class Pixel(object):
	'''
	class for raw data and attributes of a pixel.
	class attributes: (for all instances)
	-row : dimension of pixel array. 
	'''

	# shared C library for trapezoidal filter
	shp_lib = CDLL("shaper.so")
	
	# shaper functions arguments. 
	# one for multiple M (takes peaks_t struct), another for single M dataset.
	shp_lib.shaper_multi.restype = None
						  # 	     	in 			out   	 length   peaks_t struct 	baseline
	shp_lib.shaper_multi.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t), c_double]

	shp_lib.shaper_single.restype = None
						  # 	     	in 			out   	  length  	  l  	   k 		 M 	   baseline
	shp_lib.shaper_single.argtypes = [double_ptr, double_ptr, c_ulong, c_ulong, c_ulong, c_double, c_double]

	# Pixel.daq_length class attribute is defined once 
	# the data length is found in the sensor class instance.
	row = 72
	

	# these attributes are set when Sensor.analyze() is called.
	# they pertain to the analysis chain.
	
	daq_length = None
	noisethresh = None

	threshold = None
	minsep = None
	
	fit_length = None
	fudge = None

	l = None
	k = None
	M_def = None
	shaper_offset = None

	def __init__(self, number, data, avg, rms):
		'''
		Pixel attributes :
		-number : linear pixel number from 0-5183
		-loc : (x,y)=(col,row) coordinate as a tuple. top-left corner is origin.
		 rightwards increasing x, downards increasing y.
		-data : numpy array of dataset
		-avg : average voltage
		-rms : rms voltage
		-peaks : a list of 'Peak' objects
		-offset : for shifting shaper transition points +/- -> fwd/bkwd.
		 useful because it's hard to get peak locations perfect, and because we want to
		 avoid transition on peaks' rise time, which won't account for in our fit model.
		-npk : total number of peaks
		-type : descriptor for pixel. some pixels are okay, others are noisy,
		3 are for demux and some show weird basleine behavior, may be do to alpha particle impacts, or just bad pixels.
		'''

		self.number = number
		self.loc = (number % Pixel.row, int(number/Pixel.row))
		self.data = data
		self.filt = np.empty_like(data)
		self.avg = avg
		self.rms = rms
		self.peaks = []
		self.tau = None
		self.chisq = None
		self.fit_par = None
		self.fit_pts = None

		self.npk = 0
		# if (self.rms > Pixel.noisethresh) :
		# 	self.type = 'noisy'
		# else :
		# 	self.type = None

	def peak_det(self, threshold = 0.006, minsep = 40, sgwin = 15, sgorder = 4, sign = 1):
		''' 
		input:
		- data : numpy array with peaks to look for.
		- mean: average value of 'data'.
		- threshold : minimum acceptable value (volts) above the average signal level for a peak.
		- minsep : minimum number of samples between peaks.  peaks closer than minsep are discarded.
		- sign : +1 for positive data with peaks, -1 for negative data with valleys.

		*** sgwin and sgorder are filter parameters for smoothing the data with a savitsky-golay filter
		before peak detection. derivative based peak detection on noisy signals requires smoothing. ***

		- sgwin : window of savitsky-golay filter. smaller window picks up smaller features
		and vice-versa. window is number of adjacent points for polynomial fitting.
		- sgorder : order of savitsky-golay polynomial for fitting to data.
		
		return:
		- pkm : return array of peak locations given as indices of the original input array 'data'.
		- trig(optional) : namedtuple for testing each step of peak detection if desired.
		'''	

		f = savgol_scipy(self.data, sgwin, sgorder) # smooth data
		kernel = [1, 0, -1]	# calculate each derivative
		df = convolve(f, kernel, 'valid') # note: each convolution cuts length of array by len(kernel)-1
		dfn = np.sign(df) # normalize derivative to 1. 1 increasing, -1 decreasing, 0 const.
		
		# the second derivative of the normalized derivative.
		# should only have non-zero values for peaks and valleys, where the value of the derivative changes.
		ddfn = convolve(dfn, kernel, 'valid') 

		# first, find all of the positive derivative values. going up the peak.
		# this returns indices of possible candidates. we want to offset by two because
		# the convolution cuts the array length by len(kernel)-1
		if (sign == 1) :
			candidates = np.where(dfn > 0)[0] + (len(kernel)-1)
		elif (sign == -1) :
			candidates = np.where(dfn < 0)[0] + (len(kernel)-1)

		pk = sorted(set(candidates).intersection(np.where(ddfn == -sign*2)[0] + 1))
		alpha = self.avg + (sign * Pixel.threshold)

		if (sign == 1) :	# apply threshold to the raw data, not the smoothed data 
			pk = np.array(pk)[self.data[pk] > alpha]
		elif (sign == -1) :
			pk = np.array(pk)[self.data[pk] < alpha]

		# list comprehension version of toss np.array for minimum separation discrimination.
		# list is faster for 'append' operations, especially when there are many false peaks.
		# 'toss' gives the indices of the bad peaks in the array of peak locations, 
		# not the indices of the peaks in the dataset.
		toss = [i+1 for i in range(len(pk)-1) if ((pk[i+1]-pk[i]) < minsep)]
		
		# 'badpk' would give the location of thrown out peaks in the dataset, incase we want this later.
		# badpk = [pk[i+1] for i in range(len(pk)-1) if ((pk[i+1]-pk[i]) < minsep)]

		# remove peaks within the minimum separation... can do this smarter.
		#toss = np.array([])
		#for i in range(len(pk)-1) :
			# if the peaks are closer than the minimum separation and the second peak is
			# larger than the first, throw out the first peak. 
			#if ((pk[i+1]-pk[i]) < minsep) :
			#if ((peaks[i+1]-peaks[i]) < minsep) and (Y[peaks[i+1]] < Y[peaks[i]]) :
				#toss = np.append(toss, i+1)

		# pkm = peaks with minimum separation discrimination.
		# pk = peaks without minimum separation discrimination.
		pkm = np.delete(pk, toss)

		# cons = 5 # consecutively increasing values preceeding a peak
		# for j in range(len(pkm))
		# 	for k in range(cons)

		# use the 'trig' namedtuple for debugging / accessing each step of the peak detection.
		#return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=pk, toss=toss, pkm=pkm)
		
		# create a list of 'Peak' objects with these peak locations.
		self.prepeak = pk
		self.toss = toss
		self.peaks = [Peak(i) for i in pkm]
		self.npk = len(self.peaks)

	def fit_pulses(self):
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

		return: 
		- tau : array of fitted tau values. same length as number of input peaks.
		'''
		# only try to fit pulses if the list of peaks contains any peaks.
		# a list with items returns 'True'. An empty list returns 'False'.
		if self.peaks :

			# if the last peak is too close to the end of the dataset, we won't try to fit it.
			if (self.daq_length - self.peaks[self.npk-1].index < Pixel.fit_length + Pixel.fudge) :
				self.peaks.pop()
				self.npk -= 1
			# assumption: data is modeled by a 3 parameter decaying exponential.
			M = 3

			# 			   PARAMETER FORMAT: ([A, l, off]
			#				 MIN 								MAX
			bounds = ([0.0, 1.0/100, self.avg-0.3], [0.03, 1.0/10, self.avg+0.3])
			guess = [0.008, 1.0/35, self.avg]

			# include a 'fake' peak to act as endpoint for last iteration.
			self.peaks.append(Peak(Pixel.daq_length-1))

			for j in range(self.npk) :
				# if peaks are closer than fit_length, fit the space between them minus 'fudge'.
				if self.peaks[j+1].index-self.peaks[j].index < Pixel.fit_length : 
					yi = self.data[self.peaks[j].index:self.peaks[j+1].index-Pixel.fudge]
				# if peaks are farther than fit_length apart, use fit_length.
				else :								  
					yi = self.data[self.peaks[j].index:self.peaks[j].index+Pixel.fit_length]
				
				N=len(yi)
				xi = np.arange(N)

				# use this if testing with fake data / no noise.
				# if self.rms == 0 : # for fitting ideal/fake data without noise
				# 	par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
				# 		              check_finite=False, bounds=bounds, method='trf')
				# else : # for real data
					# par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
					# 	             sigma=np.ones(N)*self.rms, absolute_sigma=True, check_finite=False, \
					# 	             bounds=bounds, method='trf')

				par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=guess, \
					             sigma=np.ones(N)*self.rms, absolute_sigma=True, check_finite=False, \
					             bounds=bounds, method='trf')

				# par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, \
				# 		             sigma=np.ones(N)*rms, absolute_sigma=True, check_finite=False, \
				# 		             bounds=bounds, method='trf')


				f_xi = model_func(xi, *par)
				chisq, pval = chisquare(f_obs=yi, f_exp=f_xi, ddof=N-M)
				
				# each Pixel object has a list of peak objects,
				# determined by the number of peaks found in peak detection.
				# insert the fitted tau and chisquare to their respective peak.
				self.peaks[j].tau = 1.0/par[1]
				self.peaks[j].chisq = chisq
				
				# uncomment these attributes to do plots of fits.
				# self.peaks[j].fit_par = par
				# self.peaks[j].fit_pts = xi
				
				# Q and p values.
				# P[j] = pval
				# Q[j] = 1-pval

				# if axis object is provided, add the fit as a scatter plot to the axis.
				# if isinstance(ax, ax_obj) :
				# 	ax.scatter(xi+peaks[j], model_func(xi,*par), marker = 'o')

			# remove the 'fake' peak now that we're done iterating over
			# the fit procedure.

			self.peaks.pop()

	def filter_peaks(self):

		''' apply trapezoidal filter to data. peak locations are used to change
		filter parameters. '''

		# no peaks found, use default M. 
		
		# alternatively: to do a simple analysis do NOT run
		# peak_det or fit_pulses. this leaves npk initialized to 0 for every channel, so we just use default
		# filtering for every channel.
		if (self.npk == 0) : # no peaks found, use default M.
			Pixel.shp_lib.shaper_single(self.data, self.filt, c_ulong(Pixel.daq_length), 
				c_ulong(Pixel.l), c_ulong(Pixel.k), c_double(Pixel.M_def), c_double(self.avg))

		elif (self.npk == 1) : # one peak found.
			Pixel.shp_lib.shaper_single(self.data, self.filt, c_ulong(Pixel.daq_length), 
				c_ulong(Pixel.l), c_ulong(Pixel.k), c_double(self.peaks[0].tau), c_double(self.avg))

		elif (self.npk > 1) : # multiple peaks found.
			l_arr = np.ones(self.npk, dtype = c_ulong)*Pixel.l
			k_arr = np.ones(self.npk, dtype = c_ulong)*Pixel.k
			LR = self.pk2LR()

			# print('l', l)
			# print('k', k)
			# print('M', M)
			# print('number of peaks =', npk)
			# print('LEFT = ', LR[0])
			# print('RIGHT = ', LR[1])
			M = np.array([self.peaks[i].tau for i in range(self.npk)])
			# peaks_handle prepares a Ctypes 'Structure' object to make a C Structure.
			PEAK = peaks_handle(self.npk, LR[0], LR[1], l_arr, k_arr, M)
			# self.filt = np.empty_like(self.data)
			Pixel.shp_lib.shaper_multi(self.data, self.filt, c_ulong(len(self.data)), byref(PEAK), c_double(self.avg))

		self.filt = np.array(self.filt)

	def pk2LR(self) :
	
		''' put peak locations into arrays of left and rights for trapezoidal shaper.
		 apply desired offset. Both arrays are the same length as input 'peaks' array.

		return:
		- LEFT and RIGHT : beginning and end points for each set of 
		trapezoidal filter parameters l,k, and M.'''
	
		LEFT = np.zeros(self.npk)
		RIGHT = np.zeros(self.npk)

		for i in range(self.npk-1):
			LEFT[i]  = self.peaks[i].index   + self.shaper_offset
			RIGHT[i] = self.peaks[i+1].index + self.shaper_offset
			
		LEFT[0] = 0
		LEFT[self.npk-1] = self.peaks[self.npk-1].index + self.shaper_offset
		RIGHT[self.npk-1] = self.daq_length

		# trapezoidal filter uses size_t, or c_ulong, as its datatype
		# for left and right. they index locations possibly larger than int allows.
		LEFT = np.array(LEFT, dtype = c_ulong)
		RIGHT = np.array(RIGHT, dtype = c_ulong)

		return (LEFT, RIGHT)

	def plot(self, choose, axis = 0, lr = None) :
		''' 
		1 dimensional 'step' plot.
		inputs:
		- axis : supply a pyplot axis object to plot to. For use as a quick and dirty
		 plot in ipython, supply 0 for axis to generate a standalone fig, axis plot.
		- lr : tuple to select portion of dataset to plot. default is entire dataset.
		- choose : string to select 'data' or 'filt'
		'''
		if (lr == None) :
			lr = (0, self.daq_length-1)

		xaxis = np.arange(self.daq_length)[lr[0]:lr[1]]
		if (choose == 'd') :
			d = self.data[lr[0]:lr[1]]

		elif (choose == 'f') :
			d = self.filt[lr[0]:lr[1]]

		else :
			print('unknown input: choose "f" for filtered data, "d" for raw data.')
			d = 0

		if isinstance(axis, ax_obj) :	# if an axis is supplied, plot to it.
			axis.step(xaxis, d)

		else :	# no axis, make a quick standalone plot.
			plt.step(xaxis, d)
			plt.show()

	def text_dump(self, select, values) :
		''' output data values over channels
		onto the terminal.
		'''
		self.blob = 'fff'

	def pixel_type(self):
			self.magabob = 0

class Peak(object):

	def __init__(self, index):
		self.index = index
		# self.tau = None
		# self.chisq = None
		

class Event(object):

	def __init__(self, x, y, radius, index):
		self.radius = radius
		self.x = x
		self.y = y
		self.index = index

def get_wfm_one(infile, ch, npt, plt) :

	''' get data for a channel and calculate its average and root mean square. 
	npt defaults to max value (25890 for current dataset) for zero input or too large of input.'''
	data = pull_one(infile, ch)

	avg = np.mean(data[:npt])
	rms = np.std(data[:npt])

	if plt == True :
		print('average = ', avg, 'Volts')
		print('sigma = ', rms, 'Volts RMS')
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
		print('set calculation length to raw data array length =', length)
		npt = length
	for i in range(dead,nch) : # leave channels 0-2 with avg = 0, rms = 0.
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
	for i in range(len(pk)-1) :
		# if the peaks are closer than the minimum separation and the second peak is
		# larger than the first, throw out the first peak. 
		if ((pk[i+1]-pk[i]) < minsep) :
		#if ((peaks[i+1]-peaks[i]) < minsep) and (Y[peaks[i+1]] < Y[peaks[i]]) :
			toss = np.append(toss, i+1)

	pkm = np.delete(pk, toss)

	# cons = 5 # consecutively increasing values preceeding a peak
	# for j in range(len(pkm))
	# 	for k in range(cons)

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

	for j in range(npk) :
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

	# print('l', l)
	# print('k', k)
	# print('M', M)
	# print('number of peaks =', npk)
	# print('LEFT = ', LR[0])
	# print('RIGHT = ', LR[1])

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

	for i in range(npk-1):
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

	for i in range(nch) :
		input("press enter for next channel")
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
		axis.set_xlim(0,end) # x limits, y limits
		#axis.set_ylim()
		axis.grid(True)

	else : 			# no axis supplied, make standalone plot.
		fig = plt.figure(1)
		axis = fig.add_subplot(111)
		axis.set_title("Sensor Noise Histogram")
		axis.hist(data, nbins)
		axis.set_xlabel('Volts RMS')
		axis.set_ylabel('# of channels (5181 total)')
		axis.set_title('Sensor Noise')
		axis.set_xlim(0, end) # x limits, y limits
		#axis.set_ylim()
		axis.grid(True)
		fig.show()

def define_square(x, y, len):

	''' define a square region of interest on the 2D pixel array.
	input:
	- x, y : cartesian coordinates of upper left corner of square of interest.
	- len : length of side of the square
	
	return:
	- array of channels inside of the square region of interest.
	'''
	row = 72 # sensor is a 72x72 square
	
	return roi
def locate_pixel(x, y) :
	''' Take x and y location of square pixel array and convert to location in linear 1D array
	input:
	- x, y : cartesian coordinates of pixel in square array.
	'''
	dim = 72
	out=(y*72)+x
	print(out)
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
	if desired. 
	'''

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

	for i in range(npt) :
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
	for j in range(npt) :
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

	for i in range(len(data[channel])):
	   print(i, data[channel][i])


def demuxD(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_dbl.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def demuxF(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux_flt.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))


