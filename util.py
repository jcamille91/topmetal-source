# some useful functions

# test for slow pieces of code
import time
# print out python data objects (class instances, lists, etc)
import pprint
# for importing data from HDF5 files into Numpy arrays
import h5py 

# we use 'pickle' to save a big list of detected peaks from the analysis chain.
# we can reload them later to quickly re-filter the data without needing to
# redo the detection and fit procedure.
import pickle

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

# Ctypes and Numpy support for calling C functions defined in shared libraries.
# These functions can be slow when implemented in Python.

# datatypes ... C equivalent
from ctypes import (
	c_ulong,	# size_t
	c_int,		# int
	c_double, 	# double
	c_float		# float
	)



from ctypes import (
	CDLL,		# CDLL used to load desired shared library. 
	POINTER,	# POINTER creates equivalent of a C pointer.
	byref,		# byref() for passing pointers to function.
	Structure 	# ctypes 'Structure' class allows creation of C 'struct'.
	)

import numpy.ctypeslib as npct

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

	# let's get all of the peaks from each pixel.


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
		it may not be useful anyways because the savitsky-golay filter and trapezoidal filter will
		give different results if used on the whole dataset versus only a chunk of it.
		'''

		with h5py.File(self.infile,'r') as hf: # open file for read, then close
			d = hf.get(self.daq_event)
			data = np.array(d, dtype=np.float64)
		self.cut = cut
		
		if not cut: # do all of the data

			self.daq_length = len(data[0])
			Sensor.daq_length = len(data[0])

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


	def make_pixel(self, choose, mV_noise = 1):
		'''
		use this function to generate data for the analysis chain.

		input:
		-mV_noise : desired RMS noise voltage, in units of millivolts. 
		for reference, the typical topmetal- II- pixel  has about 2mV rms of noise.
		-choose : string for type of data to send. 'exp' for exponential pulses,
		'step' for step input. can add more options if desired.

		if using data without noise, be sure to change the RMS input in the fitting procedure,
		or  the ls fitting function will complain. This option has been commented out for later use
		in the Pixel.fit_tau() method.

		### to analyze this data, enter the created pixel index into
		### the 'select' list argument for the Sensor.analyze() function.


		'''
		# create list to store data. normally 5184 channels are loaded into this list,
		# but one or any number is fine for experimenting with fake data.
		self.pix = []

		# define number of points for dataset.
		# make sure there are enough points to accomodate 
		# whatever features are defined below.
		self.daq_length = 100000
		
		# number of points defining each pulse.
		M_npt = 1000

		# signal baseline and noise.
		baseline = 0.8

		if mV_noise == 0 :
			noise = 0
		else :
			noise = np.random.normal(loc = 0, scale = mV_noise*.001, size=self.daq_length)


		# specify M=tau, peak location, amplitude, and the number of points to define
		# the peak.
		
		# M_i = np.array([40], dtype=np.float64)
		# M_loc = np.array([5000])
		# M_a = np.array([0.015])

		if choose == 'exp' :

			# specify M=tau, peak location, amplitude, and the number of points to define the peak.
			
			#### one pulse
			# M_i = np.array([40], dtype=np.float64)
			# M_loc = np.array([50000])
			# M_a = np.array([1])

			#### a few pulses with different parameters
			M_i = np.array([70, 70, 60, 5], dtype=np.float64)
			M_loc = np.array([1000, 3000, 5000, 8000])
			# M_a = np.array([0.02, 0.015, 0.011, 0.015])
			M_a = np.array([0.008, 0.008, 0.008, 0.008])

			#### variable number (n) of pulses with same tau and amplitude
			# n = 20
			# M_i = np.ones(n, dtype=np.float64)*20
			# M_loc = np.arange(start=M_npt, stop=self.daq_length-2*M_npt, step = int((self.daq_length-4*M_npt)/n))
			# M_a = np.ones(n)*0.008

			#### three stacked pulses, 
			####let's see if the vsum calibration works if the response is no longer trapezoidal.
			# M_i = np.array([40, 40, 40, 40], dtype=np.float64)
			# M_loc = np.array([10000, 10040, 50000, 80000])
			# M_a = np.array([0.02, 0.015, 0.011, 0.015])


			npk = len(M_i)

			# build an array with exponentially decaying pulses just defined above.
			exp = np.ones(self.daq_length, dtype=np.float64)*baseline
			for i in range(npk) :
				exp[M_loc[i]:M_loc[i]+M_npt] += M_a[i]*np.exp(-(1.0/M_i[i])*np.linspace(0, M_npt-1, M_npt)) 

			signal = exp + noise

		elif choose == 'step' :

			height = 0.015
			step = np.ones(self.daq_length)*baseline
			step[int(self.daq_length/2):] += height # transition point
			
			signal = step + noise

		elif choose == 'trap' :

			# don't think this is correct
			trap = np.zeros(self.daq_length)
			A = 4.8
			l=60
			k=20
			loc = 3000
			stepsize = 4.8/20
			for i in range(k):
				trap[loc+i] = i*stepsize
			trap[loc+k:loc+l-1]=A
			for i in range(k):
				trap[loc+l-1+i] = A-stepsize*i
			signal = trap + noise

		elif choose == 'noise' :
			signal = baseline + noise


		else :
			print("unrecognized input, set 'choose' to either 'exp' or 'step' to generate the desired signal")

		self.pix.append(Pixel(0, signal, baseline, mV_noise*0.001))

	def analyze(self, read=False, simple = False, select = [],
				
				sign = 1, minsep = 50, threshold = 0.006, 								   # peak detection						
				sgwin = 15, sgorder = 4, l_s=4, k_s=2, choose='sg',		
				pd_save = 0,  			 
			    
			    fit_length = 300, fudge = 20, Mrange = [10.0,45.0],							   # lsq fit
			    bounds = ([0.0, 1.0/300, 0-0.01], [0.03, 1.0, 0+0.01]),	   			   	   # format: (min/max)
			    guess = [0.008, 1.0/20, 0],						   				   		   # amplitude, 1/tau, offset
			    
			    l = 60, k = 20, M_def = float(30), shaper_offset = -20):  				   # shaper

		''' analysis chain: 
		1. find peaks 
		2. fit peaks 
		3. apply trapezoidal shaper with fitted tau to dataset
		
		input:
		-simple : set 'True' to just filter all the channels with the default filter parameters.
		 no peak detection / exponential decay fitting.
		-select : provide list of desired channels. provide empty list to analyze all channels.
		-step : set 'True' to run filter pixel(s) with containing step input. 
		### not currently setup to work with 'step'
		-noisethresh : threshold in volts to classify channels. channels above this value are labeled with
		a string  as 'noisy'.

		(PEAK DETECTION)

		-sign : this allows peak detection algorithm to handle '+' or '-' data. 
		postitive data with peaks -> sign = 1
		negative data with valleys -> sign = -1
		-minsep : minimum separation from peaks. if two are more peaks are within minsep,
		only the first peak candidate is saved.
		-threshold : minimum value for accepted peaks. candidate peaks
		with max value smaller than the threshold are rejected.
		-l_s, k_s : these are filter parameters for the trapezoidal filter, if it's to be used for peak detection. 
		It's meant generally for optimizing SNR for signal measurement and zeroing the pixels, but it may have potential for peak detection.
		Here, the trapezoidal filter should be use small l-k and k. the skinny flat top zeros the signal well and most importantly,
		can identify constituent peaks of a larger apparent peak while also making it possible to accurately find the
		very beginning of the signal rise.
		-pd_choose : choose 'sg' savitsky-golay or 'skinny' trapezoidal filter to smooth raw data for peak det.


		*** sgwin and sgorder are filter parameters for smoothing the data with a savitsky-golay filter
		before peak detection. derivative based peak detection on noisy signals requires smoothing. ***

		-sgwin : window of savitsky-golay filter. smaller window picks up smaller features
		and vice-versa. window is number of adjacent points for polynomial fitting.
		-sgorder : order of savitsky-golay polynomial for fitting to data.
	


		(TAU FITTING)

		-fit_length : number of datapoints to fit. generally (2-3) expected tau in length.
		-fudge : if consecutive peaks are closer than fit_length, the procedure fits the 
		 smaller space between the two consecutive peaks. In this case, 'fudge' many points 
		 are removed from the end of the fit, to try and ensure the tail of the fit 
		 does not catch the rise of the next peak.
		-bounds : minimum and maximum values as possible fit values. specifying a reasonable range of
		values can speed up the fit process. the baseline for each pixel offsets the amplitude and offset fit parameters. 
		FORMAT: [(amplitude, 1/tau, offset)min, (amplitude, 1/tau, offset)max]
		-guess: inital guess for exponential fit parameters. [amplitude, 1/tau, offset]
 		-Mrange : this specifies how to split up three categories for detected peaks. lower than the
 		first argument is a noise peak. longer than the second argument is multiple peaks in close proximity.
		
		(TRAPEZOIDAL FILTER)

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

		 ### note: shaper_offset has no effect on data with 0 or 1 peaks. shaper offset is only applied to datasets
		 with multiple peaks. this is reflected in Pixel.filter_peaks() logic.

		'''
		# these 'Pixel' class attributes are used by methods for analysis.
		# we're setting the namespace (the 'Pixel' class) for all of the analysis functions here.
		# main three 'Pixel' methods are : peak_det(), fit_pulses(), filter_peaks().
		Pixel.daq_length = self.daq_length 

		Pixel.sign = sign
		Pixel.threshold = threshold
		Pixel.minsep = minsep
		Pixel.sgwin = sgwin
		Pixel.sgorder = sgorder
		Pixel.l_s = l_s
		Pixel.k_s = k_s
		Pixel.pd_choose = choose
		Pixel.pd_save = pd_save

		Pixel.fudge = fudge
		Pixel.fit_length = fit_length
		Pixel.Mrange = Mrange
		# these get offset by each pixel's average voltage (baseline) when the pixel object is initialized.
		Pixel.fit_bounds = np.array(bounds)
		Pixel.fit_guess = np.array(guess)

		Pixel.l = l # fixing l and k shaper parameters
		Pixel.k = k
		Pixel.M_def = M_def # default M is no peaks are found.
		Pixel.shaper_offset = shaper_offset # this should generally be negative

		# define a peak object for a 'no peak' channel just once, so it can be reused.
		l_arr = np.ones(1, dtype = c_ulong)*Pixel.l
		k_arr = np.ones(1, dtype = c_ulong)*Pixel.k	
		LEFT = np.array(0, dtype = c_ulong)
		RIGHT = np.array(self.daq_length, dtype = c_ulong)		
		M = np.array(M_def)
		Pixel.zero_pk = peaks_handle(1, LEFT, RIGHT, l_arr, k_arr, M)

		# define a peak object for 'skinny' trapezoidal filter parameters. this is used for peak detection.
		l_arr2 = np.ones(1, dtype = c_ulong)*Pixel.l_s
		k_arr2 = np.ones(1, dtype = c_ulong)*Pixel.k_s
		LEFT2 = np.array(0, dtype = c_ulong)
		RIGHT2 = np.array(self.daq_length, dtype = c_ulong)		
		M2 = np.array(M_def)
		Pixel.skinny_pk = peaks_handle(1, LEFT2, RIGHT2, l_arr2, k_arr2, M2)

		print('analysis parameters... \n')
		print('peak detection: sign = %i, threshold = %f, minsep = %i, \n sgwin = %i, sgorder = %i. \n' \
			% (sign, threshold, minsep, sgwin, sgorder))
		print('pulse fitting: fudge = %i, fit length = %i \n'  % (fudge, fit_length))
		print('shaping filter: l = %i, k = %i, \n default M = %i, offset = %i \n' % (l, k, M_def, shaper_offset))
		input('press enter to analyze with given parameters')


		if read : # if we have read in a saved list of peaks for each channel, 
				  # skip peak detection and fitting procedure.
			for i in self.pix :
				i.filter_peaks()

		elif simple : # do simple analysis, filter everything with default l,k,M.
			for i in self.pix :
				i.filter_peaks(simple=True)

		elif not select : # if list is empty, analyze all pixels
			
			# remove 'BAD' channels from regular analysis
			ch = [i for i in range(self.nch)]
			ch = list(set(ch) - set(Sensor.BAD))

			# regular analysis
			for i in ch:
				#print(self.pix[i].__dict__)
				print('channel no. %i' % i)
				self.pix[i].peak_det()
				print("pd done")
				#self.pix[i].__dict__

				self.pix[i].fit_pulses()
				print('fit done')
				print(self.pix[i].__dict__)
				self.pix[i].filter_peaks()
				print('filter done')


			# just default filter bad channels
			for i in Sensor.BAD :
				print('bad channel no. %i' % i)
				self.pix[i].filter_peaks()

		else :
			for i in select : # analyze only pixels in the 'select' list.
				self.pix[i].peak_det()
				self.pix[i].fit_pulses()
				self.pix[i].filter_peaks()

		
		print('finished analyzing', len(self.pix), 'pixels.') 


	def label_events(self):

		# just a quick and dirty function to manually input events of
		# interest for future reference 
		# all taken from out22.h5

		# Event methods need row to retrieve pixel selections.
		Event.row = self.row

		# the length of an event should be ~ l+k, determined by the trapezoidal
		# filter parameters.

		# circle shaped events; use select_circle(x, y, radius, index)
		self.alpha_events = [

		# events defined by circles

		Event(x=10, y=29, r=10, i=3520, shape='c'), #
		# 4800, next 320, 470, 500, 900
		Event(x=15, y=9, r=9, i=6510, shape='c'), #
		Event(x=45, y=31, r=10, i=6580, shape='c'), #
		Event(x=11, y=30, r=10, i=6845, shape='c'), #
		Event(x=18, y=45, r=15, i=7020, shape='c'), #
		Event(x=25, y=46, r=13, i=7505, shape='c'), #
		Event(x=18, y=45, r=15, i=7020, shape='c'), #
		Event(x=23, y=10, r=10, i=7865, shape='c'), #
		Event(x=19, y=10, r=10, i=8510, shape='c'), #
		Event(x=43, y=29, r=12, i=8910, shape='c'), #
		Event(x=28, y=16, r=11, i=9095, shape='c'), #
		# CHECKME Event(x=18, y=28, r=13, i=9545, shape='c'), # this one moves a little bit but is focused well
		Event(x=36, y=29, r=12, i=10610, shape='c'), #
		Event(x=32, y=10, r=10, i=10910, shape='c'), #
		Event(x=39, y=14, r=12, i=11030, shape='c'), # lopsided but focused, rectangle is better here
		Event(x=42, y=17, r=12, i=12460, shape='c'), #
		Event(x=22, y=48, r=12, i=12850, shape='c'), #
		Event(x=9, y=28, r=9, i=13060, shape='c'), #
		Event(x=42, y=17, r=12, i=12460, shape='c'), #
		Event(x=18, y=28, r=10, i=13700, shape='c'), #
		Event(x=20, y=9, r=9, i=14425, shape='c'), # elliptical and focused
		Event(x=28, y=9, r=9, i=15140, shape='c'), #
		Event(x=18, y=23, r=9, i=15475, shape='c'), #
		Event(x=40, y=42, r=14, i=16825, shape='c'), # lots of dark pixels here
		# CHECKME Event(x=12, y=7, r=7, i=17600, shape='c'), # elliptical and a little broken up
		Event(x=40, y=14, r=10, i=18345, shape='c'), # nice circle
		Event(x=42, y=30, r=10, i=18755, shape='c'), # nice circle but not isolated
		Event(x=41, y=32, r=9, i=20370, shape='c'), # moves a bit and is a little bit separated
		Event(x=40, y=30, r=11, i=20550, shape='c'), # fat elliptical circle
		# CHECK ME Event(x=24, y=42, r=11, 21140, shape='c'), # circle, lots of darker pixels
		# CHECK ME Event(x=24, y=42, r=11, 21375, shape='c'), # pretty similar in location/consistency to previous event
		Event(x=40, y=52, r=10, i=21830, shape='c'), # small ellipse
		Event(x=20, y=11, r=10, i=22270, shape='c'), # starts circular becomes elliptical
		Event(x=40, y=35, r=10, i=22675, shape='c'), # big messy circle
		Event(x=12, y=8, r=8, i=23060, shape='c'), # elliptical blip
		Event(x=7, y=25, r=7, i=23220, shape='c'), # two disjointed blobs, probably can't capture both with a circle
		Event(x=16, y=7, r=7, i=23400, shape='c'), #
		Event(x=40, y=12, r=10, i=23545, shape='c'), # circle, but not totally isolated from other signal
		Event(x=31, y=13, r=12, i=23705, shape='c'), #
		#CHECK ME Event(x=26, y=6, r=6, i=23775, shape='c'), # very small circle, not isolated
		Event(x=25, y=54, r=13, i=23890, shape='c'), # kind of scattered, not a closed circle
		Event(x=12, y=30, r=12, i=24175, shape='c'), # scattered circle, in happens very close in time to two other events
		Event(x=13, y=29, r=12, i=24940, shape='c'), # big elliptical blob, focused


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

	def save_peaks(self, outfile) :

		# get all lists of peak objects as a list with the same length as the number of channels 
		peaks = [pix.peaks for pix in self.pix]

		with open(outfile, 'wb') as output:
			pickle.dump(peaks, output)

	def read_peaks(self, infile) :

		# retrieve list of each pixel's list of peak objects
		with open(infile, 'rb') as read:
			pk =  pickle.load(read)

		if len(pk) == self.nch :

			for i in range(self.nch):
				self.pix[i].peaks = pk[i]

		else :
			print('number of lists of peaks:', len(peaks), 'does not match number of pixels:', self.nch)
			print('cannot load peaks if these numbers do not match')
		
	def vsum_all(self, nbins = 100, start=0, stop=None, use_events=True, nframe=None, axis=0):
		'''
		function to do voltage summations over event windows of all 5184 pixels.

		input:
		-start : location in dataset to start, defaults to 0.
		-nframe : number of frames to sum over the filtered data. defaults 
		to length of the trapezoidal filter response (l+k)

		''' 

		if not stop:
			stop = self.daq_length


		# use default l+k as the window, aka the filter length.
		if use_events :
			self.vsum=[]
			for ev in self.alpha_events : 
				values = np.sum(self.filt[:,ev.i:ev.i+Pixel.l+Pixel.k])
				self.vsum.append(values)


		elif not nframe :
			nframe = Pixel.l + Pixel.k
			nevent = int((stop-start)/(nframe))
			self.vsum = np.zeros(nevent)

			for e in range(nevent) :
				i = start + nframe*e
				values = np.sum(self.filt[:,i:i+nframe]) # all 5184 pixels' data points inside the event window.
				self.vsum[e] = np.sum(values)

		else :
			nevent = int((stop-start)/(nframe))
			self.vsum = np.zeros(nevent)

			for e in range(nevent) :
				i = start + nframe*e
				values = np.sum(self.filt[:,i:i+nframe]) # all 5184 pixels' data points inside the event window.
				self.vsum[e] = np.sum(values)

		print(len(self.alpha_events),'events, each summed over', Pixel.l+Pixel.k, 'frames')

		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(x=self.vsum, bins=nbins)
			axis.set_xlabel('Volts')
			axis.set_ylabel('counts')
			axis.set_title('alpha energy spectrum')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig, axis = plt.subplots(1,1)
			axis.hist(x=self.vsum, bins=nbins)
			axis.set_xlabel('Volts')
			axis.set_ylabel('counts')
			axis.set_title('alpha energy spectrum')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)
			fig.show()	 


	def vsum_select(self, show = False, nfake=100, v_window = (0,0), hist_lr = [(0, 300), (-100,100)], nbins=[20,20], axis = 0):

		'''a word on this measurement: for the experimental setup, each alpha event should deposit all of its
		energy into ionizing air. nearly all charge due to this event should be picked up by the sensor. Therefore,
		if we add all the consituent voltages of an event, we should get a measure of the alpha particle's energy
		spectrum (a "sharp" peak).

		This function makes takes selections of pixels within circles
		to define events. then, over all frames constituting an event, voltage is summed.
		
		input:
		-show : if True, step through a pixel image of each event with a circle enclosing the region of interest.
		-v_window : defaults to calculating voltage sum with l+k samples, starting at the event index. 
		+/- values move the (left,right) fwd/bkwd
		-hist_lr : sets the (left,right) bounds for the histogram
		-nbins : number of bins for histogram
		-axis : supply axis to plot to, or don't supply an axis and generate standalone plot.
		'''

		ring = [] # contains info for each as a tuple so we can easily print later.
		self.alphaE = [] # contains voltage summation for a single event.

		# get location for every event. get selection of pixels based on this location for 
		# each event. store the voltage summation for histogram.

		for ev in self.alpha_events :
			
			# retrieve a list of pixels defining the region of interest for this event.
			ev.retrieve_selection()

			# we'll use generator-expressions over list comprehensions here since we
			# don't need to store all the constituent values to be summed.
			vsum = sum(self.pix[i].filt[j] for i in ev.sel for j in range(ev.i + v_window[0], ev.i+Pixel.l+Pixel.k+v_window[1]))

			#vsum = np.sum(np.array([self.pix[i].filt[j] for i in circle for j in range(frames[0], frames[1])]))
			#valcheck = np.array([self.pix[i].filt[j] for i in circle for j in range(frames[0], frames[1])])
			
			self.alphaE.append(vsum)
			ring.append((ev.i+(Pixel.l+Pixel.k)/2, ev.x, ev.y, ev.r, vsum))


		# create some fake/random selections to observe the zero 'noise' peak in the energy spectrum.
		fake_ring  = []
		self.fake_ev=[]
		self.fake_E =[]
		for k in range(nfake) :
			x,y = np.random.randint(10,61,size=2)
			i = np.random.randint(0,self.daq_length-(Pixel.l+Pixel.l))
			self.fake_ev.append(Event(x=x, y=y, r=10, i=i, shape='c'))

					

		for f in self.fake_ev :
			
			f.retrieve_selection()
			vsum = sum(self.pix[i].filt[j] for i in f.sel for j in range(f.i + v_window[0], f.i+Pixel.l+Pixel.k+v_window[1]))			
			self.fake_E.append(vsum)
			fake_ring.append((f.i+(Pixel.l+Pixel.k)/2, f.x, f.y, f.r, vsum))



		# plot histograms for energy spectrum

		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(x=self.alphaE, bins=nbins, range=hist_lr[0])
			axis.set_xlabel('Volts, summed over event pixels and frames')
			axis.set_ylabel('counts')
			axis.set_title('alpha energy spectrum')
			#axis.set_xlim(begin, end) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig, ax = plt.subplots(1,1)
			ax.hist(x=self.alphaE, bins=nbins[0], range=hist_lr[0])
			ax.set_xlabel('Volts, summed over event pixels and frames')
			ax.set_ylabel('counts')
			ax.set_title('alpha energy peak')
			#ax.set_xlim(begin, end) # x limits, y limits
			#ax.set_ylim()
			ax.grid(True)
			fig.show()

			fig2, ax2 = plt.subplots(1,1)
			ax2.hist(x=self.fake_E, bins=nbins[1], range=hist_lr[1])
			ax2.set_xlabel('Volts, summed over event pixels and frames')
			ax2.set_ylabel('counts')
			ax2.set_title('sensor noise "zero peak"')
			#ax2.set_xlim(begin, end) # x limits, y limits
			#ax2.set_ylim()
			ax2.grid(True)
			fig2.show()		 

		# show the the middle most frame of events of interest.
		if show :

			fig3, ax3 = plt.subplots(1,1)

			# p is a tuple (midpoint, x, y, r)
			for p in ring:
				input("press enter to show next event:")	
				ax3.cla()
				ax3.set_title('frame=%i coordinate=(%i, %i) radius=%i vsum = %fV' % p)
				self.pixelate_single(sample = int(p[0]), arr=[], axis = ax3)
				# add a circle 'artist' to the event we are analyzing
				circ = plt.Circle((p[1], p[2]), p[3], color = 'r', fill=False, linewidth = 1.5, alpha=1)
				ax3.add_artist(circ)
				fig3.show()

	def fit_hist(self, nbins = [100,100,100,100], hlr=[], plr=[(0,120), (0,300), (0,400), (0,1)],  axis=0) :
		
		''' 
		generate a histogram with the range of fit parameter values we are dealing with.

		input:
		### for nbins and lr, each is a 4 element list.
		so [0]=M [1]=A [2]=OFF [3]=XSQ 
		-nbins : number of bins for histogram. 
		-plr : left and right xlimits for plot.
		-hlr : lef and right bounds for histogram bins.

		return:
		a tuple containing numpy arrays (M,A,OFF,XSQ)
		'''

		# std dev of each pixel over the entire dataset.
		#rmsvalues = [self.pix[i].rms for i in range(self.dead, self.nch)]


		#1d view of all signal values from the dataset.
		#filt1d = self.filt.ravel()

		# all detected amplitude, M, and offset values
		M = np.array([pk.tau for px in self.pix for pk in px.peaks])
		A = np.array([pk.fit_par[0] for px in self.pix for pk in px.peaks])
		OFF = np.array([pk.fit_par[2]-px.avg for px in self.pix for pk in px.peaks]) # remove baseline "px.avg" to get true offset.
		XSQ = np.array([pk.chisq for px in self.pix for pk in px.peaks])
		# algorithm to pick out M values and along with their channel.
		# self.Mbad = [(px.number, pk.index, pk.tau) for px in self.pix for pk in px.peaks if pk.tau < something and pk.tau > somethingelse]
		
		print(len(M), 'peaks detected.')


		# note: pyplot histogram function 'axis.hist' or 'pyplot.hist' 
		# calls the numpy histogram function to generate it's plot.

		if isinstance(axis, ax_obj) : # axis supplied
			axis.hist(M, nbins)
			axis.set_xlabel('M = # of time samples')
			axis.set_ylabel('# of peaks')
			axis.set_title('tau histogram')
			axis.set_xlim(plr[0]) # x limits, y limits
			#axis.set_ylim()
			axis.grid(True)

		else : 			# no axis supplied, make standalone plot.
			fig,ax = plt.subplots(1,1)
			ax.hist(M, nbins[0])
			ax.set_xlabel('M, # of time samples')
			ax.set_ylabel('# of peaks')
			ax.set_title('tau histogram')
			ax.set_xlim(plr[0]) # x limits, y limits
			#ax.set_ylim()
			ax.grid(True)

			fig2,ax2 = plt.subplots(1,1)
			ax2.hist(A*1e3, nbins[1])
			ax2.set_xlabel('milli-volts')
			ax2.set_ylabel('# of peaks')
			ax2.set_title('amplitude histogram')
			#ax2.set_xlim(plr[1]) # x limits, y limits
			#ax2.set_ylim()
			ax2.grid(True)

			fig3,ax3 = plt.subplots(1,1)
			ax3.hist(OFF*1e3, nbins[2])
			ax3.set_xlabel('milli-volts')
			ax3.set_ylabel('# of peaks')
			ax3.set_title('voltage offset histogram')
			#ax3.set_xlim(plr[2]) # x limits, y limits
			#ax3.set_ylim()
			ax3.grid(True)

			fig4,ax4 = plt.subplots(1,1)
			ax4.hist(XSQ, nbins[3])
			ax4.set_xlabel('chi-square values')
			ax4.set_ylabel('# of peaks')
			ax4.set_title('least square chi square')
			#ax4.set_xlim(plr[3]) # x limits, y limits
			#ax4.set_ylim()
			ax4.grid(True)

			fig.show()
			fig2.show()
			fig3.show()
			fig4.show()

			return (M,A,OFF,XSQ)

	def reform_data(self) :
		'''
		after analyzing the each Pixel's data and getting 
		filtered result, we can slam it all into one big 
		two-dimensional array for making 2d pixellated images.
		'''
		# self.filt = np.array([i.filt for i in self.pix])

		# always iterating over the total number of channels preserves shape of 2d array in the
		# case that we only analyze a few channels, but still want to use functions that depend on self.row
		self.filt = np.array([self.pix[i].filt for i in range(self.nch)])


	def pixelate_multi(self, start, stop, stepsize, stepthru = False, vmin = -0.001, vmax = 0.007) :

		'''
		input:
		-data : a (5184 x number of time samples) numpy array.
		-start and stop : specify desired points in time to plot.
		-stepsize : decides how many of the time points to plot. if '1',
		all of the data is plotted. if '10' for example, each successive
		10th point in time is plotted. stepsize must divide into integers. 
		-stepthru : if True, press enter to step through frame by frame. if False,
		the images stream through on their own.
		'''

		frames = np.arange(start, stop, stepsize)
		nplot = len(frames)

		d = [p.filt[start:stop:stepsize] for p in self.pix] # 5184 x nplot array
		d = np.reshape(d, (self.row, self.row, nplot)) # 72 x 72 x nplot array


		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		ax.set_title("topmetal data")

		im = ax.imshow(np.zeros((self.row,self.row)), cmap=cm.jet, vmin = vmin, vmax = vmax)
		
		fig.show()
		im.axes.figure.canvas.draw()

		# stream a series of images without user input
		if not stepthru :

			#tstart = time.time()
			for j in range(nplot) :
				#t = j*stepsize*self.tsample
				ax.set_title("Frame %i" % frames[j])
				#ax.set_title("Time elapsed: %f seconds" % t)
				im.set_data(d[:,:,j])
				im.axes.figure.canvas.draw()

		# step through series of images by keyboard input 'enter'
		else : 

			for j in range(nplot) :
				input("press enter for next frame")
				#t = j*stepsize*self.tsample
				ax.set_title("Frame %i" % frames[j])
				#ax.set_title("Time elapsed: %f seconds" % t)
				im.set_data(d[:,:,j])
				im.axes.figure.canvas.draw()

		# close the figure
		plt.close(fig)


	def pixelate_single(self, sample, values = [], vmin=-0.001, vmax=0.007, axis = 0):
		''' Take 5184 channel x 25890 data point array and plot desired points in
		time as 5184 pixel array.
		
		input:
		-sample : desired point in sample space to plot 0-25889
		-values : input pixels to show. useful to check that a selection was created properly.
		-vmin , vmax : minimum and maximum values for the colormap.
		-axis: supply an axis to impose this pixel plot to that axis.
		if no axis is supplied (default), then the function generates a standalone image.

		'''

		if isinstance(axis, ax_obj) :

			# default case, plot data by specifying sample in dataset.
			if not len(values) :
				data_2d = np.reshape(self.filt[:,sample], (self.row, -1)) # convert to square matrix
				# make value bounds for the plot and specify the color map.
				im = axis.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			
			# if 'values' is input, display these pixels
			else :
				arr = np.zeros(self.row**2)
				arr[values] = 1
				data_2d = np.reshape(arr, (self.row, -1))
				im = axis.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			axis.grid(True)

		else :	
			fig, ax = plt.subplots()

			# default case, plot data by specifying sample in dataset.
			if not len(values) :
				data_2d = np.reshape(self.filt[:,sample], (self.row, -1)) # convert to square matrix
				# make value bounds for the plot and specify the color map.
				im = ax.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
			
			# if array is input, plot this image instead.
			else :
				arr = np.zeros(self.row**2)
				arr[values] = 1
				data_2d = np.reshape(arr, (self.row, -1))
				im = ax.imshow(data_2d, cmap=cm.jet, vmin=vmin, vmax=vmax)
		
			fig.colorbar(im)
			ax.grid(True)
			fig.show()


	def plot_waveform(self, pixels=[], choose='b', avg=True, pk = True, fit = False, pd = False, lr = None, us_pt = False) :
		'''
		This function can plot a series of channels. When the user
		presses 'Enter', the next channel's plot is replaces the prior.
		User can plot raw data, filtered data, or both.

		input: 
		-pixels : 1D array/list containing Pixel indices to plot. 
		-choose : 'd' plots raw data, 'f' plots filtered data.
			'b' plots the raw data and filtered data.
		-avg : superimpose a line of the average voltage onto the raw data plot.
		-pk : pk == 'True' superimpose scatter plot of detected peaks.
		-fit : fit == 'True' superimposes scatter plot of fits 
			for a channel's peaks onto the same axis plotting the raw data.
		-lr : input a tuple to slice the plot. for example, if the dataset has
		25890 points, lr = (2000,4000) would only plot points 2000 through 4000.
		lr defaults to plotting the entire dataset. 
		-us_pt : set 'True' to show points used to calculate undershoot when using 'test' option
		inside the Pixel.iter_M() method. Use this option to see if the points used in the calculation 
		are reasonable and cover the correct portion of the trapezoidal waveform.

		'''
		if not lr :
			lr = (0, self.daq_length)

		x = np.arange(lr[0], lr[1])
		datalen = len(x)

		# setup

		# for raw data
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_title("multich plot")
		ax.set_xlim(lr)
		fig.show()

		# for filtered data

		# for peak detection filtered data
		fig3, ax3 = plt.subplots(1,1)
		ax3.set_title("multich plot")
		ax3.set_xlim(lr)
		fig3.show()

		# plot raw data and filtered data.
		if (choose == 'b') :
			
			fig2 = plt.figure()
			ax2 = fig2.add_subplot(111)
			ax2.set_xlim(lr)
			# ax2.step(x, np.zeros(datalen))
			fig2.show()
			#ax2.figure.canvas.draw()

			print('plotting raw data and filtered data')
			# loop over desired pixels
			for i in pixels :
				input("press enter for next channel")
				print('pixel number', i)
				print('avg: %f, RMS: %f' % (self.pix[i].avg, self.pix[i].rms))

				# print peaks' location, tau, and chisq
				peaks = [j for j in self.pix[i].peaks]


				try : # if M_it is defined.
					print('{0: >5s}, {1: >5s}, {2: >7s}, {3: >7s}, {4: >5s}'.format('loc', 'tau', 'chisq', 'M_it', 'type'))
					for k in peaks:
						print('{0: >5d}, {1: >5f}, {2: >5f}, {3: >5f}, {4: >5s}'.format(k.index,k.tau, k.chisq, k.M_it, k.id))
					print('\n')

				except (IndexError, TypeError) as e : # M_it isn't defined.
					print('{0: >5s}, {1: >5s}, {2: >7s}, {3: 5s}'.format('loc', 'tau', 'chisq', 'type'))
					for k in peaks:
						print('{0: >5d}, {1: >5f}, {2: >5f}, {3: >5s}'.format(k.index,k.tau, k.chisq, k.id))
					print('\n')

				# if peaks[0][3] :
				# 	# if M_it is defined
				# 	print('{0: >5s}, {1: >5s}, {2: >7s}, {3: >7s}'.format('loc', 'tau', 'chisq', 'M_it'))
				# 	for k in peaks:
				# 		print('{0: >5d}, {1: >5f}, {2: >5f}, {3: >5f}'.format(k[0],k[1], k[2], k[3]))
				# else : 
				# 	# if M_it isn't defined.
				# 	print('{0: >5s}, {1: >5s}, {2: >7s}'.format('loc', 'tau', 'chisq'))
				# 	for k in peaks:
				# 		print('{0: >5d}, {1: >5f}, {2: >5f}'.format(k[0],k[1], k[2]))

				ax.cla()
				ax.set_title('raw data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))

				# plot desired slice of data.
				ax.step(x, self.pix[i].data[lr[0]:lr[1]], color='b')

				# check for avg/pk/fit inputs to plot these features
				if avg :
					ax.step(x, np.ones(datalen)*self.pix[i].avg, color = 'y')



				if pk :
					# take peak objects only inside of the window we are plotting, determined by 'lr' tuple
					for j in peaks :
						if (j.index  < lr[1] and j.index > lr[0]):
							ax.scatter(j.index, self.pix[i].data[j.index], color='r', marker ='x')

				if fit :
					# impose a scatter plot for each fitted peak onto the raw data.
					# only plot fits that lie inside the window.
					for j in peaks :
						if (j.index + Pixel.fit_length  < lr[1] and j.index > lr[0]) :
							ax.scatter(j.fit_pts + j.index, 
								model_func(j.fit_pts, *j.fit_par), 
								marker = 'o', color = 'g')

				# now plot the filtered data.		
				ax2.cla()
				ax2.set_title('filtered data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))
				ax2.step(x, self.pix[i].filt[lr[0]:lr[1]], color='b')
				
				# superimpose peaks if desired.
				if pk :
					for j in peaks :
						if (j.index  < lr[1] and j.index > lr[0]):
							ax2.scatter(j.index, self.pix[i].filt[j.index], color='r', marker ='x')
				
				if us_pt :

					for pk in self.pix[i].peaks :

						ax2.scatter(pk.pre, self.pix[i].filt[pk.pre], marker ='x', color='olive')
						ax2.scatter(pk.flat, self.pix[i].filt[pk.flat], marker ='x', color='indigo')
						ax2.scatter(pk.post, self.pix[i].filt[pk.post], marker ='x', color='orange')

				if pd :
				
					ax3.cla()
					ax3.set_title('peak det filtered data, channel no. %i : (%i, %i)' 
					% (i, self.pix[i].loc[0], self.pix[i].loc[1]))

					# plot desired slice of data.
					ax3.step(x, self.pix[i].pd_filt[lr[0]:lr[1]], color='b')

					if pk :
						# take peak objects only inside of the window we are plotting, determined by 'lr' tuple
						for j in peaks :
							if (j.index  < lr[1] and j.index > lr[0]):
								ax3.scatter(j.index, self.pix[i].pd_filt[j.index], color='r', marker ='x')

				# ax is raw data, ax2 is filtered data.
				ax.figure.canvas.draw()
				ax2.figure.canvas.draw()
				ax3.figure.canvas.draw()
		# plot the raw data.
		elif (choose == 'd') :
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

		input('press enter to close figures')
		plt.close(fig) # raw 
		plt.close(fig2) # filt

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
	### shaper_single is broken - gives trapezoidal response but pulse height is wrong... and dependent on l,k.
	### as a temporary fix, we just implemented all of the filtering with shaper_multi, even for one and zero peaks.

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
	# they're used by the analysis chain.
	
	daq_length = None
	noisethresh = None

	threshold = None
	minsep = None
	
	fit_length = None
	fudge = None

	l = None
	k = None
	M_def = None
	M_step = -1
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
		 avoid transition on peaks' rise time, which we don't account for in our simple fit model.
		-npk : total number of peaks. -1 until set to either 0 or a positive
		integer by the peak_det() method.
		-type : descriptor for pixel. some pixels are okay, others are noisy,
		3 are for demux and some show weird basleine behavior, may be do to alpha particle impacts, or just bad pixels.
		'''

		self.number = number
		self.loc = (number % Pixel.row, int(number/Pixel.row))
		self.data = data
		self.filt = np.empty_like(data) # it's nice to initialize here, because if we only do operations on some channels,
		# having an array for every pixel will preserve the square shape of the array if we want to reform data anyways.
		
		self.avg = avg
		self.rms = rms
		self.peaks = []

	def peak_det(self):
		''' 
		input : 
		-choose: choose between using savitsky-golay filter or trapezoidal filter with tiny l-k and k. 
		relevant pixel attributes:
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

		# if Pixel.pd_choose == 'sg' : # savitsky-golay filter
		f = savgol_scipy(self.data, Pixel.sgwin, Pixel.sgorder) # smooth data
		# elif Pixel.pd_choose == 'skinny' : # trapezoidal filter with l-k and k both very small (<5)
		# 	f = np.empty_like(self.data)
		# 	Pixel.shp_lib.shaper_multi(self.data, f, c_ulong(len(self.data)), byref(Pixel.skinny_pk), c_double(self.avg))
		# 	f = np.array(f)

		kernel = [1, 0, -1]	# calculate each derivative
		df = convolve(f, kernel, 'valid') # note: each convolution cuts length of array by len(kernel)-1
		dfn = np.sign(df) # normalize derivative to 1. 1 increasing, -1 decreasing, 0 const.
		
		# the second derivative of the normalized derivative.
		# should only have non-zero values for peaks and valleys, where the value of the derivative changes.
		ddfn = convolve(dfn, kernel, 'valid') 

		# first, find all of the positive derivative values. going up the peak.
		# this returns indices of possible candidates. we want to offset by two because
		# the convolution cuts the array length by len(kernel)-1
		if (Pixel.sign == 1) :
			candidates = np.where(dfn > 0)[0] + (len(kernel)-1)
		elif (Pixel.sign == -1) :
			candidates = np.where(dfn < 0)[0] + (len(kernel)-1)

		pk = sorted(set(candidates).intersection(np.where(ddfn == -1*Pixel.sign*2)[0] + 1))
		alpha = self.avg + (Pixel.sign * Pixel.threshold)

		if (Pixel.sign == 1) :	# apply threshold to the raw data, not the smoothed data 
			pk = np.array(pk)[self.data[pk] > alpha]
		elif (Pixel.sign == -1) :
			pk = np.array(pk)[self.data[pk] < alpha]

		# list comprehension version of toss np.array for minimum separation discrimination.
		# list is faster for 'append' operations, especially when there are many false peaks.
		# 'toss' gives the indices of the bad peaks in the array of peak locations, 
		# not the indices of the peaks in the dataset.
		toss = [i+1 for i in range(len(pk)-1) if ((pk[i+1]-pk[i]) < Pixel.minsep)]
		
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
		pkm = np.delete(pk, toss) #- (len(kernel)-1)

		# cons = 5 # consecutively increasing values preceeding a peak
		# for j in range(len(pkm))
		# 	for k in range(cons)

		# use the 'trig' namedtuple for debugging / accessing each step of the peak detection.
		#return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=pk, toss=toss, pkm=pkm)

		# pull these out for debugging if desired.
		#self.prepeak = pk
		#self.toss = toss

		# create a list of 'Peak' objects with these peak locations.
		self.peaks = [Peak(i) for i in pkm]

		# option to save the filtered data. Can plot this later to look at peak det results before 
		# being imposed onto the raw data.
		# if Pixel.pd_save:
		# 	self.pd_filt = f
		# else :
		# 	pass
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
			if (Pixel.daq_length - self.peaks[-1].index < Pixel.fit_length + Pixel.fudge) :
				self.peaks.pop()

			npk = len(self.peaks)

			# assumption: data is modeled by a 3 parameter decaying exponential.
			M = 3

			# offset the fit parameters by the baseline of the signal

			self.bounds = Pixel.fit_bounds + np.array(([0,0,self.avg],[0,0,self.avg]))
			self.guess = Pixel.fit_guess + np.array([0,0,self.avg])



			# include a 'fake' peak to act as endpoint for last iteration.
			self.peaks.append(Peak(Pixel.daq_length-1))

			# for q in self.peaks:
			# 	if q+1

			for j in range(npk) :
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

				par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, p0=self.guess, \
				sigma=np.ones(N)*self.rms, absolute_sigma=True, check_finite=False, \
				bounds=self.bounds, method='trf')


				# par, cov = curve_fit(f=model_func, xdata=xi, ydata=yi, \
				# 		             sigma=np.ones(N)*rms, absolute_sigma=True, check_finite=False, \
				# 		             bounds=bounds, method='trf')


				f_xi = model_func(xi, *par)
				chisq, pval = chisquare(f_obs=yi, f_exp=f_xi, ddof=N-M)
				
				# each Pixel object has a list of peak objects,
				# determined by the number of peaks found in peak detection.
				# insert the fitted tau and chisquare to their respective peak.
				self.peaks[j].tau = 1.0/par[1]

				# assign identifying string based on tau value.

				# peaks with very small M tend to be noise spikes
				if self.peaks[j].tau < Pixel.Mrange[0] :
					self.peaks[j].id = 'noise'
				
				# peaks with large M tend to be made up of multiple constituent peaks
				elif self.peaks[j].tau > Pixel.Mrange[1] :
					self.peaks[j].id = 'multi'


				else :
					self.peaks[j].id = 'single'

				self.peaks[j].chisq = chisq
				
				# these attributes can be saved to plot fits.
				self.peaks[j].fit_par = par
				self.peaks[j].fit_pts = xi
				


				# Q and p values.
				# P[j] = pval
				# Q[j] = 1-pval

				# if axis object is provided, add the fit as a scatter plot to the axis.
				# if isinstance(ax, ax_obj) :
				# 	ax.scatter(xi+peaks[j], model_func(xi,*par), marker = 'o')

			# remove the 'fake' peak now that we're done iterating over
			# the fit procedure.
			
			self.peaks.pop()
	
	def peak_cut(self, tau_range, chisq_range) :

		''' function to remove peaks from the dataset based on some criteria.
		This is done separately from minimum separation filtering because 
		we need the least-squares fit information to do this.

		input:
		-tau_range : range of acceptable tau values. outside of this range, peaks will be removed.
		-chisq_range : range of acceptable chi square values. outside of this range, peaks will be removed.
		'''

		pass


	def filter_peaks(self, step=False, simple=False, iter=False):

		''' apply trapezoidal filter to data. peak locations are used to change
		filter parameters.
		input:
		-step : if desired, filter a single step input (sets M = -1)
		-simple : filter the data with the default M value.
		-iter : filter the data using the the newly iterated value of M
		 '''

		# only filter peaks if they have tau values in the range we're interested in.
		# print(peaks)
		# print(npk)
		pkchoose = [i for i in self.peaks if i.id == 'single']
		npk = len(pkchoose)
		for i in pkchoose :
			print(i.index)
		print('pkchoose = ', pkchoose, 'npk = ', npk)
		print(Pixel.zero_pk)
		# no peaks found or if skipping peak detection and fitting ('simple' option), use default M. 
		if (npk == 0 or simple == True) : # no peaks found, use default M.
			Pixel.shp_lib.shaper_multi(self.data, self.filt, c_ulong(len(self.data)), byref(Pixel.zero_pk), c_double(self.avg))
			print('reached this point for zero peak channel.')
		elif (npk == 1) : # one peak found.
			l_arr = np.ones(npk, dtype = c_ulong)*Pixel.l
			k_arr = np.ones(npk, dtype = c_ulong)*Pixel.k
			LEFT = np.array(0, dtype = c_ulong)
			RIGHT = np.array(Pixel.daq_length, dtype = c_ulong)
			if iter :
				M = np.array([i.M_it for i in pkchoose])
			else :			
				M = np.array([i.tau for i in pkchoose])
			PEAK = peaks_handle(npk, LEFT, RIGHT, l_arr, k_arr, M)
			Pixel.shp_lib.shaper_multi(self.data, self.filt, c_ulong(len(self.data)), byref(PEAK), c_double(self.avg))

		elif (npk > 1) : # multiple peaks found.
			l_arr = np.ones(npk, dtype = c_ulong)*Pixel.l
			k_arr = np.ones(npk, dtype = c_ulong)*Pixel.k
			LR = self.pk2LR(pkchoose, npk)

			# leave this comment in to check values later if desired
			# print('l', l)
			# print('k', k)
			# print('M', M)
			# print('number of peaks =', npk)
			# print('LEFT = ', LR[0])
			# print('RIGHT = ', LR[1])
			if iter :
				M = np.array([i.M_it for i in pkchoose])
			else :			
				M = np.array([i.tau for i in pkchoose])
			# peaks_handle prepares a Ctypes 'Structure' object to make a C Structure.
			PEAK = peaks_handle(npk, LR[0], LR[1], l_arr, k_arr, M)
			# self.filt = np.empty_like(self.data)
			Pixel.shp_lib.shaper_multi(self.data, self.filt, c_ulong(len(self.data)), byref(PEAK), c_double(self.avg))

		elif step :
			Pixel.shp_lib.shaper_single(self.data, self.filt, c_ulong(Pixel.daq_length), 
				c_ulong(Pixel.l), c_ulong(Pixel.k), c_double(Pixel.M_step), c_double(self.avg))


		self.filt = np.array(self.filt)

	def pk2LR(self, pkchoose, npk) :
	
		''' put peak locations into arrays of left and rights for trapezoidal shaper.
		 apply desired offset. Both arrays are the same length as input 'pkchoose' array.

		return:
		- LEFT and RIGHT : beginning and end points for each set of 
		trapezoidal filter parameters l,k, and M.'''
		

		LEFT = np.zeros(npk)
		RIGHT = np.zeros(npk)

		for i in range(npk-1):
			LEFT[i]  = pkchoose[i].index   + self.shaper_offset
			RIGHT[i] = pkchoose[i+1].index + self.shaper_offset
			
		LEFT[0] = 0
		LEFT[npk-1] = pkchoose[npk-1].index + self.shaper_offset
		RIGHT[npk-1] = Pixel.daq_length

		# trapezoidal filter uses size_t, or c_ulong, as its datatype
		# for left and right. they index locations possibly larger than int allows.
		LEFT = np.array(LEFT, dtype = c_ulong)
		RIGHT = np.array(RIGHT, dtype = c_ulong)

		return (LEFT, RIGHT)

	def iter_M(self, step_it = 1, max_it = 30, tail=50, midflat=0, npt_avg = 10, avg_off = -2, t_os=0.0001, t_us = 0.0001,  thresh_us=0.0001, squeeze=0, option='symm_test') :
		'''
		a function to iteratively correct the M values used for each peak.
		increases or decreases M to minimize the undershoot / overshoot.
		This is important because the trapezoidal voltage summation is sensitive to error in M.
		
		observation: 
		1. if M is bigger than the true tau, 
		there is undershoot and the flat top bulges on the left side.
		
		2. if M is smaller than the true tau, 
		there is overshoot and the flat top bulges on the right side.

		the flat-top peak is referred to as 'bump' below in the algorithm.
		
		input:

		-step_it: the stepsize for each iteration of M. smaller values can be more accurate,
		but will take longer to run.
		-max_it: maximum number of iterations, incase we keep missing the edge case.
		a reasonable number of iterations can be deduced from the range of acceptable M values
		and the initial guess.
		-option: two different ways to identify errror in M.
		-npt_avg, off_us, thresh_us : parameters for undershoot calculation. 
		'npt_avg' is the number of points to average, 'off_us' offsets the average from l+k points after the peak,
		and 'thresh_us' is for deciding when the undershoot / overshoot are too big. 
		1. 'bump' identifies characteristic bump on l/r of the flat-top section of the trapezoid.
		this method did not work well in the presence of noise on the scale we are dealing with.
		2. 'urs' calculates undershoot / overshoot 

		'''

		# initialize new M we get from iteration.
		# for pk in self.peaks :
		# 	pk.M_it = pk.tau


		# improve M by measuring symmetry of both k-ramps (the legs of the trapezoid, each has k-samples).
		if option == 'symm_test' :
			
			go = len(self.peaks)
			
			for pk in self.peaks :

				pk.pre = np.arange(pk.index + squeeze + avg_off, pk.index + Pixel.k - squeeze + avg_off)
				pk.flat= np.arange(pk.index + Pixel.k + avg_off, pk.index + Pixel.l + avg_off)
				pk.post = np.arange(pk.index + Pixel.l + squeeze + avg_off, pk.index + Pixel.l + Pixel.k - squeeze + avg_off)

				rise = np.sum(self.filt[pk.index + squeeze + avg_off : pk.index + Pixel.k - squeeze + avg_off])
				fall = np.sum(self.filt[pk.index + Pixel.l + squeeze + avg_off : pk.index + Pixel.l + Pixel.k - squeeze + avg_off])


				# if sym is positive 'blah', if it's negative then 'blah'
				sym = rise-fall

				if sym > t_us : # there is undershoot
					pk.M_it -= step_it
					go -=1
				elif sym < -1*t_os : # there is overshoot
					pk.M_it += step_it
					go -=1
				else :
					pass

				print('M=%f, M_it=%f, sym=%f, go=%i' % (pk.tau, pk.M_it, sym, go))
				self.filter_peaks(iter=True)

		if option == 'symm_test2' :
			
			go = len(self.peaks)
			
			for pk in self.peaks :

				pk.pre = np.arange(pk.index - tail + avg_off, pk.index + Pixel.k + midflat + avg_off)
				pk.flat= np.arange(pk.index + Pixel.k + avg_off, pk.index + Pixel.l + avg_off)
				pk.post = np.arange(pk.index + Pixel.l - midflat + avg_off, pk.index + Pixel.l + Pixel.k + tail + avg_off)

				rise = np.sum(self.filt[pk.index - tail + avg_off : pk.index + Pixel.k + midflat + avg_off])
				fall = np.sum(self.filt[pk.index + Pixel.l - midflat + avg_off: pk.index + Pixel.l + Pixel.k + tail + avg_off])


				# if sym is positive 'blah', if it's negative then 'blah'
				sym = rise-fall

				if sym > t_us : # there is undershoot
					pk.M_it -= step_it
					go -=1
				elif sym < -1*t_os : # there is overshoot
					pk.M_it += step_it
					go -=1
				else :
					pass

				print('M=%f, M_it=%f, sym=%f, go=%i' % (pk.tau, pk.M_it, sym, go))
				self.filter_peaks(iter=True)

		if option == 'symm_test3' :
			

			# for this option, let's try and compare to only the original values, which we will save.
			# 

			go = len(self.peaks)
			
			for pk in self.peaks :

				pk.pre = np.arange(pk.index - tail + avg_off, pk.index + Pixel.k + midflat + avg_off)
				pk.flat= np.arange(pk.index + Pixel.k + avg_off, pk.index + Pixel.l + avg_off)
				pk.post = np.arange(pk.index + Pixel.l - midflat + avg_off, pk.index + Pixel.l + Pixel.k + tail + avg_off)

				rise = np.sum(self.filt[pk.index - tail + avg_off : pk.index + Pixel.k + midflat + avg_off])
				fall = np.sum(self.filt[pk.index + Pixel.l - midflat + avg_off: pk.index + Pixel.l + Pixel.k + tail + avg_off])


				# if sym is positive 'blah', if it's negative then 'blah'
				sym = rise-fall

				if sym > t_us : # there is undershoot
					pk.M_it -= step_it
					go -=1
				elif sym < -1*t_os : # there is overshoot
					pk.M_it += step_it
					go -=1
				else :
					pass

				print('M=%f, M_it=%f, sym=%f, go=%i' % (pk.tau, pk.M_it, sym, go))
				self.filter_peaks(iter=True)	
		if option == 'symm' :

			npk = len(self.peaks)
			nit = max_it

			while(nit) :
				
				nit -= 1
				# 'squeeze' removes a specified number of points from each side of the ramp.
				# should we calculate average or just the sum...?
				for pk in self.peaks :

					rise = np.sum(self.filt[pk.index + squeeze : pk.index + Pixel.k - squeeze])
					fall = np.sum(self.filt[pk.index + Pixel.l + squeeze : pk.index + Pixel.l + Pixel.k - squeeze])


					# if sym is positive 'blah', if it's negative then 'blah'
					sym = rise-fall

					if sym > thresh_us : # there is undershoot
						pk.M_it += step_it
					elif sym < -1*thresh_us : # there is overshoot
						pk.M_it -= step_it
					print('M=%f, M_it=%f, us/os=%f, go=%i' % (pk.tau, pk.M_it, pk.shoot, go))
				self.filter_peaks(iter=True)

		# improve M by using flat top peaking
		if option == 'bump' :

			# this first loop is to initialize iteration variables

			npk = len(self.peaks) 

			for pk in self.peaks :
				bump = np.argmax(self.filt[pk.index : pk.index + Pixel.l + Pixel.k])
				
				# M is too small
				if (bump > (pk.index + Pixel.l + Pixel.k)/2) :
					pk.dir_it = 1

				# M is too big 
				else :
					pk.dir_it = -1

			self.filter_peaks(iter=True)

			# Main Loop.
			# keep going until all of the peaks are done. 
			# we'll decrement 'npk' everytime a peak is done iterating.


			# this loop gets into a trap. once a M is incorrectly made +/-, 
			# it keeps making the wrong condition filled, getting stuck in a loop forever.
			# is an undershoot calculation possible with so much noise?

			i = 0 
			while(npk) :
				i += 1
				print(i)


				for pk in self.peaks :

					bump = np.argmax(self.filt[pk.index : pk.index + Pixel.l + Pixel.k])

					# M is too small
					if (bump > (pk.index + Pixel.l + Pixel.k)/2) :

						# Keep making M bigger.
						if pk.dir_it == 1 :
							pk.M_it += step_it

						# M was too big before, but now we've made it too small. 
						# we're done.
						elif pk.dir_it == -1 :
							pk.dir_it = 0
							npk -= 1 # decrement number of remaining peaks to iterate

					# M is too big
					else :

						# Keep making M smaller.
						if pk.dir_it == -1 :
							pk.M_it -= step_it

						# M was too small before, but now we've made it too big. 
						# we're done.
						elif pk.dir_it == 1 :
							pk.dir_it = 0
							npk -= 1

					self.filter_peaks(iter=True)

		# calculate undershoot wrt baseline
		elif option == 'us_bl' :
						
			# get the region of interest for each peak.
			for pk in self.peaks :
				pk.start = pk.index + Pixel.l + Pixel.k + avg_off
				pk.stop = pk.start + npt_avg

			npk = len(self.peaks)
			go = True 
			n_it = max_it

			# loop ends if all peaks are done, or do max number of iterations.
			while(go and n_it) :
				n_it -= 1
				go = len(self.peaks)
				for pk in self.peaks : 
						
					# calculate average after trapezoid 
					pk.shoot = np.average(self.filt[pk.start:pk.stop])

					# could later make separate threshold's for undershoot / overshoot if necessary,
					# at a glance they don't appear symmetric.
					
					# excessive undershoot (M is too big)
					if pk.shoot > thresh_us :
						pk.M_it += step_it

					# excessive overshoot (M is too small)
					elif pk.shoot < -1*thresh_us : 
						pk.M_it -= step_it

					else : 
						go -= 1

					self.filter_peaks(iter=True)

# calculate undershoot wrt baseline, let's try and get this working properly.
		elif option == 'test_bsl' :
						
			# get the region of interest for each peak.
			for pk in self.peaks :
				pk.start = pk.index + Pixel.l + Pixel.k + avg_off
				pk.stop = pk.start + npt_avg

			npk = len(self.peaks)
			# go = True
			# while(go) :
			go = len(self.peaks)
			for pk in self.peaks : 
					
				# calculate average after trapezoid 
				pk.shoot = np.average(self.filt[pk.start:pk.stop])

				# could later make separate threshold's for undershoot / overshoot if necessary.
				
				# excessive overshoot (M is too big)
				if pk.shoot > thresh_us :
					pk.M_it += step_it
				
				# excessive undershoot (M is too small)
				elif pk.shoot < -1*thresh_us : 
					pk.M_it -= step_it

				else : 
					go -= 1
				print('M=%f, M_it=%f, us/os=%f, go=%i' % (pk.tau, pk.M_it, pk.shoot, go))
				
			self.filter_peaks(iter=True)

		# calculate undreshoot wrt flat-top
		elif option == 'test_ftop' : 

			for pk in self.peaks : 

				pk.pre = np.arange(pk.index - npt_avg + avg_off[0], pk.index + avg_off[0]) 
				pk.flat = np.arange(pk.index + Pixel.k + avg_off[1], pk.index + Pixel.l + avg_off[1])
				pk.post = np.arange(pk.index + Pixel.l + Pixel.k + avg_off[2],
							pk.index + Pixel.l + Pixel.k + npt_avg + avg_off[2])
				
				# calculate pre trapezoid section.
				pre=np.average(self.filt[pk.index - npt_avg + avg_off[0] : pk.index + avg_off[0]])	

				# calculate flat-top section
				flat=np.average(self.filt[pk.index + Pixel.k + avg_off[1] : pk.index + Pixel.l + avg_off[1]])

				# calculate post trapezoid section.
				post=np.average(self.filt[pk.index + Pixel.l + Pixel.k + avg_off[2] : 
							pk.index + Pixel.l + Pixel.k + npt_avg + avg_off[2]])
			
				pk.A = flat-pre
				pk.B = flat-post
				pk.shoot = pk.A-pk.B

				# M is bigger than tau, then there is undershoot.
				if pk.shoot > thresh_us :
					pk.M_it += step_it

				# M is smaller than tau, then there is overshoot.
				elif pk.shoot < -1*thresh_us:
					pk.M_it -= step_it

				print('M=%f, M_it=%f, shoot=%f, pre=%f, post=%f' % (pk.tau, pk.M_it, pk.shoot, pk.A, pk.B))
			self.filter_peaks(iter=True)	
		
		elif option == 'diff' :
			npk = len(self.peaks)

			if (len(avg_off)) == 3 :
				# loop terminates when maximum iterations completed, or all peaks are within the 
				# undershoot / overshoot threshold (when 'go' == 0).
				go = True 
				while(max_it and go):
					max_it -= 1
					go = npk
					for pk in self.peaks : 

						# calculate pre-trapezoid section.
						pre=np.average(self.filt[pk.index - npt_avg + avg_off[0] : pk.index + avg_off[0]])	

						# calculate flat-top section
						# flat=np.average(self.filt[pk.index + Pixel.k + avg_off[1] : pk.index + Pixel.l + avg_off[1]])

						# calculate post-trapezoid section.
						post=np.average(self.filt[pk.index + Pixel.l + Pixel.k + avg_off[2] : 
							pk.index + Pixel.l + Pixel.k + npt_avg + avg_off[2]])
					
						# A = flat-pre
						# B = flat-post
						# shoot = A-B

						diff = post-pre

						# M is bigger than tau, then there is undershoot.
						if diff > thresh_us :
							pk.M_it += step_it

						# M is smaller than tau, then there is overshoot.
						elif diff < -1*thresh_us :
							pk.M_it -= step_it
						# undershoot is within bounds, done.
						else :
							go -=1

					self.filter_peaks(iter=True)

			else :
				print('"iter_M not executed.')
				print('"avg_off" requires 3 elements, one for offsetting each part of the waveform:')
				print('Before the trapezoid, during the flat top, and after the trapezoid.')

	def enter_peaks(self, auto=True) :


		''' function to quickly enter M values into peaks via command line, for purposes of testing.
		'''


		if auto:
			for pk in self.peaks:
				pk.M_it = pk.tau

		else:
			i = 0
			for pk in self.peaks :
				i+=1
				print('peak %i index: %i tau: %f chisq: %f' % (i, pk.index, pk.tau, pk.chisq))

				while True :
					try :
						val = float(input('enter an M value. \n'))
						index = int(input('enter a peak location. \n'))
						break
					except ValueError :
						print('peak location and M value must have type "int" or "float". \n')

				pk.M_it = val
				pk.index = index


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
			lr = (0, Pixel.daq_length-1)

		xaxis = np.arange(Pixel.daq_length)[lr[0]:lr[1]]
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
	'''
	object for peak instances found in the dataset, 
	determined by peakdet parameters set in the analysis chain.
	'''
	def __init__(self, index) :
		
		'''
		attributes:
		-index: location in dataset of peak. (maximum value)
		-tau: decay constant extracted by least-squares fit.
		
		-id : this string identifies three types of detected peaks.
		1. 'noise' are peaks with very short time constants, they're just large and quick jitters due to noise,
		not the characteristic impulse response of the CSA.
		2. 'single' are peaks with time constants of 15-50 samples. the variation is due to different
		Rf for each pixel, and the error in fitting these pulses.
		3. 'multi' are peaks with exceptionally long fall times. so far, it looks like all of these are multiple peaks
		in close proximity. We can observe constituent peaks' locations with a very small l-k trapezoidal filter,
		if l-k is larger they aren't observable.

		-chisq : chi-square value from least-squares fit
		-shoot : undershoot / overshoot due to trapezoidal filtering.
		positive value -> overshoot (tau is too small), negative value -> undershoot (tau is too big)
		'''

		self.index = index
		self.tau = None
		self.id = None
		self.fit_pts = None
		self.fit_par = None
		self.M_it = None # this M we get from iterating to minimize the undershoot / overshoot
		self.chisq = None
		self.shoot = None

class Event(object):
	row = 72
	def __init__(self, x, y, i, r=None, a=None, b=None, angle=0, shape='c'):
		'''
		input: behave differently for different shape strings, read below.

		shape == 'c' (circle)
		-x,y : ints representing center of circle.
		-radius : radius of circular selection

		shape == 'e' (ellipse)
		-x,y : ints representing center of ellipse
		-a,b : length of major and minor axis. whichever is larger
		is automatically used as the major axis.
		-angle : angle of rotation (in degrees) from the x axis being 0 degrees.
		when angle rolls over 90, it doesn't work properly, use this method:
		0 -> 90 rotates the ellipse counter clockwise (starting lying horizontally)
		0 -> -90 rotates the ellipse clockwise (starting lying horizontally)

		shape == 'r' (rectangle)
		-x : a tuple with (left, right) points
		-y : a tuple with (top, bottom) points

		'''
		#self.npix is set after pixels from selection are retrieved.
		#self.row is set when we label all of the events, in Sensor.label_events()
		self.r = r
		self.x = x
		self.y = y
		self.a = a
		self.b = b
		self.angle = angle * np.pi/180 # convert degree to radians
		self.i = i
		self.shape = shape


	def retrieve_selection(self) :

		if self.shape == 'c' :
			self.circle()
		elif self.shape == 'r' :
			self.rectangle()
		elif self.shape == 'e' :
			self.ellipse()

	def ellipse(self) :

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
		self.angle =  self.angle-np.pi
		sin = np.sin(self.angle)
		sin2 = (np.sin(self.angle))**2
		cos = np.cos(self.angle)
		cos2 = (np.cos(self.angle))**2
		a2 = self.a**2
		b2 = self.b**2

		# make a rectangle based on the major axis
		# needs to be able to account for rotation
		if self.a > self.b :
			major = self.a
		else :
			major = self.b

		rect = [((i, j), i+j*self.row) for j in range(self.y-major, self.y+major) for i in range(self.x-major,self.x+major)]

		# cut out the ellipse. this is just an ellipse distance equation with an angle included..
		self.sel = np.array([rect[i][1] for i in range(len(rect)) if \
		(((((rect[i][0][0]-self.x)*cos+(rect[i][0][1]-self.y)*sin2)**2)/a2) + \
		((((rect[i][0][0]-self.x)*sin-(rect[i][0][1]-self.y)*cos2)**2)/b2)) < 1])

		self.npix = len(self.sel)



	def circle(self) :

		# circular selection.
		rect = [((i, j), i+j*self.row) for j in range(self.y-self.r, self.y+self.r) for i in range(self.x-self.r,self.x+self.r)]
		# linear location in 5184 element array of a 'radius' sized rectangle
		#rectangle = [((self.row)*i)+j for j in range(x_0-radius, x_0+radius) for i in range(y_0-radius, y_0+radius)]
		# list of ((x,y), i) tuples -> (xy coordinates, linear index) 
		#xy = [(self.pix[i].loc, i) for i in rectangle]
		# if the element is inside of the circle's radius, include it. We're cutting a circle out of a rectangle
		self.sel = np.array([rect[i][1] for i in range(len(rect)) if np.sqrt((rect[i][0][0]-self.x)**2 + (rect[i][0][1]-self.y)**2) < self.r])
		
		self.npix = len(self.sel)


	def rectangle(self) :

		'''make a rectangular selection of pixels.
		input:
		-lr : the left and right most points as a two element list.
		-tb : the top and bottom most points as a two element list.

		'''
		self.sel = np.array([i+j*self.row for j in range(self.y[0], self.y[1]) for i in range(self.x[0], self.x[1])])
		self.npix = len(self.sel)


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


def close_figs(): 
	'''
	closes all active matplotlib.pyplot figures
	'''

	plt.close("all")


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







