# some useful functions for testing things that are used frequently.


# for importing data from HDF5 files into Numpy arrays
import h5py 

# Ctypes and Numpy support for calling C functions defined in shared libraries.
# These are functions too slow when implemented in pure Python.
from ctypes import *
import numpy.ctypeslib as npct

# plotting tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

# math/matrices, statistics, and fittting
import numpy as np
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt, convolve, savgol_filter

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

def test_shaper(l, k, M):

	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)
	fig3, axis3 = plt.subplots(1,1)


	step = np.ones(2000, dtype=np.float32) 
	step[:1000] = 0.828 
	step[1000:] = 0.978

	exp = np.ones(2000, dtype=np.float32)*0.828
	exp[1000:] += 0.01*np.exp(-(1/40.)*np.linspace(0, 999,1000))
	
	filt = shaper_np(exp, l, k, M)
	#filt = shaper_np(step, 20, 10, -1)

	plot(exp, axis)
	#plot(step, axis2)
	plot(filt, axis3)
	

	fig.show()
	#fig2.show()
	fig3.show()

def shaper_np(data, l, k, M):
	# apply trapezoidal filter to a numpy array, return the numpy array 
	# for quick analysis / plotting.

	# import the library
	lib = CDLL("shaper.so")

	# make pointer type for 1-D arrays of floats.
	array_float = npct.ndpointer(c_float, ndim=1, flags='CONTIGUOUS')


	
	# define argument types. we just defined pointer-arrays for in/out.
	lib.trapezoid.restype = None

						  # 	  in 		  out   	 length  	l  		k 		  M
	lib.trapezoid.argtypes = [array_float, array_float, c_ulong, c_ulong, c_ulong, c_double]
 
	# allocate an array to hold output.
	filt = np.empty_like(data)
	lib.trapezoid(data, filt, c_ulong(len(data)), c_ulong(l), c_ulong(k), c_double(M))
	
	return np.array(filt)

def savgol(array, npt, order):
	
	out = savgol_filter(array, npt, order)
	return out

def fit_test(infile, pixel, pk, fit_length):

	timestep = (4*72**2)*(3.2*10**-8)
	# let's make 4 total plots. one for entire channel (25890 pts), then three
	# other plots to see how some of the fits look.
	fig1, raw_ax = plt.subplots(1,1)

	raw = pull(infile, pixel)
	plot(raw, raw_ax) # plot the raw data


	filt = savgol(raw, 15, 4)
	plot(filt, raw_ax) # plot the smoothed data

	peaks = img_derivative(filt, raw_ax) # plot peak locations
	#peaks = peakdet_cwt(filt, axis[0,0])    # plot peak locations


	# fit a single peak first
	# do a rough investigation if to some degree the peak value is associated
	# with the fall time of the pulses we are going to fit... for now we'll just fix it to 
	# a few different values and see the Q values from the chisquare test.

	# get the data for a pulse
	pulse = raw[peaks[pk]:peaks[pk]+fit_length]
	fig2, fit_ax = plt.subplots(1,1)
	baseline = np.mean(raw[:100])
	rms = (np.std(raw[:100]))*np.ones(len(pulse))

	x = np.linspace(0, fit_length-1, fit_length)
	y = pulse

	def model_func(x, amp, tau, offset):
		return amp*np.exp(-tau*(x))+offset

	par, cov = curve_fit(model_func, x, y, p0=[0.008, 1/35., baseline], sigma=rms)
	M = 1.0/par[1]
	print 'tau:', M, 'samples', M*timestep, 'seconds'
	print 'amplitude:', 1000*par[0], 'mV'
	print 'offset:', 1000*par[2], 'mV'

	plot(pulse, fit_ax)

	fit_ax.scatter(x,model_func(x,*par), marker = 'o')

	exp = model_func(x, *par)
	chisq, P = chisquare(f_obs=pulse, f_exp=exp, ddof=len(pulse)-len(par))

	print 'probability of data occuring for given parameters:', 1.0-P
	print 'Chi Square sum:', chisq

	fig3, trap_ax = plt.subplots(1,1)

	window = 100
	shp_in = raw[peaks[pk]-window:peaks[pk]+fit_length+window]
	trap = shaper_np(shp_in, 20, 10, M)

	plot(trap, trap_ax)
	### returning values with named tuples ###

	fig1.show() # plot the data
	fig2.show() # plot a pulse fit
	fig3.show() # plot trapezoidal filter of raw data

def peakdet_cwt(data, axis):
   # do a first check for peaks in the dataset. After finding peaks, should create a list of 
   # 'event' objects that will be modified as the data is further processed.
   width = np.array([1,10,20,30,40,50])
   candidates = find_peaks_cwt(data,width)	
   axis.scatter(candidates, data[candidates], marker='o', color='r', s=40)

   return candidates

def img_derivative(data, axis):
	
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



#def push(infile, pixel):

def pull(infile, pixel):
	
	# retrieve pixel signal data into 1D numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d)

	return data[pixel]

def pull_all(infile):
	
	# retrieve signal data for all 5184 pixels into 2D numpy array.
	event = 'C0' # for current dataset, only single event and sensor.
	channel = 0
	with h5py.File(infile,'r') as hf: # open file for read
		d = hf.get(event)
		data = np.array(d)

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

#def raw_plot():
