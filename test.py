from util import *

def dmux_input():
	# let's input a known signal to see if the demux is working.

	# ok... was the problem the file name? 
	# doesn't really add up but seems to work now...
	fig, ax = plt.subplots(1,1)
	infile = '../data_TM1x1/dtest.h5'
	outfile = '../data_TM1x1/stepdmuxD.h5'
	pixel = 2093
	mStart = 913
	mChLen = 4
	mNCh = 72**2
	mChOff = 1
	mChSpl = 2
	frameSize = (72**2)*4

	demux(infile, outfile, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize)
	d = pull_one(outfile, pixel)
	plot(d, ax)
	fig.show()

def dmux_comp():

	# plot the old float version + the new float/double versions

	# seems likely we used a different m-start for the old dmux file, out22_dmux.h5.
	# could compare to mstart =912 instead of =913, but the data looks reasonable/physical

	fig1, ax1 = plt.subplots(1,1)
	fig2, ax2 = plt.subplots(1,1)

	pixel = 100
	infile = '../data_TM1x1/out22.h5'
	out = '../data_TM1x1/demuxdouble.h5'
	old = '../data_TM1x1/out22_dmux.h5'
	mStart = 913 # by counting backwards should be 912.....
	mChLen = 4
	mNCh = 72**2
	mChOff = 1
	mChSpl = 2
	frameSize = (72**2)*4

	demux(infile, out, mStart, mChLen, mNCh, mChOff, mChSpl, frameSize)

	d = pull_one(out, pixel)
	o = pull_one(old, pixel)

	plot(d, ax1)
	plot(o, ax2)

	fig1.show()
	fig2.show()

def send_arr(): 

	c_ulong_p = POINTER(c_ulong)
	data = np.array([743874432,55324234555,4434235689], dtype=c_ulong)
	data_p = data.ctypes.data_as(c_ulong_p)

	lib = CDLL("shaper.so")
	lib.read_arr.argtypes = [c_ulong_p]
	lib.read_arr(data_p)

def send_struct(): 

	# make some data
	nPk = c_ulong(3)
	LEFT = (np.array([1,4,77777], dtype = c_ulong))
	RIGHT = (np.array([5,333,89543], dtype = c_ulong))
	l = (np.array([444,555,8956666], dtype = c_ulong))
	k = (np.array([11,22,3333], dtype = c_ulong))
	M = (np.array([5.09,-333.67,895.44456], dtype = c_double))

 	peak = pk_hdl(nPk, LEFT, RIGHT, l, k, M)

	lib = CDLL("shaper.so")
	lib.read_struct.restype = None
	lib.read_struct.argtypes = [POINTER(peaks_t)]

	lib.read_struct(byref(peak))

def send_peaks():

	# send an array of exponentials with expected tau and filter them with different M's using
	# shaper_peaks in the shaper.c file.
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)

	# build an array with some exponentials with known time constants
	exp = np.ones(10000, dtype=c_double)*0.828
	exp[1000:2000] += 0.01*np.exp(-(1.0/40)*np.linspace(0, 999,1000))	
	exp[5000:6000] += 0.008*np.exp(-(1.0/70)*np.linspace(0, 999,1000))
	exp[8000:9000] += 0.011*np.exp(-(1.0/100)*np.linspace(0, 999,1000))

	LEFT = (np.array([0,3000,7000], dtype = c_ulong))
	RIGHT = (np.array([3000,7000,len(exp)-1], dtype = c_ulong))
	l = (np.array([100,100,100], dtype = c_ulong))
	k = (np.array([20,20,20], dtype = c_ulong))
	M = (np.array([40.0,70.0,100.0], dtype = c_double))

	PEAK = pk_hdl(len(LEFT), LEFT, RIGHT, l, k, M)
	out = np.empty_like(exp)
	lib = CDLL("shaper.so")
	lib.shaper_peaks.restype = None
	lib.shaper_peaks.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t)]
	lib.shaper_peaks(exp, out, c_ulong(len(exp)), byref(PEAK))
	axis2.axis([0, 9999,0, 0.012])
	plot(exp, axis)
	plot(out, axis2)
	fig.show()
	fig2.show()

def test_trigger():
	# sign == 1 does positive data with peaks
	# sign == -1 does negative data with valleys

	# GSL based savgol filter is broken... values aren't ridiculously wrong... but
	# don't follow trend or baseline of the data...
	
	# GSL_Sign, Yuan's convolution by FFT
	trig = namedtuple('trig', 'mean dY S ddS cds peaks toss pkm')
	sign = 1
	threshold = 0.005
	minsep = 10
	infile = '../data_TM1x1/demuxdouble.h5'
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)
	raw = pull_one(infile, 1321)
	#data = savgol_gsl(rawdata, 4, 0, 15)
	filt = savgol_scipy(raw,15,4)
	rcut = raw[9000:9300]
	fcut = filt[9000:9300]


	# first implement for positive values

	Y = filt
	plot(raw, axis)
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
		
	axis.scatter(peaks, Y[peaks], marker='x', color='r', s=40)
	axis.scatter(peaks_minsep, Y[peaks_minsep], marker='o', color='g', s=40)
	fig.show()
	#return peaks
	return trig(mean=mean, dY=dY, S=S, ddS=ddS, cds=candidates, peaks=peaks, toss=toss, pkm=peaks_minsep)


def test(pixel, pk, fit_length, l, k):

	# temporarily fixed infile because we're just using one file. insert it into arguments
	# again if we need the flexibility.

	infile = ' ../data_TM1x1/demuxdouble.h5' # new double format file.
	#infile = '../data_TM1x1/out22_dmux.h5' # old float format file, mStart wrong by 3.

	fig1, raw_ax = plt.subplots(1,1)

	raw = pull_one(infile, pixel)
	plot(raw, raw_ax) # plot the raw data

	filt = savgol_scipy(raw, 15, 4)
	#filt = savgol_gsl(raw, 4, 0, 15)
	plot(filt, raw_ax) # plot the smoothed data

	peaks = img_der(filt, raw_ax) # plot peak locations
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
	timestep = (4*72**2)*(3.2*10**-8)

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

	# l = 20
	# k = 10
	window = 100
	shp_in = raw[peaks[pk]-window:peaks[pk]+fit_length+window]
	trap = shaper_np(shp_in, l, k, M)

	plot(trap, trap_ax)
	### returning values with named tuples ###

	fig1.show() # plot the data
	fig2.show() # plot a pulse fit
	fig3.show() # plot trapezoidal filter response to pulse

	# now let's iterate through the peaks and filter them with their appropriate time constants

def test_shaper_np(l, k, M):

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

def test_savgol():
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)
	fig3, axis3 = plt.subplots(1,1)

	data = pull_one('../data_TM1x1/out22_dmux.h5', 100)
	filt1 = savgol_scipy(data, 15, 4)
	filt2 = savgol_gsl(data, 4, 0, 15)

	plot(data, axis)
	plot(filt1, axis2)
	plot(filt2, axis3)

	fig.show()
	fig2.show()
	fig3.show()

