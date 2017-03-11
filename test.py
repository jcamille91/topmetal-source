from util import *

def sensor_noise(npt, threshold):
	
	noise = namedtuple('noise', 'data rms thresh')

	# ch 0-2 are dead channels for demux
	dead = 3
	infile = '../data_TM1x1/demuxdouble.h5'
	a = pull_all(infile)
	length = len(a[0])

	if ((npt == False) or (npt > length)) :
		print 'set calculation length to raw data array length =', length 
		npt = length

	rms = np.zeros(72**2-dead) # this array only goes from 0-5180, exlude first three pixels for calc.
	for i in np.arange(72**2-dead)+dead :
		rms[i-3] = np.std(a[i][:npt])

	thresh_i = np.where(rms > threshold)[0] + dead
	
	print 'max rms', np.amax(rms)
	print 'min rms', np.amin(rms)
	print 'avg rms', np.mean(rms)
	print 'channels above', threshold, 'volts RMS' 
	print thresh_i

	return noise(data=a, rms=rms, thresh=thresh_i)



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

def send_peaks(l, k, res, fit_length):

			### TEST THIS NEXT ###

	# test a limiting case for filters working on the edge:
	# see if the response looks good in the presence of noise / shifting baseline
	# if we switch l,k,M right at the beginning of a peak. if so, this could
	# simplify specifying different lengths of filtering for quick successive peaks.
	
	# this works doing it right on the transition.

			######################

	# send an array of exponentials with expected tau and filter them with different M's using
	# shaper_peaks in the shaper.c file.
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)

	M_i = np.array([30, 20, 10], dtype=np.float64)
	M_f = np.array([30.174071, 20.090518, 10.859501], dtype=np.float64)
	M_loc = np.array([1, 5000, 8000])
	M_a = np.array([0.01, 0.008, 0.011])
	npk = len(M_i)
	M_npt = res
	#M_npt = 100

	data_len = 10000
	noise = 0
	baseline = 0.0
	# build an array with some exponentials with known time constants
	exp = np.ones(data_len, dtype=np.float64)*baseline
	for i in np.arange(npk) :
		exp[M_loc[i]:M_loc[i]+M_npt] += M_a[i]*np.exp(-(1.0/M_i[i])*np.linspace(0, M_npt-1, M_npt)) 
	# exp[1000:2000] += 0.01*np.exp(-(1.0/M_i[0])*np.linspace(0, 999,1000))	
	# exp[5000:6000] += 0.008*np.exp(-(1.0/M_i[1])*np.linspace(0, 999,1000))
	# exp[8000:9000] += 0.011*np.exp(-(1.0/M_i[2])*np.linspace(0, 999,1000))
	peaks = get_peaks(exp, 0.828, 0, 50, 15, 4)
	tau = fit_tau(exp, peaks, 0, fit_length, 0)

	print 'peak found at:', peaks
	print 'corresponding m:', tau
	exp += noise
	LEFT = (np.array([0,5000,8000], dtype = c_ulong))
	RIGHT = (np.array([5000,8000,len(exp)], dtype = c_ulong))
	l_arr = (l*np.ones(npk, dtype = c_ulong))
	k_arr = (k*np.ones(npk, dtype = c_ulong))
	# l = (np.array([200,200,200], dtype = c_ulong))
	# k = (np.array([100,100,100], dtype = c_ulong))
	#M = (np.array(input, dtype = c_double))
	#M = (np.array([41.0,71.0,101.0], dtype = c_double))
	M = M_f
	PEAK = peaks_handle(len(LEFT), LEFT, RIGHT, l_arr, k_arr, M)
	out = np.empty_like(exp)
	lib = CDLL("shaper.so")
	lib.shaper_peaks.restype = None
	lib.shaper_peaks.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t)]
	lib.shaper_peaks(exp, out, c_ulong(len(exp)), byref(PEAK))
	#axis2.axis([0, 9999,0, 0.012])
	plot(exp, axis)
	plot(out, axis2)
	fig.show()
	fig2.show()
	return (exp , np.array(out))

def test_all() :
	
	infile = '../data_TM1x1/demuxdouble.h5'
	d = get_wfm_all(infile)


def test_one(ch) :

	infile = '../data_TM1x1/demuxdouble.h5'

	# fudge(fit_tau), minsep(get peaks), and off(shaper peaks) 
	# can cause problems, need to make them work together

	# build an array with some exponentials with known time constants
	# exp = np.ones(10000, dtype=c_double)*0.828
	# exp[1000:2000] += 0.01*np.exp(-(1.0/40)*np.linspace(0, 999,1000))	
	# exp[5000:6000] += 0.008*np.exp(-(1.0/70)*np.linspace(0, 999,1000))
	# exp[8000:9000] += 0.011*np.exp(-(1.0/100)*np.linspace(0, 999,1000))
	
	# minsep must be greater than fudge

	l = 1000
	k = 100
	threshold = 0.006
	fudge = 0
	shape_offset = c_ulong(5)
	minsep = 10
	fit_length = 100  # no more than about 3x the expected tau

	# exp = np.ones(10000, dtype=np.float64)*0.828
	# exp[1000:2000] += 0.01*np.exp(-(1.0/10)*np.arange(1000))	
	# exp[5000:6000] += 0.008*np.exp(-(1.0/30)*np.arange(1000))
	# exp[8000:9000] += 0.011*np.exp(-(1.0/45)*np.arange(1000))
 # 	noise = 0
 # 	exp = exp + noise

 	d = get_wfm(file=infile, ch=ch, npt=25889, plt=False)
 	data = d.data
 	threshold = 4*d.rms
	peaks = get_peaks(data, d.avg, threshold, minsep, 15, 4)

	if len(peaks) > 0 :
		fig, ax = plt.subplots(1,1)
		plot(data, ax)
		ax.scatter(peaks, data[peaks], marker='x', color='r', s=40)

		M = fit_tau(data, peaks, fudge, fit_length, ax)

		fig.show()
	
		filt = shaper_peaks(data, peaks, l, k, M, shape_offset)
		fig2, ax2 = plt.subplots(1,1)
		plot(filt, ax2)
		ax2.scatter(peaks, filt[peaks], marker='x', color='r', s=40)
		fig2.show()

		return (peaks, M)

	else :
		print 'No peaks found in this channel.'

	#return (peaks, M)

def trigger():

	trig = namedtuple('trig', 'mean dY S ddS cds peaks toss pkm')
	sign = 1
	threshold = 0.006
	minsep = 5
	infile = '../data_TM1x1/demuxdouble.h5'
	sgwin = 15
	sgorder = 4

	fig, axis = plt.subplots(1,1)
	d = get_wfm(infile, 1321, 25890, False)
	peaks = get_peaks(d.data, d.avg, threshold, minsep, sgwin, sgorder)

	fig, ax = plt.subplots(1,1)
	plot(d.data, ax)
	ax.scatter(peaks, d.data[peaks], marker='x', color='r', s=40)
	fig.show()


def oldtest(pixel, pk, fit_length, l, k):

	# for fitting and filtering, do a conditional:

	# fit / filter either max value, or the distance between peaks if it is less
	# than the max value. i think these two possibilities should be sufficient.

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

def shaper_np_test(l, k, M):

	### TEST THIS NEXT ###

	# test a limiting case for filters working on the edge:
	# see if the response looks good in the presence of noise / shifting baseline
	# if we switch l,k,M right at the beginning of a peak. if so, this could
	# simplify specifying different lengths of filtering for quick successive peaks.

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

def savgol_test():
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

