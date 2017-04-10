from util import *

### let's examine some busy channels. not as worried about sparse channels.

### keep track of peaks discarded

### for filtering, if peaks are closer than fit length, than filter both peaks
### with only the second peak's fitted M.

### channels for current run with instable baseline. baseline shift significant 
# compare to the desired signal. possible charge trapping inside of CSA.
BAD  = np.array([268,269,340,718,719,805,806,1047,1048,1617,1618,2037,20
    ...: 38, 2188,2189,2260,2396,2397, 2539, 2540, 2768, 2769, 4868])
EXTRA = np.array([1898, 1899, 1970])


def sensor_noise(npt):
	
	# ch 0-2 are dead channels for identifying frame start in demux algorithm.
	# npt = 25890
	nch = 72**2
	dead = 3
	infile = '../data_TM1x1/demuxdouble.h5'

	# bounds for separating channels into three categories by noise.
	x1 = 0.0012
	x2 = 0.002

	a = get_wfm_all(infile, npt) # wfm = namedtuple('wfm', 'avg rms data')
	
	quiet_ch = np.where(a.rms < med)[0]
	mid_ch = np.where((a.rms <= loud) & (a.rms >= med))[0]
	noisy_ch = np.where(a.rms > loud)[0]


	print('max rms', np.amax(a.rms[dead:nch-1]))
	print('min rms', np.amin(a.rms[dead:nch-1]))
	print('avg rms', np.mean(a.rms[dead:nch-1]))
	
	fig = plt.figure(1)
	ax = fig.add_subplot(111)

	# plot a 1D histogram of different noise levels observed (Volts RMS)
	nbin= 2000
	end_pt = 0.002 # Volts RMS
	#np.histogram(a.rms, nbin)
	hist_plot(a.rms, nbin, end_pt, ax)
	print('noisy_ch are channels above', x2, 'volts RMS')
	print('mid_ch are channels between', x1, 'and', x2, 'volts RMS')
	fig.show()
	#print(noisy_ch)

	return (a, noisy_ch, mid_ch, quiet_ch)

def compare_noise_npks():
	mvnoise = 0.002
	sn = sensor_noise(25890) #2 mV RMS threshold for noise
	noisy_ch = sn[1]

	pk = check_peaks(0.007, 40, 15, 4)
	busy_ch = pk.many

	ch = np.intersect1d(busy_ch, noisy_ch)
	print('# of noisy channels with alot of peaks:', len(ch))

	return sn[0]
	
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

def send_peaks(l, k, res, fit_length, mV_noise):

	# clear figures from previous run.
	close_figs()

	# send an array of exponentials with expected tau and filter them with different M's using
	# shaper_peaks in the shaper.c file.
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)

	# specify M=tau, peak location, amplitude, and the number of points to define
	# the peak.
	M_i = np.array([30, 20, 40], dtype=np.float64)
	M_loc = np.array([1000, 5000, 8000])
	M_a = np.array([0.02, 0.015, 0.011])
	npk = len(M_i)
	M_npt = res
	baseline = 0.8

	# parameters for peak detection.

	sgorder = 4
	sgwin = 15
	minsep = 20

	# offset for L/R of trapezoidal shaper
	shaper_offset = 0

	# make array big enough to accomodate peak location and size.
	data_len = 10000

	# apply rms noise to signal (input is in mV). if input is zero, no noise desired.
	if mV_noise == 0 :
		noise = 0
	else :
		noise = np.random.normal(loc = 0, scale = mV_noise*.001, size=data_len)

	# build an array with exponential decays just defined above.
	exp = np.ones(data_len, dtype=np.float64)*baseline
	for i in range(npk) :
		exp[M_loc[i]:M_loc[i]+M_npt] += M_a[i]*np.exp(-(1.0/M_i[i])*np.linspace(0, M_npt-1, M_npt)) 
	exp += noise

	# get the peaks. either use the locations defined above, or use get_peaks().
	peaks = M_loc
	# peaks = get_peaks(exp, baseline, mV_noise*0.001, minsep, sgwin, sgorder)

	# fit time constants for pulses
	tau = fit_tau(exp, baseline, mV_noise*0.001, peaks, 0, fit_length, axis)


	# trapezoidal filter transitions currently set to edges of pulses.
	LR = pk2LR(peaks, shaper_offset, len(exp))
	l_arr = (l*np.ones(npk, dtype = c_ulong))
	k_arr = (k*np.ones(npk, dtype = c_ulong))
	M = M_i
	#M = tau[0]

	PEAK = peaks_handle(len(LR[0]), LR[0], LR[1], l_arr, k_arr, M)
	out = np.empty_like(exp)
	lib = CDLL("shaper.so")
	lib.shaper_multi.restype = None
	lib.shaper_multi.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t), c_double]
	lib.shaper_multi(exp, out, c_ulong(len(exp)), byref(PEAK), c_double(baseline))
	
	plot(exp, axis)
	plot(out, axis2)
	fig.show()
	fig2.show()
	print('True time constants: \n')
	for i in range(npk) :
		print('Peak no.', i+1)
		print(M_i[i])
	print('\n')
	print('Fitted time constants: \n')
	for i in range(npk) :
		print('Peak no.', i+1)
		print('tau = ', tau[0][i])
		print('chi-square = ', tau[1][i])
		print('Q = ', tau[2][i])
		print('P = ', tau[3][i])
		print('\n')
	#return (exp , np.array(out))

def do_simple() :

	npix = 72**2
	infile = '../data_TM1x1/demuxdouble.h5'
	d = get_wfm_all(infile, 25890)
	data = d.data
	avg = d.avg
	rms = d.rms
	threshold = 0.007
	minsep = 40
	sgwin = 15
	sgorder = 4
	fudge = 0 
	l = 500
	k = 50
	default_M = 40
	filt = np.empty_like(data)

	for i in range(npix) :
		#filt[i] = shaper_single(data[i], l, k, default_M, data[i][0])
		filt[i] = shaper_single(data[i], l, k, default_M, avg[i])
	return filt

def check_peaks(threshold, minsep, sgwin, sgorder) :

	''' 
	inputs: 
	- threshold : minimum value to qualify a peak.
	- minsep : peaks closer than minsep are rejected.
	- sgwin : window of savitsky-golay filter. smaller window picks up smaller features
	and vice-versa. window is number of adjacent points for polynomial fitting.
	- sgorder : order of savitsky-golay polynomial for fitting to data.
	
	outputs:
	- none : no peaks found
	- few : less than 5 peaks
	- avg : 5 < peaks < 15
	- many : 15 < peaks
	'''

	dead = 3
	npix = 72**2
	infile = '../data_TM1x1/demuxdouble.h5'
	d = get_wfm_all(infile, 25890)
	data = d.data
	avg = d.avg
	rms = d.rms

	# some typical parameters for peak detection
	# threshold = 0.007
	# minsep = 40
	# sgwin = 15
	# sgorder = 4

	few_val = 5
	mid_val = 15


	none = np.array([], dtype = int)
	few = np.array([], dtype = int)
	mid = np.array([], dtype = int)
	many = np.array([], dtype = int)
	
	pkch = np.zeros(72**2, dtype=int) # array storing number of found peaks in each channel.

	for i in range(dead, npix) :

		peaks = get_peaks(data[i], avg[i], threshold, minsep, sgwin, sgorder)
		npk = len(peaks)
		pkch[i] = npk
		if npk == 0 :
			none = np.append(none, i)
		elif npk < few_val :
			few = np.append(few, i)
		elif npk < mid_val :
			mid = np.append(mid, i)
		else  :
			many = np.append(many, i)

	print('# of channels with 0 peaks:', len(none))
	print('# of channels with <', few_val, 'peaks:', len(few))
	print('# of channels with <', mid_val, 'peaks:', len(mid))
	print('# of channels with >', mid_val, 'peaks:', len(many))



	# ch_status = np.ones(72**2)*0.5
	# ch_status[busy_ch] = 1
	# ch_status[none_ch] = 0

	# pixel_status(ch_status)
	pixel_status(pkch)


	return PEAKS(none=none, few=few, mid=mid, many=many, pkch=pkch)


def do_all() :

	### let's make a plot of zero peaks, some peaks, and a lot of peaks/

	### should also do a test quantitatively comparing
	### results of sensor noise and channel status.

	dead = 3
	npix = 72**2
	infile = '../data_TM1x1/demuxdouble.h5'
	d = get_wfm_all(infile, 25890)
	data = d.data
	avg = d.avg
	rms = d.rms
	noise_thresh = 0.002 # 2mV RMS
	pk_thresh = 0.007   # 7mV amplitude
	minsep = 40
	sgwin = 15
	sgorder = 4
	fudge = 0 
	fit_length = 150
	ax = 0 # not plotting fits
	shaper_offset = 5
	l = 500
	k = 50
	default_M = 30.0

	filt = np.empty_like(data)

	max_npk = 20
	
	busy_ch = np.array([], dtype=int) # channels with number of peaks exceeding max_npk
	none_ch = np.array([], dtype=int) # channels with no found peaks
	peak_ch = np.zeros(72**2, dtype=int) # array storing number of found peaks in each channel.


	for l in range(dead) : # filter first three dead channels. no peaks to find here.

		filt[l] = shaper_single(data[l], l, k, default_M, avg[l])

	for i in range(dead,npix) :

		peaks = get_peaks(data[i], avg[i], pk_thresh, minsep, sgwin, sgorder)
		npk = len(peaks)
		peak_ch[i] = npk

		if npk > max_npk : # too many peaks
			#print('ch', i, 'has', len(peaks), 'peaks') 
			# busy_ch = np.append(busy_ch, i)
			filt[i] = shaper_single(data[i], l, k, default_M, avg[i])

		elif (npk > 1 and npk <= max_npk) : # multiple peaks on a channel
			M = fit_tau(data[i], avg[i], rms[i], peaks, fudge, fit_length, ax)
			filt[i] = shaper_multi(data[i], peaks, l, k, M, shaper_offset, avg[i])

		elif npk == 1 : # single peak on a channel
			M = fit_tau(data[i], avg[i], rms[i], peaks, fudge, fit_length, ax)
			filt[i] = shaper_single(data[i], l, k, M, avg[i])

		elif npk == 0 : # no peaks found
			#print('ch', i, 'has 0 peaks')
			# none_ch = np.append(none_ch, i)
			filt[i] = shaper_single(data[i], l, k, default_M, avg[i])

	# ch_status = np.ones(72**2)*0.5
	# ch_status[busy_ch] = 1
	# ch_status[none_ch] = 0

	# pixel_status(ch_status)
	# ps2(peak_ch)
	return (filt, peak_ch)


def do_one(ch, threshold, Msingle, fit_length) :

	close_figs()
	infile = '../data_TM1x1/demuxdouble.h5'

	# fudge(fit_tau), minsep(get peaks), and off(shaper peaks) 
	# can cause problems, need to make them work together

	# build an array with some exponentials with known time constants
	# exp = np.ones(10000, dtype=c_double)*0.828
	# exp[1000:2000] += 0.01*np.exp(-(1.0/40)*np.linspace(0, 999,1000))	
	# exp[5000:6000] += 0.008*np.exp(-(1.0/70)*np.linspace(0, 999,1000))
	# exp[8000:9000] += 0.011*np.exp(-(1.0/100)*np.linspace(0, 999,1000))
	
	# minsep must be greater than fudge

	l = 200
	k = 50
	#threshold = 
	fudge = 0
	shaper_offset = -30
	minsep = 50
	sgwin = 15
	sgorder = 4
	#fit_length = 300  # about 3x the expected tau

	# exp = np.ones(10000, dtype=np.float64)*0.828
	# exp[1000:2000] += 0.01*np.exp(-(1.0/10)*np.arange(1000))	
	# exp[5000:6000] += 0.008*np.exp(-(1.0/30)*np.arange(1000))
	# exp[8000:9000] += 0.011*np.exp(-(1.0/45)*np.arange(1000))
 # 	noise = 0
 # 	exp = exp + noise

 	# get waveform and find peaks.
	d = get_wfm_one(file=infile, ch=ch, npt=25889, plt=False)
	data = d.data
	baseline = d.avg
	rms = d.rms
	peaks = get_peaks(data, baseline, threshold, minsep, sgwin, sgorder)
	npk = len(peaks)

	print('channel #', ch, 'baseline = ', baseline, 'Vrms = ', rms)
	if npk > 0 :

		# plot peak locations
		fig, ax = plt.subplots(1,1)
		plot(data, ax)
		ax.scatter(peaks, data[peaks], marker='x', color='r', s=40)
		fig.show()
		
		# fit tau to peaks with least squares
		tau = fit_tau(data, baseline, rms, peaks, fudge, fit_length, ax)
		if npk > 1 :
			filt = shaper_multi(data, peaks, l, k, tau[0], shaper_offset, baseline)
		else :
			filt = shaper_single(data, l, k, tau[0], baseline)
		fig2, ax2 = plt.subplots(1,1)
		plot(filt, ax2)
		ax2.scatter(peaks, filt[peaks], marker='x', color='r', s=40)
		fig2.show()

		# print another plot with just a single M value used for the whole data set for comparison.
		fig3, ax3 = plt.subplots(1,1)
		filt3 = shaper_single(data, l, k, Msingle, baseline)
		plot(filt3, ax3)
		ax3.scatter(peaks, filt3[peaks], marker='x', color='r', s=40)
		fig3.show()

		print('Fitted time constants:')
		for i in range(npk) :
			print('\n')
			print('Peak no.', i+1)
			print('tau = ', tau[0][i]) 
			print('chi-square = ', tau[1][i])
			print('Q = ', tau[2][i])
			print('P = ', tau[3][i])
			

		return (peaks, tau[0])

	else :
		print('no peaks found in Channel', ch)

def trigger(channel, sgwin, threshold):
	"""
	namedtuple('trig', 'mean dY S ddS cds peaks toss pkm')
	"""


	sign = 1
	#threshold = 0.004
	minsep = 5
	infile = '../data_TM1x1/demuxdouble.h5'
	#sgwin = 11
	sgorder = 4

	fig, axis = plt.subplots(1,1)
	d = get_wfm_one(infile, channel, 25890, False)
	trig = get_peaks(d.data, d.avg, threshold, minsep, sgwin, sgorder)
	filt = savgol_scipy(d.data, sgwin, sgorder)
	fig, ax = plt.subplots(1,1)
	plot(d.data, ax)
	plot(filt, ax)
	ax.scatter(trig.pkm, d.data[trig.pkm], marker='x', color='r', s=40)
	fig.show()

	return trig

def old(pixel, pk, fit_length, l, k):

	# for fitting and filtering, do a conditional:

	# fit / filter either max value, or the distance between peaks if it is less
	# than the max value. these two possibilities might be sufficient.

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

	print('tau:', M, 'samples', M*timestep, 'seconds')
	print('amplitude:', 1000*par[0], 'mV')
	print('offset:', 1000*par[2], 'mV')

	plot(pulse, fit_ax)

	fit_ax.scatter(x,model_func(x,*par), marker = 'o')

	exp = model_func(x, *par)
	chisq, P = chisquare(f_obs=pulse, f_exp=exp, ddof=len(pulse)-len(par))

	print('probability of data occuring for given parameters:', 1.0-P)
	print('Chi Square sum:', chisq)

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

def shaper_test(l, k, M):

	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)
	fig3, axis3 = plt.subplots(1,1)


	step = np.ones(2000, dtype=np.float64) 
	step[:1000] = 0.828 
	step[1000:] = 0.978

	exp = np.ones(4000, dtype=np.float64)*0.828
	exp[1000:2000] += 0.01*np.exp(-(1/40.)*np.linspace(0, 999,1000))
	
	filt = shaper_single(exp, l, k, M)
	#filt = shaper_np(step, 20, 10, -1)

	plot(exp, axis)
	#plot(step, axis2)
	plot(filt, axis3)
	

	fig.show()
	#fig2.show()
	fig3.show()

def savgol_compare():
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

