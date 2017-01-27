import util as u

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
	 
 	# make pointers for ctypes
	LEFT_p = LEFT.ctypes.data_as(c_ulong_p)
	RIGHT_p = RIGHT.ctypes.data_as(c_ulong_p)
	l_p = l.ctypes.data_as(c_ulong_p)
	k_p = k.ctypes.data_as(c_ulong_p)
	M_p = M.ctypes.data_as(c_double_p)

 	peak = u.peaks_t(nPk, LEFT_p, RIGHT_p, l_p, k_p, M_p)

	lib = CDLL("shaper.so")
	lib.read_struct.restype = None
 	# 					  # 	  in 		out   	 length  	l  		k 		  M
	lib.read_struct.argtypes = [POINTER(peaks_t)]

	lib.read_struct(byref(peak))

def send_peaks():

	# send an array of exponentials with expected tau and filter them with different M's using
	# shaper_peaks in the shaper.c file.
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)

	# build an array with some exponentials with known time constants
	exp = np.ones(100000, dtype=c_double)*0.828
	exp[10000:20000] += 0.01*np.exp(-(1/40.)*np.linspace(0, 9999,10000))	
	exp[50000:60000] += 0.008*np.exp(-(1/70.)*np.linspace(0, 9999,10000))
	exp[80000:90000] += 0.011*np.exp(-(1/100.)*np.linspace(0, 9999,10000))

	LEFT = (np.array([0,30000,70000], dtype = c_ulong))
	RIGHT = (np.array([30000,70000,len(exp)-1], dtype = c_ulong))
	l = (np.array([25,25,25], dtype = c_ulong))
	k = (np.array([10,10,10], dtype = c_ulong))
	M = (np.array([40.0,70.0,100.0], dtype = c_double))

	PEAK = u.pk_hdl(len(LEFT), LEFT, RIGHT, l, k, M)
	out = np.empty_like(exp)
	lib = CDLL("shaper.so")
	lib.shaper_peaks.restype = None
	lib.shaper_peaks.argtypes = [double_ptr, double_ptr, c_ulong, POINTER(peaks_t)]
	lib.shaper_peaks(exp, out, c_ulong(len(exp)), byref(PEAK))
	axis2.axis([0, 99999,0, 0.012])
	u.plot(exp, axis)
	u.plot(out, axis2)
	fig.show()
	fig2.show()

def test(pixel, pk, fit_length, l, k):

	# temporarily fixed infile because we're just using one file. insert it into arguments
	# again if we need the flexibility.
	infile = '../data_TM1x1/out22_dmux.h5'

	fig1, raw_ax = plt.subplots(1,1)

	raw = u.pull_one(infile, pixel)
	u.plot(raw, raw_ax) # plot the raw data

	filt = u.savgol_scipy(raw, 15, 4)
	#filt = savgol_gsl(raw, 4, 0, 15)
	u.plot(filt, raw_ax) # plot the smoothed data

	peaks = u.img_der(filt, raw_ax) # plot peak locations
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

	u.plot(pulse, fit_ax)

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
	trap = u.shaper_np(shp_in, l, k, M)

	u.plot(trap, trap_ax)
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
	
	filt = u.shaper_np(exp, l, k, M)
	#filt = shaper_np(step, 20, 10, -1)

	u.plot(exp, axis)
	#plot(step, axis2)
	u.plot(filt, axis3)
	

	fig.show()
	#fig2.show()
	fig3.show()

def test_savgol():
	fig, axis = plt.subplots(1,1)
	fig2, axis2 = plt.subplots(1,1)
	fig3, axis3 = plt.subplots(1,1)

	data = u.pull_one('../data_TM1x1/out22_dmux.h5', 100)
	filt1 = u.savgol_scipy(data, 15, 4)
	filt2 = u.savgol_gsl(data, 4, 0, 15)

	u.plot(data, axis)
	u.plot(filt1, axis2)
	u.plot(filt2, axis3)

	fig.show()
	fig2.show()
	fig3.show()

