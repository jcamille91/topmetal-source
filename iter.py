import util as u
a=u.Sensor()
a.make_pixel(mV_noise=1, choose='exp')

# a.analyze(read = 0, select=[0], simple=0, minsep=50, threshold=0.007, 
# 	sgwin=7, sgorder=4, fudge=20, l=60, k=40, M_def=30.0, shaper_offset=-20, 
# 	bounds = ([0.008, 1.0/300, 0-0.005], [0.03, 1.0, 0+0.005]), fit_length=500, choose='sg', l_s=2, k_s=1)


a.analyze(read=0, simple=0, select=[0],
				
	sign = 1, minsep = 50, threshold = 0.006, 								# peak detection						
	sgwin = 7, sgorder = 4, l_s=4, k_s=2, choose='sg',
	pd_save=1,		  			 
    
    fit_length = 500, fudge = 20, Mrange = [10.0,45.0],						# least squares fit
    bounds = ([0.0, 1.0/300, 0-0.01], [0.03, 1.0, 0+0.01]),	   			   	# format: (min/max)
    guess = [0.008, 1.0/20, 0],						   				 		    # amplitude, 1/tau, offset
    
    l = 60, k = 20, M_def = float(30), shaper_offset = -20)					# trapezoidal filter

a.pix[0].enter_peaks(auto=1)
a.pix[0].filter_peaks(iter=True)