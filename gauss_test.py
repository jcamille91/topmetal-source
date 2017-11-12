# simple test to input a noise signal with known variance into each channel, then sum variances of multiple channels,
# and recover the correct RMS value per channel. result: it worked.

# running test with gauss_test=False, we just measure the RMS of selection of channels by summing variances.
# this gave a reasonable value of 0.8mV RMS.

### if we want to exclude some selections,they can simply be commented out of the 'label_events' method.
### if we want a very specific cut, we can implement logic in the 'rms_calc' function that discriminates
### against events with certain characteristics.

import util as u
a=u.Sensor()
a.load_pixels(npt=25890, gauss_test=False, sigma=0.7, mean=2)

# sample time is 663.552 us, so 100 points is ~66 ms to calculate the RMS value.
# in the future we could break analyze into two functions: analyze() and set_parameters()
# incase we don't always want to filter stuff and find peaks. like in this instance: we only want to 
# do stuff to the raw data, but the analysis parameters are all set in the analyze() method
# so we have to filter data to set parameters even though we don't need the data filtered..

a.analyze(read=0, simple=1, select=[],
				
	sign = 1, minsep = 50, threshold = 0.008, 								# peak detection						
	sgwin = 7, sgorder = 4, l_s=4, k_s=2, choose='sg',
	pd_save=0,		  			 
    
    fit_length = 300, fudge = 20, Mrange = [10.0,50.0],						# least squares fit
    bounds = ([0.0, 1.0/300, 0-0.01], [0.03, 1.0, 0+0.01]),	   			   	# format: (min/max)
    guess = [0.008, 1.0/20, 0],						   				 		# amplitude, 1/tau, offset
    
    l = 60, k = 20, M_def = float(20), shaper_offset = -20)

# the expected center is at sigma * sqrt(npix)

a.label_events(rms_npt=50)

a.rms_calc(nbins=50, hist_lr=[0.012,0.015], threshold=.0132)

for i in a.noisyevents :
	print('x=%i, y=%i, i=%i, rms=%f'% (i.x, i.y, i.i, i.noisyrms)) 