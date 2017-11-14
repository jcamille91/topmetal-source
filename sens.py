import util as u
a=u.Sensor()


a.load_pixels()
a.analyze(read=0, simple=1, select=[],
				
	sign = 1, minsep = 50, threshold = 0.008, 								# peak detection						
	sgwin = 7, sgorder = 4, l_s=4, k_s=2, choose='sg',
	pd_save=0,		  			 
    
    fit_length = 300, fudge = 20, Mrange = [10.0,50.0],						# least squares fit
    bounds = ([0.0, 1.0/300, 0-0.01], [0.03, 1.0, 0+0.01]),	   			   	# format: (min/max)
    guess = [0.008, 1.0/20, 0],						   				 		    # amplitude, 1/tau, offset
    
    l = 20, k = 10, M_def = float(20), shaper_offset = -20)					# trapezoidal filter




# a.make_pixel(mV_noise=1, choose='exp')

# a.analyze(read=0, simple=0, select=[0],
				
# 	sign = 1, minsep = 50, threshold = 0.006, 								# peak detection						
# 	sgwin = 7, sgorder = 4, l_s=4, k_s=2, choose='sg',
# 	pd_save=0,		  			 
    
#     fit_length = 500, fudge = 20, Mrange = [11.0,45.0],						# least squares fit
#     bounds = ([0.0, 1.0/300, 0-0.01], [0.03, 1.0, 0+0.01]),	   			   	# format: (min/max)
#     guess = [0.008, 1.0/20, 0],						   				 		    # amplitude, 1/tau, offset
    
#     l = 60, k = 20, M_def = float(30), shaper_offset = -20)	


# return breaks out of nested loop as expected, the outer loop is broken when i, j = 20
def test():
	for i in range(40):
		print(i)
		for j in range(40):
			print(j)
			if (i == 20 and j == 20):
				return (i,j)