
# this is a script to do a basic analyis.
# we start by creating a Sensor object, which contains 
# Pixel, Peak, and Event objects which we do analysis on.

# Each Pixel contains info about that specific channel, its location in the linear array, (x,y) location,
# baseline and rms values, and the raw data + filtered data associated with it.
# it also contains a list of detected Peak objects.

# Each Peak object is created during peak detection. Each peak has an associated fall time for filtering
# data. (having accurate M is very important to measurement.)

# Currently, Event objects are just labeled "by eye". Each Event object specifies when the event starts,
# and which pixels it occupies. We also use event objects to identify quiet parts of the dataset.

# we read in saved Peak data from a pickle file to save file.
# we don't have to do this, but it saves a lot of time if we don't want to repeat
# the peak detection and fitting each time (it's the slowest part of the code).
import util as u
a = u.Sensor()
a.load_pixels()
a.read_peaks('7_7_17.pkl')

# find and fit peaks

# a.analyze(read=1, M_def=20.0,threshold=0.006, l = 60, k = 20)
# a.label_events(rms_npt=50)
# a.vsum_select(show_alpha=True, show_zero=False, show_events = True, nbins=[10,40], hist_lr=[(0, 300), (-0.4,0.4)])

# use default M

# a.analyze(read=0, simple=1, M_def=20.0, l = 30, k = 20)
# a.label_events(rms_npt=50)
# a.vsum_select(show_alpha=True, show_zero=False, show_events = False, nbins=[10,40], hist_lr=[(0, 300), (-0.4,0.4)])

a.analyze(read=1, simple=0, M_def=20.0, l = 20, k = 10)
a.label_events(rms_npt=50)
a.vsum_select(v_window = (0,0), show_alpha=True, show_zero=False, \
show_events = True, nbins=[10,40], hist_lr=[(0, 300), (-0.4,0.4)])


