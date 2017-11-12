import util as u
a=u.Sensor()
a.load_pixels()
a.read_peaks('7_7_17.pkl')
a.analyze(simple=1, M_def=200.0)
a.label_events()
a.vsum_select(nfake=10000, nbins=[10,100])