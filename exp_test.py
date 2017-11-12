import util as u
a=u.Sensor()
a.make_pixel(choose='exp', mV_noise=0)
a.analyze(select=[0], simple=True, l=27, k=20, M_def=float(200) )
a.plot_waveform(pixels=[0], lr=(50000-50, 50000+2000))