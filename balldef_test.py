import util as u
a=u.Sensor()
a.make_pixel(choose='ball', mV_noise=0)
a.analyze(select=[0], simple=True, l=20, k=10, M_def=float(17.7) )
a.plot_waveform(pixels=[0], lr=(50000-2000, 50000+2000))