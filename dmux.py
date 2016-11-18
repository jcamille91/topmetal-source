import numpy as np
from ctypes import *

lib = CDLL("demux.so")
lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]

# infile outfile mStart mChLen mNCh
#lib.demux("../data3/dec_0x5f_20events.h5", "../src_SHLIB/dmuxOUT_CDLL.h5", c_ulong(6879), 4, c_ulong(5184), c_ulong(1), c_ulong(2), 20737)

#lib.demux("../data_TM1X1/out22.h5", "../data_TM1x1/out22_dmux.h5",  \
#	c_ulong(913), 4, c_ulong(5184), c_ulong(1), c_ulong(2), 20736)


lib.demux("../data_TM1X1/step.h5", "../data_TM1x1/step_dmux.h5",  \
	c_ulong(925), 4, c_ulong(5184), c_ulong(1), c_ulong(2), 20736)
