# some useful functions for testing things that are used frequently.

import numpy as np
import h5py 
from ctypes import *
import matplotlib.pyplot as plt

 # demux and pre_shaper are pre-processing functions, applied to all data.
def demux(infile, outfile, mStart, mChLen mNCh, mChOff, mChSpl, frameSize):
   lib = CDLL("demux.so")
   lib.demux.argtypes = [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
   lib.demux(c_char_p(infile), c_char_p(outfile), c_ulong(mStart), c_double(mChLen), c_ulong(mNCh), c_ulong(mChOff), c_ulong(mChSpl), c_double(frameSize))

def pre_shaper(infile, outfile, l, k, M):
   # apply preliminary trapezoidal filter to the entire dataset so we can search for peaks.
   lib = CDLL("shaper.so")
   lib.shaper.argtypes = [c_char_p, c_char_p, c_ulong, c_ulong, c_double]
   lib.shaper(c_char_p(infile), c_char_p(outfile), c_ulong(l), c_ulong(k), c_double(M))

def event_shaper(infile, pixel, evt_pt, tau):
   # apply trapezoidal filter to data on event-by-event basis
   # extract the pulse height and width.
def pull_data(infile, pixel):

   return array

def find_peaks():

    return array

def find_parameters(infile, peaks, )
def img_plot(infile, point):

def signal_plot():

def raw_plot():
