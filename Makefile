OSTYPE = $(shell uname)
ARCH   = $(shell uname -m)
##################################### Defaults ################################
CC             := gcc
CXX            := g++
INCLUDE        := -I.
CFLAGS         := -Wall -O2
CFLAGS_32      := -m32
SHLIB_CFLAGS   := -fPIC -shared
SHLIB_EXT      := .so
LIBS           := -lm
LDFLAGS        :=
############################# Library add-ons #################################
TINYSCHEME_FEATURES := -DUSE_DL=1 -DUSE_MATH=1 -DUSE_ASCII_NAMES=0
INCLUDE += -I/opt/local/include -I./tinyscheme/trunk
LIBS    += -L/opt/local/lib -lpthread -lhdf5
GSLLIBS  = $(shell gsl-config --libs)
CFLAGS  += -DH5_NO_DEPRECATED_SYMBOLS 
SHLIB_CFLAGS += -DH5_NO_DEPRECATED_SYMBOLS
LIBS    += $(GSLLIBS)
LIBS    += -lfftw3_threads -lfftw3f_threads -lfftw3 -lfftw3f
############################# OS & ARCH specifics #############################
ifneq ($(if $(filter Linux %BSD,$(OSTYPE)),OK), OK)
  ifeq ($(OSTYPE), Darwin)
    CC  = clang
	CXX = clang++
    SHLIB_CFLAGS   := -dynamiclib
    SHLIB_EXT      := .dylib
    TINYSCHEME_FEATURES += -DUSE_STRLWR=1 -D__APPLE__=1 -DOSX=1
    ifeq ($(shell sysctl -n hw.optional.x86_64), 1)
      ARCH           := x86_64
    endif
  else
    ifeq ($(OSTYPE), SunOS)
      CFLAGS         := -c -Wall -std=c99 -pedantic
    else
      # Let's assume this is win32
      SHLIB_EXT      := .dll
      TINYSCHEME_FEATURES += -DUSE_STRLWR=0
    endif
  endif
else
  TINYSCHEME_FEATURES += -DSUN_DL=1
endif

ifneq ($(ARCH), x86_64)
  CFLAGS_32 += -m32
endif

# Are all G5s ppc970s?
ifeq ($(ARCH), ppc970)
  CFLAGS += -m64
endif
############################ Define targets ###################################
#EXE_TARGETS = demux shaper
DEBUG_EXE_TARGETS = hdf5rawWaveformIo
SHLIB_TARGETS = H5SHLIB.so demux.so 
#demux_dbl.so demux_flt.so
#demux.so shaper.so filters.so
O_TARGETS = hdf5rawWaveformIo.o demux.o
#demux_dbl.o demux_flt.o
#demux.o shaper.o filters.o

ifeq ($(ARCH), x86_64) # compile a 32bit version on 64bit platforms
  # SHLIB_TARGETS += XXX_m32$(SHLIB_EXT)
endif

.PHONY: shlib_targets debug_exe_targets o_targets clean
#.PHONY: exe_targets shlib_targets debug_exe_targets o_targets clean
#exe_targets: $(EXE_TARGETS)
shlib_targets: $(SHLIB_TARGETS)
debug_exe_targets: $(DEBUG_EXE_TARGETS)
o_targets: $(O_TARGETS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@
%.o: %.C
	$(CXX) $(CFLAGS) $(ROOTCFLAGS) $(INCLUDE) -c $< -o $@
	$(CXX) $(CFLAGS) $(ROOTCFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
hdf5rawWaveformIo.o: hdf5rawWaveformIo.c hdf5rawWaveformIo.h common.h
hdf5rawWaveformIo: hdf5rawWaveformIo.c hdf5rawWaveformIo.h
	$(CC) $(CFLAGS) $(INCLUDE) -DHDF5RAWWAVEFORMIO_DEBUG_ENABLEMAIN $< $(LIBS) $(LDFLAGS) -o $@
demux_flt.o: demux_flt.c hdf5rawWaveformIo.o
demux_dbl.o: demux_dbl.c hdf5rawWaveformIo.o
demux.o: demux.c hdf5rawWaveformIo.o
demux: demux.c hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
shaper.o: shaper.c hdf5rawWaveformIo.o
shaper: shaper.c hdf5rawWaveformIo.o
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LIBS) $(LDFLAGS) -o $@
filters.o: filters.c hdf5rawWaveformIo.o
filters.so: filters.o hdf5rawWaveformIo.o
	$(CC) $^ $(SHLIB_CFLAGS) $(INCLUDE) -o filters.so $(LIBS)
H5SHLIB.so: hdf5rawWaveformIo.o
	$(CC) hdf5rawWaveformIo.o $(SHLIB_CFLAGS) $(INCLUDE) -o H5SHLIB.so $(LIBS)
demux_flt.so: demux_flt.o hdf5rawWaveformIo.o
	$(CC) $^ $(SHLIB_CFLAGS) $(INCLUDE) -o demux_flt.so $(LIBS)	
demux_dbl.so: demux_dbl.o hdf5rawWaveformIo.o
	$(CC) $^ $(SHLIB_CFLAGS) $(INCLUDE) -o demux_dbl.so $(LIBS)	
demux.so: demux.o hdf5rawWaveformIo.o
	$(CC) $^ $(SHLIB_CFLAGS) $(INCLUDE) -o demux.so $(LIBS)
shaper.so: shaper.o hdf5rawWaveformIo.o
	$(CC) $^ $(SHLIB_CFLAGS) $(INCLUDE) -o shaper.so $(LIBS)


clean:
	rm -f *.so *.dylib *.dll *.bundle
	rm -f $(SHLIB_TARGETS) $(EXE_TARGETS) $(DEBUG_EXE_TARGETS) $(O_TARGETS)
