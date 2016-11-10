# topmetal-source TM-1x1 branch (starting nov 9 2016)

The 1x8 detector ADC is faulty, so going to configure the software for some older 1x1 sensor data.

changes: 

1. this detector's ADC has a different topology, so the ADC dataformat is a signed int, as opposed to the 1x8 ADC's unsigned int. This will require small changes in the hdf5rawWaveformIo source and header files.

2. The ADC waveform attribute arrays for converting to volts from samples have 16 elements for the 1x1 file, instead of the 8 element arrays, so account for this. These attributes are yoff, yzero, and ymult.

-> change SCOPE_NCH to 16, this should fix the arrays cleanly in the hdf5 read_waveform_attribute_in_file_header function.

3. The onboard memory used for this data acquisition was different than the 1x8 detector.

-> change SCOPE_MEM_LENGTH_MAX

