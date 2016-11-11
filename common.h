/** \file common.h
 * Definitions of types, parameters and utilities commonly referenced to.
 */
#ifndef __COMMON_H__
#define __COMMON_H__

#define HDF5IO(name) hdf5io_ ## name

/* filling memory per channel:

1x8 sensor: 256 MB available. shared by 8ch. each point is 2 bytes (unsigned).
so nPt = 256*(2^20)/(8*2) = 16777216

1x16 sensor: 1GB available. shared by 16ch. each point is 2 bytes (signed)
so nPt = 1*(2^30)/(16*2) = 33554432

the 1x1 sensor uses the entire 1GB storage for the one active channel, 
but still has SCOPE_NCH set to 16. so use the nPt of the combined 16 channels for a single
channel, but still set SCOPE_NCH to get the waveform attribute format correct so we can read the
data out.

so nPt = 16 * 33554432 = 536870912
*/

#define SCOPE_NCH 16
#define SCOPE_MEM_LENGTH_MAX 33554432 // GiB memory, (16-bit x 1-ch) per point
// #define SCOPE_MEM_LENGTH_MAX 16777216 // 256MiB memory, (16-bit X 8-ch) per point
#define ANALYSIS_WAVEFORM_BASE_TYPE double
#define SCOPE_DATA_TYPE_INT int16_t
#define SCOPE_DATA_TYPE_FLOAT float
#define SCOPE_DATA_HDF5_TYPE_INT H5T_NATIVE_INT16
#define SCOPE_DATA_HDF5_TYPE_FLOAT H5T_NATIVE_FLOAT


struct waveform_attribute
{
    unsigned int chMask;
    size_t nPt;     /* number of points in each event */
    size_t nFrames; /* number of Fast Frames in each event, 0 means off */
    double dt;
    double t0;
    double ymult[SCOPE_NCH];
    double yoff[SCOPE_NCH];
    double yzero[SCOPE_NCH];
};

/* utilities */
#define bitsof(x) (8*sizeof(x))

#ifdef DEBUG
  #define debug_printf(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); \
                               } while (0)
#else
  #define debug_printf(...) ((void)0)
#endif
#define error_printf(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); \
                             } while(0)

#define MIN(x,y) (((x)>(y))?(y):(x))
#define MAX(x,y) (((x)<(y))?(y):(x))

#ifndef strlcpy
#define strlcpy(a, b, c) do {                   \
        strncpy(a, b, (c)-1);                   \
        (a)[(c)-1] = '\0';                      \
    } while (0)
#endif

char *conv16network_endian(uint16_t *buf, size_t n);
char *conv32network_endian(uint32_t *buf, size_t n);

#endif /* __COMMON_H__ */
