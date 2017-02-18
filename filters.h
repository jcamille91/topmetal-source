#ifndef __FILTERS_H__
#define __FILTERS_H__
#include <fftw3.h>
#include <gsl/gsl_wavelet.h>
#include "common.h"

#define FFTW_NTHREADS_DEFAULT 4
#define FFTW_FLAGS_DEFAULT FFTW_ESTIMATE

// types and name mangling macro for filters.c
#define FFTW(name) fftw_ ## name
#define ANALYSIS_WAVEFORM_BASE_TYPE double
#define FFT_BASE_TYPE double /* if this is float, FFTW should be fftwf_. also include -lfftw3f in the makefile. */
                            /* if this is double, FFTW should be fftw_. also include -lfftw3 in the makefile. */
#define WAVELET_BASE_TYPE double


typedef struct filters_handle 
{
    /* public for i/o */
    ANALYSIS_WAVEFORM_BASE_TYPE *inWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *outWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *respWav;
    WAVELET_BASE_TYPE *waveletWav;

    /* private for internal work */
    size_t wavLen, respLen;
    int malloced;

    /* fft */
    int fftUsed;
    size_t fftLen;
    size_t fftwNThreads;
    unsigned fftwFlags;
    FFTW(plan) fftwPlan, fftwPlan1, fftwPlan2;
    FFT_BASE_TYPE *fftwWork, *fftwWork1;
    /* window function */
    FFT_BASE_TYPE *fftwWin;
    double fftwS1, fftwS2, dt;

    /* wavelet */
    gsl_wavelet *gslDWT;
    gsl_wavelet_workspace *gslDWTWork;
} filters_t;
/* If inWav == NULL, a space is malloced and freed upon closure. 
 * If inWav is given, it is not freed.
 * n is the waveform length.  The output is also of length n */
filters_t *filters_init(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n);
 /* for convolution and fft */
filters_t *filters_init_for_convolution(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n, size_t np);
int filters_close(filters_t *fHdl);
int savgol_np(ANALYSIS_WAVEFORM_BASE_TYPE *in, ANALYSIS_WAVEFORM_BASE_TYPE *out, size_t length, int m, int ld, int window);
int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld);
int filters_raisedCosine(filters_t *fHdl);
int filters_convolute(filters_t *fHdl);
int filters_hanning_window(filters_t *fHdl);
/* compute the spectrum in Fourier space, requires init_for_convolution with np = 0 */
int filters_fft_spectrum(filters_t *fHdl);
int filters_DWT(filters_t *fHdl); /* discrete wavelet transform */
int filters_median(filters_t *fHdl, size_t n); /* median filter with moving window size n */
int filters_trapezoidal(filters_t *fHdl, size_t k, size_t l, double M);

#endif // __FILTERS_H__
