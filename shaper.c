#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "common.h"
#include "hdf5rawWaveformIo.h"

typedef SCOPE_DATA_TYPE_FLOAT IN_WFM_BASE_TYPE;
typedef SCOPE_DATA_TYPE_FLOAT OUT_WFM_BASE_TYPE;

typedef struct peaks_handle 
{
    // structure with left and rightmost points of a pulse, 
    // and the filter parameters for each pulse. might temporarily fix l,k in actual
    // actual filtering for the moment.

    size_t nPk;
    size_t *LEFT; // beginning of a peak
    size_t *RIGHT;
    size_t *l;
    size_t *k;
    double *M;

} peaks_t;

// argtype for ctypes : [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]
peaks_t *init_peaks(peaks_t *input) {

return input;

}

void read_arr(size_t *input){

    size_t *y;

    y = input;

    fprintf(stderr, "input 0: %zu\n", y[0]);
    fprintf(stderr, "input 1: %zu\n", y[1]);

}

void read_struct(peaks_t *input){

    size_t nPk; 
    size_t *LF; 
    size_t *RT; 
    size_t *l; 
    size_t *k;
    double *M;

    nPk = input->nPk;
    LF = input->LEFT;
    RT = input->RIGHT;
    l = input->l;
    k = input->k;
    M = input->M;

    fprintf(stderr, "nPk's in struct: %zu\n", nPk);
    fprintf(stderr, "left 0: %zu\n", LF[0]);
    fprintf(stderr, "right 0: %zu\n", RT[0]);
    fprintf(stderr, "left 1: %zu\n", LF[1]);
    fprintf(stderr, "right 1: %zu\n", RT[1]);


    fprintf(stderr, "l[0]: %zu\n", l[0]);
    fprintf(stderr, "k[0]: %zu\n", k[0]);
    fprintf(stderr, "M[0]: %f\n", M[0]);

    fprintf(stderr, "l[1]: %zu\n", l[1]);
    fprintf(stderr, "k[1]: %zu\n", k[1]);
    fprintf(stderr, "M[1]: %f\n", M[1]);

}

void shaper_peaks(double * input, double * filter, size_t length, peaks_t *PEAKS)
{
    /* intended to be used by python numpy arrays, to quickly test different trapezoidal
     filtering parameters on data. */

      /* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
     * rise time, l is the delay of peak, l-k is the flat-top duration, M
     * is the decay time constant (in number of samples) of the input
     * pulse.  Set M=-1.0 to deal with a step-like input function.
     */
    size_t l, k, ipk, start, stop;
    ssize_t i, j, jk, jl, jkl, idx1 = 0;
    double M, vj, vjk, vjl, vjkl, dkl, s = 0.0, pp = 0.0;
    
    fprintf(stderr, "l[0]: %zu\n", PEAKS->l[0]);
    fprintf(stderr, "k[0]: %zu\n", PEAKS->k[0]);
    fprintf(stderr, "M[0]: %f\n", PEAKS->M[0]);

    fprintf(stderr, "input[1000]: %f\n", input[1000]);
    fprintf(stderr, "input[5]: %f\n", input[5]);
    fprintf(stderr, "input[8000]: %f\n", input[8000]);


    for(ipk = 0; ipk < PEAKS->nPk; ipk++) {
        // if first peak, set "start" to 0.

        // if last peak, set "stop" to "length". or, in peak detection,
        // set final RIGHT[] value to the end of the array,
        start = PEAKS->LEFT[ipk];
        stop = PEAKS->RIGHT[ipk];
        l = PEAKS->l[ipk];
        k = PEAKS->k[ipk];
        M = PEAKS->M[ipk];

        for(i = start; i < stop; i++) {
            j=i; jk = j-k; jl = j-l; jkl = j-k-l;
    
    /*   "RESULT = CONDITION  ?     x      :     y      " 
            is a compact if-else statment, "ternary operator" */
            vj   = (j   >= 0) ? input[j]   : input[idx1];
            vjk  = (jk  >= 0) ? input[jk]  : input[idx1];
            vjl  = (jl  >= 0) ? input[jl]  : input[idx1];
            vjkl = (jkl >= 0) ? input[jkl] : input[idx1];

            dkl = vj - vjk - vjl + vjkl;
            pp = pp + dkl;


            // commented out conditional for speed. if we do testing/use steps again, uncomment.
            // for real data, we should only have positive fall times.

            s = s + pp + dkl * M;
            // if(M >= 0.0) {
            //     s = s + pp + dkl * M;
            // }
            
            // else {  //infinite decay time, so the input is a step function 
            //     s = s + dkl;
            // }
            
            filter[i] = s / (fabs(M) * (double)k);
        }
    }
}

void trapezoid(float * input, float * filter, size_t length, size_t k, size_t l, double M)
{
    /* intended to be used by python numpy arrays, to quickly test different trapezoidal
     filtering parameters on data. */

      /* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
     * rise time, l is the delay of peak, l-k is the flat-top duration, M
     * is the decay time constant (in number of samples) of the input
     * pulse.  Set M=-1.0 to deal with a step-like input function.
     */

    ssize_t i, j, jk, jl, jkl, idx1=0;
    double vj, vjk, vjl, vjkl, dkl, s = 0.0, pp = 0.0;

    /* can use these to verify that the desired input is reaching the function.
    fprintf(stderr, "length: %zu\n", length);
    fprintf(stderr, "k: %zu\n", k);
    fprintf(stderr, "l: %zu\n", l);
    fprintf(stderr, "M: %f\n", M);
    fprintf(stderr, "in[999]: %f\n", input[998]);
    fprintf(stderr, "in[1000]: %f\n", input[999]);
    fprintf(stderr, "in[1001]: %f\n", input[1000]);
    */

    for(i = 0; i < length-1; i++) {
        j=i; jk = j-k; jl = j-l; jkl = j-k-l;
/*   "RESULT = CONDITION  ?     x      :     y      " 
        is a compact if-else statment, "ternary operator" */
        vj   = (j   >= 0) ? input[j]   : input[idx1];
        vjk  = (jk  >= 0) ? input[jk]  : input[idx1];
        vjl  = (jl  >= 0) ? input[jl]  : input[idx1];
        vjkl = (jkl >= 0) ? input[jkl] : input[idx1];

        dkl = vj - vjk - vjl + vjkl;
        pp = pp + dkl;

        if(M >= 0.0) {
            s = s + pp + dkl * M;
        }
        
        else { /* infinite decay time, so the input is a step function */
            s = s + dkl;
        }
        
        filter[i] = s / (fabs(M) * (double)k);
    }
}
void shaper(char *inFileName, char *outFileName, size_t k, size_t l, double M)
{

    /* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
     * rise time, l is the delay of peak, l-k is the flat-top duration, M
     * is the decay time constant (in number of samples) of the input
     * pulse.  Set M=-1.0 to deal with a step-like input function.
     */


    /*  in this stage, we want L/R for each peak. at a higher level
    we could choose different ways to calculate L/R based on a window and peak center.
    but maybe it won't be symmetric since the rise and fall times are not.

    then we could use a struct with multiple parallel arrays.

    i guess the question is how many samples before the rise time do we need for a 
    good shaping of the pulse to estimate the height and width. and what amplitude do we expect
    from the pulse?  does it represent the actual peak of the original exponential or something else?

    in CTypes, we can use the structure functionality, with the pointer/byref() functionality to pass
    parallel arrays for peak locations.

    POINTER(c_type) is a datatype for ctypes, we can use this inside of the Structure Class.
    */
    size_t nEventsInFile;
    //size_t idx1; if this isn't signed, indices might do weird things. even though it should never be negative.
    ssize_t iCh, i, j, jk, jl, jkl, idx1;
    double vj, vjk, vjl, vjkl, dkl, s = 0.0, pp = 0.0;
    


    struct hdf5io_waveform_file *inWfmFile, *outWfmFile;
    struct waveform_attribute inWfmAttr, outWfmAttr;
    struct hdf5io_waveform_event_float inWfmEvent;
    struct hdf5io_waveform_event_float outWfmEvent;
    IN_WFM_BASE_TYPE *inWfmBuf;
    OUT_WFM_BASE_TYPE *outWfmBuf;

    inWfmFile = hdf5io_open_file_for_read(inFileName);
    hdf5io_read_waveform_attribute_in_file_header(inWfmFile, &inWfmAttr);
    fprintf(stderr, "waveform_attribute:\n"
            "     chMask  = 0x%02x\n"
            "     nPt     = %zd\n"
            "     nFrames = %zd\n"
            "     dt      = %g\n"
            "     t0      = %g\n"
            "     ymult   = %g %g %g %g %g %g %g %g\n"
            "     yoff    = %g %g %g %g %g %g %g %g\n"
            "     yzero   = %g %g %g %g %g %g %g %g\n",
            inWfmAttr.chMask, inWfmAttr.nPt, inWfmAttr.nFrames, inWfmAttr.dt,
            inWfmAttr.t0, inWfmAttr.ymult[0], inWfmAttr.ymult[1], inWfmAttr.ymult[2],
            inWfmAttr.ymult[3], inWfmAttr.ymult[4], inWfmAttr.ymult[5], inWfmAttr.ymult[6],
            inWfmAttr.ymult[7], inWfmAttr.yoff[0], inWfmAttr.yoff[1], inWfmAttr.yoff[2],
            inWfmAttr.yoff[3], inWfmAttr.yoff[4], inWfmAttr.yoff[5], inWfmAttr.yoff[6],
            inWfmAttr.yoff[7], inWfmAttr.yzero[0], inWfmAttr.yzero[1], inWfmAttr.yzero[2],
            inWfmAttr.yzero[3], inWfmAttr.yzero[4], inWfmAttr.yzero[5], inWfmAttr.yzero[6],
            inWfmAttr.yzero[7]);

    nEventsInFile = hdf5io_get_number_of_events(inWfmFile);
    fprintf(stderr, "Number of events in file: %zd\n", nEventsInFile);
    
    /*
     v = inWfmAttr.chMask;
     for(c=0; v; c++) v &= v - 1; // Brian Kernighan's way of counting bits
     chGrpLen = inWfmFile->nCh / c;
     i=0;
     for(v=0; v<SCOPE_NCH; v++)
         if((inWfmAttr.chMask >> v) & 0x01) { chGrpIdx[i] = v; i++; }
    */
  

    fprintf(stderr, "k: %zu\n", k);
    fprintf(stderr, "l: %zu\n", l);
    fprintf(stderr, "M: %f\n", M);


    /* output */
    outWfmFile = hdf5io_open_file(outFileName, inWfmFile->nWfmPerChunk, inWfmFile->nCh);
    memcpy(&outWfmAttr, &inWfmAttr, sizeof(inWfmAttr));
    outWfmAttr.nPt = inWfmAttr.nPt;
    outWfmAttr.nFrames = 0;
    hdf5io_write_waveform_attribute_in_file_header(outWfmFile, &outWfmAttr);

    inWfmBuf = (IN_WFM_BASE_TYPE*)malloc(inWfmFile->nPt * inWfmFile->nCh * sizeof(IN_WFM_BASE_TYPE));
    inWfmEvent.wavBuf = inWfmBuf;
    outWfmBuf = (OUT_WFM_BASE_TYPE*)malloc(outWfmFile->nPt * outWfmFile->nCh * sizeof(OUT_WFM_BASE_TYPE));
    outWfmEvent.wavBuf = outWfmBuf;
    
    /* this loop works for multiple sensors, the pixels repeat from 0->5184 for each
       individual sensor. e.g. nCh is (nSensor*5184) */
    for(inWfmEvent.eventId = 0; inWfmEvent.eventId < 1; inWfmEvent.eventId++) {
        hdf5io_read_event_float(inWfmFile, &inWfmEvent);
        for(iCh=0; iCh < inWfmFile->nCh; iCh++) {
            idx1 = iCh*inWfmFile->nPt;
            s = 0.0; pp = 0.0; /* reset these each pixel, because
             the algorithm is recursive. */

            // size_t k, size_t l, double M
            // double s, pp;
            // ssize_t i, j, jk, jl, jkl;
            // double vj, vjk, vjl, vjkl, dkl;
            // s = 0.0; pp = 0.0;
            for(i = 0; i < inWfmAttr.nPt; i++) {
                j=i; jk = j-k; jl = j-l; jkl = j-k-l;
                // "condition ? x : y" is a compact if-else statment, "ternary operator"
                vj   = (j   >= 0) ? inWfmBuf[idx1 + j]   : inWfmBuf[idx1];
                vjk  = (jk  >= 0) ? inWfmBuf[idx1 + jk]  : inWfmBuf[idx1];
                vjl  = (jl  >= 0) ? inWfmBuf[idx1 + jl]  : inWfmBuf[idx1];
                vjkl = (jkl >= 0) ? inWfmBuf[idx1 + jkl] : inWfmBuf[idx1];

                dkl = vj - vjk - vjl + vjkl;
                pp = pp + dkl;

                if(M >= 0.0) {
                    s = s + pp + dkl * M;
                }
                
                else { /* infinite decay time, so the input is a step function */
                    s = s + dkl;
                }
                
                outWfmBuf[idx1 + i] = s / (fabs(M) * (double)k);
            }
        }     
            outWfmEvent.eventId = inWfmEvent.eventId;
            hdf5io_write_event(outWfmFile, &outWfmEvent);
            hdf5io_flush_file(outWfmFile);
    }
        free(inWfmBuf);
        free(outWfmBuf);
        hdf5io_close_file(inWfmFile);
        hdf5io_flush_file(outWfmFile);
        hdf5io_close_file(outWfmFile);

    //return EXIT_SUCCESS;
}


int main()
{

//shaper("../data_TM1x1/out22_dmux.h5", "../data_TM1x1/out22_filter.h5", 10, 20, 40);
//shaper("../data_TM1x1/step_dmux.h5", "../data_TM1x1/step_filter.h5", 10, 20, -1);

shaper("../data_TM1x1/exp.h5", "../data_TM1x1/exp_filter.h5", 10, 20, 40);
return EXIT_SUCCESS;

}


