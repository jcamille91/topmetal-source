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

struct pulse 
{

size_t n_pulse;

/* number of pulses in the data is always different, so declare the arrays ambiguously
as pointers. Once we know the number of pulses in the data, then we will allocate the
the appropriate amount of memory to hold the data for the pulses' parameters.
*/

/* all the arrays will have the same number of points given by n_pulse */ 

/* do we want to use the absolute index of the wfmBuf, or separate the channels with
their own variables?
*/
uint *window; /* number of points to filter after the peak */
size_t *index; /* location of pulse peak in the data */
uint *M; /* time constant */
uint *k; /* rise time of ramp */
uint *l; /* delay of the peak */
         /* l-k is the flat-top duration */
};

// argtype for ctypes : [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]

void shaper_pulse(char *inFileName, char *outFileName, size_t k, size_t l, double M)
{

    /* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
     * rise time, l is the delay of peak, l-k is the flat-top duration, M
     * is the decay time constant (in number of samples) of the input
     * pulse.  Set M=-1.0 to deal with a step-like input function.
     */
    size_t nEventsInFile;
    //size_t idx1; if this isn't signed, indices might do werid things. even though it should never be negative.
    ssize_t iCh, i, j, jk, jl, jkl, idx1;
    double vj, vjk, vjl, vjkl, dkl, s = 0.0, pp = 0.0;
    

    struct pulse in_pulse;

    struct hdf5io_waveform_file *inWfmFile, *outWfmFile;
    struct waveform_attribute inWfmAttr, outWfmAttr;
    struct hdf5io_waveform_event_float inWfmEvent;
    struct hdf5io_waveform_event_float outWfmEvent;
    IN_WFM_BASE_TYPE *inWfmBuf;
    OUT_WFM_BASE_TYPE *outWfmBuf;

    /* setting the parameters here manually to see if this works */
    in_pulse.n_pulse = 3;
    win_buf = (uint*)malloc(in_pulse.n_pulse * sizeof(uint));
    in_pulse.window  = win_buf;

    inWfmBuf = (IN_WFM_BASE_TYPE*)malloc(inWfmFile->nPt * inWfmFile->nCh * sizeof(IN_WFM_BASE_TYPE));
    inWfmEvent.wavBuf = inWfmBuf;

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
    
    /* as a test that the indexing is working, set specific events and number of points 
    around those windows to some specific number. */

    /* since the filter uses memory, how many points before the peak do we begin filtering?
    what if the two peaks are close to one another ? (define a minimum separation for peaks)
    */

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

