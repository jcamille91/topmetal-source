#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "common.h"
#include "hdf5rawWaveformIo.h"

typedef SCOPE_DATA_TYPE_INT IN_WFM_BASE_TYPE;
typedef SCOPE_DATA_TYPE_FLOAT OUT_WFM_BASE_TYPE;

// argtype for ctypes : [c_char_p, c_char_p, c_ulong, c_double, c_ulong, c_ulong, c_ulong, c_double]

void shaper(char *inFileName, char *outFileName, size_t k, size_t l, double M)
{

    /* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
     * rise time, l is the delay of peak, l-k is the flat-top duration, M
     * is the decay time constant (in number of samples) of the input
     * pulse.  Set M=-1.0 to deal with a step-like input function.
     */
    
    double s, pp;
    size_t idx1;
    ssize_t i, j, jk, jl, jkl;
    double vj, vjk, vjl, vjkl, dkl;
    s = 0.0; pp = 0.0;  


    struct hdf5io_waveform_file *inWfmFile, *outWfmFile;
    struct waveform_attribute inWfmAttr, outWfmAttr;
    struct hdf5io_waveform_event_int inWfmEvent;
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
    fprintf(stderr, "M: %zu\n", M);


    /* output */
    outWfmFile = hdf5io_open_file(outFileName, inWfmFile->nWfmPerChunk, mNCh * inWfmFile->nCh);
    memcpy(&outWfmAttr, &inWfmAttr, sizeof(inWfmAttr));
    outWfmAttr.nPt = inWfmAttr.nPt;
    outWfmAttr.nFrames = 0;
    hdf5io_write_waveform_attribute_in_file_header(outWfmFile, &outWfmAttr);

    inWfmBuf = (IN_WFM_BASE_TYPE*)malloc(inWfmFile->nPt * inWfmFile->nCh * sizeof(IN_WFM_BASE_TYPE));
    inWfmEvent.wavBuf = inWfmBuf;
    outWfmBuf = (OUT_WFM_BASE_TYPE*)malloc(outWfmFile->nPt * outWfmFile->nCh * sizeof(OUT_WFM_BASE_TYPE));
    outWfmEvent.wavBuf = outWfmBuf;
    
   // number of points needs to be restricted so we don't go out of bounds... how does this recursive 
   // index??
    waveLen = inWfmAttr.nPt - (l - k) - M

    for(inWfmEvent.eventId = 0; inWfmEvent.eventId < 1; inWfmEvent.eventId++) {
        hdf5io_read_event_int(inWfmFile, &inWfmEvent);
        for(iCh=0; iCh < inWfmFile->nCh; iCh++) {
            idx1 = iCh * inWfmFile->nPt;
            // put iCh in the indices the same way we do it for the demux, so we can extend this to multiple
            // channels.

            // size_t wavLen, size_t k, size_t l, double M
            // double s, pp;
            // ssize_t i, j, jk, jl, jkl;
            // double vj, vjk, vjl, vjkl, dkl;
            // s = 0.0; pp = 0.0;
            for(i=0; i<wavLen; i++) {
                j=i; jk = j-k; jl = j-l; jkl = j-k-l;
                // "condition ? x : y" is a compact if-else statment, "ternary operator"
                vj   = (j   >= 0) ? inWfmBuf[idx1 + j]   : inWfmBuf[idx1];
                vjk  = (jk  >= 0) ? inWfmBuf[idx1 + jk]  : inWfmBuf[idx1];
                vjl  = (jl  >= 0) ? inWfmBuf[idx1 + jl]  : inWfmBuf[idx1];
                vjkl = (jkl >= 0) ? inWfmBuf[idx1 + jkl] : inWfmBuf[idx1];

                dkl = vj - vjk - vjl + vjkl;
                pp = pp + dkl;

                if(M>=0.0) {
                    s = s + pp + dkl * M;
                }
                
                else { /* infinite decay time, so the input is a step function */
                    s = s + dkl;
                }
                
                outWfmBuf[idx1 + i] = s / (fabs(M) * (double)k);

 
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

demux("dec_0x5f_20events.h5", "dmuxOUT_exe.h5", 6879, 4, 20736, 1, 2, 20737);

return EXIT_SUCCESS;
}


