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

void demux(char *inFileName, char *outFileName, size_t mStart, double mChLen, size_t mNCh, size_t mChOff, size_t mChSpl, double frameSize)

{
    ssize_t i, j, idx0, idx1, iCh, imCh, iFrame=0;
    size_t wfmOff, nEventsInFile;
    // size_t mStart, mNCh, mChOff, mChSpl; //chGrpLen;
    // double frameSize, mChLen;
    double avg = 0, /*sigma = 0,*/ val = 0;
    //unsigned int v, c;
    //size_t chGrpIdx[SCOPE_NCH] = {0};
    //char *inFileName, *outFileName;
    struct hdf5io_waveform_file *inWfmFile, *outWfmFile;
    struct waveform_attribute inWfmAttr, outWfmAttr;
    struct hdf5io_waveform_event_int inWfmEvent;
    struct hdf5io_waveform_event_float outWfmEvent;
    IN_WFM_BASE_TYPE *inWfmBuf;
    OUT_WFM_BASE_TYPE *outWfmBuf;

    /*
    if(argc<8) {
        error_printf("%s inFileName outFileName mStart, mChLen, mNCh, mChOff, mChSpl, [frameSize]\n",
                     argv[0]);
        error_printf("mStart     : starting sample in an original frame to demux\n");
        error_printf("mChLen     : length of each demux-ed channel (float, samples)\n");
        error_printf("mNCh       : number of demux-ed channels in each original channel\n");
        error_printf("mChOff     : in each demux-ed channel, only mChSpl samples starting\n");
        error_printf("mChSpl     : at mChOff sample from the channel boundary is recorded\n");
        error_printf("[frameSize]: specify when not acquired with frames\n");
        return EXIT_FAILURE;
    }



    inFileName = argv[1];
    outFileName = argv[2];
    
    mStart = atol(argv[3]);
    mChLen = atof(argv[4]);
    mNCh   = atol(argv[5]);
    mChOff = atol(argv[6]);
    mChSpl = atol(argv[7]);
    */


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
    
    // this is just some logic for handling the frameSize, which we don't currently need.
    // just input the frameSize even if it is redundant with respect to the other inputs.

    /*
    if(inWfmAttr.nFrames > 0) {
        frameSize = inWfmAttr.nPt / (double)inWfmAttr.nFrames;
    } else {
        frameSize = (double)inWfmAttr.nPt;
    }

     v = inWfmAttr.chMask;
     for(c=0; v; c++) v &= v - 1; // Brian Kernighan's way of counting bits
     chGrpLen = inWfmFile->nCh / c;
     i=0;
     for(v=0; v<SCOPE_NCH; v++)
         if((inWfmAttr.chMask >> v) & 0x01) { chGrpIdx[i] = v; i++; }
    
    wfmOff = 0;
    //if we define framesize, do this
    if(frameSize > 0) {
        frameSize = frameSize;
    */

    // wfmoff and mstart are both defined in the case that the frameSize and the length of
    // the actual data are different. this would allow you to skip over unwanted 'bad' data periodically.
    inWfmAttr.nFrames = (size_t)((inWfmAttr.nPt - mStart) / frameSize);
    wfmOff = mStart;
    mStart = 0;
    //}    

    fprintf(stderr, "mStart: %zu\n", mStart);
    fprintf(stderr, "mChLen: %f\n", mChLen);
    fprintf(stderr, "mNCh: %zu\n", mNCh);
    fprintf(stderr, "mChOff: %zu\n", mChOff);
    fprintf(stderr, "mChSpl: %zu\n", mChSpl);
    fprintf(stderr, "frameSize: %f\n", frameSize);

    /* output */
    outWfmFile = hdf5io_open_file(outFileName, inWfmFile->nWfmPerChunk, mNCh * inWfmFile->nCh);
    memcpy(&outWfmAttr, &inWfmAttr, sizeof(inWfmAttr));
    outWfmAttr.nPt = MAX(inWfmAttr.nFrames, 1);
    outWfmAttr.nFrames = 0;
    hdf5io_write_waveform_attribute_in_file_header(outWfmFile, &outWfmAttr);

    inWfmBuf = (IN_WFM_BASE_TYPE*)malloc(inWfmFile->nPt * inWfmFile->nCh * sizeof(IN_WFM_BASE_TYPE));
    inWfmEvent.wavBuf = inWfmBuf;
    outWfmBuf = (OUT_WFM_BASE_TYPE*)malloc(outWfmFile->nPt * outWfmFile->nCh * sizeof(OUT_WFM_BASE_TYPE));
    outWfmEvent.wavBuf = outWfmBuf;
    
    for(inWfmEvent.eventId = 0; inWfmEvent.eventId < 1; inWfmEvent.eventId++) {
        hdf5io_read_event_int(inWfmFile, &inWfmEvent);
        for(iCh=0; iCh < inWfmFile->nCh; iCh++) {
            for(iFrame = 0; iFrame < MAX(inWfmAttr.nFrames, 1); iFrame++) {
                for(i=mStart; i < MIN(frameSize, mStart + mNCh * mChLen); i++) {
                    imCh = (i - mStart) / mChLen;
                    j = i - mStart - imCh * mChLen - mChOff;
                    /* printf("i = %zd j = %zd iCh = %zd imCh = %zd iFrame = %zd\n", */
                    /*        i, j, iCh, imCh, iFrame); */
                    if(j >= 0 && j < mChSpl) {
                        /* idx0 = (mNCh * iCh + imCh) * outWfmFile->nPt + mChSpl * iFrame + j; */
                        idx1 = wfmOff + iCh * inWfmFile->nPt + iFrame * frameSize + i;
                        /* avg += inWfmBuf[idx1]; */
                        val = (inWfmBuf[idx1]
                            - inWfmAttr.yoff[1])
                            * inWfmAttr.ymult[1]
                            + inWfmAttr.yzero[1];
                        avg += val;
                        /* outWfmBuf[idx0] = inWfmBuf[idx1]; */
                        //[chGrpIdx[v]]
                        if (j == (mChSpl - 1)) {
                            avg /= (double)mChSpl;
                            idx0 = (mNCh * iCh + imCh) * outWfmFile->nPt + iFrame;
                            outWfmBuf[idx0] = (float)avg;
                            avg = 0; 

                        }
                    }
                }
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

demux("dec_0x5f_20events.h5", "dmuxOUT_exe.h5", 6879, 4, 20736, 1, 2, 20737);

return EXIT_SUCCESS;
}


