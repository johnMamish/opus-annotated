/* Copyright (c) 2007-2008 CSIRO
   Copyright (c) 2007-2009 Xiph.Org Foundation
   Copyright (c) 2008 Gregory Maxwell
   Written by Jean-Marc Valin and Gregory Maxwell */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef MODES_H
#define MODES_H

#include "opus_types.h"
#include "celt.h"
#include "arch.h"
#include "mdct.h"
#include "entenc.h"
#include "entdec.h"

#define MAX_PERIOD 1024

/**
 * PulseCache holds info needed to derive the per-band maximum allocation vector.
 */
typedef struct {
   int size;
   const opus_int16 *index;
   const unsigned char *bits;
   const unsigned char *caps;
} PulseCache;

/** Mode definition (opaque)
 @brief Mode definition
 */
struct OpusCustomMode {
   opus_int32 Fs;
   int          overlap;

   /**
    * eBands contains the start and stop MDCT buckets for each of the pseudo-Bark critical bands.
    * The values in eBands are determined by LM (which is determined by the number of samples in
    * the frame) as well as by whether it's a short or long frame.
    *
    * e.g. for a 2.5 ms / 120 sample frame (LM = 0), eBands would be:
    *     {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 34, 40, 48, 60, 78, 100}
    *
    * Note that the top 20 MDCT bins are left unused because - while the nominal Opus input sample
    * rate for fullband (FB) audio is 48ksps, Opus only codes up to 20kHz and throws away the rest
    * of the frequencies.
    */
   int          nbEBands;
   int          effEBands;
   opus_val16    preemph[4];
   const opus_int16   *eBands;   /**< Definition for each "pseudo-critical band" */

   int         maxLM;
   int         nbShortMdcts;

   // mode->shortMdctSize = frame_size/mode->nbShortMdcts;
   /**
    * this value holds the
    */
   int         shortMdctSize;

   /**
    * allocVectors is a table describing the number of bits in each of the 21 bands. Its values
    * are taken from table 57 in the RFC. allocVectors is populated in
    * celt/modes.c:compute_allocation_table() by interpolating a static const array. This the
    * "interpolation" that's (kind of opaquely) described in section 4.3.3 of RFC 6716.
    */
   int          nbAllocVectors; /**< Number of lines in the matrix below */
   const unsigned char   *allocVectors;   /**< Number of bits in each band for several rates */
   const opus_int16 *logN;

   const opus_val16 *window;
   mdct_lookup mdct;

   /**
    * PulseCache is defined in this file.
    */
   PulseCache cache;
};


#endif
