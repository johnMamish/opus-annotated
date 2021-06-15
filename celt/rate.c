/* Copyright (c) 2007-2008 CSIRO
   Copyright (c) 2007-2009 Xiph.Org Foundation
   Written by Jean-Marc Valin */
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "modes.h"
#include "cwrs.h"
#include "arch.h"
#include "os_support.h"

#include "entcode.h"
#include "rate.h"

static const unsigned char LOG2_FRAC_TABLE[24]={
   0,
   8,13,
  16,19,21,23,
  24,26,27,28,29,30,31,32,
  32,33,34,34,35,36,36,37,37
};

#ifdef CUSTOM_MODES

/*Determines if V(N,K) fits in a 32-bit unsigned integer.
  N and K are themselves limited to 15 bits.*/
static int fits_in32(int _n, int _k)
{
   static const opus_int16 maxN[15] = {
      32767, 32767, 32767, 1476, 283, 109,  60,  40,
       29,  24,  20,  18,  16,  14,  13};
   static const opus_int16 maxK[15] = {
      32767, 32767, 32767, 32767, 1172, 238,  95,  53,
       36,  27,  22,  18,  16,  15,  13};
   if (_n>=14)
   {
      if (_k>=14)
         return 0;
      else
         return _n <= maxN[_k];
   } else {
      return _k <= maxK[_n];
   }
}

void compute_pulse_cache(CELTMode *m, int LM)
{
   int C;
   int i;
   int j;
   int curr=0;
   int nbEntries=0;
   int entryN[100], entryK[100], entryI[100];
   const opus_int16 *eBands = m->eBands;
   PulseCache *cache = &m->cache;
   opus_int16 *cindex;
   unsigned char *bits;
   unsigned char *cap;

   cindex = (opus_int16 *)opus_alloc(sizeof(cache->index[0])*m->nbEBands*(LM+2));
   cache->index = cindex;

   /* Scan for all unique band sizes */
   for (i=0;i<=LM+1;i++)
   {
      for (j=0;j<m->nbEBands;j++)
      {
         int k;
         int N = (eBands[j+1]-eBands[j])<<i>>1;
         cindex[i*m->nbEBands+j] = -1;
         /* Find other bands that have the same size */
         for (k=0;k<=i;k++)
         {
            int n;
            for (n=0;n<m->nbEBands && (k!=i || n<j);n++)
            {
               if (N == (eBands[n+1]-eBands[n])<<k>>1)
               {
                  cindex[i*m->nbEBands+j] = cindex[k*m->nbEBands+n];
                  break;
               }
            }
         }
         if (cache->index[i*m->nbEBands+j] == -1 && N!=0)
         {
            int K;
            entryN[nbEntries] = N;
            K = 0;
            while (fits_in32(N,get_pulses(K+1)) && K<MAX_PSEUDO)
               K++;
            entryK[nbEntries] = K;
            cindex[i*m->nbEBands+j] = curr;
            entryI[nbEntries] = curr;

            curr += K+1;
            nbEntries++;
         }
      }
   }
   bits = (unsigned char *)opus_alloc(sizeof(unsigned char)*curr);
   cache->bits = bits;
   cache->size = curr;
   /* Compute the cache for all unique sizes */
   for (i=0;i<nbEntries;i++)
   {
      unsigned char *ptr = bits+entryI[i];
      opus_int16 tmp[CELT_MAX_PULSES+1];
      get_required_bits(tmp, entryN[i], get_pulses(entryK[i]), BITRES);
      for (j=1;j<=entryK[i];j++)
         ptr[j] = tmp[get_pulses(j)]-1;
      ptr[0] = entryK[i];
   }

   /* Compute the maximum rate for each band at which we'll reliably use as
       many bits as we ask for. */
   cache->caps = cap = (unsigned char *)opus_alloc(sizeof(cache->caps[0])*(LM+1)*2*m->nbEBands);
   for (i=0;i<=LM;i++)
   {
      for (C=1;C<=2;C++)
      {
         for (j=0;j<m->nbEBands;j++)
         {
            int N0;
            int max_bits;
            N0 = m->eBands[j+1]-m->eBands[j];
            /* N=1 bands only have a sign bit and fine bits. */
            if (N0<<i == 1)
               max_bits = C*(1+MAX_FINE_BITS)<<BITRES;
            else
            {
               const unsigned char *pcache;
               opus_int32           num;
               opus_int32           den;
               int                  LM0;
               int                  N;
               int                  offset;
               int                  ndof;
               int                  qb;
               int                  k;
               LM0 = 0;
               /* Even-sized bands bigger than N=2 can be split one more time.
                  As of commit 44203907 all bands >1 are even, including custom modes.*/
               if (N0 > 2)
               {
                  N0>>=1;
                  LM0--;
               }
               /* N0=1 bands can't be split down to N<2. */
               else if (N0 <= 1)
               {
                  LM0=IMIN(i,1);
                  N0<<=LM0;
               }
               /* Compute the cost for the lowest-level PVQ of a fully split
                   band. */
               pcache = bits + cindex[(LM0+1)*m->nbEBands+j];
               max_bits = pcache[pcache[0]]+1;
               /* Add in the cost of coding regular splits. */
               N = N0;
               for(k=0;k<i-LM0;k++){
                  max_bits <<= 1;
                  /* Offset the number of qtheta bits by log2(N)/2
                      + QTHETA_OFFSET compared to their "fair share" of
                      total/N */
                  offset = ((m->logN[j]+((LM0+k)<<BITRES))>>1)-QTHETA_OFFSET;
                  /* The number of qtheta bits we'll allocate if the remainder
                      is to be max_bits.
                     The average measured cost for theta is 0.89701 times qb,
                      approximated here as 459/512. */
                  num=459*(opus_int32)((2*N-1)*offset+max_bits);
                  den=((opus_int32)(2*N-1)<<9)-459;
                  qb = IMIN((num+(den>>1))/den, 57);
                  celt_assert(qb >= 0);
                  max_bits += qb;
                  N <<= 1;
               }
               /* Add in the cost of a stereo split, if necessary. */
               if (C==2)
               {
                  max_bits <<= 1;
                  offset = ((m->logN[j]+(i<<BITRES))>>1)-(N==2?QTHETA_OFFSET_TWOPHASE:QTHETA_OFFSET);
                  ndof = 2*N-1-(N==2);
                  /* The average measured cost for theta with the step PDF is
                      0.95164 times qb, approximated here as 487/512. */
                  num = (N==2?512:487)*(opus_int32)(max_bits+ndof*offset);
                  den = ((opus_int32)ndof<<9)-(N==2?512:487);
                  qb = IMIN((num+(den>>1))/den, (N==2?64:61));
                  celt_assert(qb >= 0);
                  max_bits += qb;
               }
               /* Add the fine bits we'll use. */
               /* Compensate for the extra DoF in stereo */
               ndof = C*N + ((C==2 && N>2) ? 1 : 0);
               /* Offset the number of fine bits by log2(N)/2 + FINE_OFFSET
                   compared to their "fair share" of total/N */
               offset = ((m->logN[j] + (i<<BITRES))>>1)-FINE_OFFSET;
               /* N=2 is the only point that doesn't match the curve */
               if (N==2)
                  offset += 1<<BITRES>>2;
               /* The number of fine bits we'll allocate if the remainder is
                   to be max_bits. */
               num = max_bits+ndof*offset;
               den = (ndof-1)<<BITRES;
               qb = IMIN((num+(den>>1))/den, MAX_FINE_BITS);
               celt_assert(qb >= 0);
               max_bits += C*qb<<BITRES;
            }
            max_bits = (4*max_bits/(C*((m->eBands[j+1]-m->eBands[j])<<i)))-64;
            celt_assert(max_bits >= 0);
            celt_assert(max_bits < 256);
            *cap++ = (unsigned char)max_bits;
         }
      }
   }
}

#endif /* CUSTOM_MODES */

#define ALLOC_STEPS 6

/**
 * Determines the number of relative bits to remove from each band according to the
 * "Allocation Trim". The significance of the "allocation trim" parameter is described in section
 * 4.3.3 of RFC6716: values below 5 allocate more bits for lower frequency bands, values above
 * 5 allocate more bits for higher frequency bands, and a value of 5 gives flat allocation.
 *
 * Also determines the minimum number of bits that should be used for encoding shape (PVQ). For some
 * bands, if a very small number of bits (1/8 per MDCT bin) are provided for shape encoding, it's
 * not worth it to encode shape.
 *
 * @param[in]     m          OpusCustomMode containing information about the encoding process
 * @param[in]     start      starting band (typically 0 for CELT-only mode)
 * @param[in]     end        ending band (typically 21 for full-bandwidth CELT)
 * @param[in]     alloc_trim Trim parameter on [0, 10]; values below 5 bias allocation toward
 *                           lower frequency bands, values above 5 bias allocation toward higher
 *                           freq bands, 5 is flat allocation.
 * @param[in]     C          Number of channels. Nominally 1.
 * @param[in]     LM         one of {0, 1, 2, 3} depending on which of the 4 valid CELT frame sizes
 *                           is in use. equal to log2((size of this frame) / (size of min size frame))
 * @param[out]    thresh     Returns the minimum number of bits that's worthwhile to use for shape
 *                           encoding.
 * @param[out]    trim_offset  Returns the allocation trim for each band in 1/8th (???) bits.
 */
static void calculate_trim_offset(const CELTMode* m, int start, int end, int alloc_trim, int C, int LM, int* thresh, int* trim_offset)
{
   int j;
   for (j = start; j < end;j++)
   {
      // N - number of MDCT buckets in band j
      const int N = m->eBands[j+1]-m->eBands[j];

      /* Below this threshold, we're sure not to allocate any PVQ bits */
      // Really, thresh should probably be calculated in a seperate function; there's no reason for
      // thresh and trim_offset to be calculated at the same time.
      //
      // "For each coded band, set thresh[band] to 24 times the number of MDCT bins in the band and
      //  divide by 16.  If 8 times the number of channels is greater, use that instead.  This sets
      //  the minimum allocation to one bit per channel or 48 128th bits per MDCT bin, whichever
      //  is greater."
      thresh[j] = IMAX((C)<<BITRES, (3 * N << LM << BITRES) >> 4);

      /* Tilt of the allocation curve */
      // determine the actual allocation value for band j.
      trim_offset[j] = C * N * (alloc_trim-5-LM) * (end-j-1) * (1<<(LM+BITRES)) >> 6;

      /* Giving less resolution to single-coefficient bands because they get
         more benefit from having one coarse value per coefficient*/
      // note: the number of samples in the j-th band is given by (N << LM).
      // This is described in the very last paragraph of section 4.3.3 of RFC6716; it can be found
      // by searching the document for "trim_offset"
      if ((N << LM) == 1)
         trim_offset[j] -= C<<BITRES;
   }
}

/**
 * Do a binary search over the CELT Static allocation table to determine the largest integer value
 * of 'q' (the quality factor) which doesn't exceed the remaining frame size.
 *
 * @param[in]     m          CELTMode is #defined to OpusCustomMode; defined in celt/modes.h
 * @param[in]     start      starting band (typically 0 for CELT-only mode)
 * @param[in]     end        ending band (typically 21 for full-bandwidth CELT)
 * @param[in]     offsets    Contains number of 1/8 bits to boost bands on a per-band basis
 * @param[in]     cap        Maximum possible bit allocation per band.
 * @param[in]     total      the number of 1/8th bits remaining in the frame whose allocation needs
 *                           to be spread over the bands.
 * @param[in]     C          Number of channels. Nominally 1.
 * @param[in]     LM         one of {0, 1, 2, 3} depending on which of the 4 valid CELT frame sizes
 *                           is in use. equal to log2((size of this frame) / (size of min size frame))
 * @param[in]     thresh     This holds the minimum allocation for which a PVQ vector will be
 *                           stored. "This minimum is higher than the technical limit of the PVQ
 *                           process, but very low rate allocations produce an excessively sparse
 *                           spectrum and these bands are better served by having no allocation at
 *                           all."
 * @param[in]     trim_offset  Amount of bits to add or remove on a per-band basis due to allocation
 *                             trim.
 *
 * @return Returns the largest integer value of 'q' which doesn't exceed the remaining frame budget
 */
int search_q_lo(const CELTMode* m, int start, int end, const int* offsets, const int* cap, opus_int32 total, int C, int LM, const int* thresh, const int* trim_offset)
{
   int lo, hi, j;

   lo = 1;
   hi = m->nbAllocVectors - 1;
   do
   {
      //
      int done = 0;

      int psum = 0;
      int mid = (lo+hi) >> 1;
      for (j = end; j-- > start;)
      {
         int bitsj;

         // N - number of MDCT buckets in band j
         int N = m->eBands[j+1]-m->eBands[j];

         // bitsj - number of bits in band j for quality q = mid.
         bitsj = ((C*N*m->allocVectors[mid*m->nbEBands+j]) << LM) >> 2;

         if (bitsj > 0) {
            // this IMAX is required because trim_offset is typically less than 0.
            bitsj = IMAX(0, bitsj + trim_offset[j]);
         }

         // incorporate band boost offset into bitsj
         bitsj += offsets[j];

         if ((bitsj >= thresh[j]) || done)
         {
            // If this band or any of the higher bands exceeded their minimum threshold for
            // encoding shape,
            done = 1;

            /* Don't allocate more than we can actually use */
            psum += IMIN(bitsj, cap[j]);
         } else {
            if (bitsj >= C<<BITRES)
               psum += C<<BITRES;
         }
      }

      if (psum > total) {
         // 'mid' overshot the budget
         hi = mid - 1;
      } else {
         // 'mid' undershot the budget
         lo = mid + 1;
      }
      /*printf ("lo = %d, hi = %d\n", lo, hi);*/
   } while (lo <= hi);

   return (lo - 1);
}

/**
 * Finds the exact number of bits allocated to each band via interpolation on the quality factors.
 *
 * Dear god, look at all these args.
 *
 * @param[in]     m         CELTMode is #defined to OpusCustomMode; OpusCustomMode is defined in
 *                          celt/modes.h and contains info about sample frequency, number of bands,
 *                          bit allocation vectors, and so on.
 * @param[in]     start     band number to start at, nominally 0
 * @param[in]     end       band number to end at, nominally 21
 * @param[in]     skip_start  Lowest band at which skipping can start (noninclusive). e.g. if
 *                            skip_start is 7, band 8 is the lowest one which can be skipped. This
 *                            number helps us determine if skipping bands is worth it or not.
 *                            (a band is 'skipped' if 0 bits are allocated to its shape).
 * @param[in]     bits1     (1/8th ?) Bits allocated to each band at highest integer quality q
 *                          according to table 57 and also incorporating band boosts and tilt.
 * @param[in]     bits2     (1/8th ?) Bits by which band1 will be exceeded if q is increased by 1.
 * @param[in]     thresh    Minimum number of bits for which PVQ shape will be encoded.
 * @param[in]     cap       Per-band maximum allocation vector.
 * @param[in]     total     Number of bits remaining in the frame.
 * @param[out]    _balance  Returns the number of bits remaining after initial allocation of fine
 *                          energy and shape.
 * @param[in]     skip_rsv  Number of bits to give back if no bands are skipped
 * @param[out]    intensity Returns info related to "intensity stereo" encoding, which I'll ignore
 *                          for now because I'm trying to focus on mono first.
 * @param[in]     intensity_rsv  Number of bits that are needed to signal "intensity stereo"
 *                               encoding.
 * @param[in,out] dual_stereo  In decoding mode, this returns whether the frame is in "dual stereo"
 *                             mode; in encoding mode, *dual_stereo signals whether the frame should
 *                             be encoded in "dual stereo" mode or not.
 * @param[in]     dual_stereo_rsv  Number of bits needed to signal "dual stereo mode" encoding.
 * @param[out]    bits      Bits to allocate to shape for each band. Called "pulses" in callers
 * @param[out]    ebits     Bits to allocate to energy for each band.
 * @param[out]    fine_priority  Tells which bands should be given first priority when distributing
 *                               extra bits to fine energy quantization, according to rfc6716 4.3.2.2
 * @param[in]     C         Number of channels
 * @param[in]     LM        LM is 0, 1, 2, or 3 for frames of size 120, 240, 480, and 960 respectively.
 * @param[in,out] ec        The entropy coder.
 * @param[in]     prev      Holds the previous number of bands which were not skipped (i.e. whose
 *                          shape was coded) so that some hystersis can be added to band skipping.
 * @param[in]     signalBandwidth  Allows skipping of some higher-frequency computation. This
 *                                 parameter appears to be unused and is left to set 0 in all the
 *                                 places where interp_bits2pulses() gets called, effectively
 *                                 disabling it.
 */
static OPUS_INLINE int interp_bits2pulses(const CELTMode *m, int start, int end, int skip_start,
      const int *bits1, const int *bits2, const int *thresh, const int *cap, opus_int32 total, opus_int32 *_balance,
      int skip_rsv, int *intensity, int intensity_rsv, int *dual_stereo, int dual_stereo_rsv, int *bits,
      int *ebits, int *fine_priority, int C, int LM, ec_ctx *ec, int encode, int prev, int signalBandwidth)
{
   opus_int32 psum;
   int lo, hi;
   int i, j;
   int logM;
   int stereo;
   int codedBands=-1;
   int alloc_floor;
   opus_int32 left, percoeff;
   int done;
   opus_int32 balance;
   SAVE_STACK;

   alloc_floor = C<<BITRES;
   stereo = C>1;

   logM = LM<<BITRES;

   ////////////////////////////////////////////////
   // FIND INTERPOLATION POINT
   // Binary search in increments of 1/64th for the value of q between bits1 and bits1+bits2 which
   // is as large as possible without exceeding the budget of 'total'.
   //
   // ALLOC_STEPS is determined by the precision with which interpolation between quality factors
   // in table 57 happen. RFC6716 says that it should happen in steps of 64, so the value of
   // ALLOC_STEPS is 6 and 'hi' starts at 64.
   //
   // at the end of this loop, lo contains the largest value of q (in 1/64ths) whose bit allocation
   // according to table 57 (plus allocation tilt and band boost) doesn't exceed the budget.
   lo = 0;
   hi = 1<<ALLOC_STEPS;
   for (i=0;i<ALLOC_STEPS;i++)
   {
      int mid = (lo+hi)>>1;
      psum = 0;
      done = 0;
      for (j=end;j-->start;)
      {
         // Because we're doing a binary search over 64 integers, all of our divisions can just
         // be left-shifts. Only multiply here is by mid.
         int tmp = bits1[j] + (mid*(opus_int32)bits2[j]>>ALLOC_STEPS);
         if (tmp >= thresh[j] || done)
         {
            done = 1;
            /* Don't allocate more than we can actually use */
            psum += IMIN(tmp, cap[j]);
         } else {
            if (tmp >= alloc_floor)
               psum += alloc_floor;
         }
      }
      if (psum > total)
         hi = mid;
      else
         lo = mid;
   }

   ////////////////////////////////////////////////
   // INTERPOLATE
   //
   // at the end of this, 'bits' contains the number of bits for each band.
   psum = 0;
   /*printf ("interp bisection gave %d\n", lo);*/
   done = 0;
   for (j=end;j-->start;)
   {
      // INTERPOLATION STEP
      // Note that we don't interpolate directly between values in Table 57, we first incorporate
      // band boost and tilt.
      int tmp = bits1[j] + ((opus_int32)lo*bits2[j]>>ALLOC_STEPS);

      if (tmp < thresh[j] && !done)
      {
         if (tmp >= alloc_floor)
            tmp = alloc_floor;
         else
            tmp = 0;
      } else
         done = 1;

      /* Don't allocate more than we can actually use */
      if (tmp >= cap[j]) {
          printf("Debug: tmp %i exceeds capacity %i in band %i\n", tmp, cap[j], j);
      }
      tmp = IMIN(tmp, cap[j]);
      bits[j] = tmp;
      psum += tmp;
   }

   ////////////////////////////////////////////////
   // BAND SKIPPING
   // Decide which bands to skip, working backwards from the end.
   for (codedBands=end;;codedBands--)
   {
      int band_width;
      int band_bits;
      int rem;
      j = codedBands-1;
      /* Never skip the first band, nor a band that has been boosted by
          dynalloc.
         In the first case, we'd be coding a bit to signal we're going to waste
          all the other bits.
         In the second case, we'd be coding a bit to redistribute all the bits
          we just signaled should be cocentrated in this band. */
      if (j<=skip_start)
      {
         /* Give the bit we reserved to end skipping back. */
         total += skip_rsv;
         break;
      }

      /*Figure out how many left-over bits we would be adding to this band.
        This can include bits we've stolen back from higher, skipped bands.*/
      left = total-psum;
      percoeff = celt_udiv(left, m->eBands[codedBands]-m->eBands[start]);
      left -= (m->eBands[codedBands]-m->eBands[start])*percoeff;
      rem = IMAX(left-(m->eBands[j]-m->eBands[start]),0);
      band_width = m->eBands[codedBands]-m->eBands[j];
      band_bits = (int)(bits[j] + percoeff*band_width + rem);

      /*Only code a skip decision if we're above the threshold for this band.
        Otherwise it is force-skipped.
        This ensures that we have enough bits to code the skip flag.*/
      if (band_bits >= IMAX(thresh[j], alloc_floor+(1<<BITRES)))
      {
         if (encode)
         {
            /*This if() block is the only part of the allocation function that
               is not a mandatory part of the bitstream: any bands we choose to
               skip here must be explicitly signaled.*/
            int depth_threshold;
            /*We choose a threshold with some hysteresis to keep bands from
               fluctuating in and out, but we try not to fold below a certain point. */
            if (codedBands > 17)
               depth_threshold = j<prev ? 7 : 9;
            else
               depth_threshold = 0;
#ifdef FUZZING
            if ((rand()&0x1) == 0)
#else
            if (codedBands<=start+2 || (band_bits > (depth_threshold*band_width<<LM<<BITRES)>>4 && j<=signalBandwidth))
#endif
            {
               ec_enc_bit_logp(ec, 1, 1);
               break;
            }
            ec_enc_bit_logp(ec, 0, 1);
         } else if (ec_dec_bit_logp(ec, 1)) {
            break;
         }
         /*We used a bit to skip this band.*/
         psum += 1<<BITRES;
         band_bits -= 1<<BITRES;
      }
      /*Reclaim the bits originally allocated to this band.*/
      psum -= bits[j]+intensity_rsv;
      if (intensity_rsv > 0)
         intensity_rsv = LOG2_FRAC_TABLE[j-start];
      psum += intensity_rsv;
      if (band_bits >= alloc_floor)
      {
         /*If we have enough for a fine energy bit per channel, use it.*/
         psum += alloc_floor;
         bits[j] = alloc_floor;
      } else {
         /*Otherwise this band gets nothing at all.*/
         bits[j] = 0;
      }
   }

   ////////////////////////////////////////////////
   //
   celt_assert(codedBands > start);
   /* Code the intensity and dual stereo parameters. */
   if (intensity_rsv > 0)
   {
      if (encode)
      {
         *intensity = IMIN(*intensity, codedBands);
         ec_enc_uint(ec, *intensity-start, codedBands+1-start);
      }
      else
         *intensity = start+ec_dec_uint(ec, codedBands+1-start);
   }
   else
      *intensity = 0;
   if (*intensity <= start)
   {
      total += dual_stereo_rsv;
      dual_stereo_rsv = 0;
   }
   if (dual_stereo_rsv > 0)
   {
      if (encode)
         ec_enc_bit_logp(ec, *dual_stereo, 1);
      else
         *dual_stereo = ec_dec_bit_logp(ec, 1);
   }
   else
      *dual_stereo = 0;

   ////////////////////////////////////////////////
   // ALLOCATE LEFTOVER BITS
   // Distribute leftover 1/8th bits evenly over all bands.
   left = total-psum;
   printf("leftover 1/8th bits before redistribution: %i\n", left);
   percoeff = celt_udiv(left, m->eBands[codedBands] - m->eBands[start]);
   left -= (m->eBands[codedBands] - m->eBands[start])*percoeff;

   for (j=start;j<codedBands;j++) {
      bits[j] += ((int)percoeff*(m->eBands[j+1]-m->eBands[j]));
   }

   // If there are still leftover 1/8th bits, we give them out starting at the lowest bands.
   // Every band gets some. The number of bits each band gets in excess is equal to the number
   // of MDCT samples it contains in 1/8th bits. e.g. for LM = 3 (frame size of 960 samples), the
   // lowest bucket (containing 8 samples) will get 8 1/8th bits.
   //
   // After this step, all leftover bits should be used.
   for (j=start;j<codedBands;j++)
   {
      int tmp = (int)IMIN(left, m->eBands[j+1] - m->eBands[j]);
      bits[j] += tmp;
      left -= tmp;
   }
   printf("1/8th bit allocations after step xxx?: { "); for (j=0;j<end;j++)printf("%d ", bits[j]);printf(" }\n");
   printf("left: %i\n", left);

   ////////////////////////////////////////////////
   // Split allocated bits between fine energy quantization and PVQ (???)
   //
   // Output of this loop is "ebits" and "bits".
   //
   balance = 0;
   for (j=start;j<codedBands;j++)
   {
      int N0, N, den;
      int offset;
      int NClogN;
      opus_int32 excess, bit;

      celt_assert(bits[j] >= 0);
      N0 = m->eBands[j+1]-m->eBands[j];
      N=N0<<LM;
      bit = (opus_int32)bits[j]+balance;

      if (N>1)
      {
         // what the hell is this? It's a real convoluted way to say
         // bits[j] = ((bits[j] + balance) > cap[j]) ? cap[j] : (bits[j] + balance);
         excess = MAX32(bit-cap[j],0);
         bits[j] = bit-excess;

         /* Compensate for the extra DoF in stereo */
         den = C*N;
         if ((C==2 && N>2 && !*dual_stereo && j<*intensity) ? 1 : 0) {
             den++;
         }

         NClogN = den*(m->logN[j] + logM);

         /* Offset for the number of fine bits by log2(N)/2 + FINE_OFFSET
            compared to their "fair share" of total/N */
         offset = (NClogN>>1)-den*FINE_OFFSET;

         /* N=2 is the only point that doesn't match the curve */
         if (N==2)
            offset += den<<BITRES>>2;

         /* Changing the offset for allocating the second and third
             fine energy bit */
         if (bits[j] + offset < den*2<<BITRES)
            offset += NClogN>>2;
         else if (bits[j] + offset < den*3<<BITRES)
            offset += NClogN>>3;

         /* Divide with rounding */
         ebits[j] = IMAX(0, (bits[j] + offset + (den<<(BITRES-1))));
         ebits[j] = celt_udiv(ebits[j], den)>>BITRES;

         /* Make sure not to bust */
         if (C*ebits[j] > (bits[j]>>BITRES))
            ebits[j] = bits[j] >> stereo >> BITRES;

         /* More than that is useless because that's about as far as PVQ can go */
         ebits[j] = IMIN(ebits[j], MAX_FINE_BITS);

         /* If we rounded down or capped this band, make it a candidate for the
             final fine energy pass */
         fine_priority[j] = ebits[j]*(den<<BITRES) >= bits[j]+offset;

         /* Remove the allocated fine bits; the rest are assigned to PVQ */
         bits[j] -= C*ebits[j]<<BITRES;

      } else {
         /* For N=1, all bits go to fine energy except for a single sign bit */
         // ??? That's not what it looks like this code does. I think that it should be
         //   ebits[j] = bit - excess;
         //   bits[j] = 0
         // bit idk!!! Maybe it doesn't matter because later steps will know that a band with
         // N = 1 doesn't need any shape bits.
         excess = MAX32(0,bit-(C<<BITRES));
         bits[j] = bit-excess;
         ebits[j] = 0;
         fine_priority[j] = 1;
      }

      /* Fine energy can't take advantage of the re-balancing in
          quant_all_bands().
         Instead, do the re-balancing here.*/
      if(excess > 0)
      {
         int extra_fine;
         int extra_bits;
         extra_fine = IMIN(excess>>(stereo+BITRES),MAX_FINE_BITS-ebits[j]);
         ebits[j] += extra_fine;
         extra_bits = extra_fine*C<<BITRES;
         fine_priority[j] = extra_bits >= excess-balance;
         excess -= extra_bits;
      }
      balance = excess;

      celt_assert(bits[j] >= 0);
      celt_assert(ebits[j] >= 0);
   }
   /* Save any remaining bits over the cap for the rebalancing in
       quant_all_bands(). */
   *_balance = balance;

   /* The skipped bands use all their bits for fine energy. */
   for (;j<end;j++)
   {
      ebits[j] = bits[j] >> stereo >> BITRES;
      celt_assert(C*ebits[j]<<BITRES == bits[j]);
      bits[j] = 0;
      fine_priority[j] = ebits[j]<1;
   }
   RESTORE_STACK;
   return codedBands;
}

/**
 * @param[in]     m          CELTMode is #defined to OpusCustomMode; defined in celt/modes.h
 * @param[in]     start      starting band (typically 0 for CELT-only mode)
 * @param[in]     end        ending band (typically 21 for full-bandwidth CELT)
 * @param[in]     offsets    Contains number of 1/8 bits to boost bands on a per-band basis
 *                           according to the "band boost" symbol.
 * @param[in]     cap        Maximum possible bit allocation per band.
 * @param[in]     alloc_trim Trim parameter on [0, 10]; values below 5 bias allocation toward
 *                           lower frequency bands, values above 5 bias allocation toward higher
 *                           freq bands, 5 is flat allocation.
 * @param[out]    intensity  Entropy-coded flag that determines whether there are "intensity
 *                           stereo" coded stereo channels or not.
 * @param[out]    dual_stereo  Entropy-coded flag that determines whether there are seperately-
 *                             coded stereo channels or not.
 * @param[in]     total      Not sure; I *THINK* this is just the number of 1/8th bits remaining
 *                           in the frame whose allocation needs to be spread over the bands.
 * @param[out]    balance    Remaining bits for redistribution after allocation is finished.
 * @param[out]    pulses     Bits to allocate to shape for each band. Calculated by
 *                           interp_bits2pulses()
 * @param[out]    ebits      Bits to allocate to energy for each band. Calculated by
 *                           interp_bits2pulses()
 * @param[out]    fine_priority  After all the other allocation steps are complete, any remaining
 *                               bits are distributed to the different bands for their fine energy
 *                               quantization. Starting from band 0, bands with a 'fine priority' of
 *                               0 get a bit first, followed by bands with a 'fine priority' of 1.
 *                               Leftover bits are unused. This is described in 4.3.2.2. The
 *                               details of how prioritization happens are implemented in
 *                               interp_bits2pulses.
 * @param[in]     C          Number of channels. Nominally 1.
 * @param[in]     LM         one of {0, 1, 2, 3} depending on which of the 4 valid CELT frame sizes
 *                           is in use. equal to log2((size of this frame) / (size of min size frame))
 * @param[in,out] ec         The entropy encoder we're trying to work with.
 * @param[in]     encode     1 for encode, 0 for decode
 * @param[in]     prev       Passed down to interp_bits2pulses, should contain the previous number of
 *                           bands which were not skipped (i.e. whose shape was coded) so that some
 *                           hystersis can be added to band skipping.
 * @param[in]     signalBandwidth  Allows skipping of some higher-frequency computation. This
 *                                 parameter appears to be unused and is left to set 0 in all the
 *                                 places where clt_compute_allocation() is called, effectively
 *                                 disabling it.
 *
 * A lot of this function's work happens in interp_bits2pulses().
 */
int clt_compute_allocation(const CELTMode *m, int start, int end, const int *offsets, const int *cap, int alloc_trim, int *intensity, int *dual_stereo,
      opus_int32 total, opus_int32 *balance, int *pulses, int *ebits, int *fine_priority, int C, int LM, ec_ctx *ec, int encode, int prev, int signalBandwidth)
{
   int lo, hi, len, j;
   int codedBands;
   int skip_start;
   int skip_rsv;
   int intensity_rsv;
   int dual_stereo_rsv;
   VARDECL(int, bits1);
   VARDECL(int, bits2);
   VARDECL(int, thresh);
   VARDECL(int, trim_offset);
   SAVE_STACK;


   total = IMAX(total, 0);
   len = m->nbEBands;
   skip_start = start;

   /* Reserve a bit to signal the end of manually skipped bands. */
   skip_rsv = total >= 1<<BITRES ? 1<<BITRES : 0;
   total -= skip_rsv;

   /* Reserve bits for the intensity and dual stereo parameters. */
   intensity_rsv = dual_stereo_rsv = 0;
   if (C==2)
   {
      intensity_rsv = LOG2_FRAC_TABLE[end-start];
      if (intensity_rsv>total)
         intensity_rsv = 0;
      else
      {
         total -= intensity_rsv;
         dual_stereo_rsv = total>=1<<BITRES ? 1<<BITRES : 0;
         total -= dual_stereo_rsv;
      }
   }
   ALLOC(bits1, len, int);
   ALLOC(bits2, len, int);
   ALLOC(thresh, len, int);
   ALLOC(trim_offset, len, int);

   calculate_trim_offset(m, start, end, alloc_trim, C, LM, thresh, trim_offset);

   lo = search_q_lo(m, start, end, offsets, cap, total, C, LM, thresh, trim_offset);

   /*printf ("interp between %d and %d\n", lo, hi);*/
   // This loop calculates bits1 and bits2 given trim, band-boost, and lo (the highest integer
   // value of 'q' which doesn't exceed the budget).
   //     bits1[j] - Number of bits allocated to band j by Table 57 of RFC6716 with quality set
   //                to 'lo'
   //     bits2[j] - Number of bits by which band j's allocation will exceed bits1[j] if quality
   //                is set to 'lo + 1'.
   for (j = start; j < end; j++)
   {
      int bits1j, bits2j;

      // Like, really, the innards of this loop should just be a function called "calculate bits
      // at quality factor" or something.

      int N = m->eBands[j+1]-m->eBands[j];

      // bits1j contains the number of 1/8th bits for the jth band that would be allocated if the
      // quality factor 'q' was set to 'lo'.
      bits1j = (C*N*m->allocVectors[lo*len+j] << LM) >> 2;

      // bits2j contains the number of 1/8th bits that would be allocated to the jth band if the
      // quality factor 'q' was set to 'hi'.
      // In very high bit-rate situations, 'hi' might be larger than the biggest allowed value for
      // 'q'. In that case, we just set bits2j equal to the capacity.
      bits2j = ((lo + 1) >= m->nbAllocVectors) ?
          cap[j] : (C*N*m->allocVectors[(lo + 1)*len+j] << LM) >> 2;

      // Incorporate the trim_offset
      if (bits1j > 0)
         bits1j = IMAX(0, bits1j + trim_offset[j]);
      if (bits2j > 0)
         bits2j = IMAX(0, bits2j + trim_offset[j]);

      // Incorporate band boost.
      if (lo > 0)
         bits1j += offsets[j];
      bits2j += offsets[j];

      // Looks like this searches up from the LF bands to find the band skip, but I'm not sure what
      // the rules are for determining band skip or the significance if it's set to 1.
      //
      // Seems like a frame with a single boost at a lower band (say band 3) would have skip_start
      // set to 3. Does that mean that all bands above band 3 are skipped??? That can't be right.
      // Seems like we will get more info once we look at interp_bits2pulses.
      //
      // Looking in interp_bits2pulses, it looks like band skipping starts at the top. Once the band
      // skipping loop in interp_bits2pulses reaches skip_start, it knows to stop skipping bands.
      if (offsets[j]>0)
         skip_start = j;

      //
      bits2j = IMAX(0,bits2j-bits1j);
      bits1[j] = bits1j;
      bits2[j] = bits2j;
   }

   ////////////////////////////////////////////////
   // INTERPOLATE
   codedBands = interp_bits2pulses(m, start, end, skip_start, bits1, bits2, thresh, cap,
         total, balance, skip_rsv, intensity, intensity_rsv, dual_stereo, dual_stereo_rsv,
         pulses, ebits, fine_priority, C, LM, ec, encode, prev, signalBandwidth);

   RESTORE_STACK;
   return codedBands;
}
