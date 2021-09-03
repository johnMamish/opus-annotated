/* Copyright (c) 2007 CSIRO
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

#include "laplace.h"
#include "mathops.h"

#include <stdio.h>

/* The minimum probability of an energy delta (out of 32768). */
#define LAPLACE_LOG_MINP (0)
#define LAPLACE_MINP (1<<LAPLACE_LOG_MINP)
/* The minimum number of guaranteed representable energy deltas (in one
    direction). */
#define LAPLACE_NMIN (16)

/**
 * This file implements opus's laplace-like distribution which is used to generate PDFs for range-
 * coding of Opus's coarse energy symbols.
 *
 * There's one coarse energy symbol for each band and each coarse energy symbol has two parameters
 * which implicitly define its PDF - "0 frequency" ("fs0") and "decay". The PDF that's derived
 * from these symbols consists of two geometric-like distributions, one on [1, inf) and decaying
 * with positive x and on [-1, -inf) and decaying with negative x.
 *
 *        "probability axis"
 *         ^
 *         x  <--- fs0
 *         |
 *         |
 *       x | x  <--- (1/2) * (1 - fs0 - (2 * guaranteed tail)) * (1 - decay) = freq1
 *   x x   |   x x
 * x       |       x x x x ....
 * --------|--------------------------->   "value axis"
 *         |
 * <-----|   |---------> geometric distribution
 *
 * p(n) = freq1 * (decay ^ (|n| - 1)) (for |n| > 1)
 *
 * These geometric distributions are scaled such that each one has a cumulative probability of
 *     (1/2) * (1 - fs0)
 * and each one decays with
 *     p = "decay".
 *
 * Once the distribution's probability falls below a certain value (LAPLACE_MINP), the remainder
 * of the distribution consists LAPLACE_MINP until the cumulative value of the distribution reaches
 * '1'. After that point, the probability of higher values is 0.
 *
 * The only module that uses this to encode symbols is the CELT coarse energy encoder. As a result,
 * we can have an upper bound on how many values we need above or below zero. This is defined by
 * LAPLACE_NMIN. If LAPLACE_NMIN is 16, The resulting laplace distribution is guaranteed to encode
 * values in [-16, 16]. Once the value goes above 16, the probability of encoding that symbol might
 * be 0. There should be at least 3 or 4 values before the curve hits p_min, so you should be able
 * to encode values on [-v, v] for v > 20. If you try to encode a value that's so large it "walks
 * off" the end of the curve, the encoder clamps your value to the maximum one allowed.
 *
 * For bands whose value is really expected to be close to 0 (which should be a typical case),
 * having a very aggressive (small) decay value will make the PDF largely clumped around 0,
 * decreasing the cost of encoding a 0. For bands whose value is unpredictable (like the first
 * band in an intra frame), the probability of a 0 is relatively small and the decay is large,
 * meaning that more bits need to be used to encode values because there will be more variance
 * in the distribution hence the range encoder will have more entropy to encode.
 *
 * Negative values
 * Because all range coder symbols are positive, we encode negative values by "folding over" the
 * negative half of the spectrum and interleaving negative values with positive ones. For instance,
 *     symbol value       coded value
 *     0                   0
 *     1                  -1
 *     2                   1
 *     3                  -2
 *     4                   2
 *    ...
 *     n                   floor((n - 1) / 2) for even n
 *                        -floor((n - 1) / 2) for odd n
 */

/**
 * ec_laplace_{encode,decode} approximates a discrete laplace distribution; it uses fs0 as the
 * probability of getting a 0 and then has 2 geometric distributions going to either side with
 * 'p' = decay.
 *
 * Note that the comments above 'e_prob_model' in quant_bands.c say "Laplace-Like" distribution.
 * To get a truly laplace-like distribution, f0 and decay would have to both be determined by a
 * single parameter.
 *
 * @param[in]     fs0        Q15 number representing the probability of a 0.
 * @param[in]     decay      'decay' is a Q14 number that should be positive and at most 11456.
 *                           This represents the 'p' parameter in a geometric distribution.
 */
static unsigned ec_laplace_get_freq1(unsigned fs0, int decay)
{
   unsigned ft;

   // ft = 1 - (minp * (2 * nmin)) - fs0
   // ft contains the cumulative probability of all values except 0 and the "tails"
   // on either side with laplace_minprob.
   ft = 32768 - LAPLACE_MINP * (2 * LAPLACE_NMIN) - fs0;

   // need to divide ft * (1 - decay) by 2 because ft spans in both the positive and negative
   // directions.
   return (ft * ((opus_int32)(16384 - decay)) / 2) >> 14;
}

/**
 * This function uses two geometric distributions to discretely approximate a laplace distribution.
 *
 * Once the value of
 *
 * @param[in,out] enc        Entropy encoder to add the value to
 * @param[in,out] value      Pointer to value that should be encoded.
 * @param[in]     fs         Probablity of 0 represented in Q15 format.
 * @param[in]     decay      Decay parameter, represented in Q14 format.
 */
void ec_laplace_encode(ec_enc *enc, int *value, unsigned fs, int decay)
{
   unsigned fl;
   int val = *value;
   fl = 0;

   // We need to determine fl and fs for the symbol we want to encode.
   // Inside this 'if' statement, we consider symbols one by one increasing from 1 until the laplace
   // distribution starts flooring to 0.
   //   * fl contains the cumulative probability of all buckets which have been ruled out so far.
   //   * fs contains the probability of the bucket currently under consideration.
   //   *
   if (val)
   {
      // If the value is 0,
      int i;

      // scoot
      fl = fs;
      fs = ec_laplace_get_freq1(fs, decay);

      /* Search the decaying part of the PDF.*/
      // scoot fl forward until until we find the value we want or we hit the end of the curve.
      for (i = 1; (fs >= LAPLACE_MINP) && (i < abs(val)); i++)
      {
         // we need to multiply fs by 2 because we have to consider both positive and negative
         // sides of the curve (right?)
         fl += 2 * (fs + LAPLACE_MINP);

         // fs[i] = fs * (decay^i)
         // decay is a Q0.14 number, fs is a Q0.15 number; to renormalize to a Q15 number, we need to
         // shift by 14.
         fs = (fs * (opus_int32)decay) >> 14;
      }

      if (fs < LAPLACE_MINP)
      {
         // val is in a bucket with probability LAPLACE_MINP
         // ndi_max is number of buckets over which the remaining probability is divided once
         // fs reaches a point where its probability is LAPLACE_MINP.
         int ndi_max;
         ndi_max = (32768 - fl + LAPLACE_MINP - 1) >> LAPLACE_LOG_MINP;

         if (val < 0) {
             ndi_max = (ndi_max + 1) >> 1;
         } else {
             ndi_max >>= 1;
         }


         // If the value we wanted to encode is off of the tail end of the distribution. Round it down
         // int di = IMIN(val - i, ndi_max - 1);
         // di contains the number of LAPLACE_MINP buckets that we need to scroll forward by
         int di;
         if ((ndi_max - 1) < (abs(val) - i)) {
             printf("Warning: we walked off of the end of the laplace distribution during encoding.\n");
             di = ndi_max - 1;
         } else {
             di = abs(val) - i;
         }

         if (val < 0) {
             fl += (2 * di) * LAPLACE_MINP;
             *value = -(i * di);
         } else {
             fl += (2 * di + 1) * LAPLACE_MINP;
             *value = i * di;
         }

         fs = IMIN(LAPLACE_MINP, 32768-fl);
      }
      else
      {
         // val is not in a bucket w probability LAPLACE_MINP
         // nudge fs up by
         fs += LAPLACE_MINP;

         // if val is positive, scoot fl up by 1x fs
         if (val > 0) {
             fl += fs;
         }
      }

      celt_assert((fl + fs) <= 32768);
      celt_assert(fs > 0);
   }

   printf("encoding value %i with (fl, fs) = (%i, %i); fl + fs = %i\n", *value, fl, fs, fl + fs);
   ec_encode_bin(enc, fl, fl+fs, 15);
}

int ec_laplace_decode(ec_dec *dec, unsigned fs, int decay)
{
   int val=0;
   unsigned fl;
   unsigned fm;
   fm = ec_decode_bin(dec, 15);
   fl = 0;
   if (fm >= fs)
   {
      val++;
      fl = fs;
      fs = ec_laplace_get_freq1(fs, decay)+LAPLACE_MINP;

      /* Search the decaying part of the PDF.*/
      while(fs > LAPLACE_MINP && fm >= fl+2*fs)
      {
         fs *= 2;
         fl += fs;
         fs = ((fs - (2 * LAPLACE_MINP)) * (opus_int32)decay) >> 15;
         fs += LAPLACE_MINP;
         val++;
      }

      /* Everything beyond that has probability LAPLACE_MINP. */
      if (fs <= LAPLACE_MINP)
      {
         int di;
         di = (fm-fl)>>(LAPLACE_LOG_MINP+1);
         val += di;
         fl += 2*di*LAPLACE_MINP;
      }
      if (fm < fl+fs)
         val = -val;
      else
         fl += fs;
   }
   celt_assert(fl<32768);
   celt_assert(fs>0);
   celt_assert(fl<=fm);
   celt_assert(fm<IMIN(fl+fs,32768));

   //
   ec_dec_update(dec, fl, IMIN(fl+fs,32768), 32768);
   return val;
}
