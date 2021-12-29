/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <immintrin.h>

#include "config/aom_dsp_rtcd.h"
#include "aom_dsp/x86/intrapred_x86.h"
#include "aom_dsp/x86/intrapred_utils.h"
#include "aom_dsp/x86/lpf_common_sse2.h"

// Low bit depth functions
static DECLARE_ALIGNED(16, uint8_t, Mask[2][33][16]) = {
  { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0,
      0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0,
      0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0,
      0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,
      0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff } },
  {
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0,
        0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0,
        0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0,
        0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0,
        0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0 },
      { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff },
  },
};

/* clang-format on */
static AOM_FORCE_INLINE void dr_prediction_z1_HxW_internal_sse4_1(
    int H, int W, __m128i *dst, const uint8_t *above, int upsample_above,
    int dx) {
  const int frac_bits = 6 - upsample_above;
  const int max_base_x = ((W + H) - 1) << upsample_above;

  assert(dx > 0);
  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be calculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m128i a0, a1, a32, a16;
  __m128i diff, c3f;
  __m128i a_mbase_x;

  a16 = _mm_set1_epi16(16);
  a_mbase_x = _mm_set1_epi8(above[max_base_x]);
  c3f = _mm_set1_epi16(0x3f);

  int x = dx;
  for (int r = 0; r < W; r++) {
    __m128i b, res, res1, shift;
    __m128i a0_above, a1_above;

    int base = x >> frac_bits;
    int base_max_diff = (max_base_x - base) >> upsample_above;
    if (base_max_diff <= 0) {
      for (int i = r; i < W; ++i) {
        dst[i] = a_mbase_x;  // save 4 values
      }
      return;
    }
    if (base_max_diff > H) base_max_diff = H;
    a0_above = _mm_loadu_si128((__m128i *)(above + base));
    a1_above = _mm_loadu_si128((__m128i *)(above + base + 1));

    if (upsample_above) {
      a0_above = _mm_shuffle_epi8(a0_above, *(__m128i *)EvenOddMaskx[0]);
      a1_above = _mm_srli_si128(a0_above, 8);

      shift = _mm_srli_epi16(
          _mm_and_si128(_mm_slli_epi16(_mm_set1_epi16(x), upsample_above), c3f),
          1);
    } else {
      shift = _mm_srli_epi16(_mm_and_si128(_mm_set1_epi16(x), c3f), 1);
    }
    // lower half
    a0 = _mm_cvtepu8_epi16(a0_above);
    a1 = _mm_cvtepu8_epi16(a1_above);

    diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
    a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
    a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16

    b = _mm_mullo_epi16(diff, shift);
    res = _mm_add_epi16(a32, b);
    res = _mm_srli_epi16(res, 5);

    // uppar half
    a0 = _mm_cvtepu8_epi16(_mm_srli_si128(a0_above, 8));
    a1 = _mm_cvtepu8_epi16(_mm_srli_si128(a1_above, 8));

    diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
    a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
    a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16

    b = _mm_mullo_epi16(diff, shift);
    res1 = _mm_add_epi16(a32, b);
    res1 = _mm_srli_epi16(res1, 5);

    res = _mm_packus_epi16(res, res1);

    dst[r] =
        _mm_blendv_epi8(a_mbase_x, res, *(__m128i *)Mask[0][base_max_diff]);
    x += dx;
  }
}

static AOM_FORCE_INLINE void dr_prediction_z1_32xN_internal_sse4_1(
    int N, __m128i *dstvec, __m128i *dstvec_h, const uint8_t *above,
    int upsample_above, int dx) {
  // here upsample_above is 0 by design of av1_use_intra_edge_upsample
  (void)upsample_above;
  const int frac_bits = 6;
  const int max_base_x = ((32 + N) - 1);

  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be calculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m128i a0, a1, a32, a16;
  __m128i a_mbase_x, diff, c3f;

  a16 = _mm_set1_epi16(16);
  a_mbase_x = _mm_set1_epi8(above[max_base_x]);
  c3f = _mm_set1_epi16(0x3f);

  int x = dx;
  for (int r = 0; r < N; r++) {
    __m128i b, res, res1, res16[2];
    __m128i a0_above, a1_above;

    int base = x >> frac_bits;
    int base_max_diff = (max_base_x - base);
    if (base_max_diff <= 0) {
      for (int i = r; i < N; ++i) {
        dstvec[i] = a_mbase_x;  // save 32 values
        dstvec_h[i] = a_mbase_x;
      }
      return;
    }
    if (base_max_diff > 32) base_max_diff = 32;
    __m128i shift = _mm_srli_epi16(_mm_and_si128(_mm_set1_epi16(x), c3f), 1);

    for (int j = 0, jj = 0; j < 32; j += 16, jj++) {
      int mdiff = base_max_diff - j;
      if (mdiff <= 0) {
        res16[jj] = a_mbase_x;
      } else {
        a0_above = _mm_loadu_si128((__m128i *)(above + base + j));
        a1_above = _mm_loadu_si128((__m128i *)(above + base + j + 1));

        // lower half
        a0 = _mm_cvtepu8_epi16(a0_above);
        a1 = _mm_cvtepu8_epi16(a1_above);

        diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
        a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
        a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16
        b = _mm_mullo_epi16(diff, shift);

        res = _mm_add_epi16(a32, b);
        res = _mm_srli_epi16(res, 5);

        // uppar half
        a0 = _mm_cvtepu8_epi16(_mm_srli_si128(a0_above, 8));
        a1 = _mm_cvtepu8_epi16(_mm_srli_si128(a1_above, 8));

        diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
        a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
        a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16

        b = _mm_mullo_epi16(diff, shift);
        res1 = _mm_add_epi16(a32, b);
        res1 = _mm_srli_epi16(res1, 5);

        res16[jj] = _mm_packus_epi16(res, res1);  // 16 8bit values
      }
    }

    dstvec[r] =
        _mm_blendv_epi8(a_mbase_x, res16[0],
                        *(__m128i *)Mask[0][base_max_diff]);  // 16 8bit values

    dstvec_h[r] =
        _mm_blendv_epi8(a_mbase_x, res16[1],
                        *(__m128i *)Mask[1][base_max_diff]);  // 16 8bit values
    x += dx;
  }
}

static void dr_prediction_z1_32xN_sse4_1(int N, uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *above,
                                         int upsample_above, int dx) {
  __m128i dstvec[64], dstvec_h[64];
  dr_prediction_z1_32xN_internal_sse4_1(N, dstvec, dstvec_h, above,
                                        upsample_above, dx);
  for (int i = 0; i < N; i++) {
    _mm_storeu_si128((__m128i *)(dst + stride * i), dstvec[i]);
    _mm_storeu_si128((__m128i *)(dst + stride * i + 16), dstvec_h[i]);
  }
}

static void dr_prediction_z1_64xN_sse4_1(int N, uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *above,
                                         int upsample_above, int dx) {
  // here upsample_above is 0 by design of av1_use_intra_edge_upsample
  (void)upsample_above;
  const int frac_bits = 6;
  const int max_base_x = ((64 + N) - 1);

  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be calculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m128i a0, a1, a32, a16;
  __m128i a_mbase_x, diff, c3f;
  __m128i max_base, base_inc, mask;

  a16 = _mm_set1_epi16(16);
  a_mbase_x = _mm_set1_epi8(above[max_base_x]);
  max_base = _mm_set1_epi8(max_base_x);
  c3f = _mm_set1_epi16(0x3f);

  int x = dx;
  for (int r = 0; r < N; r++, dst += stride) {
    __m128i b, res, res1;
    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        _mm_storeu_si128((__m128i *)dst, a_mbase_x);  // save 32 values
        _mm_storeu_si128((__m128i *)(dst + 16), a_mbase_x);
        _mm_storeu_si128((__m128i *)(dst + 32), a_mbase_x);
        _mm_storeu_si128((__m128i *)(dst + 48), a_mbase_x);
        dst += stride;
      }
      return;
    }

    __m128i shift =
        _mm_srli_epi16(_mm_and_si128(_mm_set1_epi16(x), c3f), 1);  // 8 element

    __m128i a0_above, a1_above, res_val;
    for (int j = 0; j < 64; j += 16) {
      int mdif = max_base_x - (base + j);
      if (mdif <= 0) {
        _mm_storeu_si128((__m128i *)(dst + j), a_mbase_x);
      } else {
        a0_above =
            _mm_loadu_si128((__m128i *)(above + base + j));  // load 16 element
        a1_above = _mm_loadu_si128((__m128i *)(above + base + 1 + j));

        // lower half
        a0 = _mm_cvtepu8_epi16(a0_above);
        a1 = _mm_cvtepu8_epi16(a1_above);

        diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
        a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
        a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16
        b = _mm_mullo_epi16(diff, shift);

        res = _mm_add_epi16(a32, b);
        res = _mm_srli_epi16(res, 5);

        // uppar half
        a0 = _mm_cvtepu8_epi16(_mm_srli_si128(a0_above, 8));
        a1 = _mm_cvtepu8_epi16(_mm_srli_si128(a1_above, 8));

        diff = _mm_sub_epi16(a1, a0);   // a[x+1] - a[x]
        a32 = _mm_slli_epi16(a0, 5);    // a[x] * 32
        a32 = _mm_add_epi16(a32, a16);  // a[x] * 32 + 16

        b = _mm_mullo_epi16(diff, shift);
        res1 = _mm_add_epi16(a32, b);
        res1 = _mm_srli_epi16(res1, 5);

        res = _mm_packus_epi16(res, res1);  // 16 8bit values

        base_inc =
            _mm_setr_epi8((uint8_t)(base + j), (uint8_t)(base + j + 1),
                          (uint8_t)(base + j + 2), (uint8_t)(base + j + 3),
                          (uint8_t)(base + j + 4), (uint8_t)(base + j + 5),
                          (uint8_t)(base + j + 6), (uint8_t)(base + j + 7),
                          (uint8_t)(base + j + 8), (uint8_t)(base + j + 9),
                          (uint8_t)(base + j + 10), (uint8_t)(base + j + 11),
                          (uint8_t)(base + j + 12), (uint8_t)(base + j + 13),
                          (uint8_t)(base + j + 14), (uint8_t)(base + j + 15));

        mask = _mm_cmpgt_epi8(_mm_subs_epu8(max_base, base_inc),
                              _mm_setzero_si128());
        res_val = _mm_blendv_epi8(a_mbase_x, res, mask);
        _mm_storeu_si128((__m128i *)(dst + j), res_val);
      }
    }
    x += dx;
  }
}

// z3 functions
static void dr_prediction_z3_4x4_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                        const uint8_t *left, int upsample_left,
                                        int dy) {
  __m128i dstvec[4], d[4];

  dr_prediction_z1_HxW_internal_sse4_1(4, 4, dstvec, left, upsample_left, dy);
  transpose4x8_8x4_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
                            &d[0], &d[1], &d[2], &d[3]);

  *(uint32_t *)(dst + stride * 0) = _mm_cvtsi128_si32(d[0]);
  *(uint32_t *)(dst + stride * 1) = _mm_cvtsi128_si32(d[1]);
  *(uint32_t *)(dst + stride * 2) = _mm_cvtsi128_si32(d[2]);
  *(uint32_t *)(dst + stride * 3) = _mm_cvtsi128_si32(d[3]);
  return;
}

static void dr_prediction_z3_8x8_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                        const uint8_t *left, int upsample_left,
                                        int dy) {
  __m128i dstvec[8], d[8];

  dr_prediction_z1_HxW_internal_sse4_1(8, 8, dstvec, left, upsample_left, dy);
  transpose8x8_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4],
                    &dstvec[5], &dstvec[6], &dstvec[7], &d[0], &d[1], &d[2],
                    &d[3]);

  _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
  _mm_storel_epi64((__m128i *)(dst + 1 * stride), _mm_srli_si128(d[0], 8));
  _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[1]);
  _mm_storel_epi64((__m128i *)(dst + 3 * stride), _mm_srli_si128(d[1], 8));
  _mm_storel_epi64((__m128i *)(dst + 4 * stride), d[2]);
  _mm_storel_epi64((__m128i *)(dst + 5 * stride), _mm_srli_si128(d[2], 8));
  _mm_storel_epi64((__m128i *)(dst + 6 * stride), d[3]);
  _mm_storel_epi64((__m128i *)(dst + 7 * stride), _mm_srli_si128(d[3], 8));
}

static void dr_prediction_z3_4x8_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                        const uint8_t *left, int upsample_left,
                                        int dy) {
  __m128i dstvec[4], d[8];

  dr_prediction_z1_HxW_internal_sse4_1(8, 4, dstvec, left, upsample_left, dy);
  transpose4x8_8x4_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &d[0],
                        &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7]);
  for (int i = 0; i < 8; i++) {
    *(uint32_t *)(dst + stride * i) = _mm_cvtsi128_si32(d[i]);
  }
}

static void dr_prediction_z3_8x4_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                        const uint8_t *left, int upsample_left,
                                        int dy) {
  __m128i dstvec[8], d[4];

  dr_prediction_z1_HxW_internal_sse4_1(4, 8, dstvec, left, upsample_left, dy);
  transpose8x8_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
                        &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7], &d[0],
                        &d[1], &d[2], &d[3]);
  _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
  _mm_storel_epi64((__m128i *)(dst + 1 * stride), d[1]);
  _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[2]);
  _mm_storel_epi64((__m128i *)(dst + 3 * stride), d[3]);
}

static void dr_prediction_z3_8x16_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[8], d[8];

  dr_prediction_z1_HxW_internal_sse4_1(16, 8, dstvec, left, upsample_left, dy);
  transpose8x16_16x8_sse2(dstvec, dstvec + 1, dstvec + 2, dstvec + 3,
                          dstvec + 4, dstvec + 5, dstvec + 6, dstvec + 7, d,
                          d + 1, d + 2, d + 3, d + 4, d + 5, d + 6, d + 7);
  for (int i = 0; i < 8; i++) {
    _mm_storel_epi64((__m128i *)(dst + i * stride), d[i]);
    _mm_storel_epi64((__m128i *)(dst + (i + 8) * stride),
                     _mm_srli_si128(d[i], 8));
  }
}

static void dr_prediction_z3_16x8_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[16], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(8, 16, dstvec, left, upsample_left, dy);
  transpose16x8_8x16_sse2(
      &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
      &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
      &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
      &d[3], &d[4], &d[5], &d[6], &d[7]);

  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
  }
}

static void dr_prediction_z3_4x16_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[4], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(16, 4, dstvec, left, upsample_left, dy);
  transpose4x16_sse2(dstvec, d);
  for (int i = 0; i < 16; i++) {
    *(uint32_t *)(dst + stride * i) = _mm_cvtsi128_si32(d[i]);
  }
}

static void dr_prediction_z3_16x4_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[16], d[8];

  dr_prediction_z1_HxW_internal_sse4_1(4, 16, dstvec, left, upsample_left, dy);
  for (int i = 4; i < 8; i++) {
    d[i] = _mm_setzero_si128();
  }
  transpose16x8_8x16_sse2(
      &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
      &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
      &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
      &d[3], &d[4], &d[5], &d[6], &d[7]);

  for (int i = 0; i < 4; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
  }
}

static void dr_prediction_z3_8x32_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[16], d[16], dstvec_h[16], d_h[16];

  dr_prediction_z1_32xN_internal_sse4_1(8, dstvec, dstvec_h, left,
                                        upsample_left, dy);
  for (int i = 8; i < 16; i++) {
    dstvec[i] = _mm_setzero_si128();
    dstvec_h[i] = _mm_setzero_si128();
  }
  transpose16x16_sse2(dstvec, d);
  transpose16x16_sse2(dstvec_h, d_h);

  for (int i = 0; i < 16; i++) {
    _mm_storel_epi64((__m128i *)(dst + i * stride), d[i]);
  }
  for (int i = 0; i < 16; i++) {
    _mm_storel_epi64((__m128i *)(dst + (i + 16) * stride), d_h[i]);
  }
}

static void dr_prediction_z3_32x8_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                         const uint8_t *left, int upsample_left,
                                         int dy) {
  __m128i dstvec[32], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(8, 32, dstvec, left, upsample_left, dy);

  transpose16x8_8x16_sse2(
      &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
      &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
      &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
      &d[3], &d[4], &d[5], &d[6], &d[7]);
  transpose16x8_8x16_sse2(
      &dstvec[0 + 16], &dstvec[1 + 16], &dstvec[2 + 16], &dstvec[3 + 16],
      &dstvec[4 + 16], &dstvec[5 + 16], &dstvec[6 + 16], &dstvec[7 + 16],
      &dstvec[8 + 16], &dstvec[9 + 16], &dstvec[10 + 16], &dstvec[11 + 16],
      &dstvec[12 + 16], &dstvec[13 + 16], &dstvec[14 + 16], &dstvec[15 + 16],
      &d[0 + 8], &d[1 + 8], &d[2 + 8], &d[3 + 8], &d[4 + 8], &d[5 + 8],
      &d[6 + 8], &d[7 + 8]);

  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    _mm_storeu_si128((__m128i *)(dst + i * stride + 16), d[i + 8]);
  }
}

static void dr_prediction_z3_16x16_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  __m128i dstvec[16], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(16, 16, dstvec, left, upsample_left, dy);
  transpose16x16_sse2(dstvec, d);

  for (int i = 0; i < 16; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
  }
}

static void dr_prediction_z3_32x32_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  __m128i dstvec[32], d[32], dstvec_h[32], d_h[32];

  dr_prediction_z1_32xN_internal_sse4_1(32, dstvec, dstvec_h, left,
                                        upsample_left, dy);
  transpose16x16_sse2(dstvec, d);
  transpose16x16_sse2(dstvec_h, d_h);
  transpose16x16_sse2(dstvec + 16, d + 16);
  transpose16x16_sse2(dstvec_h + 16, d_h + 16);
  for (int j = 0; j < 16; j++) {
    _mm_storeu_si128((__m128i *)(dst + j * stride), d[j]);
    _mm_storeu_si128((__m128i *)(dst + j * stride + 16), d[j + 16]);
  }
  for (int j = 0; j < 16; j++) {
    _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride), d_h[j]);
    _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride + 16), d_h[j + 16]);
  }
}

static void dr_prediction_z3_64x64_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  uint8_t dstT[64 * 64];
  dr_prediction_z1_64xN_sse4_1(64, dstT, 64, left, upsample_left, dy);
  transpose(dstT, 64, dst, stride, 64, 64);
}

static void dr_prediction_z3_16x32_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  __m128i dstvec[16], d[16], dstvec_h[16], d_h[16];

  dr_prediction_z1_32xN_internal_sse4_1(16, dstvec, dstvec_h, left,
                                        upsample_left, dy);
  transpose16x16_sse2(dstvec, d);
  transpose16x16_sse2(dstvec_h, d_h);
  // store
  for (int j = 0; j < 16; j++) {
    _mm_storeu_si128((__m128i *)(dst + j * stride), d[j]);
    _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride), d_h[j]);
  }
}

static void dr_prediction_z3_32x16_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  __m128i dstvec[32], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(16, 32, dstvec, left, upsample_left, dy);
  for (int i = 0; i < 32; i += 16) {
    transpose16x16_sse2((dstvec + i), d);
    for (int j = 0; j < 16; j++) {
      _mm_storeu_si128((__m128i *)(dst + j * stride + i), d[j]);
    }
  }
}

static void dr_prediction_z3_32x64_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  uint8_t dstT[64 * 32];
  dr_prediction_z1_64xN_sse4_1(32, dstT, 64, left, upsample_left, dy);
  transpose(dstT, 64, dst, stride, 32, 64);
}

static void dr_prediction_z3_64x32_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  uint8_t dstT[32 * 64];
  dr_prediction_z1_32xN_sse4_1(64, dstT, 32, left, upsample_left, dy);
  transpose(dstT, 32, dst, stride, 64, 32);
  return;
}

static void dr_prediction_z3_16x64_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  uint8_t dstT[64 * 16];
  dr_prediction_z1_64xN_sse4_1(16, dstT, 64, left, upsample_left, dy);
  transpose(dstT, 64, dst, stride, 16, 64);
}

static void dr_prediction_z3_64x16_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                          const uint8_t *left,
                                          int upsample_left, int dy) {
  __m128i dstvec[64], d[16];

  dr_prediction_z1_HxW_internal_sse4_1(16, 64, dstvec, left, upsample_left, dy);
  for (int i = 0; i < 64; i += 16) {
    transpose16x16_sse2(dstvec + i, d);
    for (int j = 0; j < 16; j++) {
      _mm_storeu_si128((__m128i *)(dst + j * stride + i), d[j]);
    }
  }
}

void av1_dr_prediction_z3_sse4_1(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                                 const uint8_t *above, const uint8_t *left,
                                 int upsample_left, int dx, int dy) {
  (void)above;
  (void)dx;
  assert(dx == 1);
  assert(dy > 0);

  if (bw == bh) {
    switch (bw) {
      case 4:
        dr_prediction_z3_4x4_sse4_1(dst, stride, left, upsample_left, dy);
        break;
      case 8:
        dr_prediction_z3_8x8_sse4_1(dst, stride, left, upsample_left, dy);
        break;
      case 16:
        dr_prediction_z3_16x16_sse4_1(dst, stride, left, upsample_left, dy);
        break;
      case 32:
        dr_prediction_z3_32x32_sse4_1(dst, stride, left, upsample_left, dy);
        break;
      case 64:
        dr_prediction_z3_64x64_sse4_1(dst, stride, left, upsample_left, dy);
        break;
      default: assert(0 && "Invalid block size");
    }
  } else {
    if (bw < bh) {
      if (bw + bw == bh) {
        switch (bw) {
          case 4:
            dr_prediction_z3_4x8_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 8:
            dr_prediction_z3_8x16_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 16:
            dr_prediction_z3_16x32_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 32:
            dr_prediction_z3_32x64_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          default: assert(0 && "Invalid block size");
        }
      } else {
        switch (bw) {
          case 4:
            dr_prediction_z3_4x16_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 8:
            dr_prediction_z3_8x32_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 16:
            dr_prediction_z3_16x64_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          default: assert(0 && "Invalid block size");
        }
      }
    } else {
      if (bh + bh == bw) {
        switch (bh) {
          case 4:
            dr_prediction_z3_8x4_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 8:
            dr_prediction_z3_16x8_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 16:
            dr_prediction_z3_32x16_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 32:
            dr_prediction_z3_64x32_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          default: assert(0 && "Invalid block size");
        }
      } else {
        switch (bh) {
          case 4:
            dr_prediction_z3_16x4_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 8:
            dr_prediction_z3_32x8_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          case 16:
            dr_prediction_z3_64x16_sse4_1(dst, stride, left, upsample_left, dy);
            break;
          default: assert(0 && "Invalid block size");
        }
      }
    }
  }
}
