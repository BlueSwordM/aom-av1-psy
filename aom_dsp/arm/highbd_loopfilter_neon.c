/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <arm_neon.h>

#include "config/aom_dsp_rtcd.h"
#include "config/aom_config.h"

#include "aom/aom_integer.h"
#include "aom_dsp/arm/transpose_neon.h"

// TODO(b/217462944): rename functions from libgav1 to libaom style after code
// is ported.

#define kBitdepth10 10

static INLINE int16x4_t Clip3S16(const int16x4_t val, const int16x4_t low,
                                 const int16x4_t high) {
  return vmin_s16(vmax_s16(val, low), high);
}

static INLINE uint16x8_t ConvertToUnsignedPixelU16(int16x8_t val,
                                                   int bitdepth) {
  const int16x8_t low = vdupq_n_s16(0);
  const uint16x8_t high = vdupq_n_u16((1 << bitdepth) - 1);

  return vminq_u16(vreinterpretq_u16_s16(vmaxq_s16(val, low)), high);
}

// (abs(p1 - p0) > thresh) || (abs(q1 - q0) > thresh)
static INLINE uint16x4_t Hev(const uint16x8_t abd_p0p1_q0q1,
                             const uint16_t thresh) {
  const uint16x8_t a = vcgtq_u16(abd_p0p1_q0q1, vdupq_n_u16(thresh));
  return vorr_u16(vget_low_u16(a), vget_high_u16(a));
}

// abs(p0 - q0) * 2 + abs(p1 - q1) / 2 <= outer_thresh
static INLINE uint16x4_t OuterThreshold(const uint16x4_t p1,
                                        const uint16x4_t p0,
                                        const uint16x4_t q0,
                                        const uint16x4_t q1,
                                        const uint16_t outer_thresh) {
  const uint16x4_t abd_p0q0 = vabd_u16(p0, q0);
  const uint16x4_t abd_p1q1 = vabd_u16(p1, q1);
  const uint16x4_t p0q0_double = vshl_n_u16(abd_p0q0, 1);
  const uint16x4_t p1q1_half = vshr_n_u16(abd_p1q1, 1);
  const uint16x4_t sum = vadd_u16(p0q0_double, p1q1_half);
  return vcle_u16(sum, vdup_n_u16(outer_thresh));
}

// abs(p1 - p0) <= inner_thresh && abs(q1 - q0) <= inner_thresh &&
//   OuterThreshold()
static INLINE uint16x4_t NeedsFilter4(const uint16x8_t abd_p0p1_q0q1,
                                      const uint16_t inner_thresh,
                                      const uint16x4_t outer_mask) {
  const uint16x8_t a = vcleq_u16(abd_p0p1_q0q1, vdupq_n_u16(inner_thresh));
  const uint16x4_t inner_mask = vand_u16(vget_low_u16(a), vget_high_u16(a));
  return vand_u16(inner_mask, outer_mask);
}

// -----------------------------------------------------------------------------
// FilterNMasks functions.

static INLINE void Filter4Masks(const uint16x8_t p0q0, const uint16x8_t p1q1,
                                const uint16_t hev_thresh,
                                const uint16x4_t outer_mask,
                                const uint16_t inner_thresh,
                                uint16x4_t *const hev_mask,
                                uint16x4_t *const needs_filter4_mask) {
  const uint16x8_t p0p1_q0q1 = vabdq_u16(p0q0, p1q1);
  // This includes cases where NeedsFilter4() is not true and so Filter2() will
  // not be applied.
  const uint16x4_t hev_tmp_mask = Hev(p0p1_q0q1, hev_thresh);

  *needs_filter4_mask = NeedsFilter4(p0p1_q0q1, inner_thresh, outer_mask);

  // Filter2() will only be applied if both NeedsFilter4() and Hev() are true.
  *hev_mask = vand_u16(hev_tmp_mask, *needs_filter4_mask);
}

// -----------------------------------------------------------------------------
// FilterN functions.

// Calculate Filter4() or Filter2() based on |hev_mask|.
static INLINE void Filter4(const uint16x8_t p0q0, const uint16x8_t p0q1,
                           const uint16x8_t p1q1, const uint16x4_t hev_mask,
                           uint16x8_t *const p1q1_result,
                           uint16x8_t *const p0q0_result) {
  const uint16x8_t q0p1 = vextq_u16(p0q0, p1q1, 4);
  // a = 3 * (q0 - p0) + Clip3(p1 - q1, min_signed_val, max_signed_val);
  // q0mp0 means "q0 minus p0".
  const int16x8_t q0mp0_p1mq1 = vreinterpretq_s16_u16(vsubq_u16(q0p1, p0q1));
  const int16x4_t q0mp0_3 = vmul_n_s16(vget_low_s16(q0mp0_p1mq1), 3);

  // If this is for Filter2() then include |p1mq1|. Otherwise zero it.
  const int16x4_t min_signed_pixel = vdup_n_s16(-(1 << (9 /*bitdepth-1*/)));
  const int16x4_t max_signed_pixel = vdup_n_s16((1 << (9 /*bitdepth-1*/)) - 1);
  const int16x4_t p1mq1 = vget_high_s16(q0mp0_p1mq1);
  const int16x4_t p1mq1_saturated =
      Clip3S16(p1mq1, min_signed_pixel, max_signed_pixel);
  const int16x4_t hev_option =
      vand_s16(vreinterpret_s16_u16(hev_mask), p1mq1_saturated);

  const int16x4_t a = vadd_s16(q0mp0_3, hev_option);

  // Need to figure out what's going on here because there are some unnecessary
  // tricks to accommodate 8x8 as smallest 8bpp vector

  // We can not shift with rounding because the clamp comes *before* the
  // shifting. a1 = Clip3(a + 4, min_signed_val, max_signed_val) >> 3; a2 =
  // Clip3(a + 3, min_signed_val, max_signed_val) >> 3;
  const int16x4_t plus_four =
      Clip3S16(vadd_s16(a, vdup_n_s16(4)), min_signed_pixel, max_signed_pixel);
  const int16x4_t plus_three =
      Clip3S16(vadd_s16(a, vdup_n_s16(3)), min_signed_pixel, max_signed_pixel);
  const int16x4_t a1 = vshr_n_s16(plus_four, 3);
  const int16x4_t a2 = vshr_n_s16(plus_three, 3);

  // a3 = (a1 + 1) >> 1;
  const int16x4_t a3 = vrshr_n_s16(a1, 1);

  const int16x8_t a3_ma3 = vcombine_s16(a3, vneg_s16(a3));
  const int16x8_t p1q1_a3 = vaddq_s16(vreinterpretq_s16_u16(p1q1), a3_ma3);

  // Need to shift the second term or we end up with a2_ma2.
  const int16x8_t a2_ma1 = vcombine_s16(a2, vneg_s16(a1));
  const int16x8_t p0q0_a = vaddq_s16(vreinterpretq_s16_u16(p0q0), a2_ma1);
  *p1q1_result = ConvertToUnsignedPixelU16(p1q1_a3, kBitdepth10);
  *p0q0_result = ConvertToUnsignedPixelU16(p0q0_a, kBitdepth10);
}

void aom_highbd_lpf_horizontal_4_neon(uint16_t *s, int pitch,
                                      const uint8_t *blimit,
                                      const uint8_t *limit,
                                      const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_horizontal_4_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  uint16_t *const dst_p1 = (uint16_t *)(s - 2 * pitch);
  uint16_t *const dst_p0 = (uint16_t *)(s - pitch);
  uint16_t *const dst_q0 = (uint16_t *)(s);
  uint16_t *const dst_q1 = (uint16_t *)(s + pitch);

  const uint16x4_t src[4] = { vld1_u16(dst_p1), vld1_u16(dst_p0),
                              vld1_u16(dst_q0), vld1_u16(dst_q1) };

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[0], src[1], src[2], src[3], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter4_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[1], src[2]);
  const uint16x8_t p1q1 = vcombine_u16(src[0], src[3]);
  Filter4Masks(p0q0, p1q1, hev_thresh, outer_mask, inner_thresh, &hev_mask,
               &needs_filter4_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter4_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  const uint64x1_t needs_filter4_mask64 =
      vreinterpret_u64_u16(needs_filter4_mask);
  if (vget_lane_u64(needs_filter4_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t needs_filter4_mask_8 =
      vcombine_u16(needs_filter4_mask, needs_filter4_mask);

  uint16x8_t f_p1q1;
  uint16x8_t f_p0q0;
  const uint16x8_t p0q1 = vcombine_u16(src[1], src[3]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f_p1q1, &f_p0q0);

  // Already integrated the Hev mask when calculating the filtered values.
  const uint16x8_t p0q0_output = vbslq_u16(needs_filter4_mask_8, f_p0q0, p0q0);

  // p1/q1 are unmodified if only Hev() is true. This works because it was and'd
  // with |needs_filter4_mask| previously.
  const uint16x8_t p1q1_mask = veorq_u16(hev_mask_8, needs_filter4_mask_8);
  const uint16x8_t p1q1_output = vbslq_u16(p1q1_mask, f_p1q1, p1q1);

  vst1_u16(dst_p1, vget_low_u16(p1q1_output));
  vst1_u16(dst_p0, vget_low_u16(p0q0_output));
  vst1_u16(dst_q0, vget_high_u16(p0q0_output));
  vst1_u16(dst_q1, vget_high_u16(p1q1_output));
}

void aom_highbd_lpf_vertical_4_neon(uint16_t *s, int pitch,
                                    const uint8_t *blimit, const uint8_t *limit,
                                    const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_vertical_4_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  // Offset by 2 uint16_t values to load from first p1 position.
  uint16_t *dst = s - 2;
  uint16_t *dst_p1 = dst;
  uint16_t *dst_p0 = dst + pitch;
  uint16_t *dst_q0 = dst + pitch * 2;
  uint16_t *dst_q1 = dst + pitch * 3;

  uint16x4_t src[4] = { vld1_u16(dst_p1), vld1_u16(dst_p0), vld1_u16(dst_q0),
                        vld1_u16(dst_q1) };
  Transpose4x4(src);

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[0], src[1], src[2], src[3], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter4_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[1], src[2]);
  const uint16x8_t p1q1 = vcombine_u16(src[0], src[3]);
  Filter4Masks(p0q0, p1q1, hev_thresh, outer_mask, inner_thresh, &hev_mask,
               &needs_filter4_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter4_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  const uint64x1_t needs_filter4_mask64 =
      vreinterpret_u64_u16(needs_filter4_mask);
  if (vget_lane_u64(needs_filter4_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t needs_filter4_mask_8 =
      vcombine_u16(needs_filter4_mask, needs_filter4_mask);

  uint16x8_t f_p1q1;
  uint16x8_t f_p0q0;
  const uint16x8_t p0q1 = vcombine_u16(src[1], src[3]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f_p1q1, &f_p0q0);

  // Already integrated the Hev mask when calculating the filtered values.
  const uint16x8_t p0q0_output = vbslq_u16(needs_filter4_mask_8, f_p0q0, p0q0);

  // p1/q1 are unmodified if only Hev() is true. This works because it was and'd
  // with |needs_filter4_mask| previously.
  const uint16x8_t p1q1_mask = veorq_u16(hev_mask_8, needs_filter4_mask_8);
  const uint16x8_t p1q1_output = vbslq_u16(p1q1_mask, f_p1q1, p1q1);

  uint16x4_t output[4] = {
    vget_low_u16(p1q1_output),
    vget_low_u16(p0q0_output),
    vget_high_u16(p0q0_output),
    vget_high_u16(p1q1_output),
  };
  Transpose4x4(output);

  vst1_u16(dst_p1, output[0]);
  vst1_u16(dst_p0, output[1]);
  vst1_u16(dst_q0, output[2]);
  vst1_u16(dst_q1, output[3]);
}
