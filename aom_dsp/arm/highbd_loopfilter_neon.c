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

// abs(p2 - p1) <= inner_thresh && abs(p1 - p0) <= inner_thresh &&
//   abs(q1 - q0) <= inner_thresh && abs(q2 - q1) <= inner_thresh &&
//   OuterThreshold()
static INLINE uint16x4_t NeedsFilter6(const uint16x8_t abd_p0p1_q0q1,
                                      const uint16x8_t abd_p1p2_q1q2,
                                      const uint16_t inner_thresh,
                                      const uint16x4_t outer_mask) {
  const uint16x8_t a = vmaxq_u16(abd_p0p1_q0q1, abd_p1p2_q1q2);
  const uint16x8_t b = vcleq_u16(a, vdupq_n_u16(inner_thresh));
  const uint16x4_t inner_mask = vand_u16(vget_low_u16(b), vget_high_u16(b));
  return vand_u16(inner_mask, outer_mask);
}

// abs(p3 - p2) <= inner_thresh && abs(p2 - p1) <= inner_thresh &&
//   abs(p1 - p0) <= inner_thresh && abs(q1 - q0) <= inner_thresh &&
//   abs(q2 - q1) <= inner_thresh && abs(q3 - q2) <= inner_thresh
//   OuterThreshold()
static INLINE uint16x4_t NeedsFilter8(const uint16x8_t abd_p0p1_q0q1,
                                      const uint16x8_t abd_p1p2_q1q2,
                                      const uint16x8_t abd_p2p3_q2q3,
                                      const uint16_t inner_thresh,
                                      const uint16x4_t outer_mask) {
  const uint16x8_t a = vmaxq_u16(abd_p0p1_q0q1, abd_p1p2_q1q2);
  const uint16x8_t b = vmaxq_u16(a, abd_p2p3_q2q3);
  const uint16x8_t c = vcleq_u16(b, vdupq_n_u16(inner_thresh));
  const uint16x4_t inner_mask = vand_u16(vget_low_u16(c), vget_high_u16(c));
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

// abs(p1 - p0) <= flat_thresh && abs(q1 - q0) <= flat_thresh &&
//   abs(p2 - p0) <= flat_thresh && abs(q2 - q0) <= flat_thresh
// |flat_thresh| == 4 for 10 bit decode.
static INLINE uint16x4_t IsFlat3(const uint16x8_t abd_p0p1_q0q1,
                                 const uint16x8_t abd_p0p2_q0q2) {
  const int flat_thresh = 1 << 2;
  const uint16x8_t a = vmaxq_u16(abd_p0p1_q0q1, abd_p0p2_q0q2);
  const uint16x8_t b = vcleq_u16(a, vdupq_n_u16(flat_thresh));
  return vand_u16(vget_low_u16(b), vget_high_u16(b));
}

static INLINE void Filter6Masks(
    const uint16x8_t p2q2, const uint16x8_t p1q1, const uint16x8_t p0q0,
    const uint16_t hev_thresh, const uint16x4_t outer_mask,
    const uint16_t inner_thresh, uint16x4_t *const needs_filter6_mask,
    uint16x4_t *const is_flat3_mask, uint16x4_t *const hev_mask) {
  const uint16x8_t abd_p0p1_q0q1 = vabdq_u16(p0q0, p1q1);
  *hev_mask = Hev(abd_p0p1_q0q1, hev_thresh);
  *is_flat3_mask = IsFlat3(abd_p0p1_q0q1, vabdq_u16(p0q0, p2q2));
  *needs_filter6_mask = NeedsFilter6(abd_p0p1_q0q1, vabdq_u16(p1q1, p2q2),
                                     inner_thresh, outer_mask);
}

// IsFlat4 uses N=1, IsFlatOuter4 uses N=4.
// abs(p[N] - p0) <= flat_thresh && abs(q[N] - q0) <= flat_thresh &&
//   abs(p[N+1] - p0) <= flat_thresh && abs(q[N+1] - q0) <= flat_thresh &&
//   abs(p[N+2] - p0) <= flat_thresh && abs(q[N+1] - q0) <= flat_thresh
// |flat_thresh| == 4 for 10 bit decode.
static INLINE uint16x4_t IsFlat4(const uint16x8_t abd_pnp0_qnq0,
                                 const uint16x8_t abd_pn1p0_qn1q0,
                                 const uint16x8_t abd_pn2p0_qn2q0) {
  const int flat_thresh = 1 << 2;
  const uint16x8_t a = vmaxq_u16(abd_pnp0_qnq0, abd_pn1p0_qn1q0);
  const uint16x8_t b = vmaxq_u16(a, abd_pn2p0_qn2q0);
  const uint16x8_t c = vcleq_u16(b, vdupq_n_u16(flat_thresh));
  return vand_u16(vget_low_u16(c), vget_high_u16(c));
}

static INLINE void Filter8Masks(const uint16x8_t p3q3, const uint16x8_t p2q2,
                                const uint16x8_t p1q1, const uint16x8_t p0q0,
                                const uint16_t hev_thresh,
                                const uint16x4_t outer_mask,
                                const uint16_t inner_thresh,
                                uint16x4_t *const needs_filter8_mask,
                                uint16x4_t *const is_flat4_mask,
                                uint16x4_t *const hev_mask) {
  const uint16x8_t abd_p0p1_q0q1 = vabdq_u16(p0q0, p1q1);
  *hev_mask = Hev(abd_p0p1_q0q1, hev_thresh);
  const uint16x4_t is_flat4 =
      IsFlat4(abd_p0p1_q0q1, vabdq_u16(p0q0, p2q2), vabdq_u16(p0q0, p3q3));
  *needs_filter8_mask =
      NeedsFilter8(abd_p0p1_q0q1, vabdq_u16(p1q1, p2q2), vabdq_u16(p2q2, p3q3),
                   inner_thresh, outer_mask);
  // |is_flat4_mask| is used to decide where to use the result of Filter8.
  // In rare cases, |is_flat4| can be true where |needs_filter8_mask| is false,
  // overriding the question of whether to use Filter8. Because Filter4 doesn't
  // apply to p2q2, |is_flat4_mask| chooses directly between Filter8 and the
  // source value. To be correct, the mask must account for this override.
  *is_flat4_mask = vand_u16(is_flat4, *needs_filter8_mask);
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

static INLINE void Filter6(const uint16x8_t p2q2, const uint16x8_t p1q1,
                           const uint16x8_t p0q0, uint16x8_t *const p1q1_output,
                           uint16x8_t *const p0q0_output) {
  // Sum p1 and q1 output from opposite directions.
  // The formula is regrouped to allow 3 doubling operations to be combined.
  //
  // p1 = (3 * p2) + (2 * p1) + (2 * p0) + q0
  //      ^^^^^^^^
  // q1 = p0 + (2 * q0) + (2 * q1) + (3 * q2)
  //                                 ^^^^^^^^
  // p1q1 = p2q2 + 2 * (p2q2 + p1q1 + p0q0) + q0p0
  //                    ^^^^^^^^^^^
  uint16x8_t sum = vaddq_u16(p2q2, p1q1);

  // p1q1 = p2q2 + 2 * (p2q2 + p1q1 + p0q0) + q0p0
  //                                ^^^^^^
  sum = vaddq_u16(sum, p0q0);

  // p1q1 = p2q2 + 2 * (p2q2 + p1q1 + p0q0) + q0p0
  //               ^^^^^
  sum = vshlq_n_u16(sum, 1);

  // p1q1 = p2q2 + 2 * (p2q2 + p1q1 + p0q0) + q0p0
  //        ^^^^^^                          ^^^^^^
  // Should dual issue with the left shift.
  const uint16x8_t q0p0 = Transpose64(p0q0);
  const uint16x8_t outer_sum = vaddq_u16(p2q2, q0p0);
  sum = vaddq_u16(sum, outer_sum);

  *p1q1_output = vrshrq_n_u16(sum, 3);

  // Convert to p0 and q0 output:
  // p0 = p1 - (2 * p2) + q0 + q1
  // q0 = q1 - (2 * q2) + p0 + p1
  // p0q0 = p1q1 - (2 * p2q2) + q0p0 + q1p1
  //                ^^^^^^^^
  const uint16x8_t p2q2_double = vshlq_n_u16(p2q2, 1);
  // p0q0 = p1q1 - (2 * p2q2) + q0p0 + q1p1
  //        ^^^^^^^^
  sum = vsubq_u16(sum, p2q2_double);
  const uint16x8_t q1p1 = Transpose64(p1q1);
  sum = vaddq_u16(sum, vaddq_u16(q0p0, q1p1));

  *p0q0_output = vrshrq_n_u16(sum, 3);
}

void aom_highbd_lpf_horizontal_6_neon(uint16_t *s, int pitch,
                                      const uint8_t *blimit,
                                      const uint8_t *limit,
                                      const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_horizontal_6_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  uint16_t *const dst_p2 = s - 3 * pitch;
  uint16_t *const dst_p1 = s - 2 * pitch;
  uint16_t *const dst_p0 = s - pitch;
  uint16_t *const dst_q0 = s;
  uint16_t *const dst_q1 = s + pitch;
  uint16_t *const dst_q2 = s + 2 * pitch;

  const uint16x4_t src[6] = { vld1_u16(dst_p2), vld1_u16(dst_p1),
                              vld1_u16(dst_p0), vld1_u16(dst_q0),
                              vld1_u16(dst_q1), vld1_u16(dst_q2) };

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[1], src[2], src[3], src[4], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter_mask;
  uint16x4_t is_flat3_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[2], src[3]);
  const uint16x8_t p1q1 = vcombine_u16(src[1], src[4]);
  const uint16x8_t p2q2 = vcombine_u16(src[0], src[5]);
  Filter6Masks(p2q2, p1q1, p0q0, hev_thresh, outer_mask, inner_thresh,
               &needs_filter_mask, &is_flat3_mask, &hev_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  // This might be faster than vaddv (latency 3) because mov to general register
  // has latency 2.
  const uint64x1_t needs_filter_mask64 =
      vreinterpret_u64_u16(needs_filter_mask);
  if (vget_lane_u64(needs_filter_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t is_flat3_mask_8 = vcombine_u16(is_flat3_mask, is_flat3_mask);
  const uint16x8_t needs_filter_mask_8 =
      vcombine_u16(needs_filter_mask, needs_filter_mask);

  uint16x8_t f4_p1q1;
  uint16x8_t f4_p0q0;
  // ZIP1 p0q0, p1q1 may perform better here.
  const uint16x8_t p0q1 = vcombine_u16(src[2], src[4]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f4_p1q1, &f4_p0q0);
  f4_p1q1 = vbslq_u16(hev_mask_8, p1q1, f4_p1q1);

  uint16x8_t p0q0_output, p1q1_output;
  // Because we did not return after testing |needs_filter_mask| we know it is
  // nonzero. |is_flat3_mask| controls whether the needed filter is Filter4 or
  // Filter6. Therefore if it is false when |needs_filter_mask| is true, Filter6
  // output is not used.
  uint16x8_t f6_p1q1, f6_p0q0;
  const uint64x1_t need_filter6 = vreinterpret_u64_u16(is_flat3_mask);
  if (vget_lane_u64(need_filter6, 0) == 0) {
    // Filter6() does not apply, but Filter4() applies to one or more values.
    p0q0_output = p0q0;
    p1q1_output = vbslq_u16(needs_filter_mask_8, f4_p1q1, p1q1);
    p0q0_output = vbslq_u16(needs_filter_mask_8, f4_p0q0, p0q0);
  } else {
    Filter6(p2q2, p1q1, p0q0, &f6_p1q1, &f6_p0q0);
    p1q1_output = vbslq_u16(is_flat3_mask_8, f6_p1q1, f4_p1q1);
    p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
    p0q0_output = vbslq_u16(is_flat3_mask_8, f6_p0q0, f4_p0q0);
    p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
  }

  vst1_u16(dst_p1, vget_low_u16(p1q1_output));
  vst1_u16(dst_p0, vget_low_u16(p0q0_output));
  vst1_u16(dst_q0, vget_high_u16(p0q0_output));
  vst1_u16(dst_q1, vget_high_u16(p1q1_output));
}

void aom_highbd_lpf_vertical_6_neon(uint16_t *s, int pitch,
                                    const uint8_t *blimit, const uint8_t *limit,
                                    const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_vertical_6_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  // Left side of the filter window.
  uint16_t *const dst = s - 3;
  uint16_t *const dst_0 = dst;
  uint16_t *const dst_1 = dst + pitch;
  uint16_t *const dst_2 = dst + 2 * pitch;
  uint16_t *const dst_3 = dst + 3 * pitch;

  // Overread by 2 values. These overreads become the high halves of src_raw[2]
  // and src_raw[3] after transpose.
  uint16x8_t src_raw[4] = { vld1q_u16(dst_0), vld1q_u16(dst_1),
                            vld1q_u16(dst_2), vld1q_u16(dst_3) };
  Transpose4x8(src_raw);
  // p2, p1, p0, q0, q1, q2
  const uint16x4_t src[6] = {
    vget_low_u16(src_raw[0]),  vget_low_u16(src_raw[1]),
    vget_low_u16(src_raw[2]),  vget_low_u16(src_raw[3]),
    vget_high_u16(src_raw[0]), vget_high_u16(src_raw[1]),
  };

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[1], src[2], src[3], src[4], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter_mask;
  uint16x4_t is_flat3_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[2], src[3]);
  const uint16x8_t p1q1 = vcombine_u16(src[1], src[4]);
  const uint16x8_t p2q2 = vcombine_u16(src[0], src[5]);
  Filter6Masks(p2q2, p1q1, p0q0, hev_thresh, outer_mask, inner_thresh,
               &needs_filter_mask, &is_flat3_mask, &hev_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  // This might be faster than vaddv (latency 3) because mov to general register
  // has latency 2.
  const uint64x1_t needs_filter_mask64 =
      vreinterpret_u64_u16(needs_filter_mask);
  if (vget_lane_u64(needs_filter_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t is_flat3_mask_8 = vcombine_u16(is_flat3_mask, is_flat3_mask);
  const uint16x8_t needs_filter_mask_8 =
      vcombine_u16(needs_filter_mask, needs_filter_mask);

  uint16x8_t f4_p1q1;
  uint16x8_t f4_p0q0;
  // ZIP1 p0q0, p1q1 may perform better here.
  const uint16x8_t p0q1 = vcombine_u16(src[2], src[4]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f4_p1q1, &f4_p0q0);
  f4_p1q1 = vbslq_u16(hev_mask_8, p1q1, f4_p1q1);

  uint16x8_t p0q0_output, p1q1_output;
  // Because we did not return after testing |needs_filter_mask| we know it is
  // nonzero. |is_flat3_mask| controls whether the needed filter is Filter4 or
  // Filter6. Therefore if it is false when |needs_filter_mask| is true, Filter6
  // output is not used.
  uint16x8_t f6_p1q1, f6_p0q0;
  const uint64x1_t need_filter6 = vreinterpret_u64_u16(is_flat3_mask);
  if (vget_lane_u64(need_filter6, 0) == 0) {
    // Filter6() does not apply, but Filter4() applies to one or more values.
    p0q0_output = p0q0;
    p1q1_output = vbslq_u16(needs_filter_mask_8, f4_p1q1, p1q1);
    p0q0_output = vbslq_u16(needs_filter_mask_8, f4_p0q0, p0q0);
  } else {
    Filter6(p2q2, p1q1, p0q0, &f6_p1q1, &f6_p0q0);
    p1q1_output = vbslq_u16(is_flat3_mask_8, f6_p1q1, f4_p1q1);
    p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
    p0q0_output = vbslq_u16(is_flat3_mask_8, f6_p0q0, f4_p0q0);
    p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
  }

  uint16x4_t output[4] = {
    vget_low_u16(p1q1_output),
    vget_low_u16(p0q0_output),
    vget_high_u16(p0q0_output),
    vget_high_u16(p1q1_output),
  };
  Transpose4x4(output);

  // dst_n starts at p2, so adjust to p1.
  vst1_u16(dst_0 + 1, output[0]);
  vst1_u16(dst_1 + 1, output[1]);
  vst1_u16(dst_2 + 1, output[2]);
  vst1_u16(dst_3 + 1, output[3]);
}

static INLINE void Filter8(const uint16x8_t p3q3, const uint16x8_t p2q2,
                           const uint16x8_t p1q1, const uint16x8_t p0q0,
                           uint16x8_t *const p2q2_output,
                           uint16x8_t *const p1q1_output,
                           uint16x8_t *const p0q0_output) {
  // Sum p2 and q2 output from opposite directions.
  // The formula is regrouped to allow 2 doubling operations to be combined.
  // p2 = (3 * p3) + (2 * p2) + p1 + p0 + q0
  //      ^^^^^^^^
  // q2 = p0 + q0 + q1 + (2 * q2) + (3 * q3)
  //                                ^^^^^^^^
  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //                    ^^^^^^^^^^^
  const uint16x8_t p23q23 = vaddq_u16(p3q3, p2q2);

  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //               ^^^^^
  uint16x8_t sum = vshlq_n_u16(p23q23, 1);

  // Add two other terms to make dual issue with shift more likely.
  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //                                   ^^^^^^^^^^^
  const uint16x8_t p01q01 = vaddq_u16(p0q0, p1q1);

  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //                                 ^^^^^^^^^^^^^
  sum = vaddq_u16(sum, p01q01);

  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //        ^^^^^^
  sum = vaddq_u16(sum, p3q3);

  // p2q2 = p3q3 + 2 * (p3q3 + p2q2) + p1q1 + p0q0 + q0p0
  //                                               ^^^^^^
  const uint16x8_t q0p0 = Transpose64(p0q0);
  sum = vaddq_u16(sum, q0p0);

  *p2q2_output = vrshrq_n_u16(sum, 3);

  // Convert to p1 and q1 output:
  // p1 = p2 - p3 - p2 + p1 + q1
  // q1 = q2 - q3 - q2 + q0 + p1
  sum = vsubq_u16(sum, p23q23);
  const uint16x8_t q1p1 = Transpose64(p1q1);
  sum = vaddq_u16(sum, vaddq_u16(p1q1, q1p1));

  *p1q1_output = vrshrq_n_u16(sum, 3);

  // Convert to p0 and q0 output:
  // p0 = p1 - p3 - p1 + p0 + q2
  // q0 = q1 - q3 - q1 + q0 + p2
  sum = vsubq_u16(sum, vaddq_u16(p3q3, p1q1));
  const uint16x8_t q2p2 = Transpose64(p2q2);
  sum = vaddq_u16(sum, vaddq_u16(p0q0, q2p2));

  *p0q0_output = vrshrq_n_u16(sum, 3);
}

void aom_highbd_lpf_horizontal_8_neon(uint16_t *s, int pitch,
                                      const uint8_t *blimit,
                                      const uint8_t *limit,
                                      const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_horizontal_8_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  uint16_t *const dst_p3 = s - 4 * pitch;
  uint16_t *const dst_p2 = s - 3 * pitch;
  uint16_t *const dst_p1 = s - 2 * pitch;
  uint16_t *const dst_p0 = s - pitch;
  uint16_t *const dst_q0 = s;
  uint16_t *const dst_q1 = s + pitch;
  uint16_t *const dst_q2 = s + 2 * pitch;
  uint16_t *const dst_q3 = s + 3 * pitch;

  const uint16x4_t src[8] = { vld1_u16(dst_p3), vld1_u16(dst_p2),
                              vld1_u16(dst_p1), vld1_u16(dst_p0),
                              vld1_u16(dst_q0), vld1_u16(dst_q1),
                              vld1_u16(dst_q2), vld1_u16(dst_q3) };

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[2], src[3], src[4], src[5], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter_mask;
  uint16x4_t is_flat4_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[3], src[4]);
  const uint16x8_t p1q1 = vcombine_u16(src[2], src[5]);
  const uint16x8_t p2q2 = vcombine_u16(src[1], src[6]);
  const uint16x8_t p3q3 = vcombine_u16(src[0], src[7]);
  Filter8Masks(p3q3, p2q2, p1q1, p0q0, hev_thresh, outer_mask, inner_thresh,
               &needs_filter_mask, &is_flat4_mask, &hev_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  // This might be faster than vaddv (latency 3) because mov to general register
  // has latency 2.
  const uint64x1_t needs_filter_mask64 =
      vreinterpret_u64_u16(needs_filter_mask);
  if (vget_lane_u64(needs_filter_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t needs_filter_mask_8 =
      vcombine_u16(needs_filter_mask, needs_filter_mask);

  uint16x8_t f4_p1q1;
  uint16x8_t f4_p0q0;
  // ZIP1 p0q0, p1q1 may perform better here.
  const uint16x8_t p0q1 = vcombine_u16(src[3], src[5]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f4_p1q1, &f4_p0q0);
  f4_p1q1 = vbslq_u16(hev_mask_8, p1q1, f4_p1q1);

  uint16x8_t p0q0_output, p1q1_output, p2q2_output;
  // Because we did not return after testing |needs_filter_mask| we know it is
  // nonzero. |is_flat4_mask| controls whether the needed filter is Filter4 or
  // Filter8. Therefore if it is false when |needs_filter_mask| is true, Filter8
  // output is not used.
  uint16x8_t f8_p2q2, f8_p1q1, f8_p0q0;
  const uint64x1_t need_filter8 = vreinterpret_u64_u16(is_flat4_mask);
  if (vget_lane_u64(need_filter8, 0) == 0) {
    // Filter8() does not apply, but Filter4() applies to one or more values.
    p2q2_output = p2q2;
    p1q1_output = vbslq_u16(needs_filter_mask_8, f4_p1q1, p1q1);
    p0q0_output = vbslq_u16(needs_filter_mask_8, f4_p0q0, p0q0);
  } else {
    const uint16x8_t is_flat4_mask_8 =
        vcombine_u16(is_flat4_mask, is_flat4_mask);
    Filter8(p3q3, p2q2, p1q1, p0q0, &f8_p2q2, &f8_p1q1, &f8_p0q0);
    p2q2_output = vbslq_u16(is_flat4_mask_8, f8_p2q2, p2q2);
    p1q1_output = vbslq_u16(is_flat4_mask_8, f8_p1q1, f4_p1q1);
    p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
    p0q0_output = vbslq_u16(is_flat4_mask_8, f8_p0q0, f4_p0q0);
    p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
  }

  vst1_u16(dst_p2, vget_low_u16(p2q2_output));
  vst1_u16(dst_p1, vget_low_u16(p1q1_output));
  vst1_u16(dst_p0, vget_low_u16(p0q0_output));
  vst1_u16(dst_q0, vget_high_u16(p0q0_output));
  vst1_u16(dst_q1, vget_high_u16(p1q1_output));
  vst1_u16(dst_q2, vget_high_u16(p2q2_output));
}

static INLINE uint16x8_t ReverseLowHalf(const uint16x8_t a) {
  return vcombine_u16(vrev64_u16(vget_low_u16(a)), vget_high_u16(a));
}

void aom_highbd_lpf_vertical_8_neon(uint16_t *s, int pitch,
                                    const uint8_t *blimit, const uint8_t *limit,
                                    const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_vertical_8_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  uint16_t *const dst = s - 4;
  uint16_t *const dst_0 = dst;
  uint16_t *const dst_1 = dst + pitch;
  uint16_t *const dst_2 = dst + 2 * pitch;
  uint16_t *const dst_3 = dst + 3 * pitch;

  // src_raw[n] contains p3, p2, p1, p0, q0, q1, q2, q3 for row n.
  // To get desired pairs after transpose, one half should be reversed.
  uint16x8_t src[4] = { vld1q_u16(dst_0), vld1q_u16(dst_1), vld1q_u16(dst_2),
                        vld1q_u16(dst_3) };

  // src[0] = p0q0
  // src[1] = p1q1
  // src[2] = p2q2
  // src[3] = p3q3
  LoopFilterTranspose4x8(src);

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask = OuterThreshold(
      vget_low_u16(src[1]), vget_low_u16(src[0]), vget_high_u16(src[0]),
      vget_high_u16(src[1]), outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter_mask;
  uint16x4_t is_flat4_mask;
  const uint16x8_t p0q0 = src[0];
  const uint16x8_t p1q1 = src[1];
  const uint16x8_t p2q2 = src[2];
  const uint16x8_t p3q3 = src[3];
  Filter8Masks(p3q3, p2q2, p1q1, p0q0, hev_thresh, outer_mask, inner_thresh,
               &needs_filter_mask, &is_flat4_mask, &hev_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  // This might be faster than vaddv (latency 3) because mov to general register
  // has latency 2.
  const uint64x1_t needs_filter_mask64 =
      vreinterpret_u64_u16(needs_filter_mask);
  if (vget_lane_u64(needs_filter_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)

  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t needs_filter_mask_8 =
      vcombine_u16(needs_filter_mask, needs_filter_mask);

  uint16x8_t f4_p1q1;
  uint16x8_t f4_p0q0;
  const uint16x8_t p0q1 = vcombine_u16(vget_low_u16(p0q0), vget_high_u16(p1q1));
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f4_p1q1, &f4_p0q0);
  f4_p1q1 = vbslq_u16(hev_mask_8, p1q1, f4_p1q1);

  uint16x8_t p0q0_output, p1q1_output, p2q2_output;
  // Because we did not return after testing |needs_filter_mask| we know it is
  // nonzero. |is_flat4_mask| controls whether the needed filter is Filter4 or
  // Filter8. Therefore if it is false when |needs_filter_mask| is true, Filter8
  // output is not used.
  const uint64x1_t need_filter8 = vreinterpret_u64_u16(is_flat4_mask);
  if (vget_lane_u64(need_filter8, 0) == 0) {
    // Filter8() does not apply, but Filter4() applies to one or more values.
    p2q2_output = p2q2;
    p1q1_output = vbslq_u16(needs_filter_mask_8, f4_p1q1, p1q1);
    p0q0_output = vbslq_u16(needs_filter_mask_8, f4_p0q0, p0q0);
  } else {
    const uint16x8_t is_flat4_mask_8 =
        vcombine_u16(is_flat4_mask, is_flat4_mask);
    uint16x8_t f8_p2q2, f8_p1q1, f8_p0q0;
    Filter8(p3q3, p2q2, p1q1, p0q0, &f8_p2q2, &f8_p1q1, &f8_p0q0);
    p2q2_output = vbslq_u16(is_flat4_mask_8, f8_p2q2, p2q2);
    p1q1_output = vbslq_u16(is_flat4_mask_8, f8_p1q1, f4_p1q1);
    p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
    p0q0_output = vbslq_u16(is_flat4_mask_8, f8_p0q0, f4_p0q0);
    p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
  }

  uint16x8_t output[4] = { p0q0_output, p1q1_output, p2q2_output, p3q3 };
  // After transpose, |output| will contain rows of the form:
  // p0 p1 p2 p3 q0 q1 q2 q3
  Transpose4x8(output);

  // Reverse p values to produce original order:
  // p3 p2 p1 p0 q0 q1 q2 q3
  vst1q_u16(dst_0, ReverseLowHalf(output[0]));
  vst1q_u16(dst_1, ReverseLowHalf(output[1]));
  vst1q_u16(dst_2, ReverseLowHalf(output[2]));
  vst1q_u16(dst_3, ReverseLowHalf(output[3]));
}

static INLINE void Filter14(
    const uint16x8_t p6q6, const uint16x8_t p5q5, const uint16x8_t p4q4,
    const uint16x8_t p3q3, const uint16x8_t p2q2, const uint16x8_t p1q1,
    const uint16x8_t p0q0, uint16x8_t *const p5q5_output,
    uint16x8_t *const p4q4_output, uint16x8_t *const p3q3_output,
    uint16x8_t *const p2q2_output, uint16x8_t *const p1q1_output,
    uint16x8_t *const p0q0_output) {
  // Sum p5 and q5 output from opposite directions.
  // p5 = (7 * p6) + (2 * p5) + (2 * p4) + p3 + p2 + p1 + p0 + q0
  //      ^^^^^^^^
  // q5 = p0 + q0 + q1 + q2 + q3 + (2 * q4) + (2 * q5) + (7 * q6)
  //                                                     ^^^^^^^^
  const uint16x8_t p6q6_x7 = vsubq_u16(vshlq_n_u16(p6q6, 3), p6q6);

  // p5 = (7 * p6) + (2 * p5) + (2 * p4) + p3 + p2 + p1 + p0 + q0
  //                 ^^^^^^^^^^^^^^^^^^^
  // q5 = p0 + q0 + q1 + q2 + q3 + (2 * q4) + (2 * q5) + (7 * q6)
  //                               ^^^^^^^^^^^^^^^^^^^
  uint16x8_t sum = vshlq_n_u16(vaddq_u16(p5q5, p4q4), 1);
  sum = vaddq_u16(sum, p6q6_x7);

  // p5 = (7 * p6) + (2 * p5) + (2 * p4) + p3 + p2 + p1 + p0 + q0
  //                                       ^^^^^^^
  // q5 = p0 + q0 + q1 + q2 + q3 + (2 * q4) + (2 * q5) + (7 * q6)
  //                     ^^^^^^^
  sum = vaddq_u16(vaddq_u16(p3q3, p2q2), sum);

  // p5 = (7 * p6) + (2 * p5) + (2 * p4) + p3 + p2 + p1 + p0 + q0
  //                                                 ^^^^^^^
  // q5 = p0 + q0 + q1 + q2 + q3 + (2 * q4) + (2 * q5) + (7 * q6)
  //           ^^^^^^^
  sum = vaddq_u16(vaddq_u16(p1q1, p0q0), sum);

  // p5 = (7 * p6) + (2 * p5) + (2 * p4) + p3 + p2 + p1 + p0 + q0
  //                                                           ^^
  // q5 = p0 + q0 + q1 + q2 + q3 + (2 * q4) + (2 * q5) + (7 * q6)
  //      ^^
  const uint16x8_t q0p0 = Transpose64(p0q0);
  sum = vaddq_u16(sum, q0p0);

  *p5q5_output = vrshrq_n_u16(sum, 4);

  // Convert to p4 and q4 output:
  // p4 = p5 - (2 * p6) + p3 + q1
  // q4 = q5 - (2 * q6) + q3 + p1
  sum = vsubq_u16(sum, vshlq_n_u16(p6q6, 1));
  const uint16x8_t q1p1 = Transpose64(p1q1);
  sum = vaddq_u16(vaddq_u16(p3q3, q1p1), sum);

  *p4q4_output = vrshrq_n_u16(sum, 4);

  // Convert to p3 and q3 output:
  // p3 = p4 - p6 - p5 + p2 + q2
  // q3 = q4 - q6 - q5 + q2 + p2
  sum = vsubq_u16(sum, vaddq_u16(p6q6, p5q5));
  const uint16x8_t q2p2 = Transpose64(p2q2);
  sum = vaddq_u16(vaddq_u16(p2q2, q2p2), sum);

  *p3q3_output = vrshrq_n_u16(sum, 4);

  // Convert to p2 and q2 output:
  // p2 = p3 - p6 - p4 + p1 + q3
  // q2 = q3 - q6 - q4 + q1 + p3
  sum = vsubq_u16(sum, vaddq_u16(p6q6, p4q4));
  const uint16x8_t q3p3 = Transpose64(p3q3);
  sum = vaddq_u16(vaddq_u16(p1q1, q3p3), sum);

  *p2q2_output = vrshrq_n_u16(sum, 4);

  // Convert to p1 and q1 output:
  // p1 = p2 - p6 - p3 + p0 + q4
  // q1 = q2 - q6 - q3 + q0 + p4
  sum = vsubq_u16(sum, vaddq_u16(p6q6, p3q3));
  const uint16x8_t q4p4 = Transpose64(p4q4);
  sum = vaddq_u16(vaddq_u16(p0q0, q4p4), sum);

  *p1q1_output = vrshrq_n_u16(sum, 4);

  // Convert to p0 and q0 output:
  // p0 = p1 - p6 - p2 + q0 + q5
  // q0 = q1 - q6 - q2 + p0 + p5
  sum = vsubq_u16(sum, vaddq_u16(p6q6, p2q2));
  const uint16x8_t q5p5 = Transpose64(p5q5);
  sum = vaddq_u16(vaddq_u16(q0p0, q5p5), sum);

  *p0q0_output = vrshrq_n_u16(sum, 4);
}

void aom_highbd_lpf_horizontal_14_neon(uint16_t *s, int pitch,
                                       const uint8_t *blimit,
                                       const uint8_t *limit,
                                       const uint8_t *thresh, int bd) {
  // TODO(b/217462944): add support for 8/12-bit.
  if (bd != 10) {
    aom_highbd_lpf_horizontal_14_c(s, pitch, blimit, limit, thresh, bd);
    return;
  }
  uint16_t *const dst_p6 = s - 7 * pitch;
  uint16_t *const dst_p5 = s - 6 * pitch;
  uint16_t *const dst_p4 = s - 5 * pitch;
  uint16_t *const dst_p3 = s - 4 * pitch;
  uint16_t *const dst_p2 = s - 3 * pitch;
  uint16_t *const dst_p1 = s - 2 * pitch;
  uint16_t *const dst_p0 = s - pitch;
  uint16_t *const dst_q0 = s;
  uint16_t *const dst_q1 = s + pitch;
  uint16_t *const dst_q2 = s + 2 * pitch;
  uint16_t *const dst_q3 = s + 3 * pitch;
  uint16_t *const dst_q4 = s + 4 * pitch;
  uint16_t *const dst_q5 = s + 5 * pitch;
  uint16_t *const dst_q6 = s + 6 * pitch;

  const uint16x4_t src[14] = {
    vld1_u16(dst_p6), vld1_u16(dst_p5), vld1_u16(dst_p4), vld1_u16(dst_p3),
    vld1_u16(dst_p2), vld1_u16(dst_p1), vld1_u16(dst_p0), vld1_u16(dst_q0),
    vld1_u16(dst_q1), vld1_u16(dst_q2), vld1_u16(dst_q3), vld1_u16(dst_q4),
    vld1_u16(dst_q5), vld1_u16(dst_q6)
  };

  // Adjust thresholds to bitdepth.
  const int outer_thresh = *blimit << 2;
  const int inner_thresh = *limit << 2;
  const int hev_thresh = *thresh << 2;
  const uint16x4_t outer_mask =
      OuterThreshold(src[5], src[6], src[7], src[8], outer_thresh);
  uint16x4_t hev_mask;
  uint16x4_t needs_filter_mask;
  uint16x4_t is_flat4_mask;
  const uint16x8_t p0q0 = vcombine_u16(src[6], src[7]);
  const uint16x8_t p1q1 = vcombine_u16(src[5], src[8]);
  const uint16x8_t p2q2 = vcombine_u16(src[4], src[9]);
  const uint16x8_t p3q3 = vcombine_u16(src[3], src[10]);
  Filter8Masks(p3q3, p2q2, p1q1, p0q0, hev_thresh, outer_mask, inner_thresh,
               &needs_filter_mask, &is_flat4_mask, &hev_mask);

#if defined(__aarch64__)
  if (vaddv_u16(needs_filter_mask) == 0) {
    // None of the values will be filtered.
    return;
  }
#else   // !defined(__aarch64__)
  // This might be faster than vaddv (latency 3) because mov to general register
  // has latency 2.
  const uint64x1_t needs_filter_mask64 =
      vreinterpret_u64_u16(needs_filter_mask);
  if (vget_lane_u64(needs_filter_mask64, 0) == 0) {
    // None of the values will be filtered.
    return;
  }
#endif  // defined(__aarch64__)
  const uint16x8_t p4q4 = vcombine_u16(src[2], src[11]);
  const uint16x8_t p5q5 = vcombine_u16(src[1], src[12]);
  const uint16x8_t p6q6 = vcombine_u16(src[0], src[13]);
  // Mask to choose between the outputs of Filter8 and Filter14.
  // As with the derivation of |is_flat4_mask|, the question of whether to use
  // Filter14 is only raised where |is_flat4_mask| is true.
  const uint16x4_t is_flat4_outer_mask = vand_u16(
      is_flat4_mask, IsFlat4(vabdq_u16(p0q0, p4q4), vabdq_u16(p0q0, p5q5),
                             vabdq_u16(p0q0, p6q6)));
  // Copy the masks to the high bits for packed comparisons later.
  const uint16x8_t hev_mask_8 = vcombine_u16(hev_mask, hev_mask);
  const uint16x8_t needs_filter_mask_8 =
      vcombine_u16(needs_filter_mask, needs_filter_mask);

  uint16x8_t f4_p1q1;
  uint16x8_t f4_p0q0;
  // ZIP1 p0q0, p1q1 may perform better here.
  const uint16x8_t p0q1 = vcombine_u16(src[6], src[8]);
  Filter4(p0q0, p0q1, p1q1, hev_mask, &f4_p1q1, &f4_p0q0);
  f4_p1q1 = vbslq_u16(hev_mask_8, p1q1, f4_p1q1);

  uint16x8_t p0q0_output, p1q1_output, p2q2_output, p3q3_output, p4q4_output,
      p5q5_output;
  // Because we did not return after testing |needs_filter_mask| we know it is
  // nonzero. |is_flat4_mask| controls whether the needed filter is Filter4 or
  // Filter8. Therefore if it is false when |needs_filter_mask| is true, Filter8
  // output is not used.
  uint16x8_t f8_p2q2, f8_p1q1, f8_p0q0;
  const uint64x1_t need_filter8 = vreinterpret_u64_u16(is_flat4_mask);
  if (vget_lane_u64(need_filter8, 0) == 0) {
    // Filter8() and Filter14() do not apply, but Filter4() applies to one or
    // more values.
    p5q5_output = p5q5;
    p4q4_output = p4q4;
    p3q3_output = p3q3;
    p2q2_output = p2q2;
    p1q1_output = vbslq_u16(needs_filter_mask_8, f4_p1q1, p1q1);
    p0q0_output = vbslq_u16(needs_filter_mask_8, f4_p0q0, p0q0);
  } else {
    const uint16x8_t use_filter8_mask =
        vcombine_u16(is_flat4_mask, is_flat4_mask);
    Filter8(p3q3, p2q2, p1q1, p0q0, &f8_p2q2, &f8_p1q1, &f8_p0q0);
    const uint64x1_t need_filter14 = vreinterpret_u64_u16(is_flat4_outer_mask);
    if (vget_lane_u64(need_filter14, 0) == 0) {
      // Filter14() does not apply, but Filter8() and Filter4() apply to one or
      // more values.
      p5q5_output = p5q5;
      p4q4_output = p4q4;
      p3q3_output = p3q3;
      p2q2_output = vbslq_u16(use_filter8_mask, f8_p2q2, p2q2);
      p1q1_output = vbslq_u16(use_filter8_mask, f8_p1q1, f4_p1q1);
      p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
      p0q0_output = vbslq_u16(use_filter8_mask, f8_p0q0, f4_p0q0);
      p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
    } else {
      // All filters may contribute values to final outputs.
      const uint16x8_t use_filter14_mask =
          vcombine_u16(is_flat4_outer_mask, is_flat4_outer_mask);
      uint16x8_t f14_p5q5, f14_p4q4, f14_p3q3, f14_p2q2, f14_p1q1, f14_p0q0;
      Filter14(p6q6, p5q5, p4q4, p3q3, p2q2, p1q1, p0q0, &f14_p5q5, &f14_p4q4,
               &f14_p3q3, &f14_p2q2, &f14_p1q1, &f14_p0q0);
      p5q5_output = vbslq_u16(use_filter14_mask, f14_p5q5, p5q5);
      p4q4_output = vbslq_u16(use_filter14_mask, f14_p4q4, p4q4);
      p3q3_output = vbslq_u16(use_filter14_mask, f14_p3q3, p3q3);
      p2q2_output = vbslq_u16(use_filter14_mask, f14_p2q2, f8_p2q2);
      p2q2_output = vbslq_u16(use_filter8_mask, p2q2_output, p2q2);
      p2q2_output = vbslq_u16(needs_filter_mask_8, p2q2_output, p2q2);
      p1q1_output = vbslq_u16(use_filter14_mask, f14_p1q1, f8_p1q1);
      p1q1_output = vbslq_u16(use_filter8_mask, p1q1_output, f4_p1q1);
      p1q1_output = vbslq_u16(needs_filter_mask_8, p1q1_output, p1q1);
      p0q0_output = vbslq_u16(use_filter14_mask, f14_p0q0, f8_p0q0);
      p0q0_output = vbslq_u16(use_filter8_mask, p0q0_output, f4_p0q0);
      p0q0_output = vbslq_u16(needs_filter_mask_8, p0q0_output, p0q0);
    }
  }

  vst1_u16(dst_p5, vget_low_u16(p5q5_output));
  vst1_u16(dst_p4, vget_low_u16(p4q4_output));
  vst1_u16(dst_p3, vget_low_u16(p3q3_output));
  vst1_u16(dst_p2, vget_low_u16(p2q2_output));
  vst1_u16(dst_p1, vget_low_u16(p1q1_output));
  vst1_u16(dst_p0, vget_low_u16(p0q0_output));
  vst1_u16(dst_q0, vget_high_u16(p0q0_output));
  vst1_u16(dst_q1, vget_high_u16(p1q1_output));
  vst1_u16(dst_q2, vget_high_u16(p2q2_output));
  vst1_u16(dst_q3, vget_high_u16(p3q3_output));
  vst1_u16(dst_q4, vget_high_u16(p4q4_output));
  vst1_u16(dst_q5, vget_high_u16(p5q5_output));
}
