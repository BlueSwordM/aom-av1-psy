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

#include "aom_dsp/arm/mem_neon.h"

#include "av1/common/quant_common.h"
#include "av1/encoder/av1_quantize.h"

static INLINE uint16x4_t quantize_4(
    const tran_low_t *coeff_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, int32x4_t v_quant_s32, int32x4_t v_dequant_s32,
    int32x4_t v_round_s32, int log_scale, int shift) {
  const int32x4_t v_coeff = vld1q_s32(coeff_ptr);
  const int32x4_t v_coeff_sign =
      vreinterpretq_s32_u32(vcltq_s32(v_coeff, vdupq_n_s32(0)));
  const int32x4_t v_abs_coeff = vabsq_s32(v_coeff);
  // ((abs_coeff << (1 + log_scale)) >= dequant_ptr[rc01])
  const int32x4_t v_abs_coeff_scaled =
      vshlq_s32(v_abs_coeff, vdupq_n_s32(1 + log_scale));
  const uint32x4_t v_mask = vcgeq_s32(v_abs_coeff_scaled, v_dequant_s32);
  // const int64_t tmp = (int64_t)abs_coeff + log_scaled_round;
  const int32x4_t v_tmp = vaddq_s32(v_abs_coeff, v_round_s32);
  //  const int abs_qcoeff = (int)((tmp * quant) >> shift);
  const int64x2_t v_abs_qcoeff_lo =
      vmull_s32(vget_low_s32(v_tmp), vget_low_s32(v_quant_s32));
#ifdef __aarch64__
  const int64x2_t v_abs_qcoeff_hi = vmull_high_s32(v_tmp, v_quant_s32);
#else
  const int64x2_t v_abs_qcoeff_hi =
      vmull_s32(vget_high_s32(v_tmp), vget_high_s32(v_quant_s32));
#endif
  // vshlq_s64 will shift right if shift value is negative.
  const int64x2_t v_shift = vdupq_n_s64(-shift);
  const int64x2_t v_abs_qcoeff_lo_sh = vshlq_s64(v_abs_qcoeff_lo, v_shift);
  const int64x2_t v_abs_qcoeff_hi_sh = vshlq_s64(v_abs_qcoeff_hi, v_shift);
  const int32x4_t v_abs_qcoeff = vcombine_s32(vmovn_s64(v_abs_qcoeff_lo_sh),
                                              vmovn_s64(v_abs_qcoeff_hi_sh));
  //  qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
  const int32x4_t v_qcoeff =
      vsubq_s32(veorq_s32(v_abs_qcoeff, v_coeff_sign), v_coeff_sign);
  //  const tran_low_t abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
  const int64x2_t v_abs_dqcoeff_lo =
      vmull_s32(vget_low_s32(v_abs_qcoeff), vget_low_s32(v_dequant_s32));
#ifdef __aarch64__
  const int64x2_t v_abs_dqcoeff_hi =
      vmull_high_s32(v_abs_qcoeff, v_dequant_s32);
#else
  const int64x2_t v_abs_dqcoeff_hi =
      vmull_s32(vget_high_s32(v_abs_qcoeff), vget_high_s32(v_dequant_s32));
#endif
  // vshlq_s64 will shift right if shift value is negative.
  const int64x2_t v_log_scale = vdupq_n_s64(-log_scale);
  const int64x2_t v_abs_dqcoeff_lo_sh =
      vshlq_s64(v_abs_dqcoeff_lo, v_log_scale);
  const int64x2_t v_abs_dqcoeff_hi_sh =
      vshlq_s64(v_abs_dqcoeff_hi, v_log_scale);
  const int32x4_t v_abs_dqcoeff = vcombine_s32(vmovn_s64(v_abs_dqcoeff_lo_sh),
                                               vmovn_s64(v_abs_dqcoeff_hi_sh));
  //  dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
  const int32x4_t v_dqcoeff =
      vsubq_s32(veorq_s32(v_abs_dqcoeff, v_coeff_sign), v_coeff_sign);

  vst1q_s32(qcoeff_ptr, vandq_s32(v_qcoeff, vreinterpretq_s32_u32(v_mask)));
  vst1q_s32(dqcoeff_ptr, vandq_s32(v_dqcoeff, vreinterpretq_s32_u32(v_mask)));

  // Used to find eob.
  return vmovn_u32(v_mask);
}

static INLINE int16x8_t get_max_lane_eob(const int16_t *iscan,
                                         int16x8_t v_eobmax,
                                         uint16x8_t v_mask) {
  const int16x8_t v_iscan = vld1q_s16(&iscan[0]);
  const int16x8_t v_iscan_plus1 = vaddq_s16(v_iscan, vdupq_n_s16(1));
  const int16x8_t v_nz_iscan = vbslq_s16(v_mask, v_iscan_plus1, vdupq_n_s16(0));
  return vmaxq_s16(v_eobmax, v_nz_iscan);
}

static INLINE uint16_t get_max_eob(int16x8_t v_eobmax) {
#ifdef __aarch64__
  return (uint16_t)vmaxvq_s16(v_eobmax);
#else
  const int16x4_t v_eobmax_3210 =
      vmax_s16(vget_low_s16(v_eobmax), vget_high_s16(v_eobmax));
  const int64x1_t v_eobmax_xx32 =
      vshr_n_s64(vreinterpret_s64_s16(v_eobmax_3210), 32);
  const int16x4_t v_eobmax_tmp =
      vmax_s16(v_eobmax_3210, vreinterpret_s16_s64(v_eobmax_xx32));
  const int64x1_t v_eobmax_xxx3 =
      vshr_n_s64(vreinterpret_s64_s16(v_eobmax_tmp), 16);
  const int16x4_t v_eobmax_final =
      vmax_s16(v_eobmax_tmp, vreinterpret_s16_s64(v_eobmax_xxx3));
  return (uint16_t)vget_lane_s16(v_eobmax_final, 0);
#endif
}

void av1_highbd_quantize_fp_neon(
    const tran_low_t *coeff_ptr, intptr_t count, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan, int log_scale) {
  (void)scan;
  (void)zbin_ptr;
  (void)quant_shift_ptr;

  const int shift = 16 - log_scale;
  const int16x4_t v_quant = vld1_s16(quant_ptr);
  const int16x4_t v_dequant = vld1_s16(dequant_ptr);
  const int16x4_t v_zero = vdup_n_s16(0);
  const uint16x4_t v_round_select = vcgt_s16(vdup_n_s16(log_scale), v_zero);
  const int16x4_t v_round_no_scale = vld1_s16(round_ptr);
  const int16x4_t v_round_log_scale =
      vqrdmulh_n_s16(v_round_no_scale, (int16_t)(1 << (15 - log_scale)));
  const int16x4_t v_round =
      vbsl_s16(v_round_select, v_round_log_scale, v_round_no_scale);
  int32x4_t v_round_s32 = vaddl_s16(v_round, v_zero);
  int32x4_t v_quant_s32 = vaddl_s16(v_quant, v_zero);
  int32x4_t v_dequant_s32 = vaddl_s16(v_dequant, v_zero);
  uint16x4_t v_mask_lo, v_mask_hi;
  int16x8_t v_eobmax = vdupq_n_s16(-1);

  // DC and first 3 AC
  v_mask_lo = quantize_4(coeff_ptr, qcoeff_ptr, dqcoeff_ptr, v_quant_s32,
                         v_dequant_s32, v_round_s32, log_scale, shift);

  // overwrite the DC round with AC round
  v_round_s32 = vdupq_lane_s32(vget_low_s32(v_round_s32), 1);
  v_quant_s32 = vdupq_lane_s32(vget_low_s32(v_quant_s32), 1);
  v_dequant_s32 = vdupq_lane_s32(vget_low_s32(v_dequant_s32), 1);

  // 4 more DC
  v_mask_hi =
      quantize_4(coeff_ptr + 4, qcoeff_ptr + 4, dqcoeff_ptr + 4, v_quant_s32,
                 v_dequant_s32, v_round_s32, log_scale, shift);

  // Find the max lane eob for the first 8 coeffs.
  v_eobmax =
      get_max_lane_eob(iscan, v_eobmax, vcombine_u16(v_mask_lo, v_mask_hi));

  count -= 8;
  do {
    coeff_ptr += 8;
    qcoeff_ptr += 8;
    dqcoeff_ptr += 8;
    iscan += 8;
    v_mask_lo = quantize_4(coeff_ptr, qcoeff_ptr, dqcoeff_ptr, v_quant_s32,
                           v_dequant_s32, v_round_s32, log_scale, shift);
    v_mask_hi =
        quantize_4(coeff_ptr + 4, qcoeff_ptr + 4, dqcoeff_ptr + 4, v_quant_s32,
                   v_dequant_s32, v_round_s32, log_scale, shift);
    // Find the max lane eob for 8 coeffs.
    v_eobmax =
        get_max_lane_eob(iscan, v_eobmax, vcombine_u16(v_mask_lo, v_mask_hi));
    count -= 8;
  } while (count);

  *eob_ptr = get_max_eob(v_eobmax);
}
