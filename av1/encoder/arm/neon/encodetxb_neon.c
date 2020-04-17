/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>

#include <math.h>

#include "av1/common/txb_common.h"

void av1_txb_init_levels_neon(const tran_low_t *const coeff, const int width,
                              const int height, uint8_t *const levels) {
  const int stride = width + TX_PAD_HOR;
  memset(levels - TX_PAD_TOP * stride, 0,
         sizeof(*levels) * TX_PAD_TOP * stride);
  memset(levels + stride * height, 0,
         sizeof(*levels) * (TX_PAD_BOTTOM * stride + TX_PAD_END));

  const int32x4_t zeros = vdupq_n_s32(0);
  int i = 0;
  uint8_t *ls = levels;
  const tran_low_t *cf = coeff;
  if (width == 4) {
    do {
      const int32x4_t coeffA = vld1q_s32(cf);
      const int32x4_t coeffB = vld1q_s32(cf + width);
      const int16x8_t coeffAB =
          vcombine_s16(vqmovn_s32(coeffA), vqmovn_s32(coeffB));
      const int16x8_t absAB = vqabsq_s16(coeffAB);
      const int8x8_t absABs = vqmovn_s16(absAB);
#if defined(__aarch64__)
      const int8x16_t absAB8 =
          vcombine_s8(absABs, vreinterpret_s8_s32(vget_low_s32(zeros)));
      const uint8x16_t lsAB =
          vreinterpretq_u8_s32(vzip1q_s32(vreinterpretq_s32_s8(absAB8), zeros));
#else
      const int32x2x2_t absAB8 =
          vzip_s32(vreinterpret_s32_s8(absABs), vget_low_s32(zeros));
      const uint8x16_t lsAB =
          vreinterpretq_u8_s32(vcombine_s32(absAB8.val[0], absAB8.val[1]));
#endif
      vst1q_u8(ls, lsAB);
      ls += (stride << 1);
      cf += (width << 1);
      i += 2;
    } while (i < height);
  } else if (width == 8) {
    do {
      const int32x4_t coeffA = vld1q_s32(cf);
      const int32x4_t coeffB = vld1q_s32(cf + 4);
      const int16x8_t coeffAB =
          vcombine_s16(vqmovn_s32(coeffA), vqmovn_s32(coeffB));
      const int16x8_t absAB = vqabsq_s16(coeffAB);
      const uint8x16_t absAB8 = vreinterpretq_u8_s8(vcombine_s8(
          vqmovn_s16(absAB), vreinterpret_s8_s32(vget_low_s32(zeros))));
      vst1q_u8(ls, absAB8);
      ls += stride;
      cf += width;
      i += 1;
    } while (i < height);
  } else {
    do {
      int j = 0;
      do {
        const int32x4_t coeffA = vld1q_s32(cf);
        const int32x4_t coeffB = vld1q_s32(cf + 4);
        const int32x4_t coeffC = vld1q_s32(cf + 8);
        const int32x4_t coeffD = vld1q_s32(cf + 12);
        const int16x8_t coeffAB =
            vcombine_s16(vqmovn_s32(coeffA), vqmovn_s32(coeffB));
        const int16x8_t coeffCD =
            vcombine_s16(vqmovn_s32(coeffC), vqmovn_s32(coeffD));
        const int16x8_t absAB = vqabsq_s16(coeffAB);
        const int16x8_t absCD = vqabsq_s16(coeffCD);
        const uint8x16_t absABCD = vreinterpretq_u8_s8(
            vcombine_s8(vqmovn_s16(absAB), vqmovn_s16(absCD)));
        vst1q_u8((ls + j), absABCD);
        j += 16;
        cf += 16;
      } while (j < width);
      *(int32_t *)(ls + width) = 0;
      ls += stride;
      i += 1;
    } while (i < height);
  }
}
