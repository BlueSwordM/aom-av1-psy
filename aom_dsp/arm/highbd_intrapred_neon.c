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

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"

// -----------------------------------------------------------------------------
// PAETH

static INLINE void highbd_paeth_4or8_x_h_neon(uint16_t *dest, ptrdiff_t stride,
                                              const uint16_t *const top_row,
                                              const uint16_t *const left_column,
                                              int width, int height) {
  const uint16x8_t top_left = vdupq_n_u16(top_row[-1]);
  const uint16x8_t top_left_x2 = vdupq_n_u16(top_row[-1] + top_row[-1]);
  uint16x8_t top;
  if (width == 4) {
    top = vcombine_u16(vld1_u16(top_row), vdup_n_u16(0));
  } else {  // width == 8
    top = vld1q_u16(top_row);
  }

  for (int y = 0; y < height; ++y) {
    const uint16x8_t left = vdupq_n_u16(left_column[y]);

    const uint16x8_t left_dist = vabdq_u16(top, top_left);
    const uint16x8_t top_dist = vabdq_u16(left, top_left);
    const uint16x8_t top_left_dist =
        vabdq_u16(vaddq_u16(top, left), top_left_x2);

    const uint16x8_t left_le_top = vcleq_u16(left_dist, top_dist);
    const uint16x8_t left_le_top_left = vcleq_u16(left_dist, top_left_dist);
    const uint16x8_t top_le_top_left = vcleq_u16(top_dist, top_left_dist);

    // if (left_dist <= top_dist && left_dist <= top_left_dist)
    const uint16x8_t left_mask = vandq_u16(left_le_top, left_le_top_left);
    //   dest[x] = left_column[y];
    // Fill all the unused spaces with 'top'. They will be overwritten when
    // the positions for top_left are known.
    uint16x8_t result = vbslq_u16(left_mask, left, top);
    // else if (top_dist <= top_left_dist)
    //   dest[x] = top_row[x];
    // Add these values to the mask. They were already set.
    const uint16x8_t left_or_top_mask = vorrq_u16(left_mask, top_le_top_left);
    // else
    //   dest[x] = top_left;
    result = vbslq_u16(left_or_top_mask, result, top_left);

    if (width == 4) {
      vst1_u16(dest, vget_low_u16(result));
    } else {  // width == 8
      vst1q_u16(dest, result);
    }
    dest += stride;
  }
}

#define HIGHBD_PAETH_NXM(W, H)                                  \
  void aom_highbd_paeth_predictor_##W##x##H##_neon(             \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,   \
      const uint16_t *left, int bd) {                           \
    (void)bd;                                                   \
    highbd_paeth_4or8_x_h_neon(dst, stride, above, left, W, H); \
  }

HIGHBD_PAETH_NXM(4, 4)
HIGHBD_PAETH_NXM(4, 8)
HIGHBD_PAETH_NXM(8, 4)
HIGHBD_PAETH_NXM(8, 8)
HIGHBD_PAETH_NXM(8, 16)

#if !CONFIG_REALTIME_ONLY
HIGHBD_PAETH_NXM(4, 16)
HIGHBD_PAETH_NXM(8, 32)
#endif

// Select the closest values and collect them.
static INLINE uint16x8_t select_paeth(const uint16x8_t top,
                                      const uint16x8_t left,
                                      const uint16x8_t top_left,
                                      const uint16x8_t left_le_top,
                                      const uint16x8_t left_le_top_left,
                                      const uint16x8_t top_le_top_left) {
  // if (left_dist <= top_dist && left_dist <= top_left_dist)
  const uint16x8_t left_mask = vandq_u16(left_le_top, left_le_top_left);
  //   dest[x] = left_column[y];
  // Fill all the unused spaces with 'top'. They will be overwritten when
  // the positions for top_left are known.
  const uint16x8_t result = vbslq_u16(left_mask, left, top);
  // else if (top_dist <= top_left_dist)
  //   dest[x] = top_row[x];
  // Add these values to the mask. They were already set.
  const uint16x8_t left_or_top_mask = vorrq_u16(left_mask, top_le_top_left);
  // else
  //   dest[x] = top_left;
  return vbslq_u16(left_or_top_mask, result, top_left);
}

#define PAETH_PREDICTOR(num)                                                  \
  do {                                                                        \
    const uint16x8_t left_dist = vabdq_u16(top[num], top_left);               \
    const uint16x8_t top_left_dist =                                          \
        vabdq_u16(vaddq_u16(top[num], left), top_left_x2);                    \
    const uint16x8_t left_le_top = vcleq_u16(left_dist, top_dist);            \
    const uint16x8_t left_le_top_left = vcleq_u16(left_dist, top_left_dist);  \
    const uint16x8_t top_le_top_left = vcleq_u16(top_dist, top_left_dist);    \
    const uint16x8_t result =                                                 \
        select_paeth(top[num], left, top_left, left_le_top, left_le_top_left, \
                     top_le_top_left);                                        \
    vst1q_u16(dest + (num * 8), result);                                      \
  } while (0)

#define LOAD_TOP_ROW(num) vld1q_u16(top_row + (num * 8))

static INLINE void highbd_paeth16_plus_x_h_neon(
    uint16_t *dest, ptrdiff_t stride, const uint16_t *const top_row,
    const uint16_t *const left_column, int width, int height) {
  const uint16x8_t top_left = vdupq_n_u16(top_row[-1]);
  const uint16x8_t top_left_x2 = vdupq_n_u16(top_row[-1] + top_row[-1]);
  uint16x8_t top[8];
  top[0] = LOAD_TOP_ROW(0);
  top[1] = LOAD_TOP_ROW(1);
  if (width > 16) {
    top[2] = LOAD_TOP_ROW(2);
    top[3] = LOAD_TOP_ROW(3);
    if (width == 64) {
      top[4] = LOAD_TOP_ROW(4);
      top[5] = LOAD_TOP_ROW(5);
      top[6] = LOAD_TOP_ROW(6);
      top[7] = LOAD_TOP_ROW(7);
    }
  }

  for (int y = 0; y < height; ++y) {
    const uint16x8_t left = vdupq_n_u16(left_column[y]);
    const uint16x8_t top_dist = vabdq_u16(left, top_left);
    PAETH_PREDICTOR(0);
    PAETH_PREDICTOR(1);
    if (width > 16) {
      PAETH_PREDICTOR(2);
      PAETH_PREDICTOR(3);
      if (width == 64) {
        PAETH_PREDICTOR(4);
        PAETH_PREDICTOR(5);
        PAETH_PREDICTOR(6);
        PAETH_PREDICTOR(7);
      }
    }
    dest += stride;
  }
}

#define HIGHBD_PAETH_NXM_WIDE(W, H)                               \
  void aom_highbd_paeth_predictor_##W##x##H##_neon(               \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,     \
      const uint16_t *left, int bd) {                             \
    (void)bd;                                                     \
    highbd_paeth16_plus_x_h_neon(dst, stride, above, left, W, H); \
  }

HIGHBD_PAETH_NXM_WIDE(16, 8)
HIGHBD_PAETH_NXM_WIDE(16, 16)
HIGHBD_PAETH_NXM_WIDE(16, 32)
HIGHBD_PAETH_NXM_WIDE(32, 16)
HIGHBD_PAETH_NXM_WIDE(32, 32)
HIGHBD_PAETH_NXM_WIDE(32, 64)
HIGHBD_PAETH_NXM_WIDE(64, 32)
HIGHBD_PAETH_NXM_WIDE(64, 64)

#if !CONFIG_REALTIME_ONLY
HIGHBD_PAETH_NXM_WIDE(16, 4)
HIGHBD_PAETH_NXM_WIDE(16, 64)
HIGHBD_PAETH_NXM_WIDE(32, 8)
HIGHBD_PAETH_NXM_WIDE(64, 16)
#endif
