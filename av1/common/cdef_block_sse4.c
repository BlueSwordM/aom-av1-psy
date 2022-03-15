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

#include "aom_dsp/aom_simd.h"
#define SIMD_FUNC(name) name##_sse4_1
#include "av1/common/cdef_block_simd.h"

void cdef_find_dir_dual_sse4_1(const uint16_t *img1, const uint16_t *img2,
                               int stride, int32_t *var_out_1st,
                               int32_t *var_out_2nd, int coeff_shift,
                               int *out_dir_1st_8x8, int *out_dir_2nd_8x8) {
  // Process first 8x8.
  *out_dir_1st_8x8 = cdef_find_dir(img1, stride, var_out_1st, coeff_shift);

  // Process second 8x8.
  *out_dir_2nd_8x8 = cdef_find_dir(img2, stride, var_out_2nd, coeff_shift);
}
