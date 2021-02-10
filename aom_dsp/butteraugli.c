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

// TODO(sdeng): update the jxl api.
#include <assert.h>
#include <jxl/encode.h>

#include "aom_dsp/butteraugli.h"

void aom_calc_butteraugli(const YV12_BUFFER_CONFIG *source,
                          const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                          float *dist_map) {
  assert(bit_depth == 8);
  assert(source->y_width == source->uv_width * 2);
  uint8_t *src_y = source->y_buffer;
  uint8_t *src_u = source->u_buffer;
  uint8_t *src_v = source->v_buffer;
  uint8_t *distorted_y = distorted->y_buffer;
  uint8_t *distorted_u = distorted->u_buffer;
  uint8_t *distorted_v = distorted->v_buffer;
  const int width = source->y_width;
  const int height = source->y_height;
  double butteraugli_diffvalue;
  JxlCalcButteraugliYuv420(width, height, src_y, source->y_stride, src_u, src_v,
                           source->uv_stride, distorted_y, distorted->y_stride,
                           distorted_u, distorted_v, distorted->uv_stride,
                           dist_map, &butteraugli_diffvalue);
  (void)bit_depth;
}
