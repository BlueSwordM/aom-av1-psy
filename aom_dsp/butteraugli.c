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

#include <assert.h>
#include <jxl/butteraugli.h>

#include "aom_dsp/butteraugli.h"
#include "aom_mem/aom_mem.h"
#include "third_party/libyuv/include/libyuv/convert_argb.h"

int aom_calc_butteraugli(const YV12_BUFFER_CONFIG *source,
                         const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                         float *dist_map) {
  (void)bit_depth;
  assert(bit_depth == 8);
  assert(source->y_width == source->uv_width * 2);
  const int width = source->y_width;
  const int height = source->y_height;

  size_t buffer_size = width * height * 3;
  uint8_t *src_rgb = (uint8_t *)aom_malloc(buffer_size);
  uint8_t *distorted_rgb = (uint8_t *)aom_malloc(buffer_size);
  if (!src_rgb || !distorted_rgb) {
    aom_free(src_rgb);
    aom_free(distorted_rgb);
    return 0;
  }

  I420ToRGB24Matrix(source->y_buffer, source->y_stride, source->u_buffer,
                    source->uv_stride, source->v_buffer, source->uv_stride,
                    src_rgb, width * 3, &kYuvH709Constants, width, height);
  I420ToRGB24Matrix(distorted->y_buffer, distorted->y_stride,
                    distorted->u_buffer, distorted->uv_stride,
                    distorted->v_buffer, distorted->uv_stride, distorted_rgb,
                    width * 3, &kYuvH709Constants, width, height);

  JxlPixelFormat pixel_format = { 3, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0 };
  JxlButteraugliApi *api = JxlButteraugliApiCreate(NULL);
  JxlButteraugliApiSetHFAsymmetry(api, 0.8f);

  JxlButteraugliResult *result = JxlButteraugliCompute(
      api, width, height, &pixel_format, src_rgb, buffer_size, &pixel_format,
      distorted_rgb, buffer_size);

  const float *distmap = NULL;
  uint32_t row_stride;
  JxlButteraugliResultGetDistmap(result, &distmap, &row_stride);
  if (distmap == NULL) {
    JxlButteraugliApiDestroy(api);
    JxlButteraugliResultDestroy(result);
    aom_free(src_rgb);
    aom_free(distorted_rgb);
    return 0;
  }

  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      dist_map[j * width + i] = distmap[j * row_stride + i];
    }
  }

  JxlButteraugliApiDestroy(api);
  JxlButteraugliResultDestroy(result);
  aom_free(src_rgb);
  aom_free(distorted_rgb);
  return 1;
}
