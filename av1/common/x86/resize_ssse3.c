/*
 *
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <tmmintrin.h>  // SSSE3

#include "av1/common/resize.h"
#include "config/av1_rtcd.h"
#include "config/aom_scale_rtcd.h"

static INLINE __m128i scale_plane_2_to_1_phase_0_kernel(
    const uint8_t *const src, const __m128i *const mask) {
  const __m128i a = _mm_loadu_si128((const __m128i *)(&src[0]));
  const __m128i b = _mm_loadu_si128((const __m128i *)(&src[16]));
  const __m128i a_and = _mm_and_si128(a, *mask);
  const __m128i b_and = _mm_and_si128(b, *mask);
  return _mm_packus_epi16(a_and, b_and);
}

static void scale_plane_2_to_1_phase_0(const uint8_t *src,
                                       const ptrdiff_t src_stride, uint8_t *dst,
                                       const ptrdiff_t dst_stride,
                                       const int dst_w, const int dst_h) {
  const int max_width = (dst_w + 15) & ~15;
  const __m128i mask = _mm_set1_epi16(0x00FF);
  int y = dst_h;

  do {
    int x = max_width;
    do {
      const __m128i d = scale_plane_2_to_1_phase_0_kernel(src, &mask);
      _mm_storeu_si128((__m128i *)dst, d);
      src += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

static INLINE __m128i scale_plane_bilinear_kernel(const __m128i *const s,
                                                  const __m128i c0c1) {
  const __m128i k_64 = _mm_set1_epi16(1 << 6);
  const __m128i t0 = _mm_maddubs_epi16(s[0], c0c1);
  const __m128i t1 = _mm_maddubs_epi16(s[1], c0c1);
  // round and shift by 7 bit each 16 bit
  const __m128i t2 = _mm_adds_epi16(t0, k_64);
  const __m128i t3 = _mm_adds_epi16(t1, k_64);
  const __m128i t4 = _mm_srai_epi16(t2, 7);
  const __m128i t5 = _mm_srai_epi16(t3, 7);
  return _mm_packus_epi16(t4, t5);
}

static void scale_plane_2_to_1_bilinear(const uint8_t *src,
                                        const ptrdiff_t src_stride,
                                        uint8_t *dst,
                                        const ptrdiff_t dst_stride,
                                        const int dst_w, const int dst_h,
                                        const __m128i c0c1) {
  const int max_width = (dst_w + 15) & ~15;
  int y = dst_h;

  do {
    int x = max_width;
    do {
      __m128i s[2], d[2];

      // Horizontal
      // Even rows
      s[0] = _mm_loadu_si128((const __m128i *)(src + 0));
      s[1] = _mm_loadu_si128((const __m128i *)(src + 16));
      d[0] = scale_plane_bilinear_kernel(s, c0c1);

      // odd rows
      s[0] = _mm_loadu_si128((const __m128i *)(src + src_stride + 0));
      s[1] = _mm_loadu_si128((const __m128i *)(src + src_stride + 16));
      d[1] = scale_plane_bilinear_kernel(s, c0c1);

      // Vertical
      s[0] = _mm_unpacklo_epi8(d[0], d[1]);
      s[1] = _mm_unpackhi_epi8(d[0], d[1]);
      d[0] = scale_plane_bilinear_kernel(s, c0c1);

      _mm_storeu_si128((__m128i *)dst, d[0]);
      src += 32;
      dst += 16;
      x -= 16;
    } while (x);
    src += 2 * (src_stride - max_width);
    dst += dst_stride - max_width;
  } while (--y);
}

void av1_resize_and_extend_frame_ssse3(const YV12_BUFFER_CONFIG *src,
                                       YV12_BUFFER_CONFIG *dst,
                                       const InterpFilter filter,
                                       const int phase, const int num_planes) {
  // We use AOMMIN(num_planes, MAX_MB_PLANE) instead of num_planes to quiet
  // the static analysis warnings.
  for (int i = 0; i < AOMMIN(num_planes, MAX_MB_PLANE); ++i) {
    const int is_uv = i > 0;
    const int src_w = src->crop_widths[is_uv];
    const int src_h = src->crop_heights[is_uv];
    const int dst_w = dst->crop_widths[is_uv];
    const int dst_h = dst->crop_heights[is_uv];

    if (2 * dst_w == src_w && 2 * dst_h == src_h) {
      if (phase == 0) {
        scale_plane_2_to_1_phase_0(src->buffers[i], src->strides[is_uv],
                                   dst->buffers[i], dst->strides[is_uv], dst_w,
                                   dst_h);
      } else if (filter == BILINEAR) {
        const int16_t c0 = av1_bilinear_filters[phase][3];
        const int16_t c1 = av1_bilinear_filters[phase][4];
        const __m128i c0c1 = _mm_set1_epi16(c0 | (c1 << 8));  // c0 and c1 >= 0
        scale_plane_2_to_1_bilinear(src->buffers[i], src->strides[is_uv],
                                    dst->buffers[i], dst->strides[is_uv], dst_w,
                                    dst_h, c0c1);
      } else {
        av1_resize_plane(src->buffers[i], src_h, src_w, src->strides[is_uv],
                         dst->buffers[i], dst_h, dst_w, dst->strides[is_uv]);
      }
    } else {
      av1_resize_plane(src->buffers[i], src_h, src_w, src->strides[is_uv],
                       dst->buffers[i], dst_h, dst_w, dst->strides[is_uv]);
    }
  }
  aom_extend_frame_borders(dst, num_planes);
}
