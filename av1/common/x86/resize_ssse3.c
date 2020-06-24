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

#include "aom_dsp/x86/mem_sse2.h"
#include "aom_dsp/x86/transpose_sse2.h"
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

static INLINE __m128i convolve8_8_ssse3(const __m128i *const s,
                                        const __m128i *const f) {
  // multiply 2 adjacent elements with the filter and add the result
  const __m128i k_64 = _mm_set1_epi16(1 << 6);
  const __m128i x0 = _mm_maddubs_epi16(s[0], f[0]);
  const __m128i x1 = _mm_maddubs_epi16(s[1], f[1]);
  const __m128i x2 = _mm_maddubs_epi16(s[2], f[2]);
  const __m128i x3 = _mm_maddubs_epi16(s[3], f[3]);
  __m128i sum1, sum2;

  // sum the results together, saturating only on the final step
  // adding x0 with x2 and x1 with x3 is the only order that prevents
  // outranges for all filters
  sum1 = _mm_add_epi16(x0, x2);
  sum2 = _mm_add_epi16(x1, x3);
  // add the rounding offset early to avoid another saturated add
  sum1 = _mm_add_epi16(sum1, k_64);
  sum1 = _mm_adds_epi16(sum1, sum2);
  // shift by 7 bit each 16 bit
  sum1 = _mm_srai_epi16(sum1, 7);
  return sum1;
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

static INLINE void shuffle_filter_ssse3(const int16_t *const filter,
                                        __m128i *const f) {
  const __m128i f_values = _mm_load_si128((const __m128i *)filter);
  // pack and duplicate the filter values
  f[0] = _mm_shuffle_epi8(f_values, _mm_set1_epi16(0x0200u));
  f[1] = _mm_shuffle_epi8(f_values, _mm_set1_epi16(0x0604u));
  f[2] = _mm_shuffle_epi8(f_values, _mm_set1_epi16(0x0a08u));
  f[3] = _mm_shuffle_epi8(f_values, _mm_set1_epi16(0x0e0cu));
}

static void scale_plane_2_to_1_general(const uint8_t *src, const int src_stride,
                                       uint8_t *dst, const int dst_stride,
                                       const int w, const int h,
                                       const int16_t *const coef,
                                       uint8_t *const temp_buffer) {
  const int width_hor = (w + 3) & ~3;
  const int width_ver = (w + 7) & ~7;
  const int height_hor = (2 * h + SUBPEL_TAPS - 2 + 7) & ~7;
  const int height_ver = (h + 3) & ~3;
  int x, y = height_hor;
  uint8_t *t = temp_buffer;
  __m128i s[11], d[4];
  __m128i f[4];

  assert(w && h);

  shuffle_filter_ssse3(coef, f);
  src -= (SUBPEL_TAPS / 2 - 1) * src_stride + SUBPEL_TAPS / 2 + 1;

  // horizontal 4x8
  do {
    load_8bit_8x8(src + 2, src_stride, s);
    // 00 01 10 11 20 21 30 31  40 41 50 51 60 61 70 71
    // 02 03 12 13 22 23 32 33  42 43 52 53 62 63 72 73
    // 04 05 14 15 24 25 34 35  44 45 54 55 64 65 74 75
    // 06 07 16 17 26 27 36 37  46 47 56 57 66 67 76 77 (overlapped)
    transpose_16bit_4x8(s, s);
    x = width_hor;

    do {
      src += 8;
      load_8bit_8x8(src, src_stride, &s[3]);
      // 06 07 16 17 26 27 36 37  46 47 56 57 66 67 76 77
      // 08 09 18 19 28 29 38 39  48 49 58 59 68 69 78 79
      // 0A 0B 1A 1B 2A 2B 3A 3B  4A 4B 5A 5B 6A 6B 7A 7B
      // 0C 0D 1C 1D 2C 2D 3C 3D  4C 4D 5C 5D 6C 6D 7C 7D
      transpose_16bit_4x8(&s[3], &s[3]);

      d[0] = convolve8_8_ssse3(&s[0], f);  // 00 10 20 30 40 50 60 70
      d[1] = convolve8_8_ssse3(&s[1], f);  // 01 11 21 31 41 51 61 71
      d[2] = convolve8_8_ssse3(&s[2], f);  // 02 12 22 32 42 52 62 72
      d[3] = convolve8_8_ssse3(&s[3], f);  // 03 13 23 33 43 53 63 73

      // 00 10 20 30 40 50 60 70  02 12 22 32 42 52 62 72
      // 01 11 21 31 41 51 61 71  03 13 23 33 43 53 63 73
      d[0] = _mm_packus_epi16(d[0], d[2]);
      d[1] = _mm_packus_epi16(d[1], d[3]);
      // 00 10 01 11 20 30 21 31  40 50 41 51 60 70 61 71
      // 02 12 03 13 22 32 23 33  42 52 43 53 62 72 63 73
      d[2] = _mm_unpacklo_epi16(d[0], d[1]);
      d[3] = _mm_unpackhi_epi16(d[0], d[1]);
      // 00 10 01 11 02 12 03 13  20 30 21 31 22 32 23 33
      // 40 50 41 51 42 52 43 53  60 70 61 71 62 72 63 73
      d[0] = _mm_unpacklo_epi32(d[2], d[3]);
      d[1] = _mm_unpackhi_epi32(d[2], d[3]);
      store_8bit_8x4_from_16x2(d, t, 2 * width_hor);

      s[0] = s[4];
      s[1] = s[5];
      s[2] = s[6];

      t += 8;
      x -= 4;
    } while (x);
    src += 8 * src_stride - 2 * width_hor;
    t += 6 * width_hor;
    y -= 8;
  } while (y);

  // vertical 8x4
  x = width_ver;
  t = temp_buffer;
  do {
    // 00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17
    // 20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37
    // 40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57
    s[0] = _mm_loadu_si128((const __m128i *)(t + 0 * width_hor));
    s[1] = _mm_loadu_si128((const __m128i *)(t + 2 * width_hor));
    s[2] = _mm_loadu_si128((const __m128i *)(t + 4 * width_hor));
    t += 6 * width_hor;
    y = height_ver;

    do {
      // 60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77
      // 80 90 81 91 82 92 83 93  84 94 85 95 86 96 87 77
      // A0 B0 A1 B1 A2 B2 A3 B3  A4 B4 A5 B5 A6 B6 A7 77
      // C0 D0 C1 D1 C2 D2 C3 D3  C4 D4 C5 D5 C6 D6 C7 77
      loadu_8bit_16x4(t, 2 * width_hor, &s[3]);
      t += 8 * width_hor;

      d[0] = convolve8_8_ssse3(&s[0], f);  // 00 01 02 03 04 05 06 07
      d[1] = convolve8_8_ssse3(&s[1], f);  // 10 11 12 13 14 15 16 17
      d[2] = convolve8_8_ssse3(&s[2], f);  // 20 21 22 23 24 25 26 27
      d[3] = convolve8_8_ssse3(&s[3], f);  // 30 31 32 33 34 35 36 37

      // 00 01 02 03 04 05 06 07  10 11 12 13 14 15 16 17
      // 20 21 22 23 24 25 26 27  30 31 32 33 34 35 36 37
      d[0] = _mm_packus_epi16(d[0], d[1]);
      d[1] = _mm_packus_epi16(d[2], d[3]);
      store_8bit_8x4_from_16x2(d, dst, dst_stride);

      s[0] = s[4];
      s[1] = s[5];
      s[2] = s[6];

      dst += 4 * dst_stride;
      y -= 4;
    } while (y);
    t -= width_hor * (2 * height_ver + 6);
    t += 16;
    dst -= height_ver * dst_stride;
    dst += 8;
    x -= 8;
  } while (x);
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
        const int buffer_stride = (dst_w + 3) & ~3;
        const int buffer_height = (2 * dst_h + SUBPEL_TAPS - 2 + 7) & ~7;
        uint8_t *const temp_buffer =
            (uint8_t *)malloc(buffer_stride * buffer_height);
        if (temp_buffer) {
          const InterpKernel *interp_kernel =
              (const InterpKernel *)av1_interp_filter_params_list[filter]
                  .filter_ptr;
          scale_plane_2_to_1_general(src->buffers[i], src->strides[is_uv],
                                     dst->buffers[i], dst->strides[is_uv],
                                     dst_w, dst_h, interp_kernel[phase],
                                     temp_buffer);
          free(temp_buffer);
        }
      }
    } else {
      av1_resize_plane(src->buffers[i], src_h, src_w, src->strides[is_uv],
                       dst->buffers[i], dst_h, dst_w, dst->strides[is_uv]);
    }
  }
  aom_extend_frame_borders(dst, num_planes);
}
