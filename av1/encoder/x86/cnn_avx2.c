/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <immintrin.h>
#include <math.h>

#include "aom_dsp/aom_dsp_common.h"
#include "av1/common/av1_common_int.h"
#include "av1/encoder/cnn.h"

// This mask rearranges source pixels in the order shown below.
// shuffle_src_layer0[0][8]: applied on source pixels 0 to 7.
// shuffle_src_layer0[1][8]: applied on source pixels 7 to 14.
// This shuffling is needed to process 3 5x5 blocks which need
// source pixels in the following order.
// 1st 5x5 block: source pixels needed are 0 to 4,
// 2nd 5x5 block: source pixels needed are 4 to 8,
// 3rd 5x5 block: source pixels needed are 8 to 12.
// Source pixels are loaded like mentioned below.
// load_src0 : 0, 1, 2, 3, 4, 5, 6, 7
// load_src1 : 7, 8, 9, 10, 11, 12, 13, 14
// After applying masks, source bytes will be in the order:
// load_src0 : 0, 1, 2, 3, 4, 4, 5, 6
//             consists 5 pixels needed for 1st 5x5 block and
//             first 3 pixels needed for 2nd 5x5 block.
// load_src1 : 7, 8, 8, 9, 10, 11, 12, x
//             consists last 2 pixels needed for 2nd 5x5 block and
//             5 pixels needed for 3rd 5x5 block.
DECLARE_ALIGNED(32, static const uint32_t,
                shuffle_src_layer0[2][8]) = { { 0, 1, 2, 3, 4, 4, 5, 6 },
                                              { 0, 1, 1, 2, 3, 4, 5, 0 } };

// This mask rearrange the weights to match shuffled source pixels order.
DECLARE_ALIGNED(32, static const uint32_t,
                shuffle_weight_layer0[2][8]) = { { 0, 1, 2, 3, 4, 0, 1, 2 },
                                                 { 3, 4, 0, 1, 2, 3, 4, 0 } };

// Load weights needed for layer 0 (for 5x5 block processing),
// and fill the registers appropriately to match source pixel mapping.
static INLINE void prepare_weights_for_5x5_concolve(
    const CNN_LAYER_CONFIG *const layer_config, int off, float weight[5][8],
    const int cstep, __m256 *shuffle_weight, const __m256i weight_mask_0,
    const __m256i weight_mask_1) {
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      weight[row][col] = layer_config->weights[off];
      off += cstep;
    }
  }
  shuffle_weight[0] = _mm256_loadu_ps(weight[0]);
  shuffle_weight[1] = _mm256_loadu_ps(weight[1]);
  shuffle_weight[2] = _mm256_loadu_ps(weight[2]);
  shuffle_weight[3] = _mm256_loadu_ps(weight[3]);
  shuffle_weight[4] = _mm256_loadu_ps(weight[4]);

  shuffle_weight[0] =
      _mm256_permutevar8x32_ps(shuffle_weight[0], weight_mask_0);
  shuffle_weight[1] =
      _mm256_permutevar8x32_ps(shuffle_weight[1], weight_mask_0);
  shuffle_weight[2] =
      _mm256_permutevar8x32_ps(shuffle_weight[2], weight_mask_0);
  shuffle_weight[3] =
      _mm256_permutevar8x32_ps(shuffle_weight[3], weight_mask_0);
  shuffle_weight[4] =
      _mm256_permutevar8x32_ps(shuffle_weight[4], weight_mask_0);
  shuffle_weight[5] =
      _mm256_permutevar8x32_ps(shuffle_weight[0], weight_mask_1);
  shuffle_weight[6] =
      _mm256_permutevar8x32_ps(shuffle_weight[1], weight_mask_1);
  shuffle_weight[7] =
      _mm256_permutevar8x32_ps(shuffle_weight[2], weight_mask_1);
  shuffle_weight[8] =
      _mm256_permutevar8x32_ps(shuffle_weight[3], weight_mask_1);
  shuffle_weight[9] =
      _mm256_permutevar8x32_ps(shuffle_weight[4], weight_mask_1);
}

// For each row, loads source pixels 0 to 7(load_src_0), 7 to 14(load_src_1) and
// arranges them appropriately to process 3 blocks.
#define PERFORM_CONVOLVE_FOR_3_5X5_BLOCKS()                            \
  {                                                                    \
    for (int row = 0; row < 5; row++) {                                \
      load_src_0 = _mm256_loadu_ps(input_ptr);                         \
      load_src_1 = _mm256_loadu_ps(input_ptr + 7);                     \
      load_src_0 = _mm256_permutevar8x32_ps(load_src_0, block0_1);     \
      load_src_1 = _mm256_permutevar8x32_ps(load_src_1, block1_2);     \
      load_src_0 = _mm256_mul_ps(load_src_0, shuffle_weight[0 + row]); \
      load_src_1 = _mm256_mul_ps(load_src_1, shuffle_weight[5 + row]); \
      accum_src_0 = _mm256_add_ps(load_src_0, accum_src_0);            \
      accum_src_1 = _mm256_add_ps(load_src_1, accum_src_1);            \
      input_ptr += in_stride;                                          \
    }                                                                  \
  }

// Do convolution of one 5x5 block.
#define PERFORM_CONVOLVE_FOR_1_5X5_BLOCK(w, accum0, in_stride)           \
  {                                                                      \
    __m128 load_src[5];                                                  \
    load_src[0] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[0][4];                      \
    input_ptr += in_stride;                                              \
    load_src[1] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[1][4];                      \
    input_ptr += in_stride;                                              \
    load_src[2] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[2][4];                      \
    input_ptr += in_stride;                                              \
    load_src[3] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[3][4];                      \
    input_ptr += in_stride;                                              \
    load_src[4] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[4][4];                      \
                                                                         \
    load_src[0] = _mm_mul_ps(load_src[0], _mm256_castps256_ps128(w[0])); \
    load_src[1] = _mm_mul_ps(load_src[1], _mm256_castps256_ps128(w[1])); \
    load_src[2] = _mm_mul_ps(load_src[2], _mm256_castps256_ps128(w[2])); \
    load_src[3] = _mm_mul_ps(load_src[3], _mm256_castps256_ps128(w[3])); \
    load_src[4] = _mm_mul_ps(load_src[4], _mm256_castps256_ps128(w[4])); \
                                                                         \
    accum0 = _mm_add_ps(load_src[0], accum0);                            \
    load_src[1] = _mm_add_ps(load_src[1], load_src[2]);                  \
    load_src[3] = _mm_add_ps(load_src[3], load_src[4]);                  \
    load_src[1] = _mm_add_ps(load_src[1], load_src[3]);                  \
    accum0 = _mm_add_ps(accum0, load_src[1]);                            \
  }

// AVX2 variant of cnn_no_maxpool_padding_valid(), when filter_width and
// filter_height are equal to 5.
// CNN convolve parsing is based on av1_intra_mode_cnn_partition_cnn_config.
// Based on the configuration set for each layer, the current encoder
// always chooses the case of no_maxpool_padding_valid.
// And also for layer 0 convolution happens at 5x5 level as the
// filter_width and filter_height are set as 5.
static void cnn_convolve_no_maxpool_padding_valid_5x5_avx2(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int channel_step) {
  assert(layer_config->filter_width == 5 && layer_config->filter_height == 5);
  assert(layer_config->skip_width == 4 && layer_config->skip_height == 4);

  // Load shuffle buffers needed for source.
  const __m256i block0_1 =
      _mm256_load_si256((const __m256i *)shuffle_src_layer0[0]);
  const __m256i block1_2 =
      _mm256_load_si256((const __m256i *)shuffle_src_layer0[1]);

  // Load shuffle buffers needed for weight.
  const __m256i weight_mask_0 =
      _mm256_load_si256((const __m256i *)shuffle_weight_layer0[0]);
  const __m256i weight_mask_1 =
      _mm256_load_si256((const __m256i *)shuffle_weight_layer0[1]);

  // Width needs to be moved to go to next iteration of processing 3 5x5 blocks.
  const int skip_width_for_next_iter = layer_config->skip_width * 3;

  // Minimum width required to process 3 5x5 blocks at a time.
  // min width (for processing 3 5x5 block) = 2*skip_width + filter_width
  // Here, skip_width specifies how much width we should move while processing
  // next block convolution and filter_width specifies for how many pixels
  // filter needs to be applied.
  const int min_width_for_3_5x5_blocks =
      (layer_config->skip_width * 2) + layer_config->filter_width;
  for (int i = start_idx; i < layer_config->out_channels; i += channel_step) {
    const float out_ch_bias = layer_config->bias[i];
    for (int k = 0; k < layer_config->in_channels; ++k) {
      __m256 shuffle_weight[10];

      // Weights needed are 5x5, for SIMD purpose made this array as 5x8.
      float weight[5][8] = { { 0 } };
      int off = k * layer_config->out_channels + i;

      // In layer 0, the convolution process happens at 5x5.
      // The weights needed for 5x5 block are same across the in-channels,
      // which is why the load of weights happens once for each in-channel.
      prepare_weights_for_5x5_concolve(layer_config, off, weight, cstep,
                                       shuffle_weight, weight_mask_0,
                                       weight_mask_1);

      for (int h = 0, u = 0; h < in_height - layer_config->filter_height + 1;
           h += layer_config->skip_height, ++u) {
        const int out_h = u * out_stride;
        int v = 0;
        int w = 0;
        int rem_width = in_width;
        // Processing 3 5x5 blocks at a time, if sufficient width is present.
        while (rem_width >= min_width_for_3_5x5_blocks) {
          __m256 load_src_0, load_src_1;
          __m256 accum_src_0 = _mm256_setzero_ps();
          __m256 accum_src_1 = _mm256_setzero_ps();
          const float *input_ptr = &input[k][h * in_stride + w];
          PERFORM_CONVOLVE_FOR_3_5X5_BLOCKS();

          // Accumulate across column.
          __m256 accum = _mm256_hadd_ps(accum_src_0, accum_src_1);
          __m128 tmp_reg_0 = _mm256_extractf128_ps(accum_src_0, 1);
          __m128 tmp_reg_1 = _mm256_extractf128_ps(accum_src_1, 1);

          __m128 accum_l = _mm256_castps256_ps128(accum);
          __m128 accum_h = _mm256_extractf128_ps(accum, 1);

          __m128 tmp_reg_2 = _mm_add_ps(accum_l, tmp_reg_0);
          __m128 tmp_reg_3 = _mm_add_ps(tmp_reg_0, accum_h);
          __m128 tmp_reg_4 = _mm_add_ps(tmp_reg_1, accum_h);

          // 1st 5x5 block output.
          output[i][out_h + v] =
              out_ch_bias + _mm_cvtss_f32(tmp_reg_2) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 1));

          // 2nd 5x5 block output.
          output[i][out_h + v + 1] =
              out_ch_bias +
              _mm_cvtss_f32(_mm_shuffle_ps(tmp_reg_3, tmp_reg_3, 1)) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 2));

          // 3rd 5x5 block output.
          output[i][out_h + v + 2] =
              out_ch_bias +
              _mm_cvtss_f32(_mm_shuffle_ps(tmp_reg_4, tmp_reg_4, 2)) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 3));

          v += 3;
          w += skip_width_for_next_iter;
          rem_width -= skip_width_for_next_iter;
        }

        // Process remaining blocks as single 5x5 block at a time.
        while (rem_width >= layer_config->filter_width) {
          float last_column_sum = 0;
          __m128 accum = _mm_setzero_ps();
          const float *input_ptr = &input[k][h * in_stride + w];
          PERFORM_CONVOLVE_FOR_1_5X5_BLOCK(shuffle_weight, accum, in_stride);

          // Accumulate across column.
          accum = _mm_hadd_ps(accum, accum);
          output[i][out_h + v] = out_ch_bias + last_column_sum +
                                 _mm_cvtss_f32(accum) +
                                 _mm_cvtss_f32(_mm_shuffle_ps(accum, accum, 1));

          v += 1;
          w += layer_config->skip_width;
          rem_width -= layer_config->skip_width;
        }
      }
    }
  }
}

// AVX2 variant of av1_cnn_convolve_no_maxpool_padding_valid_c().
void av1_cnn_convolve_no_maxpool_padding_valid_avx2(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *layer_config, float **output, int out_stride,
    int start_idx, int cstep, int channel_step) {
  if (layer_config->filter_width == 5 && layer_config->filter_height == 5) {
    cnn_convolve_no_maxpool_padding_valid_5x5_avx2(
        input, in_width, in_height, in_stride, layer_config, output, out_stride,
        start_idx, cstep, channel_step);
  } else {
    av1_cnn_convolve_no_maxpool_padding_valid_c(
        input, in_width, in_height, in_stride, layer_config, output, out_stride,
        start_idx, cstep, channel_step);
  }
}
