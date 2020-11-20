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

#ifndef AOM_AV1_ENCODER_CNN_INTERNAL_H_
#define AOM_AV1_ENCODER_CNN_INTERNAL_H_

#ifdef __cplusplus
extern "C" {
#endif

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_SAME_ZERO.
void av1_cnn_convolve_maxpool_padding_zero(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep, const int filter_width_half, const int filter_height_half);

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_SAME_REPLICATE.
void av1_cnn_convolve_maxpool_padding_replicate(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep, const int filter_width_half, const int filter_height_half);

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_VALID.
void av1_cnn_convolve_maxpool_padding_valid(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep);

// CNNConvolve specific to maxpool set as 0 with filter_height and filter_width
// equal to 1.
void av1_cnn_convolve_element_wise(const float **input, int in_width,
                                   int in_height, int in_stride,
                                   const CNN_LAYER_CONFIG *const layer_config,
                                   float **output, int out_stride,
                                   int start_idx, int step);

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_SAME_ZERO.
void av1_cnn_convolve_no_maxpool_padding_zero(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int filter_width_half,
    const int filter_height_half, const int ii_shift, const int jj_shift,
    const int channel_step);

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_SAME_REPLICATE.
void av1_cnn_convolve_no_maxpool_padding_replicate(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int ii_shift, const int jj_shift,
    const int channel_step);

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_VALID.
void av1_cnn_convolve_no_maxpool_padding_valid(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int channel_step);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_CNN_INTERNAL_H_
